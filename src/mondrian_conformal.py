import numpy as np

def compute_ece(y_true, y_prob, n_bins=10):
    y_true = np.asarray(y_true, dtype=float)
    y_prob = np.asarray(y_prob, dtype=float)
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece, n = 0.0, len(y_true)
    if n == 0:
        return 0.0
    for i in range(n_bins):
        mask = (y_prob >= bins[i]) & (y_prob < bins[i + 1] if i < n_bins - 1 else y_prob <= bins[i + 1])
        if np.any(mask):
            ece += (np.sum(mask) / n) * abs(np.mean(y_true[mask]) - np.mean(y_prob[mask]))
    return float(ece)

class AncestryMondrianConformal:
    """Non-empty Ancestry-Conditional Mondrian Split Conformal Classifier (v1.0.0)."""
    def __init__(self, alpha=0.10):
        self.alpha = alpha
        self.q_by_group = {}
        self.q_global = 1.0

    def _nonconformity(self, probs_pos, y_true):
        p_true = np.where(y_true == 1, probs_pos, 1.0 - probs_pos)
        return 1.0 - p_true

    def _finite_quantile(self, scores):
        n = len(scores)
        if n == 0:
            return 1.0
        level = min(1.0, np.ceil((n + 1) * (1.0 - self.alpha)) / n)
        return float(np.quantile(scores, level, method="higher"))

    def fit(self, probs_pos, y_true, groups):
        probs_pos = np.asarray(probs_pos, dtype=float)
        y_true = np.asarray(y_true, dtype=int)
        groups = np.asarray(groups)
        scores = self._nonconformity(probs_pos, y_true)
        self.q_global = self._finite_quantile(scores)
        self.q_by_group = {str(g): self._finite_quantile(scores[groups == g]) for g in np.unique(groups)}
        return self

    def predict_sets(self, probs_pos, groups):
        probs_pos = np.asarray(probs_pos, dtype=float)
        groups = np.asarray(groups)
        pred_sets = []
        for p, g in zip(probs_pos, groups):
            q = self.q_by_group.get(str(g), self.q_global)
            s = []
            if p <= q:          # class 0 nonconformity = 1 - (1 - p) = p
                s.append(0)
            if (1.0 - p) <= q:  # class 1 nonconformity = 1 - p
                s.append(1)
            if not s:           # v1.0.0 non-empty floor: always retain argmax class
                s.append(1 if p >= 0.5 else 0)
            pred_sets.append(s)
        return pred_sets

    def evaluate_by_group(self, probs_pos, y_true, groups):
        probs_pos = np.asarray(probs_pos, dtype=float)
        y_true = np.asarray(y_true, dtype=int)
        groups = np.asarray(groups)
        pred_sets = self.predict_sets(probs_pos, groups)
        report = {}
        for g in np.unique(groups):
            idx = np.where(groups == g)[0]
            covered = [int(y_true[i] in pred_sets[i]) for i in idx]
            sizes = [len(pred_sets[i]) for i in idx]
            report[str(g)] = {
                "n": int(len(idx)),
                "coverage": round(float(np.mean(covered)), 4),
                "avg_set_size": round(float(np.mean(sizes)), 4),
                "singleton_rate": round(float(np.mean([sz == 1 for sz in sizes])), 4),
                "empty_sets": int(sum(1 for sz in sizes if sz == 0)),
                "raw_ece": round(compute_ece(y_true[idx], probs_pos[idx]), 4),
                "q_threshold": round(self.q_by_group.get(str(g), self.q_global), 4),
            }
        return report
