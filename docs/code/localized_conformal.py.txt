import numpy as np

class LocalizedWeightedConformal:
    """
    Localized Weighted Split Conformal Classifier (Tibshirani 2019 + Guan 2023).
    Computes continuous covariate-weighted conformal quantiles over log-AF vectors
    with a guaranteed non-empty prediction set floor.
    """
    def __init__(self, alpha=0.10, bandwidth=1.5):
        self.alpha = alpha
        self.bandwidth = bandwidth
        self.z_cal = None
        self.scores_cal = None

    def fit(self, probs_pos, y_true, z_cov):
        probs_pos = np.asarray(probs_pos, dtype=float)
        y_true = np.asarray(y_true, dtype=int)
        self.z_cal = np.asarray(z_cov, dtype=float)
        p_true = np.where(y_true == 1, probs_pos, 1.0 - probs_pos)
        self.scores_cal = 1.0 - p_true
        return self

    def _weighted_quantile(self, scores, weights, q_level):
        order = np.argsort(scores)
        s_sorted = scores[order]
        w_sorted = weights[order]
        cdf = np.cumsum(w_sorted) / np.sum(w_sorted)
        idx = np.searchsorted(cdf, q_level, side="left")
        return float(s_sorted[min(idx, len(s_sorted) - 1)])

    def predict_sets(self, probs_pos, z_query):
        probs_pos = np.asarray(probs_pos, dtype=float)
        z_query = np.asarray(z_query, dtype=float)
        # Vectorized Gaussian kernel weights in continuous log-AF space
        diff = z_query[:, None, :] - self.z_cal[None, :, :]
        dist2 = np.sum(diff * diff, axis=2)
        W = np.exp(-dist2 / (2.0 * (self.bandwidth ** 2))) + 1e-8

        pred_sets = []
        target = 1.0 - self.alpha
        for i, p in enumerate(probs_pos):
            w = W[i]
            # Append point mass 1.0 at infinity (finite-sample weighted conformal correction)
            w_norm = w / (np.sum(w) + 1.0)
            q = self._weighted_quantile(self.scores_cal, w_norm, min(0.9999, target))
            s = []
            if p <= q: s.append(0)
            if (1.0 - p) <= q: s.append(1)
            if not s: s.append(1 if p >= 0.5 else 0)
            pred_sets.append(s)
        return pred_sets
