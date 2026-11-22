import os, sqlite3, sys
import numpy as np
sys.path.insert(0, os.path.abspath("src"))
from mondrian_conformal import AncestryMondrianConformal
from esm2_scorer import score_variants

def test_mondrian_ancestry_coverage():
    rng = np.random.default_rng(42)
    # Simulate European (calibrated) vs SAS/AFR (biased raw model)
    groups = np.array(["EUR"] * 400 + ["SAS"] * 400 + ["AFR"] * 400)
    y = rng.integers(0, 2, size=1200)
    probs = np.clip(y * 0.8 + (1 - y) * 0.2 + rng.normal(0, 0.12, 1200), 0.01, 0.99)
    # Inject ancestry bias into SAS & AFR
    probs[groups == "SAS"] = np.clip(probs[groups == "SAS"] * 0.65, 0.01, 0.99)
    probs[groups == "AFR"] = np.clip(probs[groups == "AFR"] * 0.60, 0.01, 0.99)

    cal_idx = np.arange(0, 1200, 2)
    test_idx = np.arange(1, 1200, 2)

    mc = AncestryMondrianConformal(alpha=0.10).fit(probs[cal_idx], y[cal_idx], groups[cal_idx])
    rep = mc.evaluate_by_group(probs[test_idx], y[test_idx], groups[test_idx])
    for g in ["EUR", "SAS", "AFR"]:
        assert rep[g]["coverage"] >= 0.85, f"Coverage too low for {g}: {rep[g]}"

def test_zero_gene_leakage():
    train_genes = {"BRCA1", "TP53", "MLH1"}
    test_genes = {"BRCA2", "MSH2", "PTEN"}
    assert len(train_genes.intersection(test_genes)) == 0

def test_sqlite_checkpoint_resume(tmp_path):
    db = str(tmp_path / "cache.sqlite")
    demo = [{"id": "V1", "seq": "ACDEF", "pos": 1, "wt": "A", "mt": "G"}]
    score_variants(demo, db_path=db, dry_run=True)
    conn = sqlite3.connect(db)
    rows = conn.execute("SELECT variant_id, llr FROM esm2_scores").fetchall()
    conn.close()
    assert len(rows) == 1 and rows[0][0] == "V1"
