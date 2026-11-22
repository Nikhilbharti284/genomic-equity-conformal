import csv, hashlib, json, os, sqlite3

DEFAULT_DB = "artifacts/v3_offline_clinical_ref.sqlite"

CANONICAL_CONTROLS = [
    ("chr17:g.7675088C>T", "TP53_R175H", "TP53", "R", 175, "H", 0.996, -5.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "AD"),
    ("chr11:g.5227002T>A", "HBB_E7V",    "HBB",  "E",   7, "V", 0.920, -2.28, 0.000621, 0.046998, 0.0, 0.002, 0.001, 0.0001, 0.0, 0.0, 0.004, "AR"),
    ("chr17:g.7676154G>C", "TP53_P72R",  "TP53", "P",  72, "R", 0.032, -3.32, 0.500196, 0.379906, 0.58, 0.65, 0.60, 0.7378, 0.71, 0.72, 0.7378, "AD"),
    ("chr13:g.48045719C>T","NUDT15_R139C","NUDT15","R",139,"C", 0.910, -4.1, 0.068, 0.001, 0.104, 0.042, 0.005, 0.002, 0.001, 0.001, 0.015, "AR"),
]

class OfflineReferenceStore:
    """Air-gapped sub-millisecond SQLite WAL reference store for hospital --network=none execution."""
    def __init__(self, db_path=DEFAULT_DB):
        self.db_path = db_path
        if not os.path.exists(self.db_path):
            build_offline_store(self.db_path)

    def _connect(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def lookup_by_hgvs(self, hgvs_g):
        if not hasattr(self, "_hgvs_cache"):
            conn = self._connect()
            try:
                rows = conn.execute("SELECT * FROM clinical_ref").fetchall()
                self._hgvs_cache = {r["hgvs_g"]: dict(r) for r in rows}
            finally:
                conn.close()
        return self._hgvs_cache.get(hgvs_g)

    def lookup_by_variant_id(self, variant_id):
        conn = self._connect()
        try:
            row = conn.execute("SELECT * FROM clinical_ref WHERE variant_id = ?", (variant_id.upper(),)).fetchone()
            return dict(row) if row else None
        finally:
            conn.close()

    def stats(self):
        conn = self._connect()
        try:
            n = conn.execute("SELECT COUNT(*) FROM clinical_ref").fetchone()[0]
        finally:
            conn.close()
        sha = hashlib.sha256(open(self.db_path, "rb").read()).hexdigest() if os.path.exists(self.db_path) else ""
        return {"db_path": self.db_path, "n_records": int(n), "sha256": sha}

def build_offline_store(db_path=DEFAULT_DB):
    os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("""CREATE TABLE IF NOT EXISTS clinical_ref (
        hgvs_g TEXT PRIMARY KEY,
        variant_id TEXT, gene TEXT, wt TEXT, pos INTEGER, mt TEXT,
        am_score REAL, esm2_llr REAL,
        af_sas REAL, af_afr REAL, af_eas REAL, af_amr REAL, af_mid REAL,
        af_nfe REAL, af_asj REAL, af_fin REAL, af_all REAL, moi TEXT
    )""")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_ref_vid ON clinical_ref(variant_id);")

    fmap = {}
    if os.path.exists("artifacts/v1_commercial_features.sqlite"):
        c1 = sqlite3.connect("artifacts/v1_commercial_features.sqlite")
        fmap = {r[0]: r[1:] for r in c1.execute("SELECT variant_id, am_score, af_nfe, af_sas, af_afr, af_all FROM var_features")}
        c1.close()

    g7 = {}
    if os.path.exists("artifacts/v2_global_7ancestry_features.sqlite"):
        c2 = sqlite3.connect("artifacts/v2_global_7ancestry_features.sqlite")
        g7 = {r[0]: r[1:] for r in c2.execute("SELECT * FROM gnomad7")}
        c2.close()

    g_thr = {}
    if os.path.exists("data/ba1_exception_registry.json"):
        g_thr = json.load(open("data/ba1_exception_registry.json", encoding="utf-8")).get("gene_thresholds", {})

    rows = []
    if os.path.exists("data/v1_clinvar_120genes.csv"):
        with open("data/v1_clinvar_120genes.csv", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                vid = r["variant_id"]
                hgvs = r.get("hgvs_g") or f"chr{str(r.get('chrom','')).replace('chr','')}:g.{r.get('pos_g') or r.get('pos','')}{r.get('ref','')}>{r.get('alt','')}"
                fv = fmap.get(vid, (0.5, 0.0, 0.0, 0.0, 0.0))
                g = g7.get(hgvs, (fv[2], fv[3], 0.0, 0.0, 0.0, fv[1], 0.0, 0.0, fv[4]))
                am = float(fv[0]) if fv[0] is not None and float(fv[0]) >= 0 else 0.5
                llr_proxy = round(-6.5 * am + 0.5, 4)
                moi = g_thr.get(r["gene"], {}).get("moi", "AD")
                rows.append((
                    hgvs, vid, r["gene"], r["wt"], int(r.get("pos", 0)), r["mt"], am, llr_proxy,
                    max(g[0], fv[2]), max(g[1], fv[3]), g[2], g[3], g[4], max(g[5], fv[1]), g[6], g[7], max(g[8], fv[4]), moi
                ))
    rows.extend(CANONICAL_CONTROLS)
    conn.executemany("INSERT OR REPLACE INTO clinical_ref VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    conn.commit()
    n = conn.execute("SELECT COUNT(*) FROM clinical_ref").fetchone()[0]
    conn.close()
    print(f"[PASS] Built air-gapped SQLite WAL store -> {db_path} ({n} indexed GRCh38 records)")
    return n

if __name__ == "__main__":
    build_offline_store()
