import csv, hashlib, io, json, math, os, re, sqlite3, urllib.request
import numpy as np
import torch
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from transformers import AutoTokenizer, AutoModelForMaskedLM
from resource_logger import ResourceLogger

AA3TO1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
HYDRO  = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
MW     = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,"G":75,"H":155,"I":131,"L":131,"K":146,"M":149,"F":165,"P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}
CHG    = {"R":1,"K":1,"H":0.1,"D":-1,"E":-1}
HGVS_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

CANDIDATES = [
    ("TP53", "urn:mavedb:00000068-a-1"),
    ("PTEN", "urn:mavedb:00000013-a-1"),
    ("BRCA1", "urn:mavedb:00000003-a-1"),
]

def fetch_seq(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
    req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/0.6"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        return json.loads(resp.read().decode("utf-8"))["results"][0]["sequence"]["value"]

def fetch_mavedb(max_variants=120):
    for gene, urn in CANDIDATES:
        url = f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores"
        req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/0.6"})
        try:
            with urllib.request.urlopen(req, timeout=25) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            seq = fetch_seq(gene)
            parsed = []
            for row in csv.DictReader(io.StringIO(text)):
                m = HGVS_RE.match((row.get("hgvs_pro") or "").strip())
                if not m: continue
                wt3, pos_s, mt3 = m.groups()
                if wt3 not in AA3TO1 or mt3 not in AA3TO1 or wt3 == mt3: continue
                wt, mt, pos = AA3TO1[wt3], AA3TO1[mt3], int(pos_s)
                if pos < 1 or pos > len(seq) or seq[pos - 1] != wt: continue
                try: sc = float(row["score"])
                except Exception: continue
                if math.isnan(sc): continue
                parsed.append({"id": f"{gene}_{wt}{pos}{mt}", "gene": gene, "wt": wt, "pos": pos, "mt": mt, "dms_score": sc})
            if len(parsed) >= 50:
                rng = np.random.default_rng(42)
                idx = rng.choice(len(parsed), size=min(max_variants, len(parsed)), replace=False)
                return gene, urn, seq, [parsed[i] for i in idx]
        except Exception as e:
            print(f"[info] MaveDB {urn} skipped: {e}")
    raise RuntimeError("Could not reach MaveDB API right now.")

def crop_window(seq, pos_1idx, max_len=1020):
    if len(seq) <= max_len: return seq, pos_1idx
    start = max(0, (pos_1idx - 1) - max_len // 2)
    end = min(len(seq), start + max_len)
    return seq[max(0, end - max_len):end], pos_1idx - max(0, end - max_len)

def main():
    logger = ResourceLogger()
    print("1/3 Downloading independent wet-lab DMS ground truth from MaveDB...")
    gene, urn, seq, dms = fetch_mavedb(max_variants=100)
    print(f"    Loaded {len(dms)} experimental wet-lab variants for {gene} ({urn})")

    print(f"2/3 Training REVEL-free model EXCLUDING {gene} (zero gene leakage, zero circularity)...")
    with open("data/pilot_variants_hg38.csv") as f:
        train_rows = [r for r in csv.DictReader(f) if r["gene"] != gene]
    conn = sqlite3.connect("artifacts/esm2_pilot.sqlite")
    cache = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT variant_id, llr, emb_norm FROM scores")}
    conn.close()

    def feat_row(wt, mt, llr, enorm):
        dh = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
        dm = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
        dc = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
        return [llr, enorm, dh, abs(dh), dm, abs(dm), dc, float(wt in ("P","G") or mt in ("P","G"))]

    X_tr, y_tr = [], []
    for r in train_rows:
        if r["variant_id"] in cache:
            X_tr.append(feat_row(r["wt"], r["mt"], *cache[r["variant_id"]]))
            y_tr.append(int(r["label"]))
    X_tr, y_tr = np.array(X_tr, dtype=float), np.array(y_tr)
    sc = StandardScaler().fit(X_tr)
    mlp = LogisticRegression(C=0.2, max_iter=500).fit(sc.transform(X_tr), y_tr)

    print(f"3/3 Scoring {len(dms)} MaveDB wet-lab variants with Frozen ESM-2...")
    tok = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod.eval()
    raw_llrs, model_probs, lab_scores = [], [], []
    for d in dms:
        wseq, wpos = crop_window(seq, d["pos"])
        inp = tok(wseq, return_tensors="pt")
        inp["input_ids"][0, wpos] = tok.mask_token_id
        with torch.no_grad():
            out = mod(**inp, output_hidden_states=True)
            lp = torch.log_softmax(out.logits[0, wpos], dim=-1)
            llr = float(lp[tok.convert_tokens_to_ids(d["mt"])] - lp[tok.convert_tokens_to_ids(d["wt"])])
            enorm = float(out.hidden_states[-1][0, wpos].norm().item())
        x = np.array([feat_row(d["wt"], d["mt"], llr, enorm)], dtype=float)
        prob = float(mlp.predict_proba(sc.transform(x))[0, 1])
        raw_llrs.append(llr)
        model_probs.append(prob)
        lab_scores.append(d["dms_score"])

    rho_esm2, p_esm2 = spearmanr(raw_llrs, lab_scores)
    rho_mod, p_mod = spearmanr(model_probs, lab_scores)
    rep = {
        "mavedb_urn": urn,
        "held_out_gene": gene,
        "n_wetlab_variants": len(dms),
        "raw_esm2_abs_spearman": round(abs(float(rho_esm2)), 4),
        "raw_esm2_pvalue": float(f"{p_esm2:.2e}"),
        "revel_free_model_abs_spearman": round(abs(float(rho_mod)), 4),
        "revel_free_model_pvalue": float(f"{p_mod:.2e}"),
        "beats_raw_esm2": bool(abs(float(rho_mod)) >= abs(float(rho_esm2))),
    }
    out_json = "artifacts/mavedb_external_validation.json"
    with open(out_json, "w") as f:
        json.dump(rep, f, indent=2)
    rep["sha256"] = logger.sha256_file(out_json)
    logger.log("mavedb_external_validation", rep)
    print(json.dumps(rep, indent=2))

if __name__ == "__main__":
    main()
