import argparse, hashlib, json, math, os, re, urllib.request
import joblib
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

HYDRO = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
MW    = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,"G":75,"H":155,"I":131,"L":131,"K":146,"M":149,"F":165,"P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}
CHG   = {"R":1,"K":1,"H":0.1,"D":-1,"E":-1}

def _logit(p):
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p / (1.0 - p))

def _num(x, default=0.0):
    if x is None: return default
    if isinstance(x, list):
        vals = [_num(v, None) for v in x]
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else default
    try:
        v = float(x)
        return default if math.isnan(v) else v
    except Exception:
        return default

def annotate_hg38_batch(hgvs_ids):
    """Auto-annotates raw CHROM/POS/REF/ALT variants with Gene, AA, gnomAD AFs, and REVEL."""
    out = {}
    url = "https://myvariant.info/v1/variant?assembly=hg38&fields=dbnsfp.genename,dbnsfp.aa,dbnsfp.revel.score,gnomad_exome.af,gnomad_genome.af,dbnsfp.gnomad_exomes"
    payload = json.dumps({"ids": hgvs_ids}).encode("utf-8")
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/0.6"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        for r in json.loads(resp.read().decode("utf-8")):
            q = r.get("query")
            if not q or r.get("notfound"): continue
            db = r.get("dbnsfp")
            if isinstance(db, list) and db: db = db[0]
            gene, wt, pos, mt, rev = "", "", 0, "", None
            if isinstance(db, dict):
                g = db.get("genename")
                gene = g[0] if isinstance(g, list) else (g or "")
                aa = db.get("aa")
                if isinstance(aa, dict):
                    wt = aa.get("ref") or ""
                    mt = aa.get("alt") or ""
                    p_raw = aa.get("pos")
                    if isinstance(p_raw, list): p_raw = p_raw[0]
                    try: pos = int(p_raw)
                    except Exception: pos = 0
                if isinstance(db.get("revel"), dict):
                    rev = _num(db["revel"].get("score"), None)
            af = {"af_nfe": 0.0, "af_sas": 0.0, "af_afr": 0.0, "af_all": 0.0}
            for k in ("gnomad_exome", "gnomad_genome"):
                sub = r.get(k)
                if isinstance(sub, list) and sub: sub = sub[0]
                if isinstance(sub, dict) and isinstance(sub.get("af"), dict):
                    a = sub["af"]
                    af = {"af_nfe": _num(a.get("af_nfe")), "af_sas": _num(a.get("af_sas")), "af_afr": _num(a.get("af_afr")), "af_all": _num(a.get("af"))}
                    if max(af.values()) > 0: break
            out[q] = {"gene": gene, "wt": wt, "pos": pos, "mt": mt, "revel": rev, **af}
    return out

def fetch_seq(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
    req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/0.6"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        res = json.loads(resp.read().decode("utf-8")).get("results", [])
    return res[0]["sequence"]["value"] if res else None

def crop_window(seq, pos_1idx, max_len=1020):
    if len(seq) <= max_len: return seq, pos_1idx
    start = max(0, (pos_1idx - 1) - max_len // 2)
    end = min(len(seq), start + max_len)
    return seq[max(0, end - max_len):end], pos_1idx - max(0, end - max_len)

def main():
    ap = argparse.ArgumentParser(description="Raw Hospital VCF 4.2 Auto-Annotating Conformal Auditor")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--default-ancestry", default="SAS", choices=["EUR", "SAS", "AFR", "ULTRA_RARE"])
    args = ap.parse_args()

    raw_rows = []
    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            cols = line.strip().split("\t")
            if len(cols) < 5: continue
            chrom = cols[0].replace("chr", "")
            pos_g, rsid, ref, alt = cols[1], cols[2], cols[3], cols[4]
            info = dict(kv.split("=", 1) for kv in cols[7].split(";") if "=" in kv) if len(cols) > 7 else {}
            raw_rows.append({"chrom": chrom, "pos_g": pos_g, "rsid": rsid, "hgvs_g": f"chr{chrom}:g.{pos_g}{ref}>{alt}", "info": info})

    print(f"1/2 Auto-annotating {len(raw_rows)} raw VCF variants via MyVariant.info (GRCh38)...")
    anns = annotate_hg38_batch([r["hgvs_g"] for r in raw_rows])

    bundle = joblib.load("artifacts/equity_model_v4.joblib")
    tok = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod.eval()

    print("2/2 Scoring with Frozen ESM-2 + Live gnomAD AF + REVEL + Mondrian Conformal...")
    results, seq_cache = [], {}
    for r in raw_rows:
        a = anns.get(r["hgvs_g"], {})
        gene = r["info"].get("GENE") or a.get("gene", "")
        anc = r["info"].get("ANC") or args.default_ancestry
        aa_info = r["info"].get("AA", "")
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", aa_info)
        if m:
            wt, pos, mt = m.group(1), int(m.group(2)), m.group(3)
        else:
            wt, pos, mt = a.get("wt", ""), int(a.get("pos", 0)), a.get("mt", "")
        if not gene or not wt or not mt or pos < 1: continue
        if gene not in seq_cache: seq_cache[gene] = fetch_seq(gene)
        seq = seq_cache[gene]
        if not seq or pos > len(seq) or seq[pos - 1] != wt: continue

        wseq, wpos = crop_window(seq, pos)
        inp = tok(wseq, return_tensors="pt")
        inp["input_ids"][0, wpos] = tok.mask_token_id
        with torch.no_grad():
            out = mod(**inp, output_hidden_states=True)
            lp = torch.log_softmax(out.logits[0, wpos], dim=-1)
            llr = float(lp[tok.convert_tokens_to_ids(mt)] - lp[tok.convert_tokens_to_ids(wt)])
            enorm = float(out.hidden_states[-1][0, wpos].norm().item())

        af_nfe, af_sas, af_afr, af_all = a.get("af_nfe", 0.0), a.get("af_sas", 0.0), a.get("af_afr", 0.0), a.get("af_all", 0.0)
        d_h = HYDRO.get(mt, 0.0) - HYDRO.get(wt, 0.0)
        d_m = (MW.get(mt, 120.0) - MW.get(wt, 120.0)) / 100.0
        d_c = abs(CHG.get(mt, 0.0) - CHG.get(wt, 0.0))
        feat = [llr, enorm, d_h, abs(d_h), d_m, abs(d_m), d_c, float(wt in ("P","G") or mt in ("P","G")),
                math.log10(af_nfe + 1e-7), math.log10(af_sas + 1e-7), math.log10(af_afr + 1e-7), math.log10(af_all + 1e-7)]
        if bundle["use_cadd"]: feat.append(bundle["med_cadd"])
        if bundle["use_rev"]:  feat.append(a["revel"] if a.get("revel") is not None else bundle["med_rev"])
        x = np.array([feat], dtype=float)
        xs = bundle["scaler"].transform(x)
        p_raw = float(0.45 * bundle["mlp"].predict_proba(xs)[0, 1] + 0.35 * bundle["hgb"].predict_proba(x)[0, 1] + 0.20 * bundle["logreg"].predict_proba(xs)[0, 1])
        p_cal = float(bundle["platt"].predict_proba(_logit(np.array([p_raw])).reshape(-1, 1))[0, 1])
        ps = bundle["mondrian"].predict_sets(np.array([p_cal]), np.array([anc]))[0]
        tier = "Pathogenic_90pct" if ps == [1] else ("Benign_90pct" if ps == [0] else "VUS_Review_Required")
        if tier == "Benign_90pct" and max(af_sas, af_afr, af_nfe) >= 0.01 and llr <= -2.0 and (a.get("revel") or 0.0) >= 0.50:
            tier = "VUS_Review_Required (High-AF Functional Conflict)" 
        results.append({
            "hgvs_g": r["hgvs_g"], "variant": f"{gene}_{wt}{pos}{mt}", "ancestry": anc,
            "live_af_sas": round(af_sas, 6), "live_af_afr": round(af_afr, 6), "live_revel": a.get("revel"),
            "esm2_llr": round(llr, 3), "platt_prob": round(p_cal, 4), "conformal_tier": tier
        })

    payload = {"vcf_source": args.vcf, "auto_annotated": True, "n_audited": len(results), "results": results}
    payload["sha256"] = hashlib.sha256(json.dumps(results, sort_keys=True).encode()).hexdigest()
    with open("artifacts/raw_vcf_audit_output.json", "w") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps(payload, indent=2))

if __name__ == "__main__":
    main()
