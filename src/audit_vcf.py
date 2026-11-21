import argparse, hashlib, json, math, os, re, sys, urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inference_engine import GenomicEquityEngine

GEC_INFO_HEADERS = [
    '##INFO=<ID=GEC_PLATT,Number=1,Type=Float,Description="Genomic Equity Conformal v1.1.0 calibrated joint probability">',
    '##INFO=<ID=GEC_ACMG,Number=1,Type=String,Description="AF-free Pejaver PP3/BP4 ACMG tier">',
    '##INFO=<ID=GEC_TIER,Number=1,Type=String,Description="Non-empty Mondrian 90pct conformal + 3-Layer BA1 final tier">',
    '##INFO=<ID=GEC_SAFEGUARD,Number=1,Type=String,Description="Balancing-selection / BA1 safeguard status">',
    '##INFO=<ID=GEC_ANC,Number=1,Type=String,Description="Patient ancestry stratum used for Mondrian calibration">',
]

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

def _vcf_safe(s):
    return re.sub(r"[^A-Za-z0-9_.\-:]+", "_", str(s or "NONE")).strip("_") or "NONE"

def annotate_hg38_batch(hgvs_ids):
    """Auto-annotates raw CHROM/POS/REF/ALT variants with Gene, AA, and gnomAD AFs (zero REVEL)."""
    if not hgvs_ids: return {}
    out = {}
    url = "https://myvariant.info/v1/variant?assembly=hg38&fields=dbnsfp.genename,dbnsfp.aa,gnomad_exome.af,gnomad_genome.af"
    payload = json.dumps({"ids": hgvs_ids}).encode("utf-8")
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json", "User-Agent": "genomic-equity-conformal/1.1"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        for r in json.loads(resp.read().decode("utf-8")):
            q = r.get("query")
            if not q or r.get("notfound"): continue
            db = r.get("dbnsfp")
            if isinstance(db, list) and db: db = db[0]
            gene, wt, pos, mt = "", "", 0, ""
            if isinstance(db, dict):
                g = db.get("genename")
                gene = g[0] if isinstance(g, list) else (g or "")
                aa = db.get("aa")
                if isinstance(aa, dict):
                    wt, mt = aa.get("ref") or "", aa.get("alt") or ""
                    p_raw = aa.get("pos")
                    if isinstance(p_raw, list): p_raw = p_raw[0]
                    try: pos = int(p_raw)
                    except Exception: pos = 0
            af = {"af_nfe": 0.0, "af_sas": 0.0, "af_afr": 0.0, "af_all": 0.0}
            for k in ("gnomad_exome", "gnomad_genome"):
                sub = r.get(k)
                if isinstance(sub, list) and sub: sub = sub[0]
                if isinstance(sub, dict) and isinstance(sub.get("af"), dict):
                    a = sub["af"]
                    af = {"af_nfe": _num(a.get("af_nfe")), "af_sas": _num(a.get("af_sas")), "af_afr": _num(a.get("af_afr")), "af_all": _num(a.get("af"))}
                    if max(af.values()) > 0: break
            out[q] = {"gene": gene, "wt": wt, "pos": pos, "mt": mt, **af}
    return out

def fetch_seq(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
    req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/1.1"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        res = json.loads(resp.read().decode("utf-8")).get("results", [])
    return res[0]["sequence"]["value"] if res else None

def crop_window(seq, pos_1idx, max_len=512):
    if len(seq) <= max_len: return seq, pos_1idx
    start = max(0, (pos_1idx - 1) - max_len // 2)
    end = min(len(seq), start + max_len)
    return seq[max(0, end - max_len):end], pos_1idx - max(0, end - max_len)

def write_annotated_vcf(in_vcf, out_vcf, card_by_hgvs):
    """Writes a standard VCF 4.2 with GEC_* INFO annotations for Nextflow/bcftools pipelines."""
    os.makedirs(os.path.dirname(out_vcf) or ".", exist_ok=True)
    with open(in_vcf, encoding="utf-8") as fin, open(out_vcf, "w", encoding="utf-8") as fout:
        inserted = False
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                if not inserted:
                    for h in GEC_INFO_HEADERS:
                        fout.write(h + "\n")
                    inserted = True
                fout.write(line)
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            while len(cols) < 8:
                cols.append(".")
            chrom = cols[0].replace("chr", "")
            hgvs_g = f"chr{chrom}:g.{cols[1]}{cols[3]}>{cols[4]}"
            card = card_by_hgvs.get(hgvs_g)
            if card:
                tags = [
                    f"GEC_PLATT={card['joint_triage_prob']:.4f}",
                    f"GEC_ACMG={_vcf_safe(card['pejaver_acmg_code'])}",
                    f"GEC_TIER={_vcf_safe(card['final_tier'])}",
                    f"GEC_SAFEGUARD={_vcf_safe(card['safeguard_triggered'] or 'NONE')}",
                    f"GEC_ANC={_vcf_safe(card['patient_ancestry'])}",
                ]
                cols[7] = ";".join(tags) if cols[7] in ("", ".") else cols[7] + ";" + ";".join(tags)
            fout.write("\t".join(cols) + "\n")

def main():
    ap = argparse.ArgumentParser(description="Genomic Equity v1.1.0 — Unified Hospital VCF 4.2 Auditor")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-vcf", default="artifacts/audited_output.vcf")
    ap.add_argument("--out-json", default="artifacts/raw_vcf_audit_output.json")
    ap.add_argument("--default-ancestry", default="SAS", choices=["EUR", "SAS", "AFR", "ULTRA_RARE"])
    args = ap.parse_args()

    raw_rows = []
    with open(args.vcf, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            cols = line.strip().split("\t")
            if len(cols) < 5: continue
            chrom = cols[0].replace("chr", "")
            pos_g, rsid, ref, alt = cols[1], cols[2], cols[3], cols[4]
            info = dict(kv.split("=", 1) for kv in cols[7].split(";") if "=" in kv) if len(cols) > 7 else {}
            raw_rows.append({"chrom": chrom, "pos_g": pos_g, "rsid": rsid, "hgvs_g": f"chr{chrom}:g.{pos_g}{ref}>{alt}", "info": info})

    print(f"1/2 Auto-annotating {len(raw_rows)} raw VCF variants via MyVariant.info (GRCh38, zero REVEL)...")
    anns = annotate_hg38_batch([r["hgvs_g"] for r in raw_rows])

    eng = GenomicEquityEngine()
    import torch
    from transformers import AutoTokenizer, AutoModelForMaskedLM
    tok = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod.eval()

    print("2/2 Scoring with Unified GenomicEquityEngine (Platt + Non-Empty Mondrian + 3-Layer BA1 Safeguard)...")
    results, card_by_hgvs, seq_cache = [], {}, {}
    for r in raw_rows:
        a = anns.get(r["hgvs_g"], {})
        gene = r["info"].get("GENE") or a.get("gene", "")
        anc = r["info"].get("ANC") or args.default_ancestry
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", r["info"].get("AA", ""))
        wt, pos, mt = (m.group(1), int(m.group(2)), m.group(3)) if m else (a.get("wt", ""), int(a.get("pos", 0)), a.get("mt", ""))
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

        card = eng.predict(
            gene=gene, wt=wt, pos=pos, mt=mt, llr=llr, enorm=enorm, ancestry=anc,
            af_nfe=a.get("af_nfe", 0.0), af_sas=a.get("af_sas", 0.0),
            af_afr=a.get("af_afr", 0.0), af_all=a.get("af_all", 0.0),
        )
        card["hgvs_g"] = r["hgvs_g"]
        card["live_af_sas"] = round(a.get("af_sas", 0.0), 6)
        card["live_af_afr"] = round(a.get("af_afr", 0.0), 6)
        results.append(card)
        card_by_hgvs[r["hgvs_g"]] = card

    write_annotated_vcf(args.vcf, args.out_vcf, card_by_hgvs)
    payload = {"engine_version": "1.1.0-unified", "commercial_license_clean": True, "vcf_source": args.vcf, "annotated_vcf": args.out_vcf, "n_audited": len(results), "results": results}
    payload["sha256"] = hashlib.sha256(json.dumps(results, sort_keys=True).encode()).hexdigest()
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps(payload, indent=2))
    print(f"[PASS] Annotated VCF 4.2 written -> {args.out_vcf}")

if __name__ == "__main__":
    main()
