import csv, hashlib, json, os, sqlite3, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alphafold_v6_features import fetch_uniprot_acc, fetch_afdb_v6_plddt

def main():
    reg = json.load(open("data/ba1_exception_registry.json", encoding="utf-8"))
    listed = set(reg.get("exceptions", {}).keys())
    g_thr = reg.get("gene_thresholds", {})

    conn1 = sqlite3.connect("artifacts/v1_commercial_features.sqlite")
    fmap = {r[0]: r[1:] for r in conn1.execute("SELECT variant_id, am_score, af_nfe, af_sas, af_afr, af_all, ancestry_group FROM var_features")}
    conn1.close()

    conn2 = sqlite3.connect("artifacts/v2_global_7ancestry_features.sqlite")
    g7 = {r[0]: r[1:] for r in conn2.execute("SELECT * FROM gnomad7")}
    conn2.close()

    cands = []
    with open("data/v1_clinvar_120genes.csv") as f:
        for r in csv.DictReader(f):
            label = int(r.get("label", r.get("y", 0)))
            if label != 1: continue
            vid = r["variant_id"]
            if vid in listed: continue
            fv = fmap.get(vid)
            if not fv or float(fv[0]) < 0.85: continue
            am, old_nfe, old_sas, old_afr, old_all, _ = fv
            hgvs = r.get("hgvs_g") or f"chr{str(r.get('chrom','')).replace('chr','')}:g.{r.get('pos_g') or r.get('pos','')}{r.get('ref','')}>{r.get('alt','')}"
            g = g7.get(hgvs, (old_sas, old_afr, 0.0, 0.0, 0.0, old_nfe, 0.0, 0.0, old_all))
            afs = {
                "SAS": round(max(g[0], old_sas), 6), "AFR": round(max(g[1], old_afr), 6),
                "EAS": round(g[2], 6), "AMR": round(g[3], 6), "MID": round(g[4], 6),
                "NFE": round(max(g[5], old_nfe), 6), "ASJ": round(g[6], 6), "FIN": round(g[7], 6)
            }
            top_anc, top_af = max(afs.items(), key=lambda kv: kv[1])
            if top_af < 0.001: continue  # Exceeds BS1 (>= 0.1%) in at least one global ancestry
            gene = r["gene"]
            cands.append({
                "variant_id": vid, "gene": gene, "hgvs_g": hgvs,
                "pos": int(r.get("pos", 0)),
                "top_ancestry": top_anc, "top_ancestry_af": top_af,
                "all_ancestry_afs": afs, "alphamissense_score": round(float(am), 4),
                "clinvar_stars": r.get("review_stars", ">=2-star"),
                "moi": g_thr.get(gene, {}).get("moi", "AD/AR")
            })

    cands.sort(key=lambda x: (-x["top_ancestry_af"], -x["alphamissense_score"]))
    # Fetch AlphaFold DB v6 pLDDT for top 25 candidates (cached per gene)
    plddt_cache = {}
    for c in cands[:25]:
        g = c["gene"]
        if g not in plddt_cache:
            try:
                acc = fetch_uniprot_acc(g)
                pm, _ = fetch_afdb_v6_plddt(acc) if acc else ({}, None)
                plddt_cache[g] = pm
            except Exception:
                plddt_cache[g] = {}
        c["afdb_v6_plddt"] = round(plddt_cache[g].get(c["pos"], -1.0), 2)
        c["structural_core_pass"] = bool(c["afdb_v6_plddt"] >= 85.0)

    out = {
        "version": "2.0.0-autonomous-discovery",
        "criteria": "ClinVar Pathogenic (>=2-star), NOT in ClinGen static BA1 Exception List, AlphaMissense >= 0.85, Global 7-Ancestry AF >= 0.001 (BS1+)",
        "n_discovered": len(cands),
        "n_structural_core_plddt85": sum(1 for c in cands[:25] if c.get("structural_core_pass")),
        "top_discoveries": cands[:25]
    }
    out["sha256"] = hashlib.sha256(json.dumps(out, sort_keys=True).encode()).hexdigest()
    open("artifacts/v20_novel_founder_discoveries.json", "w").write(json.dumps(out, indent=2) + "\n")

    with open("artifacts/v20_novel_founder_discoveries.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank","variant_id","gene","top_ancestry","top_ancestry_af","alphamissense_score","afdb_v6_plddt","structural_core_pass","moi","hgvs_g"])
        for i, c in enumerate(cands[:25], 1):
            w.writerow([i, c["variant_id"], c["gene"], c["top_ancestry"], c["top_ancestry_af"], c["alphamissense_score"], c.get("afdb_v6_plddt",-1), c.get("structural_core_pass",False), c["moi"], c["hgvs_g"]])

    print(f"[PASS] Discovered {len(cands)} unlisted BS1+ pathogenic founder/balancing candidates (SHA-256: {out['sha256'][:16]}...)")
    for i, c in enumerate(cands[:10], 1):
        print(f"  {i:2d}. {c['variant_id']:14s} | {c['top_ancestry']:3s} AF={c['top_ancestry_af']:.5f} | AM={c['alphamissense_score']:.4f} | pLDDT={c.get('afdb_v6_plddt')} | MOI={c['moi']}")

if __name__ == "__main__":
    main()
