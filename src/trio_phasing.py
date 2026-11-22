"""
Genomic Equity Conformal v4.0.0 — Global Multi-Sample Trio & Compound-Het Phasing Engine
Zero-dependency stdlib module layered directly on top of v3.0.0 (MIT License).
Detects Trans Compound Heterozygotes (AR), Cis Carriers, and De Novo candidates.
"""
from collections import defaultdict


def _parse_gt_alleles(gt: str):
    if not gt or gt in (".", "./.", ".|."):
        return None, False
    phased = "|" in gt
    sep = "|" if phased else "/"
    parts = gt.split(sep)
    if len(parts) != 2 or any(p not in ("0", "1", "2") for p in parts):
        return None, phased
    return (int(parts[0]), int(parts[1])), phased


def _is_het(alleles):
    return alleles is not None and ((alleles[0] > 0) != (alleles[1] > 0))


def _is_hom_alt(alleles):
    return alleles is not None and alleles[0] > 0 and alleles[1] > 0


def _is_hom_ref(alleles):
    return alleles is not None and alleles[0] == 0 and alleles[1] == 0


def analyze_trio_and_compound_het(records):
    """
    Input: list of dicts, each containing:
      - variant_id, gene, moi ('AR'|'AD'|'XL'), is_pathogenic_or_vus (bool)
      - proband_gt (str), optional mother_gt (str), optional father_gt (str)
    Returns: dict mapping variant_id -> trio/phasing annotation + gene-level summary.
    """
    by_gene_het = defaultdict(list)
    per_variant = {}

    for rec in records:
        vid = rec["variant_id"]
        gene = rec.get("gene", "UNKNOWN").upper()
        moi = rec.get("moi", "AR").upper()
        p_al, p_phased = _parse_gt_alleles(rec.get("proband_gt", ""))
        m_al, _ = _parse_gt_alleles(rec.get("mother_gt", ""))
        f_al, _ = _parse_gt_alleles(rec.get("father_gt", ""))

        inheritance = "Single_Sample_No_Parents"
        origin_parent = None

        if p_al is not None and m_al is not None and f_al is not None:
            m_has = not _is_hom_ref(m_al)
            f_has = not _is_hom_ref(f_al)
            if not _is_hom_ref(p_al) and not m_has and not f_has:
                inheritance = "De_Novo_Candidate"
            elif m_has and not f_has:
                inheritance = "Maternal_Inherited"
                origin_parent = "MATERNAL"
            elif f_has and not m_has:
                inheritance = "Paternal_Inherited"
                origin_parent = "PATERNAL"
            elif m_has and f_has:
                inheritance = "Biparental_Shared"
                origin_parent = "BOTH"

        per_variant[vid] = {
            "gene": gene,
            "moi": moi,
            "inheritance": inheritance,
            "origin_parent": origin_parent,
            "phased": p_phased,
            "alleles": p_al,
            "compound_het_status": "None",
        }

        if moi == "AR" and rec.get("is_pathogenic_or_vus", True) and _is_het(p_al):
            by_gene_het[gene].append(vid)

    gene_alerts = []
    for gene, vids in by_gene_het.items():
        if len(vids) < 2:
            continue
        origins = {per_variant[v]["origin_parent"] for v in vids if per_variant[v]["origin_parent"]}
        phased_haps = {per_variant[v]["alleles"] for v in vids if per_variant[v]["phased"]}

        if ("MATERNAL" in origins and "PATERNAL" in origins) or ((1, 0) in phased_haps and (0, 1) in phased_haps):
            status = "AR_Compound_Het_Trans_Affected"
        elif (len(origins) == 1 and origins <= {"MATERNAL", "PATERNAL"}) or (
            all(per_variant[v]["phased"] for v in vids) and len(phased_haps) == 1
        ):
            status = "AR_Cis_Carrier_Only"
        else:
            status = "AR_Unphased_Potential_Compound_Het"

        for v in vids:
            per_variant[v]["compound_het_status"] = status
        gene_alerts.append({"gene": gene, "variants": vids, "status": status})

    return {"per_variant": per_variant, "compound_het_alerts": gene_alerts}


if __name__ == "__main__":
    demo = [
        {"variant_id": "ATP7B_H1069Q", "gene": "ATP7B", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/1", "father_gt": "0/0", "is_pathogenic_or_vus": True},
        {"variant_id": "ATP7B_R778L",  "gene": "ATP7B", "moi": "AR", "proband_gt": "0/1", "mother_gt": "0/0", "father_gt": "0/1", "is_pathogenic_or_vus": True},
        {"variant_id": "TP53_R175H",   "gene": "TP53",  "moi": "AD", "proband_gt": "0/1", "mother_gt": "0/0", "father_gt": "0/0", "is_pathogenic_or_vus": True},
    ]
    import json
    print(json.dumps(analyze_trio_and_compound_het(demo), indent=2))
