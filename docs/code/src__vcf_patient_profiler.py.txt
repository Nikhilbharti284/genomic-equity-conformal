import math

ANCESTRIES_7 = ("SAS", "AFR", "EAS", "AMR", "MID", "NFE", "BOTTLENECK_FOUNDER")

def parse_sample_gt(format_str, sample_str, chrom=""):
    """Parses VCF 4.2 FORMAT & SAMPLE fields (GT, DP, GQ) into clinical zygosity."""
    if not format_str or not sample_str:
        return {"gt": "./.", "zygosity": "UNKNOWN", "dp": None, "gq": None}
    keys = format_str.split(":")
    vals = sample_str.split(":")
    fmap = dict(zip(keys, vals))
    gt = fmap.get("GT", "./.")
    dp = int(fmap["DP"]) if fmap.get("DP", "").isdigit() else None
    gq = int(fmap["GQ"]) if fmap.get("GQ", "").isdigit() else None
    alleles = gt.replace("|", "/").split("/")
    ch = str(chrom).upper().replace("CHR", "")
    if "." in alleles:
        zyg = "UNKNOWN"
    elif len(alleles) == 1:
        zyg = "HEMIZYGOUS" if alleles[0] != "0" and ch in ("X", "Y") else ("HOMOZYGOUS_REF" if alleles[0] == "0" else "HEMIZYGOUS")
    else:
        a1, a2 = alleles[0], alleles[1]
        if a1 == "0" and a2 == "0":
            zyg = "HOMOZYGOUS_REF"
        elif a1 != "0" and a2 != "0":
            zyg = "HOMOZYGOUS"
        else:
            zyg = "HETEROZYGOUS_CARRIER"
    return {"gt": gt, "zygosity": zyg, "dp": dp, "gq": gq}

def classify_clinical_zygosity(moi, zygosity):
    """Combines Mode of Inheritance (AD/AR/XL) with sample GT zygosity."""
    m = str(moi or "AD").upper()
    z = str(zygosity or "UNKNOWN").upper()
    if z == "HOMOZYGOUS_REF":
        return "Reference_Wildtype"
    if z == "UNKNOWN":
        return f"{m}_Zygosity_Not_Provided"
    if "AR" in m and "AD" not in m:
        if z == "HETEROZYGOUS_CARRIER":
            return "AR_Carrier_Only"
        if z in ("HOMOZYGOUS", "HEMIZYGOUS"):
            return "AR_Biallelic_Affected"
    if "XL" in m:
        if z in ("HEMIZYGOUS", "HOMOZYGOUS"):
            return "XL_Hemizygous_Affected"
        return "XL_Carrier_Monitor"
    if z == "HOMOZYGOUS":
        return "AD_Homozygous_Affected"
    return "AD_Heterozygous_Affected"

def infer_patient_ancestry_from_vcf(observed_af_dicts):
    """AIM Maximum-Likelihood Estimator across 7 global gnomAD ancestries from observed VCF variants."""
    if not observed_af_dicts:
        return "SAS", {a: round(1.0 / len(ANCESTRIES_7), 4) for a in ANCESTRIES_7}
    eps = 1e-5
    log_lik = {a: 0.0 for a in ANCESTRIES_7}
    used = 0
    for afs in observed_af_dicts:
        row = {
            "SAS": float(afs.get("af_sas", 0.0)),
            "AFR": float(afs.get("af_afr", 0.0)),
            "EAS": float(afs.get("af_eas", 0.0)),
            "AMR": float(afs.get("af_amr", 0.0)),
            "MID": float(afs.get("af_mid", 0.0)),
            "NFE": float(afs.get("af_nfe", 0.0)),
            "BOTTLENECK_FOUNDER": max(float(afs.get("af_asj", 0.0)), float(afs.get("af_fin", 0.0))),
        }
        if max(row.values()) <= 0.0:
            continue
        used += 1
        for a in ANCESTRIES_7:
            log_lik[a] += math.log(max(row[a], eps))
    if used == 0:
        return "SAS", {a: round(1.0 / len(ANCESTRIES_7), 4) for a in ANCESTRIES_7}
    m_ll = max(log_lik.values())
    exps = {a: math.exp(log_lik[a] - m_ll) for a in ANCESTRIES_7}
    tot = sum(exps.values())
    post = {a: round(exps[a] / tot, 4) for a in ANCESTRIES_7}
    best = max(post.items(), key=lambda kv: kv[1])[0]
    return best, post
