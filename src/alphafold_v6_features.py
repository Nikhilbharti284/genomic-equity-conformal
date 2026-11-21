import json, urllib.request

def fetch_uniprot_acc(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
    req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/1.1"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        res = json.loads(resp.read().decode("utf-8")).get("results", [])
    return res[0]["primaryAccession"] if res else None

def fetch_afdb_v6_plddt(uniprot_acc):
    """Dynamically resolves AlphaFold DB v6 URL (avoids 404 _v4 links) and parses CA pLDDT."""
    meta_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}"
    req = urllib.request.Request(meta_url, headers={"User-Agent": "genomic-equity-conformal/1.1"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        meta = json.loads(resp.read().decode("utf-8"))
    if not meta:
        return {}, None
    pdb_url = meta[0].get("pdbUrl")
    if not pdb_url:
        return {}, None
    req2 = urllib.request.Request(pdb_url, headers={"User-Agent": "genomic-equity-conformal/1.1"})
    plddt = {}
    with urllib.request.urlopen(req2, timeout=30) as resp:
        for raw in resp:
            line = raw.decode("utf-8", errors="ignore")
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    plddt[int(line[22:26])] = float(line[60:66])
                except Exception:
                    pass
    return plddt, pdb_url

if __name__ == "__main__":
    for gene, pos, label in [("TP53", 175, "R175H"), ("BRCA1", 1, "M1V"), ("HBB", 7, "E7V")]:
        acc = fetch_uniprot_acc(gene)
        plddt_map, url = fetch_afdb_v6_plddt(acc)
        print(f"{gene}_{label} ({acc}) -> AlphaFold DB v6 pLDDT[{pos}] = {plddt_map.get(pos)} | model={url.split('/')[-1] if url else 'N/A'}")
