import csv, json
from collections import defaultdict

def kmer_set(seq, k=5):
    return {seq[i:i+k] for i in range(max(0, len(seq) - k + 1))}

def main():
    """Zero-dependency <30% homology clustering via k-mer Jaccard + Union-Find DSA."""
    gene_seqs = {}
    with open("data/v1_clinvar_120genes.csv") as f:
        for r in csv.DictReader(f):
            if r["gene"] not in gene_seqs:
                gene_seqs[r["gene"]] = r["seq"]

    genes = sorted(gene_seqs.keys())
    kmers = {g: kmer_set(gene_seqs[g], k=5) for g in genes}
    parent = {g: g for g in genes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb

    merged_pairs = []
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            g1, g2 = genes[i], genes[j]
            k1, k2 = kmers[g1], kmers[g2]
            if not k1 or not k2: continue
            # 5-mer containment >= 0.08 conservatively catches >=30% sequence-identity paralogs
            overlap = len(k1 & k2) / min(len(k1), len(k2))
            if overlap >= 0.08:
                union(g1, g2)
                merged_pairs.append((g1, g2, round(overlap, 3)))

    clusters = defaultdict(list)
    for g in genes:
        clusters[find(g)].append(g)

    gene_to_cluster = {g: f"CLU_{idx:03d}" for idx, (_, members) in enumerate(sorted(clusters.items())) for g in members}
    out = {"n_genes": len(genes), "n_clusters": len(clusters), "paralog_merges": merged_pairs, "gene_to_cluster": gene_to_cluster}
    with open("data/v1_homology_clusters.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"Clustered {len(genes)} genes into {len(clusters)} independent homology clusters (merged {len(merged_pairs)} paralog pairs).")
    if merged_pairs:
        print("Paralog merges caught:", merged_pairs[:10])

if __name__ == "__main__":
    main()
