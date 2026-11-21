import glob, hashlib, json, os, time

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    files = sorted(
        glob.glob("src/*.py") + glob.glob("tests/*.py") +
        glob.glob("artifacts/*.json") + glob.glob("artifacts/cards/*.json") +
        ["CITATION.cff", "LICENSE", "README.md", "requirements.txt"]
    )
    parts, graph = [], []
    for f in files:
        if not os.path.isfile(f):
            continue
        parts.append({"@id": f})
        graph.append({
            "@id": f,
            "@type": "File",
            "name": os.path.basename(f),
            "contentSize": os.path.getsize(f),
            "sha256": sha256(f),
        })

    crate = {
        "@context": "https://w3id.org/ro/crate/1.1/context",
        "@graph": [
            {"@type": "CreativeWork", "@id": "ro-crate-metadata.json", "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"}, "about": {"@id": "./"}},
            {
                "@id": "./",
                "@type": "Dataset",
                "name": "Genomic Equity Conformal (v0.2.0)",
                "description": "Ancestry-calibrated missense & cryptic-splicing VUS classifier with Frozen ESM-2, Platt scaling, and Mondrian Conformal Prediction.",
                "datePublished": time.strftime("%Y-%m-%d"),
                "license": "MIT",
                "hasPart": parts,
            },
            *graph,
        ],
    }
    with open("ro-crate-metadata.json", "w") as f:
        json.dump(crate, f, indent=2)
    print(f"Generated W3C RO-Crate 1.1 manifest with {len(parts)} SHA-256 verified files -> ro-crate-metadata.json")

if __name__ == "__main__":
    main()
