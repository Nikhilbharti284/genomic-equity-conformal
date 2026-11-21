import argparse, json, os, re, urllib.request
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM
from inference_engine import GenomicEquityEngine

def fetch_seq(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
    req = urllib.request.Request(url, headers={"User-Agent": "genomic-equity-conformal/1.1"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        return json.loads(resp.read().decode("utf-8"))["results"][0]["sequence"]["value"]

def crop_window(seq, pos_1idx, max_len=512):
    if len(seq) <= max_len: return seq, pos_1idx
    start = max(0, (pos_1idx - 1) - max_len // 2)
    end = min(len(seq), start + max_len)
    return seq[max(0, end - max_len):end], pos_1idx - max(0, end - max_len)

def main():
    p = argparse.ArgumentParser(description="Genomic Equity v1.1.0 — Unified Variant Equity Card CLI")
    p.add_argument("--gene", required=True)
    p.add_argument("--variant", required=True)
    p.add_argument("--ancestry", default="SAS", choices=["EUR", "SAS", "AFR", "ULTRA_RARE"])
    p.add_argument("--af-nfe", type=float, default=0.0)
    p.add_argument("--af-sas", type=float, default=0.0)
    p.add_argument("--af-afr", type=float, default=0.0)
    args = p.parse_args()

    m = re.match(r"^([A-Z])(\d+)([A-Z])$", args.variant.upper())
    if not m: raise SystemExit("Use 1-letter variant format like R175H")
    wt, pos, mt = m.group(1), int(m.group(2)), m.group(3)
    seq = fetch_seq(args.gene.upper())
    if pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
        raise SystemExit(f"Sequence mismatch at {args.gene} pos {pos}")

    wseq, wpos = crop_window(seq, pos)
    tok = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mod.eval()
    inp = tok(wseq, return_tensors="pt")
    inp["input_ids"][0, wpos] = tok.mask_token_id
    with torch.no_grad():
        out = mod(**inp, output_hidden_states=True)
        lp = torch.log_softmax(out.logits[0, wpos], dim=-1)
        llr = float(lp[tok.convert_tokens_to_ids(mt)] - lp[tok.convert_tokens_to_ids(wt)])
        enorm = float(out.hidden_states[-1][0, wpos].norm().item())

    eng = GenomicEquityEngine()
    card = eng.predict(gene=args.gene, wt=wt, pos=pos, mt=mt, llr=llr, enorm=enorm, ancestry=args.ancestry,
                       af_nfe=args.af_nfe, af_sas=args.af_sas, af_afr=args.af_afr)
    os.makedirs("artifacts/cards", exist_ok=True)
    out_json = f"artifacts/cards/{card['variant_id']}_{args.ancestry}_v11.json"
    with open(out_json, "w") as f:
        json.dump(card, f, indent=2)
    print(json.dumps(card, indent=2))
    print(f"Saved v1.1.0 Equity Card -> {out_json}")

if __name__ == "__main__":
    main()
