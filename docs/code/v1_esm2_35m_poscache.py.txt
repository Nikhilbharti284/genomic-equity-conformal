import argparse, csv, json, os, signal, sqlite3, time
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM
from resource_logger import ResourceLogger

STOP = False
def _sig(s, f):
    global STOP
    STOP = True
signal.signal(signal.SIGTERM, _sig)
signal.signal(signal.SIGINT, _sig)

AA20 = "ACDEFGHIKLMNPQRSTVWY"

def crop_window(seq, pos_1idx, max_len=512):
    if len(seq) <= max_len: return seq, pos_1idx
    start = max(0, (pos_1idx - 1) - max_len // 2)
    end = min(len(seq), start + max_len)
    start = max(0, end - max_len)
    return seq[start:end], pos_1idx - start

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="data/v1_clinvar_120genes.csv")
    ap.add_argument("--db", default="artifacts/esm2_35m_poscache.sqlite")
    ap.add_argument("--max-wall-min", type=int, default=40)
    ap.add_argument("--limit-positions", type=int, default=0, help="0 = all unique positions")
    args = ap.parse_args()

    logger = ResourceLogger()
    with open(args.csv) as f:
        rows = list(csv.DictReader(f))

    # Deduplicate to unique (uniprot_acc, pos) — saves ~35% CPU forward passes!
    uniq_pos = {}
    for r in rows:
        k = f"{r['uniprot_acc']}:{r['pos']}"
        if k not in uniq_pos:
            uniq_pos[k] = (r["uniprot_acc"], int(r["pos"]), r["seq"])

    os.makedirs(os.path.dirname(args.db), exist_ok=True)
    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("CREATE TABLE IF NOT EXISTS pos_logits (pos_key TEXT PRIMARY KEY, logits_json TEXT, emb_norm REAL)")
    done = {r[0] for r in conn.execute("SELECT pos_key FROM pos_logits")}
    todo = [(k, *v) for k, v in uniq_pos.items() if k not in done]
    if args.limit_positions > 0:
        todo = todo[:args.limit_positions]

    print(f"Unique positions total: {len(uniq_pos)} | Already cached: {len(done)} | Scoring now: {len(todo)}")
    if not todo:
        conn.close()
        return

    model_name = "facebook/esm2_t12_35M_UR50D"
    tok = AutoTokenizer.from_pretrained(model_name)
    mod = AutoModelForMaskedLM.from_pretrained(model_name)
    mod.eval()
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA20}

    t0 = time.time()
    scored = 0
    for idx, (k, acc, pos, seq) in enumerate(todo, 1):
        if STOP or (time.time() - t0) > args.max_wall_min * 60:
            print(f"Clean checkpoint exit at {scored}/{len(todo)} positions. Re-run to resume.")
            break
        wseq, wpos = crop_window(seq, pos, max_len=512)
        inp = tok(wseq, return_tensors="pt")
        inp["input_ids"][0, wpos] = tok.mask_token_id
        with torch.no_grad():
            out = mod(**inp, output_hidden_states=True)
            lp = torch.log_softmax(out.logits[0, wpos], dim=-1)
            vec = {aa: round(float(lp[tid]), 5) for aa, tid in aa_ids.items()}
            enorm = round(float(out.hidden_states[-1][0, wpos].norm().item()), 4)
        conn.execute("INSERT OR REPLACE INTO pos_logits VALUES (?, ?, ?)", (k, json.dumps(vec), enorm))
        scored += 1
        if scored % 25 == 0:
            conn.commit()
            print(f"  Progress: {scored}/{len(todo)} unique positions cached ({round(time.time()-t0,1)}s)...")
    conn.commit()
    conn.close()
    stats = logger.log("esm2_35m_poscache", {"scored_positions": scored, "total_cached": len(done) + scored})
    print("Run summary:", stats)

if __name__ == "__main__":
    main()
