import argparse, os, signal, sqlite3, sys, time
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM
from resource_logger import ResourceLogger

STOP_REQUESTED = False
def _handle_sig(signum, frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True

signal.signal(signal.SIGTERM, _handle_sig)
signal.signal(signal.SIGINT, _handle_sig)

def init_db(db_path):
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("""
        CREATE TABLE IF NOT EXISTS esm2_scores (
            variant_id TEXT PRIMARY KEY,
            llr REAL,
            emb_norm REAL
        )
    """)
    conn.commit()
    return conn

def score_variants(variants, db_path="artifacts/esm2_cache.sqlite", max_wall_min=45, dry_run=False):
    logger = ResourceLogger()
    conn = init_db(db_path)
    done = {row[0] for row in conn.execute("SELECT variant_id FROM esm2_scores")}
    todo = [v for v in variants if v["id"] not in done]

    if not todo:
        print("All variants already cached in SQLite WAL.")
        return

    if not dry_run:
        model_name = "facebook/esm2_t6_8M_UR50D"
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModelForMaskedLM.from_pretrained(model_name)
        model.eval()

    t_start = time.time()
    scored = 0
    for v in todo:
        if STOP_REQUESTED or (time.time() - t_start) > max_wall_min * 60:
            print("Clean checkpoint exit triggered before Codespaces timeout.")
            break
        if dry_run:
            llr, emb_norm = -1.25, 10.0
        else:
            seq = v["seq"]
            pos = v["pos"]  # 1-indexed
            wt, mt = v["wt"], v["mt"]
            inputs = tokenizer(seq, return_tensors="pt")
            mask_idx = pos  # token 0 is <cls>
            inputs["input_ids"][0, mask_idx] = tokenizer.mask_token_id
            with torch.no_grad():
                out = model(**inputs, output_hidden_states=True)
                log_probs = torch.log_softmax(out.logits[0, mask_idx], dim=-1)
                wt_id = tokenizer.convert_tokens_to_ids(wt)
                mt_id = tokenizer.convert_tokens_to_ids(mt)
                llr = float(log_probs[mt_id] - log_probs[wt_id])
                emb_norm = float(out.hidden_states[-1][0, mask_idx].norm().item())
        conn.execute(
            "INSERT OR REPLACE INTO esm2_scores (variant_id, llr, emb_norm) VALUES (?, ?, ?)",
            (v["id"], llr, emb_norm),
        )
        scored += 1
        if scored % 25 == 0:
            conn.commit()
    conn.commit()
    conn.close()
    stats = logger.log("esm2_scoring", {"scored_now": scored, "dry_run": dry_run})
    print("Run summary:", stats)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--max-wall-min", type=int, default=45)
    args = p.parse_args()
    demo_variants = [
        {"id": "BRCA1_M1V", "seq": "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCSENPRDTEDVPWITLNSSIQKVNEWFSRSDELLGSDDSHDGESESNAKVADVLDVLNEVDEYSGSSEKIDLLASDPHEALICKSERVHSKSVESNIEDKIFGKTYRKKASLPNLSHVTENLIIGAFVTEPQIIQERPLTNKLKRKRRPTSGLHPEDFIKKADLAVQKTPEMINQGTNQTEQNGQVMNITNSGHENKTKGDSIQNEKNPNPIESLEKESAFKTKAEPISSSISNMELELNIHNSKAPKKNRLRRKSSTRHIHALELVVSRNLSPPNCTELQIDSCSSSEEIKKKKYNQMPVRHSRNLQLMEGKEPATGAKKSNKPNEQTSKRHDSDTFPELKLTNAPGSFTKCSNTSELKEFVNPSLPREEKEEKLETVKVSNNAEDPKDLMLSGERVLQTERSVESSSISLVPGTDYGTQESISLLEVSTLGKAKTEPNKCVSQCAAFENPKGLIHGCSKDNRNDTEGFKYPLGHEVNHSRETSIEMEESELDAQYLQNTFKVSKRQSFAPFSNPGNAEEECATFSAHSGSLKKQSPKVTFECEQKEENQGKNESNIKPVQTVNITAGFPVVGQKDKPVDNAKCSIKGGSRFCLSSQFRGNETGLITPNKHGLLQNPYRIPPLFPIKSFVKTKCKKNLLEENFEEHSMSPEREMGNENIPSTVSTISRNNIRENVFKEASSSNINEVGSSTNEVGSSINEIGSSDENIQAELGRNRGPKLNAMLRLGVLQPEVYKQSLPGSNCKHPEIKKQEYEEVVQTVNTDFSPYLISDNLEQPMGSSHASQVCSETPDDLLDDGEIKEDTSFAENDIKESSAVFSKSVQKGELSRSPSPFTHTHLAQGYRRGAKKLESSEENLSSEDEELPCFQHLLFGKVNNIPSQSTRHSTVATECLSKNTEENLLSLKNSLNDCSNQVILAKASQEHHLSEETKCSASLFSSQCSELEDLTANTNTQDPFLIGSSKQMRHQSESQGVGLSDKELVSDDEERGTGLEENNQEEQSMDSNLGEAASGCESETSVSEDCSGLSSQSDILTTQQRDTMQHNLIKLQQEMAELEAVLEQHGSQPSNSYPSIISDSSALEDLRNPEQSTSEKAVLTSQKSSEYPISQNPEGLSADKFEVSADSSTSKNKEPGVERSSPSKCPSLDDRWYMHSCSGSLQNRNYPSQEELIKVVDVEEQQLEESGPHDLTETSYLPRQDLEGTPYLESGISLFSDDPESDPSEDRAPESARVGNIPSSTSALKVPQLKVAESAQSPAAAHTTDTAGYNAMEESVSREKPELTASTERVNKRMSMVVSGLTPEEFMLVYKFARKHHITLTNLITEETTHVVMKTDAEFVCERTLKYFLGIAGGKWVVSYFWVTQSIKERKMLNEHDFEVRGDVVNGRNHQGPKRARESQDRKIFRGLEICCYGPFTNMPTDQLEWMVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWTEDNGFHAIGQMCEAPVVTREWVLDSVALYQCQELDTYLIPQIPHSHY", "pos": 1, "wt": "M", "mt": "V"},
        {"id": "TP53_R175H", "seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD", "pos": 175, "wt": "R", "mt": "H"}
    ]
    score_variants(demo_variants, max_wall_min=args.max_wall_min, dry_run=args.dry_run)
