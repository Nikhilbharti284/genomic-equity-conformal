"""
Genomic Equity Conformal v4.0.0 — Tamper-Evident SHA-256 Hash-Chain Audit Ledger
Zero-dependency compliance ledger for ISO 15189 / CAP / EU AI Act / HIPAA / GDPR / DPDP.
"""
import hashlib
import json
import os
import time

GENESIS_HASH = "0" * 64


def _canonical(obj: dict) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")


def append_audit_entry(event_type: str, payload_summary: dict, ledger_path: str = "artifacts/v4_audit_chain.jsonl") -> dict:
    os.makedirs(os.path.dirname(ledger_path) or ".", exist_ok=True)
    prev_hash = GENESIS_HASH
    seq = 1
    if os.path.exists(ledger_path):
        with open(ledger_path, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if lines:
            last = json.loads(lines[-1])
            prev_hash = last["entry_hash"]
            seq = int(last["seq"]) + 1

    body = {
        "seq": seq,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "event_type": event_type,
        "payload_sha256": hashlib.sha256(_canonical(payload_summary)).hexdigest(),
        "summary": payload_summary,
        "prev_hash": prev_hash,
    }
    entry_hash = hashlib.sha256(_canonical(body)).hexdigest()
    record = {**body, "entry_hash": entry_hash}
    with open(ledger_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True) + "\n")
    return record


def verify_audit_chain(ledger_path: str = "artifacts/v4_audit_chain.jsonl") -> dict:
    if not os.path.exists(ledger_path):
        return {"valid": True, "entries": 0, "reason": "Empty ledger"}
    prev = GENESIS_HASH
    count = 0
    with open(ledger_path, "r", encoding="utf-8") as f:
        for idx, raw in enumerate(f, 1):
            raw = raw.strip()
            if not raw:
                continue
            rec = json.loads(raw)
            expected_prev = rec.get("prev_hash")
            given_hash = rec.get("entry_hash")
            body = {k: v for k, v in rec.items() if k != "entry_hash"}
            recomputed = hashlib.sha256(_canonical(body)).hexdigest()
            if expected_prev != prev or given_hash != recomputed:
                return {"valid": False, "entries": count, "broken_at_seq": idx}
            prev = given_hash
            count += 1
    return {"valid": True, "entries": count, "head_hash": prev}


if __name__ == "__main__":
    test_log = "artifacts/v4_audit_chain_selftest.jsonl"
    if os.path.exists(test_log):
        os.remove(test_log)
    append_audit_entry("VCF_AUDIT", {"sample": "S1", "n_variants": 3}, test_log)
    append_audit_entry("FHIR_EXPORT", {"sample": "S1", "bundle": "ok"}, test_log)
    print("Audit Chain Verify:", verify_audit_chain(test_log))
    os.remove(test_log)
