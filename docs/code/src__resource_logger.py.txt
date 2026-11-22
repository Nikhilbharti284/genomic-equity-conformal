import hashlib, json, os, time, psutil

class ResourceLogger:
    def __init__(self, log_path="artifacts/resource_log.jsonl", cores=2):
        self.log_path = log_path
        self.cores = cores
        self.t0 = time.time()
        self.proc = psutil.Process(os.getpid())

    def sha256_file(self, path):
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()

    def log(self, stage, extra=None):
        wall_sec = round(time.time() - self.t0, 2)
        entry = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "stage": stage,
            "wall_sec": wall_sec,
            "core_hours_used": round(wall_sec * self.cores / 3600.0, 5),
            "peak_rss_mb": round(self.proc.memory_info().rss / (1024 * 1024), 2),
        }
        if extra:
            entry.update(extra)
        os.makedirs(os.path.dirname(self.log_path), exist_ok=True)
        with open(self.log_path, "a") as f:
            f.write(json.dumps(entry) + "\n")
        return entry
