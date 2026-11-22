import json, os

class PGxDrugSafetySidecar:
    """Zero-GPU, CPIC-aligned Global Pharmacogenomics (PGx) Drug-Toxicity Alert Sidecar (v2.0.0)."""
    def __init__(self, panel_path="data/pgx_cpic_panel.json"):
        if not os.path.exists(panel_path):
            panel_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "pgx_cpic_panel.json")
        with open(panel_path, encoding="utf-8") as f:
            self.panel = json.load(f)

    def lookup(self, variant_id=None, rsid=None):
        if variant_id:
            hit = self.panel["variants"].get(variant_id.upper())
            if hit:
                return {"matched_by": "variant_id", "key": variant_id.upper(), **hit}
        if rsid and str(rsid).lower().startswith("rs"):
            hit = self.panel["rsid_tags"].get(str(rsid).lower())
            if hit:
                return {"matched_by": "rsid_tag", "key": str(rsid).lower(), **hit}
            for vid, meta in self.panel["variants"].items():
                if meta.get("rsid", "").lower() == str(rsid).lower():
                    return {"matched_by": "rsid", "key": vid, **meta}
        return None
