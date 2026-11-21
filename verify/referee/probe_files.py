"""File locations for the historical figure-eight probes."""
import json
import os
from pathlib import Path

VERIFY = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = VERIFY.parent / "build" / "figure-eight-probes"


def work_file(name):
    return Path(os.environ.get("PENTAGON_PROBE_DIR", DEFAULT_OUTPUT)).resolve() / name


def read_probe(name):
    path = work_file(name)
    if not path.is_file():
        producer = {"ref_hits.json": "ref5_census.py", "ref10_endpoints.json": "ref10_stragglers.py"}[name]
        raise SystemExit(f"Missing {path}. Generate it with python verify/referee/{producer}.")
    return json.loads(path.read_text(encoding="utf-8"))


def write_probe(name, rows):
    path = work_file(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} records to {path}")
