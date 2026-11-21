"""Export the loose example used by the historical figure-eight illustration."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path, help="new JSON file for the loose state")
args = parser.parse_args()
states = json.loads(Path(__file__).with_name("n4_valid41.json").read_text())
theta, order = states[1]
with args.output.open("x", encoding="utf-8") as target:
    json.dump({"theta": theta, "order": order}, target)
print(f"Wrote {args.output}")
