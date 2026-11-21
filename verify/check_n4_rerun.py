"""Audit the recorded full-size search against raw checkpoints and final geometry.

This validates the numerical evidence used by the paper. It is not a proof
of cell coverage, root completeness, or full-cell local optimality.
"""
import ast
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path

from n4_tighten2 import valid_fixed2
from n4_strata2 import actual_closure_ok
from n4_walls import ORDER, all_margins

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "verify/artifacts/n4-rerun-20261023"
FINAL = RUN / "final"


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Historical source digests are immutable. The runner has new resume handling.
# Only comments changed in these three modules. Compare their complete Python
# syntax trees, including docstrings and type comments, against each recorded copy.
COMMENTS_ONLY = {"n4_census.py", "n4_walls.py", "n4_cells3.py"}


def syntax(path):
    return ast.dump(ast.parse(path.read_text(encoding="utf-8"), type_comments=True))


def replay_boundary(record):
    attempt = record["attempt"]
    start, direction = record["start"], record["direction"]
    low, high = record["bracket"]

    def point(t):
        return [a + t*v for a, v in zip(start, direction, strict=True)]

    assert len(start) == len(direction) == 7 and abs(math.hypot(*direction) - 1) < 1e-15, attempt
    assert 0 < high - low <= 0.05 / 2**45 + math.ulp(3.0), (attempt, "bisection width")
    assert point(low) == record["theta"], (attempt, "stored endpoint")

    # The last bisection straddles a floating-point threshold. Different libm
    # implementations may put its last bit on either side. Test both sides on
    # the recorded unit ray without changing the original panel/GP tolerances.
    padding = 1e-12
    assert padding < low < high < 3, attempt
    assert valid_fixed2(point(low - padding), ORDER), (attempt, "inside of bracket")
    assert not valid_fixed2(point(high + padding), ORDER), (attempt, "outside of bracket")

    margins = {k: abs(v) for k, v in all_margins(record["theta"]).items() if abs(v) < 5e-4}
    saved = {label: distance for distance, label in record["near"]}
    assert len(saved) == len(record["near"]) and saved.keys() == margins.keys(), (attempt, "contact set")
    assert all(abs(margins[k] - saved[k]) < 1e-12 for k in saved), (attempt, "contact distances")
    recorded_label = record["near"][0][1]
    assert recorded_label == min(margins, key=margins.get), (attempt, "nearest contact")
    return recorded_label


if not __debug__:
    raise SystemExit("The recorded numerical audit requires Python assertions.")

for stage in ("sweep", "walls", "descent", "strata_filter", "cells", "codewalls", "stragglers"):
    directory = FINAL / stage
    meta = read(directory / "provenance.json")
    for name, expected in meta["sources"].items():
        archived = directory / "source" / name
        current = ROOT / "verify" / name
        assert digest(archived) == expected, (stage, name, "snapshot")
        if name in COMMENTS_ONLY:
            assert syntax(current) == syntax(archived), (stage, name, "current mathematical source")
        elif name != "n4_rerun.py":
            assert digest(current) == expected, (stage, name, "current source")

sweep = read(FINAL / "sweep/result.json")
streams = [read(FINAL / f"sweep/stream-{i}.json") for i in range(8)]
assert [(s["seed"], s["samples"]) for s in streams] == [(7000+i, 1_000_000) for i in range(8)]
assert sum(s["samples"] for s in streams) == sweep["total_samples"] == 8_000_000
assert [h for s in streams for h in s["hits"]] == sweep["hits"]
assert len(sweep["hits"]) == sweep["state_layer_records"] == 42
assert len({tuple(h["theta"]) for h in sweep["hits"]}) == sweep["distinct_states"] == 21
assert all(valid_fixed2(h["theta"], h["order"]) for h in sweep["hits"])

walls = read(FINAL / "walls/result.json")
assert walls["actual_attempts"] == walls["requested_attempts"] == 4000
assert len(walls["records"]) == walls["reached_boundaries"] == walls["classified"] == 3995
assert sum(walls["histogram"].values()) == 3995 and len(walls["histogram"]) == 8
assert len(walls["records"]) + sum(walls["skipped"].values()) == 4000
attempts = [r["attempt"] for r in walls["records"]]
assert attempts == sorted(set(attempts)) and all(0 <= i < 4000 for i in attempts)
assert Counter(replay_boundary(r) for r in walls["records"]) == walls["histogram"]

descent = read(FINAL / "descent/result.json")
assert descent["n_starts"] == descent["n_finished"] == len(descent["endpoints"]) == 42
for i, row in enumerate(descent["endpoints"]):
    assert row["index"] == i and row["seed"] == 9000+i
    assert row["start"] == sweep["hits"][i]["theta"]
    assert row["order"] == sweep["hits"][i]["order"]
    assert valid_fixed2(row["theta"], row["order"])
    assert abs(sum(row["theta"][4:]) - row["ell"]) < 1e-12
assert sum(abs(r["ell"]-descent["star_core"]) < .01 for r in descent["endpoints"]) == descent["within_0_01"] == 31

raw = [read(RUN / f"strata/stratum-{i}.json") for i in range(254)]
expected = {tuple(sorted(s)) for n in range(1, 8) for s in itertools.combinations(walls["histogram"], n)}
assert {tuple(sorted(r["labels"])) for r in raw} == expected
assert sum(len(r["attempts"]) for r in raw) == 1778
assert all(len(r["attempts"]) == 7 for r in raw)
filtered = read(FINAL / "strata_filter/result.json")
assert sum(len(r["solutions"]) for r in raw) == filtered["n_solutions"] == 56
filter_meta = read(FINAL / "strata_filter/provenance.json")["config"]
assert digest(RUN / "strata/result.json") == filter_meta["harvest_sha256"]
kept = [r for r in filtered["solutions"] if actual_closure_ok(r["theta"]) and min(r["multipliers"]) >= -1e-8]
assert len(kept) == filtered["n_pass"] == 14
assert max(r["star_distance_squared"]**.5 for r in kept) < 3e-10
assert filtered["n_other_pass"] == 0

cells = read(FINAL / "cells/result.json")
assert len(cells["cells"]) == 8 and cells["restarts"] == 480
assert sum(c["shadow_hits"] > 0 for c in cells["cells"]) == 4
assert all(c["valid_hits"] == 0 for c in cells["cells"])
code = read(FINAL / "codewalls/result.json")
assert code["n_targets"] == 15 and code["n_starts"] == 4
assert code["total_trials"] == 56270
rays = read(FINAL / "stragglers/result.json")
assert rays["n_endpoints"] == 42 and rays["n_far_stragglers"] == rays["far_sampled_rays"] == 11
print("PASS: final source snapshots, 8M sweep, 3995/4000 boundaries, 42 descents")
print("PASS: 254 x 7 raw Newton attempts, 56 records, 14 repeated candidate records")
print("PASS: 8 targeted cells, 15 margin targets, 11 sampled descent-ray diagnostics")
print("Numerical evidence only; no exhaustive or full-cell optimality claim")
