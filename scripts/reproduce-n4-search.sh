#!/bin/sh
# Full-size numerical reproduction. Use a new directory for a fresh run.
set -eu
cd "$(dirname "$0")/.."
n4_run_dir=${1:-build/n4-reproduction}
n4_baseline=verify/artifacts/n4-rerun-20261023/baseline/verify
python3 - "$n4_run_dir" "$n4_baseline" <<'CHECK_OUTPUT'
from pathlib import Path
import sys

output, baseline = map(Path, sys.argv[1:])
if output.resolve().is_relative_to(baseline.resolve().parents[1]):
    raise SystemExit("Choose a new output directory outside the archived run.")
CHECK_OUTPUT
mkdir -p "$n4_run_dir/baseline/verify"
for n4_input in n4_valid41.json n4_census.json n4_wall_inventory.json; do
  if [ -e "$n4_run_dir/baseline/verify/$n4_input" ]; then
    cmp "$n4_baseline/$n4_input" "$n4_run_dir/baseline/verify/$n4_input"
  else
    cp -p "$n4_baseline/$n4_input" "$n4_run_dir/baseline/verify/$n4_input"
  fi
done
for n4_stage in sweep walls descent strata cells; do
  uv run python -u verify/n4_rerun.py --stage "$n4_stage" --workers 4 --output "$n4_run_dir" > "$n4_run_dir/$n4_stage.log" 2>&1
  tail -n 2 "$n4_run_dir/$n4_stage.log"
done
for n4_stage in codewalls stragglers; do
  uv run python -u verify/n4_rerun_diagnostics.py --stage "$n4_stage" --output "$n4_run_dir" > "$n4_run_dir/$n4_stage.log" 2>&1
  tail -n 2 "$n4_run_dir/$n4_stage.log"
done
printf 'Numerical reproduction complete: %s\n' "$n4_run_dir"
