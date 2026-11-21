"""Run each original certificate and its mathematical mutation in both Python modes."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

VERIFY_DIR = Path(__file__).resolve().parent


def _require(cond, msg="guard check failed"):
    if not cond:
        raise RuntimeError(msg)


MUTATIONS = [
    (
        "exact_global_core.py",
        "-6*s*q**2+(s**4-4*s**2-1)*q",
        "-5*s*q**2+(s**4-4*s**2-1)*q",
        "wrong global length envelope",
        "reduced length polynomial",
    ),
    (
        "n4_exact_cell.py",
        "(4*v_value-2*u_value)/root",
        "(4*v_value+2*u_value)/root",
        "wrong coupled cap orientation",
        "upper tail cap",
    ),
    (
        "exact_slice.py",
        "cos108 = -(sqrt(5) - 1) / 4",
        "cos108 = 0",
        "wrong pentagon interior angle",
        "One or more standard-slice exact identities failed",
    ),
    (
        "exact_unique_h.py",
        "hf = sp.Rational(1, 12)",
        "hf = 1 + sp.Rational(1, 12)",
        "perturbed h-dagger polynomial",
        "h-linear relation failed mod p_x",
    ),
    (
        "exact_trunc_nbhd.py",
        '_require(red(uu - sp.Rational(1, 2)) == 0, f"{tag} on {nm}: u != 1/2")',
        '_require(red(uu - sp.Rational(1, 3)) == 0, f"{tag} on {nm}: u != 1/3")',
        "wrong truncation-cut midpoint",
        "outgoing D on EH: u != 1/3",
    ),
    (
        "exact_kauffman2.py",
        "field.convert(28) + field.convert(12) * S5",
        "field.convert(29) + field.convert(12) * S5",
        "corrupted Farkas multiplier lambda_w4",
        "Farkas residual nonzero",
    ),
    (
        "exact_kkt_dagger.py",
        "R(17, 96) * x ** 11",
        "R(19, 96) * x ** 11",
        "corrupted KKT multiplier lambda_w6",
        "KKT stationarity residual nonzero",
    ),
    (
        "exact_opening.py",
        "dstar = (R(-93, 1000), R(2197, 1000), R(1021, 500), R(509, 250), R(419, 200))",
        "dstar = (R(93, 1000), R(-2197, 1000), R(-1021, 500), R(-509, 250), R(-419, 200))",
        "negated opening direction at P*",
        "Pentagon opening direction has the wrong core slope",
    ),
    (
        "exact_reach.py",
        "Aell = sp.expand(",
        "Aell = -sp.expand(",
        "reversed core-derivative numerator",
        "Aell: wrong interior sign",
    ),
    (
        "exact_branch.py",
        'COND["AC:s>0"]=sAC;',
        'COND["AC:s>0"]=sAC*(1000*x-323)**2;',
        "injected interior root at x=0.323 in branch preservation",
        "One or more branch-preservation resultants have a root",
    ),
    (
        "exact_census_roots.py",
        'threshold = arb("-0.315")',
        'threshold = arb("-100.0")',
        "impossible geometric sliver threshold in census root certification",
        "Saved centre row 0",
    ),
    (
        "n4_exact_kkt.py",
        "(21, 3, 21, 13, 4, 17, 5)",
        "(23, 3, 21, 13, 4, 17, 5)",
        "corrupted figure-eight KKT multiplier",
        "Multiplier mismatch for g1",
    ),
]


def run_guards():
    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)
    env.update(PYTHONPATH=str(VERIFY_DIR), PYTHONDONTWRITEBYTECODE="1")
    with tempfile.TemporaryDirectory(prefix="pentagon-guard-") as tmp:
        for filename, old, new, label, diagnostic in MUTATIONS:
            source = (VERIFY_DIR / filename).read_text(encoding="utf-8")
            _require(source.count(old) == 1, f"Expected one mutation site in {filename}: {old}")
            path = Path(tmp) / filename
            for changed, code in ((False, source), (True, source.replace(old, new, 1))):
                path.write_text(code, encoding="utf-8")
                for flags in ([], ["-O"]):
                    mode = " ".join(flags) or "default"
                    run = subprocess.run([sys.executable, *flags, str(path)], cwd=VERIFY_DIR.parent,
                                         env=env, capture_output=True, text=True, timeout=900)
                    if changed:
                        last = run.stderr.rstrip().splitlines()[-1:] or [""]
                        passed = run.returncode == 1 and last[0].startswith("RuntimeError: " + diagnostic)
                    else:
                        passed = run.returncode == 0
                    _require(passed, f"{filename} ({label}, {mode}, mutated={changed}): "
                             f"exit {run.returncode}\n{run.stdout}\n{run.stderr}")
            print(f"[guard] {filename}: original passed and {label} rejected in both modes", flush=True)
    print(f"PASS: all {len(MUTATIONS)} mathematical mutations rejected in default and optimized Python")


if __name__ == "__main__":
    run_guards()
