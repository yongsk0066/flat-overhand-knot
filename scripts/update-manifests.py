"""Refresh integrity manifests, including nested numerical-run evidence."""
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def eligible(path):
    return not any(part in {".venv", "__pycache__", ".git"} for part in path.parts)


groups = {
    "hashes-verify.txt": [
        *ROOT.joinpath("verify").rglob("*.py"),
        *ROOT.joinpath("scripts").glob("*.py"),
        *ROOT.joinpath("scripts").glob("*.sh"),
    ],
    "hashes-logs.txt": [
        *ROOT.joinpath("logs").rglob("*.log"),
        *ROOT.joinpath("verify/artifacts").rglob("*.log"),
    ],
    "hashes-artifacts.txt": [
        ROOT / "paper/main.pdf", ROOT / "paper/main.bbl",
        *ROOT.joinpath("verify").rglob("*.json"),
        *ROOT.joinpath("verify").rglob("*.jsonl"),
        *ROOT.joinpath("paper/figs").glob("*.svg"),
        *ROOT.joinpath("paper/figs").glob("*.pdf"),
        *ROOT.joinpath("verify/artifacts/figure-history").glob("*.svg"),
        *ROOT.joinpath("verify/artifacts/figure-history").glob("*.pdf"),
    ],
}

for name, candidates in groups.items():
    paths = sorted({p for p in candidates if p.is_file() and eligible(p)})
    lines = [f"{sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}\n"
             for p in paths]
    ROOT.joinpath("logs", name).write_text("".join(lines), encoding="utf-8")
    print(f"{name}: {len(paths)} files")
