"""Build and compare the standalone manuscript and verification archive."""
import argparse
import gzip
import os
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "ancillary.tar.gz"
REQUIRED = (
    "README.md", "REPRODUCE.md", "LICENSE", "CITATION.cff", ".python-version",
    "pyproject.toml", "uv.lock", "paper/main.pdf", "paper/main.bbl", "paper/references.bib",
)
DISTRIBUTION_DIRS = ("scripts", "verify", "formal", "logs", "paper/figs")
SKIP = {".git", ".lake", ".venv", "__pycache__", ".DS_Store", "build"}


def release_files():
    files = {ROOT / name for name in REQUIRED}
    files.update((ROOT / "paper").glob("*.tex"))
    for directory in DISTRIBUTION_DIRS:
        for folder, dirs, names in os.walk(ROOT / directory):
            dirs[:] = [part for part in dirs if part not in SKIP]
            for name in names:
                if name not in SKIP and not name.endswith((".pyc", ".pyo")):
                    files.add(Path(folder) / name)
    return {path.relative_to(ROOT).as_posix(): path for path in sorted(files)}


def check_archive(files):
    seen = set()
    with tarfile.open(ARCHIVE, "r:gz") as archive:
        for member in archive:
            name = member.name
            if not member.isfile() or name in seen or name not in files:
                raise SystemExit(f"Unexpected archive member: {name}")
            if archive.extractfile(member).read() != files[name].read_bytes():
                raise SystemExit(f"Archive content differs from working file: {name}")
            if bool(member.mode & 0o111) != os.access(files[name], os.X_OK):
                raise SystemExit(f"Archive executable bit differs: {name}")
            seen.add(name)
    if missing := files.keys() - seen:
        raise SystemExit(f"Archive is missing: {', '.join(sorted(missing))}")
    print(f"PASS: {len(seen)} distribution files match the working tree")


def build_archive(files):
    with tempfile.TemporaryDirectory(dir=ROOT, prefix=".ancillary-") as tmp:
        output = Path(tmp) / ARCHIVE.name
        with output.open("wb") as target, gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=target, mtime=0) as stream:
            with tarfile.open(mode="w|", fileobj=stream, format=tarfile.PAX_FORMAT) as archive:
                for name, path in files.items():
                    member = tarfile.TarInfo(name)
                    member.size = path.stat().st_size
                    member.mode = 0o755 if os.access(path, os.X_OK) else 0o644
                    with path.open("rb") as source:
                        archive.addfile(member, source)
        output.replace(ARCHIVE)
    print(f"Wrote {ARCHIVE.name} with {len(files)} files")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare with the archive without writing")
    args = parser.parse_args()
    files = release_files()
    if args.check:
        check_archive(files)
    else:
        build_archive(files)
