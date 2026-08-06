#!/bin/sh
# Build ancillary.tar.gz — the self-contained archive deposited alongside the paper.
# Regenerates the SHA-256 manifests first; excludes large build artifacts
# (verify/.venv, formal/.lake, __pycache__).
set -e
cd "$(dirname "$0")/.."
shasum -a 256 verify/*.py verify/referee/*.py > logs/hashes-verify.txt
shasum -a 256 logs/*.log > logs/hashes-logs.txt
tar --exclude='verify/.venv' --exclude='formal/.lake' --exclude='__pycache__' \
    -czf ancillary.tar.gz \
    README.md REPRODUCE.md LICENSE CITATION.cff \
    verify formal logs paper/*.tex paper/references.bib paper/figs
echo "ancillary.tar.gz: $(du -h ancillary.tar.gz | cut -f1)"
