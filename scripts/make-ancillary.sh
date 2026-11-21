#!/bin/sh
# Assemble the manuscript and verification inputs, or check an existing archive.
set -eu
cd "$(dirname "$0")/.."
if [ "$#" -eq 0 ]; then
    python3 scripts/update-manifests.py
fi
exec python3 scripts/package-ancillary.py "$@"
