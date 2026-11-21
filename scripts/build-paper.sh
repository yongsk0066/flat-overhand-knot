#!/bin/sh
# Rebuild the distributed manuscript with the publication's TeX Live image.
set -eu
cd "$(dirname "$0")/.."
tex_image=ghcr.io/xu-cheng/texlive-alpine@sha256:c06416ebbc69c8e90e07d1f1a6f7a3f9b334e773daa3543ce879ea533317c3a7
docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -v "$PWD/paper:/paper" -w /paper "$tex_image" sh -ec '
        mkdir -p /tmp/manuscript
        latexmk -pdf -file-line-error -interaction=nonstopmode -halt-on-error \
            -outdir=/tmp/manuscript main.tex
        cp /tmp/manuscript/main.pdf /tmp/manuscript/main.bbl .
    '
