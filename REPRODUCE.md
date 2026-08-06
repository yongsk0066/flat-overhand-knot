# Reproducibility guide

Every numerical, exact-algebraic, and machine-checked claim in the paper can be
re-run from this repository. All commands are executed from the repository root;
the Python environment is managed by [uv](https://docs.astral.sh/uv/) from
`pyproject.toml` + `uv.lock` (run scripts with `uv run python ...`).

## Environment

| Component | Pinned version |
|---|---|
| Python | 3.12+ (developed on 3.14; managed by uv) |
| sympy / mpmath | 1.14.0 / 1.3.0 |
| numpy / scipy / shapely | 2.5.1 / 1.18.0 / 2.1.2 |
| python-flint (Arb) | 0.9.0 |
| Lean 4 | 4.33.0-rc1 (elan; Lake 5.0.0), Mathlib pinned in `formal/lake-manifest.json` |
| Figures | `rsvg-convert` (librsvg) |
| Paper | TeX Live (`pdflatex` + `bibtex`) |

```bash
uv sync            # creates .venv from the committed uv.lock
# pip fallback:  python3 -m venv .venv && .venv/bin/pip install -r verify/requirements.txt
```

## Quick verification (~3 minutes)

```bash
uv run python verify/n4_exact_valid.py  # all 39 contact margins of the figure-eight tight state, exactly in Q(sqrt 15)
uv run python verify/exact_taco.py      # taco entanglement identity, endpoint exclusions, family certificates
```

## Trefoil suite

**Model and slice.** `check_slice.py`, `check_slice3.py`, `check_topology.py`, `exact_slice.py`.

**Walls and the pentagonal junction.** `check_layers.py`, `check_final.py`,
`exact_farkas.py`, `exact_farkas2.py`, `check_taco.py`, `exact_taco.py` (~2 min).

**The tight state θ†.** `check_edge.py`, `family_trace.py`, `check_branch.py`,
`exact_elim.py` (~1 min), `certify.py` (1–2 min), `exact_kkt_dagger.py`,
`exact_second.py`, `exact_opening.py`, `exact_notch.py`, `exact_outline.py` (~1 min),
`exact_reach.py` (~1 min; monotone feasible descent from the pentagon),
`exact_branch.py` (~15 min; branch preservation along the second family),
`exact_unique_h.py` (uniqueness of the critical height over x†).

**Truncation length, realization, fold number.** `check_kauffman_compare.py`,
`exact_kauffman.py`, `exact_kauffman2.py` (~2 min), `check_trunc.py`,
`exact_trunc_nbhd.py`, `check_realize.py`, `check_global.py`, `check_gauss.py`, `check_gp.py`.

**Blocked-state census** (order matters):

```bash
uv run python verify/uniq_setup.py      # rational wall system, cross-validated
uv run python verify/uniq_poly.py       # polynomialization (0/1600 mismatches)
uv run python verify/uniq_strata.py     # KKT harvest, implementation 1 (215 strata x 30 seeds, ~1 min on 8 cores)
uv run python verify/uniq_census.py     # harvest, implementation 2 (100 seeds, ~4 min)
uv run python verify/uniq_certify2.py   # exclusion certificates (tens of minutes)
```

Expected outcome: every script passes its assertions; θ† is the unique candidate
passing all filters; the pentagon appears only with a negative multiplier; all 72
regular candidates are excluded by interval certificates.

**Independent re-verification in Arb.** `arb_recheck.py` (~15 min): the root
bracket of p_x, the reduced Hessian, the KKT multipliers, and the 67 sliver-sign
census exclusions, re-evaluated in ball arithmetic independently of mpmath.

## Figure-eight suite

```bash
uv run python verify/n4_tighten2.py --test      # Fox 5-coloring oracle unit tests
uv run python verify/n4_cells3.py               # cell enumeration by exact determinants
uv run python verify/n4_tighten2.py --harvest   # 3M-sample harvest (~2 min on 8 cores)
uv run python verify/n4_tighten2.py --tighten2  # exact-axiom tightening
uv run python verify/n4_reduce2.py              # closed form via 60-digit Newton + integer relations
uv run python verify/n4_certify.py              # exact certification in Q(sqrt 15), incl. the hexagon tiling 24 = 24
uv run python verify/n4_exact_kkt.py            # exact KKT multipliers, L_trunc = 32 sqrt(15)/15, opening direction
uv run python verify/n4_exact_valid.py          # all 39 candidate contact margins, exactly
uv run python verify/n4_census.py               # 8M-sample sweep + 42 descents (~12 min on 8 cores)
uv run python verify/n4_walls.py                # wall classification, 3995 boundary probes (~6 min)
uv run python verify/n4_codewall_probe.py       # code/length wall margins >= 0.525 (~4 min)
uv run python verify/n4_chart.py                # tangent half-angle chart; 8 wall polynomials vanish at θ*
uv run python verify/n4_krawczyk.py             # existence/uniqueness enclosure (radius 1e-25)
uv run python verify/n4_strata2.py              # 254-strata KKT harvest (~8 min on 8 cores)
uv run python verify/n4_cells_probe.py          # exclusion probes for the 8 abstract cells (~15 min)
uv run python verify/n4_straggler.py --stalls   # Farkas rejection of descent stalls
```

Expected outcome: cos τ₁ = −1/4 exactly, L* = (6,8,6)/√15, core length 4√15/3,
truncation 32√15/15 with ratio 8/5; KKT residual identically zero; the Krawczyk
box isolates θ*; all strata survivors coincide with θ*.

## Implementation discipline

- `sympy.lambdify` is always called with `cse=True` (raw-tree evaluation is ~570× slower).
- Expanded polynomials are never evaluated in floating point (catastrophic cancellation);
  numerics use structured forms, exact algebra uses expanded forms at ≥40 digits.
- Census exclusion certificates prefer the unconditionally scaled walls
  (w1, w4, w6, w7, w8); conditionally scaled walls (w5, w9, w10) are used only
  with a co-certified scale factor.

## Historical-record scripts

`check_walls.py`, `check_descent.py`, `n4_tighten.py`, `n4_cells2.py`, and
`n4_kkt.py` are retained deliberately: they reproduce defects that were found and
repaired during the internal adversarial-referee cycles described in the paper's
methods appendix. Their outputs are superseded by the scripts listed above.
`verify/referee/` contains the sixteen adversarial-referee scripts themselves.

## Certificate logs and hashes

`logs/` contains the recorded output of every exact and interval certificate.
Integrity check:

```bash
sha256sum -c logs/hashes-verify.txt && sha256sum -c logs/hashes-logs.txt
```

## Figures (SVG → PDF)

```bash
uv run python verify/tref_render.py        # fig1-pentagon.svg, theta-dagger-heptagon.svg
uv run python verify/n4_pair.py            # fig41-loose.svg, fig41-tight.svg
uv run python verify/make_fig_notation.py  # fig0-notation.svg
for f in fig0-notation fig1-pentagon fig3-mountain-pass fig4-two-families fig5-slab-ramp fig41-code fig41-loose fig41-tight; do
  rsvg-convert -f pdf paper/figs/$f.svg -o paper/figs/$f.pdf; done
rsvg-convert -f pdf paper/figs/theta-dagger-heptagon.svg -o paper/figs/fig2-heptagon.pdf
```

## Building the paper

```bash
cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Lean formalization

```bash
cd formal && lake exe cache get && lake build Formal.Pentagon Formal.ThetaDagger
```

## Packaging the archived version

`scripts/make-ancillary.sh` regenerates the hash manifests and produces
`ancillary.tar.gz`, the archive deposited alongside the paper.
