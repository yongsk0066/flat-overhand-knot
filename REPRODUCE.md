# Reproducibility guide

The exact certificates and the recorded numerical search procedures can be
re-run from this repository. Numerical searches are not completeness proofs. All commands are executed from the repository root;
the Python environment is managed by [uv](https://docs.astral.sh/uv/) from
`pyproject.toml` + `uv.lock` (run scripts with `uv run python ...`).

## Environment

| Component | Pinned version |
|---|---|
| Python | 3.13.18 (pinned by `.python-version`; managed by uv) |
| uv | 0.13.9 (tested with managed Python and used in CI) |
| sympy / mpmath | 1.14.0 / 1.3.0 |
| numpy / scipy / shapely | 2.5.5 / 1.18.1 / 2.2.0 |
| python-flint (Arb) | 0.9.0 |
| z3-solver | 5.4.1.0 |
| Lean 4 | 4.36.0 (elan; Lake 5.0.0), Mathlib pinned in `formal/lake-manifest.json` |
| Figures | `rsvg-convert` (librsvg) |
| Paper | TeX Live (`pdflatex` + `bibtex`) |

```bash
uv sync --locked   # creates .venv from the committed uv.lock
# pip fallback:  python3 -m venv .venv && .venv/bin/pip install -r verify/requirements.txt
```

## Current regression certificates

```bash
sh scripts/verify-core.sh
```

A successful run ends with `PASS: core certificates, mutation controls and figure checks`.
All Python filenames below are in [verify/](verify/). The core command runs the
entries marked **core**, then checks the shared geometry export and the paper's SVG
generators. The Python CI workflow also compiles the scripts and checks the three
integrity manifests. CI also rebuilds the archive and compares its bytes with the published file.

### Theorem-to-certificate guide

| Theorem and claim | Command or Python entry point | Method and expected result |
|---|---|---|
| A · Three-crease crossing rigidity and layer order | `check_gauss.py`, `check_gp.py` (extended), `exact_common_overlap.py` (core) | Diagram enumeration gives the alternating trefoil. The global order is B > D > A > C. The rational loose example realizes cyclic local orders. |
| B · Slice, walls and chart | `exact_slice.py`, `exact_core_descent.py`, `exact_common_overlap.py` (core), `exact_taco.py` (extended) | Exact slice identities and the affine core gradient (1/2 in all five coordinates). Four-panel overlap fixes the order along the specified slices. |
| C · Junction and reachability | `exact_core_descent.py`, `exact_reach.py`, `exact_branch.py` (core) | The inward ray has core slope −1/20. Five root-free resultants, endpoint identities and interval signs certify monotonicity. Twenty-two further conditions preserve the crossing branch. |
| D · Algebraic state and notched outline | `exact_field.py`, `exact_unique_h.py`, `exact_kkt_dagger.py`, `exact_second.py`, `exact_outline.py`, `exact_notch.py` (core) | Degree 12, associated sextic group S₆, four positive multipliers, positive reduced curvature and a seven-edge outline. |
| E · Kauffman truncation | `exact_kauffman2.py`, `exact_trunc_nbhd.py` (core) | Exact Farkas and flush-edge certificates give a strict local minimum at the regular pentagon. |
| F · Realization on the globally ordered branch | `exact_opening.py`, `exact_reach_geometry.py`, `exact_common_overlap.py` (core), `check_realize.py` (numerical illustrations) | Fixed-angle dilation gives positive panel and crease clearances at every valid state in the chart for small positive dilation. Open states have strong slab realizations, and zero-clearance states have limit realizations. Explicit openings at the two optima and the loose cyclic-order example have separate certificates. |
| G · Unique global core minimum in the positive chart | `exact_global_core.py` (core) | Necessary panel inequalities, strict averaging, a supporting double root in the dagger field and an independent discriminant/Sturm lower bound. |
| G · Reachability from every slack valid state in the positive chart | `exact_reach_geometry.py`, `exact_reach_transport.py`, `check_reach_paths.py` (core) | Midpoint geometry verifies the six contact functions and crossing order. Finite-step identities and 33 unsatisfiable exact real counterexample queries check transport and the algebraic endpoint. A separate rational clipping check samples 4,310 panel configurations across the three balancing cases and both label orientations. |
| G · Saved census candidates | `exact_census_roots.py` (core), `recertify_census_pending.py` (four-row extended check) | 104 distinct saved candidates, 72 isolating boxes each violating a true geometric wall below −0.315, and 32 unconfirmed centres classified 4 / 6 / 18 / 4. |
| H · Figure-eight state | `n4_exact_valid.py`, `n4_exact_opening.py`, `n4_exact_cell.py`, `n4_cell_crosscheck.py`, `n4_exact_kkt.py`, `n4_exact_local.py`, `n4_krawczyk.py`, `n4_certify.py` (core) | Valid oblique panels, an opening and a strict local minimum in the geometric cell, including both double-corner contacts. Core length 4√15/3 and truncation 32√15/15. |
| H · Recorded search | `check_n4_geometry.py`, `check_n4_rerun.py` (core), `n4_cells3.py` (extended) | Checkpoints, source records, finite search counts and determinant-based enumeration. |
| Regression controls | `test_guards.py` (core) | Fourteen original certificate scripts pass, and fourteen specific mathematical faults raise their expected diagnostic in both default and optimized Python. |

The [Lean guide](formal/README.md) identifies the 16 algebraic statements. Near the figure-eight hexagon, the eight explicit contacts are equivalent to nonpenetration against the actual crease caps. The global classification of figure-eight cells remains open.

## Further trefoil calculations

- `exact_elim.py` derives the field equations. `certify.py` isolates the stationary point,
  and `arb_recheck.py` independently checks the root bracket, Hessian, multipliers
  and 71 saved wall-verdict rows using Arb ball arithmetic. `family_trace.py`
  and `check_branch.py` are numerical traces of the same family.
- `check_slice.py`, `check_slice3.py`, `check_topology.py`, `check_layers.py` and
  `check_final.py` illustrate slice and wall geometry. `check_taco.py` and
  `check_trunc.py` are numerical probes. `exact_taco.py` and
  `exact_trunc_nbhd.py` provide the corresponding algebraic calculations.

### Census inputs and results

Read the saved harvest (`verify/uniq_census.json`) and dispositions
(`verify/uniq_certified.json`) in place with:

```sh
uv run python verify/exact_census_roots.py
uv run python verify/recertify_census_pending.py
```

The latter checks four specified W7 exclusions. Its optional `--write` flag
updates their verdicts in the local copy. The 32 inconclusive centres have
4 cases with `x < 0`, 6 with `|q + x| < 10⁻⁶`, 18 with `|1 − qx| < 10⁻⁶`,
and 4 further cases with negative `h` and `a`. All 104 saved centres have a
radius-10⁻¹² box violating a true sliver wall. The 72 regular roots have
separate certified isolating boxes. Root existence for the other 32 centres and exhaustiveness of the multistart search remain open. The length dilation in the global-core proof establishes MFCQ at every valid point of the positive trefoil chart.

To **regenerate** the harvest, extract the ancillary archive into a new directory,
change into that directory, and run the following commands there. They write
JSON results under that copy's `verify/` directory.

```sh
uv sync --locked
uv run python verify/uniq_setup.py
uv run python verify/uniq_poly.py
uv run python verify/uniq_strata.py
uv run python verify/uniq_census.py
uv run python verify/uniq_certify2.py
```

The two harvests cover 215 strata using their documented start sets.
`uniq_poly.py` compares 1,600 wall evaluations. The root-box checks in the
core suite refer to the saved candidates in this distribution.

## Figure-eight suite

The full-size numerical protocol was rerun on 2026-10-23 with the actual
oblique crease caps, unbounded tails, layer axiom and centreline GP checks.
The archived search retains its original environment records. The current Python lock reproduces the eight seeded streams, all 42 descent
endpoints, and the coordinates, multipliers and residuals of the 56 recorded
Newton solutions. Wall-clock fields depend on the machine.
The final evidence is in `verify/artifacts/n4-rerun-20261023/final/`:

| Stage | Recorded result | Status |
|---|---|---|
| Blind sweep | 8,000,000 configurations; 21 distinct accepted states, 42 state/layer records | numerical |
| Boundary probes | 4,000 attempts; 3,995 boundaries, all assigned to 8 contact labels; 5 invalid starts | numerical |
| Descent | All 42 starts completed; 31 endpoints within 0.01 of the exact core value | numerical proximity, not a convergence proof |
| KKT harvest | 254 subsets × 7 seeds = 1,778 attempts; 56 convergence records, 14 nonnegative-multiplier records at the same candidate | numerical; repeated roots are not distinct solutions |
| Other abstract cells | 8 cells × 60 restarts; four with shadow hits, none with valid-state hits | finite search, not an emptiness proof |
| Code/length margins | 15 targets × 4 starts, 56,270 trials; observed minimum 0.5252598788 | not a global lower bound |
| Remaining descent endpoints | Sampled decreasing rays at all 11 endpoints outside the 0.01 threshold | 101 checks per ray, not a continuous-interval certificate |

Each boundary record contains the starting point, a unit direction, the final
45-step bisection bracket and its contact distances. At double precision the last
bit of a bracket endpoint can depend on the platform's trigonometric library.
The replay reconstructs the saved endpoint, verifies validity on the inner side
and invalidity on the outer side within a parameter padding of 10⁻¹² along that
ray, then recomputes the near contacts and their histogram. The geometric
distance, crease-parameter and general-position tolerances are those of the
recorded search.

Each stage retains input hashes, source snapshots, tolerances and checkpoints.
The top-level `sweep/`, `walls/` and `descent/` are explicitly marked interim
runs from before the endpoint-GP correction. The geometry-independent
`strata/` contains all 1,778 raw Newton attempts; `final/strata_filter/` applies
the final predicate to every raw convergence record. Original historical JSON
files are preserved in `baseline/` and are not relabelled as new results.

Audit the stored evidence without repeating the large searches:

```bash
uv run python verify/check_n4_rerun.py
```

Repeat the entire numerical protocol in a **new** output directory:

```bash
sh scripts/reproduce-n4-search.sh build/n4-reproduction
```

This copies the archived starting inputs, uses the original eight seeds
7000–7007 with one million samples each, runs all 42 descents and all 254 × 7
KKT attempts, and then runs the supplementary cell and margin probes. It uses
four workers. Sweep, descent and KKT checkpoints can be resumed with unchanged
inputs and settings. Existing receipts reject source/configuration mismatches;
use a new directory after changing the inputs. The final recorded sweep used
two workers and the final descent used four; the random streams are independent
of worker scheduling.

The pointwise and local mathematical certificates are separate:

```bash
uv run python verify/n4_tighten2.py --test      # Fox colouring truth-table checks
uv run python verify/n4_cells3.py              # exact combinatorial enumeration
uv run python verify/n4_reduce2.py             # algebraic candidate identification
uv run python verify/n4_certify.py             # exact hexagon union, with overlap
uv run python verify/n4_exact_kkt.py            # exact multipliers and opening direction
uv run python verify/n4_exact_valid.py          # all actual-panel forbidden pairs
uv run python verify/n4_exact_opening.py        # continuous valid opening interval
uv run python verify/n4_exact_local.py          # first-order growth on the contact cone
uv run python verify/n4_exact_cell.py           # actual two-corner cap and side conditions
uv run python verify/n4_cell_crosscheck.py      # independent rational cap and KKT data
uv run python verify/n4_chart.py               # rational half-angle contact equations
uv run python verify/n4_krawczyk.py             # existence/uniqueness in the root box
```

Expected exact values: cos τ₁ = −1/4, L* = (6,8,6)/√15, core length 4√15/3,
truncation 32√15/15 with ratio 8/5. The geometric-cell proof combines the two-corner cap data with the independent rank and positive-multiplier certificate. Exhaustive root and cell coverage remains open.

## Implementation discipline

- The census KKT evaluators use `sympy.lambdify(..., cse=True)` to reuse common subexpressions in the multistart loop. Smaller symbolic checks also use the default evaluator.
- The numerical census uses structured expressions to limit cancellation. Quotient-field identities are decided by exact reduction.
- Census exclusion certificates prefer the unconditionally scaled walls
  (w1, w4, w6, w7, w8); conditionally scaled walls (w5, w9, w10) are used only
  with a co-certified scale factor.

## Historical-record scripts

The following names remain available to reproduce earlier calculations.

| Historical entry | Interpretation and current counterpart |
|---|---|
| `check_walls.py`, `check_descent.py` | Sampled-interior probes can miss a sliver near a crease endpoint. Use the analytic walls and `exact_core_descent.py` for the junction. |
| `exact_farkas.py` | Six-wall multiplier search. The current eight-wall descent certificate is `exact_core_descent.py`. The command `exact_farkas2.py` is a compatibility entry point to that same function. |
| `n4_cells.py`, `n4_cells2.py` | Historical determinant-one mask pair. Use `n4_cells3.py` for the figure-eight (determinant five). |
| `n4_tighten.py`, `n4_kkt.py`, `n4_straggler.py` | Exploratory optimization and stalled endpoints. Use `n4_exact_local.py` for the contact cone, `n4_exact_cell.py` for the capped panels and `check_n4_rerun.py` for the recorded search. |
| `n4_fig.py`, `n4_fig_loose.py` | Historical single-panel SVG and loose-state JSON. Each takes a new output filename. `n4_pair.py` owns the two paper illustrations. |
| [Historical figure-eight probes](verify/referee/README.md) | Data-generation order, inputs, output location and the corrected clipping example. |

## Certificate logs and hashes

`logs/` contains selected recorded certificate outputs. Superseded logs for the reachability, census-root and mutation checks are retained under `logs/history/`. Source snapshots in `verify/artifacts/n4-rerun-20261023/` are checked byte for byte against their original receipts. The current search modules also match those files, with a syntax-tree comparison for three comment-only edits. The runner has separately documented archive and resume handling.
The ten- and twelve-case mutation receipts are `logs/test_guards.log` and `logs/test_guards_global_cell.log`. The fourteen-case receipt is `logs/test_guards_slack_reach.log`.
The three new path checks record their output in `logs/exact_reach_geometry.log`, `logs/exact_reach_transport.log` and `logs/check_reach_paths.log`.
Integrity check:

```bash
shasum -a 256 -c logs/hashes-verify.txt
shasum -a 256 -c logs/hashes-logs.txt
shasum -a 256 -c logs/hashes-artifacts.txt
```

## Figures and exported geometry data

`verify/ribbon_geometry.py` is the common width-one geometric source used by
the trefoil renderer, exact outline/notch certificates, and JSON geometry
exporter. The `e_B` edge is the slanted line through `F`, not the horizontal
line `y = 1`, and the flush-cut silhouette includes all four panels.

```bash
uv run python verify/export_figure_geometry.py
uv run python verify/export_figure_geometry.py --check
```

`verify/artifacts/figure_geometry.json` exports the two distinguished states
(`P*` and `θ†`), notch measurements, and sampled standard/second-family curves
for external interactive visualizations. To refresh a downstream copy under its own filename, pass
`--output PATH` to the exporter, then add `--check` to verify that copy. These arrays sample the two families in Figure 4. The scripts `exact_reach.py` and `exact_branch.py` certify those particular curves. Theorem 10.4 constructs a path from an arbitrary slack start using the contact and transport formulas.

## Figures (SVG → PDF)

Figure 9 uses the same algebraic crease endpoints as the exact figure-eight
validity check. Converting those endpoints to drawing coordinates preserves
the coincident horizontal edges when the hidden-line view is clipped. The
[earlier SVG and PDF](verify/artifacts/figure-history/) preserve the rendering
that used rounded trigonometric coordinates.

```bash
export SOURCE_DATE_EPOCH=0
uv run python verify/tref_render.py        # fig1-pentagon.svg, theta-dagger-heptagon.svg
uv run python verify/render_families.py --pdf  # computed family curves and vertical zoom
uv run python verify/n4_pair.py --pdf      # fig41-loose/tight SVG and PDF
uv run python verify/make_fig_notation.py --pdf  # fig0-notation SVG and PDF
uv run python verify/render_junction.py --pdf   # fig3-mountain-pass SVG and PDF
for f in fig0-notation fig1-pentagon fig3-mountain-pass fig4-two-families fig5-slab-ramp fig41-code fig41-loose fig41-tight; do
  rsvg-convert -f pdf paper/figs/$f.svg -o paper/figs/$f.pdf; done
rsvg-convert -f pdf paper/figs/theta-dagger-heptagon.svg -o paper/figs/fig2-heptagon.pdf
```

## Building the paper

For the distributed PDF and bibliography, run the same pinned TeX Live image
as the paper CI job. The command requires Docker and writes `paper/main.pdf`
and `paper/main.bbl`. Intermediates stay inside the container.

```sh
sh scripts/build-paper.sh
```

With a local TeX installation, the manuscript also builds with the usual sequence.
Different TeX Live package versions can change glyph mappings and PDF bytes.

```bash
cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Lean formalization

```bash
cd formal && lake exe cache get && lake build Formal
```

## Packaging the archived version

`scripts/make-ancillary.sh` regenerates the hash manifests and produces
`ancillary.tar.gz` with the PDF, bibliography, figure sources and computational evidence. It can be unpacked without Git. Rebuilding from identical inputs preserves the archive bytes.
