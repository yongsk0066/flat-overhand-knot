# Historical figure-eight probes

These programs preserve exploratory calculations on knot type, contact sets and truncation. Their numbering identifies the original files. The present mathematical claims and the commands that certify them are listed in the [reproducibility guide](../../REPRODUCE.md).

| Programs | Calculation |
|---|---|
| `ref1` through `ref4` | Candidate geometry, multipliers, the linearized cone and inward directions |
| `ref5`, `ref6`, `ref10`, `ref16` | A separate seeded harvest, integer-determinant and finite-tail comparisons, endpoint descents and boundary tolerances |
| `ref7` through `ref9` | Enumeration and diagram checks for the historical masks and the determinant-five pair |
| `ref11`, `ref12`, `ref13` | Truncation at the outline versus at a centreline crossing. The signed segment test in `ref11` accepts line intersections outside the segment. `ref12` uses convex clipping. |
| `ref14`, `ref15` | Numerical escape directions and contact-cone probes at the independent descent endpoints |

To check the historical mask example, run `uv run python verify/referee/ref8_assign.py`. The two masks in `n4_cells2.py` have determinant 1. The current enumeration in `n4_cells3.py` selects determinant 5. For the outline convention, `uv run python verify/referee/ref12_trunc2.py` computes the length 32√15/15.

The independent harvest and descent generate their inputs in sequence from the repository root:

```sh
uv run python verify/referee/ref5_census.py
uv run python verify/referee/ref10_stragglers.py
uv run python verify/referee/ref13_final.py
```

The first command writes `ref_hits.json` in `build/figure-eight-probes`. The second reads those hits and writes `ref10_endpoints.json` there. Probes 6 and 16 read the hits, and probes 11 and 13–15 read the endpoints. Set the same `PENTAGON_PROBE_DIR` for all commands to choose another directory. These samples use seeds starting at 424242. The main paper records a separate eight-stream search starting at 7000.
