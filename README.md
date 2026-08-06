# Two notions of tightness for the flat overhand knot

[![verify](https://github.com/yongsk0066/flat-overhand-knot/actions/workflows/verify.yml/badge.svg)](https://github.com/yongsk0066/flat-overhand-knot/actions/workflows/verify.yml)
[![lean](https://github.com/yongsk0066/flat-overhand-knot/actions/workflows/lean.yml/badge.svg)](https://github.com/yongsk0066/flat-overhand-knot/actions/workflows/lean.yml)
[![paper](https://github.com/yongsk0066/flat-overhand-knot/actions/workflows/paper.yml/badge.svg)](https://github.com/yongsk0066/flat-overhand-knot/actions/workflows/paper.yml)

Tie an overhand knot in a paper strip and pull it flat: folklore says you get a regular
pentagon. This repository contains a paper, and its complete verification suite, showing
that *what "tight" means decides what you get*. On the same five-dimensional moduli space
of three-crease flat folded states, with classified non-penetration walls:

- **Kauffman's truncation length** selects the regular pentagon as a strict local
  minimizer (an exact golden-ratio Farkas certificate with linear growth);
- **the core length that pulling extracts** is *not* minimized there: a certified
  pulling path crosses the pentagon and terminates at a strict local minimizer
  θ† whose coordinate field has degree 12 with Galois group S₆ — no radical
  expression exists; its outline is a heptagon with an exact V-notch;
- one crossing higher, the **figure-eight branch** tightens to an exactly solvable
  state in ℚ(√15) — core length 4√15/3, a centrally symmetric hexagonal outline,
  and truncation length 32/√15, exactly Kauffman's conjectured figure-eight value.

Everything is certified: exact arithmetic in quotient fields, resultant–Sturm root
isolation, Krawczyk interval enclosures, an independent Arb re-verification, and a
Lean 4 formalization of the algebraic core.

## Layout

| Path | Contents |
|---|---|
| `paper/` | LaTeX source and the compiled paper (`main.pdf`) |
| `paper/figs/` | Figure SVG sources and the PDFs the paper includes |
| `verify/` | The certificate suite: numerical probes, exact symbolic derivations, interval certificates, and the internal adversarial-referee scripts (`verify/referee/`) |
| `formal/` | Lean 4 project (16 machine-checked statements against Mathlib) |
| `logs/` | Certificate output logs and SHA-256 manifests |
| `REPRODUCE.md` | Step-by-step reproduction commands |

## Quick start

```bash
uv sync                                    # locked environment from pyproject.toml + uv.lock
uv run python verify/n4_exact_valid.py     # ~1 min: exact validity at the figure-eight tight state
uv run python verify/exact_taco.py         # ~2 min: taco entanglement identities and certificates
```

See `REPRODUCE.md` for the full suite (some certificates run for 10–20 minutes) and for
figure regeneration. The Lean project builds with `cd formal && lake exe cache get && lake build`.

## Citation

See `CITATION.cff`. An arXiv identifier and a Zenodo DOI for the archived version will
be added at submission.

## License

Code (`verify/`, `formal/`, figure scripts): MIT — see `LICENSE`.
Paper text and figures (`paper/`): © 2026 Yongseok Jang, licensed CC BY 4.0.
