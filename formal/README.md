# Lean algebra certificates

Run from this directory:

```sh
lake exe cache get
lake build Formal
```

`lean-toolchain` and `lake-manifest.json` pin Lean and Mathlib.
There are 16 named theorems: nine in `Formal/Pentagon.lean` and seven in
`Formal/ThetaDagger.lean`. They prove polynomial identities, sign and root
properties, the slice core-length formula and its monotonicity, and the
algebraic separation of the pentagon parameter from the dagger parameter.

`side_eq_of_contact` proves the equality of two squared-length expressions,
not a theorem about all five polygon sides. The fifth-root-of-unity lemma
is an algebraic implication, not a formal construction of the polygon.
The coordinate-to-geometry bridge, layer axioms, non-penetration, complete
outline, Galois group, optimization and realization theorems are **not**
formalized here. Their human proofs and symbolic/interval certificates are
in the paper and `../verify/`.
