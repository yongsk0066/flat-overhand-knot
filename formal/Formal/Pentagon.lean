/-
Copyright (c) 2026 Yongseok Jang. All rights reserved.
Released under the MIT license as described in the file LICENSE.
Authors: Yongseok Jang
-/
module

public import Mathlib

/-!
# Pentagon — the algebraic core of the standard-slice analysis

Machine-checked algebra supporting the folklore half of the pentagonal-junction theorem
(Theorem C of the paper). Every `linear_combination` certificate below was
pre-verified in exact arithmetic by `verify/exact_slice.py` (certificates
E1–E4); coordinates follow the standard-slice normalization of the paper.

The nine theorems in this file are

1. `delta_formula` — closed form of the gap function Δ
2. `deriv_numerator_pos` — positivity of the numerator of Δ′ (⟹ Δ strictly decreasing)
3. `p_pos` — p(y) > 0 for y > 0, by an SOS decomposition
4. `x0_root` — the tight parameter is a root of the contact polynomial
5. `root_unique` — uniqueness of that root in the standard range
6. `side_eq_of_contact` — contact ⟹ |EH|² = |EF|² (a P(x) identity)
7. `closure_root_of_unity` — side-vector closure ⟹ z⁵ = 1 (the 72° argument)
8. `ell_core_closed_form` — ℓ_core = 1/x − x on the slice
9. `ell_core_strict_anti` — ℓ_core strictly decreasing (discrete form)
-/

public section

namespace Pentagon

/-- **Closed form of the gap function**:
`k(x) − λ(x) = (5x⁴−10x²+1) / (4x(1−x²))`. -/
theorem delta_formula (x : ℝ) (hx : x ≠ 0) (h1 : 1 - x ^ 2 ≠ 0) :
    (1 - 3 * x ^ 2) / (2 * x) - (1 + x ^ 2) ^ 2 / (4 * x * (1 - x ^ 2))
      = (5 * x ^ 4 - 10 * x ^ 2 + 1) / (4 * x * (1 - x ^ 2)) := by
  field_simp
  ring

/-- **Positivity of the numerator of Δ′**: for every real x,
`5x⁶ − 5x⁴ + 7x² + 1 > 0`, by the SOS decomposition
`5x²(x²−1/2)² + (23/4)x² + 1` (degree-6 form of certificate E4). -/
theorem deriv_numerator_pos (x : ℝ) :
    0 < 5 * x ^ 6 - 5 * x ^ 4 + 7 * x ^ 2 + 1 := by
  nlinarith [mul_nonneg (sq_nonneg x) (sq_nonneg (x ^ 2 - 1 / 2)),
             sq_nonneg x, sq_nonneg (x ^ 2 - 1 / 2)]

/-- `p(y) = 5y³−5y²+7y+1 > 0` for y > 0 by the sum-of-squares decomposition
`5y(y−1/2)² + (23/4)y + 1` (certificate E4). -/
theorem p_pos (y : ℝ) (hy : 0 < y) : 0 < 5 * y ^ 3 - 5 * y ^ 2 + 7 * y + 1 := by
  nlinarith [mul_nonneg hy.le (sq_nonneg (y - 1 / 2)), hy]

/-- **The tight parameter is a root of the contact polynomial**:
if `s² = 5` and `s·x² = s − 2` (the denominator-cleared form of
x² = 1 − 2/√5), then `5x⁴ − 10x² + 1 = 0`. Certificate E1. -/
theorem x0_root (s x : ℝ) (hs : s ^ 2 = 5) (hx : s * x ^ 2 = s - 2) :
    5 * x ^ 4 - 10 * x ^ 2 + 1 = 0 := by
  linear_combination (s * x ^ 2 - s - 2) * hx - (x ^ 2 - 1) ^ 2 * hs

/-- **Uniqueness of the root in the standard range**: if `y = x²`
satisfies `5y² − 10y + 1 = 0` and `y < 1/3` (the range `0<x<1/√3`),
then `s·y = s − 2`, i.e. only the branch `y = 1 − 2/√5` survives. -/
theorem root_unique (s y : ℝ) (hs : s ^ 2 = 5) (hs0 : 0 < s)
    (hy : 5 * y ^ 2 - 10 * y + 1 = 0) (hy3 : y < 1 / 3) :
    s * y = s - 2 := by
  have hfac : (s * y - s + 2) * (s * y - s - 2) = 0 := by
    linear_combination hy + (y - 1) ^ 2 * hs
  rcases mul_eq_zero.mp hfac with h' | h'
  · linarith
  · exfalso
    have hmul : s * y < s * (1 / 3) := by
      exact mul_lt_mul_of_pos_left hy3 hs0
    nlinarith

/-- **Contact ⟹ equal sides**: if the contact polynomial vanishes, then
`|EH|² = ((1−3x²)/(2x))²` equals `|EF|² = 1 + x²`. Certificate E2. -/
theorem side_eq_of_contact (x : ℝ) (hx : x ≠ 0)
    (h : 5 * x ^ 4 - 10 * x ^ 2 + 1 = 0) :
    ((1 - 3 * x ^ 2) / (2 * x)) ^ 2 = 1 + x ^ 2 := by
  field_simp
  linear_combination h

/-- **Pentagon closure ⟹ fifth root of unity**: if the side vectors sum
to zero then `z⁵ = 1` — the source of the 72° exterior angle. Certificate E3. -/
theorem closure_root_of_unity (z : ℂ) (h : 1 + z + z ^ 2 + z ^ 3 + z ^ 4 = 0) :
    z ^ 5 = 1 := by
  linear_combination (z - 1) * h

/-- **Closed form of ℓ_core on the slice**: the two core segments sum
to `1/x − x`. -/
theorem ell_core_closed_form (x : ℝ) (hx : x ≠ 0) :
    (1 - x ^ 2) / (2 * x) + (1 - x ^ 2) / (2 * x) = 1 / x - x := by
  field_simp
  ring

/-- **ℓ_core is strictly decreasing**: if `0 < a < b` then
`ℓ_core(b) < ℓ_core(a)` — the formalization of pulling as increasing x. -/
theorem ell_core_strict_anti (a b : ℝ) (ha : 0 < a) (hab : a < b) :
    1 / b - b < 1 / a - a := by
  have hb : 0 < b := ha.trans hab
  have h1 : 1 / b < 1 / a := one_div_lt_one_div_of_lt ha hab
  linarith

end Pentagon
