import Mathlib

/-!
# ThetaDagger — the algebra of the tight state θ†

`p x = x¹² + 4x¹⁰ − 7x⁴ + 20x² − 2` is the minimal polynomial of the
x-coordinate of θ† (irreducibility and the S₆ Galois group are sympy
certificates, `verify/exact_elim.py`). Machine-checked here:

1. `ell_eq_h_add_b` — the field-representation identity ℓ = h + b (exact over ℚ[x])
2. `deriv_pos` — p′ > 0 on (0,1]
3. `root_exists` — existence of a root in (0.32, 0.33) (intermediate value theorem)
4. `p_strictMonoOn` / `root_unique` — uniqueness of the root in (0,1)
5. `pentagon_not_thetaDagger` — no x satisfying the pentagon relations
   (s²=5, s·x²=s−2) is a root of p: the machine-checked **algebraic
   separation of the pentagonal junction from θ†**. The
   `linear_combination` certificate was pre-verified in exact algebra
   by `verify/certify.py` (certificate C4).
-/

namespace ThetaDagger

noncomputable def p (x : ℝ) : ℝ := x ^ 12 + 4 * x ^ 10 - 7 * x ^ 4 + 20 * x ^ 2 - 2

/-- Field-representation identity: ℓ† = h† + b† holds as an exact
identity at the level of ℚ[x] (the coordinate form of ℓ = h + b). -/
theorem ell_eq_h_add_b (x : ℝ) :
    (5 * x ^ 11 + 20 * x ^ 9 - x ^ 7 - 4 * x ^ 5 - 34 * x ^ 3 + 106 * x) / 12
      = (2 * x ^ 11 + 8 * x ^ 9 - x ^ 7 - 4 * x ^ 5 - 13 * x ^ 3 + 40 * x) / 12
        + (x ^ 11 + 4 * x ^ 9 - 7 * x ^ 3 + 22 * x) / 4 := by
  ring

/-- The derivative of p is positive on (0,1]: 12x¹¹ + 40x⁹ − 28x³ + 40x > 0,
since 40x − 28x³ = 12x + 28x(1−x²) ≥ 12x > 0. -/
theorem deriv_pos (x : ℝ) (hx : 0 < x) (hx1 : x ≤ 1) :
    0 < 12 * x ^ 11 + 40 * x ^ 9 - 28 * x ^ 3 + 40 * x := by
  have h2 : x ^ 2 ≤ 1 := by nlinarith
  nlinarith [pow_pos hx 11, pow_pos hx 9,
             mul_nonneg hx.le (sub_nonneg.mpr h2)]

theorem p_hasDerivAt (x : ℝ) :
    HasDerivAt p (12 * x ^ 11 + 40 * x ^ 9 - 28 * x ^ 3 + 40 * x) x := by
  have h :=
    ((((hasDerivAt_pow 12 x).add ((hasDerivAt_pow 10 x).const_mul 4)).sub
        ((hasDerivAt_pow 4 x).const_mul 7)).add
      ((hasDerivAt_pow 2 x).const_mul 20)).sub_const 2
  convert h using 1
  · rfl
  · rfl
  · funext y
    simp [p]
  · norm_num
    ring

/-- p is strictly increasing on [0,1]. -/
theorem p_strictMonoOn : StrictMonoOn p (Set.Icc (0 : ℝ) 1) := by
  apply strictMonoOn_of_deriv_pos (convex_Icc 0 1)
  · exact fun x _ => ((p_hasDerivAt x).continuousAt).continuousWithinAt
  · intro x hx
    rw [interior_Icc] at hx
    rw [(p_hasDerivAt x).deriv]
    exact deriv_pos x hx.1 hx.2.le

/-- Existence: p has a root in (0.32, 0.33), by the intermediate value theorem. -/
theorem root_exists : ∃ x ∈ Set.Ioo (0.32 : ℝ) 0.33, p x = 0 := by
  have hc : ContinuousOn p (Set.Icc (0.32 : ℝ) 0.33) := by
    intro x _
    exact ((p_hasDerivAt x).continuousAt).continuousWithinAt
  have hmem : (0 : ℝ) ∈ Set.Ioo (p 0.32) (p 0.33) := by
    constructor
    · norm_num [p]
    · norm_num [p]
  have h := intermediate_value_Ioo (by norm_num : (0.32 : ℝ) ≤ 0.33) hc hmem
  obtain ⟨x, hx, hpx⟩ := h
  exact ⟨x, hx, hpx⟩

/-- Uniqueness: p has at most one root in (0,1). -/
theorem root_unique {a b : ℝ} (ha : a ∈ Set.Ioo (0 : ℝ) 1) (hb : b ∈ Set.Ioo (0 : ℝ) 1)
    (hpa : p a = 0) (hpb : p b = 0) : a = b := by
  rcases lt_trichotomy a b with h | h | h
  · have := p_strictMonoOn (Set.mem_Icc_of_Ioo ha) (Set.mem_Icc_of_Ioo hb) h
    rw [hpa, hpb] at this
    exact absurd this (lt_irrefl 0)
  · exact h
  · have := p_strictMonoOn (Set.mem_Icc_of_Ioo hb) (Set.mem_Icc_of_Ioo ha) h
    rw [hpa, hpb] at this
    exact absurd this (lt_irrefl 0)

/-- **Algebraic separation of the junction from the tight state**: no x
satisfying the pentagon relations s² = 5, s·x² = s − 2 is a root of p.
Certificate (`verify/certify.py`, C4):
125·p(x) = (9664 − 4320s) + α·(sx² − s + 2) + β·(s² − 5). -/
theorem pentagon_not_thetaDagger (s x : ℝ) (hs : s ^ 2 = 5) (hx : s * x ^ 2 = s - 2) :
    p x ≠ 0 := by
  intro h0
  have key : (125 : ℝ) * p x = 9664 - 4320 * s := by
    unfold p
    linear_combination
      (s * (25 * x ^ 10 + 125 * x ^ 8 + 145 * x ^ 6 + 265 * x ^ 4 + 326 * x ^ 2 + 1274)
        - 50 * x ^ 8 - 300 * x ^ 6 - 590 * x ^ 4 - 1120 * x ^ 2 - 1772) * hx
      + (-25 * x ^ 12 - 100 * x ^ 10 - 20 * x ^ 8 - 120 * x ^ 6 - 61 * x ^ 4
        - 948 * x ^ 2 + 1274) * hs
  rw [h0, mul_zero] at key
  have h1 : (4320 : ℝ) * s = 9664 := by linarith
  have h2 : (4320 : ℝ) ^ 2 * s ^ 2 = 9664 ^ 2 := by nlinarith [h1]
  rw [hs] at h2
  norm_num at h2

end ThetaDagger
