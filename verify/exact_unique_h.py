"""Uniqueness of the critical ordinate h-dagger above x-dagger via subresultants.

Extracts the h-linear polynomial A(x)*h + B(x) from the subresultant chain of
(N5, K) with respect to h and certifies that gcd(A, p_x) = 1 and
A(x)*h_dagger(x) + B(x) == 0 (mod p_x).
"""
import sympy as sp

def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)


def certify():
    x, h = sp.symbols("x h")
    N5 = (
        8 * h ** 3 * x ** 3
        + 4 * h ** 2 * x ** 4
        + 4 * h ** 2 * x ** 2
        + 2 * h * x ** 5
        - 4 * h * x ** 3
        - 6 * h * x
        + x ** 6
        + 3 * x ** 4
        + 3 * x ** 2
        + 1
    )
    K = sp.expand(2 * x ** 2 * sp.diff(N5, x) + (1 - x ** 2) * sp.diff(N5, h))
    sub = sp.subresultants(N5, K, h)
    linear = [p for p in sub if sp.degree(p, h) == 1]
    _require(linear, "No h-linear subresultant found")
    A_, B_ = sp.Poly(linear[0], h).all_coeffs()

    px = x ** 12 + 4 * x ** 10 - 7 * x ** 4 + 20 * x ** 2 - 2
    hf = sp.Rational(1, 12) * (2 * x ** 11 + 8 * x ** 9 - x ** 7 - 4 * x ** 5 - 13 * x ** 3 + 40 * x)

    rem_val = sp.rem(sp.expand(A_ * hf + B_), px, x)
    gcd_val = sp.gcd(sp.Poly(sp.expand(A_), x, domain="QQ"), sp.Poly(px, x, domain="QQ"))
    print(f"h-linear subresultant: deg A = {sp.degree(A_, x)}, deg B = {sp.degree(B_, x)}")
    print(f"A*h_dagger + B mod p_x = {rem_val}")
    print(f"gcd(A, p_x) = {gcd_val.as_expr()}")
    _require(rem_val == 0, f"h-linear relation failed mod p_x: {rem_val}")
    _require(gcd_val.is_one, f"A shares a factor with p_x: {gcd_val}")
    print("PASS: critical ordinate h-dagger is unique above x-dagger")


if __name__ == "__main__":
    certify()
