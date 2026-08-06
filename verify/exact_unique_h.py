# exact_unique_h.py — <N5,K>의 h-선형 subresultant: x† 위 임계 h의 유일성 투명화 (리뷰03 §7-1) — h-linear subresultant of <N5,K>: makes transparent the uniqueness of the critical h over x† (review03 §7-1)
import sympy as sp
x,h=sp.symbols('x h')
N5 = 8*h**3*x**3 + 4*h**2*x**4 + 4*h**2*x**2 + 2*h*x**5 - 4*h*x**3 - 6*h*x + x**6 + 3*x**4 + 3*x**2 + 1
K  = sp.expand(2*x**2*sp.diff(N5,x) + (1-x**2)*sp.diff(N5,h))
sub=sp.subresultants(N5,K,h)
L=[p for p in sub if sp.degree(p,h)==1][0]
A_,B_=sp.Poly(L,h).all_coeffs()
px=x**12+4*x**10-7*x**4+20*x**2-2
hf=sp.Rational(1,12)*(2*x**11+8*x**9-x**7-4*x**5-13*x**3+40*x)
print("h-선형 subresultant: deg A =",sp.degree(A_,x)," deg B =",sp.degree(B_,x))
print("A·h†+B mod p_x =", sp.rem(sp.expand(A_*hf+B_),px,x), " (0 ⟺ 유일성)")
print("gcd(A,p_x) =", sp.gcd(sp.Poly(sp.expand(A_),x),sp.Poly(px,x)), " (1 ⟺ A(x†)≠0)")
