# Eliminate W1, W5 and W7 on a census stratum. When Q0=F, the equations are dependent and leave a curve.
import sys, sympy as sp
sys.path.insert(0,'.')
import uniq_poly as P
import uniq_setup as U0
x,g,h,a,b = U0.TH
W = {n: sp.expand(P.POLY[n]) for n in ["W1","W5","W7","W9","W10","W4"]}
g_s = sp.solve(W["W1"], g)[0]
a_s = sp.solve(W["W5"], a)[0].subs(g, g_s)
b_s = sp.solve(W["W7"], b)[0].subs(g, g_s)
sub = {g: g_s, a: a_s, b: b_s}
E9  = sp.numer(sp.cancel(W["W9"].subs(sub)))
E10 = sp.numer(sp.cancel(W["W10"].subs(sub)))
R = sp.resultant(E9, E10, h)
print("[층 {W1,W5,W7,W9,W10}] res_h(E9,E10) =", sp.factor(R))
assert R == 0
print("  ⇒ 0차원 아님: F-다양체 흡수 (명제 16.7'의 대수적 재확인) — 이 층의 차단 후보는")
print("     더 작은 |S| 층이 커버 (Caratheodory) [기록 완료]")
g9 = sp.gcd(E9, E10)
print("  공통 인수 gcd(E9,E10) 차수(h):", sp.degree(g9, h), " — 공통 성분의 정의 다항")

# ---- 결과 2 실행: taco-1 층 {W1,W5,W7,W9} + KKT(λ-소거 det) ---- — run Result 2: taco-1 stratum {W1,W5,W7,W9} + KKT (λ-elimination det) ----
if __name__ == "__main__" and "--taco1" in sys.argv:
    import uniq_census  # noqa: 준비된 dELL 사용 대신 직접 구성 — built directly instead of using the prepared dELL
if "--taco1" in sys.argv:
    from sympy import Matrix, together, numer
    Np2 = sp.expand(2*P.Np); Dp = sp.expand(P.Dp)
    dNp2 = [sp.expand(sp.diff(Np2, v)) for v in U0.TH]
    dDp = [sp.expand(sp.diff(Dp, v)) for v in U0.TH]
    LIN2 = [1,1,1,0,0]
    gradL = [sp.expand(LIN2[c]*Dp**2 + dNp2[c]*Dp - Np2*dDp[c]) for c in range(5)]  # = 2Dp²∇ℓ
    S = ["W1","W5","W7","W9"]
    dW = {n: [sp.diff(W[n], v) for v in U0.TH] for n in S}
    Mx = Matrix([[dW[n][c] for n in S] + [gradL[c]] for c in range(5)])
    import time; t0=time.time()
    DET = Mx.det(method="berkowitz")
    print(f"[taco-1 {{W1,W5,W7,W9}}] det 구성 {time.time()-t0:.0f}s; 소거 대입 중...")
    t0=time.time()
    DETs = sp.cancel(sp.together(DET.subs({g: g_s, a: a_s, b: b_s})))
    E9s  = sp.cancel(W["W9"].subs({g: g_s, a: a_s, b: b_s}))
    print(f"  대입 {time.time()-t0:.0f}s; res_h 계산...")
    t0=time.time()
    Rt = sp.resultant(sp.numer(DETs), sp.numer(E9s), h)
    print(f"  res_h 완료 {time.time()-t0:.0f}s; 인수분해...")
    F = sp.factor_list(Rt)
    print("  res_h 인수 개수:", len(F[1]), "; x-실근 (0,1):")
    found=[]
    for fac,_m in F[1]:
        if fac.has(x):
            try:
                for r in sp.real_roots(sp.Poly(fac, x)):
                    rv = float(r.evalf(30))
                    if 0 < rv < 1: found.append(round(rv,8))
            except Exception as e: print("   skip:", str(e)[:50])
    print("  후보 x:", sorted(set(found)))
