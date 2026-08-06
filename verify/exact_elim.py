"""
09 검증: theta† 의 정확한 대수적 정체 — 최소다항식 확정.
09 verification: the exact algebraic identity of theta† — pinning down the minimal polynomial.

결과 (전부 차수 12, QQ 위 기약, 짝수 거듭제곱만 = 미러 대칭 x->-x 의 흔적):
results (all of degree 12, irreducible over QQ, only even powers = a trace of the mirror symmetry x->-x):

  p_x(x)   = x^12 + 4x^10 - 7x^4 + 20x^2 - 2
  p_h = p_a: 32h^12 - 80h^10 + 67h^8 - 36h^6 + 40h^4 - 20h^2 - 4
  p_b(b)   = 32b^12 - 144b^10 + 155b^8 - 8b^6 - 52b^4 + 20b^2 - 4
  p_g(g)   = 32g^12 + 48g^10 - 85g^8 - 188g^6 - 106g^4 - 16g^2 - 1
  p_ell(L) = 8L^12 - 88L^10 + 267L^8 - 456L^6 + 88L^4 - 80L^2 - 16

  y = x†^2 는 차수 6:  y^6 + 4y^5 - 7y^2 + 20y - 2  (Q(sqrt5) 위에서도 기약
  y = x†^2 has degree 6:  y^6 + 4y^5 - 7y^2 + 20y - 2  (irreducible over Q(sqrt5) as well
  => sqrt5 는 Q(theta†) 에 없음 — theta† 는 오각형의 황금비 체를 떠난다).
  => sqrt5 is not in Q(theta†) — theta† leaves the golden-ratio field of the pentagon).
  갈루아 군 Gal(q) = S6 (위수 720, 비가해) => theta† 는 근호로 표현 불가 —
  Galois group Gal(q) = S6 (order 720, not solvable) => theta† is not expressible by radicals —
  x0 = sqrt(1-2/sqrt5) (차수 4, 근호식) 와의 최대 대비.
  the sharpest possible contrast with x0 = sqrt(1-2/sqrt5) (degree 4, a radical expression).

  체 구조: Q(theta†) = Q(x†), 차수 12. 명시적으로 (mod p_x 기호 검증됨):
  field structure: Q(theta†) = Q(x†), of degree 12. Explicitly (symbolically verified mod p_x):
    h = a = (2x^11 + 8x^9 - x^7 - 4x^5 - 13x^3 + 40x)/12
    b     = (x^11 + 4x^9 - 7x^3 + 22x)/4
    ell   = (5x^11 + 20x^9 - x^7 - 4x^5 - 34x^3 + 106x)/12

트랙 B (기호 소거, 발견 경로):
Track B (symbolic elimination, the path of discovery):
  변수 (x,h,b);  g=(1-x^2)/(2x) [w1=0],  a=sI(x,h) [w5=0] 대입.
  variables (x,h,b);  substitute g=(1-x^2)/(2x) [w1=0] and a=sI(x,h) [w5=0].
  N4=numer(w4), N6=numer(w6);  제2족 = V(N4,N6) 의 theta† 가지.
  N4=numer(w4), N6=numer(w6);  the second family = the theta† branch of V(N4,N6).
  (1) R1=res_b(N4,N6) 의 기약인수 중 (x†,h†)에서 소멸하는 것이 정확히
  (1) among the irreducible factors of R1=res_b(N4,N6), the one vanishing at (x†,h†) is exactly
      N5=numer(sI-h)  =>  관찰 8.1 (족 위 a=h) 이 기호적으로 증명됨.
      N5=numer(sI-h)  =>  Observation 8.1 (a=h along the family) is thereby proved symbolically.
  (2) 라그랑주 임계: L = det[D*grad ell; grad N4; grad N6] (D=분모 소거),
  (2) Lagrange criticality: L = det[D*grad ell; grad N4; grad N6] (D = clearing the denominator),
      유관 기약인수 Lrel 선택.
      then select the relevant irreducible factor Lrel.
  (3) C2 = res_b(N4, Lrel) 의 소멸 인수;  p_x = res_h(N5, C2) 의 소멸 인수.
  (3) C2 = the vanishing factor of res_b(N4, Lrel);  p_x = the vanishing factor of res_h(N5, C2).
  (4) 나머지 좌표: res_x(N5,p_x) -> p_h;  res_h(N4,N5) -> C_b(b는 1차!)
  (4) the remaining coordinates: res_x(N5,p_x) -> p_h;  res_h(N4,N5) -> C_b (linear in b!)
      -> p_b;  numer(ell-t) -> res_b,res_h,res_x 사슬 -> p_ell;  res_x(g-u 분자,
      -> p_b;  numer(ell-t) -> the res_b, res_h, res_x chain -> p_ell;  res_x(numerator of g-u,
      p_x) -> p_g.

족 보조정리 (기호 증명, 이 스크립트의 부수 발견):
Family lemmas (symbolic proofs, a by-product discovery of this script):
  (L1) 가지 위에서  b = (x^2+1)/(2x) = g + x        [C_b 가 정확히 -2bx+x^2+1]
  (L1) on the branch,  b = (x^2+1)/(2x) = g + x        [C_b is exactly -2bx+x^2+1]
  (L2) 가지 위에서  ell = h + b                      [numer(ell-h-b)|_{b=(x^2+1)/2x}
  (L2) on the branch,  ell = h + b                      [numer(ell-h-b)|_{b=(x^2+1)/2x}
                                                      가 QQ(x)[h] 에서 N5 로 나눠떨어짐]
                                                      is divisible by N5 in QQ(x)[h]]
  => 제2족은 평면곡선 N5(x,h)=0 으로 완전 축약:  ell(x) = h(x) + (x^2+1)/(2x).
  => the second family reduces completely to the plane curve N5(x,h)=0:  ell(x) = h(x) + (x^2+1)/(2x).
  제3 독립 유도: 임계조건 h' = (1-x^2)/(2x^2), h' = -N5_x/N5_h 에서
  a third independent derivation: from the criticality conditions h' = (1-x^2)/(2x^2) and h' = -N5_x/N5_h,
  K := 2x^2 N5_x + (1-x^2) N5_h;  res_h(N5, K) = x^12 (x^2+1)^2 p_x  — 라그랑주
  K := 2x^2 N5_x + (1-x^2) N5_h;  res_h(N5, K) = x^12 (x^2+1)^2 p_x  — the same p_x as
  경로와 동일한 p_x. (전 계수 일치 assert.)
  obtained by the Lagrange route. (all coefficients asserted to agree.)

트랙 A (고정밀 수치, 독립 확인):
Track A (high-precision numerics, independent confirmation):
  (N4,N6,L)=0 을 정확 야코비안 뉴턴으로 dps 360 에서 풂 (350자리; 40자리
  (N4,N6,L)=0 is solved by Newton with an exact Jacobian at dps 360 (350 digits; from a 40-digit
  초기값에서 4회 제곱수렴). PSLQ(findpoly) 로 독립 재발견 — 일치.
  initial value it converges quadratically in 4 steps). Independently rediscovered by PSLQ (findpoly) — they agree.
  주의: mpmath findpoly 기본 maxsteps=100 은 차수 12 관계를 못 찾는다.
  caution: with mpmath findpoly's default maxsteps=100 a degree-12 relation cannot be found.
  08 의 "PSLQ 소득 없음"의 원인 = 차수 상한 8 + maxsteps 기본값.
  the cause of the "PSLQ yielded nothing" in 08 = the degree cap of 8 + the default maxsteps.

논리 지위: 소거 사슬은 [정확대수]. "theta† 가 어느 기약인수의 근인가"의
logical status: the elimination chain is [exact algebra]. Only the identification of "which
식별만 350자리 수치(인수값 ~1e-334 vs 타 인수 ~0.3)에 의존. p_x 기약이므로
irreducible factor theta† is a root of" relies on 350-digit numerics (that factor's value ~1e-334 vs ~0.3 for the others). Since p_x is irreducible,
p_x 가 x† 의 최소다항식이다. mod p_x 항등식 검증으로 h,b,ell,g 도 폐쇄.
p_x is the minimal polynomial of x†. Verifying the identities mod p_x closes h, b, ell, g as well.

실행:  verify/.venv/bin/python verify/exact_elim.py     (~1-2분)
Run:  verify/.venv/bin/python verify/exact_elim.py     (~1-2 min)
출력:  verify/exact_elim.json
Output:  verify/exact_elim.json
"""
import json, os, time
import mpmath as mp
from sympy import (symbols, S, Rational, cancel, together, fraction, diff,
                   Matrix, Poly, expand, numer, factor_list, resultant,
                   lambdify, lcm, div, rem, factor, sqrt as ssqrt)

T0 = time.time()
def tick(msg): print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

x, h, b, t, u, y = symbols('x h b t u y', real=True)

# 40자리 시드 (08-verdict.md) — 40-digit seed (08-verdict.md)
SEED = dict(x='0.3221274973112271792397142057801235136978',
            h='1.036386300430797104445252761398636908419',
            b='1.713244187064195830219559769216923187227')


# ---------------------------------------------------------------- 기호계 — the symbolic system
def build_system():
    """제2족의 기호계: N4, N6, N5, 라그랑주 L, ell, a=sI. 전부 (x,h,b) 유리식. — Symbolic system of the second family: N4, N6, N5, the Lagrange determinant L, ell, a=sI. All of them rational functions in (x,h,b)."""
    g = (1 - x**2)/(2*x)

    def refl_lin(d):
        dx, dy = d
        n = dx**2 + dy**2
        return lambda v: (((dx**2-dy**2)*v[0] + 2*dx*dy*v[1])/n,
                          (2*dx*dy*v[0] + (dy**2-dx**2)*v[1])/n)

    def crs(p, q): return p[0]*q[1] - p[1]*q[0]
    def sub2(p, q): return (p[0]-q[0], p[1]-q[1])
    def add2(p, q): return (p[0]+q[0], p[1]+q[1])

    E = (S(0), S(0)); F = (-x, S(1)); G = (g, S(1)); H = (h, S(0))
    d2 = refl_lin(sub2(F, E))((S(1), S(0)))
    w0 = refl_lin(sub2(G, H))((S(1), S(0)))
    w0 = (cancel(w0[0]), cancel(w0[1]))
    a = cancel(crs(sub2(F, G), d2) / crs(w0, d2))          # w5=0 => a = sI(x,h)
    Q1 = (cancel(g + a*w0[0]), cancel(1 + a*w0[1]))
    Q0 = (h + b*w0[0], b*w0[1])
    dQ = (cancel(Q0[0]-Q1[0]), cancel(Q0[1]-Q1[1]))
    rl3 = refl_lin(dQ)
    rho3 = lambda P: add2(Q1, rl3(sub2(P, Q1)))
    dD = rl3(w0)
    e1p = rho3(G)
    tick("기하 구성")
    w4 = cancel(together(crs(sub2(H, e1p), dD)))
    w6 = cancel(together(crs(sub2(E, rho3(H)), dD)))
    M2 = ((g+h)/2, Rational(1, 2))
    s23 = cancel(crs(sub2(Q1, M2), sub2(Q0, Q1)) / crs(w0, sub2(Q0, Q1)))
    ell = cancel(together((g + h + x)/2 + s23))
    tick("w4, w6, ell 소거·정리")

    N4 = expand(numer(w4)); N6 = expand(numer(w6))
    N5 = expand(numer(cancel(a - h)))                       # a=h 관찰의 (x,h) 다항식 — the (x,h) polynomial behind the a=h observation

    # 라그랑주: grad ell 행은 분모 D 를 곱해 다항식화 (det 는 D 배 — 영점 동일) — Lagrange: the grad ell row is made polynomial by multiplying through by the denominator D (the det is scaled by D — same zero set)
    fr = [fraction(cancel(diff(ell, v))) for v in (x, h, b)]
    Dl = lcm([f[1] for f in fr])
    row_ell = [expand(f[0]*cancel(Dl/f[1])) for f in fr]
    gN4 = [diff(N4, v) for v in (x, h, b)]
    gN6 = [diff(N6, v) for v in (x, h, b)]
    L = expand(Matrix([row_ell, gN4, gN6]).det(method='berkowitz'))
    tick(f"라그랑주 det L: {len(Poly(L, x, h, b).terms())}항")
    return dict(N4=N4, N6=N6, N5=N5, L=L, ell=ell, a=a, g=g)


# ---------------------------------------------------------------- 트랙 A-1: 뉴턴 — Track A-1: Newton
def newton_theta(sysd, dps=360):
    """(N4, N6, L) = 0 뉴턴, 정확 기호 야코비안. 350자리 theta†. — Newton on (N4, N6, L) = 0 with an exact symbolic Jacobian. Gives theta† to 350 digits."""
    Fs = [sysd['N4'], sysd['N6'], sysd['L']]
    fF = [lambdify((x, h, b), e, modules='mpmath') for e in Fs]
    fJ = [[lambdify((x, h, b), diff(e, v), modules='mpmath') for v in (x, h, b)]
          for e in Fs]
    mp.mp.dps = dps
    X = [mp.mpf(SEED['x']), mp.mpf(SEED['h']), mp.mpf(SEED['b'])]
    hist = []
    for it in range(40):
        Fv = mp.matrix([f(*X) for f in fF])
        Jv = mp.matrix(3, 3)
        for i in range(3):
            for j in range(3):
                Jv[i, j] = fJ[i][j](*X)
        dX = mp.lu_solve(Jv, Fv)
        X = [X[i] - dX[i] for i in range(3)]
        step = max(abs(v) for v in dX)
        hist.append(mp.nstr(step, 4))
        if step < mp.mpf(10)**(-(dps-30)):
            break
    xd, hd, bd = X
    resid = max(abs(f(*X)) for f in fF)
    fell = lambdify((x, h, b), sysd['ell'], modules='mpmath')
    fa = lambdify((x, h), sysd['a'], modules='mpmath')
    vals = dict(x=xd, h=hd, b=bd, ell=fell(xd, hd, bd), a=fa(xd, hd),
                g=(1 - xd*xd)/(2*xd))
    tick(f"뉴턴 완료: {len(hist)}회, |스텝|={hist[-1]}, |잔차|={mp.nstr(resid, 4)}")
    print(f"    a†-h† = {mp.nstr(vals['a']-vals['h'], 5)}  (관찰 8.1, {dps-10}자리)")
    return vals, dict(iters=len(hist), steps=hist, residual=mp.nstr(resid, 5))


# ---------------------------------------------------------------- 인수 선택 — factor selection
def pick_vanishing(expr, varlist, vals, name, tol='1e-150'):
    """기약인수 중 고정밀 theta† 에서 소멸하는 것 — 정확히 하나여야 한다. — The irreducible factor that vanishes at the high-precision theta† — there must be exactly one."""
    got, log = [], []
    for fac, m in factor_list(expr)[1]:
        f = lambdify(varlist, fac, modules='mpmath')
        v = f(*vals)
        P = Poly(fac, *varlist)
        log.append((tuple(P.degree(w) for w in varlist), m, mp.nstr(abs(v), 3)))
        if abs(v) < mp.mpf(tol):
            got.append((fac, v))
    assert len(got) == 1, f"{name}: 소멸 인수 {len(got)}개 (유일해야 함) — {log}"
    fac, v = got[0]
    print(f"    {name}: 인수 {len(log)}개 중 유일 소멸 "
          f"deg={tuple(int(d) for d in Poly(fac, *varlist).degree_list())} "
          f"|값|={mp.nstr(abs(v), 3)}")
    return fac, log, v


# ---------------------------------------------------------------- 트랙 B — Track B
def track_b(sysd, vals):
    N4, N6, N5, L = sysd['N4'], sysd['N6'], sysd['N5'], sysd['L']
    xd, hd, bd, ed, gd = vals['x'], vals['h'], vals['b'], vals['ell'], vals['g']
    out = {}

    # (1) a=h 의 기호 증명: R1 = res_b(N4,N6) 의 유일 소멸 인수 == N5 — (1) symbolic proof of a=h: the unique vanishing factor of R1 = res_b(N4,N6) equals N5
    R1 = resultant(N4, N6, b)
    tick("R1 = res_b(N4, N6)")
    fac1, log1, _ = pick_vanishing(R1, (x, h), (xd, hd), "R1 소멸 인수")
    a_eq_h = expand(fac1 - N5) == 0 or expand(fac1 + N5) == 0
    qq, rr = div(R1, N5, x, h)
    out['a_eq_h'] = dict(
        vanishing_factor_is_N5=bool(a_eq_h),
        N5_divides_R1=(rr == 0),
        R1_factor_degrees=[list(map(int, d)) for d, m, v in log1],
        numeric_a_minus_h=mp.nstr(vals['a']-vals['h'], 5))
    print(f"    a=h 증명: 소멸인수==N5 {a_eq_h}, N5|R1 {rr == 0}")

    # (2) 라그랑주 유관 인수 — (2) the relevant Lagrange factor
    Lrel, logL, _ = pick_vanishing(L, (x, h, b), (xd, hd, bd), "L 유관 인수")

    # (3) x† 최소다항식 — (3) the minimal polynomial of x†
    R2 = resultant(N4, Lrel, b)
    tick("R2 = res_b(N4, Lrel)")
    C2, _, _ = pick_vanishing(R2, (x, h), (xd, hd), "C2")
    Rx = resultant(N5, C2, h)
    tick("Rx = res_h(N5, C2)")
    px, logx, resx = pick_vanishing(Rx, (x,), (xd,), "p_x")
    Ppx = Poly(px, x)
    if Ppx.LC() < 0:
        px = expand(-px); Ppx = Poly(px, x)
    assert Ppx.is_irreducible
    out['px_search'] = dict(Rx_factor_degrees=[list(map(int, d)) for d, m, v in logx])

    # (4) 나머지 좌표의 최소다항식 — (4) minimal polynomials of the remaining coordinates
    def minpoly_from(rel, var, val, name):
        p, _, res = pick_vanishing(rel, (var,), (val,), name)
        P = Poly(p, var)
        if P.LC() < 0:
            p = expand(-p); P = Poly(p, var)
        assert P.is_irreducible
        return p, P, res

    ph, Pph, resh = minpoly_from(resultant(N5, px, x), h, hd, "p_h")
    pg, Ppg, resg = minpoly_from(resultant(-x**2 - 2*u*x + 1, px, x), u, gd, "p_g")
    B1 = resultant(N4, N5, h)
    Cb, _, _ = pick_vanishing(B1, (x, b), (xd, bd), "C_b")
    assert Poly(Cb, b).degree() == 1          # b 는 x 의 유리함수! — b is a rational function of x!
    pb, Ppb, resb = minpoly_from(resultant(Cb, px, x), b, bd, "p_b")
    E1 = expand(numer(cancel(together(sysd['ell'] - t))))
    S1 = resultant(E1, N4, b)
    S1r, _, _ = pick_vanishing(S1, (x, h, t), (xd, hd, ed), "S1r")
    S2 = resultant(S1r, N5, h)
    S2r, _, _ = pick_vanishing(S2, (x, t), (xd, ed), "S2r")
    pell, Ppell, rese = minpoly_from(resultant(S2r, px, x), t, ed, "p_ell")
    tick("좌표별 최소다항식 완료")

    def rec(P, res, sym):
        cs = [int(c) for c in P.all_coeffs()]
        assert all(c == 0 for c in cs[1::2]), "짝수 거듭제곱 구조 위반"
        return dict(var=sym, degree=int(P.degree()), coeffs_desc=cs,
                    irreducible=True, even_powers_only=True,
                    abs_residual_at_350digits=mp.nstr(abs(res), 4))

    out['min_polys'] = dict(
        x=rec(Ppx, resx, 'x'), h=rec(Pph, resh, 'h'), a=rec(Pph, resh, 'a'),
        b=rec(Ppb, resb, 'b'), g=rec(Ppg, resg, 'g'), ell=rec(Ppell, rese, 'ell'))

    # x†^2 의 차수 6 다항식 — the degree-6 polynomial of x†^2
    qy = Poly(Ppx.all_coeffs()[::2], y)
    fq = lambdify(y, qy.as_expr(), modules='mpmath')
    assert qy.is_irreducible
    out['x_squared'] = dict(degree=6, coeffs_desc=[int(c) for c in qy.all_coeffs()],
                            irreducible=True,
                            abs_residual=mp.nstr(abs(fq(xd*xd)), 4))

    # 갈루아 군 (차수 6 은 sympy 지원 범위): S6 => 비가해 => 근호 표현 불가 — Galois group (degree 6 is within sympy's supported range): S6 => not solvable => not expressible by radicals
    try:
        from sympy.polys.numberfields.galoisgroups import galois_group
        Ggal, _ = galois_group(qy)
        out['x_squared']['galois'] = dict(
            order=int(Ggal.order()), is_solvable=bool(Ggal.is_solvable),
            name=("S6" if int(Ggal.order()) == 720 else str(Ggal.order())),
            radical_expressible=bool(Ggal.is_solvable))
        print(f"    Galois(x†^2): 위수 {Ggal.order()}, 가해 {Ggal.is_solvable} "
              f"=> 근호 표현 {'가능' if Ggal.is_solvable else '불가'}")
    except Exception as ex:                          # sympy 버전 편차 대비 — guards against differences between sympy versions
        out['x_squared']['galois'] = dict(error=str(ex))

    # sqrt5 와의 관계: Q(sqrt5) 위에서도 기약 => sqrt5 not in Q(theta†) — relation to sqrt5: irreducible over Q(sqrt5) as well => sqrt5 not in Q(theta†)
    f5x = factor_list(px, extension=ssqrt(5))
    f5y = factor_list(qy.as_expr(), extension=ssqrt(5))
    out['golden_field'] = dict(
        px_irreducible_over_Q_sqrt5=(len(f5x[1]) == 1),
        qy_irreducible_over_Q_sqrt5=(len(f5y[1]) == 1),
        conclusion="sqrt5 not in Q(theta†) — theta† 는 오각형의 체 Q(sqrt5) 밖")

    # p_x 실근: 정확히 2개 (±x†) — x† 는 유일 양의 실근 — real roots of p_x: exactly 2 of them (±x†) — x† is the unique positive real root
    from sympy import real_roots
    rts = real_roots(Ppx)
    out['px_real_roots'] = dict(count=len(rts),
                                roots_40d=[mp.nstr(mp.mpf(str(r.evalf(45))), 40)
                                           for r in rts])
    tick("트랙 B 완료")
    return out, px, Ppx


# ------------------------------------------------- 체 구조: Q(theta†) = Q(x†) — field structure: Q(theta†) = Q(x†)
def field_structure(sysd, vals, px):
    """h,b,ell ∈ Q(x†) 의 명시 표현을 PSLQ 로 찾고 mod p_x 로 기호 검증. — Find explicit expressions for h,b,ell ∈ Q(x†) by PSLQ, then verify them symbolically mod p_x."""
    mp.mp.dps = 300
    xd = vals['x']
    powers = [xd**i for i in range(12)]
    exprs, sym_exprs = {}, {}
    for nm, val in (('h', vals['h']), ('b', vals['b']), ('ell', vals['ell'])):
        r = mp.pslq([val] + powers, maxcoeff=10**14, maxsteps=200000)
        assert r and r[0] != 0, f"{nm}: Q(x†) 표현 PSLQ 실패"
        c0 = -r[0]
        e = sum(Rational(int(r[i+1]), int(c0))*x**i for i in range(12))
        resid = mp.nstr(abs(sum(mp.mpf(c)*v for c, v in zip(r, [val]+powers))), 4)
        exprs[nm] = dict(pslq_vector=[int(v) for v in r], expr=str(e),
                         abs_residual_at_300digits=resid)
        sym_exprs[nm] = e
    # h_expr 은 족 방정식 N5 도 만족해야 (mod p_x) — h_expr must also satisfy the family equation N5 (mod p_x)
    n5ok = rem(expand(sysd['N5'].subs(h, sym_exprs['h'])), px, x) == 0
    return exprs, sym_exprs, n5ok


def verify_membership_mod_px(sym_exprs, minp, px):
    """p_h(h_expr) ≡ p_b(b_expr) ≡ p_ell(ell_expr) ≡ 0  (mod p_x)  [정확대수]. — p_h(h_expr) ≡ p_b(b_expr) ≡ p_ell(ell_expr) ≡ 0  (mod p_x)  [exact algebra]."""
    out = {}
    for nm, var in (('h', h), ('b', b), ('ell', t)):
        P = Poly(minp['min_polys'][nm]['coeffs_desc'], var).as_expr()
        r = rem(expand(P.subs(var, sym_exprs[nm])), px, x)
        out[nm] = (r == 0)
    # 부수 항등식: ell_expr == h_expr + b_expr  (정확히, 표현 수준에서) — side identity: ell_expr == h_expr + b_expr  (exactly, at the level of the expressions)
    out['ell_eq_h_plus_b_expr'] = \
        expand(sym_exprs['ell'] - sym_exprs['h'] - sym_exprs['b']) == 0
    return out


# ------------------------------------------------- 족 보조정리 + 제3 유도 — family lemmas + the third derivation
def family_lemmas(sysd, vals, px, Ppx):
    """(L1) b=(x^2+1)/(2x)=g+x, (L2) ell=h+b [기호];  제3 유도로 p_x 재확인. — (L1) b=(x^2+1)/(2x)=g+x, (L2) ell=h+b [symbolic];  p_x reconfirmed by the third derivation."""
    from sympy import QQ, solve
    N4, N5, N6 = sysd['N4'], sysd['N5'], sysd['N6']
    out = {}

    # (L1) C_b = res_h(N4,N5) 의 유관 인수 — 정확히 -2bx + x^2 + 1 — (L1) C_b = the relevant factor of res_h(N4,N5) — exactly -2bx + x^2 + 1
    B1 = resultant(N4, N5, h)
    Cb, _, _ = pick_vanishing(B1, (x, b), (vals['x'], vals['b']), "C_b(보조정리)")
    bsol = solve(Cb, b)[0]
    out['L1_b_closed_form'] = dict(C_b=str(Cb), b_on_branch=str(bsol),
                                   equals_g_plus_x=(cancel(bsol - sysd['g'] - x) == 0))

    # (L2) Z = numer(ell-h-b);  Z|_{b=bsol} ≡ 0 (mod N5) in QQ(x)[h]
    Z = expand(numer(cancel(together(sysd['ell'] - h - b))))
    Zc = expand(numer(cancel(together(Z.subs(b, bsol)))))
    PZc = Poly(Zc, h, domain=QQ.frac_field(x))
    PN5 = Poly(N5, h, domain=QQ.frac_field(x))
    _, r = PZc.div(PN5)
    out['L2_ell_eq_h_plus_b'] = dict(Z_terms=len(Poly(Z, x, h, b).terms()),
                                     remainder_zero=bool(r.is_zero))
    assert r.is_zero, "족 항등식 ell=h+b 기호 증명 실패"

    # 제3 독립 유도: K = 2x^2 N5_x + (1-x^2) N5_h;  res_h(N5,K) 의 소멸 인수 == p_x — third independent derivation: K = 2x^2 N5_x + (1-x^2) N5_h;  the vanishing factor of res_h(N5,K) equals p_x
    K = expand(2*x**2*diff(N5, x) + (1 - x**2)*diff(N5, h))
    Rx2 = resultant(N5, K, h)
    px2, _, _ = pick_vanishing(Rx2, (x,), (vals['x'],), "p_x(제3 유도)")
    P2 = Poly(px2, x)
    if P2.LC() < 0:
        P2 = Poly(expand(-px2), x)
    same = P2.all_coeffs() == Ppx.all_coeffs()
    out['third_derivation'] = dict(K=str(K), px_matches_lagrange_route=bool(same))
    assert same, "제3 유도의 p_x 가 라그랑주 경로와 불일치"
    print(f"    족 보조정리: b=g+x {out['L1_b_closed_form']['equals_g_plus_x']}, "
          f"ell=h+b {r.is_zero}, 제3 유도 p_x 일치 {same}")
    return out


# ---------------------------------------------------------------- 트랙 A-2: PSLQ — Track A-2: PSLQ
def track_a_pslq(vals, Ppx, minp):
    """독립 재발견: findpoly 를 차수/계수한도 격자로. (발견 시 350자리 검증) — Independent rediscovery: findpoly run over a grid of degree / coefficient bounds. (Anything found is verified to 350 digits.)"""
    mp.mp.dps = 250
    targets = dict(x=vals['x'], h=vals['h'], b=vals['b'], ell=vals['ell'])
    known = {nm: minp['min_polys'][nm]['coeffs_desc'] for nm in targets}
    out = {}
    for nm, val in targets.items():
        found, log = None, []
        for H in (10**6, 10**12, 10**18):
            for deg in (2, 3, 4, 6, 8, 10, 12, 16, 20, 24):
                if (deg+1)*mp.log(H, 10) > 220:      # dps 250 검출 한계 — the detection limit at dps 250
                    log.append(dict(degree=deg, maxcoeff=float(H), skipped=True))
                    continue
                r = mp.findpoly(val, deg, maxcoeff=H, maxsteps=20000)
                if r is None:
                    log.append(dict(degree=deg, maxcoeff=float(H), found=False))
                    continue
                # 350자리 검증 — verification to 350 digits
                mp.mp.dps = 360
                resid = abs(sum(mp.mpf(int(c))*val**(len(r)-1-i)
                                for i, c in enumerate(r)))
                mp.mp.dps = 250
                if resid < mp.mpf('1e-250'):
                    found = dict(degree=len(r)-1, maxcoeff=float(H),
                                 coeffs_desc=[int(c) for c in r],
                                 abs_residual_at_350digits=mp.nstr(resid, 4))
                    break
                log.append(dict(degree=deg, maxcoeff=float(H), found=True,
                                spurious=True))
            if found:
                break
        # 트랙 B 와 대조 (부호/상수배 무시: 원시부 비교) — compare against Track B (ignoring sign/constant multiples: the primitive parts are compared)
        agree = None
        if found:
            c = found['coeffs_desc']
            k = known[nm]
            from math import gcd
            from functools import reduce
            cc = [v//reduce(gcd, map(abs, c)) for v in c]
            if cc[0] < 0:
                cc = [-v for v in cc]
            agree = (cc == k)
        out[nm] = dict(found=found, agrees_with_track_B=agree, search_log=log)
        print(f"    PSLQ {nm}: " +
              (f"차수 {found['degree']} 발견, 트랙B 일치={agree}" if found else "미발견"))
    out['note'] = ("mpmath findpoly 기본 maxsteps=100 으로는 차수 12 관계 탐지 불가. "
                   "08 의 'PSLQ 소득 없음' = 차수<=8 상한 + maxsteps 기본값의 결과.")
    tick("트랙 A PSLQ 완료")
    return out


# ---------------------------------------------------------------- main
def main():
    sysd = build_system()
    P4 = Poly(sysd['N4'], x, h, b); P6 = Poly(sysd['N6'], x, h, b)
    PL = Poly(sysd['L'], x, h, b)
    vals, newton_info = newton_theta(sysd)

    trackB, px, Ppx = track_b(sysd, vals)
    lemmas = family_lemmas(sysd, vals, px, Ppx)
    exprs, sym_exprs, n5ok = field_structure(sysd, vals, px)
    memb = verify_membership_mod_px(sym_exprs, trackB, px)
    print(f"    mod p_x 기호 검증: {memb}, N5(h_expr)≡0: {n5ok}")
    assert all(memb.values()) and n5ok, "mod p_x 검증 실패"
    pslq = track_a_pslq(vals, Ppx, trackB)

    mp.mp.dps = 310
    result = dict(
        script="verify/exact_elim.py",
        summary=dict(
            x_minpoly="x^12 + 4x^10 - 7x^4 + 20x^2 - 2",
            all_coordinates_degree=12,
            field="Q(theta†) = Q(x†), [Q(x†):Q] = 12; x†^2 는 차수 6",
            sqrt5_in_field=False,
            identities=["족 전체: ell = h + b [기호 증명]",
                        "가지 위: b = (x^2+1)/(2x) = g + x [기호 증명]",
                        "족 전체: a = h [기호 증명]"]),
        system=dict(
            N4_degrees=[int(P4.degree(v)) for v in (x, h, b)],
            N4_terms=len(P4.terms()),
            N6_degrees=[int(P6.degree(v)) for v in (x, h, b)],
            N6_terms=len(P6.terms()),
            L_degrees=[int(PL.degree(v)) for v in (x, h, b)],
            L_terms=len(PL.terms()),
            N5=str(sysd['N5'])),
        newton=newton_info,
        theta_dagger_300digits={k: mp.nstr(v, 300) for k, v in vals.items()},
        track_B=trackB,
        family_lemmas=lemmas,
        field_structure=dict(expressions=exprs,
                             membership_verified_mod_px=memb,
                             h_expr_satisfies_N5_mod_px=n5ok),
        track_A_pslq=pslq)

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'exact_elim.json')
    with open(path, 'w') as f:
        json.dump(result, f, ensure_ascii=False, indent=1)
    tick(f"저장: {path}")

    print("\n=== 최소다항식 (계수 내림차순) ===")
    for nm in ('x', 'h', 'b', 'g', 'ell'):
        d = trackB['min_polys'][nm]
        print(f"  p_{nm}: deg {d['degree']}, {d['coeffs_desc']}, "
              f"|잔차|={d['abs_residual_at_350digits']}")
    print(f"  x†^2: deg 6, {trackB['x_squared']['coeffs_desc']}")
    print("\n=== Q(x†) 표현 ===")
    for nm in ('h', 'b', 'ell'):
        print(f"  {nm} = {exprs[nm]['expr']}")


if __name__ == '__main__':
    main()
