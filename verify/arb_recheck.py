# arb_recheck.py — 핵심 구간 계산의 Arb(FLINT) 독립 재인증 [리뷰 04 #11] — independent Arb (FLINT) re-certification of the core interval computations [review 04 #11]
# mpmath와 독립인 ball-arithmetic 백엔드로: (i) p_x 근 브래킷+유일성 (ii) ℓ''(x†)>0 (iii) KKT λ>0 — using a ball-arithmetic backend independent of mpmath: (i) root bracket + uniqueness for p_x, (ii) ℓ''(x†)>0, (iii) KKT λ>0
from flint import arb, ctx
from fractions import Fraction as F
ctx.dps = 60

def poly_eval(coeffs, X):
    # coeffs: 유리 계수 리스트 (최고차 먼저) — Horner (arb) — coeffs: list of rational coefficients (highest degree first)
    v = arb(0)
    for c in coeffs:
        v = v*X + arb(c.numerator)/arb(c.denominator)
    return v
def to_coeffs(pairs, deg):
    cs=[F(0)]*(deg+1)
    for d,c in pairs: cs[deg-d]=c
    return cs

# p_x = x^12+4x^10-7x^4+20x^2-2
PX = to_coeffs([(12,F(1)),(10,F(4)),(4,F(-7)),(2,F(20)),(0,F(-2))],12)
DPX= to_coeffs([(11,F(12)),(9,F(40)),(3,F(-28)),(1,F(40))],11)
lo="0.32212749731122717"; hi="0.32212749731122719"
XL,XH=arb(lo),arb(hi)
HULL = XL.union(XH) if hasattr(XL,'union') else (XL+XH)/2 + (XH-XL)/2*arb("0",rad=1) if False else None
# hull: mid±rad 구성 — hull: built as mid±rad
mid=(F(32212749731122718),)  # 사용 안함 — 안전한 구성: — unused; the safe construction follows:
X_hull = arb("0.32212749731122718", "1e-17")
vlo=poly_eval(PX,XL); vhi=poly_eval(PX,XH); dv=poly_eval(DPX,X_hull)
print("[Arb] p_x(lo) =",vlo," <0 :",vlo<0)
print("[Arb] p_x(hi) =",vhi," >0 :",vhi>0)
print("[Arb] p_x' on hull =",dv," >0 :",dv>0)
print("  ⟹ 브래킷 내 근 존재·유일 (부호 변화 + 순증)" if (vlo<0 and vhi>0 and dv>0) else "  판정 실패")

# ℓ''(x†) = 2483x^11/24+849x^9/2+175x^7/4+53x^5/12−17377x^3/24+7977x/4
H2 = to_coeffs([(11,F(2483,24)),(9,F(849,2)),(7,F(175,4)),(5,F(53,12)),(3,F(-17377,24)),(1,F(7977,4))],11)
v=poly_eval(H2,X_hull)
print("[Arb] ℓ''(x†) =",v," >0 :",v>0)

# KKT 승수 (θ†): λ_w1=λ_w4, λ_w5, λ_w6 — KKT multipliers at θ†
L14 = to_coeffs([(11,F(7,48)),(9,F(7,12)),(7,F(-1,48)),(5,F(-1,16)),(3,F(-23,24)),(1,F(53,16))],11)
L5  = to_coeffs([(10,F(1,24)),(8,F(1,6)),(4,F(-1,24)),(2,F(-3,8)),(0,F(25,24))],10)
L6  = to_coeffs([(11,F(17,96)),(9,F(35,48)),(7,F(1,16)),(5,F(-1,12)),(3,F(-119,96)),(1,F(161,48))],11)
for nm,cs in (("λ_w1=λ_w4",L14),("λ_w5",L5),("λ_w6",L6)):
    v=poly_eval(cs,X_hull)
    print(f"[Arb] {nm} =",v," >0 :",v>0)
print("\n결론: mpmath와 독립인 Arb ball 산술로 세 핵심 양수성·근 유일성 재인증")

# ---- [리뷰 05 #8.2] census 벽-부호 제외(67건)의 Arb 독립 재인증 ---- — [review 05 #8.2] independent Arb re-certification of the census wall-sign exclusions (67 cases)
# uniq_certified.json의 "W_k<0 구간 인증" 판정 후보들: 해당 벽의 부호 확정 유리식을 — the candidates in uniq_certified.json whose verdict is "W_k<0 interval certificate": reduce the sign-determining rational expression of that wall
# 정수계수 분자/분모 다항으로 환원해 반경 1e-12 arb ball에서 재평가 (mpmath 독립). — to integer-coefficient numerator/denominator polynomials and re-evaluate on an arb ball of radius 1e-12 (independent of mpmath)
import sys, os, json, io, contextlib
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    import uniq_setup as US

def arb_eval_poly(P, vals):
    """정수계수 sympy Poly(5변수)를 arb ball 점에서 평가.
    Evaluate an integer-coefficient sympy Poly (5 variables) at an arb ball point."""
    tot=arb(0)
    for mono,coef in P.terms():
        term=arb(int(coef))
        for v,e in zip(vals,mono):
            if e: term=term*(v**e)
        tot=tot+term
    return tot

C=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'uniq_certified.json')))
signed={nm: US.sign_fix[nm]*raw for nm,raw in zip(US.NAMES,US.RAWS)}
npass=0; nskip=0; nfail=[]
for it in C:
    v=str(it['verdict'])
    import re
    m=re.match(r'(W\d+)<0', v)
    if not m: nskip+=1; continue
    nm=m.group(1)
    expr=sp.cancel(sp.together(signed[nm]))
    num,den=sp.fraction(expr)
    def lcm_den(e):
        ds=[sp.denom(c) for c in sp.Poly(e,US.TH).coeffs()]
        out=1
        for dd in ds: out=sp.ilcm(out,dd)
        return out
    l1=lcm_den(num); l2=lcm_den(den)
    Pn=sp.Poly(sp.expand(num*l1),US.TH); Pd=sp.Poly(sp.expand(den*l2),US.TH)
    vals=[arb(str(t),"1e-12") for t in it['theta']]
    vn=arb_eval_poly(Pn,vals); vd=arb_eval_poly(Pd,vals)
    if vd>0 and vn<0: npass+=1
    elif vd<0 and vn>0: npass+=1
    else: nfail.append((it['S'],nm))
print(f"\n[Arb] census 벽-부호 제외 재인증: 통과 {npass} / 건너뜀(비벽 판정: Krawczyk·미결·λ) {nskip} / 실패 {len(nfail)} {nfail if nfail else ''}")
