"""Kauffman 절단 길이 규약 검증과 오각형/theta† 비교 (11 §2).
   Verification of the Kauffman truncation-length convention, and pentagon/theta† comparison (11 §2).
   규약: L_trunc = 레이인 스텁 + ell_core + 레이아웃 스텁 (중심선, 외곽 flush 절단).
   Convention: L_trunc = ray-in stub + ell_core + ray-out stub (centerline, flush truncation at the outline).
   오각형에서 4/sqrt(5-2 sqrt5)와 10자리 일치 확인.
   Confirms 10-digit agreement with 4/sqrt(5-2 sqrt5) for the pentagon."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_edge import config, slice_theta

def sub(a,b): return (a[0]-b[0],a[1]-b[1])
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,s): return (a[0]*s,a[1]*s)
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def norm(a): return math.hypot(*a)

def meetpt(P1,d1,P2,d2):
    det=d1[0]*(-d2[1])-d1[1]*(-d2[0]); rhs=sub(P2,P1)
    s=(rhs[0]*(-d2[1])-rhs[1]*(-d2[0]))/det
    return add(P1,mul(d1,s)), s

def L_trunc(theta):
    c = config(theta)
    M1,M2,M3 = c['M1'],c['M2'],c['M3']
    ell = norm(sub(M2,M1)) + norm(sub(M3,M2))
    _, s_in = meetpt(M1, mul(c['d2'],-1.0), c['G'], sub(c['Q1'],c['G']))
    _, t_out = meetpt(M3, c['v3'], (0.0,0.0), (1.0,0.0))
    return s_in + ell + t_out, ell, s_in, t_out

x0 = math.sqrt(1-2/math.sqrt(5))
Lp, ep, sp, tp = L_trunc(slice_theta(x0))
K = 4/math.sqrt(5-2*math.sqrt(5))
print(f"오각형: L_trunc = {Lp:.10f}  (Kauffman 4/sqrt(5-2sqrt5) = {K:.10f}, 차 {Lp-K:.1e})")
print(f"        = 스텁 {sp:.7f} + core {ep:.7f} + 스텁 {tp:.7f}")
thd = (0.3221274973112272, 1.3911166897529687, 1.0363863004307971, 1.0363863004307971, 1.7132441870641958)
Ld, ed, sd, td = L_trunc(thd)
print(f"theta†: L_trunc = {Ld:.10f}  = 스텁 {sd:.7f} + core {ed:.7f} + 스텁 {td:.7f}")
print(f"\n이분법: core는 theta†가 승 ({ed:.7f} < {ep:.7f}), L_trunc는 오각형이 승 ({Lp:.7f} < {Ld:.7f})")
