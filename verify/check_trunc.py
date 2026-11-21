# Probe the stability of the outline truncation lines in Definition def:trunc.
#
# 주장: P* 근방의 유효 상태에서, 인-꼬리(M1에서 −d2 방향 후방 연장)가 코어 외곽을
# 마지막으로 가로지르는 모서리는 ℓ2 위 변 [G,Q1]이고, 아웃-꼬리(M3에서 v3 방향)는
# {y=0} 위 변 [E,H]다 — 게다가 교차는 변의 상대 내부·횡단(각도 하한)이다.
# 이들이 열린 조건이므로 정의 7.1은 P* 근방 전체에서 Kauffman의 flush 절단과 일치한다.
# Claim: for valid states in a neighborhood of P*, the edge that the in-tail (the backward
# extension from M1 along −d2) crosses last on the core outline is the edge [G,Q1] on ℓ2, and
# for the out-tail (from M3 along v3) it is the edge [E,H] on {y=0} — moreover the crossings are
# in the relative interior of the edge and transversal (angle lower bound).
# Since these are open conditions, Definition 7.1 agrees with Kauffman's flush truncation
# throughout a neighborhood of P*.

import math, random
import uniq_setup as U0
import uniq_strata as ST

def _meet_param(P1, d1, A, B):
    """ray P1+s·d1 (s>0) 와 선분 [A,B]: (s, u) — u∈(0,1)이면 변 내부 교차.
    ray P1+s·d1 (s>0) against segment [A,B]: (s, u) — u∈(0,1) means the crossing is interior to the edge."""
    d2 = (B[0]-A[0], B[1]-A[1])
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    if abs(det) < 1e-13: return None
    rhs = (A[0]-P1[0], A[1]-P1[1])
    s = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    u = (d1[0]*rhs[1] - d1[1]*rhs[0])/det
    ang = math.degrees(math.asin(min(1.0, abs(det)/math.hypot(*d1)/math.hypot(*d2))))
    return s, u, ang

def audit(th):
    c = U0._config(th)
    E, F, G, H, Q1, Q0 = c['E'], c['F'], c['G'], c['H'], c['Q1'], c['Q0']
    edges = {
        "crease1[E,F]": (E, F), "crease2[H,G]": (H, G), "crease3[Q1,Q0]": (Q1, Q0),
        "ℓ2[G,Q1]": (G, Q1), "y=0[E,H]": (E, H), "mouth[Q0,F]": (Q0, F),
    }
    out = {}
    for tail, (P0, d) in [("in", (c['M1'], (-c['d2'][0], -c['d2'][1]))),
                          ("out", (c['M3'], c['v3']))]:
        hits = []
        for name, (A, B) in edges.items():
            r = _meet_param(P0, d, A, B)
            if r is None: continue
            s, u, ang = r
            if s > 1e-9 and 0 < u < 1:
                hits.append((s, name, u, ang))
        hits.sort()
        out[tail] = hits
    return out

X0 = math.sqrt(1 - 2/math.sqrt(5))
P = U0._slice_theta(X0)
r = audit(P)
print("[P*] 인-꼬리 외곽 교차(파라미터 순):", [(f"{s:.4f}", n, f"u={u:.3f}", f"{a:.0f}°") for s, n, u, a in r['in']])
print("[P*] 아웃-꼬리 외곽 교차:", [(f"{s:.4f}", n, f"u={u:.3f}", f"{a:.0f}°") for s, n, u, a in r['out']])
last_in = r['in'][-1]; last_out = r['out'][-1]
assert last_in[1] == "ℓ2[G,Q1]" and last_out[1] == "y=0[E,H]", (last_in, last_out)
print(f"[P*] 마지막 교차: in → {last_in[1]} (s={last_in[0]:.7f}), out → {last_out[1]} (t={last_out[0]:.7f})  ✓ 정의 7.1의 두 직선")
assert abs(last_in[0] - 1.3763819) < 1e-6 and abs(last_out[0] - 1.3763819) < 1e-6

rnd = random.Random(7)
n_ok = 0; mins = dict(u_in=1, u_out=1, a_in=90, a_out=90)
tries = 0
while n_ok < 300 and tries < 100000:
    tries += 1
    th = tuple(v + rnd.uniform(-0.02, 0.02) for v in P)
    try:
        ws, c = U0._walls8f(th)
    except Exception:
        continue
    if min(ws[:6]) < -1e-9 or max(ws[6], ws[7]) < -1e-9: continue
    ok, _ = ST.valid_full(c)
    if not ok: continue
    rr = audit(th)
    if not rr['in'] or not rr['out']: continue
    li, lo = rr['in'][-1], rr['out'][-1]
    if li[1] != "ℓ2[G,Q1]" or lo[1] != "y=0[E,H]":
        print("  ✗ 예외 상태:", th, li, lo); raise SystemExit(1)
    n_ok += 1
    mins['u_in'] = min(mins['u_in'], min(li[2], 1 - li[2]))
    mins['u_out'] = min(mins['u_out'], min(lo[2], 1 - lo[2]))
    mins['a_in'] = min(mins['a_in'], li[3]); mins['a_out'] = min(mins['a_out'], lo[3])
print(f"\n[근방 감사] P* 근방(‖δθ‖∞≤0.02) 유효 상태 {n_ok}개 (시도 {tries}):")
print(f"  마지막 외곽 교차가 전부 (in → ℓ2[G,Q1], out → y=0[E,H])  ✓")
print(f"  변-내부 여유 min(u,1−u): in ≥ {mins['u_in']:.3f}, out ≥ {mins['u_out']:.3f};  교차각: in ≥ {mins['a_in']:.1f}°, out ≥ {mins['a_out']:.1f}°")
print("\n결론: 정의 7.1의 절단선(ℓ2, y=0)은 P* 근방 전체에서 '마지막 flush 모서리'와 일치 —")
print("      내부성·횡단성·최종성이 전부 열린 조건이므로 근방 명제로 승격 가능. [수치 감사 + 개방성 논증]")
