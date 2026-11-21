# check_gauss.py — combinatorial verification for Theorem A (fold number and crossing rigidity).
#
# G1. Gauss 홀짝 조건: n=3 통과열 후보 4개 중 평면 실현 가능한 것은 ABCABC뿐. — Gauss parity condition: of the 4 candidate n=3 passage words, only ABCABC is realizable in the plane.
# G2. ≤2 교차 코드: 실현 가능한 것은 전부 인접 반복(R1 컬) 보유 → 축약 → 자명. — codes with ≤2 crossings: every realizable one carries an adjacent repeat (an R1 curl) → reduce → unknot.
# G3. 8개 상하 배정 = sigma_1^{s1,s2,s3}의 폐포: Kauffman 괄호 → Jones. — the 8 over/under assignments = the closure of sigma_1^{s1,s2,s3}: Kauffman bracket → Jones.
#     혼합 부호 6개 → V=1 (R2 소거와 정합), 상수 부호 2개 → 세잎 다항식 (비자명). — the 6 mixed-sign cases → V=1 (consistent with R2 cancellation), the 2 constant-sign cases → the trefoil polynomial (nontrivial).
# G4. 삼색성(mod 3): 상수 부호 → 9 (비자명 인증, 초등), 혼합 → 3. — tricolorability (mod 3): constant sign → 9 (an elementary certificate of nontriviality), mixed → 3.

import itertools
import sympy as sp

A, t = sp.symbols('A t')

# ---------- G1 + G2: Gauss 홀짝 ---------- — Gauss parity
def parity_ok(code):
    """실현 가능성의 필요조건: 각 기호의 두 출현 사이 기호 수가 짝수.
    A necessary condition for realizability: the number of symbols between the two occurrences of each symbol is even."""
    for ch in set(code):
        i, j = code.index(ch), code.rindex(ch)
        if (j - i - 1) % 2 == 1:
            return False
    return True

def has_adjacent_repeat(code):
    n = len(code)
    return any(code[i] == code[(i + 1) % n] for i in range(n))  # 폐곡선: 순환 인접 — closed curve: cyclic adjacency

print("[G1] n=3 통과열 후보 (p0 내부 {a1,b1} 순서 × p3 내부 {b2,c2} 순서):")
cases = {
    "a1<b1, b2<c2": "ABCABC",
    "b1<a1, b2<c2": "BACABC",
    "a1<b1, c2<b2": "ABCACB",
    "b1<a1, c2<b2": "BACACB",
}
survivors = []
for k, code in cases.items():
    ok = parity_ok(code)
    print(f"     {k}: {code} → 홀짝 {'통과' if ok else '위반 (실현 불가)'}")
    if ok:
        survivors.append(code)
assert survivors == ["ABCABC"], survivors
print("     ⇒ 유일 생존 = ABCABC (전 쌍 연결(linked) 코드)\n")

print("[G2] ≤2 교차 폐코드 전수:")
small = ["XX", "XXYY", "XYYX", "XYXY"]
for code in small:
    ok = parity_ok(code)
    rep = has_adjacent_repeat(code)
    note = "실현 불가" if not ok else ("인접 반복 → R1 컬 축약 가능" if rep else "!! 반례 !!")
    print(f"     {code}: 홀짝 {'통과' if ok else '위반'} — {note}")
    assert (not ok) or rep
print("     ⇒ 실현 가능한 ≤2교차 다이어그램은 전부 R1 축약 → 자명형 (c(K)≥3)\n")

# ---------- G3: Kauffman 괄호 / Jones ----------
delta = sp.expand(-A**2 - A**-2)

def bracket_2braid(signs):
    """sigma_1^{s1} sigma_1^{s2} sigma_1^{s3} (2-가닥) 폐포의 괄호.
    Bracket of the closure of sigma_1^{s1} sigma_1^{s2} sigma_1^{s3} (2 strands).
    상태합: 각 교차를 identity(|) 또는 e(cup-cap)로 스무딩.
    State sum: smooth each crossing as identity(|) or e(cup-cap).
    TL_2: e^k = delta^(k-1) e, closure(1)=2 루프, closure(e)=1 루프 ⇒ 루프수 = 2(k=0) / k(k≥1).
    TL_2: e^k = delta^(k-1) e, closure(1)=2 loops, closure(e)=1 loop ⇒ loop count = 2 (k=0) / k (k≥1)."""
    tot = sp.Integer(0)
    for sm in itertools.product([0, 1], repeat=3):  # 1 = e-스무딩 개수 — 1 = the number of e-smoothings
        k = sum(sm)
        loops = 2 if k == 0 else k
        w = sp.Integer(1)
        for s, m in zip(signs, sm):
            if s == +1:
                w *= A if m == 0 else A**-1
            else:
                w *= A**-1 if m == 0 else A
        tot += w * delta**(loops - 1)
    return sp.expand(tot)

def jones_2braid(signs):
    wr = sum(signs)                     # 폐포는 평행 방향 ⇒ 교차 부호 = s_i — the closure is parallel-oriented ⇒ crossing sign = s_i
    f = sp.expand((-A**3)**(-wr) * bracket_2braid(signs))
    return sp.expand(sp.powsimp(f.subs(A, t**sp.Rational(-1, 4))))

TRE_L = sp.expand(-t**-4 + t**-3 + t**-1)   # 세잎 (한쪽 손지기) — trefoil (one handedness)
TRE_R = sp.expand(-t**4 + t**3 + t)         # 거울상 — the mirror image

print("[G3] 8개 상하 배정 (sigma_1 부호열) → Jones 다항식:")
for signs in itertools.product([+1, -1], repeat=3):
    V = jones_2braid(signs)
    if V == 1:
        tag = "자명 (V=1)"
    elif V == TRE_L:
        tag = "세잎 V(3_1)"
    elif V == TRE_R:
        tag = "세잎 V(3_1 거울)"
    else:
        tag = f"?? {V}"
    print(f"     s={signs}: {tag}")
    mixed = len(set(signs)) > 1
    assert (V == 1) == mixed, (signs, V)     # 혼합 ⟺ V=1, 상수 ⟺ 세잎 — mixed ⟺ V=1, constant ⟺ trefoil
    if not mixed:
        assert V in (TRE_L, TRE_R)
print("     ⇒ 혼합 부호 6개 전부 V=1 (R2 소거와 정합), 상수 부호 2개 = 좌·우 세잎\n")

# ---------- G4: 삼색성 ---------- — tricolorability
def tricolor_count(signs):
    """폐포 곡선의 6 통과 사건: 첫 하강 j=0,1,2 (위치 1,2,1 진입),
    복귀 j=0,1,2 (위치 2,1,2 진입). sigma_1: 위치-1 가닥이 위 (규약).
    s=+1이면 위치1 진입이 over, s=-1이면 위치2 진입이 over.
    The 6 passage events of the closure curve: the first descent j=0,1,2 (entering at positions 1,2,1),
    the return j=0,1,2 (entering at positions 2,1,2). sigma_1: the position-1 strand is on top (convention).
    If s=+1 the position-1 entry is the over-pass; if s=-1 the position-2 entry is the over-pass."""
    events = []                          # (crossing, is_over) 순환열 — cyclic sequence of (crossing, is_over)
    pos = 1
    for j in range(3):
        events.append((j, (signs[j] == +1) == (pos == 1)))
        pos = 3 - pos
    for j in range(3):
        events.append((j, (signs[j] == +1) == (pos == 1)))
        pos = 3 - pos
    unders = [i for i, (_, o) in enumerate(events) if not o]
    assert len(unders) == 3
    arc_of = {}                          # 사건 → 호 번호 (U 사건에서 절단; U 자체는 경계) — event → arc index (cut at the U events; U itself is a boundary)
    n = len(events)
    for a_idx in range(3):
        start = unders[a_idx]
        end = unders[(a_idx + 1) % 3]
        i = (start + 1) % n
        while True:
            arc_of[i] = a_idx            # start 직후부터 end까지의 O 사건들 + end(U 진입) — the O events from just after start through end, plus end (the U entry)
            if i == end:
                break
            i = (i + 1) % n
    cnt = 0
    for cols in itertools.product(range(3), repeat=3):
        good = True
        for a_idx in range(3):
            u_ev = unders[a_idx]         # 이 U 사건: 들어오는 호 = arc_of[u_ev], 나가는 호 = 다음 run — this U event: incoming arc = arc_of[u_ev], outgoing arc = the next run
            j = events[u_ev][0]
            over_evs = [i for i, (jj, o) in enumerate(events) if jj == j and o]
            c_over = cols[arc_of[over_evs[0]]]
            c_in = cols[arc_of[u_ev]]
            c_out = cols[arc_of[(u_ev + 1) % n]]   # U 직후 사건이 속한 호 = 나가는 호 — the arc containing the event just after U = the outgoing arc
            if (2 * c_over - c_in - c_out) % 3 != 0:
                good = False
                break
        if good:
            cnt += 1
    return cnt

print("[G4] 삼색성 (mod 3 해 개수; 자명형=3):")
for signs in itertools.product([+1, -1], repeat=3):
    c = tricolor_count(signs)
    mixed = len(set(signs)) > 1
    print(f"     s={signs}: {c}")
    assert c == (3 if mixed else 9), (signs, c)
print("     ⇒ 상수 부호 = 9 > 3: 비자명 (초등 인증); 혼합 = 3 (자명과 정합)")

print("\n결론: 통과열 ABCABC 강제 [G1], c(K)≥3 [G2],")
print("      상하 배정은 교대(상수 부호) 2개만 묶임 — 좌·우 세잎, 나머지 6개 자명 [G3,G4].")

# ---------- Exhaustive chord-pairing enumeration and machine certification of shadow uniqueness ----------

def all_chord_codes(c):
    """{1..2c}의 c-현 완전 짝짓기 전수 → 순환 회전 + 재명명 정규형으로 분류.
    All c-chord perfect matchings of {1..2c} → classified by the normal form under cyclic rotation + relabeling."""
    import itertools as it
    pts = list(range(2*c))
    def pairings(rem):
        if not rem: yield []
        else:
            a0 = rem[0]
            for j in range(1, len(rem)):
                for rest in pairings(rem[1:j] + rem[j+1:]):
                    yield [(a0, rem[j])] + rest
    def code_of(pairs):
        lab = {}
        s = []
        nxt = 0
        for p in pts:
            ch = next(i for i, pr in enumerate(pairs) if p in pr)
            if ch not in lab:
                lab[ch] = chr(ord('A') + nxt); nxt += 1
            s.append(lab[ch])
        return "".join(s)
    def canon(code):
        n = len(code)
        best = None
        for r in range(n):
            rot = code[r:] + code[:r]
            lab = {}; nxt = 0; out = []
            for ch in rot:
                if ch not in lab:
                    lab[ch] = chr(ord('A') + nxt); nxt += 1
                out.append(lab[ch])
            cand = "".join(out)
            best = cand if best is None or cand < best else best
        return best
    classes = {}
    for pr in pairings(pts):
        classes.setdefault(canon(code_of(pr)), 0)
        classes[canon(code_of(pr))] += 1
    return sorted(classes)

print("\n[G1'] c=3 순환 코드 전수 (짝짓기 15개 → 정규류):")
cls3 = all_chord_codes(3)
print("     류:", cls3)
ok3 = [c0 for c0 in cls3 if parity_ok(c0) and not has_adjacent_repeat(c0)]
print("     홀짝 통과 ∧ 인접반복 없음:", ok3)
assert ok3 == ["ABCABC"], ok3
print("     ⇒ 실현 가능·비축약 3교차 순환 코드는 ABCABC 유일 (전수)")

print("\n[G2'] c=2 순환 코드 전수:")
cls2 = all_chord_codes(2)
print("     류:", cls2)
for c0 in cls2:
    ok = parity_ok(c0); rep = has_adjacent_repeat(c0)
    assert (not ok) or rep, c0
    print(f"     {c0}: 홀짝 {'통과' if ok else '위반'} / 인접반복 {'유' if rep else '무'}")
print("     ⇒ 실현 가능한 c≤2 코드는 전부 R1 컬 보유 (전수)")

print("\n[G5] ABCABC shadow uniqueness — all 8 rotation systems (Lemma 3.12):")
# 폐곡선 통과열 p0..p5 = A1 B1 C1 A2 B2 C2; 변 e_i: p_i → p_{i+1 (mod 6)} — the passage word of the closed curve p0..p5 = A1 B1 C1 A2 B2 C2; edge e_i: p_i → p_{i+1 (mod 6)}
# 다트: (i,+1)=정방향, (i,-1)=역방향. 정점별 out-다트 4개, 횡단성 ⇒ 회전 2택. — darts: (i,+1)=forward, (i,-1)=backward. 4 out-darts per vertex; transversality ⇒ 2 rotation choices.
VIS = {'A': (0, 3), 'B': (1, 4), 'C': (2, 5)}     # 각 교차의 두 통과 위치 — the two passage positions of each crossing
def out_darts(vname):
    i1, i2 = VIS[vname]
    # 통과 i의 진입 변 = e_{i-1}(역다트로 나감), 진출 변 = e_i(정다트) — the entering edge of passage i = e_{i-1} (leaving as a backward dart), the exiting edge = e_i (a forward dart)
    return dict(out1=(i1, +1), back1=((i1 - 1) % 6, -1),
                out2=(i2, +1), back2=((i2 - 1) % 6, -1))
def rotations(choice):
    """choice[v] ∈ {0,1}: 정점 v의 횡단 회전 2택.
    choice[v] ∈ {0,1}: the 2 transversal rotation choices at vertex v.
    가닥 교대 순서: (out1, X, back1, Y), {X,Y}={out2, back2}.
    Alternating strand order: (out1, X, back1, Y), {X,Y}={out2, back2}."""
    rot = {}
    for v in "ABC":
        d = out_darts(v)
        X, Y = ((d['out2'], d['back2']) if choice[v] == 0 else (d['back2'], d['out2']))
        rot[v] = [d['out1'], X, d['back1'], Y]
    return rot
def vertex_of_dart_head(dart):
    i, s = dart
    p = (i + 1) % 6 if s == +1 else i          # 다트가 도착하는 통과 위치 — the passage position the dart arrives at
    return next(v for v, (a, b) in VIS.items() if p in (a, b))
def faces(rot):
    darts = [(i, s) for i in range(6) for s in (+1, -1)]
    nxt = {}
    for v, cyc in rot.items():
        for k, d in enumerate(cyc):
            nxt[d] = cyc[(k + 1) % 4]
    def phi(d):
        rev = (d[0], -d[1])
        return nxt[rev]
    seen = set(); F = []
    for d0 in darts:
        if d0 in seen: continue
        cyc = []; d = d0
        while d not in seen:
            seen.add(d); cyc.append(d); d = phi(d)
        F.append(len(cyc))
    return sorted(F)
import itertools as _it
planar = []
for bits in _it.product((0, 1), repeat=3):
    choice = dict(zip("ABC", bits))
    Fv = faces(rotations(choice))
    eul = 3 - 6 + len(Fv)
    print(f"     회전택 {bits}: 면 {len(Fv)}개 {Fv}, 오일러 {eul} {'← 평면' if eul == 2 else ''}")
    if eul == 2:
        planar.append((bits, Fv))
assert len(planar) == 2, planar
assert all(Fv == [2, 2, 2, 3, 3] for _, Fv in planar), planar
b0, b1 = planar[0][0], planar[1][0]
assert all(x + y == 1 for x, y in zip(b0, b1)), (b0, b1)   # 전 정점 반전 = 거울쌍 — flipping every vertex = the mirror pair
print("     ⇒ 평면 회전계 정확히 2개, 둘 다 면벡터 (2,2,2,3,3), 서로 전 정점 반전(거울쌍)")
print("     => Lemma 3.12 (shadow uniqueness) machine-checked")
