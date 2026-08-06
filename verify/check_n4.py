# check_n4.py — 프로그램 2단계 (17): n=4 묶음 데이터 열거 1차 — 블록 코드 전수 + 홀짝 필터 — program stage 2 (17): first-pass enumeration of the n=4 tying data — exhaustive block codes + parity filter
# 조각 P0..P4 (P0,P4 = 반직선), 비인접 쌍 6개, 직선성 ⇒ 쌍당 교차 ≤ 1. — pieces P0..P4 (P0,P4 = rays), 6 non-adjacent pairs, straightness ⇒ at most 1 crossing per pair.
import itertools
PAIRS = [(0,2),(0,3),(0,4),(1,3),(1,4),(2,4)]
def parity_ok(code):
    for ch in set(code):
        i, j = code.index(ch), code.rindex(ch)
        if (j - i - 1) % 2 == 1: return False
    return True
def canon(code):
    n=len(code); best=None
    for r in range(n):
        rot=code[r:]+code[:r]
        lab={}; out=[]
        for ch in rot:
            if ch not in lab: lab[ch]=chr(ord('A')+len(lab))
            out.append(lab[ch])
        cand="".join(out)
        best=cand if best is None or cand<best else best
    return best
from collections import defaultdict
summary=defaultdict(set)
for k in range(3,7):
    for sub in itertools.combinations(range(6), k):
        pairs=[PAIRS[i] for i in sub]
        marks=defaultdict(list)          # piece -> [chord ids]
        for cid,(p,q) in enumerate(pairs):
            marks[p].append(cid); marks[q].append(cid)
        piece_orders=[list(itertools.permutations(marks[p])) for p in range(5)]
        for combo in itertools.product(*piece_orders):
            seq=[]
            for p in range(5): seq.extend(combo[p])
            code="".join(chr(ord('A')+c) for c in seq)
            if parity_ok(code):
                summary[k].add(canon(code))
print("n=4: 홀짝-실현가능 순환 코드류 (블록 구조 전수):")
for k in sorted(summary):
    print(f"  교차 {k}개: {len(summary[k])}류")
    if k<=4:
        for c in sorted(summary[k]): print("     ", c)
