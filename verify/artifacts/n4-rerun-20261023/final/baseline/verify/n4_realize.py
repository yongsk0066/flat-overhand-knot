import math, random, sys
def canon(code):
    n=len(code); best=None
    for r in range(n):
        rot=code[r:]+code[:r]; lab={}; out=[]
        for ch in rot:
            if ch not in lab: lab[ch]=chr(ord('A')+len(lab))
            out.append(lab[ch])
        cand="".join(out)
        best=cand if best is None or cand<best else best
    return best
def diagram(tau,L,ray=80.0):
    cs=[sum(tau[:j]) for j in range(5)]
    d=[(math.cos(c),math.sin(c)) for c in cs]
    M=[(0.0,0.0)]
    for i in range(3): M.append((M[-1][0]+L[i]*d[i+1][0], M[-1][1]+L[i]*d[i+1][1]))
    segs={0:(M[0],(-d[0][0],-d[0][1]),ray),1:(M[0],d[1],L[0]),2:(M[1],d[2],L[1]),3:(M[2],d[3],L[2]),4:(M[3],d[4],ray)}
    X=[]
    for i in range(5):
        for j in range(i+2,5):
            P1,d1,l1=segs[i]; P2,d2,l2=segs[j]
            det=d1[0]*d2[1]-d1[1]*d2[0]
            if abs(det)<1e-12: continue
            rx,ry=P2[0]-P1[0],P2[1]-P1[1]
            s=(rx*d2[1]-ry*d2[0])/det; u=(rx*d1[1]-ry*d1[0])/det
            if 1e-9<s<l1-1e-9 and 1e-9<u<l2-1e-9: X.append(((i,j),s,u))
    return X
def code_of(X):
    if len(X)!=4: return None
    lab={p:chr(ord('A')+k) for k,(p,_,_) in enumerate(X)}
    ev=[]
    for p,s,u in X:
        i,j=p
        ev.append((i, -s if i==0 else s, lab[p]))   # P0: 재료순 = s 내림차순 ⟹ -s 오름 — P0: material order = decreasing s ⟹ increasing -s
        ev.append((j, -u if False else u, lab[p]) if j!=0 else (0,-u,lab[p]))
        if j==4: pass
    ev.sort()
    return "".join(t[2] for t in ev)
rnd=random.Random(77)
hits=[]
N=1200000
for k in range(N):
    rho=rnd.choice([1,-1])
    tau=[rho*rnd.uniform(0.6,3.05) for _ in range(3)]+[rnd.choice([1,-1])*rnd.uniform(0.3,3.05)]
    L=[rnd.uniform(0.2,3.5) for _ in range(3)]
    X=diagram(tau,L)
    if len(X)!=4: continue
    c=code_of(X)
    if c and canon(c)=="ABCADCBD":
        hits.append((tau,L,[p for p,_,_ in X],c))
        if len(hits)>=5: break
print(f"발견 {len(hits)} / 시도 {k+1}")
for tau,L,ps,c in hits[:3]:
    print(" τ=",[round(t,3) for t in tau]," L=",[round(l,3) for l in L]," pairs=",ps," code=",c)
