"""Supplementary n=4 numerical diagnostics on the corrected actual panels.

Code/length margin minima are observed values, not global bounds. Sampled
short descent rays are numerical witnesses, not interval-certified paths.
No absent witness is interpreted as a proof of blocking.
"""
import argparse,hashlib,json,math,random,time
from pathlib import Path
from n4_rerun import DEFAULT,stage_setup,write_json
from n4_tighten2 import config,pieces,crossings,valid_fixed2,ell
from n4_walls import ORDER


def xmargins(th):
    c=config(th)
    if c is None:return None
    segs=pieces(c);out={}
    for (i,j),s,u in crossings(c):
        out[f'({i},{j})s']=s;out[f'({i},{j})u']=u
        if i not in (0,4):out[f'({i},{j})S']=segs[i][2]-s
        if j not in (0,4):out[f'({i},{j})U']=segs[j][2]-u
    for k in (4,5,6):out[f'L{k-3}']=th[k]
    return out


def tuples(value):
    return tuple(map(tuples,value)) if isinstance(value,list) else value


def codewalls(out):
    path=out/'baseline/verify/n4_valid41.json'
    starts=[th for th,o in json.loads(path.read_text()) if th[0]>0][:4]
    assert len(starts)==4 and all(valid_fixed2(th,ORDER) for th in starts)
    keys=sorted(xmargins(starts[0]))
    dest=stage_setup(out,'codewalls',['n4_rerun_diagnostics.py','n4_rerun.py','n4_tighten2.py','n4_geometry.py'],
                     dict(seed=21,n_starts=4,starts=starts,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          max_iterations=3000,trials_per_iteration=20,step_start=.05,step_stop=1e-6,
                          keys=keys,scope='observed numerical minima; not global lower bounds'))
    rnd=random.Random(21);results=[]
    for index,key in enumerate(keys):
        saved=dest/f'margin-{index}.json'
        if saved.exists():
            row=json.loads(saved.read_text());rnd.setstate(tuples(row['rng_state_after']));results.append(row);continue
        started=time.monotonic();runs=[]
        for si,th0 in enumerate(starts):
            th=list(th0);step=.05;cur=xmargins(th)[key];trials=0;accepted=0
            for it in range(3000):
                adv=False
                for _ in range(20):
                    trials+=1;cand=[th[q]+step*rnd.gauss(0,1) for q in range(7)]
                    if not valid_fixed2(cand,ORDER):continue
                    m=xmargins(cand)
                    if m and m[key]<cur-1e-12:
                        th=cand;cur=m[key];adv=True;accepted+=1;break
                if not adv:
                    step*=.6
                    if step<1e-6:break
            assert valid_fixed2(th,ORDER)
            runs.append(dict(start_index=si,theta=th,minimum=cur,iterations=it+1,trials=trials,accepted=accepted,final_step=step))
        row=dict(key=key,minimum=min(r['minimum'] for r in runs),runs=runs,seconds=time.monotonic()-started,rng_state_after=rnd.getstate())
        write_json(saved,row);results.append(row)
        print(f"codewall {index+1}/{len(keys)} {key}: observed_min={row['minimum']:.9f}, {row['seconds']:.1f}s",flush=True)
    output=dict(n_targets=len(keys),n_starts=4,total_trials=sum(r['trials'] for row in results for r in row['runs']),
                minima={r['key']:r['minimum'] for r in results},overall_observed_min=min(r['minimum'] for r in results),
                scope='sampled targeted minimization; no certified lower bound')
    write_json(dest/'result.json',output)
    print(f"COMPLETE codewalls: {len(keys)} targets x4starts; observed minimum={output['overall_observed_min']:.9f}",flush=True)


def stragglers(out):
    import numpy as np
    from scipy.optimize import linprog
    from n4_strata2 import margin_fn,grad,GRAD_ELL
    source=out/'descent/result.json';data=json.loads(source.read_text())
    inv=json.loads((out/'walls/result.json').read_text())['histogram']
    labels=[k for k,v in inv.items() if v>=2 and k.startswith('a:')]
    dest=stage_setup(out,'stragglers',['n4_rerun_diagnostics.py','n4_rerun.py','n4_strata2.py','n4_walls.py','n4_tighten2.py','n4_geometry.py'],
                     dict(input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),labels=labels,
                          active_tolerance=1e-3,lp_bounds=[-1,1],trial_lengths=[1e-3,1e-4,1e-5,1e-6,1e-7],
                          samples_per_ray=101,scope='sampled numerical rays only; no negative Farkas completeness claim'))
    rows=[]
    for endpoint in data['endpoints']:
        th=endpoint['theta'];order=endpoint['order']
        # Top-bottom reversal leaves panel admissibility unchanged. The
        # lateral functions use the positive-turn reference at this point.
        assert valid_fixed2(th,order)
        fns=[(lb,margin_fn(lb,th)) for lb in labels]
        active=[(lb,fn) for lb,fn in fns if abs(fn(th))<1e-3]
        G=np.array([grad(fn,th) for _,fn in active])
        solution=linprog(GRAD_ELL,A_ub=-G if len(G) else None,
                         b_ub=np.zeros(len(G)) if len(G) else None,bounds=[(-1,1)]*7,method='highs')
        row=dict(index=endpoint['index'],core_gap=endpoint['ell']-data['star_core'],
                 active=[lb for lb,_ in active],lp_success=bool(solution.success),status='unresolved_by_this_probe')
        if solution.success and solution.fun < -1e-7:
            direction=solution.x.tolist()
            for length in [1e-3,1e-4,1e-5,1e-6,1e-7]:
                samples=[[a+(length*j/100)*b for a,b in zip(th,direction)] for j in range(101)]
                if all(valid_fixed2(point,order) for point in samples) and ell(samples[-1])<ell(th)-1e-12:
                    row.update(status='sampled_descending_ray_found',direction=direction,length=length,
                               derivative=float(solution.fun),end=samples[-1],core_decrease=ell(th)-ell(samples[-1]),samples=101)
                    break
        rows.append(row)
    far=[r for r in rows if r['core_gap']>=.01]
    output=dict(n_endpoints=len(rows),n_far_stragglers=len(far),
                far_sampled_rays=sum(r['status']=='sampled_descending_ray_found' for r in far),
                all_sampled_rays=sum(r['status']=='sampled_descending_ray_found' for r in rows),rows=rows,
                scope='finite ray sampling does not certify a continuous feasible path; unresolved is not blocked')
    write_json(dest/'result.json',output)
    print(f"COMPLETE stragglers: sampled rays {output['far_sampled_rays']}/{len(far)} far endpoints; no blocking claim for remainder",flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['codewalls','stragglers'],required=True)
    parser.add_argument('--output',type=Path,default=DEFAULT);args=parser.parse_args()
    globals()[args.stage](args.output)

if __name__=='__main__':main()
