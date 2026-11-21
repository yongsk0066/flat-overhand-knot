"""Reproducible full-size n=4 reruns, with separate provenance/checkpoints.

Historical JSON files are never overwritten. Eight original random streams
(seed7000..7007, one million configurations each) are scheduling independent.
All numerical outcomes remain numerical evidence, not exhaustive proofs.
"""
import argparse, collections, datetime, hashlib, json, math, multiprocessing as mp
import platform, shutil, subprocess, sys, time
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
ARCHIVED = ROOT / 'verify/artifacts/n4-rerun-20261023'
DEFAULT = ROOT / 'build/n4-reproduction'


def write_json(path,data):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    tmp.replace(path)


def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)


def _source_revision():
    if not (ROOT / '.git').exists():
        return 'archive-no-git'
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, stderr=subprocess.DEVNULL, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return 'archive-no-git'


def _lock_digest():
    lock = ROOT / 'uv.lock'
    return hashlib.sha256(lock.read_bytes()).hexdigest() if lock.exists() else None


def stage_setup(out,stage,names,config):
    dest = (out / stage).resolve()
    _require(not dest.is_relative_to(ARCHIVED.resolve()), "choose a new output directory outside the archived run")
    dest.mkdir(parents=True, exist_ok=True)
    # Normalize tuples before comparing a resumed JSON receipt.
    config=json.loads(json.dumps(config))
    sources={name:hashlib.sha256((ROOT/'verify'/name).read_bytes()).hexdigest() for name in names}
    lock_sha=_lock_digest()
    meta=dict(stage=stage,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              python=sys.version,executable=sys.executable,platform=platform.platform(),
              git_head=_source_revision(),uv_lock_sha256=lock_sha,
              sources=sources,config=config)
    old=dest/'provenance.json'
    if old.exists():
        previous=json.loads(old.read_text())
        _require(previous['sources']==sources and previous['config']==config, 'resume input mismatch')
        _require(previous.get('uv_lock_sha256') is not None and previous.get('uv_lock_sha256')==lock_sha,
                 'resume dependency-lock mismatch (use a fresh --output directory)')
        _require(previous.get('python')==sys.version,
                 'resume python-version mismatch (use a fresh --output directory)')
    else:
        _require(not any(dest.iterdir()), "orphan stage files without provenance: use a fresh output directory")
        for name in names:
            copy=dest/'source'/name;copy.parent.mkdir(exist_ok=True);shutil.copy2(ROOT/'verify'/name,copy)
        write_json(old,meta)
    return dest


def sweep_worker(job):
    from n4_census import _hworker
    index,seed,count=job; started=time.monotonic()
    shadow,cells,hits=_hworker((seed,count))
    return dict(index=index,seed=seed,samples=count,seconds=time.monotonic()-started,
                shadows=[dict(pairs=p,count=n) for p,n in sorted(shadow.items())],
                cells=[dict(pairs=p,order=o,count=n) for (p,o),n in sorted(cells.items())],
                hits=[dict(theta=th,order=o) for th,o in hits])


def sweep(out,workers):
    jobs=[(i,7000+i,1_000_000) for i in range(8)]
    dest=stage_setup(out,'sweep',['n4_census.py','n4_tighten2.py','n4_geometry.py'],
                     dict(total_samples=8_000_000,jobs=jobs,workers=workers,gp_tolerance=1e-9,
                          distance_tolerance=1e-12,parameter_tolerance=1e-11))
    todo=[j for j in jobs if not (dest/f'stream-{j[0]}.json').exists()]
    with mp.get_context('fork').Pool(workers) as pool:
        for result in pool.imap_unordered(sweep_worker,todo):
            write_json(dest/f'stream-{result["index"]}.json',result)
            print(f"sweep stream {result['index']}: {result['samples']:,} samples, {len(result['hits'])} state/layer records, {result['seconds']:.1f}s",flush=True)
    streams=[json.loads((dest/f'stream-{i}.json').read_text()) for i in range(8)]
    assert sum(r['samples'] for r in streams)==8_000_000
    shadow=collections.Counter();cells=collections.Counter();hits=[]
    for row in streams:
        for r in row['shadows']: shadow[tuple(map(tuple,r['pairs']))]+=r['count']
        for r in row['cells']: cells[(tuple(map(tuple,r['pairs'])),tuple(r['order']))]+=r['count']
        hits.extend(row['hits'])
    unique=len({tuple(h['theta']) for h in hits})
    result=dict(total_samples=8_000_000,distinct_states=unique,state_layer_records=len(hits),
                shadows=[dict(pairs=p,count=n) for p,n in sorted(shadow.items())],
                cells=[dict(pairs=p,order=o,count=n) for (p,o),n in sorted(cells.items())],hits=hits)
    write_json(dest/'result.json',result)
    print(f"COMPLETE sweep: 8,000,000 samples; {unique} distinct states; {len(hits)} state/layer records; {len(cells)} cells with layer order",flush=True)


def descent_worker(job):
    from n4_census import _dworker
    from n4_tighten2 import valid_fixed2
    i,hit=job;started=time.monotonic()
    assert valid_fixed2(hit['theta'],hit['order'])
    ell,theta,order=_dworker((9000+i,tuple(hit['theta']),hit['order']))
    assert valid_fixed2(theta,order)
    return dict(index=i,seed=9000+i,theta=theta,ell=ell,order=order,
                start=hit['theta'],seconds=time.monotonic()-started)


def descent(out,workers):
    source=out/'sweep/result.json';data=json.loads(source.read_text());hits=data['hits']
    dest=stage_setup(out,'descent',['n4_census.py','n4_tighten2.py','n4_geometry.py'],
                     dict(inputs_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                          n_starts=len(hits),all_starts=True,rounds=[8,10],steps_per_round=2500,
                          trials_per_step=30,seed_start=9000,workers=workers))
    jobs=[(i,h) for i,h in enumerate(hits) if not (dest/f'endpoint-{i}.json').exists()]
    with mp.get_context('fork').Pool(workers) as pool:
        for r in pool.imap_unordered(descent_worker,jobs):
            write_json(dest/f'endpoint-{r["index"]}.json',r)
            print(f"descent {r['index']+1}/{len(hits)}: core={r['ell']:.9f}; {r['seconds']:.1f}s",flush=True)
    results=[json.loads((dest/f'endpoint-{i}.json').read_text()) for i in range(len(hits))]
    star=4*math.sqrt(15)/3
    output=dict(n_starts=len(hits),n_finished=len(results),star_core=star,
                within_0_01=sum(abs(r['ell']-star)<.01 for r in results),
                clusters=dict(collections.Counter(f"{r['ell']:.2f}" for r in results)),endpoints=results)
    write_json(dest/'result.json',output)
    print(f"COMPLETE descent: {len(results)}/{len(hits)}; within0.01 of star={output['within_0_01']}",flush=True)


def walls(out,workers):
    from n4_walls import boundary_probe
    baseline=out/'baseline/verify'
    starts=[tuple(th) for th,o in json.loads((baseline/'n4_valid41.json').read_text())]
    starts += [tuple(r['theta']) for r in json.loads((baseline/'n4_census.json').read_text())['endpoints']]
    input_hashes={p:hashlib.sha256((baseline/p).read_bytes()).hexdigest() for p in ['n4_valid41.json','n4_census.json']}
    dest=stage_setup(out,'walls',['n4_walls.py','n4_tighten2.py','n4_geometry.py'],
                     dict(requested_attempts=4000,original_seed=11,inputs=input_hashes,
                          starts=starts,bisections=45,direction_limit=3.,near_tolerance=5e-4))
    started=time.monotonic()
    hits,co,records,skipped=boundary_probe(4000,11,starts=starts,return_records=True)
    output=dict(requested_attempts=4000,actual_attempts=4000,reached_boundaries=len(records),
                classified=sum(n for k,n in hits.items() if k!='(미분류)'),histogram=dict(hits),
                simultaneous=[dict(labels=k,count=n) for k,n in co.most_common()],
                skipped=skipped,records=records,seconds=time.monotonic()-started)
    write_json(dest/'result.json',output)
    print(f"COMPLETE walls: {len(records)}/4000 boundaries; {output['classified']} classified; {len(hits)} labels; skipped={skipped}",flush=True)
    print(json.dumps(dict(hits),ensure_ascii=False),flush=True)
    assert len(records)>=3995, 'Below historical3995 reached boundaries; supplement required explicitly'


def strata_worker(job):
    from n4_strata2 import solve_stratum,actual_closure_ok
    index,labels,seeds,ref=job;started=time.monotonic();solutions=[];attempts=[]
    for si,seed in enumerate(seeds):
        t=time.monotonic();r=solve_stratum(labels,seed,ref)
        attempts.append(dict(seed_index=si,converged=r is not None,seconds=time.monotonic()-t))
        if r is None:continue
        theta,lam,residual=r;valid=actual_closure_ok(theta)
        solutions.append(dict(seed_index=si,theta=theta,multipliers=lam,residual=residual,
                              actual_approximate_closure=valid,nonnegative_multipliers=min(lam)>=-1e-8,
                              star_distance_squared=sum((a-b)**2 for a,b in zip(theta,ref))))
    return dict(index=index,labels=labels,attempts=attempts,solutions=solutions,seconds=time.monotonic()-started)


def strata(out,workers):
    import itertools
    from n4_strata2 import actual_closure_ok
    source=out/'walls/result.json';inventory=json.loads(source.read_text())['histogram']
    labels=[k for k,v in sorted(inventory.items(),key=lambda p:-p[1]) if v>=2 and k!='(미분류)']
    old=json.loads((out/'baseline/verify/n4_wall_inventory.json').read_text())
    # Original8-family subset is tracked separately even if new faces occur.
    original=[k for k,v in sorted(old.items(),key=lambda p:-p[1]) if v>=2 and k!='(미분류)']
    labels=original+[k for k in labels if k not in original]
    t=math.acos(-.25);ref=[t,2*math.pi-2*t,-(2*math.pi-2*t),-t,6/math.sqrt(15),8/math.sqrt(15),6/math.sqrt(15)]
    archived=json.loads((out/'baseline/verify/n4_valid41.json').read_text())
    seeds=[ref]+[th for th,o in archived if th[0]>0][:6]
    combinations=[list(s) for r in range(1,8) for s in itertools.combinations(labels,r)]
    assert len([s for s in combinations if set(s)<=set(original)])==254
    dest=stage_setup(out,'strata',['n4_strata2.py','n4_walls.py','n4_tighten2.py','n4_geometry.py'],
                     dict(walls=labels,original_walls=original,n_strata=len(combinations),seeds=seeds,
                          seed_actual_closure=[actual_closure_ok(s) for s in seeds],workers=workers,
                          newton_iterations=60,residual_tolerance=5e-8,geometry_tolerance=1e-6,
                          multiplier_tolerance=1e-8,scope='sampled smooth contact systems; no root completeness claim'))
    jobs=[(i,lb,seeds,ref) for i,lb in enumerate(combinations) if not (dest/f'stratum-{i}.json').exists()]
    with mp.get_context('fork').Pool(workers) as pool:
        for count,r in enumerate(pool.imap_unordered(strata_worker,jobs),1):
            write_json(dest/f'stratum-{r["index"]}.json',r)
            if count%20==0 or count==len(jobs):
                print(f"strata checkpoint {count}/{len(jobs)} remaining jobs",flush=True)
    rows=[json.loads((dest/f'stratum-{i}.json').read_text()) for i in range(len(combinations))]
    solutions=[dict(labels=row['labels'],**sol) for row in rows for sol in row['solutions']]
    passed=[r for r in solutions if r['actual_approximate_closure'] and r['nonnegative_multipliers']]
    rejected_geometry=[r for r in solutions if not r['actual_approximate_closure']]
    negative=[r for r in solutions if r['actual_approximate_closure'] and not r['nonnegative_multipliers']]
    others=[r for r in passed if r['star_distance_squared']>=1e-6]
    output=dict(n_strata=len(rows),n_attempts=sum(len(r['attempts']) for r in rows),n_seeds=len(seeds),
                n_solutions=len(solutions),n_pass=len(passed),n_rejected_geometry=len(rejected_geometry),
                n_negative_multipliers=len(negative),n_other_pass=len(others),others=others,
                solutions=solutions,walls=labels,scope='numerical KKT harvest; no exhaustive root or full-cell optimality proof')
    write_json(dest/'result.json',output)
    print(f"COMPLETE strata: {len(rows)} systems x {len(seeds)} seeds; roots={len(solutions)}, actual-closure+lambda-pass={len(passed)}, other-pass={len(others)}",flush=True)


def cells(out,workers):
    import contextlib,io,runpy
    from n4_cells_probe import CELLS,_work
    dest=stage_setup(out,'cells',['n4_cells3.py','n4_fig8.py','n4_cells_probe.py','n4_tighten2.py','n4_geometry.py'],
                     dict(seeds=list(range(1000,1008)),restarts_per_cell=60,steps_per_restart=4000,workers=workers))
    for script in ['n4_cells3.py','n4_fig8.py']:
        buf=io.StringIO()
        with contextlib.redirect_stdout(buf):runpy.run_path(str(ROOT/'verify'/script))
        (dest/(script+'.log')).write_text(buf.getvalue())
    jobs=[(cell,1000+i,60,4000) for i,cell in enumerate(CELLS)]
    results=[]
    with mp.get_context('fork').Pool(workers) as pool:
        for cell,best,sh,vh,theta in pool.imap_unordered(_work,jobs):
            row=dict(cell=cell,best_mismatch=best,shadow_hits=sh,valid_hits=vh,theta=theta)
            results.append(row);write_json(dest/f'cell-{len(results)}.json',row)
            print(f'cell {cell}: shadow={sh}, valid={vh}',flush=True)
    write_json(dest/'result.json',dict(cells=results,restarts=480,max_steps_per_restart=4000))
    print('COMPLETE cells: all8 targetedcells,60restarts each',flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['sweep','descent','walls','strata','cells'],required=True)
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--output',type=Path,default=DEFAULT)
    args=parser.parse_args();assert 1<=args.workers<=4
    globals()[args.stage](args.output,args.workers)

if __name__=='__main__': main()
