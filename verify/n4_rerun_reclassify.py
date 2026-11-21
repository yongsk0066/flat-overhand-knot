"""Apply the final actual-panel/GP predicate to completed raw KKT roots.

The Newton solves do not use geometric validity, so source changes confined
to that predicate need no duplicate numerical solve. Both source snapshots
and the unfiltered harvest are retained with separate provenance.
"""
import argparse,hashlib,json
from pathlib import Path
from n4_rerun import DEFAULT,stage_setup,write_json
from n4_strata2 import actual_closure_ok
from n4_tighten2 import valid_fixed2
from n4_walls import ORDER


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--harvest',type=Path,default=DEFAULT/'strata')
    parser.add_argument('--output',type=Path,default=DEFAULT/'final');args=parser.parse_args()
    source=args.harvest/'result.json';data=json.loads(source.read_text())
    meta=args.harvest/'provenance.json';harvest_meta=json.loads(meta.read_text())
    assert data['n_strata']==254 and data['n_attempts']==254*7
    dest=stage_setup(args.output,'strata_filter',['n4_rerun_reclassify.py','n4_rerun.py','n4_strata2.py','n4_walls.py','n4_tighten2.py','n4_geometry.py'],
                     dict(harvest_path=str(args.harvest.resolve()),harvest_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                          harvest_provenance_sha256=hashlib.sha256(meta.read_bytes()).hexdigest(),
                          n_strata=data['n_strata'],n_attempts=data['n_attempts'],geometry_tolerance=1e-6,
                          parameter_tolerance=1e-8,gp_tolerance=1e-9,multiplier_tolerance=1e-8,
                          justification='Newton equation solve and seeds are independent of geometry filter'))
    rows=[]
    for source_row in data['solutions']:
        row=dict(source_row);row['previous_approximate_closure']=row['actual_approximate_closure']
        row['actual_approximate_closure']=actual_closure_ok(row['theta'])
        row['strict_numeric_oracle_pass']=valid_fixed2(row['theta'],ORDER) is not None
        rows.append(row)
    passed=[r for r in rows if r['actual_approximate_closure'] and r['nonnegative_multipliers']]
    others=[r for r in passed if r['star_distance_squared']>=1e-6]
    output=dict(n_strata=data['n_strata'],n_attempts=data['n_attempts'],n_seeds=7,n_solutions=len(rows),
                n_pass=len(passed),n_other_pass=len(others),
                n_rejected_geometry=sum(not r['actual_approximate_closure'] for r in rows),
                n_negative_multipliers=sum(r['actual_approximate_closure'] and not r['nonnegative_multipliers'] for r in rows),
                solutions=rows,others=others,harvest_sources=harvest_meta['sources'],
                scope='numerical harvest and approximate-closure filter; no root completeness or full-cell optimality proof')
    write_json(dest/'result.json',output)
    print(f"COMPLETE final strata filter: {data['n_strata']}x7 solves, roots={len(rows)}, pass={len(passed)}, other-pass={len(others)}")

if __name__=='__main__':main()
