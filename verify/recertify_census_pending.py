"""Reproduce the four formerly unresolved census dispositions.

Only the four matching records' verdict fields are updated with --write;
no harvest or unrelated certificate is regenerated. Re-running without
--write checks the same records after their status has been updated.
"""
import argparse
import copy
import json
from pathlib import Path
import uniq_certify2 as certificate

TARGETS = {
    ('W1','W4','W10'),
    ('W1','W4','W5','W10'),
    ('W1','W4','W6','W10'),
    ('W1','W4','W8','W10'),
}
VERDICT = 'W7<0 참-벽 구간 인증 (재인증)'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    path = Path(__file__).with_name('uniq_certified.json')
    recorded = json.loads(path.read_text())
    original = copy.deepcopy(recorded)
    census = json.loads(path.with_name('uniq_census.json').read_text())['census']
    checked = []
    for index, row in enumerate(recorded):
        if tuple(row['S']) not in TARGETS:
            continue
        if row['verdict'] != VERDICT and not row['verdict'].startswith('미결'):
            continue
        matches = [r for r in census if r['S'] == row['S'] and r['theta'] == row['theta']]
        assert len(matches) == 1, row
        r = matches[0]
        active = tuple(certificate.NAMES.index(name) for name in r['S'])
        box = certificate.krawczyk(active, r['theta']+r['lam'])
        assert box is not None, ('Krawczyk failed', row)
        assert certificate.true_wall_sign_on_box('W7', box) == '-', row
        row['verdict'] = VERDICT
        checked.append(index)
        print(f'PASS: {r["S"]}, x={r["theta"][0]:.17g}: isolated root and true W7<0', flush=True)
    assert len(checked) == 4
    assert {tuple(recorded[i]['S']) for i in checked} == TARGETS
    for i, (old, new) in enumerate(zip(original, recorded)):
        if i not in checked:
            assert old == new
        else:
            assert {k:v for k,v in old.items() if k != 'verdict'} == {k:v for k,v in new.items() if k != 'verdict'}
    if args.write:
        path.write_text(json.dumps(recorded, indent=1)+'\n')
    print(f'PASS: all four certificates reproduced; updated={args.write}')


if __name__ == '__main__':
    main()
