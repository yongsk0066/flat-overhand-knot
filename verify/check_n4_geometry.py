"""Focused regression tests for actual crease caps and the numerical oracle."""
import json
import math
from pathlib import Path
from n4_geometry import halfplanes, open_chord_interval
from n4_tighten2 import config, crossings, general_position_ok, valid_fixed2

sign = lambda v: int(v>0)-int(v<0)
# The incoming pass is y in (-.5,.5), x<y, not the rectangle x<0.
fixture = dict(d=[(1.,0.)]*5, M=[(0.,0.)]*4,
               CH=[((-.5,-.5),(.5,.5))]*4, L=[1.,1.,1.])
planes = halfplanes(fixture, 0, sign=sign)
assert open_chord_interval(((.2,.3),(.3,.4)), planes, sign=sign) is not None
assert open_chord_interval(((-.3,-.4),(-.2,-.3)), planes, sign=sign) is None
ray_fixture = dict(d=[(1.,0.),(1.,0.),(1.,0.),(1.,0.),(0.,-1.)],
                   M=[(0.,0.),(0.,1.),(0.,2.),(-100.,1.)], L=[1.,1.,1.])
assert (0,4) in [pair for pair, _, _ in crossings(ray_fixture)]

t = math.acos(-.25)
star = [t, 2*math.pi-2*t, -(2*math.pi-2*t), -t,
        6/math.sqrt(15), 8/math.sqrt(15), 6/math.sqrt(15)]
order = [1,4,2,0,3]
assert valid_fixed2(star, order)
cstar=config(star);xstar=crossings(cstar)
assert general_position_ok(cstar,xstar)
assert not general_position_ok(cstar,xstar+[xstar[0]])
# Independent reviewer witness: M2=(-1,0) lies on the nonadjacent incoming
# ray. crossings() omits that endpoint event; GP must still reject it.
fold_contact=(math.pi/4,2.819842099193151,1.1071487177940904,math.pi/2,
              math.sqrt(2),math.sqrt(5),1.)
cfold=config(fold_contact)
assert not crossings(cfold)
assert not general_position_ok(cfold,crossings(cfold))
v = [0,0,0,0,8*math.sqrt(15)/15,4*math.sqrt(15)/5,8*math.sqrt(15)/15]
for step in [1e-3,1e-4,1e-6]:
    assert valid_fixed2([a+step*b for a,b in zip(star,v)], order)
# The archived harvest includes sampled-oracle starts; report rather than
# relabel it as a fresh exhaustive census.
path = Path(__file__).with_name('n4_valid41.json')
starts = json.loads(path.read_text())
accepted = sum(bool(valid_fixed2(theta, layering)) for theta,layering in starts)
print('PASS: oblique cap counterexamples, exact-state numerical oracle, opening-ray probes')
print(f'Archived sampled starts accepted by actual-panel oracle: {accepted}/{len(starts)} (not a new census)')
