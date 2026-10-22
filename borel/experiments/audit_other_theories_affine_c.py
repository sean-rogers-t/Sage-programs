"""Independent affine type C divisor contents via actual real-coroot Chevalley.

This does not use Schur P coproducts or the weighted-multinomial recurrence.
Only degree-bounded affine Weyl groups are enumerated.
"""
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from sage.all import RootSystem

ROOT=Path(__file__).resolve().parents[1]


def run(n,maximum):
    R=RootSystem(['C',n,1]);lattice=R.root_lattice();W=lattice.weyl_group()
    simple={i:W.simple_reflections()[i] for i in range(n+1)};one=W.one();levels=[[one]];seen={one}
    for d in range(maximum):
        target=[]
        for w in levels[-1]:
            for s in simple.values():
                v=w*s
                if v not in seen:seen.add(v);target.append(v)
        assert all(w.length()==d+1 for w in target)
        levels.append(target)
    coroots=R.coroot_lattice().simple_roots();reflections={}
    for w in seen:
        word=w.reduced_word()
        for i,s in simple.items():
            a=coroots[i]
            for j in reversed(word):a=a.simple_reflection(j)
            coefficients=[int(a.monomial_coefficients().get(j,0)) for j in range(n+1)]
            if max(coefficients)<=0:coefficients=[-x for x in coefficients]
            assert min(coefficients)>=0
            reflection=w*s*w.inverse()
            if reflection in reflections:assert reflections[reflection]==coefficients
            else:reflections[reflection]=coefficients
    grass=[[w for w in level if all((w*simple[i]).length()>w.length() for i in range(1,n+1))] for level in levels]
    coefficients={one:1};rows=[]
    for d in range(1,maximum+1):
        target={}
        for v in grass[d]:
            total=0
            for w,c in coefficients.items():
                a=reflections.get(w.inverse()*v)
                if a is not None:total+=c*a[0]
            if total:target[v]=total
        content=0
        for c in target.values():content=gcd(content,c)
        coefficients=target
        rows.append(dict(degree=d,Schubert_rows=len(grass[d]),nonzero_divisor_rows=len(target),
                         content=str(content)))
    return dict(rank=n,degree_through=maximum,affine_Weyl_elements_enumerated=len(seen),degrees=rows)


def main():
    primary=json.loads((ROOT/'experiments/other-theories-affine-pilot.json').read_text())
    cases=[]
    for n,maximum in ((2,20),(3,12)):
        case=run(n,maximum);p=next(c for c in primary['type_C_content_cases'] if c['rank']==n)
        assert [row['content'] for row in case['degrees']]==p['initial_contents'][1:maximum+1]
        cases.append(case);print('affine C',n,'Chevalley audit through degree',maximum,'passed',flush=True)
    report=dict(result='passed',date='2026-10-22',cases=cases,
                independent_method='Real-coroot affine Chevalley multiplication on Grassmannian minimal cosets; no Schur P coproduct or factorial formula used',
                source_sha256={'experiments/audit_other_theories_affine_c.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()},
                raw_sha256={'experiments/other-theories-affine-pilot.json':sha256((ROOT/'experiments/other-theories-affine-pilot.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-affine-c-audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
