"""Check the closed SL3 literal K-image Smith formulas from actual characters.

The formulas use the source homothety and cubic trace index. The audit
instead forms complete Laurent Demazure image matrices on every proper
SL3 flag for each selected actual representation.
"""
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from sage.all import QQ,ZZ,matrix
from sage_borel import BorelMap,SpringerCharacterMap
from other_theories_k_pilot import run,euler

ROOT=Path(__file__).resolve().parents[1]


def classify(highest):
    ch=SpringerCharacterMap('A2',highest);B=ch.linear_basis_matrix
    s=0
    for x in B.list():s=gcd(s,int(x))
    kind='P' if abs(B.det())==s*s else 'Q'
    assert abs(B.det())==(s*s if kind=='P' else 3*s*s)
    I=sum(m*(p[0]+2*p[1])**2 for p,m in ch.weight_multiplicities)/QQ(6)
    A=-sum(m*(p[0]+2*p[1])**3 for p,m in ch.weight_multiplicities)/QQ(6)
    c=s*(1-A/I)/2
    assert c.denominator()==1
    c=int(c);g=gcd(s,c)
    expected=[1,g,g,s**3//g,s**3//g,s**3] if kind=='P' else [1,g,3*g,3*s**3//g,3*s**3//g,3*s**3]
    actual,_=run('A2',highest)
    assert actual['K_smith_factors']==expected
    projective=[1,g,s**3//g] if kind=='P' else [1,3*g,3*s**3//g]
    roots=[(2,-1),(-1,2)];x,y=ch.ring.gens()
    for node in (1,2):
        b=BorelMap('A2',omitted_nodes=node)
        words=[tuple(w.reduced_word()) for d in range(3) for w in b.basis(d)]
        h=x if node==1 else y
        q=(2*x-y)*(x+y) if node==1 else (2*y-x)*(x+y)
        fs=[ch.ring.one(),s*h,s*s*h*h] if kind=='P' else [ch.ring.one(),3*s*h,s*s*q]
        columns=[]
        for f in fs:
            terms={tuple(map(int,p)):v for p,v in ch.pullback(f).dict().items()}
            assert all(v.denominator()==1 for v in terms.values())
            columns.append([euler(terms,w,roots) for w in words])
        M=matrix(ZZ,columns).transpose()
        assert [int(a) for a in M.elementary_divisors()]==projective
    return dict(highest_weights=[list(w) for w in highest],dimension=ch.dimension,
                source_form=kind,source_scale=s,trace_index=str(I),cubic_index=str(A),
                correction_parameter=c,gcd_parameter=g,
                full_flag_smith_factors=expected,projective_smith_factors=projective,
                actual_character_Demazure_matrices_agree=True)


def main():
    highest=[((i,j),) for i in range(6) for j in range(6) if i+j]
    highest += [((1,0),)*4+((2,0),),((1,0),)*4+((0,2),),
                ((1,0),(0,1)),((1,0),(1,1)),((2,0),(0,2)),((1,1),(1,1))]
    cases=[classify(h) for h in highest]
    report=dict(result='passed',date='2026-10-22',cases=cases,
                representation_cases=len(cases),proper_flag_images_checked=3*len(cases),
                scope='Finite complete actual-character image checks corroborate the separately written all-representation SL3 proof; cubic indices come from moments independently of cotangent recursion.',
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in ('experiments/other_theories_sl3_classification.py','experiments/other_theories_k_pilot.py')})
    (ROOT/'experiments/other-theories-sl3-classification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED',len(cases),'representations on all three proper SL3 flags',flush=True)


if __name__=='__main__':main()
