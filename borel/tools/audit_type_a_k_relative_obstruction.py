"""Independent polynomial and dual-functional replay of exact order five."""
from functools import lru_cache
from hashlib import sha256
from math import comb,prod
from pathlib import Path
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(None)
def powers(weights,n):
    if not weights:
        return ((),) if n==0 else ()
    return tuple((i,)+tail for i in range(n//weights[0]+1)
                 for tail in powers(weights[1:],n-i*weights[0]))


def main():
    name='experiments/type-a-k-relative-five-obstruction.json'
    data=json.loads((ROOT/name).read_text())
    d,k=data['weight'],data['coefficient_degree']
    assert (d,k,data['order'])==(12,3,5)
    M,L,*aa=s.symbols('M L '+' '.join('a%d'%i for i in range(1,d+1)))
    q=[M]+[M*L**j+sum(comb(j,i)*aa[i-1]*L**(j-i) for i in range(1,j+1))
           for j in range(1,d+1)]
    T=s.expand(sum(c*(q[i]*q[j]*q[l]-M**2*q[d]) for c,i,j,l in data['source_triples']))
    identity=0
    for item in data['ideal_witness']:
        e=item['word']
        assert len(e)==2*d and sum(e)<=k and any(e[:d])
        assert sum(i*(e[i-1]+e[d+i-1]) for i in range(1,d+1))==d
        identity+=item['coefficient']*M**(k-sum(e))*prod(g**t for g,t in zip(aa+q[1:],e))
    assert s.expand(identity-5*T)==0
    functional={tuple(item['powers']):s.Rational(item['numerator'],item['denominator'])
                for item in data['functional']}
    def evaluate(poly):
        terms=s.Poly(poly.subs(M,1),L,*aa).terms()
        return sum(c*functional.get(e,0) for e,c in terms)
    count=0
    for e in powers(tuple(list(range(1,d+1))*2),d):
        if sum(e)>k or not any(e[:d]):
            continue
        value=evaluate(prod(g**t for g,t in zip(aa+q[1:],e)))
        assert s.denom(value)==1
        count+=1
    value=evaluate(T)
    assert s.denom(value)==5 and value==s.Rational(data['functional_value'])
    names=['TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md','tools/audit_type_a_k_relative_obstruction.py',name]
    report=dict(result='passed',exact_relative_order=5,weight=d,coefficient_degree=k,
                literal_integral_ideal_identity=True,all_ideal_generators_checked=count,
                functional_value=str(value),sage_imported=False,smith_form_used=False,
                actual_elementary_claim=False,
                source_sha256={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for n in names})
    (ROOT/'type-a-k-relative-obstruction-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key!='source_sha256'}))


if __name__=='__main__':
    main()
