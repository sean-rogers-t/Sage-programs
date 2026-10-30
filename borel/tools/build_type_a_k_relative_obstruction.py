"""Exact order-five negative control for overbroad relative saturation claims.

Exports an integral ideal identity and a rational separating functional.
The independent replay uses neither Smith form nor the source matrix.
"""
from hashlib import sha256
from math import comb
from pathlib import Path
import json
from sage.all import ZZ, QQ, PolynomialRing, matrix, vector
from verify_type_a_k_stable_multiplicities import exponents

ROOT = Path(__file__).resolve().parents[1]


def main():
    d,k = 12,3
    R=PolynomialRing(ZZ,names=['L']+['a%d'%i for i in range(1,d+1)])
    L,*aa=R.gens()
    q=[R(1)]+[L**j+sum(comb(j,i)*aa[i-1]*L**(j-i)
                       for i in range(1,j+1)) for j in range(1,d+1)]
    basis=[e for e in exponents([1]+list(range(1,d+1)),d) if sum(e[1:])<=k]
    index={e:i for i,e in enumerate(basis)}
    def row(f):
        out=[ZZ(0)]*len(basis)
        for powers,c in f.dict().items():
            out[index[tuple(powers)]]=c
        return vector(ZZ,out)
    words=[e for e in exponents(list(range(1,d+1))*2,d)
           if sum(e)<=k and any(e[:d])]
    C=matrix(ZZ,[row(R.prod(g**t for g,t in zip(aa+q[1:],e))) for e in words])
    terms=[(1,1,2,9),(12,2,3,7),(21,2,4,6),(-63,2,5,5),
           (-84,3,3,6),(252,3,4,5),(-140,4,4,4)]
    T=row(sum(c*(q[i]*q[j]*q[l]-q[d]) for c,i,j,l in terms))
    D,U,V=C.smith_form()
    y=T*V
    rank=C.rank()
    assert all(x==0 for x in y[rank:])
    ratios=[QQ(y[i])/D[i,i] for i in range(rank)]
    assert all((5*x).denominator()==1 for x in ratios)
    pivot=next(i for i,x in enumerate(ratios) if x.denominator()==5)
    z=vector(ZZ,[ZZ(5*x) for x in ratios]+[0]*(C.nrows()-rank))
    witness=z*U
    assert witness*C==5*T
    functional=V.column(pivot)/QQ(D[pivot,pivot])
    assert all(x.denominator()==1 for x in C*functional)
    assert (T*functional).denominator()==5
    certificate=dict(weight=d,coefficient_degree=k,order=5,
        source_triples=terms,
        ideal_witness=[dict(coefficient=int(c),word=list(e))
                       for c,e in zip(witness,words) if c],
        functional=[dict(powers=list(e),numerator=int(c.numerator()),denominator=int(c.denominator()))
                    for e,c in zip(basis,functional) if c],
        functional_value=str(T*functional),
        homogenization='M power = coefficient degree minus number of a and q factors')
    path=ROOT/'experiments/type-a-k-relative-five-obstruction.json'
    path.write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8')
    names=['TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md','tools/build_type_a_k_relative_obstruction.py',
           'tools/verify_type_a_k_stable_multiplicities.py',str(path.relative_to(ROOT)).replace('\\','/')]
    report=dict(result='passed',weight=d,coefficient_degree=k,exact_relative_order=5,
                ideal_generators_checked=C.nrows(),ideal_identity_terms=len(certificate['ideal_witness']),
                functional_value=certificate['functional_value'],
                all_degree_saturation_claim_refuted=True,actual_elementary_claim=False,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in names})
    (ROOT/'type-a-k-relative-obstruction-build.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key!='source_sha256'}))


if __name__=='__main__':
    main()
