"""Direct low-degree BGG checks, separate from trace-invariant quotients.

Reconstruct literal elementary functions on the six A5 weights and compare
all three parabolics. Source completeness is supplied by the separately
replayed coefficient certificates, not inferred from the BGG agreement.
"""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from time import perf_counter
import json,gzip
from sage.all import QQ,ZZ,matrix,prod
from sage_borel import BorelMap

ROOT=Path(__file__).resolve().parents[1]


def main():
    start=perf_counter();cases=[];archives={}
    for kind,omitted,bound in (('A5',(2,),6),('A2A2A1',(4,),6),('A2A2',(2,4),5)):
        data=json.loads((ROOT/('experiments/e6-remaining-'+('a5' if kind=='A5' else 'a2a2a1')+'.json')).read_text())
        if 'certificate_archive' in data:
            archive=ROOT/data['certificate_archive']
            digest=sha256(archive.read_bytes()).hexdigest()
            assert digest==data['certificate_archive_sha256']
            archives[data['certificate_archive']]=digest
            data=json.loads(gzip.decompress(archive.read_bytes()))
        b=BorelMap('E6',omitted);w=b.weights
        xs=[sum(QQ(c)*v for c,v in zip(row,w)) for row in data['branching']['x_coordinates']]
        z=sum(QQ(c)*v for c,v in zip(data['branching']['z_coordinates'],w))
        def elementary(alphabet,k):return sum(prod(f for f in subset) for subset in combinations(alphabet,k))
        if kind=='A5':
            gs=[elementary(xs,k) for k in range(2,7)];cw=1
        else:
            gs=[sum(xs[:3]),elementary(xs[:3],2),prod(xs[:3]),elementary(xs[3:],2),prod(xs[3:])]
            cw=2 if kind=='A2A2A1' else 1
        for g in gs+[z**cw]:b._validate_invariance(g)
        for d in range(bound+1):
            labels=b.basis(d);nc=[];mc=[]
            for k in range(d//cw+1):
                r=data['sources'][d-cw*k];exps=r['exponents'];n=len(exps)
                B=matrix(QQ,r['basis']) if 'basis' in r else matrix(QQ,n,n,{(i,j):QQ(c) for i,j,c in r['basis_sparse']})
                columns=[]
                for e in exps:
                    f=z**(cw*k)*prod(g**a for g,a in zip(gs,e));expanded=b.expand(f)
                    columns.append([expanded.get(v,0) for v in labels])
                A=matrix(QQ,columns).transpose() if columns else matrix(QQ,len(labels),0)
                nc.append(A);mc.append(A*B)
            N=nc[0];M=mc[0]
            for A in nc[1:]:N=N.augment(A)
            for A in mc[1:]:M=M.augment(A)
            den=N.denominator()
            H=matrix(ZZ,den*N).transpose().hermite_form(include_zero_rows=False).transpose()/den
            assert H.is_square() and H.nrows()==len(labels)
            V=H.inverse()*M;assert V.denominator()%3
            vals=[int(a.valuation(3)) for a in matrix(ZZ,V.denominator()*V).elementary_divisors()]
            expected=list(data['degrees'][d]['smith_exponents'])
            if kind=='A2A2' and d:expected+=data['degrees'][d-1]['smith_exponents']
            assert vals==sorted(expected),(kind,d,vals,expected)
            cases.append(dict(case=kind,degree=d,rank=len(labels),
                              multiplicities=[vals.count(a) for a in (1,2,3)]))
            print('BGG',kind,d,cases[-1]['multiplicities'],round(perf_counter()-start,2),flush=True)
    paths=[Path(__file__).resolve(),ROOT/'experiments/e6-remaining-a5.json',ROOT/'experiments/e6-remaining-a2a2a1.json']
    result=dict(result='passed',cases=cases,seconds=round(perf_counter()-start,3),raw_sha256=archives,
                scope='20 direct low-degree BGG comparisons across all three remaining parabolics, not all degrees.',
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'experiments/e6-remaining-bgg-check.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
