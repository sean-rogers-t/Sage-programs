"""Full integral image closure, not a source-degree cutoff.

Run with sage -python. Uses Singular strong standard bases over ZZ.
The symmetric coefficient ring base-changes freely to the actual torus.
"""
from pathlib import Path
from math import comb
from hashlib import sha256
import json
from sage.all import ZZ, PolynomialRing, matrix, vector, singular

ROOT=Path(__file__).resolve().parents[1]

def run(n,order='dp'):
    S=PolynomialRing(ZZ,names=['c'+str(i) for i in range(1,n)])
    cs=[S(1)]+list(S.gens())+[S(1)]
    P=PolynomialRing(S,'L'); L=P.gen()
    f=sum((-1)**j*cs[j]*L**(n-j) for j in range(n+1))
    fs=[sum((-1)**k*comb(n-k,j-k)*cs[k]*L**(j-k)
            for k in range(j+1)) for j in range(1,n)]
    cols=[vector(S,[(g%f)[i] for i in range(n)]) for g in [P(1)]+fs]
    mult=[matrix(S,n,n,lambda i,k: ((g*L**k)%f)[i]) for g in fs]
    singular.eval('ring r=integer,('+','.join(S.variable_names())+'),'+order+';')
    def literal(M):
        return ','.join(str(a).replace('**','^') for a in M.list())
    initial=matrix(S,cols).transpose()
    singular.eval('matrix M[%d][%d]=%s; module G=std(module(M));'%(n,n,literal(initial)))
    for j,A in enumerate(mult):
        singular.eval('matrix T%d[%d][%d]=%s;'%(j,n,n,literal(A)))
    rounds=0
    while True:
        rounds+=1
        singular.eval('module extra;')
        for j in range(n-1):
            singular.eval('extra=extra,reduce(module(T%d*matrix(G)),G);'%j)
        remains=int(singular.eval('size(extra);'))
        if remains==0:break
        singular.eval('G=std(G,extra); kill extra;')
        if rounds>=30:raise RuntimeError('No closure certificate within 30 rounds')
    B=singular('matrix(G)').sage()
    # Independent exact multiplication-closure remainder certificate in Singular.
    for j in range(n-1):
        assert int(singular.eval('size(reduce(module(T%d*matrix(G)),G));'%j))==0
    singular.eval('module relations=syz(G); matrix rel=matrix(relations);')
    rel=singular('rel').sage()
    assert B*rel==0
    vals={S.gen(i-1):ZZ(comb(n,i)) for i in range(1,n)}
    ev=lambda a: ZZ(S(a).subs(vals))
    B0=matrix(ZZ,B.nrows(),B.ncols(),[ev(x) for x in B.list()])
    R0=matrix(ZZ,rel.nrows(),rel.ncols(),[ev(x) for x in rel.list()])
    elementary=list(R0.elementary_divisors())
    torsion=[int(abs(x)) for x in elementary if abs(x)>1]
    assert B0.rank()==n
    assert R0.rank()==B0.ncols()-n
    ordinary=[int(x) for x in B0.elementary_divisors() if x]
    # Direct ordinary source closure is a separate integer calculation.
    ds=[ZZ(1)]
    from sage.all import gcd
    for d in range(1,n):ds.append(gcd([comb(n,j)*ds[d-j] for j in range(1,d+1)]))
    want=matrix.diagonal(ZZ,ds).elementary_divisors()
    assert ordinary==list(want)
    return dict(result='passed',n=n,order=order,closure_rounds=rounds,image_generators=B.ncols(),
                relations=rel.ncols(),specialization_kernel_factors=torsion,
                ordinary_smith=ordinary,full_module_closure=True,
                presentation=[[str(x) for x in row] for row in B.rows()],
                syzygies=[[str(x) for x in row] for row in rel.rows()],
                source_sha256={'experiments/type_a_k_specialization.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})

if __name__=='__main__':
    import sys
    ns=list(map(int,sys.argv[1:])) or list(range(2,7))
    for n in ns:
        row=run(n)
        alternate=run(n,'Dp')
        assert row['specialization_kernel_factors']==alternate['specialization_kernel_factors']
        assert row['ordinary_smith']==alternate['ordinary_smith']
        row['alternate_order_check']={'order':'Dp','image_generators':alternate['image_generators'],'relations':alternate['relations'],'same_kernel_and_ordinary_factors':True}
        (ROOT/'experiments'/('type-a-k-specialization-%d.json'%n)).write_text(json.dumps(row,indent=2)+'\n')
        print(json.dumps({k:v for k,v in row.items() if k not in ('presentation','syzygies')}),flush=True)
