"""One-parameter image defects: quotients of the full defect, not upper bounds."""
from pathlib import Path
from math import comb
import json
from hashlib import sha256
from sage.all import ZZ, PolynomialRing, matrix, vector, singular

ROOT=Path(__file__).resolve().parents[1]

def run(n,jvary=1,additive_trace=False,order='dp'):
    S=PolynomialRing(ZZ,['a']);a=S.gen()
    cs=[S(comb(n,j))+(a*comb(n-1,j-1) if additive_trace and 1<=j<n else (a if not additive_trace and j==jvary else 0)) for j in range(n+1)]
    P=PolynomialRing(S,'L');L=P.gen()
    f=sum((-1)**j*cs[j]*L**(n-j) for j in range(n+1))
    fs=[sum((-1)**k*comb(n-k,j-k)*cs[k]*L**(j-k) for k in range(j+1)) for j in range(1,n)]
    cols=[vector(S,[(g%f)[i] for i in range(n)]) for g in [P(1)]+fs]
    mult=[matrix(S,n,n,lambda i,k:((g*L**k)%f)[i]) for g in fs]
    singular.eval('ring r=integer,(a),'+order+';')
    literal=lambda M:','.join(str(x).replace('**','^') for x in M.list())
    initial=matrix(S,cols).transpose()
    singular.eval('matrix M[%d][%d]=%s;module G=std(module(M));'%(n,n,literal(initial)))
    for j,A in enumerate(mult):singular.eval('matrix T%d[%d][%d]=%s;'%(j,n,n,literal(A)))
    rounds=0
    while True:
        rounds+=1;singular.eval('module extra;')
        for j in range(n-1):singular.eval('extra=extra,reduce(module(T%d*matrix(G)),G);'%j)
        if int(singular.eval('size(extra);'))==0:break
        singular.eval('G=std(G,extra);kill extra;')
        if rounds>=30:raise RuntimeError('No closure certificate within limit')
    B=singular('matrix(G)').sage()
    singular.eval('module sy=syz(G);matrix syzM=matrix(sy);')
    rel=singular('syzM').sage()
    assert B*rel==0
    R0=matrix(ZZ,rel.nrows(),rel.ncols(),[ZZ(x.subs(a=0)) for x in rel.list()])
    factors=[int(abs(x)) for x in R0.elementary_divisors() if abs(x)>1]
    B0=matrix(ZZ,B.nrows(),B.ncols(),[ZZ(x.subs(a=0)) for x in B.list()])
    assert B0.rank()==n and R0.rank()==B0.ncols()-n
    return dict(result='passed',n=n,varied_coefficient=jvary if not additive_trace else 'additive_trace',order=order,rounds=rounds,generators=B.ncols(),
                kernel_factors=factors,annihilated_by_n=all(n%x==0 for x in factors),
                presentation=[[str(x) for x in row] for row in B.rows()],
                syzygies=[[str(x) for x in row] for row in rel.rows()],
                source_sha256={'experiments/type_a_k_slices.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})

if __name__=='__main__':
    import sys
    trace='--trace' in sys.argv
    ns=list(map(int,[x for x in sys.argv[1:] if x!='--trace'])) or ([3,5,7] if trace else [7,8])
    for n in ns:
        row=run(n,additive_trace=trace)
        check=run(n,additive_trace=trace,order='Dp')
        assert row['kernel_factors']==check['kernel_factors']
        row['alternate_order_verified']=True
        if trace:
            fac=list(ZZ(n).factor());assert len(fac)==1 and n>=3
            p,e=fac[0]
            want=sorted(int(p**(e-ZZ(j).valuation(p))) for j in range(1,n))
            assert row['kernel_factors']==want
        (ROOT/'experiments'/('type-a-k-%s-%d.json'%('trace-slice' if trace else 'slice',n))).write_text(json.dumps(row,indent=2)+'\n')
        print(json.dumps({k:v for k,v in row.items() if k not in ('presentation','syzygies')}),flush=True)
