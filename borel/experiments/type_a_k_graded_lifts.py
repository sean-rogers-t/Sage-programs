"""Homogeneous full-source membership certificates, not coefficient slices.

The forward integer matrices certify every auxiliary row is in the source.
No completeness assertion about a computed standard basis is needed to
verify the resulting identities. Run using sage -python, e.g. with 3 4 7.
"""
from pathlib import Path
from math import comb
from hashlib import sha256
import gzip,json
from sage.all import ZZ,QQ,PolynomialRing,matrix

ROOT=Path(__file__).resolve().parents[1]
def exponents(weights,d):
    if not weights:
        return [()] if d==0 else []
    return [(k,)+tail for k in range(d//weights[0]+1)
            for tail in exponents(weights[1:],d-k*weights[0])]
def sparse_rows(M):
    return [[[j,int(x)] for j,x in enumerate(row) if x] for row in M.rows()]

def run(m):
    fac=list(ZZ(m).factor());assert len(fac)==1
    p,s=fac[0]
    R=PolynomialRing(ZZ,names=['L']+['c%d'%i for i in range(1,m)])
    L=R.gen(0);cs=[R(1)]+list(R.gens()[1:])
    cs.append((-1)**(m+1)*sum((-1)**k*cs[k]*L**(m-k) for k in range(m)))
    fs=[R(1)]+[sum((-1)**k*comb(m-k,j-k)*cs[k]*L**(j-k) for k in range(j+1)) for j in range(1,m)]
    gens=cs[1:]+fs[1:];degrees=list(range(1,m+1))+list(range(1,m))
    weights=[1]+list(range(1,m));bases={0:[(0,)*m]};Js={0:matrix(ZZ,[[1]])};records=[]
    for d in range(1,2*m-1):
        basis=exponents(weights,d);bases[d]=basis;index={e:i for i,e in enumerate(basis)}
        blocks=[];cblocks=[]
        for g,w,gi in zip(gens,degrees,range(len(gens))):
            if w>d:continue
            B=bases[d-w]
            op=matrix(ZZ,len(B),len(basis),sparse=True)
            for a,e in enumerate(B):
                for z,coefficient in g.dict().items():
                    op[a,index[tuple(x+y for x,y in zip(e,z))]]+=coefficient
            block=Js[d-w]*op
            blocks.append(block)
            if gi<m:cblocks.append(block)
        M=matrix(ZZ,[row for block in blocks for row in block.rows()])
        H,U=M.hermite_form(transformation=True,include_zero_rows=False)
        assert U*M==H and H.nrows()==len(basis)
        Js[d]=H
        C=matrix(ZZ,[row for block in cblocks for row in block.rows()])
        IH,IU=C.hermite_form(transformation=True,include_zero_rows=False)
        target_rows=[];pairs=[]
        for i in range(1,m):
            j=d-i
            if not i<=j<m:continue
            if d<m:
                ratio=QQ(comb(m,i)*comb(m,j))/comb(m,d)
                denominator=ratio.denominator();numerator=ratio.numerator()
                assert denominator%p
                poly=m*(denominator*fs[i]*fs[j]-numerator*fs[d])
            else:
                denominator=1;numerator=0;poly=m*fs[i]*fs[j]
            row=[ZZ(0)]*len(basis)
            for e,x in poly.dict().items():row[index[tuple(e)]]=x
            target_rows.append(row);pairs.append(dict(i=i,j=j,denominator=int(denominator),numerator=int(numerator)))
        T=matrix(ZZ,target_rows) if target_rows else matrix(ZZ,0,len(basis))
        if target_rows:
            X=IH.transpose().solve_right(T.transpose()).transpose()
            assert all(x.denominator()==1 for x in X.list()),(m,d,'nonintegral lift')
            W=matrix(ZZ,X)*IU
            assert W*C==T
        else:W=matrix(ZZ,0,C.nrows())
        records.append(dict(degree=d,basis=[list(e) for e in basis],J=sparse_rows(H),
                            construction=sparse_rows(U),candidate_count=M.nrows(),
                            coefficient_candidate_count=C.nrows(),pairs=pairs,lifts=sparse_rows(W)))
        print(json.dumps(dict(m=m,degree=d,ambient_rank=len(basis),targets=len(pairs),result='passed')),flush=True)
    data=dict(result='passed',m=m,prime=int(p),records=records,
              scope='Full universal homogeneous source; all quadratic relations certified in this rank. Translation and nonzerodivisor descent are justified in the proof note.',
              source_sha256={'experiments/type_a_k_graded_lifts.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    path=ROOT/'experiments'/('type-a-k-graded-lifts-%d.json.gz'%m)
    with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as f:f.write(json.dumps(data,separators=(',',':')).encode())
    return data
if __name__=='__main__':
    import sys
    for m in map(int,sys.argv[1:] or [3,4,7]):run(m)
