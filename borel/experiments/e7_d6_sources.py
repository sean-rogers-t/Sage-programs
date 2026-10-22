"""Exact full graded D6 invariant sources for E7/P1 at the prime two.

S6 permutation coordinates give a polynomial coefficient ring for each
actual lattice. Every recorded component is saturated independently;
no all-degree generator-completeness assertion is made by this script.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from functools import lru_cache
import argparse,gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,matrix,vector,prod,binomial,gcd
ROOT=Path(__file__).resolve().parents[1]
WEIGHTS=(1,1,2,3,4,5,6)
CLASSICAL=(1,2,4,6,6,8,10)

def powers(ws,d):
 if d<0:return
 if not ws:
  if not d:yield ()
  return
 for j in range(d//ws[0]+1):
  for e in powers(ws[1:],d-j*ws[0]):yield (j,)+e

def context(kind):
 P=PolynomialRing(QQ,7,names=('z','c1','c2','c3','c4','c5','c6'),order=TermOrder('wdegrevlex',WEIGHTS));z,*cs=P.gens();es=[P.one()]+cs
 h=cs[0]-2*z if kind=='N' else cs[0]+2*z
 shift=-h if kind=='N' else -cs[0]-z
 ps=[P(6)]
 for d in range(1,11):ps.append(P(sum((-1)**(j-1)*es[j]*ps[d-j] for j in range(1,min(d,7)))+((-1)**(d-1)*d*es[d] if d<=6 else 0)))
 xp={d:P(sum(binomial(d,j)*2**j*shift**(d-j)*ps[j] for j in range(d+1))) for d in (2,4,6,8,10)}
 sq=[P.one()]
 for d in range(1,6):sq.append(P(sum((-1)**(j-1)*sq[d-j]*xp[2*j] for j in range(1,d+1))/d))
 eu=P(sum(2**j*shift**(6-j)*es[j] for j in range(7)))
 raw=[h,sq[1],sq[2],sq[3],eu,sq[4],sq[5]];divs=[ZZ(gcd(f.coefficients())) for f in raw]
 return P,[P(f/k) for f,k in zip(raw,divs)],divs

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=12);args=ap.parse_args();start=perf_counter();certificate={};reports={}
 for kind in ('N','M'):
  P,orig,divs=context(kind);basis={0:[P.one()]};gs=[];rows=[];records=[]
  @lru_cache(None)
  def classical(e):
   if not any(e):return P.one()
   j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
   return orig[j]*classical(tuple(u))
  for d in range(1,args.max_degree+1):
   exps=list(powers(CLASSICAL,d));support=list(powers(WEIGHTS,d));fs=[classical(e) for e in exps];dicts=[{tuple(e):c for e,c in f.dict().items()} for f in fs]
   F=matrix(QQ,[[f.get(e,0) for f in dicts] for e in support]);assert F.denominator()==1
   # Track the basis change during saturation instead of reconstructing it
   # by a large rational inverse. The starting columns are a rational basis;
   # RREF kernel pivots make each replacement an exact 2-local enlargement.
   A=matrix(ZZ,F);C=matrix(QQ,len(exps),len(exps),1);stages=0
   while True:
    K=matrix(GF(2),A).right_kernel().basis_matrix().echelon_form();piv=list(K.pivots())
    if not piv:break
    changes=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in changes.list());changes=matrix(ZZ,changes.nrows(),changes.ncols(),[c//2 for c in changes.list()])
    coordinates=C*matrix(QQ,K.transpose())/2
    for j,k in enumerate(piv):A.set_column(k,changes.column(j));C.set_column(k,coordinates.column(j))
    stages+=1
   assert A.ncols()==len(exps) and matrix(GF(2),A).rank()==len(exps)
   assert F*C==A
   polys=[P({support[i]:c for i,c in col.dict().items()}) for col in A.columns()];basis[d]=polys
   # Find minimal algebra generators modulo two, using complete lower components.
   products=[P(f*b) for f,k in gs for b in basis.get(d-k,[])];old=matrix(GF(2),len(support),0)
   if products:
    dictionaries=[{tuple(e):c for e,c in f.dict().items()} for f in products];old=matrix(GF(2),[[f.get(e,0) for f in dictionaries] for e in support])
   pivots=old.augment(matrix(GF(2),A)).pivots();new=[j-old.ncols() for j in pivots if j>=old.ncols()]
   for j in new:gs.append((polys[j],d))
   assert len(pivots)==len(exps)
   record=dict(degree=d,rank=len(exps),coefficient_rows=len(support),saturation_rounds=stages,new_generator_indices=new)
   records.append(dict(**record,classical_exponents=[list(e) for e in exps],rational_coordinates=[[str(c) for c in row] for row in C],basis_polynomials=list(map(str,polys))))
   rows.append(dict(**record,new_generators=[str(polys[j]) for j in new]));print(kind,'source',d,'rank',len(exps),'rows',len(support),'new',len(new),'seconds',round(perf_counter()-start,2),flush=True)
  certificate[kind]=dict(original_generators=list(map(str,orig)),original_divisors=list(map(int,divs)),degrees=records)
  reports[kind]=dict(generator_degrees=[d for f,d in gs],generators=[[str(f),d] for f,d in gs],degrees=rows)
 path=ROOT/'experiments/e7-d6-sources-certificates.json.gz';path.write_bytes(gzip.compress(json.dumps(certificate,separators=(',',':')).encode(),mtime=0))
 report=dict(result='passed',scope='Full 2-local source components only through the recorded bound; no all-degree completeness claim.',max_degree=args.max_degree,coefficient_weights=list(WEIGHTS),classical_degrees=list(CLASSICAL),sources=reports,seconds=round(perf_counter()-start,3),certificate_archive=path.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(path.read_bytes()).hexdigest(),source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
 (ROOT/'experiments/e7-d6-sources.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
