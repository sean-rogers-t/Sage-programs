"""Full bounded root invariant components for E7/(D5 A1), localized at two.

The full A4 x A1 coefficient ring on Q is polynomial in h,b,c2..c5,r,
with sum(u)=b, r=(w-z)(h-w-z) and t=r+hb/2-b^2/4. Saturate the full
rational D5 x A1 invariant component, retaining exact preimage columns.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import argparse,gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,matrix,vector,binomial,gcd
ROOT=Path(__file__).resolve().parents[1]
WEIGHTS=(1,1,2,3,4,5,2)
CLASSICAL=(1,2,2,4,5,6,8)

def powers(ws,d):
 if d<0:return
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def context(kind='M'):
 P=PolynomialRing(QQ,7,names=('h','z','c2','c3','c4','c5','r'),order=TermOrder('wdegrevlex',WEIGHTS));h,z,c2,c3,c4,c5,r=P.gens();es=[P.one(),z if kind=='M' else 2*z,c2,c3,c4,c5];t=r+h*z/2-z*z/4 if kind=='M' else r
 ps=[P(5)]
 for k in range(1,9):ps.append(sum((-1)**(j-1)*es[j]*ps[k-j] for j in range(1,min(k,6)))+((-1)**(k-1)*k*es[k] if k<=5 else 0))
 xp={k:sum(binomial(k,j)*2**j*(-h)**(k-j)*ps[j] for j in range(k+1)) for k in (2,4,6,8)};sq=[P.one()]
 for k in range(1,5):sq.append(sum((-1)**(j-1)*sq[k-j]*xp[2*j] for j in range(1,k+1))/k)
 pf=sum(2**j*(-h)**(5-j)*es[j] for j in range(6));raw=[h,t,sq[1],sq[2],pf,sq[3],sq[4]];factors=[QQ(gcd(f.coefficients())) for f in raw]
 return P,[P(f/k) for f,k in zip(raw,factors)],factors

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=16);args=ap.parse_args();start=perf_counter();P,orig,factors=context();P2=P.change_ring(GF(2));basis={0:[P.one()]};basis2={0:[P2.one()]};coordinates={0:matrix(QQ,1,1,[1])};exponents={0:[(0,)*7]};gs=[];rows=[]
 @lru_cache(None)
 def mon(e):
  if not any(e):return P.one()
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return orig[j]*mon(tuple(u))
 def coefficients(fs,support,ring):
  ix={e:i for i,e in enumerate(support)};return matrix(ring,len(support),len(fs),{(ix[tuple(e)],j):c for j,f in enumerate(fs) for e,c in f.dict().items()})
 for d in range(1,args.max_degree+1):
  es=list(powers(CLASSICAL,d));support=list(powers(WEIGHTS,d));idx={e:i for i,e in enumerate(es)};n=len(es);polys=[orig[0]*b for b in basis[d-1]];cols=[]
  for col in coordinates[d-1].columns():
   v=vector(QQ,n)
   for e,c in zip(exponents[d-1],col):v[idx[(e[0]+1,)+e[1:]]]=c
   cols.append(v)
  for e in es:
   if not e[0]:polys.append(mon(e));cols.append(vector(QQ,n,{idx[e]:1}))
  A=coefficients(polys,support,ZZ);C=matrix(QQ,cols).transpose();assert A.ncols()==n;stages=0
  while True:
   K=matrix(GF(2),A).right_kernel().basis_matrix().echelon_form();piv=list(K.pivots())
   if not piv:break
   B=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in B.list());B=matrix(ZZ,B.nrows(),B.ncols(),[c//2 for c in B.list()]);D=C*matrix(QQ,K.transpose())/2
   for j,k in enumerate(piv):A.set_column(k,B.column(j));C.set_column(k,D.column(j))
   stages+=1
  polys=[P({support[i]:c for i,c in col.dict().items()}) for col in A.columns()];basis[d]=polys;basis2[d]=[P2(f) for f in polys];coordinates[d]=C;exponents[d]=es
  products=[P2(g)*b for g,k in gs for b in basis2[d-k]];old=coefficients(products,support,GF(2));piv=old.augment(matrix(GF(2),A)).pivots();assert len(piv)==n;new=[j-old.ncols() for j in piv if j>=old.ncols()]
  for j in new:gs.append((polys[j],d))
  rows.append(dict(degree=d,rank=n,coefficient_rows=len(support),saturation_rounds=stages,classical_exponents=[list(e) for e in es],rational_coordinates=[[str(c) for c in row] for row in C],new_generator_indices=new,new_generator_polynomials=[str(polys[j]) for j in new]))
  print('M source',d,'rank',n,'new',len(new),'generator degrees',[k for f,k in gs],'seconds',round(perf_counter()-start,2),flush=True)
 cert=dict(original_generators=list(map(str,orig)),normalization_factors=list(map(str,factors)),degrees=rows);archive=ROOT/'experiments/e7-d5a1-sources-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0))
 report=dict(result='passed',max_degree=args.max_degree,coefficient_weights=list(WEIGHTS),classical_degrees=list(CLASSICAL),generator_degrees=[d for f,d in gs],generators=[[str(f),d] for f,d in gs],degrees=[{k:v for k,v in row.items() if k not in ('classical_exponents','rational_coordinates','new_generator_polynomials')} for row in rows],seconds=round(perf_counter()-start,3),scope='Full root invariant components through the stated bound only. No unbounded algebra completeness or image claim.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()});(ROOT/'experiments/e7-d5a1-sources.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
