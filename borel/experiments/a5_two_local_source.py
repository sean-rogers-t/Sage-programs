"""Full bounded A5 root invariant components at two.

The S5 subgroup permutes d_i=x_i-x6, an actual root-lattice basis.
Its exact coefficient ring is R[c1,..,c5]. Rational S6 invariants are
the centered six-variable e2,..,e6. Weight-lattice invariants are the
ordinary six-variable elementary classes with e1=0.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from functools import lru_cache
import argparse,gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,matrix,gcd,binomial,prod
ROOT=Path(__file__).resolve().parents[1]
CW=(1,2,3,4,5);DW=(2,3,4,5,6)

def powers(ws,d):
 if d<0:return
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=24);args=ap.parse_args();start=perf_counter();P=PolynomialRing(QQ,5,names=('c1','c2','c3','c4','c5'),order=TermOrder('wdegrevlex',CW));cs=[P.one()]+list(P.gens())+[P.zero()];mu=cs[1]/6;raw=[sum(binomial(6-j,k-j)*(-mu)**(k-j)*cs[j] for j in range(k+1)) for k in range(2,7)];factors=[QQ(gcd(f.coefficients())) for f in raw];orig=[f/k for f,k in zip(raw,factors)];P2=P.change_ring(GF(2));basis2={0:[P2.one()]};rows=[];gs=[]
 @lru_cache(None)
 def mon(e):return P(prod(g**k for g,k in zip(orig,e)))
 def coeff(fs,support,ring):
  ix={e:i for i,e in enumerate(support)};return matrix(ring,len(support),len(fs),{(ix[tuple(e)],j):c for j,f in enumerate(fs) for e,c in f.dict().items()})
 for d in range(1,args.max_degree+1):
  es=list(powers(DW,d));support=list(powers(CW,d));n=len(es)
  if not n:basis2[d]=[];continue
  A=coeff([mon(e) for e in es],support,ZZ);C=matrix.identity(QQ,n);stages=0
  while True:
   K=matrix(GF(2),A).right_kernel().basis_matrix().echelon_form();piv=list(K.pivots())
   if not piv:break
   B=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in B.list());B=matrix(ZZ,B.nrows(),B.ncols(),[c//2 for c in B.list()]);D=C*matrix(QQ,K.transpose())/2
   for j,k in enumerate(piv):A.set_column(k,B.column(j));C.set_column(k,D.column(j))
   stages+=1
  polys=[P({support[i]:c for i,c in col.dict().items()}) for col in A.columns()];basis2[d]=[P2(f) for f in polys];products=[P2(f)*b for f,k in gs for b in basis2[d-k]];old=coeff(products,support,GF(2));piv=old.augment(matrix(GF(2),A)).pivots();assert len(piv)==n;new=[j-old.ncols() for j in piv if j>=old.ncols()]
  for j in new:gs.append((polys[j],d))
  rows.append(dict(degree=d,rank=n,classical_exponents=[list(e) for e in es],rational_coordinates=[[str(c) for c in row] for row in C],new_generator_indices=new));print('A5 root source',d,n,'generators',[k for f,k in gs],round(perf_counter()-start,2),flush=True)
 cert=dict(original_generators=list(map(str,orig)),normalization_factors=list(map(str,factors)),degrees=rows);archive=ROOT/'experiments/a5-two-local-source-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));report=dict(result='passed',max_degree=args.max_degree,generator_degrees=[k for f,k in gs],generators=[[str(f),k] for f,k in gs],scope='Full bounded source components only. An all-degree algebra certificate is a separate requirement.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),seconds=round(perf_counter()-start,3),source_sha256={Path(__file__).resolve().relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()});(ROOT/'experiments/a5-two-local-source.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
