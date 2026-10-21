"""Construct the full 2-local D4 invariant source by coefficient saturation."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json,argparse
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,prod,binomial,gcd,identity_matrix,TermOrder
ROOT=Path(__file__).resolve().parents[1]
WEIGHTS=(1,1,1,1,2,3)
def powers(ws,d):
 if d<0:return
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def context():
 P=PolynomialRing(QQ,6,names=('h','y','z','v','c2','c3'),order=TermOrder('wdegrevlex',WEIGHTS));h,y,z,v,c2,c3=P.gens();es=[P.one(),2*z-v-y,c2,c3]
 ps=[P(3)]
 for d in range(1,7):ps.append(P(sum((-1)**(j-1)*es[j]*ps[d-j] for j in range(1,min(d,4)))+((-1)**(d-1)*d*es[d] if d<=3 else 0)))
 xp=[P(4)]+[P(sum(binomial(d,j)*2**j*(-h)**(d-j)*ps[j] for j in range(d+1))+(2*v-h)**d) for d in range(1,7)]
 sq=[P.one()]
 for d in range(1,4):sq.append(P(sum((-1)**(j-1)*sq[d-j]*xp[2*j] for j in range(1,d+1))/d))
 eu=P((2*v-h)*sum(2**j*(-h)**(3-j)*es[j] for j in range(4)))
 gs=[(h,1),(y,1),(sq[1],2),(sq[2],4),(eu,4),(sq[3],6)]
 gs=[(P(f/gcd(f.coefficients())),d) for f,d in gs]
 return P,gs

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=12);args=ap.parse_args();start=perf_counter();P,original=context();gs=[];basis={0:[P.one()]};record=[]
 for d in range(1,args.max_degree+1):
  gs.extend((f,k) for f,k in original if k==d)
  fs=[P(g*b) for g,k in gs if k<=d for b in basis[d-k]]
  support=list(powers(WEIGHTS,d));ds=[{tuple(e):c for e,c in f.dict().items()} for f in fs];A=matrix(QQ,[[f.get(e,0) for f in ds] for e in support]);den=A.denominator();A=matrix(ZZ,den*A).transpose().hermite_form(include_zero_rows=False).transpose();assert den==1
  before=matrix(ZZ,A);stages=0
  while True:
   K=matrix(GF(2),A).right_kernel().basis_matrix().echelon_form();piv=list(K.pivots())
   if not piv:break
   changes=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in changes.list());changes=matrix(ZZ,changes.nrows(),changes.ncols(),[c//2 for c in changes.list()])
   for j,k in enumerate(piv):A.set_column(k,changes.column(j))
   stages+=1
  # Coefficient columns are independent modulo two: whole rational invariant component.
  es=list(powers((1,1,2,4,4,6),d));assert A.ncols()==len(es),(d,A,es,fs,support,ds)
  rows=list(matrix(GF(2),A).transpose().pivots());T=A.matrix_from_rows(rows).change_ring(QQ).inverse()*before.matrix_from_rows(rows);assert T.denominator()%2
  old=matrix(GF(2),T);rank=old.rank();new=[]
  for j in range(A.ncols()):
   unit=vector(GF(2),A.ncols(),{j:1});more=old.augment(unit.column())
   if more.rank()>rank:
    f=P(sum(c*P({e:1}) for e,c in zip(support,A.column(j))));gs.append((f,d));new.append(str(f));old=more;rank+=1
  assert rank==A.ncols()
  basis[d]=[P(sum(c*P({e:1}) for e,c in zip(support,col))) for col in A.columns()]
  row=dict(degree=d,rank=A.ncols(),coefficient_rows=len(support),saturation_rounds=stages,new_generators=new);record.append(row)
  print('SOURCE',d,A.dimensions(),'rounds',stages,'new',len(new),'total',len(gs),round(perf_counter()-start,2),flush=True)
  for f in new:print(' NEW terms',len(P(f).dict()),flush=True)
 report=dict(result='bounded_source_components_passed',degrees=record,generators=[[str(f),d] for f,d in gs],source_sha256=sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest())
 (ROOT/'experiments/e7-d4-source.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
