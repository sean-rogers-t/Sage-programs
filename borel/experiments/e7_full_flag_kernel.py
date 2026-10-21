"""Exact 2-saturated full-flag characteristic kernels in bounded E7 degrees.

The rational Borel kernel is reduced by the common monic quadratic.
Integral kernel lattices are certified by mod-2 saturation, before images.
No rational generator subalgebra is presumed to be the integral kernel.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from math import comb
import argparse,gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,CartanMatrix,matrix,vector,gcd,prod,identity_matrix
ROOT=Path(__file__).resolve().parents[1]

def compositions(n,k=6):
 if n<0:return
 if k==1:yield (n,);return
 for i in range(n+1):
  for rest in compositions(n-i,k-1):yield (i,)+rest

def module_exponents(d):
 return [e+(parity,) for parity in (0,1) for e in compositions(d-parity)]

def encoded(f):return [[list(e),str(c)] for e,c in sorted(f.dict().items())]
def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=12);args=ap.parse_args()
 assert 10<=args.max_degree<=13
 start=perf_counter();C=matrix(QQ,CartanMatrix('E7'));XB=C.matrix_from_columns(range(6)).augment(vector(QQ,[0]*6+[1]).column())
 R=PolynomialRing(QQ,7,names=['a'+str(i) for i in range(1,7)]+['w']);a=R.gens()[:6];w=R.gen(6)
 B=R(-QQ(3)/2*(vector(R,a)*C[:6,:6].inverse()*vector(R,a)));quad=w*w-B
 Bpowers=[R.one()]
 for i in range(1,7):Bpowers.append(Bpowers[-1]*B)
 def reduced(f):
  out=R.zero()
  for exp,c in R(f).dict().items():
   e=tuple(exp);out+=c*R({e[:6]+(e[6]%2,):1})*Bpowers[e[6]//2]
  return R(out)
 orbit={(0,0,0,0,0,0,1)};pending=list(orbit)
 while pending:
  t=pending.pop()
  for i in range(7):
   u=tuple(t[j]-t[i]*C[i,j] for j in range(7))
   if u not in orbit:orbit.add(u);pending.append(u)
 assert len(orbit)==56
 inverse=XB.inverse();forms=[sum(c*x for c,x in zip(inverse*vector(QQ,t),R.gens())) for t in sorted(orbit)]
 assert reduced(sum(f*f for f in forms))==0
 relations={}
 for degree in (6,8,10,12):
  if degree>args.max_degree:continue
  f=reduced(sum(l**degree for l in forms));den=f.denominator();g=gcd((den*f).coefficients());f=R(den*f/g)
  assert f.denominator()==1;relations[degree]=f
  print('fresh moment',degree,len(f.dict()),round(perf_counter()-start,2),flush=True)
 prior=json.loads((ROOT/'experiments/e7-refinement-relations.json').read_text())
 for d,key in ((6,'sextic'),(8,'octic')):
  assert relations[d]==R({tuple(e):QQ(c) for e,c in prior[key]})
 q=PolynomialRing(ZZ,'q').gen();P=prod(sum(q**i for i in range(d)) for d in (2,6,8,10,12,14,18));half=prod(sum(q**i for i in range(d)) for d in (6,8,10,12,14,18))
 certificate=dict(quadratic=encoded(quad),relations={str(d):encoded(f) for d,f in relations.items()},degrees=[])
 report=dict(result='running',coefficient_ring='Z_(2)',retained_nodes=[],degrees=[])
 for degree in range(10,args.max_degree+1):
  support=module_exponents(degree);index={e:i for i,e in enumerate(support)}
  labels=[(d,e) for d in relations for e in module_exponents(degree-d)]
  expected=len(support)-int(P[degree]);assert len(labels)==expected,(degree,len(labels),expected)
  entries={}
  for col,(d,e) in enumerate(labels):
   f=reduced(relations[d]*R({e:1}))
   for u,c in f.dict().items():entries[index[tuple(u)],col]=ZZ(c)
  A=matrix(ZZ,len(support),len(labels),entries,sparse=False);del entries
  stages=[]
  for iteration in range(24):
   mod=matrix(GF(2),A);K=mod.right_kernel().basis_matrix();nullity=K.nrows()
   print('kernel',degree,'iteration',iteration,'defect',nullity,round(perf_counter()-start,2),flush=True)
   if not nullity:break
   K=K.echelon_form();pivots=list(K.pivots());assert K.matrix_from_columns(pivots)==identity_matrix(GF(2),nullity)
   changes=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in changes.list())
   changes=matrix(ZZ,changes.nrows(),changes.ncols(),[c//2 for c in changes.list()])
   stages.append(dict(replace_columns=list(map(int,pivots)),combinations=[[j for j,c in enumerate(row) if c] for row in K]))
   for i,j in enumerate(pivots):A.set_column(j,changes.column(i))
  else:raise RuntimeError('saturation iteration bound')
  # There are no ordinary relation syzygies below degree 14. Full mod-2
  # rank therefore also certifies the required rational kernel dimension.
  assert mod.ncols()-K.nrows()==expected
  oddrows=[i for i,e in enumerate(support) if e[6]==1]
  oddrank=mod.matrix_from_rows(oddrows).rank();m=comb(degree+4,5)-int(oddrank)
  prediction=int(half[degree-1]);row=dict(degree=degree,coefficient_rows=len(support),kernel_rank=expected,saturation_rounds=len(stages),saturation_index_valuation=sum(len(s['replace_columns']) for s in stages),odd_relation_rank=int(oddrank),multiplicity_2=m,half_poincare_prediction=prediction,prediction_agrees=m==prediction)
  report['degrees'].append(row);certificate['degrees'].append(dict(degree=degree,initial_columns=[[d,list(e)] for d,e in labels],saturation=stages))
  print('EXACT',row,round(perf_counter()-start,2),flush=True)
 archive=ROOT/'experiments/e7-full-flag-kernel-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(certificate,separators=(',',':')).encode(),mtime=0))
 report.update(result='passed',seconds=round(perf_counter()-start,3),certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),scope='Exact full-flag relative image multiplicities only in the recorded degrees. Agreement with q P_W/(1+q) is evidence, not an all-degree theorem.',source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in (Path(__file__).resolve(),ROOT/'experiments/e7-refinement-relations.json')})
 (ROOT/'experiments/e7-full-flag-kernel.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED',report['seconds'],flush=True)
if __name__=='__main__':main()
