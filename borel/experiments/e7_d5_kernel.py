"""Exact 2-local D5 invariant Borel kernels; integral coefficient saturation."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json,gzip,argparse
from sage.all import QQ,ZZ,GF,CartanMatrix,PolynomialRing,matrix,vector,prod,gcd,identity_matrix
from e7_full_flag_kernel import encoded,compositions,digest

ROOT=Path(__file__).resolve().parents[1]
DEGREES=(1,2,4,5,6,8)
def powers(ws,n):
 if n<0:return
 if not ws:
  if not n:yield ()
  return
 for k in range(n//ws[0]+1):
  for e in powers(ws[1:],n-k*ws[0]):yield (k,)+e

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=13);args=ap.parse_args()
 assert 10<=args.max_degree<=13
 start=perf_counter();C=matrix(QQ,CartanMatrix('E7'))
 R=PolynomialRing(QQ,7,names=['a'+str(i) for i in range(1,7)]+['w']);a=R.gens()[:6];w=R.gen(6)
 B=R(-QQ(3)/2*(vector(R,a)*C[:6,:6].inverse()*vector(R,a)))
 bp=[R.one()]
 for k in range(args.max_degree//2):bp.append(bp[-1]*B)
 def reduced(f):
  out=R.zero()
  for e,c in R(f).dict().items():out+=c*R({tuple(e[:6])+(e[6]%2,):1})*bp[e[6]//2]
  return R(out)
 priorpath=ROOT/'experiments/e7-full-flag-kernel.json';prior=json.loads(priorpath.read_text());archive=ROOT/prior['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==prior['certificate_archive_sha256']
 old=json.loads(gzip.decompress(archive.read_bytes()));relations={int(d):R({tuple(e):QQ(c) for e,c in f}) for d,f in old['relations'].items()}
 om=C[:6,:6].inverse()*vector(R,a);h=om[0]
 xs=[-h-2*om[1]+2*om[2],-h+2*om[1]+2*om[2]-2*om[3],-h+2*om[3]-2*om[4],-h+2*om[4]-2*om[5],-h+2*om[5]]
 from itertools import combinations
 elem=lambda j:sum(prod(xs[i]**2 for i in inds) for inds in combinations(range(5),j))
 gs=[h,elem(1),elem(2),prod(xs),elem(3),elem(4)]
 # Clear individual rational denominators: rational spanning suffices before saturation.
 gs=[R(g.denominator()*g) for g in gs]
 for node in (2,3,4,5,6):
  new=[R.gen(j)-C[j,node-1]*R.gen(node-1) if j<6 else w for j in range(7)]
  assert all(g(*new)==g for g in gs)
 inv={d:[prod(g**k for g,k in zip(gs,e)) for e in powers(DEGREES,d)] for d in range(args.max_degree+1)}
 q=PolynomialRing(ZZ,'q').gen();WG=prod(sum(q**i for i in range(d)) for d in (2,6,8,10,12,14,18));WI=prod(sum(q**i for i in range(d)) for d in (2,4,5,6,8));P=WG//WI;pred=(q*P)//(1+q)
 report=dict(result='running',retained_nodes=[2,3,4,5,6],degrees=[])
 cert=dict(generators=[encoded(g) for g in gs],quadratic=encoded(w*w-B),relations={str(d):encoded(f) for d,f in relations.items()},degrees=[])
 for degree in range(10,args.max_degree+1):
  fs=[];labels=[]
  for rd,f in relations.items():
   for par in (0,1):
    for ex,g in zip(powers(DEGREES,degree-rd-par),inv.get(degree-rd-par,[])):
     z=reduced(f*g*w**par);den=z.denominator();z=R(den*z/gcd((den*z).coefficients()));fs.append(z);labels.append([rd,par,list(ex)])
  expected=len(inv[degree])+len(inv[degree-1])-int(P[degree]);assert len(fs)==expected
  support=sorted(set().union(*(set(f.dict()) for f in fs)));index={tuple(e):i for i,e in enumerate(support)}
  entries={(index[tuple(e)],j):ZZ(c) for j,f in enumerate(fs) for e,c in f.dict().items()}
  A=matrix(ZZ,len(support),len(fs),entries,sparse=False);stages=[]
  print('built',degree,A.dimensions(),round(perf_counter()-start,2),flush=True)
  for iteration in range(40):
   mod=matrix(GF(2),A);K=mod.right_kernel().basis_matrix().echelon_form();piv=list(K.pivots())
   if not piv:break
   assert K.matrix_from_columns(piv)==identity_matrix(GF(2),len(piv))
   changes=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in changes.list());changes=matrix(ZZ,changes.nrows(),changes.ncols(),[c//2 for c in changes.list()])
   stages.append(dict(replace_columns=list(map(int,piv)),combinations=[[j for j,c in enumerate(row) if c] for row in K]))
   for j,k in enumerate(piv):A.set_column(k,changes.column(j))
  else:raise RuntimeError('bound')
  assert mod.rank()==expected
  odd=[i for i,e in enumerate(support) if e[6]==1];r=int(mod.matrix_from_rows(odd).rank());m=len(inv[degree-1])-r
  row=dict(degree=degree,invariant_rank=len(inv[degree-1]),kernel_rank=expected,odd_relation_rank=r,multiplicity_2=m,half_poincare_prediction=int(pred[degree]),prediction_agrees=m==pred[degree],saturation_rounds=len(stages),saturation_index_valuation=sum(len(z['replace_columns']) for z in stages));report['degrees'].append(row)
  cert['degrees'].append(dict(degree=degree,initial_columns=labels,initial_polynomials=[encoded(f) for f in fs],saturation=stages))
  print('EXACT',row,round(perf_counter()-start,2),flush=True)
 archive=ROOT/'experiments/e7-d5-kernel-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0))
 report.update(result='passed',seconds=round(perf_counter()-start,3),certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in (Path(__file__).resolve(),priorpath,ROOT/'experiments/e7_full_flag_kernel.py')})
 (ROOT/'experiments/e7-d5-kernel.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
