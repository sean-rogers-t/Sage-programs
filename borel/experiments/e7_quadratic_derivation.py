"""Test the proposed central derivation on the *saturated* E7 kernel.

Fresh Weyl-orbit moments, no Smith tables, and no dimension prediction for
the relative image. The classical rational Borel dimensions certify the
integral kernel after coefficient saturation. This is bounded, not an
all-degree structural proof.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,CartanMatrix,matrix,vector,gcd,prod
from e7_full_flag_kernel import compositions,module_exponents,encoded
ROOT=Path(__file__).resolve().parents[1]

def main():
 start=perf_counter();C=matrix(QQ,CartanMatrix('E7'))
 XB=C.matrix_from_columns(range(6)).augment(vector(QQ,[0]*6+[1]).column())
 R=PolynomialRing(QQ,7,names=['a'+str(i) for i in range(1,7)]+['w'])
 a=R.gens()[:6];w=R.gen(6);B=R(-QQ(3)/2*(vector(R,a)*C[:6,:6].inverse()*vector(R,a)))
 bp=[R.one()]
 for j in range(6):bp.append(bp[-1]*B)
 def reduced(f):
  return R(sum(c*R({tuple(e[:6])+(e[6]%2,):1})*bp[e[6]//2] for e,c in f.dict().items()))
 orbit={(0,0,0,0,0,0,1)};pending=list(orbit)
 while pending:
  t=pending.pop()
  for i in range(7):
   u=tuple(t[j]-t[i]*C[i,j] for j in range(7))
   if u not in orbit:orbit.add(u);pending.append(u)
 assert len(orbit)==56
 inv=XB.inverse();forms=[sum(c*x for c,x in zip(inv*vector(QQ,t),R.gens())) for t in sorted(orbit)]
 assert reduced(sum(f*f for f in forms))==0
 relations={}
 for d in (6,8,10,12):
  f=reduced(sum(l**d for l in forms));den=f.denominator();f=R(den*f/gcd((den*f).coefficients()))
  assert f.denominator()==1;relations[d]=f
  print('fresh moment',d,round(perf_counter()-start,2),flush=True)
 q=PolynomialRing(ZZ,'q').gen();P=prod(sum(q**j for j in range(k)) for k in (2,6,8,10,12,14,18))
 certificate=dict(quadratic=encoded(w*w-B),relations={str(d):encoded(f) for d,f in relations.items()},degrees=[])
 kernels={};supports={};rows=[]
 for d in range(14):
  support=module_exponents(d);ix={e:i for i,e in enumerate(support)};supports[d]=support
  labels=[(k,e) for k in relations for e in module_exponents(d-k)]
  A=matrix(ZZ,len(support),len(labels),sparse=False)
  for j,(k,e) in enumerate(labels):
   for u,c in reduced(relations[k]*R({e:1})).dict().items():A[ix[tuple(u)],j]=ZZ(c)
  assert len(labels)==len(support)-int(P[d]) # no ordinary syzygy before 14
  stages=[]
  while True:
   mod=matrix(GF(2),A);K=mod.right_kernel().basis_matrix().echelon_form()
   if not K.nrows():break
   positions=list(K.pivots());changes=A*matrix(ZZ,K.transpose())
   assert all(c%2==0 for c in changes.list())
   changes=matrix(ZZ,changes.nrows(),changes.ncols(),[c//2 for c in changes.list()])
   stages.append(dict(replace_columns=list(map(int,positions)),combinations=[[j for j,c in enumerate(row) if c] for row in K]))
   for j,pos in enumerate(positions):A.set_column(pos,changes.column(j))
   assert len(stages)<=8
  kernels[d]=mod;assert mod.rank()==len(labels)
  # d/dw on the quadratic quotient sends a+w*b to b. Its square is zero.
  # Check the FULL primitive kernel, including new half-sum relations.
  delta=matrix(GF(2),len(supports[d-1]) if d else 0,len(labels))
  if d:
   previous={e:i for i,e in enumerate(supports[d-1])}
   for i,e in enumerate(support):
    if e[-1]:delta.set_row(previous[e[:-1]+(0,)],mod.row(i))
   assert kernels[d-1].augment(delta).rank()==kernels[d-1].rank()
  even=[i for i,e in enumerate(support) if not e[-1]];odd=[i for i,e in enumerate(support) if e[-1]]
  er=mod.matrix_from_rows(even).rank();orr=mod.matrix_from_rows(odd).rank()
  assert er+orr==mod.rank() # independently checks parity splitting
  rows.append(dict(degree=d,kernel_rank=int(mod.rank()),saturation_index_valuation=sum(len(s['replace_columns']) for s in stages),even_projection_rank=int(er),odd_projection_rank=int(orr),derivative_preserves_full_kernel=True,parity_split=True))
  certificate['degrees'].append(dict(degree=d,initial_columns=[[k,list(e)] for k,e in labels],saturation=stages))
  print('kernel derivation',d,rows[-1],round(perf_counter()-start,2),flush=True)
 archive=ROOT/'experiments/e7-quadratic-derivation-certificates.json.gz'
 archive.write_bytes(gzip.compress(json.dumps(certificate,separators=(',',':')).encode(),mtime=0))
 out=dict(result='passed',maximum_degree=13,full_flag_only=True,orbit_size=56,degrees=rows,seconds=round(perf_counter()-start,3),scope='Derivation and parity splitting on the actual integral full-flag kernel through degree 13, without relative Smith counts. No all-degree or bad-Levi invariant-lifting conclusion.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7_full_flag_kernel.py')})
 (ROOT/'experiments/e7-quadratic-derivation.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASSED',out['seconds'],flush=True)
if __name__=='__main__':main()
