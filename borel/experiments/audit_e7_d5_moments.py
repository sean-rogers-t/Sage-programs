"""Fresh-weight, exact unisolvent audit of all six D5 quotient relations.

No Sage or project imports. Equality is proved by full-rank evaluation
matrices over F_101 and exact rational identities on the same sample set.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from random import Random
from math import comb,prod
from time import perf_counter
import json,gzip
import sympy as s
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
WS=(1,1,4,5,6,8)
def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for i in range(d//ws[0]+1):
  for e in powers(ws[1:],d-i*ws[0]):yield (i,)+e

def main():
 start=perf_counter();path=ROOT/'experiments/e7-d5-images.json';primary=json.loads(path.read_text());archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 C=s.Matrix([[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]])
 ci=C[:6,:6].inv();alpha=[C.col(i) for i in range(6)];om=[sum((ci[k,j]*alpha[j] for j in range(6)),s.zeros(7,1)) for k in range(6)];h=om[0]
 xs=[-h-2*om[1]+2*om[2],-h+2*om[1]+2*om[2]-2*om[3],-h+2*om[3]-2*om[4],-h+2*om[4]-2*om[5],-h+2*om[5]];w=s.eye(7).col(6);T=s.Matrix.hstack(h,*xs,w);TI=T.inv()
 orbit={(0,0,0,0,0,0,1)};todo=list(orbit)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=tuple(v[j]-v[i]*int(C[i,j]) for j in range(7))
   if u not in orbit:orbit.add(u);todo.append(u)
 assert len(orbit)==56;forms=[tuple(mpq(int(c.p),int(c.q)) for c in TI*s.Matrix(v)) for v in sorted(orbit)]
 # A nonzero seven-by-seven trace Jacobian proves these invariants are
 # algebraically independent. Their degrees are exactly the published E7
 # basic degrees, so their polynomial subalgebra has the full Hilbert series.
 jacobian_degrees=(2,6,8,10,12,14,18);point=(2,3,5,7,11,13,17)
 jacobian=[[sum(d*prod_value**(d-1)*form[j] for form in forms for prod_value in [sum(c*x for c,x in zip(form,point))]) for j in range(7)] for d in jacobian_degrees]
 A=[list(row) for row in jacobian];det=mpq(1)
 for k in range(7):
  pivot=next(j for j in range(k,7) if A[j][k])
  if pivot!=k:A[k],A[pivot]=A[pivot],A[k];det=-det
  v=A[k][k];det*=v
  for j in range(k+1,7):
   c=A[j][k]/v;A[j]=[a-c*b for a,b in zip(A[j],A[k])]
 assert det
 # Independently establish the 27-weight E6 branching behind the primary formula.
 E6={(1,0,0,0,0,0)};todo=list(E6)
 while todo:
  v=todo.pop()
  for i in range(6):
   u=tuple(v[j]-v[i]*int(C[i,j]) for j in range(6))
   if u not in E6:E6.add(u);todo.append(u)
 proposed={tuple(h)}
 for x in xs:
  proposed.add(tuple((-h+x)/2));proposed.add(tuple((-h-x)/2))
 for mask in range(32):
  if mask.bit_count()%2==0:proposed.add(tuple((h+sum(((-1 if mask>>i&1 else 1)*xs[i] for i in range(5)),s.zeros(7,1)))/4))
 actual={tuple(sum((v[i]*om[i] for i in range(6)),s.zeros(7,1))) for v in E6};assert len(E6)==27 and actual==proposed
 e7from27={tuple(w),tuple(-w)}
 for lam in actual:
  t=s.Matrix(lam)-w/3;e7from27.add(tuple(t));e7from27.add(tuple(-t))
 assert e7from27=={tuple(s.Matrix(v)) for v in orbit}
 variables=s.symbols('h w v4 v5 v6 v8');env=dict(zip(map(str,variables),variables));degrees=(6,8,10,12,14,18)
 polynomials={d:{tuple(e):mpq(int(c.p),int(c.q)) for e,c in s.Poly(s.sympify(text.replace('^','**'),locals=env),*variables).terms()} for d,text in zip(degrees,cert['relations'])}
 exps={d:list(powers(WS,d)) for d in degrees};pivots={d:{} for d in degrees};rng=Random(314159265);samples=0
 def finite(c):return (int(c.numerator)%101)*pow(int(c.denominator)%101,-1,101)%101
 def insert(d,row):
  row=list(map(finite,row))
  for j,old in sorted(pivots[d].items()):
   if row[j]:
    scale=row[j];row=[(x-scale*y)%101 for x,y in zip(row,old)]
  pivot=next((j for j,x in enumerate(row) if x),None)
  if pivot is not None:
   inv=pow(row[pivot],-1,101);pivots[d][pivot]=[(x*inv)%101 for x in row]
 while any(len(pivots[d])<len(exps[d]) for d in degrees):
  samples+=1;assert samples<=512
  h0=rng.randrange(-19,20);xx=[rng.randrange(-23,24) for i in range(5)];sq=[v*v for v in xx]
  e1=sum(sq);e2=sum(prod(sq[i] for i in inds) for inds in combinations(range(5),2));e3=sum(prod(sq[i] for i in inds) for inds in combinations(range(5),3));e4=sum(prod(sq[i] for i in inds) for inds in combinations(range(5),4));p=prod(xx);B=-mpq(3,8)*(e1+3*h0*h0)
  values=(mpq(h0),mpq(1),mpq(e2),mpq(p),mpq(e3),mpq(e4))
  def value(e):return prod(v**k for v,k in zip(values,e))*B**(e[1]//2)
  for d in degrees:
   expected=[mpq(0),mpq(0)]
   for e,c in polynomials[d].items():expected[e[1]%2]+=c*value(e)
   actual=[mpq(0),mpq(0)]
   for form in forms:
    a=form[0]*h0+sum(c*v for c,v in zip(form[1:6],xx));b=form[6]
    for j in range(d+1):actual[j%2]+=comb(d,j)*a**(d-j)*b**j*B**(j//2)
   assert actual==expected,(samples,d)
   for parity in (0,1):insert(d,[value(e) if e[1]%2==parity else mpq(0) for e in exps[d]])
  if samples%25==0:print('exact unisolvent samples',samples,[(d,len(pivots[d])) for d in degrees],round(perf_counter()-start,2),flush=True)
 report=dict(result='passed',engine='Fresh E7 and E6 Weyl orbits, explicit rational torus change, exact rational moments and independent mod-101 Gaussian ranks. No Sage or project imports.',weight_orbit=56,E6_branching_orbit=27,samples=samples,trace_jacobian_point=list(point),trace_jacobian_determinant=str(det),basic_invariant_algebra_complete=True,identities=[dict(degree=d,rational_coefficient_dimension=len(exps[d]),evaluation_rank_mod101=len(pivots[d]),exact_identity=True) for d in degrees],scope='All six primitive characteristic-zero moment identities. Combined with the separate Buchberger/dimension audit, these prove the complete rational Borel ideal.',seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-d5-moments-audit.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED independent moment identities',report['seconds'],flush=True)
if __name__=='__main__':main()
