"""Independent E7/D6 lattice, branching and exact moment audit; no Sage."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations,product
from math import prod
from random import Random
from time import perf_counter
import gzip,json
import sympy as s
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
DEGREES=(2,6,8,10,12,14,18)
WS=(1,2,4,6,6,8,10)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def main():
 start=perf_counter();path=ROOT/'experiments/e7-d6-images.json';primary=json.loads(path.read_text());archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 C=s.Matrix([[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]])
 om=[s.eye(7).col(j) for j in range(7)];h=om[0];z=om[1]
 us=[om[6],om[5]-om[6],om[4]-om[5],om[3]-om[4],om[1]+om[2]-om[3],om[0]+om[1]-om[2]]
 assert sum(us,s.zeros(7,1))==h+2*z
 assert abs(s.Matrix.hstack(z,*us).det())==1
 vs=[u-z for u in us];QB=s.Matrix.hstack(2*z,*vs);root_coordinates=C.inv()*QB;assert abs(QB.det())==2 and root_coordinates.det() in (-1,1) and all(c.q==1 for c in root_coordinates);assert all(c.q==1 for c in C.inv()*h)
 for i in (2,3,4,5,6):
  assert z[i]==0
  assert {tuple(u-u[i]*C.col(i)) for u in us}==set(map(tuple,us))
 ys=[2*u-h for u in us]
 for i in range(1,7):
  reflected=[y-y[i]*C.col(i) for y in ys];matches=[next((j,t) for j,x in enumerate(ys) for t in (-1,1) if y==t*x) for y in reflected]
  assert len({j for j,t in matches})==6 and prod(t for j,t in matches)==1 and h[i]==0
 T=s.Matrix.hstack(h,*ys);TI=T.inv();orbit={tuple(om[6])};todo=list(orbit)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=tuple(v[j]-v[i]*int(C[j,i]) for j in range(7))
   if u not in orbit:orbit.add(u);todo.append(u)
 assert len(orbit)==56
 forms={tuple(mpq(int(c.p),int(c.q)) for c in TI*s.Matrix(v)) for v in orbit}
 proposed=set()
 for i in range(6):
  for a,b in product((-1,1),repeat=2):proposed.add((mpq(a,2),)+tuple(mpq(b,2) if j==i else mpq(0) for j in range(6)))
 for signs in product((-1,1),repeat=6):
  if prod(signs)==1:proposed.add((mpq(0),)+tuple(mpq(t,4) for t in signs))
 assert forms==proposed
 point=(2,3,5,7,11,13,17)
 J=s.Matrix([[s.Rational(str(sum(d*sum(a*x for a,x in zip(f,point))**(d-1)*f[j] for f in forms))) for j in range(7)] for d in DEGREES]);det=J.det();assert det
 variables=s.symbols('h v2 v4 v6 p6 v8 v10');env=dict(zip(map(str,variables),variables))
 fs={d:{tuple(e):mpq(int(c.p),int(c.q)) for e,c in s.Poly(s.sympify(cert['full_moments'][str(d)].replace('^','**'),locals=env),*variables).terms()} for d in DEGREES}
 full={d:s.sympify(cert['full_moments'][str(d)].replace('^','**'),locals=env) for d in DEGREES};hv,v2=variables[:2]
 assert s.expand(full[2]-3*(v2+2*hv**2))==0
 for d,text in zip(DEGREES[1:],cert['relations']):assert s.expand(full[d].subs(v2,-2*hv**2)-s.sympify(text.replace('^','**'),locals=env))==0
 assert all(sum(a*b for a,b in zip(e,WS))==d for d in DEGREES for e in fs[d])
 exps={d:list(powers(WS,d)) for d in DEGREES};pivots={d:{} for d in DEGREES};rng=Random(2718281828);samples=0
 while any(len(pivots[d])<len(exps[d]) for d in DEGREES):
  samples+=1;assert samples<1000
  hv=rng.randrange(-31,32);yy=[rng.randrange(-37,38) for j in range(6)];sq=[y*y for y in yy]
  elem={k:sum(prod(sq[i] for i in ix) for ix in combinations(range(6),k)) for k in range(1,6)}
  values=(hv,elem[1],elem[2],elem[3],prod(yy),elem[4],elem[5]);weightvalues=[sum(a*x for a,x in zip(f,(hv,*yy))) for f in forms]
  for d in DEGREES:
   assert sum(c*prod(v**k for v,k in zip(values,e)) for e,c in fs[d].items())==sum(v**d for v in weightvalues),(samples,d)
   row=[prod(pow(v%101,k,101) for v,k in zip(values,e))%101 for e in exps[d]]
   for j,old in sorted(pivots[d].items()):
    if row[j]:
     a=row[j];row=[(x-a*y)%101 for x,y in zip(row,old)]
   j=next((j for j,x in enumerate(row) if x),None)
   if j is not None:
    a=pow(row[j],-1,101);pivots[d][j]=[a*x%101 for x in row]
 print('exact fresh-weight identities:',samples,'unisolvent samples',flush=True)
 out=dict(result='passed',engine='Native Python/SymPy/gmpy2; fresh Weyl orbit, exact rational evaluations and independent finite-field ranks. No Sage or project imports.',weight_orbit=56,weight_lattice_basis_determinant=1,root_lattice_basis_index=2,permutation_subgroup='A5 = S6',levi_action='even signed permutations of six coordinates',trace_jacobian_point=list(point),trace_jacobian_determinant=str(det),basic_invariant_algebra_complete=True,samples=samples,identities=[dict(degree=d,coefficient_space_dimension=len(exps[d]),evaluation_rank_mod101=len(pivots[d])) for d in DEGREES],seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-d6-moments-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],flush=True)

if __name__=='__main__':main()
