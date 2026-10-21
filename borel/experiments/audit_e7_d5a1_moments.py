"""Independent E7/P6 coordinates and exact fresh-weight trace identities.

Native SymPy/gmpy2 only. The full weight source A[t] reuses the audited
D5 algebra A; this checks the actual integral lattice and identical action.
"""
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
WS=(1,2,2,4,5,6,8)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def main():
 start=perf_counter();path=ROOT/'experiments/e7-d5a1-images.json';primary=json.loads(path.read_text());archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 C=2*s.eye(7)
 for i,j in ((1,3),(3,4),(2,4),(4,5),(5,6),(6,7)):C[i-1,j-1]=C[j-1,i-1]=-1
 om=[s.eye(7).col(j) for j in range(7)];h,w,z=om[5],om[6],om[4];us=[om[0],om[2]-om[0],om[3]-om[2],om[1]+om[4]-om[3],om[4]-om[1]];ys=[2*u-h for u in us]
 assert sum(us,s.zeros(7,1))==2*z;assert abs(s.Matrix.hstack(h,w,z,*us[:4]).det())==1
 QB=s.Matrix.hstack(h,2*z,w-z,*us[:4]);RC=C.inv()*QB;assert abs(RC.det())==1 and all(c.q==1 for c in RC)
 assert abs(s.Matrix.hstack(h,w,z,*(u-us[4] for u in us[:4])).det())==5
 for i in (0,1,2,3):
  assert {tuple(u-u[i]*C.col(i)) for u in us}==set(map(tuple,us));assert h[i]==w[i]==z[i]==0
 for i in range(5):
  matches=[next((j,t) for j,y in enumerate(ys) for t in (-1,1) if v-v[i]*C.col(i)==t*y) for v in ys];assert len({j for j,t in matches})==5 and prod(t for j,t in matches)==1;assert h[i]==w[i]==0
 assert all(v[6]==0 for v in [h,z,*us,*ys]);assert w-w[6]*C.col(6)==h-w
 p=w-z;assert p-p[6]*C.col(6)==h-2*z-p
 T=s.Matrix.hstack(h,*ys,w);TI=T.inv();orbit={tuple(om[6])};todo=list(orbit)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=tuple(v[j]-v[i]*int(C[j,i]) for j in range(7))
   if u not in orbit:orbit.add(u);todo.append(u)
 assert len(orbit)==56;forms={tuple(mpq(int(c.p),int(c.q)) for c in TI*s.Matrix(v)) for v in orbit};proposed=set()
 for a in (-1,1):
  proposed.add((mpq(0),)+(mpq(0),)*5+(mpq(a),));proposed.add((mpq(a),)+(mpq(0),)*5+(mpq(-a),))
 for i in range(5):
  for a,b in product((-1,1),repeat=2):proposed.add((mpq(-a,2),)+tuple(mpq(b,2) if j==i else mpq(0) for j in range(5))+(mpq(a),))
 for signs in product((-1,1),repeat=5):
  if prod(signs)==1:
   for a in (-1,1):proposed.add((mpq(a,4),)+tuple(mpq(a*t,4) for t in signs)+(mpq(0),))
 assert forms==proposed
 point=(2,3,5,7,11,13,17);J=s.Matrix([[s.Rational(str(sum(d*sum(a*x for a,x in zip(f,point))**(d-1)*f[j] for f in forms))) for j in range(7)] for d in DEGREES]);det=J.det();assert det
 variables=s.symbols('h t v2 v4 p5 v6 v8');env=dict(zip(map(str,variables),variables));full={d:s.sympify(cert['full_moments'][str(d)].replace('^','**'),locals=env) for d in DEGREES};fs={d:{tuple(e):mpq(int(c.p),int(c.q)) for e,c in s.Poly(full[d],*variables).terms()} for d in DEGREES};hv,tv,v2=variables[:3]
 assert s.expand(full[2]-3*(v2+3*hv**2-8*tv))==0
 for d,text in zip(DEGREES[1:],cert['relations']):assert s.expand(full[d].subs(v2,8*tv-3*hv**2)-s.sympify(text.replace('^','**'),locals=env))==0
 assert all(sum(a*b for a,b in zip(e,WS))==d for d in DEGREES for e in fs[d])
 exps={d:list(powers(WS,d)) for d in DEGREES};pivots={d:{} for d in DEGREES};rng=Random(6180339887);samples=0
 while any(len(pivots[d])<len(exps[d]) for d in DEGREES):
  samples+=1;assert samples<1500;hv=rng.randrange(-31,32);ww=rng.randrange(-31,32);yy=[rng.randrange(-37,38) for j in range(5)];sq=[y*y for y in yy];elem={k:sum(prod(sq[i] for i in ix) for ix in combinations(range(5),k)) for k in range(1,5)};values=(hv,ww*(hv-ww),elem[1],elem[2],prod(yy),elem[3],elem[4]);weightvalues=[sum(a*x for a,x in zip(f,(hv,*yy,ww))) for f in forms]
  for d in DEGREES:
   assert sum(c*prod(v**k for v,k in zip(values,e)) for e,c in fs[d].items())==sum(v**d for v in weightvalues),(samples,d)
   row=[prod(pow(v%101,k,101) for v,k in zip(values,e))%101 for e in exps[d]]
   for j,old in sorted(pivots[d].items()):
    if row[j]:
     a=row[j];row=[(x-a*y)%101 for x,y in zip(row,old)]
   j=next((j for j,x in enumerate(row) if x),None)
   if j is not None:
    a=pow(row[j],-1,101);pivots[d][j]=[a*x%101 for x in row]
 print('fresh E7/P6 moments:',samples,'unisolvent exact samples',flush=True)
 shared=[ROOT/'experiments/e7-d5-source-algebra.json',ROOT/'experiments/e7-d5-source-audit.json'];assert all(json.loads(p.read_text())['result']=='passed' for p in shared)
 out=dict(result='passed',engine='Native SymPy/gmpy2: actual lattice/action, fresh 56-weight Weyl orbit, exact unisolvent evaluations and independent ranks; no Sage or project imports.',weight_orbit=56,weight_basis_determinant=1,root_basis_determinant_in_roots=1,centered_A4_basis_index=5,coefficient_rings=dict(N='R[h,z,c2,c3,c4,c5,t], e1(u)=2z',M='R[h,b,c2,c3,c4,c5,r], e1(u)=b'),full_weight_source='Previously audited D5 algebra A with a polynomial variable t=w(h-w); D5 fixes w, and s7 fixes every D5 coefficient and sends w to h-w.',trace_jacobian_point=list(point),trace_jacobian_determinant=str(det),basic_invariant_algebra_complete=True,samples=samples,identities=[dict(degree=d,coefficient_space_dimension=len(exps[d]),evaluation_rank_mod101=len(pivots[d])) for d in DEGREES],seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path,*shared)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()});(ROOT/'experiments/e7-d5a1-moments-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],flush=True)

if __name__=='__main__':main()
