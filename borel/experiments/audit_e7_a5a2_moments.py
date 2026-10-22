"""Fresh E7 orbit, odd-index lattice and exact quotient-moment audit."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from functools import lru_cache
from random import Random
from math import prod
from time import perf_counter
import gzip,json
import sympy as s
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
DS=(2,6,8,10,12,14,18);WS=(2,3,3,4,5,6)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def main():
 start=perf_counter();ip=ROOT/'experiments/e7-a5a2-images.json';ir=json.loads(ip.read_text());ia=ROOT/ir['certificate_archive'];assert sha256(ia.read_bytes()).hexdigest()==ir['certificate_archive_sha256'];ic=json.loads(gzip.decompress(ia.read_bytes()));C=2*s.eye(7)
 for i,j in ((1,3),(3,4),(2,4),(4,5),(5,6),(6,7)):C[i-1,j-1]=C[j-1,i-1]=-1
 chain=(1,3,4,5,6);A=C.extract(chain,chain);roots=[C.col(i) for i in chain];theta=C*s.Matrix([2,2,3,4,3,2,1]);a2=[C.col(0),-theta];eta=[sum((s.Matrix([[2,-1],[-1,2]]).inv()[i,j]*a2[j] for j in range(2)),s.zeros(7,1)) for i in range(2)];omega=[sum((A.inv()[i,j]*roots[j] for j in range(5)),s.zeros(7,1)) for i in range(5)];xx=[omega[0]]+[omega[i]-omega[i-1] for i in range(1,5)]+[-omega[4]];zz=[eta[0],eta[1]-eta[0],-eta[1]]
 assert abs((C.inv()*s.Matrix.hstack(*roots,*a2)).det())==3;Q=C.inv()*s.Matrix.hstack(*roots,*eta);assert all(int(c.q)%2 for c in Q) and int(Q.det().p)%2 and int(Q.det().q)%2;X=s.Matrix.hstack(*xx[:5],*eta);assert int(X.det().p)%2 and int(X.det().q)%2
 orbit={(0,0,0,0,0,0,1)};todo=list(orbit)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=tuple(v[j]-v[i]*C[j,i] for j in range(7))
   if u not in orbit:orbit.add(u);todo.append(u)
 ext={tuple(sum((xx[i] for i in ids),s.zeros(7,1))) for ids in combinations(range(6),3)};branch=ir['branch_sign'];assert orbit==ext|{tuple(eps*(x+branch*z)) for x in xx for z in zz for eps in (-1,1)} and len(orbit)==56
 forms=[tuple(mpq(int(c.p),int(c.q)) for c in X.inv()*s.Matrix(v)) for v in orbit];point=(2,3,5,7,11,13,17);J=s.Matrix([[s.Rational(str(sum(d*sum(a*x for a,x in zip(f,point))**(d-1)*f[j] for f in forms))) for j in range(7)] for d in DS]);det=J.det();assert det
 vars=s.symbols('q2 q3 e3 e4 e5 e6');env=dict(zip(map(str,vars),vars));fs={d:{tuple(e):mpq(int(c.p),int(c.q)) for e,c in s.Poly(s.sympify(text.replace('^','**'),locals=env),*vars).terms()} for d,text in zip(DS[1:],ic['relations'])};exps={d:list(powers(WS,d)) for d in DS[1:]};ranks={d:{} for d in exps};rng=Random(1618033988);samples=0
 while any(len(ranks[d])<len(exps[d]) for d in exps):
  samples+=1;assert samples<1000;xs=[rng.randrange(-19,20) for _ in range(5)];xs.append(-sum(xs));h=rng.randrange(-19,20);es={k:sum(prod(xs[j] for j in ids) for ids in combinations(range(6),k)) for k in range(2,7)};t=h*h-es[2];values=(-es[2],-h*t,es[3],es[4],es[5],es[6]);powersum={d:[mpq(0),mpq(0)] for d in DS}
  for form in forms:
   a=sum(c*v for c,v in zip(form[:5],xs[:5]))+form[6]*h;b=form[5];v0=mpq(1);v1=mpq(0)
   for d in range(1,19):
    v0,v1=a*v0-b*t*v1,b*v0+(a+b*h)*v1
    if d in powersum:powersum[d][0]+=v0;powersum[d][1]+=v1
  assert powersum[2]==[0,0]
  for d in exps:
   expected=sum(c*prod(v**k for v,k in zip(values,e)) for e,c in fs[d].items());assert powersum[d]==[expected,0],(samples,d);row=[prod(pow(v%101,k,101) for v,k in zip(values,e))%101 for e in exps[d]]
   for j,old in sorted(ranks[d].items()):
    if row[j]:
     a=row[j];row=[(x-a*y)%101 for x,y in zip(row,old)]
   j=next((j for j,c in enumerate(row) if c),None)
   if j is not None:
    a=pow(row[j],-1,101);ranks[d][j]=[a*c%101 for c in row]
 out=dict(result='passed',engine='Native SymPy/gmpy2; fresh E7 orbit and subsystem lattices, independent exact quadratic-algebra evaluations and finite-field unisolvence. No Sage or primary project imports.',subsystem_index=3,root_and_weight_splittings_at_two=True,weight_orbit=56,branch_sign=branch,trace_jacobian_determinant=str(det),basic_global_invariants_complete=True,samples=samples,identities=[dict(degree=d,coefficient_space_dimension=len(exps[d]),rank_mod101=len(ranks[d])) for d in exps],A2_common_fiber_bases=dict(A1=['1','h','h^2'],trivial=['1','h','h^2','u','u*h','u*h^2']),A2_relations=['h^3+q2*h+q3=0','u^2-h*u+q2+h^2=0'],seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ip)},raw_sha256={ia.relative_to(ROOT).as_posix():sha256(ia.read_bytes()).hexdigest()});(ROOT/'experiments/e7-a5a2-moments-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED A5+A2 moments',samples,out['seconds'],flush=True)

if __name__=='__main__':main()
