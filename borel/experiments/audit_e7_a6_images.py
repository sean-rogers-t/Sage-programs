"""Plain Python audit of complete E7/A6 local image bases.

No Sage or project imports. Uses a sparse rational weighted reducer,
Buchberger checks, Gaussian inverses, and mod-2 bit ranks. The separate
adjoint-moment audit checks the rational ideal via a different input.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from functools import lru_cache
from time import perf_counter
import gzip,json
import sympy as s
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
WS=(1,3,4,5,6,7)
ZERO=(0,)*6

def rank2(rows):
 pivots={}
 for row0 in rows:
  row=sum((int(c.numerator)%2)<<j for j,c in enumerate(row0))
  while row:
   k=row.bit_length()-1
   if k in pivots:row^=pivots[k]
   else:pivots[k]=row;break
 return len(pivots)

def inverse(A):
 n=len(A);B=[list(row)+[mpq(int(i==j)) for j in range(n)] for i,row in enumerate(A)]
 for k in range(n):
  pivot=next(j for j in range(k,n) if B[j][k]);B[k],B[pivot]=B[pivot],B[k]
  c=B[k][k];B[k]=[v/c for v in B[k]]
  for j in range(n):
   if j!=k and B[j][k]:
    c=B[j][k];B[j]=[v-c*u for v,u in zip(B[j],B[k])]
 return [row[n:] for row in B]

def matmul(A,B):
 return [[sum((a*b for a,b in zip(row,col)),mpq(0)) for col in zip(*B)] for row in A]

def local(A):return all(int(c.denominator)%2 for row in A for c in row)

def main():
 start=perf_counter();primary_path=ROOT/'experiments/e7-a6-images.json';ideal_path=ROOT/'experiments/e7-a6-ideal-audit.json'
 primary=json.loads(primary_path.read_text());ideal=json.loads(ideal_path.read_text());assert primary['result']==ideal['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 variables=s.symbols('c1 c3 c4 c5 c6 c7');env=dict(zip(('c1','c3','c4','c5','c6','c7'),variables))
 def parse(text):
  f=s.Poly(s.sympify(text.replace('^','**'),locals=env),*variables,domain=s.QQ)
  return {tuple(e):mpq(int(c.p),int(c.q)) for e,c in f.terms()}
 def key(e):return (sum(w*j for w,j in zip(WS,e)),tuple(-j for j in reversed(e)))
 G=[parse(f) for f in cert['groebner_basis']];leads=[max(f,key=key) for f in G]
 assert all(len({key(e)[0] for e in f})==1 for f in G)
 reducers=[(e,{u:-c/f[e] for u,c in f.items() if u!=e}) for e,f in zip(leads,G)]
 @lru_cache(None)
 def normal_monomial(e):
  # Every monomial above degree 42 lies in the leading ideal, as verified
  # below by exhaustive enumeration of its 576 standard monomials.
  if key(e)[0]>42:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     exponent=tuple(a+b for a,b in zip(u,shift))
     for v,d in normal_monomial(exponent).items():out[v]=out.get(v,mpq(0))+c*d
    return {u:c for u,c in out.items() if c}
  return {e:mpq(1)}
 def normal(f):
  out={}
  for e,c in f.items():
   for u,d in normal_monomial(e).items():out[u]=out.get(u,mpq(0))+c*d
  return {u:c for u,c in out.items() if c}
 # Leading-ideal enumeration is independent of the recorded monomial list.
 standard={ZERO};pending=[ZERO]
 while pending:
  e=pending.pop()
  for j in range(6):
   v=list(e);v[j]+=1;v=tuple(v)
   if v not in standard and not any(all(a>=b for a,b in zip(v,u)) for u in leads):standard.add(v);pending.append(v)
 assert len(standard)==576 and max(key(e)[0] for e in standard)==42
 bydegree={d:[tuple(e) for e in cert['normal_exponents'][str(d)]] for d in range(43)}
 assert set().union(*(set(es) for es in bydegree.values()))==standard
 assert all(key(e)[0]==d for d,es in bydegree.items() for e in es)
 # Buchberger's criterion, skipping relatively prime leading monomials.
 pairs=0
 for i in range(len(G)):
  for j in range(i):
   a,b=leads[i],leads[j]
   if not any(x and y for x,y in zip(a,b)):continue
   lcm=tuple(max(x,y) for x,y in zip(a,b))
   if key(lcm)[0]>42:continue
   shifts=[tuple(x-y for x,y in zip(lcm,lead)) for lead in (a,b)];poly={}
   for f,lead,shift,sign in ((G[i],a,shifts[0],1),(G[j],b,shifts[1],-1)):
    for e,c in f.items():
     u=tuple(x+y for x,y in zip(e,shift));poly[u]=poly.get(u,mpq(0))+sign*c/f[lead]
   assert not normal(poly),(i,j)
   pairs+=1
 print('independent Buchberger/standard basis',pairs,'pairs',round(perf_counter()-start,2),flush=True)
 assert all(not normal(parse(f)) for f in cert['relations'])
 # Independently check the actual torus/root lattices and the weight sets.
 C=s.Matrix([[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]])
 T=s.Matrix([[-1,1,0,0,0,0,0],[1,1,-1,0,0,0,0],[0,1,1,-1,0,0,0],[0,0,0,1,-1,0,0],[0,0,0,0,1,-1,0],[0,0,0,0,0,1,-1],[0,0,0,0,0,0,1]])
 assert abs(T.det())==3
 ti=[T.row(j).T for j in range(7)];center=sum(ti,s.zeros(7,1))/7;V=[ti[j]-ti[6] for j in range(6)]
 XB=s.Matrix.hstack(*V,center);QB=C.inv()*s.Matrix.hstack(*V,2*center)
 assert abs(XB.det())==abs(QB.det())==3 and all(s.denom(v)%2 for v in list(XB)+list(QB))
 trace_center=sum(ti,s.zeros(7,1))/3
 minweights=set();roots=set()
 def put(S,v):S.add(tuple(map(int,v)));S.add(tuple(map(int,-v)))
 for v in ti:put(minweights,v);put(roots,v-trace_center)
 for u,v in combinations(ti,2):put(minweights,u+v-trace_center);put(roots,u-v)
 for u,v,w in combinations(ti,3):put(roots,u+v+w-trace_center)
 assert len(minweights)==56 and len(roots)==126
 def orbit_check(weights,start):
  seen={start};todo=[start]
  while todo:
   v=todo.pop()
   for i in range(7):
    u=tuple(v[j]-v[i]*int(C[i,j]) for j in range(7));assert u in weights
    if u not in seen:seen.add(u);todo.append(u)
  assert seen==weights
 orbit_check(minweights,(0,0,0,0,0,0,1));orbit_check(roots,tuple(map(int,C.row(0))))
 # Rebuild both source generator lists directly from elementary identities.
 c1,c3,c4,c5,c6,c7=variables;es=[s.Integer(1),c1,s.Rational(4,9)*c1*c1,c3,c4,c5,c6,c7];mu=c1/7
 source_expr={'N':[(es[i],i) for i in range(1,8)],'M':[(2*mu,1)]+[(sum(s.binomial(7-j,k-j)*(-mu)**(k-j)*es[j] for j in range(k+1)),k) for k in range(2,8)]}
 sources={kind:[(parse(str(f)),d) for f,d in gs] for kind,gs in source_expr.items()}
 assert sources=={kind:[(parse(f),d) for f,d in gs] for kind,gs in cert['source_generators'].items()}
 proposed={kind:{0:[[mpq(1)]]} for kind in sources};inverses={kind:{0:[[mpq(1)]]} for kind in sources};records=[]
 @lru_cache(None)
 def multiplication(kind,index,degree):
  f,k=sources[kind][index];target=bydegree[degree];ix={e:i for i,e in enumerate(target)};columns=[]
  for e in bydegree[degree-k]:
   out=normal({tuple(a+b for a,b in zip(u,e)):c for u,c in f.items()});col=[mpq(0)]*len(target)
   for u,c in out.items():col[ix[u]]=c
   columns.append(col)
  return list(map(list,zip(*columns)))
 for entry in cert['image_degrees']:
  degree=entry['degree'];n=len(bydegree[degree])
  for kind in sources:
   H=[[mpq(c) for c in row] for row in entry['weight_basis' if kind=='N' else 'root_basis']];HI=inverse(H);proposed[kind][degree]=H;inverses[kind][degree]=HI
   products=[]
   for i,(f,k) in enumerate(sources[kind]):
    if k>degree:continue
    values=matmul(multiplication(kind,i,degree),proposed[kind][degree-k]);coords=matmul(HI,values);assert local(coords)
    products.extend(list(col) for col in zip(*coords))
   assert rank2(list(map(list,zip(*products))))==n
  transition=matmul(inverses['N'][degree],proposed['M'][degree]);assert local(transition)
  twice_inverse=[[2*c for c in row] for row in inverse(transition)];assert local(twice_inverse)
  m=n-rank2(transition);assert m==entry['multiplicity_2']==primary['multiplicity_coefficients'][degree]
  records.append(dict(degree=degree,rank=n,multiplicity_2=m,both_image_closures_and_primitive_generation=True,twice_weight_image_in_root_image=True))
  if degree%5==0 or degree==42:print('independent complete images degree',degree,round(perf_counter()-start,2),flush=True)
 output=dict(result='passed',engine='Plain Python/gmpy2 sparse weighted reduction and Buchberger checks; own rational inverses and bit ranks. SymPy parses inputs and verifies small Cartan identities. No Sage or project imports.',degree_rows=records,total_rank=576,dimension=42,buchberger_pairs_checked=pairs,weight_basis_index=3,root_basis_index=3,verified_weight_orbit=56,verified_root_orbit=126,seconds=round(perf_counter()-start,3),scope='Full source lattice split, rational standard basis, every proposed image basis by closure and primitive generation, and every local comparison. The rational-ideal bridge is separately checked with fresh adjoint moments; that bridge shares Singular arithmetic with primary.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path,ideal_path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-a6-images-audit.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED',output['seconds'],flush=True)
if __name__=='__main__':main()
