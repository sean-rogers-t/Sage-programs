"""Independent full D4 source audit: no Sage and no project imports.

Checks the coefficient coordinate lattice, rational invariant identities,
GF(2) graph ideal by substitution and Buchberger, and its Hilbert numerator.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb
from time import perf_counter
import gzip,json,sys
import sympy as s
from gmpy2 import mpq
sys.setrecursionlimit(20000)
ROOT=Path(__file__).resolve().parents[1]
Z=(0,)*6

def add(f,g,c=1):
 out=dict(f)
 for e,v in g.items():out[e]=out.get(e,mpq(0))+c*v
 return {e:v for e,v in out.items() if v}
def mul(f,g):
 out={}
 for e,a in f.items():
  for u,b in g.items():
   v=tuple(x+y for x,y in zip(e,u));out[v]=out.get(v,mpq(0))+a*b
 return {e:v for e,v in out.items() if v}
def scale(f,c):return {e:v*c for e,v in f.items() if v*c}
def power(f,k):
 out={Z:mpq(1)}
 while k:
  if k%2:out=mul(out,f)
  f=mul(f,f);k//=2
 return out
def parse(text,variables,mod2=False):
 f=s.Poly(s.sympify(text.replace('^','**'),locals=dict(zip(map(str,variables),variables))),*variables,domain=s.QQ)
 if mod2:return frozenset(tuple(e) for e,c in f.terms() if int(c)%2)
 return {tuple(e):mpq(int(c.p),int(c.q)) for e,c in f.terms()}
def bitmul(f,g):
 out=set()
 for e in f:
  for u in g:
   v=tuple(x+y for x,y in zip(e,u))
   if v in out:out.remove(v)
   else:out.add(v)
 return frozenset(out)
def bitpower(f,k,n):
 out=frozenset([(0,)*n])
 while k:
  if k%2:out=bitmul(out,f)
  f=frozenset(tuple(2*x for x in e) for e in f);k//=2
 return out

def main():
 start=perf_counter();path=ROOT/'experiments/e7-d4-source-algebra.json';source=json.loads(path.read_text());assert source['result']=='passed'
 towerpath=ROOT/'experiments/e7-d4-tower-certificates.json.gz';tower=json.loads(gzip.decompress(towerpath.read_bytes()))
 C=s.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
 om=[s.eye(6).row(i) for i in range(6)];hh=om[0];zz=om[2];us=[-om[1]+om[2],om[1]+om[2]-om[3],om[3]-om[4],om[4]-om[5],om[5]];yy=us[4];vv=us[3]
 assert C.det()==3 and sum(us,s.zeros(1,6))==2*zz
 assert abs(s.Matrix.vstack(hh,yy,zz,vv,us[0],us[1]).det())==1
 assert abs(s.Matrix.vstack(hh,yy,zz,vv,us[0]-us[2],us[1]-us[2]).det())==3
 for i in (1,3):
  assert set(map(tuple,[u-u[i]*C.row(i) for u in us[:3]]))==set(map(tuple,us[:3]))
  assert all(a[i]==0 for a in (hh,yy,zz,vv))
 xx=[2*u-hh for u in us[:4]]
 for i in (1,2,3,4):
  matches=[next((j,k) for j,x in enumerate(xx) for k in (-1,1) if a-a[i]*C.row(i)==k*x) for a in xx]
  assert len(set(j for j,k in matches))==4 and __import__('math').prod(k for j,k in matches)==1
  assert hh[i]==yy[i]==0
 variables=s.symbols('h y z v c2 c3');units=[{tuple(int(i==j) for i in range(6)):mpq(1)} for j in range(6)];h,y,z,v,c2,c3=units;one={Z:mpq(1)}
 es=[one,add(add(scale(z,2),v,-1),y,-1),c2,c3];ps=[scale(one,3)]
 for k in range(1,7):
  f={}
  for j in range(1,min(k,4)):f=add(f,mul(es[j],ps[k-j]),(-1)**(j-1))
  if k<=3:f=add(f,es[k],(-1)**(k-1)*k)
  ps.append(f)
 xp={}
 for k in (2,4,6):
  f=power(add(scale(v,2),h,-1),k)
  for j in range(k+1):f=add(f,mul(power(h,k-j),ps[j]),comb(k,j)*2**j*(-1)**(k-j))
  xp[k]=f
 sq=[one]
 for k in range(1,4):
  f={}
  for j in range(1,k+1):f=add(f,mul(sq[k-j],xp[2*j]),(-1)**(j-1))
  sq.append(scale(f,mpq(1,k)))
 eu={}
 for j in range(4):eu=add(eu,mul(power(h,3-j),es[j]),2**j*(-1)**(3-j))
 eu=mul(add(scale(v,2),h,-1),eu)
 original=[h,y,sq[1],sq[2],eu,sq[3]];divisors=[]
 for i,f in enumerate(original):
  g=__import__('functools').reduce(__import__('math').gcd,(abs(int(c)) for c in f.values()));divisors.append(g);original[i]=scale(f,mpq(1,g))
 @lru_cache(None)
 def monomial(e):
  if not any(e):return one
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return mul(original[j],monomial(tuple(u)))
 for item,recorded in zip(tower['source_preimages'],source['integral_generators']):
  f={}
  for e,c in zip(item['exponents'],item['coefficients']):f=add(f,monomial(tuple(e)),mpq(c))
  assert f==parse(recorded,variables)==parse(item['source'],variables)
  assert all(c.denominator==1 for c in f.values())
 assert len(tower['source_preimages'])==10
 print('independent D4 lattice and ten rational invariant identities',round(perf_counter()-start,2),flush=True)
 vs=s.symbols(' '.join(source['elimination_variables']));ws=source['elimination_grading'];n=len(vs)
 G=[parse(f,vs,True) for f in source['elimination_groebner']];relations=[parse(f,vs,True) for f in source['elimination_relations']]
 def key(e):return (sum(a*b for a,b in zip(e[:3],ws[:3])),tuple(-x for x in reversed(e[:3])),sum(a*b for a,b in zip(e[3:],ws[3:])),tuple(-x for x in reversed(e[3:])))
 leads=[max(f,key=key) for f in G]
 @lru_cache(None)
 def nm(e):
  for f,l in zip(G,leads):
   if all(a>=b for a,b in zip(e,l)):
    shift=tuple(a-b for a,b in zip(e,l));out=frozenset()
    for u in f:
     if u!=l:out=out.symmetric_difference(nm(tuple(a+b for a,b in zip(u,shift))))
    return out
  return frozenset([e])
 def normal(f):
  out=frozenset()
  for e in f:out=out.symmetric_difference(nm(e))
  return out
 assert all(not normal(f) for f in relations)
 pairs=0;chain_pairs=0;resolved=set();queue=[]
 for i in range(len(G)):
  for j in range(i):
   lcm=tuple(max(x,y) for x,y in zip(leads[i],leads[j]));queue.append((key(lcm),i,j,lcm))
 for _,i,j,lcm in sorted(queue):
  a,b=leads[i],leads[j];pair=(i,j)
  if not any(x and y for x,y in zip(a,b)):resolved.add(pair);continue
  witness=next((k for k,c in enumerate(leads) if k not in (i,j) and all(x<=y for x,y in zip(c,lcm)) and tuple(max(x,y) for x,y in zip(a,c))!=lcm and tuple(max(x,y) for x,y in zip(b,c))!=lcm and (max(i,k),min(i,k)) in resolved and (max(j,k),min(j,k)) in resolved),None)
  if witness is not None:chain_pairs+=1;resolved.add(pair);continue
  sa=tuple(x-y for x,y in zip(lcm,a));sb=tuple(x-y for x,y in zip(lcm,b))
  f=frozenset(tuple(x+y for x,y in zip(e,sa)) for e in G[i]).symmetric_difference(tuple(x+y for x,y in zip(e,sb)) for e in G[j]);assert not normal(f),(i,j);pairs+=1;resolved.add(pair)
 dependent=list(range(6,13));free=list(range(6));replacements={i:frozenset([tuple(int(j==i) for j in range(6))]) for i in free}
 for f in relations:
  candidates=[i for i in dependent if tuple(int(i==j) for j in range(n)) in f];assert len(candidates)==1;i=candidates[0];tail=f-{tuple(int(i==j) for j in range(n))};assert all(not any(e[j] for j in dependent) for e in tail);replacements[i]=frozenset(e[:6] for e in tail)
 def evaluate(f,inputs):
  out=frozenset()
  for e in f:
   term=frozenset([Z])
   for a,k in zip(inputs,e):
    if k:term=bitmul(term,bitpower(a,k,6))
   out=out.symmetric_difference(term)
  return out
 for f in G:assert not evaluate(f,[replacements[i] for i in range(n)])
 z0,v0,b0,h0,y0,a0=[replacements[i] for i in range(6)]
 c20=a0.symmetric_difference(bitmul(h0,z0)).symmetric_difference(bitpower(v0,2,6)).symmetric_difference(bitmul(v0,y0))
 c30=b0.symmetric_difference(bitmul(c20,v0)).symmetric_difference(bitmul(bitpower(h0,2,6),z0)).symmetric_difference(bitmul(bitmul(h0,y0),z0)).symmetric_difference(bitmul(h0,bitpower(z0,2,6)))
 for i,f in enumerate(source['integral_generators']):assert evaluate(parse(f,variables,True),[h0,y0,z0,v0,c20,c30])==replacements[3+i]
 pure=[f for f in G if all(not any(e[:3]) for e in f)];assert len(pure)==37
 assert {frozenset(e[3:] for e in f) for f in pure}=={parse(f,vs[3:],True) for f in source['kernel']}
 numerator={(0,)*10:1}
 for f in pure:
  l=max(f,key=key)[3:];new=dict(numerator)
  for e,c in numerator.items():u=tuple(max(a,b) for a,b in zip(e,l));new[u]=new.get(u,0)-c
  numerator={e:c for e,c in new.items() if c}
 univariate={}
 for e,c in numerator.items():d=sum(a*b for a,b in zip(e,ws[3:]));univariate[d]=univariate.get(d,0)+c
 univariate={d:c for d,c in univariate.items() if c};expected={0:1}
 for d in (5,8,8,8):
  nxt=dict(expected)
  for k,c in expected.items():nxt[k+d]=nxt.get(k+d,0)-c
  expected={k:c for k,c in nxt.items() if c}
 assert univariate==expected
 print('independent graph Buchberger and full-source Hilbert series',pairs,round(perf_counter()-start,2),flush=True)
 # Rebuild D5 embedding from the five elementary symmetric functions.
 esbits=[frozenset([Z]),z0.symmetric_difference(z0).symmetric_difference(v0).symmetric_difference(y0),c20,c30];tail=[frozenset([Z]),v0.symmetric_difference(y0),bitmul(v0,y0)];full=[]
 for k in range(6):
  f=frozenset()
  for j in range(4):
   if 0<=k-j<=2:f=f.symmetric_difference(bitmul(esbits[j],tail[k-j]))
  full.append(f)
 previous=json.loads((ROOT/'experiments/e7-d5-source-algebra.json').read_text());dv=s.symbols('h z c2 c3 c4 c5')
 embedded=[parse(f,vs[3:],True) for f,d in source['relative_D5_generators']]
 for f,actual in zip(previous['integral_generators'],embedded):assert evaluate(parse(f,dv,True),[h0,z0]+full[2:])==evaluate(actual,[replacements[i] for i in range(3,n)])
 # Independently derive the relative ideal by Buchberger, starting from
 # the 37 kernel relations and the ten freshly checked D5 restrictions.
 qkey=lambda e:(sum(a*b for a,b in zip(e,ws[3:])),tuple(-x for x in reversed(e)))
 basis=[];bl=[]
 def reduce_bits(f):
  f=set(f);out=set()
  while f:
   e=max(f,key=qkey);found=False
   for g,l in zip(basis,bl):
    if all(a>=b for a,b in zip(e,l)):
     shift=tuple(a-b for a,b in zip(e,l));f.symmetric_difference_update(tuple(a+b for a,b in zip(u,shift)) for u in g);found=True;break
   if not found:f.remove(e);out.add(e)
  return frozenset(out)
 queue=[]
 def append(f):
  f=reduce_bits(f)
  if not f:return
  i=len(basis);l=max(f,key=qkey);basis.append(f);bl.append(l)
  for j,b in enumerate(bl[:-1]):
   if any(x and y for x,y in zip(l,b)):queue.append((i,j))
 for f in sorted([frozenset(e[3:] for e in f) for f in pure]+embedded,key=lambda f:qkey(max(f,key=qkey))):append(f)
 qpairs=0
 while queue:
  i,j=queue.pop(0);a,b=bl[i],bl[j];l=tuple(max(x,y) for x,y in zip(a,b));sa=tuple(x-y for x,y in zip(l,a));sb=tuple(x-y for x,y in zip(l,b));f=frozenset(tuple(x+y for x,y in zip(e,sa)) for e in basis[i]).symmetric_difference(tuple(x+y for x,y in zip(e,sb)) for e in basis[j]);append(f);qpairs+=1
 QG=[parse(f,vs[3:],True) for f in source['relative_quotient_groebner']];assert all(not reduce_bits(f) for f in QG)
 qleads=[max(f,key=qkey) for f in QG];zero=(0,)*10;found={zero};todo=[zero]
 while todo:
  e=todo.pop()
  for j in range(10):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   if u not in found and not any(all(a<=b for a,b in zip(l,u)) for l in qleads):found.add(u);todo.append(u)
 assert len(found)==34 and found==set(map(tuple,tower['module_exponents']))
 # This also verifies that the recorded quotient GB spans the same ideal:
 # each derived GB element reduces by QG, and its S-pairs vanish because
 # the independently derived GB is a GB and QG has the same leading ideal.
 own_min={l for l in bl if not any(k!=l and all(a<=b for a,b in zip(k,l)) for k in bl)}
 q_min={l for l in qleads if not any(k!=l and all(a<=b for a,b in zip(k,l)) for k in qleads)};assert own_min==q_min
 report=dict(result='passed',engine='No Sage or project imports: native sparse rational/GF(2) arithmetic, own Buchberger and Hilbert numerator; SymPy parses text and checks six-dimensional determinants.',full_source_complete_all_degrees=True,integral_invariant_generators=10,generator_degrees=source['generator_degrees'],rational_generator_divisors=divisors,graph_buchberger_pairs=pairs,graph_buchberger_chain_pairs=chain_pairs,kernel_relations=37,relative_buchberger_pairs=qpairs,relative_module_generators=34,relative_module_degrees=sorted(sum(a*b for a,b in zip(e,ws[3:])) for e in found),hilbert_numerator=univariate,seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path,ROOT/'experiments/e7-d5-source-algebra.json')},raw_sha256={towerpath.relative_to(ROOT).as_posix():sha256(towerpath.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-d4-source-audit.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED independent complete D4 source',report['seconds'],flush=True)

if __name__=='__main__':main()

