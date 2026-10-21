"""Independent full D5 source audit: no Sage and no project imports.

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
 start=perf_counter();srcpath=ROOT/'experiments/e7-d5-source-algebra.json';source=json.loads(srcpath.read_text());assert source['result']=='passed'
 imagepath=ROOT/'experiments/e7-d5-images.json';image=json.loads(imagepath.read_text());archive=ROOT/image['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==image['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 # Exact E6 weight lattice and the A4 permutation compression.
 C=s.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
 om=[s.eye(6).row(i) for i in range(6)];h0=om[0];p=om[2];us=[-om[1]+om[2],om[1]+om[2]-om[3],om[3]-om[4],om[4]-om[5],om[5]]
 assert sum(us,s.zeros(1,6))==2*p
 assert abs(s.Matrix.vstack(h0,p,*us[:4]).det())==1
 assert abs(s.Matrix.vstack(h0,p,*[u-us[4] for u in us[:4]]).det())==5
 for i in (1,3,4,5):
  reflected=[u-u[i]*C.row(i) for u in us];assert set(map(tuple,reflected))==set(map(tuple,us))
 xx=[2*u-h0 for u in us]
 for i in (1,2,3,4,5):
  reflected=[x-x[i]*C.row(i) for x in xx];signs=[];permutation=[]
  for x in reflected:
   match=next((j,sign) for j,y in enumerate(xx) for sign in (-1,1) if x==sign*y);permutation.append(match[0]);signs.append(match[1])
  assert len(set(permutation))==5 and __import__('math').prod(signs)==1
 assert C.det()==3
 variables=s.symbols('h z c2 c3 c4 c5');units=[{tuple(int(i==j) for i in range(6)):mpq(1)} for j in range(6)];h,z,c2,c3,c4,c5=units;one={Z:mpq(1)};es=[one,scale(z,2),c2,c3,c4,c5]
 ps=[scale(one,5)]
 for k in range(1,9):
  f={}
  for j in range(1,min(k,6)):f=add(f,mul(es[j],ps[k-j]),(-1)**(j-1))
  if k<=5:f=add(f,es[k],(-1)**(k-1)*k)
  ps.append(f)
 xp={0:scale(one,5)}
 for k in (2,4,6,8):
  out={}
  for j in range(k+1):out=add(out,mul(power(h,k-j),ps[j]),comb(k,j)*2**j*(-1)**(k-j))
  xp[k]=out
 sq=[one]
 for k in range(1,5):
  out={}
  for j in range(1,k+1):out=add(out,mul(sq[k-j],xp[2*j]),(-1)**(j-1))
  sq.append(scale(out,mpq(1,k)))
 eu={}
 for j in range(6):eu=add(eu,mul(power(h,5-j),es[j]),2**j*(-1)**(5-j))
 original=[h,sq[1],scale(sq[2],mpq(1,2)),eu,scale(sq[3],mpq(1,2)),sq[4]]
 @lru_cache(None)
 def monomial(e):
  if not any(e):return one
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return mul(original[j],monomial(tuple(u)))
 for item,recorded in zip(cert['preimages'],source['integral_generators']):
  actual={}
  for e,c in zip(item['rational_exponents'],item['rational_coefficients']):actual=add(actual,monomial(tuple(e)),mpq(c))
  assert actual==parse(recorded,variables)==parse(item['polynomial'],variables),(item['degree'],add(actual,parse(recorded,variables),-1))
  assert all(c.denominator==1 for c in actual.values())
 # Check the actual source-to-rational-quotient substitution as well.
 quotient_variables=s.symbols('h w v4 v5 v6 v8');unitq=[{tuple(int(i==j) for i in range(6)):mpq(1)} for j in range(6)]
 e1q=add(scale(power(unitq[0],2),-3),scale(power(unitq[1],2),-mpq(8,3)))
 qgens=[unitq[0],e1q,scale(unitq[2],mpq(1,2)),unitq[3],scale(unitq[4],mpq(1,2)),unitq[5]]
 @lru_cache(None)
 def quotient_monomial(e):
  if not any(e):return one
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return mul(qgens[j],quotient_monomial(tuple(u)))
 for i,item in enumerate(cert['preimages']):
  out={}
  for e,c in zip(item['rational_exponents'],item['rational_coefficients']):out=add(out,quotient_monomial(tuple(e)),mpq(c))
  assert out==parse(cert['source_generators']['N'][i+1][0],quotient_variables)
 assert len(cert['preimages'])==10
 print('independent lattice and ten integral invariant identities',round(perf_counter()-start,2),flush=True)
 # Independent GF(2) graph-ideal arithmetic.
 vs=s.symbols(' '.join(source['elimination_variables']));ws=source['elimination_grading'];n=len(vs)
 G=[parse(f,vs,True) for f in source['elimination_groebner']];relations=[parse(f,vs,True) for f in source['elimination_relations']]
 def key(e):return (sum(a*b for a,b in zip(e[:3],ws[:3])),tuple(-x for x in reversed(e[:3])),sum(a*b for a,b in zip(e[3:],ws[3:])),tuple(-x for x in reversed(e[3:])))
 leads=[max(f,key=key) for f in G]
 @lru_cache(None)
 def normal_monomial(e):
  for f,lead in zip(G,leads):
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out=frozenset()
    for u in f:
     if u!=lead:out=out.symmetric_difference(normal_monomial(tuple(a+b for a,b in zip(u,shift))))
    return out
  return frozenset([e])
 def normal(f):
  out=frozenset()
  for e in f:out=out.symmetric_difference(normal_monomial(e))
  return out
 assert all(not normal(f) for f in relations)
 pairs=0
 for i in range(len(G)):
  for j in range(i):
   a,b=leads[i],leads[j]
   if not any(x and y for x,y in zip(a,b)):continue
   lcm=tuple(max(x,y) for x,y in zip(a,b));shifts=[tuple(x-y for x,y in zip(lcm,l)) for l in (a,b)]
   S=frozenset(tuple(x+y for x,y in zip(e,shifts[0])) for e in G[i]).symmetric_difference(tuple(x+y for x,y in zip(e,shifts[1])) for e in G[j])
   assert not normal(S),(i,j);pairs+=1
 # The seven graph equations solve seven dependent variables over the
 # six free variables z,b,c4,h,b1,b3. Substitution checks reverse containment.
 dependent=[5,7,8,9,10,11,12];free=[i for i in range(n) if i not in dependent];assert len(free)==6
 replacements={i:frozenset([tuple(int(j==k) for j in range(6))]) for k,i in enumerate(free)}
 for f in relations:
  candidates=[i for i in dependent if tuple(int(i==j) for j in range(n)) in f];assert len(candidates)==1;i=candidates[0];unit=tuple(int(i==j) for j in range(n));tail=f-{unit}
  assert all(not any(e[j] for j in dependent) for e in tail)
  replacements[i]=frozenset(tuple(e[j] for j in free) for e in tail)
 # Verify that this graph really describes the ten source polynomials,
 # by rebuilding its triangular coordinate substitution independently.
 def evaluate_bits(f,inputs):
  out=frozenset()
  for e in f:
   term=frozenset([Z])
   for a,k in zip(inputs,e):
    if k:term=bitmul(term,bitpower(a,k,6))
   out=out.symmetric_difference(term)
  return out
 hbit=replacements[3];zbit=replacements[0];abit=replacements[4];bbit=replacements[1];c4bit=replacements[2];tbit=replacements[6]
 c2bit=abit.symmetric_difference(bitmul(hbit,zbit));c3bit=bbit.symmetric_difference(bitmul(bitpower(hbit,2,6),zbit)).symmetric_difference(bitmul(hbit,bitpower(zbit,2,6)))
 source_bits=[parse(f,variables,True) for f in source['generators']]
 c5bit=tbit.symmetric_difference(evaluate_bits(source_bits[3],[hbit,zbit,c2bit,c3bit,c4bit,frozenset()]))
 actual_inputs=[hbit,zbit,c2bit,c3bit,c4bit,c5bit]
 for i,f in enumerate(source_bits):assert evaluate_bits(f,actual_inputs)==replacements[3+i]
 @lru_cache(None)
 def substitute_monomial(e):
  out=frozenset([Z])
  for i,k in enumerate(e):
   if k:out=bitmul(out,bitpower(replacements[i],k,6))
  return out
 for f in G:
  out=frozenset()
  for e in f:out=out.symmetric_difference(substitute_monomial(e))
  assert not out
 elimination=[f for f in G if all(not any(e[:3]) for e in f)];assert len(elimination)==20
 # Taylor inclusion-exclusion on the leading monomial ideal gives its
 # Hilbert numerator, independently of Sage/Singular Hilbert-series code.
 numerator={(0,)*10:1}
 for f in elimination:
  lead=max(f,key=key)[3:];new=dict(numerator)
  for e,c in numerator.items():
   u=tuple(max(a,b) for a,b in zip(e,lead));new[u]=new.get(u,0)-c
  numerator={e:c for e,c in new.items() if c}
 univariate={}
 for e,c in numerator.items():
  d=sum(a*b for a,b in zip(e,ws[3:]));univariate[d]=univariate.get(d,0)+c
 univariate={d:c for d,c in univariate.items() if c};expected={0:1}
 for d in (8,9,12,16):
  nxt=dict(expected)
  for k,c in expected.items():nxt[k+d]=nxt.get(k+d,0)-c
  expected={k:c for k,c in nxt.items() if c}
 assert univariate==expected
 report=dict(result='passed',engine='Plain Python sparse rational and GF(2) polynomials; own Buchberger reductions, graph substitution, and monomial-ideal Hilbert numerator. SymPy only parses polynomials and checks six-dimensional lattice determinants.',integral_invariant_generators=10,source_generator_degrees=source['generator_degrees'],elimination_groebner_size=len(G),buchberger_pairs_checked=pairs,elimination_kernel_generators=20,hilbert_numerator=univariate,full_source_complete_all_degrees=True,seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),srcpath,imagepath)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-d5-source-audit.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED independent complete source',report['seconds'],flush=True)
if __name__=='__main__':main()
