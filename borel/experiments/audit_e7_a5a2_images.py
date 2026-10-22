"""Plain Python audit of E7/(A5+A2) full-source image bases.

No Sage or project imports. Uses a sparse rational weighted reducer,
Buchberger checks, Gaussian inverses, and mod-2 bit ranks. The separate
fresh-weight moment audit checks the rational ideal via a different input.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from functools import lru_cache
from time import perf_counter
import gzip,json
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
WS=(2,3,3,4,5,6)
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
 start=perf_counter();primary_path=ROOT/'experiments/e7-a5a2-images.json'
 primary=json.loads(primary_path.read_text());assert primary['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 variables=('q2','q3','e3','e4','e5','e6');ix={v:i for i,v in enumerate(variables)}
 def parse(text):
  out={}
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   coefficient=mpq(-1 if term.startswith('-') else 1);term=term.lstrip('-');exponents=[0]*6
   for factor in term.split('*'):
    name,_,power=factor.partition('^')
    if name in ix:exponents[ix[name]]+=int(power) if power else 1
    else:coefficient*=mpq(factor)
   e=tuple(exponents);out[e]=out.get(e,mpq(0))+coefficient
  return {e:c for e,c in out.items() if c}
 def key(e):return (sum(w*j for w,j in zip(WS,e)),tuple(-j for j in reversed(e)))
 G=[parse(f) for f in cert['groebner_basis']];leads=[max(f,key=key) for f in G]
 assert all(len({key(e)[0] for e in f})==1 for f in G)
 reducers=[(e,{u:-c/f[e] for u,c in f.items() if u!=e}) for e,f in zip(leads,G)]
 @lru_cache(None)
 def normal_monomial(e):
  # The leading ideal has no standard monomials above degree 45, as
  # verified by exhaustive enumeration of its 672 standard monomials.
  if key(e)[0]>45:return {}
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
 assert len(standard)==672 and max(key(e)[0] for e in standard)==45
 bydegree={d:[tuple(e) for e in cert['normal_exponents'][str(d)]] for d in range(46)}
 assert set().union(*(set(es) for es in bydegree.values()))==standard
 assert all(key(e)[0]==d for d,es in bydegree.items() for e in es)
 # Buchberger's criterion, skipping relatively prime leading monomials.
 pairs=0;chain_pairs=0;product_pairs=0;high_pairs=0;resolved=set();queue=[]
 for i in range(len(G)):
  for j in range(i):
   a,b=leads[i],leads[j];lcm=tuple(max(x,y) for x,y in zip(a,b));queue.append((key(lcm)[0],i,j,lcm))
 for degree,i,j,lcm in sorted(queue):
  a,b=leads[i],leads[j];pair=(i,j)
  if not any(x and y for x,y in zip(a,b)):product_pairs+=1;resolved.add(pair);continue
  if degree>45:high_pairs+=1;resolved.add(pair);continue
  # Buchberger chain criterion: both smaller-lcm pairs already have
  # certified standard representations, hence so does this S-pair.
  witness=next((k for k,c in enumerate(leads) if k not in (i,j)
     and all(x<=y for x,y in zip(c,lcm))
     and tuple(max(x,y) for x,y in zip(a,c))!=lcm
     and tuple(max(x,y) for x,y in zip(b,c))!=lcm
     and (max(i,k),min(i,k)) in resolved and (max(j,k),min(j,k)) in resolved),None)
  if witness is not None:chain_pairs+=1;resolved.add(pair);continue
  shifts=[tuple(x-y for x,y in zip(lcm,lead)) for lead in (a,b)];poly={}
  for f,lead,shift,sign in ((G[i],a,shifts[0],1),(G[j],b,shifts[1],-1)):
   for e,c in f.items():
    u=tuple(x+y for x,y in zip(e,shift));poly[u]=poly.get(u,mpq(0))+sign*c/f[lead]
  assert not normal(poly),(i,j)
  pairs+=1;resolved.add(pair)
 print('independent Buchberger/standard basis',pairs,'pairs',round(perf_counter()-start,2),flush=True)
 assert all(not normal(parse(f)) for f in cert['relations'])
 # Reconstruct the complete root algebra from the independently checked
 # coefficient identities, rather than trusting source-image strings.
 sourcepath=ROOT/'experiments/a5-two-local-source.json';source=json.loads(sourcepath.read_text());sourcearchive=ROOT/source['certificate_archive'];assert sha256(sourcearchive.read_bytes()).hexdigest()==source['certificate_archive_sha256'];sources=json.loads(gzip.decompress(sourcearchive.read_bytes()))
 def multiply(f,g):
  out={}
  for e,c in f.items():
   for u,b in g.items():
    v=tuple(a+b for a,b in zip(e,u));out[v]=out.get(v,mpq(0))+c*b
  return {e:c for e,c in out.items() if c}
 target=[parse(f) for f in ('-q2','e3','e4','e5','e6')];target=[{e:c/mpq(k) for e,c in f.items()} for f,k in zip(target,sources['normalization_factors'])]
 @lru_cache(None)
 def mon(e):
  if not any(e):return {ZERO:mpq(1)}
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return normal(multiply(target[j],mon(tuple(u))))
 generators=dict(N=[(parse(f),d) for f,d in zip(variables,WS)],M=[(parse('q2'),2),(parse('q3'),3)])
 for row in sources['degrees']:
  mons=[mon(tuple(e)) for e in row['classical_exponents']]
  for j in row['new_generator_indices']:
   out={}
   for g,line in zip(mons,row['rational_coordinates']):
    c=mpq(line[j])
    for e,b in g.items():out[e]=out.get(e,mpq(0))+c*b
   generators['M'].append(({e:c for e,c in out.items() if c},row['degree']))
 assert all(gs==[(parse(f),d) for f,d in cert['generators'][kind]] for kind,gs in generators.items());history={k:{0:[[mpq(1)]]} for k in generators};records=[]
 for row in cert['degrees']:
  d=row['degree'];n=len(bydegree[d]);bases={};inverses={}
  for kind,gs in generators.items():
   cols=[]
   for g,k in gs:
    if k<=d:
     for col in zip(*history[kind][d-k]):
      f=normal(multiply(g,{e:c for e,c in zip(bydegree[d-k],col) if c}));cols.append([f.get(e,mpq(0)) for e in bydegree[d]])
   A=list(map(list,zip(*cols)));assert A==[[mpq(c) for c in line] for line in row['image_columns'][kind]],(kind,d,'generator products');H=[[mpq(c) for c in line] for line in row['image_bases'][kind]];HI=inverse(H);C=matmul(HI,A);assert local(C) and rank2(C)==n,(kind,d,'primitive full image');history[kind][d]=H;bases[kind]=H;inverses[kind]=HI
  T=matmul(inverses['N'],bases['M']);assert local(T);assert local([[2*c for c in line] for line in inverse(T)]),(d,'exponent two');m=n-rank2(T);assert row['smith_multiplicities']==({'1':m} if m else {});records.append(dict(degree=d,rank=n,multiplicity_2=m,full_source_image_verified=True,twice_weight_in_root=True))
  if d%5==0 or d==45:print('independent A5+A2 images',d,round(perf_counter()-start,2),flush=True)
 ap=ROOT/'experiments/a5-two-local-audit.json';assert json.loads(ap.read_text())['full_source_complete_all_degrees'];out=dict(result='passed',complete=True,subsystem='A5+A2',top_polynomial_degree=45,total_rank=672,degrees=records,buchberger_pairs=pairs,buchberger_chain_pairs=chain_pairs,buchberger_product_pairs=product_pairs,buchberger_above_top_degree_pairs=high_pairs,seconds=round(perf_counter()-start,3),engine='Native Python/gmpy2: separate sparse rational normal forms, Buchberger criteria, source reconstruction, exact inverses and mod-two primitive ranks. No Sage or primary project imports.',scope='Complete embedded full-source Borel images for the odd-index maximal-rank subgroup. Lattice, fresh-weight invariants and unbounded source completeness are audited separately.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path,sourcepath,ap)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (archive,sourcearchive)});(ROOT/'experiments/e7-a5a2-images-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED A5+A2 images',out['seconds'],flush=True)

if __name__=='__main__':main()
