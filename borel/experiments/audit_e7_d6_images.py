"""Plain Python audit of E7/D6 full-component and sandwich image bases.

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
WS=(1,4,6,6,8,10)
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
 start=perf_counter();primary_path=ROOT/'experiments/e7-d6-images.json'
 primary=json.loads(primary_path.read_text());assert primary['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 variables=('h','v4','v6','p6','v8','v10');ix={v:i for i,v in enumerate(variables)}
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
  # The leading ideal has no standard monomials above degree 33, as
  # verified by exhaustive enumeration of its 126 standard monomials.
  if key(e)[0]>33:return {}
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
 assert len(standard)==126 and max(key(e)[0] for e in standard)==33
 bydegree={d:[tuple(e) for e in cert['normal_exponents'][str(d)]] for d in range(34)}
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
  if degree>33:high_pairs+=1;resolved.add(pair);continue
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
 # Reconstruct every source image from the independently audited complete
 # source component, rather than assuming an algebra generating set.
 sourcepath=ROOT/'experiments/e7-d6-sources.json';source=json.loads(sourcepath.read_text());sourcearchive=ROOT/source['certificate_archive'];assert sha256(sourcearchive.read_bytes()).hexdigest()==source['certificate_archive_sha256'];sources=json.loads(gzip.decompress(sourcearchive.read_bytes()))
 target=[parse(f) for f in ('h','-2*h^2','v4','v6','p6','v8','v10')]
 def multiply(f,g):
  out={}
  for e,a in f.items():
   for u,b in g.items():
    v=tuple(x+y for x,y in zip(e,u));out[v]=out.get(v,mpq(0))+a*b
  return {e:c for e,c in out.items() if c}
 @lru_cache(None)
 def rational_monomial(kind,e):
  if not any(e):return {ZERO:mpq(1)}
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;div=sources[kind]['original_divisors'][j]
  return normal(multiply({v:c/div for v,c in target[j].items()},rational_monomial(kind,tuple(u))))
 root_generators=[]
 for row in sources['M']['degrees']:
  mons=[rational_monomial('M',tuple(e)) for e in row['classical_exponents']]
  for j in row['new_generator_indices']:
   f={}
   for g,line in zip(mons,row['rational_coordinates']):
    c=mpq(line[j])
    for e,a in g.items():f[e]=f.get(e,mpq(0))+c*a
   root_generators.append(({e:c for e,c in f.items() if c},row['degree']))
 assert root_generators==[(parse(f),k) for f,k in cert['root_generators']]
 records=[];history={0:[[mpq(1)]]}
 for entry in cert['degrees']:
  degree=entry['degree'];n=len(bydegree[degree])
  bases={};inverses={}
  for kind in sources:
   if degree<=len(sources[kind]['degrees']):
    row=sources[kind]['degrees'][degree-1];mons=[rational_monomial(kind,tuple(e)) for e in row['classical_exponents']];F=[[f.get(e,mpq(0)) for f in mons] for e in bydegree[degree]];C=[[mpq(c) for c in line] for line in row['rational_coordinates']];A=matmul(F,C)
   else:
    assert kind=='M';columns=[]
    for g,k in root_generators:
     if k<=degree:
      for col in zip(*history[degree-k]):
       f=normal(multiply(g,{e:c for e,c in zip(bydegree[degree-k],col) if c}));columns.append([f.get(e,mpq(0)) for e in bydegree[degree]])
    A=list(map(list,zip(*columns)))
   assert A==[[mpq(c) for c in line] for line in entry['source_image_columns'][kind]]
   H=[[mpq(c) for c in line] for line in entry['image_bases'][kind]];HI=inverse(H);coords=matmul(HI,A);assert local(coords) and rank2(coords)==n
   bases[kind]=H;inverses[kind]=HI
   if kind=='M':history[degree]=H
  transition=matmul(inverses['N'],bases['M']);assert local(transition)
  twice_inverse=[[2*c for c in row] for row in inverse(transition)];assert local(twice_inverse)
  m=n-rank2(transition);assert ({'1':m} if m else {})==entry['smith_multiplicities']
  if degree==8:
   w=cert['first_nonroot_weight_image'];j=w['source_column'];assert w==primary['first_nonroot_weight_image']
   assert w['integral_weight_invariant']==sources['N']['degrees'][7]['basis_polynomials'][j]
   col=[mpq(line[j]) for line in entry['source_image_columns']['N']];f={e:c for e,c in zip(bydegree[8],col) if c};assert parse(w['borel_normal_form'])==f
   v=matmul(inverses['M'],[[c] for c in col]);assert [str(line[0]) for line in v]==w['root_basis_coordinates'];assert not local(v) and local([[2*line[0]] for line in v]);assert w['nonzero_relative_class_of_order']==2
  if degree>len(sources['M']['degrees']):assert m==0
  records.append(dict(degree=degree,rank=n,multiplicity_2=m,both_image_containments_and_primitive_generation=True,root_image_full_by_weight_equality=degree>len(sources['M']['degrees']),twice_weight_image_in_root_image=True))
  if degree%5==0 or degree==33:print('independent complete images degree',degree,round(perf_counter()-start,2),flush=True)
 output=dict(result='passed',engine='Plain Python/gmpy2 sparse weighted reduction and Buchberger checks; independent reconstruction of every full-component image, own rational inverses and bit ranks. No Sage or project imports.',degree_rows=records,max_degree=primary['max_degree'],total_rank=126,dimension=33,buchberger_pairs_checked=pairs,buchberger_chain_pairs=chain_pairs,buchberger_product_pairs=product_pairs,buchberger_above_top_degree_pairs=high_pairs,seconds=round(perf_counter()-start,3),scope='Rational standard basis, every proposed bounded image basis by containment and primitive generation, and every local comparison. Full-source completeness and the fresh-weight ideal identity are checked in separate audits.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path,sourcepath)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (archive,sourcearchive)})
 (ROOT/'experiments/e7-d6-images-audit.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED',output['seconds'],flush=True)
if __name__=='__main__':main()
