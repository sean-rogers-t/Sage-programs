"""Native A5 source audit: fresh rational identities and own GF(2) algebra."""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb
from time import perf_counter
import gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mpoly_ctx
from gmpy2 import mpq
sys.setrecursionlimit(20000)
ROOT=Path(__file__).resolve().parents[1]

def parse(text,names):
 ix={v:i for i,v in enumerate(names)};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=mpq(-1 if term.startswith('-') else 1);e=[0]*len(names)
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in ix:e[ix[name]]+=int(power) if power else 1
   else:c*=mpq(f)
  e=tuple(e);out[e]=out.get(e,mpq(0))+c
 return {e:c for e,c in out.items() if c}

def bitmul(a,b):
 out=set()
 for e in a:
  for f in b:
   u=tuple(x+y for x,y in zip(e,f))
   if u in out:out.remove(u)
   else:out.add(u)
 return frozenset(out)

def bitpower(a,k,n):
 out=frozenset([(0,)*n])
 while k:
  if k&1:out=bitmul(out,a)
  a=frozenset(tuple(2*x for x in e) for e in a);k//=2
 return out

def main():
 start=perf_counter();sp=ROOT/'experiments/a5-two-local-source.json';source=json.loads(sp.read_text());sa=ROOT/source['certificate_archive'];assert sha256(sa.read_bytes()).hexdigest()==source['certificate_archive_sha256'];sc=json.loads(gzip.decompress(sa.read_bytes()));ap=ROOT/'experiments/a5-two-local-algebra.json';algebra=json.loads(ap.read_text());aa=ROOT/algebra['certificate_archive'];assert sha256(aa.read_bytes()).hexdigest()==algebra['certificate_archive_sha256'];ac=json.loads(gzip.decompress(aa.read_bytes()));assert algebra['full_source_complete_all_degrees']
 names=('c1','c2','c3','c4','c5');P=fmpq_mpoly_ctx.get(names);cs=[P.constant(1)]+list(P.gens())+[P.constant(0)];mu=cs[1]/6;raw=[sum((comb(6-j,k-j)*(-mu)**(k-j)*cs[j] for j in range(k+1)),P.constant(0)) for k in range(2,7)];factors=list(map(fmpq,sc['normalization_factors']));orig=[f/k for f,k in zip(raw,factors)];assert orig==[P.from_dict({e:fmpq(str(c)) for e,c in parse(f,names).items()}) for f in sc['original_generators']]
 @lru_cache(None)
 def mon(e):
  if not any(e):return P.constant(1)
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return orig[j]*mon(tuple(u))
 generated=[]
 for row in sc['degrees']:
  for j in row['new_generator_indices']:
   actual=sum((mon(tuple(e))*fmpq(line[j]) for e,line in zip(row['classical_exponents'],row['rational_coordinates'])),P.constant(0));recorded,d=source['generators'][len(generated)];assert d==row['degree'];assert actual==P.from_dict({e:fmpq(str(c)) for e,c in parse(recorded,names).items()});assert all(c.denominator==1 for c in actual.to_dict().values());generated.append(actual)
 assert len(generated)==9
 # d_i=x_i-x6 are an exact Q(A5) basis. The first five permutations
 # permute this basis, so symmetric coefficients are the full test ring.
 D=[[int(i==j)+1 for j in range(5)] for i in range(5)];from flint import fmpz_mat
 assert fmpz_mat(D).det()==6
 vs=ac['elimination_variables'];ws=ac['elimination_grading'];n=len(vs)
 def bits(text,nn):return frozenset(e for e,c in parse(text,nn).items() if int(c)&1)
 G=[bits(f,vs) for f in ac['elimination_groebner']];relations=[bits(f,vs) for f in ac['elimination_relations']]
 def key(e):return (sum(a*b for a,b in zip(e[:3],ws[:3])),tuple(-a for a in reversed(e[:3])),sum(a*b for a,b in zip(e[3:],ws[3:])),tuple(-a for a in reversed(e[3:])))
 leads=[max(g,key=key) for g in G]
 dependent=[3+i for i in range(9) if i not in (1,3)];free=[i for i in range(n) if i not in dependent];assert len(free)==5;unit=lambda j:frozenset([tuple(int(i==j) for i in range(5))]);rep={j:unit(k) for k,j in enumerate(free)}
 for g in relations:
  found=[i for i in dependent if tuple(int(i==j) for j in range(n)) in g];assert len(found)==1;i=found[0];tail=g-{tuple(int(i==j) for j in range(n))};assert all(not any(e[j] for j in dependent) for e in tail);rep[i]=frozenset(tuple(e[j] for j in free) for e in tail)
 def evaluate(f,inputs):
  out=frozenset()
  for e in f:
   term=frozenset([(0,)*5])
   for a,k in zip(inputs,e):
    if k:term=bitmul(term,bitpower(a,k,5))
   out=out.symmetric_difference(term)
  return out
 c1,c2,c4=rep[0],rep[1],rep[2];c3=rep[4].symmetric_difference(bitpower(c1,3,5));c5=rep[6].symmetric_difference(bitmul(bitpower(c1,2,5),c3)).symmetric_difference(bitmul(c2,c3)).symmetric_difference(bitmul(c1,c4));inputs=[c1,c2,c3,c4,c5]
 for i,(f,d) in enumerate(source['generators']):assert evaluate(bits(f,names),inputs)==rep[3+i]
 for g in G:assert not evaluate(g,[rep[i] for i in range(n)])
 # Independent exact Hilbert certificates replace expensive S-pair replay.
 # Since every recorded polynomial is in the true graph ideal I, its
 # leading ideal L is contained in in(I). Equality of their known graded
 # Hilbert series proves L=in(I), hence the recorded list is a GB of I.
 def clean(gens):
  out=[]
  for e in sorted(set(gens),key=lambda e:(sum(e),e)):
   if not any(all(a<=b for a,b in zip(f,e)) for f in out):out.append(e)
  return tuple(sorted(out))
 def padd(a,b,shift=0,scale=1):
  out=dict(a)
  for d,c in b.items():out[d+shift]=out.get(d+shift,0)+scale*c
  return {d:c for d,c in out.items() if c}
 @lru_cache(None)
 def hnum(weights,gens):
  if not gens:return {0:1}
  if not any(gens[0]):return {}
  n=len(weights);maxima=[max(e[j] for e in gens) for j in range(n)]
  j=min(range(n),key=lambda j:(maxima[j]/max(1,sum(e[j]>0 for e in gens)),maxima[j]))
  m=maxima[j];ww=weights[:j]+weights[j+1:];out={}
  for k in range(m+1):
   gg=clean([e[:j]+e[j+1:] for e in gens if e[j]<=k]);v=hnum(ww,gg)
   if k<m:out=padd(padd(out,v,k*weights[j]),v,(k+1)*weights[j],-1)
   else:out=padd(out,v,k*weights[j])
  return out
 def numerator(degrees):
  out={0:1}
  for d in degrees:out=padd(out,out,d,-1)
  return out
 graphnum=hnum(tuple(ws),clean(leads));assert graphnum==numerator([d for i,d in enumerate(source['generator_degrees']) if i not in (1,3)])
 elimination=[g for g in G if all(not any(e[:3]) for e in g)];assert len(elimination)==algebra['kernel_relations'];ll=clean([max(g,key=key)[3:] for g in elimination]);univ=hnum(tuple(ws[3:]),ll);assert univ==numerator((8,9,10,12))
 print('A5 exact graph and kernel Hilbert certificates',hnum.cache_info(),round(perf_counter()-start,2),flush=True)
 out=dict(result='passed',full_source_complete_all_degrees=True,generator_degrees=source['generator_degrees'],root_basis_index_in_weights=6,groebner_certified_by_exact_graph_hilbert_series=True,graph_groebner_size=len(G),kernel_relations=len(elimination),hilbert_numerator=univ,seconds=round(perf_counter()-start,3),engine='Native FLINT rational identities and integer lattice determinant; own sparse GF(2) substitution and recursive monomial Hilbert numerators, certifying the graph Groebner basis by exact Hilbert equality. No Sage or primary project imports.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),sp,ap)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (sa,aa)});(ROOT/'experiments/a5-two-local-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED A5 full source',out['seconds'],flush=True)

if __name__=='__main__':main()
