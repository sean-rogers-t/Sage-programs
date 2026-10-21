"""Independent exact full-component audit for E7/P6, using native FLINT.

The invariant rational basis is reconstructed from Newton identities. Its
recorded rational changes of basis must give integral, primitive columns.
No Sage or project modules are imported; no algebra-generation assertion.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb,gcd,lcm
from time import perf_counter
import argparse,gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mpoly_ctx,fmpz_mod_mpoly_ctx,fmpz_mat
ROOT=Path(__file__).resolve().parents[1]
WS=(1,1,2,3,4,5,2)
DS=(1,2,2,4,5,6,8)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def parse(text,ctx):
 indices={str(v):j for j,v in enumerate(ctx.names())};terms={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(indices)
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in indices:e[indices[name]]+=int(power) if power else 1
   else:c*=fmpq(f)
  e=tuple(e);terms[e]=terms.get(e,fmpq(0))+c
 return ctx.from_dict(terms)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--degree-bound',type=int,default=20);args=ap.parse_args();stem='e7-d5a1-sources'+('' if args.degree_bound==20 else '-'+str(args.degree_bound));start=perf_counter();path=ROOT/('experiments/'+stem+'.json');source=json.loads(path.read_text());archive=ROOT/source['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==source['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()));records=[]
 ctx=fmpq_mpoly_ctx.get(('h','z','c2','c3','c4','c5','r'));h,z,c2,c3,c4,c5,r=ctx.gens();one=ctx.constant(1);zero=ctx.constant(0);es=[one,z,c2,c3,c4,c5];ps=[5*one]
 for k in range(1,9):ps.append(sum(((-1)**(j-1)*es[j]*ps[k-j] for j in range(1,min(k,6))),zero)+((-1)**(k-1)*k*es[k] if k<=5 else zero))
 xp={k:sum((comb(k,j)*2**j*(-h)**(k-j)*ps[j] for j in range(k+1)),zero) for k in (2,4,6,8)};sq=[one]
 for k in range(1,5):sq.append(sum(((-1)**(j-1)*sq[k-j]*xp[2*j] for j in range(1,k+1)),zero)/k)
 pf=sum((2**j*(-h)**(5-j)*es[j] for j in range(6)),zero);raw=[h,r+h*z/2-z*z/4,sq[1],sq[2],pf,sq[3],sq[4]]
 factors=[fmpq(gcd(*(int(c.numerator) for c in f.to_dict().values())),lcm(*(int(c.denominator) for c in f.to_dict().values()))) for f in raw];orig=[f/k for f,k in zip(raw,factors)]
 assert factors==list(map(fmpq,cert['normalization_factors']));assert orig==[parse(f,ctx) for f in cert['original_generators']]
 @lru_cache(None)
 def mon(e):
  if not any(e):return one
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return orig[j]*mon(tuple(u))
 assert [row['degree'] for row in cert['degrees']]==list(range(1,source['max_degree']+1))
 generated=[];modctx=fmpz_mod_mpoly_ctx.get(ctx.names(),modulus=2);modbases={0:[modctx.constant(1)]};modgens=[];generation=[]
 for row in cert['degrees']:
  d=row['degree'];exps=list(powers(DS,d));support=list(powers(WS,d));n=len(exps);assert [list(e) for e in exps]==row['classical_exponents'];F=[mon(e).to_dict() for e in exps];assert all(c.denominator==1 for f in F for c in f.values())
  C=[[fmpq(c) for c in line] for line in row['rational_coordinates']];D=lcm(*(int(c.denominator) for line in C for c in line));assert D>0 and D&(D-1)==0;B=fmpz_mat([[int(c*D) for c in line] for line in C]);assert B.nrows()==B.ncols()==n;bits=[0]*n
  proposed={j:parse(f,ctx).to_dict() for j,f in zip(row['new_generator_indices'],row['new_generator_polynomials'])}
  for k in range(0,len(support),128):
   ss=support[k:k+128];left=fmpz_mat([[int(f.get(e,0)) for f in F] for e in ss])*B
   for i,e in enumerate(ss):
    for j in range(n):
     a=int(left[i,j]);assert a%D==0,(d,k+i,j);v=a//D;bits[j]|=(v&1)<<(k+i)
     if j in proposed:assert v==proposed[j].get(e,0),(d,j,e)
  piv={}
  for a in bits:
   while a:
    k=a.bit_length()-1
    if k in piv:a^=piv[k]
    else:piv[k]=a;break
  assert len(piv)==n
  modbases[d]=[modctx.from_dict({support[i]:1 for i in range(len(support)) if (a>>i)&1}) for a in bits]
  if d>20:
   ix={e:i for i,e in enumerate(support)};piv={}
   def insert(a):
    while a:
     k=a.bit_length()-1
     if k in piv:a^=piv[k]
     else:piv[k]=a;return True
    return False
   for g,k in modgens:
    for b in modbases[d-k]:insert(sum(1<<ix[e] for e,c in (g*b).to_dict().items()))
   oldrank=len(piv);newrank=0
   for j in row['new_generator_indices']:newrank+=int(insert(bits[j]))
   assert len(piv)==n and newrank==len(row['new_generator_indices'])
   generation.append(dict(degree=d,full_component_rank=n,products_of_lower_generators_rank=oldrank,necessary_new_generators=newrank));print('independent indecomposables',generation[-1],flush=True)
  for j in row['new_generator_indices']:modgens.append((modbases[d][j],d))
  generated.extend([[f,d] for f in row['new_generator_polynomials']]);records.append(dict(degree=d,rank=n,integral_invariant_columns=True,primitive_mod2=True));print('independent full root source',d,n,round(perf_counter()-start,2),flush=True)
 assert generated==source['generators']
 out=dict(result='passed',engine='Native Python/python-flint: fresh Newton identities, exact integer matrix products, independent mod-two ranks and products. No Sage or project imports.',max_degree=source['max_degree'],degrees=records,indecomposable_checks=generation,seconds=round(perf_counter()-start,3),scope='Full localized root invariant components through the stated bound only. The coefficient ring and Weyl action are checked in the coordinate/moment audit. Additional mod-two product ranks prove the displayed generator truncations fail; no unbounded generation claim.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()});(ROOT/('experiments/'+stem+'-audit.json')).write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],flush=True)

if __name__=='__main__':main()
