"""Independent fixed-transcript audit of the extra degree-23 source generator.

Start from the coefficient-complete degree-22 seed. Multiplication by h
preserves the echelon ordering, so only the h-free rational complement
needs new saturation at each degree. Track the Borel images alongside
integer FLINT polynomials; no dense source basis-change matrix is built.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb
from time import perf_counter
import argparse,gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mpoly_ctx,fmpz_mpoly_ctx,fmpz_mat
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
CW=(1,1,2,3,4,5,2);DW=(1,2,2,4,5,6,8);IW=(1,2,4,5,6,8)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def parse(text,ctx):
 ix={str(v):i for i,v in enumerate(ctx.names())};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(ix)
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in ix:e[ix[name]]+=int(power) if power else 1
   else:c*=fmpq(f)
  e=tuple(e);out[e]=out.get(e,fmpq(0))+c
 return ctx.from_dict(out)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=27);args=ap.parse_args();start=perf_counter();sp=ROOT/'experiments/e7-d5a1-sources-22.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];assert sha256(sa.read_bytes()).hexdigest()==sr['certificate_archive_sha256'];sc=json.loads(gzip.decompress(sa.read_bytes()));ip=ROOT/'experiments/e7-d5a1-images.json';ir=json.loads(ip.read_text());ia=ROOT/ir['certificate_archive'];assert sha256(ia.read_bytes()).hexdigest()==ir['certificate_archive_sha256'];ic=json.loads(gzip.decompress(ia.read_bytes()))
 Q=fmpq_mpoly_ctx.get(('h','z','c2','c3','c4','c5','r'));Z=fmpz_mpoly_ctx.get(Q.names());origq=[parse(f,Q) for f in sc['original_generators']];orig=[Z.from_dict({e:int(c) for e,c in f.to_dict().items()}) for f in origq];hh=Z.gen(0);one=Z.constant(1)
 def poly_mon(e):return __import__('functools').reduce(lambda a,b:a*b,(orig[j]**k for j,k in enumerate(e) if k),one)
 ivars=('h','t','v4','p5','v6','v8');imctx=fmpq_mpoly_ctx.get(ivars);ix={n:i for i,n in enumerate(ivars)}
 def iparse(text):return {tuple(e):mpq(str(c)) for e,c in parse(text,imctx).to_dict().items()}
 def key(e):return (sum(w*k for w,k in zip(IW,e)),tuple(-k for k in reversed(e)))
 GG=[iparse(f) for f in ic['groebner_basis']];reducers=[]
 for g in GG:
  e=max(g,key=key);reducers.append((e,{u:-c/g[e] for u,c in g.items() if u!=e}))
 @lru_cache(None)
 def normal_mon(e):
  if key(e)[0]>42:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,b in normal_mon(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,mpq(0))+c*b
    return {e:c for e,c in out.items() if c}
  return {e:mpq(1)}
 def inormal(f):
  out={}
  for e,c in f.items():
   for u,b in normal_mon(e).items():out[u]=out.get(u,mpq(0))+c*b
  return {e:c for e,c in out.items() if c}
 def imul(a,b):
  out={}
  for e,c in a.items():
   for f,v in b.items():
    u=tuple(x+y for x,y in zip(e,f));out[u]=out.get(u,mpq(0))+c*v
  return inormal(out)
 norms=list(map(mpq,sc['normalization_factors']));target=[iparse(f) for f in ('h','t','8*t-3*h^2','v4','p5','v6','v8')];target=[{e:c/k for e,c in f.items()} for f,k in zip(target,norms)]
 @lru_cache(None)
 def image_mon(e):
  if not any(e):return {(0,)*6:mpq(1)}
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return imul(target[j],image_mon(tuple(u)))
 def image_vec(f,d):return [f.get(tuple(e),mpq(0)) for e in ic['normal_exponents'][str(d)]]
 support=list(powers(CW,22));exps=[tuple(e) for e in sc['degrees'][-1]['classical_exponents']];C=sc['degrees'][-1]['rational_coordinates'];den=max(int(mpq(c).denominator) for line in C for c in line);assert den&(den-1)==0;B=fmpz_mat([[int(mpq(c)*den) for c in line] for line in C]);F=[poly_mon(e).to_dict() for e in exps];cols=[{} for _ in exps]
 for k in range(0,len(support),128):
  ss=support[k:k+128];A=fmpz_mat([[int(f.get(e,0)) for f in F] for e in ss])*B
  for i,e in enumerate(ss):
   for j in range(len(exps)):
    v=int(A[i,j]);assert v%den==0
    if v:cols[j][e]=v//den
 seedpolys=[Z.from_dict(f) for f in cols];seedimages=[[mpq(line[j]) for line in ic['degrees'][21]['image_columns']['M']] for j in range(len(exps))];del C,B,F,cols
 # A fixed transcript is checked, rather than discovering an echelon basis.
 existing=ROOT/'experiments/e7-d5a1-extra23-certificates.json.gz'
 if existing.exists():
  packaged=json.loads(gzip.decompress(existing.read_bytes()));transcript=dict(steps=packaged['seed22_and_degree23_operations']);expected_polynomial=packaged['source_polynomial'];expected_image=list(map(mpq,packaged['borel_normal_coordinates']))
 else:
  xp=ROOT/'experiments/e7-d5a1-sparse-23.json';xr=json.loads(xp.read_text());xa=ROOT/xr['certificate_archive'];assert sha256(xa.read_bytes()).hexdigest()==xr['certificate_archive_sha256'];transcript=json.loads(gzip.decompress(xa.read_bytes()));expected_polynomial=transcript['final_source_polynomials'][780];expected_image=[mpq(row[780]) for row in transcript['image_columns']]
 seedstep,step=transcript['steps'];assert seedstep['degree']==22 and step['degree']==23;basis=[];images=[]
 for row in seedstep['columns']:
  f=seedpolys[row['seed_column']];v=seedimages[row['seed_column']][:]
  for j in row['operations']:
   assert j>=0;f-=basis[j];v=[a-b for a,b in zip(v,images[j])]
  basis.append(f);images.append(v)
 previous=basis;oldimages=images;basis=[hh*f for f in previous];images=[image_vec(imul(iparse('h'),dict(zip(map(tuple,ic['normal_exponents']['22']),v))),23) for v in oldimages]
 for row in step['columns']:
  e=tuple(row['classical_exponent']);assert not e[0];f=poly_mon(e);v=image_vec(image_mon(e),23)
  for j in row['operations']:
   if j<0:
    assert j==-1 and all(int(c)%2==0 for c in f.to_dict().values());f=f//2;v=[a/2 for a in v]
   else:f-=basis[j];v=[a-b for a,b in zip(v,images[j])]
  basis.append(f);images.append(v)
 assert len(basis)==781;chosen=780;actual=parse(expected_polynomial,Q);assert actual==Q.from_dict({e:int(c) for e,c in basis[chosen].to_dict().items()});assert images[chosen]==expected_image
 from flint import fmpq_mat
 HM=fmpq_mat([[fmpq(c) for c in row] for row in ic['degrees'][22]['image_bases']['M']]);V=fmpq_mat([[fmpq(str(c))] for c in images[chosen]]);coordinates=HM.inv()*V;assert all(int(c.denominator)%2 for c in coordinates.entries());print('extra degree-23 generator: exact source and Borel transcript passed',round(perf_counter()-start,3),flush=True)
 # Package only the required witness and compact transcript, not all
 # 781 expanded source polynomials from the exploratory checkpoint.
 extra=dict(degree=23,source_polynomial=str(basis[chosen]),borel_normal_coordinates=list(map(str,images[chosen])),root_lower_image_coordinates=[str(c) for c in coordinates.entries()],seed22_and_degree23_operations=transcript['steps']);dest=ROOT/'experiments/e7-d5a1-extra23-certificates.json.gz';dest.write_bytes(gzip.compress(json.dumps(extra,separators=(',',':')).encode(),mtime=0));report=dict(result='passed',degree=23,source_generator_coefficient_integral=True,rational_invariance_preserved_by_exact_operations=True,borel_image_in_constructed_root_image=True,source_terms=len(basis[chosen].to_dict()),seconds=round(perf_counter()-start,3),engine='Independent replay of the fixed source transcript with exact integer polynomials and rational Borel vectors, starting at the independently verified full degree-22 coefficient basis. Each halving is checked coefficientwise. The primary rank-discovery algorithm is not called; native FLINT is shared generic arithmetic.',certificate_archive=dest.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(dest.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),sp,ip,ROOT/'experiments/e7-d5a1-sources-22-audit.json',ROOT/'experiments/e7-d5a1-images-audit.json')},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (sa,ia)});(ROOT/'experiments/e7-d5a1-extra23-audit.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
