"""Complete two-local D4A1 root source via an initial-monoid Hilbert certificate.

Finite integral-product witnesses give a monomial Hilbert lower bound.
An exact all-degree Hilbert identity proves source generation by saturation
and degreewise Nakayama, without an unbounded generator-cutoff assumption.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import argparse,gzip,json,os,sys
import numpy as np
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq_mpoly_ctx
ROOT=Path(__file__).resolve().parents[1]
CW=(1,1,1,2,3,4,2);DW=(1,1,2,2,4,4,6)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def monoid_hilbert(sg,periods,weights,target_weights):
 from heapq import heappush,heappop
 from itertools import product
 from functools import reduce
 residues=list(product(*(range(p) for p in periods)));ideals={r:[] for r in residues};zero=(0,)*7;ideals[zero]=[zero];witness={(zero,zero):(0,)*len(sg)};queue=[(0,zero,zero)];ww=tuple(a*b for a,b in zip(periods,weights))
 while queue:
  degree,r,q=heappop(queue)
  if q not in ideals[r]:continue
  for gi,(g,k) in enumerate(sg):
   rr=tuple((a+b)%p for a,b,p in zip(r,g,periods));qq=tuple(a+(b+c)//p for a,b,c,p in zip(q,r,g,periods));old=ideals[rr]
   if any(all(a<=b for a,b in zip(v,qq)) for v in old):continue
   ideals[rr]=[v for v in old if not all(a<=b for a,b in zip(qq,v))]+[qq];cc=list(witness[(r,q)]);cc[gi]+=1;witness[(rr,qq)]=tuple(cc);heappush(queue,(sum(a*b for a,b in zip(rr,weights))+sum(a*b for a,b in zip(qq,ww)),rr,qq))
 def add(a,b,shift=0,scale=1):
  out=dict(a)
  for d,c in b.items():out[d+shift]=out.get(d+shift,0)+scale*c
  return {d:c for d,c in out.items() if c}
 def product_poly(ds):
  out={0:1}
  for d in ds:out=add(out,out,d,-1)
  return out
 numerator={};rows=[]
 for r in residues:
  if not ideals[r]:continue
  # Exact lcm-lattice inclusion/exclusion for a monomial ideal.
  terms={zero:1}
  for g in ideals[r]:
   nxt=dict(terms)
   for e,c in terms.items():
    u=tuple(max(a,b) for a,b in zip(e,g));nxt[u]=nxt.get(u,0)-c
   terms={e:c for e,c in nxt.items() if c}
  hq={}
  for e,c in terms.items():
   d=sum(a*b for a,b in zip(e,ww));hq[d]=hq.get(d,0)+c
  idealnum=add({0:1},hq,scale=-1);numerator=add(numerator,idealnum,sum(a*b for a,b in zip(r,weights)));rows.append(dict(residue=list(r),minimal_exponents=[list(q) for q in ideals[r]],monoid_witnesses=[list(witness[(r,q)]) for q in ideals[r]]))
 lhs={}
 for d,c in numerator.items():lhs=add(lhs,product_poly(target_weights),d,c)
 rhs=product_poly(ww);return dict(periods=list(periods),polynomial_parameter_weights=list(ww),residues=len(residues),module_minima=sum(map(len,ideals.values())),module_ideals=rows,hilbert_numerator={str(d):c for d,c in sorted(numerator.items())},target_hilbert_series_equal=lhs==rhs)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=32);args=ap.parse_args();start=perf_counter();sp=ROOT/'experiments/e7-d4a1-sources.json';sr=json.loads(sp.read_text());P=fmpq_mpoly_ctx.get(('h','y','b','c2','c3','c4','r'));names=dict(zip(P.names(),range(7)))
 def parse(text):
  out={}
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   c=-1 if term.startswith('-') else 1;e=[0]*7
   for f in term.lstrip('-').split('*'):
    n,_,k=f.partition('^')
    if n in names:e[names[n]]+=int(k) if k else 1
    else:c*=int(f)
   if c&1:out[tuple(e)]=out.get(tuple(e),0)^1
  return frozenset(e for e,c in out.items() if c)
 gens=[dict(degree=d,terms=parse(f),witness=dict(kind='bounded_root_source_16',index=i)) for i,(f,d) in enumerate(sr['generators'])]
 supports={d:list(powers(CW,d)) for d in range(args.max_degree+1)};index={d:{e:i for i,e in enumerate(ee)} for d,ee in supports.items()};basis={0:[1]};provenance={0:[dict(unit=True)]};records=[];extra_witnesses=[];complete=True;semigroup={0:{(0,)*7}};sg=[];sgpoly=[];source_complete=False;monoid_certificate=None
 for d in range(1,args.max_degree+1):
  rows=len(supports[d]);hzero=sum(not e[0] for e in supports[d]);target=sum(1 for _ in powers(DW,d));piv={};selected=[];source=[];tested=0;maps={};arrays={}
  def insert(bits,witness):
   nonlocal tested
   tested+=1;b=bits
   while b:
    k=(b&-b).bit_length()-1
    if k in piv:b^=piv[k]
    else:piv[k]=b;selected.append(bits);source.append(witness);return True
   return False
  def multiply(g,prev,k):
   if len(g)==1 and next(iter(g))==(1,0,0,0,0,0,0):return prev<<hzero
   if k not in arrays:
    size=len(supports[d-k]);arrays[k]={};arrays[k]['size']=size
   size=arrays[k]['size'];bits=np.unpackbits(np.frombuffer(prev.to_bytes((size+7)//8,'little'),dtype=np.uint8),bitorder='little')[:size];out=np.zeros(rows,dtype=np.uint8)
   for e in g:
    if e not in maps:
     maps[e]=np.fromiter((index[d][tuple(a+b for a,b in zip(u,e))] for u in supports[d-k]),dtype=np.int32,count=size)
    out[maps[e]]^=bits
   return int.from_bytes(np.packbits(out,bitorder='little').tobytes(),'little')
  for gi,g in sorted(enumerate(gens),key=lambda it:(it[1]['degree'],it[0])):
   k=g['degree']
   if k>d:continue
   for j,bb in enumerate(basis[d-k]):
    insert(multiply(g['terms'],bb,k),dict(generator=gi,previous_degree=d-k,previous_column=j))
    if len(selected)==target:break
   if len(selected)==target:break
  rank=len(selected);full=rank==target;record=dict(degree=d,rational_source_rank=target,generated_mod_two_rank=rank,full_source_component=full,products_tested=tested,seconds=round(perf_counter()-start,3));records.append(record);print('integral product source',record,flush=True)
  if not full:complete=False;break
  basis[d]=selected;provenance[d]=source
  reachable=set()
  for ee,kk in sg:
   if kk<=d:reachable.update(tuple(a+b for a,b in zip(ee,u)) for u in semigroup[d-kk])
  pivot_exponents={supports[d][i] for i in piv};assert reachable<=pivot_exponents
  for i in sorted(piv):
   ee=supports[d][i]
   if ee not in reachable:sg.append((ee,d));sgpoly.append(dict(degree=d,leading_exponent=list(ee),polynomial_parity_hex=hex(piv[i])));reachable.add(ee)
  semigroup[d]=reachable;assert len(reachable)==target;print('initial semigroup generators',d,len(sg),flush=True)
  periods=[]
  for j in range(7):
   pp=[ee[j] for ee,kk in sg if ee[j] and sum(bool(a) for a in ee)==1];periods.append(min(pp) if pp else None)
  if all(periods):
   monoid_certificate=monoid_hilbert(sg,periods,CW,DW);source_complete=monoid_certificate['target_hilbert_series_equal'];print('monomial Hilbert lower bound',d,source_complete,monoid_certificate['residues'],monoid_certificate['module_minima'],flush=True)
   if source_complete:break
  del maps,arrays,piv
 archive=ROOT/'experiments/e7-d4a1-initial-source-certificates.json.gz';certificate=dict(generators=[dict(degree=g['degree'],coefficient_exponents=[list(e) for e in sorted(g['terms'])],witness=g['witness']) for g in gens],selected_products={str(d):rows for d,rows in provenance.items()},component_parity_sha256={str(d):sha256(b''.join(bb.to_bytes((len(supports[d])+7)//8,'little') for bb in cols)).hexdigest() for d,cols in basis.items()},final_parity_columns_hex=[hex(bb) for bb in basis[max(basis)]],initial_semigroup_generators=sgpoly,monoid_hilbert_certificate=monoid_certificate);archive.write_bytes(gzip.compress(json.dumps(certificate,separators=(',',':')).encode(),mtime=0));out=dict(result='passed' if source_complete else 'incomplete',full_source_complete_all_degrees=source_complete,all_required_source_components_complete=complete,max_complete_degree=max(basis),first_missing_degree=None if complete else records[-1]['degree'],degrees=records,extra_source_generators=extra_witnesses,seconds=round(perf_counter()-start,3),scope=('Actual integral invariant product witnesses give an initial-monoid Hilbert lower bound equal to the full rational invariant Hilbert series. Saturation and degreewise Nakayama prove integral source generation in all degrees.' if source_complete else 'Finite integral-product ranks verify only the recorded full source components; the all-degree Hilbert identity has not been established.'),certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),sp)},raw_sha256={});(ROOT/'experiments/e7-d4a1-initial-source.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
