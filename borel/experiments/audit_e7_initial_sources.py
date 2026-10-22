"""Independent full-source proof via a finite initial-monoid certificate.

Verify actual invariant products, leading-polynomial membership, a monoid
module by integer witnesses and closure, and its entire Hilbert series by
an independent recursive monomial-ideal calculation. No degree cutoff is
used as an unbounded-generation hypothesis.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import argparse,gzip,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--levi',choices=['d4a1','d5a1'],required=True);args=ap.parse_args();kind=args.levi;start=perf_counter();CW=(1,1,1,2,3,4,2) if kind=='d4a1' else (1,1,2,3,4,5,2);DW=(1,1,2,2,4,4,6) if kind=='d4a1' else (1,2,2,4,5,6,8);names=('h','y','b','c2','c3','c4','r') if kind=='d4a1' else ('h','z','c2','c3','c4','c5','r');rp=ROOT/('experiments/e7-'+kind+'-initial-source.json');r=json.loads(rp.read_text());assert r['result']=='passed' and r['full_source_complete_all_degrees'];archive=ROOT/r['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==r['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()));sp=ROOT/('experiments/e7-'+kind+'-sources'+('-22' if kind=='d5a1' else '')+'.json');sr=json.loads(sp.read_text());dependencies=[rp,sp];raws=[archive]
 def parse(text):
  ans=set()
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   c=-1 if term.startswith('-') else 1;e=[0]*7
   for factor in term.lstrip('-').split('*'):
    name,_,power=factor.partition('^')
    if name in names:e[names.index(name)]+=int(power or 1)
    else:c*=int(factor)
   if c%2:
    ee=tuple(e)
    if ee in ans:ans.remove(ee)
    else:ans.add(ee)
  return sorted(ans)
 expected=[(d,parse(f)) for f,d in sr['generators']]
 if kind=='d5a1':
  cp=ROOT/'experiments/e7-d5a1-norm.json';cr=json.loads(cp.read_text());ca=ROOT/cr['certificate_archive'];cc=json.loads(gzip.decompress(ca.read_bytes()));cap=ROOT/'experiments/e7-d5a1-norm-audit.json';sap=ROOT/'experiments/e7-d5a1-sources-22-audit.json';xp=ROOT/'experiments/e7-d5a1-extra23-audit.json';xr=json.loads(xp.read_text());xa=ROOT/xr['certificate_archive'];extra=json.loads(gzip.decompress(xa.read_bytes()));assert all(json.loads(p.read_text())['result']=='passed' for p in (cap,sap,xp));expected +=[(32,parse(cc['norm_mod2'])),(23,parse(extra['source_polynomial']))];dependencies +=[cp,cap,sap,xp];raws +=[ca,xa]
 else:
  sap=ROOT/'experiments/e7-d4a1-generators-audit.json';assert json.loads(sap.read_text())['result']=='passed';dependencies.append(sap)
 assert expected==[(g['degree'],list(map(tuple,g['coefficient_exponents']))) for g in cert['generators']]
 top=r['max_complete_degree'];ranks=[0]*(top+1);ranks[0]=1
 for weight in DW:
  for d in range(weight,top+1):ranks[d]+=ranks[d-weight]
 def mons(d,j=0):
  if j==7:
   if not d:yield ()
   return
  for k in range(d//CW[j],-1,-1):
   for tail in mons(d-k*CW[j],j+1):yield (k,)+tail
 supports={d:list(mons(d)) for d in range(top+1)};indices={d:{e:i for i,e in enumerate(es)} for d,es in supports.items()};basis={0:[1]};table=bytes(int(format(i,'08b')[::-1],2) for i in range(256));verified_polys=[];necessary_source_corrections=[]
 def reverse(bits,size):
  n=(size+7)//8;return int.from_bytes(bits.to_bytes(n,'little').translate(table),'big')>>(8*n-size)
 for d in range(1,top+1):
  m=len(supports[d]);cols=[];piv={};maps={};digest=sha256()
  def productbits(gg,old,bd):
   if gg==[(1,0,0,0,0,0,0)]:value=old
   else:
    size=len(supports[bd]);bits=np.unpackbits(np.frombuffer(old.to_bytes((size+7)//8,'little'),dtype=np.uint8),bitorder='little')[:size];out=np.zeros(m,dtype=np.uint8)
    for shift in gg:
     if shift not in maps:maps[shift]=np.asarray([indices[d][tuple(e[i]+shift[i] for i in range(7))] for e in supports[bd]],dtype=np.int32)
     np.bitwise_xor.at(out,maps[shift],bits)
    value=int.from_bytes(np.packbits(out,bitorder='little').tobytes(),'little')
   return value
  for row in cert['selected_products'][str(d)]:
   k,gg=expected[row['generator']];bd=d-k;old=basis[bd][row['previous_column']];assert row['previous_degree']==bd;value=productbits(gg,old,bd)
   cols.append(value);vv=value
   while vv:
    pivot=(vv&-vv).bit_length()-1
    if pivot in piv:vv^=piv[pivot]
    else:piv[pivot]=vv;break
   digest.update(reverse(value,m).to_bytes((m+7)//8,'little'))
  assert len(piv)==len(cols)==ranks[d];assert digest.hexdigest()==cert['component_parity_sha256'][str(d)];basis[d]=cols
  for row in cert['initial_semigroup_generators']:
   if row['degree']!=d:continue
   vv=reverse(int(row['polynomial_parity_hex'],16),m);lead=supports[d][vv.bit_length()-1];assert list(lead)==row['leading_exponent'];original=vv
   while vv:
    k=(vv&-vv).bit_length()-1;assert k in piv;vv^=piv[k]
   verified_polys.append((tuple(lead),d))
  if kind=='d5a1' and d==23:
   oldpiv={}
   for k,gg in expected[:len(sr['generators'])]:
    for old in basis[d-k]:
     vv=productbits(gg,old,d-k)
     while vv:
      kk=(vv&-vv).bit_length()-1
      if kk in oldpiv:vv^=oldpiv[kk]
      else:oldpiv[kk]=vv;break
   assert len(oldpiv)==780 and ranks[d]==781;necessary_source_corrections.append(dict(degree=23,lower_generator_products_rank=780,full_source_rank=781,necessary_new_generators=1))
  print('independent leading-source polynomials',kind,d,len(verified_polys),round(perf_counter()-start,2),flush=True)
 assert verified_polys==[(tuple(row['leading_exponent']),row['degree']) for row in cert['initial_semigroup_generators']]
 mc=cert['monoid_hilbert_certificate'];pp=tuple(mc['periods']);ww=tuple(a*b for a,b in zip(pp,CW));sg=[e for e,d in verified_polys];zero=(0,)*7
 assert all(tuple(p if i==j else 0 for i in range(7)) in sg for j,p in enumerate(pp));ideals={tuple(row['residue']):[tuple(q) for q in row['minimal_exponents']] for row in mc['module_ideals']};assert zero in ideals and zero in ideals[zero]
 for row in mc['module_ideals']:
  rr=tuple(row['residue']);assert all(0<=a<p for a,p in zip(rr,pp));assert len(row['minimal_exponents'])==len(row['monoid_witnesses'])
  for qq,counts in zip(row['minimal_exponents'],row['monoid_witnesses']):
   assert all(k>=0 for k in counts) and len(counts)==len(sg);assert tuple(a+p*b for a,p,b in zip(rr,pp,qq))==tuple(sum(c*g[i] for c,g in zip(counts,sg)) for i in range(7))
   for g in sg:
    nr=tuple((a+b)%p for a,b,p in zip(rr,g,pp));nq=tuple(a+(b+c)//p for a,b,c,p in zip(qq,rr,g,pp));assert nr in ideals and any(all(a<=b for a,b in zip(v,nq)) for v in ideals[nr])
 # Witnesses prove every module minimum belongs to the generated monoid;
 # closure and the unit prove the opposite inclusion. No normality of
 # the monoid is presumed or used.
 def clean(gens):
  ans=[]
  for e in sorted(set(gens),key=lambda e:(sum(e),e)):
   if not any(all(a<=b for a,b in zip(v,e)) for v in ans):ans.append(e)
  return tuple(sorted(ans))
 def add(a,b,shift=0,scale=1):
  out=dict(a)
  for d,c in b.items():out[d+shift]=out.get(d+shift,0)+scale*c
  return {d:c for d,c in out.items() if c}
 @lru_cache(None)
 def quotient_numerator(weights,gens):
  if not gens:return {0:1}
  if not any(gens[0]):return {}
  # Slice by a different variable heuristic from the discovery code's
  # lcm-lattice inclusion/exclusion.
  maxima=[max(e[j] for e in gens) for j in range(len(weights))];j=min(range(len(weights)),key=lambda j:(maxima[j],-sum(e[j]>0 for e in gens)));m=maxima[j];other=weights[:j]+weights[j+1:];out={}
  for k in range(m+1):
   gg=clean([e[:j]+e[j+1:] for e in gens if e[j]<=k]);hh=quotient_numerator(other,gg);out=add(out,hh,k*weights[j])
   if k<m:out=add(out,hh,(k+1)*weights[j],-1)
  return out
 numerator={}
 for rr,gg in ideals.items():numerator=add(numerator,add({0:1},quotient_numerator(ww,clean(gg)),scale=-1),sum(a*b for a,b in zip(rr,CW)))
 def denominator(ds):
  out={0:1}
  for d in ds:out=add(out,out,d,-1)
  return out
 lhs={}
 for d,c in numerator.items():lhs=add(lhs,denominator(DW),d,c)
 assert lhs==denominator(ww);assert numerator=={int(d):c for d,c in mc['hilbert_numerator'].items()};out=dict(result='passed',full_source_complete_all_degrees=True,levi=kind,generator_degrees=[d for d,g in expected],finite_polynomial_witness_degree=top,initial_monoid_generators=len(sg),periods=list(pp),polynomial_parameter_weights=list(ww),nonempty_residue_modules=len(ideals),module_minima=sum(map(len,ideals.values())),necessary_source_corrections=necessary_source_corrections,all_monoid_module_witnesses_and_closure_verified=True,entire_hilbert_series_verified=True,seconds=round(perf_counter()-start,3),engine='Native reverse-order integral-product reconstruction and opposite Gaussian pivots; exact nonnegative monoid witnesses and module closure; independent recursive monomial-ideal Hilbert numerators. No Sage or primary project imports.',scope='An all-degree Hilbert lower bound equals the rational invariant Hilbert upper bound. Saturation of the full invariant lattice and degreewise Nakayama prove unbounded integral source generation. Finite polynomial witnesses are not a conjectural generator cutoff.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+dependencies},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in raws});(ROOT/('experiments/e7-'+kind+'-initial-source-audit.json')).write_text(json.dumps(out,indent=2)+'\n');print('PASSED complete all-degree source proof',kind,out['seconds'],flush=True)

if __name__=='__main__':main()
