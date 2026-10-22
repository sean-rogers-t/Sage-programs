"""Independent full characteristic-polynomial and Borel containment audit.

Use the already independently verified spin elementary coefficients. A
Horner calculation in the quadratic algebra computes S16(w-U), followed
by its conjugate product. This checks every coefficient, its integral
root expansion, and its membership in the constructed root Borel image.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb,gcd
from time import perf_counter
import gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mat,fmpq_mpoly_ctx
ROOT=Path(__file__).resolve().parents[1]
IW=(1,2,4,5,6,8)

def main():
 start=perf_counter();cp=ROOT/'experiments/e7-d5a1-orbit-chern.json';cr=json.loads(cp.read_text());ca=ROOT/cr['certificate_archive'];assert sha256(ca.read_bytes()).hexdigest()==cr['certificate_archive_sha256'];cc=json.loads(gzip.decompress(ca.read_bytes()));np=ROOT/'experiments/e7-d5a1-norm.json';nr=json.loads(np.read_text());na=ROOT/nr['certificate_archive'];nc=json.loads(gzip.decompress(na.read_bytes()));nap=ROOT/'experiments/e7-d5a1-norm-audit.json';assert json.loads(nap.read_text())['result']=='passed';bp=ROOT/'experiments/e7-d5a1-images.json';br=json.loads(bp.read_text());ba=ROOT/br['certificate_archive'];bc=json.loads(gzip.decompress(ba.read_bytes()));sp=ROOT/'experiments/e7-d5a1-sources-22.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];sc=json.loads(gzip.decompress(sa.read_bytes()))
 def parse(text,ctx):
  names=dict(zip(ctx.names(),range(len(ctx.names()))));out={}
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(names)
   for f in term.lstrip('-').split('*'):
    name,_,k=f.partition('^')
    if name in names:e[names[name]]+=int(k) if k else 1
    else:c*=fmpq(f)
   e=tuple(e);out[e]=out.get(e,fmpq(0))+c
  return ctx.from_dict(out)
 # U is a genuine polynomial indeterminate; Horner avoids the primary
 # binomial-coefficient expansion of each separate remainder coefficient.
 P=fmpq_mpoly_ctx.get(('h','t','v2','v4','p5','v6','v8','U'));h,t,E1,E2,p5,E3,E4,U=P.gens();es=[parse(f,P) for f in nc['spin_elementary']];a=P.constant(0);b=P.constant(1)
 for k in range(1,17):a,b=b+(h-U)*a,-U*b-t*a+(-1)**k*es[k]
 F=b*b+h*a*b+t*a*a;recorded=P.constant(0)
 for k,text in enumerate(cc['characteristic_coefficients']):recorded+=parse(text,P)*U**k
 assert F==recorded
 R=fmpq_mpoly_ctx.get(('h','z','c2','c3','c4','c5','r'));raw=[parse(f,R)*fmpq(k) for f,k in zip(sc['original_generators'],sc['normalization_factors'])];I=fmpq_mpoly_ctx.get(('h','t','v4','p5','v6','v8'));ih,it,iv4,ip5,iv6,iv8=I.gens();reducers=[]
 def key(e):return (sum(a*b for a,b in zip(e,IW)),tuple(-a for a in reversed(e)))
 for text in bc['groebner_basis']:
  f=parse(text,I).to_dict();e=max(f,key=key);reducers.append((e,{u:-c/f[e] for u,c in f.items() if u!=e}))
 @lru_cache(None)
 def normal_mon(e):
  if key(e)[0]>42:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,b0 in normal_mon(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,fmpq(0))+c*b0
    return {e:c for e,c in out.items() if c}
  return {e:fmpq(1)}
 def image(f):
  out={}
  for e,c in f.to_dict().items():
   assert not e[7]
   for j in range(e[2]+1):
    v=(e[0]+2*(e[2]-j),e[1]+j,e[3],e[4],e[5],e[6]);scale=c*comb(e[2],j)*8**j*(-3)**(e[2]-j)
    for u,b0 in normal_mon(v).items():out[u]=out.get(u,fmpq(0))+scale*b0
  return out
 records=[]
 for row in cc['generators']:
  if row['orbit_shift_by_h']:continue
  d=row['degree'];f=parse(row['classical_polynomial'],P);poly=F.to_dict();coeff=P.from_dict({e[:-1]+(0,):c*(-1)**d for e,c in poly.items() if e[-1]==32-d});assert coeff/fmpq(row['primitive_divisor'])==f
  # Compose in the actual polynomial coefficient ring; denominators are
  # checked before reducing modulo two or applying Borel.
  expanded=f.compose(*raw,R.constant(0));assert all(c.denominator==1 for c in expanded.to_dict().values());assert gcd(*(abs(int(c)) for c in expanded.to_dict().values()))==1;actual={e:int(c)&1 for e,c in expanded.to_dict().items() if int(c)&1};assert actual=={e:int(c) for e,c in parse(row['root_polynomial_mod2'],R).to_dict().items()}
  if d>22:
   im=image(f);norm=list(map(tuple,bc['normal_exponents'][str(d)]));H=fmpq_mat([[fmpq(c) for c in rr] for rr in bc['degrees'][d-1]['image_bases']['M']]);coords=H.inv()*fmpq_mat([[im.get(e,0)] for e in norm]);assert all(int(c.denominator)%2 for c in coords.entries()),(d,'new root image');records.append(dict(degree=d,primitive_integral_source_verified=True,contained_in_constructed_root_image=True))
  print('independent Chern coefficient',d,round(perf_counter()-start,2),flush=True)
 out=dict(result='passed',all_characteristic_coefficients_verified=True,primitive_integral_root_generators_through=32,all_additional_generators_in_constructed_root_image=True,generators=records,seconds=round(perf_counter()-start,3),engine='Native FLINT Horner calculation in a quadratic algebra, independently audited spin moments, actual root coefficient expansions and separate Borel normal forms with exact local membership.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),cp,np,nap,bp,sp)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (ca,na,ba,sa)});(ROOT/'experiments/e7-d5a1-orbit-chern-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED independent source and Borel Chern audit',out['seconds'],flush=True)

if __name__=='__main__':main()
