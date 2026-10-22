"""Independent native audit of the D5-refinement image upper bound."""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mat,fmpq_mpoly_ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
WS=(1,2,4,5,6,8)

def main():
 start=perf_counter();boundp=ROOT/'experiments/e7-d5a1-refinement-bound.json';br=json.loads(boundp.read_text());ba=ROOT/br['certificate_archive'];assert sha256(ba.read_bytes()).hexdigest()==br['certificate_archive_sha256'];bc=json.loads(gzip.decompress(ba.read_bytes()))
 def read(name):
  p=ROOT/'experiments'/name;r=json.loads(p.read_text());a=ROOT/r['certificate_archive'];assert sha256(a.read_bytes()).hexdigest()==r['certificate_archive_sha256'];return p,r,a,json.loads(gzip.decompress(a.read_bytes()))
 xp,xr,xa,xc=read('e7-d5a1-images.json');fp,fr,fa,fc=read('e7-d5-images.json');P=fmpq_mpoly_ctx.get(('h','t','v4','p5','v6','v8'));h,t,v4,p5,v6,v8=P.gens()
 def parse(text,ctx=P):
  out={};names=dict(zip(ctx.names(),range(ctx.nvars())))
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(names)
   for f in term.lstrip('-').split('*'):
    name,_,power=f.partition('^')
    if name in names:e[names[name]]+=int(power) if power else 1
    else:c*=fmpq(f)
   e=tuple(e);out[e]=out.get(e,fmpq(0))+c
  return ctx.from_dict(out)
 def key(e):return (sum(a*b for a,b in zip(e,WS)),tuple(-a for a in reversed(e)))
 reducers=[]
 for text in xc['groebner_basis']:
  f=parse(text).to_dict();e=max(f,key=key);reducers.append((e,{u:-c/f[e] for u,c in f.items() if u!=e}))
 @lru_cache(None)
 def normal_mon(e):
  if key(e)[0]>42:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,b in normal_mon(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,fmpq(0))+c*b
    return {e:c for e,c in out.items() if c}
  return {e:fmpq(1)}
 def red(f):
  out={}
  for e,c in f.to_dict().items():
   for u,b in normal_mon(e).items():out[u]=out.get(u,fmpq(0))+c*b
  return P.from_dict(out)
 def mul(a,b):return (red(a[0]*b[0]-t*a[1]*b[1]),red(a[0]*b[1]+a[1]*b[0]+h*a[1]*b[1]))
 targets=[(-h,P.constant(fmpq(4,3))),(P.constant(0),P.constant(1)),(v4,P.constant(0)),(-p5,P.constant(0)),(v6,P.constant(0)),(v8,P.constant(0))]
 @lru_cache(None)
 def mon(e):
  if not any(e):return (P.constant(1),P.constant(0))
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return mul(targets[j],mon(tuple(u)))
 O=fmpq_mpoly_ctx.get(('h','w','v4','v5','v6','v8'))
 for text in fc['relations']:
  out=[P.constant(0),P.constant(0)]
  for e,c in parse(text,O).to_dict().items():
   a,b=mon(e);out[0]+=c*a;out[1]+=c*b
  assert all(f.is_zero() for f in out)
 def mat(rows):return fmpq_mat([[fmpq(c) for c in row] for row in rows])
 def local(A):return all(int(c.denominator)%2 for c in A.entries())
 def rankcols(cc):
  piv={}
  for v in cc:
   while v:
    k=v.bit_length()-1
    if k in piv:v^=piv[k]
    else:piv[k]=v;break
  return len(piv)
 norm={d:[tuple(e) for e in xc['normal_exponents'].get(str(d),[])] for d in range(44)};NH={0:fmpq_mat([[1]]),43:fmpq_mat(0,0)}
 for row in xc['degrees']:NH[row['degree']]=mat(row['image_bases']['N'])
 records=[]
 for row,cr in zip(fc['image_degrees'],bc['degrees']):
  d=row['degree'];n=len(norm[d]);prev=len(norm[d-1]);oldmons=[mon(tuple(e)) for e in fc['normal_exponents'][str(d)]];T=fmpq_mat(n+prev,len(oldmons))
  for j,(a,b) in enumerate(oldmons):
   assert str(a)==str(parse(cr['old_normal_images'][j][0])) and str(b)==str(parse(cr['old_normal_images'][j][1]))
   ad=a.to_dict();bd=b.to_dict()
   for i,e in enumerate(norm[d]):T[i,j]=ad.get(e,0)
   for i,e in enumerate(norm[d-1]):T[n+i,j]=bd.get(e,0)
  H=fmpq_mat(n+prev,n+prev)
  for off,A in ((0,NH[d]),(n,NH[d-1])):
   for i in range(A.nrows()):
    for j in range(A.ncols()):H[off+i,off+j]=A[i,j]
  W=H.inv()*T*mat(row['weight_basis']);B=H.inv()*T*mat(row['root_basis']);assert local(W) and local(W.inv()) and local(B) and local(2*B.inv());assert W==mat(cr['weight_transition']) and B==mat(cr['root_transition']);cc=[sum((int(B[i,j].numerator)&1)<<i for i in range(B.nrows())) for j in range(B.ncols())];r=rankcols(cc);upper=rankcols(cc+[1<<i for i in range(n)])-r;assert upper==cr['refinement_intersection_cokernel_multiplicity'];records.append(dict(degree=d,intersection_cokernel_multiplicity=upper,matching_lower_image=cr['full_image_proved']))
 C=S.Matrix([[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]]);E=C[:6,:6];om=[sum((E.inv()[i,j]*C.col(j) for j in range(6)),S.zeros(7,1)) for i in range(6)];w=S.eye(7).col(6);B=S.Matrix.hstack(*om,w);perm=[5,1,4,3,2,0];w0=S.Matrix.hstack(*[-om[i] for i in perm],w)*B.inv();assert all(S.denom(c)==1 for c in C.inv()*w0*C) and (C.inv()*w0*C).det() in (1,-1);assert w0*om[0]==-S.eye(7).col(5)+4*w/3 and w0*w==w
 out=dict(result='passed',degrees=records,all_recorded_upper_bounds_verified=True,full_image_components_by_sandwich=[r['degree'] for r in records if r['matching_lower_image']],engine='Native FLINT rational polynomial identities, separate weighted reducer and block matrix inverses; own mod-two augmented column ranks. Fresh SymPy opposition lattice check. No Sage or primary experiment imports.',seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),boundp,xp,fp)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (ba,xa,fa)});(ROOT/'experiments/e7-d5a1-refinement-bound-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED independent refinement bound',out['seconds'],flush=True)

if __name__=='__main__':main()
