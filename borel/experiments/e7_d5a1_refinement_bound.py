"""Bound the full E7/P6 root image by the known D5 refinement image.

For J={1,..,5} subset I={1,..,5,7}, pullback gives
L_I subset M_I subset N_I intersect M_J. Equality of the two finite
indices proves the full image, without exactness of taking invariants
after the Borel quotient. The known D5 pair is a shared audited input.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import gzip,json
from sage.all import QQ,GF,PolynomialRing,TermOrder,matrix,block_diagonal_matrix
ROOT=Path(__file__).resolve().parents[1]
WS=(1,2,4,5,6,8)

def parse(text,P):
 ix={v:i for i,v in enumerate(P.variable_names())};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=QQ(-1 if term.startswith('-') else 1);e=[0]*P.ngens()
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in ix:e[ix[name]]+=int(power) if power else 1
   else:c*=QQ(f)
  e=tuple(e);out[e]=out.get(e,0)+c
 return P(out)

def read(name):
 p=ROOT/'experiments'/name;r=json.loads(p.read_text());assert r['result']=='passed';a=ROOT/r['certificate_archive'];assert sha256(a.read_bytes()).hexdigest()==r['certificate_archive_sha256'];return p,r,a,json.loads(gzip.decompress(a.read_bytes()))

def main():
 start=perf_counter();xp,x,xa,xc=read('e7-d5a1-images.json');fp,f,fa,fc=read('e7-d5-images.json');P=PolynomialRing(QQ,6,names=('h','t','v4','p5','v6','v8'),order=TermOrder('wdegrevlex',WS));h,t,v4,p5,v6,v8=P.gens();G=[parse(g,P) for g in xc['groebner_basis']];norm={d:[tuple(e) for e in xc['normal_exponents'].get(str(d),[])] for d in range(44)};bases={0:matrix(QQ,1,1,[1]),43:matrix(QQ,0,0)}
 for row in xc['degrees']:bases[row['degree']]=matrix(QQ,row['image_bases']['N'])
 # The longest E6 element sends the original D5 to nodes 1,..,5:
 # old h -> -Omega6=-h+4w/3, old y_i -> -y_(6-i), old w -> w.
 targets=[(-h,P(QQ(4)/3)),(P.zero(),P.one()),(v4,P.zero()),(-p5,P.zero()),(v6,P.zero()),(v8,P.zero())]
 def multiply(a,b):return ((a[0]*b[0]-t*a[1]*b[1]).reduce(G),(a[0]*b[1]+a[1]*b[0]+h*a[1]*b[1]).reduce(G))
 @lru_cache(None)
 def mon(e):
  if not any(e):return (P.one(),P.zero())
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return multiply(targets[j],mon(tuple(u)))
 def evaluate(text):
  O=PolynomialRing(QQ,6,names=('h','w','v4','v5','v6','v8'));out=[P.zero(),P.zero()]
  for e,c in parse(text,O).dict().items():
   a,b=mon(tuple(e));out[0]+=c*a;out[1]+=c*b
  return tuple(out)
 assert all(evaluate(g)==(P.zero(),P.zero()) for g in fc['relations']);print('all refinement relations vanish',round(perf_counter()-start,2),flush=True)
 records=[];certrows=[]
 for row in fc['image_degrees']:
  d=row['degree'];oldmons=[mon(tuple(e)) for e in fc['normal_exponents'][str(d)]];nn=len(norm[d]);prev=len(norm[d-1]);T=matrix(QQ,nn+prev,len(oldmons),{(i,j):a.monomial_coefficient(P({e:1})) for j,(a,b) in enumerate(oldmons) for i,e in enumerate(norm[d])})
  for j,(a,b) in enumerate(oldmons):
   for i,e in enumerate(norm[d-1]):T[nn+i,j]=b.monomial_coefficient(P({e:1}))
  assert T.is_square() and T.is_invertible();H=block_diagonal_matrix(bases[d],bases[d-1]);W=H.inverse()*T*matrix(QQ,row['weight_basis']);assert W.denominator()%2 and W.inverse().denominator()%2
  B=H.inverse()*T*matrix(QQ,row['root_basis']);assert B.denominator()%2 and (2*B.inverse()).denominator()%2;F=matrix(GF(2),B);K=F.transpose().right_kernel().basis_matrix();assert K.nrows()==row['multiplicity_2'];upper=K.matrix_from_columns(range(nn)).rank();lower=next((r['smith_multiplicities'].get('1',0) for r in xc['degrees'] if r['degree']==d),0);assert upper<=lower
  exact=upper==lower;r=dict(degree=d,rank=nn,refinement_rank=nn+prev,constructed_root_cokernel_multiplicity=lower,refinement_intersection_cokernel_multiplicity=upper,full_image_proved=exact);records.append(r);certrows.append(dict(**r,old_normal_images=[[str(a),str(b)] for a,b in oldmons],weight_transition=[[str(c) for c in line] for line in W],root_transition=[[str(c) for c in line] for line in B],annihilator_mod2=[[int(c) for c in line] for line in K]));print('refinement bound',r,round(perf_counter()-start,2),flush=True)
 complete=all(r['full_image_proved'] for r in records);archive=ROOT/'experiments/e7-d5a1-refinement-bound-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(dict(targets=[[str(a),str(b)] for a,b in targets],degrees=certrows),separators=(',',':')).encode(),mtime=0));report=dict(result='passed',complete=complete,all_degree_multiplicities_proved=complete,retained_nodes=[1,2,3,4,5,7],refinement_nodes=[1,2,3,4,5],degrees=records,seconds=round(perf_counter()-start,3),scope='An unconditional upper image from a known finer parabolic and an independently constructed lower image. Equality proves each recorded full image; inequalities remain unresolved. No invariant-lifting exactness is assumed.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),xp,fp)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (xa,fa)});(ROOT/'experiments/e7-d5a1-refinement-bound.json').write_text(json.dumps(report,indent=2)+'\n');print('COMPLETE',complete,flush=True)

if __name__=='__main__':main()
