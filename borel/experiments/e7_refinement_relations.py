"""Sage proof certificate: common quadratic and first higher E7 relation.

Independent input route: 56 minuscule weights, rather than E6 psi images.
"""
from sage.all import QQ,ZZ,GF,PolynomialRing,CartanMatrix,matrix,vector,gcd,prod
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
ROOT=Path(__file__).resolve().parents[1]
def main():
 start=perf_counter();C=matrix(QQ,CartanMatrix('E7'));XB=C.matrix_from_columns(range(6)).augment(vector(QQ,[0]*6+[1]).column())
 assert XB.det()==3
 R=PolynomialRing(QQ,7,names=['a'+str(i) for i in range(1,7)]+['w']);a=R.gens()[:6];w=R.gen(6)
 B=-QQ(3)/2*(vector(R,a)*C[:6,:6].inverse()*vector(R,a))
 quad=w*w-B
 # Reduce as a polynomial in w explicitly; avoid dependence on GB term order.
 def reduce_quad(f):
  result=R.zero()
  for e,c in R(f).dict().items():
   e=tuple(e);k=e[6];result+=c*R({tuple(e[:-1])+(k%2,):1})*B**(k//2)
  return R(result)
 orbit={(0,0,0,0,0,0,1)};pending=list(orbit)
 while pending:
  t=pending.pop()
  for i in range(7):
   u=tuple(t[j]-t[i]*C[i,j] for j in range(7))
   if u not in orbit:orbit.add(u);pending.append(u)
 assert len(orbit)==56
 forms=[sum(c*x for c,x in zip(XB.inverse()*vector(QQ,t),R.gens())) for t in sorted(orbit)]
 assert reduce_quad(sum(f*f for f in forms))==0, (B,reduce_quad(sum(f*f for f in forms)))
 moment=reduce_quad(sum(f**6 for f in forms));den=moment.denominator();content=gcd((den*moment).coefficients());F=R(den*moment/content)
 even=F.subs({w:0});odd=R(sum(c*R({tuple(e[:-1])+(0,):1}) for e,c in F.dict().items() if e[6]==1))
 assert F==even+w*odd and all(c.denominator()==1 for c in F.coefficients())
 assert min(ZZ(c).valuation(2) for c in even.coefficients())==0
 odd_val=min(ZZ(c).valuation(2) for c in odd.coefficients());assert odd_val==2
 # First test the octic relation for an extension of the truncation.
 octic=reduce_quad(sum(f**8 for f in forms));den8=octic.denominator();g8=gcd((den8*octic).coefficients());H=R(den8*octic/g8)
 for step in range(8):
  modring=PolynomialRing(GF(2),7,R.variable_names())
  if modring(H)==modring(B*F):H=(H+B*F)/2
  else:break
 even8=H.subs({w:0});odd8=R(sum(c*R({tuple(e[:6])+(0,):1}) for e,c in H.dict().items() if e[6]==1))
 val8=min(ZZ(c).valuation(2) for c in odd8.coefficients())
 assert val8>=1 and H.denominator()==1
 def powers(n,k=6):
  if k==1:yield (n,);return
  for j in range(n+1):
   for rest in powers(n-j,k-1):yield (j,)+rest
 def module_monomials(d):
  if d<0:return []
  return [prod(x**j for x,j in zip(a,e)) for e in powers(d)]+([w*prod(x**j for x,j in zip(a,e)) for e in powers(d-1)] if d else [])
 kernel_checks=[]
 for degree in range(6,10):
  polys=[reduce_quad(F*f) for f in module_monomials(degree-6)]+[reduce_quad(H*f) for f in module_monomials(degree-8)]
  maps=[f.dict() for f in polys];support=sorted(set().union(*(set(t) for t in maps)))
  mat=matrix(GF(2),[[t.get(e,0) for t in maps] for e in support]);assert mat.rank()==len(polys)
  kernel_checks.append(dict(degree=degree,columns=len(polys),primitive_mod2_rank=int(mat.rank())))
 # E6 reflections preserve the primitive relation in the quadratic quotient.
 for i in range(6):
  new=[R.gen(j)-C[j,i]*R.gen(i) if j<6 else w for j in range(7)]
  assert F(*new)==F and H(*new)==H and B(*new)==B
 F2=PolynomialRing(GF(2),7,R.variable_names());bar=F2(quad)
 # Nonzero mixed terms in B exclude a square, so w^2-B is irreducible.
 cross=[e for e,c in F2(B).dict().items() if sum(k%2 for k in e)]
 assert cross and F2(odd)==0 and F2(even)!=0
 # Mod-2 polar matrix of B is nonsingular on the six-dimensional V.
 polar=matrix(GF(2),6,6,lambda i,j:B.derivative(a[i]).derivative(a[j]))
 assert polar.rank()==6
 def encoded(f):return [[list(e),str(c)] for e,c in sorted(f.dict().items())]
 out=dict(result='passed',method='56-weight Weyl orbit, sixth/eighth moments, monic quadratic reduction, primitive integer coefficient normalization',orbit_size=56,common_complement_index=3,quadratic_polar_rank_mod2=6,first_higher_relation_degree=6,odd_coefficient_2_valuation=int(odd_val),quadratic=encoded(quad),sextic=encoded(F),octic=encoded(H),octic_odd_valuation=int(val8),kernel_primitivity=kernel_checks,seconds=round(perf_counter()-start,3),scope='Coefficient certificate for the uniform proof through degree 9 on all 64 refinements, not an all-degree multiplicity calculation.',source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
 (ROOT/'experiments/e7-refinement-relations.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],'primitive kernel degrees 6..9',flush=True)
if __name__=='__main__':main()
