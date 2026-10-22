"""Native polynomial proof of the Spin7 obstruction and classical identities.

SymPy polynomial identities and monic division replace the primary BGG
image and integer-kernel computation. No Sage, runtime or Smith tables.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main():
 u,v,z=s.symbols('u v z');eps=(u,v-u,2*z-v)
 q=s.expand(sum(x*x for x in eps)/2)
 r4=s.expand(sum(eps[i]**2*eps[j]**2 for i in range(3) for j in range(i+1,3)))
 f=-u**4+2*u**3*v-2*u**2*v**2+u*v**3+u**2*v*z-u*v**2*z+z**4
 t=-5*u**2+5*u*v-v**2+2*v*z+2*z*z
 assert s.expand(4*f-r4-q*t)==0
 assert q==u*u-u*v+v*v-2*v*z+2*z*z
 derivative=s.diff(f,z);g=u*u*v+u*v*v
 assert s.Poly(derivative-g,u,v,z,modulus=2).is_zero
 quotient,remainder=s.reduced(g,[q],u,v,z,modulus=2)
 assert quotient==[v] and remainder==v**3
 roots=s.Matrix([[2,-1,0],[-1,2,-1],[0,-2,2]])
 assert roots.det()==2 and all(int(roots[2,j])%2==0 for j in range(3))
 # A monic primitive quadratic times arbitrary linear forms is the entire
 # degree-three rational-kernel intersection. Its mod-two column rank is 3.
 mons=[u**i*v**j*z**(3-i-j) for i in range(4) for j in range(4-i)]
 columns=[]
 for x in (u,v,z):
  f3=s.Poly(q*x,u,v,z,modulus=2);columns.append(sum((int(f3.coeff_monomial(m))%2)<<i for i,m in enumerate(mons)))
 piv={}
 for col in columns:
  while col:
   j=col.bit_length()-1
   if j in piv:col^=piv[j]
   else:piv[j]=col;break
 assert len(piv)==3
 # Direct checks of the formulas used in the arbitrary-rank written proofs.
 rows=[]
 for n in range(1,7):
  xs=s.symbols('x:'+str(n))
  for k in range(1,n+1):
   e=sum(s.prod(c) for c in combinations(xs,k));delta=sum(s.diff(e,x) for x in xs)
   previous=sum(s.prod(c) for c in combinations(xs,k-1))
   assert s.expand(delta-(n-k+1)*previous)==0
   square_e=sum(s.prod(x*x for x in c) for c in combinations(xs,k))
   assert s.Poly(sum(s.diff(square_e,x) for x in xs),*xs,modulus=2).is_zero
  rows.append(dict(number_of_variables=n,elementary_translation_formula=True,square_invariants_killed_mod2=True))
 primary_path=ROOT/'experiments/central-parity-kernel.json';primary=json.loads(primary_path.read_text());assert primary['result']=='passed'
 case=primary['Spin7'];assert case['degrees'][4]['actual_integral_kernel_preserved'] is False
 def decode(terms):return sum(s.Rational(c)*u**e[0]*v**e[1]*z**e[2] for e,c in terms)
 assert s.expand(decode(case['witness']['polynomial'])-f)==0
 assert s.expand(decode(case['witness']['derivative'])-derivative)==0
 out=dict(result='passed',Spin7=dict(root_lattice_index=2,fundamental_weight_parity=[0,0,1],kernel_relation=str(f),identity='4F = e2(epsilon^2) + q*(-5u^2+5uv-v^2+2vz+2z^2)',derivative_mod2=str(g),primitive_cubic_kernel_rank=3,monic_division_remainder='v^3',integral_kernel_preservation=False),classical_identity_checks=rows,scope='Exact polynomial proof of the actual Spin7/B kernel failure, independent of BGG. Classical generator derivative identities checked in up to six variables; all-rank type A/C proofs are in the note. No complete classification or novelty claim.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path)})
 (ROOT/'experiments/central-parity-kernel-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED: Spin7 counterexample and classical generator identities')
if __name__=='__main__':main()
