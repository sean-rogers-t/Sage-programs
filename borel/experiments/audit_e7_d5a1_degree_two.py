"""Independent hand-basis proof and native audit of the E7/P6 first gap."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations_with_replacement
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def rank2(A):
 piv={}
 for col in zip(*A):
  b=sum((int(c.p)%2)<<j for j,c in enumerate(col))
  while b:
   k=b.bit_length()-1
   if k in piv:b^=piv[k]
   else:piv[k]=b;break
 return len(piv)

def local(A):return all(c.q%2 for c in A)

def main():
 path=ROOT/'experiments/e7-d5a1-degree-two.json';primary=json.loads(path.read_text());assert primary['result']=='passed'
 C=2*s.eye(7)
 for i,j in ((1,3),(3,4),(2,4),(4,5),(5,6),(6,7)):C[i-1,j-1]=C[j-1,i-1]=-1
 om=[s.eye(7).col(j) for j in range(7)];h=om[5];w=om[6];z=om[4];us=[om[0],om[2]-om[0],om[3]-om[2],om[1]+om[4]-om[3],om[4]-om[1]];assert sum(us,s.zeros(7,1))==2*z
 NB=s.Matrix.hstack(h,w,z,*us[:4]);QB=s.Matrix.hstack(h,2*z,w-z,*us[:4]);assert abs(NB.det())==1;RC=C.inv()*QB;assert abs(RC.det())==1 and all(c.q==1 for c in RC)
 ys=[2*u-h for u in us]
 for i in range(5):
  matches=[next((j,t) for j,y in enumerate(ys) for t in (-1,1) if v-v[i]*C.col(i)==t*y) for v in ys];assert len({j for j,t in matches})==5 and __import__('math').prod(t for j,t in matches)==1
 assert all(y[6]==0 for y in ys) and w-w[6]*C.col(6)==h-w
 omega=s.symbols('w1:8');alpha=s.symbols('a1:8');env=dict(zip(map(str,omega+alpha),omega+alpha));lin=lambda v:sum(a*x for a,x in zip(v,omega));hh,ww,zz=lin(h),lin(w),lin(z);E=sum(lin(y)**2 for y in ys);t=ww*(hh-ww)
 N=[hh**2,s.expand((E-5*hh**2)/8),t];M=[hh**2,s.expand((E-5*hh**2)/4),4*t];sub=dict(zip(omega,C.inv()*s.Matrix(alpha)))
 exps=[tuple(int(k==i)+int(k==j) for k in range(7)) for i,j in combinations_with_replacement(range(7),2)]
 def coefficients(fs,variables):
  dictionaries=[s.Poly(f,*variables).as_dict() for f in fs];return s.Matrix([[f.get(e,s.Rational(0)) for f in dictionaries] for e in exps])
 analytic=dict(N=coefficients(N,omega),M=coefficients([s.expand(f.subs(sub,simultaneous=True)) for f in M],alpha))
 for kind in ('N','M'):
  A=analytic[kind];assert all(c.q==1 for c in A) and rank2(A.tolist())==3
  vs=omega if kind=='N' else alpha;recorded=coefficients([s.sympify(f.replace('^','**'),locals=env) for f in primary['sources'][kind]['invariant_polynomials']],vs);B=recorded.gauss_jordan_solve(A)[0];assert recorded*B==A and local(B) and local(B.inv())
 # Fresh minimal-weight quadratic, independent of the primary adjoint roots.
 orbit={tuple(om[6])};todo=list(orbit)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=tuple(v[j]-v[i]*C[j,i] for j in range(7))
   if u not in orbit:orbit.add(u);todo.append(u)
 assert len(orbit)==56;moment=s.expand(sum(sum(a*x for a,x in zip(v,omega))**2 for v in orbit));relation=s.expand(E+3*hh**2-8*t);assert s.expand(moment-3*relation)==0
 quad=s.sympify(primary['primitive_global_quadratic'].replace('^','**'),locals=env);ratio=s.Poly(moment,*omega).LC()/s.Poly(quad,*omega).LC();assert s.expand(moment-ratio*quad)==0
 # In the rational Borel quotient E=8t-3h^2, the candidate bases have
 # images (h^2,t-h^2,t) and (h^2,2t-2h^2,4t), respectively.
 assert s.Matrix([[1,-1,0],[0,1,1]]).rank()==s.Matrix([[1,-2,0],[0,2,4]]).rank()==2
 # Reconstruct and verify both independently proposed local image bases.
 p=next(i for i,e in enumerate(exps) if s.Poly(quad,*omega).coeff_monomial(e));qc=coefficients([quad],omega);cols={};bases={}
 for kind in ('N','M'):
  fs=[s.sympify(f.replace('^','**'),locals=env) for f in primary['sources'][kind]['invariant_polynomials']]
  if kind=='M':fs=[s.expand(f.subs(dict(zip(alpha,C*s.Matrix(omega))),simultaneous=True)) for f in fs]
  A=coefficients(fs,omega);A=A-qc*A[p,:]/qc[p,0];A=A[primary['target_pivot_rows'],:];assert A==s.Matrix(primary['image_columns'][kind]);H=s.Matrix(primary['image_bases'][kind]);D=H.inv()*A;assert local(D) and rank2(D.tolist())==2;bases[kind]=H
 T=bases['N'].inv()*bases['M'];assert local(T) and local(2*T.inv()) and rank2(T.tolist())==1
 out=dict(result='passed',retained_nodes=[1,2,3,4,5,7],first_nonzero_relative_degree=2,cohomological_degree=4,degree_two_relative_group='Z/2',degree_two_full_source_bases=dict(N=['h^2','(E-5*h^2)/8','t'],M=['h^2','(E-5*h^2)/4','4*t']),degree_two_image_bases=dict(N=['h^2','t'],M=['h^2','2*t']),quadratic_relation='E+3*h^2-8*t=0',weight_orbit=56,engine='Native SymPy: independent integral lattice bases, analytic full-source bases with primitive coefficient ranks, fresh 56-weight quadratic, and complete reconstruction of the Sage image certificate. No Sage or project imports.',scope='All degrees 0,1,2 only: fixed characters agree in degree one; the first relative class is Z/2 in degree two. Higher multiplicities remain unresolved.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path)})
 (ROOT/'experiments/e7-d5a1-degree-two-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED analytic primitive bases and independent image audit: first gap degree two, Z/2')

if __name__=='__main__':main()
