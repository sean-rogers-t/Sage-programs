"""Direct full invariant-kernel/Borel comparison in E7/P6 degree two."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,prod,pari,inverse_mod,gcd
ROOT=Path(__file__).resolve().parents[1]
C=matrix(ZZ,[[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]])
I=(0,1,2,3,4,6)

def main():
 X=PolynomialRing(QQ,7,names=['w'+str(i+1) for i in range(7)]);Q=PolynomialRing(QQ,7,names=['a'+str(i+1) for i in range(7)]);omega=X.gens();alpha=[sum(C[i,j]*omega[j] for j in range(7)) for i in range(7)];embed=Q.hom(alpha,X);exps=[tuple(int(k==i)+int(k==j) for k in range(7)) for i in range(7) for j in range(i,7)]
 source={};polys={}
 for kind,P in (('N',X),('M',Q)):
  vs=P.gens();mons=[P({e:1}) for e in exps];blocks=[]
  for i in I:
   reflected=[vs[j]-(alpha[i] if kind=='N' else C[i,j]*vs[i]) if kind=='M' or j==i else vs[j] for j in range(7)];phi=P.hom(reflected,P);blocks.append(matrix(ZZ,28,28,{(exps.index(tuple(e)),j):c for j,f in enumerate(mons) for e,c in (phi(f)-f).dict().items()}))
  A=blocks[0]
  for block in blocks[1:]:A=A.stack(block)
  K=A.right_kernel().basis_matrix().transpose();assert K.ncols()==3 and A*K==0 and matrix(GF(2),K).rank()==3
  fs=[P({e:c for e,c in zip(exps,col) if c}) for col in K.columns()];polys[kind]=[embed(f) if kind=='M' else f for f in fs];source[kind]=dict(invariant_polynomials=list(map(str,fs)),kernel_columns=[[int(c) for c in row] for row in K])
 roots={tuple(C.row(i)) for i in range(7)};todo=list(roots)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=tuple(v[j]-v[i]*C[i,j] for j in range(7))
   if u not in roots:roots.add(u);todo.append(u)
 assert len(roots)==126;quad=sum(sum(c*w for c,w in zip(v,omega))**2 for v in roots);quad/=gcd(quad.coefficients());qc=[quad.monomial_coefficient(X({e:1})) for e in exps];p=next(j for j,c in enumerate(qc) if c);images={}
 for kind,fs in polys.items():
  reduced=[f-f.monomial_coefficient(X({exps[p]:1}))/qc[p]*quad for f in fs];images[kind]=matrix(QQ,28,3,{(exps.index(tuple(e)),j):c for j,f in enumerate(reduced) for e,c in f.dict().items()})
 piv=list(images['N'].transpose().pivots());assert len(piv)==2
 def basis(A):
  k=max([0]+[int(c.denominator().valuation(2)) for c in A.list()]);mod=ZZ(2)**64;Z=matrix(ZZ,A.nrows(),A.ncols(),[ZZ(c.numerator()*inverse_mod(c.denominator(),mod))%mod for c in (2**k*A).list()]);H=matrix(QQ,pari(Z).mathnfmodid(mod))/2**k;D=H.inverse()*A;assert D.denominator()%2 and matrix(GF(2),D).rank()==2;return H
 images={kind:A.matrix_from_rows(piv) for kind,A in images.items()};bases={kind:basis(A) for kind,A in images.items()};T=bases['N'].inverse()*bases['M'];assert T.denominator()%2;vals=[int(c.valuation(2)) for c in matrix(ZZ,T.denominator()*T).elementary_divisors()];assert vals==[0,1]
 out=dict(result='passed',retained_nodes=[1,2,3,4,5,7],degree=2,full_source_kernel_ranks=dict(N=3,M=3),cohomology_rank=2,relative_smith_2_valuations=vals,sources=source,primitive_global_quadratic=str(quad),target_pivot_rows=piv,image_columns={kind:[[str(c) for c in row] for row in A] for kind,A in images.items()},image_bases={kind:[[str(c) for c in row] for row in H] for kind,H in bases.items()},scope='Full degree-two comparison only; higher degrees remain unresolved.',engine='Sage full integral kernels of all six reflection equations in each actual lattice, adjoint-root quadratic and exact local image bases.',source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
 (ROOT/'experiments/e7-d5a1-degree-two.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED full degree-two kernels: relative factors (1,2)')

if __name__=='__main__':main()
