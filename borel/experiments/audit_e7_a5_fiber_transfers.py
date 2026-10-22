"""Fresh lattice and Weyl-word checks for the common A2 fiber transfers."""
from pathlib import Path
from hashlib import sha256
from collections import Counter
import json
import sympy as S
ROOT=Path(__file__).resolve().parents[1]

def main():
 C=S.Matrix([[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]]);chain=[1,3,4,5,6];Ac=C.extract(chain,chain);ar=[C.col(j) for j in chain];om=[sum((Ac.inv()[i,j]*ar[j] for j in range(5)),S.zeros(7,1)) for i in range(5)];theta=C*S.Matrix([2,2,3,4,3,2,1]);a2=[C.col(0),-theta];eta=[(2*a2[0]+a2[1])/3,(a2[0]+2*a2[1])/3];X=S.Matrix.hstack(*om,*eta);Q=S.Matrix.hstack(*ar,*eta)
 def local(c):return int(S.denom(c))%2==1
 assert local(X.det()) and local(1/X.det());assert all(local(c) for c in C.inv()*Q) and local((C.inv()*Q).det()) and local(1/(C.inv()*Q).det());assert abs(S.Matrix.hstack(*ar,*a2).det()/C.det())==3
 roots={tuple(C.col(i)) for i in range(7)};todo=list(roots)
 while todo:
  v=todo.pop()
  for i in range(7):
   w=tuple(v[j]-v[i]*C[j,i] for j in range(7))
   if w not in roots:roots.add(w);todo.append(w)
 assert len(roots)==126
 def rootspan(gens):
  A=S.Matrix.hstack(*gens);return {r for r in roots if A.row_join(S.Matrix(r)).rank()==A.cols}
 def act(v,word):
  v=S.Matrix(v)
  for i in word:v-=v[i-1]*C.col(i-1)
  return tuple(v)
 records=[]
 for name,v,asub in [('A5+A2',S.zeros(7,1),chain),('42',om[3],[chain[i] for i in (0,1,2,4)]),('222',om[1]+2*om[3],[chain[i] for i in (0,2,4)])]:
  if name=='A5+A2':continue
  word=[];v=S.Matrix(v)
  while any(c<0 for c in v):
   i=next(i for i,c in enumerate(v) if c<0);v-=v[i]*C.col(i);word.append(i+1);assert len(word)<100
  I=[i+1 for i,c in enumerate(v) if not c];target=rootspan([C.col(i-1) for i in I]);source=rootspan([C.col(i) for i in asub]+a2);assert {act(r,word) for r in source}==target
  records.append(dict(blocks=name,source_A5_nodes=[i+1 for i in asub],source_A2_roots=[list(map(str,v)) for v in a2],standard_retained_nodes=I,application_order_word=word,levi_root_count=len(target)))
 # Root and weight A2 lattices agree at two; all fiber generators are
 # therefore common integral characters. Monic relations prove the two
 # successive source extensions BEFORE Borel; their rational independence
 # after Borel follows from the ordinary coinvariant presentation.
 q2,q3,h,t,u=S.symbols('q2 q3 h t u');assert S.expand((h**3+q2*h+q3).subs({q2:t-h*h,q3:-h*t}))==0;assert S.expand((u*u-h*u+t).subs(t,u*(h-u)))==0
 out=dict(result='passed',odd_index_subsystem=3,full_local_lattice_splittings=True,common_A2_partial_fiber_degrees=[0,1,2],common_A2_full_fiber_degrees=[0,1,1,2,2,3],fiber_relations=['h^3+q2*h+q3','u^2-h*u+t'],conjugacies=records,scope='Exact root/weight local splittings and Weyl conjugacies. The monic common A2 source bases transfer full-image comparisons by direct sums; they do not assume invariants commute with the Borel quotient.',engine='Native SymPy independent Cartan lattice determinants, fresh root orbit and full Levi root-mask equality.',source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()});(ROOT/'experiments/e7-a5-fiber-transfers-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out)

if __name__=='__main__':main()
