"""Integer-only D4 polynomial, overlap and completed-family inventory audit."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from fractions import Fraction
import json
ROOT=Path(__file__).resolve().parents[1]

def multiply(f,g):
 h=[0]*(len(f)+len(g)-1)
 for i,a in enumerate(f):
  for j,b in enumerate(g):h[i+j]+=a*b
 return h

def quotient(f,g):
 assert g[0]==1;out=[0]*len(f)
 for i,a in enumerate(f):out[i]=a-sum(g[j]*out[i-j] for j in range(1,min(i+1,len(g))))
 while out and not out[-1]:out.pop()
 assert multiply(out,g)==f
 return out

def determinant(A):
 A=[list(row) for row in A];n=len(A);previous=1;sign=1
 for k in range(n-1):
  if not A[k][k]:
   p=next(i for i in range(k+1,n) if A[i][k]);A[k],A[p]=A[p],A[k];sign=-sign
  pivot=A[k][k];assert pivot
  for i in range(k+1,n):
   for j in range(k+1,n):value=pivot*A[i][j]-A[i][k]*A[k][j];assert value%previous==0;A[i][j]=value//previous
  previous=pivot
 return sign*A[-1][-1]

def inverse(A):
 n=len(A);B=[[Fraction(x) for x in row]+[Fraction(i==j) for j in range(n)] for i,row in enumerate(A)]
 for k in range(n):
  p=next(i for i in range(k,n) if B[i][k]);B[p],B[k]=B[k],B[p];c=B[k][k];B[k]=[x/c for x in B[k]]
  for i in range(n):
   if i!=k:
    c=B[i][k];B[i]=[x-c*y for x,y in zip(B[i],B[k])]
 out=[row[n:] for row in B];assert all(sum(A[i][k]*out[k][j] for k in range(n))==int(i==j) for i in range(n) for j in range(n));return out

def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]

def main():
 names=('e7-d4-images.json','e7-d4-images-audit.json','e7-d4-source-audit.json','e7-d5-summary.json','e7-refinement-independent.json','e7-full-flag-transfer.json');data={n:json.loads((ROOT/'experiments'/n).read_text()) for n in names};assert all(d['result']=='passed' for d in data.values())
 primary,audit,source,d5,older,transfer=[data[n] for n in names];assert primary['complete'] and source['full_source_complete_all_degrees']
 fiber=multiply([1]*5,[1,0,0,0,1]);f=multiply(d5['multiplicity_coefficients'],fiber)
 assert len(f)==52 and sum(f)==7560 and f[-1]==1 and max(f)==356
 assert f==[r['multiplicity_2'] for r in primary['degrees']]
 assert len(primary['degrees'])==len(audit['degree_rows'])==52
 assert all(r['prediction_agrees'] and r['indices']['N']==r['indices']['M'] for r in primary['degrees'])
 assert [r['indices'] for r in primary['degrees']]==[r['indices'] for r in audit['degree_rows']]
 ranks=[f[i]+(f[i+1] if i<51 else 0) for i in range(52)]
 assert ranks==list(reversed(ranks)) and sum(ranks)==15120
 assert ranks==[r['rank'] for r in primary['degrees']]==[r['rank'] for r in audit['degree_rows']]
 C=[a-b for a,b in zip(ranks,f)]
 while C and not C[-1]:C.pop()
 assert multiply(C,[1,1])==ranks and [0]+C==f
 row=next(r for r in older['refinements'] if r['retained_nodes']==[2,3,4,5]);assert row['multiplicity_2']==f[:10]
 PE7=[1]
 for d in (2,6,8,10,12,14,18):PE7=multiply(PE7,[1]*d)
 PD4=[1]
 for d in (2,4,4,6):PD4=multiply(PD4,[1]*d)
 assert quotient(quotient([0]+PE7,PD4),[1,1])==f
 settled={tuple(r['retained_nodes']) for r in transfer['good_levi_parabolics']};settled.update([(1,2,3,4,5,6),(1,2,3,4,5),(2,3,4,5,6),(2,3,4,5)])
 family=[]
 for row in older['refinements']:
  I=tuple(row['retained_nodes']);assert I in settled;PI=[1]
  for d in row['invariant_degrees']:PI=multiply(PI,[1]*d)
  flag=quotient(PE7,PI);multiplicity=quotient([0]+flag,[1,1]);assert all(x>=0 for x in multiplicity)
  assert row['multiplicity_2']==(multiplicity+[0]*10)[:10]
  family.append(dict(retained_nodes=list(I),multiplicity_coefficients=multiplicity,total_relative_factors_2=sum(multiplicity),cohomology_rank=sum(flag),dimension=len(flag)-1))
 assert len(family)==64 and len(settled)==112
 remaining=[list(I) for n in range(7) for I in combinations(range(1,8),n) if I not in settled];assert len(remaining)==15 and all(7 in I for I in remaining)
 center_nodes={2,5,7};central_scope={I for n in range(7) for I in combinations(range(1,8),n) if not center_nodes<=set(I)};assert settled==central_scope
 cartan=[[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]]
 assert all(sum(row[j-1] for j in center_nodes)%2==0 for row in cartan)
 assert determinant(cartan)==2 and all(center_nodes<=set(I) for I in remaining)
 ci=inverse(cartan);lifts=[]
 for node,n0 in ((2,7),(5,15),(7,3)):
  K=[j for j in range(7) if j!=node-1];CK=[[cartan[i][j] for j in K] for i in K];n=determinant(CK);assert n==n0 and ci[node-1][node-1]==Fraction(n,2);CKi=inverse(CK)
  assert all(c.denominator%2 for row in CKi for c in row) and all((CKi[j][j]/2).denominator%2 for j in range(6))
  basis=[cartan[j] for j in K]+[[int(j==node-1) for j in range(7)]];assert abs(determinant(basis))==n
  gram_inverse=[[CKi[i][j] if i<6 and j<6 else Fraction(2,n) if i==j==6 else Fraction(0) for j in range(7)] for i in range(7)]
  assert mm(mm(list(map(list,zip(*basis))),gram_inverse),basis)==cartan
  for j in range(6):
   S=[[Fraction(int(a==b))-CK[a][j]*int(b==j) for b in range(6)] for a in range(6)]
   assert mm(mm(list(map(list,zip(*S))),CKi),S)==CKi
  assert all((-Fraction(n,2)*CKi[j][j]).denominator==1 for j in range(6)) and all((-n*CKi[i][j]).denominator==1 for i in range(6) for j in range(i))
  lifts.append(dict(omitted_node=node,weight_and_root_split_odd_index=n,quadratic_omega_norm=str(Fraction(n,2)),quadratic_lift_diagonal_coefficients=[str(-Fraction(n,2)*CKi[j][j]) for j in range(6)],quadratic_lift_off_diagonal_coefficients=[str(-n*CKi[i][j]) for i in range(6) for j in range(i)],all_coefficients_integral=True,all_coefficients_2_local=True,all_simple_reflections_checked=True,global_quadratic_identity_checked=True))
 report=dict(result='passed',formula='q P_W(E7)/((1+q) P_W(D4))',multiplicity_coefficients=f,cohomology_ranks=ranks,mod_2_root_subimage_ranks=C,quadratic_algebra_corollary='On all 112 cases, choose any omitted i in {2,5,7}. N/2N is free on (1,beta(omega_i)) over C=image(M -> N/2N), with beta(omega_i)^2=beta(B_i). Formal consequence of the common source split, exponent bound, and certified Hilbert dimensions; C is not M/2M.',quadratic_lifts=lifts,total_relative_factors_2=7560,relative_module_generators=34,all_correction_indices_agree=True,proper_parabolics_settled=112,P7_refinements_settled=64,P7_refinements_remaining=[],other_proper_parabolics_remaining=remaining,central_character_scope='Exactly the 112 retained sets I that do not contain all of {2,5,7}; equivalently X^WI -> X/Q is onto. The other 15 have N^1=M^1, so the same polynomial cannot apply unchanged.',completed_P7_family=family,scope='D4 completes all 64 P7 refinements using prior P7/A6/D5 proofs plus the new full-source module/index proof. All-degree on this finite family, not a uniform all-types theorem. Other 15 proper E7 parabolics, ambient H/N and priority remain separate.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+[ROOT/'experiments'/n for n in names]})
 (ROOT/'experiments/e7-d4-summary.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED D4: 7560 factors; 112 proper E7 parabolics and all 64 P7 refinements')

if __name__=='__main__':main()
