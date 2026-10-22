"""Small maximal-rank A5+A2 intermediary for E7 integral Borel images."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from functools import lru_cache
from math import comb
import gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,matrix,vector,prod,pari,inverse_mod
ROOT=Path(__file__).resolve().parents[1]
WS=(2,3,3,4,5,6)

def main():
 start=perf_counter();C=matrix(QQ,7,7,[2,0,-1,0,0,0,0,0,2,0,-1,0,0,0,-1,0,2,-1,0,0,0,0,-1,-1,2,-1,0,0,0,0,0,-1,2,-1,0,0,0,0,0,-1,2,-1,0,0,0,0,0,-1,2]);chain=[1,3,4,5,6];A5=C.matrix_from_rows_and_columns(chain,chain);theta=C*vector(QQ,[2,2,3,4,3,2,1]);roots=[C.column(j) for j in chain];om=[sum((A5.inverse()[i,j]*roots[j] for j in range(5)),vector(QQ,7)) for i in range(5)];a2=[C.column(0),-theta];CI=matrix(QQ,[[2,-1],[-1,2]]).inverse();eta=[sum((CI[i,j]*a2[j] for j in range(2)),vector(QQ,7)) for i in range(2)];assert abs(matrix(QQ,roots+a2).transpose().det()/C.det())==3;XB=matrix(QQ,om+eta).transpose();QB=matrix(QQ,roots+eta).transpose();assert XB.det().valuation(2)==0 and (C.inverse()*QB).det().valuation(2)==0 and all(c.denominator()%2 for c in C.inverse()*QB)
 xx=[om[0]]+[om[i]-om[i-1] for i in range(1,5)]+[-om[4]];zz=[eta[0],eta[1]-eta[0],-eta[1]];orbit={(0,0,0,0,0,0,1)};todo=list(orbit)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=tuple(v[j]-v[i]*C[j,i] for j in range(7))
   if u not in orbit:orbit.add(u);todo.append(u)
 from itertools import combinations
 ext={tuple(sum((xx[i] for i in ids),vector(QQ,7))) for ids in combinations(range(6),3)};branch=None
 for sign in (1,-1):
  proposed=ext|{tuple(eps*(x+sign*z)) for x in xx for z in zz for eps in (-1,1)}
  if proposed==orbit:branch=sign;break
 assert branch is not None and len(orbit)==56;print('A5+A2 branching verified',branch,round(perf_counter()-start,2),flush=True)
 P=PolynomialRing(QQ,6,names=('q2','q3','e3','e4','e5','e6'),order=TermOrder('wdegrevlex',WS));q2,q3,e3,e4,e5,e6=P.gens();es=[P.one(),P.zero(),-q2,e3,e4,e5,e6];ps=[P(6)]
 for k in range(1,19):ps.append(P(sum((-1)**(j-1)*es[j]*ps[k-j] for j in range(1,min(k,7)))+((-1)**(k-1)*k*es[k] if k<=6 else 0)))
 zes=[P.one(),P.zero(),q2,q3];zp=[P(3)]
 for k in range(1,19):zp.append(P(sum((-1)**(j-1)*zes[j]*zp[k-j] for j in range(1,min(k,4)))+((-1)**(k-1)*k*zes[k] if k<=3 else 0)))
 moments={}
 for d in (2,6,8,10,12,14,18):
  triple=sum(comb(d,j)*comb(d-j,k)*ps[j]*ps[k]*ps[d-j-k] for j in range(d+1) for k in range(d-j+1));double=sum(comb(d,j)*2**j*ps[j]*ps[d-j] for j in range(d+1));extmoment=(triple-3*double+2*3**d*ps[d])/6;moments[d]=P(extmoment+2*sum(comb(d,j)*branch**(d-j)*ps[j]*zp[d-j] for j in range(d+1)))
 assert not moments[2];relations=[moments[d] for d in (6,8,10,12,14,18)];G=P.ideal(relations).groebner_basis();assert P.ideal(G).vector_space_dimension()==672;leads=[tuple(g.lm().exponents()[0]) for g in G];found={(0,)*6};todo=list(found)
 while todo:
  e=todo.pop()
  for j in range(6):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   if u not in found and not any(all(a<=b for a,b in zip(l,u)) for l in leads):found.add(u);todo.append(u)
 norm={d:sorted(e for e in found if sum(a*b for a,b in zip(e,WS))==d) for d in range(46)};assert len(found)==672 and sum(map(len,norm.values()))==672;print('A5+A1 quotient',len(G),'GB elements',round(perf_counter()-start,2),flush=True)
 sp=ROOT/'experiments/a5-two-local-source.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];assert sha256(sa.read_bytes()).hexdigest()==sr['certificate_archive_sha256'];sc=json.loads(gzip.decompress(sa.read_bytes()));algebra=ROOT/'experiments/a5-two-local-algebra.json';assert json.loads(algebra.read_text())['full_source_complete_all_degrees'];factors=list(map(QQ,sc['normalization_factors']));target=[-q2,e3,e4,e5,e6];gens=dict(N=[(f,d) for f,d in zip(P.gens(),WS)],M=[(q2,2),(q3,3)])
 for row in sc['degrees']:
  mons=[prod((f/k)**a for f,k,a in zip(target,factors,e)).reduce(G) for e in row['classical_exponents']];CC=matrix(QQ,row['rational_coordinates'])
  for j in row['new_generator_indices']:gens['M'].append((P(sum(f*CC[i,j] for i,f in enumerate(mons))).reduce(G),row['degree']))
 bases={k:{0:matrix(QQ,1,1,[1])} for k in gens};records=[];certrows=[]
 @lru_cache(None)
 def multiplication(kind,i,d):
  f,k=gens[kind][i];ix={e:j for j,e in enumerate(norm[d])};A=matrix(QQ,len(norm[d]),len(norm[d-k]))
  for j,e in enumerate(norm[d-k]):
   for u,c in (f*P({e:1})).reduce(G).dict().items():A[ix[tuple(u)],j]=c
  return A
 def local_basis(A):
  k=max([0]+[int(c.denominator().valuation(2)) for c in A.list()]);scaled=2**k*A
  for precision in (64,128,256):
   mod=ZZ(2)**precision;Z=matrix(ZZ,A.nrows(),A.ncols(),[ZZ(c.numerator()*inverse_mod(c.denominator(),mod))%mod for c in scaled.list()]);H=matrix(QQ,pari(Z).mathnfmodid(mod))/2**k;B=H.inverse()*A
   if B.denominator()%2 and matrix(GF(2),B).rank()==H.nrows():return H
  raise AssertionError('local precision')
 for d in range(1,46):
  columns={}
  for kind,gs in gens.items():
   blocks=[multiplication(kind,i,d)*bases[kind][d-k] for i,(f,k) in enumerate(gs) if k<=d];A=blocks[0] if blocks else matrix(QQ,len(norm[d]),0)
   for block in blocks[1:]:A=A.augment(block)
   bases[kind][d]=local_basis(A) if norm[d] else matrix(QQ,0,0);columns[kind]=A
  B=bases['N'][d].inverse()*bases['M'][d];assert B.denominator()%2;vals=[int(c.valuation(2)) for c in matrix(ZZ,B.denominator()*B).elementary_divisors()];counts={str(k):vals.count(k) for k in sorted(set(vals)) if k};row=dict(degree=d,rank=len(norm[d]),smith_multiplicities=counts);records.append(row);certrows.append(dict(**row,image_columns={k:[[str(c) for c in line] for line in A] for k,A in columns.items()},image_bases={k:[[str(c) for c in line] for line in bases[k][d]] for k in gens}));print('A5+A1 image',row,round(perf_counter()-start,2),flush=True)
 cert=dict(branch_sign=branch,subsystem_A5_roots=[list(map(str,v)) for v in roots],subsystem_A2_roots=[list(map(str,v)) for v in a2],A5_fundamental_weights=[list(map(str,v)) for v in om],A2_fundamental_weights=[list(map(str,v)) for v in eta],relations=list(map(str,relations)),groebner_basis=list(map(str,G)),normal_exponents={str(d):[list(e) for e in ee] for d,ee in norm.items()},generators={k:[[str(f),d] for f,d in gs] for k,gs in gens.items()},degrees=certrows);archive=ROOT/'experiments/e7-a5a2-images-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));report=dict(result='passed',complete=True,subsystem='A5+A2',top_polynomial_degree=45,total_rank=672,subsystem_index=3,branch_sign=branch,degrees=records,seconds=round(perf_counter()-start,3),scope='Complete primary Borel images of full invariants for the odd-index maximal-rank A5+A2 subgroup, embedded in full-flag cohomology. It is an intermediary, not a standard parabolic. Common integral A2 fibers transfer the comparison to A5+A1 and A5 parabolics. Independent audits are required separately.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),sp,algebra)},raw_sha256={sa.relative_to(ROOT).as_posix():sha256(sa.read_bytes()).hexdigest()});(ROOT/'experiments/e7-a5a2-images.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
