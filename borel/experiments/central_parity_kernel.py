"""Actual integral kernel checks for central parity: Spin7 failure and type C.

Reflection kernels and BGG coefficient kernels are computed over ZZ before
reduction mod two. This avoids substituting the kernel of H -> H/2H for
the reduced integral kernel, which is precisely the Spin7 obstruction.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
from sage.all import ZZ,GF,matrix,vector,identity_matrix
from sage_borel import BorelMap
ROOT=Path(__file__).resolve().parents[1]

def exponents(n,d):
 if n==1:yield (d,);return
 for k in range(d+1):
  for e in exponents(n-1,d-k):yield (k,)+e

def encode(f):return [[list(e),str(c)] for e,c in sorted(f.dict().items())]
def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def check(typ,rank,I,maximum):
 allmap=BorelMap([typ,rank]);b=BorelMap([typ,rank],[j for j in range(1,rank+1) if j not in I]);P=allmap.ring
 chi=[int(i==rank) for i in range(1,rank+1)] if typ=='B' else [i%2 for i in range(1,rank+1)]
 previous=None;previousK=None;previous_ix=None;rows=[];witness=None
 for d in range(maximum+1):
  es=list(exponents(rank,d));mons=[P({e:1}) for e in es];ix={e:j for j,e in enumerate(es)}
  R=matrix(ZZ,len(I)*len(es),len(es))
  for k,node in enumerate(I):
   for j,f in enumerate(mons):
    for e,c in (allmap.reflect(f,node)-f).dict().items():R[k*len(es)+ix[tuple(e)],j]=ZZ(c)
  U=R.right_kernel().basis_matrix().transpose() if I else identity_matrix(ZZ,len(es))
  polys=[sum(c*m for c,m in zip(v,mons)) for v in U.columns()]
  words=b.basis(d);E=matrix(ZZ,len(words),len(polys))
  for j,f in enumerate(polys):
   out=b.expand(f)
   for i,w in enumerate(words):E[i,j]=ZZ(out.coefficient(w))
  K=E.right_kernel().basis_matrix().transpose();actualK=U*K
  row=dict(degree=d,source_rank=int(U.ncols()),integral_kernel_rank=int(K.ncols()))
  if d:
   delta=matrix(GF(2),previous.nrows(),len(es))
   for j,e in enumerate(es):
    for t,power in enumerate(e):
     if power%2 and chi[t]:
      ee=list(e);ee[t]-=1;delta[previous_ix[tuple(ee)],j]+=1
   prev=matrix(GF(2),previous);prevK=matrix(GF(2),previousK)
   sourceok=prev.augment(delta*matrix(GF(2),U)).rank()==prev.ncols()
   kernelok=prevK.augment(delta*matrix(GF(2),actualK)).rank()==prevK.ncols()
   row.update(actual_source_preserved=bool(sourceok),actual_integral_kernel_preserved=bool(kernelok))
   if not kernelok:
    for j in range(actualK.ncols()):
     image=delta*matrix(GF(2),actualK.column(j).column())
     if prevK.augment(image).rank()>prevK.ncols():
      f=P(sum(c*m for c,m in zip(actualK.column(j),mons)))
      df=sum(c*f.derivative(x) for c,x in zip(chi,P.gens()))
      assert not b.expand(f)
      witness=dict(polynomial=encode(f),derivative=encode(df),derivative_Schubert_terms=b.expand(df).to_terms(),previous_integral_kernel=[encode(sum(c*P({e:1}) for c,e in zip(v,previous_es))) for v in previousK.columns()])
      break
   if typ=='C':assert sourceok and kernelok
  rows.append(row)
  previous=U;previousK=actualK;previous_ix=ix;previous_es=es
 return dict(type=typ,rank=rank,retained_nodes=list(I),maximum_degree=maximum,degrees=rows,witness=witness)

def main():
 start=perf_counter();spin=check('B',3,(),4)
 assert spin['degrees'][4]['actual_integral_kernel_preserved'] is False
 assert all(r.get('actual_integral_kernel_preserved',True) for r in spin['degrees'][:4])
 c=[check('C',3,I,d) for I,d in (((),10),((2,3),6),((1,2),7))]
 out=dict(result='passed',Spin7=spin,type_C_checks=c,seconds=round(perf_counter()-start,3),scope='Complete actual ZZ reflection/BGG kernel tests in the listed degrees, including a Spin7/B counterexample in degree four. Type C all-rank and all-parabolic stability is proved from the classical presentation in the note; finite checks are corroboration.',source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in (Path(__file__).resolve(),ROOT/'src/sage_borel/core.py',ROOT/'src/sage_borel/expansion.py')})
 (ROOT/'experiments/central-parity-kernel.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],flush=True)
if __name__=='__main__':main()
