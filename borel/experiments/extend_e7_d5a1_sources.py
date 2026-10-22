"""Extend full E7/P6 root components from a rechecked compact seed.

Output has a distinct stem so a longer experimental run does not silently
replace the bounded certificate used by existing image/audit reports.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import argparse,gzip,json
from sage.all import QQ,ZZ,GF,matrix,vector
from e7_d5a1_sources import context,powers,WEIGHTS,CLASSICAL
ROOT=Path(__file__).resolve().parents[1]

def parse(text,P):
 ix={name:i for i,name in enumerate(P.variable_names())};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=QQ(-1 if term.startswith('-') else 1);e=[0]*P.ngens()
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in ix:e[ix[name]]+=int(power) if power else 1
   else:c*=QQ(f)
  e=tuple(e);out[e]=out.get(e,0)+c
 return P(out)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=22);args=ap.parse_args();start=perf_counter();seedpath=ROOT/'experiments/e7-d5a1-sources.json';seed=json.loads(seedpath.read_text());seedarchive=ROOT/seed['certificate_archive'];assert sha256(seedarchive.read_bytes()).hexdigest()==seed['certificate_archive_sha256'];data=json.loads(gzip.decompress(seedarchive.read_bytes()));P,orig,factors=context();assert list(map(str,orig))==data['original_generators'];P2=P.change_ring(GF(2));basis={0:[P.one()]};basis2={0:[P2.one()]};coordinates={0:matrix(QQ,1,1,[1])};exponents={0:[(0,)*7]};rows=data['degrees'][:];gs=[(parse(f,P),d) for f,d in seed['generators']]
 @lru_cache(None)
 def mon(e):
  if not any(e):return P.one()
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return orig[j]*mon(tuple(u))
 def coefficients(fs,support,ring):
  ix={e:i for i,e in enumerate(support)};return matrix(ring,len(support),len(fs),{(ix[tuple(e)],j):c for j,f in enumerate(fs) for e,c in f.dict().items()})
 for row in rows:
  d=row['degree'];es=[tuple(e) for e in row['classical_exponents']];support=list(powers(WEIGHTS,d));C=matrix(QQ,row['rational_coordinates']);A=coefficients([mon(e) for e in es],support,QQ)*C;assert A.denominator()==1 and matrix(GF(2),A).rank()==len(es);polys=[P({support[i]:c for i,c in col.dict().items()}) for col in A.columns()];basis[d]=polys;basis2[d]=[P2(f) for f in polys];coordinates[d]=C;exponents[d]=es
 print('reconstructed seed',seed['max_degree'],round(perf_counter()-start,2),flush=True)
 for d in range(seed['max_degree']+1,args.max_degree+1):
  es=list(powers(CLASSICAL,d));support=list(powers(WEIGHTS,d));idx={e:i for i,e in enumerate(es)};n=len(es);polys=[orig[0]*b for b in basis[d-1]];cols=[]
  for col in coordinates[d-1].columns():
   v=vector(QQ,n)
   for e,c in zip(exponents[d-1],col):v[idx[(e[0]+1,)+e[1:]]]=c
   cols.append(v)
  for e in es:
   if not e[0]:polys.append(mon(e));cols.append(vector(QQ,n,{idx[e]:1}))
  A=coefficients(polys,support,ZZ);C=matrix(QQ,cols).transpose();assert A.ncols()==n;stages=0
  while True:
   K=matrix(GF(2),A).right_kernel().basis_matrix().echelon_form();piv=list(K.pivots())
   if not piv:break
   B=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in B.list());B=matrix(ZZ,B.nrows(),B.ncols(),[c//2 for c in B.list()]);D=C*matrix(QQ,K.transpose())/2
   for j,k in enumerate(piv):A.set_column(k,B.column(j));C.set_column(k,D.column(j))
   stages+=1
  polys=[P({support[i]:c for i,c in col.dict().items()}) for col in A.columns()];basis[d]=polys;basis2[d]=[P2(f) for f in polys];coordinates[d]=C;exponents[d]=es
  products=[P2(g)*b for g,k in gs for b in basis2[d-k]];old=coefficients(products,support,GF(2));piv=old.augment(matrix(GF(2),A)).pivots();assert len(piv)==n;new=[j-old.ncols() for j in piv if j>=old.ncols()]
  for j in new:gs.append((polys[j],d))
  rows.append(dict(degree=d,rank=n,coefficient_rows=len(support),saturation_rounds=stages,classical_exponents=[list(e) for e in es],rational_coordinates=[[str(c) for c in row] for row in C],new_generator_indices=new,new_generator_polynomials=[str(polys[j]) for j in new]));print('extended M',d,'new',len(new),'degrees',[k for f,k in gs],round(perf_counter()-start,2),flush=True)
 stem='e7-d5a1-sources-'+str(args.max_degree);cert=dict(original_generators=list(map(str,orig)),normalization_factors=list(map(str,factors)),degrees=rows);archive=ROOT/('experiments/'+stem+'-certificates.json.gz');archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));report=dict(result='passed',max_degree=args.max_degree,coefficient_weights=list(WEIGHTS),classical_degrees=list(CLASSICAL),generator_degrees=[d for f,d in gs],generators=[[str(f),d] for f,d in gs],degrees=[{k:v for k,v in row.items() if k not in ('classical_exponents','rational_coordinates','new_generator_polynomials')} for row in rows],seconds=round(perf_counter()-start,3),scope='Full root invariant components through the stated bound only. No all-degree source or image claim.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7_d5a1_sources.py',seedpath)},raw_sha256={seedarchive.relative_to(ROOT).as_posix():sha256(seedarchive.read_bytes()).hexdigest()});(ROOT/('experiments/'+stem+'.json')).write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
