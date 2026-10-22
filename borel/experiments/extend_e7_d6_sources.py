"""Full N components to degree 33 and full M components to degree 27.

The seed consists of exact polynomial/coordinate certificates, independently
rechecked by audit_e7_d6_sources.py. It is data, not a generator-completeness
assumption. Higher N components are proved complete by primitive products.
Higher M components are saturated with a complete lower-component basis.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import gzip,json
from sage.all import QQ,ZZ,GF,matrix,vector
from e7_d6_sources import context,powers,WEIGHTS,CLASSICAL
ROOT=Path(__file__).resolve().parents[1]

def main():
 start=perf_counter();seed=ROOT/'experiments/e7-d6-source-seed20.json.gz';cert=json.loads(gzip.decompress(seed.read_bytes()));reports={}
 for kind,limit in (('N',33),('M',27)):
  P,orig,divs=context(kind);P2=P.change_ring(GF(2));rows=cert[kind]['degrees'];basis={0:[P.one()]};coords={0:matrix(QQ,1,1,[1])};exponents={0:[(0,)*7]};gs=[];gc=[];basis2={0:[P2.one()]};zero=(0,)*7
  for row in rows:
   d=row['degree'];basis[d]=[P(f) for f in row['basis_polynomials']];basis2[d]=[P2(f) for f in basis[d]];coords[d]=matrix(QQ,row['rational_coordinates']);exponents[d]=[tuple(e) for e in row['classical_exponents']]
   for j in row['new_generator_indices']:gs.append((basis[d][j],d));gc.append({e:c for e,c in zip(exponents[d],coords[d].column(j)) if c})
  @lru_cache(None)
  def classical(e):
   if not any(e):return P.one()
   j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
   return orig[j]*classical(tuple(u))
  def coeff_matrix(polys,support,ring):
   indices={e:i for i,e in enumerate(support)}
   return matrix(ring,len(support),len(polys),{(indices[tuple(e)],j):c for j,f in enumerate(polys) for e,c in f.dict().items()})
  for d in range(21,limit+1):
   exps=list(powers(CLASSICAL,d));support=list(powers(WEIGHTS,d));indices={e:i for i,e in enumerate(exps)};r=len(exps);stages=0
   if kind=='N':
    products=[];descriptors=[]
    for i,(g,k) in enumerate(gs):
     for j,b in enumerate(basis2[d-k]):products.append(P2(g)*b);descriptors.append((i,j))
    pivots=coeff_matrix(products,support,GF(2)).pivots();assert len(pivots)==r
    polys=[];columns=[]
    for p in pivots:
     i,j=descriptors[p];g,k=gs[i];polys.append(g*basis[d-k][j]);col=vector(QQ,r)
     for a,c in gc[i].items():
      for b,v in zip(exponents[d-k],coords[d-k].column(j)):
       if v:col[indices[tuple(x+y for x,y in zip(a,b))]]+=c*v
     columns.append(col)
    C=matrix(QQ,columns).transpose();A=coeff_matrix(polys,support,ZZ);new=[]
   else:
    # Multiplication by the primitive invariant h identifies odd components
    # exactly. In even degree it supplies a complete lower-degree summand.
    polys=[orig[0]*b for b in basis[d-1]];columns=[]
    for col in coords[d-1].columns():
     v=vector(QQ,r)
     for e,c in zip(exponents[d-1],col):v[indices[(e[0]+1,)+e[1:]]]=c
     columns.append(v)
    for e in exps:
     if e[0]==0:polys.append(classical(e));columns.append(vector(QQ,r,{indices[e]:1}))
    C=matrix(QQ,columns).transpose();A=coeff_matrix(polys,support,ZZ);assert A.ncols()==r
    while True:
     K=matrix(GF(2),A).right_kernel().basis_matrix().echelon_form();piv=list(K.pivots())
     if not piv:break
     D=A*matrix(ZZ,K.transpose());assert all(c%2==0 for c in D.list());D=matrix(ZZ,D.nrows(),D.ncols(),[c//2 for c in D.list()]);V=C*matrix(QQ,K.transpose())/2
     for j,k in enumerate(piv):A.set_column(k,D.column(j));C.set_column(k,V.column(j))
     stages+=1
    polys=[P({support[i]:c for i,c in col.dict().items()}) for col in A.columns()]
    products=[P2(g)*b for g,k in gs for b in basis2[d-k]];old=coeff_matrix(products,support,GF(2));pivots=old.augment(matrix(GF(2),A)).pivots();new=[j-old.ncols() for j in pivots if j>=old.ncols()]
    for j in new:gs.append((polys[j],d));gc.append({e:c for e,c in zip(exps,C.column(j)) if c})
   assert matrix(GF(2),A).rank()==r
   basis[d]=polys;basis2[d]=[P2(f) for f in polys];coords[d]=C;exponents[d]=exps
   record=dict(degree=d,rank=r,coefficient_rows=len(support),saturation_rounds=stages,new_generator_indices=new,method='primitive products of certified generators' if kind=='N' else 'saturation extending the complete previous component')
   rows.append(dict(**record,classical_exponents=[list(e) for e in exps],rational_coordinates=[[str(c) for c in line] for line in C],basis_polynomials=list(map(str,polys))))
   print(kind,'source',d,'rank',r,'new',len(new),'seconds',round(perf_counter()-start,2),flush=True)
  cert[kind]['degrees']=rows;reports[kind]=dict(generator_degrees=[d for f,d in gs],generators=[[str(f),d] for f,d in gs],degrees=[{k:v for k,v in row.items() if k not in ('rational_coordinates','basis_polynomials','classical_exponents')} for row in rows])
 archive=ROOT/'experiments/e7-d6-sources-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0))
 report=dict(result='passed',scope='Complete full invariant-source components for N through 33 and M through 27. No unbounded algebra-generation assertion.',max_degree=33,source_degree_bounds=dict(N=33,M=27),coefficient_weights=list(WEIGHTS),classical_degrees=list(CLASSICAL),sources=reports,seconds=round(perf_counter()-start,3),certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),seed_certificate_archive=seed.relative_to(ROOT).as_posix(),seed_archive_sha256=sha256(seed.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7_d6_sources.py')})
 (ROOT/'experiments/e7-d6-sources.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED',report['seconds'],flush=True)

if __name__=='__main__':main()
