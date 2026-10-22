"""Full A5 root sources for even Levi blocks (4,2) and (2,2,2).

Both exact coefficient rings reduce to the same involution
(u,v,a,b,c) -> (u,v,-a,b-u*a,c-v*a). Its invariant algebra has a
five-generator module with one relation over R[u,v,a^2,N(b),N(c)].
The normal-form coefficient proof, also valid modulo two, establishes
unbounded completeness; no generator-degree cutoff is used.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from math import comb
import gzip,json
from sage.all import QQ,GF,PolynomialRing,TermOrder,matrix,prod
ROOT=Path(__file__).resolve().parents[1]

def main():
 start=perf_counter();sourcepath=ROOT/'experiments/a5-two-local-source.json';source=json.loads(sourcepath.read_text());allcases=[]
 for kind,dv,dc in (('42',3,4),('222',1,2)):
  weights=[1,dv,2,4,2*dc,2,dc,2+dv,2+dc];R=PolynomialRing(QQ,9,names=('u','v','A','Nb','Nc','B','C','D','E'),order=TermOrder('wdegrevlex',weights));u,v,A,Nb,Nc,B,C,D,E=R.gens();relations=[v*B-u*C-2*D,B*B-u*u*A-4*Nb,C*C-v*v*A-4*Nc,B*C-u*v*A-2*E,B*D-2*v*Nb+u*E,C*D-v*E+2*u*Nc,B*E-2*Nb*C+u*A*D,C*E-2*Nc*B-v*A*D,D*D-v*v*Nb-u*u*Nc+u*v*E,D*E-v*Nb*C+u*Nc*B+u*v*A*D,E*E-4*Nb*Nc-A*(v*v*Nb+u*u*Nc-u*v*E)]
  P=PolynomialRing(QQ,5,names=('u','v','a','b','c'));uu,vv,a,b,c=P.gens();invariants=[uu,vv,a*a,b*b-uu*a*b,c*c-vv*a*c,2*b-uu*a,2*c-vv*a,vv*b-uu*c,2*b*c-a*(vv*b+uu*c)];sigma=P.hom([uu,vv,-a,b-uu*a,c-vv*a],P);assert all(sigma(f)==f for f in invariants);phi=R.hom(invariants,P);assert all(not phi(f) for f in relations)
  # Exact free ambient module over the five invariant parameters.
  S=PolynomialRing(QQ,8,names=('a','b','c','u','v','A','Nb','Nc'),order=TermOrder('wdegrevlex',(1,2,dc))+TermOrder('wdegrevlex',(1,dv,2,4,2*dc)));a,b,c,uu,vv,AA,bb,cc=S.gens();ambient_relations=[a*a-AA,b*b-uu*a*b-bb,c*c-vv*a*c-cc];AG=S.ideal(ambient_relations).groebner_basis();assert {tuple(g.lm().exponents()[0]) for g in AG}=={(2,0,0,0,0,0,0,0),(0,2,0,0,0,0,0,0),(0,0,2,0,0,0,0,0)};sub=P.hom([uu,vv,a,b,c],S)
  if kind=='42':
   # s=-u/9, c1(first four)=2a-3s; v=c3-a*c2+a^3,
   # c=c4-a^2*c2/3. These are triangular 2-local coordinates.
   s=-uu/9;c1=2*a-3*s;c2=b;c3=vv+a*b-a**3;c4=c+a*a*b/3;full=[S.one(),c1+a,c2+a*c1,c3+a*c2,c4+a*c3,a*c4]
  else:
   # First two pairs have sums u+a,v+a and products b,c; the fifth
   # root coordinate is a. Multiply their characteristic polynomials.
   full=[S.one()]
   for coeff in ([S.one(),uu+a,b],[S.one(),vv+a,c],[S.one(),a]):
    nxt=[S.zero()]*(len(full)+len(coeff)-1)
    for i,f in enumerate(full):
     for j,g in enumerate(coeff):nxt[i+j]+=f*g
    full=nxt
  # Convert a full-S6 invariant to the five module generators using its
  # unique eight-term ambient normal form and an explicit mod-two syzygy.
  K=PolynomialRing(QQ,5,names=('u','v','A','Nb','Nc'));ku,kv,kA,kb,kc=K.gens();K2=K.change_ring(GF(2));toR=K.hom([u,v,A,Nb,Nc],R);mapped=[];old=PolynomialRing(QQ,5,names=('c1','c2','c3','c4','c5'));oldmap=old.hom(full[1:],S)
  for text,d in source['generators']:
   f=oldmap(old(text)).reduce(AG);coeff={e:K.zero() for e in [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]}
   for e,x in f.dict().items():coeff[tuple(e[:3])]+=K({tuple(e[3:]):x})
   f0=coeff[(0,0,0)];fE=coeff[(0,1,1)]/2;fb=coeff[(0,1,0)];fc=coeff[(0,0,1)];assert not coeff[(1,1,1)];assert fE.denominator()%2;quo,rem=K2(fb).quo_rem(K2.gen(1));assert not rem;fD=K({tuple(e):int(x) for e,x in quo.dict().items()});fB=(fb-kv*fD)/2;fC=(fc+ku*fD)/2;assert all(g.denominator()%2 for g in (f0,fB,fC,fD,fE));rr=toR(f0)+toR(fB)*B+toR(fC)*C+toR(fD)*D+toR(fE)*E;assert sub(phi(rr)).reduce(AG)==f;mapped.append(dict(degree=d,polynomial=str(rr)))
  R2=R.change_ring(GF(2));J=R2.ideal([R2(f) for f in relations]+[R2(R(row['polynomial'])) for row in mapped]);G=J.groebner_basis();assert J.dimension()==0;leads=[tuple(g.lm().exponents()[0]) for g in G];found={(0,)*9};todo=list(found)
  while todo:
   e=todo.pop()
   for j in range(9):
    w=tuple(a+int(i==j) for i,a in enumerate(e))
    if w not in found and not any(all(a<=b for a,b in zip(l,w)) for l in leads):found.add(w);todo.append(w)
  basis=sorted(found,key=lambda e:(sum(a*b for a,b in zip(e,weights)),e));row=dict(blocks=kind,source_complete_all_degrees=True,generator_names=list(R.variable_names()),generator_degrees=weights,coefficient_invariants=list(map(str,invariants)),relations=list(map(str,relations)),ambient_variables=list(S.variable_names()),ambient_grading=[1,2,dc,1,dv,2,4,2*dc],ambient_relations=list(map(str,ambient_relations)),full_A5_root_generators=mapped,module_quotient_groebner=list(map(str,G)),module_generator_exponents=[list(e) for e in basis],module_generator_degrees=[sum(a*b for a,b in zip(e,weights)) for e in basis],module_generators=len(basis),rational_module_rank=15 if kind=='42' else 90);allcases.append(row);print('A5 even blocks',kind,'full source 9 generators; module',len(basis),'max degree',max(row['module_generator_degrees']),round(perf_counter()-start,2),flush=True)
 archive=ROOT/'experiments/a5-even-block-sources-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(dict(cases=allcases),separators=(',',':')).encode(),mtime=0));report=dict(result='passed',full_sources_complete_all_degrees=True,cases=[{k:r[k] for k in ('blocks','generator_degrees','module_generators','rational_module_rank')} for r in allcases],seconds=round(perf_counter()-start,3),scope='Uniform involution presentation gives full integral sources; zero-dimensional mod-two source quotients give complete finite modules over the full A5 invariant source by graded Nakayama. Independent audit remains required.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),sourcepath)});(ROOT/'experiments/a5-even-block-sources.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
