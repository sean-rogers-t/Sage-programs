"""Full D4 source and its relative quotient over the proved D5 source."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
from sage.all import GF,PolynomialRing,TermOrder,prod
ROOT=Path(__file__).resolve().parents[1]
def main():
 start=perf_counter();data=json.loads((ROOT/'experiments/e7-d4-source.json').read_text());ints=[('h',1),('y',1)]+[(f,r['degree']) for r in data['degrees'] for f in r['new_generators']];ds=[d for f,d in ints];print('source degrees',ds,flush=True)
 P=PolynomialRing(GF(2),6,names=('h','y','z','v','c2','c3'));gs=[(P(f),d) for f,d in ints]
 names=('z','v','b','h','y')+tuple('g'+str(i) for i in range(2,len(gs)))
 S=PolynomialRing(GF(2),len(names),names=names,order=TermOrder('wdegrevlex',(1,1,3))+TermOrder('wdegrevlex',ds));z,v,b=S.gens()[:3];ss=S.gens()[3:];h,y,a=ss[:3]
 c2=a+h*z+v*v+v*y;c3=b+c2*v+h*h*z+h*y*z+h*z*z;phi=P.hom([h,y,z,v,c2,c3],S)
 I=S.ideal([ss[i]-phi(f) for i,(f,d) in enumerate(gs) if i>=3]);G=I.groebner_basis(algorithm='singular:slimgb');pure=[f for f in G if all(not any(e[:3]) for e in f.exponents())]
 print('graph',len(G),'pure',len(pure),'seconds',perf_counter()-start,flush=True)
 T=PolynomialRing(GF(2),len(ds),names=names[3:],order=TermOrder('wdegrevlex',ds));rho=S.hom([T.zero()]*3+list(T.gens()),T);J=T.ideal([rho(f) for f in pure]);H=J.hilbert_series(grading=ds);q=H.parent().gen();complete=H==1/prod(1-q**d for d in (1,1,2,4,4,6));print('SOURCE COMPLETE',complete,flush=True)
 # D5 signed-permutation subgroup acts on all five u's; D4 fixes u5=y.
 D=PolynomialRing(GF(2),6,names=('h','z','c2','c3','c4','c5'));es=[S.one(),2*z-v-y,c2,c3]
 full=[sum(es[j]*(1 if k-j==0 else v+y if k-j==1 else v*y) for j in range(4) if 0<=k-j<=2) for k in range(6)]
 substitute=D.hom([h,z]+full[2:],S);previous=json.loads((ROOT/'experiments/e7-d5-source-algebra.json').read_text());embedded=[]
 for f,d in zip(previous['generators'],previous['generator_degrees']):
  p=substitute(D(f)).reduce(G);assert all(not any(e[:3]) for e in p.exponents());embedded.append((rho(p),d))
 quotient=J+T.ideal([f for f,d in embedded]);QG=quotient.groebner_basis();dim=quotient.vector_space_dimension();relH=quotient.hilbert_series(grading=ds)
 print('RELATIVE QUOTIENT',dim,relH,'seconds',perf_counter()-start,flush=True)
 report=dict(result='passed' if complete else 'incomplete',generator_degrees=ds,integral_generators=[f for f,d in ints],generators=[str(f) for f,d in gs],elimination_variables=list(names),elimination_grading=[1,1,3]+ds,elimination_relations=[str(f) for f in I.gens()],elimination_groebner=[str(f) for f in G],kernel=[str(rho(f)) for f in pure],hilbert_series=str(H),relative_D5_generators=[[str(f),d] for f,d in embedded],relative_quotient_groebner=[str(f) for f in QG],relative_quotient_dimension=int(dim),relative_hilbert_series=str(relH),seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7-d4-source.json',ROOT/'experiments/e7-d5-source-algebra.json')})
 (ROOT/'experiments/e7-d4-source-algebra.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
