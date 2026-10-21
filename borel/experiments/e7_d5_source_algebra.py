"""Certify the mod-two source-algebra Hilbert series by elimination."""
from pathlib import Path
import json,time
from hashlib import sha256
from sage.all import GF,PolynomialRing,TermOrder,prod
ROOT=Path(__file__).resolve().parents[1]
def main():
 start=time.time();data=json.loads((ROOT/'experiments/e7-d5-source.json').read_text());P=PolynomialRing(GF(2),6,names=('h','z','c2','c3','c4','c5'))
 gs=[(P('h'),1)]+[(P(f),row['degree']) for row in data['degrees'] for f in row['new_generators']]
 ds=[d for f,d in gs];print('degrees',ds,flush=True)
 names=('z','b','c4','h')+tuple('b'+str(i) for i in range(1,len(gs)))
 S=PolynomialRing(GF(2),len(names),names=names,order=TermOrder('wdegrevlex',(1,3,4))+TermOrder('wdegrevlex',ds));z,b,c4=S.gens()[:3];vv=S.gens()[3:];h=vv[0];a=vv[1];t=vv[3]
 c2=a+h*z;c3=b+h*h*z+h*z*z
 partial=P.hom([h,z,c2,c3,c4,S.zero()],S);c5=t+partial(gs[3][0]);phi=P.hom([h,z,c2,c3,c4,c5],S)
 I=S.ideal([vv[i]-phi(f) for i,(f,d) in enumerate(gs) if i not in (0,1,3)]);print('elimination start',len(I.gens()),flush=True)
 G=I.groebner_basis(algorithm='singular:slimgb');out=[f for f in G if not any(any(e[:3]) for e in f.exponents())];print('GB',len(G),'kernel',len(out),'seconds',time.time()-start,flush=True)
 T=PolynomialRing(GF(2),len(gs),names=names[3:],order=TermOrder('wdegrevlex',ds));rho=S.hom([T.zero()]*3+list(T.gens()),T);J=T.ideal([rho(f) for f in out]);H=J.hilbert_series(grading=ds);print('HILBERT',H,flush=True)
 q=H.parent().gen();target=1/prod(1-q**d for d in (1,2,4,5,6,8));print('complete',H==target,flush=True)
 report=dict(result='passed' if H==target else 'incomplete',generator_degrees=ds,integral_generators=['h']+[f for row in data['degrees'] for f in row['new_generators']],generators=[str(f) for f,d in gs],kernel=[str(rho(f)) for f in out],hilbert_series=str(H),seconds=time.time()-start,elimination_variables=list(names),elimination_grading=[1,3,4]+ds,elimination_relations=[str(f) for f in I.gens()],elimination_groebner=[str(f) for f in G],source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7-d5-source.json')})
 (ROOT/'experiments/e7-d5-source-algebra.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
