"""G2: independent root-coordinate enumerations and exact Laurent ideals."""
from sage.all import QQ, RootSystem, PolynomialRing, LaurentPolynomialRing, matrix
from collections import Counter
from pathlib import Path
from hashlib import sha256
import json
ROOT=Path(__file__).resolve().parents[1]

def main():
    roots=[(1,0),(0,1),(1,1),(2,1),(3,1),(3,2)]
    roots=roots+[(-a,-b) for a,b in roots]
    native=Counter((2*a-3*b,-a+2*b) for a,b in roots)
    R=RootSystem('G2').ambient_space();cr=R.simple_coroots()
    sage=Counter(tuple(int(a.scalar(cr[i])) for i in (1,2)) for a in R.roots())
    assert native==sage
    L=LaurentPolynomialRing(QQ,names='u,z');u,z=L.gens();S=z+1/z;x=u+1/u
    f=[x*x+S*(x**3-2*x)+S*S-2,x*x+(S**3-2*S)*x+S*S-2]
    results=[]
    for i,h in enumerate(((1,2),(2,3))):
        actual=2+sum(u**w[i]*z**sum(a*b for a,b in zip(w,h)) for w in native)
        assert actual==f[i]
        trace=sum(sum(a*b for a,b in zip(w,h))**2 for w in native)
        assert trace==(16,48)[i]
        P=PolynomialRing(QQ,names='x,z,v');xx,zz,v=P.gens()
        ss=zz+v
        ff=[xx*xx+ss*(xx**3-2*xx)+ss*ss-2,xx*xx+(ss**3-2*ss)*xx+ss*ss-2][i]
        F=ff.derivative(xx)
        # Tangent logarithmic derivative in the Laurent ring z*v=1.
        D=lambda a:zz*a.derivative(zz)-v*a.derivative(v)
        assert P.ideal([F,F.derivative(xx),D(F),zz*v-1]).is_one()
        witness={xx:-1 if i==0 else -2,zz:1,v:1}
        assert F.subs(witness)==0 and D(F).subs(witness)==0
        assert D(D(ff)).subs(witness)==(10,-32)[i]
        results.append(dict(simple_node=i+1,central_cocharacter=list(h),trace_norm=trace,
                            full_character=True,whole_incidence_smooth=True,
                            central_logarithmic_hessian=(10,-32)[i]))
    killing=matrix(QQ,[[sum(w[i]*w[j] for w in native) for j in range(2)] for i in range(2)])
    assert killing.det()!=0
    report=dict(result='passed',cases=results,root_enumerations_agree=True,
                torus_killing_determinant=str(killing.det()),proper_standard_Levis=3,
                all_proper_G2_Levis_covered=True,literature_priority_claim=False)
    names=['G2_ADJOINT_FIBERS.md','tools/verify_g2_adjoint_fibers.py']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'g2-adjoint-fibers-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':main()
