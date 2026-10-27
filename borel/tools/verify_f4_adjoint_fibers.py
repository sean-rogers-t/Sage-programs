"""Full F4 root characters, five cubic reductions and all block boundaries."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ,PolynomialRing, matrix, vector,prod
from explore_f4_completion import restriction
from verify_affine_branching_classification import subsets
ROOT=Path(__file__).resolve().parents[1]

def main():
    rows=[restriction(I) for I in subsets(4)]
    assert len(rows)==15 and sum(r['degree']<=2 for r in rows)==10
    checked=[]
    for row in rows:
        if row['degree']<=2:continue
        n=len(row['nodes']);s=4-n
        P=PolynomialRing(QQ,['x'+str(i) for i in range(n)]+['t'+str(i) for i in range(s)])
        F=P.fraction_field();xs=P.gens()[:n];ts=P.gens()[n:]
        f=sum(v['coefficient']*prod(x**e for x,e in zip(xs,v['exponent']))*
              prod(F(t)**e for t,e in zip(ts,v['charge'])) for v in row['terms'])
        if n==2:
            x,y=xs;p,q=ts
            u,v={(1,3):(p,p*q),(1,4):(p*q*q,p*q),(2,4):(p*q*q,q)}[tuple(row['nodes'])]
            A=u+1/F(u);B=v+1/F(v)
            ell=A*y*y+B*y+A*(B*B-3)
            assert f==x*x+x*ell+(B*B-1)*y*y+A*B*y+A*A-B*B
            reduced=f.subs({x:-ell/2})
            J=2*reduced.derivative(y)
            assert J==-2*A*A*y**3-3*A*B*y*y+(-2*A*A*B*B+6*A*A+3*B*B-4)*y+A*B*(5-B*B)
            # Direct Laurent numerator check includes all trace ramification.
            num=P(J.numerator());common=num.gcd(num.derivative(y))
            assert len(common.dict())==1 and all(e==0 for e in common.exponents()[0][:n])
            checked.append(dict(nodes=row['nodes'],case='binary',full_formula=True,squarefree=True))
        elif row['nodes']==[1,2,4]:
            a,b,x=xs;t=ts[0];S=t**3+1/F(t**3)
            fa=f.derivative(a)-b;fb=f.derivative(b)-a;graph={a:-fb,b:-fa}
            reduced=f.subs(graph)
            M=4*x**3+3*S*x*x+2*(S*S-4)*x-S
            assert reduced.derivative(x)==-M
            num=P(M.numerator());assert num.gcd(num.derivative(x)).is_constant()
            checked.append(dict(nodes=row['nodes'],case='paired_A2',full_formula=True,squarefree=True))
        else:
            x,y,z=xs;t=ts[0];S=t**3+1/F(t**3)
            H=matrix(F,[[2*(x/t+t*t),1],[1,2*(x*t+1/F(t*t))]])
            ell=vector(F,[-x*t-1/F(t*t),-x/t-t*t]);v=vector(F,[y,z])
            Q=x*x+S*x-1
            assert f==Q+(v*H*v)/2+ell*v
            D=H.det();assert D==4*x*x+4*S*x+3
            R=Q-(ell*H.adjugate()*ell)/(2*D)
            N=P(t**9*D**2*R.derivative(x))
            assert N.degree(x)==5 and N.coefficient({x:5})==32*t**9
            assert N.gcd(N.derivative(x)).is_constant()
            k=P(t**3*D)
            c1=P((H.adjugate()*ell)[0].numerator())
            # Nonzero resultant bounds the entire consistent base by finitely
            # many central values. No specialization is used for this bound.
            resultant=k.resultant(c1,x);assert resultant and resultant.degree(x)==0
            B=PolynomialRing(QQ,names='X,T,V');X,T,V=B.gens()
            hom=P.hom([X,B(0),B(0),T],B)
            c2=P((H.adjugate()*ell)[1].numerator())
            assert B.ideal([hom(k),hom(k.derivative(x)),hom(c1),hom(c2),T*V-1]).is_one()
            checked.append(dict(nodes=row['nodes'],case='quadratic_block',full_formula=True,
                                squarefree=True,consistent_boundary_base_finite=True,
                                double_determinant_consistency_empty=True))
    report=dict(result='passed',all_proper_F4_Levis_covered=True,proper_standard_Levis=15,
                quadratic_cases=10,cubic_cases=checked,characters=rows,
                generic_nonemptiness_classification_claim=False,literature_priority_claim=False)
    names=['F4_ADJOINT_FIBERS.md','tools/explore_f4_completion.py','tools/verify_f4_adjoint_fibers.py']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'f4-adjoint-fibers-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','characters')},indent=2))
if __name__=='__main__':main()
