"""Exact root-string and Schur checks; no new Levi coverage or smoothness claim."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from collections import defaultdict
import json
import sympy as s
from audit_affine_branching_certificates import cartan, root_data, inverse

ROOT=Path(__file__).resolve().parents[1]
KINDS=('A3','B3','C3','D4','G2','F4','E6','E7','E8')


def matrix_for(kind):
    if kind[0] in 'EFG':return cartan(kind)
    n=int(kind[1:]); C=[[2*int(i==j) for j in range(n)] for i in range(n)]
    edges=[(i,i+1) for i in range(n-1)]
    if kind[0]=='D':edges=[(i,i+1) for i in range(n-2)]+[(n-3,n-1)]
    for i,j in edges:C[i][j]=C[j][i]=-1
    if kind[0]=='B':C[-1][-2]=-2
    if kind[0]=='C':C[-2][-1]=-2
    return tuple(map(tuple,C))


def root_identity(kind):
    C=matrix_for(kind);n=len(C);metric,positive=root_data(C)
    roots=tuple(positive)+tuple(tuple(-x for x in a) for a in positive)
    def dot0(a,b):return sum(a[i]*metric[i][j]*b[j] for i in range(n) for j in range(n))
    scale=F(2)/max(dot0(a,a) for a in roots)
    dot=lambda a,b:scale*dot0(a,b)
    ai=inverse(C)
    coordinates=lambda a:tuple(sum(x*y for x,y in zip(row,a)) for row in ai)
    highest=max(positive,key=lambda a:sum(coordinates(a)))
    rho=tuple(sum(a[i] for a in positive)/F(2) for i in range(n))
    eigen=dot(highest,highest)+2*dot(highest,rho)
    assert eigen==2*{'A3':4,'B3':5,'C3':4,'D4':6,'G2':4,'F4':9,'E6':12,'E7':18,'E8':30}[kind]
    total=defaultdict(F)
    # Pair beta with s_alpha(beta), then divide by 1-e^(-alpha)
    # by a finite geometric series, before summing positive roots.
    for alpha in positive:
        for beta in roots:
            p=dot(beta,alpha)
            if p<=0:continue
            length=2*p/dot(alpha,alpha);assert length.denominator==1
            for k in range(int(length)):
                for shift in (k,k+1):
                    weight=tuple(b-shift*a for a,b in zip(alpha,beta))
                    total[weight]+=p
    total={a:q for a,q in total.items() if q}
    expected={a:eigen-dot(a,a) for a in roots};expected[(0,)*n]=eigen*n
    assert total==expected
    report=dict(root_count=len(roots),casimir_eigenvalue=str(eigen),
        root_length_squares=sorted({str(dot(a,a)) for a in roots}),
        coefficientwise_radial_identity=True)
    if kind=='E6':
        central=(2,3,4,6,4,2);I=(0,1,2,4,5)
        added={a for a in roots if tuple(a[i] for i in I) in ((0,1,0,0,0),(0,-1,0,0,0))
               and abs(sum(x*y for x,y in zip(a,central)))==3}
        assert len(added)==4
        monomial=lambda a:F(2)**int(sum(coordinates(a)))
        def residual(extra):
            m={a:1+int(extra and a in added) for a in roots}
            flat=sum(dot(a,a)*m[a]*monomial(a) for a in roots)
            drift=F(0)
            for alpha in positive:
                e=monomial(alpha)
                drift+=(e+1)/(e-1)*sum(dot(a,alpha)*m[a]*monomial(a) for a in roots)
            return flat+drift-eigen*(n+sum(m[a]*monomial(a) for a in roots))
        assert residual(False)==0
        defect=residual(True);assert defect!=0
        report['altered_Levi_radial_defect']=str(defect)
    return report


def schur_checks():
    torus,trace=s.symbols('torus trace')
    assert s.cancel((torus**2+1)/(torus**2-1)*(torus-1/torus)-(torus+1/torus))==0
    assert s.expand(trace*s.diff(trace**2-1,trace)-2*(trace**2-1)-2)==0
    x,u,v,y,t,r,a,b=s.symbols('x u v y t r a b')
    source=ROOT/'adjoint-nilpotent-source-screen.json'
    f=s.sympify(json.loads(source.read_text())['potential'].replace('^','**'))
    f=s.expand(f.subs({u:a/t,v:t*b},simultaneous=True))
    H=s.hessian(f,(a,b));L=s.Matrix([s.diff(f,q).subs({a:0,b:0}) for q in (a,b)])
    Q=f.subs({a:0,b:0});D=s.expand(H.det());N=s.expand(Q*D-(L.T*H.adjugate()*L)[0]/2)
    point={t:2,r:3,x:1,y:0}
    wedge=s.cancel((s.diff(N,x)*s.diff(D,y)-s.diff(N,y)*s.diff(D,x)).subs(point)/D.subs(point))
    assert D.subs(point)==-s.Rational(9,4) and N.subs(point)==s.Rational(27,4)
    assert wedge==s.Rational(515,6)
    scalar_wedge=(s.diff(Q,x)*s.diff(D,y)-s.diff(Q,y)*s.diff(D,x)).subs(point)
    assert scalar_wedge==-s.Rational(525,2)
    pole={t:2,r:3,x:-s.Rational(1,2),y:0}
    assert D.subs(pole)==0 and N.subs(pole)==12 and s.diff(D,x).subs(pole)==6
    assert H.subs(pole).rank()==1 and H.subs(pole).row_join(L.subs(pole)).rank()==2
    # Frechet derivative of the Schur scalar, with all coefficients independent.
    h,j,k,l,m,q=s.symbols('h j k l m q')
    dh,dj,dk,dl,dm,dq=s.symbols('dh dj dk dl dm dq')
    M=s.Matrix([[h,j],[j,k]]);ell=s.Matrix([l,m]);star=-M.inv()*ell
    scalar=q-(ell.T*M.inv()*ell)[0]/2
    derivative=sum(s.diff(scalar,z)*dz for z,dz in zip((h,j,k,l,m,q),(dh,dj,dk,dl,dm,dq)))
    expected=dq+(star.T*s.Matrix([dl,dm]))[0]+(star.T*s.Matrix([[dh,dj],[dj,dk]])*star)[0]/2
    assert s.cancel(derivative-expected)==0
    return dict(E8_Q_D_wedge=str(scalar_wedge),E8_R_D_wedge=str(wedge),
        E8_regular_block_point={'t':2,'r':3,'x':1,'y':0,'D':'-9/4','N':'27/4'},
        E8_pole_point={'t':2,'r':3,'x':'-1/2','y':0,'D':0,'N':12,'D_x':6},
        E8_pole_point_quadratic_equations_inconsistent=True,
        E8_R_D_algebraically_independent_over_central_function_field=True,
        universal_Schur_envelope_derivative_verified=True,
        SL2_radial_derivation_regular_on_quotient_discriminant=True,
        any_alternative_quadratic_block_excluded=False,
        main_conjecture_counterexample=False)


def main():
    records={kind:root_identity(kind) for kind in KINDS}
    report=dict(result='passed',engine='Python fractions, root strings and SymPy',
        root_identities=records,schur=schur_checks(),
        radial_character_identity_is_classical=True,
        general_singular_image_bound_proved=False,new_Levi_coverage=False,
        proof_scope='Written all-type identity and Schur transport; finite checks are corroboration, not a new smoothness theorem.',
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'tools/audit_affine_branching_certificates.py',
                      ROOT/'adjoint-nilpotent-source-screen.json',ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-radial-schur-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
