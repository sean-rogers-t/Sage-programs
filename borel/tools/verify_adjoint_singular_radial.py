"""Independent Sage/Singular verification; no native audit functions imported."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, matrix, vector

ROOT=Path(__file__).resolve().parents[1]


def check():
    P=PolynomialRing(QQ,names='q,y,z,k,c');q,y,z,k,c=P.gens();K=P.fraction_field()
    S=K(z)+1/z;U=K(z)-1/z
    f=S*(1+y*q*q)-c/k
    D=lambda g:z*g.derivative(z)
    b=k*S/U
    assert b*D(f)==k*f+c
    assert f(z=1/z)==f
    E=U**6
    assert D(E)/E==b(k=6) and E(z=1/z)==E
    F=f.derivative(q);R=f.derivative(y);mu=D(f)
    # Clear only the explicitly inverted coefficients z,S,2. These are units
    # on z*(z^2-1)*(z^2+1) != 0, unlike the block coefficient y.
    def exact_polynomial(g):
        quotient,remainder=g.numerator().quo_rem(g.denominator())
        assert remainder==0
        return quotient
    I=P.ideal([exact_polynomial(F/(2*S)),exact_polynomial(R/S)])
    assert I==P.ideal(q*q,q*y)
    assert P.ideal(q).intersection(P.ideal(q*q,y))==I
    assert q not in I and q*q in I
    assert I.radical()==P.ideal(q)
    assert I.dimension()==4  # y,z,k,c are the support coordinates.
    assert mu-U==U*y*q*q
    assert f.derivative(q,q)(y=0)==0
    # On the bad support all derived critical-equation differentials vanish.
    J=matrix(K,[[g.derivative(v) for v in (q,y,z)] for g in (F,R)])
    assert J.apply_map(lambda g:g(q=0,y=0))==matrix(K,2,3,0)
    assert D(U)==S
    assert U*(k-D(b))==b*D(U)==k*S*S/U
    assert f(y=0).derivative(q)==0 and f.derivative(y)(y=0)==S*q*q
    # Target generic fiber: over QQ(h), impose z^2-hz-1 and the ACTUAL ideal.
    T=PolynomialRing(QQ,'h');h=T.gen();B=T.fraction_field()
    R=PolynomialRing(B,names='q,y,z');qq,yy,zz=R.gens()
    generic=R.ideal(qq*qq,qq*yy,zz*zz-h*zz-1)
    assert generic.dimension()==1
    assert qq not in generic and qq*qq in generic
    # The localization z*S*U removes none of this generic fiber.
    for unit in (zz,zz*zz-1,zz*zz+1):
        assert generic+R.ideal(unit)==R.ideal(1)
    singular_support=R.ideal(qq,yy,zz*zz-h*zz-1)
    assert singular_support.vector_space_dimension()==2
    # Independent polynomial jet calculation: all coefficients of H,L,Q
    # and V are algebraically independent at a point.
    P=PolynomialRing(QQ,names='q,y,z,h0,hy,hz,l0,ly,lz,Q0,Qy,Qz,a0,aq,ay,az,d0,dq,dy,dz,b0,bq,by,bz,k,c')
    (q,y,z,h0,hy,hz,l0,ly,lz,Q0,Qy,Qz,a0,aq,ay,az,d0,dq,dy,dz,b0,bq,by,bz,k,c)=P.gens()
    hh=h0+hy*y+hz*z;ll=l0+ly*y+lz*z;Q=Q0+Qy*y+Qz*z
    a=a0+aq*q+ay*y+az*z;d=d0+dq*q+dy*y+dz*z;b=b0+bq*q+by*y+bz*z
    f=hh*q*q/2+ll*q+Q;F=f.derivative(q);E=f.derivative(y);mu=z*f.derivative(z)
    defect=a*F+d*E+b*mu-k*f-c
    assert defect.derivative(q)==hh*a+F.derivative(y)*d+z*F.derivative(z)*b+aq*F+dq*E+bq*mu-k*F
    xs=(q,y);A=matrix(P,[[f.derivative(i,j) for j in xs] for i in xs])
    mixed=vector(P,[mu.derivative(x) for x in xs])
    da=matrix(P,[[g.derivative(x) for g in (a,d)] for x in xs])
    db=vector(P,[b.derivative(x) for x in xs])
    assert A*vector(P,[a,d])+b*mixed+da*vector(P,[F,E])+mu*db-k*vector(P,[F,E])==vector(P,[defect.derivative(x) for x in xs])
    return dict(radial_identity=True,block_and_full_derivative_identities=True,
        independent_polynomial_jet_check=True,primary_intersection_verified=True,
        actual_generic_fiber_dimension=1,nonzero_generic_fiber_nilpotent=True,
        generic_bad_support_length=2,localization_keeps_generic_fiber=True,
        differentiated_rhs_nonzero=True,normal_equations_must_be_retained=True)


def main():
    report=dict(result='passed',engine='Sage exact fraction fields and Singular ideal operations',
        checks=check(),native_audit_imported=False,
        abstract_radial_identity_suffices=False,uniform_adjoint_conjecture_proved=False,
        adjoint_counterexample=False,new_Levi_coverage=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-singular-radial-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
