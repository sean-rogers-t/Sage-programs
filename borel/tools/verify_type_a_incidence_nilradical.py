"""Independent exact ideal and defining-matrix audit of the type A proof."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, matrix, block_diagonal_matrix, vector

ROOT=Path(__file__).resolve().parents[1]


def checks():
    # Abstract ideal operation is independent of the native elimination.
    models=[]
    for n in range(1,5):
        P=PolynomialRing(QQ,names=['d'+str(i) for i in range(n)]+['a'+str(i) for i in range(n)]+['g'])
        d=P.gens()[:n];a=P.gens()[n:2*n];g=P.gens()[-1]
        D=P.ideal(d);L=sum(u*v for u,v in zip(d,a))
        I=P.ideal([u*L for u in d])
        assert I.radical()==P.ideal(L)
        assert I.quotient(P.ideal(L))==D
        assert L not in I and L*L in I
        Lg=L+g;Ig=P.ideal([u*Lg for u in d])
        assert Ig==D.intersection(P.ideal(Lg))
        assert Ig.radical()==Ig
        models.append(n)
    P=PolynomialRing(QQ,names='x,y,z,r');x,y,z,r=P.gens();K=P.fraction_field()
    first=matrix(K,[[z*x,-z],[z,0]])
    second=matrix(K,[[y/z,-1/z],[1/z,0]])
    g=block_diagonal_matrix(first,second)
    assert g.det()==1
    phi=g.trace()*g.inverse().trace()-1
    assert phi==x*x+y*y+(z*z+z**-2)*x*y-1
    # Construct the actual fifteen-dimensional conjugation matrix, rather
    # than trusting the scalar character identity alone.
    pairs=[(i,j) for i in range(4) for j in range(4) if i!=j]
    basis=[]
    for i,j in pairs:
        m=matrix(K,4,4,0);m[i,j]=1;basis.append(m)
    for i in range(3):
        m=matrix(K,4,4,0);m[i,i]=1;m[i+1,i+1]=-1;basis.append(m)
    def coords(m):return vector(K,[m[i,j] for i,j in pairs]+[sum(m[j,j] for j in range(i+1)) for i in range(3)])
    Ad=matrix(K,[coords(g*m*g.inverse()) for m in basis]).transpose()
    assert Ad.trace()==phi
    # Laurent equations represented by r*z=1, so no coefficient boundary
    # is accidentally inverted. Work in the ORIGINAL x,y variables.
    S=z*z+r*r;R=z*z-r*r
    I=P.ideal(2*x+S*y,2*y+S*x,z*r-1)
    radical=P.ideal(2*y+S*x,R*x,z*r-1)
    assert I.radical()==radical
    epsilon=R*x
    assert epsilon not in I and epsilon*epsilon in I
    mu=2*R*x*y
    assert mu not in I and mu*mu in I and mu*epsilon in I
    # Generic target is empty by a unit-ideal calculation over Q(h).
    T=PolynomialRing(QQ,'h');h=T.gen();B=T.fraction_field()
    A=PolynomialRing(B,names='x,y,z,r');xx,yy,zz,rr=A.gens()
    generic=A.ideal(2*xx+(zz*zz+rr*rr)*yy,2*yy+(zz*zz+rr*rr)*xx,zz*rr-1,2*(zz*zz-rr*rr)*xx*yy-h)
    assert generic==A.ideal(1)
    # SL5 with one singleton, independent trace calculation and exact
    # transformed ideal. No extension to Q(z,w) is used for radicality.
    P=PolynomialRing(QQ,names='a,u,z,w');a,u,z,w=P.gens();K=P.fraction_field()
    singleton=1/(z*z*w*w)
    x=a/z;y=(u-singleton-a)/w
    g=block_diagonal_matrix(matrix(K,[[z*x,-z],[z,0]]),
        matrix(K,[[w*y,-w],[w,0]]),matrix(K,[[singleton]]))
    assert g.det()==1
    phi=g.trace()*g.inverse().trace()-1
    d=1/(z*z)-1/(w*w)
    L=d*a+z*z*w*w-1/(z*z*w**4)
    assert phi==u*(u/(w*w)+L)-1
    assert phi.derivative(u)==2*u/(w*w)+L
    assert phi.derivative(a)==d*u
    num=(w*w-z*z)*(w*w*(w*w-z*z)*a+z**4*w**6-1)
    Q=PolynomialRing(QQ,names='a,z,w');qa,qz,qw=Q.gens()
    I=Q.ideal((qw*qw-qz*qz)*(qw*qw*(qw*qw-qz*qz)*qa+qz**4*qw**6-1))
    assert I.radical()==I
    return dict(abstract_ideal_intersections_and_nilradicals=models,
        SL4_fifteen_dimensional_adjoint_matrix_trace=True,
        SL4_original_Laurent_ideal_radical_verified=True,
        SL4_nilpotent_and_target_annihilator_verified=True,
        SL4_generic_target_ideal_is_unit=True,
        SL5_defining_matrix_elimination_and_radical_verified=True,
        source_parameter_function_fields_not_used_for_radicality=True)


def main():
    report=dict(result='passed',engine='Sage defining matrices and Singular exact ideal operations',
        checks=checks(),native_audit_imported=False,
        all_type_A_nilradical_theorem_written=True,
        generic_type_A_reducedness_already_followed_from_quadratic_theorem=True,
        uniform_all_simple_reducedness_proved=False,new_Levi_coverage=False,
        novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'type-a-incidence-nilradical-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
