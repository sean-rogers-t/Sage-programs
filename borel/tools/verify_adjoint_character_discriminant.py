"""Independent defining-matrix adjoint check, using Sage over Q(x,z)."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, matrix, vector, identity_matrix

ROOT=Path(__file__).resolve().parents[1]


def checks():
    R=PolynomialRing(QQ,names='x,z');x,z=R.gens();K=R.fraction_field()
    # SL2 companion matrix, so x=+/-2 is allowed. No diagonalizing cover
    # or divided torus derivative is used for the quotient calculation.
    g=matrix(K,[[z*x,-z,0],[z,0,0],[0,0,z**-2]])
    assert g.det()==1
    inverse=g.inverse()
    pairs=[(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]
    basis=[]
    for i,j in pairs:
        m=matrix(K,3,3,0);m[i,j]=1;basis.append(m)
    basis.extend([matrix(K,[[1,0,0],[0,-1,0],[0,0,0]]),
                  matrix(K,[[0,0,0],[0,1,0],[0,0,-1]])])
    def coords(m):
        assert m.trace()==0
        return vector(K,[m[i,j] for i,j in pairs]+[m[0,0],-m[2,2]])
    Ad=matrix(K,[coords(g*m*inverse) for m in basis]).transpose()
    Phi=Ad.trace()
    S=z**3+z**-3;U=z**3-z**-3
    assert Phi==x*x+S*x
    # Determinant on the external four-dimensional root-space module.
    ext=[1,3,4,5]
    V=PolynomialRing(K,'v');v=V.gen()
    Pext=(identity_matrix(V,4)-v*Ad.matrix_from_rows_and_columns(ext,ext)).det()
    P=(identity_matrix(V,8)-v*Ad).det()
    assert P==(1-v)**2*(1-v*(x*x-2)+v*v)*Pext
    assert Pext(1)==(x-S)**2
    c1=g.trace();c2=inverse.trace()
    J=matrix(K,[[c.derivative(x),z*c.derivative(z)] for c in (c1,c2)])
    assert J.det()**2==4*Pext(1)
    assert J.det()(x=2,z=2)==QQ(49)/4
    phi2=(Ad*Ad).trace()
    critical=-S/2
    remainder=phi2.derivative(x)(x=critical)
    assert remainder==-QQ(3)/2*S*U*U
    mu=z*Phi.derivative(z)
    assert mu(x=critical)==-QQ(3)/2*S*U
    assert remainder(z=2)==-QQ(773955)/1024
    assert mu(x=-QQ(65)/16,z=2)==-QQ(12285)/128
    # Work directly over Q(z): original critical equation and Adams
    # derivative generate the unit ideal. Adding the latter deletes the
    # whole generic critical incidence, not just a stray special point.
    T=PolynomialRing(QQ,'t');t=T.gen();L=T.fraction_field()
    A=PolynomialRing(L,'w');w=A.gen();ss=t**3+t**-3
    original=A.ideal(2*w+ss)
    moment_derivative=2*w*(2*w*w+ss*ss-6)
    assert moment_derivative not in original
    assert original+A.ideal(moment_derivative)==A.ideal(1)
    # Full spectrum's first FOUR power traces, independently from the
    # native root-product/Newton implementation.
    for m in range(1,5):
        assert (Ad**m).trace()==(g**m).trace()*(inverse**m).trace()-1
    return dict(A2_adjoint_constructed_on_eight_matrix_basis_elements=True,
        full_and_external_matrix_determinants_verified=True,
        Jacobian_square_identity=True,internal_discriminant_retained=True,
        independent_power_trace_checks=4,
        second_Adams_gradient_at_specialization='-773955/1024',
        generic_critical_plus_second_Adams_gradient_is_unit_ideal=True,
        criticality_propagates_to_Adams_characters=False)


def main():
    report=dict(result='passed',engine='Sage defining matrices, determinants and exact ideals',
        checks=checks(),native_audit_imported=False,
        uniform_adjoint_conjecture_proved=False,new_Levi_coverage=False,
        adjoint_counterexample=False,novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-character-discriminant-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
