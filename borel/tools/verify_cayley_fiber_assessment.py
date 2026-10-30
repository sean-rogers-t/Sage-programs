"""One exact obstruction check supporting the theorem/literature assessment.

This script does not certify the literature audit or replace the uniform
proofs. It uses the original trace representation, no Borel recurrence.
"""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import (QQ, PolynomialRing, NumberField, matrix, identity_matrix,
                      block_diagonal_matrix, vector)
from sage.env import SAGE_VERSION

ROOT=Path(__file__).resolve().parents[1]


def main():
    Q=PolynomialRing(QQ,'z'); z=Q.gen()
    K=NumberField(z*z+1,'ii'); ii=K.gen()
    R=PolynomialRing(K,names=['a','b','c','d','e','f'])
    a,b,c,d,e,f=R.gens()
    N=matrix(R,[[0,a,d,f],[0,0,b,e],[0,0,0,c],[0,0,0,0]])
    I=identity_matrix(R,4)
    U=I+N; Uinv=I-N+N**2-N**3
    g=ii*U; ginv=-ii*Uinv
    rho=block_diagonal_matrix(g,ginv.transpose())
    basis=[]
    for j in range(3):
        H=matrix(QQ,4);H[j,j]=1;H[j+1,j+1]=-1;basis.append(H)
    for j in range(4):
        for k in range(4):
            if j!=k:
                E=matrix(QQ,4);E[j,k]=1;basis.append(E)
    drho=[block_diagonal_matrix(X,-X.transpose()) for X in basis]
    gram=matrix(QQ,[[(X*Y).trace() for Y in drho] for X in drho])
    rhs=vector(R,[(rho*X).trace() for X in drho])
    coefficients=gram.change_ring(R.fraction_field()).solve_right(rhs)
    actual=sum((coeff*X for coeff,X in zip(coefficients,basis)),matrix(R,4))
    expected=ii/2*N**2*Uinv
    assert actual==expected==(g-ginv)/2-ii*I
    Y=N-N**2/2+3*N**3/8
    assert expected==ii/2*Y**2
    recovered=Y+Y**2/2+Y**3/8
    assert recovered==N
    inverse=N+N**2/2+N**3/8
    assert inverse-inverse**2/2+3*inverse**3/8==N
    assert all(p.total_degree()>=2 for p in expected.list() if p)
    # Work over QQ for the ideal decomposition, which holds after extension.
    S=PolynomialRing(QQ,names=['a','b','c','d','e','f'])
    a,b,c,d,e,f=S.gens()
    ideal=S.ideal(a*b,b*c,a*e+d*c)
    P=S.ideal(a,c); Qideal=S.ideal(b,a*e+d*c)
    assert ideal==P.intersection(Qideal)
    assert P.is_prime() and Qideal.is_prime()
    assert ideal.radical()==ideal
    assert P.dimension()==Qideal.dimension()==4
    assert (P+Qideal).dimension()==3
    # Direct trace moment degree two: J superdiagonal, test E31 pairs
    # to J^2, with B_rho=2tr, so choose E31/2 as the dual test element.
    J=matrix(QQ,4,{(j,j+1):1 for j in range(3)})
    test=matrix(QQ,4,{(2,0):QQ(1)/2})
    dJ=block_diagonal_matrix(J,-J.transpose())
    dtest=block_diagonal_matrix(test,-test.transpose())
    rhot=block_diagonal_matrix(ii*identity_matrix(K,4),-ii*identity_matrix(K,4))
    assert (rhot*dJ*dtest).trace()==0
    assert (rhot*dJ**2*dtest).trace()/2==ii/2
    paths=['CAYLEY_FIBER_ASSESSMENT.md','tools/verify_cayley_fiber_assessment.py']
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                actual_representation='SL4: standard direct sum dual; central t=iI',
                full_trace_dimension=8,lie_gram_dimension=15,
                symbolic_trace_projection_checked=True,
                square_map_coordinate_change_and_inverse_checked=True,
                zero_linear_nonzero_quadratic_moment_checked=True,
                zero_fiber_reduced=True,zero_fiber_components=2,
                component_dimensions=[4,4],intersection_dimension=3,
                uniform_proofs='Written in the assessment; finite calculations do not prove their all-rank scope.',
                literature_scope='Known translated maps, Springer isomorphism structure and square-zero geometry identified; no exhaustive search or novelty certificate.',
                source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'cayley-fiber-assessment-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    main()
