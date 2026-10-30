"""Independent full exterior-representation checks of the realization proof."""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import (QQ, CyclotomicField, PolynomialRing, matrix,
                      identity_matrix, block_diagonal_matrix, binomial)
from sage.env import SAGE_VERSION
from verify_type_a_character_degeneration import wedge_lie, compound

ROOT=Path(__file__).resolve().parents[1]
K=CyclotomicField(4)
ii=K.gen()


def spectrum(m, p):
    D=binomial(2*m-2,m-1)
    R=PolynomialRing(K,'z'); z=R.gen()
    E=[K(0)]+[D*ii**(-m)*(-1)**(m+h-1)*
               sum(binomial(h-1,k-1)*p[m-k-1] for k in range(1,h+1))
               for h in range(1,m)]
    c=K(1)
    while True:
        Q=1+(-z)**m+c*sum(E[h]*z**h for h in range(1,m))
        if Q(-1) and Q.gcd(Q.derivative()).degree()==0:
            break
        c+=1
    e=Q.list()
    # charpoly of the gamma block: x^m Q(-1/x).
    F=sum((-1)**h*e[h]*z**(m-h) for h in range(m+1))
    A=matrix(K,m)
    for j in range(m-1):
        A[j+1,j]=1
    for j in range(m):
        A[j,m-1]=-F[j]
    assert A.charpoly().list()==F.list()
    assert F.gcd(F.derivative()).degree()==0
    assert (A-identity_matrix(K,m)).det()!=0
    t=ii*block_diagonal_matrix(identity_matrix(K,m),A)
    assert t.det()==1
    return c,t,Q


def ambient_basis(n):
    basis=[]
    for j in range(n-1):
        H=matrix(QQ,n); H[j,j]=1; H[j+1,j+1]=-1; basis.append(H)
    for j in range(n):
        for k in range(n):
            if j!=k:
                basis.append(matrix(QQ,n,{(j,k):1}))
    return basis


def pair(A,B):
    return sum(value*A[j,i] for (i,j),value in B.dict().items())


def main():
    # Independent symbolic inversion, before numeric constructions.
    for m in range(2,15):
        R=PolynomialRing(QQ,names=['e%s'%h for h in range(1,m)])
        e=list(R.gens())
        S=[sum(binomial(h-1,k-1)*e[k-1] for k in range(1,h+1))
           for h in range(1,m)]
        recovered=[sum((-1)**(h-k)*binomial(h-1,k-1)*S[k-1]
                       for k in range(1,h+1)) for h in range(1,m)]
        assert recovered==e

    checks=[]; pairings=0
    for m in (2,3,4):
        basis=ambient_basis(2*m)
        lie=[wedge_lie(Y,m) for Y in basis]
        D=binomial(2*m-2,m-1)
        assert (lie[0]**2).trace()==D*(basis[0]**2).trace()
        polynomials=[('power_%s'%r,[K(int(j==r)) for j in range(1,m)])
                     for r in range(1,m)]
        polynomials += [('mixed',[K((-1)**j*(j+1)) for j in range(1,m)]),
                        ('zero',[K(0)]*(m-1))]
        for label,p in polynomials:
            c,t,Q=spectrum(m,p)
            N=matrix(K,m,{(j,k):j+2*k+1 for j in range(m) for k in range(j+1,m)})
            u=block_diagonal_matrix(identity_matrix(K,m)+N,identity_matrix(K,m))
            # Independent minors in the actual irreducible representation.
            difference=compound(t*u,m)-compound(t,m)
            target=c*sum((p[j-1]*N**j for j in range(1,m)),matrix(K,m))
            expected=block_diagonal_matrix(target,matrix(K,m))
            for Y,dY in zip(basis,lie):
                assert pair(difference,dY)==D*(expected*Y).trace(),(m,label,Y)
                pairings+=1
            if label.startswith('power_'):
                r=int(label.split('_')[1])
                assert target==c*N**r and target!=0
            if label=='zero':
                assert difference!=matrix(K,difference.nrows())
                assert expected==0
            checks.append(dict(m=m,ambient_n=2*m,label=label,scalar=str(c),
                               exterior_dimension=int(binomial(2*m,m)),
                               spectrum_polynomial=str(Q),semisimple=True,
                               centralizer='S(GL_%s x (C*)^%s)'%(m,m)))
            print('Full ambient trace projection passed:',m,label,flush=True)

    # A fully symbolic nilpotent block, including a mixed polynomial.
    m=3; p=[K(2),K(-3)]; c,t,Q=spectrum(m,p)
    R=PolynomialRing(K,names=['a','b','d']); a,b,d=R.gens()
    N=matrix(R,[[0,a,d],[0,0,b],[0,0,0]])
    u=block_diagonal_matrix(identity_matrix(R,m)+N,identity_matrix(R,m))
    difference=compound(t.change_ring(R)*u,m)-compound(t,m).change_ring(R)
    expected=block_diagonal_matrix(c*(2*N-3*N**2),matrix(R,m))
    D=binomial(4,2)
    symbolic_pairings=0
    for Y in ambient_basis(6):
        assert pair(difference,wedge_lie(Y,3))==D*(expected*Y).trace()
        symbolic_pairings+=1

    paths=['IRREDUCIBLE_CENTRALIZER_REALIZATION.md',
           'tools/verify_irreducible_centralizer_realization.py',
           'tools/verify_type_a_character_degeneration.py']
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                symbolic_pascal_inverse_ranks=list(range(2,15)),
                actual_irreducible_matrix_cases=checks,
                exact_ambient_tangent_pairings=pairings,
                symbolic_SL6_tangent_pairings=symbolic_pairings,
                group_determinants_and_exact_centralizers_checked=True,
                zero_map_does_not_mean_constant_representation_matrix=True,
                uses_character_or_borel_runtime=False,
                scope='Finite exact companion/compound/wedge matrix corroboration of the separate arbitrary-rank proof. Shared Sage arithmetic and earlier wedge helpers; not an external audit or priority certificate.',
                source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'irreducible-centralizer-realization-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','actual_irreducible_matrix_cases')},indent=2))


if __name__=='__main__':
    main()
