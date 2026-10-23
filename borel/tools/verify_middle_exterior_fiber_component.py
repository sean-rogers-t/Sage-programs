"""Exact independent checks for the full-fiber component and parity proof."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import (QQ, CyclotomicField, PolynomialRing, LaurentPolynomialRing,
                      matrix, identity_matrix, block_diagonal_matrix, binomial)
from sage.env import SAGE_VERSION
from verify_irreducible_centralizer_realization import ambient_basis, pair
from verify_type_a_character_degeneration import wedge_lie, compound

ROOT=Path(__file__).resolve().parents[1]
K=CyclotomicField(4); alpha=K.gen()


def companion(coefficients):
    """Monic charpoly with elementary coefficients a_1,...,a_m."""
    m=len(coefficients); R=coefficients[0].parent()
    A=matrix(R,m)
    for j in range(m-1): A[j+1,j]=1
    for j in range(m): A[j,m-1]=-(-1)**(m-j)*coefficients[m-j-1]
    return A


def elementary(g):
    p=g.charpoly(); n=g.nrows()
    return [(-1)**j*p[n-j] for j in range(n+1)]


def raw(g,p):
    e=elementary(g)
    return sum(((-1)**(k-1)*e[p-k]*g**k for k in range(1,p+1)),
               matrix(g.base_ring(),g.nrows()))


def theta(g,m):
    F=raw(g,m); n=g.nrows()
    return (F-F.trace()/n*identity_matrix(g.base_ring(),n))/binomial(2*m-2,m-1)


def jacobian(g,m):
    n=g.nrows(); e=elementary(g); I=identity_matrix(K,n)
    powers=[I]+[g**k for k in range(1,m+1)]
    lower=[matrix(K,n)]+[raw(g,l) for l in range(1,m)]
    columns=[]
    for X in ambient_basis(n):
        dg=g*X
        dF=matrix(K,n)
        for k in range(1,m+1):
            de=(lower[m-k]*X).trace() if m-k else 0
            dp=sum((powers[r]*dg*powers[k-1-r] for r in range(k)),matrix(K,n))
            dF+=(-1)**(k-1)*(de*powers[k]+e[m-k]*dp)
        dF=(dF-dF.trace()/n*I)/binomial(2*m-2,m-1)
        columns.append(dF.list())
    return matrix(K,columns).transpose()


def main():
    ranks=[]; trace_checks=[]; curves=[]
    for m in range(2,7):
        n=2*m; D=binomial(n-2,m-1)
        B=companion([K(0)]*(m-1)+[(-1)**m*alpha**m])
        J=matrix(K,m,{(j,j+1):1 for j in range(m-1)})
        A=alpha*(identity_matrix(K,m)+J)
        g=block_diagonal_matrix(A,B)
        t=block_diagonal_matrix(alpha*identity_matrix(K,m),B)
        assert g.det()==t.det()==1 and theta(g,m)==theta(t,m)
        h=theta(t,m)
        commutator=matrix(K,[(h*X-X*h).list() for X in ambient_basis(n)]).transpose()
        centralizer_dim=n*n-1-commutator.rank()
        expected_K=m*m+m-1+(m-1 if m%2 else 0)
        assert centralizer_dim==expected_K
        derivative=jacobian(g,m)
        tangent=n*n-1-derivative.rank()
        dimension=m*m-1 if m%2 else m*m-m
        assert tangent==dimension+(1 if m%2 else 0),(m,tangent)
        ranks.append(dict(m=m,component_dimension=dimension,tangent_dimension=tangent,
                          differential_rank=derivative.rank(),target_centralizer_dimension=centralizer_dim))
        if m<=4:
            # Original trace construction in the actual irreducible exterior matrices.
            basis=ambient_basis(n); lie=[wedge_lie(X,m) for X in basis]
            rep=compound(g,m)
            rows=[]
            for dX in lie:
                product=rep*dX
                rows.append([pair(product,dY) for dY in lie])
            trace_matrix=matrix(K,rows).transpose()
            assert trace_matrix.rank()==derivative.rank()
            polynomial_pairings=matrix(K,[[D*(matrix(K,n,n,derivative.column(j))*Y).trace()
                                          for j in range(n*n-1)] for Y in basis])
            assert trace_matrix==polynomial_pairings
            trace_checks.append(dict(m=m,exterior_dimension=int(binomial(n,m)),
                                     pairings=(n*n-1)**2,rank=trace_matrix.rank()))
        print('Full differential checked:',ranks[-1],flush=True)

    for m in range(2,9):
        R=LaurentPolynomialRing(K,'s'); s=R.gen(); n=2*m
        coefficients=[binomial(m,l)*R(alpha)**l*s**(l-m) for l in range(1,m)]
        coefficients += [R(alpha)**m*s**(-m)]
        A=companion(coefficients)
        B=companion([R(0)]*(m-1)+[(-1)**m*R(alpha)**m])
        g=block_diagonal_matrix(A,s*B)
        assert g.det()==1
        base=matrix(K,n,n,[x(K(1)) for x in g.list()])
        epsilon=(-1)**(m+1)
        residual=1+epsilon*s**(2*m)-(1+epsilon)*s**m
        expected=R(alpha)**m/(2*binomial(2*m-2,m-1))*s**(-m)*residual*block_diagonal_matrix(
            identity_matrix(R,m),-identity_matrix(R,m))
        assert theta(g,m)-theta(base,m).change_ring(R)==expected
        assert residual(K(1))==0
        assert (residual.derivative()(K(1))==0)==bool(m%2)
        if m%2: assert residual.derivative().derivative()(K(1))==2*m*m
        else: assert residual.derivative()(K(1))==-2*m
        if m<=3:
            difference=compound(g,m)-compound(base,m).change_ring(R)
            for Y in ambient_basis(n):
                assert pair(difference,wedge_lie(Y,m))==binomial(2*m-2,m-1)*(expected*Y).trace()
        curves.append(dict(m=m,residual=str(residual),order_at_one=2 if m%2 else 1,
                           direct_exterior_trace_checked=m<=3))

    # Exact roots of unity verify all target collisions without numerical comparisons.
    collision_checks=[]
    for m in range(2,17):
        F=CyclotomicField(2*m); z=F.gen()
        gamma=[z**(2*j+1) for j in range(m)]
        w=[(1-x)**m for x in gamma]
        assert all(w)
        pairs=[]
        for j in range(m):
            for k in range(j+1,m):
                equal=w[j]==w[k]
                assert equal==(bool(m%2) and gamma[j]*gamma[k]==1)
                if equal: pairs.append([j,k])
        collision_checks.append(dict(m=m,equal_pairs=len(pairs)))

    paths=['MIDDLE_EXTERIOR_FIBER_COMPONENT.md','tools/verify_middle_exterior_fiber_component.py',
           'tools/verify_irreducible_centralizer_realization.py','tools/verify_type_a_character_degeneration.py']
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                differential_checks=ranks,independent_exterior_derivatives=trace_checks,
                exact_transverse_curve_checks=curves,cyclotomic_collision_checks=collision_checks,
                scheme_multiplicity_scope='The double multiplicity uses the written local elimination proof; tangent ranks alone do not prove it.',
                scope='Exact finite corroboration, original wedge traces versus polynomial differentiation, and exact curve identities. Separate arbitrary-rank local scheme proof; shared Sage arithmetic, no external audit or novelty claim.',
                source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'middle-exterior-fiber-component-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
