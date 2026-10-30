"""Independent wedge/compound matrix checks of the written all-rank proof.

No runtime character or Borel API is imported. Finite corroboration only.
"""
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from sage.all import (QQ, PolynomialRing, CyclotomicField, matrix,
                      identity_matrix, binomial, factorial)
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def wedge_lie(X, p):
    """Derivation action on ordered exterior basis, by direct substitution."""
    subsets = list(combinations(range(X.nrows()), p))
    index = {v:i for i,v in enumerate(subsets)}
    A = matrix(X.base_ring(), len(subsets), sparse=True)
    for col, subset in enumerate(subsets):
        for pos, old in enumerate(subset):
            for new in range(X.nrows()):
                if not X[new, old]:
                    continue
                terms = list(subset)
                terms[pos] = new
                if len(set(terms)) != p:
                    continue
                inversions = sum(terms[i] > terms[j]
                                 for i in range(p) for j in range(i+1,p))
                A[index[tuple(sorted(terms))], col] += (-1)**inversions * X[new,old]
    return A


def compound(g, p):
    subsets = list(combinations(range(g.nrows()), p))
    return matrix(g.base_ring(), [[g.matrix_from_rows_and_columns(list(I),list(J)).det()
                                  for J in subsets] for I in subsets])


def index_formula(m,p,d):
    return sum((-1)**(l-1)*l**(d-1)*binomial(m,p-l) for l in range(1,p+1))


def coefficients(m, a, zeta):
    D = sum(a[p]*binomial(m-2,p-1) for p in range(1,m))
    return [(-1)**(j-1)*sum(a[p]*zeta**p*binomial(m-j-1,p-j)
                            for p in range(j,m))/D for j in range(1,m)]


def main():
    moments = {}
    count = 0
    for m in range(2,9):
        J = matrix(QQ,m,{(i,i+1):1 for i in range(m-1)})
        for p in range(1,m):
            A = wedge_lie(J,p)
            power = identity_matrix(QQ,A.nrows(),sparse=True)
            for k in range(1,m):
                power *= A
                Y = matrix(QQ,m,{(k,0):1})
                moment = (power*wedge_lie(Y,p)).trace()
                assert moment == index_formula(m,p,k+1), (m,p,k,moment)
                moments[m,p,k] = moment
                count += 1
    print('Independent wedge trace moments passed:',count,flush=True)

    central_checks = 0
    order_checks = []
    for m in range(2,9):
        K = CyclotomicField(m)
        zeta = K.gen()
        J = matrix(K,m,{(i,i+1):1 for i in range(m-1)})
        N = sum((J**k/factorial(k) for k in range(1,m)),matrix(K,m))
        # All central points, all single exterior powers, and one mixed input.
        inputs = [[int(p==q) for p in range(m+1)] for q in range(1,m)]
        inputs += [[p*p+1 for p in range(m+1)]]
        for h in range(m):
            for a in inputs:
                D = sum(a[p]*binomial(m-2,p-1) for p in range(1,m))
                actual = sum((sum(a[p]*zeta**(h*p)*moments[m,p,k]
                                   for p in range(1,m))/(D*factorial(k))*J**k
                              for k in range(1,m)),matrix(K,m))
                b = coefficients(m,a,zeta**h)
                expected = sum((b[k-1]*N**k for k in range(1,m)),matrix(K,m))
                assert actual == expected
                assert any(b)
                central_checks += 1
        if m%2 == 0:
            for s in range(m-1):
                a = [p**s for p in range(m+1)]
                b = coefficients(m,a,QQ(-1))
                D = sum(a[p]*binomial(m-2,p-1) for p in range(1,m))
                r = next(j for j in range(1,m) if b[j-1])
                assert r == m-1-s
                assert b[r-1] == (-1)**(s+1)*factorial(s)/D
                independent = [sum(a[p]*(-1)**p*moments[m,p,k]
                                   for p in range(1,m))/(D*factorial(k))
                               for k in range(1,m)]
                assert next(k for k,c in enumerate(independent,1) if c) == r
                assert independent[r-1] == b[r-1]
                order_checks.append(dict(m=m,s=s,r=r,leading=str(b[r-1])))
            # An arbitrary nonnegative list: interpolation, not chosen monomial.
            S = PolynomialRing(QQ,'x')
            a = [((p+2)**3 % 7)+1 for p in range(m+1)]
            P = S.lagrange_polynomial([(p,a[p]) for p in range(1,m)])
            b = coefficients(m,a,QQ(-1))
            r = next(j for j in range(1,m) if b[j-1])
            assert r == m-1-P.degree()
    print('Central comparisons passed:',central_checks,flush=True)

    # Full symbolic SL4 trace projection, all tangent directions, no index formula.
    R = PolynomialRing(QQ,names=['z','a','b','c','d','e','f'])
    z,a,b,c,d,e,f = R.gens()
    N = matrix(R,[[0,a,d,f],[0,0,b,e],[0,0,0,c],[0,0,0,0]])
    U = identity_matrix(R,4)+N
    basis = []
    for i in range(3):
        H = matrix(QQ,4); H[i,i]=1; H[i+1,i+1]=-1; basis.append(H)
    for i in range(4):
        for j in range(4):
            if i!=j:
                basis.append(matrix(QQ,4,{(i,j):1}))
    projections = 0
    for p in range(1,4):
        rep = z**p*compound(U,p)
        expected = sum(((-1)**(j-1)*z**p*binomial(3-j,p-j)*N**j
                        for j in range(1,p+1)),matrix(R,4))
        for Y in basis:
            assert (rep*wedge_lie(Y,p)).trace() == (expected*Y).trace()
            projections += 1
        H = basis[0]
        assert (wedge_lie(H,p)**2).trace() == binomial(2,p-1)*(H**2).trace()
    print('Symbolic SL4 tangent pairings passed:',projections,flush=True)

    # Explicit E tensor E for SL2: E=1 + V + 1, dimension 4.
    from sage.all import block_diagonal_matrix
    S = PolynomialRing(QQ,'u'); u = S.gen()
    g = matrix(S,[[-1,-u],[0,-1]])
    repE = block_diagonal_matrix(matrix(S,[[1]]),g,matrix(S,[[1]]))
    rep = repE.tensor_product(repE)
    basis2 = [matrix(QQ,[[1,0],[0,-1]]),matrix(QQ,[[0,1],[0,0]]),
              matrix(QQ,[[0,0],[1,0]])]
    for Y in basis2:
        dE = block_diagonal_matrix(matrix(QQ,[[0]]),Y,matrix(QQ,[[0]]))
        drep = dE.tensor_product(identity_matrix(QQ,4))+identity_matrix(QQ,4).tensor_product(dE)
        assert (rep*drep).trace() == 0
        assert (drep**2).trace() == 8*(Y**2).trace()

    paths = ['TYPE_A_CHARACTER_DEGENERATION.md','tools/verify_type_a_character_degeneration.py']
    report = dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                  exterior_matrix_trace_moments=count,matrix_sizes_m=list(range(2,9)),
                  exact_cyclotomic_central_checks=central_checks,
                  arbitrary_interpolation_checks=4,
                  all_order_witnesses=order_checks,
                  full_symbolic_SL4_tangent_pairings=projections,
                  explicit_tensor_square_zero_map_checked=True,
                  independent_of_character_and_borel_runtime=True,
                  scope='Finite exact matrix corroboration of separate written arbitrary-rank proofs; shared Sage arithmetic; no priority certification.',
                  source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'type-a-character-degeneration-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    main()
