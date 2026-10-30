"""Finite exact corroboration of the simple-group branching deductions.

Weyl weight multiplicities are independent of the characteristic-polynomial
identities. Symbolic derivatives check the block argument; they do not
replace the all-rank proof or test singular quotient points.
"""
from pathlib import Path
from hashlib import sha256
import json
from math import prod, comb, factorial
from itertools import permutations
from sage.all import (QQ, PolynomialRing, WeylCharacterRing, SymmetricFunctions,
                      matrix, identity_matrix, Partition, Partitions, CyclotomicField,
                      block_diagonal_matrix, diagonal_matrix, companion_matrix, vector)
from sage.env import SAGE_VERSION
from verify_type_a_character_degeneration import compound, wedge_lie

ROOT = Path(__file__).resolve().parents[1]


def elementary(values, k):
    if k < 0:
        return QQ(0)
    a = [QQ(1)] + [QQ(0)]*k
    for v in values:
        for j in range(k, 0, -1):
            a[j] += v*a[j-1]
    return a[k]


def weyl_checks():
    rows = []
    for family, ranks in [('B', range(2, 5)), ('C', range(2, 5)), ('D', range(4, 6))]:
        for n in ranks:
            W = WeylCharacterRing(family+str(n), style='coroots')
            top = n if family == 'C' else n-1 if family == 'B' else n-2
            for k in range(1, top+1):
                hw = tuple(int(j == k-1) for j in range(n))
                weights = W(*hw).weight_multiplicities()
                for seed in (2, 5):
                    t = [QQ(seed+j) for j in range(n)]
                    actual = sum(mult * prod(t[j]**int(w[j]) for j in range(n))
                                 for w, mult in weights.items())
                    alphabet = t + [1/x for x in t] + ([QQ(1)] if family == 'B' else [])
                    expected = elementary(alphabet, k)
                    if family == 'C':
                        expected -= elementary(alphabet, k-2)
                    assert actual == expected, (family, n, k, seed)
                rows.append(dict(type=family, rank=n, fundamental_node=k,
                                 representation_dimension=int(sum(weights.values())), points=2))
    return rows


def block_checks():
    rows = []
    for family in ('A', 'B', 'C', 'D'):
        for m in range(1, 9):
            top = m-1 if family in ('A', 'B') else m if family == 'C' else m-2
            if top < 1:
                continue
            R = PolynomialRing(QQ, names=tuple('a'+str(j) for j in range(1, m+1))
                               + tuple('b'+str(j) for j in range(1, 9)))
            aa = [R(1)] + list(R.gens()[:m]); bb = [R(1)] + list(R.gens()[m:])
            Z = PolynomialRing(R, 'z'); z = Z.gen()
            if family == 'A':
                A = sum(aa[j]*z**j for j in range(m+1))
            else:
                degree = 2*m + (family == 'B')
                A = sum(aa[min(j, degree-j)]*z**j for j in range(degree+1))
            E = sum(bb[j]*z**j for j in range(len(bb)))
            C = A*E*(1-z*z if family == 'C' else 1)
            for k in range(1, top+1):
                assert C[k].derivative(aa[k]) == 1, (family, m, k)
                assert all(A[j].derivative(aa[k]) == 0 for j in range(k))
                rows.append(dict(type=family, zero_or_single_block_rank=m, k=k, derivative=1))
    return rows


def hook_checks():
    SF = SymmetricFunctions(QQ); s = SF.schur(); e = SF.elementary()
    rows = []
    for d in range(1, 13):
        for a in range(1, d+1):
            b = d-a
            # Sage basis conversion, independent of the hook expansion proof.
            expr = e(s([a]+[1]*b))
            assert expr.monomial_coefficients().get(Partition([d]), 0) == (-1)**(a-1)
            assert all(max(part) < d or list(part) == [d] for part in expr.support())
            rows.append(dict(size=d, arm=a, leg=b, top_coefficient=(-1)**(a-1)))
    return rows


def endpoint_controls():
    # Original 3-dimensional Sym^2 matrices for SL2, not a quotient gradient.
    R = PolynomialRing(QQ, names=('p', 'q', 'r')); F = R.fraction_field()
    p, q, r = map(F, R.gens()); s = (1+q*r)/p
    g = matrix(F, [[p, q], [r, s]])
    sym = matrix(F, [[p*p, p*q, q*q], [2*p*r, p*s+q*r, 2*q*s], [r*r, r*s, s*s]])
    actions = [matrix(F, [[0, 1, 0], [0, 0, 2], [0, 0, 0]]),
               matrix(F, [[0, 0, 0], [2, 0, 0], [0, 1, 0]]),
               matrix(F, [[2, 0, 0], [0, 0, 0], [0, 0, -2]])]
    basis = [matrix(F, [[0, 1], [0, 0]]), matrix(F, [[0, 0], [1, 0]]),
             matrix(F, [[1, 0], [0, -1]])]
    a = g.trace()
    assert sym.trace() == a*a-1
    for act, X in zip(actions, basis):
        assert (sym*act).trace() == 2*a*(g*X).trace()
    # tr(g)=0 with determinant 1 is regular semisimple; Sym^2 theta vanishes.
    g0 = matrix(QQ, [[0, -1], [1, 0]])
    assert g0.det() == 1 and g0.trace() == 0 and g0 != identity_matrix(QQ, 2)
    # D2 zero block: Spin4 = SL2 x SL2, vector character u*v, not a free coordinate.
    R2 = PolynomialRing(QQ, names=('u', 'v')); u, v = R2.gens()
    assert R2.ideal([(u*v).derivative(x) for x in R2.gens()]) == R2.ideal([u, v])
    # An actual SO4 involution is regular: centralizer dimension 1+1=rank 2.
    d = matrix(QQ, 4, 4, lambda i, j: (1 if i < 2 else -1) if i == j else 0)
    skew = []
    for i in range(4):
        for j in range(i+1, 4):
            X = matrix(QQ, 4); X[i, j] = 1; X[j, i] = -1; skew.append(X)
    assert all((d*X).trace() == 0 for X in skew)
    commutators = matrix(QQ, [(d*X-X*d).list() for X in skew])
    assert 6-commutators.rank() == 2
    return dict(SL2_Sym2_trace_pairings=3, strict_hook_degree_endpoint=True,
                B1_vector_endpoint=True, D2_vector_regular_zero=True)


def middle_exterior_checks():
    rows = []
    for k in (2, 3):
        n = 2*k; I = QQ(comb(2*k-2, k-1)); K = CyclotomicField(4*k); zeta = K.gen()
        for partition, double in [(part, dbl) for part in Partitions(k)
                                  for dbl in (False, True)]:
            # D=2^k generically, D=i at the double root; lambda fixed by the quadratic.
            D = zeta**k if double else K(2**k)
            lam = -(D-1/D)/(2*I)
            exponent = (1 if double else 0) + (2 if k % 2 == 0 else 0)
            beta = [(1 if double else 2)*zeta**(exponent+4*j) for j in range(k)]
            assert prod(beta) == D
            assert all(elementary(beta, j) == 0 for j in range(1, k))
            offsets = [QQ(j+1) for j, size in enumerate(partition) for unused in range(size)]
            mean = sum(offsets)/k
            offsets = [x-mean for x in offsets]
            targets = [-lam+c for c in offsets]
            assert len(set(targets)) == len(partition) and all(x != lam for x in targets)
            M = matrix(K, k, k-1, lambda i, j:
                       beta[i]*elementary(beta[:i]+beta[i+1:], k-(j+1)-1))
            a = M.solve_right(vector(K, [I*c for c in offsets]))
            R = PolynomialRing(K, 'u'); u = R.gen()
            coeff = [K(1)] + list(a) + [1/D]
            poly = sum((-1)**j*coeff[j]*u**(k-j) for j in range(k+1))
            A = companion_matrix(poly, format='right')
            assert A.charpoly('u') == poly and A.det() == 1/D
            g = block_diagonal_matrix([A, diagonal_matrix(K, beta)])
            h = diagonal_matrix(K, [lam]*k+targets)
            assert g.det() == 1 and h.trace() == 0
            rep = compound(g, k); count = 0
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    X = matrix(K, n); X[i, j] = 1
                    assert (rep*wedge_lie(X, k)).trace() == I*(h*X).trace()
                    count += 1
            for i in range(n-1):
                X = matrix(K, n); X[i, i] = 1; X[-1, -1] = -1
                action = wedge_lie(X, k)
                assert (rep*action).trace() == I*(h*X).trace()
                assert (action*action).trace() == I*(X*X).trace()
                count += 1
            # Regularity independently: cyclic A and no spectrum shared with beta.
            assert all(poly(b) != 0 for b in beta)
            rows.append(dict(k=k, remaining_target_multiplicities=list(map(int, partition)),
                             double_component=double, trace_pairings=count,
                             regular_source=True, index_checked=True))
    residual = []
    K = CyclotomicField(4); ii = K.gen()
    R = PolynomialRing(K, names=('D', 'lam')); D, lam = R.gens()
    for k in range(2, 15):
        I = comb(2*k-2, k-1); f = D*D+2*I*lam*D-1
        assert f.subs({lam: 0}) == (D-1)*(D+1)
        assert f.subs({lam: -ii/I}) == (D-ii)**2
        assert f.subs({lam: ii/I}) == (D+ii)**2
        residual.append(dict(k=k,
                             generic_multiplicity=1, exceptional_multiplicity=2))
        if k <= 7:
            # Independently enumerate unordered root assignments for every
            # partition in these ranks; do not store untested larger tables.
            counts = [dict(multiplicities=list(map(int, part)),
                           weighted_components=2*factorial(k)//prod(factorial(int(j)) for j in part))
                      for part in Partitions(k)]
            for record in counts:
                sizes = record['multiplicities']; orbits = set()
                for perm in permutations(range(k)):
                    start = 0; blocks = []
                    for size in sizes:
                        blocks.append(tuple(sorted(perm[start:start+size]))); start += size
                    orbits.add(tuple(blocks))
                assert 2*len(orbits) == record['weighted_components']
            residual[-1]['partition_counts'] = counts
            residual[-1]['root_assignments_independently_enumerated'] = True
    return rows, residual


def main():
    original, residual = middle_exterior_checks()
    report = dict(result='passed', date='2026-10-23', sage_version=SAGE_VERSION,
                  weyl_character_comparisons=weyl_checks(), block_derivatives=block_checks(),
                  hook_elementary_coefficients=hook_checks(), endpoint_controls=endpoint_controls(),
                  middle_exterior_original_matrices=original, middle_exterior_residuals=residual,
                  all_rank_proof_is_written_not_inferred_from_tests=True,
                  exceptional_Levi_classification_claimed=False)
    paths = [Path(__file__), ROOT/'SIMPLE_GROUP_CAYLEY_BRANCHING.md',
             ROOT/'tools/verify_type_a_character_degeneration.py']
    report['source_sha256'] = {p.relative_to(ROOT).as_posix():
        sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in paths}
    (ROOT/'simple-group-cayley-branching-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(result=report['result'],
                         weyl_representations=len(report['weyl_character_comparisons']),
                         symbolic_blocks=len(report['block_derivatives']),
                         hooks=len(report['hook_elementary_coefficients']),
                         middle_exterior_pairings=sum(x['trace_pairings'] for x in original),
                         controls=report['endpoint_controls']), indent=2))


if __name__ == '__main__':
    main()
