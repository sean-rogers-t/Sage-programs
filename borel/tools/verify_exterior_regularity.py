"""Independent exact finite checks of the uniform exterior-power proof.

No Sage or runtime API is used. Projection is reconstructed from the
actual representation eigenvalues and coroot trace equations, rather
than the closed formula in the note. The all-rank argument is the note.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, factorial, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def elementary(values, k):
    if k < 0 or k > len(values):
        return Q(0)
    coefficients = [Q(1)] + [Q(0)] * k
    for value in values:
        for j in range(k, 0, -1):
            coefficients[j] += value * coefficients[j-1]
    return coefficients[k]


def solve(matrix, vector):
    rows = [[Q(x) for x in row] + [Q(y)] for row, y in zip(matrix, vector)]
    n = len(rows)
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [x / scale for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [a-scale*b for a, b in zip(rows[i], rows[j])]
    return [row[-1] for row in rows]


def projection_from_weights(z, k):
    n = len(z)
    subsets = list(combinations(range(n), k))
    # A simple coroot has diagonal (0,...,1,-1,...,0).
    weights = [[int(j in S)-int(j+1 in S) for j in range(n-1)] for S in subsets]
    eigenvalues = [prod(z[i] for i in S) for S in subsets]
    gram = [[sum(w[i]*w[j] for w in weights) for j in range(n-1)] for i in range(n-1)]
    rhs = [sum(w[j]*value for w, value in zip(weights, eigenvalues)) for j in range(n-1)]
    coroot_coordinates = solve(gram, rhs)
    output = [sum(c*(int(i == j)-int(i == j+1))
                  for j, c in enumerate(coroot_coordinates)) for i in range(n)]
    scale = comb(n-2, k-1)
    assert gram == [[scale*(2*int(i == j)-int(abs(i-j) == 1))
                     for j in range(n-1)] for i in range(n-1)]
    return output


def check_point(z, k):
    n = len(z)
    assert all(z) and prod(z) == 1 and len(set(z)) == n
    output = projection_from_weights(z, k)
    scale = comb(n-2, k-1)
    expected = [(z[i]*elementary(z[:i]+z[i+1:], k-1)
                 - Q(k, n)*elementary(z, k))/scale for i in range(n)]
    assert output == expected and sum(output) == 0
    for i, j in combinations(range(n), 2):
        rest = [z[a] for a in range(n) if a not in (i, j)]
        assert output[i]-output[j] == (z[i]-z[j])*elementary(rest, k-1)/scale
    if k in (1, n-1):
        assert len(set(output)) == n
    return output


def witness(n, k):
    m, r = n-2, k-1
    u = list(map(Q, range(1, m)))
    c = elementary(u, r)/elementary(u, r-1)
    F = factorial(m-1)
    a = m+1 + 1/(c*c*F)
    b = -1/(a*c*F)
    assert c > 0 and a > m-1 and -c < b < 0
    z = [a, b] + u + [-c]
    assert elementary(z[2:], r) == 0
    return z


def monomial_cancellation(n, k):
    """Compare polynomial coefficient dictionaries, without numeric inputs."""
    left, right = {}, {}
    for S in combinations(range(n), k):
        coefficient = int(0 in S)-int(1 in S)
        if coefficient:
            exponent = tuple(int(i in S) for i in range(n))
            left[exponent] = coefficient
    for S in combinations(range(2, n), k-1):
        for first, coefficient in ((0, 1), (1, -1)):
            exponent = tuple(int(i == first or i in S) for i in range(n))
            right[exponent] = right.get(exponent, 0)+coefficient
    assert left == right


def multiply(A, B):
    return [[sum(a*b for a, b in zip(row, column)) for column in zip(*B)] for row in A]


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def transpose(A):
    return list(map(list, zip(*A)))


def inverse(A):
    n = len(A)
    columns = [solve(A, [int(i == j) for i in range(n)]) for j in range(n)]
    return transpose(columns)


def determinant(A):
    n = len(A)
    total = Q(0)
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        total += (-1)**inversions*prod(A[i][p[i]] for i in range(n))
    return total


def compound(A, k):
    subsets = list(combinations(range(len(A)), k))
    return [[determinant([[A[i][j] for j in J] for i in I])
             for J in subsets] for I in subsets]


def exterior_derivative(X, k):
    n = len(X)
    subsets = list(combinations(range(n), k))
    positions = {S: i for i, S in enumerate(subsets)}
    result = [[Q(0) for _ in subsets] for _ in subsets]
    # Differentiate each factor in an ordered wedge basis.
    for column, S in enumerate(subsets):
        for position, j in enumerate(S):
            for i in range(n):
                if not X[i][j]:
                    continue
                word = list(S)
                word[position] = i
                if len(set(word)) != k:
                    continue
                inversions = sum(word[a] > word[b] for a in range(k) for b in range(a+1, k))
                result[positions[tuple(sorted(word))]][column] += (-1)**inversions*X[i][j]
    return result


def full_matrix_endpoint_checks():
    checked = 0
    for n in (3, 4):
        for diagonal in ([Q(1)]*n, [Q(2)]*(n-1)+[Q(1, 2**(n-1))]):
            g = [[diagonal[i] if i == j else Q(int(j == i+1)) for j in range(n)] for i in range(n)]
            assert determinant(g) == 1
            inv = inverse(g)
            # The full Jordan chain makes g cyclic, hence regular; these
            # matrices are nonsemisimple (the second has a repeated block).
            e_last = [[int(i == n-1)] for i in range(n)]
            columns = []
            v = e_last
            for _ in range(n):
                columns.append([row[0] for row in v])
                v = multiply(g, v)
            assert determinant(transpose(columns)) != 0
            for k in sorted({1, n-1}):
                H = g if k == 1 else [[-v for v in row] for row in inv]
                theta = [[H[i][j]-Q(trace(H), n)*int(i == j) for j in range(n)] for i in range(n)]
                rho_g = compound(g, k)
                drho_theta = exterior_derivative(theta, k)
                # Every off-diagonal matrix unit plus all simple coroots.
                basis = [[[Q(int(a == i and b == j)) for b in range(n)] for a in range(n)]
                         for i in range(n) for j in range(n) if i != j]
                basis += [[[Q(int(a == b == i)-int(a == b == i+1)) for b in range(n)]
                           for a in range(n)] for i in range(n-1)]
                for X in basis:
                    D = exterior_derivative(X, k)
                    assert trace(multiply(drho_theta, D)) == trace(multiply(rho_g, D))
                assert multiply(g, theta) == multiply(theta, g)
                checked += 1
    return checked


def main():
    generic_count, symbolic_count = 0, 0
    witnesses = []
    for n in range(2, 13):
        generic = list(map(Q, range(2, n+1))) + [Q(1, factorial(n))]
        for k in range(1, n):
            check_point(generic, k)
            generic_count += 1
            monomial_cancellation(n, k)
            symbolic_count += 1
            if 2 <= k <= n-2:
                z = witness(n, k)
                output = check_point(z, k)
                assert output[0] == output[1]
                if (n, k) == (4, 2):
                    assert z == [4, Q(-1, 4), 1, -1]
                    assert output == [0, 0, Q(15, 8), Q(-15, 8)]
                witnesses.append(dict(n=n, k=k, torus=list(map(str, z)),
                                      image_diagonal=list(map(str, output))))
    endpoint_count = full_matrix_endpoint_checks()
    paths = ['EXTERIOR_REGULARITY.md', 'tools/verify_exterior_regularity.py']
    report = dict(result='passed', date='2026-10-22', coefficient_field='QQ',
        group='SL_n(C)', representation='exterior^k(C^n)',
        uniform_theorem='Torus and global regularity preservation iff k=1 or k=n-1',
        arbitrary_rank_proof='EXTERIOR_REGULARITY.md, Sections 1-3; finite checks do not replace it',
        finite_n_range=[2, 12], regular_test_point_projection_checks=generic_count,
        sparse_polynomial_pair_identities=symbolic_count,
        rational_interior_counterexamples=len(witnesses),
        nonsemisimple_full_matrix_trace_checks=endpoint_count,
        independent_method='Actual subset eigenvalues, simple-coroot Gram equations and rational Gaussian elimination; direct exterior matrices at endpoints',
        witnesses=witnesses, runtime_API_changed=False, novelty_claim=False,
        source_sha256={p: sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in paths})
    (ROOT/'exterior-regularity-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('witnesses', 'source_sha256')}, indent=2))


if __name__ == '__main__':
    main()
