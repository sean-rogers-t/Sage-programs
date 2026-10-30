"""Independent exact checks of the all-rank type C fundamental obstruction.

Uses only the standard library. Kernel multiplicities come from actual
symplectic contraction matrices, not character subtraction. Algebraic
witnesses are verified in free quadratic algebras without approximate roots.
"""
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import comb, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def binomial(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def elementary(values, degree):
    if degree < 0 or degree > len(values):
        return Q(0)
    a = [Q(1)]+[Q(0)]*degree
    for x in values:
        for j in range(degree, 0, -1):
            a[j] += x*a[j-1]
    return a[degree]


def fundamental_x(r, s, values):
    assert len(values) == r and 0 <= s <= r
    return sum((binomial(r-s+2*j, j)-binomial(r-s+2*j, j-1))
               * elementary(values, s-2*j) for j in range(s//2+1))


def rank(matrix):
    if not matrix:
        return 0
    a = [list(map(Q, row)) for row in matrix]
    row = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        leading = a[row][column]
        for j in range(column, len(a[0])):
            a[row][j] /= leading
        for i in range(row+1, len(a)):
            leading = a[i][column]
            if leading:
                for j in range(column, len(a[0])):
                    a[i][j] -= leading*a[row][j]
        row += 1
        if row == len(a):
            break
    return row


def exterior_weight(S, n):
    return tuple(int(i in S)-int(n+i in S) for i in range(n))


def primitive_kernel_weights(n, k):
    domain = defaultdict(list)
    for S in combinations(range(2*n), k):
        domain[exterior_weight(S, n)].append(S)
    target = defaultdict(list)
    if k >= 2:
        for S in combinations(range(2*n), k-2):
            target[exterior_weight(S, n)].append(S)
    result, blocks, max_columns = {}, 0, 0
    for weight, columns in domain.items():
        rows = {S: j for j, S in enumerate(target.get(weight, []))}
        matrix = [[0]*len(columns) for _ in rows]
        for column, S in enumerate(columns):
            positions = {letter: i for i, letter in enumerate(S)}
            for i in range(n):
                if i in positions and n+i in positions:
                    a, b = positions[i], positions[n+i]
                    lower = tuple(letter for letter in S if letter not in (i, n+i))
                    matrix[rows[lower]][column] += (-1)**(a+b-1)
        multiplicity = len(columns)-rank(matrix)
        assert multiplicity >= 0
        if multiplicity:
            result[weight] = multiplicity
        blocks += 1
        max_columns = max(max_columns, len(columns))
    assert sum(result.values()) == binomial(2*n, k)-binomial(2*n, k-2)
    return result, blocks, max_columns


def determinant(matrix):
    a = [list(map(Q, row)) for row in matrix]
    answer, n = Q(1), len(a)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            answer = -answer
        leading = a[j][j]
        answer *= leading
        for i in range(j+1, n):
            c = a[i][j]/leading
            for b in range(j+1, n):
                a[i][b] -= c*a[j][b]
    return answer


def weyl_character(n, k, z):
    numerator = [[x**(n-i+int(i < k))-x**(-(n-i+int(i < k)))
                  for x in z] for i in range(n)]
    denominator = [[x**(n-i)-x**(-(n-i)) for x in z] for i in range(n)]
    return determinant(numerator)/determinant(denominator)


def regular(z):
    pair_sums = [x+1/x for x in z]
    return all(x not in (2, -2) for x in pair_sums) and len(set(pair_sums)) == len(z)


def trace_coordinates(weights, z, index):
    values = {w: prod(x**a for x, a in zip(z, w)) for w in weights}
    return [sum(m*w[i]*values[w] for w, m in weights.items())/(2*index)
            for i in range(len(z))]


def generic_check(n, k, weights):
    z = list(map(Q, range(2, n+2)))
    assert regular(z)
    index = binomial(2*n-2, k-1)-binomial(2*n-2, k-3)
    gram = [[sum(m*w[i]*w[j] for w, m in weights.items())
             for j in range(n)] for i in range(n)]
    assert gram == [[2*index*int(i == j) for j in range(n)] for i in range(n)]
    alphabet = [v for x in z for v in (x, 1/x)]
    character = sum(m*prod(x**a for x, a in zip(z, w)) for w, m in weights.items())
    assert character == elementary(alphabet, k)-elementary(alphabet, k-2)
    assert character == weyl_character(n, k, z)
    assert character == fundamental_x(n, k, [x+1/x for x in z])
    y = trace_coordinates(weights, z, index)
    for i in range(n):
        rest = [v for j, x in enumerate(z) if j != i for v in (x, 1/x)]
        factor = elementary(rest, k-1)-elementary(rest, k-3)
        assert y[i] == (z[i]-1/z[i])*factor/(2*index)
    for i, j in combinations(range(n), 2):
        rest = [v for a, x in enumerate(z) if a not in (i, j) for v in (x, 1/x)]
        # Coefficients of R(u)=(1-u^2) product of the other paired factors.
        R = lambda degree: elementary(rest, degree)-elementary(rest, degree-2)
        A, B = R(k-1)+R(k-3), R(k-2)
        a, b = z[i], z[j]
        assert 2*index*(y[i]-y[j]) == (a-b)*((1+1/(a*b))*A+2*(a+b)/(a*b)*B)
        assert 2*index*(y[i]+y[j]) == (a*b-1)*((a+b)/(a*b)*A+2*(a*b+1)/(a*b)*B)
    return index


def zero_pair_sums(r, s):
    trials = []
    for t in range(2, s+3):
        other = [Q(t*j) for j in range(2, r+1)]
        B = fundamental_x(r, s, [Q(0)]+other)
        A = fundamental_x(r, s, [Q(1)]+other)-B
        assert A > 0 and B >= 0
        trials.append(dict(t=t, A=str(A), B=str(B)))
        if B != 2*A:
            values = [-B/A]+other
            assert len(set(values)) == r
            assert all(x not in (2, -2, Q(5, 2)) for x in values)
            assert fundamental_x(r, s, values) == 0
            return values, trials
    raise AssertionError((r, s, 'The bounded construction failed'))


def direct_quadratic_trace(weights, pair_sums):
    """First torus coordinate is 2; the remaining ones are quadratic.

    Multiplication occurs in the free algebra with each w_i satisfying
    w_i^2-x_i*w_i+1=0. A weight has only exponents -1,0,1, so each
    inverse can be expanded as x_i-w_i without further products.
    """
    total = defaultdict(Q)
    for weight, multiplicity in weights.items():
        if not weight[0]:
            continue
        terms = {0: Q(multiplicity*weight[0])*Q(2)**weight[0]}
        for i, exponent in enumerate(weight[1:]):
            bit = 1 << i
            assert exponent in (-1, 0, 1)
            if exponent == 1:
                terms = {mask | bit: value for mask, value in terms.items()}
            elif exponent == -1:
                updated = defaultdict(Q)
                for mask, value in terms.items():
                    updated[mask] += pair_sums[i]*value
                    updated[mask | bit] -= value
                terms = dict(updated)
        for mask, value in terms.items():
            total[mask] += value
    assert all(value == 0 for value in total.values())


def symbolic_ballot_check(r):
    polynomial = {(0, 0): 1}
    for i in range(r):
        updated = defaultdict(int)
        for (degree, mask), coefficient in polynomial.items():
            updated[(degree, mask)] += coefficient
            updated[(degree+1, mask | (1 << i))] += coefficient
            updated[(degree+2, mask)] += coefficient
        polynomial = dict(updated)
    for s in range(1, r+1):
        for mask in range(1 << r):
            actual = polynomial.get((s, mask), 0)-polynomial.get((s-2, mask), 0)
            a = mask.bit_count()
            expected = 0
            if a <= s and (s-a) % 2 == 0:
                j = (s-a)//2
                expected = binomial(r-a, j)-binomial(r-a, j-1)
                assert expected > 0
            assert actual == expected


def main():
    records, total_blocks, max_block, short_roots = [], 0, 0, 0
    for n in range(1, 9):
        for k in range(1, n+1):
            weights, blocks, largest = primitive_kernel_weights(n, k)
            total_blocks += blocks
            max_block = max(max_block, largest)
            index = generic_check(n, k, weights)
            short_roots += n*(n-1)//2
            if n == 1:
                continue
            if k == 1:
                z = [Q(2), Q(-1, 2)]+list(map(Q, range(3, n+1)))
                assert regular(z)
                y = trace_coordinates(weights, z, index)
                assert y[0] == y[1] == Q(3, 4)
                records.append(dict(n=n, k=k, trace_index=index,
                    method='rational defining collision', torus=list(map(str, z)),
                    image=list(map(str, y)), root='epsilon_1-epsilon_2'))
            else:
                values, trials = zero_pair_sums(n-1, k-1)
                direct_quadratic_trace(weights, values)
                records.append(dict(n=n, k=k, trace_index=index,
                    method='exact primitive-kernel trace in a free quadratic algebra',
                    first_torus_entry='2', remaining_pair_sums=list(map(str, values)),
                    trials=trials, vanished_coordinate='y_1',
                    all_branches_regular=True, algebra_dimension=2**(n-1)))
    for r in range(1, 11):
        symbolic_ballot_check(r)
    auxiliary = []
    for r in range(1, 21):
        for s in range(1, r+1):
            values, trials = zero_pair_sums(r, s)
            auxiliary.append(dict(r=r, s=s, chosen_t=trials[-1]['t'], trials=len(trials)))
    # The Sp4/omega2 example uses i^2=-1, with the second coordinate
    # reconstructed independently from the actual kernel weights.
    weights, _, _ = primitive_kernel_weights(2, 2)
    imaginary_coefficient = Q(0)
    for (a, b), multiplicity in weights.items():
        if b:
            imaginary_coefficient += multiplicity*Q(2)**a*b*b
    assert imaginary_coefficient/4 == Q(5, 4)
    paths = ['TYPE_C_FUNDAMENTAL_REGULARITY.md',
             'tools/verify_type_c_fundamental_regularity.py']
    report = dict(result='passed', date='2026-10-22', group='Sp_(2n)(C)',
        representation='primitive exterior^k of the defining representation, 1<=k<=n',
        theorem='No fundamental map preserves torus/global regularity for n>=2; rank-one defining case is positive',
        all_rank_proof='TYPE_C_FUNDAMENTAL_REGULARITY.md, Sections 2-5',
        finite_n_range=[1, 8], primitive_kernel_cases=36,
        exact_contraction_weight_blocks=total_blocks, largest_weight_block_columns=max_block,
        independent_Weyl_alternant_characters=36, trace_Gram_and_coordinate_checks=36,
        short_root_pair_formula_checks=short_roots,
        rational_defining_counterexamples=7, quadratic_algebra_counterexamples=28,
        symbolic_ballot_characters=55, symbolic_ballot_maximum_rank=10,
        auxiliary_zero_lemma_cases=len(auxiliary), auxiliary_maximum_rank=20,
        maximum_bounded_search_trials=max(row['trials'] for row in auxiliary),
        auxiliary_trials=auxiliary, witnesses=records,
        Sp4_omega2_example=dict(pair_sums=['5/2', '0'], torus=['2', 'i'], image=['0', '5*i/4']),
        scope='Independent exact finite corroboration of the written arbitrary-rank proof, not a novelty certificate or classification of all type C irreducibles.',
        runtime_API_changed=False, novelty_claim=False,
        source_sha256={p: sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in paths})
    (ROOT/'type-c-fundamental-regularity-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('auxiliary_trials', 'witnesses', 'source_sha256')}, indent=2))


if __name__ == '__main__':
    main()
