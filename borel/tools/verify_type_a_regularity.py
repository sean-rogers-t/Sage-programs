"""Exact finite corroboration of the irreducible SL_n regularity proof.

Uses only the Python standard library. Characters from Gelfand--Tsetlin
weights are compared with Jacobi--Trudi determinants and derivatives.
An algebraic witness is specified by its polynomial, not a floating root.
The arbitrary-rank argument is TYPE_A_REGULARITY.md, not these tests.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
from math import comb, factorial, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def trim(p):
    p = list(p)
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p or [0]


def add(p, q):
    return trim([(p[i] if i < len(p) else 0)+(q[i] if i < len(q) else 0)
                 for i in range(max(len(p), len(q)))])


def scale(p, c):
    return trim([c*x for x in p])


def mul(p, q):
    result = [0]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            result[i+j] += a*b
    return trim(result)


def evaluate(p, x):
    result = 0
    for c in reversed(p):
        result = result*x+c
    return result


def derivative(p):
    return trim([i*p[i] for i in range(1, len(p))])


def complete_polynomial(alphabet, degree):
    coefficients = [[1]]+[[0] for _ in range(degree)]
    for variable in alphabet:
        for j in range(1, degree+1):
            coefficients[j] = add(coefficients[j], mul(variable, coefficients[j-1]))
    return coefficients[degree]


def divide_one_minus_x(p):
    assert len(p) > 1
    q, previous = [], 0
    for c in p[:-1]:
        previous += c
        q.append(previous)
    assert p[-1] == -previous
    q = trim(q)
    assert mul(q, [1, -1]) == trim(p)
    return q


def dual_add(a, b):
    return add(a[0], b[0]), add(a[1], b[1])


def dual_mul(a, b):
    return mul(a[0], b[0]), add(mul(a[0], b[1]), mul(a[1], b[0]))


def dual_scale(a, c):
    return scale(a[0], c), scale(a[1], c)


def jt_character_and_derivative(p):
    """At z=(1,2,...,n-1,x); epsilon differentiates z_1 at 1.

    Uses h_r and a determinant, with no Gelfand--Tsetlin weights.
    """
    n = len(p)
    shape = tuple(x for x in p if x)
    length = len(shape)
    max_index = shape[0]+length-1
    h = [([1], [0])] + [([0], [0]) for _ in range(max_index)]
    alphabet = [([1], [1])] + [([j], [0]) for j in range(2, n)] + [([0, 1], [0])]
    for variable in alphabet:
        for r in range(1, max_index+1):
            h[r] = dual_add(h[r], dual_mul(variable, h[r-1]))
    result = ([0], [0])
    for permutation in permutations(range(length)):
        inversions = sum(permutation[i] > permutation[j] for i in range(length) for j in range(i+1, length))
        term = ([1], [0])
        for i, j in enumerate(permutation):
            index = shape[i]-i+j
            term = dual_mul(term, h[index] if index >= 0 else ([0], [0]))
        result = dual_add(result, dual_scale(term, (-1)**inversions))
    return result


@lru_cache(None)
def gt_weights(p):
    n = len(p)
    if n == 1:
        return {(p[0],): 1}
    result = defaultdict(int)
    for lower in product(*(range(p[i+1], p[i]+1) for i in range(n-1))):
        last = sum(p)-sum(lower)
        for weight, multiplicity in gt_weights(tuple(lower)).items():
            result[weight+(last,)] += multiplicity
    return dict(result)


def dimension(p):
    return prod(Q(p[i]-p[j]+j-i, j-i)
                for i in range(len(p)) for j in range(i+1, len(p)))


def solve(A, b):
    rows = [[Q(x) for x in row]+[Q(value)] for row, value in zip(A, b)]
    n = len(rows)
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        leading = rows[j][j]
        rows[j] = [a/leading for a in rows[j]]
        for i in range(n):
            if i != j:
                leading = rows[i][j]
                rows[i] = [a-leading*c for a, c in zip(rows[i], rows[j])]
    return [row[-1] for row in rows]


def divided_root(weights):
    result = defaultdict(int)
    numerator = {}
    for a, multiplicity in weights.items():
        reflected = (a[-1],)+a[1:-1]+(a[0],)
        assert weights[reflected] == multiplicity
        if a[0] != a[-1]:
            numerator[a] = multiplicity*(a[0]-a[-1])
        if a[0] > a[-1]:
            gap = a[0]-a[-1]
            for s in range(gap):
                monomial = (a[0]-1-s,)+a[1:-1]+(a[-1]+s,)
                result[monomial] += multiplicity*gap
    reconstructed = defaultdict(int)
    for a, c in result.items():
        reconstructed[(a[0]+1,)+a[1:]] += c
        reconstructed[a[:-1]+(a[-1]+1,)] -= c
    assert {a: c for a, c in reconstructed.items() if c} == numerator
    assert result and all(c > 0 for c in result.values())
    return dict(result)


def trace_projection_check(p, weights):
    n, D = len(p), sum(p)
    point = list(map(Q, range(2, n+1)))+[Q(1, factorial(n))]
    assert prod(point) == 1 and len(set(point)) == n
    weight_list = list(weights)
    multiplicities = [weights[a] for a in weight_list]
    coroot_weights = [[a[i]-a[i+1] for i in range(n-1)] for a in weight_list]
    eigenvalues = [prod(z**a for z, a in zip(point, weight)) for weight in weight_list]
    gram = [[sum(m*w[i]*w[j] for m, w in zip(multiplicities, coroot_weights))
             for j in range(n-1)] for i in range(n-1)]
    rhs = [sum(m*w[i]*ev for m, w, ev in zip(multiplicities, coroot_weights, eigenvalues))
           for i in range(n-1)]
    c = solve(gram, rhs)
    y = [sum(v*(int(i == j)-int(i == j+1)) for j, v in enumerate(c)) for i in range(n)]
    index = Q(sum(m*(a[0]-a[-1])**2 for a, m in weights.items()), 2)
    assert index > 0
    assert gram == [[index*(2*int(i == j)-int(abs(i-j) == 1))
                     for j in range(n-1)] for i in range(n-1)]
    character = sum(m*ev for m, ev in zip(multiplicities, eigenvalues))
    logarithmic = [sum(m*a[i]*ev for m, a, ev in zip(multiplicities, weight_list, eigenvalues))
                   for i in range(n)]
    assert y == [(v-Q(D, n)*character)/index for v in logarithmic]
    if p[0] == 1 and sum(p) in (1, n-1):
        assert len(set(y)) == n
    return index


def partitions(total, max_part=None):
    if not total:
        yield ()
        return
    if max_part is None:
        max_part = total
    for first in range(min(total, max_part), 0, -1):
        for rest in partitions(total-first, first):
            yield (first,)+rest


def main():
    cases = set()
    for n in range(2, 7):
        for D in range(1, 8):
            for shape in partitions(D):
                if len(shape) < n:
                    cases.add(shape+(0,)*(n-len(shape)))
        for d in range(8, 11):
            cases.add((d,)+(0,)*(n-1))
    witnesses = []
    endpoints, interior_fundamental, symmetric, central_quotient_cases = 0, 0, 0, 0
    for p in sorted(cases, key=lambda p: (len(p), sum(p), p)):
        n, D, L = len(p), sum(p), p[0]
        weights = gt_weights(p)
        assert weights[p] == 1 and all(sum(a) == D for a in weights)
        assert sum(weights.values()) == dimension(p)
        chi, numerator = [0], [0]
        for a, multiplicity in weights.items():
            value = multiplicity*prod(j**a[j-1] for j in range(1, n))
            chi = add(chi, [0]*a[-1]+[value])
            numerator = add(numerator, [0]*a[-1]+[value*(a[0]-a[-1])])
        jt_chi, jt_dz1 = jt_character_and_derivative(p)
        assert chi == jt_chi
        assert numerator == add(jt_dz1, scale(mul([0, 1], derivative(jt_chi)), -1))
        q = divide_one_minus_x(numerator)
        full_Q = divided_root(weights)
        specialized = [0]
        for a, coefficient in full_Q.items():
            value = coefficient*prod(j**a[j-1] for j in range(1, n))
            specialized = add(specialized, [0]*a[-1]+[value])
        assert q == specialized
        endpoint = L == 1 and D in (1, n-1)
        assert (len(full_Q) == 1) == endpoint
        index = trace_projection_check(p, weights)
        if len([a for a in p if a]) == 1:
            symmetric += 1
            assert index == comb(n+D, n+1)
            alphabet = [[j] for j in range(1, n)]+[[0, 1], [1], [0, 1]]
            assert q == complete_polynomial(alphabet, D-1)
        for central_order in range(2, n+1):
            if n % central_order == 0 and D % central_order == 0:
                assert not endpoint
                central_quotient_cases += 1
        if endpoint:
            endpoints += 1
        elif L == 1:
            interior_fundamental += 1
        else:
            # Exact algebraic witness: any root of q, then common scaling.
            # Positivity excludes 0 and all the other positive coordinates.
            assert len(q) >= L and q[0] > 0 and all(c >= 0 for c in q)
            assert all(evaluate(q, j) > 0 for j in range(1, n))
            assert numerator == mul([1, -1], q)
            witnesses.append(dict(n=n, partition=list(p), highest_weight=[p[i]-p[i+1] for i in range(n-1)],
                dimension=int(dimension(p)), trace_index=str(index),
                root_polynomial_coefficients_ascending=q,
                remaining_coordinates=list(range(1, n)),
                determinant_normalization='c^n = 1/((n-1)! x), z=c(1,2,...,n-1,x)',
                image_collision=[1, n]))
    selected = next(row for row in witnesses if row['partition'] == [2, 1, 0])
    assert selected['root_polynomial_coefficients_ascending'] == [8, 5]
    assert next(row for row in witnesses if row['partition'] == [2, 0, 0])['root_polynomial_coefficients_ascending'] == [4, 2]
    raw_point = (1, 2, -2)
    example_weights = gt_weights((2, 0, 0))
    chi_example = sum(m*prod(z**a for z, a in zip(raw_point, weight))
                      for weight, m in example_weights.items())
    raw_y = [(sum(m*weight[i]*prod(z**a for z, a in zip(raw_point, weight))
                   for weight, m in example_weights.items())-Q(2, 3)*chi_example)/5
             for i in range(3)]
    assert raw_y == [Q(-4, 15), Q(8, 15), Q(-4, 15)]
    paths = ['TYPE_A_REGULARITY.md', 'tools/verify_type_a_regularity.py',
             'EXTERIOR_REGULARITY.md', 'exterior-regularity-verification.json']
    report = dict(result='passed', date='2026-10-22', group='SL_n(C)',
        uniform_theorem='Nontrivial irreducible torus/global regularity preservation iff highest weight omega_1 or omega_(n-1)',
        proof='TYPE_A_REGULARITY.md, Sections 2-5; not extrapolation from finite checks',
        finite_n_range=[2, 6], all_normalized_partitions_through_size=7,
        additional_symmetric_powers=[8, 9, 10], character_cases=len(cases),
        independent_Jacobi_Trudi_character_and_derivative_checks=len(cases),
        exact_simple_coroot_trace_solves=len(cases),
        algebraic_witness_polynomials=len(witnesses), endpoint_cases=endpoints,
        interior_fundamental_cases=interior_fundamental,
        symmetric_power_cases=symmetric,
        symmetric_power_complete_factor_checks=symmetric,
        nontrivial_central_quotient_descent_cases=central_quotient_cases,
        SL3_symmetric_square_example=dict(raw_torus=list(raw_point),
            raw_image=list(map(str, raw_y)), determinant_scaling='c^3=-1/4, image scales by c^2'),
        inherited_exterior_proof='EXTERIOR_REGULARITY.md',
        root_existence='Nonconstant polynomial over QQ; q(0)>0 and q(x)>0 for x>0. No approximate roots are used.',
        witnesses=witnesses, runtime_API_changed=False, novelty_claim=False,
        source_sha256={p: sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in paths})
    (ROOT/'type-a-regularity-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('witnesses', 'source_sha256')}, indent=2))


if __name__ == '__main__':
    main()
