"""Exact finite corroboration of the all-highest-weight type C proof.

Run with Sage Python. Weyl alternant/Chebyshev slices are compared with
Freudenthal weights, their rank-one strings and actual trace coordinates.
Symmetric monomials provide a third character calculation in one family.
"""
from collections import defaultdict
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json

from sage.all import QQ, PolynomialRing, QuadraticField, WeylCharacterRing, matrix, prod
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]
R = PolynomialRing(QQ, 'X')
X = R.gen()


def partitions(total, bound=None):
    if not total:
        yield ()
        return
    if bound is None:
        bound = total
    for first in range(min(bound, total), 0, -1):
        for tail in partitions(total-first, first):
            yield (first,)+tail


def sl2_characters(maximum):
    a = [R(1)]
    if maximum:
        a.append(X)
    while len(a) <= maximum:
        a.append(X*a[-1]-a[-2])
    return a


def character_slice(p, other_pair_sums):
    """One-variable type C Weyl alternant, without weight multiplicities."""
    n = len(p)
    assert len(other_pair_sums) == n-1
    P = sl2_characters(p[0]+n-1)
    exponents = [p[i]+n-i-1 for i in range(n)]
    numerator = R(0)
    for i, degree in enumerate(exponents):
        minor = matrix(QQ, [[P[exponents[j]](v) for v in other_pair_sums]
                            for j in range(n) if j != i])
        numerator += (-1)**i*minor.det()*P[degree]
    denominator = prod(X-v for v in other_pair_sums)*prod(
        other_pair_sums[i]-other_pair_sums[j]
        for i in range(n-1) for j in range(i+1, n-1))
    quotient, remainder = numerator.quo_rem(R(denominator))
    assert not remainder and quotient.degree() == p[0]
    return quotient


def weight_oracle(p, rings):
    n = len(p)
    if n not in rings:
        rings[n] = WeylCharacterRing('C'+str(n), style='coroots')
    dynkin = tuple(p[i]-(p[i+1] if i+1 < n else 0) for i in range(n))
    character = rings[n](dynkin)
    result = {}
    for w, multiplicity in character.weight_multiplicities().items():
        coordinates = tuple(w.to_vector())
        assert len(coordinates) == n and all(a == int(a) for a in coordinates)
        result[tuple(int(a) for a in coordinates)] = int(multiplicity)
    assert result[p] == 1 and sum(result.values()) == character.degree()
    return result


def rank_one_strings(weights, L):
    rest = {w[1:] for w in weights}
    modules = {}
    reconstructed = defaultdict(int)
    for a in range(L+1):
        module = {}
        for v in rest:
            m = weights.get((a,)+v, 0)-weights.get((a+2,)+v, 0)
            assert m >= 0
            if m:
                module[v] = m
                for first in range(-a, a+1, 2):
                    reconstructed[(first,)+v] += m
        modules[a] = module
    assert dict(reconstructed) == weights and modules[L]
    return modules


def branch_slice(modules, z_other, L):
    P = sl2_characters(L)
    result, coefficients = R(0), {}
    for a, module in modules.items():
        value = sum(QQ(m)*prod(QQ(z)**w for z, w in zip(z_other, v))
                    for v, m in module.items())
        assert value >= 0 and (value > 0) == bool(module)
        result += value*P[a]
        coefficients[a] = value
    assert coefficients[L] > 0
    return result, coefficients


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for tail in compositions(total-first, length-1):
                yield (first,)+tail


def symmetric_weights(n, degree):
    result = defaultdict(int)
    for a in compositions(degree, 2*n):
        result[tuple(a[i]-a[n+i] for i in range(n))] += 1
    return dict(result)


def symmetric_derivative(degree, pair_sums):
    maximum = degree-1
    P = sl2_characters(maximum)
    series = [R(1)]+[R(0)]*maximum
    for pair_sum in [X, X]+list(pair_sums):
        factor = [polynomial(pair_sum) for polynomial in P]
        series = [sum(series[j]*factor[d-j] for j in range(d+1))
                  for d in range(maximum+1)]
    return series[maximum]


def multiply_quadratic(a, b, q):
    f, g = a
    h, j = b
    return ((f*h-g*j) % q, (f*j+g*h+X*g*j) % q)


def exact_witness_trace(weights, z_other, q, L):
    """Direct trace in Q[X]/q plus a root of z^2-Xz+1, all exact."""
    powers = {0: (R(1), R(0))}
    for exponent in range(1, L+1):
        powers[exponent] = multiply_quadratic(powers[exponent-1], (R(0), R(1)), q)
        powers[-exponent] = multiply_quadratic(powers[-exponent+1], (X, R(-1)), q)
    total = [R(0), R(0)]
    for weight, multiplicity in weights.items():
        c = multiplicity*weight[0]*prod(QQ(z)**a for z, a in zip(z_other, weight[1:]))
        f, g = powers[weight[0]]
        total[0] += c*f
        total[1] += c*g
    assert not total[0] % q and not total[1] % q


def generic_trace_check(p, weights, modules):
    n, L = len(p), p[0]
    z = [QQ(i) for i in range(2, n+2)]
    pair_sums = [v+1/v for v in z]
    index = sum(QQ(m)*w[0]**2 for w, m in weights.items())/2
    assert index > 0
    gram = [[sum(m*w[i]*w[j] for w, m in weights.items())
             for j in range(n)] for i in range(n)]
    assert gram == [[2*index*int(i == j) for j in range(n)] for i in range(n)]
    for i in range(n):
        other = pair_sums[:i]+pair_sums[i+1:]
        f = character_slice(p, other)
        actual_character = sum(m*prod(v**a for v, a in zip(z, w)) for w, m in weights.items())
        assert actual_character == f(pair_sums[i])
        rhs = sum(m*w[i]*prod(v**a for v, a in zip(z, w)) for w, m in weights.items())
        assert rhs == (z[i]-1/z[i])*f.derivative()(pair_sums[i])
    f = character_slice(p, pair_sums[1:])
    branched, coefficients = branch_slice(modules, z[1:], L)
    assert f == branched
    negative = character_slice(p, [-v for v in pair_sums[1:]])
    assert negative(-X) == (-1)**sum(p)*f
    assert f.derivative()(2) > 0
    if p[1:] == (0,)*(n-1):
        assert index == comb(2*n+L, 2*n+1)
        assert f.derivative() == symmetric_derivative(L, pair_sums[1:])
    return index


def universal_witness_check(p, weights, index):
    n = len(p)
    K = QuadraticField(-1, 'i')
    imaginary = K.gen()
    z = [K(value) for a in range(2, n//2+2) for value in (a, -a)]
    if n % 2:
        z.append(imaginary)
    assert len(z) == n and n >= 2
    assert all(a*a != 1 for a in z)
    assert all(z[a] != z[b] and z[a]*z[b] != 1
               for a in range(n) for b in range(a+1, n))
    transformed = list(z)
    for a in range(0, n-1, 2):
        transformed[a], transformed[a+1] = z[a+1], z[a]
    if n % 2:
        transformed[-1] = 1/z[-1]
    assert transformed == [-a for a in z]
    y = [sum(m*w[a]*prod(value**exponent for value, exponent in zip(z, w))
             for w, m in weights.items())/(2*index) for a in range(n)]
    wy = list(y)
    for a in range(0, n-1, 2):
        wy[a], wy[a+1] = y[a+1], y[a]
    if n % 2:
        wy[-1] = -y[-1]
    central_scalar = (-1)**sum(p)
    assert all((sum(w)-sum(p)) % 2 == 0 for w in weights)
    assert wy == [central_scalar*a for a in y]
    assert y[0] == central_scalar*y[1]
    return dict(n=n, partition=list(p), central_scalar=central_scalar,
                torus=list(map(str, z)), image=list(map(str, y)),
                vanishing_root='epsilon_1-epsilon_2' if central_scalar == 1 else 'epsilon_1+epsilon_2')


def main():
    cases = set()
    for n in range(1, 5):
        for degree in range(1, 7):
            for shape in partitions(degree):
                if len(shape) <= n:
                    cases.add(shape+(0,)*(n-len(shape)))
    for degree in range(1, 7):
        cases.add((degree, 0, 0, 0, 0))
    cases.update({(2, 1, 0, 0, 0), (2, 1, 1, 0, 0), (2, 2, 0, 0, 0), (3, 1, 0, 0, 0)})
    rings, witnesses, universal_checks = {}, [], []
    symmetric_count, inherited_count, coordinate_checks, maximum_trials = 0, 0, 0, 0
    total_distinct_weights = 0
    for p in sorted(cases, key=lambda p: (len(p), sum(p), p)):
        n, L = len(p), p[0]
        weights = weight_oracle(p, rings)
        total_distinct_weights += len(weights)
        assert max(w[0] for w in weights) == L
        assert all(max(abs(v) for v in w) <= L for w in weights)
        modules = rank_one_strings(weights, L)
        index = generic_trace_check(p, weights, modules)
        coordinate_checks += n
        if n >= 2:
            universal_checks.append(universal_witness_check(p, weights, index))
        if p[1:] == (0,)*(n-1):
            symmetric_count += 1
            assert symmetric_weights(n, L) == weights
        if L == 1:
            inherited_count += 1
            continue
        grids = [range(2+j*(L+1), 2+(j+1)*(L+1)) for j in range(n-1)]
        trials = 0
        for z_other in product(*grids):
            trials += 1
            xs = [QQ(z)+1/QQ(z) for z in z_other]
            f = character_slice(p, xs)
            q = f.derivative()
            if q(-2):
                break
        else:
            raise AssertionError((p, 'nonzero polynomial vanished on the entire grid'))
        assert trials <= (L+1)**(n-1)
        maximum_trials = max(maximum_trials, trials)
        branched, coefficients = branch_slice(modules, z_other, L)
        assert f == branched
        assert q.degree() == L-1 and q(2) > 0 and q(-2)
        assert all(q(x) > 0 and x > 2 for x in xs)
        assert len(set(xs)) == len(xs)
        exact_witness_trace(weights, z_other, q, L)
        witnesses.append(dict(n=n, partition=list(p), dimension=sum(weights.values()),
            trace_index=str(index), other_torus_coordinates=list(z_other),
            other_pair_sums=list(map(str, xs)), trials=trials,
            character_slice_coefficients_ascending=list(map(str, f.list())),
            witness_polynomial_coefficients_ascending=list(map(str, q.list())),
            witness_polynomial_at_minus_two=str(q(-2)),
            witness_polynomial_at_two=str(q(2)),
            quadratic_algebra_dimension=2*int(q.degree()), image_collision='y_1=0'))
    example = next(w for w in witnesses if w['partition'] == [2, 0])
    assert example['other_torus_coordinates'] == [2]
    assert example['witness_polynomial_coefficients_ascending'] == ['5/2', '2']
    w2 = weight_oracle((2, 0), rings)
    # In Q[i sqrt(39)] the first pair sum is -5/4; the second trace
    # coordinate follows directly from the first-coordinate SL2 strings.
    # A direct Laurent reduction in the first coordinate confirms the first.
    q = X+QQ(5)/4
    exact_witness_trace(w2, (2,), q, 2)
    second = QQ(0)
    for (a, b), m in w2.items():
        # Sum inversion-paired terms via the recurrence z^a+z^-a.
        if a == 0:
            second += m*b*QQ(2)**b
        elif a > 0:
            assert w2.get((-a, b)) == m
            # General recurrence avoids any floating algebraic root.
            T0, T1 = QQ(2), QQ(-5)/4
            for _ in range(2, a+1):
                T0, T1 = T1, QQ(-5)/4*T1-T0
            second += m*b*QQ(2)**b*T1
    assert second/12 == QQ(15)/32
    same_parity_sums = []
    for n in range(2, 6):
        for degrees in ((1, 3), (2, 4)):
            combined = defaultdict(int)
            for degree in degrees:
                for w, m in weight_oracle((degree,)+(0,)*(n-1), rings).items():
                    combined[w] += m
            if degrees[0] == 2:
                combined[(0,)*n] += 3
            index = sum(QQ(m)*w[0]**2 for w, m in combined.items())/2
            row = universal_witness_check((degrees[-1],)+(0,)*(n-1), combined, index)
            row.pop('partition')
            row['representation'] = 'Sym^'+str(degrees[0])+' + Sym^'+str(degrees[1])+(' + 3 trivial' if degrees[0] == 2 else '')
            same_parity_sums.append(row)
    paths = ['TYPE_C_REGULARITY.md', 'tools/verify_type_c_regularity.py',
             'TYPE_C_FUNDAMENTAL_REGULARITY.md', 'type-c-fundamental-regularity-verification.json']
    report = dict(result='passed', date='2026-10-22', sage_version=SAGE_VERSION,
        group='Sp_(2n)(C)', theorem='No nontrivial irreducible preserves torus/global regularity for n>=2; n=1 only defining is positive',
        arbitrary_rank_proof='TYPE_C_REGULARITY.md, Section 1a (central--Weyl) and Sections 2-5 (branching)',
        finite_all_partitions_n_range=[1, 4], finite_all_partitions_maximum_size=6,
        additional_rank_five_symmetric_degrees=list(range(1, 7)),
        additional_rank_five_mixed_partitions=4, character_cases=len(cases),
        Freudenthal_weight_and_alternant_comparisons=len(cases),
        rank_one_string_reconstructions=len(cases), direct_trace_coordinate_checks=coordinate_checks,
        symmetric_monomial_character_comparisons=symmetric_count,
        inherited_fundamental_cases=inherited_count,
        nonfundamental_exact_algebraic_witnesses=len(witnesses),
        total_distinct_weights_in_cases=total_distinct_weights,
        maximum_grid_trials_used=maximum_trials,
        universal_Qi_trace_checks=len(universal_checks),
        universal_witnesses=universal_checks,
        same_central_parity_reducible_checks=len(same_parity_sums),
        same_central_parity_reducible_witnesses=same_parity_sums,
        witnesses=witnesses,
        Sp4_adjoint_example=dict(first_pair_sum='-5/4', second_torus_coordinate='2',
            second_image_coordinate='15/32'),
        runtime_API_changed=False, novelty_claim=False,
        scope='Exact finite independent corroboration; not an extrapolation proof, external review, degree or integral-image classification.',
        source_sha256={p: sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in paths})
    (ROOT/'type-c-regularity-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('witnesses', 'universal_witnesses', 'same_central_parity_reducible_witnesses', 'source_sha256')}, indent=2))


if __name__ == '__main__':
    main()
