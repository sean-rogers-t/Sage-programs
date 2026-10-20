"""Independent F4 certificate audit, Python standard library only.

Reconstructs the torus, all polynomial source lattices, direct low-degree
Laurent checks, and image-lattice arithmetic. The separate Sage experiment
compares every Borel column by localization and root/coroot Chevalley products.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
from time import perf_counter
import json

from audit_e6_character_form import inverse, matmul, combine, multiply

ROOT = Path(__file__).resolve().parents[1]
ZERO = (0, 0, 0, 0)
ONE = {ZERO: 1}
C = ((2, -1, 0, 0), (-1, 2, -1, 0), (0, -2, 2, -1), (0, 0, -1, 2))
ROOTS = tuple(tuple(C[j][i] for j in range(4)) for i in range(4))


def determinant(A):
    """Fraction-free Bareiss determinant, independent of Sage's algorithms."""
    A = [list(map(int, r)) for r in A]
    if not A:
        return 1
    sign, previous = 1, 1
    for k in range(len(A)-1):
        i = next((i for i in range(k, len(A)) if A[i][k]), None)
        if i is None:
            return 0
        if i != k:
            A[k], A[i] = A[i], A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1, len(A)):
            lower = A[i][k]
            for j in range(k+1, len(A)):
                numerator = A[i][j]*pivot-lower*A[k][j]
                assert numerator % previous == 0
                A[i][j] = numerator//previous
            A[i][k] = 0
        previous = pivot
    return sign*A[-1][-1]


def power(f, n):
    out = ONE
    for _ in range(n):
        out = multiply(out, f)
    return out


def orbit(w):
    found, pending = {w}, [w]
    while pending:
        w = pending.pop()
        for i, a in enumerate(ROOTS):
            v = tuple(w[j]-w[i]*a[j] for j in range(4))
            if v not in found:
                found.add(v)
                pending.append(v)
    return found


def substitute(f, values):
    out = {}
    for e, c in f.items():
        term = ONE
        for v, n in zip(values, e):
            term = multiply(term, power(v, n))
        out = combine([out, term], [1, c])
    return out


def exact_integer_columns(polys, B):
    """Combine rational polynomial columns using integer arithmetic."""
    denominators = [lcm(*(F(c).denominator for c in f.values())) for f in polys]
    integers = [{e: int(c*d) for e, c in f.items()} for f, d in zip(polys, denominators)]
    output = []
    for col in zip(*B):
        coefficients = [F(c)/d for c, d in zip(col, denominators)]
        den = lcm(*(c.denominator for c in coefficients))
        p = combine(integers, [int(c*den) for c in coefficients])
        assert all(c % den == 0 for c in p.values())
        output.append({e: c//den for e, c in p.items() if c})
    return output


def minor_certificate(polys, certificate, row_exponents):
    divisor = 0
    for minor in certificate['maximal_minors']:
        selected = [row_exponents[i] for i in minor['rows']]
        det = abs(determinant([[f.get(e, 0) for f in polys] for e in selected]))
        assert det == int(minor['determinant'])
        divisor = gcd(divisor, det)
    assert divisor == certificate['minor_gcd'] == 1


def smith_factors(H):
    previous, previous_factor = 1, 1
    result = []
    for k in range(1, len(H)+1):
        divisor = 0
        # Every k-minor is divisible by this bound once earlier factors
        # have been certified. Equality permits an exact early exit.
        lower_bound = previous*previous_factor
        stop = False
        for rs in combinations(range(len(H)), k):
            for cs in combinations(range(len(H)), k):
                divisor = gcd(divisor, abs(determinant([[H[i][j] for j in cs] for i in rs])))
                if divisor == lower_bound:
                    stop = True
                    break
            if stop:
                break
        assert divisor
        previous_factor = divisor//previous
        result.append(previous_factor)
        previous = divisor
    return result


def main():
    started = perf_counter()
    path = ROOT/'experiments/f4-character-form-results.json'
    data = json.loads(path.read_text())
    assert data['result'] == 'passed' and data['cartan'] == list(map(list, C))
    weights = orbit((0, 0, 0, 1))
    assert len(weights) == 24
    assert {tuple(w): m for w, m in data['weights']} == dict(
        [(w, 1) for w in weights]+[(ZERO, 2)])
    gram = [[sum(w[i]*w[j] for w in weights) for j in range(4)] for i in range(4)]
    assert gram == [[int(c) for c in r] for r in data['trace_gram']]
    theta = [combine([{w: w[i] for w in weights if w[i]} for i in range(4)], r)
             for r in inverse(gram)]
    assert all(sum(f.values()) == 0 for f in theta)
    assert [[sum(c*e[j] for e, c in f.items()) for j in range(4)] for f in theta] == [
        [int(i == j) for j in range(4)] for i in range(4)]
    bs = [combine([f], [c]) for f, c in zip(theta, (6, 6, 12, 12))]
    assert all(F(c).denominator == 1 for f in bs for c in f.values())
    unit_exponents = list(map(tuple, data['linear_unit_minor_exponents']))
    assert abs(determinant(unit_exponents)) == 1  # weighted weights span X(T)
    assert abs(determinant([[f.get(e, 0) for f in bs] for e in unit_exponents])) == 1
    assert data['linear_basis'] == [[str(c) for c in r] for r in
        ((6, 0, 0, 0), (0, 6, 0, 0), (0, 0, 12, 0), (0, 0, 0, 12))]
    cert = data['torus_certificate']
    U = cert['basis_change']
    assert abs(determinant(U)) == 1
    transformed = [combine(bs, r) for r in U]
    exponents = [max(f) for f in transformed]
    assert [f[e] for f, e in zip(transformed, exponents)] == [1]*4
    assert determinant(exponents)
    assert cert['certified'] and cert['exceptional_primes'] == []
    assert cert['cases'] == [[0, list(map(list, exponents)), [1]*4, 4]]
    # Weyl-stability plus a basis of the indicated lattice identifies Q_long.
    long = orbit(ROOTS[0])
    assert len(long) == 24
    K = tuple(tuple(C[i][j]*(1 if j < 2 else 2) for j in range(4)) for i in range(4))
    inverseK = inverse(K)
    assert all(F(c).denominator == 1 for w in long
               for row in matmul(inverseK, [[a] for a in w]) for c in row)
    assert all(tuple(K[i][j] for i in range(4)) in long or
               tuple(K[i][j]//2 for i in range(4)) in weights for j in range(4))
    # 2*short roots belong to the long lattice: 2 alpha3/4 are sums of
    # two long roots, checked exhaustively without a root-length convention.
    sums = {tuple(a+b for a, b in zip(u, v)) for u in long for v in long}
    assert all(tuple(K[i][j] for i in range(4)) in long | sums for j in range(4))
    transition = [[F(c) for c in r] for r in data['long_root_transition']]
    assert transition == [[F(K[i][j], 1 if i < 2 else 2) for j in range(4)] for i in range(4)]
    assert abs(determinant(transition)) == 1
    W = [{tuple(int(j == i) for j in range(4)): 1} for i in range(4)]
    source_rows, source_columns, image_rows, direct_rows = 0, 0, 0, 0
    for case in data['cases']:
        node = case['parabolic']
        generators = [{tuple(e): F(c) for e, c in terms} for terms in case['generator_terms']]
        for f in generators:
            for i in range(4):
                if i+1 == node:
                    continue
                values = list(W)
                values[i] = combine([W[i], combine(W, ROOTS[i])], [1, -1])
                assert substitute(f, values) == f
        normalized = [{e: c/F(2)**(e[2]+e[3]) for e, c in f.items()} for f in generators]
        @lru_cache(None)
        def monomial(p):
            if not any(p):
                return ONE
            j = next(j for j, e in enumerate(p) if e)
            q = list(p)
            q[j] -= 1
            return multiply(normalized[j], monomial(tuple(q)))
        for row in case['degrees']:
            d = row['degree']
            powers = list(map(tuple, row['monomial_exponents']))
            # Enumerate the complete rational invariant monomial basis.
            def enumerate_powers(weights, remaining, prefix=()):
                if not weights:
                    return [prefix] if remaining == 0 else []
                return [p for e in range(remaining//weights[0]+1)
                        for p in enumerate_powers(weights[1:], remaining-e*weights[0], prefix+(e,))]
            assert powers == enumerate_powers(case['generator_degrees'], d)
            polys = [monomial(p) for p in powers]
            support = sorted(set().union(*(set(f) for f in polys)))
            # Runtime merges duplicate coefficient rows then orders them
            # lexicographically. Reconstruct that ordering independently.
            row_map = {tuple(f.get(e, 0) for f in polys): e for e in support}
            row_exponents = [row_map[r] for r in sorted(row_map)]
            assert len(row_map) == row['coefficient_rows']
            B0 = [list(map(F, r)) for r in row['normalized_source_basis']]
            integral_polys = exact_integer_columns(polys, B0)
            minor_certificate(integral_polys, row['completeness_certificate'], row_exponents)
            source_rows += 1
            source_columns += len(powers)
            B = [[6**d*c for c in r] for r in B0]
            image = [[int(c) for c in r] for r in row['image_matrix']]
            assert matmul([list(map(F, r)) for r in row['rational_image_matrix']], B) == image
            H = [[int(c) for c in r] for r in row['image_basis']]
            coordinates = matmul(inverse(H), image)
            assert all(c.denominator == 1 for r in coordinates for c in r)
            minor_polys = [{(j,): int(c) for j, c in enumerate(r) if c} for r in coordinates]
            minor_certificate(minor_polys, row['image_completeness_certificate'],
                              [(j,) for j in range(len(powers))])
            assert smith_factors(H) == row['smith_factors']
            for comparison in row['comparisons']:
                oldH = comparison['old_basis']
                old_inverse = inverse(oldH)
                new_inverse = inverse(H)
                assert all(c.denominator == 1 for r in matmul(old_inverse, H) for c in r) == comparison['old_contains_character']
                assert all(c.denominator == 1 for r in matmul(new_inverse, oldH) for c in r) == comparison['character_contains_old']
                S, I = comparison['sum_basis'], comparison['intersection_basis']
                assert abs(determinant(I)*determinant(S)) == abs(determinant(H)*determinant(oldH))
                # Containment holds in every saved F4 case, so the actual
                # sum must equal the old lattice and the intersection the
                # character lattice. Check equal indices, not merely the
                # product identity (which alone would be insufficient).
                assert comparison['old_contains_character']
                assert abs(determinant(S)) == abs(determinant(oldH))
                assert abs(determinant(I)) == abs(determinant(H))
                for A, expected in (
                    (matmul(inverse(S), H), comparison['sum_mod_character']),
                    (matmul(inverse(S), oldH), comparison['sum_mod_old']),
                    (matmul(old_inverse, I), comparison['old_mod_intersection']),
                    (matmul(new_inverse, I), comparison['character_mod_intersection'])):
                    assert all(c.denominator == 1 for r in A for c in r)
                    assert smith_factors(A) == expected
            image_rows += 1
        pulled = [substitute(f, theta) for f in generators]
        for direct in data['direct_laurent_rows']:
            if direct['parabolic'] != node:
                continue
            polys = []
            blocks = []
            degrees = range(direct['degree']+1) if direct['mode'] == 'filtered' else (direct['degree'],)
            labels = []
            for d in degrees:
                r = case['degrees'][d]
                blocks.append([[6**d*F(c) for c in row] for row in r['normalized_source_basis']])
                for p in r['monomial_exponents']:
                    f = ONE
                    for g, e in zip(pulled, p):
                        f = multiply(f, power(g, e))
                    polys.append(f)
                    labels.append([d, p])
            assert labels == direct['labels']
            B = [[F(0)]*len(polys) for _ in polys]
            offset = 0
            for block in blocks:
                for i, r in enumerate(block):
                    for j, c in enumerate(r):
                        B[offset+i][offset+j] = c
                offset += len(block)
            support = sorted(set().union(*(set(f) for f in polys)))
            assert len(support) == direct['support_rows']
            integral = exact_integer_columns(polys, B)
            minor_certificate(integral, direct['certificate'], support)
            direct_rows += 1
        print("Audited P", node, "elapsed", round(perf_counter()-started, 2), flush=True)
        monomial.cache_clear()
    assert source_rows == image_rows == data['degree_rows'] == 74
    assert source_columns == data['independent_columns']
    assert direct_rows == 32
    files = [path, Path(__file__).resolve(), ROOT/'experiments/audit_e6_character_form.py']
    result = dict(result='passed', polynomial_source_lattices=source_rows,
                  source_basis_columns=source_columns, image_lattices=image_rows,
                  direct_laurent_lattices=direct_rows,
                  seconds=round(perf_counter()-started, 3),
                  scope='Python standard library: independently reconstructed torus and all source coefficient lattices; exact integral combinations and gcd-one minor certificates; all image-lattice Smith factors and inclusion/quotient arithmetic. Borel columns independently checked by two algorithms in the Sage experiment.',
                  source_hash_semantics='SHA256 with CRLF normalized to LF.',
                  source_sha256={p.relative_to(ROOT).as_posix():
                    sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/f4-character-form-audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
