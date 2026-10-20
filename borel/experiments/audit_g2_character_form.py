"""Independent exact G2 audit: Python standard library, no Sage/library import."""
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
from time import perf_counter
import json
from audit_e6_character_form import determinant, inverse, multiply, combine

ROOT = Path(__file__).resolve().parents[1]
ONE = {(0, 0): F(1)}
W1, W2 = {(1, 0): F(1)}, {(0, 1): F(1)}
ROOTS = ((2, -1), (-3, 2))
SIMPLE = (((-1, 0), (1, 1)), ((1, 3), (0, -1)))


def power(f, n):
    out = ONE
    for _ in range(n):
        out = multiply(out, f)
    return out


def matmul(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def word_matrix(word):
    A = ((1, 0), (0, 1))
    for i in word:
        A = matmul(A, SIMPLE[i-1])
    return A


def substitute(f, values):
    terms = [multiply(power(values[0], a), power(values[1], b)) for a, b in f]
    return combine(terms, list(f.values())) if terms else {}


def divided_difference(f, node):
    i = node-1
    alpha = combine([W1, W2], ROOTS[i])
    variables = [W1, W2]
    variables[i] = combine([variables[i], alpha], [1, -1])
    numerator = combine([f, substitute(f, variables)], [1, -1])
    quotient = {}
    lead = max(alpha)
    while numerator:
        e = max(numerator)
        q = tuple(e[j]-lead[j] for j in range(2))
        assert min(q) >= 0
        c = numerator[e]/alpha[lead]
        quotient[q] = quotient.get(q, 0)+c
        numerator = combine([numerator, multiply({q: c}, alpha)], [1, -1])
    return {e: c for e, c in quotient.items() if c}


def coefficient(f, word):
    for node in reversed(word):
        f = divided_difference(f, node)
    assert all(e == (0, 0) for e in f)
    return f.get((0, 0), F(0))


def orbit(start):
    found, pending = {start}, [start]
    while pending:
        w = pending.pop()
        for i, alpha in enumerate(ROOTS):
            v = tuple(w[j]-w[i]*alpha[j] for j in range(2))
            if v not in found:
                found.add(v)
                pending.append(v)
    return found


def certify_leading(fs, certificate):
    transformed = [combine(fs, row) for row in certificate['basis_change']]
    assert abs(determinant(certificate['basis_change'])) == 1
    exceptional = set()
    for f in transformed:
        n = abs(int(f[max(f)]))
        p = 2
        while p*p <= n:
            if n % p == 0:
                exceptional.add(p)
                while n % p == 0:
                    n //= p
            p += 1
        if n > 1:
            exceptional.add(n)
    assert sorted(exceptional) == certificate['exceptional_primes']
    assert [c[0] for c in certificate['cases']] == [0]+sorted(exceptional)
    for p, exponents, coefficients, rank in certificate['cases']:
        leading = [max(e for e, c in f.items() if not p or c % p) for f in transformed]
        assert [list(e) for e in leading] == exponents
        assert [f[e] for f, e in zip(transformed, leading)] == coefficients
        assert rank == 2 and determinant(leading)
    assert certificate['certified']


def main():
    start = perf_counter()
    path = ROOT/'experiments/g2-character-form-results.json'
    data = json.loads(path.read_text())
    assert data['result'] == 'passed'
    weights = orbit((1, 0)) | {(0, 0)}
    assert len(weights) == 7 and sorted(map(list, weights)) == data['torus_weights']
    gram = [[sum(w[i]*w[j] for w in weights) for j in range(2)] for i in range(2)]
    assert gram == data['trace_gram'] == [[12, -6], [-6, 4]]
    moments = [{w: F(w[i]) for w in weights if w[i]} for i in range(2)]
    theta = [combine(moments, row) for row in inverse(gram)]
    b1, b2 = [combine([f], [c]) for f, c in zip(theta, (6, 2))]
    assert all(c.denominator == 1 for f in (b1, b2) for c in f.values())
    # These two coefficient rows give a unit minor: exact degree-one lattice.
    assert determinant([[f.get(e, 0) for f in (b1, b2)]
                        for e in ((2, -1), (1, 0))]) == 1
    assert data['linear_basis'] == [[6, 0], [0, 2]]
    certify_leading((b1, b2), data['torus_certificate'])
    adj_weights = orbit(ROOTS[0]) | orbit(ROOTS[1])
    assert len(adj_weights) == 12
    adj_gram = [[sum(w[i]*w[j] for w in adj_weights) for j in range(2)] for i in range(2)]
    assert adj_gram == [[4*x for x in row] for row in gram]
    adj_moments = [{w: F(w[i]) for w in adj_weights if w[i]} for i in range(2)]
    adj_theta = [combine(adj_moments, row) for row in inverse(adj_gram)]
    adj_pullbacks = [combine([f], [c]) for f, c in zip(adj_theta, (24, 8))]
    assert all(c.denominator == 1 for f in adj_pullbacks for c in f.values())
    # Coefficients at (2,-1) and (1,0) also give a unit minor here.
    assert determinant([[f.get(e, 0) for f in adj_pullbacks]
                        for e in ((2, -1), (1, 0))]) == 1
    certify_leading(adj_pullbacks, data['adjoint_certificate'])
    pulled = {1: (b1, multiply(b2, combine([b1, b2], [1, -1]))),
              2: (b2, multiply(b1, combine([b2, b1], [3, -1])))}
    positions = 0
    assert len(data['direct_laurent_rows']) == 24
    for row in data['direct_laurent_rows']:
        d = row['degree']
        expected = {(a, b) for a in range(d+1) for b in range(d//2+1)
                    if a+2*b == d or (row['mode'] == 'filtered' and a+2*b < d)}
        assert set(map(tuple, row['exponents'])) == expected
        fs = pulled[row['parabolic']]
        expansions = [multiply(power(fs[0], a), power(fs[1], b)) for a, b in row['exponents']]
        support = sorted(set().union(*(set(f) for f in expansions)))
        N = [[f.get(e, F(0)) for f in expansions] for e in support]
        assert all(c.denominator == 1 for r in N for c in r)
        divisor = 0
        for minor in row['certificate']['maximal_minors']:
            det = abs(determinant([N[i] for i in minor['rows']]))
            assert det == int(minor['determinant'])
            divisor = gcd(divisor, int(det))
        assert divisor == 1
        positions += len(support)*len(expansions)
    assert positions == data['coefficient_positions']
    # Build the 12-element Weyl group without Sage, with shortest-word lengths.
    identity = ((1, 0), (0, 1))
    group, pending = {identity: ()}, [identity]
    while pending:
        A = pending.pop(0)
        for i, reflection in enumerate(SIMPLE, 1):
            B = matmul(A, reflection)
            if B not in group:
                group[B] = group[A]+(i,)
                pending.append(B)
    assert len(group) == 12
    oldpath = ROOT/'experiments/exceptional-integral-results.json'
    old = json.loads(oldpath.read_text())
    x, y = combine([W1], [6]), combine([W2], [2])
    inputs = {1: (x, multiply(y, combine([x, y], [1, -1]))),
              2: (y, multiply(x, combine([y, x], [3, -1])))}
    column_count, comparisons = 0, 0
    for case in data['cases']:
        node = case['parabolic']
        parabolic = {A: w for A, w in group.items()
                     if len(group[matmul(A, SIMPLE[2-node])]) > len(w)}
        assert sorted(map(len, parabolic.values())) == list(range(6))
        first, second = inputs[node]
        for row in case['degrees']:
            d = row['degree']
            assert row['exponents'] == [[d-2*j, j] for j in range(d//2+1)]
            if d < 6:
                word, = row['schubert_words']
                assert word_matrix(word) in parabolic and len(word) == d
            else:
                assert not row['schubert_words']
            values = []
            for (a, b), column in zip(row['exponents'], row['columns']):
                f = multiply(power(first, a), power(second, b))
                if d < 6:
                    c = coefficient(f, word)
                    assert column == [c] and c.denominator == 1
                    values.append(int(c))
                else:
                    for top in group.values():
                        if len(top) == 6:
                            assert coefficient(f, top) == 0
                    assert not column
                column_count += 1
            factor = abs(gcd(*values)) if values else 0
            assert row['smith_factors'] == ([factor] if d < 6 else [])
            if d < 6:
                assert factor > 0
                assert row['h_power_coefficients'] == [coefficient(power((W1, W2)[node-1], d), word)]
            for comparison in row['comparisons']:
                prior_case = next(c for c in old['cases'] if c['kind'] == 'G2'
                                  and c['parabolic'] == node and c['form'] == comparison['form'])
                prior = prior_case['degrees'][d]
                old_factor = abs(gcd(*prior['matrix'][0]))
                assert comparison['old_smith_factors'] == [old_factor]
                assert comparison['character_contained'] and factor % old_factor == 0
                assert comparison['old_mod_character_order'] == factor//old_factor
                comparisons += 1
        # The stated image-ring relation in degree two.
        relation = combine([power(first, 2), second], [1, -3]) if node == 1 else \
                   combine([second, power(first, 2)], [1, -3])
        for word in group.values():
            if len(word) == 2:
                assert coefficient(relation, word) == 0
    assert column_count == data['independent_chevalley_columns'] == 32
    assert comparisons == 24
    files = [path, oldpath, Path(__file__).resolve(), ROOT/'experiments/audit_e6_character_form.py']
    out = dict(result='passed', engine='Python standard library only',
               torus_certificates=2, linear_lattice_unit_minors=2,
               direct_laurent_lattices=24, coefficient_positions=positions,
               independent_bgg_columns=column_count, embedded_comparisons=comparisons,
               image_ring_relations=True, seconds=round(perf_counter()-start, 3),
               source_hash_semantics='SHA256 with CRLF normalized to LF.',
               source_sha256={p.relative_to(ROOT).as_posix():
                   sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/g2-character-form-audit.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
