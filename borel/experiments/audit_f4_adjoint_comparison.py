"""Independent standard-library audit of the F4 adjoint scaling proof."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json

from audit_f4_character_form import (
    ZERO, ONE, ROOTS, orbit, determinant, inverse, combine, multiply,
    substitute, power, exact_integer_columns, minor_certificate, smith_factors)

ROOT = Path(__file__).resolve().parents[1]


def main():
    start = perf_counter()
    path = ROOT/'experiments/f4-adjoint-comparison-results.json'
    data = json.loads(path.read_text())
    base_path = ROOT/'experiments/f4-character-form-results.json'
    base = json.loads(base_path.read_text())
    audit_path = ROOT/'experiments/f4-character-form-audit.json'
    previous = json.loads(audit_path.read_text())
    assert data['result'] == base['result'] == previous['result'] == 'passed'
    assert sha256(base_path.read_bytes().replace(b'\r\n', b'\n')).hexdigest() == previous['source_sha256'][base_path.relative_to(ROOT).as_posix()]
    short, long = orbit((0, 0, 0, 1)), orbit(ROOTS[0])
    assert len(short) == len(long) == 24 and not short.intersection(long)
    weights = short | long
    assert {tuple(w): m for w, m in data['weights']} == dict(
        [(w, 1) for w in weights]+[(ZERO, 4)])
    assert sum(m for _, m in data['weights']) == data['adjoint_dimension'] == 52
    gram = [[sum(w[i]*w[j] for w in weights) for j in range(4)] for i in range(4)]
    assert gram == [[int(c) for c in r] for r in data['trace_gram']]
    assert gram == [[3*int(c) for c in r] for r in base['trace_gram']]
    theta = [combine([{w: w[i] for w in weights if w[i]} for i in range(4)], r)
             for r in inverse(gram)]
    assert all(sum(f.values()) == 0 for f in theta)
    assert [[sum(c*e[j] for e, c in f.items()) for j in range(4)] for f in theta] == [
        [int(i == j) for j in range(4)] for i in range(4)]
    scalars = (18, 18, 36, 36)
    assert data['linear_basis'] == [[str(scalars[i] if i == j else 0) for j in range(4)]
                                    for i in range(4)]
    bs = [combine([f], [c]) for f, c in zip(theta, scalars)]
    assert all(c.denominator == 1 for f in bs for c in f.values())
    unit_exponents = list(map(tuple, data['linear_unit_minor_exponents']))
    assert abs(determinant([[f.get(e, 0) for f in bs] for e in unit_exponents])) == 1
    cert = data['certificate']
    assert abs(determinant(cert['basis_change'])) == 1
    transformed = [combine(bs, r) for r in cert['basis_change']]
    assert cert['certified'] and cert['exceptional_primes'] == [2]
    assert [row[0] for row in cert['cases']] == [0, 2]
    leading_determinants = {}
    for p, es, cs, rank in cert['cases']:
        actual = [max(e for e, c in f.items() if not p or c % p) for f in transformed]
        assert list(map(list, actual)) == es
        assert [f[e] for f, e in zip(transformed, actual)] == cs
        assert cs == ([1, 1, 2, 1] if not p else [1, 1, 1, 1])
        assert rank == 4 and determinant(actual)
        leading_determinants[str(p)] = determinant(actual)
    witness = data['map_distinction']
    e = tuple(witness['exponent'])
    assert e in long and e not in short
    assert [f.get(e, 0) for f in theta] == list(map(F, witness['adjoint_coefficients']))
    assert any(F(c) for c in witness['adjoint_coefficients'])
    assert witness['minimal_coefficients'] == ['0']*4
    image_rows, direct_rows = 0, 0
    for case, oldcase in zip(data['cases'], base['cases']):
        assert case['parabolic'] == oldcase['parabolic']
        assert len(case['degrees']) == len(oldcase['degrees'])
        for row, old in zip(case['degrees'], oldcase['degrees']):
            d = row['degree']
            assert d == old['degree']
            q = 3**d
            H = [[int(c) for c in r] for r in old['image_basis']]
            newH = [[int(c) for c in r] for r in row['adjoint_image_basis']]
            assert row['rank'] == len(H) and int(row['degree_scale']) == q
            assert newH == [[q*c for c in r] for r in H]
            # Independent determinantal-divisor reconstruction, not a
            # second call to the Sage Smith algorithm.
            oldfactors = smith_factors(H)
            assert oldfactors == row['minimal_smith_factors'] == old['smith_factors']
            assert smith_factors(newH) == row['adjoint_smith_factors'] == [q*f for f in oldfactors]
            assert row['minimal_mod_adjoint_factors'] == ([q]*len(H) if d else [])
            image_rows += 1
        generators = [{tuple(e): F(c) for e, c in terms} for terms in oldcase['generator_terms']]
        pulled = [substitute(f, theta) for f in generators]
        for direct in data['direct_laurent_rows']:
            if direct['parabolic'] != case['parabolic']:
                continue
            ds = range(direct['degree']+1) if direct['mode'] == 'filtered' else (direct['degree'],)
            polys, blocks, labels = [], [], []
            for d in ds:
                r = oldcase['degrees'][d]
                blocks.append([[18**d*F(c) for c in row] for row in r['normalized_source_basis']])
                for p in r['monomial_exponents']:
                    f = ONE
                    for g, n in zip(pulled, p):
                        f = multiply(f, power(g, n))
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
        print('Audited P', case['parabolic'], 'elapsed', round(perf_counter()-start, 2), flush=True)
    assert image_rows == 74 and direct_rows == 32
    assert len(data['fresh_api_rows']) == data['fresh_api_degrees'] == 20
    for r in data['fresh_api_rows']:
        old = base['cases'][r['parabolic']-1]['degrees'][r['degree']]
        assert [list(map(F, row)) for row in r['source_basis']] == [
            [18**r['degree']*F(c) for c in row] for row in old['normalized_source_basis']]
    files = [Path(__file__).resolve(), path, base_path, audit_path,
             ROOT/'experiments/audit_f4_character_form.py',
             ROOT/'experiments/audit_e6_character_form.py']
    result = dict(result='passed', leading_exponent_determinants=leading_determinants,
                  image_scaling_rows=image_rows, direct_laurent_lattices=direct_rows,
                  scope='Independent root orbits, trace projection, linear unit minor, both characteristic cases, all 74 derived embedded image scalings and determinantal-divisor Smith lists, and 32 fresh Laurent completeness certificates. Earlier minimal Borel columns are reused, not claimed as freshly recomputed.',
                  seconds=round(perf_counter()-start, 3),
                  source_hash_semantics='SHA256 with CRLF normalized to LF.',
                  source_sha256={p.relative_to(ROOT).as_posix():
                    sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/f4-adjoint-comparison-audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
