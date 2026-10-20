"""Check the finite identities used in the structural F4 duality argument."""
from hashlib import sha256
import json
from pathlib import Path
import sys
from sage.all import PolynomialRing, QQ
from sage_borel import exceptional_borel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.test_exceptional import RootChevalleyOracle


def main():
    oracle = RootChevalleyOracle("F4")
    identities = 0
    reflections = 0
    def zero(f):
        nonlocal identities
        assert not oracle.expand(f)
        identities += 1
    maps = {i: exceptional_borel("F4", i) for i in range(1, 5)}
    reduced = {}
    for i, b in maps.items():
        polynomials = {g.name: b.generator_polynomial(g.name) for g in b.generators()}
        if i in (1, 4):
            h, u = polynomials['y4'], polynomials['e2']
            zero(polynomials['e1'] + h*h)
            zero(polynomials['e3'] + h*h*u)
            reduced[i] = (h, u)
        else:
            h = next(polynomials[g.name] for g in b.generators() if g.polynomial_degree == 1)
            second_h = [polynomials[g.name] for g in b.generators() if g.polynomial_degree == 1][1]
            A, B = [polynomials[g.name] for g in b.generators() if g.polynomial_degree == 2]
            Q = A+B
            zero(second_h-h)
            zero(A - (-11*h*h+9*Q if i == 2 else 11*h*h-8*Q))
            zero(B - (11*h*h-8*Q if i == 2 else -11*h*h+9*Q))
            reduced[i] = (h, Q, polynomials['e3'])
    for i, b in maps.items():
        target = maps[5-i]
        substitution = (2*target.weights[3], 2*target.weights[2],
                        target.weights[1], target.weights[0])
        for node in b.nodes:
            for weight in b.weights:
                assert b.reflect(weight, node)(*substitution) == target.reflect(weight(*substitution), 5-node)
                reflections += 1
        for f, g in zip(reduced[i], reduced[5-i]):
            scale = 2**f.total_degree() if i in (1, 2) else 1
            zero(target.ring(f(*substitution))-scale*g)
    path = ROOT / 'experiments/exceptional-integral-results.json'
    data = json.loads(path.read_text())
    assert data['result'] == 'passed'
    cases = {c['parabolic']: c for c in data['cases']
             if c['kind'] == 'F4' and c['form'] == 'image_denominators_cleared'}
    entries = 0
    degree_matrices = 0
    for source, target in ((4, 1), (3, 2)):
        a, b = cases[source], cases[target]
        ra = PolynomialRing(QQ, len(a['retained_generators']), names=a['retained_generators'])
        rb = PolynomialRing(QQ, len(b['retained_generators']), names=b['retained_generators'])
        for x, y in zip(a['degrees'], b['degrees']):
            target_rows = {oracle.W.from_reduced_word(w): j for j, w in enumerate(y['schubert_words'])}
            target_columns = {tuple(rb(m).exponents()[0]): j for j, m in enumerate(y['monomials'])}
            assert x['rank'] == y['rank']
            for row, word in enumerate(x['schubert_words']):
                target_row = target_rows[oracle.W.from_reduced_word([5-node for node in word])]
                scalar = 2**sum(node in (1,2) for node in word)
                for col, monomial in enumerate(x['monomials']):
                    target_col = target_columns[tuple(ra(monomial).exponents()[0])]
                    assert y['matrix'][target_row][target_col] == scalar*x['matrix'][row][col]
                    entries += 1
            def odd_part(n):
                while n % 2 == 0:
                    n //= 2
                return n
            assert [odd_part(n) for n in x['smith_factors']] == [odd_part(n) for n in y['smith_factors']]
            degree_matrices += 1
    result = {'result': 'passed', 'root_intertwining_identities': reflections,
              'independent_chevalley_image_identities': identities,
              'paired_degree_matrices': degree_matrices, 'row_scaling_entries': entries,
              'scope': 'Finite identities supporting the written root-duality proof, and complete original-form paired matrix comparison.',
              'source_hash_semantics': 'SHA256 of UTF-8 file bytes with CRLF normalized to LF, matching Git text blobs.',
              'source_sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                                for p in (path, Path(__file__).resolve(), ROOT/'tests/test_exceptional.py',
                                          ROOT/'src/sage_borel/exceptional.py',
                                          ROOT/'src/sage_borel/data/exceptional_generators.json')}}
    (ROOT/'experiments/exceptional-f4-duality.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
