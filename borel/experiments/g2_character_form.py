"""Full G2 character sources/images and exact independent checks.

Run: PYTHONPATH=src sage -python experiments/g2_character_form.py
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
import sys
from sage.all import QQ, ZZ, identity_matrix, matrix, prod, gcd
from sage_borel import springer_character_map, g2_character_algebra
from e6_character_form import integral_preimage, certify_preimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.test_exceptional import RootChevalleyOracle


def powers(degree):
    return [(degree-2*j, j) for j in range(degree//2+1)]


def main():
    start = perf_counter()
    source = springer_character_map("G2", (1, 0))
    # Independent reflection orbit with alpha1=(2,-1), alpha2=(-3,2).
    roots = ((2, -1), (-3, 2))
    weights, pending = {(1, 0)}, [(1, 0)]
    while pending:
        w = pending.pop()
        for i, alpha in enumerate(roots):
            v = tuple(w[j]-w[i]*alpha[j] for j in range(2))
            if v not in weights:
                weights.add(v)
                pending.append(v)
    assert len(weights) == 6
    weights.add((0, 0))
    assert dict(source.weight_multiplicities) == {w: 1 for w in weights}
    assert source.trace_gram == matrix(QQ, [[12, -6], [-6, 4]])
    assert source.linear_basis_matrix == matrix(QQ, [[6, 0], [0, 2]])
    cert = source.certify_polynomial_source()
    assert cert.certified and not cert.exceptional_primes
    t, u = source.torus_ring.gens()
    # Published Part III torus coordinates are (t, u/t).
    a, b = t, u/t
    published = ((2*a-2/a+b-1/b+a/b-b/a)/6, (a-1/a-1/b+b)/2)
    assert source.theta == published
    assert matrix(QQ, 2, 2, lambda i, j:
                  source.theta[i].derivative(source.torus_ring.gen(j))(1, 1)) == identity_matrix(2)
    adjoint = springer_character_map("G2", (0, 1))
    adj_cert = adjoint.certify_polynomial_source([[1, -3], [0, 1]])
    assert adj_cert.certified and adjoint.linear_basis_matrix == 4*source.linear_basis_matrix
    report = dict(scope="Full seven-dimensional G2 character forms and both maximal-parabolic images; adjoint torus source additionally certified.",
                  torus_weights=sorted(map(list, weights)),
                  trace_gram=[list(map(int, row)) for row in source.trace_gram],
                  linear_basis=[list(map(int, row)) for row in source.linear_basis_matrix],
                  torus_certificate=cert.__dict__, adjoint_certificate=adj_cert.__dict__,
                  direct_laurent_rows=[], cases=[])
    old_path = ROOT/'experiments/exceptional-integral-results.json'
    old = json.loads(old_path.read_text())
    oracle = RootChevalleyOracle("G2")
    checked_columns = 0
    for node in (1, 2):
        algebra = g2_character_algebra(node)
        generators = algebra.generator_polynomials
        borel = algebra.borel_map
        assert all(borel.reflect(f, 3-node) == f for f in generators)
        pulled = [source.pullback(f) for f in generators]
        for d in range(6):
            for mode in ("homogeneous", "filtered"):
                exponents = [p for e in range(d+1) if mode == "filtered" or e == d
                             for p in powers(e)]
                expansions = [prod(f**e for f, e in zip(pulled, p)) for p in exponents]
                dictionaries = [f.dict() for f in expansions]
                support = sorted(set().union(*(set(D) for D in dictionaries)))
                A = matrix(QQ, [[D.get(e, 0) for D in dictionaries] for e in support])
                B = integral_preimage(A)
                assert B.denominator() == 1 and abs(B.det()) == 1
                certificate = certify_preimage(A, identity_matrix(len(exponents)))
                report['direct_laurent_rows'].append(dict(
                    parabolic=node, degree=d, mode=mode, exponents=exponents,
                    support_rows=len(support), coefficient_positions=A.nrows()*A.ncols(),
                    certificate=certificate))
        case = dict(parabolic=node, generators=list(map(str, generators)), degrees=[])
        for d in range(7):
            row = algebra.integral_degree(d)
            exponents = powers(d)
            columns = []
            for p in exponents:
                f = prod(g**e for g, e in zip(generators, p))
                actual = borel.expand(f)
                independent = oracle.expand(f)
                independent = {borel.W.from_reduced_word(w.reduced_word()): c
                               for w, c in independent.items()}
                assert dict(actual) == independent
                columns.append([int(actual.get(w, 0)) for w in row.schubert_indices])
                checked_columns += 1
            # Match columns by exponents, independently of monomial order.
            for j, monomial in enumerate(row.monomials):
                p = tuple(map(int, monomial.exponents()[0]))
                assert list(row.matrix.column(j)) == columns[exponents.index(p)]
            factors = [] if d == 6 else [int(gcd([c[0] for c in columns]))]
            assert factors == list(row.smith_factors)
            comparisons = []
            for previous in old['cases']:
                if previous['kind'] != "G2" or previous['parabolic'] != node:
                    continue
                if d >= len(previous['degrees']):
                    assert d == 6 and not factors
                    continue
                prior = previous['degrees'][d]
                if d < 6:
                    old_factor = int(gcd(prior['matrix'][0]))
                    assert [old_factor] == prior['smith_factors']
                    assert factors[0] % old_factor == 0
                    quotient = factors[0]//old_factor
                else:
                    assert not factors and not prior['smith_factors']
                    quotient = 1
                comparisons.append(dict(form=previous['form'],
                    old_smith_factors=prior['smith_factors'], character_contained=True,
                    old_mod_character_order=quotient))
            h = borel.weights[node-1]
            hp = borel.expand(h**d)
            case['degrees'].append(dict(degree=d, exponents=exponents,
                schubert_words=[list(w.reduced_word()) for w in row.schubert_indices],
                columns=columns, smith_factors=factors,
                h_power_coefficients=[int(hp.get(w, 0)) for w in row.schubert_indices],
                comparisons=comparisons))
        first, second = generators
        if node == 1:
            assert not borel.expand(first**2-3*second)
            assert not borel.expand(second**3)
            assert borel.expand(second) == borel.expand(12*borel.weights[0]**2)
        else:
            assert not borel.expand(second-3*first**2)
            assert not borel.expand(first**6)
        report['cases'].append(case)
    files = [Path(__file__).resolve(), ROOT/'experiments/e6_character_form.py',
             old_path, ROOT/'tests/test_exceptional.py', ROOT/'tests/test_characters.py']
    files += sorted((ROOT/'src/sage_borel').glob('*.py'))
    report.update(result='passed', independent_chevalley_columns=checked_columns,
                  direct_laurent_lattices=len(report['direct_laurent_rows']),
                  coefficient_positions=sum(r['coefficient_positions'] for r in report['direct_laurent_rows']),
                  seconds=round(perf_counter()-start, 3),
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():
                     sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/g2-character-form-results.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('result', 'independent_chevalley_columns',
          'direct_laurent_lattices', 'coefficient_positions', 'seconds')}, indent=2))


if __name__ == "__main__":
    main()
