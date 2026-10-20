"""Full F4 character sources and all four maximal-parabolic integral images."""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json
import sys

from sage.all import QQ, ZZ, GF, CartanMatrix, identity_matrix, matrix, prod, vector
from sage_borel import f4_character_source
from e6_character_form import certify_preimage
from e7_character_form import lattice_compare
from exceptional_localization import FixedPointBorel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from exceptional_chevalley_audit import ChevalleyAction


def rows(M):
    return [[str(c) for c in row] for row in M]


def terms(f):
    return [[list(map(int, e)), str(c)] for e, c in sorted(f.dict().items())]


def main():
    start = perf_counter()
    first = f4_character_source(1)
    s = first.springer_map
    C = CartanMatrix("F4")
    weights, pending = {(0, 0, 0, 1)}, [(0, 0, 0, 1)]
    while pending:
        w = pending.pop()
        for i in range(4):
            v = tuple(w[j]-w[i]*C[j, i] for j in range(4))
            if v not in weights:
                weights.add(v)
                pending.append(v)
    assert len(weights) == 24
    expected = {tuple(map(int, w)): 1 for w in weights}
    expected[(0, 0, 0, 0)] = 2
    assert dict(s.weight_multiplicities) == expected
    assert s.trace_gram == matrix(QQ, [[12, -6, 0, 0], [-6, 12, -12, 0],
                                     [0, -12, 24, -12], [0, 0, -12, 24]])
    assert s.linear_basis_matrix == matrix.diagonal(QQ, [6, 6, 12, 12])
    assert matrix(QQ, 4, 4, lambda i, j:
                  s.theta[i].derivative(s.torus_ring.gen(j))(1, 1, 1, 1)) == identity_matrix(4)
    # The first two simple roots are long. alpha1,alpha2,2alpha3,2alpha4
    # span the long-root lattice; its Weyl orbit supplies every long root.
    K = C*matrix.diagonal(QQ, [1, 1, 2, 2])
    T = (s.linear_basis_matrix/6).inverse()*K
    assert T.denominator() == 1 and abs(T.det()) == 1
    support = sorted(set().union(*(set(f.dict()) for f in s.linear_pullbacks)))
    linear = matrix(ZZ, [[f.dict().get(e, 0) for f in s.linear_pullbacks] for e in support])
    pivots = linear.change_ring(GF(2)).transpose().pivots()
    assert abs(linear.matrix_from_rows(pivots).det()) == 1
    report = dict(scope="Full torus source and all four maximal-Levi sources; complete images in every nonzero cohomological degree.",
                  cartan=[list(map(int, r)) for r in C],
                  weights=[[list(w), m] for w, m in s.weight_multiplicities],
                  trace_gram=rows(s.trace_gram), linear_basis=rows(s.linear_basis_matrix),
                  linear_unit_minor_exponents=[list(map(int, support[i])) for i in pivots],
                  torus_certificate=first.torus_certificate.__dict__,
                  long_root_transition=rows(T), cases=[], direct_laurent_rows=[])
    old = json.loads((ROOT/'experiments/exceptional-integral-results.json').read_text())
    oracle = ChevalleyAction("F4")
    column_count = 0
    for node in range(1, 5):
        source = first if node == 1 else f4_character_source(node)
        a = source.rational_algebra
        b = a.borel_map
        loc = FixedPointBorel("F4", node)
        fs = a.generator_polynomials
        # Independent generator degrees and Jacobian certify rational fullness
        # for reflection invariants: product(degrees) = order(W_L).
        assert prod(a.generator_degrees) == (48 if node in (1, 4) else 12)
        assert matrix([[f.derivative(x) for x in b.weights] for f in fs]).det()
        assert all(b.reflect(f, i) == f for f in fs for i in range(1, 5) if i != node)
        vals = [loc.polynomial_values(f) for f in fs]

        @lru_cache(None)
        def chevalley(powers):
            if not any(powers):
                return {oracle.W.one(): QQ.one()}
            i = next(i for i, e in enumerate(powers) if e)
            lower = list(powers)
            lower[i] -= 1
            return oracle.multiply_polynomial(fs[i], chevalley(tuple(lower)))

        case = dict(parabolic=node, complex_dimension=loc.top,
                    names=list(a.generator_names), generator_degrees=list(a.generator_degrees),
                    generator_terms=[terms(f) for f in fs], degrees=[])
        print("P", node, "dimension", loc.top, flush=True)
        for d in range(loc.top+1):
            r = source.degree(d)
            powers = [tuple(map(int, m.exponents()[0])) for m in r.monomials]
            polynomials = [prod(g**e for g, e in zip(source._inputs, p)) for p in powers]
            dictionaries = [f.dict() for f in polynomials]
            exponents = sorted(set().union(*(set(D) for D in dictionaries)))
            coefficients = sorted(set(tuple(D.get(e, 0) for D in dictionaries) for e in exponents))
            M = matrix(QQ, coefficients)
            B0 = r.basis_matrix/6**d
            certificate = certify_preimage(M, B0)
            target = [w for w in loc.basis if w.length() == d]
            columns = []
            for p in powers:
                values = [prod(vals[j][i]**e for j, e in enumerate(p)) for i in range(len(loc.basis))]
                actual = loc.coefficients(values, d)
                expected = {loc.W.from_reduced_word(w.reduced_word()): c
                            for w, c in chevalley(p).items()}
                assert actual == expected, (node, d, p)
                columns.append([actual.get(w, 0) for w in target])
                column_count += 1
            rational = matrix(QQ, columns).transpose()
            image = matrix(ZZ, rational*r.basis_matrix)
            H = image.transpose().hermite_form(include_zero_rows=False).transpose()
            assert H.nrows() == H.ncols() == len(target)
            image_certificate = certify_preimage(
                (H.change_ring(QQ).inverse()*image).transpose(), identity_matrix(len(target)))
            factors = list(map(int, H.elementary_divisors()))
            comparisons = []
            for previous in old['cases']:
                if previous['kind'] != "F4" or previous['parabolic'] != node:
                    continue
                prior = previous['degrees'][d]
                oldM = matrix(ZZ, [prior['matrix'][prior['schubert_words'].index(list(w.reduced_word()))]
                                  for w in target])
                oldH = oldM.transpose().hermite_form(include_zero_rows=False).transpose()
                comparisons.append(lattice_compare(H, oldH, previous['form'], prior['smith_factors']))
            case['degrees'].append(dict(
                degree=d, monomial_exponents=[list(p) for p in powers],
                normalized_source_basis=rows(B0), coefficient_rows=len(coefficients),
                completeness_certificate=certificate,
                schubert_words=[list(w.reduced_word()) for w in target],
                rational_image_matrix=rows(rational),
                image_matrix=rows(image), image_basis=rows(H),
                image_completeness_certificate=image_certificate,
                smith_factors=factors, comparisons=comparisons))
            # Direct Laurent checks include all source generators occurring
            # through degree 3, and compare full filtered intersections.
            if d <= 3:
                for mode in ("homogeneous", "filtered"):
                    ds = range(d+1) if mode == "filtered" else (d,)
                    labels = [(e, m) for e in ds for m in a.monomials(e)]
                    expansions = [s.pullback(a.substitute(m)) for _, m in labels]
                    DD = [f.dict() for f in expansions]
                    supp = sorted(set().union(*(set(D) for D in DD)))
                    A = matrix(QQ, [[D.get(e, 0) for D in DD] for e in supp])
                    from sage.all import block_diagonal_matrix
                    B = block_diagonal_matrix([source.degree(e).basis_matrix for e in ds],
                                              subdivide=False)
                    cert = certify_preimage(A, B)
                    report['direct_laurent_rows'].append(dict(
                        parabolic=node, degree=d, mode=mode,
                        labels=[[e, list(map(int, m.exponents()[0]))] for e, m in labels],
                        support_rows=len(supp), certificate=cert))
            print("P", node, "degree", d, "rank", len(powers), "Smith", factors,
                  "seconds", round(perf_counter()-start, 2), flush=True)
        case['above_top_image_zero'] = not loc.coefficients(
            [v**(loc.top+1) for v in vals[0]], a.generator_degrees[0]*(loc.top+1))
        report['cases'].append(case)
        chevalley.cache_clear()
        (ROOT/'experiments/f4-character-form-results.json').write_text(json.dumps(report, indent=2)+'\n')
    files = [Path(__file__).resolve(), ROOT/'experiments/e6_character_form.py',
             ROOT/'experiments/e7_character_form.py', ROOT/'experiments/exceptional_localization.py',
             ROOT/'experiments/exceptional_chevalley_audit.py',
             ROOT/'experiments/exceptional-integral-results.json', ROOT/'tests/test_exceptional.py',
             ROOT/'src/sage_borel/data/exceptional_generators.json']
    files += sorted((ROOT/'src/sage_borel').glob('*.py'))
    report.update(result='passed', independent_columns=column_count,
                  degree_rows=sum(len(c['degrees']) for c in report['cases']),
                  seconds=round(perf_counter()-start, 3),
                  source_hash_semantics='SHA256 with CRLF normalized to LF.',
                  source_sha256={p.relative_to(ROOT).as_posix():
                    sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/f4-character-form-results.json').write_text(json.dumps(report, indent=2)+'\n')
    print("PASSED", report['degree_rows'], "rows", column_count, "independent columns", flush=True)


if __name__ == "__main__":
    main()
