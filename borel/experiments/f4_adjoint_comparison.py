"""F4 adjoint all-primes certificate and degreewise comparison with V(omega4).

All 74 image comparisons are deduced from the proved scaling and the audited
minimal-representation tables. Fresh API checks cover degrees 0..4 and direct
Laurent checks cover homogeneous/filtered degrees 0..3, for all four Levis.
"""
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json

from sage.all import (QQ, ZZ, CartanMatrix, block_diagonal_matrix,
                      identity_matrix, matrix)
from sage_borel import f4_character_source
from e6_character_form import certify_preimage

ROOT = Path(__file__).resolve().parents[1]


def rows(M):
    return [[str(c) for c in row] for row in M]


def main():
    start = perf_counter()
    base_path = ROOT/'experiments/f4-character-form-results.json'
    baseline = json.loads(base_path.read_text())
    baseline_audit_path = ROOT/'experiments/f4-character-form-audit.json'
    baseline_audit = json.loads(baseline_audit_path.read_text())
    assert baseline['result'] == baseline_audit['result'] == 'passed'
    assert sha256(base_path.read_bytes().replace(b'\r\n', b'\n')).hexdigest() == baseline_audit['source_sha256'][base_path.relative_to(ROOT).as_posix()]
    minimal = f4_character_source(1)
    adjoint = f4_character_source(1, representation='adjoint')
    s, t = minimal.springer_map, adjoint.springer_map
    C = CartanMatrix('F4')
    found, pending = {(1, 0, 0, 0)}, [(1, 0, 0, 0)]
    while pending:
        w = pending.pop()
        for i in range(4):
            v = tuple(w[j]-w[i]*C[j, i] for j in range(4))
            if v not in found:
                found.add(v)
                pending.append(v)
    short = {w for w, _ in s.weight_multiplicities if any(w)}
    assert len(found) == len(short) == 24 and not found.intersection(short)
    weights = {tuple(map(int, w)): 1 for w in found | short}
    weights[(0, 0, 0, 0)] = 4
    assert dict(t.weight_multiplicities) == weights and t.dimension == 52
    assert t.trace_gram == 3*s.trace_gram
    assert t.linear_basis_matrix == 3*s.linear_basis_matrix
    assert matrix(QQ, 4, 4, lambda i, j:
                  t.theta[i].derivative(t.torus_ring.gen(j))(1, 1, 1, 1)) == identity_matrix(4)
    cert = adjoint.torus_certificate
    assert cert.certified and cert.exceptional_primes == (2,)
    exponents = baseline['linear_unit_minor_exponents']
    pullbacks = [{tuple(e): c for e, c in f.dict().items()} for f in t.linear_pullbacks]
    unit = matrix(ZZ, [[D.get(tuple(e), 0) for D in pullbacks] for e in exponents])
    assert abs(unit.det()) == 1
    long_exponent = (1, 0, 0, 0)
    coefficient = lambda f: {tuple(e): c for e, c in f.dict().items()}.get(long_exponent, 0)
    long_coefficients = [coefficient(f) for f in t.theta]
    assert any(long_coefficients) and all(not coefficient(f) for f in s.theta)
    report = dict(
        scope='All-degree adjoint torus/Levi source and image scaling; all 74 maximal-parabolic image rows derived from the audited minimal tables.',
        linear_scale=3, adjoint_dimension=t.dimension,
        trace_gram=rows(t.trace_gram), linear_basis=rows(t.linear_basis_matrix),
        weights=[[list(w), m] for w, m in t.weight_multiplicities],
        certificate=cert.__dict__, linear_unit_minor_exponents=exponents,
        map_distinction=dict(exponent=list(long_exponent),
                             minimal_coefficients=['0']*4,
                             adjoint_coefficients=list(map(str, long_coefficients))),
        fresh_api_rows=[], direct_laurent_rows=[], cases=[])
    for oldcase in baseline['cases']:
        node = oldcase['parabolic']
        a = f4_character_source(node)
        b = f4_character_source(node, representation='adjoint')
        assert a._inputs == b._inputs and b._scale == 3*a._scale
        algebra = b.rational_algebra
        case = dict(parabolic=node, complex_dimension=oldcase['complex_dimension'], degrees=[])
        for old in oldcase['degrees']:
            d = old['degree']
            q = ZZ(3)**d
            H = matrix(ZZ, old['image_basis'])
            newH = q*H
            factors = list(map(int, newH.elementary_divisors()))
            assert factors == [q*c for c in old['smith_factors']]
            transition = H.change_ring(QQ).inverse()*newH
            assert transition == q*identity_matrix(H.nrows())
            case['degrees'].append(dict(
                degree=d, rank=H.nrows(), degree_scale=str(q),
                minimal_smith_factors=old['smith_factors'],
                adjoint_smith_factors=factors, adjoint_image_basis=rows(newH),
                minimal_mod_adjoint_factors=[int(q)]*H.nrows() if d else []))
            if d <= 4:
                ra, rb = a.degree(d), b.degree(d)
                assert rb.basis_matrix == q*ra.basis_matrix
                expected_source = 18**d*matrix(QQ, old['normalized_source_basis'])
                assert rb.basis_matrix == expected_source
                im = rb.image()
                assert list(im.smith_factors) == factors
                old_lattice = ra.image()
                comparison = im.compare(old_lattice)
                assert comparison.left_in_right and comparison.right_in_left == (d == 0)
                quotient = old_lattice.quotient_by(im)
                assert quotient.free_rank == 0
                assert list(quotient.torsion_factors) == case['degrees'][-1]['minimal_mod_adjoint_factors']
                report['fresh_api_rows'].append(dict(
                    parabolic=node, degree=d, source_basis=rows(rb.basis_matrix),
                    smith_factors=factors,
                    minimal_mod_adjoint_factors=list(map(int, quotient.torsion_factors))))
        for d in range(4):
            for mode in ('homogeneous', 'filtered'):
                degrees = range(d+1) if mode == 'filtered' else (d,)
                labels = [(e, m) for e in degrees for m in algebra.monomials(e)]
                fs = [t.pullback(algebra.substitute(m)) for _, m in labels]
                dictionaries = [f.dict() for f in fs]
                support = sorted(set().union(*(set(D) for D in dictionaries)))
                A = matrix(QQ, [[D.get(e, 0) for D in dictionaries] for e in support])
                B = block_diagonal_matrix([b.degree(e).basis_matrix for e in degrees],
                                          subdivide=False)
                certificate = certify_preimage(A, B)
                report['direct_laurent_rows'].append(dict(
                    parabolic=node, degree=d, mode=mode,
                    labels=[[e, list(map(int, m.exponents()[0]))] for e, m in labels],
                    support_rows=len(support), certificate=certificate))
        report['cases'].append(case)
        print('P', node, 'checked; elapsed', round(perf_counter()-start, 2), flush=True)
    files = [Path(__file__).resolve(), base_path, baseline_audit_path,
             ROOT/'experiments/e6_character_form.py']
    files += sorted((ROOT/'src/sage_borel').glob('*.py'))
    report.update(result='passed', derived_image_rows=74, fresh_api_degrees=20,
                  direct_laurent_lattices=32, seconds=round(perf_counter()-start, 3),
                  source_hash_semantics='SHA256 with CRLF normalized to LF.',
                  source_sha256={p.relative_to(ROOT).as_posix():
                    sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/f4-adjoint-comparison-results.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASSED', report['seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
