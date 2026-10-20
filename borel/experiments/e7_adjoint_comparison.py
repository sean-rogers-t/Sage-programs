"""E7 adjoint full source, compressed weight model and complete P7 images."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from unittest.mock import patch
from itertools import combinations
import json
from sage.all import (QQ, ZZ, CartanMatrix, PolynomialRing, matrix, vector,
                      prod, identity_matrix, block_diagonal_matrix)
from sage_borel import e7_character_source, IntegralImageLattice
from sage_borel.character_invariants import _integral_preimage
from e6_character_form import certify_preimage
from exceptional_published_rings import PublishedMinusculeRing

ROOT = Path(__file__).resolve().parents[1]


def rows(M):
    return [[str(c) for c in r] for r in M]


def compression_certificate(source, C):
    A = PolynomialRing(QQ, ["w"+str(i) for i in range(1, 8)])
    w = vector(A, A.gens())
    a = C*w
    chain = [0, 2, 3, 4, 5]
    vs = [a[6]+sum(a[j] for j in chain[i:]) for i in range(6)]
    u = a[1]+sum(vs[:3])
    z = (u-sum(vs))/2
    T = matrix(ZZ, [[f.monomial_coefficient(x) for x in w] for f in vs+[z]])
    assert abs(T.det()) == 1
    assert sum(vs)+3*z == w[6]
    permutations = []
    for i in chain:
        reflect = lambda f: f(*[w[j]-int(i == j)*a[i] for j in range(7)])
        perm = [vs.index(reflect(v)) for v in vs]
        assert reflect(z) == z
        permutations.append(dict(node=i+1, permutation=perm))
    definition = json.loads((ROOT/"src/sage_borel/data/exceptional_generators.json").read_text())["conventions"]["psi_e7"]
    for v, r in zip(vs, definition["x_simple_root_coefficients"]):
        assert 6*v-sum(vs) == sum(c*x for c, x in zip(r, a))
    assert -3*sum(vs)-12*z == sum(c*x for c, x in zip(definition["z_simple_root_coefficients"], a))
    V = PolynomialRing(QQ, ["v"+str(i) for i in range(1, 7)]+["z"])
    v, zz = V.gens()[:6], V.gen(6)
    E = sum(v)
    es = [sum(prod(v[i] for i in I) for I in combinations(range(6), k)) for k in range(1, 7)]
    xx = [6*x-E for x in v]
    Z = -3*E-12*zz
    assert source._inputs[0](*(es+[zz])) == E+3*zz
    for d, f in zip((2, 5, 6, 8, 9, 12), source._inputs[1:]):
        direct = sum((x+Z)**d+(x-Z)**d for x in xx)
        direct += sum((-xx[i]-xx[j])**d for i in range(6) for j in range(i+1, 6))
        assert f(*(es+[zz])) == direct
    return dict(weight_basis_change=rows(T), determinant=int(T.det()),
                permutation_checks=permutations,
                compressed_inputs=[[[list(map(int, e)), str(c)] for e, c in sorted(f.dict().items())]
                                   for f in source._inputs])


def main():
    start = perf_counter()
    output = ROOT/"experiments/e7-adjoint-comparison-results.json"
    baseline_path = ROOT/"experiments/minuscule-character-api-results.json"
    baseline_audit = ROOT/"experiments/e6-adjoint-api-audit.json"
    base, prior = json.loads(baseline_path.read_text()), json.loads(baseline_audit.read_text())
    assert base["result"] == prior["result"] == "passed"
    assert sha256(baseline_path.read_bytes().replace(b"\r\n", b"\n")).hexdigest() == prior["source_sha256"][baseline_path.relative_to(ROOT).as_posix()]
    minimal_rows = next(c["degrees"] for c in base["cases"] if c["kind"] == "E7")
    source = e7_character_source(7, representation="adjoint")
    s, algebra = source.springer_map, source.rational_algebra
    C = CartanMatrix("E7")
    roots, pending = {tuple(C.column(0))}, [tuple(C.column(0))]
    while pending:
        w = pending.pop()
        for i in range(7):
            v = tuple(w[j]-w[i]*C[j, i] for j in range(7))
            if v not in roots:
                roots.add(v)
                pending.append(v)
    assert len(roots) == 126
    weights = {tuple(map(int, w)): 1 for w in roots}
    weights[(0,)*7] = 7
    assert dict(s.weight_multiplicities) == weights and s.dimension == 133
    assert s.trace_gram == 36*C and s.linear_basis_matrix == 36*identity_matrix(7)
    assert matrix(QQ, [[f.derivative(t)(*([1]*7)) for t in s.torus_ring.gens()] for f in s.theta]) == identity_matrix(7)
    unit_exponents = [tuple(map(int, C.column(i))) for i in range(7)]
    bs = [{tuple(e): c for e, c in f.dict().items()} for f in s.linear_pullbacks]
    assert matrix(QQ, [[D.get(e, 0) for D in bs] for e in unit_exponents]) == identity_matrix(7)
    report = dict(result="running", weights=[[list(w), m] for w, m in s.weight_multiplicities],
                  trace_gram=rows(s.trace_gram), linear_basis=rows(s.linear_basis_matrix),
                  certificate=source.torus_certificate.__dict__, unit_exponents=unit_exponents,
                  compression_certificate=compression_certificate(source, C),
                  direct_laurent_rows=[], degrees=[])
    for d in range(4):
        for mode in ("homogeneous", "filtered"):
            ds = range(d+1) if mode == "filtered" else (d,)
            labels = [(e, m) for e in ds for m in algebra.monomials(e)]
            fs = [s.pullback(algebra.substitute(m)) for e, m in labels]
            dictionaries = [{tuple(e): c for e, c in f.dict().items()} for f in fs]
            support = sorted(set().union(*(set(D) for D in dictionaries)))
            M = matrix(QQ, [[D.get(e, 0) for D in dictionaries] for e in support])
            B = block_diagonal_matrix([source.degree(e).basis_matrix for e in ds], subdivide=False)
            report["direct_laurent_rows"].append(dict(
                degree=d, mode=mode, monomial_exponents=[list(map(int, m.exponents()[0])) for e, m in labels],
                source_basis=rows(B), support_rows=len(support), completeness_certificate=certify_preimage(M, B)))
    # Clear bounded source caches so every full-row calculation captures its
    # actual coefficient constraints, including degrees used by Laurent checks.
    source._cache.clear()
    output.write_text(json.dumps(report, indent=2)+"\n")
    print("Torus, compression and direct Laurent checks passed", round(perf_counter()-start, 2), flush=True)
    pub = PublishedMinusculeRing("E7")
    minimal = e7_character_source(7)
    for old in minimal_rows:
        d = old["degree"]
        captured = []
        def preimage(M):
            captured.append(M)
            return _integral_preimage(M)
        with patch("sage_borel.character_invariants._integral_preimage", side_effect=preimage):
            row = source.degree(d, max_matrix_entries=10000000)
        assert len(captured) == 1
        M = captured[0]
        cert = certify_preimage(M, row.basis_matrix/(36**d))
        im = row.image()
        rational = algebra.degree(d)
        words = [list(w.reduced_word()) for w in rational.schubert_indices]
        powers = [list(map(int, m.exponents()[0])) for m in row.monomials]
        assert powers == old["monomial_exponents"]
        columns = []
        for p in powers:
            polynomial = prod(pub.generator_expressions[n]**e for n, e in zip(algebra.generator_names, p))
            check = {algebra.borel_map.W.from_reduced_word(w.reduced_word()): c
                     for w, c in pub.expand(polynomial, d).items()}
            columns.append([check.get(w, 0) for w in rational.schubert_indices])
        assert rational.matrix == matrix(QQ, columns).transpose()
        prior_R = matrix(QQ, [old["rational_image_matrix"][old["schubert_words"].index(w)] for w in words])
        assert prior_R == rational.matrix
        prior_H = matrix(ZZ, [old["image_basis"][old["schubert_words"].index(w)] for w in words])
        a = IntegralImageLattice._from_columns(d, rational.schubert_indices, prior_H, algebra.borel_map)
        comparison = a.compare(im)
        if d <= 5:
            fresh = minimal.degree(d).image().compare(a)
            assert fresh.left_in_right and fresh.right_in_left
        transition = matrix(QQ, comparison.intersection.basis).inverse()*(3**d*a.basis)
        assert transition.denominator() == 1 and abs(transition.det()) == 1
        transition = matrix(QQ, comparison.sum.basis).inverse()*(im.basis/(3**d))
        assert transition.denominator() == 1 and abs(transition.det()) == 1
        entry = dict(degree=d, source_rank=len(powers), monomial_exponents=powers,
            source_basis=rows(row.basis_matrix), coefficient_rows=M.nrows(), completeness_certificate=cert,
            schubert_words=words, rational_image_matrix=rows(rational.matrix),
            image_matrix=rows(rational.matrix*row.basis_matrix), image_basis=rows(im.basis),
            smith_factors=list(map(int, im.smith_factors)), minimal_image_basis=rows(a.basis),
            minimal_smith_factors=list(map(int, a.smith_factors)),
            minimal_in_adjoint=comparison.left_in_right, adjoint_in_minimal=comparison.right_in_left,
            sum_basis=rows(comparison.sum.basis), intersection_basis=rows(comparison.intersection.basis),
            sum_mod_minimal=comparison.sum.quotient_by(a).as_dict(),
            sum_mod_adjoint=comparison.sum.quotient_by(im).as_dict(),
            minimal_mod_intersection=a.quotient_by(comparison.intersection).as_dict(),
            adjoint_mod_intersection=im.quotient_by(comparison.intersection).as_dict())
        report["degrees"].append(entry)
        output.write_text(json.dumps(report, indent=2)+"\n")
        print("degree", d, "rank", len(powers), "quotient", entry["sum_mod_minimal"]["torsion_factors"],
              "elapsed", round(perf_counter()-start, 2), flush=True)
    files = [Path(__file__).resolve(), baseline_path, baseline_audit,
             ROOT/"experiments/e6_character_form.py", ROOT/"experiments/exceptional_published_rings.py",
             ROOT/"tests/test_exceptional.py", ROOT/"src/sage_borel/data/exceptional_generators.json"]
    files += sorted((ROOT/"src/sage_borel").glob("*.py"))
    assert len(report["degrees"]) == 28 and sum(r["source_rank"] for r in report["degrees"]) == 1967
    report.update(result="passed", fresh_adjoint_degrees=28, independently_checked_borel_columns=1967,
                  fresh_minimal_api_degrees=6, seconds=round(perf_counter()-start, 3),
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest() for p in files})
    output.write_text(json.dumps(report, indent=2)+"\n")
    print("PASSED", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
