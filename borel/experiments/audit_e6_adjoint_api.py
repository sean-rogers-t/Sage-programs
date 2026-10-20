"""Independent Python-only E6 adjoint proof and E6/E7 API data audit.

Reconstructs the adjoint roots, trace map and all-primes certificate, all
adjoint source coefficient lattices, eight direct Laurent intersections,
and image/comparison arithmetic. Borel columns are checked separately by
the Sage ordinary-ring oracle. Historical minimal source proofs are reused.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path
from time import perf_counter, sleep
import json
from audit_e6_character_form import (
    C, Y, inverse, matmul, combine, multiply, smith_square, check_image_basis)
from audit_f4_character_form import determinant, exact_integer_columns

ROOT = Path(__file__).resolve().parents[1]
ZERO = (0,)*6
ONE = {ZERO: 1}


def orbit(w):
    found, pending = {w}, [w]
    while pending:
        w = pending.pop()
        for i in range(6):
            v = tuple(w[j]-w[i]*C[j][i] for j in range(6))
            if v not in found:
                found.add(v)
                pending.append(v)
    return found


def integer_matrix(A):
    assert all(F(x).denominator == 1 for r in A for x in r)
    return [[int(x) for x in r] for r in A]


def same_lattice(A, B):
    if A == B:
        return
    T = integer_matrix(matmul(inverse(A), B))
    assert abs(determinant(T)) == 1


def contained(A, B):
    return all(F(x).denominator == 1 for r in matmul(inverse(B), A) for x in r)


def quotient(A, B):
    T = integer_matrix(matmul(inverse(A), B))
    fs = [x for x in smith_square(T) if x > 1]
    return dict(free_rank=0, torsion_factors=fs, index=abs(determinant(T)))


def generators(variables, low=False):
    ys = [combine(variables, r) for r in Y]
    squares = [multiply(y, y) for y in ys[1:]]
    def elementary(k):
        terms = []
        for subset in combinations(squares, k):
            p = ONE
            for f in subset:
                p = multiply(p, f)
            terms.append(p)
        return combine(terms, [1]*len(terms))
    if low:
        return [ys[0], elementary(1)]
    euler = ONE
    for y in ys[1:]:
        euler = multiply(euler, y)
    return [ys[0], elementary(1), elementary(2), euler, elementary(3), elementary(4)]


def monomials_for(inputs):
    @lru_cache(maxsize=512)
    def monomial(p):
        if not any(p):
            return ONE
        i = next(i for i, e in enumerate(p) if e)
        q = list(p)
        q[i] -= 1
        return multiply(inputs[i], monomial(tuple(q)))
    return monomial


def check_minors(N, certificate):
    divisor = 0
    for minor in certificate["maximal_minors"]:
        det = abs(determinant([N[i] for i in minor["rows"]]))
        assert det == int(minor["determinant"])
        divisor = gcd(divisor, det)
    assert divisor == certificate["minor_gcd"] == 1


def completed_cases(path):
    """Audit finished cases while the independent Sage run finishes E7.

    Final success requires the complete passed report. A partial checkpoint
    never produces a full audit result; stop after 30 minutes if unfinished.
    """
    seen, deadline = 0, perf_counter()+1800
    while True:
        try:
            data = json.loads(path.read_text())
        except json.JSONDecodeError:  # a checkpoint write may be in progress
            data = {"result": "running", "cases": []}
        assert data["result"] in ("running", "passed")
        for case in data["cases"][seen:]:
            yield case
            seen += 1
        if data["result"] == "passed":
            assert seen == 3
            return
        if perf_counter() > deadline:
            raise TimeoutError("the complete passed API report is not available")
        sleep(1)


def main():
    start = perf_counter()
    torus_path = ROOT/"experiments/e6-adjoint-certificate-results.json"
    path = ROOT/"experiments/minuscule-character-api-results.json"
    torus = json.loads(torus_path.read_text())
    assert torus["result"] == "passed"
    roots = orbit(tuple(C[i][0] for i in range(6)))
    assert len(roots) == 72
    assert determinant(C) == 3
    inverse_C = inverse(C)
    integer_matrix([[3*c for c in r] for r in inverse_C])
    assert any(r[0].denominator == 3 for r in inverse_C)
    assert {tuple(w): m for w, m in torus["weights"]} == dict(
        [(w, 1) for w in roots]+[(ZERO, 6)])
    gram = [[sum(w[i]*w[j] for w in roots) for j in range(6)] for i in range(6)]
    assert gram == [[24*c for c in r] for r in C]
    assert gram == [[int(c) for c in r] for r in torus["trace_gram"]]
    moments = [{w: w[i] for w in roots if w[i]} for i in range(6)]
    theta = [combine(moments, r) for r in inverse(gram)]
    assert all(sum(f.values()) == 0 for f in theta)
    assert [[sum(c*e[j] for e, c in f.items()) for j in range(6)] for f in theta] == [
        [int(i == j) for j in range(6)] for i in range(6)]
    bs = [combine([f], [24]) for f in theta]
    assert all(F(c).denominator == 1 for f in bs for c in f.values())
    unit = [[f.get(tuple(e), 0) for f in bs] for e in torus["unit_exponents"]]
    assert unit == [[int(i == j) for j in range(6)] for i in range(6)]
    cert = torus["certificate"]
    assert abs(determinant(cert["basis_change"])) == 1
    fs = [combine(bs, r) for r in cert["basis_change"]]
    leading = [max(f) for f in fs]
    assert [list(e) for e in leading] == cert["cases"][0][1]
    assert [f[e] for f, e in zip(fs, leading)] == [1]*6 == cert["cases"][0][2]
    assert cert["certified"] and cert["exceptional_primes"] == []
    leading_det = determinant(leading)
    assert leading_det
    laurent_monomial = monomials_for(generators(theta, low=True))
    for row in torus["direct_laurent_rows"]:
        polys = [laurent_monomial(tuple(p)) for p in row["monomial_exponents"]]
        support = sorted(set().union(*(set(p) for p in polys)))
        assert len(support) == row["support_rows"]
        B = [[F(c) for c in r] for r in row["source_basis"]]
        integral = exact_integer_columns(polys, B)
        N = [[f.get(e, 0) for f in integral] for e in support]
        check_minors(N, row["completeness_certificate"])
    print("Torus and 8 direct Laurent certificates passed", round(perf_counter()-start, 2), flush=True)
    variables = [{tuple(int(i == j) for i in range(6)): 1} for j in range(6)]
    polynomial_monomial = monomials_for(generators(variables))
    minimal_images = {}
    counts = {"image_rows": 0, "minimal_source_lattices": 0, "adjoint_source_lattices": 0}
    files = [path, torus_path, Path(__file__).resolve(),
             ROOT/"experiments/audit_e6_character_form.py",
             ROOT/"experiments/audit_f4_character_form.py"]
    for case in completed_cases(path):
        kind, representation = case["kind"], case["representation"]
        baseline_path = ROOT/("experiments/"+kind.lower()+"-character-form-results.json")
        baseline_audit = ROOT/("experiments/"+kind.lower()+"-character-form-audit.json")
        old, previous = json.loads(baseline_path.read_text()), json.loads(baseline_audit.read_text())
        assert old["result"] == previous["result"] == "passed"
        assert sha256(baseline_path.read_bytes().replace(b"\r\n", b"\n")).hexdigest() == previous["source_sha256"][baseline_path.relative_to(ROOT).as_posix()]
        files += [baseline_path, baseline_audit]
        for row in case["degrees"]:
            d = row["degree"]
            B = [[F(c) for c in r] for r in row["source_basis"]]
            R = [[F(c) for c in r] for r in row["rational_image_matrix"]]
            image = integer_matrix(matmul(R, B))
            assert image == [[int(c) for c in r] for r in row["image_matrix"]]
            H = [[int(c) for c in r] for r in row["image_basis"]]
            check_image_basis(image, H)
            assert smith_square(H) == row["smith_factors"]
            counts["image_rows"] += 1
            if representation == "minimal":
                prior = old["degrees"][d]
                assert row["monomial_exponents"] == prior["monomial_exponents"]
                same_lattice([[F(c) for c in r] for r in prior["source_basis_matrix"]], B)
                oldH = [prior["image_basis"][prior["schubert_words"].index(w)] for w in row["schubert_words"]]
                same_lattice(oldH, H)
                if kind == "E6":
                    minimal_images[d] = H
                counts["minimal_source_lattices"] += 1
            else:
                polys = [polynomial_monomial(tuple(p)) for p in row["monomial_exponents"]]
                B0 = [[c/(24**d) for c in r] for r in B]
                integral = exact_integer_columns(polys, B0)
                support = set().union(*(set(p) for p in polys))
                coefficient_rows = sorted(set(tuple(p.get(e, 0) for p in polys) for e in support))
                assert len(coefficient_rows) == row["coefficient_rows"]
                # All coefficients are integral by exact_integer_columns;
                # recompute the selected minors in the experiment's row order.
                for minor in row["completeness_certificate"]["maximal_minors"]:
                    selected = [coefficient_rows[i] for i in minor["rows"]]
                    N = integer_matrix(matmul(selected, B0))
                    assert abs(determinant(N)) == int(minor["determinant"])
                assert gcd(*(int(m["determinant"]) for m in row["completeness_certificate"]["maximal_minors"])) == 1
                a = minimal_images[d]
                c = row["comparison"]
                S = [[int(x) for x in r] for r in c["sum_basis"]]
                I = [[int(x) for x in r] for r in c["intersection_basis"]]
                assert contained(a, H) == c["minimal_in_adjoint"]
                assert contained(H, a) == c["adjoint_in_minimal"]
                check_image_basis([ra+rb for ra, rb in zip(a, H)], S)
                assert contained(I, a) and contained(I, H)
                assert abs(determinant(I)*determinant(S)) == abs(determinant(a)*determinant(H))
                same_lattice([[4**d*x for x in r] for r in a], I)
                same_lattice([[F(x, 4**d) for x in r] for r in H], S)
                for name, large, small in (("sum_mod_minimal", S, a), ("sum_mod_adjoint", S, H),
                    ("minimal_mod_intersection", a, I), ("adjoint_mod_intersection", H, I)):
                    assert quotient(large, small) == c[name]
                # Consequences of Q_(2)=P_(2) and 4 being a 3-local unit.
                old_factors = smith_square(a)
                def valuation(n, p):
                    v = 0
                    while n % p == 0:
                        n //= p
                        v += 1
                    return v
                assert [valuation(x, 2) for x in row["smith_factors"]] == [
                    valuation(x, 2)+2*d for x in old_factors]
                assert all(x % 3 != 0 for x in c["sum_mod_adjoint"]["torsion_factors"])
                counts["adjoint_source_lattices"] += 1
                print("Adjoint degree", d, "passed", round(perf_counter()-start, 2), flush=True)
    assert counts == dict(image_rows=62, minimal_source_lattices=45, adjoint_source_lattices=17)
    result = dict(result="passed", engine="Python standard library only",
                  leading_exponent_determinant=leading_det, direct_laurent_lattices=8,
                  **counts, seconds=round(perf_counter()-start, 3),
                  scope="Fresh adjoint torus and source coefficients; all image and comparison arithmetic. Minimal source bases compared with historically audited tables. Independent Borel evaluation is the separate Sage published-ring check, with its documented BGG anchor dependencies.",
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest() for p in files})
    (ROOT/"experiments/e6-adjoint-api-audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k != "source_sha256"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
