"""Python-only audit of E7 adjoint sources and complete image comparisons.

No Sage or runtime-library imports. Symmetric power sums are reconstructed
by the closed partition formula rather than the runtime Newton recurrence.
The Sage companion independently checks Borel columns by the ordinary ring.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, factorial, gcd
from pathlib import Path
from time import perf_counter, sleep
import json
from audit_e6_character_form import inverse, matmul, combine, multiply, smith_square, check_image_basis
from audit_f4_character_form import determinant, exact_integer_columns
from audit_e7_character_form import C, X, Z

ROOT = Path(__file__).resolve().parents[1]
ZERO = (0,)*7
ONE = {ZERO: 1}


def integer_matrix(A):
    assert all(F(c).denominator == 1 for r in A for c in r)
    return [[int(c) for c in r] for r in A]


def same_lattice(A, B):
    T = integer_matrix(matmul(inverse(A), B))
    assert abs(determinant(T)) == 1


def contained(A, B):
    return all(F(c).denominator == 1 for r in matmul(inverse(B), A) for c in r)


def quotient(A, B):
    T = integer_matrix(matmul(inverse(A), B))
    return dict(free_rank=0, torsion_factors=[x for x in smith_square(T) if x > 1],
                index=abs(determinant(T)))


def power(f, n):
    out = ONE
    for _ in range(n):
        out = multiply(out, f)
    return out


def variables():
    return [{tuple(int(i == j) for i in range(7)): 1} for j in range(7)]


def symmetric_inputs():
    """Girard partition sums and direct expansions of the 27 linear forms."""
    es = variables()
    E, z = es[0], es[6]
    def partitions(k, j=1, prefix=()):
        if j == 7:
            if k == 0:
                yield prefix
            return
        for m in range(k//j+1):
            yield from partitions(k-j*m, j+1, prefix+(m,))
    ps = [{ZERO: 6}]
    for k in range(1, 13):
        p = {}
        for m in partitions(k):
            length = sum(m)
            denominator = 1
            for a in m:
                denominator *= factorial(a)
            numerator = (-1)**(k+length)*k*factorial(length-1)
            assert numerator % denominator == 0
            p[m+(0,)] = numerator//denominator
        ps.append(p)
    plus = combine([E, z], [-4, -12])
    minus = combine([E, z], [2, 12])
    pair_shift = combine([E], [2])
    inputs = [combine([E, z], [1, 3])]
    for d in (2, 5, 6, 8, 9, 12):
        terms, factors = [], []
        for k in range(d+1):
            shifted = combine([power(plus, d-k), power(minus, d-k)], [1, 1])
            terms.append(multiply(ps[k], shifted))
            factors.append(comb(d, k)*6**k)
            pairs = combine([multiply(ps[j], ps[k-j]) for j in range(k+1)]+[ps[k]],
                            [comb(k, j) for j in range(k+1)]+[-2**k])
            assert all(c % 2 == 0 for c in pairs.values())
            pairs = {e: c//2 for e, c in pairs.items()}
            terms.append(multiply(pairs, power(pair_shift, d-k)))
            factors.append(comb(d, k)*(-6)**k)
        inputs.append(combine(terms, factors))
    return inputs


def monomial_function(inputs, size=400):
    @lru_cache(maxsize=size)
    def monomial(p):
        if not any(p):
            return ONE
        i = next(i for i, e in enumerate(p) if e)
        q = list(p)
        q[i] -= 1
        return multiply(inputs[i], monomial(tuple(q)))
    return monomial


def read_report(path):
    for _ in range(1800):
        try:
            return json.loads(path.read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            sleep(1)
    raise TimeoutError("Sage report is not available")


def completed_rows(path):
    seen, deadline = 0, perf_counter()+1800
    while True:
        data = read_report(path)
        assert data["result"] in ("running", "passed")
        for row in data["degrees"][seen:]:
            yield row
            seen += 1
        if data["result"] == "passed":
            assert seen == 28
            return
        if perf_counter() > deadline:
            raise TimeoutError("complete passed Sage report not available")
        sleep(1)


def main():
    start = perf_counter()
    path = ROOT/"experiments/e7-adjoint-comparison-results.json"
    data = read_report(path)
    roots, pending = {tuple(C[i][0] for i in range(7))}, [tuple(C[i][0] for i in range(7))]
    while pending:
        w = pending.pop()
        for i in range(7):
            v = tuple(w[j]-w[i]*C[j][i] for j in range(7))
            if v not in roots:
                roots.add(v)
                pending.append(v)
    assert len(roots) == 126 and determinant(C) == 2
    invC = inverse(C)
    integer_matrix([[2*c for c in r] for r in invC])
    assert any(c.denominator == 2 for c in invC[6])
    assert {tuple(w): m for w, m in data["weights"]} == dict([(w, 1) for w in roots]+[(ZERO, 7)])
    gram = [[sum(w[i]*w[j] for w in roots) for j in range(7)] for i in range(7)]
    assert gram == [[36*c for c in r] for r in C]
    assert gram == [[int(c) for c in r] for r in data["trace_gram"]]
    moments = [{w: w[i] for w in roots if w[i]} for i in range(7)]
    theta = [combine(moments, r) for r in inverse(gram)]
    assert all(sum(f.values()) == 0 for f in theta)
    assert [[sum(c*e[j] for e, c in f.items()) for j in range(7)] for f in theta] == [
        [int(i == j) for j in range(7)] for i in range(7)]
    bs = [combine([f], [36]) for f in theta]
    assert all(F(c).denominator == 1 for f in bs for c in f.values())
    assert [[f.get(tuple(e), 0) for f in bs] for e in data["unit_exponents"]] == [
        [int(i == j) for j in range(7)] for i in range(7)]
    cert = data["certificate"]
    assert abs(determinant(cert["basis_change"])) == 1
    transformed = [combine(bs, r) for r in cert["basis_change"]]
    leading = [max(f) for f in transformed]
    assert [list(e) for e in leading] == cert["cases"][0][1]
    assert [f[e] for f, e in zip(transformed, leading)] == [1]*7 == cert["cases"][0][2]
    leading_det = determinant(leading)
    assert leading_det and cert["certified"] and not cert["exceptional_primes"]
    # Verify the integral P basis and A5 permutation action in fundamental weights.
    chain = [0, 2, 3, 4, 5]
    v = [[C[6][j]+sum(C[k][j] for k in chain[i:]) for j in range(7)] for i in range(6)]
    u = [C[1][j]+sum(v[i][j] for i in range(3)) for j in range(7)]
    E = [sum(r[j] for r in v) for j in range(7)]
    z = [F(u[j]-E[j], 2) for j in range(7)]
    T = integer_matrix(v+[z])
    compression = data["compression_certificate"]
    assert T == [[int(c) for c in r] for r in compression["weight_basis_change"]]
    assert abs(determinant(T)) == 1
    assert [E[j]+3*z[j] for j in range(7)] == [0]*6+[1]
    for check in compression["permutation_checks"]:
        i = check["node"]-1
        def reflect(r):
            return [r[j]-r[i]*C[i][j] for j in range(7)]
        assert [v.index(reflect(r)) for r in v] == check["permutation"]
        assert reflect(z) == z
    assert [[6*r[j]-E[j] for j in range(7)] for r in v] == matmul(X, C)
    assert [-3*E[j]-12*z[j] for j in range(7)] == matmul([Z], C)[0]
    inputs = symmetric_inputs()
    assert inputs == [{tuple(e): F(c) for e, c in f} for f in compression["compressed_inputs"]]
    # Fresh direct Laurent coefficients, bounded through degree three.
    xs = [combine(moments, [F(c, 36) for c in r]) for r in X]
    zz = combine(moments, [F(c, 36) for c in Z])
    linears = [combine([x, zz], [1, sign]) for x in xs for sign in (1, -1)]
    linears += [combine([xs[i], xs[j]], [-1, -1]) for i in range(6) for j in range(i+1, 6)]
    quadratic = combine([multiply(f, f) for f in linears], [1]*27)
    laurent = monomial_function([theta[6], quadratic])
    for row in data["direct_laurent_rows"]:
        polys = [laurent(tuple(p)) for p in row["monomial_exponents"]]
        support = sorted(set().union(*(set(f) for f in polys)))
        assert len(support) == row["support_rows"]
        integral = exact_integer_columns(polys, [[F(c) for c in r] for r in row["source_basis"]])
        divisor = 0
        for minor in row["completeness_certificate"]["maximal_minors"]:
            det = abs(determinant([[f.get(support[i], 0) for f in integral] for i in minor["rows"]]))
            assert det == int(minor["determinant"])
            divisor = gcd(divisor, det)
        assert divisor == 1
    print("Torus, P coordinates, symmetric inputs and eight Laurent certificates passed", round(perf_counter()-start, 2), flush=True)
    basepath = ROOT/"experiments/minuscule-character-api-results.json"
    auditpath = ROOT/"experiments/e6-adjoint-api-audit.json"
    base, prior = json.loads(basepath.read_text()), json.loads(auditpath.read_text())
    assert base["result"] == prior["result"] == "passed"
    assert sha256(basepath.read_bytes().replace(b"\r\n", b"\n")).hexdigest() == prior["source_sha256"][basepath.relative_to(ROOT).as_posix()]
    oldrows = next(c["degrees"] for c in base["cases"] if c["kind"] == "E7")
    monomial = monomial_function(inputs)
    count = 0
    for row in completed_rows(path):
        d = row["degree"]
        B = [[F(c) for c in r] for r in row["source_basis"]]
        B0 = [[c/(36**d) for c in r] for r in B]
        polys = [monomial(tuple(p)) for p in row["monomial_exponents"]]
        integral = exact_integer_columns(polys, B0)
        support = set().union(*(set(f) for f in polys))
        representatives = {tuple(f.get(e, 0) for f in polys): e for e in support}
        coefficient_rows = sorted(representatives)
        assert len(coefficient_rows) == row["coefficient_rows"]
        divisor = 0
        for minor in row["completeness_certificate"]["maximal_minors"]:
            exponents = [representatives[coefficient_rows[i]] for i in minor["rows"]]
            N = [[f.get(e, 0) for f in integral] for e in exponents]
            det = abs(determinant(N))
            assert det == int(minor["determinant"])
            divisor = gcd(divisor, det)
        assert divisor == 1
        R = [[F(c) for c in r] for r in row["rational_image_matrix"]]
        image = integer_matrix(matmul(R, B))
        assert image == [[int(c) for c in r] for r in row["image_matrix"]]
        H = [[int(c) for c in r] for r in row["image_basis"]]
        a = [[int(c) for c in r] for r in row["minimal_image_basis"]]
        check_image_basis(image, H)
        assert smith_square(H) == row["smith_factors"]
        old = oldrows[d]
        oldH = [[int(c) for c in old["image_basis"][old["schubert_words"].index(w)]]
                for w in row["schubert_words"]]
        same_lattice(oldH, a)
        assert smith_square(a) == row["minimal_smith_factors"]
        S = [[int(c) for c in r] for r in row["sum_basis"]]
        I = [[int(c) for c in r] for r in row["intersection_basis"]]
        assert contained(a, H) == row["minimal_in_adjoint"]
        assert contained(H, a) == row["adjoint_in_minimal"]
        check_image_basis([ra+rb for ra, rb in zip(a, H)], S)
        assert contained(I, a) and contained(I, H)
        assert abs(determinant(I)*determinant(S)) == abs(determinant(a)*determinant(H))
        same_lattice([[3**d*c for c in r] for r in a], I)
        same_lattice([[F(c, 3**d) for c in r] for r in H], S)
        for name, large, small in (("sum_mod_minimal", S, a), ("sum_mod_adjoint", S, H),
                                  ("minimal_mod_intersection", a, I), ("adjoint_mod_intersection", H, I)):
            assert quotient(large, small) == row[name]
        assert row["sum_mod_adjoint"]["torsion_factors"] == ([3**d]*len(H) if d else [])
        assert all(2**d % c == 0 for c in row["sum_mod_minimal"]["torsion_factors"])
        count += 1
        print("Audited degree", d, "elapsed", round(perf_counter()-start, 2), flush=True)
    files = [path, basepath, auditpath, Path(__file__).resolve(),
             ROOT/"experiments/audit_e6_character_form.py",
             ROOT/"experiments/audit_e7_character_form.py", ROOT/"experiments/audit_f4_character_form.py"]
    result = dict(result="passed", engine="Python standard library only", root_count=126,
                  leading_exponent_determinant=leading_det, compressed_generators=7,
                  direct_laurent_lattices=8, adjoint_source_degrees=count, image_comparison_degrees=count,
                  scope="Fresh roots, trace, all-primes proof, integral P coordinates, Girard-formula symmetric inputs, all source coefficient lattices and image arithmetic. Borel columns independently checked by the separate Sage ordinary-ring calculation; historical minimal image tables reused.",
                  seconds=round(perf_counter()-start, 3),
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest() for p in files})
    (ROOT/"experiments/e7-adjoint-comparison-audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k != "source_sha256"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
