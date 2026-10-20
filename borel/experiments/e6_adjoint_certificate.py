"""Fresh adjoint torus certificate and direct Laurent intersections."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
from sage.all import QQ, ZZ, CartanMatrix, matrix, identity_matrix, block_diagonal_matrix
from sage_borel import e6_character_source
from e6_character_form import certify_preimage

ROOT = Path(__file__).resolve().parents[1]


def rows(M):
    return [[str(c) for c in r] for r in M]


def main():
    started = perf_counter()
    source = e6_character_source(representation="adjoint")
    s = source.springer_map
    C = CartanMatrix("E6")
    roots, pending = {tuple(C.column(0))}, [tuple(C.column(0))]
    while pending:
        w = pending.pop()
        for i in range(6):
            v = tuple(w[j]-w[i]*C[j, i] for j in range(6))
            if v not in roots:
                roots.add(v)
                pending.append(v)
    assert len(roots) == 72
    weights = {tuple(map(int, w)): 1 for w in roots}
    weights[(0,)*6] = 6
    assert dict(s.weight_multiplicities) == weights
    assert s.trace_gram == 24*C
    assert matrix(QQ, [[f.derivative(t)(*([1]*6)) for t in s.torus_ring.gens()]
                       for f in s.theta]) == identity_matrix(6)
    pullbacks = [{tuple(e): c for e, c in f.dict().items()} for f in s.linear_pullbacks]
    unit_exponents = [tuple(map(int, C.column(i))) for i in range(6)]
    assert matrix(QQ, [[D.get(e, 0) for D in pullbacks] for e in unit_exponents]) == identity_matrix(6)
    report = dict(result="running", weights=[[list(w), m] for w, m in s.weight_multiplicities],
                  trace_gram=rows(s.trace_gram), linear_basis=rows(s.linear_basis_matrix),
                  unit_exponents=unit_exponents, certificate=source.torus_certificate.__dict__,
                  direct_laurent_rows=[])
    algebra = source.rational_algebra
    for d in range(4):
        for mode in ("homogeneous", "filtered"):
            degrees = range(d+1) if mode == "filtered" else (d,)
            labels = [(e, m) for e in degrees for m in algebra.monomials(e)]
            fs = [s.pullback(algebra.substitute(m)) for e, m in labels]
            dictionaries = [{tuple(e): c for e, c in f.dict().items()} for f in fs]
            support = sorted(set().union(*(set(D) for D in dictionaries)))
            M = matrix(QQ, [[D.get(e, 0) for D in dictionaries] for e in support])
            B = block_diagonal_matrix([source.degree(e).basis_matrix for e in degrees],
                                      subdivide=False)
            cert = certify_preimage(M, B)
            report["direct_laurent_rows"].append(dict(
                degree=d, mode=mode, monomial_exponents=[list(map(int, m.exponents()[0])) for e, m in labels],
                source_basis=rows(B), support_rows=len(support), completeness_certificate=cert))
    files = [Path(__file__).resolve(), ROOT/"experiments/e6_character_form.py"]
    files += sorted((ROOT/"src/sage_borel").glob("*.py"))
    report.update(result="passed", seconds=round(perf_counter()-started, 3),
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest() for p in files})
    (ROOT/"experiments/e6-adjoint-certificate-results.json").write_text(json.dumps(report, indent=2)+"\n")
    print("PASSED", report["seconds"], source.torus_certificate, flush=True)


if __name__ == "__main__":
    main()
