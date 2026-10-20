"""Further representations and finite certificates for the trace-source theorem."""
from dataclasses import asdict
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from time import perf_counter

from sage.all import CartanMatrix, GF, QQ, ZZ, matrix
from sage_borel import SpringerCharacterMap, modular_character_relation

ROOT = Path(__file__).resolve().parents[1]
PRIMES = (2, 3, 5, 7, 11, 13)


def cases():
    out = []
    def add(kind, weights):
        out.extend((kind, (tuple(w),)) for w in weights)
    add("A1", [(k,) for k in range(1, 9)])
    add("A2", [(1,0), (0,1), (2,0), (1,1), (3,0), (2,1), (2,2)])
    add("B2", [(1,0), (0,1), (2,0), (1,1)])
    add("A3", [(1,1,0), (0,2,0)])
    add("G2", [(a,b) for a in range(4) for b in range(4) if 1 <= a+b <= 3]
               +[(4,0), (2,2)])
    for kind, n in (("F4",4), ("E6",6)):
        add(kind, [tuple(int(i == j) for i in range(n)) for j in range(n)])
    add("F4", [(0,0,0,2), (1,0,0,1)])
    add("E6", [(2,0,0,0,0,0), (1,0,0,0,0,1)])
    add("E7", [tuple(int(i == j) for i in range(7)) for j in (0,1,5,6)]
               +[(0,0,0,0,0,0,2), (1,0,0,0,0,0,1)])
    add("E8", [(0,0,0,0,0,0,0,1), (1,0,0,0,0,0,0,0)])
    out += [("G2", ((1,0),(0,1))), ("G2", ((1,0),)*3),
            ("G2", ((1,0),(0,0))),
            ("E6", ((1,0,0,0,0,0),(0,0,0,0,0,1))),
            ("E7", ((1,0,0,0,0,0,0),(0,0,0,0,0,0,1)))]
    return out


def rows(M):
    return [[str(c) for c in r] for r in M.rows()]


def saturation_certificate(C):
    """Gcd-one maximal minors certify a primitive column lattice."""
    selected = C.transpose().pivots()
    det = abs(ZZ(C.matrix_from_rows(selected).det()))
    minors = [dict(rows=list(selected), determinant=str(det), selection_prime=0)]
    divisor = det
    while divisor != 1:
        p = divisor.prime_divisors()[0]
        selected = C.change_ring(GF(p)).transpose().pivots()
        assert len(selected) == C.ncols()
        det = abs(ZZ(C.matrix_from_rows(selected).det()))
        assert det % p
        minors.append(dict(rows=list(selected), determinant=str(det), selection_prime=int(p)))
        divisor = ZZ(gcd(divisor, det))
    return dict(minor_gcd=1, maximal_minors=minors)


def main():
    started = perf_counter()
    report = dict(result="running", cases=[], bounded_searches=[])
    path = ROOT/"experiments/trace-source-survey-results.json"
    for kind, highest in cases():
        s = SpringerCharacterMap(kind, highest)
        cert = s.certify_trace_source()
        assert cert.certified
        dictionaries = [f.dict() for f in s.linear_pullbacks]
        support = sorted(set().union(*(set(f) for f in dictionaries)), reverse=True)
        C = matrix(ZZ, [[f.get(e,0) for f in dictionaries] for e in support])
        checks = []
        for p in PRIMES:
            A = C.transpose().change_ring(GF(p))
            selected = A.pivots()
            assert len(selected) == s.rank
            U = A.matrix_from_columns(selected).inverse()
            transformed = U*A
            leading = [next(i for i,c in enumerate(r) if c) for r in transformed.rows()]
            assert leading == list(selected)
            E = matrix(ZZ, [support[i] for i in leading])
            assert E.det()
            checks.append(dict(prime=p, pivots=list(selected), change=rows(U),
                               leading_exponents=[list(support[i]) for i in selected],
                               exponent_determinant=str(E.det())))
        report["cases"].append(dict(cartan_type=kind, highest_weights=highest,
            dimension=s.dimension, cartan_matrix=rows(CartanMatrix(kind)),
            weights=[[w,m] for w,m in s.weight_multiplicities],
            trace_gram=rows(s.trace_gram), linear_basis=rows(s.linear_basis_matrix),
            support=[list(e) for e in support], coefficients=rows(C),
            certificate=asdict(cert), saturation_certificate=saturation_certificate(C),
            prime_checks=checks))
        # Search actual Laurent expansions, separately from the structural test.
        if (kind, highest) in [("G2", ((2,0),)), ("G2", ((1,1),)),
                               ("G2", ((0,2),)), ("G2", ((2,1),)),
                               ("A2", ((2,1),)), ("G2", ((1,0),(0,1)))]:
            for p in (2,3,5):
                assert s.modular_obstruction(p, max_degree=3) is None
                report["bounded_searches"].append(dict(case=len(report["cases"])-1,
                    prime=p, max_degree=3, relation=None))
        path.write_text(json.dumps(report, indent=2)+"\n")
        print(kind, highest, "dimension", s.dimension, "passed", round(perf_counter()-started,2), flush=True)
    assert len(report["cases"]) == 59 and len(report["bounded_searches"]) == 18
    files = [Path(__file__).resolve(), ROOT/"src/sage_borel/characters.py",
             ROOT/"src/sage_borel/__init__.py", ROOT/"src/sage_borel/character_invariants.py",
             ROOT/"tests/test_trace_source.py", ROOT/"TRACE_SOURCE_THEOREM.md"]
    report.update(result="passed", representation_cases=59,
                  prime_adapted_checks=59*len(PRIMES), bounded_relation_searches=18,
                  seconds=round(perf_counter()-started,3),
                  scope="Finite tests support the displayed all-degree proof; they do not replace it.",
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n",b"\n")).hexdigest() for p in files})
    path.write_text(json.dumps(report, indent=2)+"\n")
    print("PASSED",report["seconds"],flush=True)


if __name__ == "__main__":
    main()
