"""Finite exact B/D row comparisons on R plus e_H times the smaller rectangle.

Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_extended_chern_comparison.py
Here R=H*(Gr(H,H+K);Z) and the second source is e_H times
H*(Gr(H,H+K-1);Z). This is the full D source after swapping the two
coordinate blocks, when the original even tail has size H. Equality
with the B coefficient-row lattice is checked in finite cases; the general
comparison is proved separately in TYPE_D_EVEN_GLUING_THEOREM.md.
The comparison on arbitrary surviving Chern coordinates is false; an
explicit six-column counterexample is retained as a scope guard.
"""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import signal
from time import perf_counter

from sage.all import (Partition, Partitions, PolynomialRing, QQ, ZZ,
                      diagonal_matrix, matrix, prod)
from sage.env import SAGE_VERSION
from sage_borel import classical_borel
from sage_borel.operations import elementary_polynomial

from type_d_chern_frontier import BinaryChernModel, SCHUR
from type_b_cartesian_frontier import CartesianCFrontier
from type_d_cartesian_frontier import CartesianDFrontier
from verify_type_d_three_row_smith import integral_inverse_product
from verify_type_d_even_gluing import matrix_exponents
from verify_type_d_square_comparison import schur_polynomials


def coefficient_matrix(model, degree, kind, source):
    """Actual integral Schubert coefficients of arbitrary specified Chern columns."""
    coordinates = model.coordinates(degree, kind)
    shapes = [tuple(partition) for partition in
              Partitions(degree, max_part=model.maximum-(kind == "D"), max_length=model.a)
              if all(list(partition).count(part) <= 1
                     for part in partition if part > model.k)]
    values = [model.decompose(model.theta(partition, kind)) for partition in shapes]
    change = matrix(ZZ, len(coordinates), len(shapes), lambda i, j:
                    values[j].get(coordinates[i][0], SCHUR.zero()).coefficient(Partition(coordinates[i][1])))
    assert change.is_square() and abs(change.det()) == 1
    assert all(value in coordinates for value in source)
    inputs = matrix(ZZ, len(coordinates), len(source), lambda i, j: int(coordinates[i] == source[j]))
    coefficients = integral_inverse_product(change, inputs)
    actual = diagonal_matrix(ZZ, [2**sum(part > model.k for part in partition)
                                 for partition in shapes])*coefficients
    assert actual.base_ring() is ZZ
    return actual


def extended_source(model, degree):
    H, K = model.a, model.k
    zero, top = (0,)*H, (0,)*(H-1)+(1,)
    source = [(zero, partition) for partition in model.rectangle_partitions(degree//2)]
    if degree >= H and (degree-H) % 2 == 0:
        source += [(top, tuple(partition)) for partition in
                   Partitions((degree-H)//2, max_part=K-1, max_length=H)]
    return source


def cartesian_extended_matrices(H, K, degree, cache):
    """Independent Cartesian D BGG and Cartesian C BGG rescaled to B."""
    if (H, K) not in cache:
        n = H+K
        c, d = classical_borel("C", n, K), classical_borel("D", n, K)
        ring = PolynomialRing(QQ, n, names=["x%d" % i for i in range(n)])
        cache[H, K] = CartesianCFrontier(c), CartesianDFrontier(d), ring
    c_frontier, d_frontier, ring = cache[H, K]
    x = ring.gens()
    polynomials = []
    for sector in (0, 1):
        remaining = degree-sector*H
        if remaining < 0 or remaining % 2:
            continue
        formal_ring, _partitions, formal = schur_polynomials(K-sector, H, remaining//2)
        substitution = formal_ring.hom([
            elementary_polynomial(ring, [value*value for value in x[:H]], i)
            for i in range(1, H+1)], ring)
        factor = prod(x[:H]) if sector else ring.one()
        polynomials.extend(factor*substitution(value) for value in formal)
    columns = [[], []]
    for polynomial in polynomials:
        c_value = c_frontier.expand(polynomial)
        columns[0].append({word: 2**sum(entry < 0 for entry in word.matrix().list())*value
                           for word, value in c_value.items()})
        columns[1].append(d_frontier.expand(polynomial))
    answers = []
    for sector in columns:
        words = tuple({word for column in sector for word in column})
        answers.append(matrix(ZZ, len(words), len(sector),
                              lambda i, j: sector[j].get(words[i], 0)))
    return answers


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, here.with_name("type_d_chern_frontier.py"),
              here.with_name("type_b_cartesian_frontier.py"), here.with_name("type_d_cartesian_frontier.py"),
              here.with_name("verify_type_d_three_row_smith.py"), here.with_name("verify_type_d_even_gluing.py"),
              here.with_name("verify_type_d_candidate.py"),
              here.with_name("verify_type_d_square_comparison.py"),
              borel/"src"/"sage_borel"/"core.py", borel/"src"/"sage_borel"/"operations.py"]
    hashes = {str(path.relative_to(borel)): sha256(path.read_bytes()).hexdigest() for path in inputs}
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "source_sha256": hashes,
              "scope": "Finite exact ZZ evidence for the special source R+e_H R0. The general comparison is proved separately in TYPE_D_EVEN_GLUING_THEOREM.md; these cases do not replace the proof.",
              "resource_limits": {"component_seconds": 90, "total_seconds": 600},
              "completed_components": [], "excluded_components": []}
    output = here.with_name("type_d_extended_chern_comparison_verification.json")

    def save():
        report["seconds"] = round(perf_counter()-started, 6)
        payload = json.dumps(report, indent=2)+"\n"
        try:
            output.write_text(payload)
        except OSError:
            Path("/tmp/type_d_extended_chern_comparison_verification.json").write_text(payload)

    # A necessary guard against replacing this source by every Chern class.
    guard_model = BinaryChernModel(3, 2)
    guard_source = guard_model.coordinates(8, "D")
    guard_B = coefficient_matrix(guard_model, 8, "C", guard_source)
    guard_D = coefficient_matrix(guard_model, 8, "D", guard_source)
    assert guard_B.row_module() != guard_D.row_module()
    assert matrix_exponents(guard_B) == [1, 1, 1, 1, 1, 2]
    assert matrix_exponents(guard_D) == [1, 1, 1, 1, 2, 2]
    report["arbitrary_surviving_Chern_source_counterexample"] = {
        "H": 3, "K": 2, "polynomial_degree": 8,
        "source_binary_Chern_exponents_and_square_Schur_partitions": [
            [list(map(int, basis)), list(map(int, partition))] for basis, partition in guard_source],
        "B_exponents": matrix_exponents(guard_B), "D_exponents": matrix_exponents(guard_D),
        "B_Hermite_row_basis": [list(map(int, row)) for row in guard_B.row_module().basis_matrix()],
        "D_Hermite_row_basis": [list(map(int, row)) for row in guard_D.row_module().basis_matrix()],
        "ZZ_row_modules_differ": True}
    save()

    def timeout(_signal, _frame):
        raise TimeoutError("component wall-clock bound reached")

    signal.signal(signal.SIGALRM, timeout)
    cases = [(H, K, degree) for H, K in ((4, 3), (4, 4), (4, 5), (6, 2))
             for degree in range(0, 2*H*K+1, 2)]
    cases += [(6, 3, degree) for degree in range(0, 25, 2)]
    cases += [(6, 4, 20), (6, 4, 24)]
    models, cartesian_cache = {}, {}
    for H, K, degree in cases:
        stamp = perf_counter()
        remaining = int(600-(stamp-started))
        if remaining <= 0:
            report["excluded_components"].append({"H": H, "K": K, "polynomial_degree": degree,
                                                   "reason": "overall wall-clock bound"})
            continue
        signal.alarm(min(90, remaining))
        try:
            if (H, K) not in models:
                models[H, K] = BinaryChernModel(H, K)
            model = models[H, K]
            source = extended_source(model, degree)
            B, D = [coefficient_matrix(model, degree, kind, source) for kind in ("C", "D")]
            equal = B.row_module() == D.row_module()
            cross = (H, K, degree) in ((4, 3, 8), (4, 3, 12), (4, 4, 12))
            if cross:
                actual_B, actual_D = cartesian_extended_matrices(H, K, degree, cartesian_cache)
                assert actual_B.row_module() == B.row_module(), (H, K, degree, "B")
                assert actual_D.row_module() == D.row_module(), (H, K, degree, "D")
            row = {"H": H, "K": K, "original_D_head": K, "original_D_even_tail": H,
                   "polynomial_degree": degree, "source_columns": len(source),
                   "exact_ZZ_B_D_row_module_equality": equal,
                   "B_exponents": matrix_exponents(B), "D_exponents": matrix_exponents(D),
                   "independent_Cartesian_B_and_D_row_module_crosschecks": cross,
                   "seconds": round(perf_counter()-stamp, 6), "base_ring": "ZZ"}
            report["completed_components"].append(row)
            print("H=%d,K=%d,d=%d: exact ZZ B/D rows equal=%s" % (H, K, degree, equal), flush=True)
        except TimeoutError as error:
            report["excluded_components"].append({"H": H, "K": K, "polynomial_degree": degree,
                                                   "reason": str(error), "seconds": round(perf_counter()-stamp, 6)})
        finally:
            signal.alarm(0)
        save()
    assert all(sha256(path.read_bytes()).hexdigest() == hashes[str(path.relative_to(borel))] for path in inputs)
    report["proof_snapshot_sha256_at_finalization"] = {
        name: sha256((borel/name).read_bytes()).hexdigest()
        for name in ("TYPE_D_CANDIDATE.md", "TYPE_D_EVEN_GLUING_THEOREM.md",
                     "TYPE_D_RELATIVE_EULER_EXTENSION.md", "TYPE_D_HEAD_EULER_B_COMPARISON.md")}
    rows = report["completed_components"]
    report["summary"] = {
        "exact_ZZ_B_D_extended_source_comparisons": len(rows),
        "unequal_ZZ_row_modules": sum(not row["exact_ZZ_B_D_row_module_equality"] for row in rows),
        "independent_Cartesian_B_and_D_row_module_crosscheck_pairs": sum(
            row["independent_Cartesian_B_and_D_row_module_crosschecks"] for row in rows),
        "excluded_components": len(report["excluded_components"]),
        "arbitrary_Chern_source_false_generalization_guard_passed": True}
    save()
    print(report["summary"], flush=True)


if __name__ == "__main__":
    main()
