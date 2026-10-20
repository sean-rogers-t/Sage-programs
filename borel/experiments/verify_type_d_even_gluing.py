"""Exact finite checks of full even-tail type-D images and saturated gluing.

Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_even_gluing.py
The full gluing dimension is computed from THREE integer Smith indices;
it is not replaced by the overlap of the original reductions modulo two.
Finite agreement with unit-gluing exhaustion does not prove that conjecture.
"""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import signal
from time import perf_counter

from sage.all import GF, PolynomialRing, QQ, ZZ, diagonal_matrix, gcd, matrix, prod
from sage.env import SAGE_VERSION
from sage_borel import classical_borel
from sage_borel.operations import elementary_polynomial
from type_d_cartesian_frontier import CartesianDFrontier
from verify_type_d_three_row_smith import ThreeRootChernModel, integral_inverse_product, smith_factors
from verify_type_d_two_row_smith import TwoRootChernModel
from verify_type_d_square_comparison import schur_polynomials
from verify_type_d_candidate import type_b_counts, expanded, unit_intersection, full_candidate


def three_head_square(K, degree, model):
    """Partition-indexed actual coefficients on the primitive Y_lambda basis."""
    partitions = model.isotropic_partitions(degree, model.maximum-1)
    coordinates = model.coordinates(degree)
    ideal = model.coordinates(degree, ideal=True)
    coordinates = [item for item in coordinates if item not in ideal]
    source = [(basis, p) for basis, p in coordinates if not any(basis)]
    values = [model.decompose(model.theta(partition, "D")) for partition in partitions]
    change = matrix(ZZ, [model.coordinate_vector(value, coordinates)
                         for value in values]).transpose()
    assert change.is_square() and abs(change.det()) == 1
    inputs = matrix(ZZ, len(coordinates), len(source),
                    lambda i, j: int(coordinates[i] == source[j]))
    coefficients = integral_inverse_product(change, inputs)
    actual = diagonal_matrix(ZZ, [2**sum(part > K for part in partition)
                                  for partition in partitions])*coefficients
    return partitions, actual, [partition for _basis, partition in source]


def two_head_euler(K, degree, model):
    """B(2,K) coefficients; Proposition 2 identifies these with D differences."""
    if degree < 0 or degree % 2:
        return [], matrix(ZZ, 0, 0), []
    m = degree//2
    source = model.rectangle_partitions(m)
    e_source = model.rectangle_partitions(m-1)
    partitions = model.isotropic_partitions(degree)
    values = [model.theta(partition) for partition in partitions]
    coordinates = [(0, p) for p in source]+[(1, p) for p in e_source]
    change = matrix(ZZ, len(coordinates), len(partitions),
                    lambda i, j: values[j][coordinates[i][0]].coefficient(coordinates[i][1]))
    assert change.is_square() and abs(change.det()) == 1
    inputs = matrix(ZZ, len(coordinates), len(source), lambda i, j: int(i == j))
    coefficients = integral_inverse_product(change, inputs)
    actual = diagonal_matrix(ZZ, [2**sum(part > K for part in partition)
                                  for partition in partitions])*coefficients
    return [tuple(p) for p in partitions], actual, [tuple(p) for p in source]


def three_head_typed_matrices(K, degree, models=None):
    """Return square/Euler/full matrices, including partition-indexed ingredients.

    Columns are s_gamma(head^2) and p*s_delta(head^2), with gamma in K^3,
    delta in K^2. Rows for a positive-type shape lambda are the two typed
    Schubert rows [a_lambda,+b_(lambda-K)] and [a_lambda,-b_(lambda-K)].
    """
    square_model, euler_model = models or (ThreeRootChernModel(K), TwoRootChernModel(K))
    shapes, square, plus_basis = three_head_square(K, degree, square_model)
    smaller, euler, minus_basis = two_head_euler(K, degree-K, euler_model)
    indices = {partition: i for i, partition in enumerate(smaller)}
    plus_rows, minus_rows, full_rows, typed_rows = [], [], [], []
    for i, partition in enumerate(shapes):
        a, b = list(square.row(i)), [ZZ(0)]*euler.ncols()
        if K in partition:
            mu = list(partition)
            mu.remove(K)
            if tuple(mu) in indices:
                b = list(euler.row(indices[tuple(mu)]))
            for kind, sign in ((1, 1), (2, -1)):
                signed = [sign*entry for entry in b]
                plus_rows.append(a)
                minus_rows.append(signed)
                full_rows.append(a+signed)
                typed_rows.append((partition, kind))
        else:
            plus_rows.append(a)
            minus_rows.append(b)
            full_rows.append(a+b)
            typed_rows.append((partition, 0))
    return {"plus": matrix(ZZ, plus_rows), "minus": matrix(ZZ, minus_rows),
            "full": matrix(ZZ, full_rows), "square_partitions": plus_basis,
            "Euler_partitions": minus_basis, "typed_rows": typed_rows,
            "Y_partitions": shapes, "Y_square_coefficients": square,
            "Euler_B_partitions": smaller, "Euler_B_coefficients": euler}


def matrix_exponents(coefficients):
    assert coefficients.base_ring() is ZZ
    result = []
    for factor in smith_factors(coefficients):
        exponent = int(factor.valuation(2))
        assert factor == 2**exponent
        result.append(exponent)
    return result


def analyze_matrices(a, k, degree, data):
    plus, minus, full = data["plus"], data["minus"], data["full"]
    assert all(value.base_ring() is ZZ for value in (plus, minus, full))
    exponents = [matrix_exponents(value) for value in (plus, minus, full)]
    assert len(exponents[0]) == plus.ncols() and len(exponents[1]) == minus.ncols()
    assert len(exponents[2]) == plus.ncols()+minus.ncols()
    assert exponents[0] == expanded(type_b_counts(a, k, degree))
    assert exponents[1] == expanded(type_b_counts(a-1, k, degree-k))
    # [S:L] = [S:S_plus+S_minus] * [S_plus:L_plus]*[S_minus:L_minus].
    full_gluing = sum(exponents[2])-sum(exponents[0])-sum(exponents[1])
    assert 0 <= full_gluing <= min(plus.ncols(), minus.ncols())
    mod2 = [value.change_ring(GF(2)).rank() for value in (plus, minus, full)]
    original_overlap = mod2[0]+mod2[1]-mod2[2]
    gaussian_overlap = unit_intersection(a, k, degree)
    assert original_overlap == gaussian_overlap
    assert full_gluing >= original_overlap
    expected = expanded(full_candidate(a, k, degree))
    row = {"a": a, "k": k, "n": a+k, "polynomial_degree": degree,
           "coordinate_rows": full.nrows(),
           "square_exponents": exponents[0], "Euler_exponents": exponents[1],
           "full_exponents": exponents[2], "candidate_exponents": expected,
           "sector_indices_2_adic": [sum(exponents[0]), sum(exponents[1])],
           "full_index_2_adic": sum(exponents[2]),
           "FULL_saturated_gluing_dimension": full_gluing,
           "original_image_mod2_overlap": original_overlap,
           "Gaussian_unit_overlap": gaussian_overlap,
           "FULL_gluing_exhausted_by_unit_overlap": full_gluing == original_overlap,
           "full_candidate_matches": exponents[2] == expected,
           "base_ring": "ZZ"}
    if full_gluing == original_overlap:
        row["gluing_flag_pairs"] = {"0,0": full_gluing}
    else:
        row["gluing_flag_pairs"] = "Positive-level flags require a separate calculation."
    return row


def explicit_three_head_unit_class(K, degree, data):
    """Verify the known head representatives gamma=(h,h,j), delta=(h,j)."""
    h = K//2
    j = degree//2-K
    if not 0 <= j <= h:
        return None
    gamma = tuple(value for value in (h, h, j) if value)
    delta = tuple(value for value in (h, j) if value)
    plus_i = data["square_partitions"].index(gamma)
    minus_i = data["Euler_partitions"].index(delta)
    u, v = data["plus"].column(plus_i), data["minus"].column(minus_i)
    assert gcd(list(u)) == gcd(list(v)) == 1
    assert all(entry % 2 == 0 for entry in u+v)
    half = (u+v)/2
    assert all(entry.denominator() == 1 for entry in half)
    return {"square_head_partition": list(map(int, gamma)),
            "Euler_head_partition": list(map(int, delta)),
            "both_sector_classes_primitive": True,
            "half_sum_integral": True,
            "meaning": "A nonzero unit-gluing class (f_gamma+p*f_delta)/2."}


def cartesian_matrices(a, k, degree, cache=None, force_head=False):
    """Independent actual D BGG matrices, with the smaller polynomial alphabet."""
    if cache is None:
        cache = {}
    if (a, k) not in cache:
        borel = classical_borel("D", a+k, k)
        ring = PolynomialRing(QQ, a+k, names=["x%d" % i for i in range(a+k)])
        cache[a, k] = borel, CartesianDFrontier(borel), ring, ring.gens()
    _borel, frontier, ring, x = cache[a, k]
    columns, bases = [], []
    for sector in (0, 1):
        remaining = degree-sector*k
        if remaining < 0 or remaining % 2:
            columns.append([])
            bases.append([])
            continue
        if k < a and not force_head:
            width, height, alphabet = a-sector, k, x[a:]
        else:
            width, height, alphabet = k, a-sector, x[:a]
        formal_ring, partitions, formal = schur_polynomials(width, height, remaining//2)
        substitution = formal_ring.hom([
            elementary_polynomial(ring, [value*value for value in alphabet], i)
            for i in range(1, height+1)], ring)
        # In the Euler head alphabet, omitted top elementary terms are killed
        # by p*prod(head roots)=0, so only a-1 formal elementary inputs remain.
        factor = prod(x[a:]) if sector else ring.one()
        columns.append([frontier.expand(factor*substitution(f)) for f in formal])
        bases.append(partitions)
    words = tuple({word for sector in columns for column in sector for word in column})
    coefficients = [matrix(ZZ, len(words), len(sector),
                           lambda i, j: sector[j].get(words[i], 0)) for sector in columns]
    return {"plus": coefficients[0], "minus": coefficients[1],
            "full": coefficients[0].augment(coefficients[1]),
            "square_partitions": bases[0], "Euler_partitions": bases[1]}


def timeout_handler(_signal, _frame):
    raise TimeoutError("component wall-clock limit reached")


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, here.with_name("type_d_cartesian_frontier.py"),
              here.with_name("verify_type_d_three_row_smith.py"),
              here.with_name("verify_type_d_two_row_smith.py"),
              here.with_name("verify_type_d_square_comparison.py"),
              here.with_name("verify_type_d_candidate.py"),
              borel/"src"/"sage_borel"/"core.py", borel/"src"/"sage_borel"/"operations.py"]
    hashes = {str(path.relative_to(borel)): sha256(path.read_bytes()).hexdigest() for path in inputs}
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION,
              "scope": "Finite full-D evidence. The general theorem is proved separately in TYPE_D_EVEN_GLUING_THEOREM.md; these computations do not replace that proof.",
              "normalization": "q_i=e_i(x_tail^2), p=prod(x_tail), polynomial degree d means cohomological 2d",
              "source_sha256": hashes, "three_head_components": [], "higher_Cartesian_components": [],
              "excluded_components": [],
              "resource_limits": {"Cartesian_component_seconds": 90, "total_seconds": 480}}
    output = here.with_name("type_d_even_gluing_verification.json")

    def save():
        report["seconds"] = round(perf_counter()-started, 6)
        payload = json.dumps(report, indent=2)+"\n"
        try:
            output.write_text(payload)
        except OSError:
            fallback = Path("/tmp/type_d_even_gluing_verification.json")
            fallback.write_text(payload)
            print("Report checkpoint fallback: %s" % fallback, flush=True)

    cartesian_cache = {}
    for K in (2, 4, 6, 8, 10, 12):
        models = ThreeRootChernModel(K), TwoRootChernModel(K)
        for degree in range(0, 6*K+1, 2):
            stamp = perf_counter()
            data = three_head_typed_matrices(K, degree, models)
            row = analyze_matrices(3, K, degree, data)
            row["method"] = "Chern/typed-Schubert assembly"
            row["explicit_unit_gluing_generator"] = explicit_three_head_unit_class(K, degree, data)
            cross = K == 2 or (K == 4 and degree in (8, 10, 12, 14, 16, 18, 20, 22, 24))
            if cross:
                independent = cartesian_matrices(3, K, degree, cartesian_cache, force_head=True)
                assert independent["full"].row_module() == data["full"].row_module(), (K, degree)
            row["independent_FULL_Cartesian_row_module_crosscheck"] = cross
            row["seconds"] = round(perf_counter()-stamp, 6)
            report["three_head_components"].append(row)
            save()
            if not row["FULL_gluing_exhausted_by_unit_overlap"]:
                print("COUNTEREXAMPLE: %s" % row, flush=True)
                return
        print("a=3,k=%d: all %d degrees, FULL W and original overlap computed separately" %
              (K, 3*K+1), flush=True)
    signal.signal(signal.SIGALRM, timeout_handler)
    higher = ((5, 4, 14), (5, 4, 16), (6, 2, 22),
              (4, 4, 18), (4, 4, 24), (5, 4, 22))
    for a, k, degree in higher:
        remaining = int(480-(perf_counter()-started))
        stamp = perf_counter()
        if remaining <= 0:
            report["excluded_components"].append({"a": a, "k": k, "polynomial_degree": degree,
                                                   "reason": "overall wall-clock limit"})
            continue
        signal.alarm(min(90, remaining))
        try:
            data = cartesian_matrices(a, k, degree, cartesian_cache)
            row = analyze_matrices(a, k, degree, data)
            row["method"] = "independent Cartesian D divided differences"
            row["supersedes_historical_fundamental_weight_timeout"] = (a, k, degree) in ((5, 4, 14), (5, 4, 16))
            row["seconds"] = round(perf_counter()-stamp, 6)
            report["higher_Cartesian_components"].append(row)
            print("Cartesian a=%d,k=%d,d=%d: W=%d, original overlap=%d" %
                  (a, k, degree, row["FULL_saturated_gluing_dimension"], row["original_image_mod2_overlap"]), flush=True)
        except TimeoutError as error:
            report["excluded_components"].append({"a": a, "k": k, "polynomial_degree": degree,
                                                   "reason": str(error), "seconds": round(perf_counter()-stamp, 6)})
        finally:
            signal.alarm(0)
        save()
    assert all(sha256(path.read_bytes()).hexdigest() == hashes[str(path.relative_to(borel))] for path in inputs)
    proof_files = ("TYPE_D_CANDIDATE.md", "TYPE_D_FROBENIUS.md", "TYPE_D_SMITH.md", "TYPE_D_SQUARE_THEOREM.md",
                   "TYPE_D_EVEN_GLUING_THEOREM.md", "TYPE_D_RELATIVE_EULER_EXTENSION.md",
                   "TYPE_D_HEAD_EULER_B_COMPARISON.md")
    report["proof_snapshot_sha256_at_finalization"] = {
        name: sha256((borel/name).read_bytes()).hexdigest() for name in proof_files}
    all_rows = report["three_head_components"]+report["higher_Cartesian_components"]
    report["summary"] = {
        "three_head_full_components": len(report["three_head_components"]),
        "independent_FULL_Cartesian_crosschecks": sum(row["independent_FULL_Cartesian_row_module_crosscheck"] for row in report["three_head_components"]),
        "higher_Cartesian_full_components": len(report["higher_Cartesian_components"]),
        "FULL_saturated_gluing_dimensions_computed": len(all_rows),
        "explicit_primitive_unit_gluing_generators": sum(row["explicit_unit_gluing_generator"] is not None for row in report["three_head_components"]),
        "old_timeouts_superseded": sum(row["supersedes_historical_fundamental_weight_timeout"] for row in report["higher_Cartesian_components"]),
        "excluded_components": len(report["excluded_components"]),
        "positive_level_gluing_counterexamples": sum(not row["FULL_gluing_exhausted_by_unit_overlap"] for row in all_rows),
        "full_candidate_mismatches": sum(not row["full_candidate_matches"] for row in all_rows)}
    save()
    print(report["summary"], flush=True)


if __name__ == "__main__":
    main()
