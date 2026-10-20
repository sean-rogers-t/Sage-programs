"""Exact typed-Chern checks of the two-head even-tail D gluing theorem.

Run with Sage Python from borel/experiments:
    sage -python verify_type_d_two_head_gluing.py

The previously audited TwoRootChernModel supplies ordinary integral
Schur/LR arithmetic and the symplectic two-root presentation. No Borel,
Weyl, divided-difference, or candidate-verifier code is imported.
This checks the new extraction lattice and computes actual full typed
Schubert Smith lists, separately comparing the claimed formula.
"""

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import Partition, QQ, ZZ, diagonal_matrix, matrix, vector
from sage.env import SAGE_VERSION

from verify_type_d_two_row_smith import TwoRootChernModel, expected_factors


def partitions(size, width):
    if size < 0:
        return []
    if not size:
        return [Partition([])]
    return [Partition([size-second, second])
            for second in range(size//2+1) if size-second <= width]


def tau(partition, half_tail):
    first = partition[0] if partition else 0
    second = partition[1] if len(partition) > 1 else 0
    return 0 if first <= half_tail-1 else 2 if second >= half_tail+1 else 1


def nonzero_smith(value):
    return tuple(abs(ZZ(factor)) for factor in value.elementary_divisors() if factor)


def parameter_checks(k):
    assert k >= 2 and k % 2 == 0
    h = k//2
    model = TwoRootChernModel(k)
    checks = []
    for half_degree in range(2*k+1):
        degree = 2*half_degree
        source_parts = partitions(half_degree, k)
        e_parts = partitions(half_degree-1, k)
        quotient_parts = partitions(half_degree-1, k-1)
        schubert_parts = model.isotropic_partitions(degree)
        values = [model.theta(partition) for partition in schubert_parts]
        coordinates = [(0, part) for part in source_parts]
        coordinates += [(1, part) for part in e_parts]
        change = matrix(ZZ, len(coordinates), len(schubert_parts),
                        lambda i, j: values[j][coordinates[i][0]].coefficient(coordinates[i][1]))
        assert change.is_square() and abs(change.det()) == 1, (k, degree)
        source = matrix(ZZ, len(coordinates), len(source_parts),
                        lambda i, j: int(i == j))
        coefficients = change.inverse()*source
        assert all(value.denominator() == 1 for value in coefficients.list())
        coefficients = matrix(ZZ, coefficients)
        retained = [i for i, part in enumerate(schubert_parts)
                    if not part or part[0] < model.M]
        weights = [sum(part > k for part in shape) for shape in schubert_parts]

        # Four times the exact rational extraction coordinates of the
        # integral Y basis; comparison is equality of integer row modules.
        extraction = matrix(QQ, len(retained), len(quotient_parts),
                            lambda i, j: values[retained[i]][1].coefficient(quotient_parts[j])
                            / ZZ(2)**weights[retained[i]])
        scaled = matrix(ZZ, 4*extraction)
        expected_diagonal = [ZZ(4)//ZZ(2)**tau(part, h) for part in quotient_parts]
        expected_lattice = diagonal_matrix(ZZ, expected_diagonal)
        assert scaled.row_module() == expected_lattice.row_module(), (k, degree)

        euler_index = None
        euler_level = None
        typed_shape = None
        glue_test = None
        if degree >= k and (degree-k) % 2 == 0:
            euler_index = (degree-k)//2
            if not 0 <= euler_index <= k:
                euler_index = None
        if euler_index is not None:
            r = 2*euler_index
            typed_shape = tuple(sorted([k, r], reverse=True))
            if not typed_shape[-1]:
                typed_shape = typed_shape[:-1]
            typed_row = retained.index(schubert_parts.index(typed_shape))
            extraction_row = scaled.row(typed_row)
            glue_test = extraction_row in (2*expected_lattice).row_module()
            assert glue_test == (euler_index == h), (k, degree, typed_shape, extraction_row)
            euler_level = 0 if r <= k else 1

        # Actual full D matrix in the individual typed Schubert basis.
        # Equal square coefficients go on both members of a typed pair;
        # the Euler column has opposite signs on its unique support pair.
        rows = []
        euler_support_count = 0
        for i in retained:
            shape = schubert_parts[i]
            square = [ZZ(2)**weights[i]*value for value in coefficients.row(i)]
            typed = k in shape
            if typed:
                for sign in [1, -1]:
                    row = list(square)
                    if euler_index is not None:
                        entry = sign*ZZ(2)**euler_level if shape == typed_shape else ZZ(0)
                        if entry:
                            euler_support_count += 1
                        row.append(entry)
                    rows.append(row)
            else:
                row = list(square)
                if euler_index is not None:
                    row.append(ZZ(0))
                rows.append(row)
        column_count = len(source_parts)+int(euler_index is not None)
        full_matrix = matrix(ZZ, len(rows), column_count, [value for row in rows for value in row])
        assert full_matrix.rank() == column_count, (k, degree)
        if euler_index is not None:
            assert euler_support_count == 2

        actual = nonzero_smith(full_matrix)
        expected = list(expected_factors(source_parts, k))
        if euler_index is not None:
            expected.append(ZZ(2)**euler_level)
        if degree == 2*k:
            expected.remove(ZZ(1))
            expected.remove(ZZ(1))
            expected.extend([ZZ(1), ZZ(2)])
        expected = tuple(sorted(expected))
        assert actual == expected, (k, degree, actual, expected)
        assert all(value.divides(ZZ(4)) for value in actual)
        checks.append({
            "head_dimension": 2, "tail_dimension": k, "polynomial_degree": degree,
            "source_square_partitions": [list(map(int, part)) for part in source_parts],
            "Chern_change_unimodular": True,
            "extraction_quotient_partitions": [list(map(int, part)) for part in quotient_parts],
            "extraction_denominator_exponents": [tau(part, h) for part in quotient_parts],
            "extraction_scaled_integer_row_modules_equal": True,
            "euler_index": euler_index,
            "euler_sector_exponent": euler_level,
            "typed_partition": list(typed_shape) if typed_shape is not None else None,
            "typed_extraction_in_2T": glue_test,
            "gluing_dimension": int(glue_test) if glue_test is not None else 0,
            "full_smith_factors": list(map(str, actual)),
            "full_closed_formula_matches": True,
            "full_factors_divide_four": True,
        })
    print("a=2 k=%d: %d components, extraction lattice and full typed Smith lists pass"
          % (k, len(checks)), flush=True)
    return checks


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, borel/"TYPE_D_TWO_HEAD_GLUING.md",
              here.with_name("verify_type_d_two_row_smith.py"),
              borel/"TYPE_D_TWO_ROW_SMITH.md"]
    hashes = {str(path.relative_to(borel)): sha256(path.read_bytes()).hexdigest() for path in inputs}
    checks = []
    for k in range(2, 25, 2):
        checks.extend(parameter_checks(k))
    summary = {
        "even_tail_dimensions_checked": 12,
        "all_even_degree_components_checked": len(checks),
        "exact_extraction_lattice_equalities": len(checks),
        "Euler_line_gluing_tests": sum(row["euler_index"] is not None for row in checks),
        "nonzero_gluing_components": sum(row["gluing_dimension"] for row in checks),
        "full_typed_Smith_lists_checked": len(checks),
        "counterexamples": 0,
    }
    report = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sage_version": SAGE_VERSION,
        "status": "passed",
        "scope": "Finite exact typed-Chern corroboration of the all-even-tail two-head proof; no Borel/Weyl engine.",
        "normalization": "normalized full tail inputs; polynomial degree is half cohomological degree",
        "source_sha256": hashes,
        "summary": summary,
        "degree_checks": checks,
        "seconds": round(perf_counter()-started, 6),
    }
    assert all(sha256(path.read_bytes()).hexdigest() == hashes[str(path.relative_to(borel))]
               for path in inputs)
    output = here.with_name("type_d_two_head_gluing_verification.json")
    output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": "passed", **summary, "seconds": report["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
