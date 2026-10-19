"""Exact Chern/LR checks for the type-D two-row square comparison.

Run from borel/: sage -python experiments/verify_type_d_two_row_smith.py

This uses integral Schur Littlewood--Richardson multiplication and the
two-variable Chern presentation. It imports no Borel, Weyl-group, or
divided-difference code. The finite checks corroborate the written proof;
the all-rank theorem is not inferred from a finite range of parameters.
"""

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import Partition, SymmetricFunctions, ZZ, diagonal_matrix, matrix, vector
from sage.env import SAGE_VERSION


SCHUR = SymmetricFunctions(ZZ).schur()


class TwoRootChernModel:
    """C = R[e,h]/(e^2-t, h^2-s-2e), R = H*(Gr(2,K+2);Z).

    An even class is represented by a pair (A,B) for A+e*B, where A,B
    are truncated ordinary Schur expansions in the two squared roots.
    An odd Chern class is stored as h times such a pair.
    """

    def __init__(self, tail_dimension):
        self.K = tail_dimension
        self.M = 2*tail_dimension+2

    def truncate(self, polynomial):
        return SCHUR._from_dict({
            partition: coefficient for partition, coefficient in polynomial
            if len(partition) <= 2 and (not partition or partition[0] <= self.K)
        })

    def product(self, left, right):
        return self.truncate(left*right)

    def complete(self, degree):
        return SCHUR[degree] if 0 <= degree <= self.K else SCHUR.zero()

    def chern(self, index):
        if index < 0 or index > self.M:
            return 0, (SCHUR.zero(), SCHUR.zero())
        if index % 2:
            return 1, (self.complete(index//2), SCHUR.zero())
        return 0, (self.complete(index//2), self.complete(index//2-1))

    def even_product(self, left, right):
        a, b = left
        c, d = right
        return (self.product(a, c)+self.product(SCHUR[1, 1], self.product(b, d)),
                self.product(a, d)+self.product(b, c))

    def chern_product(self, first, second):
        if first < 0 or first > self.M or second < 0 or second > self.M:
            return SCHUR.zero(), SCHUR.zero()
        first_parity, left = self.chern(first)
        second_parity, right = self.chern(second)
        assert (first_parity+second_parity) % 2 == 0
        result = self.even_product(left, right)
        if first_parity and second_parity:
            result = self.even_product(result, (SCHUR[1], 2*SCHUR.one()))
        return result

    @staticmethod
    def add(left, right, scalar=1):
        return left[0]+scalar*right[0], left[1]+scalar*right[1]

    def theta(self, partition):
        """Published type-C two-row theta Giambelli class in C."""
        if not partition:
            return SCHUR.one(), SCHUR.zero()
        if len(partition) == 1:
            return self.chern(partition[0])[1]
        first, second = partition
        result = self.chern_product(first, second)
        if first+second < 2*self.K+2:
            return self.add(result, self.chern_product(first+1, second-1), -1)
        for shift in range(1, second+1):
            result = self.add(result, self.chern_product(first+shift, second-shift),
                              2*(-1)**shift)
        return result

    def rectangle_partitions(self, size):
        if size < 0:
            return []
        return [Partition([size-second, second]) for second in range(size//2+1)
                if size-second <= self.K]

    def isotropic_partitions(self, size):
        if not size:
            return [()]
        return [(size-second, second) if second else (size,)
                for second in range(size//2+1)
                if size-second <= self.M
                and not (size-second == second and second > self.K)]


def expected_factors(partitions, K):
    upper, lower = (K+1)//2, K//2
    result = []
    for partition in partitions:
        first = partition[0] if partition else 0
        second = partition[1] if len(partition) > 1 else 0
        exponent = (0 if first <= upper and second <= lower
                    else 2 if second > upper else 1)
        result.append(ZZ(2)**exponent)
    return tuple(sorted(result))


def parameter_checks(K):
    model = TwoRootChernModel(K)
    checks = []
    for half_degree in range(2*K+1):
        degree = 2*half_degree
        source_partitions = model.rectangle_partitions(half_degree)
        e_partitions = model.rectangle_partitions(half_degree-1)
        schubert_partitions = model.isotropic_partitions(degree)
        values = [model.theta(partition) for partition in schubert_partitions]
        coordinates = [(0, partition) for partition in source_partitions]
        coordinates += [(1, partition) for partition in e_partitions]
        change = matrix(ZZ, len(coordinates), len(schubert_partitions),
                        lambda i, j: values[j][coordinates[i][0]].coefficient(coordinates[i][1]))
        assert change.is_square() and abs(change.det()) == 1, (K, degree)
        source = matrix(ZZ, len(coordinates), len(source_partitions),
                        lambda i, j: int(i == j))
        rational_coefficients = change.inverse()*source
        assert all(value.denominator() == 1 for value in rational_coefficients.list())
        coefficients = matrix(ZZ, rational_coefficients)
        weights = [sum(part > K for part in partition) for partition in schubert_partitions]
        b_matrix = diagonal_matrix(ZZ, [2**weight for weight in weights])*coefficients
        retained = [index for index, partition in enumerate(schubert_partitions)
                    if not partition or partition[0] < model.M]
        missing = [index for index in range(len(schubert_partitions)) if index not in retained]
        d_matrix = b_matrix.matrix_from_rows(retained)

        relations, parity_checks = [], 0
        for missing_index in missing:
            partition = schubert_partitions[missing_index]
            second = partition[1] if len(partition) > 1 else 0
            assert second % 2 == 0
            s = second//2
            target = Partition([K, s])
            assert values[missing_index] == (SCHUR.zero(), SCHUR[target])
            extraction = [value[1].coefficient(target) for value in values]
            assert all(extraction[index] == int(index == missing_index) for index in missing)
            unweighted_relation = coefficients.row(missing_index)
            for index in retained:
                unweighted_relation += extraction[index]*coefficients.row(index)
            assert unweighted_relation.is_zero(), (K, degree, partition)

            weighted_relation = vector(ZZ, len(source_partitions))
            relation_coefficients = []
            for index in retained:
                if weights[missing_index] == 1 and weights[index] == 2:
                    assert extraction[index] % 2 == 0, (K, degree, partition, schubert_partitions[index])
                    parity_checks += 1
                scalar = -ZZ(2)**(weights[missing_index]-weights[index])*extraction[index]
                assert scalar.denominator() == 1
                scalar = ZZ(scalar)
                weighted_relation += scalar*b_matrix.row(index)
                if scalar:
                    relation_coefficients.append({"retained_partition": list(schubert_partitions[index]),
                                                  "coefficient": str(scalar)})
            assert weighted_relation == b_matrix.row(missing_index)
            relations.append({"missing_partition": list(partition),
                              "weight_exponent": weights[missing_index],
                              "weighted_relation_coefficients": relation_coefficients,
                              "exact_integral_row_relation_verified": True})

        assert b_matrix.row_module() == d_matrix.row_module(), (K, degree)
        b_factors = tuple(abs(ZZ(value)) for value in b_matrix.elementary_divisors() if value)
        d_factors = tuple(abs(ZZ(value)) for value in d_matrix.elementary_divisors() if value)
        expected = expected_factors(source_partitions, K)
        assert b_factors == d_factors == expected, (K, degree, b_factors, d_factors, expected)
        assert all(factor.divides(ZZ(2)**min(2, K)) for factor in d_factors)
        checks.append({"a": 2, "k": K, "polynomial_degree": degree,
                       "source_partitions": [list(map(int, partition)) for partition in source_partitions],
                       "rank": len(source_partitions),
                       "C_schubert_to_R_plus_eR_unimodular": True,
                       "B_smith_factors": list(map(str, b_factors)),
                       "D_square_smith_factors": list(map(str, d_factors)),
                       "closed_partition_formula_matches": True,
                       "weighted_integer_row_modules_equal": True,
                       "missing_row_relations": relations,
                       "weight_one_vs_weight_two_parity_checks": parity_checks})
    print("a=2 k=%d: %d degree components, every weighted row relation verified" % (K, len(checks)), flush=True)
    return checks


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, borel/"TYPE_D_TWO_ROW_SMITH.md", borel/"TYPE_B_TWO_ROW_SMITH.md"]
    hashes = {str(path.relative_to(borel)): sha256(path.read_bytes()).hexdigest() for path in inputs}
    checks = []
    for K in range(1, 13):
        checks.extend(parameter_checks(K))
    # The proof uses the already-proved block duality for the other order;
    # no second finite calculation is counted for exchanging the blocks.
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "status": "passed",
              "scope": "Finite integral Chern/LR corroboration of the all-rank two-row proof; no Borel or Weyl engine.",
              "normalization": "normalized square-only inputs, polynomial degree is half cohomological degree",
              "source_sha256": hashes, "degree_checks": checks,
              "summary": {"tail_dimensions_checked": 12,
                          "degree_components_checked": len(checks),
                          "weighted_integer_row_module_equalities": len(checks),
                          "exact_missing_row_relations": sum(len(row["missing_row_relations"]) for row in checks),
                          "weight_one_vs_weight_two_parity_checks": sum(row["weight_one_vs_weight_two_parity_checks"] for row in checks),
                          "closed_partition_formula_matches": len(checks),
                          "counterexamples": 0},
              "seconds": round(perf_counter()-started, 6)}
    assert all(sha256(path.read_bytes()).hexdigest() == hashes[str(path.relative_to(borel))] for path in inputs)
    output = here.with_name("type_d_two_row_smith_verification.json")
    output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], **report["summary"], "seconds": report["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
