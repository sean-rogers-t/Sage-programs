"""Integral Chern/raising-operator checks for three head roots in type D.

Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_three_row_smith.py

This finite experiment computes actual integer row modules, not just rational
row spaces. It uses the binary Chern model and both published Giambelli
thresholds, with independent Cartesian divided-difference cross-checks in
specified small cases. No conclusion for untested ranks follows from it.
"""

from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from time import perf_counter

from sage.all import (Partition, Partitions, PolynomialRing, QQ,
                      SymmetricFunctions, ZZ, diagonal_matrix, matrix, prod,
                      vector)
from sage.env import SAGE_VERSION

from sage_borel import classical_borel
from type_d_cartesian_frontier import CartesianDFrontier


SCHUR = SymmetricFunctions(ZZ).schur()
POLYNOMIALS = PolynomialRing(ZZ, 3, names=["e1", "e2", "e3"])
E = POLYNOMIALS.gens()
Q_RING = PolynomialRing(ZZ, "q")
q = Q_RING.gen()


@lru_cache(None)
def gaussian(n, k):
    if k < 0 or k > n:
        return Q_RING.zero()
    if k == 0 or k == n:
        return Q_RING.one()
    return gaussian(n-1, k)+q**(n-k)*gaussian(n-1, k-1)


def closed_multiplicities(K):
    """The G,F,H,J formulas of TYPE_B_THREE_ROW_SMITH.md."""
    u, v, t = (K+2)//2, (K+1)//2, (K+2)//2
    G = gaussian(K+3, 3)
    F = sum((q**(r+s+w) for w in range(v+1)
             for s in range(w, u+1) for r in range(s, K+1)), Q_RING.zero())
    H = q**(3*(t+1))*gaussian(K-t+2, 3)
    b = K//2
    if K % 2:
        J = gaussian(b+3, 3)+q**(b+1)*gaussian(b+2, 2)
    else:
        J = (gaussian(b+3, 3)
             + q**(b+1)*(gaussian(b+2, 2)-q**(2*b))
             - q**(2*b+1)*sum((q**i for i in range(b)), Q_RING.zero()))
    n = K+3
    r, c = n//2, (n+1)//2
    positive_J = gaussian(c, 3)+q**(c-2)*gaussian(r, 1)*gaussian(c-1, 1)
    assert J == positive_J
    answers = (J, F-J, G-H-F, H)
    positive_answers = (
        positive_J,
        q**c*gaussian(r, 1)*(gaussian(c-1, 2)+q**(c-2)*gaussian(r-1, 1)),
        q**(2*c)*gaussian(r, 2)*gaussian(c-2, 1),
        q**(3*c)*gaussian(r, 3),
    )
    assert answers == positive_answers
    assert all(coefficient >= 0 for answer in answers for coefficient in answer)
    assert sum(answers, Q_RING.zero()) == G
    return answers


class ThreeRootChernModel:
    r"""C = R[e1,e2,e3] with square relations and R=Gr(3,K+3).

    The relations are e1^2=s1+2e2, e2^2=s2+2e1e3, e3^2=s3,
    where si is elementary in the three squared roots. Binary e-monomials
    form the relative integral basis. In even degrees the basis is
    1,e2,e1e3,e1e2e3. The D Chern quotient removes the coordinates
    e1e3*s_(K,mu) and e1e2e3*s_(K,mu), namely the even ideal of e3*h_K.
    """

    def __init__(self, K):
        self.K = K
        self.maximum = 2*K+3

    def truncate(self, polynomial):
        return SCHUR._from_dict({
            partition: coefficient for partition, coefficient in polynomial
            if len(partition) <= 3 and (not partition or partition[0] <= self.K)
        })

    def multiply(self, left, right):
        return self.truncate(left*right)

    @lru_cache(None)
    def reduce_binary(self, exponent):
        if max(exponent, default=0) <= 1:
            return {exponent: SCHUR.one()}
        position = next(i for i in range(3) if exponent[i] >= 2)
        r = position+1
        base = list(exponent)
        base[position] -= 2
        result = {basis: self.multiply(SCHUR[tuple([1]*r)], value)
                  for basis, value in self.reduce_binary(tuple(base)).items()}
        for shift in range(1, r+1):
            if r+shift > 3:
                continue
            following = list(base)
            following[r+shift-1] += 1
            if r-shift > 0:
                following[r-shift-1] += 1
            for basis, value in self.reduce_binary(tuple(following)).items():
                result[basis] = result.get(basis, SCHUR.zero())-2*(-1)**shift*value
        return {basis: value for basis, value in result.items() if value}

    def decompose(self, polynomial):
        result = {}
        for exponent, coefficient in polynomial.dict().items():
            for basis, value in self.reduce_binary(tuple(exponent)).items():
                result[basis] = result.get(basis, SCHUR.zero())+coefficient*value
        return {basis: value for basis, value in result.items() if value}

    @lru_cache(None)
    def chern(self, index):
        if index < 0 or index > self.maximum:
            return POLYNOMIALS.zero()
        if index == 0:
            return POLYNOMIALS.one()
        return sum((-1)**(j-1)*E[j-1]*self.chern(index-j)
                   for j in range(1, min(3, index)+1))

    @lru_cache(None)
    def theta(self, partition, kind="C"):
        """Unnormalized C or D Giambelli numerator, including all 3 pairs."""
        if not partition:
            return POLYNOMIALS.one()
        padded = list(partition)+[0]*(3-len(partition))
        pairs = ((0, 1), (0, 2), (1, 2))
        critical = {
            (i, j): (padded[i]+padded[j] > 2*self.K+j-i if kind == "C"
                     else padded[i]+padded[j] >= 2*self.K+j-i)
            for i, j in pairs
        }

        def coefficient(i, j, shift):
            if shift == 0:
                return 1
            if critical[i, j]:
                return 2*(-1)**shift
            return -1 if shift == 1 else 0

        result = POLYNOMIALS.zero()
        # A transfer 2->1 can be fed by a later 3->2 transfer. Bounds
        # therefore refer to the FINAL indices, not intermediate positivity.
        for r23 in range(padded[2]+1):
            a23 = coefficient(1, 2, r23)
            if not a23:
                continue
            for r13 in range(padded[2]-r23+1):
                a13 = coefficient(0, 2, r13)
                if not a13:
                    continue
                for r12 in range(padded[1]+r23+1):
                    a12 = coefficient(0, 1, r12)
                    if not a12:
                        continue
                    indices = (padded[0]+r12+r13,
                               padded[1]-r12+r23, padded[2]-r13-r23)
                    result += a12*a13*a23*prod(self.chern(index) for index in indices)
        return result

    def isotropic_partitions(self, degree, width):
        return [tuple(partition) for partition in
                Partitions(degree, max_part=width, max_length=3)
                if all(list(partition).count(part) <= 1
                       for part in partition if part > self.K)]

    def rectangle_partitions(self, size):
        return ([tuple(partition) for partition in
                 Partitions(size, max_part=self.K, max_length=3)]
                if size >= 0 else [])

    def coordinates(self, degree, ideal=False):
        result = []
        for basis in product(range(2), repeat=3):
            basis_degree = sum((i+1)*basis[i] for i in range(3))
            if (degree-basis_degree) % 2 or (ideal and not basis[2]):
                continue
            for partition in self.rectangle_partitions((degree-basis_degree)//2):
                if ideal and (not partition or partition[0] != self.K):
                    continue
                result.append((basis, partition))
        return result

    @staticmethod
    def coordinate_vector(value, coordinates):
        return vector(ZZ, [value.get(basis, SCHUR.zero()).coefficient(Partition(partition))
                           for basis, partition in coordinates])


def integral_inverse_product(change, source):
    rational = change.inverse()*source
    assert all(value.denominator() == 1 for value in rational.list())
    answer = matrix(ZZ, rational)
    assert answer.base_ring() is ZZ
    return answer


def smith_factors(coefficients):
    assert coefficients.base_ring() is ZZ
    return tuple(abs(ZZ(value)) for value in coefficients.elementary_divisors() if value)


def cartesian_head_matrix(K, degree, partitions, cache):
    """An independently calculated D matrix in the same head Schur basis."""
    if K not in cache:
        borel = classical_borel("D", K+3, K)
        frontier = CartesianDFrontier(borel)
        ring = PolynomialRing(QQ, K+3, names=["x%d" % i for i in range(K+3)])
        x = ring.gens()[:3]
        u = tuple(value*value for value in x)
        elementary = (u[0]+u[1]+u[2],
                      u[0]*u[1]+u[0]*u[2]+u[1]*u[2], prod(u))
        complete = [ring.one()]
        for index in range(1, K+3):
            complete.append(sum((-1)**(j-1)*elementary[j-1]*complete[index-j]
                                for j in range(1, min(3, index)+1)))
        cache[K] = borel, frontier, ring, complete
    borel, frontier, ring, complete = cache[K]
    columns = []
    for partition in partitions:
        length = len(partition)
        polynomial = matrix(ring, length, length,
                            lambda i, j: complete[partition[i]-i+j]
                            if partition[i]-i+j >= 0 else ring.zero()).det()
        columns.append(frontier.expand(polynomial))
    words = tuple({word for column in columns for word in column})
    return matrix(ZZ, len(words), len(columns),
                  lambda i, j: columns[j].get(words[i], 0))


def relative_quadric_checks():
    r"""A nontrivial-base test of the compensated parity projection.

    For B=CP^2 and its rank-two quotient W, c1(W)=h,c2(W)=h^2.
    The two rulings of Q(W+W*) give the projective-bundle algebra
    Z[h,alpha,beta]/(h^3,alpha^2-h alpha,beta^2-h beta+h^2).
    The formulas below can also be checked at its four relative torus
    fixed points: z has weights +/-w1,+/-w2 and zp=w1*w2.
    """
    ring = PolynomialRing(ZZ, 3, names=["alpha", "beta", "h"])
    alpha, beta, h = ring.gens()

    @lru_cache(None)
    def reduce_term(a, b, c):
        if c >= 3:
            return ring.zero()
        if a >= 2:
            return reduce_term(a-1, b, c+1)
        if b >= 2:
            return reduce_term(a, b-1, c+1)-reduce_term(a, b-2, c+2)
        return alpha**a*beta**b*h**c

    def normal(polynomial):
        return sum((coefficient*reduce_term(*map(int, exponent))
                    for exponent, coefficient in ring(polynomial).dict().items()), ring.zero())

    z, p = alpha+beta-h, alpha-beta
    u, euler = z*z, h*h
    compensated = ring.hom([h-alpha, h-beta, h], ring)
    assert normal(euler) != 0
    assert normal(z*p-euler) == 0
    assert normal(p*p+z*z+h*h) == 0
    assert normal(u*u+h*h*u) == 0
    assert normal(compensated(z)+z) == normal(compensated(p)+p) == 0
    denominator_identities = [
        (2*alpha, z+p+h), (2*beta, z-p+h),
        (2*alpha*beta, z*z+h*z+h*h),
    ]
    assert all(normal(left-right) == 0 for left, right in denominator_identities)
    projection_checks = 0
    for base_degree in range(3):
        for head_degree in range(7):
            for euler_power in range(2):
                f = h**base_degree*z**head_degree*p**euler_power
                summed = normal(f+compensated(f))
                assert all(coefficient % 2 == 0 for coefficient in summed.coefficients())
                averaged = ring(summed/2)
                if not euler_power and head_degree % 2 == 0:
                    expected = h**base_degree*u**(head_degree//2)
                elif euler_power and head_degree % 2:
                    expected = h**(base_degree+2)*u**((head_degree-1)//2)
                else:
                    expected = ring.zero()
                assert normal(averaged-expected) == 0
                projection_checks += 1
    # Averaging on the whole integral cohomology need not be integral.
    assert normal(alpha+compensated(alpha)) == h
    # An uncompensated deletion of odd head powers sends zp-h^2 to -h^2.
    assert normal(-euler) != 0
    degree_checks = []
    for degree in range(5):
        target = [(a, b, c) for a, b in product(range(2), repeat=2)
                  for c in range(3) if a+b+c == degree]
        source = [h**i*u**j for i in range(3) for j in range(2) if i+2*j == degree]
        coefficients = matrix(ZZ, len(target), len(source),
                              lambda i, j: normal(source[j]).monomial_coefficient(
                                  alpha**target[i][0]*beta**target[i][1]*h**target[i][2]))
        factors = smith_factors(coefficients)
        assert len(factors) == len(source)
        assert all(factor.divides(ZZ(2)) for factor in factors)
        if degree == 4:
            assert factors == (ZZ(2),)
        degree_checks.append({"polynomial_degree": degree, "source_rank": len(source),
                              "square_ring_smith_factors": list(map(str, factors))})
    return {"base": "CP^2", "tautological_quotient_rank": 2,
            "head_size": 1, "tail_size": 1,
            "nonzero_base_Euler_class": "h^2",
            "integral_polynomial_projection_checks": projection_checks,
            "twice_cohomology_in_Chern_Euler_subring": True,
            "compensated_projection_preserves_nonzero_Euler_relation": True,
            "uncompensated_head_projection_fails": True,
            "averaging_all_integral_cohomology_fails": "rho(alpha)=h/2",
            "quaternionic_square_relation": "u^2+h^2*u=0",
            "square_ring_injection_and_denominator_bound_checks": degree_checks}


def parameter_checks(K, cartesian_cache):
    model = ThreeRootChernModel(K)
    multiplicities = closed_multiplicities(K)
    checks = []
    for half_degree in range(3*K+1):
        degree = 2*half_degree
        partitions = model.isotropic_partitions(degree, model.maximum)
        coordinates = model.coordinates(degree)
        ideal_coordinates = model.coordinates(degree, ideal=True)
        source_coordinates = [(basis, p) for basis, p in coordinates if not any(basis)]
        values = [model.decompose(model.theta(partition)) for partition in partitions]
        c_change = matrix(ZZ, [model.coordinate_vector(value, coordinates)
                              for value in values]).transpose()
        assert c_change.is_square() and abs(c_change.det()) == 1
        source = matrix(ZZ, len(coordinates), len(source_coordinates),
                        lambda i, j: int(coordinates[i] == source_coordinates[j]))
        c_coefficients = integral_inverse_product(c_change, source)
        retained = [i for i, p in enumerate(partitions) if not p or p[0] < model.maximum]
        missing = [i for i in range(len(partitions)) if i not in retained]
        weights = [sum(part > K for part in partition) for partition in partitions]
        d_coordinates = [coordinate for coordinate in coordinates if coordinate not in ideal_coordinates]
        d_values = [model.decompose(model.theta(partitions[i], "D")) for i in retained]
        d_change = matrix(ZZ, [model.coordinate_vector(value, d_coordinates)
                              for value in d_values]).transpose()
        assert d_change.is_square() and abs(d_change.det()) == 1
        d_source = matrix(ZZ, len(d_coordinates), len(source_coordinates),
                          lambda i, j: int(d_coordinates[i] == source_coordinates[j]))
        d_coefficients = integral_inverse_product(d_change, d_source)
        prefixes = []
        cumulative = 0
        for level in range(4):
            c_prefix = c_coefficients.matrix_from_rows([
                i for i in range(len(partitions)) if weights[i] <= level]).row_module()
            d_prefix = d_coefficients.matrix_from_rows([
                j for j, i in enumerate(retained) if weights[i] <= level]).row_module()
            assert c_prefix.base_ring() is ZZ and d_prefix.base_ring() is ZZ
            assert c_prefix == d_prefix, (K, degree, level)
            assert all(value in (0, 1, -1)
                       for value in c_prefix.basis_matrix().elementary_divisors())
            cumulative += ZZ(multiplicities[level][half_degree])
            assert c_prefix.rank() == d_prefix.rank() == cumulative
            prefixes.append({"level": level, "base_ring": "ZZ", "rank": int(cumulative),
                             "exact_integer_modules_equal": True, "primitive": True})
        b_matrix = diagonal_matrix(ZZ, [2**weight for weight in weights])*c_coefficients
        d_matrix = diagonal_matrix(ZZ, [2**weights[i] for i in retained])*d_coefficients
        assert b_matrix.base_ring() is ZZ and d_matrix.base_ring() is ZZ
        assert b_matrix.row_module() == d_matrix.row_module()
        expected = tuple(ZZ(2)**level for level in range(4)
                         for _ in range(int(multiplicities[level][half_degree])))
        assert smith_factors(b_matrix) == smith_factors(d_matrix) == expected
        assert all(factor.divides(ZZ(8)) for factor in expected)
        independently_checked = K == 3 or (K == 4 and degree in (10, 14, 18))
        if independently_checked:
            direct = cartesian_head_matrix(K, degree,
                                           [p for _basis, p in source_coordinates], cartesian_cache)
            assert direct.base_ring() is ZZ
            assert direct.row_module() == d_matrix.row_module(), (K, degree)
        regression = None
        if K == 3 and degree == 10:
            projection = matrix(ZZ, [model.coordinate_vector(value, ideal_coordinates)
                                     for value in values]).transpose()
            boundary = projection.matrix_from_columns(missing)
            assert boundary.is_square() and abs(boundary.det()) == 1
            relation = boundary.inverse()*projection
            bi = partitions.index((9, 1))
            ri = partitions.index((5, 4, 1))
            entry = relation[missing.index(bi), ri]
            weighted_entry = ZZ(2)**(weights[bi]-weights[ri])*entry
            assert entry == -1 and weighted_entry == -QQ(1)/2
            regression = {"missing_C_partition": [9, 1], "retained_C_partition": [5, 4, 1],
                          "canonical_weighted_coefficient": str(weighted_entry),
                          "conclusion": "Deleting C boundary rows alone does not give the correct weighted D model."}
        checks.append({"a": 3, "k": K, "polynomial_degree": degree,
                       "source_partitions": [list(map(int, p)) for _b, p in source_coordinates],
                       "rank": len(source_coordinates),
                       "C_binary_Schubert_change_unimodular": True,
                       "D_quotient_raising_change_unimodular": True,
                       "prefixes": prefixes,
                       "weighted_integer_row_modules_equal": True,
                       "B_smith_factors": list(map(str, expected)),
                       "D_square_smith_factors": list(map(str, expected)),
                       "B_three_row_closed_formula_matches": True,
                       "four_positive_Gaussian_formulas_match": True,
                       "D_exponents_at_most_three": True,
                       "independent_Cartesian_D_row_module_crosscheck": independently_checked,
                       "canonical_boundary_regression": regression})
    print("K=%d: %d components; all integer prefixes and weighted row modules agree" %
          (K, len(checks)), flush=True)
    return checks


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, here.with_name("type_d_cartesian_frontier.py"),
              borel/"src"/"sage_borel"/"core.py",
              borel/"src"/"sage_borel"/"operations.py",
              borel/"TYPE_D_SQUARE_THEOREM.md"]
    hashes = {str(path.relative_to(borel)): sha256(path.read_bytes()).hexdigest()
              for path in inputs}
    checks, cartesian_cache = [], {}
    for K in range(3, 9):
        checks.extend(parameter_checks(K, cartesian_cache))
    relative_check = relative_quadric_checks()
    assert all(sha256(path.read_bytes()).hexdigest() == hashes[str(path.relative_to(borel))]
               for path in inputs)
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "status": "passed",
              "scope": "Finite exact integral Chern/raising corroboration only; no untested-rank conclusion.",
              "normalization": "q_i=e_i(x_tail^2); head Schur input basis; polynomial d means cohomological 2d",
              "coefficient_row_modules": "All compared over ZZ, with explicit base-ring assertions after inverse changes of basis.",
              "independence": "Schur LR and raising-operator model, cross-checked against Cartesian D divided differences in 13 components.",
              "source_sha256": hashes, "degree_checks": checks,
              "relative_twisted_quadric": relative_check,
              "summary": {"tail_dimensions_checked": 6, "degree_components_checked": len(checks),
                          "exact_integer_prefix_module_equalities": 4*len(checks),
                          "primitive_prefix_checks": 4*len(checks),
                          "weighted_integer_row_module_equalities": len(checks),
                          "closed_formula_matches": len(checks),
                          "four_positive_Gaussian_formula_matches": len(checks),
                          "independent_Cartesian_D_components": sum(row["independent_Cartesian_D_row_module_crosscheck"] for row in checks),
                          "canonical_boundary_half_counterexamples": sum(row["canonical_boundary_regression"] is not None for row in checks),
                          "relative_compensated_projection_checks": relative_check["integral_polynomial_projection_checks"],
                          "all_checks_passed": True},
              "seconds": round(perf_counter()-started, 6)}
    output = here.with_name("type_d_three_row_smith_verification.json")
    output.write_text(json.dumps(report, indent=2)+"\n")
    print(report["summary"], flush=True)


if __name__ == "__main__":
    main()
