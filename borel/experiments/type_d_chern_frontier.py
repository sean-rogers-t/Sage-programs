"""Small exact Chern/raising-operator models for even-tail D image checks.

This experimental integral coefficient engine corroborates the separate
proof in TYPE_D_EVEN_GLUING_THEOREM.md. It generalizes the independently
checked two- and three-root models. Its finite checks include full Cartesian D BGG
comparisons, which do not use this presentation or its raising operators.

The head-square ring is R=H*(Gr(a,a+k);Z), in its Schur rectangle basis.
Adjoining elementary head classes e_i gives the binary integral basis via
  e_r^2 = s_(1^r) - 2*sum_j (-1)^j e_(r-j)e_(r+j).
For D the Chern kernel e_a*h_k removes exactly the binary coordinates
with e_a present and Schur first part k. Complete head classes are the
quotient Chern classes. The C and D unnormalized Giambelli numerators
use, respectively, lambda_i+lambda_j > and >= 2k+j-i. Their normalized
Schubert coefficients multiply row lambda by 2^#{lambda_i>k}.

For positive-type D shapes, the primitive Y class is the sum of the two
typed Schubert classes. The smaller B class indexed by mu maps to their
difference at shape lambda=mu union (k). Thus the full typed rows are
[A_lambda,+B_mu] and [A_lambda,-B_mu]. These are integral coefficient
matrices; division by powers of two occurs only when explicitly forming
the unweighted filtration in the separate verifier.

Sources and conventions: BKT Giambelli/eta formulas and Euler comparison,
as cited in TYPE_D_SMITH.md and TYPE_D_SQUARE_THEOREM.md. In particular,
the D equality threshold is intentional. Polynomial degree d means
cohomological degree 2d. The public image API below is for even d and k.
"""
from functools import lru_cache
from itertools import product

from sage.all import (Partition, Partitions, PolynomialRing,
                      SymmetricFunctions, ZZ, diagonal_matrix, matrix, prod)

from verify_type_d_three_row_smith import integral_inverse_product


SCHUR = SymmetricFunctions(ZZ).schur()


class BinaryChernModel:
    """The finite R-relative binary Chern presentation for a head of size a."""

    def __init__(self, a, k):
        if a < 1 or k < 2:
            raise ValueError("Require a>=1 and k>=2")
        self.a, self.k = int(a), int(k)
        self.polynomials = PolynomialRing(ZZ, a, names=["e%d" % i for i in range(1, a+1)])
        self.elementary = self.polynomials.gens()
        self.maximum = 2*k+a

    def truncate(self, value):
        return SCHUR._from_dict({
            partition: coefficient for partition, coefficient in value
            if len(partition) <= self.a and (not partition or partition[0] <= self.k)})

    def multiply(self, left, right):
        return self.truncate(left*right)

    @lru_cache(None)
    def reduce_binary(self, exponent):
        if max(exponent, default=0) <= 1:
            return {exponent: SCHUR.one()}
        position = next(i for i in range(self.a) if exponent[i] >= 2)
        r = position+1
        base = list(exponent)
        base[position] -= 2
        answer = {basis: self.multiply(SCHUR[tuple([1]*r)], value)
                  for basis, value in self.reduce_binary(tuple(base)).items()}
        for shift in range(1, r+1):
            if r+shift > self.a:
                continue
            following = list(base)
            following[r+shift-1] += 1
            if r-shift:
                following[r-shift-1] += 1
            for basis, value in self.reduce_binary(tuple(following)).items():
                answer[basis] = answer.get(basis, SCHUR.zero())-2*(-1)**shift*value
        return {basis: value for basis, value in answer.items() if value}

    def decompose(self, polynomial):
        answer = {}
        for exponent, coefficient in polynomial.dict().items():
            for basis, value in self.reduce_binary(tuple(exponent)).items():
                answer[basis] = answer.get(basis, SCHUR.zero())+coefficient*value
        return {basis: value for basis, value in answer.items() if value}

    @lru_cache(None)
    def chern(self, index):
        if index < 0 or index > self.maximum:
            return self.polynomials.zero()
        if index == 0:
            return self.polynomials.one()
        return sum((-1)**(j-1)*self.elementary[j-1]*self.chern(index-j)
                   for j in range(1, min(self.a, index)+1))

    @lru_cache(None)
    def chern_product(self, indices):
        return prod(self.chern(index) for index in indices)

    @lru_cache(None)
    def theta(self, partition, kind):
        if kind not in ("C", "D"):
            raise ValueError("kind must be C or D")
        padded = tuple(partition)+(0,)*(self.a-len(partition))
        critical = {
            (i, j): (padded[i]+padded[j] > 2*self.k+j-i if kind == "C"
                     else padded[i]+padded[j] >= 2*self.k+j-i)
            for j in range(1, self.a) for i in range(j)}
        states = {padded: ZZ(1)}
        # Process larger source indices first. Transfers into j have already
        # happened when its outgoing transfers are bounded by its current size.
        for j in range(self.a-1, 0, -1):
            for i in range(j-1, -1, -1):
                following = {}
                for indices, coefficient in states.items():
                    bound = indices[j] if critical[i, j] else min(indices[j], 1)
                    for shift in range(bound+1):
                        factor = (1 if shift == 0 else
                                  2*(-1)**shift if critical[i, j] else -1)
                        moved = list(indices)
                        moved[j] -= shift
                        moved[i] += shift
                        moved = tuple(moved)
                        following[moved] = following.get(moved, ZZ(0))+coefficient*factor
                states = {indices: coefficient for indices, coefficient in following.items()
                          if coefficient}
        # Chern factors commute, so combine permutations before expanding.
        combined = {}
        for indices, coefficient in states.items():
            key = tuple(sorted(indices))
            combined[key] = combined.get(key, ZZ(0))+coefficient
        return sum((coefficient*self.chern_product(indices)
                    for indices, coefficient in combined.items() if coefficient),
                   self.polynomials.zero())

    def rectangle_partitions(self, size):
        return ([tuple(partition) for partition in
                 Partitions(size, max_part=self.k, max_length=self.a)] if size >= 0 else [])

    def coordinates(self, degree, kind):
        answer = []
        for basis in product(range(2), repeat=self.a):
            basis_degree = sum((i+1)*basis[i] for i in range(self.a))
            if (degree-basis_degree) % 2:
                continue
            for partition in self.rectangle_partitions((degree-basis_degree)//2):
                if kind == "D" and basis[-1] and partition and partition[0] == self.k:
                    continue
                answer.append((basis, partition))
        return answer

    def coefficient_matrix(self, degree, kind):
        """Return (shapes, actual integral rows, head Schur source partitions)."""
        if degree < 0:
            return [], matrix(ZZ, 0, 0), []
        if degree % 2:
            raise ValueError("The square input source has even polynomial degree")
        coordinates = self.coordinates(degree, kind)
        width = self.maximum-(kind == "D")
        shapes = [tuple(partition) for partition in
                  Partitions(degree, max_part=width, max_length=self.a)
                  if all(list(partition).count(part) <= 1
                         for part in partition if part > self.k)]
        values = [self.decompose(self.theta(partition, kind)) for partition in shapes]
        change = matrix(ZZ, len(coordinates), len(shapes), lambda i, j:
                        values[j].get(coordinates[i][0], SCHUR.zero()).coefficient(
                            Partition(coordinates[i][1])))
        assert change.is_square() and abs(change.det()) == 1, (self.a, self.k, degree, kind)
        source = [(basis, partition) for basis, partition in coordinates if not any(basis)]
        inputs = matrix(ZZ, len(coordinates), len(source), lambda i, j:
                        int(coordinates[i] == source[j]))
        coefficients = integral_inverse_product(change, inputs)
        actual = diagonal_matrix(ZZ, [2**sum(part > self.k for part in partition)
                                     for partition in shapes])*coefficients
        assert actual.base_ring() is ZZ
        return shapes, actual, [partition for _basis, partition in source]


def typed_matrices(a, k, degree, models=None):
    """Full actual D matrix in head Schur-square and p times smaller head bases."""
    if a < 2 or k % 2 or degree % 2:
        raise ValueError("Require a>=2 and even k,degree")
    square_model, euler_model = models or (BinaryChernModel(a, k), BinaryChernModel(a-1, k))
    shapes, square, plus_source = square_model.coefficient_matrix(degree, "D")
    smaller, euler, minus_source = euler_model.coefficient_matrix(degree-k, "C")
    indices = {partition: i for i, partition in enumerate(smaller)}
    plus_rows, minus_rows, full_rows, typed_rows = [], [], [], []
    for i, partition in enumerate(shapes):
        left, right = list(square.row(i)), [ZZ(0)]*euler.ncols()
        if k in partition:
            mu = list(partition)
            mu.remove(k)
            if tuple(mu) in indices:
                right = list(euler.row(indices[tuple(mu)]))
            for kind, sign in ((1, 1), (2, -1)):
                signed = [sign*entry for entry in right]
                plus_rows.append(left)
                minus_rows.append(signed)
                full_rows.append(left+signed)
                typed_rows.append((partition, kind))
        else:
            plus_rows.append(left)
            minus_rows.append(right)
            full_rows.append(left+right)
            typed_rows.append((partition, 0))
    return {"plus": matrix(ZZ, plus_rows), "minus": matrix(ZZ, minus_rows),
            "full": matrix(ZZ, full_rows), "square_partitions": plus_source,
            "Euler_partitions": minus_source, "typed_rows": typed_rows,
            "Y_partitions": shapes, "Y_square_coefficients": square,
            "Euler_B_partitions": smaller, "Euler_B_coefficients": euler}
