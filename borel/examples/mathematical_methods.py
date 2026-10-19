"""Mathematical examples of bounded images and complete kernel presentations.

Run from borel/: PYTHONPATH=src sage -python examples/mathematical_methods.py

All grading degrees printed here are polynomial degrees; double them for
cohomological degrees. These are executable examples accompanying METHODS.md.
"""

from sage_borel import BorelMap, GeneratorAlgebra, two_step_symplectic


def projective_plane():
    """H*(P^2; QQ) = QQ[h]/(h^3), with h in cohomological degree two."""
    borel = BorelMap("A2", omitted_nodes=1)
    algebra = GeneratorAlgebra(borel, {"h": borel.weights[0]})
    h, = algebra.domain.gens()

    # Bounded reports use Schubert-coordinate columns of a linear map.
    degree_two = algebra.degree(2)
    assert degree_two.monomials == (h**2,)
    assert degree_two.matrix.ncols() == 1
    assert degree_two.rank == 1
    assert degree_two.kernel_dimension == 0
    assert algebra.degree(3).relations == (h**3,)

    # A presentation proves the full ideal equality, beyond a degree cutoff.
    presentation = algebra.presentation()
    assert presentation.relations == (h**3,)
    assert presentation.standard_monomials == (1, h, h**2)
    assert presentation.dimension == 3
    assert presentation.reduce(h**4 + 2*h**2) == 2*h**2
    assert presentation.verify().certified
    assert presentation.verify(backend="bgg").certified
    print("P^2: degree-2 image matrix {}, rank {}; kernel {}; Hilbert {}".format(
        degree_two.matrix, degree_two.rank, presentation.relations,
        presentation.hilbert_series))


def two_step_symplectic_image():
    """The chosen middle-plus-tail image is smaller than all flag cohomology."""
    family = two_step_symplectic(n=3, b=1, k=1)
    algebra = family.joint_algebra()
    c, e = algebra.domain.gens()
    presentation = algebra.presentation()

    assert presentation.generator_degrees == (1, 2)
    assert set(presentation.relations) == {c**4 + c**2*e + e**2, e**3}
    assert presentation.dimension == 12
    assert set(presentation.standard_monomials) == {
        c**i * e**j for i in range(4) for j in range(3)}
    assert presentation.hilbert_coefficients == (1, 1, 2, 2, 2, 2, 1, 1)
    assert presentation.reduce(c**4) == -c**2*e - e**2
    assert presentation.reduce(e**3) == 0
    assert presentation.verify().certified
    assert presentation.verify(backend="bgg").certified
    assert (presentation.kernel_ideal()
            == algebra.presentation(backend="bgg").kernel_ideal())
    print("C3 two-step: weights {}; image dimension {} (ambient dimension 24); "
          "kernel {}; Hilbert {}".format(
              presentation.generator_degrees, presentation.dimension,
              presentation.relations, presentation.hilbert_series))


def rank_nine_tail():
    """Relations in the invariant subring can be computed in the full quotient."""
    borel = BorelMap("C9", omitted_nodes=8)
    # In type C, x9 = omega9 - omega8, so q = x9^2 has degree two.
    x9 = borel.weights[8] - borel.weights[7]
    assert borel.tail_polynomial(1, 1) == x9**2
    algebra = GeneratorAlgebra(borel, {"q": x9**2})
    q, = algebra.domain.gens()
    presentation = algebra.presentation()

    assert presentation.backend == "coinvariant"
    assert presentation.generator_degrees == (2,)
    assert presentation.relations == (q**9,)
    assert presentation.standard_monomials == tuple(q**i for i in range(9))
    assert presentation.dimension == 9
    assert presentation.hilbert_coefficients == tuple(
        int(degree % 2 == 0) for degree in range(17))
    # This defaults to a fresh coinvariant certificate; do not run BGG at rank 9.
    assert presentation.verify().certified
    print("C9/P8 tail: backend {}; kernel {}; dimension {}; Hilbert {}".format(
        presentation.backend, presentation.relations, presentation.dimension,
        presentation.hilbert_series))


def orthogonal_euler_generator():
    """The type-D Euler product is an additional, necessary invariant."""
    borel = BorelMap("D4", omitted_nodes=2)
    _, w2, w3, w4 = borel.weights
    x3, x4 = w3 + w4 - w2, w4 - w3
    q_input, u_input = x3**2 + x4**2, x3*x4
    assert q_input == borel.tail_polynomial(2, 1)
    assert u_input**2 == borel.tail_polynomial(2, 2)
    algebra = GeneratorAlgebra(borel, {"q": q_input, "u": u_input})
    q, u = algebra.domain.gens()
    # q and u both have polynomial degree two, and remain independent images.
    assert algebra.generator_degrees == (2, 2)
    assert algebra.degree(2).rank == 2
    presentation = algebra.presentation()
    assert presentation.reduce(u) != 0
    assert presentation.reduce(u**2) != 0
    assert presentation.reduce(q - u) != 0
    assert presentation.verify().certified
    assert presentation.verify(backend="bgg").certified
    print("D4/P2 tail: q=e1(x3^2,x4^2), u=x3*x4; u^2=e2(x3^2,x4^2); "
          "degree-2 rank 2; kernel {}; dimension {}; Hilbert {}".format(
              presentation.relations, presentation.dimension,
              presentation.hilbert_series))


def main():
    projective_plane()
    two_step_symplectic_image()
    rank_nine_tail()
    orthogonal_euler_generator()


if __name__ == "__main__":
    main()
