"""Run from borel/: PYTHONPATH=src sage -python examples/integral_images.py."""
from sage_borel import BorelMap, GeneratorAlgebra, classical_borel, exceptional_borel


def main():
    b = BorelMap("A3", omitted_nodes=1)
    h = b.schubert([1])
    A = GeneratorAlgebra(b, {"x": 2*b.weights[0]})
    x, = A.domain.gens()
    D = A.integral_degree(2, include_ambient=True)
    assert D.smith_factors == (4,)
    assert not D.contains(h**2) and D.contains(4*h**2)
    assert D.preimage(4*h**2) == x**2
    print("P3, generator 2h, degree two:", D.smith_factors)
    print("Integer preimage of 4h^2:", D.preimage(4*h**2))

    c = classical_borel("C", 3, 1)
    q = c.tail_polynomial(1, 1)
    normalized = GeneratorAlgebra(c, {"q": q})
    literal = GeneratorAlgebra(c, {"Q": 4*q})
    assert normalized.integral_degree(4).smith_factors == (1,)
    assert literal.integral_degree(4).smith_factors == (16,)
    print("C3 tail q, degree four:", normalized.integral_degree(4).smith_factors)
    print("Literal polynomial-character form 4q:", literal.integral_degree(4).smith_factors)

    d = classical_borel("D", 4, 2)
    AD = GeneratorAlgebra(d, {"q": d.tail_polynomial(2, 1),
                             "u": d.euler_polynomial(2)})
    DD = AD.integral_degree(4, include_ambient=True)
    assert DD.smith_factors == (1, 2, 2)
    assert DD.saturation_index == 4 and DD.free_cokernel_rank == 1
    print("D4/P2 tail plus Euler, degree four:", DD.smith_factors)
    print("Index in saturation / free cokernel rank:",
          DD.saturation_index, DD.free_cokernel_rank)

    f = exceptional_borel("F4", 4)
    AF = GeneratorAlgebra(f, {info.name: f.generator_polynomial(info.name)
                              for info in f.generators()})
    assert AF.integral_degree(4).smith_factors == (1, 12)
    print("F4/P4 catalogue, degree four:", AF.integral_degree(4).smith_factors)

    g = exceptional_borel("G2", 1)
    AG = GeneratorAlgebra(g, {info.name: g.generator_polynomial(info.name)
                              for info in g.generators()})
    try:
        AG.integral_degree(0)
    except ValueError as error:
        assert "nonintegral Schubert coefficient" in str(error)
        print("Original G2/P1 catalogue form correctly rejected:", error)
    else:
        raise AssertionError("nonintegral generator form was accepted")
    print("All integral-image examples passed.")


if __name__ == "__main__":
    main()
