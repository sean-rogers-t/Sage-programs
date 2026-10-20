"""Run from borel/: PYTHONPATH=src sage -python examples/schubert_calculus.py.

Exact rational cup products, polynomial lifts and intersection numbers.
Polynomial degree d corresponds to cohomological degree 2*d throughout.
"""

from sage_borel import BorelMap


def main():
    # A4/P2 is Gr(2,5), of complex dimension 6.  The two head coordinates
    # describe the rank-two block; symmetric polynomials in them are literal
    # Levi invariants and hence valid inputs to its Borel map.
    grass = BorelMap("A4", omitted_nodes=2)
    x1, x2 = grass.classical_coordinates()[:2]
    e2 = x1*x2
    h2 = x1**2 + x1*x2 + x2**2
    unit = grass.schubert(())
    hyperplane = grass.schubert((2,))
    row_two = grass.expand(h2)
    assert len(row_two) == 1

    # This uses two nondivisor degree-two classes; the product has degree
    # four, strictly below the top degree six.  The Littlewood-Richardson
    # identity s_(2)^2=s_(4)+s_(3,1)+s_(2,2) loses s_(4), which lies outside
    # the 2-by-3 rectangle.  In two variables s_(3,1)=e2*h2 and s_(2,2)=e2^2.
    product = grass.multiply(row_two, row_two)
    assert len(product) == 2
    assert product == grass.expand(e2*h2 + e2**2)
    index = next(iter(row_two))
    assert product == grass.schubert_product(index, index) == row_two**2
    print("Gr(2,5), sigma_(2):", row_two)
    print("sigma_(2)^2 = sigma_(3,1) + sigma_(2,2):", product)

    # A lift is a chosen polynomial representative, not a unique inverse.
    # Its characteristic-map image is the requested class exactly.
    lifted = grass.schubert_polynomial(index)
    assert grass.expand(lifted) == row_two
    print("Chosen polynomial lift of sigma_(2):", lifted)

    mixed = unit + 3*hyperplane + 2*row_two
    mixed_product = grass.multiply(mixed, unit + row_two)
    assert mixed_product == mixed * (unit + row_two)
    assert grass.expand(grass.lift(mixed)) == mixed
    assert mixed_product == grass.expand(grass.lift(mixed)
                                        * grass.lift(unit + row_two))
    print("Mixed-class product:", mixed_product)

    # The Pluecker degree of Gr(2,5) is 5: h^6=5[point].  Powers use cup
    # multiplication, while the final two-factor integral uses duality.
    top_power = hyperplane**6
    assert grass.integrate(top_power) == 5
    assert grass.intersection_number(hyperplane**2, hyperplane**2,
                                     hyperplane**2) == 5
    assert grass.integrate(unit + top_power) == 5
    print("h^6:", top_power)
    print("Integral of h^6:", grass.integrate(top_power))

    # Low-degree inverse calculations also work in exceptional types.
    # The degreewise linear method avoids expanding E7's global discriminant.
    exceptional = BorelMap("E7", omitted_nodes=7)
    e7_index = exceptional.basis(2)[0]
    e7_class = exceptional.schubert(e7_index)
    e7_lift = exceptional.schubert_polynomial(e7_index, method="linear")
    assert exceptional.expand(e7_lift) == e7_class
    print("E7/P7 degree-two Schubert class:", e7_class)
    print("Its polynomial lift:", e7_lift)
    print("All Schubert-calculus examples passed.")


if __name__ == "__main__":
    main()
