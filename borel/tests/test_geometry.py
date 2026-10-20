"""Exact quotient geometry, root directions, and generic rather than special fibers."""
import unittest
from sage.all import QQ, matrix, identity_matrix
from sage_borel import (SpringerGeometry, SpringerCharacterMap,
                       springer_character_map, g2_adjoint_quotient,
                       quadratic_character_pullback)


class GeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = springer_character_map("G2", (1, 0))
        cls.g = SpringerGeometry(cls.s)

    def test_minimal_g2_compact_quotient(self):
        quotient = g2_adjoint_quotient()
        x, y = quotient.character_polynomials[0].parent().gens()
        A = 2*x*x+13*x+y+13
        D = x*x+2*x-4*y-7
        self.assertEqual(quotient.degrees, (2, 6))
        self.assertEqual(quotient.character_polynomials,
                         ((x*x-4*x-y-7)/12, D*A*A/46656))
        self.assertEqual(quotient.character_jacobian().det(),
                         -(y-2)*(x+1)*A/10368)
        self.assertEqual(quotient.generic_degree(), 6)

    def test_adjoint_g2_degree_and_normalization(self):
        quotient = g2_adjoint_quotient(representation="adjoint")
        self.assertEqual(quotient.generic_degree(), 18)
        s = springer_character_map("G2", (0, 1))
        g = SpringerGeometry(s)
        psi = g.character_coordinates(g.cayley_determinant())
        self.assertEqual(psi(7, 14), 1)
        x, y = psi.parent().gens()
        H = -9*x**3-5*x*x+26*x*y+21*x+30*y+25
        L = x**3+3*x*x-4*x*y+2*y*y-5*x+4*y+1
        S = 9*x**3-27*x*x*y-11*x*x-18*x*y+18*y*y+23*x+17*y-13
        self.assertEqual(psi, -H*L*L*S*S/9618527719784448)

    def test_full_minimal_determinant_and_root_division(self):
        s, g = self.s, self.g
        J = g.torus_jacobian()
        self.assertEqual(matrix(QQ, [[f(1, 1) for f in row] for row in J]), identity_matrix(2))
        x, y = g.character_ring.gens()
        A = 2*x*x+13*x+y+13
        self.assertEqual(g.character_coordinates(J.det()), (y-2)/12)
        psi = g.cayley_determinant()
        self.assertEqual(g.character_coordinates(psi), (y-2)*(x+1)**2*A*A/35831808)
        self.assertEqual(psi(1, 1), 1)
        self.assertNotEqual(psi, J.det())
        for a, r in g.root_factors():
            from sage.all import prod
            e = prod(t**k for t, k in zip(s.torus_ring.gens(), a))
            self.assertEqual((e-1)*r, sum(k*f for k, f in zip(a, s.theta)))

    def test_sl2_and_type_a_two(self):
        s = springer_character_map("A1", (1,))
        g = SpringerGeometry(s)
        t = s.torus_ring.gen()
        self.assertEqual(g.cayley_determinant(), (t+t**-1)/2)
        self.assertEqual(g.adjoint_quotient((s.ring.gen()**2,)).generic_degree(), 2)
        s = springer_character_map("A2", (1, 0))
        g = SpringerGeometry(s)
        w, v = s.ring.gens()
        quotient = g.adjoint_quotient((w*w-w*v+v*v, w*v*(w-v)))
        self.assertEqual(quotient.generic_degree(), 3)
        self.assertEqual(g.cayley_determinant()(1, 1), 1)

    def test_character_coordinates_and_repeated_representation(self):
        x, y = self.g.character_ring.gens()
        self.assertEqual(self.g.character_coordinates(self.s.torus_ring(3), max_steps=1), 3)
        self.assertEqual(self.g.character_coordinates(self.s.torus_ring(0)), 0)
        doubled = SpringerGeometry(SpringerCharacterMap("G2", [(1, 0), (1, 0)]))
        self.assertEqual(doubled.cayley_determinant(), self.g.cayley_determinant())

    def test_uniform_quadratic_casimir_identity(self):
        from sage.all import WeylCharacterRing
        for kind, weight in (("A1", (1,)), ("G2", (1, 0)), ("G2", (0, 1)),
                             ("F4", (0, 0, 0, 1)), ("E6", (1, 0, 0, 0, 0, 0)),
                             ("E7", (0, 0, 0, 0, 0, 0, 1))):
            s = springer_character_map(kind, weight)
            from sage_borel import irreducible_character_lattice
            info = irreducible_character_lattice(kind, weight)
            from sage.all import vector
            w = vector(s.ring.gens())
            q = (w*(info.trace_gram/info.trace_index)*w)/2
            actual = SpringerGeometry(s).character_coordinates(s.pullback(q))
            WCR = WeylCharacterRing(kind, base_ring=QQ, style="coroots")
            fundamentals = [WCR(tuple(int(i == j) for j in range(s.rank))) for i in range(s.rank)]
            self.assertEqual(actual(*fundamentals), quadratic_character_pullback(kind, weight))
        with self.assertRaises(ValueError):
            quadratic_character_pullback("E7", (1, 0, 0, 0, 0, 0, 0), max_dimension=100)

    def test_validation_and_budgets(self):
        g, s = self.g, self.s
        with self.assertRaises(TypeError):
            SpringerGeometry("G2")
        with self.assertRaises(ValueError):
            g.character_coordinates(s.torus_ring.gen(0))
        for invariants in ((s.ring.gen(0),), (s.ring.gen(0)**2, s.ring.gen(1)**6)):
            with self.assertRaises(ValueError):
                g.adjoint_quotient(invariants)
        for kwargs in ({"max_terms": 1}, {"max_rank": 1}, {"max_rank": True}):
            with self.assertRaises(ValueError):
                g.cayley_determinant(**kwargs)
        f = g2_adjoint_quotient().torus_pullbacks[0]
        for kwargs in ({"max_terms": 1}, {"max_steps": 1}, {"max_character_dimension": 1}):
            with self.assertRaises(ValueError):
                g.character_coordinates(f, **kwargs)
        with self.assertRaises(ValueError):
            g2_adjoint_quotient().generic_degree(max_rank=1)
        with self.assertRaises(ValueError):
            g2_adjoint_quotient(representation="unknown")


if __name__ == "__main__":
    unittest.main()
