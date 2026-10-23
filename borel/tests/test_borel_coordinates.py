"""Independent matrix/trace checks and exact fixed-t fiber tests."""
import unittest
from sage.all import QQ, matrix, identity_matrix, diagonal_matrix, vector
from sage_borel import borel_springer_map


def unit(n, i, j):
    result = matrix(QQ, n)
    result[i, j] = 1
    return result


def independent_data(group):
    """Construct group coordinates and full Lie basis without API matrices."""
    if group == "SL3":
        n = 3
        roots = [unit(n, 0, 1), unit(n, 1, 2), unit(n, 0, 2)]
        basis = [unit(n, i, j) for i in range(n) for j in range(n) if i != j]
        basis += [unit(n, i, i)-unit(n, i+1, i+1) for i in range(n-1)]
        return roots, basis, None
    n = 4
    J = matrix(QQ, [[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]])
    roots = [unit(n, 0, 1)-unit(n, 2, 3), unit(n, 1, 2),
             unit(n, 0, 2)+unit(n, 1, 3), unit(n, 0, 3)]
    # X belongs to sp(J) iff J X is symmetric: no closed projection formula used.
    basis = [J.inverse()*(unit(n, i, j)+(unit(n, j, i) if i != j else 0))
             for i in range(n) for j in range(i, n)]
    return roots, basis, J


class BorelCoordinateTests(unittest.TestCase):
    def test_symbolic_trace_projection(self):
        for group in ("SL3", "Sp4"):
            api = borel_springer_map(group)
            x, y = api.torus_parameters
            roots, basis, J = independent_data(group)
            n = roots[0].nrows()
            u = identity_matrix(api.field, n)
            for parameter, root in zip(api.unipotent_parameters, roots):
                u *= identity_matrix(api.field, n)+parameter*root
            diagonal = [x, y, 1/(x*y)] if n == 3 else [x, y, 1/y, 1/x]
            g = diagonal_matrix(diagonal)*u
            self.assertEqual(g, api.group_matrix())
            self.assertEqual(g.det(), 1)
            if J is not None:
                self.assertEqual(g.transpose()*J*g, J)
                for e in basis:
                    self.assertTrue((e.transpose()*J+J*e).is_zero())
            gram = matrix(QQ, [[(e*f).trace() for f in basis] for e in basis])
            coefficients = gram.solve_right(vector(api.field, [(g*e).trace() for e in basis]))
            projected = sum((c*e for c, e in zip(coefficients, basis)), matrix(api.field, n))
            self.assertEqual(projected, api.lie_matrix())
            direct = g-g.trace()/3*identity_matrix(api.field, n) if n == 3 else (g-g.inverse())/2
            self.assertEqual(projected, direct)

    def test_jacobian_and_root_factors(self):
        for group in ("SL3", "Sp4"):
            api = borel_springer_map(group)
            x, y = api.torus_parameters
            h = api.toral_coordinates()
            if group == "SL3":
                characters = [x/y, x*y*y, x*x*y]
                alpha_h = [h[0]-h[1], h[1]-h[2], h[0]-h[2]]
                expected = x*x*y
            else:
                characters = [x/y, y*y, x*y, x*x]
                alpha_h = [h[0]-h[1], 2*h[1], h[0]+h[1], 2*h[0]]
                expected = x*(x*y+1)*(x+y)/4
            jac = api.vertical_jacobian()
            self.assertEqual(jac.det(), expected)
            for i, (character, value) in enumerate(zip(characters, alpha_h)):
                self.assertEqual(jac[i, i], character*value/(character-1))
            self.assertEqual(api.vertical_rank((1, 1), [0]*api.dimension), api.dimension)
            self.assertEqual(api.vertical_rank((1, 1), [2]*api.dimension), api.dimension)

    def test_open_stratum_inverse(self):
        for group in ("SL3", "Sp4"):
            api = borel_springer_map(group)
            for torus in ((1, 1), (2, 3), (-2, -3)):
                for source in ([0]*api.dimension, list(range(1, api.dimension+1))):
                    output = api.coordinates(torus, source)
                    self.assertEqual(api.inverse(torus, output), tuple(source))

    def test_sp4_all_collapse_strata_and_fibers(self):
        api = borel_springer_map("Sp4")
        # Each fiber is verified identically in its free polynomial parameters.
        cases = [((2, QQ(-1)/2), (0, 3, 5, 7), 1),
                 ((2, -2), (3, 4, -8, 7), 1),
                 ((1, -1), (0, 3, 5, 7), 1),
                 ((-1, 1), (0, 3, 5, 7), 1),
                 ((1, -1), (0, 0, 0, 7), 2),
                 ((-1, 1), (0, 0, 0, 7), 2)]
        formulas = api.coordinates()
        for torus, target, dimension in cases:
            fiber = api.fiber(torus, target)
            self.assertIsNotNone(fiber)
            self.assertEqual(fiber['dimension'], dimension)
            base = fiber['parameters'][0].parent()
            values = tuple(map(base, torus))+tuple(map(base, fiber['coordinates']))
            evaluate = api.ring.hom(values, base)
            for f, v in zip(formulas, target):
                self.assertEqual(evaluate(f.numerator())/evaluate(f.denominator()), v)
            # Rank over the free-parameter field, not a numeric sample.
            J = api.vertical_jacobian()
            J = matrix(base.fraction_field(), [[evaluate(f.numerator())/evaluate(f.denominator())
                                                for f in row] for row in J.rows()])
            self.assertEqual(J.rank(), 4-dimension)
        for torus, target in [((2, QQ(-1)/2), (1, 0, 0, 0)),
                              ((2, -2), (3, 4, 0, 7)),
                              ((1, -1), (0, 0, 1, 0)),
                              ((-1, 1), (0, 0, 1, 0))]:
            self.assertIsNone(api.fiber(torus, target))

    def test_symbolic_stratum_ranks(self):
        api = borel_springer_map("Sp4")
        x, y = api.torus_parameters
        a, b, c, d = api.unipotent_parameters
        J = api.vertical_jacobian()
        def specialize(subs):
            return matrix(api.field, [[f.subs(subs) for f in row] for row in J.rows()])
        self.assertEqual(specialize({y: -1/x}).rank(), 3)
        self.assertEqual(specialize({y: -x}).rank(), 3)
        for sign in (-1, 1):
            crossing = specialize({x: sign, y: -sign})
            self.assertEqual(crossing.rank(), 3)
            self.assertEqual(matrix(api.field, [[f.subs({b: 0}) for f in row]
                                               for row in crossing.rows()]).rank(), 2)

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            borel_springer_map("G2")
        api = borel_springer_map("Sp4")
        for torus, u in [((0, 1), (0, 0, 0, 0)), ((1,), (0, 0, 0, 0)),
                          ((1, 1), (0, 0, 0))]:
            with self.assertRaises(ValueError):
                api.coordinates(torus, u)
        with self.assertRaises(ValueError):
            api.inverse((1, -1), (0, 0, 0, 0))
        with self.assertRaises(ValueError):
            api.coordinates((1, 1))


if __name__ == "__main__":
    unittest.main()
