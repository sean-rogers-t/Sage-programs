"""Independent representation matrices versus torus/commutator reconstruction."""
import unittest
from sage.all import QQ, PolynomialRing, QuadraticField, matrix, vector, identity_matrix, diagonal_matrix, prod
from sage_borel import borel_springer_map, reconstruct_borel_map, reconstruct_root_coordinates
from test_borel_coordinates import independent_data, unit


def tensor_square(g, alternating=False):
    """Direct action on e_i e_j or e_i wedge e_j, by expanding columns."""
    n = g.nrows()
    pairs = [(i, j) for i in range(n) for j in range(i+int(alternating), n)]
    def entry(k, l, i, j):
        if alternating:
            return g[k, i]*g[l, j]-g[k, j]*g[l, i]
        return g[k, i]*g[l, j]+(g[l, i]*g[k, j] if k != l else 0)
    return matrix(g.base_ring(), [[entry(k, l, i, j) for i, j in pairs] for k, l in pairs])


def representation_model(g, basis, kind):
    """Return rho(g) and d rho of a full Lie basis, without characters."""
    if kind == 'defining':
        return g, basis
    if kind == 'dual':
        return g.inverse().transpose(), [-e.transpose() for e in basis]
    if kind == 'adjoint':
        flattened = matrix(QQ, [e.list() for e in basis]).transpose()
        rows = flattened.transpose().pivots()
        inverse = flattened.matrix_from_rows(rows).inverse()
        def coordinates(e):
            return inverse*vector(e.base_ring(), [e.list()[i] for i in rows])
        adjoint = matrix(g.base_ring(), [coordinates(g*e*g.inverse()) for e in basis]).transpose()
        differentials = [matrix(QQ, [coordinates(e*f-f*e) for f in basis]).transpose() for e in basis]
        return adjoint, differentials
    if kind in ('symmetric_square', 'exterior_square'):
        alternating = kind == 'exterior_square'
        R = PolynomialRing(QQ, 'v')
        v = R.gen()
        differentials = []
        for e in basis:
            lifted = tensor_square(identity_matrix(R, g.nrows())+v*e, alternating)
            differentials.append(matrix(QQ, [[f.derivative(v)(0) for f in row] for row in lifted.rows()]))
        return tensor_square(g, alternating), differentials
    raise ValueError(kind)


def independent_projection(api, kinds):
    roots, basis, J = independent_data(api.group)
    n = roots[0].nrows()
    x, y = api.torus_parameters
    u = identity_matrix(api.field, n)
    for parameter, root in zip(api.unipotent_parameters, roots):
        u *= identity_matrix(api.field, n)+parameter*root
    g = diagonal_matrix([x, y, 1/(x*y)] if n == 3 else [x, y, 1/y, 1/x])*u
    gram = matrix(QQ, len(basis))
    rhs = vector(api.field, [0]*len(basis))
    for kind in kinds:
        rho, differentials = representation_model(g, basis, kind)
        gram += matrix(QQ, [[(e*f).trace() for f in differentials] for e in differentials])
        rhs += vector(api.field, [(rho*e).trace() for e in differentials])
    coefficients = gram.solve_right(rhs)
    return sum((c*e for c, e in zip(coefficients, basis)), matrix(api.field, n))


class BorelReconstructionTests(unittest.TestCase):
    def test_low_level_other_root_data(self):
        # Independent A1 and A3 input, without either ready coordinate backend.
        for n in (2, 4):
            pairs = sorted(((i, j) for i in range(n) for j in range(i+1, n)),
                           key=lambda pair: (pair[1]-pair[0], pair[0]))
            roots = [tuple(int(i <= k < j) for k in range(n-1)) for i, j in pairs]
            names = tuple('z%d' % i for i in range(n-1))+tuple('u%d' % i for i in range(len(pairs)))
            R = PolynomialRing(QQ, names=names)
            K = R.fraction_field()
            generators = tuple(map(K, R.gens()))
            diagonal = generators[:n-1]+(1/prod(generators[:n-1]),)
            parameters = generators[n-1:]
            matrices = [unit(n, i, j) for i, j in pairs]
            brackets = {}
            for i in range(len(pairs)):
                for j in range(i+1, len(pairs)):
                    bracket = matrices[i]*matrices[j]-matrices[j]*matrices[i]
                    if not bracket.is_zero():
                        nonzero = [(k, l, bracket[k, l]) for k in range(n) for l in range(n)
                                   if bracket[k, l] != 0]
                        self.assertEqual(len(nonzero), 1)
                        brackets[i, j] = nonzero[0][2]
            result = reconstruct_root_coordinates(roots, brackets,
                [diagonal[i]/diagonal[j] for i, j in pairs],
                [diagonal[i]-diagonal[i+1] for i in range(n-1)], parameters)
            g = diagonal_matrix(diagonal)
            for u, e in zip(parameters, matrices):
                g *= identity_matrix(K, n)+u*e
            self.assertEqual(result, tuple(g[i, j] for i, j in pairs))

    def test_defining_recovery_and_commutators(self):
        for group in ('SL3', 'Sp4'):
            old = borel_springer_map(group)
            new = reconstruct_borel_map(group)
            self.assertEqual(old.lie_matrix(), new.lie_matrix())
            roots, _, _ = independent_data(group)
            expected = {(0, 1): (2, 1)}
            if group == 'Sp4':
                expected[0, 2] = (3, 2)
            for i in range(len(roots)):
                for j in range(i+1, len(roots)):
                    bracket = roots[i]*roots[j]-roots[j]*roots[i]
                    if (i, j) in expected:
                        k, c = expected[i, j]
                        self.assertEqual(bracket, c*roots[k])
                    else:
                        self.assertTrue(bracket.is_zero())
            custom = reconstruct_borel_map(group, torus_map=lambda x, y: old.toral_coordinates()[:2])
            self.assertEqual(custom.coordinates(), old.coordinates())

    def test_all_parameter_independent_representation_projections(self):
        cases = [('SL3', r, [r]) for r in ('defining', 'dual', 'symmetric_square', 'adjoint')]
        cases += [('Sp4', 'defining', ['defining']),
                  ('Sp4', 'second_fundamental', ['exterior_square']),
                  ('Sp4', 'adjoint', ['adjoint'])]
        for group, name, models in cases:
            with self.subTest(group=group, representation=name):
                api = reconstruct_borel_map(group, name)
                self.assertEqual(independent_projection(api, models), api.lie_matrix())
                self.assertTrue(api.denominators_cancelled)
        # Exterior square of Sp4 is omega_2 plus a trivial summand. The
        # latter has zero differential, so it does not change trace projection.
        for group in ('SL3', 'Sp4'):
            d = reconstruct_borel_map(group)
            a = reconstruct_borel_map(group, 'adjoint')
            api = reconstruct_borel_map(group, character=lambda x, y: d.character+a.character)
            self.assertEqual(independent_projection(api, ['defining', 'adjoint']), api.lie_matrix())

    def test_regular_extension_and_normalization(self):
        for group, names in [('SL3', ('defining', 'dual', 'symmetric_square', 'adjoint')),
                             ('Sp4', ('defining', 'second_fundamental', 'adjoint', 'symmetric_square'))]:
            for name in names:
                api = reconstruct_borel_map(group, name)
                self.assertEqual(api.toral_coordinates((1, 1)), (0,)*(3 if group == 'SL3' else 4))
                self.assertEqual(api.vertical_jacobian((1, 1), [0]*api.dimension),
                                 identity_matrix(QQ, api.dimension))
                # Includes e^alpha(t)=1 and the old Sp4 crossing torus points.
                for t in ((1, 1), (2, 2), (1, -1), (-1, 1)):
                    coordinates = api.coordinates(t, [2]*api.dimension)
                    self.assertTrue(all(c in QQ for c in coordinates))
                self.assertEqual(api.vertical_jacobian().det(),
                                 prod(api.root_factors()))

    def test_triangular_inverse_other_representations(self):
        for group in ('SL3', 'Sp4'):
            for name in ('defining', 'dual', 'symmetric_square', 'adjoint'):
                api = reconstruct_borel_map(group, name)
                point = tuple(range(1, api.dimension+1))
                for t in ((1, 1), (2, 3)):
                    self.assertEqual(api.inverse(t, api.coordinates(t, point)), point)

    def test_changed_collapse_loci(self):
        # SL3 adjoint has new collapse even though the defining map has none.
        api = reconstruct_borel_map('SL3', 'adjoint')
        x, y = api.torus_parameters
        a, b, c = api.unipotent_parameters
        self.assertEqual(tuple(f.subs({x: 1, y: -1}) for f in api.coordinates()),
                         (0, b/3, -a*b/6))
        self.assertEqual(api.vertical_rank((1, -1), (0, 0, 0)), 1)
        self.assertEqual(api.vertical_rank((1, -1), (2, 3, 4)), 2)
        self.assertEqual(borel_springer_map('SL3').vertical_rank((1, -1), (0, 0, 0)), 3)
        with self.assertRaises(ValueError):
            api.inverse((1, -1), (0, 0, 0))
        # Sp4 omega_2 has additional long-root zeros at z1^2=-1.
        api = reconstruct_borel_map('Sp4', 'second_fundamental')
        K = QuadraticField(-1, 'i')
        for obj, rank in ((api, 3), (borel_springer_map('Sp4'), 4)):
            evaluate = obj.ring.hom([K.gen(), K(2), K(0), K(0), K(0), K(0)], K)
            jac = matrix(K, [[evaluate(f.numerator())/evaluate(f.denominator())
                              for f in row] for row in obj.vertical_jacobian().rows()])
            self.assertEqual(jac.rank(), rank)

    def test_character_scaling_and_trivial_summands(self):
        for group in ('SL3', 'Sp4'):
            original = reconstruct_borel_map(group, 'adjoint')
            repeated = reconstruct_borel_map(group, character=lambda x, y: 2*original.character+7)
            self.assertEqual(repeated.coordinates(), original.coordinates())
            self.assertEqual(repeated.toral_coordinates(), original.toral_coordinates())

    def test_input_guards(self):
        for kwargs in [dict(character=lambda x, y: 3), dict(character=lambda x, y: x),
                       dict(torus_map=lambda x, y: (1, 0)),
                       dict(character=lambda x, y: x+y, torus_map=lambda x, y: (0, 0)),
                       dict(representation='adjoint', character=lambda x, y: x)]:
            with self.assertRaises(ValueError):
                reconstruct_borel_map('SL3', **kwargs)
        with self.assertRaises(ValueError):
            reconstruct_borel_map('Sp4', 'unknown')
        api = reconstruct_borel_map('SL3')
        with self.assertRaises(ValueError):
            reconstruct_root_coordinates(((1, 0),), {}, (1,), (1, 2), api.unipotent_parameters[:1])


if __name__ == '__main__':
    unittest.main()
