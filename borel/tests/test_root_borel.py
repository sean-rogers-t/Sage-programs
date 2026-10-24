"""Automatic root/character adapter: independent matrix and trace routes."""
from functools import lru_cache
import unittest
from sage.all import (QQ, matrix, vector, diagonal_matrix, identity_matrix,
                      RootSystem, prod)
from sage_borel import root_borel_map, SpringerCharacterMap, borel_springer_map
from test_borel_coordinates import unit, independent_data
from test_borel_reconstruction import representation_model


def exponential_nilpotent(X):
    result = identity_matrix(X.base_ring(), X.nrows(), sparse=True)
    term = result
    for k in range(1, X.nrows()+1):
        term = term*X/k
        if term.is_zero():
            return result
        result += term
    raise AssertionError('expected a nilpotent matrix')


def classical_model(api, torus):
    """Standard A/B/C/D matrices; signs of nonsimple vectors aligned explicitly."""
    family, n = api.cartan_type.type(), api.rank
    size = n+1 if family == 'A' else 2*n+int(family == 'B')
    E, F = [], []
    if family == 'A':
        E = [unit(size, i, i+1) for i in range(n)]
        F = [e.transpose() for e in E]
        diagonal = [torus[0]]+[torus[i]/torus[i-1] for i in range(1, n)]+[1/torus[-1]]
    else:
        offset = n+int(family == 'B')
        for i in range(n-1):
            e = unit(size, i, i+1)-unit(size, offset+i+1, offset+i)
            E.append(e)
            F.append(e.transpose())
        if family == 'C':
            E.append(unit(size, n-1, 2*n-1))
            F.append(E[-1].transpose())
            eps = [torus[0]]+[torus[i]/torus[i-1] for i in range(1, n)]
        elif family == 'B':
            E.append(unit(size, n-1, n)-unit(size, n, 2*n))
            F.append(2*E[-1].transpose())
            eps = [torus[0]]+[torus[i]/torus[i-1] for i in range(1, n-1)]
            eps.append(torus[-1]**2/torus[-2])
        else:
            E.append(unit(size, n-2, 2*n-1)-unit(size, n-1, 2*n-2))
            F.append(E[-1].transpose())
            eps = [torus[0]]+[torus[i]/torus[i-1] for i in range(1, n-2)]
            eps += [torus[-2]*torus[-1]/torus[-3], torus[-1]/torus[-2]]
        diagonal = eps+([QQ.one()] if family == 'B' else [])+[1/v for v in eps]
    H = [e*f-f*e for e, f in zip(E, F)]
    index = {r: i for i, r in enumerate(api.positive_roots)}
    roots = []
    for alpha in api.positive_roots:
        if sum(alpha) == 1:
            roots.append(E[alpha.index(1)])
            continue
        for (i, j), c in api.brackets.items():
            if tuple(a+b for a, b in zip(api.positive_roots[i], api.positive_roots[j])) == alpha:
                roots.append((roots[i]*roots[j]-roots[j]*roots[i])/c)
                break
        else:
            raise AssertionError('positive nonsimple root has no decomposition')
    # Every bracket is checked, not just the decomposition used to align signs.
    for i, a in enumerate(api.positive_roots):
        for j in range(i+1, len(roots)):
            bracket = roots[i]*roots[j]-roots[j]*roots[i]
            k = index.get(tuple(x+y for x, y in zip(a, api.positive_roots[j])))
            expected = roots[k]*api.brackets[i, j] if k is not None else matrix(QQ, size)
            assert bracket == expected
    for i, h in enumerate(H):
        for j, e in enumerate(roots):
            assert h*e-e*h == api.root_weights[j][i]*e
    return diagonal_matrix(diagonal), roots, H


@lru_cache(maxsize=8)
def adjoint_model(kind):
    """Full Lie-algebra matrices, separate from the Borel recurrence/arrows."""
    from sage.all import LieAlgebra
    L = LieAlgebra(QQ, cartan_type=list(kind))
    keys = tuple(L.basis().keys())
    index = {key: i for i, key in enumerate(keys)}
    basis = L.basis()
    def ad(x):
        entries = {}
        for j, key in enumerate(keys):
            for out, c in x.bracket(basis[key]).monomial_coefficients().items():
                entries[index[out], j] = c
        return matrix(QQ, len(keys), len(keys), entries, sparse=True)
    matrices = {key: ad(basis[key]) for key in keys}
    H = [ad(L.e(i).bracket(L.f(i))) for i in range(1, kind[1]+1)]
    return keys, matrices, H


def trace_pair(A, B):
    return sum(v*B[j, i] for (i, j), v in A.dict().items())


def adjoint_data(api, torus):
    keys, matrices, H = adjoint_model((api.cartan_type.type(), api.rank))
    R = RootSystem(api.cartan_type).root_lattice()
    simple = R.simple_roots()
    roots = [sum(c*simple[i+1] for i, c in enumerate(a)) for a in api.positive_roots]
    weights = {r: w for r, w in zip(roots, api.root_weights)}
    weights.update({-r: tuple(-v for v in w) for r, w in zip(roots, api.root_weights)})
    diagonal = [prod(t**e for t, e in zip(torus, weights[key])) if key in weights else QQ.one()
                for key in keys]
    return diagonal_matrix(diagonal, sparse=True), [matrices[r] for r in roots], [matrices[-r] for r in roots], H


def product_unipotent(roots, unipotent):
    result = identity_matrix(QQ, roots[0].nrows(), sparse=True)
    for e, u in zip(roots, unipotent):
        if u:
            result *= exponential_nilpotent(u*e)
    return result


class RootBorelTests(unittest.TestCase):
    def test_previous_coordinate_conventions(self):
        for kind, group in [(('A', 2), 'SL3'), (('C', 2), 'Sp4')]:
            api = root_borel_map(kind, (1, 0))
            torus = (QQ(2), QQ(3))
            T, roots, H = classical_model(api, torus)
            if group == 'Sp4':
                P = identity_matrix(QQ, 4)
                P.swap_rows(2, 3)
                roots = [P*e*P for e in roots]
                H = [P*h*P for h in H]
            old_roots, _, _ = independent_data(group)
            scales = []
            for new, old in zip(roots, old_roots):
                ij = next(ij for ij, c in old.dict().items() if c)
                c = new[ij]/old[ij]
                self.assertEqual(new, c*old)
                scales.append(c)
            u = tuple(map(QQ, range(1, api.dimension+1)))
            actual = sum((c*h for c, h in zip(api.toral_coordinates(torus), H)), matrix(QQ, T.nrows()))
            actual += sum(c*e for c, e in zip(api.coordinates(torus, u), roots))
            expected = borel_springer_map(group).lie_matrix((torus[0], torus[1]/torus[0]),
                                                           tuple(c*v for c, v in zip(scales, u)))
            self.assertEqual(actual, expected)

    def test_classical_independent_matrices(self):
        for kind in (('A', 3), ('B', 3), ('C', 3), ('D', 4)):
            api = root_borel_map(kind, (1,)+(0,)*(kind[1]-1))
            for torus in (tuple(map(QQ, [2, 3, 5, 7][:kind[1]])),
                          (QQ.one(),)*kind[1], (QQ(1), QQ(-1))+(QQ.one(),)*(kind[1]-2)):
                T, roots, H = classical_model(api, torus)
                u = tuple(QQ((i % 3)-1) for i in range(api.dimension))
                g = T*product_unipotent(roots, u)
                expected = g-g.trace()/g.nrows()*identity_matrix(QQ, g.nrows()) if kind[0] == 'A' else (g-g.inverse())/2
                observed = sum((c*h for c, h in zip(api.toral_coordinates(torus), H)), matrix(QQ, g.nrows()))
                observed += sum(c*e for c, e in zip(api.coordinates(torus, u), roots))
                self.assertEqual(observed, expected)

    def test_larger_symmetric_square_trace(self):
        api = root_borel_map(('A', 3), (2, 0, 0))
        for t in ((QQ(2), QQ(3), QQ(5)), (QQ.one(),)*3):
            T, roots, H = classical_model(api, t)
            u = tuple(map(QQ, range(1, api.dimension+1)))
            g = T*product_unipotent(roots, u)
            basis = H+roots+[e.transpose() for e in roots]
            rho, differentials = representation_model(g, basis, 'symmetric_square')
            gram = matrix(QQ, [[trace_pair(e, f) for f in differentials] for e in differentials])
            rhs = vector(QQ, [trace_pair(rho, e) for e in differentials])
            expected = sum((c*e for c, e in zip(gram.solve_right(rhs), basis)), matrix(QQ, g.nrows()))
            actual = sum((c*h for c, h in zip(api.toral_coordinates(t), H)), matrix(QQ, g.nrows()))
            actual += sum(c*e for c, e in zip(api.coordinates(t, u), roots))
            self.assertEqual(actual, expected)

    def test_exceptional_adjoint_full_trace(self):
        for kind in (('G', 2), ('F', 4), ('E', 6), ('E', 7), ('E', 8)):
            R = RootSystem(list(kind))
            highest = max(R.root_lattice().positive_roots(), key=lambda r: sum(r.to_vector()))
            weight = tuple(R.cartan_type().cartan_matrix()*highest.to_vector())
            api = root_borel_map(kind, weight)
            u = tuple(QQ(i+1) if i < api.rank else QQ(int(i == api.dimension-1))
                      for i in range(api.dimension))
            for t in ((QQ.one(),)*api.rank,
                      (QQ(1), QQ(2))+(QQ.one(),)*(api.rank-2),
                      tuple(map(QQ, [2, 3, 5, 7, 11, 13, 17, 19][:api.rank]))):
                T, positive, negative, H = adjoint_data(api, t)
                g = T*product_unipotent(positive, u)
                expected = tuple(trace_pair(g, f)/trace_pair(e, f) for e, f in zip(positive, negative))
                gram = matrix(QQ, [[trace_pair(h, k) for k in H] for h in H])
                toral = tuple(gram.solve_right(vector(QQ, [trace_pair(g, h) for h in H])))
                self.assertEqual(api.coordinates(t, u), expected, kind)
                self.assertEqual(api.toral_coordinates(t), toral, kind)
                self.assertTrue(all(trace_pair(g, e) == 0 for e in positive))

    def test_exceptional_minimal_conjugation(self):
        for kind, weight in [(('G', 2), (1, 0)), (('F', 4), (0, 0, 0, 1)),
                             (('E', 6), (1, 0, 0, 0, 0, 0)), (('E', 7), (0, 0, 0, 0, 0, 0, 1))]:
            api = root_borel_map(kind, weight)
            t = tuple(map(QQ, [2, 3, 5, 7, 11, 13, 17][:api.rank]))
            T, roots, _, H = adjoint_data(api, t)
            # Build v t v^-1, collect t^-1 v t v^-1 into ordered root coordinates,
            # and compare with Ad(v)theta(t); no fixed-vector recurrence here.
            parameters = [QQ.zero()]*api.dimension
            i, j = next(iter(api.brackets))
            parameters[i], parameters[j], parameters[-1] = QQ(2), QQ(-1), QQ(1)/2
            v = product_unipotent(roots, parameters)
            remaining = T.inverse()*v*T*v.inverse()
            collected = []
            for e in roots:
                # At this height remaining-I has u_i ad(E_i), but an arbitrary
                # matrix entry can contain higher weights. Use a Cartan column.
                keys, _, _ = adjoint_model((kind[0], kind[1]))
                R = RootSystem(list(kind)).root_lattice()
                root_keys = set(R.roots())
                cartan_columns = [k for k, key in enumerate(keys) if key not in root_keys]
                entry = next(ij for ij, c in e.dict().items() if c and ij[1] in cartan_columns)
                a = (remaining[entry]-int(entry[0] == entry[1]))/e[entry]
                collected.append(a)
                if a:
                    remaining = exponential_nilpotent(-a*e)*remaining
            self.assertEqual(remaining, identity_matrix(QQ, remaining.nrows()))
            h = sum((c*m for c, m in zip(api.toral_coordinates(t), H)), matrix(QQ, T.nrows(), sparse=True))
            expected = v*h*v.inverse()
            actual = h+sum(c*e for c, e in zip(api.coordinates(t, collected), roots))
            self.assertEqual(actual, expected, kind)

    def test_symbolic_jets_and_api_connection(self):
        for kind, weight in [(('A', 3), (1, 0, 0)), (('C', 3), (1, 0, 0)), (('G', 2), (1, 0))]:
            source = SpringerCharacterMap(kind, [weight, weight])
            api = source.borel_reconstruction()
            fs = api.symbolic_coordinates()
            for t in ((QQ.one(),)*api.rank, tuple(map(QQ, [2, 3, 5][:api.rank]))):
                u = tuple(map(QQ, range(1, api.dimension+1)))
                evaluate = api.symbolic_ring.hom(t+u, QQ)
                self.assertEqual(api.coordinates(t, u),
                    tuple(evaluate(f.numerator())/evaluate(f.denominator()) for f in fs))
                factors = tuple(f.derivative(v) for f, v in zip(fs, api.unipotent_parameters))
                self.assertEqual(api.root_factors(t),
                    tuple(evaluate(f.numerator())/evaluate(f.denominator()) for f in factors))
                self.assertEqual(api.root_factors(t), tuple(f(*t) for f in api.root_factors()))
                self.assertEqual(api.vertical_determinant(t), prod(api.root_factors(t)))
            self.assertEqual(api.root_factors((1,)*api.rank), (1,)*api.dimension)
            self.assertIs(api.springer_map, source)

    def test_root_strings_and_mixed_characters(self):
        for kind, weight in [(('G', 2), (1, 0)), (('F', 4), (0, 0, 0, 1)),
                             (('E', 6), (1, 0, 0, 0, 0, 0)), (('E', 7), (0, 0, 0, 0, 0, 0, 1)),
                             (('E', 8), (0, 0, 0, 0, 0, 0, 0, 1))]:
            api = root_borel_map(kind, weight)
            roots = set(api.positive_roots) | {tuple(-c for c in r) for r in api.positive_roots}
            for (i, j), c in api.brackets.items():
                a, b = api.positive_roots[i], api.positive_roots[j]
                p = 0
                while tuple(y-(p+1)*x for x, y in zip(a, b)) in roots:
                    p += 1
                self.assertEqual(abs(c), p+1)
        first = root_borel_map(('G', 2), (1, 0))
        second = root_borel_map(('G', 2), (0, 1))
        mixed = SpringerCharacterMap(('G', 2), [(1, 0), (0, 1)]).borel_reconstruction()
        c = second.springer_map.trace_gram[0, 0]/first.springer_map.trace_gram[0, 0]
        for t in ((1, 1), (1, -1), (2, 3)):
            u = tuple(range(1, 7))
            expected = tuple((a+c*b)/(1+c) for a, b in zip(first.coordinates(t, u), second.coordinates(t, u)))
            self.assertEqual(mixed.coordinates(t, u), expected)

    def test_validation_and_budget(self):
        api = root_borel_map(('F', 4), (0, 0, 0, 1))
        with self.assertRaises(ValueError):
            api.symbolic_coordinates()
        for t, u in [((0, 1, 1, 1), [0]*24), ((1, 1, 1), [0]*24), ((1,)*4, [0]*23)]:
            with self.assertRaises(ValueError):
                api.coordinates(t, u)
        with self.assertRaises(ValueError):
            root_borel_map(api.springer_map, (1, 0, 0, 0))
        with self.assertRaises(ValueError):
            root_borel_map(('A', 2))


if __name__ == '__main__':
    unittest.main()
