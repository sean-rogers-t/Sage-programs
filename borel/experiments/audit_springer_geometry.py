"""Independent G2 matrix-representation differential and generic-fiber audit.

Build the seven-dimensional representation from its sl2 strings, then the
adjoint matrices from commutators. No WeylCharacterRing or root-factor formula
is used in this independent matrix construction.
"""
from pathlib import Path
import json
from sage.all import (QQ, PolynomialRing, diagonal_matrix, identity_matrix,
                      matrix, vector, prod, Polyhedron, NumberField)
from sage_borel import SpringerGeometry, springer_character_map, g2_adjoint_quotient


def representation():
    weights = [(1, 0), (-1, 1), (2, -1), (0, 0), (-2, 1), (1, -1), (-1, 0)]
    roots = [(2, -1), (-3, 2)]
    H = [diagonal_matrix(QQ, [w[i] for w in weights]) for i in range(2)]
    E, F = [], []
    for i, alpha in enumerate(roots):
        e, f = matrix(QQ, 7), matrix(QQ, 7)
        for j, w in enumerate(weights):
            lower = tuple(w[k]-alpha[k] for k in range(2))
            if lower in weights:
                f[weights.index(lower), j] = 1
            upper = tuple(w[k]+alpha[k] for k in range(2))
            if upper in weights:
                top, depth = w, 0
                while tuple(top[k]+alpha[k] for k in range(2)) in weights:
                    top = tuple(top[k]+alpha[k] for k in range(2))
                    depth += 1
                e[weights.index(upper), j] = depth*(top[i]-depth+1)
        E.append(e)
        F.append(f)
    C = matrix(QQ, [[2, -3], [-1, 2]])
    comm = lambda a, b: a*b-b*a
    for i in range(2):
        for j in range(2):
            assert comm(H[i], E[j]) == C[i, j]*E[j]
            assert comm(H[i], F[j]) == -C[i, j]*F[j]
            assert comm(E[i], F[j]) == (H[i] if i == j else matrix(QQ, 7))
            if i != j:
                for generators in (E, F):
                    x = generators[j]
                    for _ in range(1-int(C[i, j])):
                        x = comm(generators[i], x)
                    assert not x
    basis = H+E+F
    lie_weights = [(0, 0), (0, 0)]+roots+[tuple(-v for v in a) for a in roots]
    columns = lambda xs: matrix(QQ, [x.list() for x in xs]).transpose()
    changed = True
    while changed:
        changed = False
        for i in range(len(basis)):
            for j in range(len(basis)):
                x = comm(basis[i], basis[j])
                if columns(basis+[x]).rank() > len(basis):
                    lie_weights.append(tuple(a+b for a, b in zip(lie_weights[i], lie_weights[j])))
                    basis.append(x)
                    changed = True
    assert len(basis) == 14
    coordinate_matrix = columns(basis)
    adjoint = [matrix(QQ, [coordinate_matrix.solve_right(vector(comm(x, y).list()))
                          for y in basis]).transpose() for x in basis]
    return basis, weights, adjoint, lie_weights


def trace_differential(basis, weights, torus_ring):
    n, d = len(basis), len(weights)
    ts = torus_ring.gens()
    T = [prod(t**w for t, w in zip(ts, weight)) for weight in weights]
    A = matrix(torus_ring, n, n, lambda i, j:
               sum(T[k]*sum(basis[j][k, l]*basis[i][l, k] for l in range(d))
                   for k in range(d)))
    B = matrix(QQ, n, n, lambda i, j: (basis[j]*basis[i]).trace())
    J = B.inverse()*A
    assert matrix(QQ, [[f(1, 1) for f in row] for row in J]) == identity_matrix(n)
    # The direct trace matrix has one Cartan block and twelve scalar blocks.
    assert all(J[i, j] == 0 for i in range(n) for j in range(n)
               if i != j and (i >= 2 or j >= 2))
    return J[:2, :2].det()*prod(J[i, i] for i in range(2, n)), J[:2, :2]


def main():
    minimal, weights, adjoint, adj_weights = representation()
    report = []
    for name, highest, basis, ws in (("minimal", (1, 0), minimal, weights),
                                      ("adjoint", (0, 1), adjoint, adj_weights)):
        s = springer_character_map("G2", highest)
        g = SpringerGeometry(s)
        independent, J = trace_differential(basis, ws, s.torus_ring)
        assert independent == g.cayley_determinant()
        assert J == g.torus_jacobian()
        if name == "minimal":
            point = (QQ(-1), QQ(2))
            vanished_root = (-3, 2)
        else:
            Z = PolynomialRing(QQ, 'z')
            z = Z.gen()
            field = NumberField(3*z*z+2*z+3, 'eta')
            point = (field(-1), field.gen())
            vanished_root = (2, -1)
        # Every root character differs from 1: the group element is regular
        # semisimple. Its Lie image has a zero root, hence is nonregular.
        assert all(prod(t**k for t, k in zip(point, root)) != 1
                   for root in adj_weights if any(root))
        h = tuple(f(*point) for f in s.theta)
        assert sum(k*v for k, v in zip(vanished_root, h)) == 0
        assert g.cayley_determinant()(*point) == 0
        quotient = g2_adjoint_quotient(representation=name)
        K = PolynomialRing(QQ, names=("a", "b")).fraction_field()
        a, b = K.gens()
        R = PolynomialRing(K, 2, names=("x", "y"))
        x, y = R.gens()
        equations = [R(f(x, y))-target for f, target in zip(quotient.character_polynomials, (a, b))]
        I = R.ideal(equations)
        degree = int(I.vector_space_dimension())
        print(name, "generic quotient degree", degree, flush=True)
        # Independent Newton polygon check: every exposed edge has exactly
        # two weights, and the two coefficient columns are independent.
        support = sorted(set(w for w in ws if any(w)))
        polygon = Polyhedron(vertices=support)
        edges = []
        for face in polygon.faces(1):
            vertices = [tuple(map(int, v.vector())) for v in face.vertices()]
            p, q = vertices
            on_edge = [w for w in support
                       if matrix(QQ, [vector(q)-vector(p), vector(w)-vector(p)]).det() == 0
                       and all(min(p[i], q[i]) <= w[i] <= max(p[i], q[i]) for i in range(2))]
            assert len(on_edge) == 2
            determinant = int(matrix(QQ, vertices).det())
            assert determinant
            edges.append(dict(vertices=vertices, determinant=determinant))
        normalized_area = int(2*polygon.volume())
        assert degree == normalized_area == (6 if name == "minimal" else 18)
        report.append(dict(representation=name, generic_degree=degree,
                           normalized_newton_area=normalized_area, edges=edges,
                           full_trace_matrix_size=14, coefficientwise_determinant_agreement=True,
                           torus_jacobian_agreement=True,
                           regular_to_nonregular_witness=dict(point=list(map(str, point)),
                               image=list(map(str, h)), vanishing_root=vanished_root),
                           groebner_basis=[str(f) for f in I.groebner_basis()]))
    Path(__file__).with_name("springer-geometry-audit.json").write_text(
        json.dumps(dict(result="passed", cases=report), indent=2)+"\n")


if __name__ == "__main__":
    main()
