"""Defining Springer maps on B=T U for SL3 and Sp4, in characteristic zero.

This concerns a Borel subgroup, not the cohomological BorelMap API.
Torus coordinates are epsilon characters, not fundamental-weight coordinates.
Symbolic formulas live over QQ; point/fiber methods accept rational points.
"""
from sage.all import QQ, PolynomialRing, matrix


class BorelSpringerMap:
    """Ordered positive-root coordinates for the defining representation.

    SL3 uses E12,E23,E13; Sp4 uses E12-E34,E23,E13+E24,E14
    in the symplectic basis (e1,e2,f2,f1). Elements are t times the
    ordered product of root elements. Outputs are *Lie* root coefficients.
    All ranks and fibers below hold with the torus point fixed.
    """

    def __init__(self, group):
        if group not in ("SL3", "Sp4"):
            raise ValueError("supported groups are 'SL3' and 'Sp4' (defining representation)")
        self.group = group
        self.representation = "defining"
        self.root_order = ((1, 0), (0, 1), (1, 1)) if group == "SL3" else (
            (1, 0), (0, 1), (1, 1), (2, 1))
        self.dimension = len(self.root_order)
        names = ("z1", "z2", "a", "b", "c") + (() if group == "SL3" else ("d",))
        self.ring = PolynomialRing(QQ, names=names)
        self.field = self.ring.fraction_field()
        generators = tuple(map(self.field, self.ring.gens()))
        self.torus_parameters = generators[:2]
        self.unipotent_parameters = generators[2:]

    def _point(self, torus, unipotent):
        t = tuple(QQ(x) for x in torus)
        u = tuple(QQ(x) for x in unipotent)
        if len(t) != 2 or any(x == 0 for x in t):
            raise ValueError("torus must contain two nonzero rational numbers")
        if len(u) != self.dimension:
            raise ValueError("wrong number of unipotent/root coordinates")
        return t, u

    def _arguments(self, torus, unipotent):
        if torus is None and unipotent is None:
            return self.torus_parameters, self.unipotent_parameters
        if torus is None or unipotent is None:
            raise ValueError("supply both torus and unipotent coordinates, or neither")
        return self._point(torus, unipotent)

    def coordinates(self, torus=None, unipotent=None):
        """Positive Lie-root coefficients of theta(tu); no args gives formulas."""
        (x, y), u = self._arguments(torus, unipotent)
        a, b, c = u[:3]
        if self.group == "SL3":
            return (x*a, y*b, x*(c+a*b))
        d = u[3]
        return ((x+1/y)*a/2, y*b, (x+y)*c/2+x*a*b/2, x*(d+a*c))

    def toral_coordinates(self, torus=None):
        """The full diagonal of theta(t), in the stated vector basis."""
        x, y = self.torus_parameters if torus is None else self._point(
            torus, [0]*self.dimension)[0]
        if self.group == "SL3":
            z = 1/(x*y)
            mean = (x+y+z)/3
            return (x-mean, y-mean, z-mean)
        h1, h2 = (x-1/x)/2, (y-1/y)/2
        return (h1, h2, -h2, -h1)

    def group_matrix(self, torus=None, unipotent=None):
        """The source matrix tu, useful for independent calculations."""
        (x, y), u = self._arguments(torus, unipotent)
        a, b, c = u[:3]
        if self.group == "SL3":
            return matrix([[x, x*a, x*(a*b+c)], [0, y, y*b], [0, 0, 1/(x*y)]])
        d = u[3]
        return matrix([[x, x*a, x*(a*b+c), x*(d+a*c)],
                       [0, y, y*b, y*c], [0, 0, 1/y, -a/y], [0, 0, 0, 1/x]])

    def lie_matrix(self, torus=None, unipotent=None):
        """The complete matrix theta(tu), including its toral part."""
        self._arguments(torus, unipotent)
        h = self.toral_coordinates(torus)
        v = self.coordinates(torus, unipotent)
        if self.group == "SL3":
            p, q, r = v
            return matrix([[h[0], p, r], [0, h[1], q], [0, 0, h[2]]])
        p, q, r, s = v
        return matrix([[h[0], p, r, s], [0, h[1], q, r],
                       [0, 0, h[2], -p], [0, 0, 0, h[3]]])

    def vertical_jacobian(self, torus=None, unipotent=None):
        """Derivative only along U, with t fixed (not the full B derivative)."""
        result = matrix(self.field, [[f.derivative(v) for v in self.unipotent_parameters]
                                     for f in self.coordinates()])
        if torus is None and unipotent is None:
            return result
        t, u = self._arguments(torus, unipotent)
        substitution = dict(zip(self.ring.gens(), t+u))
        return matrix(QQ, [[f.subs(substitution) for f in row] for row in result.rows()])

    def vertical_rank(self, torus, unipotent):
        """Exact vertical rank at a rational point."""
        return self.vertical_jacobian(torus, unipotent).rank()

    def fiber(self, torus, target):
        """Exact fiber over positive Lie-root coefficients, at fixed rational t.

        Return None for an empty fiber. Otherwise return dimension,
        parameters, and source coordinates parametrizing the fiber as affine
        space. The target's diagonal is implicitly theta(t). Free parameters
        are polynomial indeterminates over QQ. This also describes the fiber
        after base change to any characteristic-zero extension of QQ.
        """
        (x, y), v = self._point(torus, target)
        p, q, r = v[:3]
        if self.group == "SL3":
            a, b = p/x, q/y
            return dict(dimension=0, parameters=(), coordinates=(a, b, r/x-a*b))
        s = v[3]
        A, C = (x+1/y)/2, (x+y)/2
        b = q/y
        if A != 0 and C != 0:
            a = p/A
            c = (r-x*a*b/2)/C
            return dict(dimension=0, parameters=(), coordinates=(a, b, c, s/x-a*c))
        if A == 0 and p != 0:
            return None
        if A != 0:  # C=0: c is free, with a fixed by p.
            a = p/A
            if r != x*a*b/2:
                return None
            c = PolynomialRing(QQ, "v").gen()
            parameters = (c,)
        elif C != 0:  # a is free; c is solved by r.
            a = PolynomialRing(QQ, "v").gen()
            parameters = (a,)
            c = (r-x*a*b/2)/C
        elif b != 0:  # crossing: the mixed ab term still determines a.
            a = 2*r/(x*b)
            c = PolynomialRing(QQ, "v").gen()
            parameters = (c,)
        else:
            if r != 0:
                return None
            a, c = PolynomialRing(QQ, names=("v", "w")).gens()
            parameters = (a, c)
        return dict(dimension=len(parameters), parameters=parameters,
                    coordinates=(a, b, c, s/x-a*c))

    def inverse(self, torus, target):
        """Invert the fixed-t map only on torus strata where it is an isomorphism."""
        (x, y), v = self._point(torus, target)
        if self.group == "Sp4" and (x+1/y)*(x+y) == 0:
            raise ValueError("the fixed-t map is not an isomorphism; use fiber()")
        return self.fiber((x, y), v)["coordinates"]


def borel_springer_map(group):
    """Build Borel-subgroup coordinates for defining SL3 or Sp4."""
    return BorelSpringerMap(group)
