"""Automatic finite-type Borel reconstruction in the character API's bases.

Torus variables are fundamental-weight characters; toral outputs are
simple-coroot coefficients. Positive root vectors use Sage's Chevalley
basis, with a public height order. Point arithmetic is exact over QQ.
"""
from functools import lru_cache
from types import MappingProxyType
from sage.all import (QQ, ZZ, LieAlgebra, RootSystem, PolynomialRing,
                      PowerSeriesRing, matrix, vector, prod, lcm)
from .characters import SpringerCharacterMap, springer_character_map, _integer
from .borel_reconstruction import reconstruct_root_coordinates


@lru_cache(maxsize=32)
def _root_data(kind):
    lattice = RootSystem(list(kind)).root_lattice()
    roots = sorted(lattice.positive_roots(), key=lambda r: (sum(r.to_vector()),
                                                          tuple(-c for c in r.to_vector())))
    coefficients = tuple(tuple(map(int, r.to_vector())) for r in roots)
    indices = {r: i for i, r in enumerate(roots)}
    algebra = LieAlgebra(QQ, cartan_type=list(kind))
    basis = algebra.basis()
    cartan = matrix(QQ, algebra.cartan_type().cartan_matrix())
    weights = tuple(tuple(map(int, cartan*vector(QQ, r))) for r in coefficients)
    brackets = {}
    for i, alpha in enumerate(roots):
        for j in range(i+1, len(roots)):
            beta = roots[j]
            if alpha+beta in indices:
                result = basis[alpha].bracket(basis[beta]).monomial_coefficients()
                if set(result) != {alpha+beta}:
                    raise ArithmeticError("Chevalley bracket did not have the expected root support")
                brackets[i, j] = QQ(result[alpha+beta])
    return coefficients, weights, tuple(brackets.items()), algebra


class RootBorelMap:
    """Borel subgroup reconstruction attached to a SpringerCharacterMap.

    Supports the same standard finite irreducible crystallographic types
    and simply connected character lattice as the source map. Evaluation
    is lazy and works on root hypertori using finite exact torus jets.
    Full symbolic expansion is separately budgeted, not done on creation.
    """

    def __init__(self, springer_map):
        if not isinstance(springer_map, SpringerCharacterMap):
            raise TypeError("supply a SpringerCharacterMap")
        self.springer_map = springer_map
        self.cartan_type = springer_map.cartan_type
        self.rank = springer_map.rank
        data = _root_data((self.cartan_type.type(), self.rank))
        self.positive_roots, self.root_weights, constants, self.lie_algebra = data
        self.brackets = MappingProxyType(dict(constants))
        self.dimension = len(self.positive_roots)
        self.max_height = max(map(sum, self.positive_roots))
        self.basis_convention = 'Sage Chevalley root vectors; fundamental-weight torus; simple-coroot toral output'
        self._cartan = matrix(QQ, self.cartan_type.cartan_matrix())
        rho = self._cartan.transpose().inverse()*vector(QQ, [1]*self.rank)
        self._height_scale = lcm(c.denominator() for c in rho)
        self.deformation_cocharacter = tuple(ZZ(self._height_scale*c) for c in rho)
        # Sparse arrows of ad(E_i) on (H_1,...,H_r,E_1,...,E_N).
        arrows = [dict() for _ in range(self.dimension)]
        for i, weight in enumerate(self.root_weights):
            for j, c in enumerate(weight):
                if c:
                    arrows[i][j] = (self.rank+i, QQ(-c))
        index = {r: i for i, r in enumerate(self.positive_roots)}
        for (i, j), c in self.brackets.items():
            k = index[tuple(a+b for a, b in zip(self.positive_roots[i], self.positive_roots[j]))]
            arrows[i][self.rank+j] = (self.rank+k, c)
            arrows[j][self.rank+i] = (self.rank+k, -c)
        self._arrows = tuple(arrows)
        self._symbolic = None
        self._root_factors = None

    def _point(self, torus, unipotent=None):
        t = tuple(QQ(v) for v in torus)
        if len(t) != self.rank or any(v == 0 for v in t):
            raise ValueError("supply one nonzero rational fundamental-weight coordinate per node")
        if unipotent is None:
            return t
        u = tuple(QQ(v) for v in unipotent)
        if len(u) != self.dimension:
            raise ValueError("unipotent coordinates must match positive_roots")
        return t, u

    def _action_columns(self, unipotent):
        """Ad(u) from sparse root-string operations, with exact rational input."""
        columns = [{j: QQ.one()} for j in range(self.rank+self.dimension)]
        for i in reversed(range(self.dimension)):
            parameter = unipotent[i]
            if parameter == 0:
                continue
            arrows = self._arrows[i]
            for j, column in enumerate(columns):
                total, term = dict(column), column
                for degree in range(1, self.max_height+2):
                    nxt = {}
                    for source, value in term.items():
                        if source in arrows:
                            target, c = arrows[source]
                            nxt[target] = nxt.get(target, 0)+parameter*c*value/degree
                    term = {k: v for k, v in nxt.items() if v}
                    if not term:
                        break
                    for k, v in term.items():
                        total[k] = total.get(k, 0)+v
                else:
                    raise ArithmeticError("positive root action failed to terminate")
                columns[j] = {k: v for k, v in total.items() if v}
        return columns

    def _toral_values(self, torus, series=False, precision=None):
        if series:
            precision = self.max_height+1 if precision is None else precision
            ring = PowerSeriesRing(QQ, 's', default_prec=precision)
            base = (ring(1)+ring.gen()).add_bigoh(precision)
            values = tuple(ring(t)*base**c for t, c in zip(torus, self.deformation_cocharacter))
        else:
            values = torus
        moments = [0]*self.rank
        for weight, multiplicity in self.springer_map.weight_multiplicities:
            monomial = multiplicity*prod(t**e for t, e in zip(values, weight))
            for j, c in enumerate(weight):
                if c:
                    moments[j] += c*monomial
        inverse = self.springer_map.trace_gram.inverse()
        h = tuple(sum(inverse[i, j]*moments[j] for j in range(self.rank)) for i in range(self.rank))
        return values, h

    def toral_coordinates(self, torus=None):
        """Simple-coroot coefficients, matching SpringerCharacterMap.theta."""
        return self.springer_map.theta if torus is None else self._toral_values(self._point(torus))[1]

    def coordinates(self, torus=None, unipotent=None):
        """Positive-root coefficients; no arguments requests budgeted symbolic output."""
        if torus is None and unipotent is None:
            return self.symbolic_coordinates()
        if torus is None or unipotent is None:
            raise ValueError("supply both torus and unipotent coordinates")
        t, u = self._point(torus, unipotent)
        root_values = tuple(prod(tj**e for tj, e in zip(t, w)) for w in self.root_weights)
        singular_torus = any(r == 1 for r in root_values)
        values, h = self._toral_values(t, series=singular_torus)
        if singular_torus:
            root_values = tuple(prod(tj**e for tj, e in zip(values, w)) for w in self.root_weights)
        columns = self._action_columns(u)
        result = []
        for i, r in enumerate(root_values):
            row = self.rank+i
            lower = sum(columns[j][row]*h[j] for j in range(self.rank) if row in columns[j])
            lower += sum(columns[self.rank+j][row]*result[j] for j in range(i)
                         if row in columns[self.rank+j])
            if lower == 0 and not hasattr(lower, 'precision_absolute'):
                lower = h[0].parent()(0) if singular_torus else QQ.zero()
            f = -r*lower/(r-1)
            if singular_torus and (f.valuation() < 0 or f.precision_absolute() < 1):
                raise ArithmeticError("insufficient jet precision or an uncancelled pole")
            result.append(f)
        return tuple(QQ(f[0]) for f in result) if singular_torus else tuple(map(QQ, result))

    def root_factors(self, torus=None):
        """Regular leading factors; no argument returns Laurent formulas.

        Formula construction uses the existing SpringerGeometry root
        division API and never expands the full unipotent coordinates.
        """
        if torus is None:
            if self._root_factors is None:
                from .geometry import SpringerGeometry
                ts = self.springer_map.torus_ring.gens()
                factors = {w: f*prod(t**e for t, e in zip(ts, w))
                           for w, f in SpringerGeometry(self.springer_map).root_factors()}
                self._root_factors = tuple(factors[w] for w in self.root_weights)
            return self._root_factors
        t = self._point(torus)
        roots = tuple(prod(tj**e for tj, e in zip(t, w)) for w in self.root_weights)
        use_jets = any(r == 1 for r in roots)
        _, h = self._toral_values(t, series=use_jets, precision=2)
        factors = []
        for alpha, weight, r in zip(self.positive_roots, self.root_weights, roots):
            value = sum(c*v for c, v in zip(weight, h))
            if use_jets:
                factors.append(QQ(value[1])/(self._height_scale*sum(alpha)) if r == 1
                               else r*QQ(value[0])/(r-1))
            else:
                factors.append(r*value/(r-1))
        return tuple(factors)

    def vertical_determinant(self, torus):
        """Fixed-t determinant along U, not the full group differential."""
        return prod(self.root_factors(torus))

    def symbolic_coordinates(self, *, max_roots=12, max_terms=50000):
        """Expand all coordinates, with an up-front positive-root budget.

        Defaults include A3, C3, D4 and G2. Larger types should normally
        use exact point evaluation. max_terms is a post-expansion output
        guard, not a bound on intermediate expression size or running time.
        """
        limit = _integer(max_roots, 'max_roots', 1)
        terms = _integer(max_terms, 'max_terms', 1)
        if self.dimension > limit:
            raise ValueError("symbolic root budget exceeded; use coordinates(torus, unipotent) or raise max_roots")
        if self._symbolic is None:
            names = tuple('t%d' % (i+1) for i in range(self.rank))
            names += tuple('u%d' % (i+1) for i in range(self.dimension))
            self.symbolic_ring = PolynomialRing(QQ, names=names)
            field = self.symbolic_ring.fraction_field()
            gs = tuple(map(field, self.symbolic_ring.gens()))
            self.torus_parameters, self.unipotent_parameters = gs[:self.rank], gs[self.rank:]
            h = tuple(sum(field(c)*prod(t**e for t, e in zip(self.torus_parameters, exponent))
                          for exponent, c in f.dict().items()) for f in self.springer_map.theta)
            roots = tuple(prod(t**e for t, e in zip(self.torus_parameters, w)) for w in self.root_weights)
            simple_values = tuple(self._cartan.transpose()*vector(field, h))
            fs = reconstruct_root_coordinates(self.positive_roots, self.brackets, roots,
                                               simple_values, self.unipotent_parameters)
            for f in fs:
                if (len(f.denominator().monomials()) != 1 or
                        any(v in f.denominator().variables() for v in self.symbolic_ring.gens()[self.rank:])):
                    raise ArithmeticError("reconstructed Springer coordinate has an uncancelled pole")
            self._symbolic = fs
        if any(len(f.numerator().monomials()) > terms for f in self._symbolic):
            raise ValueError("symbolic output term budget exceeded")
        return self._symbolic


def root_borel_map(source, highest_weight=None):
    """Attach to an existing character map, or construct from type and weight."""
    if isinstance(source, SpringerCharacterMap):
        if highest_weight is not None:
            raise ValueError("an existing character map already specifies its representation")
        return RootBorelMap(source)
    if highest_weight is None:
        raise ValueError("provide a highest weight or an existing SpringerCharacterMap")
    return RootBorelMap(springer_character_map(source, highest_weight))
