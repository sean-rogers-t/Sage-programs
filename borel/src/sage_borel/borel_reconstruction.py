"""Recover Borel Springer coordinates from a torus map and positive brackets.

The low-level recurrence is root-data based. The ready-to-use coordinate
backends are SL3 and Sp4, in the conventions of borel_coordinates.py.
No matrices of the input representation are constructed by this module.
"""
from sage.all import QQ, matrix, identity_matrix, vector
from .borel_coordinates import BorelSpringerMap


def _reduced(f):
    """Cancel polynomial factors and normalize the denominator's leading unit."""
    f.reduce()
    scale = f.denominator().lc()
    return f.parent()(f.numerator()/scale)/f.parent()(f.denominator()/scale)


def reconstruct_root_coordinates(positive_roots, brackets, root_characters,
                                 simple_root_values, parameters):
    """Solve Ad(tu)(h+n)=h+n using exp(ad E) and height-triangular substitution.

    Roots are nonnegative simple-root coefficient tuples, in height order.
    ``brackets[(i,j)] = c`` for i<j means [E_i,E_j]=c E_(alpha_i+alpha_j).
    Omitted brackets vanish. ``simple_root_values`` are alpha_i(h), while
    ``root_characters`` are e^alpha(t). Parameters and values must coerce
    to the same characteristic-zero rational function field. The caller supplies consistent
    root/bracket data; this routine checks shape, not a root-datum axiom set.
    The result is rational: regular extension is checked by the high-level API.
    """
    roots = tuple(tuple(QQ(c) for c in r) for r in positive_roots)
    n, rank = len(roots), len(simple_root_values)
    if (not n or len(parameters) != n or len(root_characters) != n
            or len(set(roots)) != n or any(len(r) != rank for r in roots)
            or any(any(c < 0 or c.denominator() != 1 for c in r) or sum(r) == 0 for r in roots)
            or list(map(sum, roots)) != sorted(map(sum, roots))):
        raise ValueError("supply distinct positive integral roots in height order and matching data")
    field = parameters[0].parent()
    if (field.characteristic() != 0 or not field.is_field()
            or not hasattr(field(1).denominator(), 'lc')):
        raise ValueError("parameters must belong to a characteristic-zero rational function field")
    index = {root: i for i, root in enumerate(roots)}
    constants = {}
    for (i, j), c in brackets.items():
        if not (0 <= i < j < n):
            raise ValueError("bracket keys must have 0 <= i < j < number of roots")
        root = tuple(a+b for a, b in zip(roots[i], roots[j]))
        if root not in index:
            raise ValueError("a supplied nonzero bracket must have a root as its sum")
        constants[i, j] = (index[root], field(c))
        constants[j, i] = (index[root], -field(c))
    # Basis (h,E_1,...,E_N). Treat h as one formal vector with its given
    # root evaluations. This is only an N+1 dimensional calculation.
    action = identity_matrix(field, n+1)
    for i, parameter in enumerate(parameters):
        ad = matrix(field, n+1)
        ad[i+1, 0] = -sum(c*v for c, v in zip(roots[i], simple_root_values))
        for j in range(n):
            if (i, j) in constants:
                k, c = constants[i, j]
                ad[k+1, j+1] = c
        exponential = identity_matrix(field, n+1)
        term = identity_matrix(field, n+1)
        for k in range(1, n+2):
            term = term*(parameter*ad)/k
            if term.is_zero():
                break
            exponential += term
        else:
            raise ValueError("root action is not nilpotent with the supplied data")
        action *= exponential
    result = []
    for i, character in enumerate(root_characters):
        r = field(character)
        if r == 1:
            raise ValueError("reconstruct symbolically before specializing root characters to one")
        lower = action[i+1, 0]+sum(action[i+1, j+1]*result[j] for j in range(i))
        result.append(_reduced(-r*lower/(r-1)))
    return tuple(result)


class ReconstructedBorelMap:
    """SL3/Sp4 Borel coordinates from a character or a supplied torus map.

    A torus_map callback receives (z1,z2), returning the first two diagonal
    entries of theta(t). A character callback receives the same torus
    coordinates. Custom characters are checked for Weyl invariance and
    trace nondegeneracy, not positivity in the representation ring.
    """

    def __init__(self, group, representation="defining", *, character=None, torus_map=None):
        self._source = BorelSpringerMap(group)
        for name in ("group", "root_order", "dimension", "ring", "field",
                     "torus_parameters", "unipotent_parameters"):
            setattr(self, name, getattr(self._source, name))
        if character is not None and torus_map is not None:
            raise ValueError("supply a character or a torus map, not both")
        if (character is not None or torus_map is not None) and representation != "defining":
            raise ValueError("do not combine a named representation with custom input")
        x, y = self.torus_parameters
        self.trace_form = None
        self.character = None
        if torus_map is not None:
            self.representation = "custom_torus_map"
            h = tuple(self.field(v) for v in torus_map(x, y))
            if len(h) != 2:
                raise ValueError("torus_map must return the first two diagonal entries")
        else:
            self.representation = representation if character is None else "custom_character"
            chi = self.field(character(x, y) if character is not None
                             else self._named_character(representation))
            self._check_laurent(chi, torus_only=True)
            substitutions = ({x: y, y: x},
                             {y: 1/(x*y)} if group == "SL3" else {y: 1/y})
            if any(chi.subs(subs) != chi for subs in substitutions):
                raise ValueError("character must be Weyl invariant")
            self.character = chi
            gradient = vector(self.field, [z*chi.derivative(z) for z in (x, y)])
            self.trace_form = matrix(QQ, [[(w*f.derivative(w)).subs({x: 1, y: 1})
                                           for w in (x, y)] for f in gradient])
            if self.trace_form.det() == 0:
                raise ValueError("character trace form is degenerate")
            h = tuple(_reduced(f) for f in self.trace_form.solve_right(gradient))
        for f in h:
            self._check_laurent(f, torus_only=True)
        h1, h2 = h
        if group == "SL3":
            self._toral = (h1, h2, -h1-h2)
            simple_values = (h1-h2, h1+2*h2)
            characters = (x/y, x*y*y, x*x*y)
            brackets = {(0, 1): 1}
        else:
            self._toral = (h1, h2, -h2, -h1)
            simple_values = (h1-h2, 2*h2)
            characters = (x/y, y*y, x*y, x*x)
            brackets = {(0, 1): 1, (0, 2): 2}
        self._coordinates = reconstruct_root_coordinates(
            self.root_order, brackets, characters, simple_values, self.unipotent_parameters)
        for f in self._coordinates:
            self._check_laurent(f)
        self.denominators_cancelled = True

    def _named_character(self, representation):
        x, y = self.torus_parameters
        weights = (x, y, 1/(x*y)) if self.group == "SL3" else (x, y, 1/y, 1/x)
        defining = sum(weights)
        if representation == "defining":
            return defining
        if representation == "dual":
            return sum(1/z for z in weights)
        if representation == "symmetric_square":
            return sum(weights[i]*weights[j] for i in range(len(weights)) for j in range(i, len(weights)))
        if representation == "adjoint":
            if self.group == "SL3":
                return defining*sum(1/z for z in weights)-1
            return self._named_character("symmetric_square")
        if representation == "second_fundamental" and self.group == "Sp4":
            return sum(weights[i]*weights[j] for i in range(4) for j in range(i+1, 4))-1
        raise ValueError("unsupported named representation for this group")

    def _check_laurent(self, f, torus_only=False):
        denominator = f.denominator()
        if (len(denominator.monomials()) != 1
                or any(v in denominator.variables() for v in self.ring.gens()[2:])):
            raise ValueError("reconstruction has uncancelled poles; input does not give a regular Borel map")
        if torus_only and any(v in f.numerator().variables() for v in self.ring.gens()[2:]):
            raise ValueError("torus input must not involve unipotent parameters")

    def _evaluate(self, formulas, torus, unipotent):
        if torus is None and unipotent is None:
            return formulas
        t, u = self._source._arguments(torus, unipotent)
        substitution = dict(zip(self.ring.gens(), t+u))
        return tuple(QQ(f.subs(substitution)) for f in formulas)

    def coordinates(self, torus=None, unipotent=None):
        return self._evaluate(self._coordinates, torus, unipotent)

    def toral_coordinates(self, torus=None):
        return self._toral if torus is None else self._evaluate(self._toral, torus, [0]*self.dimension)

    def group_matrix(self, torus=None, unipotent=None):
        return self._source.group_matrix(torus, unipotent)

    def lie_matrix(self, torus=None, unipotent=None):
        return BorelSpringerMap.lie_matrix(self, torus, unipotent)

    def _arguments(self, torus, unipotent):
        return self._source._arguments(torus, unipotent)

    def vertical_jacobian(self, torus=None, unipotent=None):
        return BorelSpringerMap.vertical_jacobian(self, torus, unipotent)

    def vertical_rank(self, torus, unipotent):
        return self.vertical_jacobian(torus, unipotent).rank()

    def root_factors(self):
        """Regular Laurent leading factors in the ordered positive roots."""
        return tuple(_reduced(f.derivative(v)) for f, v in zip(self._coordinates, self.unipotent_parameters))

    def inverse(self, torus, target):
        """Triangular inverse at rational t where all leading factors are nonzero.

        General singular fibers are not classified by this API.
        """
        t, target = self._source._point(torus, target)
        substitution = dict(zip(self.ring.gens()[:2], t))
        answer = []
        for f, q, v, desired in zip(self._coordinates, self.root_factors(),
                                   self.unipotent_parameters, target):
            leading = QQ(q.subs(substitution))
            if leading == 0:
                raise ValueError("a root factor vanishes; singular fibers require separate analysis")
            tail = f-q*v
            value = (desired-QQ(tail.subs(substitution)))/leading
            substitution[v] = value
            answer.append(value)
        return tuple(answer)


def reconstruct_borel_map(group, representation="defining", *, character=None, torus_map=None):
    """Reconstruct from a named character, custom character, or torus map."""
    return ReconstructedBorelMap(group, representation, character=character, torus_map=torus_map)
