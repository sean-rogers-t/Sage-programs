"""Small, explicitly normalized families of partial flag varieties.

All degrees in this module are polynomial degrees, half the cohomological
degree.  The two-step symplectic family retains a GL(b) factor and an Sp(2k)
factor while the complementary GL(a) block grows.
"""

from .core import BorelMap
from .operations import _integer, elementary_polynomial


class TwoStepSymplectic(BorelMap):
    r"""The isotropic flag variety ``U < V < C^(2n)``.

    Set ``a = n-b-k``.  The dimensions of ``U`` and ``V`` are respectively
    ``a`` and ``a+b=n-k``.  Its Levi is
    ``GL(a) x GL(b) x Sp(2k)`` and its omitted Dynkin nodes are ``(a,a+b)``.
    Require ``b,k >= 1`` and ``n >= b+k+1`` so both flag steps are genuine.

    ``middle_generators()`` returns ``c_j=e_j(x_(a+1),...,x_(a+b))``;
    ``tail_generators()`` returns ``e_i=e_i(x_(n-k+1)^2,...,x_n^2)``.
    Thus the respective polynomial degrees are ``j`` and ``2*i``.

    These are literal Levi-invariant input polynomials, not their Schubert
    expansions.  Use an algebra helper or ``expand`` to compute their images.
    """

    def __init__(self, n, b, k):
        n = _integer(n, "n", minimum=1)
        b = _integer(b, "b", minimum=1)
        k = _integer(k, "k", minimum=1)
        if n < b + k + 1:
            raise ValueError("n must be at least b+k+1 for two genuine flag steps")
        self._family_n = n
        self._family_b = b
        self._family_k = k
        self._family_a = n - b - k
        super().__init__(["C", n], omitted_nodes=(self.a, self.a + b))
        self._family_middle = None
        self._family_tail = None
        self._family_algebras = {}

    @property
    def n(self):
        """Symplectic rank, fixed at construction."""
        return self._family_n

    @property
    def a(self):
        """Dimension of the growing head block, fixed at construction."""
        return self._family_a

    @property
    def b(self):
        """Dimension of the middle block, fixed at construction."""
        return self._family_b

    @property
    def k(self):
        """Rank of the tail symplectic factor, fixed at construction."""
        return self._family_k

    def middle_generators(self):
        """Return an ordered dictionary of the middle GL(b) generators."""
        if self._family_middle is None:
            coordinates = self.classical_coordinates()[self.a:self.a + self.b]
            self._family_middle = {
                "c{}".format(j): elementary_polynomial(self.ring, coordinates, j)
                for j in range(1, self.b + 1)
            }
        return dict(self._family_middle)

    def tail_generators(self):
        """Return an ordered dictionary of the tail Sp(2k) generators."""
        if self._family_tail is None:
            coordinates = tuple(x ** 2 for x in self.classical_coordinates()[-self.k:])
            self._family_tail = {
                "e{}".format(i): elementary_polynomial(self.ring, coordinates, i)
                for i in range(1, self.k + 1)
            }
        return dict(self._family_tail)

    def joint_generators(self):
        """Return middle then tail generators, with names c1,...,e1,...."""
        generators = self.middle_generators()
        generators.update(self.tail_generators())
        return generators

    def _family_algebra(self, kind):
        # Local import keeps the family usable without introducing a circular
        # dependency between the generic image calculator and BorelMap.
        from .algebra import GeneratorAlgebra

        if kind not in self._family_algebras:
            generators = {
                "joint": self.joint_generators,
                "middle": self.middle_generators,
                "tail": self.tail_generators,
            }[kind]()
            self._family_algebras[kind] = GeneratorAlgebra(self, generators)
        return self._family_algebras[kind]

    def joint_algebra(self):
        """Analyze the joint middle/tail image, including mixed relations."""
        return self._family_algebra("joint")

    def middle_algebra(self):
        """Analyze just the middle GL(b) factor."""
        return self._family_algebra("middle")

    def tail_algebra(self):
        """Analyze just the tail Sp(2k) factor."""
        return self._family_algebra("tail")

    def clear_cache(self):
        """Clear Borel data and all previously created factor-algebra caches."""
        super().clear_cache()
        for algebra in getattr(self, "_family_algebras", {}).values():
            algebra.clear_cache()

    def predicted_relations(self):
        r"""Return ``h_(a+1)(Z),...,h_n(Z)`` in the joint generator ring.

        Here ``Z`` is the union of the squared middle coordinates and squared
        tail coordinates.  These arise by eliminating the complementary head
        block from the ordinary Borel presentation.  This construction uses
        symmetric-function recurrences, independently of Schubert expansions.
        Each relation ``h_r(Z)`` has polynomial degree ``2*r``.
        """
        algebra = self.joint_algebra()
        ring = algebra.domain
        formal = dict(zip(algebra.generator_names, ring.gens()))
        middle = [ring.one()] + [formal["c{}".format(j)]
                                  for j in range(1, self.b + 1)]
        # If C(u)=prod(1+x_j*u), then C(u)C(-u)=prod(1-x_j^2*u^2).
        middle_squared = [ring.one()]
        for r in range(1, self.b + 1):
            middle_squared.append((-1) ** r * sum(
                (-1) ** ell * middle[j] * middle[ell]
                for j in range(self.b + 1)
                for ell in (2 * r - j,)
                if 0 <= ell <= self.b))
        tail = [ring.one()] + [formal["e{}".format(i)]
                                for i in range(1, self.k + 1)]
        elementary = [ring(sum(middle_squared[j] * tail[r - j]
                               for j in range(self.b + 1)
                               if 0 <= r - j <= self.k))
                      for r in range(self.b + self.k + 1)]
        complete = [ring.one()]
        for r in range(1, self.n + 1):
            complete.append(ring(sum(
                (-1) ** (j + 1) * elementary[j] * complete[r - j]
                for j in range(1, min(r, self.b + self.k) + 1))))
        return tuple(complete[self.a + 1:self.n + 1])

    def predicted_hilbert_prefix(self, max_degree):
        r"""Coefficients through a bound of the presentation's Hilbert series.

        The series is ``prod_(r=a+1)^n(1-t^(2r))`` divided by
        ``prod_(j=1)^b(1-t^j) prod_(i=1)^k(1-t^(2i))``.  It is calculated
        without constructing a Weyl group basis or evaluating any polynomial.
        """
        max_degree = _integer(max_degree, "max_degree")
        coefficients = [1] + [0] * max_degree
        weights = list(range(1, self.b + 1)) + [2 * i for i in range(1, self.k + 1)]
        for weight in weights:
            for degree in range(weight, max_degree + 1):
                coefficients[degree] += coefficients[degree - weight]
        for r in range(self.a + 1, self.n + 1):
            weight = 2 * r
            coefficients = [value - (coefficients[degree - weight]
                                      if degree >= weight else 0)
                            for degree, value in enumerate(coefficients)]
        return tuple(coefficients)

    def family_parameters(self):
        """Return JSON-compatible dimensions, nodes and generator degrees."""
        return {
            "type": "C",
            "n": self.n,
            "a": self.a,
            "b": self.b,
            "k": self.k,
            "flag_dimensions": [self.a, self.a + self.b],
            "ambient_dimension": 2 * self.n,
            "omitted_nodes": list(self.omitted_nodes),
            "middle_generator_degrees": list(range(1, self.b + 1)),
            "tail_generator_degrees": [2 * i for i in range(1, self.k + 1)],
        }


def two_step_symplectic(n, b, k):
    """Construct :class:`TwoStepSymplectic`; see its parameter conventions."""
    return TwoStepSymplectic(n, b, k)
