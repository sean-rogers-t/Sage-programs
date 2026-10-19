"""Exact products and integration of rational Schubert classes.

The characteristic map is a ring homomorphism: if ``beta(p_u)=sigma_u`` and
``beta(p_v)=sigma_v``, their cup product is ``beta(p_u*p_v)``.  Polynomial
lifts are therefore sufficient for arbitrary multiplication; the result does
not depend on which lifts were chosen, since positive-degree global Weyl
invariants generate the Borel ideal.  Divisor products instead use Chevalley's
formula, and complementary degrees use Poincare duality directly.
"""

from sage.all import QQ

from .expansion import Expansion


class ProductsMixin:
    """Cup products and the normalized top-degree integral on ``G/P``.

    All public arguments representing cohomology classes are immutable
    :class:`Expansion` objects in this map's flag-variety context.  Weyl
    words are accepted only by :meth:`schubert_product`.  Polynomial degree
    is half cohomological degree, and the top Schubert class integrates to 1.
    """

    def _product_expansion(self, expansion):
        if not isinstance(expansion, Expansion):
            raise TypeError("cohomology class must be a Schubert Expansion")
        self._expansion()._check_context(expansion)
        return expansion

    def _product_element(self, word):
        element = self._coerce_weyl(word)
        if not self._is_minimal(element):
            raise ValueError("Schubert index is not a minimal right-coset representative")
        return element

    def _basis_product(self, left, right):
        """Multiply validated basis indices, caching unordered pairs."""
        if not hasattr(self, "_product_basis_cache"):
            self._product_basis_cache = {}
        # Weyl elements do not have a useful portable ordering.  Reduced words
        # give a deterministic commutative key without enumerating the group.
        if (left.length(), self._word(left)) > (right.length(), self._word(right)):
            left, right = right, left
        key = (left, right)
        if key in self._product_basis_cache:
            return self._product_basis_cache[key]

        degree_left = int(left.length())
        degree_right = int(right.length())
        degree = degree_left + degree_right
        top = self.longest_representative()
        top_degree = int(top.length())
        if degree > top_degree:
            result = self._expansion()
        elif not degree_left:
            result = self._expansion({right: QQ.one()})
        elif degree == top_degree:
            # Integral(sigma_u*sigma_v) is delta_(v, dual(u)).  As the
            # top-degree group is one-dimensional this is the entire product.
            result = (self._expansion({top: QQ.one()})
                      if self.poincare_dual(left) == right else self._expansion())
        elif degree_left == 1:
            node, = self._word(left)
            result = self.chevalley(node, right)
        else:
            result = self.expand(self.schubert_polynomial(left)
                                 * self.schubert_polynomial(right))
        self._product_basis_cache[key] = result
        return result

    def schubert_product(self, left_word, right_word):
        r"""Return ``sigma_u cup sigma_v`` as a Schubert expansion.

        Each index is a Weyl element in this map or a Dynkin word, and must
        be a minimal right-coset representative for ``G/P``.  Polynomial lifts
        and the Borel ring homomorphism handle arbitrary classes.  Units,
        products above top degree, divisors and complementary degrees have
        direct shortcuts; repeated unordered pairs are cached.
        """
        left = self._product_element(left_word)
        right = self._product_element(right_word)
        return self._basis_product(left, right)

    def multiply(self, left, right):
        r"""Return the cup product of two Schubert ``Expansion`` objects.

        Bilinearity handles arbitrary rational, possibly inhomogeneous
        classes.  For homogeneous basis factors, ``deg(sigma_w)=length(w)``
        in polynomial degree.  Contributions above ``dim_C(G/P)`` vanish.
        Both inputs must belong to this map's flag-variety context.
        """
        left = self._product_expansion(left)
        right = self._product_expansion(right)
        if not left or not right:
            return self._expansion()
        identity = self.W.one()
        if len(left) == 1 and identity in left:
            return right * left[identity]
        if len(right) == 1 and identity in right:
            return left * right[identity]
        terms = {}
        for u, a in left.items():
            for v, b in right.items():
                coefficient = a * b
                for w, factor in self._basis_product(u, v).items():
                    terms[w] = terms.get(w, QQ.zero()) + coefficient * factor
                    if not terms[w]:
                        del terms[w]
        return self._expansion(terms)

    def integrate(self, expansion):
        r"""Return the integral over ``G/P``, with ``integral(sigma_top)=1``.

        Only the component in cohomological degree ``2*dim_C(G/P)``
        contributes; lower-degree terms of an inhomogeneous class are ignored.
        The result is an exact rational number.  On the point ``G/G`` the
        unit is already the top class and integrates to 1.
        """
        expansion = self._product_expansion(expansion)
        return QQ(expansion[self.longest_representative()])

    def _product_pairing(self, left, right):
        """Integrate two validated expansions by sparse Poincare duality."""
        if len(left) > len(right):
            left, right = right, left
        return sum((coefficient * right[self.poincare_dual(word)]
                    for word, coefficient in left.items()), QQ.zero())

    def intersection_number(self, *expansions):
        r"""Integrate the cup product of any number of Schubert expansions.

        The final multiplication is replaced by the pairing
        ``integral(sigma_u*sigma_v)=delta_(v, w_0*u*w_(0,P))``.  Thus no
        Schubert product or polynomial lift is needed for two factors.
        Inhomogeneous classes are allowed.  With no factors this returns the
        integral of the unit: 1 for ``G/G`` and 0 otherwise.
        """
        expansions = tuple(self._product_expansion(value) for value in expansions)
        if not expansions:
            return QQ.one() if not self.longest_representative().length() else QQ.zero()
        if any(not value for value in expansions):
            return QQ.zero()
        if len(expansions) == 1:
            return self.integrate(expansions[0])
        accumulated = expansions[0]
        for factor in expansions[1:-1]:
            accumulated = self.multiply(accumulated, factor)
            if not accumulated:
                return QQ.zero()
        return self._product_pairing(accumulated, expansions[-1])

    def _clear_product_caches(self):
        """Release cached arbitrary basis products; retain no input classes."""
        if hasattr(self, "_product_basis_cache"):
            self._product_basis_cache.clear()
