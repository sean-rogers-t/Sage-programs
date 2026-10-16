"""Classical conventions and signed-permutation embeddings."""

from sage.all import CartanType, WeylGroup

from .core import BorelMap
from .operations import _integer


def classical_borel(family, n, k):
    """Construct the maximal parabolic used in the historical programs.

    In type A, ``n`` is the GL(n) dimension, so the Dynkin type is A_(n-1).
    In B/C/D it is the Dynkin rank.  In each case the omitted node is n-k.
    The supported historical ranges are A: 1 <= k < n; B/C: 0 <= k < n;
    D: 2 <= k < n.  Use BorelMap directly for spin or other parabolics.
    """
    if not isinstance(family, str) or family.upper() not in ("A", "B", "C", "D"):
        raise ValueError("family must be A, B, C or D")
    family = family.upper()
    n = _integer(n, "n", minimum=4 if family == "D" else 2)
    k = _integer(k, "k", minimum=2 if family == "D" else (1 if family == "A" else 0))
    if k >= n:
        raise ValueError("k must be smaller than n")
    return BorelMap([family, n - 1 if family == "A" else n], omitted_nodes=n - k)


def to_type_a(element):
    """Embed a B/C/D Weyl element into A_(2n-1) as a permutation.

    Signed coordinate order is (e1,...,en,-en,...,-e1).  Reduced words come
    from Sage's structural API; printed generator names are never parsed.
    The type D fork uses s_n -> r_n r_(n+1) r_(n-1) r_n.
    """
    if not hasattr(element, "parent") or not hasattr(element, "reduced_word"):
        raise TypeError("element must be a Sage Weyl-group element")
    try:
        kind = element.parent().cartan_type()
    except AttributeError as error:
        raise TypeError("element must have a Weyl-group parent") from error
    if not kind.is_finite() or not kind.is_irreducible() or kind.type() not in ("B", "C", "D"):
        raise ValueError("the signed embedding requires finite irreducible type B, C or D")
    n = int(kind.rank())
    if kind != CartanType([kind.type(), n]):
        raise ValueError("the signed embedding requires standard Sage Dynkin numbering")
    target = WeylGroup(["A", 2 * n - 1], prefix="r")
    generators = target.simple_reflections()
    result = target.one()
    for node in element.reduced_word():
        node = int(node)
        if node < n:
            indices = (node, 2 * n - node)
        elif kind.type() == "D":
            indices = (n, n + 1, n - 1, n)
        else:
            indices = (n,)
        for index in indices:
            result = result * generators[index]
    return result
