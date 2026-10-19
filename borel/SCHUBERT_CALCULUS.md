# Schubert multiplication and inverse Borel calculations

The library can multiply arbitrary rational Schubert classes on `G/P`, produce
polynomial representatives of those classes, and compute their intersection
numbers. These operations use the same Borel map and Schubert normalization as
the generator and kernel calculations. They apply to a general standard
parabolic, including partial flags with several omitted nodes.

Run the examples with SageMath, for example
`PYTHONPATH=src sage -python your_example.py` from `borel/`.
[METHODS.md](METHODS.md) describes the surrounding API;
[MATHEMATICS.md](MATHEMATICS.md) explains its Borel/Springer interpretation.

## The map and its chosen inverse section

Let `W_P` be the Levi Weyl group, and set

$$
S=\mathbb Q[\omega_1,\ldots,\omega_r],\qquad
I_W=S^W_+S.
$$

The rational Borel map is the surjective ring homomorphism

$$
\beta^P:S^{W_P}\longrightarrow H^*(G/P;\mathbb Q),\qquad
H^*(G/P;\mathbb Q)\cong S^{W_P}/(I_W\cap S^{W_P}).
$$

The indices `w` of the Schubert basis are minimal **right** coset
representatives in `W/W_P`. The class `sigma_w` has polynomial degree
`length(w)`, or cohomological degree `2*length(w)`. A word `[i1,...,id]` means
`s_i1*...*s_id`, with the rightmost divided-difference operator acting first.
An omitted Dynkin node `i` gives the divisor `beta(omega_i)=sigma_(s_i)`.

`b.schubert_polynomial(w)` chooses a homogeneous polynomial `p_w` in
`S^{W_P}` such that `beta^P(p_w)=sigma_w`. `b.lift(x)` extends these choices
linearly to an arbitrary Schubert expansion `x`. Consequently,

$$
\beta^P(\operatorname{lift}(x))=x.
$$

This is a chosen rational linear section of a surjection. Polynomial lifts are
not unique: their differences lie in the Borel kernel. Different methods can
return different polynomials representing the same class. The section need
not preserve multiplication, so `lift(x*y)` need not equal `lift(x)*lift(y)`
as a polynomial, although their Borel images agree. It also does not invert
an arbitrary character or an arbitrary Springer pullback.

The returned polynomials are **literally invariant** under every retained
simple reflection, before passing to the coinvariant quotient. They can be
passed directly to the same map's `expand` and `coefficient` methods.

## Public methods

| Method | Result |
| --- | --- |
| `b.schubert(word)` | The single class `sigma_w` as an immutable `Expansion`; `()` gives the unit |
| `b.schubert_polynomial(word, method="auto")` | A homogeneous Levi-invariant polynomial lifting that class |
| `b.lift(expansion, method="auto")` | A Levi-invariant polynomial lifting a possibly inhomogeneous rational class |
| `b.schubert_product(u,v)` | The Schubert expansion of `sigma_u*sigma_v`; inputs are Weyl elements or words |
| `b.multiply(x,y)` or `x*y` | The cup product of two `Expansion` objects |
| `x**n` | A nonnegative integral cup power; exponent zero gives the unit |
| `b.integrate(x)` | The coefficient of the top Schubert class, normalized to integrate to one |
| `b.intersection_number(x1,...,xm)` | The integral of the product of the supplied expansions |

Weyl indices must be minimal right-coset representatives. Expansion inputs
must have the same root datum, parabolic, and Weyl-group parent as the map.
Scalar multiplication of expansions remains available: use `2*x` or
`(QQ(1)/2)*x` for rational scaling. The methods
`lift`, `multiply`, and `integrate` accept expansions, not bare polynomial
inputs; first apply `expand` when starting with a polynomial.

## An arbitrary nondivisor product on Gr(2,5)

Take type `A4` with omitted node `2`. This is the Grassmannian of
two-dimensional subspaces of a five-dimensional space. In our convention its
partition indices fit in a two-by-three rectangle. The following classes are
created directly from their Dynkin words:

```python
from sage_borel import BorelMap

b = BorelMap("A4", omitted_nodes=2)
sigma2 = b.schubert([3, 2])
sigma31 = b.schubert([4, 3, 1, 2])
sigma22 = b.schubert([2, 3, 1, 2])

assert sigma2 * sigma2 == sigma31 + sigma22
assert b.multiply(sigma2, sigma2) == sigma31 + sigma22
assert b.schubert_product([3, 2], [3, 2]) == sigma31 + sigma22
assert b.schubert(()) * sigma2 == sigma2
```

This computes the Littlewood--Richardson identity

$$
\sigma_2^2=\sigma_{31}+\sigma_{22}.
$$

The missing `sigma_4` term is outside the rectangle and is zero in the
Grassmannian cohomology ring. Both factors have polynomial degree two, and
the result has polynomial degree four. This calculation goes beyond the
divisor-only Chevalley API.

For a connection with Cartan polynomials, write `h=omega_2` and
`c2=omega_1*(omega_2-omega_1)`. These are the elementary symmetric functions
of the first two coordinate weights. The Schur polynomials give

$$
\sigma_2=\beta(h^2-c_2),\qquad
\sigma_{31}=\beta(c_2(h^2-c_2)),\qquad
\sigma_{22}=\beta(c_2^2).
$$

The following block continues the example and checks the inverse operations:

```python
from sage.all import QQ

w1, w2, w3, w4 = b.weights
h, c2 = w2, w1 * (w2 - w1)
assert b.expand(h**2 - c2) == sigma2
assert b.expand(c2 * (h**2 - c2)) == sigma31
assert b.expand(c2**2) == sigma22

p_linear = b.schubert_polynomial([3, 2], method="linear")
p_bgg = b.schubert_polynomial([3, 2], method="bgg")
assert b.expand(p_linear) == sigma2
assert b.expand(p_bgg) == sigma2
assert b.expand(p_linear - p_bgg) == b.expand(0)

x = QQ(3)/7 * b.schubert(()) + 2 * sigma2 - sigma31
assert b.expand(b.lift(x)) == x
assert b.expand(b.lift(sigma2 * sigma2) - b.lift(sigma2)**2) == b.expand(0)
for node in b.retained_nodes:
    assert b.reflect(p_linear, node) == p_linear
    assert b.reflect(p_bgg, node) == p_bgg
```

## The degreewise linear inverse

The `linear` method works in the requested polynomial degree `d` without
constructing a representative of the global top class.

1. Enumerate the degree-`d` monomials in the fundamental weights. Their number
   is `binomial(r+d-1,d)`.
2. Stack the coefficient matrices of `1-s_i` for the retained nodes. The
   right kernel is exactly the homogeneous invariant space `(S^{W_P})_d`.
3. Apply `expand` to a basis `f_1,...,f_m` of that invariant space. This gives
   the matrix

   $$
   C_{w,j}=(\partial_w f_j)(0),\qquad w\in W^P,\quad\ell(w)=d.
   $$

4. Check that `C` has full row rank and solve `C*v=e_w` over `QQ`. The
   polynomial `sum_j v_j*f_j` lifts `sigma_w`.

Surjectivity of the homogeneous Borel map proves that these equations have
solutions. The implementation checks rank, checks the solved matrix equations
exactly, and checks literal Levi invariance of the resulting polynomials.
It caches the lifts of every class in that degree after a single matrix solve.

This construction depends on deterministic basis and linear-solve choices.
It does not assert a preferred geometric polynomial representative, nor a
minimal number of terms. It enumerates the Schubert basis **in that degree**,
so degreewise inverse calculation can cost more than evaluating a sparse
polynomial with `expand`.

## The normalized-discriminant BGG inverse

Let `w0` be the longest element of the full Weyl group and put

$$
\Delta=\prod_{\alpha\in\Phi^+}\alpha,\qquad
\partial_i(f)=\frac{f-s_i(f)}{\alpha_i}.
$$

With the library's indexing and operator order, the `bgg` method returns

$$
p_w=\partial_{w^{-1}w_0}\left(\frac{\Delta}{|W|}\right).
$$

The normalization is fixed by `partial_(w0)(Delta)=|W|`. The nil-Coxeter
relations and the degree-matching coefficient rule imply

$$
(\partial_u p_w)(0)=\delta_{u,w}\qquad\text{when }\ell(u)=\ell(w).
$$

Thus `beta(p_w)=sigma_w`, and `p_w` is homogeneous of degree `length(w)`.
The rightmost simple operator in a reduced word is applied first, as in the
existing `coefficient` method.

The same formula supplies partial-flag lifts without averaging. If `i` is
retained and `w` is minimal in its right coset, then `length(w*s_i)` is
`length(w)+1`. Consequently

$$
\ell(s_iw^{-1}w_0)=\ell(w^{-1}w_0)-1.
$$

The nil-Coxeter relation gives `partial_i*partial_(w^-1*w0)=0`, and therefore
`partial_i(p_w)=0`. Since `f-s_i(f)=alpha_i*partial_i(f)` in the polynomial
ring, this proves `s_i(p_w)=p_w` as an actual polynomial identity.

The code constructs the positive roots in the root lattice and converts them
to the map's fundamental-weight coordinates. It does not assume that roots
have equal lengths. Both inverse methods also handle relabeled Dynkin nodes
supported by `BorelMap`.

These are applications of the classical Borel and BGG framework, rather than
new existence theorems. See [BGG, *Schubert cells and cohomology of the spaces
G/P*](https://doi.org/10.1070/RM1973v028n03ABEH001557), and
[Billey, *Kostant polynomials and the cohomology ring for G/B*, Proposition 2.1
and Section 6](https://sites.math.washington.edu/~billey/papers/kostant.polys.pdf)
for divided differences, the top discriminant and cup-product calculations.
That paper's sign and indexing conventions require translation to the
normalization `beta(omega_i)=sigma_(s_i)` and rightmost-first operator order
used here. Its orbit-value matrix algorithm is distinct from the degreewise
invariant matrix construction implemented above.

## Why multiplication is independent of the lifts

If `beta(p)=x` and `beta(q)=y`, the ring homomorphism property gives

$$
xy=\beta(pq).
$$

Changing either representative changes `pq` by an element of the Borel ideal,
so the resulting Schubert expansion is independent of the inverse method.
The implementation handles general expansions by bilinearity, and caches
unordered pairs of basis indices. For each pair it uses these shortcuts:

- Products above `dim_C(G/P)` in polynomial degree are zero.
- Multiplication by the unit returns the other factor.
- Complementary-degree basis products are determined by Poincare duality.
- Multiplication by an omitted-node divisor uses the existing Chevalley rule.
- Other basis products multiply the chosen polynomial lifts and apply `expand`.

The library does not need to enumerate the whole Weyl group for the direct
shortcuts. General products can require degreewise basis generation or a large
discriminant through the inverse calculation. These product methods compute
ordinary rational cup products; equivariant and quantum parameters are absent.

## Integration and intersection numbers

Put `w_top=w0*w0P` and `w_dual=w0*w*w0P`. The normalization is

$$
\int_{G/P}\sigma_{w_{\rm top}}=1,\qquad
\int_{G/P}\sigma_u\sigma_v=\delta_{v,w_0uw_{0,P}}.
$$

`integrate(x)` extracts the top coefficient. Lower-degree components of an
inhomogeneous input make no contribution. For two expansions the intersection
method pairs their sparse coefficients directly using the dual-index formula;
it does not construct polynomial lifts. With more factors it multiplies all
but the last and then uses this pairing.

Continuing on `Gr(2,5)`, whose complex dimension is six:

```python
H = b.schubert([2])
point = b.schubert(b.longest_representative())
assert b.integrate(point) == 1
assert b.integrate(H) == 0
assert b.intersection_number(H, H, H, H, H, H) == 5
assert b.integrate(b.expand(h**6)) == 5

u = next(iter(sigma2))
dual_class = b.schubert(b.poincare_dual(u))
assert b.intersection_number(sigma2, dual_class) == 1
assert b.integrate(sigma2 * dual_class) == 1
```

The number five is the Plucker degree of `Gr(2,5)`. On the point `G/G` the
unit is already the top class and integrates to one. A zero-factor
`intersection_number()` means the integral of the unit, so it is zero for a
positive-dimensional flag variety and one for the point.

## Choosing a method and understanding costs

`method="linear"` and `method="bgg"` explicitly request the constructions
above. `method="auto"` uses the following size estimates, with `N=length(w0)`:

- Use the degreewise linear inverse if `binomial(r+d-1,d) <= 150`.
- Otherwise use BGG if `binomial(r+N-1,N) <= 20000`, an upper bound for the
  number of monomials in the full discriminant.
- Otherwise use the degreewise linear inverse.

The unit returns `1`, and a divisor returns its fundamental weight directly.
The heuristic estimates polynomial-space sizes; it does not guarantee a time
or memory bound. A high-degree inverse in a large root system may be expensive
by either method. Low degrees in exceptional types can remain small: degree
two in rank seven has only 28 weight monomials, without expanding a degree-63
E7 discriminant.

`b.clear_cache()` releases the inverse polynomials, degreewise inverse data,
normalized discriminant, and arbitrary basis-product cache, together with the
map's existing caches. Repeated calculations on the same map reuse those
objects until they are cleared.

All coefficients and linear equations use exact rational arithmetic.
Integral Schubert structure constants can emerge from rational polynomial
lifts; denominators in a lift do not themselves indicate an error. This API
does not certify an integral image lattice or identify torsion primes.
