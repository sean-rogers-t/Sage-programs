# Conventions and correctness of the sparse calculation

Work over QQ, using fundamental-weight coordinates and Sage's root datum. The
chosen full-flag normalization is `beta(omega_i) = sigma_(s_i)`; on G/P this
divisor statement applies to omitted nodes. A word `[i1,...,id]`
means `s_i1 * ... * s_id`. For a minimal right-coset representative `w` of
length `d`, the coefficient of `sigma_w` in `beta(f)` is the constant term of

\[
\partial_{i_1}\cdots\partial_{i_d}f,
\qquad \partial_i f=\frac{f-s_i f}{\alpha_i}.
\]

Operators on the right act first. This convention agrees with the independent
discriminant-derived Schubert representatives in the tests. In fundamental-weight
coordinates a reflection changes only `omega_i`, replacing it by `omega_i-alpha_i`.
The code performs exact polynomial division and checks that the remainder is zero.

## Nonzero frontiers

For each homogeneous component of degree `d`, start with `{identity: f}`. If a
frontier contains `w -> partial_w(f)`, apply `partial_i` and prepend `s_i` when
`s_i*w` has greater length. A zero divided difference has no nonzero descendants,
so it is discarded immediately. At depth `d`, the surviving polynomials are
constant and give the coefficients.

Different reduced paths reaching the same Weyl element give the same operator
by the braid relations. The engine checks equality and keeps one value. Adding
the values would incorrectly multiply the coefficient by the number of paths.

For G/P, only minimal right-coset representatives are needed. Every suffix of a
reduced expression of a minimal right-coset representative is also minimal: a
right descent in the Levi for that suffix would also give a reduced expression
ending in that Levi generator for the whole word. Thus restricting each frontier
to minimal representatives cannot lose any desired coefficient. The input is
checked for Levi invariance before this pruning is used.

Mixed-degree inputs are split into homogeneous components. A component of
degree `d` contributes only to length `d`. Constants give the identity class.
An individual coefficient can instead be found by applying just one reversed
reduced word and taking its constant term.

## Divisor multiplication

The independent second route uses the ordinary Chevalley formula

\[
\sigma_{s_i}\sigma_w=
\sum_{\substack{\alpha>0\\\ell(ws_\alpha)=\ell(w)+1}}
\langle\omega_i,\alpha^\vee\rangle\,\sigma_{ws_\alpha}.
\]

Root lengths and the D fork are obtained from Sage's associated coroots.
Positive-root reflections are cached once per map; covers are cached once per
encountered class. Coefficients are combined after each product. On G/P the
allowed divisor indices are omitted nodes, for which this agrees with the
full-flag pullback of the product.

## Complementary elementary functions

Split the coordinate alphabet into head `X` and tail `Y`. In type A, the positive
degree coefficients of `E_(X union Y)(t)` vanish in the coinvariant quotient.
In B/C/D the same holds for the squared alphabet. Therefore

\[
E_Y(t)=E_X(t)^{-1},\qquad
e_i(Y)=(-1)^i h_i(X),
\]

with `X,Y` replaced by their squares in B/C/D. These are equalities of
characteristic-map images, not literal polynomial identities. The library uses
the smaller of the two monomial-count estimates and checks that the Levi
preserves the alphabet partition before applying the shortcut. Literal input
polynomials remain available through `tail_polynomial`.

The type D Euler generator is the product of the tail coordinates, with the
sign determined by the stated coordinate order. It is evaluated directly.
Poincare duality on minimal right-coset representatives uses `w0*w*w0P`.

## Scope of the evidence

The independent low-rank oracle explicitly implements orthonormal coordinate
swaps and sign changes. It starts with the product of positive roots divided by
`|W|`, obtains Schubert representatives by BGG operators, and multiplies those
polynomials. This checks the Chevalley route independently of its implementation.
The benchmarks separately verify the literal coordinate polynomial against every
class in the relevant degree, including zero coefficients, without the tail
shortcut. Exhaustive exceptional Weyl-group enumeration is not part of the test
suite; those tests verify divisors, quadratic invariants and API behavior.

Sage's underlying APIs are documented in its
[Weyl-group reference](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/weyl_group.html)
and [Coxeter-group reference](https://doc.sagemath.org/html/en/reference/categories/sage/categories/coxeter_groups.html).
