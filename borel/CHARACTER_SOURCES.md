# Integral character sources: criterion, certificates and API

The Springer/Borel construction has two distinct steps: pullback by a
representation's normalized Springer map, followed by the Borel map on
Levi-invariant polynomials. This guide concerns the integral source of the
first step. The [G2 application](G2_CHARACTER_FORM.md) identifies that source
and its entire Borel image for both maximal parabolics.

The general criterion below is an elementary saturation argument, recorded
with proof. It is not presented as a new theorem in the literature. Its
applications to representations require their own certificates. The API
does not assert that every representation satisfies the criterion.

## Definition and hypotheses

Let G be simply connected, semisimple and of finite irreducible type,
T a maximal torus, X=X(T) its full weight lattice, and rho a nontrivial
finite-dimensional representation over C. Direct sums are allowed.
Use the normalized trace-projection map theta_rho:G -> g, so
d(theta_rho)_1 is the identity. Let t_Q^*=X tensor Q.

The torus source and its degree-one lattice are

\[
\mathcal C_\rho=\{f\in S_{\mathbb Q}(\mathfrak t^*):
                   \theta_\rho^*f\in\mathbb Z[X]\},\qquad
\Lambda_\rho=\mathcal C_\rho\cap\mathfrak t_{\mathbb Q}^*.
\]

No grading assumption is made about C_rho. Integral characters mean integral
virtual characters; positivity is not required.

Restricting the coefficients of f to Q loses no integral characters from
the complex polynomial image. For any fixed finite polynomial-degree bound,
pullback has an injective rational Laurent coefficient matrix. An integral
right-hand side therefore has rational coordinates in its image.

The trace projection and its representation-change identities are published
in [Kostant--Michor, Sections 2.1, 2.4 and Theorem 2.8](https://arxiv.org/pdf/math/0109066).
The polynomial-character subring and its transported Borel map are published
in [Kumar, Definition 4 and Theorem 5](https://kumar.math.unc.edu/papers/kumar76.pdf).
These are the underlying constructions, rather than consequences of our
integral calculations.

Theta has rational Laurent coefficients in these coordinates and is dominant:
its differential at 1 is the identity. Hence theta* is injective. Clearing
the finitely many linear-coordinate denominators shows that Lambda_rho has
full rank. It is a lattice, obtained by an exact integral coefficient
preimage, and is Weyl-stable.

In fact Lambda_rho is contained in X. For a linear input ell, differentiate
its integral Laurent pullback in each torus direction at 1. Those derivatives
are integers, and d(theta)_1=id identifies them with the fundamental-weight
coefficients of ell. This observation does not assert integrality of every
higher-degree input in C_rho.

## The linear lattice directly from weights

Write a weight mu as an integer column in the fundamental-weight basis,
with multiplicity m_mu, and put

\[
M_\rho=\sum_\mu m_\mu\mu\mu^t,\qquad
D_\rho=\sum_\mu\mathbb Z(m_\mu\mu)\subset\mathbb Z^r.
\]

Let D_rho^#={v in Q^r : v^t D_rho is integral}, the dual lattice for the
coordinate pairing. Then

\[
\boxed{\Lambda_\rho=M_\rho D_\rho^\#.}
\]

Indeed, the Laurent coefficient of theta*(ell) at e^mu is
m_mu ell^t M_rho^(-1) mu. Requiring all of these coefficients to be
integers is exactly M_rho^(-1)ell in D_rho^#, since M_rho is symmetric.
This is a direct consequence of trace projection, not an additional
all-degree theorem. Multiplicities must be included in D_rho.
For example, repeating a representation m times multiplies both M_rho
and D_rho by m and leaves Lambda_rho unchanged.

In the minimal G2, F4, E6 and E7 cases the weighted coefficient lattice
D_rho is the full coordinate lattice. The trace Gram matrix then directly
gives Lambda_rho as its column lattice. The all-primes certificates are
still needed to prove that this linear lattice generates the complete
polynomial-character intersection.

## Exact criterion

Choose a Z-basis ell_1,...,ell_r of Lambda_rho and put b_i=theta*(ell_i).
Then the following are equivalent:

1. C_rho = Sym_Z(Lambda_rho).
2. For every prime p, the reductions of b_1,...,b_r are algebraically
   independent in F_p[X].

**Proof of 2 implies 1.** Since the ell_i are a rational basis, write any
element of C_rho uniquely as F(ell)/D, with F in Z[u1,...,ur] and D a
least positive common denominator. If a prime p divides D, some coefficient
of F is not divisible by p. Integrality of F(b)/D makes F(b) zero modulo p.
This is a nonzero polynomial relation, contradicting 2. Therefore D=1.
The opposite inclusion is immediate.

**Proof of 1 implies 2.** Lift a nonzero modular relation to an integer
polynomial F with at least one coefficient not divisible by p. Every Laurent
coefficient of F(b) is divisible by p, so F(ell)/p lies in C_rho. It does
not lie in Sym_Z(Lambda_rho), by uniqueness of polynomial coordinates in
the rational basis ell. This contradicts 1.

The modular relation may be inhomogeneous. This is why testing only
homogeneous coefficient blocks would not prove the general criterion.
An obstruction at the torus level need not be invariant under a given Levi.

Once the criterion holds, the source is graded. For any standard Levi L,
equivariance gives the all-degree identification

\[
\{f\in S_{\mathbb Q}(\mathfrak t^*)^{W_L}:
        \theta^*f\in\mathbb Z[X]^{W_L}\}
=\operatorname{Sym}_{\mathbb Z}(\Lambda_\rho)^{W_L}.
\]

No division by |W_L| or integral averaging is used. The formula determines
an invariant ring, not necessarily a polynomial algebra over Z. Since
Lambda_rho is contained in X, its Borel image is integral: on G/B this
follows from integral weight classes, and pullback from G/P is primitive
in the integral Schubert basis.

## Finite sufficient certificate

Order Laurent exponents lexicographically. If the leading exponent vectors
of b_1,...,b_r are linearly independent over Q and their coefficients
survive modulo p, distinct monomials in the b_i have distinct leading
exponents modulo p. They are therefore algebraically independent there.

Only primes dividing a characteristic-zero leading coefficient can change
that leading term. The implementation checks those finitely many primes
separately. A unimodular integer row change on the b_i is allowed and does
not change the generated integral algebra.

Exponent rank is computed over Q, even for a positive-characteristic case.
For example t1^2,t2^2 are algebraically independent in characteristic 2.
A nonzero Jacobian modulo every prime is a stronger condition and is not
required. Finite Laurent supports and an additive lexicographic order
suffice; a well-order on all integer exponents is unnecessary.

A failed leading test is **inconclusive**. A search finding no modular
relation through degree d proves absence only through that degree. Neither
is evidence of the opposite all-degree conclusion.

For an abstract example of a genuine obstruction, b1=t and b2=t^2+2v
have a primitive linear coefficient lattice. Nevertheless
(b2-b1^2)/2=v is integral. Modulo 2 the relation is u2-u1^2.
This example illustrates the criterion; it is not claimed to arise from
a representation or to be a Springer-map counterexample.

## API

Version 0.10 uses Sage's standard Dynkin numbering and the simply connected
weight lattice. Affine, reducible and relabeled types are rejected here.
Nontrivial dominant integral highest weights specify genuine representations.

```python
from sage_borel import springer_character_map, g2_character_algebra

source = springer_character_map("G2", (1, 0))
assert source.dimension == 7
assert tuple(map(str, source.linear_generators)) == ("6*w1", "2*w2")
assert source.certify_polynomial_source().certified
assert not source.is_integral(source.ring.gen(0))
assert source.is_integral(source.linear_generators[0])

image = g2_character_algebra(1)
assert image.integral_degree(3).smith_factors == (144,)
```

The returned algebra's `presentation()` retains the existing rational-kernel
API. Complete integral kernels for these two cases are proved separately in
the G2 note; this release does not add a general integral-presentation solver.

`trace_gram` is the trace pairing in the simple-coroot basis.
`theta` contains the pullbacks of the fundamental weights as rational
Laurent polynomials in the fundamental-weight torus characters.
`linear_basis_matrix` has columns equal to the fundamental-weight
coordinates of a complete Z-basis of Lambda_rho.
`linear_generators` and `linear_pullbacks` are the corresponding tuples.

The linear lattice is computed by the row Hermite lattice of the full
Laurent coefficient matrix, accounting for its common denominator.
It is not obtained by independently clearing coordinate denominators.
`pullback(f)` and `is_integral(f)` evaluate actual finite Laurent coefficients.
They impose no Levi-invariance condition.

```python
from sage_borel import SpringerCharacterMap

twice = SpringerCharacterMap("G2", [(1, 0), (1, 0)])
assert twice.theta == source.theta
assert twice.linear_basis_matrix == source.linear_basis_matrix
assert twice.trace_gram == 2 * source.trace_gram

adjoint = springer_character_map("G2", (0, 1))
assert not adjoint.certify_polynomial_source().certified  # inconclusive
certificate = adjoint.certify_polynomial_source([[1, -3], [0, 1]])
assert certificate.certified
assert certificate.exceptional_primes == (3,)
assert adjoint.linear_basis_matrix == 4 * source.linear_basis_matrix
```

The repetition identity is an essential normalization check: the raw trace
index changes under repetition, but the map and its integral source do not.

`modular_obstruction(p, max_degree=d)` searches all monomials through d.
It returns `None` for no relation within that bound, or a `CharacterObstruction`
with the prime, a lifted relation in u1,...,ur, its polynomial witness in
fundamental weights, and its integral Laurent pullback. A returned witness
is outside Sym_Z(Lambda_rho); it need not be homogeneous or Levi-invariant.

```python
assert source.modular_obstruction(2, max_degree=3) is None
```

For other explicitly given integral Laurent inputs, use
`polynomial_source_certificate(pullbacks, basis_change=None)` or
`modular_character_relation(pullbacks, p, max_degree=3)`.
The latter returns only a lifted modular relation. An assertion that its
divided lift enlarges the polynomial algebra additionally requires
characteristic-zero algebraic independence.

The bounded relation search defaults to at most 1,000 monomials, 100,000
Laurent support exponents and 2,000,000 coefficient-matrix entries.
These can be raised explicitly with `max_monomials`, `max_terms` and
`max_matrix_entries`. Support bounds are checked after each expansion;
they are not a strict memory or time guarantee. Highest-weight character
construction can also be expensive for large representations.

## Verification and next cases

The [tests](tests/test_characters.py) include the direct-sum and duality
identities, rank-one normalization, the adjoint characteristic-three case,
a true abstract higher-degree obstruction, and an inseparable example
whose Jacobian test would give the wrong conclusion.
The [release verifier](tools/verify_character_release.py) runs the full
library suite and saves its result and source hashes in
[character-verification.json](character-verification.json).

The [G2 proof and independent audit](G2_CHARACTER_FORM.md) establish the
first complete application. The [F4 proof and API](F4_CHARACTER_FORM.md)
add Sym_Z(6 Q_long), all four maximal-parabolic image tables, and the
general `IntegralCharacterSource` degreewise intersection interface.
The [F4 representation comparison](F4_REPRESENTATION_COMPARISON.md)
certifies the adjoint source Sym_Z(18 Q_long) and proves degreewise
3^d scaling of every Levi source and image. Its general scaling argument
requires full-source certificates for both representations, not just
proportional linear lattices.
The [E6/E7 factories](MINUSCULE_CHARACTER_API.md) now package the proved
minimal sources and complete images. The [E6 adjoint comparison](E6_REPRESENTATION_COMPARISON.md)
proves Sym_Z(24P), which is not homothetic to the minimal Sym_Z(6Q).
It gives all-Levi sum/intersection formulas and complete incomparable
images on E6/P1.
A classification of representations
satisfying the criterion remains open work in this repository. No claim
that those questions are open in the literature, or explicitly posed by
Kostant--Michor, is made.
