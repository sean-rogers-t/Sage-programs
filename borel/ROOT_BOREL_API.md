# Automatic finite-type Borel reconstruction

## Outcome and conventions

The [root reconstruction implementation](src/sage_borel/root_borel.py)
connects the existing `SpringerCharacterMap` directly to automatic
Chevalley root data. It supports the standard finite irreducible
crystallographic types accepted by that character API, including all
classical families and G2, F4, E6, E7, E8. The representation is an actual
nontrivial highest-weight module or a specified direct sum of such modules.
The trace form and its normalization are inherited unchanged.

The input lattice is the **full weight lattice of the simply connected
group**. In particular vector representations in B/D may have central
kernels; their actual torus coordinates are still the Spin weight-lattice
coordinates. Central-quotient descent is not inferred or substituted.
The coefficient field of the point API is Q; the proofs apply in
characteristic zero and the formulas define maps over C.

Three basis conventions are explicit:

* The torus coordinate t_i is e^(omega_i)(t), a fundamental-weight character.
* The toral output h_i is the coefficient of the simple coroot H_i.
* Positive-root outputs use Sage's Chevalley vectors E_alpha, with public
  root tuples in simple-root coordinates. They are sorted by height, then
  reverse lexicographic order, so the simple roots occur in node order.

If A is Sage's Cartan matrix and a is a root's simple-root column, its
fundamental-weight column is w=A a. Thus e^alpha(t)=product t_i^w_i and
alpha(theta(t))=sum w_i h_i. Transposing this conversion incorrectly
interchanges root and coroot conventions in non-simply-laced types.
The tests check the conversion against explicit Cartan actions.

Sage's nonsimple root vectors may differ in sign from those of the
[earlier SL3/Sp4 matrix formulas](BOREL_COORDINATES.md). The old APIs remain
unchanged. Exact basis conversions are checked in the new tests; their
raw root coefficients must not be compared without that conversion.

## Mathematics: sparse recurrence and regular evaluation

The [reconstruction theorem](BOREL_RECONSTRUCTION.md) supplies the
uniform formula. Let N be the number of positive roots and write

\[
 u=\prod_{i=1}^N\exp(u_iE_{\alpha_i}),\qquad
 \theta(tu)=h(t)+\sum_i F_iE_{\alpha_i}.
\]

The fixed-vector equation Ad(tu)theta(tu)=theta(tu), together with its
prescribed toral part h(t), gives the height-triangular recurrence

\[
 F_i=-\frac{r_i}{r_i-1}
 \left(v_i+\sum_{j<i}M_{ij}F_j\right),\qquad
 r_i=e^{\alpha_i}(t),
\]

where Ad(u)h=h+sum v_i E_alpha_i and Ad(u)E_alpha_j=sum M_ij E_alpha_i.
The new implementation constructs the positive-root brackets from
Sage's finite Lie algebra, not hand-written tables. It applies the
terminating exp(u_i ad E_alpha_i) operators as sparse root-string
operations on the Borel Lie algebra. It never constructs matrices of
the user's representation in order to evaluate the map. The representation
enters through the character API's torus coefficients and trace matrix.

The root factors are the same regular Laurent functions already provided
by the geometry API, with its left-translation factor included:

\[
 q_\alpha(t)=\frac{e^\alpha(t)\alpha(h(t))}{e^\alpha(t)-1},\qquad
 \det D_U\theta(tu)=\prod_{\alpha>0}q_\alpha(t).
\]

The determinant is for t held fixed. Its zero set tells where the
unipotent derivative is singular; it does not give the rank on that set
or the full B/G critical image and target branch/non-finite loci.

### A finite-jet evaluation lemma

The rational recurrence is directly evaluable when no r_i equals one.
On root hypertori, including t=1, apparent 0/0 expressions need regular
extension. Expanding the entire multivariable map first is unnecessary.

Let H be the maximum positive-root height. Choose a positive integer L
so that the cocharacter c=L A^(-transpose)(1,...,1)^transpose has integer
coordinates in the simple-coroot basis. Deform the torus by

\[
 t_j(s)=t_j(1+s)^{c_j}.
\]

Then

\[
 e^\alpha(t(s))=e^\alpha(t)(1+s)^{L\operatorname{ht}(\alpha)}.
\]

In characteristic zero each vanishing denominator r_alpha(s)-1 has
order exactly one. All other denominators are units. **Toral jets modulo
s^(H+1) suffice to recover the exact constant term of every F_alpha.**

Proof. The initial h(t(s)) is known modulo s^(H+1). A nonzero M_alpha,beta
requires ht(beta)<ht(alpha). Inductively F_beta is known to precision
at least H+1-ht(beta). Before the next division the precision is therefore
at least H+2-ht(alpha); division by a denominator of order at most one
leaves at least H+1-ht(alpha), which is at least one. Since the original
Springer map is regular, no negative Laurent coefficient survives. The
constant term is exactly its value at s=0. The program rejects a negative
valuation or insufficient remaining precision. Zero matrix arrows are
omitted so that unrelated truncated series do not cause artificial
precision loss.

This is exact truncated-series arithmetic, not numerical approximation
or an inference from a nearby regular point. At t=1 every root denominator
vanishes and the same argument still applies. A pure root subgroup is
sent to its Lie root line with derivative one. For general u the nonlinear
terms are retained; the method does not replace theta(u) by log(u).

For root factors alone, a first-order torus jet suffices. If r_alpha(t)=1,
the regular value is the derivative of alpha(h(t(s))) divided by
L ht(alpha). If r_alpha(t)!=1, use the ordinary quotient. This gives
exact leading factors on all rational torus points without full expansion.

## API examples

```python
from sage_borel import springer_character_map, root_borel_map, SpringerCharacterMap

source = springer_character_map(("E", 6), (1, 0, 0, 0, 0, 0))
b = source.borel_reconstruction()
assert b.dimension == 36  # number of positive roots, not dimension of V
assert b.springer_map is source

t = (2, 3, 5, 7, 11, 13)  # fundamental-weight characters
u = [0]*b.dimension
u[0] = 2
root_coefficients = b.coordinates(t, u)
coroot_coefficients = b.toral_coordinates(t)
assert len(root_coefficients) == 36
assert len(coroot_coefficients) == 6
assert b.root_factors((1,)*6) == (1,)*36

# Same factory with a highest weight supplied directly.
e7 = root_borel_map(("E", 7), (0, 0, 0, 0, 0, 0, 1))
assert len(e7.coordinates((1,)*7, [0]*63)) == 63

# An actual mixed representation; each highest weight is one summand.
mixed = SpringerCharacterMap(("G", 2), [(1, 0), (0, 1)])
g2 = root_borel_map(mixed)
assert g2.springer_map is mixed
assert g2.vertical_determinant((1, 1)) == 1
```

`positive_roots`, `root_weights`, `brackets`, `rank`, `dimension`,
`max_height`, `deformation_cocharacter` and `basis_convention` expose the
root conventions. The bracket table is read-only. `dimension` counts
positive roots; use `springer_map.dimension` for representation dimension.

`root_factors()` with no point returns regular Laurent formulas in the
source's torus ring. `root_factors(t)` evaluates them, including on root
hypertori, through the exact jet method. `vertical_determinant(t)` takes
their product without expanding all unipotent coordinates.

```python
from sage.all import QQ
g2 = root_borel_map(("G", 2), (1, 0))
formulas = g2.symbolic_coordinates()
assert len(formulas) == 6
assert len(g2.root_factors()) == 6
assert g2.coordinates((1, 1), (1, 2, 3, 4, 5, 6)) == (1, 2, 2, QQ(4)/3, QQ(1)/2, QQ(39)/2)
```

Full symbolic expansion is a separate operation. Its default budget is
12 positive roots, covering A3, C3, D4 and G2. E6/E7/E8 evaluation does
not require full symbolic expansion. Users may explicitly raise
`max_roots`; `max_terms` checks output numerator size after expansion
and is not an intermediate memory/time bound. The character API itself
enumerates representation weights, so very large highest weights can
still be costly. This is not a universal complexity guarantee.

## Independent checks and practical scope

[Tests](tests/test_root_borel.py) use several distinct constructions:

| Check | Scope | Independence and shared inputs |
| --- | --- | --- |
| Standard matrix projection | A3, B3, C3, D4 at three torus points each | Explicit classical root matrices; every bracket checked after aligning signs; direct g-tr(g)I/n or (g-g^(-1))/2 projection |
| Larger representation trace solve | SL4 Sym^2 at two points | Direct symmetric-square matrices and full trace Gram solve |
| Full adjoint trace projection | G2, F4, E6, E7, E8 at three points each | Full Lie-algebra matrices and trace pairing; no Borel recurrence; shares Sage Chevalley constants |
| Minimal-representation conjugation | G2/F4/E6/E7 | Build vtv^(-1), independently collect root coordinates and compare with Ad(v)theta(t); shares the audited torus input and Chevalley constants |
| Symbolic versus jet evaluation | A3, C3, G2 | Exact multivariable cancellation versus finite-series specialization; includes t=1 |
| Root-string magnitudes | All five exceptional types | Check absolute bracket constants using root-string lengths |
| Previous coordinate conventions | SL3 and Sp4 | Explicit basis/sign conversion to the earlier matrix formulas |
| Dense E8 identity-torus case | All 120 unipotent parameters nonzero | Exact jet evaluation compared with a separate 248-dimensional adjoint trace calculation |

The adjoint checks include the identity, a partially singular torus point
(1,2,1,...,1), and a regular point of distinct prime coordinates.
Classical checks also include a noncentral sign point. A G2 mixed-character
check verifies trace-form-weighted direct-sum normalization at regular
and singular torus points. Input and expansion-budget failures are tested.

The [verification report](root-borel-verification.json) records exact
counts, source hashes and a dense E8 timing. That one timing is a
benchmark on this machine, not a general runtime claim. These are
independent constructions within Sage. Exceptional checks share its
Chevalley structure constants and are not an external or separate-CAS
verification. The uniform recurrence and finite-jet proofs are separate
from the finite test cases.

## Published inputs, deductions, and next step

Published inputs are the generalized Cayley trace projection and
equivariance of Kostant--Michor,
[The generalized Cayley map](https://arxiv.org/abs/math/0109066), and the
weight/trace-matrix torus formula in Rogers,
[An Explicit Determination of the Springer Morphism](https://arxiv.org/abs/1701.01538).
Automatic structure constants use Sage's documented
[Chevalley-basis Lie algebra](https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/lie_algebras/classical_lie_algebra.html).
The height recurrence, finite-jet bound and API are repository deductions
and implementations. No novelty or exhaustive literature-priority claim
is made. The documentation was checked against the installed Sage 9.5
API and the official reference on 2026-10-22.

This completes the automatic finite-type adapter milestone. It does not
classify singular fibers, full-group critical images, branch loci, other
isogeny forms, affine groups or Kac--Moody groups. The next useful goal is
to compute vertical rank strata from root factors and the sparse mixed
terms, beginning with G2 minimal versus adjoint. The purpose is a reusable
criterion for when mixed terms recover directions after several leading
factors vanish, with exact fiber checks on a small exceptional benchmark.
Adding larger expanded formula tables alone would not settle that question.
