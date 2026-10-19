# Mathematics of the Borel and kernel calculations

This guide explains the mathematical maps computed by `sage-borel`, the
conventions that determine signs and indices, and why its two exact evaluation
methods agree. [METHODS.md](METHODS.md) connects these ideas to runnable API
examples. [PRESENTATIONS.md](PRESENTATIONS.md) gives the full kernel algorithm
and the proof behind its completeness certificate.

The implementation supports finite irreducible crystallographic root systems.
All calculations here use rational coefficients. Polynomial degree `d` means
cohomological degree `2*d`. A polynomial in squared coordinate weights of
elementary index `i` therefore lies in cohomological degree `4*i`.

## 1. The geometric map and its polynomial source

Let `G` be a connected complex semisimple group, `T` a maximal torus, `B` a
Borel subgroup, and `P` a standard parabolic containing `B`. Let `L` be its
standard Levi subgroup, `W` the Weyl group of `G`, and `W_P` the Weyl group of
`L`. Work in the rational weight space and put

$$
S=\operatorname{Sym}_{\mathbb Q}(\mathfrak t^*_{\mathbb Q}),
\qquad I=S^W_+S,
$$

where `S^W_+` is the ideal of positive-degree invariants in `S^W`. The rational
Borel characteristic map gives the presentations

$$
H^*(G/B;\mathbb Q)\cong S/I,\qquad
H^*(G/P;\mathbb Q)\cong
S^{W_P}/(S^W_+S^{W_P}).
$$

Thus `b.expand(f)` computes the class `beta(f)` for a polynomial
`f` in `S^{W_P}`, expressed in a Schubert basis. It does not evaluate a
numerical function at a point. Weyl-invariant polynomials of positive degree
map to zero, and two input polynomials can represent the same class.

The variables `b.weights` are fundamental weights `omega_i`. The chosen
normalization is

$$
\beta(\omega_i)=\sigma_{s_i}\quad\text{on }G/B.
$$

On `G/P`, this formula applies to omitted Dynkin nodes, whose fundamental
weights are `W_P`-invariant. The degree `2*length(w)` classes `sigma_w` are
the Schubert basis denoted `epsilon_w` in the research papers. This fixes the
convention without relying on a printed name for a Weyl-group generator.

### Parabolics and their Schubert indices

`omitted_nodes` specifies the nodes removed from the Dynkin diagram to obtain
the semisimple Levi diagram. The **retained** nodes generate `W_P`. Hence a
single omitted node gives a maximal parabolic, several give a general partial
flag variety, every node omitted gives `G/B`, and none gives `G/G`.

The Schubert basis of `G/P` is indexed by minimal **right** coset
representatives `W^P` in `W/W_P`. Equivalently, `w` has no right descent at
any retained node. A word `[i1, ..., id]` means
`s_i1 * ... * s_id`, in that multiplication order. The code rejects inputs
that fail invariance under a retained simple reflection.

If `w0` and `w0P` are the longest elements of `W` and `W_P`, the top polynomial
degree is `N = length(w0*w0P) = dim_C(G/P)`. The Poincare-dual index of `w`
in this convention is `w0*w*w0P`. Every positive-degree class is nilpotent,
and any homogeneous class above degree `N` is zero.

### Why full-flag coordinates suffice for a partial flag

The projection `pi: G/B -> G/P` induces an injective rational pullback. In
the presentations above it comes from including `S^{W_P}` in `S`. Algebraically,
average over `W_P` to obtain the Reynolds operator

$$
\rho(f)=\frac{1}{|W_P|}\sum_{w\in W_P}w(f).
$$

It is an `S^W`-linear retraction onto `S^{W_P}`. Applying it to an expression
in `I` proves

$$
I\cap S^{W_P}=S^W_+S^{W_P}.
$$

Consequently, vanishing, linear dependence, and multiplication of our
Levi-invariant inputs can all be tested in `S/I`. This justifies the direct
coinvariant backend: passing to the full flag loses no relations among the
selected partial-flag classes.

## 2. How this relates to the Springer-induced maps

For an almost faithful irreducible representation `V(lambda)` (one with finite
kernel), restrict its Springer morphism `G -> g` to `Theta_lambda: T -> t`.
The papers use its pullback to
identify a polynomial part of the complexified
representation ring with invariant Cartan polynomials:

$$
\Theta_\lambda^*:
\operatorname{Sym}_{\mathbb C}(\mathfrak t^*)^{W_P}
\xrightarrow{\ \sim\ }
\operatorname{Rep}_{\lambda\text{-poly}}^{\mathbb C}(L).
$$

The corresponding map into flag cohomology is the transported Borel map

$$
\xi_\lambda^P=\beta_{\mathbb C}^P\circ(\Theta_\lambda^*)^{-1}.
$$

Here the Borel map is complexified, with target `H*(G/P;C)`.

If a representation-ring function has already been written as
`Theta_lambda^*(q)`, pass `q` in fundamental-weight coordinates to the
library. The generic engine computes the characteristic-map stage. It does
not take an arbitrary character, decide membership in this polynomial
subring, or automatically invert an arbitrary Springer pullback. The current
implementation requires rational polynomial coefficients, even though the
research construction above is stated over `C`.

The [exceptional catalogue](EXCEPTIONAL.md) supplies specific named inputs
from the audited papers. Explicit exceptional Springer formulas and their
separate checks are in [articles/calculations](../articles/calculations/).
Changing the representation can change its embedded polynomial representation
subalgebra; it does not create a different abstract Borel map after this
identification.

## 3. BGG evaluation in Schubert coordinates

For a simple root `alpha_i`, the divided-difference operator is

$$
\partial_i(f)=\frac{f-s_i(f)}{\alpha_i}.
$$

It lowers polynomial degree by one. A simple reflection fixes every
fundamental weight except `omega_i`, which it replaces by
`omega_i-alpha_i`. The numerator is divisible by `alpha_i`; the code performs
exact division and checks that the remainder is zero.

For a reduced expression `w=s_i1 ... s_id`, put
`partial_w=partial_i1 ... partial_id`. Operators on the **right act first**.
For homogeneous `f` of degree `d` and `w` of length `d`,

$$
[\sigma_w]\beta(f)=\partial_w(f),
$$

whose value is a constant. Equivalently, it is its constant term. Mixed-degree
inputs are split into homogeneous components: degree `d` contributes only
to indices of length `d`. Constants contribute to the unit class.

The operators satisfy braid relations, so `partial_w` does not depend on the
chosen reduced expression. The sparse engine begins with `identity -> f`.
From an entry `w -> partial_w(f)`, it applies `partial_i` and prepends `s_i`
when length increases. Two paths reaching the same Weyl element give the same
polynomial; their values are checked for equality and retained **once**.
Adding them would count reduced expressions instead of Schubert coefficients.

A zero divided difference has only zero descendants, so that branch can be
discarded. For a minimal right-coset representative every suffix of a reduced
expression is also minimal: a suffix with a Levi right descent would give a
reduced expression for the whole element ending in that Levi generator.
Therefore pruning nonminimal elements cannot remove a required coefficient.
This avoids enumerating the entire Weyl group, although a dense answer can
still have large support.

## 4. Divisor products and complementary symmetric functions

The independent divisor route uses the ordinary Chevalley formula

$$
\sigma_{s_i}\sigma_w=
\sum_{\substack{\alpha>0\\\ell(ws_\alpha)=\ell(w)+1}}
\langle\omega_i,\alpha^\vee\rangle\,\sigma_{ws_\alpha}.
$$

The coefficient pairs the weight with the **coroot**. In nonsimply laced
types, replacing the coroot by the root gives wrong factors. Sage supplies
associated coroots, including the type D fork conventions. The library caches
positive-root reflections and the length-one covers encountered in products.
For a partial flag, `i` must be omitted; terms with zero pairing vanish, and
the result agrees with pullback from `G/P` to `G/B`. Nonminimal full-flag cover
terms have zero coefficient in this omitted-node case, so the surviving
indices belong to `W^P`.

For a coordinate alphabet `X`, define

$$
E_X(t)=\prod_{x\in X}(1+xt),\qquad
H_X(t)=\prod_{x\in X}(1-xt)^{-1}.
$$

Their coefficients are elementary and complete symmetric polynomials. Split
the full alphabet into head `X` and tail `Y`. In type A, the positive-degree
coefficients of `E_(X union Y)` vanish in the coinvariant quotient. In B/C/D
the analogous assertion holds for the squared alphabet. Thus, in cohomology,

$$
E_Y(t)=E_X(t)^{-1}=H_X(-t),\qquad
\beta(e_i(Y))=(-1)^i\beta(h_i(X)),
$$

with `X,Y` replaced by `X^2,Y^2` in B/C/D. These are quotient identities,
not literal equalities of input polynomials. `tail_elementary` uses the side
with the smaller monomial-count estimate when the Levi preserves the split.
`tail_polynomial` always provides the literal input. Elementary and complete
polynomials themselves use coefficient recurrences, without listing all
subsets or multisets.

## 5. Classical coordinates and the coinvariant ideal

The public ring uses fundamental weights. The optimized backend substitutes
them into standard coordinate variables and reduces modulo the full
Weyl-invariant ideal. For type `A_r` put `m=r+1`; otherwise put `m=r`.

| Type | Fundamental weights in the coordinate realization |
| --- | --- |
| `A_r` | `omega_i = x1+...+xi - (i/m)*(x1+...+xm)` |
| `B_r` | `omega_i = x1+...+xi` for `i<r`; `omega_r = (x1+...+xr)/2` |
| `C_r` | `omega_i = x1+...+xi` for all `i` |
| `D_r` | Prefix sums for `i<=r-2`; `omega_(r-1)=(x1+...+x_(r-1)-xr)/2`; `omega_r=(x1+...+xr)/2` |

Type A uses `m` variables with the trace relation included. The trace
correction makes the substitution an actual weight-space realization before
reduction. The public `classical_coordinates()` are trace-zero coordinate
weights. The central `GL(m)` character has zero Borel image. The half factors
in B and D are essential, even for calculations over `QQ`.

The coordinate invariant ideals are

$$
\begin{array}{ll}
A_r:&(e_1(x),\ldots,e_m(x)),\\
B_r,C_r:&(e_1(x^2),\ldots,e_r(x^2)),\\
D_r:&(e_1(x^2),\ldots,e_{r-1}(x^2),x_1\cdots x_r).
\end{array}
$$

The type D Euler product cannot be omitted. Although
`e_r(x^2)=(x1*...*xr)^2`, killing that square alone does not kill the Euler
class. `euler_polynomial` returns the selected tail Euler product; its
invariance is checked when used as an input generator.

For A the backend uses the triangular Groebner basis
`g_i=h_i(x_i,...,x_m)`, `1<=i<=m`. For B/C it uses
`g_i=h_i(x_i^2,...,x_r^2)`, `1<=i<=r`. Elementary/complete identities give an
invertible triangular change from the invariant generators above. In degree
reverse lexicographic order the leading monomials are respectively `x_i^i`
and `x_i^(2*i)`. They are pairwise coprime, so Buchberger's product criterion
proves that these are already Groebner bases. Type D computes a Groebner
basis from its squared invariants and Euler product instead.

The normal form in this quotient is unique. Multiplying normal forms and
reducing again gives the normal form of the original product. The evaluator
caches reduced generator images and reduced products, avoiding raw high
powers in fundamental weights. It does not construct the entire ambient
monomial or Schubert basis. Standard-numbered A/B/C/D types use this backend
automatically for presentations; other supported types use BGG evaluation.
Either method supplies exact coordinates for the same images.

## 6. Chosen image algebras and complete kernels

Choose homogeneous Levi-invariant polynomials `f_1,...,f_s` of positive degrees
`d_i`. `GeneratorAlgebra` creates formal variables `y_i` with these weights and
the map

$$
\phi:\mathbb Q[y_1,\ldots,y_s]\longrightarrow H^*(G/P;\mathbb Q),
\qquad y_i\longmapsto\beta(f_i).
$$

Its quotient by `ker(phi)` is the **selected image algebra**. Unless the
chosen classes generate everything, it is smaller than the full cohomology
ring. A relation is a polynomial in the formal `y_i` whose substituted
cohomology class vanishes; it need not vanish before taking the Borel quotient.

`A.degree(d)` finds the vector space of relations of one weighted degree.
`A.presentation()` finds generators of the entire kernel ideal, a standard
monomial basis, and the finite Hilbert series. It discovers independent
monomial images, records dependencies at their border, and continues until
that border is exhausted. Finite-dimensional cohomology ensures termination;
no guessed degree cutoff is used. The basis size is the selected image
dimension, not necessarily `|W/W_P|`.

The completeness certificate checks vanishing relations, independence of
basis images, and closure of the monomial border, along with output
consistency. These finite checks imply that every kernel polynomial reduces
to zero; [PRESENTATIONS.md](PRESENTATIONS.md) gives the proof. `P.verify()`
rechecks with its selected evaluator. Explicitly choosing the other backend
gives a separate coordinate check when available. Certification is exact
rational computation, not a machine-checked formal proof or an automatic
assertion about integral cohomology.

For the rank-nine example `C9/P8`, choose `q=x9^2`. The result is

$$
\operatorname{im}(\phi)\cong\mathbb Q[q]/(q^9),\qquad
\deg(q)=2,\qquad \operatorname{Hilb}(t)=1+t^2+\cdots+t^{16}.
$$

This records nine independent classes, in cohomological degrees `0,4,...,32`.
The complete-symmetric relation comes from eliminating the eight head squared
generators. The computation supplies the same answer directly from invariant
normal forms; it does not assume this family formula. See the
[benchmark record](benchmarks/README.md) and the separate
[two-step derivation](experiments/TWO_STEP_DERIVATION.md).

## References and scope of verification

These implement established characteristic-map, divided-difference, and
finite-dimensional algebra methods; no novelty is claimed for those results.

* A. Borel, *Sur la cohomologie des espaces fibres principaux et des espaces
  homogenes de groupes de Lie compacts*, Annals of Mathematics 57 (1953),
  115-207, especially Theorem 26.1.
  [Original paper](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Borel-Sur.pdf).
* I. N. Bernstein, I. M. Gelfand and S. I. Gelfand, *Schubert cells and
  cohomology of the spaces G/P*, Russian Mathematical Surveys 28 (1973), 1-26.
  [Publisher and DOI](https://doi.org/10.1070/RM1973v028n03ABEH001557).
* S. Kumar, *Representation ring of Levi subgroups versus cohomology ring of
  flag varieties*, and S. Kumar--S. Rogers, Part II, for the Springer transport
  and its compatibility with flag projections.
  [Part I](https://kumar.math.unc.edu/papers/kumar76.pdf),
  [Part II](https://kumar.math.unc.edu/papers/kumar81.pdf).
* J.-C. Faugere, P. Gianni, D. Lazard and T. Mora, *Efficient Computation of
  Zero-dimensional Groebner Bases by Change of Ordering*, Journal of Symbolic
  Computation 16 (1993), 329-344, for the FGLM linear-algebra viewpoint.
  [Publisher and DOI](https://www.sciencedirect.com/science/article/pii/S0747717183710515).

The low-rank test oracle uses explicit coordinate swaps/sign changes and
discriminant-derived Schubert representatives rather than the Chevalley
implementation. Presentation tests compare complete kernels with independent
elimination calculations and cross-check the two evaluators. These checks
support the implementation; they do not exhaust all ranks or exceptional
cases. See [PRESENTATIONS.md](PRESENTATIONS.md) for the verification records.
