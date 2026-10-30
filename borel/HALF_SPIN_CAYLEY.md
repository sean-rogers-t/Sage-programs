# Individual half-spin Springer maps: a Pfaffian term and the zero fiber

Date: 22 October 2026. Status: uniform written deductions with exact
independent corroboration. External review and publication priority are
not established. This is a first geometric result outside the scalar
full-spin family, not a complete degree or finiteness classification.

Follow-up: the [boundary theorem](HALF_SPIN_BOUNDARY.md) now completes
the underlying global zero fiber and proves its dimension is
r(2r-1)-12 for r>=4. The lower-bound and open-boundary statements below
describe the scope of this preceding chart calculation.

## 1. Representation, cover and normalization

Work over C with G=Spin_(2r), r>=3. Let S_sigma, sigma=+1 or -1,
be the half-spin representation with weights

\[
 \tfrac12(\epsilon_1,\ldots,\epsilon_r),\qquad
 \epsilon_i\in\{1,-1\},\quad\prod_i\epsilon_i=\sigma.
\]

Thus dim S_sigma=2^(r-1). Define theta_sigma by the actual trace form:
B_sigma(theta_sigma(g),H)=tr(S_sigma(g)dS_sigma(H)). Identify the
Lie algebra with so_(2r) via the differential of the standard cover
pi:Spin_(2r) -> SO_(2r). In epsilon coordinates the Cartan trace Gram
matrix is 2^(r-3) I. Rank two is excluded because one half-spin of
Spin_4 is not infinitesimally faithful.

Neither half-spin descends through this standard pi: its kernel
element -1 acts as -I. If r is odd, half-spin is faithful. If r is
even, its kernel has order two, a different central subgroup from
ker pi. These cover choices must not be silently exchanged.

For example, use actual Spin torus coordinates z_i,b with
b^2=product z_i. The weight epsilon/2 has character
b product z_i^((epsilon_i-1)/2). At a central lift of -I in SO,
all z_i=-1, b^2=(-1)^r, and the half-spin scalar is sigma b.
Together with the two lifts of I this proves the kernel assertions.

## 2. Formula beyond scalar Cayley rescaling

Retain the full-spin scalar
c=chi_(S_+ direct-sum S_-)/2^r, and on c!=0 set

\[
 X=(\pi(g)-I)(\pi(g)+I)^{-1},\qquad c^2\det(I-X)=1.
\]

Choose the invariant Pfaffian polynomial p(X) whose value on the split
Cartan diag(x_1,...,x_r,-x_1,...,-x_r) is product x_i. In particular
p(X)^2=(-1)^r det X. Let J(X) be the polynomial Pfaffian cofactor map,
characterized on invertible X by

\[
 J(X)=p(X)X^{-1},\qquad XJ(X)=p(X)I.
\]

It has degree r-1 and extends polynomially to singular matrices.
Changing the sign of p interchanges the labels of the half-spins.

**Theorem 1.** On this entire Cayley chart,

\[
 \boxed{\theta_\sigma(g)=2c\bigl(X+\sigma J(X)\bigr).}     \tag{1}
\]

This is generally not a scalar multiple of X. The additional term
has a different direction and can cancel X without any scalar
character factor vanishing.

**Proof.** Temporarily put z_i=a_i^2 and b=product a_i, and define
A=product(a_i+a_i^(-1)), B=product(a_i-a_i^(-1)). Then
chi_sigma=(A+sigma B)/2, c=A/2^r, and
x_i=(a_i-a_i^(-1))/(a_i+a_i^(-1)). Direct logarithmic
differentiation of this character, divided by the Gram scalar,
gives

\[
 y_i=2c\left(x_i+\sigma\prod_{j\ne i}x_j\right).
\]

No division by x_i is needed. This is (1) on the Cartan. Both sides
are regular equivariant maps on the chart, so density of regular
semisimple elements proves (1) everywhere, including singular X. QED.

## 3. Exact zero fiber inside the Cayley chart

**Theorem 2.** For either half-spin, the scheme theta_sigma^(-1)(0)
intersect {c!=0} is smooth and is the disjoint union of:

- the two reduced points 1,-1 in ker pi;
- for each nontrivial (r-2)-nd root of unity u, two closed semisimple
  G-conjugacy classes, characterized by

\[
 X^2=uI,\qquad p(X)=-\sigma u,\qquad
 c=\pm(1-u)^{-r/2}.                                      \tag{2}
\]

The last notation means the two square roots of (1-u)^(-r).
Each noncentral class has dimension r(r-1), and the full differential
of theta_sigma has rank r^2 there. Thus there are 2(r-3) such classes
for each chosen half-spin. The half-spin character on each is
2^(r-1)c(1-u), which is nonzero.

**Proof of completeness.** By (1), a zero satisfies X+sigma J(X)=0.
If p=0, a skew matrix has rank at most 2r-2. Its Pfaffian cofactor
has rank two when rank X=2r-2 and is zero when rank X<=2r-4.
Equality of the ranks of X and J(X), with r>=3, therefore forces X=0.
The cover equation then gives c=+/-1 and hence g=+/-1.

If p!=0, multiply the equation by X to obtain X^2=-sigma p I.
Put u=-sigma p. Then u!=0, X is semisimple, and p^2=u^r while
p^2=u^2. Consequently u^(r-2)=1. The Cayley chart excludes u=1.
Conversely these conditions imply J=pX^(-1)=-sigma X, so they give
exactly the noncentral zeros.

The eigenvalues are +/-sqrt(u), each with multiplicity r. Such
matrices have two SO orbits distinguished by p; condition (2) selects
one. Their SO centralizer is GL_r, so the orbit dimension is
r(2r-1)-r^2=r(r-1). The two Spin lifts have distinct nonzero values
of the invariant c and cannot be conjugate. Each gives one G orbit.
This proves the asserted set and class count.

**Differential and scheme structure.** At a noncentral zero the
factor X+sigma J(X) vanishes, so its differential is simply multiplied
by 2c. Differentiating J=pX^(-1), using
dp(H)=p tr(X^(-1)H)/2, gives

\[
 D(X+\sigma J)_X(H)
 =H+XHX/u-\frac{\operatorname{tr}(XH)}{2u}X.               \tag{3}
\]

The tangent space splits into the conjugacy-orbit directions and the
centralizer gl_r. The differential is zero on the orbit directions.
On a centralizer block H=diag(A,-A^T), after choosing the +/-sqrt(u)
eigenspaces, it is represented by 2A-tr(A)I. Its eigenvalues are 2
on sl_r and 2-r on the scalar line; both are nonzero for r>=3.
Thus its rank is r^2 and its kernel is precisely the orbit tangent
space. The fiber has the same local dimension as that tangent space,
so its local ring is regular and the fiber scheme is smooth there.
At 1 and -1 the left-coordinate derivatives are I and -I, respectively.
Those two points are reduced. This proves the scheme assertion. QED.

For Spin_10, r=5, this gives **four 20-dimensional semisimple
conjugacy classes**, besides 1 and -1, in the chart zero fiber of
each faithful 16-dimensional half-spin representation. The full
differential has rank 25 on each class.

## 4. A large additional zero family on the boundary

Theorem 2 does not describe the whole zero fiber. In fact the boundary
contains a considerably larger family.

**Proposition.** For r>=4 the global zero fiber contains the Spin
lifts of the semisimple SO family with a four-dimensional -1
eigenspace, a four-dimensional +1 eigenspace, and r-4 remaining
distinct reciprocal eigenvalue pairs different from +/-1. Its
conjugacy-saturated family has dimension

\[
 \boxed{r(2r-1)-12.}                                     \tag{4}
\]

**Proof.** In the torus take z_1=z_2=-1 and z_3=z_4=1, with
the other z_i generic. The product A in the character formula has
two vanishing factors, as does B. Hence each first logarithmic
derivative of A and B vanishes. The trace Gram matrix is invertible,
so theta_sigma=0. Equivariance proves the statement on every conjugate.
The generic SO centralizer has dimension 6+6+(r-4)=r+8, and the
remaining eigenvalue pairs supply r-4 parameters modulo finite Weyl
identifications. The total dimension is r(2r-1)-(r+8)+(r-4), giving
(4). The finite Spin cover preserves dimension. QED.

In Spin_10 this is a **33-dimensional boundary family** mapping to
zero, compared with the 20-dimensional chart classes. Its closure
also lies in the zero fiber. This is a lower bound on the dimension
of the global zero fiber, not a classification of all boundary points
or proof that this family is an irreducible component of maximal size.
In particular the individual half-spin maps are not globally finite
for any r>=4.

## 5. Low ranks and what is already explained by known results

For r=3, Spin_6=SL_4 and the half-spins are defining and dual
defining. Their generalized Cayley degree is four by the known SL_4
case. Theorem 2 has no noncentral chart classes, as required; the
other two central zeros lie outside this chart.

For r=4, triality identifies each half-spin with the vector
representation composed with a group automorphism. Thus Spin_8 is
a control, not the first case beyond classical identifications.
The vector trace map on SO_8 is (a-a^(-1))/2. Over a generic Cartan
target each torus coordinate solves z_i^2-2y_i z_i-1=0, giving 16
points. Pullback along Spin_8 -> SO_8 gives degree 32, and triality
gives the same degree for each half-spin map on Spin_8. The order-two
cover is essential. This is a deduction from known triality and the
elementary vector formula, not a new degree claim.

Spin_10 is the first case in this sequence without those reductions.
Its generic half-spin degree and maximal finite target locus are
**not determined here**. Nor do we claim that the full-spin finite
locus transfers to either individual half-spin.

## 6. Independent checks and route forward

The [verifier](tools/verify_half_spin_cayley.py) compares sign-cube
weights with Sage's highest-weight character calculation, checks
normalized moments, and computes actual chiral Clifford traces.
Exact cyclotomic calculations check the zero classes and the complete
differential ranks, rather than merely substituting into (1).
Separate weight checks verify the boundary family, cover kernels,
and the low-rank degree controls. The
[report](half-spin-cayley-verification.json) records finite ranges,
methods and source hashes. These checks support the written proof;
the all-rank statement does not follow just from finite calculations.

The current advance is an explicit exceptional-fiber mechanism:
the Pfaffian cofactor cancels the linear term on specified conjugacy
classes, while a different character-vanishing mechanism creates a
larger boundary family. The previous scalar argument cannot be
transferred without handling these two terms.

The subsequent [boundary theorem](HALF_SPIN_BOUNDARY.md) determines the
remaining underlying zero fiber, including nonsemisimple points, and
proves that the 33-dimensional Spin_10 family is maximal in dimension.
Generic degree remains the next target. On the dense
chart c p!=0 the equations reduce to coupled quadratics

\[
 x_i^2-\frac{y_i}{2c}x_i+\sigma p=0,\qquad
 p=\prod_i x_i,\qquad c^2\prod_i(1-x_i^2)=1.
\]

They give an exact starting system, not a degree count: choices of
quadratic roots, cover conditions and excluded divisors still need
to be accounted for.

## 7. Sources and status

Known inputs are the half-spin weights and Clifford realization,
the normalized trace-projection construction, and the full-spin
Cayley chart in [Kostant--Michor](https://arxiv.org/pdf/math/0109066),
Sections 7.4--7.15. Individual chiral representations are different
from their full-spin convention in even dimension.
For the triality control, see
[Knus, Notes on Triality, Proposition 32](https://people.math.ethz.ch/~knus/papers/lens.pdf)
and [Bryant's spinor notes](https://sites.math.duke.edu/~bryant/Spinors.pdf).

The Pfaffian formula, chart zero-fiber classification and boundary
dimension bound are written deductions here from these standard
inputs. The selected source checks do not establish literature
priority. There is no claim of external review, a complete half-spin
degree/finite-locus classification, or a new runtime API in this step.
