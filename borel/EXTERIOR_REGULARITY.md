# Regularity of the exterior-power Springer maps for SL_n

Date: 22 October 2026. Status: a written uniform deduction from the
published trace construction and weight formula, with independent finite
checks. No claim of literature priority or external review is made.
This concerns characteristic zero and complex regular elements, not
invertibility of the torus Jacobian or an integral Borel image.

## Statement and conventions

Let G=SL_n(C), n>=2, T its diagonal determinant-one torus, and
rho_k=exterior^k(C^n), 1<=k<=n-1. These representations have faithful
differential; their possible finite central kernels are allowed. Write
theta_k:G -> sl_n for the normalized trace Springer morphism:

\[
 \operatorname{tr}(d\rho_k(\theta_k(g))d\rho_k(X))
 =\operatorname{tr}(\rho_k(g)d\rho_k(X))\quad(X\in\mathfrak{sl}_n).
\]

Here regular means centralizer dimension n-1. On T this means that all
diagonal entries are pairwise distinct; on the diagonal Lie algebra it
has the same meaning, without a determinant constraint.

**Theorem.** The following are equivalent:

1. theta_k(T_reg) is contained in t_reg.
2. theta_k(G_reg) is contained in (sl_n)_reg.
3. k=1 or k=n-1.

For every interior case 2<=k<=n-2, a rational diagonal, regular
semisimple element of SL_n maps to a nonregular element. In particular,
the negative result is not limited to numerical examples or bounded n.

The standard positive case is already in Kostant--Michor Section 6.2.
The dual positive case is its elementary inverse-matrix analogue. The
interior obstruction and the resulting all-n equivalence are deductions
proved here; this audit has not established that they are unpublished.

## 1. Published inputs and the specialized weight formula

The trace projection is the construction of
[Kostant--Michor](https://arxiv.org/pdf/math/0109066), Section 2.1,
also used in the Kumar papers. The explicit general weight formula is
from [Rogers](https://arxiv.org/pdf/1701.01538v1), Theorem 1. Neither is
a new theorem of this note. The exterior-power character is the standard
elementary symmetric polynomial e_k.

Put z=(z_1,...,z_n), product_i z_i=1, and

\[
 I_{n,k}=\binom{n-2}{k-1}.
\]

For trace-zero diagonal X=diag(x_i), Y=diag(x'_i), counting subsets gives

\[
 \operatorname{tr}(d\rho_k(X)d\rho_k(Y))
 =\left(\binom{n-1}{k-1}-\binom{n-2}{k-2}\right)
       \sum_i x_i x'_i
 =I_{n,k}\operatorname{tr}(XY).
\]

The binomial with lower index -1 is zero. By invariance, this is the
trace form on all sl_n. Set v_i=z_i e_{k-1}(z with z_i omitted).
The representation eigenvalues are product_{i in S} z_i for |S|=k,
and their derivative against X gives sum_i v_i x_i. Since
sum_i v_i=k e_k(z), trace projection gives

\[
 \boxed{\quad y_i=\frac{z_i e_{k-1}(z_{\widehat i})
                -\frac{k}{n}e_k(z)}{I_{n,k}},\qquad
 \theta_k(\operatorname{diag}z)=\operatorname{diag}y.\quad}
\]

Here a hat means omission, and e_0=1. Separating subsets containing j
from those not containing j gives the crucial cancellation

\[
 \boxed{\quad y_i-y_j=
 \frac{z_i-z_j}{I_{n,k}}
 e_{k-1}(z_{\widehat i\widehat j}).\quad}                 \tag{1}
\]

Consequently, on T_reg, the failure locus is exactly

\[
 \bigcup_{i<j}\{z:e_{k-1}(z_{\widehat i\widehat j})=0\}.
\]

This states image regularity. It does not state that the torus derivative
is singular at exactly those points.

## 2. A uniform rational obstruction for every interior power

Assume 2<=k<=n-2. Put m=n-2 and r=k-1, so 1<=r<=m-1.

There is first a short simultaneous obstruction for all interior powers
at a fixed n. Let w consist of all m-th roots of unity, and set

\[
 z=(2,(-1)^{m-1}/2,w_1,...,w_m).
\]

The w_i are distinct, their product is (-1)^(m-1), and their intermediate
elementary symmetric functions vanish, because product_i(x-w_i)=x^m-1.
Thus z has product one. Its first two entries have absolute values 2
and 1/2, so they are distinct from each other and from the w_i, which
have absolute value one. Equation (1) now gives y_1=y_2 for every
2<=k<=n-2 at the same torus element. This alone proves all the negative
cases over C (indeed over a cyclotomic number field).

For independently checkable rational witnesses, the following construction
avoids number-field arithmetic and also proves the stronger rational claim.
For u=(1,2,...,m-1), define the positive rational number

\[
 c=\frac{e_r(u)}{e_{r-1}(u)}>0,\qquad
 w=(1,2,...,m-1,-c).
\]

Both symmetric functions in the ratio are strictly positive, since
r<=m-1. Multiaffinity proves

\[
 e_r(w)=e_r(u)-c e_{r-1}(u)=0.
\]

All entries of w are distinct and nonzero. Let F=(m-1)! and set

\[
 a=m+1+\frac1{c^2F},\qquad b=-\frac1{acF},\qquad
 z=(a,b,1,2,...,m-1,-c).
\]

Then product_i z_i=ab(-cF)=1. Also a>m-1, and
-c<b<0 because ac^2F>1. Thus all n entries are distinct and nonzero:
a is larger than every positive entry of w, while b lies strictly
between its sole negative entry and zero. Hence diag(z) belongs to T_reg.
But (1) gives y_1-y_2=(a-b)e_r(w)/I_{n,k}=0. The image is nonregular.
This proves failure in every interior case.

For example, n=4,k=2 gives c=1 and

\[
 z=(4,-1/4,1,-1),\qquad
 y=(0,0,15/8,-15/8).
\]

This is a different rational witness from the SL4 example in the
[literature audit](CAYLEY_LITERATURE_AUDIT.md); both are checked by the
trace definition.

## 3. Endpoints, including nonsemisimple elements

At k=1,

\[
 \theta_1(g)=g-\frac{\operatorname{tr}(g)}n I.
\]

Adding a scalar matrix leaves its centralizer unchanged. Thus
Z_G(theta_1(g))=Z_G(g) for every g, proving global regularity preservation.
This argument and conclusion are already in KM Section 6.2.

For k=n-1, the SL_n-equivariant identification exterior^(n-1)(C^n)
with the dual defining representation gives rho(g)=g^(-t),
d rho(X)=-X^t and the ordinary trace form. Therefore

\[
 \theta_{n-1}(g)=-g^{-1}
                    +\frac{\operatorname{tr}(g^{-1})}n I.
\]

Since Z_G(g)=Z_G(g^{-1}), this also preserves the entire centralizer
for every g, including nonsemisimple elements. Combined with the
interior torus counterexamples, this proves all three equivalent
statements of the theorem. For n=2 the endpoints coincide, and for
n=3 there are no interior cases; both edge cases are included.

## 4. Consequence for harmonic separation and limitations

KM Theorem 3.2 is a published equivalence, not a new harmonic theorem:
the multiplication map

\[
 \mathbb C[G]^G\otimes\theta_k^*\operatorname{Harm}(\mathfrak{sl}_n)
 \longrightarrow\mathbb C[G]
\]

is an isomorphism exactly when theta_k preserves all regular elements.
Thus among these fundamental representations this particular harmonic
space works without localization precisely at k=1,n-1. The conclusion
does not rule out a different harmonic complement for an interior power.

The proof does not classify arbitrary highest weights, direct sums,
representations of other groups, or positive-characteristic maps. It
does not determine any generic degree or Smith factor. Standard trace
normalization and the group SL_n, rather than an unspecified isogenous
image, are used throughout.

## 5. Verification and literature boundary

The [independent check](tools/verify_exterior_regularity.py) enumerates the
actual exterior-representation eigenvalues and simple-coroot weights,
builds their full trace Gram matrix and solves the projection equations
using exact rational elimination. This reconstructs theta without using
the displayed closed coordinate formula. It compares every coordinate
and every root difference at deterministic regular rational test points and at all
interior witnesses for 2<=n<=12, 1<=k<=n-1. A separate sparse-monomial
check verifies the pair-difference cancellation for the same finite range.
Direct exterior matrices on nonsemisimple SL3 and SL4 elements check the
dual endpoint formula as well. The report is
[exterior-regularity-verification.json](exterior-regularity-verification.json).

The arbitrary-rank proof is Sections 1--3, not extrapolation from these
finite checks. This is an internal argument audit, not a separate human
or agent review. Targeted searches on 22 October 2026 for Springer/Cayley
regularity and exterior/fundamental representations did not locate this
exact classification. The negative search is not evidence of originality.
The [broader source audit](CAYLEY_LITERATURE_AUDIT.md) remains applicable.

The next question at this stage was arbitrary irreducible SL_n highest
weights. The subsequent [positive-character proof](TYPE_A_REGULARITY.md)
now settles it: only the same defining and dual endpoints preserve
regularity. Other root systems and reducible representations remain
separate questions; literature priority remains unverified.
