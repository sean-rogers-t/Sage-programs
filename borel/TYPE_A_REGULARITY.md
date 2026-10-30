# Regularity classification for irreducible SL_n Springer maps

Date: 22 October 2026. Status: a uniform written deduction from the trace
construction, positive character multiplicities and the exterior-power
argument. Internal exact corroboration is provided separately. Literature
priority and external review remain unverified. This note concerns C and
ordinary regular elements, not positive characteristic or integral images.

## 1. The statement

Let G=SL_n(C), n>=2, and let V_lambda be a nontrivial finite-dimensional
irreducible rational representation with dominant highest weight

\[
 \lambda=\sum_{i=1}^{n-1}a_i\omega_i,\qquad a_i\in\mathbb Z_{\ge0}.
\]

Its differential is faithful, although its group kernel may be finite.
Let theta_lambda:G -> sl_n be the normalized trace Springer map. Let T
be the diagonal determinant-one torus. Regular means centralizer dimension
n-1; on T, it means that all diagonal entries are distinct.

**Theorem.** The following statements are equivalent:

1. theta_lambda(T_reg) is contained in t_reg.
2. theta_lambda(G_reg) is contained in (sl_n)_reg.
3. lambda is omega_1 or omega_(n-1).

At n=2 these two weights coincide. In every other irreducible case there
is a regular semisimple torus element with algebraic entries whose image
has a repeated diagonal entry. Thus this is a classification for all n
and all irreducible highest weights, not just symmetric powers or finite
rank calculations. The symmetric-power corollary is

\[
 \theta_{\operatorname{Sym}^d(\mathbb C^n)}
 \text{ preserves torus/global regularity}\quad\Longleftrightarrow\quad d=1
 \qquad(d\ge1).
\]

The same holds for the duals of these symmetric powers. The trivial
representation is excluded: its trace form does not define this projection.
Reducible representations are outside the theorem's scope.

## 2. Published inputs and normalization

The trace construction is from
[Kostant--Michor](https://arxiv.org/pdf/math/0109066), Section 2.1;
the explicit general weight formula is from
[Rogers](https://arxiv.org/pdf/1701.01538v1), Theorem 1. These are published
inputs, not new claims of this note. The standard global regularity case
is already in KM Section 6.2.

Use the unique polynomial GL_n lift indexed by the partition

\[
 p=(p_1,\ldots,p_n),\qquad
 p_j=\sum_{i=j}^{n-1}a_i,\quad p_n=0.
\]

Its character chi=s_p has nonnegative integer monomial coefficients,
is symmetric in its n arguments, is homogeneous of degree D=|p|,
and contains the monomial z^p with coefficient one. These are classical
highest-weight/Schur-character facts; for the polynomial-character and
tableau description see [Stanley, GL(n,C) for combinatorialists](https://math.mit.edu/~rstan/pubs/pubfiles/57.pdf),
Section 1, especially Theorem 1.2. The highest tableau has every entry
in row i equal to i.

There is a positive trace index I such that

\[
 \operatorname{tr}(d\rho(X)d\rho(Y))=I\operatorname{tr}(XY)
 \quad(X,Y\in\mathfrak{sl}_n).
\]

Indeed invariance and simplicity make the forms proportional, and for
H=diag(1,0,...,0,-1), the weight formula gives
2I=sum_mu m_mu(mu_1-mu_n)^2>0. Nontriviality ensures positivity.
Consequently, for product_i z_i=1, the diagonal output y is

\[
 y_i=I^{-1}\left(z_i\partial_i\chi(z)-\frac Dn\chi(z)\right),
 \qquad
 y_i-y_j=I^{-1}(z_i-z_j)Q_{ij}(z),
\]

where

\[
 Q_{ij}(z)=\frac{z_i\partial_i\chi(z)-z_j\partial_j\chi(z)}{z_i-z_j}.
                                                               \tag{1}
\]

Although theta is defined on SL_n, the displayed polynomials can also
be evaluated on any nonzero tuple before imposing determinant one.
They are homogeneous: y has degree D and Q has degree D-1. Multiplying
the GL_n lift by a determinant character multiplies Q by a nonvanishing
monomial and has no effect on the torus zero test.

## 3. A positive-coefficient divided-root lemma

**Lemma.** If chi is a symmetric polynomial with nonnegative real
monomial coefficients, then Q_ij in (1) is a polynomial with nonnegative
real coefficients.

**Proof.** Fix the exponents in the remaining variables. Symmetry pairs
terms c z_i^a z_j^b and c z_i^b z_j^a, where a>b and c>=0 (times the
fixed remaining monomial). Their contribution to Q is

\[
 c(a-b)(z_i z_j)^b\sum_{s=0}^{a-b-1}z_i^{a-b-1-s}z_j^s.
                                                               \tag{2}
\]

All coefficients are nonnegative. Terms with a=b contribute zero.
Summing the pairs proves the assertion, including divisibility by z_i-z_j.
This elementary lemma is a deduction, not a new positive-character theorem.

## 4. Every nonfundamental highest weight has a torus obstruction

Suppose L=p_1=sum_i a_i>=2. The highest monomial and its transposition
in coordinates 1,n contribute to Q_(1,n) the polynomial

\[
 L\left(\prod_{j=2}^{n-1}z_j^{p_j}\right)
       \sum_{s=0}^{L-1}z_1^{L-1-s}z_n^s.                       \tag{3}
\]

The character contains both monomials with coefficient one. By the
lemma no other contribution can cancel any term in (3).
Specialize z_j=j for 1<=j<=n-1 and put

\[
 q(x)=Q_{1n}(1,2,\ldots,n-1,x).
\]

This polynomial has nonnegative rational coefficients, a strictly positive
constant term, and degree at least L-1>=1. It therefore has a complex
algebraic root x_0. No such root is zero or positive real: q(0)>0, and
q(x)>0 for x>0. In particular x_0 differs from 1,2,...,n-1.
The tuple v=(1,2,...,n-1,x_0) has distinct nonzero entries and
Q_(1,n)(v)=0.

Choose any algebraic c with

\[
 c^n=\frac1{(n-1)!x_0},\qquad z=cv.
\]

Then product_i z_i=1, its entries are still distinct, and homogeneity
gives Q_(1,n)(z)=c^(D-1) Q_(1,n)(v)=0. Equation (1) proves y_1=y_n.
Thus diag(z) is regular semisimple and its image is nonregular. This
argument covers every highest weight with sum_i a_i>=2, including all
higher symmetric powers. The use of an algebraic root and scaling is
essential to the stated field: a rational witness is not asserted here.

## 5. Fundamental weights and the positive global cases

The only remaining nonzero highest weights have sum_i a_i=1, hence
lambda=omega_k and V_lambda=exterior^k(C^n). The
[exterior-power proof](EXTERIOR_REGULARITY.md) gives

\[
 Q_{ij}=e_{k-1}(z_{\widehat i\widehat j}).
\]

Every interior k=2,...,n-2 has a distinct-entry determinant-one torus
counterexample. A simultaneous complex witness uses the n-2 roots of
unity together with 2 and (-1)^(n-3)/2; a rational construction is also
given in that note. This proves all the remaining negative cases.

For the endpoints, for every g in G,

\[
 \theta_{\omega_1}(g)=g-\frac{\operatorname{tr}g}{n}I_n,
 \qquad
 \theta_{\omega_{n-1}}(g)=-g^{-1}
             +\frac{\operatorname{tr}(g^{-1})}{n}I_n.
\]

Scalar shifts and inversion preserve the entire centralizer. These
formulas therefore prove global regularity preservation, including
nonsemisimple elements. Global preservation implies torus preservation,
so all three statements of the theorem are now proved.

## 6. Symmetric powers explicitly

For Sym^d(C^n), chi=h_d is the complete homogeneous symmetric polynomial,
and the trace index is

\[
 I_{n,d}=\binom{n+d}{n+1}.
\]

One direct count uses weak compositions a of d. The generating function
for sum_a(a_1^2-a_1 a_2) is t/(1-t)^(n+2), giving this index.
The character generating function product_i(1-z_i t)^(-1) gives

\[
 Q_{ij}=h_{d-1}(z_1,\ldots,z_n,z_i,z_j).                     \tag{4}
\]

Here z_i and z_j occur twice each in the alphabet. Thus a concrete
witness polynomial for d>=2 is

\[
 q_{n,d}(x)=h_{d-1}(1,2,\ldots,n-1,x,1,x).
\]

Its degree is d-1, its constant term is positive, and every coefficient
is positive. Any root, followed by the determinant normalization of
Section 4, provides a regular-torus counterexample. For n=2,d=2,
q=2+2x, giving diag(i,-i) in SL2 with theta=0. This rank-one case is
also consistent with KM's already published all-SL2 formulas in Section 6.1;
it is not a new rank-one construction.

For a higher-rank example, n=3,d=2 gives q=4+2x. With c^3=-1/4,

\[
 t=c\operatorname{diag}(1,2,-2)\in SL_3,\qquad
 \theta_{\operatorname{Sym}^2}(t)
 =c^2\operatorname{diag}(-4/15,8/15,-4/15).
\]

The input has three distinct eigenvalues and the output has a repeated
eigenvalue. The trace index is five; a direct monomial-weight trace check
in the verifier confirms the displayed output.

For any representation with nondegenerate trace form, directly from
the trace identity,

\[
 \theta_{\rho^*}(g)=-\theta_\rho(g^{-1}).
\]

Thus regularity preservation is also invariant under taking the dual.
This is a general elementary consequence, not a separate classification.

## 7. Consequences, checks and limits

**Central-quotient corollary.** If C is a nontrivial subgroup of the
center of SL_n and G'=SL_n/C, then no nontrivial irreducible representation
of G' has a trace Springer map preserving torus regularity (or global
regularity). To see this, pull the representation back along the central
isogeny pi:SL_n -> G'. Under the differential Lie-algebra isomorphism,
trace projection commutes with pi. The adjoint actions also correspond,
so Lie centralizer dimensions, and therefore regularity, agree at lifts.
The classification would force the pullback to be defining or dual,
but neither is trivial on any nontrivial subgroup C of the center.
Thus neither descends, a contradiction. This is an elementary corollary
of the theorem, not a new isogeny or representation-descent theorem.

KM Theorem 3.2 already identifies global regularity preservation with
the multiplication isomorphism

\[
 \mathbb C[G]^G\otimes\theta_\lambda^*\operatorname{Harm}(\mathfrak{sl}_n)
 \longrightarrow\mathbb C[G].
\]

For irreducible SL_n representations, this particular harmonic space
therefore works without localization exactly at omega_1,omega_(n-1).
This is an application of a published equivalence. It does not deny
existence of other harmonic complements, alter the polynomial-character
or Borel maps in the Kumar papers, or compute a Cayley degree or an
integral image. Beyond the stated type A central quotients, no classification
of other root systems or reducible
representations is asserted.

The [exact verifier](tools/verify_type_a_regularity.py) reconstructs
characters by Gelfand--Tsetlin interlacing weights, and independently
evaluates Jacobi--Trudi determinants using polynomial dual numbers.
It compares characters and logarithmic character derivatives, checks
the positive divided-root expansion coefficient by coefficient, and
checks witness-polynomial divisibility using actual weight traces.
Full simple-coroot Gram solves additionally check trace normalization
at distinct-entry rational determinant-one points. The checked finite
scope is recorded in [the report](type-a-regularity-verification.json).
The arbitrary-rank proof is Sections 2--5, not an inference from tests.
These are independent computations within an internal audit, not a
separate human review or a formal proof certificate.

Targeted searches on 22 October 2026 for Springer/Cayley regularity,
irreducible representations and symmetric powers did not locate this
exact classification. The published construction, character facts,
standard case and harmonic criterion are explicitly attributed above;
the absence of a located classification is not a novelty claim.
The [literature audit](CAYLEY_LITERATURE_AUDIT.md) records the broader boundary.

Next question: how much of this classification extends to the other
simple root systems? The type A proof uses the polynomial Schur-character
lift, permutation symmetry, and a common scalar determinant normalization.
These steps must be replaced or justified before applying the argument
to another group; it does not settle that extension automatically.
