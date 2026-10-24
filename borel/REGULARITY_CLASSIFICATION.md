# Regularity table for all complex simple groups and actual representations

Date: 22 October 2026. Status: a written uniform deduction with an internal
argument audit and exact finite corroboration. Publication priority and
external review are not established. The known trace construction, positive
weight multiplicities, elementary Laurent factorization and standard root
systems are inputs. This note completes the previously unfinished table;
it does not claim a new published theorem or a formal proof certificate.

## 1. Scope, notation and the completed table

Let G be a connected simply connected complex simple algebraic group,
T a maximal torus, P=X^*(T) its weight lattice and Phi its reduced root
system. Let rho be an actual, nontrivial finite-dimensional rational
representation, allowing arbitrary direct sums and trivial summands.
Its differential is injective and its trace form is nondegenerate. Define
the normalized trace-projection Springer morphism by

\[
 B_\rho(\theta_\rho(g),H)=\operatorname{tr}(\rho(g)d\rho(H)),
 \qquad B_\rho(X,H)=\operatorname{tr}(d\rho(X)d\rho(H)).
\]

Torus preservation means theta(T_reg) is contained in t_reg; global
preservation means theta(G_reg) is contained in g_reg. Here regular means
centralizer dimension rank G. The first tests only regular semisimple
inputs; the second includes nonsemisimple regular elements. Put
D(g)=det(I+g) when a defining symplectic representation is being used.

**Theorem.** The following table is complete for these hypotheses.

| Simply connected group | Representations preserving torus regularity | Representations preserving global regularity |
| --- | --- | --- |
| SL_2=Sp_2 | Actual characters A D^m+B, A,m positive integers | Copies of the defining representation plus trivial summands |
| SL_N, N>=3 | Copies of the defining representation, or copies of its dual, plus trivial summands; no mixture of the two | Exactly the same |
| Sp_(2n), n>=2 | Actual characters A D^m+B, A,m positive integers | None |
| Spin_(2n+1), n>=3 | None | None |
| Spin_(2n), n>=4 | None | None |
| G2, F4, E6, E7, E8 | None | None |

There is at least one nontrivial copy in each positive entry. In A D^m+B,
B is an integer for which the expression is an actual character; it may
remove available trivial summands. Tensor powers of the total exterior
algebra of the defining representation give D^m. Low-rank identifications
are essential: Spin3=SL2, Spin5=Sp4, Spin6=SL4. Spin4 is not simple and
is outside this table.

**Corollary.** Among connected complex simple groups of arbitrary isogeny
form, a globally regularity-preserving representation exists exactly for
SL_N itself, and the positive representations are those just listed.
None of the positive torus-preserving representations above descends
to a nontrivial central quotient. The table therefore also determines
all torus-preserving cases for the other simple isogeny forms.

Virtual characters with negative nontrivial multiplicities, characteristic
p, real groups and general reductive/product groups are not classified.
Neither a critical point nor failure of an integral Borel-image computation
is being used as the definition of regularity failure.

## 2. Positivity along every real root direction

Fix a positive definite Weyl-invariant inner product on the real weight
space. The trace form is a positive scalar c_rho times that form on the
Cartan. For a root alpha put f_alpha(t)=alpha(theta_rho(t)). The known
weight formula gives

\[
 f_\alpha(t)=c_\rho^{-1}\sum_\mu m_\mu(\alpha,\mu)e^\mu(t),
 \qquad c_\rho>0.
\]

For a real Cartan element H, pair a weight mu with s_alpha mu. If
a=(alpha,mu)>0, their contribution to f_alpha(exp H) is

\[
 c_\rho^{-1}m_\mu a
     \left(e^{\mu(H)}-e^{s_\alpha\mu(H)}\right),\qquad
 (\mu-s_\alpha\mu)(H)=\frac{2a}{(\alpha,\alpha)}\alpha(H).
\]

All these contributions have the sign of alpha(H); one is nonzero
when alpha(H)!=0, since the differential is faithful. Consequently

\[
 \boxed{\quad f_\alpha(\exp H)=0\ \Longleftrightarrow\ \alpha(H)=0.
 \quad}                                                     \tag{1}
\]

This holds on every real root wall, not just in a chosen dominant chamber.
No positivity of arbitrary character coefficients in other coordinates
is assumed. Actual nonnegative weight multiplicities are essential.
Also d theta_1=id, directly from the defining trace identity.

## 3. Primitive roots force binomial root coordinates

Suppose every root is primitive in P. Then C[P] is a Laurent polynomial
unique-factorization domain and each e^beta-1 is irreducible. If theta
preserves torus regularity, the zeros of f_alpha lie in the union of
the root hypertori e^beta=1. Unique factorization and the Nullstellensatz
therefore give a monomial unit times a product of these root factors.

No factor with beta not parallel to alpha can occur: choose a real H
with beta(H)=0 and alpha(H)!=0. At exp H that factor vanishes, whereas
(1) says f_alpha does not. Since the root system is reduced, only the
factor e^alpha-1 remains, up to a monomial unit. Thus

\[
 f_\alpha=c e^\lambda(e^\alpha-1)^k.
\]

At identity the derivative of f_alpha is alpha because d theta_1=id.
Hence k=1 and c=1. We have proved the decisive necessary condition

\[
 \boxed{\quad f_\alpha=e^{u_\alpha}-e^{v_\alpha},\qquad
 u_\alpha-v_\alpha=\alpha.\quad}                           \tag{2}
\]

In particular every root coordinate has exactly two monomials with
coefficients +1 and -1. This conclusion concerns actual characters of
the given group, not numerical evaluations at a chosen torus point.

There is a second useful consequence. A monomial e^mu in the weight
formula has coefficient c_rho^(-1)m_mu(alpha,mu). Contributions from
different summands of rho at the same weight have the same sign, so
they cannot cancel. Therefore the union of all weights not orthogonal
to alpha must consist of exactly the two exponents in (2). This is why
mixing representations cannot evade the primitive-root argument.

## 4. Two local root configurations exclude all non-A primitive cases

We use only linearity f_(alpha+beta)=f_alpha+f_beta when alpha+beta
is a root. If two binomials with exponent differences independent alpha
and beta sum to a binomial of the form (2), exactly one opposite-sign
endpoint must cancel. Equal-sign coincidences would give coefficient
2 or -2, and a second opposite-sign cancellation would make the roots
dependent. Thus the two oriented edges form a path of three distinct
monomials.

### 4a. A root string of length at least three

Suppose alpha,beta,alpha+beta,2alpha+beta are roots and alpha,beta
are independent. Write f_alpha=p-q, where p,q denote distinct monomials.
The path condition forces either f_beta=q-r or f_beta=r-p, with p,q,r
distinct. But then

\[
 2f_\alpha+f_\beta=2p-q-r\quad\text{or}\quad p-2q+r,
\]

which cannot be a binomial of the required form. This excludes all
primitive-root types containing a B2 or G2 root string, in particular
B_n for n>=3, F4 and G2.

### 4b. A three-armed simply laced fork

Let alpha_0 be the central simple root and alpha_1,alpha_2,alpha_3
three neighboring simple roots. Each alpha_0+alpha_i is a root, and
each alpha_0+alpha_i+alpha_j, i!=j, is a root. Fix f_0=p-q. Each leaf
is either q-r_i or r_i-p. Among three leaves two have the same
orientation. Their sum with f_0 is then

\[
 p+q-r_i-r_j\quad\text{or}\quad r_i+r_j-p-q.
\]

All four monomials are distinct: an equality would make a leaf zero,
equal to another leaf, or equal to plus/minus alpha_0. The result has
four terms, contradicting (2) for alpha_0+alpha_i+alpha_j. This excludes
D_n, n>=4, and E6,E7,E8 by their Dynkin forks.

### 4c. Primitivity and the precise exceptions

In simply laced rank at least two, a simple root is 2omega_i minus
the sum of its neighboring fundamental weights. A coefficient -1 makes
it primitive. For B_n, n>=3, and F4, each simple-root column of the
Cartan matrix also has a coefficient -1; for G2 the columns have
gcd(2,1)=1 and gcd(2,3)=1. Weyl transformations preserve primitivity,
so every root is primitive in all the types excluded above and in A.

The exceptions are A1 and C_n: their long roots are twice primitive
weights. B2 is C2 under the low-rank identification. Their extra
nonidentity components e^beta=-1 are exactly why the factor argument
must not be applied to them. Section 6 treats them separately.

Thus all type D and exceptional rows are now closed uniformly for
arbitrary actual representations, not by testing finitely many highest
weights. The type B conclusion also follows and agrees with the earlier
independent mixed-character proof in [the B/C note](MIXED_REGULARITY.md).

## 5. Arbitrary representations of SL_N, N>=3

The roots are primitive, so (2) applies. Decompose rho into irreducibles
and let lambda=sum a_i omega_i be any nontrivial summand. For the
highest root alpha=epsilon_1-epsilon_N, put L=sum a_i. Its highest
weight and reflected weight are lambda and lambda-L alpha. Both have
nonzero pairing with alpha and occur in rho. If L>=2 they cannot be
the two exponents separated by alpha in (2). No other summand can
cancel their coefficients. Therefore every summand has L=1 and is a
fundamental representation exterior^k C^N.

In exterior degree k the number of distinct weights with nonzero
alpha pairing is 2 binomial(N-2,k-1). This equals two only at k=1,N-1.
Hence every nontrivial summand is defining or dual defining.

For N>=3, mixing those two produces four distinct alpha-nonorthogonal
weights: mu_1,mu_N,-mu_1,-mu_N, where
mu_i=epsilon_i-(epsilon_1+...+epsilon_N)/N. They do not coincide,
so (2) prohibits the mixture. We conclude that rho is either

\[
 a\mathbb C^N\oplus b\mathbb C
 \quad\text{or}\quad
 a(\mathbb C^N)^*\oplus b\mathbb C,\qquad a\ge1,\ b\ge0.
\]

These representations have the same normalized maps as one defining
copy or its dual. Their formulas are

\[
 \theta_V(g)=g-\frac{\operatorname{tr}g}{N}I,
 \quad
 \theta_{V^*}(g)=-g^{-1}+\frac{\operatorname{tr}g^{-1}}{N}I.
\]

Both preserve the entire centralizer, including nonsemisimple inputs,
and hence preserve global regularity. The defining case is already
published in Kostant--Michor Section 6.2. The present support argument
extends the previous [irreducible A proof](TYPE_A_REGULARITY.md) to all
actual sums and does not re-prove that published positive example.

## 6. The nonprimitive cases: Sp and rank one

For Sp_(2n), n>=2, the [mixed-representation proof](MIXED_REGULARITY.md)
already shows that torus preservation is equivalent to chi=A D^m+B.
No such map preserves every regular element: at a regular element with
one J_2(-1) block the image is rank-one nilpotent for m=1 and zero
for m>=2. All other representations already fail on the torus.

For SL2=Sp2 use x=z+z^(-1), chi=F(x). The known weight formula gives
y=(z-z^(-1))F'(x)/S with S>0. Torus regularity means y!=0 whenever
x!=2,-2. The one-variable factor lemma and positivity F'(x)>0 on x>=2
force

\[
 F'(x)=c(x+2)^{m-1},\qquad F(x)=A(x+2)^m+B,
 \qquad A\in\mathbb Z_{>0},\ m\ge1.
\]

Integrality of A follows from the highest weight m, whose monomial z^m
has coefficient A. These are actual exterior/tensor characters with an
allowed trivial adjustment. They all preserve torus regularity.

The normalized full map is the same determinant formula as for Sp:

\[
 \theta_m(g)=\frac{2}{4^m}D(g)^{m-1}(g-I)\operatorname{adj}(g+I).
\]

At m=1 it is the defining SL2 map and preserves global regularity.
At m>=2 it vanishes on the regular Jordan block J_2(-1), so it fails
global regularity. The m=1 actual characters are exactly copies of
the defining representation plus trivial summands. This rank-one
mixture phenomenon was outside the earlier irreducible A theorem.

## 7. Other simple isogeny forms

Pulling a representation back along a central isogeny identifies the
normalized trace maps through the differential. Lie centralizer dimensions
and torus regularity agree at lifts. Thus positive cases on a quotient
must occur in the simply connected table and descend through its center.

For SL_N, N>=3, the defining and dual representations have faithful
central character, and neither positive sum can descend to a nontrivial
central quotient. For SL2, the coefficient of z in (z+2+z^(-1))^m is
positive for every m>=1 and is unaffected by B, so the character does
not descend to PSL2. For Sp_(2n), n>=2, D^m contains the weights
(m,...,m) and (m,...,m,m-1), both nonzero; they have opposite central
parity. At least one is odd, remains after any trivial adjustment, and
prevents descent to PSp. All other simply connected cases were already
negative. This proves the isogeny corollary with actual character lattices.

## 8. Verification, known inputs and remaining work

The [exact verifier](tools/verify_regularity_classification.py) and
[report](regularity-classification-verification.json) separately audit:
root primitivity in the full weight lattice; the root strings and fork
roots used in the proof; all endpoint-cancellation patterns; actual
character weights in the missing classical/exceptional cases; exact
algebraic regular-to-nonregular witnesses; the A mixture obstruction;
and rank-one torus-positive/global-negative controls. Independent
Gelfand--Tsetlin versus Freudenthal weights check the A support argument;
root orbits and root-system adjoint weights check exceptional characters.
Finite checks corroborate the written arbitrary-rank argument and are
not an exhaustion of all representations.

Published inputs are the trace construction and weight differentiation,
Weyl invariance and positive multiplicities of actual characters, the
classical root systems, and elementary factorization. The root-coordinate
formula is compatible with [Rogers Theorem 1](https://arxiv.org/pdf/1701.01538v1).
[Kostant--Michor](https://arxiv.org/pdf/math/0109066), Sections 2 and 3,
provide trace equivariance, tensor scaling and the harmonic criterion;
their Section 6.2 already proves the defining SL_N positive case.
[Lu--Ye](https://arxiv.org/pdf/2411.02071) classify applicability of the
classical matrix Cayley transform in a representation, a different
property. In particular their Spin8 half-spin examples do not imply
regularity preservation for these trace-projection morphisms.

The targeted literature searches used here did not locate this exact
all-representation table. That limited search is not evidence sufficient
to assert originality. The primitive-root and local-support proofs are
written deductions; external proof review and priority remain outstanding.
For this scope, the regularity table itself has no unclassified rows.
General reductive/product groups, virtual characters and characteristic p
remain outside the theorem. Minimal harmonic localization for arbitrary
negative maps also remains a separate question; the known-corollary
[symplectic determinant localization](SYMPLECTIC_HARMONIC_LOCALIZATION.md)
does not settle those other localizations. No PDF/LaTeX or runtime API
is changed by this note.
