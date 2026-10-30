# Polynomial spin characters: degrees, descent and collapsed fibers

Date: 22 October 2026. Status: uniform written deduction extending the
tensor-power calculation, with independent finite exact checks. Published
inputs, the deduction, and its scope are separated below. External review
and literature priority are not established.

## 1. Exact hypotheses and statement

Work over C. Let G=Spin_N, N>=3, r=floor(N/2), and K=2^r.
Let S be irreducible spin for odd N and full spin S+ direct-sum S-
for even N. Write u=chi_S and c=u/K. Identify Lie(G) with so_N by
the differential of pi:G -> SO_N. Spin_4 is included as a semisimple
group; a single half-spin representation is excluded throughout.

Let F in Z[u] have degree q>=1 and assume that F(chi_S) is the
character of an **actual** rational representation rho_F. This includes
nonnegative polynomial mixtures of tensor powers, but also some
polynomials with negative coefficients. It does not mean arbitrary
virtual characters, nor all representations of Spin_N. Nonconstant F
has a nonconstant character and hence rho_F is nontrivial.

Use the trace normalization

\[
 B_F(A,H)=\operatorname{tr}(d\rho_F(A)d\rho_F(H)),\qquad
 B_F(\theta_F(g),H)=\operatorname{tr}(\rho_F(g)d\rho_F(H)).
\]

**Theorem.** Under these hypotheses,

\[
 \boxed{\deg\theta_F=2r\deg F.}                              \tag{1}
\]

The torus and adjoint-quotient generic degrees are the same. The map
descends through pi to SO_N if and only if F is even. In that case
rho_F itself descends and the SO map has degree r deg F.
There is an explicit nonempty open target over which theta_F is finite
flat of rank 2r deg F, described in Section 3.

Within this character family, for nontrivial actual rho_F and rho_H,

\[
 \deg\theta_{\rho_F\oplus\rho_H}
     =\max(\deg\theta_F,\deg\theta_H),\qquad
 \deg\theta_{\rho_F\otimes\rho_H}
     =\deg\theta_F+\deg\theta_H.                            \tag{2}
\]

Multiplicities of a representation and added trivial summands do not
affect its map. The minimum degree within this specified nontrivial
polynomial-character family is 2r, attained exactly when deg F=1.
This is not a minimum over all representations or all rational
equivariant maps. There is no parabolic or cohomology lattice here.

## 2. The character derivative and normalization

Since G is semisimple, du at identity is zero. Differentiating the
character twice gives B_F=F'(K)B_S. The trace form of S is nondegenerate;
actuality and compact-form positivity imply F'(K)>0. Indeed, if it
vanished the actual representation would have zero trace form and be
trivial on the connected semisimple group, contradicting nonconstant F.
This argument also covers Spin_4.

Differentiating at an arbitrary g gives the global identity

\[
 \theta_F(g)=\frac{F'(u(g))}{F'(K)}\theta_S(g).              \tag{3}
\]

In the [Spin double-cover Cayley chart](SPIN_CAYLEY_DEGREES.md),

\[
 X=(\pi(g)-I)(\pi(g)+I)^{-1},\quad
 c^2\det(I-X)=1,\quad c\ne0,\quad \theta_S=2cX.
\]

Define

\[
 h(c)=\frac{cF'(Kc)}{F'(K)}.
\]

Then h is a polynomial of degree q, h(0)=0, h(1)=1, and

\[
 \theta_F=2h(c)X.                                         \tag{4}
\]

Equation (3) is a character-chain-rule consequence of the trace
construction, consistent with Kostant--Michor's direct-sum and tensor
identities. It is not a new definition of the Springer map.

## 3. Fiber polynomial and absence of degree cancellation

Put epsilon=N-2r and define the even characteristic polynomial Q_Y by

\[
 \det(tI-Y)=t^\epsilon Q_Y(t),\qquad
 Q_Y(t)=t^{2r}+a_1(Y)t^{2r-2}+\cdots+a_r(Y).
\]

For h(c)!=0, a fiber over Y has X=Y/(2h(c)). Substitution in the
double-cover equation gives c^2 Q_Y(2h(c))=(2h(c))^(2r).
Because h is divisible by c, division by c^2 gives the polynomial

\[
 \boxed{P_{F,Y}(c)=Q_Y(2h(c))
                -4^r c^{2r-2}\bigl(h(c)/c\bigr)^{2r}.}    \tag{5}
\]

Here h(c)/c denotes a polynomial, including at c=0. If ell is the
leading coefficient of h, the two terms have respective degrees
2rq and 2rq-2. Thus the leading coefficient of P is 4^r ell^(2r),
a nonzero constant independent of Y. **There can be no cancellation
of its leading term**, even if F has negative coefficients.
The unavoidable factor c^2 in the uncleared equation has already been
removed. Counting that uncleared polynomial would overcount by two.

Its constant term is

\[
 P_{F,Y}(0)=
 \begin{cases}a_r(Y),&r>1,\\
 a_1(Y)-4h'(0)^2,&r=1.
 \end{cases}                                              \tag{6}
\]

Choose the nonempty target open

\[
 U_F=\{Y:a_r(Y)P_{F,Y}(0)\ne0\}.                          \tag{7}
\]

**Proposition.** Over U_F the entire map has the presentation

\[
 \theta_F^{-1}(U_F)\simeq
 \operatorname{Spec}_{U_F}
       \bigl(\mathcal O_{U_F}[c]/(P_{F,Y}(c))\bigr).        \tag{8}
\]

**Proof.** In the right-hand ring c is a unit because P(0) is a unit.
Moreover h(c) is a unit: modulo h, equation (5), after inverting c,
reduces to a_r(Y)=0, impossible over U_F. Equivalently the ideals
(P,h) and (P,c) are the unit ideal in the localized polynomial ring.
We may therefore reconstruct X=Y/(2h(c)); equation (5) gives
det(I-X)=c^(-2). The double-cover chart reconstructs a unique Spin
element with these X and c, regularly and with regular inverse.

No source points outside this chart were lost. The identity
P_{F,theta_F(g)}(c(g))=0 holds on the dense open c h(c)!=0 and hence
on all of G. If c(g)=0 it forces P_{F,theta_F(g)}(0)=0. If c(g)!=0
but h(c(g))=0, (4) gives theta_F(g)=0, so a_r(theta_F(g))=0.
Both cases are excluded by (7). This proves (8). QED.

Dividing P by its constant leading coefficient gives a monic polynomial
of degree 2rq. Hence (8) is finite flat with basis 1,c,...,c^(2rq-1).
It is a nonempty open subset of the irreducible group G. Its function
field has degree 2rq over C(g), with c a primitive generator. In
characteristic zero the extension is separable, giving (1) and 2rq
distinct points in a generic fiber. This proof does not infer generic
degree from any finite collection of specialized fiber counts.

For generic regular semisimple Y, the reconstruction stays in the
preimage of its SO maximal torus. Thus the torus degree agrees with
the full degree. The quotient-degree equality is also the published
Kostant--Michor Corollary 3.3.

## 4. Descent and operations on representations

The deck involution fixes X and negates c. Consequently theta_F
descends to SO_N exactly when h(-c)=h(c). A nonconstant character c
has dense image in A^1, so this is a polynomial identity, equivalent
to F having no odd powers. Conversely, for actual rho_F, an even F
gives chi_rho_F(-g)=chi_rho_F(g); evaluating at identity shows that
the central involution has no negative eigenspace. Thus rho_F itself
descends. The degree on SO_N is half the Spin degree.

To prove (2), write F=a_q u^q+lower terms. The weight
(q/2,...,q/2) has multiplicity one in S^(tensor q) and is absent
from every lower tensor power. Its multiplicity in F(chi_S) is a_q.
Actuality forces a_q>0. Thus two actual polynomial characters cannot
cancel their highest coefficient when added, while their product has
the sum of their degrees. Apply (1). This argument does not assume
all coefficients of F are nonnegative.

Actual examples with negative coefficients include u^2-1,
(u^2-1)(u+1), u^3-u and (u^2-1)^2. Indeed S is self-dual (including
full spin in even dimension), so S tensor S contains a trivial
summand. Remove one such summand and then use actual sums and tensor
products. Nonnegative weight multiplicities alone would not certify
that an arbitrary polynomial is an actual character.

## 5. Mixing changes special fibers even when the degree law persists

Inside the open c!=0, (4) gives the exact set-theoretic zero fiber

\[
 \theta_F^{-1}(0)\cap G_{c\ne0}
   =\{1,-1\}\ \cup\!
     \bigcup_{\substack{a\ne0\\F'(a)=0}}\{g: \chi_S(g)=a\}.  \tag{9}
\]

Here 1,-1 are the two elements over the identity of SO_N. They
correspond to X=0,c=+/-1. Equation (9) asserts no scheme multiplicities
and does not classify the extra fiber points at c=0.

For example F(u)=u+u^2 is the actual representation S direct-sum
S^(tensor 2). Its degree is 4r, while the entire level chi_S=-1/2
maps to zero. This level is nonempty: setting all but one torus
coordinate to identity gives c=(t+t^(-1))/2, which attains every
complex value. Since the character is nonconstant on the smooth
irreducible group, each nonempty level is a hypersurface. Thus the
zero fiber contains a divisor and this mixed map is not globally
finite. Pure powers have no nonzero root of F' and no such additional
zero-fiber levels within this chart.

This is the distinction the generic degree alone misses: mixtures
preserve the degree rule but can collapse whole character levels.

## 6. Evidence, published inputs and remaining questions

The [verifier](tools/verify_spin_polynomial_degree.py) constructs weights
of F(S) by tensor convolution and addition, then computes trace moments
on the actual Spin torus b^2=product z_i. It uses the parity of each
individual weight; a mixture has no single tensor-power parity.
It checks reconstructed fibers simultaneously in finite QQ algebras,
computes original Laurent Groebner fibers independently of (5), checks
SO descent, and tests the collapsed-level example directly from weights.
The [report](spin-polynomial-cayley-degrees-verification.json) gives
all finite ranges and source hashes. No runtime API was changed.

Published inputs: the trace-projection construction, character/trace
identities, the Spin double-cover Cayley chart and base degree, and
quotient-degree equality in [Kostant--Michor](https://arxiv.org/pdf/math/0109066),
Theorem 2.8, Corollary 3.3 and Sections 7.4--7.15. We rechecked those
passages; this is not an exhaustive literature-priority search. The
present theorem follows by character differentiation and elimination
in their chart. In rank one it agrees with their SL_2 calculations.

The new repository deduction is the polynomial-character formula,
the explicit open fiber presentation, its descent/operation laws, and
the stated special-fiber description. It is not claimed as a novel
published theorem or a resolution of a published open problem.

The subsequent [geometry note](SPIN_CAYLEY_GEOMETRY.md) determines the
maximal finite target locus of the full map, gives its global Jacobian,
and identifies genuine ramification on the finite locus. It also shows
why extending this scalar discriminant across its chart can introduce
a false branch component in Spin_4. Adjoint-quotient finiteness remains
a separate question; equality of generic degrees does not settle it.
