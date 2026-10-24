# Tensor-power degrees for the Spin generalized Cayley map

Date: 22 October 2026. Status: uniform written deduction from published
Kostant--Michor formulas, with independent exact finite corroboration.
Neither external review nor publication priority is established.

## 1. Representation and theorem

Work over C. Let N>=3, r=floor(N/2), G=Spin_N, and let
pi:G -> SO_N be the standard double cover. Identify Lie(G) with so_N
using d pi. Choose the following representation S:

- For N=2r+1, S is the irreducible spin representation, of dimension 2^r.
- For N=2r, S=S^+ direct-sum S^- is the **full spin representation**,
  of dimension 2^r. A single half-spin representation is not covered.

Kostant--Michor use two copies of the irreducible spin in odd dimension.
Taking repeated copies changes neither the normalized character nor the
trace-projection map, including after tensor powers. Our convention gives
exactly their normalized map. Spin_4 is semisimple rather than simple;
it is included, with both half-spin summands essential to infinitesimal
faithfulness. Spin_2 and its torus are excluded.

For m>=1 put rho_m=S^(tensor m), and define theta_m by

\[
 B_m(A,H)=\operatorname{tr}(d\rho_m(A)d\rho_m(H)),\qquad
 B_m(\theta_m(g),H)=\operatorname{tr}(\rho_m(g)d\rho_m(H)).
\]

The trace form is nondegenerate and d theta_m at identity is the identity.
These are full group-to-Lie-algebra maps; no parabolic, Borel quotient or
integral cohomology lattice is involved.

**Theorem.** For every N>=3 and m>=1,

\[
 \boxed{\deg(\theta_m)=2rm
   =\begin{cases}Nm,&N\text{ even},\\(N-1)m,&N\text{ odd}.\end{cases}}       \tag{1}
\]

This is also the generic degree on a maximal torus and on the adjoint
quotient. For even m the map descends through pi to a map SO_N -> so_N
of degree rm. For odd m it does not descend through pi.

An explicit scalar polynomial below describes the entire map over a
specified nonempty open target and proves that it is finite flat of
degree 2rm there. Global finiteness or surjectivity is not asserted.
This is the degree of the specified representation-induced morphism,
not the minimum Cayley degree among rational equivariant maps.

## 2. Published inputs and the double cover

Put c(g)=chi_S(g)/2^r. Kostant--Michor, Sections 7.5, 7.10, 7.14 and
Theorem 7.15, give, on c!=0,

\[
 a=\pi(g),\quad X=(a-I)(a+I)^{-1},\quad
 \delta(X)=\det(I-X),\qquad c^2\delta(X)=1,\qquad
 \theta_1(g)=2cX.                                             \tag{2}
\]

In their notation Gamma(a)=-X. Their identity
c^2=2^(-N)det(I+a) gives the stated relation because
det(I+a)=2^N/delta(X). Their exterior-exponential inverse identifies
the *whole* open Spin variety with

\[
 \{(X,c)\in\mathfrak{so}_N\times\mathbb G_m:
                         c^2\det(I-X)=1\}.                  \tag{3}
\]

In particular, X alone specifies an SO element, whereas (X,c) specifies
a unique Spin lift. The deck involution is (X,c) -> (X,-c).
Replacing (3) by the SO Cayley chart would lose a factor of two.

The tensor formula in Kostant--Michor Theorem 2.8(3) gives

\[
 \theta_m=c^{m-1}\theta_1=2c^mX.                              \tag{4}
\]

One can also derive it directly: semisimplicity gives tr(dS(H))=0;
the tensor trace form is m(2^r)^(m-1) B_1 and the character derivative
is m chi_S^(m-1) times the original derivative. This works for Spin_4
as well. Equations (2)--(4) are published inputs or immediate rewritings,
not new representation formulas.

## 3. A scalar presentation of the fibers

For Y in so_N set epsilon=N-2r and define the even polynomial Q_Y by

\[
 \det(tI-Y)=t^\epsilon Q_Y(t),\qquad
 Q_Y(t)=t^{2r}+a_1(Y)t^{2r-2}+\cdots+a_r(Y).
\]

The coefficients are invariant polynomial functions of Y. On the Cartan,
with eigenvalues +/-y_1,...,+/-y_r (and an extra zero when N is odd),
Q_Y(t)=product_i(t^2-y_i^2). This identity also holds at nonsemisimple
Y by polynomial continuation. The forced zero eigenvalue for odd N has
already been removed; it must not be counted as a root of a fiber equation.

For a fiber theta_m(g)=Y in the chart (3), X=Y/(2c^m). Thus

\[
 \boxed{P_{m,Y}(c)=Q_Y(2c^m)-4^r c^{2rm-2}=0.}              \tag{5}
\]

This has degree 2rm and leading coefficient 4^r. In diagonal coordinates,

\[
 P_{m,Y}(c)=\prod_{i=1}^r(4c^{2m}-y_i^2)
                                      -4^r c^{2rm-2}.        \tag{6}
\]

Let U_m be the target open defined by P_{m,Y}(0)!=0. Explicitly,

\[
 U_m=\begin{cases}
 \{a_r(Y)\ne0\},&rm>1,\\
 \{a_1(Y)-4\ne0\},&r=m=1.
 \end{cases}                                                \tag{7}
\]

The rank-one exception matters: the second term of (5) is constant
when rm=1. Both opens are nonempty.

**Fiber-presentation proposition.** There is an isomorphism over U_m

\[
 \theta_m^{-1}(U_m)\simeq
 \operatorname{Spec}_{U_m}\bigl(\mathcal O_{U_m}[c]/(P_{m,Y}(c))\bigr).
                                                               \tag{8}
\]

**Proof.** On (3), the substitution X=Y/(2c^m) changes its defining
equation into (5), after multiplication by a unit power of c.
Conversely, any invertible c satisfying (5) gives X=Y/(2c^m),
delta(X)=c^(-2), then a=(I+X)(I-X)^(-1) and its unique lift with
normalized spin character c. These constructions are regular inverse
maps, so the assertion is scheme-theoretic, not just a root count.

There are no extra points outside (3) above U_m. Indeed,
P_{m,theta_m(g)}(c(g))=0 holds on the dense chart and hence everywhere
on the irreducible group G. At c(g)=0 it implies
P_{m,theta_m(g)}(0)=0, which excludes U_m. Finally, in the right-hand
ring of (8), c is already a unit because the constant term of P is a
unit on U_m. This proves (8). QED.

Dividing P by its constant leading coefficient 4^r makes it monic.
Consequently 1,c,...,c^(2rm-1) are a free basis in (8). This proves
finite flatness and degree (1). Characteristic zero makes the finite
function-field extension separable, so a generic fiber has 2rm distinct
points; no assumption about distinct roots at every target is needed.
The same presentation proves

\[
 \mathbb C(G)=\theta_m^*\mathbb C(\mathfrak g)(c),
\]

with (5), up to its constant leading coefficient, the minimal polynomial
over that embedded field. Irreducibility follows from the irreducibility
of the nonempty open subset of G in (8).

## 4. Central parity, torus and quotient

The nontrivial element -1 of ker(pi) acts on S as -I, so
c(-g)=-c(g) and theta_m(-g)=(-1)^m theta_m(g). For even m,
rho_m itself descends to SO_N and theta_m is its trace-projection map
under d pi. Multiplicativity of generic degrees in this degree-two
cover gives deg(theta_m on SO_N)=rm. For odd m, dominance prevents
theta_m=-theta_m identically, so there is no descent through pi.
This statement concerns the kernel of the Spin-to-SO cover, not a
classification of every quotient by the possibly larger center.

For generic regular semisimple Y in a chosen Cartan, the reconstruction
X=Y/(2c^m) and a=(I+X)/(I-X) lies in its SO torus, and both Spin lifts
lie in the preimage torus. Thus all generic full fibers are counted by
the torus map. The adjoint-quotient equality is Kostant--Michor
Corollary 3.3 (or follows by taking these regular semisimple fibers).
Passing to the Weyl quotient does not introduce another factor |W|.

## 5. Independent checks and limits

The [exact Sage verifier](tools/verify_spin_cayley_degree.py) builds S
from all weights (epsilon_1/2,...,epsilon_r/2), epsilon_i=+/-1, and
constructs tensor powers by weight convolution. It computes the trace
Gram matrix and character derivatives directly, without using (4) to
construct those moments. In doubled-weight coordinates k it verifies

\[
 B_m|_{\mathfrak t}=m\,2^{rm-2} I_r.
\]

The actual Spin torus has coordinates z_1,...,z_r,b, all invertible,
with b^2=product_i z_i. Its character lattice is
Z^r+Z(1/2,...,1/2), and b represents the half-sum weight. A weight
k/2 in S^(tensor m), with p=m mod 2, gives the monomial
b^p product_i z_i^((k_i-p)/2). The verifier uses this torus, not
the larger auxiliary torus obtained by independently adjoining every
square root of z_i.

Independent checks consist of direct weight moments, simultaneous
reconstruction in exact finite fiber algebras, central-deck parity,
and small original Laurent Groebner fibers, saturated only by torus
coordinates. Even-power SO fibers are counted separately. The
[machine-readable report](spin-cayley-degrees-verification.json) records
the finite ranges and source hashes. These checks corroborate the
written uniform proof; finite specializations alone would not prove it.

## 6. Published versus deduced, and next question

Published: the trace map, tensor-scaling identity, Spin double-cover
chart and scalar Cayley formula, m=1 degree 2r, and equality with the
quotient degree. Source: [Kostant--Michor, The generalized Cayley map
from an algebraic group to its Lie algebra](https://arxiv.org/pdf/math/0109066),
Theorem 2.8(3), Corollary 3.3 and Sections 7.4--7.15.

Deduced here: applying those formulas uniformly to tensor powers yields
degree 2rm, the explicit presentation (5)--(8), and the parity-dependent
SO degree rm. The proof deliberately follows KM's m=1 method. The
rank-one degree law is also covered by their SL_2 computations. This
is a useful extension and benchmark, not a claim to solve a published
open problem or to establish a new literature result.

The [symplectic determinant family](SYMPLECTIC_CAYLEY_DEGREES.md) has
a comparable linear tensor-power law. The subsequent
[polynomial spin-character proof](SPIN_POLYNOMIAL_CAYLEY_DEGREES.md)
extends this calculation to every actual character F(chi_S): degree
2r deg(F), with descent to SO precisely for even F. It retains the
double cover and proves that the leading term cannot cancel. It also
identifies additional character levels that mixtures can collapse to
zero; this special-fiber behavior is invisible in the generic degree.
