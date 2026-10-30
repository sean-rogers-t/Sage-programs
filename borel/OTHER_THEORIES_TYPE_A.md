# The defining type A source is a change of Chern coordinates

Exploration, 22 October 2026. This is a useful general baseline: the
defining type A representation has much less new oriented image
arithmetic than a general representation. The ambient flag-bundle
presentation is published; the Springer-image statement is its direct
application. No originality claim is made.

## 1. Statement for every partial flag

Let G=SL_m(C), let rho be its defining representation, let P be any
standard parabolic, and use the literal integral source
A=Sym_Z(Q)^(W_L). Let E be a multiplicative complex-oriented theory
with even coefficients R; R need not be torsion-free. On G/B put

\[
          y_i=-c_1^E(\mathcal L_{-\epsilon_i}).
\]

Use the trace-coordinate lift specified in
[the finite comparison](OTHER_THEORIES_FINITE.md#4-a-uniform-oriented-lift-with-explicit-hypotheses).
There is an integral filtered ring isomorphism

\[
 \Phi_{E,P}:H^{2*}(G/P;\mathbb Z)\otimes R\xrightarrow{\sim}E^*(G/P)
                                                                    \tag{1}
\]

under which the lifted source map is exactly Phi_(E,P) composed with
the base change of the old Borel map. On the full flag Phi sends the
ordinary Chern roots epsilon_i to y_i. On a partial flag it sends the
elementary symmetric functions in a Levi block to the corresponding
signed Chern classes of its dual tautological subquotient.

In particular, for ordinary integral K-theory (after the usual Bott
identification), **all ambient total Smith factors agree with the old
ordinary ones for every P**. The same statement holds for any explicitly
specified subalgebra of this integral source. For characteristic-p
coefficients, its embedded image is the old ordinary matrix reduced
modulo p in these coordinates. Its rank is independent of formal-group
height. This statement concerns the defining rho; it is false for
arbitrary rho, as the SL3 family in the finite note shows.

## 2. Proof and published input

The associated defining dual G-bundle is trivial. Its splitting roots
on G/B are c1^E(L_(-epsilon_i)), so the Whitney formula gives
e_j(y_1,...,y_m)=0 for every j>0. The known type A flag-bundle
presentation therefore gives the map (1). Equivalently one can first
use these relations on the full flag: its map has identity first jet,
and the finite even-cell filtrations make it an isomorphism integrally
over R. Block Chern classes define the partial-flag version; the same
first-jet argument proves it is an isomorphism there. This uses the
integral block Chern presentation, not an assumption that taking Weyl
invariants commutes with modular reduction.

The full-flag cobordism presentation is explicitly
[Hornbostel--Kiritchenko, Theorems 1.1 and 2.6](https://arxiv.org/abs/0903.3936).
Their cellular comparison also connects algebraic and topological
cobordism over C. The partial-flag construction follows the standard
flag-bundle Chern presentation and projective bundle theorem.

For this normalized Springer morphism,

\[
 \theta^*(\epsilon_i-\epsilon_j)=e^{\epsilon_i}-e^{\epsilon_j}.
\]

The trace-coordinate lift consequently sends every root epsilon_i-epsilon_j
to y_i-y_j. These roots generate Q integrally, so polynomial substitution
proves the asserted commutative comparison on A. In K-theory y_i=L_i-1,
and the ambient change of integral bases is unimodular. QED.

## 3. Projective modular images

On P1, put h=c1^E(O(1)) and y=-iota_F(h), a strict coordinate with
first term h. The generators E_j=e_j(epsilon_1-epsilon_2,...,
epsilon_1-epsilon_m) have images binomial(m,j)y^j. This follows directly
from the trivial defining dual bundle and the elementary symmetric
identities, without invoking torsion-free coefficients.

If R is a field of characteristic p and a=v_p(m), Lucas' binomial
criterion shows that every nonzero generator has degree divisible by
p^a, and the first possible one is p^a. If p^a<m its coefficient is a
unit. Therefore the complete image is

\[
       R[y^{p^a}]\ \subset R[y]/y^m,\qquad
       \operatorname{rank}_R\operatorname{im}=m/p^a.               \tag{2}
\]

If m=p^a, (2) means constants only. If p does not divide m, the image
is the whole ring. This is also immediate by reducing the old integer
image factors m/gcd(m,r). It is a base-change deduction, not a new
ambient computation of generalized cohomology of projective space.

Compare (2) with the symplectic defining representation on the same
odd-dimensional projective spaces: there the image rank is
ceil(m/2^s) at characteristic-two height s. The difference comes from
the representation and its actual characters, not from a change in the
underlying space.

## 4. Independent finite check

The [native K-theory audit](experiments/audit_other_theories_k.py)
constructs the integral ambient change of bases for SL3/B and SL4/B
using a full weight Artin basis. Its ordinary coordinates come from
additive BGG operators and its K coordinates from Todd/BGG operators.
Both square basis matrices have determinant plus or minus one, and
their quotient carries every actual defining Springer source column
to its K column. This verifies the full embedded comparison in these
finite examples, independently of simply comparing two Smith lists.
