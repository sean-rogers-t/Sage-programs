# Equivariant K-images and the root part of the Cayley determinant

Exploration, 22 October 2026. This is a characteristic-zero deduction
from the published Borel/Pittie--Steinberg and localization frameworks,
applied to the explicit Springer morphism. It relates the earlier
quotient-coordinate and Cayley-determinant calculations directly to a
new image problem. It does not give all integral Fitting ideals or a
new construction of ambient equivariant K-theory. No priority claim is
made for this deduction.

## 1. The source family uses the quotient coordinates already computed

Let G be finite, simply connected and simple; rho has faithful
differential and theta=theta_rho has the normalization theta(1)=0,
d theta_1=id. Fix P with Weyl group W_L. Put
S=Sym_Q(X_Q), A=S^(W_L), B=S^W and R=Q[X]=R(T)_Q.
The representation-dependent coefficient map is

\[
 B\longrightarrow R,\qquad f\longmapsto f(\theta(t)).
\]

For basic invariants these are exactly the higher adjoint-quotient
coordinates from the earlier geometry project. Define the finite free
R-algebra

\[
       \mathcal A_{\rho,T}^P=R\otimes_{B,\,f\mapsto f(\theta(t))}A.
                                                                    \tag{1}
\]

Its rank is n=|W/W_L|. The equivariant associated-bundle map gives

\[
 \eta_{\rho,T}^P:\mathcal A_{\rho,T}^P\longrightarrow K_T^0(G/P)_\mathbb Q,
                    \quad f\longmapsto\kappa_T(\theta^*f).         \tag{2}
\]

Global characters now become their restrictions to T, rather than
their dimensions, so the balancing in (1) is essential. At t=1 the
source becomes the rational ordinary Borel source quotient. This
statement is over Q; a chosen integral source form can have additional
arithmetic issues and need not be free over its global invariant ring.

Published inputs for the target and restriction methods are
[Griffeth--Ram, Sections 1--2](https://arxiv.org/abs/math/0405333)
and [Willems's equivariant restriction formulas](https://arxiv.org/abs/math/0204265).

## 2. The root quotients

For each root-length orbit a, let Phi_a^+ be its positive roots and put

\[
 \Delta_a(x)=\prod_{\alpha\in\Phi_a^+}\alpha(x),\qquad
 \rho_a=\tfrac12\sum_{\alpha\in\Phi_a^+}\alpha,\qquad
 \delta_a(t)=e^{\rho_a}\prod_{\alpha\in\Phi_a^+}(1-e^{-\alpha}),
 \qquad
 \mathcal J_{\rho,a}(t)=\frac{\Delta_a(\theta(t))}{\delta_a(t)}.     \tag{3}
\]

For a simply connected group rho_a is a weight: its simple-coroot
pairings are one on the simple roots of that length and zero on the
others. The numerator and denominator have the same Weyl sign character.
The numerator vanishes on every root hypertorus in that orbit, because
theta is W-equivariant. Thus their quotient is a regular invariant
Laurent polynomial over Q. Its value at 1 is one by the identity first
jet. Write J_rho for the product over root lengths.

## 3. A formula for every parabolic

Let N_a=|Phi_a^+| and N_(L,a)=|Phi_L intersect Phi_a^+|. Put

\[
               m_{P,a}=\frac n2\left(1-\frac{N_{L,a}}{N_a}\right).
                                                                    \tag{4}
\]

These are nonnegative integers: m_(P,a) is half the number of cosets
in W/W_L moved by a reflection of length a.

**Theorem.** The map (2) is injective, its image is finite free of rank
n over R, and for its torsion cokernel C the zeroth Fitting ideal is

\[
             \boxed{\operatorname{Fitt}_0(C)=
                     \left(\prod_a\mathcal J_{\rho,a}^{m_{P,a}}\right).}
                                                                    \tag{5}
\]

The generator is specified up to a unit of R. Consequently the map is
an isomorphism after inverting the displayed root quotients. This is
a statement about the determinant ideal, not all elementary divisors
over the non-PID ring R.

### Proof

The rational invariant ring A is free over B with a homogeneous basis
f_1,...,f_n lifting any rational cohomology basis of G/P. The fixed-point
evaluation matrix of (2) has entries f_j(w^(-1)theta(t)), w in W/W_L.
Whenever a reflection fixes theta(t), each moved pair of cosets gives
one coincident pair of rows. The polynomial evaluation determinant
therefore contains alpha(theta(t))^(m_(P,a)) for each positive root
alpha of length a. Its total degree is

\[
 \sum_j\deg f_j=n\dim(G/P)/2=\sum_a N_a m_{P,a}.
\]

The generic orbit is free on W/W_L, so that determinant is nonzero.
Degree comparison makes it a nonzero rational constant times the product
of precisely those root powers. In the target, the Schubert restriction
matrix is triangular, with tangent Euler factors on its diagonal. Its
determinant is, up to an R-unit, the product of delta_a^(m_(P,a)).
Taking their ratio proves (5). This determinant is nonzero, proving
injectivity between the finite free modules and the localization claim.

For completeness, a coset wW_L is fixed by s_alpha exactly when
w^(-1)alpha belongs to Phi_L. Transitivity on roots of length a counts
the fixed cosets as n*N_(L,a)/N_a, proving (4). QED.

## 4. Connection to the full Cayley determinant

For the full flag all m_(B,a)=|W|/2. In the existing normalization the
root-direction part of the full Springer differential determinant is
exactly J_rho^2. Thus

\[
 D_\rho|_T=\det(d\theta|_{\mathfrak t})\,\mathcal J_\rho^2,
 \qquad
 \operatorname{Fitt}_0(C_{\rho,T}^{B})=(\mathcal J_\rho^{|W|/2}).    \tag{6}
\]

See the already established [root-factor identity](SPRINGER_GEOMETRY.md)
and [Kostant--Michor's generalized Cayley framework](https://arxiv.org/abs/math/0109066).
Formula (6) says exactly which part of that determinant the equivariant
K-image detects. The Cartan Jacobian is separate. A zero of a root
quotient is not automatically the entire critical divisor or a target
branch locus.

In the following simply laced cases, (4) gives

| Partial flag | n | dim(G/P) | m_P | Determinant ideal |
| --- | ---: | ---: | ---: | --- |
| E6/P1 | 27 | 16 | 6 | (J_rho^6) |
| E7/P7 | 56 | 27 | 12 | (J_rho^12) |
| E8/P8 | 240 | 57 | 57 | (J_rho^57) |

E8/P8 is the adjoint variety, not minuscule; its row uses the same
simply laced formula without that extra hypothesis. For SL_m/P_k the
exponent is binomial(m-2,k-1). For G2/P1 the short/long exponents are
(3,2), and for G2/P2 they are (2,3). Thus the formula is uniform and
does not depend on obtaining a new giant Schubert matrix in each case.

## 5. A first integral example and a representation test

For SL2 with the adjoint representation, the literal source is
Z[4omega]. Put s=t+t^(-1), let L be the right tautological character,
and use K_T(P1)=Z[t^+-1,L^+-1]/((L-t)(L-t^(-1))). The source generator
maps to b=L^2-L^(-2)=2sL-s^2 and b^2 is a coefficient-ring element.
The image has basis 1,b and its ambient quotient is

\[
              \mathbb Z[t^{\pm1}]/(2(t+t^{-1})).                   \tag{7}
\]

Here J_adj=s/2. Over Q the horizontal divisor is s=0; integrally (7)
also retains the factor two. Setting t=1 recovers the old factor four.
This example separates rational support from integral source scaling.

Another useful uniform test is rho tensor rho. Its normalized trace
formula is theta_(rho tensor rho)=(chi_rho/dim rho)*theta_rho. Hence
J_(rho tensor rho,a)=(chi_rho/dim rho)^(N_a)*J_(rho,a). The normalized
cotangent data, and therefore the nonequivariant rational K-coordinate
maps, are unchanged. Equivariant images can nevertheless acquire the
additional character divisor. This connects the
[primitive K-data](OTHER_THEORIES_COTANGENT.md) with information that is
lost under augmentation.

## 6. What remains

The [independent exact pilot](experiments/other_theories_equivariant_k_pilot.py)
checks minimal/adjoint full-flag determinants in A2, C2 and G2 at three
rational torus points each. One method uses Weyl fixed-point evaluation
and a Steinberg character basis; the other uses Laurent Demazure Euler
operators. Their determinant ratios agree with (5). A separate small
weight-orbit enumeration checks the reflection exponents for both G2
maximal flags, the F4/P1 and P4 pair, and E6/P1, E7/P7 and E8/P8. It
does not construct the enormous exceptional Weyl groups or claim to
check their entire K-image matrices. [Report](experiments/other-theories-equivariant-k-pilot.json).

The source family (1), determinant ideal (5) and root-sector interpretation
are characteristic-zero statements proved above. Integral source freeness,
higher Fitting ideals, faithful specialization and a full API remain
separate questions. The existing E6/E7 quotient coordinates and root
determinants are direct inputs to this project, rather than unrelated
background. A targeted computation should begin with (7), G2 and one
minuscule exceptional flag, retaining both arithmetic and horizontal
degeneration data.

## 7. An integral symplectic projective family

There is also a complete integral model in arbitrary rank. Take
G=Sp_(2r), rho defining, P=P1, and use its literal source
Z[z,q_1,...,q_(r-1)] from the [height note](OTHER_THEORIES_HEIGHT.md).
Put R=Z[t_1^+-1,...,t_r^+-1], s_i=t_i+t_i^(-1), and

\[
 K_T(P^{2r-1})=R[L]/\prod_{i=1}^r(L-t_i)(L-t_i^{-1}),
 \quad b=L-L^{-1},\quad B_i=t_i-t_i^{-1},\quad E_j=e_j(B_1^2,...,B_r^2).
\]

The full image is exactly R[b], free with basis 1,b,...,b^(2r-1).
The complete source kernel is equivalently given by
q_0=1, q_r=0 and

\[
              q_j+z^2q_{j-1}=E_j\quad(1\le j\le r).
\]

Eliminating the q_j leaves the monic relation
prod_i(z^2-B_i^2)=0. These relations hold because the left sides are
the global signed-symmetric Springer characters. Conversely the quotient
is free of rank 2r, and its generic fixed-point values are the 2r distinct
numbers plus/minus B_i. Thus it embeds and the kernel is complete.

**Integral determinant formula.** For the actual embedded image,

\[
 \boxed{\operatorname{Fitt}_0(K_T/R[b])=
           \left(2^r\prod_{i<j}(s_i+s_j)^2\right).}               \tag{8}
\]

Indeed take u_1,...,u_(2r)=t_1,t_1^(-1),...,t_r,t_r^(-1). Comparing
the Vandermonde matrices for L and b gives the determinant
prod_(a<b)(1+1/(u_a u_b)). Each reciprocal pair contributes two.
The four cross pairs for i<j contribute (s_i+s_j)^2, proving (8).

At t_i=1 this becomes 2^(r(2r-1)), the old total index. The source
image specializes faithfully to Z[b] in ordinary K-theory: after
specialization its square inclusion matrix still has nonzero integer
determinant, so its kernel, equivalently the Tor_1 specialization defect,
is zero. This does not make its Smith factors the ordinary cohomology
ones. They are computed, in part uniformly, in the height note.

As a check on the rational theorem, the defining symplectic root
quotients are J_long=1 and
J_short=4^(-binomial(r,2))*prod_(i<j)(s_i+s_j). The parabolic exponents
are m_long=1 and m_short=2 for r>=2. Multiplying J_short^2 by the
literal-source index recovers (8). This family therefore retains both
the horizontal character divisor and the vertical integer factors.
The [independent symbolic pilot](experiments/other_theories_symplectic_equivariant_pilot.py)
constructs the actual projective-bundle K multiplication matrices for
r=1,2,3, verifies (8) as a polynomial identity and computes the full
characteristic polynomial of b without using root factors.
[Report](experiments/other-theories-symplectic-equivariant-pilot.json).

## 8. Two exact local Smith patterns in characteristic zero

The determinant does not by itself determine higher Fitting ideals.
There are, however, two structural local patterns. Work over the DVR
R_D at a height-one prime D of Q[X], with uniformizer pi. These are
horizontal local factors in characteristic zero, not the integer-prime
Smith factors of ordinary cohomology.

**A single reflection wall.** Suppose t is regular modulo D, precisely
one positive root alpha of theta(t) vanishes modulo D, and it has
valuation e. Then the local cokernel of (2) is

\[
                        (R_D/(\pi^e))^{m_{P,a}},                  \tag{9}
\]

where a is the length of alpha. To prove this, pass to a strict henselian
completion, which preserves the Smith valuations. The only stabilizer
of the residual theta is the order-two reflection. Its fixed cosets
give unramified rank-one factors; each moved pair gives the rank-two
reflection slice. Since two is a unit, this slice is
R_D[z]/(z^2-c^2) -> R_D direct-sum R_D, z -> (c,-c), with v_D(c)=e.
Its matrix has Smith factors (1,pi^e). There are exactly m_(P,a)
moved pairs. This is the usual tame reflection local model applied to
the specified map; it requires all the stated regularity hypotheses.
One can obtain that model directly by completing the finite Weyl
quotient: the completed ring at an orbit is the completed ring at one
point invariant under its stabilizer, since averaging commutes with
completion in characteristic zero. A reflection negates one coordinate
and fixes the others. Partial-quotient orbit factors are consequently
rank one or the displayed rank-two model.

**Scalar collapse.** Suppose theta_rho=a(t)*theta_sigma for an invariant
Laurent function a, v_D(a)=e>0, and all the root quotients J_(sigma,b)
are units at D. Then the complete local Smith valuations are

\[
                    \boxed{\{e\ell(w):w\in W/W_L\},}             \tag{10}
\]

with minimal coset lengths and their usual Betti multiplicities. Indeed
choose the homogeneous source basis f_w of degree ell(w). Its image
matrix is exactly M_sigma*diag(a^(ell(w))). The reference matrix is
invertible at D by (5), so elementary operations give (10). No Schubert
matrix of exceptional size is required for this assertion.

For defining SL_m, J_sigma=1, and tensor squaring gives
a=chi_V/m. Hence (10) applies unconditionally at each characteristic-zero
component of the character divisor chi_V=0, with its actual valuation e,
for every parabolic. On SL3/P1 the pattern is (1,pi^e,pi^(2e)). On the
full SL3 flag it is (1,pi^e,pi^e,pi^(2e),pi^(2e),pi^(3e)). This shows
why a general claim that every horizontal factor occurs to exponent one
would be false, even though it holds at a transverse single-wall point.

In fact this family has a **global** answer over R=Q[X]. For the defining
SL_m module V and any tensor power V^(tensor k), k>=1,

\[
 \boxed{C_{V^{\otimes k},T}^P\cong
       \bigoplus_{w\in W/W_L}R/(\chi_V^{(k-1)\ell(w)}).}           \tag{11}
\]

Terms of exponent zero are zero modules. To prove this, the defining
Springer map has J_V=1, so its equivariant matrix is invertible over R
globally. Equivalently it is the translation of the right Chern roots
L_i by the coefficient chi_V(t)/m in the type A flag presentation.
Tensor differentiation gives
theta_(V tensor k)=(chi_V/m)^(k-1)*theta_V. The preceding column-scaling
argument now gives a global diagonal presentation, proving (11).

If d_1<=...<=d_n are the minimal coset lengths, its higher Fitting
ideals are explicitly
Fitt_i(C)=(chi_V^((k-1)*sum_(j=1)^(n-i) d_j)) for 0<=i<n;
Fitt_i is the unit ideal for i>=n. This is a special diagonal
presentation over a non-PID, not a claim of general Smith normal form
over R. It is a rational coefficient theorem: actual integral source
scaling and integer-prime factors are separate.
The [exact Laurent/Demazure check](experiments/check_other_theories_scalar.py)
verifies the global column identity and unit reference matrix on the
full SL3 flag. [Report](experiments/other-theories-scalar-check.json).

The reflection slice is established local invariant theory; the scalar
statement is an elementary homogeneous-basis deduction. Their image
interpretation is recorded here without a priority claim. Higher-rank
stabilizers, multiple walls, common factors of a and the reference root
quotients, integral coefficient primes and global higher Fitting ideals
still require separate work.
