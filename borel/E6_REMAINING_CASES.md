# The three exceptional E6 root-to-weight comparisons

Status, 20 October 2026: complete computer-assisted deduction. Both
base calculations and their independent all-degree certificate replays
have passed: 22 degrees for A5 and 30 for A2+A2+A1. The third case is
proved by the common-basis degree-shift argument in Section 6.
This note extends the [60-parabolic theorem](E6_ROOT_WEIGHT_COMPARISON.md).
Polynomial degree d means cohomological degree 2d; coefficients are
R=Z_(3). Literature priority is separate from validity.

## 1. Precisely what is being completed

G is simply connected complex E6, Q is its root lattice, X its weight
lattice, and I is a retained-node set. Keep the same maps and integral
images as in the preceding theorem:

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dQ)^{W_I}),\qquad
 N_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dX)^{W_I}).
\]

The three retained sets are

\[
 I_A=\{1,3,4,5,6\},\quad
 I_B=\{1,2,3,5,6\},\quad I_C=\{1,3,5,6\}.
\]

Their derived Levi types are A5, A2+A2+A1 and A2+A2. The respective
dimensions of G/P_I are 21,29,30, with total Schubert ranks 72,720,1440.
The comparison quotients are finite and 3-primary; localization loses
none of them. We compare full invariant sources before applying beta,
including every polynomial in the rational Borel kernel.

## 2. Multiplicity statement

**Theorem.** For the three retained sets in Section 1, write
T_(I,a)(q)=sum_d m_(I,a,d) q^d. In every degree,

\[
 N_I^d/M_I^d\cong\bigoplus_{a=1}^3
 (\mathbb Z/3^a)^{m_{I,a,d}}.
\tag{1}
\]

For I_A the exact polynomials are

\[
\begin{aligned}
 T_{I_A,1}&=\sum_{d=10}^{15}q^d+2q^{16}+2q^{17}
                       +\sum_{d=18}^{21}q^d,\\
 T_{I_A,2}&=q^5+2q^6+3\sum_{d=7}^{14}q^d
                       +2q^{15}+q^{16}+q^{17}+q^{18},\\
 T_{I_A,3}&=q^3+q^4.
\end{aligned}
\tag{2}
\]

For I_B the coefficients m_(I_B,a,d) are specified degree by degree in
the [complete tables](E6_REMAINING_TABLES.md). The third case follows exactly from

\[
 \boxed{T_{I_C,a}(q)=(1+q)T_{I_B,a}(q).}
\tag{3}
\]

Section 6 proves (3) independently of any computed table. The unit
factors in each degree are the ambient Betti rank minus the three
displayed multiplicities.

Together with the 60-case theorem, (1)--(3) specify the relative
comparison on every proper standard E6 parabolic. This is an all-degree
computer-assisted theorem for a fixed root system, not a formula for
arbitrary exceptional types. The G/G comparison is zero.

For example, in polynomial degree 3 the A5 quotient is Z/27 and the
A2+A2+A1 quotient is (Z/3)^2 plus Z/27. The A2+A2 quotient in the same
degree is (Z/3)^3 plus Z/27. These describe the difference between M
and N, not either image's cokernel in the whole cohomology group.

The [stepwise proof/source audit](E6_REMAINING_AUDIT.md) distinguishes
published ingredients, exact finite certificates and deductions.

## 3. Full-rank A5+A1 coordinates

Let alpha_i and omega_i use the previous E6 numbering. For I_A take
the A5 chain (alpha1,alpha3,alpha4,alpha5,alpha6) and z=omega2,
the highest E6 root. For I_B and I_C use the A5 chain

\[
 (\alpha_1,\alpha_3,\gamma,\alpha_5,\alpha_6),\qquad
 \gamma=-\alpha_1-\alpha_2-2\alpha_3-2\alpha_4-2\alpha_5-\alpha_6,
 \quad z=\alpha_2.
\]

The vector gamma has fundamental-weight coordinates (0,0,-1,1,-1,0).
It is an actual E6 root. In both cases z is a root perpendicular to
the A5 chain. Define x1,...,x6 by x_i-x_(i+1) equal to the five chain
roots and sum x_i=0. Direct lattice determinants give

\[
 Q_R=Q_{A_5,R}\oplus Rz,\qquad
 X_R=X_{A_5,R}\oplus Rz.
\tag{4}
\]

These identities are 3-local: the discarded index is a power of 2.
Root differences r_i=x_i-x6, i=1,...,5, form a basis of the A5 root
lattice; x1,...,x5 form its weight-lattice basis. Membership and unit
determinants are checked in the actual E6 Q and X, not only over Q.

This is the structural reason for using a full-rank subsystem: at 3
the exceptional root/weight defect becomes the classical A5 isogeny
defect, with a common independent A1 coordinate. The E6 geometry is
then imposed through its rational Borel relations. To reuse this method
in another type one must check both actual lattices and the parabolic
action; the existence of a full-rank subsystem alone is insufficient.

The Weyl actions are S6 for I_A, S3 on each triple together with
z -> -z for I_B, and just the two S3 groups for I_C. Therefore the
complete weight-invariant sources are

\[
\begin{array}{ll}
 I_A:&R[z,c_2,c_3,c_4,c_5,c_6],\quad c_j=e_j(x_1,\ldots,x_6),\\
 I_B:&R[a_1,a_2,a_3,b_2,b_3,z^2],\\
 I_C:&R[a_1,a_2,a_3,b_2,b_3,z],
\end{array}
\tag{5}
\]

where a_j and b_j are elementary functions of the first and second
triples and b1=-a1. This is integral symmetric-polynomial theory: start
with independent permutation variables and quotient by sum x_i=0.
Surjectivity on invariants follows from H^1 of each permutation-monomial
module being zero over the torsion-free ring R (Shapiro reduces it to
Hom of a finite stabilizer into R). The relation is invariant and a
non-zero-divisor. Thus no additional integral weight invariants are missed.

## 4. Exact rational Borel presentation

The 27 minimal weights, in either coordinate system, are

\[
 \{x_i+x_j:i<j\}\ \cup\ \{-x_i+z/2,-x_i-z/2:1\le i\le6\},
\tag{6}
\]

up to replacing every weight by its negative. Equality with the
27-element E6 reflection orbit is checked exactly. Let f_k be the sum
of their kth powers, for k=2,5,6,8,9,12. These are W_E6-invariant.
The quotient by these six polynomials is finite; its dimension in the
I_A invariant ring is 72. Hence the f_k are algebraically independent.
Their degrees are exactly the published fundamental invariant degrees
of E6, so equality of Hilbert series proves that they generate the
rational Weyl-invariant ring. The classical Borel theorem then gives
the desired rational characteristic-map kernel.

This argument selects and verifies a convenient basic invariant system;
it does not assert a new Borel presentation theorem. The trace-power
description is also independent of which Springer representation is
being compared: it is used here only to describe the rational Borel kernel.

In particular f2=3(z^2-4c2). After eliminating z^2, the full-rank
A5+A1 invariant quotient is a 36-dimensional algebra in c2,...,c6,
with relations f5,f6,f8,f9,f12. For I_A adjoining z with z^2=4c2
gives dimension 72 and top degree 21.

For I_B use the universal Gr(3,6) relations. With c_j the total
elementary functions, put

\[
 b_2=c_2-a_2+a_1^2,\qquad
 b_3=c_3-a_3-a_1c_2+2a_1a_2-a_1^3.
\]

The three additional relations are

\[
 a_1b_3+a_2b_2-a_1a_3=c_4,\quad
 a_2b_3+a_3b_2=c_5,\quad a_3b_3=c_6.
\tag{7}
\]

Their relative Groebner basis is monic, with leading monomials all
degree-four monomials in (a1,a2,a3), where degree here counts factors.
It follows that a1^i a2^j a3^k, i+j+k<=3, give a 20-element basis
over the total symmetric ring. The relative leading monomials and the
36-dimensional base leading monomials involve disjoint variables.
Their union is a Groebner basis by the coprime-leading-term criterion.
This proves dimension 720 and top degree 29, without an expensive
six-variable Groebner calculation or a guessed degree cutoff.

Equations (7) and the symmetric/Schur basis are classical Grassmannian
identities. Their use here is an exact computational reduction.

## 5. Full root sources and finite proof certificates

For I_A, work in the S5-invariant polynomial ring on the r_i. If
E_j=e_j(r1,...,r5), E0=1, E6=0 and mu=E1/6, then

\[
 c_k=\sum_{j=0}^k\binom{6-j}{k-j}(-\mu)^{k-j}E_j.
\tag{8}
\]

For I_B and I_C, the subgroup S3 x S2 permutes r1,r2,r3 and r4,r5.
Write Uj=e_j(r1,r2,r3), Vj=e_j(r4,r5), mu=(U1+V1)/6. The five
weight-source generators become

\[
\begin{aligned}
 a_1&=(U_1-V_1)/2,\\
 a_2&=U_2-2\mu U_1+3\mu^2,&
 a_3&=U_3-\mu U_2+\mu^2U_1-\mu^3,\\
 b_2&=V_2-2\mu V_1+3\mu^2,&
 b_3&=-\mu V_2+\mu^2V_1-\mu^3.
\end{aligned}
\tag{9}
\]

These smaller permutation-invariant rings detect root integrality
exactly. A polynomial already invariant under the full Levi belongs
to Sym_R(Q) if and only if its coefficients in E, or in U,V, belong
to R. This avoids expanding in five individual difference variables.

In each degree take every weighted monomial in the rational generators
in (5), retaining the kernel. Let C be its coefficient matrix after
(8) or (9), and let B have columns spanning

\[
 \{v\in\mathbb Q^n:Cv\in R^m\}.
\tag{10}
\]

The computer certificate checks CB is R-integral and its reduction
modulo 3 has column rank n. These two properties *prove completeness*:
if v=Bu had a coefficient with a 3-denominator, clearing its least
denominator and reducing modulo 3 would give a nonzero vector in the
kernel of CB. Thus every allowed root invariant occurs, not only a
convenient denominator-cleared subset.

The source lattices are calculated by modular Hermite reduction of
the dual coefficient rows. Adding 3^k times the identity imposes
weight-integral coefficients; reduction modulo 3^k preserves all local
conditions. The independent primitive-column check makes the result
independent of this particular choice of lattice algorithm.

The z coordinate is independent of the A5 lattice. Consequently the
degree-d source is the direct sum of these base coefficient lattices
multiplied by z^k, or by z^(2k). Their rational Borel images are formed
only *after* the complete sources have been certified. In each target
degree, image-basis certificates again check local membership and full
rank modulo 3. Integer Smith factors then give (1).

All certificates passed, giving an exact finite, computer-assisted
proof. The rational presentations prove that target degrees above
21 or 29 vanish, so the checks exhaust the two base spaces.

For a hand-checkable example, consider I_A in degree 3. Formula (8) gives

\[
 c_2=E_2-\frac5{12}E_1^2,\qquad
 c_3=E_3-\frac23E_1E_2+\frac5{27}E_1^3.
\]

The weight source has basis (z^3,zc2,c3). Since z is independent, its
root-integral sublattice has basis (z^3,3zc2,27c3): the displayed
coefficients force, and suffice for, these valuations. In the rational
target z^3=4zc2, and there are no further relations in degree 3.
Thus N has basis (zc2,c3) and M has basis (zc2,27c3), since 4 is a
unit at 3. This proves N^3/M^3=Z/27 directly. It also illustrates why
root integrality must be imposed before using the Borel relation.

The certificate argument itself works over Z_(p) for any prime: for an
injective rational coefficient matrix C and a nonsingular rational B,
integrality of CB and full column rank modulo p imply
C^(-1)(Z_(p)^m)=B Z_(p)^n. This is the standard primitive-lattice
criterion expressed in polynomial coefficient coordinates. Its use here
provides a reusable way to certify full sources without an integral
presentation of the root-invariant algebra.

## 6. Why the third case is a degree shift

The reflection s2 commutes with W_(I_C), fixes the two A2 blocks and
sends z=alpha2 to -z. At 3 both root and weight lattices have the
same direct A1 summand Rz. Thus their I_C-invariant sources decompose
as the respective I_B sources plus z times those sources.

The rational Borel relations involve only z^2. After f2, the I_C
target is the monic extension of the I_B target by z^2=4c2. It is free
with basis (1,z). The two embedded image decompositions are therefore
direct in a *common* target basis, and

\[
 (N_{I_C}/M_{I_C})^d\cong
 (N_{I_B}/M_{I_B})^d\oplus(N_{I_B}/M_{I_B})^{d-1}.
\]

Equivalently, G/P_(I_C) -> G/P_(I_B) has P1 fiber, with z restricting
to twice its divisor generator, a unit at 3. This is the familiar
projective-bundle/Leray--Hirsch argument, not a new cohomology theorem.
It proves (3) without extrapolating any finite data.

## 7. Attribution and verification scope

Published inputs include the rational Borel characteristic theorem,
E6 invariant degrees, symmetric-polynomial and Newton identities,
Grassmannian Chern relations and Schur bases, and Hermite/Smith theory.
The uniform exponent bound 27 follows from the published adjoint
E6 torsion index 54; see [Totaro, Theorem 6.1](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf).
The integral E6/P2 presentation in [Duan--Zhao, Theorem 3](https://arxiv.org/abs/math/0511332)
is established background, not a ring presentation claimed here as new.
The characteristic-map framework is explained in
[Calmes--Petrov--Zainoulline](https://arxiv.org/abs/0905.1341).

The deductions recorded here concern the *full root-to-weight image
quotients*, their exact degreewise factors and the completion of the
E6 comparison. Their interpretation as Springer-image comparisons uses
the earlier literal-source formulas J27=6^d M and J78=24^d N.
No originality or literature-priority claim follows from the calculations.

The separate audit imports neither sage_borel nor the main experiment.
Plain Python/SymPy reconstructs the trace invariants from the 27 literal
branching weights. SymPy uses a different Groebner order; a separately
written sparse reducer evaluates every image column. Newton sums replace
the main coefficient-substitution formulas. Every proposed source and
image basis is checked for membership and local primitivity. An independent
Python diagonalization modulo 81 checks the image Smith multiplicities;
the bound 27 is also checked directly on the image lattices. Sage exact
polynomial/matrix arithmetic is shared, so this is not an independent CAS
for every arithmetic operation or external peer review.

A further direct BGG calculation checks the first seven degrees of
each base case and the first six degrees of the P1 extension. It does
not recompute every high-degree Schubert column. It checks a different
mathematical algorithm from the rational trace-quotient calculation.
