# E7/P7: an explicit basis for the 2-primary comparison

Status, 20 October 2026: complete computer-assisted deduction. The
primary calculation and all three independent checks have passed.
The basis argument and its finite verification are given below;
no claim of literature priority is made. See the
[proof/source audit](E7_ROOT_WEIGHT_AUDIT.md).

## 1. Maps, lattices and the exact statement

Let G be simply connected complex E7. Let X be its weight lattice,
Q its root lattice and K={1,2,3,4,5,6}. The retained-node parabolic
P_K is the maximal parabolic P7; its derived Levi is E6. Use the
Bourbaki/Sage edges (1,3),(3,4),(2,4),(4,5),(5,6),(6,7).
For any retained set I put

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dQ)^{W_I}),\qquad
 N_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dX)^{W_I}).
\]

These are the full, unscaled characteristic images, as additive lattices
inside integral cohomology. They are not selected-generator algebras.
Since X/Q=Z/2, the finite group N_I^d/M_I^d is 2-primary. We can therefore
work over R=Z_(2) without losing any relative factors. Polynomial degree
d means cohomological degree 2d.

Let h,s,t be the published Schubert generators on G/P7 of degrees 1,5,9:
h=sigma_[7], s=sigma_[2,4,5,6,7], t=sigma_[1,5,4,2,3,4,5,6,7].
Put S=2s, T=2t and epsilon_i=1 for even i, 2 for odd i.

**Theorem.** The following are respective R-bases of the two full images:

\[
\begin{array}{ll}
 N_K\otimes R:&
 h^i\ (0\le i\le27),\quad h^iS\ (0\le i\le13),\quad
 h^iT\ (0\le i\le13),\\[2mm]
 M_K\otimes R:&
 \epsilon_i h^i\ (0\le i\le27),\quad
 \epsilon_i h^iS\ (0\le i\le13),\quad
 \epsilon_i h^iT\ (0\le i\le13).
\end{array}
\tag{1}
\]

In particular, only the exponent of h controls the relative multiplier;
it is not the parity of the total degree. Define [r]_(q^2)=sum_(j=0)^(r-1) q^(2j).
Then in every degree

\[
 \boxed{N_K^d/M_K^d\cong(\mathbb Z/2)^{m_d},\qquad
 \sum_d m_d q^d=q[14]_{q^2}+(q^6+q^{10})[7]_{q^2}.}
\tag{2}
\]

There are no relative factors 4 or 8. The total additive quotient has
28 factors 2. For example m5=m6=1, m10=m14=2, m24=m26=0 and m27=1.
Both images and all multiplicities vanish above degree 27.
The displayed bases are bases of N and M at 2, not integral Schubert
bases of the whole cohomology ring.

An additional consequence, proved in Section 6 without a degreewise
computation, is

\[
 \boxed{2N_I^d\subseteq M_I^d\quad\text{for every }I\subseteq K
 \text{ and every }d.}
\tag{3}
\]

Thus all 64 standard parabolics refining P7 have elementary 2-primary
relative quotients. Formula (2) gives the multiplicities only for P7;
it does not determine those of the other 63 refinements, or address
the other proper retained sets containing node 7.

## 2. Published cohomology input

[Duan--Zhao, Theorem 5](https://arxiv.org/abs/math/0511332) gives

\[
 H^*(G/P_7,\mathbb Z)=\mathbb Z[h,s,t]/
 (s^2-2ht,\ 2st-9h^4s^2+6h^9s-h^{14},\
 t^2+10h^3s^3-9h^8s^2+2h^{13}s).
\tag{4}
\]

Its dimension is 27 and its Schubert rank is 56. This ring, these
generators, and their Schubert meanings are published results. The
Weyl-group Poincare polynomial is

\[
 P_{W(E_7)/W(E_6)}(q)
 =\frac{[10]_q[14]_q[18]_q}{[5]_q[9]_q}
 =[28]_q+(q^5+q^9)[14]_q.
\tag{5}
\]

The rational invariant generators omega7, psi2, psi5, psi6, psi8, psi9,
psi12 and their Borel images use the unchanged convention of
[Kumar--Xie, Part III, Theorem 8.2(g)](https://kumar.math.unc.edu/papers/kumar86.pdf).
In particular, the inputs needed for the explicit lifts below are

\[
\begin{aligned}
 \beta(\omega_7)&=h,&\beta(\psi_2)&=-144h^2,\\
 \beta(\psi_5/5760)&=70h^5-162s,\\
 \beta(\psi_9/387072)&=-3332h^9-486h^4s+39366t.
\end{aligned}
\tag{6}
\]

The source forms and transport through the normalized Springer maps
come from the [earlier source proofs](E7_REPRESENTATION_COMPARISON.md).
Formula (2) concerns their integral lattice comparison; the complex
surjectivity statements in the Kumar papers do not by themselves
determine these finite groups.

## 3. The basis reduction in the published ring

Substitution S=2s, T=2t into (4) gives

\[
\begin{aligned}
 S^2&=4hT,\\
 ST&=2h^{14}-6h^9S+18h^5T,\\
 T^2&=-40h^{18}+116h^{13}S-324h^9T.
\end{aligned}
\tag{7}
\]

Comparing S(ST) with S^2T gives
98h19-285h14S+798h10T=0. Comparing (ST)T with ST^2, and then
eliminating h18S, gives the following useful relations over R:

\[
 h^{14}S=\frac{98}{285}h^{19}+\frac{14}{5}h^{10}T,\qquad
 h^{14}T=\frac{238}{1311}h^{23},\qquad h^{28}=0.
\tag{8}
\]

All denominators are odd. The last identity also follows from the top
degree. The [plain Python/SymPy check](experiments/e7-two-primary-ring-check.json)
checks all six identities directly in (4).

Let L_N be the R-span of the 56 proposed N-basis elements in (1).
Relations (7) reduce products involving S and T to three towers.
Relations (8) reduce the excess S and T exponents, and h28 kills the
remaining excess h exponents. Hence L_N is exactly R[h,S,T]. Over Q,
these generators generate all of (4). Its dimension is 56, so the 56
spanning elements are independent. This proves an explicit R-basis of
the candidate image algebra, without assuming an integral Schubert basis.

Let L_M be the span of the proposed M-basis. Every displayed element is
a monomial in 2h,h2,S,T. Conversely this span is closed under those four
generators. Multiplication by h2 preserves the parity of the h exponent;
multiplication by 2h supplies the required factor 2 when parity changes.
In (7), every term changing that parity has an even coefficient.
The coefficients in (8) that change parity also have even numerators
and odd denominators. Thus the same property holds after reduction.
This proves

\[
 L_N=R[h,2s,2t],\qquad L_M=R[2h,h^2,2s,2t].
\tag{9}
\]

This is a ring calculation from a published presentation. The assertion
that these two rings are the *full* characteristic images still needs
the next two steps.

## 4. Explicit source lifts: inclusion of the candidate rings

The integral permutation coordinates v1,...,v6,u for Q from
[the source proof](E7_CHARACTER_FORM.md) give
Sym_Z(Q)^(W_A5)=Z[E1,...,E6,u], where Ej=e_j(v1,...,v6).
Their change of coordinates has determinant -1. The weight-lattice
version replaces u by E1+2z and has determinant -1 in fundamental weights.
Thus coefficients in these symmetric coordinates detect actual root or
weight integrality, including at 2.

The three source polynomials

\[
 B=-\psi_2/144,\qquad F_5=\psi_5/5760,\qquad
 F_9=\psi_9/387072
\]

have integer root coefficients; their compressed coefficient lists are
reconstructed and checked by both implementations. B also belongs to
Sym_Z(V)^(W_E6), where V is the span of alpha1,...,alpha6.
For example B is an integer quadratic in these six simple roots and
has no alpha7 term. It maps to h2 by (6).

Put

\[
\begin{aligned}
 A&=2\omega_7,\\
 C&=(35AB^2-F_5)/81,\\
 D&=(F_9+1666AB^4+243B^2C)/19683.
\end{aligned}
\tag{10}
\]

A is an E7 root-lattice element, and all the displayed denominators
in C,D are odd. Therefore A,B,C,D are integral root-invariant
polynomials over R. Equation (6) proves their Borel images are
2h,h2,2s,2t. This gives L_M contained in M_K tensor R. The weight
source additionally contains omega7, proving L_N contained in N_K tensor R.

## 5. Full-source containment: the computer-assisted step

We must exclude extra image elements coming from invariant polynomial
combinations, especially combinations involving the rational Borel
kernel. Checking the four lifts alone would not do this.

In each degree d<=27, take **every** weighted monomial in
omega7,psi2,psi5,psi6,psi8,psi9,psi12, with weights
(1,2,5,6,8,9,12). They are a basis of the rational invariant source.
There are 1,967 such monomials over all these degrees.
Expand them in the root or weight permutation coordinates, giving the
injective rational coefficient matrix C_d.

The historical full-source bases are used as proposed bases B_d after
removing 12^d or 36^d. They are not trusted because an earlier Smith
table passed. For both sources in every degree, the new computation checks

\[
 C_dB_d\text{ is R-integral},\qquad
 \operatorname{rank}_{\mathbb F_2}(C_dB_d\bmod2)
 =\text{the full source rank}.
\tag{11}
\]

These conditions prove that B_d is the entire local source lattice:
if an allowed rational vector B_d u had a coefficient with a least
2-denominator, multiply u by that denominator and reduce modulo 2.
One obtains a nonzero kernel vector of C_dB_d modulo 2, a contradiction.
The A5 symmetric-coordinate model makes coefficient integrality equivalent
to actual root or weight integrality. No averaging by a Weyl order is used.

Only after this certification do we apply beta and express every image
column in the explicit candidate basis. All coordinates are R-integral,
and the resulting columns have full row rank modulo 2. Thus their R-span
equals the candidate module. This proves the reverse inclusions and
completes (1). In particular all source-kernel combinations are retained
before integrality and before the image is taken.

The primary [report](experiments/e7-two-primary-basis.json) records all 28
degrees for both sources and 392 module-closure checks. The
[independent replay](experiments/e7-two-primary-basis-audit.json) reconstructs
every coefficient condition and image column using different formulas
and a different rational normal form. A [separate Schubert check](experiments/e7-two-primary-schubert-check.json)
uses fresh BGG calculations for all seven generators and Chevalley
multiplication for every proposed basis class. See the audit for the
shared arithmetic dependencies and limits of independence.

Since the target vanishes above degree 27, this is an exhaustive finite
proof of the full P7 image identification, not a bounded-degree conjecture.
Finally, the simultaneous bases immediately give (2): count odd exponents
i in each of the three towers. No Smith reduction is needed for that step.

## 6. A structural exponent bound on all refinements of P7

Let V_R be the R-span of the first six simple roots and let w=omega7.
The actual lattice computations give

\[
 X_R=V_R\oplus Rw,\qquad Q_R=V_R\oplus R(2w).
\tag{12}
\]

The basis matrices have index 3 in X and Q respectively, a unit at 2.
These are local identities, not integral decompositions at 3. Every
W_I with I contained in K preserves V_R and fixes w. Consequently, with
A_I=(Sym_R V_R)^(W_I), taking invariants coefficient by coefficient gives

\[
 (\operatorname{Sym}_R X_R)^{W_I}=A_I[w],\qquad
 (\operatorname{Sym}_R Q_R)^{W_I}=A_I[2w].
\tag{13}
\]

This does not require the E6 torsion index to be invertible. The root
polynomial B in Section 4 is invariant under W_K and maps to h2, so its
pullback proves h2 is in M_I for every refinement. Reducing powers of h
therefore gives N_I=M_I+hM_I over R, while 2h is already in M_I. Thus
2N_I is contained in M_I over R. The integral quotient is 2-primary, so
the same inclusion holds for the original integer images, proving (3).

The quadratic lift can also be proved directly from classical invariant
forms, independently of the numerical catalogue. Normalize root lengths
to 2, put a=(alpha1,...,alpha6)^t and write C6 for the E6 Cartan matrix.
In the orthogonal basis consisting of these roots and w, the E7 invariant
quadratic is

\[
 F_7=F+\frac23w^2,\qquad F=a^t C_6^{-1}a.
\]

The matrix 3C6^(-1) is integral with even diagonal. Consequently
B=-(3/2)F is an integer root polynomial invariant under W_E6.
The rational Borel map kills the positive-degree full Weyl invariant F7,
so beta(B)=h2. Direct expansion also identifies this B with -psi2/144.
The plain Python check verifies both quadratic identities and the matrix
integrality. Thus (3) has a hand proof using Cartan data and the classical
Borel relation; it does not depend on the high-degree source certificates.

More generally, suppose at a prime p the actual source lattices split
as X_R=V_R plus Rw and Q_R=V_R plus R(p^e w), with a common invariant
complement and w fixed by the Levi Weyl group. If h=beta(w) satisfies
h^r in M for a positive integer r, then polynomial reduction gives

\[
 N=\sum_{j=0}^{r-1}h^jM,\qquad p^{e(r-1)}N\subseteq M.
\tag{14}
\]

Indeed p^(ej)h^j=(p^e h)^j belongs to M. This is an elementary adjoining-
generator lemma, not a new Leray--Hirsch theorem. Its hypotheses include
the actual stable lattice splitting and the root-image lift of h^r;
an index condition on its own is insufficient. The E7 application has
e=1,r=2. It gives an exponent bound, not the finer-flag multiplicities.

The [refinement extension](E7_REFINEMENT_MULTIPLICITIES.md) now supplies
a uniform multiplicity formula through degree nine for all 64 cases.
It also reduces the remaining all-degree task to three base polynomials
and exhibits an obstruction to a polynomial common-source-basis shortcut.

## 7. Connection to the Springer maps and attribution

For the normalized minimal and adjoint Springer morphisms, the earlier
source theorem gives J56,I^d=12^d M_I^d and J133,I^d=36^d N_I^d.
Their previously proved sum/intersection identities imply
(J56+J133)/J56 is isomorphic, as an additive finite group, to N/M.
Thus (2) gives its exact factors on P7, and (3) bounds it by exponent 2
on all refinements. It does not compute an inverse Springer map, a
generic map degree, or either image's whole ambient cohomology cokernel.

Published inputs are the characteristic-map framework, the rational
Levi invariant generators and their images, the integral cohomology
presentation (4), Weyl-group degrees and Schubert multiplication.
[Totaro, Theorem 6.1](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf)
already gives adjoint E7 torsion index 24; with divided differences and
X/Q=Z/2 it supplies the general relative exponent bound 8.
Our sharper bound 2 under (12)--(13), the explicit simultaneous bases,
and formula (2) are deductions recorded here. The factor lists agree
with the earlier finite E7/P7 table; they are not new numerical data.
The additional work is their basis explanation and the structural
refinement consequence. No claim that the cited papers state these exact
conclusions, or that these deductions are original, is made.

## Reproduction

From borel/:

```sh
PYTHONPATH=src sage -python experiments/e7_two_primary_basis.py
sage -python experiments/audit_e7_two_primary_basis.py
PYTHONPATH=src sage -python experiments/e7_two_primary_schubert_check.py
python experiments/e7_two_primary_ring_check.py
```

The source proposal files are preserved unchanged. Runtime and the
installed library API remain version 0.20.0; these are research proof
experiments with hash-bound reports.
