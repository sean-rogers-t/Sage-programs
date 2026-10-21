# E7 full-flag comparison through an A6 Levi

Status, 20 October 2026: complete computer-assisted deduction, with all
independent checks passed. This resolves the all-degree full-flag target
from the [refinement note](E7_REFINEMENT_MULTIPLICITIES.md). Classical
common-fiber transfer gives 108 proper standard parabolics; the separate
P7 basis theorem adds one more. Literature priority is not established.

## 1. Objects and theorem

Let G be simply connected complex E7, X its weight lattice and Q its
root lattice, in Bourbaki/Sage numbering. For retained nodes I define

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dQ)^{W_I}),\qquad
 N_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dX)^{W_I}).
\]

These are the full unscaled integral characteristic images. Since X/Q
has order two, N_I/M_I is 2-primary, and R=Z_(2) retains the entire
relative group. Polynomial degree d means cohomological degree 2d.
Write [r]_q=1+q+...+q^(r-1) and

\[
 P_W(q)=[2]_q[6]_q[8]_q[10]_q[12]_q[14]_q[18]_q.
\]

**Full-flag theorem.** In every degree,

\[
 N_{\varnothing}^d/M_{\varnothing}^d
       \simeq(\mathbb Z/2)^{m_d},\qquad
 \boxed{\sum_d m_dq^d=\frac{q}{1+q}P_W(q)
       =q[6]_q[8]_q[10]_q[12]_q[14]_q[18]_q.}
 \tag{1}
\]

The first coefficients, starting at degree zero, are

\[
 0,1,6,21,56,126,252,461,786,1265,1940,2855,4054,5578.
\]

There are 1,451,520 nonunit factors in total, all equal to 2. The top
factor in degree 63 is 2. That top factor already follows from the
published torsion-index ratio 24/12; it is not a new top-degree index.

**Good-Levi corollary.** Suppose W_I is a product of type A Weyl groups
and the Cartan submatrix C[:,I] has full column rank modulo 2. Then

\[
 \boxed{N_I^d/M_I^d\simeq(\mathbb Z/2)^{m_{I,d}},\qquad
 \sum_d m_{I,d}q^d=\frac{qP_W(q)}{(1+q)P_{W_I}(q)}.}
 \tag{2}
\]

The hypotheses select 108 proper standard parabolics. An independent
[transfer table](experiments/e7-full-flag-transfer.json) enumerates them.
The separately proved [P7 formula](E7_ROOT_WEIGHT_COMPARISON.md) also
agrees with (2), adding a 109th proper parabolic. No assertion that (2)
holds on the other 18 proper standard parabolics is made.

Among the 64 retained subsets of K={1,...,6}, this proves every degree
for 61 cases. The remaining three are I={2,3,4,5} and the two equivalent
D5 cases {1,2,3,4,5}, {2,3,4,5,6}. Thus two independent all-degree
polynomials remain in that family. Their earlier bounded-degree theorem
remains valid; they are not settled by the present transfer argument.
The subsequent [D5 source proof](E7_D5_REFINEMENT.md) settles the two
D5 cases and leaves only D4 in that family.

For the normalized minimal and adjoint Springer maps, the earlier source
theorem gives J56,I^d=12^d M_I^d and J133,I^d=36^d N_I^d. The previously
proved sum comparison identifies (J56+J133)/J56 with N_I/M_I. Equations
(1)--(2) therefore determine those Springer/Borel image-sum factors.
They do not classify the ambient cohomology cokernels H/N or H/M.

The [irreducible source classification](IRREDUCIBLE_CHARACTER_LATTICES.md)
therefore makes these results usable for any two nontrivial irreducible
E7 representations on the covered parabolics. Each source is a scalar
multiple of Q or X. At 2 the relative local exponents are 0 or 1 with
multiplicities supplied by (1)--(2); scalar powers add their degree times
the Dynkin-index valuations. Sums and intersections then use minimum and
maximum exponents. At odd primes Q and X agree locally. This is elementary
Smith-basis bookkeeping combined with the earlier source theorem; it is
not a classification for arbitrary direct sums or for ambient H/J.

## 2. A common A6 complement at 2

Let A={1,3,4,5,6,7}, so G/P_A is E7/P2 with Levi A6. Put

\[
\begin{split}
 t_1&=-\omega_1+\omega_2,&t_2&=\omega_1+\omega_2-\omega_3,\\
 t_3&=\omega_2+\omega_3-\omega_4,&t_4&=\omega_4-\omega_5,\\
 t_5&=\omega_5-\omega_6,&t_6&=\omega_6-\omega_7,\qquad t_7=\omega_7.
\end{split}
\]

Their determinant in X is 3, a unit in R. W_A permutes the seven t_i.
All seven have the same nontrivial class in X/Q. Set

\[
 \mu=(t_1+\cdots+t_7)/7,\qquad y_i=t_i-\mu,
 \qquad V_A=\sum_{i=1}^6R(t_i-t_7).
\]

Because 7 is a unit, all y_i lie in V_A and sum to zero. The actual
lattice equalities are

\[
 X_R=V_A\oplus R\mu,\qquad Q_R=V_A\oplus R(2\mu).
 \tag{3}
\]

The [plain Python audit](experiments/e7-a6-images-audit.json) checks
both basis matrices have local index a unit, using the full E7 Cartan
matrix. These are equalities at 2, not claimed as global Z splittings.

Write a_j=e_j(y1,...,y7), with a1=0. Elementary symmetric polynomial
theory, together with 7 invertible, gives the full invariant sources

\[
 (\operatorname{Sym}_R X)^{W_A}=R[\mu,a_2,\ldots,a_7],\qquad
 (\operatorname{Sym}_R Q)^{W_A}=R[2\mu,a_2,\ldots,a_7].
 \tag{4}
\]

Here is a direct justification that does not average by 7!. The symmetric
ring in the t_i is R[e1,...,e7]. The change from t_i to centered y_i
rewrites it as R[mu,a2,...,a7], since e1=7mu. Taking the coefficient
independent of mu proves (Sym_R V_A)^(W_A)=R[a2,...,a7]. Coefficientwise
invariance and (3) then give (4).

The obstruction to a six-generator integral E6 invariant ring recorded
in the refinement note remains valid. The intermediary here is A6.

## 3. Rational quotient and exact integral images

Put c_j=e_j(t1,...,t7), with c1=3omega2. The common quadratic Borel
relation is

\[
 c_2=4c_1^2/9.
 \tag{5}
\]

After imposing it, work in the rational polynomial ring

\[
 \mathbb Q[c_1,c_3,c_4,c_5,c_6,c_7],\qquad
 \deg(c_1,c_3,c_4,c_5,c_6,c_7)=(1,3,4,5,6,7).
\]

The 56 minuscule weights in these coordinates are

\[
 \{\pm t_i\}_{i=1}^7\ \cup\
 \{\pm(t_i+t_j-c_1/3):i<j\}.
 \tag{6}
\]

Their power sums of degrees 6,8,10,12,14,18 give the remaining rational
Borel relations. Newton recurrence compresses them into the c_j. Their
ideal has 221 weighted Groebner basis elements and 576 standard
monomials, with degree profile P_W/P_(W_A) and top degree 42. Since the
relations are full-Weyl invariant, they lie in the classical Borel
kernel; equality of quotient dimensions proves they generate that
rational kernel. No integrality of this rational presentation is assumed.

A second construction uses the 126 adjoint roots

\[
 \{t_i-t_j:i\ne j\}\ \cup\
 \{\pm(t_i-c_1/3)\}_i\ \cup\
 \{\pm(t_i+t_j+t_k-c_1/3):i<j<k\}.
\]

Girard partitions replace the primary Newton recurrence. The independently
rebuilt adjoint ideal equals the minuscule ideal in both directions and
has the same dimension. This [ideal audit](experiments/e7-a6-ideal-audit.json)
shares Sage/Singular exact polynomial arithmetic with the primary
calculation; its representation input and construction are different.
The plain Python audit separately verifies both weight sets by Weyl
reflection closure and checks the 221-element basis by Buchberger's
criterion (9,731 nontrivial pairs in the relevant degrees). Its own
leading-ideal enumeration gives exactly the same 576 monomials.

In this rational quotient the full sources in (4) have images generated by

\[
 N_A: c_1,c_2,\ldots,c_7,\qquad
 M_A:2c_1/7,a_2,\ldots,a_7,
\]

where

\[
 a_k=\sum_{j=0}^k {7-j\choose k-j}(-c_1/7)^{k-j}c_j.
 \tag{7}
\]

In particular a2=c1^2/63, so beta(mu)^2=(9/7) beta(a2) is in the root
image. This also gives the exponent-two bound directly on P2.

For each source, begin with 1 in degree zero. In degree d, multiply each
source generator by the already obtained image bases in degree d minus
its degree; take the exact column lattice of all these products. This
induction computes the entire source image. It does not select a smaller
generator algebra because (4) has already identified the full sources.
Hermite bases are taken in the common rational standard-monomial basis;
no integral Schubert normalization is needed for a relative comparison.

The [primary certificate](experiments/e7-a6-images.json) records bases in
all 42 positive degrees. The independent plain Python audit verifies,
for both sources and in every degree:

- every generator product has coordinates in R in the proposed basis;
- those product coordinates have full row rank modulo 2;
- the root basis lies in the weight basis over R;
- twice the weight basis lies in the root basis over R.

The first two checks inductively prove equality with the complete image
lattices. The last two show that the relative group has only factors 2;
its multiplicity is rank minus the mod-2 rank of the inclusion matrix.
The audit uses sparse rational reduction, its own rational matrix inverses
and Python bit elimination, rather than the primary Hermite/Smith route.
All target degrees have been exhausted, since the target vanishes above
42. The resulting exact identity is

\[
 T_A(q)=\frac{q}{1+q}\frac{P_W(q)}{P_{W_A}(q)}.
 \tag{8}
\]

It has 288 nonunit factors in total. This finite coefficient/basis
verification is the computer-assisted part of the theorem.

## 4. Transferring the result to the full flag

The symmetric polynomial ring of the centered y_i is free over
R[a2,...,a7]. One common Artin basis is

\[
 \mathcal B=\{y_1^{b_1}\cdots y_6^{b_6}:
                  0\le b_i\le7-i\}.
\]

It has 7!=5040 elements and degree polynomial P_(W_A)=[2]...[7]. Every
basis element lies in Sym_R V_A, hence in both actual source rings in
(3). The same basis supplies their free-module decompositions over the
two respective invariant sources (4).

These classes restrict to an R-basis of the A6 flag fiber: at 2 the A6
root and weight lattices agree because their index is 7, and the
ordinary type A characteristic map gives the Artin coinvariant basis.
Leray--Hirsch on G/B -> G/P_A then makes the corresponding Borel image
decompositions embedded direct sums. Consequently

\[
 T_{\varnothing}(q)=P_{W_A}(q)T_A(q).
\]

Combining with (8) proves (1), including all multiplicities. This transfer
is a consequence of classical symmetric polynomial and flag-fiber
results; no division by the even Weyl-group order is involved.

More generally, if I satisfies the hypotheses of (2), the actual adjoint
Levi has a 2-local simply connected derived group and type A flag factors.
The usual common polynomial fiber lifts therefore give
T_empty=P_(W_I)T_I, exactly as in the
[E6 good-Levi argument](E6_ROOT_WEIGHT_COMPARISON.md#6-good-levi-deconvolution-60-proper-parabolics).
Division proves (2). The transfer checker independently enumerates all
128 retained sets, verifies the mod-2 Cartan rank and type-A conditions,
and checks all polynomial divisions and nonnegative multiplicities.
The P7 agreement uses its separate completed basis proof, not an
unavailable type E6 non-torsion fiber argument at 2.

## 5. Direct full-flag corroboration and the missing integral relations

A separate [full-flag kernel calculation](experiments/e7-full-flag-kernel.json)
recomputes degrees 10--13 without the A6 transfer. In the earlier E6
complement V plus Rw, it reduces by the monic quadratic and generates
the rational kernel from fresh 56-weight moments. There are no ordinary
relation syzygies below degree 14. Its integer coefficient column
lattices require the following 2-saturations:

| degree | kernel rank | initial mod-2 defect | exact multiplicity of 2 |
| --- | --- | --- | --- |
| 10 | 210 | 0 | 1940 |
| 11 | 462 | 1 | 2855 |
| 12 | 924 | 8 | 4054 |
| 13 | 1716 | 35 | 5578 |

Each defect is removed by one round of exact integral half-sums, after
which the coefficient lattice has full column rank modulo 2 and is
therefore the entire local kernel. Odd-w coefficient ranks then compute
the relative factors. These match (1). The degree-11 correction confirms
that the rational relation multiples alone do not give the full integral
kernel, even though the final multiplicity prediction survives.

The [independent kernel audit](experiments/e7-full-flag-kernel-audit.json)
reconstructs the initial moments from 27 E6 forms, independently verifies
the resulting 56-weight set, and replays every saturation witness using
integer coefficients modulo 2^64. It explicitly tracks precision lost
at each halving and retains at least 63 valid binary digits. Its own bit
ranks check both the complete kernel and odd coefficient images. These
are exact finite congruence certificates, not numerical approximations.
They corroborate the separate all-degree proof; they are not its basis.

## 6. Attribution, limitations and next step

Published ingredients are rational Borel presentations, elementary
symmetric/Artin bases, integral type A characteristic maps and
Leray--Hirsch. The full lattice split and finite complete-image verification
are the applications recorded here. Their exact literature priority is
unresolved. Sources for the framework include
[Demazure](https://doi.org/10.1007/BF01418790),
[Calmes--Petrov--Zainoulline](https://arxiv.org/abs/0905.1341),
[Kumar Part I](https://kumar.math.unc.edu/papers/kumar76.pdf) and
[Kostant--Michor](https://arxiv.org/abs/math/0109066).
[Kaji--Nakagawa](https://arxiv.org/abs/0709.3702) gives the torus coordinates
and related integral flag presentations. The independent trace construction
above does not assume their printed high-degree integral relations.

[Totaro](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf)
gives the E7 adjoint/simply-connected torsion indices 24/12. The full-flag
top factor is thus known background. Related mod-2 characteristic-image
and canonical-dimension data appear in
[Zainoulline](https://www.math.uni-bielefeld.de/LAG/man/199.pdf).
Those concern the image after reduction in ambient cohomology. That image
is different from N/2N when N is nonprimitive; such mod-2 data alone do
not supply the relative integral factors proved here. The source checks
performed do not establish a novelty claim for (1)--(2).

The subsequent [D5 completion](E7_D5_REFINEMENT.md) settles both D5
refinements using their full integral invariant source, bringing the
totals to 111 proper E7 parabolics and 63/64 P7 refinements. D4 remains.
The torsion-prime fiber prevents applying the common unit-index
argument as stated. The six-generator E6 source obstruction remains
relevant to those routes. Other excluded E7 parabolics, ambient H/N
factors, direct-sum representations and equivariant refinements remain
separate questions. Runtime and the installed library API stay at 0.20.0.

## Reproduction

From borel/:

```sh
sage -python experiments/e7_a6_images.py
sage -python experiments/audit_e7_a6_ideal.py
python experiments/audit_e7_a6_images.py
sage -python experiments/e7_full_flag_kernel.py --max-degree 13
python experiments/audit_e7_full_flag_kernel.py
python experiments/e7_full_flag_transfer_check.py
```
