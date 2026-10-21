# Sharp 3-primary root-to-weight images in E6

Status, 20 October 2026: a written deduction for the Cayley plane and 60
proper standard parabolics. A common relative Schubert basis first
transfers the exact factors to minuscule refinements, then a good-Levi
deconvolution covers 12 more cases. The [three-case completion](E6_REMAINING_CASES.md)
now settles the other proper parabolics, with [complete tables](E6_REMAINING_TABLES.md)
and a separate [proof/source audit](E6_REMAINING_AUDIT.md).
The space G/G has zero comparison quotient. No claim of literature
priority is made. See the [proof/source audit](E6_ROOT_WEIGHT_AUDIT.md).

## 1. The objects compared

Let G be simply connected complex E6, with Bourbaki/Sage nodes and edges
(1,3), (3,4), (2,4), (4,5), (5,6). Let X be its weight lattice and Q its
root lattice, so X/Q is cyclic of order 3. For retained nodes I, let P_I
be the standard parabolic with Levi Weyl group W_I, and put

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^d Q)^{W_I}),\qquad
 N_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^d X)^{W_I}).
\]

Both are embedded integral lattices in H^(2d)(G/P_I,Z), of full rank,
and M_I^d is contained in N_I^d. Since 3X is contained in Q, the quotient
N_I^d/M_I^d is finite and 3-primary. It suffices to work over R=Z_(3).
All degrees below are polynomial degrees, half cohomological degrees.

These are unscaled characteristic images. For the minimal and adjoint
Springer representations the earlier source theorems give

\[
 J_{27,I}^d=6^d M_I^d,\qquad J_{78,I}^d=24^d N_I^d.
\]

Their [previous all-Levi comparison](E6_REPRESENTATION_COMPARISON.md)
therefore identifies (J27+J78)/J27 with N/M in each degree.

## 2. Exact theorem and multiplicity polynomials

Let K={2,3,4,5,6}, so P_K is the minuscule P1 and G/P_K is the Cayley
plane. Define three polynomials

\[
\begin{aligned}
 C_1(q)&=q+q^2+q^7+q^8+q^9+2q^{10}+q^{11}+q^{12}+q^{14},\\
 C_2(q)&=q^3+q^4+q^6+q^7+q^9+q^{11}+q^{12}+q^{13}+q^{15}+q^{16},\\
 C_3(q)&=q^5+q^8.
\end{aligned}
\tag{1}
\]

Write [r]_q=1+q+...+q^(r-1), put
D(q)=[2]_q[4]_q[5]_q[6]_q[8]_q, and let P_(W_I)(q) be the Weyl length
polynomial of the retained Levi. The variable q just records degree.

**Theorem.** Let I be a proper retained set that does not contain
{1,3,5,6}. Then C_a(q)D(q)/P_(W_I)(q) is a polynomial with nonnegative
integer coefficients, and

\[
 N_I^d/M_I^d\simeq\bigoplus_{a=1}^3
 (\mathbb Z/3^a)^{[q^d]C_a(q)D(q)/P_{W_I}(q)}.
\tag{2}
\]

This covers 60 of the 63 proper standard parabolics. The excluded
retained sets are {1,3,5,6}, {1,2,3,5,6}, and {1,3,4,5,6}, of types
A2+A2, A2+A2+A1, and A5. The separate G/G comparison is zero.
The good-Levi argument proving this assertion is in Section 6. The
excluded cases have different multiplicities, proved in the linked
completion; formula (2) is asserted only under its stated hypothesis.

**Minuscule-refinement form.** Suppose I is contained in K. Let

\[
 B_{K/I}(q)=\sum_{v\in W_K/W_I\text{ minimal}}q^{\ell(v)}
 \quad\text{and}\quad m_{a,I,d}=[q^d]C_a(q)B_{K/I}(q).
\]

Then in every degree

\[
\boxed{
 N_I^d/M_I^d\simeq
 (\mathbb Z/3)^{m_{1,I,d}}\oplus
 (\mathbb Z/9)^{m_{2,I,d}}\oplus
 (\mathbb Z/27)^{m_{3,I,d}}.}
\tag{2a}
\]

The same assertion holds for I contained in {1,2,3,4,5}, by the E6
diagram involution exchanging nodes 1 and 6, and 3 and 5. These two sets
cover 32+32-16=48 standard parabolics. In the overlap either computation
describes the same quotient.

The unit Smith factors have multiplicity
rank H^(2d)(G/P_I,Z) minus the three multiplicities in (2).
For G/B the fiber polynomial is the D5 Weyl-group length polynomial,
with fundamental degrees 2,4,5,6,8. Thus (2) also gives every full-flag
factor without constructing its 51,840-dimensional total cohomology.

For the two-node omission {1,6}, the fiber is the 8-dimensional quadric:

\[
 B_{K/I}(q)=1+q+q^2+q^3+2q^4+q^5+q^6+q^7+q^8.
\]

## 3. Source splitting and classical invariant theory

Write h0=omega1 and use the following five D5 coordinates:

\[
\begin{aligned}
 x_1&=-\omega_1-2\omega_2+2\omega_3,&
 x_2&=-\omega_1+2\omega_2+2\omega_3-2\omega_4,\\
 x_3&=-\omega_1+2\omega_4-2\omega_5,&
 x_4&=-\omega_1+2\omega_5-2\omega_6,\\
 x_5&=-\omega_1+2\omega_6.
\end{aligned}
\]

The vectors (3h0,x1,...,x5) lie in Q with index 64; (h0,x1,...,x5)
lie in X with index 64. Thus, over R,

\[
 Q_R=R(3h0)\oplus V,\qquad X_R=R h0\oplus V,
 \qquad V=\bigoplus_{j=1}^5 Rx_j.
\tag{3}
\]

W_K fixes h0 and acts on the x's by even signed permutations. Since
2 is invertible, its invariant ring is

\[
 R[V]^{W_K}=R[e_1,e_2,p,e_3,e_4],\qquad
 e_j=e_j(x_1^2,\ldots,x_5^2),\quad p=x_1\cdots x_5.
\]

Consequently the two complete K-invariant sources are
R[3h0,e1,e2,p,e3,e4] and R[h0,e1,e2,p,e3,e4]. This is a classical D5
invariant-ring calculation, combined with the explicit lattice splitting.

## 4. Cayley-plane images and their exact local bases

Use h=beta(omega1), and let s be the degree-four Schubert class with
word (2,4,3,1). The published Duan--Zhao presentation, transported from
P6 to P1, is

\[
 H^{2*}(G/P_1,\mathbb Z)=
 \mathbb Z[h,s]/(2h^9-6h^5s+3hs^2,\ h^{12}-6h^4s^2+s^3).
\tag{4}
\]

See [Duan--Zhao, Theorem 4](https://arxiv.org/abs/math/0511332).
The invariant-to-Schubert calculation is classical in method; see
[Iliev--Manivel, Proposition 6.1](https://arxiv.org/abs/math/0306329)
for a published D5-invariant calculation in their coordinates.
Our chosen coordinates give the following six root-source images:

\[
\begin{array}{ll}
 A=3h,&B=-3h^2,\\
 C=42h^4-96s,&D=-21h^5+48hs,\\
 E=282h^6-672h^2s,&
 F=-1371h^8+3552h^4s-768s^2.
\end{array}
\tag{5}
\]

They are the images of (3h0,e1,e2,p,e3,e4), respectively.
The coefficient identities in (5) have been recomputed by Sage BGG and
independently by sparse binomial divided differences in ordinary Python.

It follows that

\[
 M_K\otimes R=R[A,B,C,D,E,F],\qquad
 N_K\otimes R=R[h,3s,3s^2].
\tag{6}
\]

For the second equality, C and h supply 3s because 32 is a unit, and F
then supplies 3s^2 because 256 is a unit. The other generators lie in
that ring. This calculation uses all Levi invariants, not a selected
proper subalgebra.

Reducing (4) modulo 3 gives F3[h,s]/(h^9,s^3). Since integral cohomology
is torsion-free, the 27 monomials h^i s^j, 0<=i<=8 and 0<=j<=2,
are an R-basis. In this basis N has coefficients 1 for j=0 and 3 for
j=1,2. This description is closed under h,3s,3s^2 using (4).

Here are local bases for M and the complete Smith exponents relative
to N. Exponent zero records a unit factor.

| d | R-basis for M_K^d | Smith exponents |
| --- | --- | --- |
| 0 | 1 | 0 |
| 1 | A | 1 |
| 2 | B | 1 |
| 3 | AB | 2 |
| 4 | C, B^2 | 0,2 |
| 5 | D, AB^2 | 0,3 |
| 6 | E, BC | 0,2 |
| 7 | BD, AE | 1,2 |
| 8 | F, BE, B^2C | 0,1,3 |
| 9 | CD, B^2D | 1,2 |
| 10 | D^2, CE | 1,1 |
| 11 | DE, BCD | 1,2 |
| 12 | E^2, CF | 1,2 |
| 13 | BDE | 2 |
| 14 | EF | 1 |
| 15 | CDE | 2 |
| 16 | F^2 | 2 |

This is an exact finite algebra proof. Each displayed basis element is
a source monomial. Their R-span contains 1 and is closed under all six
generators, using only (4); all higher target degrees vanish. Therefore
it equals the entire generated image ring, rather than just a lattice
inside it. The determinant/minor valuations of its degree matrices give
the last column. The [certificate](experiments/e6-root-weight-comparison.json)
records every matrix, selected monomial basis and two-sided witness;
193 local module closure checks are performed. The
[independent audit](experiments/e6-root-weight-independent.json)
rebuilds all 17 image matrices and their determinantal valuations.
Reading the nonzero exponents gives precisely (1).

For example, degree five uses the N-basis (h^5,3hs). The two M columns
D and AB^2 give the matrix ((-21,27),(16,0)). Its entries have a unit
gcd and its determinant has 3-adic valuation 3. Hence the quotient is
Z/27, not merely bounded by 3^5.

## 5. Common relative Schubert basis: the transfer argument

This step explains why the factors transfer and why independently
choosing bases for the two source rings would be insufficient.

Let I be contained in K, and put A0=Sym_R(V), B0=A0^(W_K).
The D5 torsion index t is a power of 2, hence a unit in R. Choose a
degree-20 polynomial F0 in A0 with partial_(w0)(F0)=t and define

\[
 b_w=t^{-1}\partial_{w^{-1}w_0}F_0\qquad(w\in W_K).
\]

The classical full-flag characteristic-map theorem identifies
A0/(B0,+) with H*(D5/B,R), and the b_w map to its Schubert basis.
Graded Nakayama gives generation over B0; the rational invariant Hilbert
series gives independence. Thus the b_w form a B0-basis of A0.

The nil-Hecke relations give exact polynomial identities

\[
 \partial_i b_w=
 \begin{cases}b_{ws_i},&\ell(ws_i)<\ell(w),\\
 0,&\ell(ws_i)>\ell(w).
 \end{cases}
\]

Coefficients in B0 are invariant, so partial_i is B0-linear. Writing
f=sum_w c_w b_w, the condition partial_i(f)=0 forces c_w=0 for every
w with a right descent at i: their images are distinct basis elements.
Since the polynomial ring is a domain, partial_i(f)=0 is equivalent
to s_i(f)=f. Intersecting these conditions for i in I proves

\[
 A_0^{W_I}=\bigoplus_{v\in W_K/W_I\text{ minimal}} B_0 b_v.
\]

Each b_v has degree length(v), is W_I-invariant, and lifts the fiber
Schubert class. No division by |W_I| is involved. The classical inputs
are discussed in [Demazure](https://doi.org/10.1007/BF01418790)
and [Calmes--Petrov--Zainoulline, Sections 5--6](https://arxiv.org/abs/0905.1341).

Crucially, b_v involves only the common V in (3). Therefore the *same*
polynomials give both decompositions

\[
 (\operatorname{Sym}_R Q_R)^{W_I}
   =\bigoplus_v b_v R[3h0,B_0],\qquad
 (\operatorname{Sym}_R X_R)^{W_I}
   =\bigoplus_v b_v R[h0,B_0].
\tag{7}
\]

Consider the projection pi:G/P_I -> G/P_K. The classes beta(b_v)
restrict to a basis of its D5/P_I fiber. Leray--Hirsch gives an
H*(G/P_K,R)-basis of H*(G/P_I,R). Applying beta to (7) consequently
gives embedded direct sums

\[
 M_I\otimes R=\bigoplus_v\beta(b_v)\pi^*(M_K\otimes R),\qquad
 N_I\otimes R=\bigoplus_v\beta(b_v)\pi^*(N_K\otimes R).
\tag{8}
\]

There is no hidden extension between these summands: their independence
holds in the ambient cohomology by Leray--Hirsch. Taking the quotient of
(8) proves (2), including all multiplicities. The original quotient is
3-primary, so its localization recovers the same finite abelian group.

The mechanism is general: a common complementary lattice, a fiber with
torsion index invertible at p, and common invariant Schubert lifts yield
the same transfer rule. It is a formal consequence of these classical
inputs, rather than a new Leray--Hirsch theorem.

## 6. Good-Levi deconvolution: 60 proper parabolics

Let L_I^ad be the Levi in the adjoint E6 group with retained set I,
and L_I^sc its inverse image in the simply connected group. Both have
the same flag fiber F_I=L_I/B_I and Weyl group W_I. Suppose the
3-part of the torsion index of L_I^ad is one. Then there are homogeneous
polynomials b_v in Sym_R(Q_R) whose restrictions to F_I form its
Schubert basis. They restrict to the same basis for L_I^sc, because
central isogenies identify the flag fibers.

Here is the general common-basis argument without the splitting (3).
Classical restriction at a non-torsion prime identifies H*(BL_I^ad,R)
with (Sym_R Q_R)^(W_I), and likewise for L_I^sc and X. Leray--Hirsch
on F_I -> BT -> BL_I gives free polynomial-module decompositions
over these respective invariant rings with the same b_v. Applying
the characteristic maps and Leray--Hirsch again, now to
G/B -> G/P_I, gives

\[
 M_B\otimes R=\bigoplus_v\beta(b_v)\pi^*(M_I\otimes R),\qquad
 N_B\otimes R=\bigoplus_v\beta(b_v)\pi^*(N_I\otimes R).
\tag{9}
\]

No relation between the two invariant rings is lost, since the same
basis is used in both source modules and in ambient cohomology. Thus,
if T_(I,a)(q) counts the cyclic factors of order 3^a, then

\[
 T_{B,a}(q)=P_{W_I}(q)T_{I,a}(q).
\]

Section 5 already gives T_(B,a)=C_a D from the D5 minuscule projection.
Dividing therefore proves (2). Polynomiality and nonnegative coefficients
follow from (9), not from an assumption about a rational function.

For completeness, the required Levi hypothesis can be checked without
enumerating any representations. With C the E6 Cartan matrix, the
restriction of Q to the simply connected derived Levi weight lattice
has coefficient matrix C[:,I]. Its reduction modulo 3 has full column
rank exactly when I does not contain {1,3,5,6}: the kernel of C modulo
3 is generated by (1,0,2,0,1,2). Full column rank says that the derived
central isogeny has degree prime to 3. Every proper E6 derived Levi
is a product of types A and D; its simply connected torsion index is
a power of 2. Consequently the adjoint-form Levi also has torsion
index a unit at 3 under this rank condition. This supplies the b_v
used in (9).

Conversely the three exceptional proper retained sets contain the
displayed support; the common-fiber-basis hypothesis does not follow.
In these three cases C_a D/P_(W_I) has a nonzero polynomial-division
remainder for every a. Thus this particular deconvolution rule cannot
be asserted for them. Their actual image multiplicities are determined
by the separate [three-case completion](E6_REMAINING_CASES.md), not by
this deconvolution rule. For published restriction/divided-difference background, see
[Holm--Sjamaar, Section 2](https://arxiv.org/abs/math/0607069), as well
as the characteristic-map references in Section 5.

## 7. Consequences and precise limits

For every E6 parabolic, the quotient N/M is killed by 27. Indeed Totaro's
published t(E6 adjoint)=54 and integral divided-difference lifts give
54H^(2d)(G/P,Z) inside the root characteristic image; at 3 this is
27H. This is a classical consequence of the adjoint torsion index.
Degree five of the Cayley plane attains 27, so the uniform exponent
bound is sharp. See [Totaro, Theorem 6.1](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf).

In full-flag top degree the relative quotient is Z/9. This agrees with
the already known ratio t(E6 adjoint)/t(E6 simply connected)=54/6=9.
The top-degree number is therefore not a new torsion-index calculation.
The all-degree factors and common-basis transfer are the deductions
recorded here; their exact literature priority remains unresolved.

The theorem also removes the coprimality restriction for scaled root
versus weight image comparisons in its 60-parabolic scope. If M in N
has local diagonal exponents s_i, the images a^d M and b^d N have,
in a common local N-basis, exponents

\[
 u_i=d v_3(a)+s_i,\qquad v_i=d v_3(b).
\]

Their sum uses min(u_i,v_i), and their intersection uses max(u_i,v_i).
The quotient of the sum by a^d M has exponents max(u_i-v_i,0), while
the quotient by b^d N has exponents max(v_i-u_i,0). This is elementary
Smith-basis arithmetic, now with the s_i supplied by (2). The previous
[irreducible source classification](IRREDUCIBLE_CHARACTER_LATTICES.md)
reduces every irreducible E6 source to a scalar multiple of Q or X.
It does not give the same shortcut for arbitrary direct sums.

The three proper retained sets listed in Section 2 are now settled by
the [full-source certificate proof](E6_REMAINING_CASES.md), including the
maximal A5 and A2+A2+A1 Levis. That proof uses an A5+A1 subsystem at 3;
the bound 27 alone would not determine their multiplicities.
The A4+A1 maximal Levis are covered by the good-Levi argument. Direct sums,
equivariant specialization kernels and ambient H/N factors are separate
questions.

## 8. Reproduction

From borel/ run the following in an environment with Sage, and run the
last command with Python containing SymPy:

```sh
PYTHONPATH=src sage -python experiments/e6_root_weight_comparison.py
PYTHONPATH=src sage -python experiments/e6_root_weight_refinement_check.py
python experiments/verify_e6_root_weight_independent.py
```

The first run reconstructs the six BGG images and all finite certificates.
The second directly recomputes full-flag degrees 0--4, the D4-refined
flag degrees 0--6, the retained A1+A1 flag {1,6} in degrees 0--4,
and maximal P3 in degrees 0--4. The last uses full integral invariant
reflection kernels. None uses the transfer model to supply its columns.
The last uses neither Sage nor sage_borel: it checks the lattice indices,
all 27 minuscule labels, 11 BGG coefficients, all 17 base lattices and
all 32 D5 fiber polynomials by even signed permutations. It also
enumerates the 63 proper E6 Levi Weyl groups on the 27-weight orbit,
checks the local Cartan-rank criterion and all 60 deconvolution formulas.
Finite
checks corroborate the written argument; they do not replace the
classical geometry or certify literature priority.
