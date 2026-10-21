# E7 D4 refinement and completion of the P7 family

This note completes the remaining D4 case by an exhaustive, computer-assisted
integral image calculation. It builds on the independently audited
[D5 theorem](E7_D5_REFINEMENT.md), retaining all integral source generators.
The new arithmetic has its own native-Python audit. This is a finite-type
deduction, with no claim of established literature priority.

## 1. Statement and connection to the Springer maps

Let G be simply connected of type E7 over C. Use Bourbaki numbering, with
edges 1--3, 3--4, 2--4, 4--5, 5--6, 6--7. Let X be the weight lattice,
Q the root lattice, W the Weyl group, and R=Z_(2). A retained set I denotes
the simple roots of the Levi of P_I. Put

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^d Q)^{W_I}),\qquad
 N_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^d X)^{W_I}).
\]

These are the full unscaled embedded integral Borel images. Polynomial
degree d corresponds to cohomological degree 2d.

**D4 theorem.** For I={2,3,4,5}, in every degree,

\[
 N_I^d/M_I^d\simeq(\mathbb Z/2)^{m_d},\qquad
 \sum_{d=0}^{51}m_dq^d
 =\frac{qP_{W(E7)}(q)}{(1+q)P_{W(D4)}(q)}.                 \tag{1}
\]

Here P_(W(D4))=[2]_q[4]_q^2[6]_q,
P_(W(E7))=[2][6][8][10][12][14][18], and [a]_q=1+...+q^(a-1).
An explicit polynomial form is

\[
 q(1+q^5)(1+q^4+q^8)(1+q^4)[5]_q[7]_{q^2}[18]_q.
\]

There are 7,560 nonunit relative factors in total. The flag variety has
complex dimension 51 and 15,120 Schubert classes. Both images vanish
above degree 51.

**Completed-family corollary.** Together with the
[P7 basis theorem](E7_ROOT_WEIGHT_COMPARISON.md), the
[A6/full-flag transfer](E7_FULL_FLAG_COMPARISON.md) and the D5 theorem,
for **every** I contained in {1,...,6},

\[
 N_I^d/M_I^d\simeq(\mathbb Z/2)^{[q^d]
     qP_{W(E7)}(q)/((1+q)P_{W_I}(q))}.                    \tag{2}
\]

Thus all 64 standard refinements of E7/P7 are covered in all degrees.
Across the existing transfer scope this makes 112 proper E7 parabolics;
the other 15 proper parabolics remain outside this proof.

**Scope criterion.** Formula (2) holds more generally for exactly the
112 retained sets I that omit at least one of nodes **2,5,7**. Equivalently,
the restriction of parabolic characters X^(W_I) to X/Q is surjective.
This follows by identifying the existing 108 good-Levi cases and the
four additional P7/D5/D4 cases with that set of 112 diagrams, as checked
in the summary certificate. It packages the case proofs into a statement
with an intrinsic hypothesis, rather than replacing them with a new
uniform central-cover argument.

Indeed the E7 Cartan determinant is 2 and its mod-2 null vector has
support {2,5,7}. Thus sum(a_j omega_j) lies in Q exactly when
a2+a5+a7 is even. Since X^(W_I) is spanned by the fundamental weights
whose nodes are omitted, the criterion follows. Every one of the other
15 proper retained sets contains all three nodes. For those sets
N_I^1=M_I^1. Formula (2), which always has degree-one multiplicity 1,
therefore **cannot** extend unchanged to them. Their higher-degree
relative images remain a separate problem.

For all 112 cases there is also an algebra interpretation. Choose any
omitted node i in {2,5,7}.
Put N_I=direct-sum_d N_I^d and
let C_I be the **image** of M_I in N_I/2N_I. It is not M_I/2M_I.
Let w=beta(omega_i), and let B_i be the common root-source quadratic
constructed in Section 2. Then

\[
 N_I/2N_I\simeq C_I[t]/(t^2-\overline{\beta(B_i)}),\qquad
 \deg t=1,\quad t\longmapsto\bar w.                    \tag{3}
\]

In particular the actual weight image modulo 2 is free over the image
of the root source modulo 2, with basis (1,bar(w)). To prove this
interpretation, work over R and write A_image for the image of the common
invariant source A_i defined in Section 2. The source split gives N_I=A_image+w A_image
and M_I=A_image+2w A_image. It therefore gives the displayed quadratic
surjection. The prior inclusion 2N_I contained in M_I gives an exact
sequence 0 -> C_I -> N_I/2N_I -> N_I/M_I -> 0. The multiplicity formula says
Hilb(C_I)=P_(G/P_I)/(1+q). Both sides of the surjection have the same
graded dimension, proving (3). This is a formal consequence of (2) and
the existing quadratic identity, rather than a separate unverified
freeness claim about integral invariant sources.
Consequently N_I/M_I is the C_I-module bar(w) C_I, isomorphic to C_I
with its grading shifted up by one. The factors 2 record exactly these
additional weight-image classes.

For the normalized Springer morphisms of the minimal 56-dimensional and
adjoint 133-dimensional representations, the established
[source classification](IRREDUCIBLE_CHARACTER_LATTICES.md) gives

\[
 \mathcal J_{56,I}^d=12^dM_I^d,\qquad
 \mathcal J_{133,I}^d=36^dN_I^d.
\]

Equations (1)--(2) also describe
(J_56,I^d+J_133,I^d)/J_56,I^d. They quantify the relative integral
effect of changing the Springer representation before applying the
Kumar map beta composed with (theta_rho^*)^-1. They do not determine
the ambient cohomology quotients H/N or arbitrary direct-sum sources.

## 2. Common quadratic lifts and the full integral D4 source

First the common-complement argument works uniformly for the three
possible omitted nodes. For i=2,5,7 let K_i be the other six nodes,
C_i their Cartan matrix, n_i=det(C_i), V_i=sum_(j in K_i) R alpha_j,
and eta_i=omega_i. The three n_i are respectively 7,15,3, all odd.
The determinants of the bases (alpha_(K_i),eta_i) and
(alpha_(K_i),2eta_i), relative to X and Q respectively, both have
absolute value n_i. Consequently

\[
 X_R=V_i\oplus R\eta_i,\qquad Q_R=V_i\oplus R(2\eta_i).
\]

For I contained in K_i the full sources are A_i[eta_i] and A_i[2eta_i],
where A_i=(Sym_R V_i)^(W_I). The Gram matrix on (alpha_(K_i),eta_i)
is diag(C_i,n_i/2). The invariant quadratic tensor and the rational
Borel relation therefore give

\[
 B_i=-\frac{n_i}{2}\alpha_{K_i}^{\mathsf T}C_i^{-1}\alpha_{K_i},
 \qquad \beta(\eta_i)^2=\beta(B_i).
\]

This is an integral root polynomial: n_i C_i^-1 is integral, and its
diagonal is even because C_i is invertible alternating modulo 2.
For completeness, the inverse of an invertible alternating matrix is
alternating: substitute v=C_i u in v^T C_i^-1 v=u^T C_i u=0 modulo 2.
The off-diagonal terms in the displayed quadratic already contain the
factor 2. Hence B_i lies in Sym_Z^2(sum_(j in K_i) Z alpha_j), and is
W_(K_i)-invariant. This proves N_I=A_image+w A_image,
M_I=A_image+2w A_image and 2N_I contained in M_I at 2, uniformly on
the 112 cases. The summary audit checks the three determinants, the
quadratic coefficients and all 18 simple-reflection identities using
separate exact rational arithmetic. The mod-2 algebra assertion (3)
then follows from the multiplicity formula and the dimension argument
above. These quadratic/complement facts are elementary consequences of
known root-lattice and Borel theory; they do not supply the multiplicities
without the separate image proofs.

For the D4 calculation specialize to i=7. Put K={1,...,6},
V=R alpha1+...+R alpha6 and w=omega7. The W_K-stable split is
X_R=V direct-sum Rw and Q_R=V direct-sum R(2w). Hence the actual sources
are A[w] and A[2w], where A=(Sym_R V)^(W_I).
The common quadratic lift from the P7 proof already implies 2N_I is
contained in M_I for every I contained in K. Consequently only factors
2 can occur; it remains to determine their multiplicities.

Write Omega=C6^-1 alpha for the E6 fundamental weights and set

\[
\begin{split}
 h&=\Omega_1,& z&=\Omega_3,&y&=\Omega_6,&v&=\Omega_5-\Omega_6,\\
 u_1&=-\Omega_2+\Omega_3,&u_2&=\Omega_2+\Omega_3-\Omega_4,
 &u_3&=\Omega_4-\Omega_5.
\end{split}
\]

Then u1+u2+u3=2z-v-y. The six vectors (h,y,z,v,u1,u2) have determinant
of absolute value 1 in the E6 weight lattice; replacing u1,u2 by
u1-u3,u2-u3 gives index 3. Also det(C6)=3. All these indices are units
in R. The subgroup on nodes {2,4} permutes u1,u2,u3 and fixes h,y,z,v.
The symmetric-polynomial theorem, after centering by (2z-v-y)/3, gives
its **full** invariant ring

\[
 P_R=R[h,y,z,v,c_2,c_3],\qquad c_j=e_j(u_1,u_2,u_3).       \tag{4}
\]

Thus coefficientwise 2-local integrality in (4) tests actual source
integrality. D4 acts by even signed permutations on
x_i=2u_i-h (i=1,...,4), with u4=v, and fixes h,y. Over Q the source is
polynomial on h,y,e1(x^2),e2(x^2),prod(x_i),e3(x^2), of degrees
(1,1,2,4,4,6). This rational statement does not imply integral generation.

The [source construction](experiments/e7_d4_source.py) finds ten integral
generators of degrees

\[
                    (1,1,2,4,4,5,6,8,8,8).               \tag{5}
\]

The [exact source report](experiments/e7-d4-source-algebra.json) gives every
polynomial. Its mod-2 graph has a certified 58-element Groebner basis;
eliminating three auxiliary variables gives 37 kernel relations.
The Hilbert numerator relative to the ten generator variables is
(1-q5)(1-q8)^3, so the generated mod-2 subalgebra has Hilbert series

\[
             \frac{1}{(1-q)^2(1-q^2)(1-q^4)^2(1-q^6)}.
\]

This equals the rational invariant Hilbert series. The actual invariant
source is saturated in P_R: if 2f is invariant then f is invariant.
It therefore injects modulo 2 into P_R/2. Equality of each finite
graded dimension, followed by Nakayama, proves (5) generates **all** A,
in every degree. The bounded generator search alone is not the proof.

## 3. Why a free D5 transfer fails, and the replacement

Let J={2,3,4,5,6}, and B=(Sym_R V)^(W_J), the full D5 source already
proved. Here A is an algebra over B. Embed B explicitly by appending
u4=v and u5=y to the three u variables in (4). If s_j=e_j(u1,u2,u3), then
the five-variable elementary functions are
s_j+(v+y)s_(j-1)+vy s_(j-2).

The rational D5/D4 fiber has rank 10, with degrees
(0,1,2,3,4,4,5,6,7,8). However the exact quotient A/(2,B_+A) has dimension
**34**, not 10. Thus an integral free rank-10 invariant-source transfer
cannot be used. This obstruction concerns the invariant source, not the
Schubert cohomology of the quadric fiber.

The relative mod-2 ideal has a certified 20-element Groebner basis.
Its 34 standard monomials lift to a complete graded B-module generating
set b_i of A, by graded induction and Nakayama in each degree. Their
degrees range from 0 to 19. All are retained in the image calculation.

For the rational tower put t=x5=2y-h, r=x1*x2*x3*x4 and
E_j=e_j(x1^2,...,x5^2), while E5 denotes prod(x1,...,x5) as in the D5
files. Thus the labels E1,E2,E3,E4 are square elementary functions,
whereas E5 is the degree-five product. The relations are

\[
\begin{split}
 tr&=E_5,\\
 r^2&=E_4-E_3t^2+E_2t^4-E_1t^6+t^8,\\
 t^9&=E_5r-E_4t+E_3t^3-E_2t^5+E_1t^7.                \tag{6}
\end{split}
\]

With weights deg(t)=1, deg(r)=4 these form a monic Groebner basis with
leads tr,r^2,t^9. The two nontrivial S-pairs reduce to zero; the third
has coprime leads. Hence (1,t,...,t8,r) is a free basis over the rational
D5 source, and remains a basis after the rational E7 Borel quotient.
This also proves injectivity of the base quotient in this tower.

Choose the integral quartic

\[
 f_4=-h^3z-h^2yz+h^2z^2-hvc_2-hc_3+2vc_3.
\]

Its coefficient of r in (6) is 1/8. Therefore
(1,y,...,y8,f4), with the same fiber degrees, is another rational basis.
Every element of this basis belongs to the common integral source A.
No assertion that it is primitive in the ambient integral fiber is needed.

## 4. Index comparison lemma and exhaustive certificate

Here is the reusable lattice argument. Suppose E is a rational algebra
free over a rational base D with a homogeneous basis F_j in a common
integral source. Suppose the full finer source is generated over the
base integral source by finitely many b_i. For each degree form the
reference lattices

\[
 L_N^d=\bigoplus_j F_jN_J^{d-\deg F_j},\qquad
 L_M^d=\bigoplus_j F_jM_J^{d-\deg F_j}.
\]

These embed by rational freeness, and lie in the actual finer images
N_I,M_I. Those images are generated by all b_i N_J and b_i M_J,
respectively. If a_N(d)=log2[N_I^d:L_N^d] and
a_M(d)=log2[M_I^d:L_M^d], then index multiplicativity gives

\[
 \log_2[N_I^d:M_I^d]
 =\sum_j m_{J,d-\deg F_j}+a_N(d)-a_M(d).                \tag{7}
\]

This lemma is elementary lattice algebra. Equality of the two correction
indices is a separate fact to check, not a consequence of the lemma.

The [primary image calculation](experiments/e7_d4_images.py) expresses all
ten generators in rational signed invariants and all 34 b_i in the
reference basis above. It uses the previously audited D5 rational Borel
quotient, where E1=-3h^2-(8/3)w^2, and its complete M_J,N_J image bases.
For each degree and each image, it constructs the coefficient matrix C
of every b_i times every base image basis vector, in L_M or L_N coordinates.

If 2^k clears all even denominators, the extension index is computed
exactly from the integer column lattice generated by 2^k C and 2^k times
the identity, after reducing odd denominators modulo 2^k. This finite
modulus suffices because the reference identity lattice is included.
PARI modular HNF gives the primary index, with full residue matrices saved.

The [complete image report](experiments/e7-d4-images.json) exhausts degrees
0..51 and proves a_N(d)=a_M(d) in each. The reference lattices are not
the actual images: already a_N(4)=a_M(4)=1, and at degree 26 both indices
have 2-adic valuation 420. Equality, rather than their vanishing, proves
that (7) is the D5 polynomial times (1+q4)[5]_q, giving (1).
Vanishing above the top dimension and the prior exponent-2 theorem
finish the all-degree assertion. Since X/Q has order 2, localization
at every odd prime identifies the two sources; the result is integral.

## 5. Independent verification and attribution

The [native source audit](experiments/audit_e7_d4_source.py) uses no Sage
or project imports. It independently checks the lattice indices and
reflections, all ten rational invariant identities, both containments
of the graph ideal, Buchberger reductions and the monomial-ideal Hilbert
numerator. It rebuilds the D5 embedding and independently derives the
relative ideal by Buchberger, verifying its 34 standard monomials.

The [native image audit](experiments/audit_e7_d4_images.py) independently
reconstructs every new substitution and all 34 tower coefficient lists.
It rebuilds each image extension matrix with its own sparse rational
reducer and FLINT rational block arithmetic. Its index algorithm uses
unit-pivot Smith elimination over Z/(2^k), with exact bounded integers,
independently of PARI HNF. It checks every degree and both images.
Rational matrix arithmetic can share FLINT with Sage; the independently
implemented parts are the polynomial reconstruction and the index algorithm.
The D5 base relations and image bases are shared, previously independently
audited inputs; this is not a fresh construction of the E7 weight system.
The [integer-only summary](experiments/e7-d4-summary.json) checks the
closed polynomial, all ranks, the older degree-nine overlap and the
completed parabolic inventory.

Published foundations include the rational Borel presentation, Weyl
invariant theory, symmetric polynomials and characteristic maps; see the
context in [Kumar Part I](https://kumar.math.unc.edu/papers/kumar76.pdf) and
[Kumar--Xie Part III](https://kumar.math.unc.edu/papers/kumar86.pdf).
The source-completeness criterion uses the standard Nakayama lemma.
The rational quadric relations (6) are elementary signed-permutation
identities, explicitly proved here by monic reduction.

Integral spin invariant theory is already studied; for example
[Duan, Sections 7--8](https://arxiv.org/abs/1810.03799) describes integral
spin Weyl invariants. Its lattice must be matched with this actual E7
Levi lattice before importing integral generators. The full ten-generator
source certificate, the index equality for this D4 flag, and the completed
64-refinement formula are deductions recorded here. The inspected sources
are not being credited with this exact relative-image formula. This is
not a complete novelty search and does not establish priority.

The next structural question is why the correction indices in (7) agree,
without checking them degree by degree, and how to formulate a comparison
for the 15 remaining E7 parabolics with zero degree-one gap. The same
polynomial already fails for them in degree one. The present proof supplies a
complete benchmark and a finite-module method for testing such a theorem.
It adds no runtime API; library version remains 0.20.0.

## Reproduction

From borel/, with Sage for the primary calculations and native Python
with gmpy2, SymPy, NumPy and python-flint 0.9.0 for the audits:

```sh
sage -python experiments/e7_d4_source.py --max-degree 12
sage -python experiments/e7_d4_source_algebra.py
sage -python experiments/e7_d4_images.py
python experiments/audit_e7_d4_source.py
python experiments/audit_e7_d4_images.py
python experiments/e7_d4_summary_check.py
```

An isolated python-flint installation can be selected with
PYTHON_FLINT_PATH. The image archive contains every extension matrix;
the tower archive contains all rational preimages and module substitutions.
Text hashes normalize CRLF to LF; archive hashes use raw bytes.
