# The remaining E7 parabolics: reductions and the D6 calculation

This note continues the [112-case result](E7_D4_REFINEMENT.md). It separates
classical inputs, exact repository deductions and unresolved cases. No
literature-priority claim is made. The preceding D4 work was merged into
master through [PR 9](https://github.com/sean-rogers-t/Sage-programs/pull/9),
merge commit `2cee8aa`.

## Objects and scope

Let G be simply connected of type E7 over C, with Bourbaki edges
1--3, 3--4, 2--4, 4--5, 5--6, 6--7. Write X for its weight lattice,
Q for its root lattice, W for its Weyl group and R=Z_(2). A retained
set I specifies the simple roots of the Levi of P_I. Define the full images

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^d Q)^{W_I}),
 \qquad
 N_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^d X)^{W_I}).
\]

Polynomial degree d means cohomological degree 2d. We compare N_I/M_I,
not the ambient cohomology/image quotient. Since X/Q=Z/2, the images
coincide after localization at any odd prime. Thus a complete 2-local
comparison determines this integral relative cokernel.
Here [n]_q=1+q+...+q^(n-1), and P_H(q) is the Coxeter length polynomial
of a finite Weyl group H.

For the normalized minimal 56-dimensional and adjoint 133-dimensional
Springer maps, the earlier source classification gives

\[
 J_{56,I}^d=12^d M_I^d,\qquad J_{133,I}^d=36^d N_I^d.
\]

Consequently (J_56+J_133)/J_56 is isomorphic degreewise to N_I/M_I:
after dividing the common 12^d, multiplication by 3^d acts invertibly
on the 2-primary group N_I/M_I. These calculations concern the integral
polynomial-character/Borel images of our Springer maps.

## Why the old polynomial stops, and four conjugacy reductions

The omitted fundamental weights generate X^(W_I). The inverse E7 Cartan
matrix shows that precisely omega_2, omega_5 and omega_7 have nontrivial
class in X/Q. Therefore the 15 proper retained sets containing {2,5,7}
have X^(W_I)=Q^(W_I), and hence N_I^1=M_I^1. The previous polynomial
q P_W/((1+q)P_WI) has coefficient one in degree one and cannot apply.
Degree-one equality does not imply equality in higher degrees.

Weyl conjugacy of Levi root systems identifies the graded relative image
comparisons. Indeed W acts integrally on H*(G/B;Z), as realized through
G/T, and the characteristic map is W-equivariant. Its action sends both
the Q-invariant and X-invariant images for I to those for the conjugate
system. This is a formal use of the classical characteristic map and
Weyl action, not a new geometric theorem.

The following words identify four pairs. Entries are simple reflections
in **application order**: apply each reflection to the current root vector
from left to right.

| First retained set | Conjugate retained set | Word |
| --- | --- | --- |
| {1,2,5,7} | {2,3,5,7} | (3,1) |
| {2,4,5,7} | {2,5,6,7} | (6,5,4,2,7,6,5,4) |
| {1,2,4,5,7} | {1,2,5,6,7} | (6,5,4,2,7,6,5,4) |
| {1,2,4,5,7} | {2,3,5,6,7} | (6,5,4,2,7,6,5,4,3,1) |

The [root-mask search](experiments/e7-remaining-conjugacy.json) constructs
all 126 roots and checks these identities on entire Levi root systems.
The [independent audit](experiments/e7-remaining-conjugacy-audit.json) uses
root-basis reflection matrices instead and checks unimodular equality of
the two Levi root lattices. This reduces 15 diagrams to **11** cases.
Their different root-system types prevent further Weyl identifications
between these representatives.

| Representative I | Levi type | Complex dimension | Schubert rank |
| --- | --- | --- | --- |
| {2,5,7} | A1^3 | 60 | 362,880 |
| {1,2,5,7} | A1^4 | 59 | 181,440 |
| {2,4,5,7} | A3+A1 | 56 | 60,480 |
| {1,2,3,5,7} | A2+A1^3 | 57 | 60,480 |
| {1,2,4,5,7} | A3+A1^2 | 55 | 30,240 |
| {2,3,4,5,7} | D4+A1 | 50 | 7,560 |
| {2,4,5,6,7} | A5 | 48 | 4,032 |
| {1,2,3,4,5,7} | D5+A1 | 42 | 756 |
| {1,2,3,5,6,7} | A3+A2+A1 | 53 | 10,080 |
| {1,2,4,5,6,7} | A5+A1 | 47 | 2,016 |
| {2,3,4,5,6,7} | D6 | 33 | 126 |

Conjugacy alone does not calculate any of their higher multiplicities.

## D6 theorem: complete relative multiplicities

For I={2,3,4,5,6,7}, in every degree,

\[
 N_I^d/M_I^d\simeq(\mathbb Z/2)^{m_d},\qquad
 \sum_d m_dq^d
 =q^8(1+q)\bigl([10]_{q^2}+q^4[6]_{q^2}+q^{12}\bigr).       \tag{1}
\]

There are **34** nonunit factors in total. Explicitly m_d is 1 in
degrees 8--11 and 24--27, 2 in degrees 12--19 and 22--23, 3 in
degrees 20--21, and zero otherwise. This is a computer-assisted all-degree
deduction for this one flag, with separate native arithmetic checks.
The new scope is **113** proper E7 parabolics. The remaining **14**
diagrams reduce to **ten** distinct cases under the conjugacies above;
their complete formulas are not proved here.

The sources are certified componentwise for N through degree 33 and M
through degree 27. Thirteen integral N generators have primitive products
spanning every required component through degree 33. This does not assert
their unbounded algebra completeness. In degrees 28--33, products of actual
root invariants and already certified lower root images fill the complete
weight image. If L_d is this computed lower root image, then
L_d contained in M_I^d contained in N_I^d and L_d=N_I^d prove equality.
Above degree 33 the cohomology vanishes, so all image degrees are exhausted.

## D6: actual integral coordinate lattices

Take I={2,3,4,5,6,7}, so G/P_I is conventionally E7/P1. Put h=omega_1,
z=omega_2 and

\[
\begin{split}
 u_1&=\omega_7,&u_2&=\omega_6-\omega_7,
 &u_3&=\omega_5-\omega_6,\\
 u_4&=\omega_4-\omega_5,&u_5&=\omega_2+\omega_3-\omega_4,
 &u_6&=\omega_1+\omega_2-\omega_3.
\end{split}
\]

Then sum u_i=h+2z, and (z,u_1,...,u_6) is an integral basis of X.
The A5 subgroup on nodes {3,4,5,6,7} fixes z and permutes the u_i.
For v_i=u_i-z and a=2z, (a,v_1,...,v_6) is an integral basis of Q,
with sum v_i=h-2a. Thus the full A5 invariant coefficient rings are

\[
 (\operatorname{Sym}X)^{S_6}=\mathbb Z[z,e_1(u),\ldots,e_6(u)],
 \qquad
 (\operatorname{Sym}Q)^{S_6}=\mathbb Z[a,e_1(v),\ldots,e_6(v)].
\]

This uses the integral symmetric-polynomial theorem on an actual permutation
lattice. No splitting by the even index six is assumed.

Set y_i=2u_i-h=2v_i-a-e_1(v). The full D6 group acts by even signed
permutations of the six y_i and fixes h. Its rational invariant algebra
has generators h, e_1(y^2), e_2(y^2), e_3(y^2), product y_i,
e_4(y^2), e_5(y^2), in degrees (1,2,4,6,6,8,10).
These rational generators are not assumed to generate either integral source.

In each degree, form their complete rational monomial basis in the relevant
integral S6 coefficient ring. An integral polynomial matrix A with
F C=A, det C nonzero, and full column rank modulo two gives the entire
2-local invariant component: its columns are invariant, have the full
rational rank, and are primitive in the coefficient lattice. This is an
elementary lattice-saturation criterion. It certifies the **full** source,
rather than a selected subalgebra.

## D6: rational quotient and finite verification method

A fresh orbit of omega_7 gives the 56 weights in the following coordinates:
32 half-spin forms sum epsilon_i y_i/4 with product epsilon_i=1, and
24 forms (+/-h +/-y_i)/2. Their power sums in degrees
2,6,8,10,12,14,18 have nonzero Jacobian at an explicitly recorded rational
point. These are the published basic invariant degrees of E7. Therefore
the power sums generate the rational W(E7) invariant algebra: algebraic
independence and the matching Hilbert series give equality.

The degree-two relation gives e_1(y^2)=-2h^2. The other six moments give
a quotient in six generators of weights (1,4,6,6,8,10). Its standard
monomials have total dimension 126 and top degree 33, with Hilbert series

\[
 P_{E7/P1}(q)=\frac{[12]_q[14]_q[18]_q}{[4]_q[6]_q}
 =(1+q^4+q^8)(1+q^6+q^{12})[14]_q.
\]

This rational Borel presentation and Weyl-degree dimension formula are
classical inputs. In particular, a rational presentation in generators
of degrees 1,4,6 with relations of degrees 12,14,18 is already available
by setting the quantum parameter to zero in Chaput--Perrin, Proposition 5.6.

For any proposed local image basis H, the verification checks that H^-1 A
has odd denominators and full row rank modulo two. This proves equality
with the image lattice by an invertible minor over R. For T=H_N^-1 H_M,
integrality of T and 2T^-1 proves that the relative factors are at most 2;
their multiplicity is dim(H^2d)-rank(T mod 2). Every step uses exact arithmetic.

The [summary](experiments/e7-d6-summary.json) binds all checks and verifies
formula (1) coefficientwise. The degree-eight certificate also provides an
explicit integral weight invariant whose Borel class has order two modulo
the root image. In the rational quotient coordinates its normal form is

\[
 \frac{25}{256}v_4^2+\frac{1871}{768}h^2v_6
 -\frac{83}{6}h^2p_6+\frac{285}{256}v_8.
\]

Here v_4=e_2(y^2), v_6=e_3(y^2), p_6=product y_i and v_8=e_4(y^2).
The denominators belong to this rational coordinate presentation; its
integral X-polynomial lift and the coordinates proving nonmembership in M
and membership of twice the class are recorded and independently checked.

## D5+A1: an explicit first-gap proof

For I={1,2,3,4,5,7}, conventionally E7/P6, the first nonzero relative
degree is **two**, with N_I^2/M_I^2=Z/2. This is a partial result for this
flag: its higher multiplicities remain unresolved here. In particular the
condition X^(W_I) contained in Q alone cannot distinguish the first-gap
degrees two and eight in these two examples.

Here is a direct proof, independently corroborated by full invariant kernels.
Put h=omega_6, w=omega_7, z=omega_5 and

\[
 (u_1,\ldots,u_5)=(\omega_1,\omega_3-\omega_1,
 \omega_4-\omega_3,\omega_2+\omega_5-\omega_4,
 \omega_5-\omega_2).
\]

Their sum is 2z. The tuples (h,w,z,u_1,...,u_4) and
(h,b=2z,p=w-z,u_1,...,u_4) are integral bases of X and Q, respectively.
The D5 factor acts by even signed permutations of y_i=2u_i-h and
fixes h,w; s_7 fixes h,y_i,z and sends w to h-w. Thus, with
E=sum y_i^2 and t=w(h-w), the rational degree-two invariant space
has basis (h^2,E,t). Let c_2=e_2(u) and r=p(h-b-p). Then

\[
 E=16z^2-8c_2-8hz+5h^2
   =4b^2-8c_2-4hb+5h^2,\qquad
 t=r+hb/2-b^2/4.
\]

The following are bases of the **full** degree-two invariant sources over R:

\[
 A_X^2=R\{h^2,(E-5h^2)/8,t\},\qquad
 A_Q^2=R\{h^2,(E-5h^2)/4,4t\}.                         \tag{2}
\]

All three candidates in each case are integral in the displayed actual
lattice basis. Their coefficient columns are primitive modulo two. For X,
the second candidate is 2z^2-c_2-hz, whose u_1^2 coefficient distinguishes
it from h^2 and t; t has a w^2 term. For Q, the latter two candidates
reduce to b^2+hb and b^2, and are independent of h^2. They have the full
rational rank three, so primitivity proves completeness of (2).

The unique rational E7 degree-two relation is

\[
 E+3h^2-8t=0.                                          \tag{3}
\]

A fresh 56-weight trace is exactly three times its left side. Substituting
(3) in (2) gives full images R{h^2,t} and R{h^2,2t}. The two remaining
rational classes are independent because the only degree-two relation
is (3). Therefore the relative factors are (1,2). Constants and degree-one
characters agree, proving the assertion about the first degree.

The [Sage calculation](experiments/e7-d5a1-degree-two.json) derives full
integral kernels from all six reflection equations in both lattices and
uses an independently constructed 126-root quadratic. The
[native audit](experiments/e7-d5a1-degree-two-audit.json) verifies the
explicit bases (2), their primitive coefficient ranks, the fresh 56-weight
relation (3), and every image basis and relative factor. This is a deduction
from classical invariant theory and Borel's presentation, with no novelty
claim and no assertion about degrees above two.

## Independent checks, reproduction and literature

The primary computation uses Sage/PARI. The independent source verifier uses
native Python/FLINT, reconstructs the Newton polynomials and checks every
exact rational invariant identity, full-rank rational basis change, and primitive
mod-two coefficient rank. The independent moment verifier uses a fresh E7
orbit, a different torus-basis implementation and 120 unisolvent exact
evaluation samples; independent mod-101 ranks prove those samples determine
each full polynomial identity. Its Jacobian check proves completeness of the
rational basic invariants. The image verifier uses native Python/gmpy2,
a separate parser, sparse rational reduction, Buchberger checks, its own
matrix inverses and bit ranks. It reconstructs every image column from the
full-source certificates. No independent verifier imports Sage or project code.
FLINT is shared as generic exact arithmetic by the source verifier, not as
a project invariant-source implementation.

Reports and exact certificates:

- [Full source components](experiments/e7-d6-sources.json),
  [independent source audit](experiments/e7-d6-sources-audit.json).
- [Image comparisons](experiments/e7-d6-images.json),
  [independent image audit](experiments/e7-d6-images-audit.json).
- [Fresh-weight/moment audit](experiments/e7-d6-moments-audit.json).
- [Exact source seed through degree 20](experiments/e7-d6-source-seed20.json.gz),
  [full source certificates](experiments/e7-d6-sources-certificates.json.gz),
  [image/quotient certificates](experiments/e7-d6-images-certificates.json.gz).

From borel/, the primary commands are:

```sh
sage -python experiments/extend_e7_d6_sources.py
python experiments/pack_e7_d6_sources.py
sage -python experiments/e7_d6_images.py --max-degree 33
python experiments/audit_e7_d6_sources.py
python experiments/audit_e7_d6_moments.py
python experiments/audit_e7_d6_images.py
python experiments/e7_remaining_conjugacy.py
python experiments/audit_e7_remaining_conjugacy.py
sage -python experiments/e7_d5a1_degree_two.py
python experiments/audit_e7_d5a1_degree_two.py
```

The source audit requires python-flint; PYTHON_FLINT_PATH can point to an
isolated installation. The other independent audits require SymPy/gmpy2.
The slower experiments/e7_d6_sources.py can regenerate an arbitrary finite
source bound from scratch. The faster extension uses a certificate seed,
but the native source audit rechecks every seed identity and completeness
condition. Higher expanded source polynomials are omitted from the compact
certificate because F C specifies them uniquely. The native audit reconstructs
every coefficient, checks integrality and verifies full primitive rank.

The literature comparison on 21 October 2026 has the following scope:

- [Kumar--Rogers, Part I](https://kumar.math.unc.edu/papers/kumar76.pdf)
  and the existing [source classification](IRREDUCIBLE_CHARACTER_LATTICES.md)
  supply the previously audited Springer/Borel construction and lattices;
  the earlier paper audit is reused here.
- [Chaput--Perrin, *On the quantum cohomology of adjoint varieties*,
  Proposition 5.6](https://arxiv.org/pdf/0904.4824): the rational cohomology
  presentation of E7/P1 is known. It is not being claimed here as new.
- [Duan, *The characteristic classes and Weyl invariants of Spinor groups*](https://arxiv.org/abs/1810.03799):
  integral spin Weyl-invariant rings are established background. The actual
  E7 root/weight Levi lattices and their Borel images still need to be matched;
  the computation above verifies them directly rather than assuming that
  abstract Levi type alone determines the integral source.

The relative root/weight Smith multiplicities are the finite deductions
being investigated here. No exact matching formula was identified in this
targeted check. That does not establish novelty or replace a comprehensive
literature audit. A structural explanation and formulas for the unresolved
representatives remain separate goals.

## D5+A1 continuation

The [subsequent D5+A1 note](E7_D5_A1_COMPARISON.md) extends the degree-two
argument above. It proves 2N^d subset M^d in every degree on E7/P6 and
determines full relative multiplicities through degree 22. Degree 41 also
has equality by containment. The full weight source is the previously
audited D5 algebra A with the variable t=w(h-w). Complete bounded root
components, a fresh 56-weight quotient calculation and separate native
arithmetic audits prove the stated scope. Necessary new root generators
occur in degrees 21 and 22; a degree-20 cutoff is not an unbounded source
presentation. Higher multiplicities are still bounds. Thus 113 parabolics
remain complete, and 14 diagrams still lack full formulas.

## Next proof targets

The next target is to finish the multiplicity polynomial for D5+A1,
retained {1,2,3,4,5,7}, with 756 Schubert classes. The continuation proves
the exponent and degrees through 22. A complete higher root source, or a
proof that a smaller root source has the same Borel image, is still needed.
Invariants and quotient images need not commute at two; averaging by two
is unavailable here.

After that, A5+A1 (rank 2,016) and A5 (rank 4,032) provide smaller
intermediaries for the other cases. A common integral relative basis would
transfer formulas to finer parabolics, but such a basis must be proved;
rational Leray--Hirsch alone does not compare these integral image lattices.
The conjugacy reductions above can be used unconditionally.

The general questions remain: does 2N_I lie in M_I for all proper E7
parabolics, and what diagram/lattice data controls the first nonzero
relative degree and its multiplicities when X^(W_I) lies in Q? The
112-case theorem, the D6 computation and the D5+A1 exponent theorem do not
prove an answer for all remaining diagrams. No uniform conjectural
polynomial is asserted from these base cases.
