# The first E7 half-sum: integral invariant lifting in degree eleven

Status, 21 October 2026: the degree-eleven step is proved. A rank-one
saturation lemma and exact quadratic-form identities extend the structural
argument to every parabolic in the 112-case family through degree 11.
Two independent implementations check the coefficient inputs. The
uniform all-degree structural proof remains incomplete. This does not
change the already established all-degree multiplicity formulas or claim
literature priority.

## 1. Objects and theorem

Use the [quadratic-derivation note](E7_QUADRATIC_DERIVATION.md): G is
simply connected complex E7, X=X(T), Q its root lattice, R=Z_(2), and
I denotes retained Bourbaki nodes. Put

\[
 M_I=\beta^{P_I}((\operatorname{Sym}_R Q_R)^{W_I}),\qquad
 N_I=\beta^{P_I}((\operatorname{Sym}_R X_R)^{W_I}).
\]

Polynomial degree d means cohomological degree 2d. The 112-case family
consists of I omitting at least one of nodes 2,5,7. For an omitted such
node i, write h=beta(omega_i) and
C_I=im(M_I -> N_I/2N_I). The common quadratic lift already proves
2N_I contained in M_I; in particular C_I=M_I/2N_I, not M_I/2M_I.

**Theorem.** The intrinsic source derivative

\[
 \partial(\textstyle\sum a_j\omega_j)=a_2+a_5+a_7\pmod2
\]

respects the actual integral invariant-source Borel kernel through
degree **11**, for all 112 cases. Consequently, in those degrees,

\[
 N_I/2N_I=C_I\oplus\overline h C_I,\qquad
 N_I/M_I\simeq C_I[-1],\qquad
 m_{I,d}=[q^d]\,\frac{qP_{W(E_7)}(q)}{(1+q)P_{W_I}(q)}.
 \tag{1}
\]

Here C_I[-1]_d=C_I^(d-1). The derivative and the indicated multiplication
give this splitting without importing the known relative Smith counts.
On the 108 good type A Levi cases, classical common-fiber descent also
transfers the existing full-flag derivative check through degree **13**.
The other four cases are E6, the two D5 diagrams, and D4 inside the E6
Levi. No uniform claim beyond degree 11 is made for these four.

For the Springer/Borel composite in the Kumar papers,
xi_rho^(P_I)=beta^(P_I) after (theta_rho^*)^(-1), the earlier literal-source
theorem gives J_(56,I)^d=12^d M_I^d and J_(133,I)^d=36^d N_I^d. Thus
(1) supplies an explanation of the relative arithmetic of those same
maps. It is not a new Springer morphism or a determination of H/N_I.

## 2. The cubic and quintic come from the quadratic form

Let V=sum_(j=1)^6 R alpha_j, w=omega7 and C6 the E6 Cartan matrix.
Then X_R=V plus Rw, Q_R=V plus R(2w), and the common root quadratic is

\[
 B=-\frac32 a^tC_6^{-1}a,\qquad
 \mathcal D=\operatorname{Sym}_R(X_R)/(w^2-B)
             =\operatorname{Sym}_R(V)\oplus w\operatorname{Sym}_R(V).
\]

Write a_j=alpha_j as polynomial variables. From bar(B) alone define
two polynomials over F2:

\[
 G_3=\sum_{j=1}^6 a_j^2\frac{\partial\overline B}{\partial a_j},
 \qquad
 G_5=\sum_{j=1}^6 a_j^4\frac{\partial\overline B}{\partial a_j}.     \tag{2}
\]

These formal Frobenius derivations commute with linear changes of
variables over F2. Since B is W(E6)-invariant, G3 and G5 are modular
invariants. That observation does **not** assert they lift to integral
invariants. Both have ten monomials. Denote their coefficient-0-or-1
integral polynomial lifts by g3,g5.

Fresh sixth, eighth and tenth moments of the 56 weights, reduced in
D and divided by coefficient contents, give primitive full-Weyl-invariant
classes Fj=Aj+w C_(j-1), j=6,8,10. Exact coefficient identities are

\[
 \overline F_6=G_3^2,\qquad \overline F_8=G_3G_5,\qquad
 C_5\in4\operatorname{Sym}_R^5 V,\quad
 C_7\in2\operatorname{Sym}_R^7 V,\quad C_9\in4\operatorname{Sym}_R^9 V,
 \tag{3}
\]
\[
 \overline{C_5/4}=G_5+\overline B G_3,\qquad
 \overline{C_7/2}=\overline B G_5+\overline B^2G_3.             \tag{4}
\]

They explain the first modular coefficient syzygy and its half-sum:

\[
 F_6g_5+F_8g_3\equiv0\pmod2,\qquad
 H_{11}=\frac{F_6g_5+F_8g_3}{2}\in\mathcal D_{11}.           \tag{5}
\]

Let K_d be the full integral characteristic kernel in D_d and put

\[
 L_{11}=F_6\mathcal D_5\oplus F_8\mathcal D_3\oplus F_{10}\mathcal D_1.
\]

The 462 initial columns have rank 461 modulo two. Appending H11 raises
that rank to 462, the rational Borel-kernel dimension. All columns are in
the rational kernel, and 2H11 lies in L11. Thus the augmented span is
primitive and is the **entire integral kernel**:

\[
 K_{11}=L_{11}+R H_{11},\qquad K_{11}/L_{11}\simeq R/2R.      \tag{6}
\]

The initial 462 columns are rationally independent: H11 is already in
their rational span, whose rank is at least the augmented mod-two rank.
Consequently the coefficient expression in L11 is unique. Its sole
mod-two coefficient syzygy is (G5,G3,0). These statements are proved
from primitive kernel columns, not a relative image count.

## 3. A rank-one saturation lifting lemma

The following argument applies to any graded algebra of this form,
with the invariance, uniqueness, saturation and coefficient identities
just stated. It is not an exactness assertion for modular invariants.

Let Gamma preserve V, fix w, preserve B, and fix the classes F6,F8,F10.
Set A=(Sym_R V)^Gamma, so D^Gamma=A plus w A. Assume the kernel is
Gamma-stable. The normal-form operator p(a+wb)=b commutes with Gamma
over R; modulo two it is the derivation partial.

**Rank-one lemma.** Under (3)--(6), partial respects the integral
invariant kernel through degree 11, provided it does so in the preceding
degrees. No lifting of all modular Gamma-invariants is required.

**Proof.** For f in K11^Gamma, 2f has a unique expression

\[
                 2f=F_6u_6+F_8u_8+F_{10}u_{10}.             \tag{7}
\]

Each u_j is Gamma-invariant, by uniqueness and invariance of Fj.
If f belongs to L11, its unique coefficients are already integral
invariants; (3) and the mod-two Leibniz rule prove the assertion.

Otherwise (6) shows that the coefficients in (7) reduce to (G5,G3,0).
Write u_j=x_j+w y_j with x_j,y_j in A. These reductions have no w
component, so every y_j is even in A; also x_10 is even in A. Subtract
the integral invariant kernel element
sum_j Fj*w(y_j/2)+F10*(x_10/2). It suffices to handle the remaining
integral invariant element

\[
 f'=\frac{F_6x_6+F_8x_8}{2},\qquad
 \overline x_6=G_5,\quad\overline x_8=G_3.                   \tag{8}
\]

These particular integral invariant lifts exist because of f; the
argument has not assumed a general modular lifting theorem. Taking
the odd-w coefficient and using (3)--(4) gives

\[
 \overline{p(f')}=
   \overline{C_5/2}\,G_5+\overline{C_7/2}\,G_3
   =\overline B^2\,\overline F_6+\overline B\,\overline F_8.  \tag{9}
\]

The right side has the integral invariant kernel lift B^2 F6+B F8.
Thus p(f') minus that lift is divisible by two in D. Both terms are
Gamma-invariant, so division by two stays in D^Gamma. The previously
subtracted invariant kernel element has derivative in the lower kernel
modulo two by (3). This proves

\[
 p(K_{11}^{\Gamma})\subseteq K_{10}^{\Gamma}+2\mathcal D_{10}^{\Gamma}.
\]

This is precisely the kernel condition in the quadratic-derivation
criterion. The proof uses neither averaging nor invariants commuting
with reduction modulo two.

Apply the lemma to every subgroup Gamma of W(E6), in particular all
64 retained subsets I contained in {1,...,6}. The preceding degree-10
argument supplies the lower-degree hypothesis. This proves their new
degree-eleven step.

## 4. Why the half-sum is absent on the four exceptional Levis

There is a second useful explanation for the four diagrams that the
good-Levi transfer misses. If K11^WI contained an element outside L11,
equation (8) would give an integral invariant cubic in Sym_R(V) reducing
to G3. Such a lift does not exist in these four cases.

For E6, there is no rational invariant cubic. For D5 and D4, the
nontrivial reflection factors also have no basic invariant of degree
three. Consequently every rational invariant cubic belongs to the ideal
generated by the fixed linear forms. This assertion includes products
of a fixed linear form with the invariant quadratic.

The full fixed linear lattice is primitive in V_R: E6 root and weight
lattices agree at two since det C6=3. It is spanned by the E6 fundamental
weights at the omitted nodes. Therefore its polynomial ideal is
saturated. The reduction of any integral invariant cubic must belong
to the corresponding mod-two linear ideal. Direct restriction gives:

| retained I | Levi | reduced fixed linear forms | restriction of G3 |
| --- | --- | --- | --- |
| 1,2,3,4,5,6 | E6 | none | nonzero |
| 1,2,3,4,5 | D5 | a2+a5 | nonzero |
| 2,3,4,5,6 | D5 | a2+a3 | nonzero |
| 2,3,4,5 | D4 | a2+a3, a2+a5 | a4*a5*(a4+a5), nonzero |

Hence G3 does not lift, and in these four cases

\[
                  K_{11}^{W_I}=L_{11}^{W_I}.                \tag{10}
\]

Thus the first full-flag half-sum does not contribute a new invariant
kernel element on these four Levis. This is an invariant-lifting
obstruction that helps the proof, rather than a missing assumption.
The rank-one lemma above also handles cases where such a lift does exist.

## 5. Descent to the other parabolics

Here is the second reusable lemma. Suppose the actual weight image N_I
embeds into the full-flag image N_empty and the quotient is R-torsion-free.
Suppose the source derivative preserves the reduced actual invariant
source. If it descends to N_empty/2N_empty through degree d, it descends
to N_I/2N_I through that degree: the latter injects into the former, and
the diagram with source differentiation commutes. In particular a
derivative of an invariant kernel element vanishes in N_I/2N_I. This
argument concerns actual images, not their images in ambient H/2H.

For a product-of-type-A Levi whose Cartan column matrix C[:,I] has full
column rank modulo two, the previously proved common polynomial fiber
basis gives N_empty as a direct sum of shifted copies of N_I, with the
unit fiber vector giving the stated inclusion. The quotient is therefore
torsion-free. This is the classical common-fiber/Leray--Hirsch argument in
the [full-flag proof](E7_FULL_FLAG_COMPARISON.md#4-transferring-the-result-to-the-full-flag),
used here to transfer the derivative, not numerical multiplicities.

An independent graph and bit-rank enumeration finds 108 such proper
retained sets. They all belong to the 112-case family; its other four
are exactly the diagrams in the table above. Their source derivative
preserves actual invariants by the fixed-odd-weight split. The full-flag
degree-eleven assertion is already a special case of the rank-one lemma;
common-fiber descent proves it for all 108. Combining with either section
3 or section 4 proves (1) for all 112.

The earlier independent full-flag certificate reaches degree 13. The same
descent proves the derivative criterion through 13 on those 108 cases.
No image table or exactness of a modular invariant functor is needed.

## 6. Independent checks, attribution and the next obstruction

The [primary calculation](experiments/e7-degree11-lifting.json) starts
from a fresh 56-weight Weyl orbit. The cubic and quintic are constructed
from B by (2), not extracted from the previous saturation witness. It
checks (3)--(5), the 461-to-462 primitive kernel rank increase, the single
coefficient syzygy, the four cubic obstructions and the 108/112 inventory.

The [independent audit](experiments/e7-degree11-lifting-audit.json)
constructs the weights from 27 E6 forms and integer multinomial moments.
It rebuilds B,G3,G5 and H11 exactly over the integers. Its own bit ranks
verify the primitive kernel; its own linear-substitution algorithm checks
the four obstructions, independently of the primary Sage Groebner
reducer. SymPy is used only for the six-dimensional rational matrix
inverse. No relative Smith data or old halving witness are read. The
Cartan matrix, rational Borel degree list and algebraic definitions are
shared mathematical inputs. The rank-one and descent lemmas are the
written arguments just given, not assertions certified by a matrix test.

Published ingredients remain rational Borel presentations, rational Weyl
invariant degrees and the type A common-fiber construction. The new
repository deduction explains the first integral correction and proves
the structural criterion in the stated range. No new literature-priority
search was undertaken for these lemmas or identities in this step.

The previous all-degree E7 Smith formulas remain established by their
case/transfer proofs. They are not used to prove this extension. The
earlier presentation of C_I using only F6,F8,F10 remains truncated at
degree 10; the new kernel can add constant-sector relations in degree 11.
Quadratic splitting does not assert that this older three-relation
presentation holds unchanged.

The degree-by-degree continuation is now deferred in favor of the
[general root-data question](CENTRAL_PARITY_KERNEL.md), which has an
all-rank type C theorem and an actual Spin7 counterexample.

**Deferred E7 target:** determine the degree-12 correction module and its Levi
invariant behavior. The full-flag saturation index valuation is eight
there, so the rank-one lemma alone does not cover it. An invariant-compatible
description of the new relations is required before extending the proof
for E6, the two D5 cases and D4. The uniform all-degree structural theorem
and a closed presentation of C_I remain unresolved here.

## Reproduction

From borel/:

```sh
PYTHONPATH=src sage -python experiments/e7_degree11_lifting.py
python experiments/audit_e7_degree11_lifting.py
```

The native audit uses NumPy and SymPy, not Sage or the library runtime.
Both reports bind their sources and the explicit coefficient archive.
