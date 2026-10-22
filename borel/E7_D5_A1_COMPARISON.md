# E7/P6: the complete integral root-to-weight image comparison

This continues the [remaining-parabolic comparison](E7_REMAINING_PARABOLICS.md).
The full multiplicity formula is now proved. The
[complete local-source argument](E7_LOCAL_SOURCES.md) supplies an audited
25-generator root source in all degrees. Its degree-23 correction has
Borel image in the earlier constructed image, which therefore equals the
full root image. The earlier bounded calculations are retained below as
inputs to this completion.

## Statement and relation to the Springer maps

Let G be simply connected of type E7 over C. Use Bourbaki numbering, with
edges 1--3, 3--4, 2--4, 4--5, 5--6, 6--7. Let X be the weight lattice,
Q the root lattice, and I={1,2,3,4,5,7} the retained Levi nodes. Thus
G/P_I=E7/P6 has Levi type D5+A1, complex dimension 42, and 756 Schubert
classes. Put R=Z_(2) and define the full, unscaled Borel images

\[
 M^d=\beta^{P_I}((\operatorname{Sym}_R^d Q_R)^{W_I}),\qquad
 N^d=\beta^{P_I}((\operatorname{Sym}_R^d X_R)^{W_I}).
\]

The target is rational cohomology, with the integral Borel-image lattices
embedded in it. Polynomial degree d means cohomological degree 2d.

**Theorem (computer-assisted).** For every degree d,

\[
                  2N^d\subseteq M^d\subseteq N^d.
\]

More precisely, N^d/M^d=(Z/2)^(m_d), where m_d is given by the complete
table below. Its sum is 366. The first nonzero quotient is
N^2/M^2=Z/2, and all m_d vanish above degree 42. The
[complete summary](experiments/e7-d5a1-complete-summary.json) binds the
source proof, degree-23 witness, orbit norm and independently audited
embedded image calculation.

Since Q and X become equal at odd primes, this also describes the possible
nonunit factors of the integral relative comparison. It does not compute
the absolute cokernel H^(2d)(G/P_I;Z)/N^d.

For the previously established normalized Springer sources,

\[
 \mathcal J_{56}^d=12^dM^d,\qquad
 \mathcal J_{133}^d=36^dN^d.
\]

Multiplication by 3 is invertible on N^d/M^d. Thus the same relative group
describes (J_56^d+J_133^d)/J_56^d. These are consequences of the
[earlier character-lattice classification](IRREDUCIBLE_CHARACTER_LATTICES.md),
not new determinations of the Springer morphisms themselves.

## Actual lattice coordinates and the full weight source

Set h=omega6, w=omega7, z=omega5, and

\[
 (u_1,\ldots,u_5)=(\omega_1,\omega_3-\omega_1,
 \omega_4-\omega_3,\omega_2+\omega_5-\omega_4,
 \omega_5-\omega_2),\qquad y_i=2u_i-h.
\]

Then sum(u_i)=2z. The tuples (h,w,z,u1,...,u4) and
(h,b=2z,p=w-z,u1,...,u4) are integral bases of X and Q respectively.
The D5 factor acts by even signed permutations of the y_i and fixes h,w.
The commuting reflection s7 fixes h,z,u_i and sends w to h-w.

The A4 subgroup on nodes {1,2,3,4} permutes the u_i. On the weight
lattice, replacing the four variable u_i by u_i-u5 has odd index 5.
Centering by 2z/5 therefore proves that the full permutation-invariant
coefficient ring over R is R[h,w,z,c2,c3,c4,c5], with c_j=e_j(u).
On the root lattice, (h,p,u1,...,u5) is itself a basis, and the ordinary
symmetric-polynomial theorem gives R[h,p,b,c2,c3,c4,c5].

For any torsion-free coefficient ring D, the invariants of D[v] under
v -> a-v are D[v(a-v)]. Indeed an invariant of odd leading degree has
leading coefficient equal to its negative, so that coefficient is zero.
For even degree, subtract a multiple of a power of v(a-v) and induct.
No averaging by two is used. Applying this observation gives the exact
coefficient rings after A4 and A1 invariants:

\[
 \begin{array}{ll}
 X:&R[h,z,c_2,c_3,c_4,c_5,t],\quad t=w(h-w),\ e_1(u)=2z,\\
 Q:&R[h,b,c_2,c_3,c_4,c_5,r],\quad r=p(h-b-p),\ e_1(u)=b,
 \end{array}\qquad t=r+hb/2-b^2/4.
\]

These are coefficient integrality tests for rational D5+A1 invariants;
they are not yet their full D5+A1 invariant rings.

There is an additional simplification for the weight source. The
[previously audited D5 ring](E7_D5_REFINEMENT.md) A has exactly the same
formal coefficient ring R[h,z,c2,...,c5], the same y_i=2u_i-h action,
and ten algebra generators of degrees 1,2,4,5,6,8,8,9,12,16.
Since D5 fixes w, its full invariants are A[w]. Since s7 fixes A,

\[
                 (\operatorname{Sym}_R X_R)^{W_I}=A[t].
\]

This reuses the earlier all-degree D5 source theorem as a shared input.
It does not infer integral source completeness from rational generators.

## Rational quotient and independent check

Write E_j=e_j(y_1^2,...,y_5^2) and p5=product(y_i). The rational Levi
invariant algebra is polynomial on h,t,E1,E2,p5,E3,E4, with degrees
1,2,2,4,5,6,8. This is the classical even-signed-permutation invariant
theorem and the commuting A1 norm calculation.

A fresh Weyl orbit of omega7 has 56 weights. In these coordinates they are

\[
 \begin{split}
 &\pm w,\quad\pm(h-w),\\
 &\pm(w-h/2)\pm y_i/2\quad(1\leq i\leq5),\\
 &\pm(h+\textstyle\sum_i\epsilon_i y_i)/4,
       \qquad\prod_i\epsilon_i=1.
 \end{split}
\]

Their power sums in degrees 2,6,8,10,12,14,18 are independent basic E7
invariants: the exact Jacobian is nonzero at the recorded point. The
quadratic sum is 3(E1+3h^2-8t), so eliminate E1=8t-3h^2. Borel's
rational presentation then gives the six-relation quotient in h,t,E2,p5,E3,E4.

The [independent moment audit](experiments/e7-d5a1-moments-audit.json)
derives the actual lattice action and 56-weight orbit afresh. It verifies
every full trace identity by exact evaluation on 333 samples whose
evaluation matrices have full rank modulo 101. Thus these are polynomial
identities, not merely agreement at selected points. A separate sparse
Buchberger audit certifies the 756 standard monomials and top degree 42.

## Full root components and the exponent proof

In each bounded degree the rational classical invariant monomials are
expanded in the exact root coefficient ring above. Repeatedly halve every
integral mod-two relation among their coefficient columns. Stop when the
columns are independent modulo two. Such columns span the full local
invariant component: a primitive lattice of full rational rank equals the
intersection of that rational subspace with the coefficient lattice.
In the root-source JSON certificates the variable named `z` represents
b=2omega5; the weight-source notation in this note uses z=omega5.

The [native FLINT audit](experiments/e7-d5a1-sources-audit.json) reconstructs
the classical polynomials using fresh Newton identities. It checks every
coefficient of the recorded rational basis change, integrality, and full
primitive mod-two rank. The [extension audit](experiments/e7-d5a1-sources-22-audit.json)
verifies full components through degree 22, including the new indecomposables
in degrees 21 and 22.
It does not assert an unbounded algebra presentation.

Let L^d be the image of the explicitly constructed root invariants and
their products; in the bounded source range use the entire full component.
Every such generator is integral and invariant, so L^d is contained in
M^d. The full weight source A[t] gives N^d in every degree. For each
d=1,...,42 the [independent image audit](experiments/e7-d5a1-images-audit.json)
reconstructs the proposed local bases and verifies

\[
       L^d\subseteq N^d,\qquad 2N^d\subseteq L^d.
\]

Both basis-containment directions and primitive generation are checked
over R; a nonzero determinant or rational rank alone would not suffice.
In particular, if T expresses an L-basis in an N-basis, T and 2T^(-1)
must both have odd denominators. This proves the theorem by
2N^d subset L^d subset M^d subset N^d. Degree zero is trivial and all
images vanish above 42.

The original component calculation gave exact multiplicities through 22
and bounds afterwards. The all-degree source theorem now proves L=M in
every degree: all 25 source generators map to L, including the additional
degree-23 generator. Thus m_d=rank(N^d)-rank(T mod 2) is exact throughout.

The complete formula is as follows. Each entry m_d means (Z/2)^(m_d), not
the index in ambient cohomology.

| Polynomial degrees d | Multiplicities m_d in the same order |
| --- | --- |
| 0 through 7 | 0, 0, 1, 1, 1, 2, 3, 4 |
| 8 through 15 | 4, 5, 7, 9, 9, 11, 12, 15 |
| 16 through 22 | 14, 15, 17, 19, 17, 17, 19 |
| 23 through 30 | 19, 17, 14, 17, 15, 14, 10, 12 |
| 31 through 38 | 9, 9, 5, 7, 4, 4, 2, 3 |
| 39 through 42 | 1, 1, 0, 1 |

The new root generators in degrees 21 and 22 have Borel images already
in the previously generated root image in their own degrees. The previous
upper bounds 17 and 19 are therefore the exact multiplicities. This
illustrates why source generation and generation after the Borel map must
be checked separately.

## The root-orbit norm and completion of the source proof

The actual root-lattice character p=w-z has a 32-element D5+A1 orbit.
If z_epsilon=(5h+sum epsilon_i y_i)/4, product epsilon_i=1, its norm is

\[
       F_{32}=\prod_{\epsilon}(t-hz_\epsilon+z_\epsilon^2).
\]

It is integral because it is the product of actual root characters.
Newton identities for the 16-element orbit of z compute its characteristic
polynomial S16(Z). If S16(Z)=A Z+B modulo Z^2-hZ+t, then

\[
                      F_{32}=A^2t+ABh+B^2.
\]

The [independent norm audit](experiments/e7-d5a1-norm-audit.json) checks
the fresh root orbit, spin moment identities, Newton recurrence, quadratic
resultant and all root-coefficient integrality. Modulo two its restriction
to the r-axis is r^16. Together with the degree-23 correction and the
generators through 22, it gives the complete 25-generator root source.

The degree-20 generating list is incomplete:
the component extension finds another necessary generator in each of degrees
21 and 22. Independent mod-two product ranks are 547 versus full rank 548
in degree 21, and 697 versus 698 in degree 22. Thus the degree-20 list,
even with a degree-32 norm added, cannot generate the full source. The
[failed completeness report](experiments/e7-d5a1-source-algebra.json) records
this obstruction. A separate finite-module attempt obtains fiber rank 32
but the wrong degree numerator; it is
[recorded as incomplete](experiments/e7-d5a1-free-source.json), with no
inference that the full source is nonfree. These earlier failed attempts
are superseded by the [initial-monoid proof](E7_LOCAL_SOURCES.md): 25
leading-polynomial witnesses through degree 32 give a monoid whose entire
Hilbert series equals the rational invariant Hilbert series. Explicit
residue ideals retain its holes. Independent product, monoid and Hilbert
checks then prove source completeness by saturation and degreewise
Nakayama. The degree-23 generator's exact source/Borel transcript proves
that its image was already in L. This closes every higher multiplicity.

## Published inputs, deductions and reproduction

Borel's rational presentation, classical D5 invariants, the degree list
for E7 basic invariants, the symmetric-polynomial theorem, and elementary
commutative algebra are established inputs. The full D5 integral algebra
and the normalized Springer source classification are earlier audited
repository deductions reused here. The all-degree exponent comparison on
this particular parabolic, the complete two-local root source and the
full multiplicity formula are computer-assisted deductions. No new literature search or priority claim
is made in this continuation. See the [earlier targeted literature check](E7_REMAINING_PARABOLICS.md#independent-checks-reproduction-and-literature).

From borel/, the primary and independent commands are:

```sh
sage -python experiments/e7_d5a1_sources.py --max-degree 20
python experiments/audit_e7_d5a1_sources.py
sage -python experiments/extend_e7_d5a1_sources.py --max-degree 22
python experiments/audit_e7_d5a1_sources.py --degree-bound 22
sage -python experiments/e7_d5a1_norm.py
python experiments/audit_e7_d5a1_norm.py
sage -python experiments/e7_d5a1_images.py
python experiments/audit_e7_d5a1_moments.py
python experiments/audit_e7_d5a1_images.py
python experiments/audit_e7_d5a1_extra23.py
python experiments/e7_d5a1_initial_source.py
python experiments/audit_e7_initial_sources.py --levi d5a1
python experiments/summarize_e7_d5a1_complete.py
```

The source/norm audits require python-flint (optionally through
PYTHON_FLINT_PATH); the remaining native audits use SymPy and/or gmpy2.
All independent audits use separate arithmetic and source reconstruction
and import no Sage or primary project implementation. The earlier complete
D5 theorem remains an explicitly shared mathematical input.
