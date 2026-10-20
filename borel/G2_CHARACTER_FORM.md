# G2: literal character forms and complete maximal-parabolic images

For the seven-dimensional representation V(omega1), the complete torus
character source is

\[
\boxed{\mathcal C_7=\mathbb Z[6\omega_1,2\omega_2].}
\]

This is an all-degree identification, including inhomogeneous polynomials.
Put x=6omega1 and y=2omega2. For the standard maximal parabolics, omitting
node 1 or 2 respectively, the integral Levi-character forms are

\[
\boxed{\mathcal A_{1,\mathbb Z}=\mathbb Z[x,q],\quad q=y(x-y),\qquad
       \mathcal A_{2,\mathbb Z}=\mathbb Z[y,p],\quad p=x(3y-x).}
\]

These source algebras are polynomial rings with generator degrees 1 and 2.
The [general criterion and API](CHARACTER_SOURCES.md) explain why the
character requirement differs from simply making Borel images integral.

## Conventions and published inputs

G is the complex group of type G2 (whose root and weight lattices coincide).
The representation has highest weight omega1 and dimension 7. Simple roots
in fundamental-weight coordinates are alpha1=2omega1-omega2 and
alpha2=-3omega1+2omega2. Polynomial degree d means H^(2d).

The normalized Springer map is the trace projection, with derivative the
identity at 1. The character condition means an integral virtual character
of the actual Levi, not a positivity requirement.

The Springer/Borel construction is published in
[Kumar, Part I, Definition 4 and Theorem 5](https://kumar.math.unc.edu/papers/kumar76.pdf).
The G2 Springer formulas and complex generator images are published in
[Kumar--Xie, Part III, Lemma 5.1 and Theorem 5.2](https://kumar.math.unc.edu/papers/kumar86.pdf).
Our torus coordinates are t1=e^(omega1), t2=e^(omega2); their Section 5
coordinates are (t1,t2/t1). The calculation checks agreement after this
change, rather than identifying the two coordinate conventions.

The integral intersections, source presentations, image presentations and
tables below are deductions recorded here with exact checks. The general
criterion is an elementary saturation argument. Neither its formulation
nor these deductions are claimed to be new in the literature. A focused
source check does not establish priority.

## Full torus intersection

The weights are

\[
(1,0),\ (-1,1),\ (2,-1),\ (0,0),\
(-2,1),\ (1,-1),\ (-1,0),
\]

each with multiplicity one. Their trace Gram matrix in the simple-coroot
basis is

\[
M=\begin{pmatrix}12&-6\\-6&4\end{pmatrix}.
\]

If Theta_i=theta*(omega_i), trace projection gives M Theta=b, where
b_i=sum_mu mu_i e^mu. In particular the integral pullbacks of x,y are

\[
\begin{aligned}
B_1=\theta^*x={}&t_1^2/t_2+2t_1-t_1/t_2+t_2/t_1-2/t_1-t_2/t_1^2,\\
B_2=\theta^*y={}&t_1-t_1/t_2+t_2/t_1-1/t_1.
\end{aligned}
\]

The coefficient rows at Laurent exponents (2,-1) and (1,0) are (1,0)
and (2,1). Their determinant is 1. Hence a rational linear combination
a B1+b B2 is integral only if a and b are integers. The exact linear
lattice is therefore Lambda7=6Z omega1 + 2Z omega2.

In lexicographic exponent order, the leading exponents of B1,B2 are
(2,-1),(1,0), with coefficients 1,1 and determinant 1. This survives
over every prime field. Distinct polynomial monomials have distinct leading
Laurent exponents, proving algebraic independence modulo every prime.
The [denominator criterion](CHARACTER_SOURCES.md#exact-criterion) now proves
C7=Z[x,y] in every degree.

This lattice can also be written
6Z alpha1 + 2Z alpha2. Indeed their coefficient columns form M, with
determinant 12, the same index as 6omega1,2omega2. Thus this example has
different factors for the two root lengths; a scalar multiple of the root
lattice is not the correct general pattern.

## Integral Levi invariants

The simple reflections act by

\[
s_1(x)=3y-x,\quad s_1(y)=y,\qquad
s_2(x)=x,\quad s_2(y)=x-y.
\]

For P1 the retained reflection is s2, and for P2 it is s1.
The formulas for A1 and A2 follow from the elementary identity

\[
\mathbb Z[a,b]^{\,b\mapsto ca-b}
=\mathbb Z[a,b(ca-b)]\qquad(c\in\mathbb Z).
\]

For completeness, write an invariant as a polynomial in b over Z[a].
Its highest b-degree is even: an odd highest degree would force twice
its leading coefficient to be zero. Subtract a suitable integer multiple
of (b(ca-b))^k to cancel that highest term, and induct. The leading
b-degrees also prove algebraic independence of the displayed generators.
This establishes equality over Z, without averaging or division by 2.

Equivariance and injectivity of theta* identify these invariant rings with
the literal integral Levi-character forms. Their grading is a conclusion
of the torus proof.

## Image rings, integral kernels and complete tables

For each space let h=beta(omega_i) be its positive degree-one Schubert
class. The root/coroot Chevalley calculation gives, in the unique integral
Schubert basis sigma_d in each degree 0,...,5,

\[
\begin{array}{c|rrrrrr}
d&0&1&2&3&4&5\\ \hline
h^d\text{ on }G_2/P_1&1&1&1&2&2&2\\
h^d\text{ on }G_2/P_2&1&1&3&6&18&18
\end{array}
\]

The source generators have images

\[
\begin{array}{c|cc}
&\text{linear generator}&\text{quadratic generator}\\ \hline
G_2/P_1&\beta(x)=6h&\beta(q)=12h^2\\
G_2/P_2&\beta(y)=2h&\beta(p)=12h^2=3\beta(y)^2.
\end{array}
\]

Consequently the entire image rings have integral presentations

\[
\boxed{
J^{\rm char}_{1}=\mathbb Z[6h,12h^2]
 \cong\mathbb Z[X,Q]/(X^2-3Q,Q^3),\qquad
J^{\rm char}_{2}=\mathbb Z[2h]\cong\mathbb Z[Y]/(Y^6).
}
\]

Here deg X=deg Y=1 and deg Q=2. The kernels on the full source rings
are (x^2-3q,q^3) and (p-3y^2,y^6), respectively.

These are complete integral kernels: the first quotient is free with basis
1,X,Q,XQ,Q^2,XQ^2; the second has basis 1,Y,...,Y^5. Their images are
nonzero in different degrees. Thus no further additive or polynomial
relation can occur. This argument also excludes hidden torsion in the
presented rings; a rational kernel calculation alone would not suffice.

Each cohomology group in these degrees has rank one. Its Smith factor is
therefore the positive generator of the embedded image:

| d | Character G2/P1 | Character G2/P2 | Earlier cleared/primitive form, either space |
| --- | --- | --- | --- |
| 0 | 1 | 1 | 1 |
| 1 | 6 | 2 | 1 |
| 2 | 12 | 12 | 1 |
| 3 | 144 | 48 | 2 |
| 4 | 288 | 288 | 2 |
| 5 | 1728 | 576 | 2 |

All groups above degree 5 vanish. Both character images are properly
contained in both earlier images in every positive nonzero degree. The
quotient orders are obtained by dividing the character entry by the earlier
entry. The earlier P1 form clears the catalogue quadratic's image
denominator; it is not the unmodified fractional journal input.
See the [earlier source conventions](EXCEPTIONAL_INTEGRAL.md).

## Adjoint representation: a second certified source

The same calculation for the 14-dimensional adjoint representation gives

\[
\Lambda_{14}=24\mathbb Z\omega_1+8\mathbb Z\omega_2=4\Lambda_7.
\]

Let C1,C2 be the pullbacks of 24omega1,8omega2. Replace them by
(C1-3C2,C2), a unimodular change. Outside characteristic 3 their
leading exponents are (3,-2),(3,-1), with coefficients 3,1. In
characteristic 3 they become (2,-1),(3,-1), with coefficients 1,1.
Both exponent determinants are nonzero. Coefficient rows at (2,-1)
and (1,0) give a unit minor before the change, proving the linear
lattice is complete.

Thus the adjoint torus form is Z[24omega1,8omega2]. For any fixed
standard Levi, its homogeneous source and Borel image in degree d are
4^d times those of the seven-dimensional representation. This is a
deduction from these two certified sources and homogeneous scaling; it
does not identify the two Springer maps. The public G2 image helper
currently selects the seven-dimensional representation.

## Independent verification

The [Sage experiment](experiments/g2_character_form.py) independently builds
the seven weights by reflections and compares them with WeylCharacterRing.
It checks the published coordinate formula, trace matrix, derivative,
leading certificates and both invariant generators.

It computes 24 direct Laurent intersections (homogeneous and filtered
through degree 5, for both Levis), using 4,940 coefficient positions.
Hermite preimages and gcd-one maximal-minor certificates agree with the
proposed polynomial bases. All 32 source monomial columns through degree 6,
including vanishing above the top degree, agree between BGG and the
separate root/coroot Chevalley implementation.

The [standard-library Python audit](experiments/audit_g2_character_form.py)
imports no Sage or runtime library code. It reconstructs both weight
systems, trace and leading certificates, the two unit minors, all direct
Laurent matrices and their completeness minors. It builds the 12-element
Weyl group and applies independently implemented polynomial divided
differences to all 32 columns. It verifies 24 embedded comparisons and
the image-ring relations.

The [dataset](experiments/g2-character-form-results.json) and
[audit report](experiments/g2-character-form-audit.json) retain source
hashes with CRLF normalized to LF. Earlier experiment reports remain
historical snapshots of their recorded source versions; the version-0.10
package export change does not rewrite their hashes.

```sh
PYTHONPATH=src sage -python experiments/g2_character_form.py
python experiments/audit_g2_character_form.py
PYTHONPATH=src sage -python -m unittest tests.test_characters -v
```
