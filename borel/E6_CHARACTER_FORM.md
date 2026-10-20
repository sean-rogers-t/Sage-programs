# E6/P1: the integral form defined by Levi characters

For simply connected E6, the normalized Springer map attached to the
27-dimensional representation of highest weight omega1 gives the following
exact integral form. Let Q be the E6 root lattice and let W_L be generated
by s2,...,s6, the Weyl group of the Levi of P1. Inside the rational symmetric
algebra on the Cartan dual, the character integral form is

\[
\boxed{\mathcal A_{\mathbb Z}=
       \operatorname{Sym}_{\mathbb Z}(6Q)^{W_L}.}
\]

This determines the form in all degrees. It is an invariant-ring
description, not a claim that a particular list of six invariant
polynomials freely generates it over Z. Its image under the Borel map has
been calculated in every degree of E6/P1 and compared with both previously
chosen forms. See the [complete comparison table](E6_CHARACTER_TABLES.md)
and [exact dataset](experiments/e6-character-form-results.json).

Already in polynomial degree 1 the character image is 18 Z h, whereas
the original catalogue image is 3 Z h and the individually primitive image
is Z h. In degree 4 the three Smith lists are respectively (3888,23328),
(3,288), and (1,16). The literal character image is contained in both
previous images in every degree. This is established by the recorded
integer containment matrices, not by comparing Smith lists alone.

## Definition and published inputs

Let P denote the weight lattice of simply connected E6, so the torus
character ring is Z[P]. Write L=L1 for the Levi of P1, theta=theta_omega1,
and beta for the Borel characteristic map. Define

\[
\mathcal A_{\mathbb Z}
=\{f\in\operatorname{Sym}_{\mathbb Q}(P\otimes\mathbb Q)^{W_L}:
          \theta^*f\in\mathbb Z[P]^{W_L}\}.
\]

Then theta* identifies this ring with the rational Springer-polynomial
subalgebra intersected with Rep(L). Here Rep(L) means integral **virtual**
characters; positivity is not required. The identification
Rep(L)=Z[P]^(W_L) follows from the unit triangular highest-weight expansions
of irreducible characters in Weyl-orbit sums. Rational coefficients cause
no loss relative to the complex intersection: the pullback has rational
coefficients, is injective, and every polynomial has finite degree, so a
rational linear system with integral output has rational unique input.
The torus lattice is the full E6 weight lattice. Replacing this Levi by
an abstract SO(10) times a one-dimensional torus would change the character
lattice and is not part of this calculation.

[Kumar, Part I, Definition 4 and Theorem 5](https://kumar.math.unc.edu/papers/kumar76.pdf)
provide the polynomial-character construction and its transported Borel
map. The rational Levi invariants used below are those of
[Kumar--Xie, Part III, Theorem 7.2(a)](https://kumar.math.unc.edu/papers/kumar86.pdf).
The six Springer coordinates retain the corrected denominator 18 from the
earlier [normalization audit](EXCEPTIONAL.md), not the printed denominator 6.

The all-degree integral-intersection argument below is a deduction supplied
here. The character definition is published; this note does not claim that
the exact E6 integral form or its Smith tables are new in the literature.

## Proof of the full torus intersection

Use fundamental-weight torus coordinates t_i=e^(omega_i) and simple
coroots H_i. Let Omega be the 27 weights of V(omega1), each with multiplicity
one. For mu in Omega put mu_i=mu(H_i), and set

\[
b_i=\sum_{\mu\in\Omega}\mu_i e^\mu.
\]

The trace Gram matrix is

\[
M_{ij}=\sum_{\mu\in\Omega}\mu_i\mu_j=6C_{ij},
\]

where C is the symmetric E6 Cartan matrix with edges (1,3), (3,4),
(2,4), (4,5), (5,6). The trace-projection construction gives M Theta=b.
Since alpha_i=sum_j C_ij omega_j, it follows that

\[
\theta^*(6\alpha_i)=b_i\in\mathbb Z[P].
\]

Thus Z[6alpha1,...,6alpha6] maps into Z[P]. To prove that no further
rational polynomials become integral characters, it suffices to prove
that the six b_i remain algebraically independent over every prime field.

Make the following integer change of coordinates, c=Ub:

\[
U=\begin{pmatrix}
0&-2&-1&-3&-2&-1\\
1&-1&1&0&0&0\\
0&-2&-1&-2&-2&-1\\
0&-2&-1&-2&-1&-1\\
0&-2&-1&-2&-1&0\\
0&-3&-1&-3&-2&-1
\end{pmatrix},\qquad \det U=-1.
\]

Order Laurent exponents lexicographically, comparing the t1 exponent
first. The leading exponents of c1,...,c6 in characteristic other than 3
are the rows of

\[
E=\begin{pmatrix}
1&1&0&-1&0&0\\
1&0&0&0&0&0\\
1&0&-1&1&-1&0\\
1&0&-1&0&1&-1\\
1&0&-1&0&0&1\\
1&-1&0&0&0&0
\end{pmatrix}.
\]

Their coefficients are (1,1,1,1,1,3), and det E=-3. In characteristic 3,
the leading exponent of c6 becomes (0,1,0,0,0,-1), with coefficient -2;
the other five stay the same. The resulting exponent matrix has
determinant 2.

In either case the six exponent vectors are linearly independent over Q.
Hence distinct monomials in the c_i have distinct leading Laurent
exponents, and no nonzero polynomial in them can vanish. This uses
independence of integer exponent vectors, not invertibility of their
determinant in the coefficient field. Lexicographic order is compatible
with addition on Z^6; finite Laurent supports suffice, so a well-order on
all Laurent monomials is unnecessary. This proves algebraic independence
over every prime field.

Now write a rational polynomial in the b_i as F(b)/D, with F in Z[X1,...,X6]
and D a least positive common denominator. If it is an integral Laurent
polynomial and a prime p divides D, then F(b) is zero modulo p. The least
denominator implies that F is nonzero modulo p, contradicting the
independence just proved. Therefore D=1. Consequently

\[
\mathbb Q[b_1,\ldots,b_6]\cap\mathbb Z[P]
=\mathbb Z[b_1,\ldots,b_6],
\qquad
(\theta^*)^{-1}\mathbb Z[P]=\operatorname{Sym}_{\mathbb Z}(6Q).
\]

The root lattice is W_L-stable and theta* is injective and equivariant.
Taking W_L-invariants therefore proves the boxed description. In particular
the integral form is graded. Mixed-degree cancellation cannot enlarge it;
this conclusion follows from the all-degree proof, rather than being
assumed in the degreewise calculations.

The full-torus argument also gives Sym_Z(6Q)^(W_L) for every standard Levi
in this same simply connected E6 with this same representation. Only the
P1 cohomology image is calculated here. Integrality of the Borel images
follows from integral weight polynomials and the primitive Schubert-basis
inclusion H*(G/P1,Z) into H*(G/B,Z).

## Exact finite-degree calculation

Write y=y1 and let e_j denote the elementary symmetric polynomial in
y2^2,...,y6^2; put p=y2*y3*y4*y5*y6. The rational invariant ring is freely
generated by y,e1,e2,p,e3,e4, of degrees 1,2,4,5,6,8. The redundant catalogue
generator e5 equals p^2. Thus their weighted monomials give a complete
rational basis in each degree.

For degree d, substitute omega=C^(-1) alpha in that monomial basis and
form its coefficient matrix M_d in ordinary simple-root monomials. The
full invariant lattice for Sym_Z(Q) is exactly

\[
\{v\in\mathbb Q^{n_d}:M_dv\in\mathbb Z^{N_d}\}.
\]

The character-form lattice is 6^d times this lattice. No invariant is
discarded just because its Borel image is zero: such terms can change
which combinations are integral characters.

The implementation computes a lattice basis B by Hermite reduction and
then checks completeness independently. For N=M_d B it verifies integral
entries and selects several nonzero maximal minors using modular pivots.
Their exact integer determinants have gcd 1. Thus the gcd of all maximal minors is 1,
so N has primitive column lattice. This proves that B contains every
allowed rational coefficient vector, rather than just producing some
integral characters.

The first homogeneous character-form bases, expressed in catalogue
polynomials, are:

| d | A Z-basis of the source component |
| --- | --- |
| 0 | 1 |
| 1 | 6y |
| 2 | 36e1; (9/2)(y^2-5e1) |
| 3 | 216ye1; 27y^3-135ye1 |
| 4 | 648e2; 162e1^2-405e2; 162y^2e1-81e2; (81/4)y^4-(81/2)y^2e1-(567/4)e1^2+243e2 |

These fractional catalogue expressions are integral virtual characters
after Springer pullback. Merely multiplying each catalogue generator by
its individual Laurent denominator would miss some of them.

## Independence of the checks

The [Sage calculation](experiments/e6_character_form.py) reconstructs the
27 weights by reflection and compares them with WeylCharacterRing. It
checks the trace normalization, identity derivative, and the leading-term
certificate. Through degree 4 it also computes the literal Laurent
coefficient intersection directly, both homogeneously and allowing all
lower degrees. Those lattices agree with the root-polynomial method by
unimodular basis changes.

A [separate standard-library Python audit](experiments/audit_e6_character_form.py)
uses no Sage or sage-borel imports. It reconstructs the weights, trace
matrix, Springer coordinates and leading terms with exact fractions. It
also verifies all 36 reflection identities for the moment coordinates and
independently reconstructs the ten homogeneous/filtered Laurent lattices
through degree 4. Both integrality and completeness are checked, using
integer determinants and modular elimination. Its
[report](experiments/e6-character-form-audit.json) records the exact scope.
It checks 256047 Laurent coefficient positions, and separately verifies
all 51 image lattices/Smith lists and all 34 embedded containments in the
three-way comparison, using rational inverses and gcds of integer minors.

For every degree through 16, the Borel monomial images are separately
calculated by fixed-point localization and by the published integral
cohomology ring with ordinary Chevalley multiplication, as in the
[previous E6/E7 verification](EXCEPTIONAL_MINUSCULE.md). The latter retains
its documented BGG starting-image dependencies. Smith factors are checked
by determinantal divisors of an integer lattice basis. Containment in
both previous images is checked in common Schubert coordinates.
The calculation uses 382 source-basis monomials across the 17 degrees.

Because the source form is proved graded and H^(2d)(E6/P1)=0 for d>16,
these finite tables determine its entire cohomology image. They do not
assert that the listed low-degree source elements generate the full
invariant ring over Z. The root-lattice description already determines
that full ring without requiring a minimal integral presentation.

## Reproduction

From `borel/`:

```sh
PYTHONPATH=src sage -python experiments/e6_character_form.py
python experiments/audit_e6_character_form.py
python tools/render_e6_character_tables.py
```

The result files include source hashes with CRLF normalized to LF. The
calculation is currently a research experiment; the public exceptional
catalogue and its existing normalization are unchanged.
