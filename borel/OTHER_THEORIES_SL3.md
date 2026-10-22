# Complete literal K-image factors for every SL3 representation

Exploration, 22 October 2026. This is a complete low-rank classification
in terms of the source lattice and the established quadratic/cubic
representation indices. Its purpose is to test the general K-image
framework on every representation, rather than only a selected pair.
The ambient K-ring and index theory are known. The image formulas are
deductions here, with no literature-priority claim.

## 1. Statement

Let rho be any finite-dimensional SL3 representation with nonzero
differential, including direct sums, and use the normalized trace
Springer morphism and full literal integral source. Write P for the
weight lattice and Q for the root lattice. There is a unique positive
integer s such that Lambda_rho is either sP or sQ.

Normalize the trace index I by I_V=1 and the cubic index A by A_V=1
for the defining three-dimensional V; then A_(V*)=-1. Put

\[
                c=\frac s2(1-A/I),\qquad g=\gcd(s,c)>0.            \tag{1}
\]

The number c is an integer for these actual literal sources. The
complete ambient K^0-image Smith factors, including units, are

| Full source lattice | SL3/B | Either SL3/P1 or SL3/P2 |
| --- | --- | --- |
| sP | (1,g,g,s^3/g,s^3/g,s^3) | (1,g,s^3/g) |
| sQ | (1,g,3g,3s^3/g,3s^3/g,3s^3) | (1,3g,3s^3/g) |

For G/G the image is Z. Thus this covers every standard parabolic.
The theorem does not assert that every pair of integers (s,c) is realized
by an actual representation.
The factors in each row are already in divisibility order. Equal factors
do not assert equality of embedded images. In particular dualizing rho
preserves these lists, since it replaces c by s-c.

## 2. The source and the only higher parameter

The old source theorem gives a full W-stable lattice Lambda contained
in P. Let s be the gcd of its fundamental-weight coordinates. Stability
gives lambda-s_i lambda=lambda_i alpha_i in Lambda. Transitivity on
the A2 roots and Bezout therefore give sQ contained in Lambda contained
in sP. Since P/Q has prime order three, there are exactly the two forms
in the theorem. This does not assume rho irreducible.

Write the reduced character differential in the two fundamental
characters as a_1 e_1+a_2 e_2. The quadratic and cubic index functionals
give a_1+a_2=I and a_1-a_2=A. Hence c=s*a_2/I. This uses the
[normalized primitive data](OTHER_THEORIES_COTANGENT.md), rather than
claiming a new definition of generalized Dynkin indices.

Use the integral type A ambient K-coordinates
x=L_1-1, y=L_2-1. They have the ordinary coinvariant relations
x^2+xy+y^2=0, x^3=y^3=0, and integral basis
1,x,y,x^2,xy,x^2y. The defining representation sends
omega_1 to x and omega_2 to x+y. The dual defining representation sends
them to x-x^2 and x+y+xy. The character-differential formula therefore
gives the exact normalized rational map

\[
 \omega_1\longmapsto x-(c/s)x^2,\qquad
 \omega_2\longmapsto x+y+(c/s)xy.                                \tag{2}
\]

This map fixes every class of source degree at least two. For sP the
integrality of (2) immediately gives c in Z. For sQ, evaluate (2) on
s(2omega_1-omega_2) and s(-omega_1+2omega_2); their quadratic coefficient
vectors are c*(-2,-1) and c*(1,2), so integrality again forces c in Z.

## 3. Smith reduction on the full flag

For sP use the image basis
1,U,V,U^2,UV,U^2V, U=s omega_1, V=s omega_2. It is a complete integral
basis by the ordinary weight coinvariant presentation and the common
K-kernel theorem. The top-degree factor is s^3. In the remaining four
rows, (2) gives a block matrix

\[
       \begin{pmatrix}sI_2&0\\cB&s^2I_2\end{pmatrix},\qquad
                  B=\begin{pmatrix}-1&-1\\0&1\end{pmatrix}.
\]

Since B is unimodular, this is equivalent to two copies of
[[s,0],[c,s^2]], each with factors (g,s^3/g). This proves the sP row.

For sQ put a=alpha_1=x-y and b=alpha_2=x+2y. The ordinary integral
root image has basis 1,a,b,a^2,ab,a^2b: its quadratic relation gives
b^2=-a^2-ab, its cubic products give a^3=b^3=-2a^2b and ab^2=a^2b,
and all degree-four products vanish. Their target images are independent,
with a^2=-3xy, ab=3(x^2+xy), a^2b=3x^2y, proving completeness of the
listed basis without assuming an unsaturated invariant-ideal presentation.
Scale by s in each source degree. The degree-one/two inclusion block is

\[
 \begin{pmatrix}sC&0\\cE&3s^2D\end{pmatrix},\quad
 C=\begin{pmatrix}1&1\\-1&2\end{pmatrix},\quad
 E=\begin{pmatrix}-2&1\\-1&2\end{pmatrix},\quad
 D=\begin{pmatrix}0&1\\-1&1\end{pmatrix}.
\]

Integral elementary operations make it the direct sum
[[s,0],[-c,3s^2]] and 3*[[s,0],[c,s^2]]. Their factors are
(g,3s^3/g) and (3g,3s^3/g). The top class contributes 3s^3, giving
exactly the sQ list. No rational row operation is used here.

## 4. Projective flags and verification

For sP the complete source image on P1 has basis 1,z,z^2, z=s omega_1.
Its two positive K-columns in the basis v=L-1, v^2 are (s,-c) and
(0,s^2). For sQ use the complete image basis
1,E_1,E_2, E_1=3s omega_1,
E_2=s^2(2omega_1-omega_2)(omega_1+omega_2). The two positive columns
are (3s,-3c) and (0,3s^2). Their gcds and determinants give the two
projective lists. P2 replaces c by s-c, which has the same gcd with s.
This proves all assertions.

The [actual-character audit](experiments/other_theories_sl3_classification.py)
independently computes complete Laurent Demazure matrices on all three
proper flags for 35 nontrivial irreducibles with highest weights
0<=a,b<=5 and six reducible examples. It extracts I and A directly
from second and third weight moments, rather than the cotangent
recursion, and obtains every displayed factor. The written proof,
rather than this bounded set, establishes the all-representation theorem.
[Report](experiments/other-theories-sl3-classification.json).

The higher-rank problem is to replace this single cubic correction by
the full primitive character vector and classify the resulting integral
extensions. The uniform SL3 answer is a completed benchmark for that
problem, not evidence that higher rank depends on only one gcd.
