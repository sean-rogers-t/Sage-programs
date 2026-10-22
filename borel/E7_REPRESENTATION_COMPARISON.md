# E7: minimal and adjoint integral Springer/Borel images

Update, 20 October 2026: the [2-primary basis theorem](E7_ROOT_WEIGHT_COMPARISON.md)
now explains the P7 relative factors by explicit simultaneous bases.
Only factors 2 occur, with polynomial q[14]_(q2)+(q6+q10)[7]_(q2).
A stable common lattice complement and a root-image quadratic also
prove exponent 2 on all 64 retained subsets of the E6 Levi. Their
multiplicities are not supplied by that bound. See the separate
[proof/source audit](E7_ROOT_WEIGHT_AUDIT.md).

The subsequent [complete E7 theorem](E7_ALL_PARABOLICS.md) now determines
the relative full root/weight factors on all 127 proper standard
parabolics, supplying the corresponding scaled image comparisons.
The [A6 intermediary theorem](E7_FULL_FLAG_COMPARISON.md) was the first
109-case step. The completion retains the source normalization below;
every nonunit relative root/weight image factor is two.

For simply connected E7, with root lattice Q and weight lattice P,
the full integral polynomial-character sources are

\[
\boxed{\mathcal C_{56}=\operatorname{Sym}_{\mathbb Z}(12Q),
\qquad \mathcal C_{133}=\operatorname{Sym}_{\mathbb Z}(36P).}
\]

The representations are V(omega7), dimension 56, and the adjoint
V(omega1), dimension 133. Here the character lattice of the maximal
torus T is P, Z[P] is its integral Laurent character ring, and

\[
\mathcal C_\rho=
\{f\in\operatorname{Sym}_{\mathbb Q}(P\otimes\mathbb Q):
(\theta_\rho|_T)^*f\in\mathbb Z[P]\}.
\]

The Springer morphism theta_rho uses trace projection normalized to
have identity differential. For every standard Levi L, let W_L be
its Weyl group, beta the rational Borel characteristic map, and
A_(rho,L)=C_rho^(W_L) and J_(rho,L)=beta(A_(rho,L)).
Then in every polynomial degree d,

\[
\boxed{J_{56,L}^{d}+J_{133,L}^{d}=3^{-d}J_{133,L}^{d},
\qquad J_{56,L}^{d}\cap J_{133,L}^{d}=3^dJ_{56,L}^{d}.}
\]

These formulas compare embedded lattices in common cohomology coordinates.
They also hold for the source lattices in the common rational polynomial
ring. Polynomial degree d corresponds to H^(2d). The cohomology images
are integral because both source lattices consist of integral weight
polynomials and integral Schubert pullback from G/P_L to G/B is primitive.

Already on E7/P7 the divisor images 24 Z h and 36 Z h are incomparable,
so the two complete image rings are incomparable. Some individual higher
degrees do have containment, as the complete table below shows.

## Published inputs and status

We use the normalized trace-projection construction of
[Kostant--Michor, Section 2.1](https://arxiv.org/pdf/math/0109066),
and the Springer-polynomial/Borel construction of
[Kumar, Part I](https://kumar.math.unc.edu/papers/kumar76.pdf).
The rational E7/P7 invariant generators omega7, psi2, psi5, psi6, psi8,
psi9 and psi12 come from
[Kumar--Xie, Part III, Theorem 8.2(g)](https://kumar.math.unc.edu/papers/kumar86.pdf).
Their catalogue normalization is unchanged.

The [earlier minimal-source proof](E7_CHARACTER_FORM.md) establishes
C56=Sym_Z(12Q). The adjoint source certificate and integral comparison
below are checked deductions in this repository, not claims attributed
to those papers. No originality or comprehensive literature-priority
claim is made. A general representation classification is not asserted.

We keep the simply connected group throughout. The adjoint representation
has a finite central kernel and faithful differential. Its Laurent
pullbacks have support in Q inside P; integrality in Z[Q] and in Z[P]
agrees for those functions. The invariant polynomial sources are
transported through their respective Springer maps, whose Laurent
functions are different.

## The full adjoint torus source

Let C be the symmetric Cartan matrix with edges (1,3), (3,4), (2,4),
(4,5), (5,6), (6,7). In fundamental-weight coordinates the adjoint
weights are the 126 roots with multiplicity one and zero with multiplicity
seven. Set

\[
b_i=\sum_{\alpha\in\Phi}\alpha_i e^\alpha,\qquad
M=\sum_{\alpha\in\Phi}\alpha\alpha^t=36C.
\]

Trace projection gives theta_133^*(omega)=(36C)^(-1)b.
The pullbacks B_i=theta_133^*(36omega_i) have coefficient vector
C^(-1)alpha at e^alpha, hence integer coefficients. At the seven
simple-root exponents this coefficient matrix is the identity. Thus the
exact linear lattice is 36P.

The following determinant-one row change gives an all-primes certificate:

\[
U=\begin{pmatrix}
1&-1&-1&1&0&0&0\\
0&0&0&1&-1&0&0\\
0&0&0&0&1&-1&0\\
0&0&0&0&0&1&-1\\
0&1&0&-1&0&0&1\\
0&1&0&-1&0&0&0\\
0&0&1&-1&0&0&0
\end{pmatrix}.
\]

In descending lexicographic order the leading exponents of UB are the
rows of

\[
E=\begin{pmatrix}
2&0&-1&0&0&0&0\\
1&1&0&0&-1&0&0\\
1&1&0&-1&1&-1&0\\
1&1&0&-1&0&1&-1\\
1&1&0&-1&0&0&1\\
1&1&-1&0&0&0&-1\\
1&0&1&-1&0&0&0
\end{pmatrix}.
\]

All leading coefficients are 1, and det(E)=2. The exponent vectors
are independent over Q, so distinct monomials in UB have distinct
leading Laurent exponents in every characteristic. Independence over
the coefficient field is not required; in particular determinant 2
causes no failure in characteristic two.

The [least-denominator argument](CHARACTER_SOURCES.md#exact-criterion)
therefore proves C133=Sym_Z(36P), including inhomogeneous polynomials.
Taking W_L invariants proves the source statement for every standard
Levi. This all-degree argument is separate from the finite cohomology
table computation.

## Integral weight coordinates for efficient source calculation

Use the integral Q basis v1,...,v6,u from the
[minimal E7 proof](E7_CHARACTER_FORM.md#integral-symmetric-coordinates-make-the-computation-feasible):

\[
v_i=\alpha_7+\sum_{j\in(1,3,4,5,6)\text{ from position }i}\alpha_j
\quad(1\leq i\leq6),\qquad
u=\alpha_2+v_1+v_2+v_3.
\]

The sum is empty for i=6, so v6=alpha7. The A5 subgroup generated by
s1,s3,s4,s5,s6 permutes the six v's and fixes u. Put E1=sum(v_i) and

\[
z=(u-E_1)/2=\omega_7-u.
\]

Since P/Q has order 2 and omega7 represents its nonzero class,
v1,...,v6,z is an integral basis of P; its determinant in fundamental
weights is -1. A5 fixes z. Consequently

\[
\operatorname{Sym}_{\mathbb Z}(P)^{W_{A_5}}
=\mathbb Z[E_1,E_2,E_3,E_4,E_5,E_6,z],
\]

where Ej=e_j(v1,...,v6). Every full E6-Levi invariant lies in this
A5-invariant ring. Its integrality in Sym_Z(P) is equivalent to integer
coefficients in these symmetric coordinates.

The old root-coordinate formulas become

\[
u=2z+E_1,\qquad \omega_7=E_1+3z,\qquad
x_i=6v_i-E_1,\qquad Z=-3E_1-12z.
\]

The catalogue psi_d is the sum of the d-th powers of the 27 forms
x_i+Z, x_i-Z and -x_i-x_j (i<j). Substituting u=2z+E1 into the
previous Newton formulas gives its compressed polynomial. The Sage
check verifies the coordinate determinant, all five permutation actions,
the catalogue's defining linear forms and every polynomial identity.

Independently, the Python audit reconstructs the power sums by the closed
partition formula

\[
p_k(v)=\sum_{\sum j m_j=k}
(-1)^{k+\sum m_j}
\frac{k(\sum m_j-1)!}{\prod m_j!}\prod E_j^{m_j},
\]

and expands the 27 linear forms by binomial identities. It does not
import the runtime Newton implementation.

## Why the sums and intersections follow

Fix L and d and write

\[
M_d=J_{56,L}^d,\qquad
N_d=12^d\beta((\operatorname{Sym}_{\mathbb Z}^dP)^{W_L}).
\]

The inclusions Q subset P and 2P subset Q give
2^d N_d subset M_d subset N_d, while J133,L^d=3^d N_d.
Multiplication by 3^d is an automorphism of the finite 2-group N_d/M_d.
It follows that M_d+3^dN_d=N_d and M_d intersect 3^dN_d=3^dM_d.
This proves the displayed formulas directly for image lattices.

If r_d is the common image rank, then

\[
(J_{56,L}^d+J_{133,L}^d)/J_{133,L}^d
\cong(\mathbb Z/3^d\mathbb Z)^{r_d}.
\]

The quotient of the sum by J56,L^d is a 2-group killed by 2^d.
At prime 3 the adjoint image is 3^d times the minimal image; at prime 2
the minimal image is contained in the adjoint image; at all other primes
the images agree. In particular each ordered 3-adic Smith valuation
increases by d.

The additive sum over all degrees is the image of Sym_Z(12P)^(W_L).
It is already a graded subring, so mixed products add nothing beyond
the additive sum.

More generally, the same elementary argument applies to certified sources
Sym_Z(a Lambda) and Sym_Z(ac Lambda'), when Lambda subset Lambda',
the quotient has exponent e, and the positive integers c and e are
coprime. Their degree-d image sum is the image from a Lambda', and
their intersection is c^d times the first image. This requires full-source
identifications, not merely proportional trace matrices.
The [E6 comparison](E6_REPRESENTATION_COMPARISON.md) has (c,e)=(4,3);
E7 has (c,e)=(3,2).

## Divisors and complete E7/P7 images

Since omega7 has order 2 in P/Q and h=beta(omega7) is primitive,

\[
J_{56}^1=24\mathbb Zh,\qquad J_{133}^1=36\mathbb Zh,\qquad
J_{56}^1+J_{133}^1=12\mathbb Zh,\qquad
J_{56}^1\cap J_{133}^1=72\mathbb Zh.
\]

In degree two the Smith factors are 144 and 1296. Their ratio is 9,
whereas squaring the degree-one ratio gives 9/4. Thus the image
comparison is not a uniform degreewise scalar rescaling.

The complete table is linked in the
[comparison dataset](experiments/e7-adjoint-comparison-results.json).
It retains every source basis, image lattice, Smith list, sum, intersection
and quotient, in all 28 degrees. The target vanishes above degree 27,
so these tables determine the entire E7/P7 cohomology images.

In this finite calculation the quotient of the sum by the minimal image
is elementary abelian:

\[
(J_{56}^d+J_{133}^d)/J_{56}^d\cong(\mathbb Z/2\mathbb Z)^{c_d}.
\]

| Polynomial degrees d | c_d |
| --- | --- |
| 0, 2, 4, 24, 26 | 0 |
| 10, 12, 14, 16, 18 | 2 |
| All other degrees from 1 through 27 | 1 |

Thus J133 is properly contained in J56 in degrees 2, 4, 24 and 26.
The lattices are incomparable in every other positive degree through
27, and agree in degree zero. The elementary-abelian assertion here
is the fully checked E7/P7 table; the general all-Levi argument above
only gives exponent dividing 2^d.

The two P7 source lattices themselves are incomparable in every positive
degree. Indeed, (36omega7)^d belongs to the adjoint source, but after
dividing by 12^d its v1^d coefficient in integral Q coordinates is
(-3/2)^d, so it does not belong to the minimal source. Conversely,
(24omega7)^d belongs to the minimal source since 2omega7 lies in Q,
but division by 36^d gives the nonintegral weight coefficient (2/3)^d.
Thus the four image containments above arise after passage through the
Borel kernel; they are not containments of the character sources.

## API and independent verification

```python
from sage_borel import e7_character_source

minimal = e7_character_source(7)
adjoint = e7_character_source(7, representation="adjoint")
a = minimal.degree(1).image()
b = adjoint.degree(1).image()
c = a.compare(b)
assert a.smith_factors == (24,)
assert b.smith_factors == (36,)
assert c.sum.smith_factors == (12,)
assert c.intersection.smith_factors == (72,)
assert c.sum.quotient_by(a).torsion_factors == (2,)
assert c.sum.quotient_by(b).torsion_factors == (3,)
```

Only P7 is packaged; the minimal representation remains the default.
High degrees may require an explicit larger coefficient budget, such as
`adjoint.degree(27, max_matrix_entries=10000000)`.
See the [API guide](MINUSCULE_CHARACTER_API.md).

The [Sage experiment](experiments/e7_adjoint_comparison.py) compares the
adjoint root orbit with the Weyl character calculation, checks the
identity differential and all-primes certificate, and verifies the
integral P model. Eight direct Laurent-intersection checks cover
homogeneous and filtered degrees 0..3. All 28 source coefficient lattices
have integrality and gcd-one maximal-minor completeness certificates.
The experiment captures the coefficient matrix actually passed to the
runtime Hermite solver, then certifies its output.

Every one of the 1,967 rational Borel columns is freshly compared with
the separate published ordinary-ring/Chevalley calculation. That engine
retains its documented BGG generator-image dependencies. The minimal
tables are reused from their previous full audit, with six fresh
minimal API degree checks.

The [Python-only audit](experiments/audit_e7_adjoint_comparison.py)
independently reconstructs roots, trace projection, the torus certificate,
the integral P basis, compressed generators, every source coefficient
constraint and image/comparison arithmetic. It checks the sum/intersection
identities in every degree. Its [saved result](experiments/e7-adjoint-comparison-audit.json)
binds the exact dataset and historical tables used.
