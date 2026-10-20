# E6: minimal and adjoint character sources have incomparable images

For simply connected E6, changing from V(omega1), dimension 27, to the
adjoint V(omega2), dimension 78, does not give the uniform scaling found
in G2 and F4. The full torus sources are

\[
\boxed{\mathcal C_{27}=\operatorname{Sym}_{\mathbb Z}(6Q),
\qquad \mathcal C_{78}=\operatorname{Sym}_{\mathbb Z}(24P),}
\]

where Q is the root lattice and P the weight lattice, with P/Q of order 3.
For every standard Levi L these give
A_(rho,L)=C_rho^(W_L), and J_(rho,L)=beta(A_(rho,L)).
All sources are compared inside the common rational Cartan polynomial
ring, and all images in the common integral cohomology ring of G/P_L.
Polynomial degree d corresponds to cohomological degree 2d.

Already on E6/P1 the divisor images are 18 Z h and 24 Z h, respectively,
where h=beta(omega1). Thus the two full integral image rings are
incomparable. The complete finite calculation also finds incomparable
lattices in every positive degree through the top degree 16.

## Published inputs and scope

The trace projection is the normalized Springer/generalized Cayley map
of [Kostant--Michor, Section 2.1](https://arxiv.org/pdf/math/0109066).
The polynomial-character/Borel construction is from
[Kumar, Part I](https://kumar.math.unc.edu/papers/kumar76.pdf).
The rational E6/P1 Levi generators are the Part III inputs recorded in
the [catalogue](EXCEPTIONAL.md), from
[Kumar--Xie, Section 7](https://kumar.math.unc.edu/papers/kumar86.pdf).
The earlier [minimal-character proof](E6_CHARACTER_FORM.md) establishes
C27=Sym_Z(6Q).

The adjoint identification and comparison below are deductions recorded
here. They are not attributed to those papers, and no originality or
literature-priority claim is made. This is not a classification of E6
representations. The all-Levi source statement is mathematical; the new
public image factory currently implements only E6/P1.

We use the simply connected group for both representations. The adjoint
representation has a central kernel, but its differential is faithful.
Its Laurent pullbacks lie in Z[Q] inside Z[P]; testing their integrality
in either group ring gives the same answer. We do not silently replace
the group or the flag variety when changing representation.

## Adjoint torus proof

In fundamental-weight coordinates let C be the symmetric E6 Cartan
matrix, with edges (1,3), (3,4), (2,4), (4,5), (5,6).
The adjoint weights are the 72 roots, with multiplicity one, and zero
with multiplicity six. Writing

\[
b_i=\sum_{\alpha\in\Phi}\alpha_i e^\alpha,
\qquad M=\sum_{\alpha\in\Phi}\alpha\alpha^t,
\]

direct root enumeration gives M=24C. Consequently theta_78^*(omega_i)
are the entries of (24C)^(-1)b.

Put x_i=24omega_i and B_i=theta_78^*(x_i). The coefficient vector at
e^alpha in B is C^(-1)alpha, the integral simple-root coordinates of
alpha. Thus all B_i are integral. At the six simple-root exponents
their coefficient matrix is the identity. This proves that the complete
linear lattice is exactly 24P.

To prove the full polynomial intersection, use the unimodular matrix

\[
U=\begin{pmatrix}
1&-1&-1&1&0&0\\
0&0&0&1&-1&0\\
0&0&0&0&1&-1\\
0&1&0&-1&0&1\\
0&1&0&-1&0&0\\
0&0&1&-1&0&0
\end{pmatrix}.
\]

The six leading exponents of UB in descending lexicographic order are
the rows of

\[
E=\begin{pmatrix}
2&0&-1&0&0&0\\
1&1&0&0&-1&0\\
1&1&0&-1&1&-1\\
1&1&0&-1&0&1\\
1&1&-1&0&0&-1\\
1&0&1&-1&0&0
\end{pmatrix}.
\]

Here det(U)=1, every leading coefficient is 1, and det(E)=3. These leading
terms survive in every prime characteristic. Independence of their
exponents over Q means distinct monomials in UB have distinct leading
Laurent exponents. Therefore B1,...,B6 are algebraically independent
modulo every prime.
The determinant 3 does not invalidate the characteristic-three case:
the exponent vectors need to be independent over Q, not over the
coefficient field, since their integer exponent sums distinguish monomials.

By the [least-denominator argument](CHARACTER_SOURCES.md#exact-criterion),
an integral Laurent polynomial in Q[B1,...,B6] must belong to
Z[B1,...,B6]. This proves C78=Sym_Z(24P), including inhomogeneous
polynomials. Taking W_L invariants proves the stated Levi source for
every L. The independent audit reconstructs the roots, unit coefficient
minor and leading terms using integer arithmetic without Sage.

## All-degree sums and intersections

Fix a standard Levi and degree d, and abbreviate

\[
M_d=J_{27,L}^d,\qquad
N_d=6^d\beta\bigl((\operatorname{Sym}_{\mathbb Z}^dP)^{W_L}\bigr).
\]

Since Q is contained in P and 3P is contained in Q,

\[
3^dN_d\subseteq M_d\subseteq N_d,\qquad
J_{78,L}^d=4^dN_d.
\]

The quotient N_d/M_d is killed by 3^d. Multiplication by 4^d is an
automorphism of that quotient. Hence M_d+4^dN_d=N_d, and
4^dn in M_d implies n in M_d. This proves, for every Levi and every d,

\[
\boxed{
J_{27,L}^d+J_{78,L}^d=N_d=4^{-d}J_{78,L}^d,\qquad
J_{27,L}^d\cap J_{78,L}^d=4^dJ_{27,L}^d.}
\]

The divided lattice on the right is integral by its displayed definition;
this is not informal denominator clearing. The same proof applies to
the source lattices before beta.

The additive sum over all degrees is the Borel image of
Sym_Z(6P)^(W_L), so it is already a graded subring. Allowing mixed
products of minimal and adjoint image classes adds nothing beyond this
additive sum.

If r_d is their common image rank, it follows that

\[
\frac{J_{27,L}^d+J_{78,L}^d}{J_{78,L}^d}
\cong(\mathbb Z/4^d\mathbb Z)^{r_d}.
\]

The quotient of the sum by the minimal image is a finite 3-group of
exponent at most 3^d. Its multiplicities are computed for E6/P1 below;
no general formula for those multiplicities for arbitrary Levis is asserted.

Equivalently, at prime 2 the adjoint image is 4^d times the minimal
image; at prime 3 the minimal image is contained in the adjoint image;
and at every other prime the two images agree. Thus each ordered
2-adic Smith valuation increases by 2d. There is no uniform global
scaling rule, because the root-to-weight lattice change affects prime 3.

## Divisors and the complete E6/P1 comparison

The W_L-fixed linear space for P1 is Q omega1.
The class of omega1 has order 3 in P/Q, so the integral linear sources
on this line are 18 Z omega1 and 24 Z omega1. Therefore

\[
J_{27}^1=18\mathbb Zh,\quad J_{78}^1=24\mathbb Zh,\quad
J_{27}^1+J_{78}^1=6\mathbb Zh,\quad
J_{27}^1\cap J_{78}^1=72\mathbb Zh.
\]

The degree-two Smith factors are 108 and 576. Their ratio is 16/3,
whereas squaring the divisor ratio 24/18 gives 16/9. Thus even the
embedded images on this one flag variety cannot be related by a single
degreewise scalar c^d.

The following records the nonunit invariant factors of
(J27^d+J78^d)/J27^d. The quotient by J78^d always has r_d factors 4^d
in positive degree. This small table captures the additional 3-primary
information alongside the general theorem.

| d | r_d | Sum modulo minimal |
|---|---|---|
| 0 | 1 | zero |
| 1 | 1 | 3 |
| 2 | 1 | 3 |
| 3 | 1 | 9 |
| 4 | 2 | 9 |
| 5 | 2 | 27 |
| 6 | 2 | 9 |
| 7 | 2 | 3, 9 |
| 8 | 3 | 3, 27 |
| 9 | 2 | 3, 9 |
| 10 | 2 | 3, 3 |
| 11 | 2 | 3, 9 |
| 12 | 2 | 3, 9 |
| 13 | 1 | 9 |
| 14 | 1 | 3 |
| 15 | 1 | 9 |
| 16 | 1 | 9 |

All positive rows have nonzero quotients in both directions, proving
degreewise incomparability in the complete E6/P1 table.
The [dataset](experiments/minuscule-character-api-results.json) includes
all adjoint Smith lists and embedded image, sum and intersection bases.
The sources retain their kernels; for example the degree-two source
has rank two while the image has rank one.

## API and verification

Use `e6_character_source(1, representation="adjoint")`; see the
[API guide](MINUSCULE_CHARACTER_API.md) for examples and resource limits.
The [torus experiment](experiments/e6_adjoint_certificate.py) checks the
root orbit against the Weyl character calculation, the identity
differential, the unit minor, the all-primes certificate and eight direct
homogeneous/filtered Laurent intersections through degree three.

The [full API experiment](experiments/minuscule_character_api.py) computes
all 17 adjoint source lattices with completeness certificates, and checks
every image column against the separate published-ring/Chevalley engine.
The [Python-only audit](experiments/audit_e6_adjoint_api.py) reconstructs
all adjoint source coefficients and checks every image and comparison by
exact arithmetic. Borel-method independence and the ordinary-ring
engine's shared BGG anchor dependencies are stated in the API guide.

The finite tables establish all degrees of E6/P1 cohomology because the
target vanishes above degree 16. The all-degree, all-Levi source and
sum/intersection assertions follow from the proofs above, independently
of the finite range of the tables.
