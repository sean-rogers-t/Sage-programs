# The F4 minimal Springer map has generic degree 96

2026-10-22. This is a characteristic-zero boundary proof, with independent
weight and fiber checks. The degree is a deduction recorded here; novelty
and external review are not established. The finite-field calculations
corroborate the written proof and are not used to infer its conclusion.

## 1. Theorem and the maps involved

Let G=F4(C) and rho=V(omega_4), the 26-dimensional minimal representation,
with the numbering used in the Kumar-paper notes. Define theta_rho:G->g
by trace projection, with no additional rescaling. Then

\[
 \boxed{\deg\theta_\rho=\deg(\theta_\rho|_T)
       =\deg(T/W\longrightarrow\mathfrak t/W)=96.}
\]

The torus restriction and induced adjoint-quotient map are, in fact,
finite flat morphisms of rank 96. Consequently all their scheme-theoretic
fibers have length 96. General fibers are reduced. The full group map
is not finite; Section 6 gives an explicit positive-dimensional zero-fiber
witness. Generic degree and global finiteness must be distinguished.

This is the degree of this particular Springer/generalized Cayley map,
not the minimum Cayley degree of F4. No parabolic or cohomology coefficient
ring is involved. F4 has trivial center, so there is no central cover
factor to restore.

## 2. Weights, character lattice and normalization

In orthonormal epsilon coordinates the nonzero weights are the 24 short
roots, each with multiplicity one:

\[
 \{\pm e_i:1\le i\le4\}
 \ \cup\
 \{\tfrac12(\pm e_1\pm e_2\pm e_3\pm e_4)\}.
\]

Zero has multiplicity two and contributes no derivative. The actual
character lattice is

\[
 M=P(F4)=Q(F4)=\mathbb Z^4+\mathbb Z\delta,
 \qquad\delta=(1/2,1/2,1/2,1/2).
\]

Its Euclidean covolume is 1/2. An integral basis is (e_1,e_2,e_3,delta);
the exponent vector of w in this basis is

\[
 (w_1-w_4,w_2-w_4,w_3-w_4,2w_4).
\]

The weights generate M. Computing the trace moment matrix in epsilon
coordinates gives sum_w ww^T=6I_4: the eight coordinate weights contribute
2I_4 and the sixteen half-sum weights contribute 4I_4. Therefore

\[
 \theta_\rho(t)=\frac16\sum_w w\,e^w(t).                 \tag{1}
\]

The 1/6 is an invertible target scaling and does not change degree.
The highest-weight character calculation independently recovers the
24 nonzero weights, their multiplicities, and the two zero weights.

## 3. Every boundary face is nondegenerate

Let Delta be the convex hull of the 24 nonzero weights. It is the
24-cell. Its 24 facets have supporting hyperplanes

\[
 \pm x_i\pm x_j=1\qquad(i<j).
\]

These normals are primitive in Hom(M,Z), so every facet has height one.
For example x_1+x_2=1 has the six vertices

\[
 e_1,e_2,\quad(1/2,1/2,\pm1/2,\pm1/2).
\]

Its center is (1/2,1/2,0,0), and its opposite-pair directions are

\[
 (1/2,-1/2,0,0),\quad(0,0,1/2,1/2),\quad(0,0,1/2,-1/2).
\]

Thus each facet is a three-dimensional crosspolytope, or octahedron.
All lower-dimensional faces are simplices. There are no other weights
on a facet: the only remaining weight, zero, is strictly interior.

After a monomial translation and a finite torus cover, the facet
character is

\[
 f=z_1+z_1^{-1}+z_2+z_2^{-1}+z_3+z_3^{-1}.              \tag{2}
\]

At a critical point, z_i-z_i^-1=0, so every z_i is +/-1. The possible
critical values of f are -6,-2,2,6. None is zero. There is therefore no
solution of f=df=0. Multiplying by a monomial preserves the singular-zero
equations, and a finite torus cover is etale in characteristic zero,
so this exclusion holds on the actual face torus as well.

A simplex face has no singular zero: its affine exponent independence
makes its polynomial and tangent derivatives an invertible linear
system on nonzero monomial values. This checks every proper face,
including all lower-dimensional ones, rather than only the facets.

By the [general boundary criterion](GENERAL_CAYLEY_DEGREE.md), the degree
equals the normalized volume in M. There is no boundary correction.

## 4. Volume 96, with an explicit determinant calculation

Triangulate each octahedral facet along one opposite pair. The four
tetrahedra give four four-dimensional cones from the origin. In the
example above, take the opposite pair e_1,e_2 and choose one vertex from
each of the other two opposite pairs. Each resulting four-by-four
determinant has absolute value 1/2 in epsilon coordinates.

Dividing by the lattice covolume 1/2 gives normalized cone volume one.
Thus the cone over each facet has volume four, and

\[
 \operatorname{Vol}_M(\Delta)=24\cdot4=96.               \tag{3}
\]

Here normalized volume means 4! times ordinary lattice volume. The
Euclidean volume of Delta is two, consistently giving 4!*2/(1/2)=96.
Using the Euclidean integer lattice by mistake would lose a factor two.
This factor is already included in (3); it is not a central-isogeny factor.

The verifier triangulates all 24 facets in the integral M basis. It
also computes the polytope volume independently from the Weyl orbit
of omega_4 in the fundamental-weight basis. Both give 96. Its face
vector is (24,96,96,24), and the weight-lattice Smith index is one.

The group/torus/quotient degree equality is the standard equivariant
degree reduction, consistent with
[Kostant--Michor, Corollary 3.3](https://arxiv.org/pdf/math/0109066).
The trace construction and that equality are published inputs;
the F4 boundary check above supplies the particular degree calculation.

## 5. What generalizes

**Odd-crosspolytope criterion.** Let chi be a Laurent polynomial on a
rank-r torus with actual character lattice M. Suppose its Newton polytope
is full dimensional with origin interior, and every facet polynomial is
either supported on simplex vertices, or is, up to monomial and finite
cover, a nonzero common multiple of

\[
 \sum_{i=1}^d(z_i+z_i^{-1})\quad\text{with d odd}.
\]

Suppose there are no additional facet monomials. Then every proper face
has nonsingular zero set and the gradient degree equals normalized volume.
The proof is exactly the critical-value and simplex argument above.
This includes both the F4 minimal case (d=3) and our earlier E6 minimal
case (d=5, with simplex facets also present).

The coefficients matter. For a general crosspolytope polynomial
sum_i(a_i z_i+b_i z_i^-1), the critical values are
2 sum_i epsilon_i sqrt(a_i b_i), with epsilon_i=+/-1. An odd number of
terms alone does not prevent vanishing. For example

\[
 z_1+z_1^{-1}+z_2+z_2^{-1}+2(z_3+z_3^{-1})
\]

has a singular zero at (1,1,-1). Likewise an extra monomial at the facet
center shifts the critical values and requires a new check. Thus this
criterion is a sufficient condition with explicit coefficient hypotheses,
not a classification of all exceptional representations.

**Finite-flat strengthening of the no-boundary criterion.** Under the
no-boundary hypothesis of the general degree theorem, the logarithmic
gradient is finite flat onto affine r-space. To prove it, consider the
cleared equations D_i chi-y_i as sections on X_Delta times A^r. On every
boundary orbit the y_i terms vanish because zero is interior; absence of
face singular zeros therefore excludes the boundary for every target,
not only general targets. Their common zero scheme is the graph of the
gradient, closed in this projective compactification. Projection to A^r
is proper. It is also affine, hence finite
([Stacks Project, Tag 01WN](https://stacks.math.columbia.edu/tag/01WN)).
Its nonempty source has dimension r, so its closed image has dimension r
and is all of A^r. A torus is smooth and Cohen--Macaulay; the finite
equal-dimensional map to regular A^r is flat by
[miracle flatness, Tag 00R4](https://stacks.math.columbia.edu/tag/00R4).
Its rank is the normalized-volume generic degree.

For F4 the Weyl-equivariant map consequently induces a finite map
T/W -> t/W: C[T] is finite over C[t]^W, and its invariant subring
C[T]^W is a finite submodule. Both quotients are affine four-space;
on the source this is the fundamental-character theorem for a simply
connected semisimple group
([Steinberg, Regular elements of semi-simple algebraic groups](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf)).
Miracle flatness again gives rank 96. The same reasoning applies to any
previous no-boundary example with a smooth source quotient, including
the simply connected E6 minimal example. These are deductions from
standard algebraic geometry, with no novelty claim.

## 6. Why this does not make the full group map finite

Take a torus point with z_i=e^{e_i}(t)=-1 for all four i, and
b=e^delta(t)=1. These coordinates satisfy b^2=product z_i=1, so the
point belongs to the actual F4 torus. It is not the identity, and F4 has
trivial center. Formula (1), or the product formula below, gives theta(t)=0.
Equivariance sends its entire nontrivial conjugacy class to zero.

Exactly 32 roots centralize this point, so its centralizer has dimension
4+32=36 and its conjugacy class has dimension 52-36=16. The verifier
checks the root evaluations and zero image directly. Thus G->g has a
positive-dimensional special fiber even though T->t and T/W->t/W are
finite flat. None of the degree or flatness assertions claims that every
group fiber consists of 96 points.

## 7. Independent original-fiber calculations

The first model enumerates the 24 exponents in the integral basis
(e_1,e_2,e_3,delta) and differentiates their monomials. It clears Laurent
denominators and saturates by the product of torus coordinates.

The second model uses the independently factored character

\[
 \chi=2+\sum_{i=1}^4(z_i+z_i^{-1})+b^{-1}F(z),
 \qquad F(z)=\prod_{i=1}^4(1+z_i),\quad b^2=\prod_i z_i.
\]

Its normalized epsilon-coordinate equations at target y are

\[
 2b(z_i^2-1)+z_i(2z_i\partial_iF-F)-12y_i b z_i=0,
\]

together with b^2-product z_i=0 and the saturation equation
v b product z_i=1. The factor 12 comes from doubling the half-weight
derivative and the trace normalization 6I_4. The relation for b is part
of the actual torus presentation; it must not be counted as an extra
cover after the calculation.

At normalized epsilon target (1,2,3,5), each model is checked over
F_32003 and F_32009. All four saturated fibers have length 96. An
independent integer staircase count agrees with the Sage normal basis.
Multiplication by b+2z_1+3z_2 has squarefree characteristic polynomial
in every case, certifying reducedness. For each prime its entire
characteristic polynomial agrees across the two models, so the check
compares the same target and coordinate function, not only dimensions.

Both implementations use Sage/Singular arithmetic; independence is in
their defining equations, lattice description and the geometric proof.
Specialized finite-field lengths alone would not establish the
characteristic-zero generic degree or global flatness.

```bash
PYTHONPATH=src sage -python tools/verify_f4_minimal_degree.py
```

See the [verifier](tools/verify_f4_minimal_degree.py) and
[exact report](f4-minimal-degree-verification.json). No runtime API was
changed by this research calculation.

## 8. Literature scope and next goal

The construction, equivariant degree reduction, toric intersection theory,
affine-proper criterion and flatness inputs are established results.
The odd-crosspolytope test is an elementary deduction from those inputs
and the displayed character calculation. The number 96 is a repository
deduction with the proof above. Targeted searches for F4 generalized
Cayley/Springer degrees did not locate a direct published statement of
this value; the results were mostly unrelated uses of the word Cayley.
That limited search does not establish priority or settle publication status.

**Follow-up completed:** the [minimal-adjoint family proof](F4_MIXTURE_DEGREES.md)
now gives the degree for every a V(omega_4) direct sum b V(omega_1),
with a,b nonnegative integers and a+b positive. For b>0 the degree is
384 except at a=b (degree 240) and a=5b (degree 336). The two ratios
have primitive-height-two boundary corrections, checked against twenty
reduced fibers in two models. The next goal is a coefficient-discriminant
criterion for predicting such jumps from facet characters and normal
jets; it aims to make these calculations reusable across representation
families rather than merely extending a table of isolated degrees.
