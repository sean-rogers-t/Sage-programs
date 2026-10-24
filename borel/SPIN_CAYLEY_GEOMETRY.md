# Spin polynomial Cayley maps: finiteness, ramification and boundary

Date: 22 October 2026. Status: uniform written deductions with internal
exact corroboration; not externally reviewed and no novelty claim.
The purpose is to distinguish lost inverse branches from ramification,
and to determine whether the scalar family explains more than a degree.

## 1. Scope and structural theorem

Use the hypotheses and trace normalization of the
[polynomial-character theorem](SPIN_POLYNOMIAL_CAYLEY_DEGREES.md):
G=Spin_N(C), N>=3, r=floor(N/2), S irreducible spin for odd N and
full spin S+ direct-sum S- for even N, K=2^r, and an actual character
F(chi_S) with q=deg F>=1. Put

\[
 d=N(N-1)/2,\quad c=\chi_S/K,\quad
 k(c)=F'(Kc)/F'(K),\quad h(c)=c k(c),\quad a=\pi(g).
\]

Write det(tI-Y)=t^(N-2r) Q_Y(t), with Q_Y(t) even of degree 2r,
and let a_r(Y)=Q_Y(0). Let D_0 be the reduced zero set of a_r.
For even N, D_0 is the singular skew-matrix hypersurface (its equation
a_r is the square of a Pfaffian in orthonormal coordinates). For odd N,
all skew matrices are singular; D_0 instead means that the forced zero
eigenvalue has algebraic multiplicity at least three. It is not the
same condition as having matrix rank at most N-3.

By the **finite target locus** we mean the union of target open subsets
V such that theta_F^(-1)(V) -> V is finite. Its complement is called
the non-finite target locus here. Since these maps are affine, finite
and proper over a target open are equivalent. This definition concerns
an entire neighborhood and all its preimages, not just one fiber.

**Theorem A (maximal finite locus).**

| Group and polynomial | Finite target locus | Degree there |
| --- | --- | --- |
| N>=5, any q>=1 | so_N minus D_0 | 2rq |
| N=4, q>=2 | so_4 minus D_0 | 4q |
| N=3, q>=2 | so_3 minus D_0 | 2q |
| N=3 or 4, q=1 | the whole Lie algebra | 2r |

The map is finite flat on each listed locus. Thus, for N>=5, the
non-finite target locus is independent of F within this family.
This theorem concerns the full group-to-Lie-algebra morphism. Equality
of generic degrees with the adjoint-quotient map does not prove equality
of their finite loci; the quotient finite locus is not classified here.

**Theorem B (full normalized Jacobian).** For N>=4, globally on G,

\[
 \boxed{\Psi_F(g)=c^{(N-1)(N-4)/2} k(c)^{d-1}
 \left[k(c)+(k(c)+c k'(c))\frac{\operatorname{tr}a-N}{4}\right].} \tag{1}
\]

Here the derivative is taken in left-invariant source coordinates and
the same Lie-algebra target basis, so Psi_F(1)=1. In rank one, N=3,
the regular polynomial formula is instead

\[
 \boxed{\Psi_F=k(c)^2\bigl[c k(c)+(c^2-1)k'(c)\bigr].}    \tag{2}
\]

These formulas describe the entire critical locus as a zero set,
including the source boundary c=0. They are not an irreducible
factorization claim; displayed factors can share components.

On c h(c)!=0, the critical equation is

\[
 L=h+c h' R=0,\qquad R=(\operatorname{tr}a-N)/4.             \tag{3}
\]

Every such critical point has differential corank exactly one. The
source and image have equal Lie-centralizer dimensions there, so
regularity is preserved there even when the differential is singular.

## 2. Differential proof and geometric meaning

The double-cover Cayley chart has coordinates (X,c), with
c^2 det(I-X)=1 and theta_F=2h(c)X. Projection to X is etale. For
a tangent matrix H,

\[
 dc(H)=\frac c2\operatorname{tr}((I-X)^{-1}H),\qquad
 dc(X)=cR,\qquad R=\tfrac12\operatorname{tr}((I-X)^{-1}X).
\]

Differentiating theta_F in these chart coordinates gives

\[
 D\theta_F=2h I+2h' X\otimes dc,\qquad
 \det D\theta_F=2^d h^{d-1}(h+c h'R).                       \tag{4}
\]

The second summand has rank at most one. If h!=0, all differential
failure is therefore in the radial direction: at L=0 the kernel is
the line spanned by X. Neither X nor h' can vanish at such a point.

For a left-invariant input H in so_N,

\[
 dX(H)=\tfrac12(I+X)H(I-X).
\]

In orthonormal coordinates I-X=(I+X)^T. The action A H A^T on
skew matrices is the exterior-square action, with determinant
det(A)^(N-1). Consequently

\[
 \det(dX)=2^{-d}\det(I-X)^{N-1}=2^{-d}c^{-2(N-1)}.
\]

Multiplying this by (4), substituting h=c k and
tr(a)=N+4R, proves (1) on the dense chart. Its right side is
regular for N>=4, so the identity holds globally. For N=3 use
tr(a)=4c^2-1 to cancel the apparent negative power of c; this gives
(2). No boundary value is obtained by substituting into a pole.

Finally, for c h!=0, a and X determine one another by rational
matrix functions and Y=2hX is a nonzero scalar multiple of X.
Their Lie centralizers agree. The finite Spin cover does not change
Lie-centralizer dimensions. This proves the regularity assertion.

For pure powers F(u)=u^m, equation (1) reduces on c!=0 to

\[
 \Psi_m=c^{md-2(N-1)}
            \left(1+\frac m4(\operatorname{tr}a-N)\right). \tag{5}
\]

For every r>=2 and m>=1 there are regular semisimple critical points
in the finite target locus. To see this, choose distinct
s_i=z_i+z_i^(-1), all different from +/-2, with
sum s_i=2r-4/m. This is possible on a positive-dimensional affine
hyperplane, avoiding finitely many proper subvarieties. Choose roots
z_i and either Spin lift. Then c!=0, the source and image are regular
semisimple, a_r(Y)!=0, but the last factor of (5) vanishes.
Thus criticality is not synonymous with failure of regularity.

## 3. The boundary rank bound

The essential boundary fact is

\[
 c(g)=0\quad\Longrightarrow\quad
                 \operatorname{rank}(\theta_F(g))\le2.   \tag{6}
\]

Here rank is in the defining N-dimensional orthogonal matrix model.
This is a consequence of Kostant--Michor's Clifford construction.
Write g_j for its exterior-degree-j component in the Clifford algebra.
Their exterior exponential on the Cayley chart gives

\[
 c g_4=\tfrac12 g_2\wedge g_2.
\]

Both sides are polynomial on Spin, so this identity extends across
c=0. There g_2 wedge g_2=0. A bivector with square zero is decomposable
or zero, hence has rank at most two under the Lie identification.
The base Springer map is g_2, and theta_F=k(c)theta_S. This proves
(6). If k(0)=0 the boundary image is zero.

## 4. Why the finite target locus is maximal

For N>=4, the earlier monic scalar presentation over a_r!=0 already
proves finite flatness of degree 2rq. Let

\[
 P_Y(c)=Q_Y(2h(c))-4^r h(c)^{2r}/c^2,
\]

where the quotient is polynomial. It is the minimal polynomial of c
over the embedded target function field, up to a nonzero constant.

We use a standard integral-extension argument. If theta_F were finite
over an affine target neighborhood V, its coordinate ring B would be
finite over A=O(V). The subring A[c] in B is A[T]/(P_Y(T)), using the
constant leading coefficient and minimality over Frac(A). Every element
of B is integral over A and hence over A[c]. Lying over then forces
every scalar root (Y,c_0) of P to lift to a genuine g over Y with
c(g)=c_0. A scalar root with no possible lift rules out finiteness on
any such neighborhood. This does not count specialized roots as a
proof of generic degree.

For N>=5 and Y in D_0 with rank Y>=4, c=0 is a scalar root but (6)
excludes every lift. Such Y are dense in every component of D_0:
the rank-at-most-two skew-matrix variety has dimension 2N-3, hence
codimension (N-2)(N-3)/2>1 in so_N. It cannot contain a hypersurface
component. The non-finite target locus is closed by its definition
as the complement of the union of finite target opens. It therefore
contains all of D_0. This proves the N>=5 row of Theorem A.

For N=4, q>=2, k is nonconstant. If k(0)=0, the same c=0 scalar
root has no lift over any nonzero Y in D_0 because theta_F vanishes
at every c=0 source point. If k(0)!=0, choose a nonzero root c_0 of k.
Then P_Y(c_0)=a_r(Y)=0, but theta_F is zero on c=c_0, again excluding
a lift for nonzero Y. These targets are dense in D_0. This proves
the N=4, q>=2 row.

For N=3 use Spin_3=SL_2 directly. Write g=cI+Z, tr Z=0, and
Z^2=(c^2-1)I. In these Lie coordinates theta_F=k(c)Z. For target W
the equation is

\[
 k(c)^2(c^2-1)+\det W=0.                                  \tag{7}
\]

For det W!=0 this is monic up to a constant, has degree 2q, and
forces k(c) invertible. It reconstructs g=cI+W/k(c) and proves finite
flatness. If q>=2, any root of k is a scalar root over det W=0,
but every source with that c maps to zero. Nonzero nilpotent targets
cannot lift this root. They are dense in the nilpotent cone D_0.

Finally q=1 gives k=1. For SL_2 the global presentation
c^2=1-det W makes the map finite flat of degree two. For
Spin_4=SL_2 times SL_2, full spin is the sum of the two defining
representations, so theta_S is the product of these two maps,
globally finite flat of degree four. Added copies and trivial
summands change neither map. This proves the remaining rows.

## 5. Branching on the finite locus, and a discriminant trap

On the monogenic finite open with c h!=0, direct differentiation of
P, using the source equation, gives

\[
 P'_Y(c)=2^{2r+1}h(c)^{2r-1}c^{-3}L.                       \tag{8}
\]

The coefficient multiplying L is a unit there. Thus the relative
differential module is generated by dc with relation L dc=0, and
the target branch support is exactly the zero set of disc_c(P_Y).
Its complement is the finite etale locus. This is the standard
monogenic discriminant criterion, not a new theorem about arbitrary
Springer maps. For N=3, presentation (7) applies on all of det W!=0,
including source c=0, and its derivative gives the same Jacobian
criterion. For q=1,N=3 or 4 use the global presentations above.

The restriction to the correct presentation matters. In Spin_4's
base case write p=-det W_1 and v=-det W_2. The two SL_2 equations
are x^2=1+p and y^2=1+v, and c=(x+y)/2. The scalar equation for c is

\[
 c^4-\tfrac12(p+v+2)c^2+(p-v)^2/16=0.
\]

Its discriminant is (p-v)^2(p+1)^2(v+1)^2. However the actual product
map branches only over p=-1 or v=-1. The factor p-v arises because c
fails to generate the entire coordinate algebra there; it does not
describe ramification of the original map at generic p=v. The map
remains finite there. Thus blindly extending the scalar discriminant
across its excluded chart creates a false branch component.

## 6. Assessment, checks and next step

This family supplies a structural result worth retaining: an exact
maximal finite target locus, a full Jacobian, a boundary rank bound,
and regular semisimple critical points. It separates three issues:
ramification of a finite cover, failure of finiteness, and loss of
regularity. None of these notions may be substituted for another.
The non-finite locus assertion does not classify every exceptional
fiber or every escaping curve individually.

The [verifier](tools/verify_spin_cayley_geometry.py) checks the full
Jacobian from original weight moments and root directions, compares
actual Clifford-matrix derivatives including boundary points,
checks (8) in exact fiber algebras, and verifies the Spin_4 discriminant
artifact symbolically. Rational invariant tuples also corroborate the
regular semisimple critical-point construction. The
[report](spin-cayley-geometry-verification.json)
records the finite scope. These computations corroborate the separate
uniform proofs; they are not an external proof review.

Published inputs: [Kostant--Michor](https://arxiv.org/pdf/math/0109066),
Sections 7.7 and 7.14--7.15, for the Clifford projection and exterior
exponential. Standard algebra inputs include
[integral and finite morphisms](https://stacks.math.columbia.edu/tag/01WG)
and the [finite locally free discriminant criterion](https://stacks.math.columbia.edu/tag/0BVH).
The differential determinant lemma, exterior-square determinant and
rank-two skew-matrix dimension are elementary linear algebra.
The above combination is a repository deduction; no exhaustive
priority search or resolution of a published open problem is claimed.

Recommendation: stop expanding this scalar family into more tables.
Use its finiteness and boundary results as benchmarks for an individual
half-spin representation. That changes the representation family and
tests whether the scalar model captures the next problem at all. A
separate structural question is which non-finite behavior disappears
after taking adjoint quotients; it must not be inferred from degree equality.
