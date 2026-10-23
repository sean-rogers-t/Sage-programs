# Reconstructing Borel Springer formulas from torus data

## The theorem and its scope

Let G be a connected reductive complex group, B=T U a Borel subgroup,
and theta:G -> mathfrak g a conjugation-equivariant regular map.
In the Springer application, theta=theta_rho is trace projection for a
representation rho with nondegenerate trace form on mathfrak g. Let
h(t)=theta(t), which belongs to mathfrak t. Choose positive root vectors
E_alpha and x_alpha(s)=exp(s E_alpha), ordered by increasing root height.
Write

\[
 u=\prod_{i=1}^N x_{\alpha_i}(u_i),\qquad
 \theta(tu)=h(t)+\sum_{i=1}^N F_i(t,u)E_{\alpha_i}.
\]

**Reconstruction theorem (written deduction).** The functions F_i are
uniquely recoverable from h(t), the root characters r_i=e^{alpha_i}(t),
and the positive-root brackets. Define

\[
 M=\prod_{i=1}^N\exp(u_i\operatorname{ad}E_{\alpha_i}),\qquad
 M h=h+\sum_i v_iE_{\alpha_i},\qquad
 M E_{\alpha_j}=\sum_i M_{ij}E_{\alpha_i}.
\]

Then the exact recursive formula is

\[
 \boxed{F_i=-\frac{r_i}{r_i-1}
       \left(v_i+\sum_{j<i}M_{ij}F_j\right).}
\]

Compute this in the rational function field before specializing t.
For a genuine equivariant regular theta, all apparent poles along
r_i=1 cancel, and the result belongs to C[T][u_1,...,u_N]. The calculation
does not need matrices of rho. Only the adjoint action on mathfrak b is
involved; for a fixed h it can be computed on the formal basis
(h,E_alpha1,...,E_alphaN), of size N+1.

Proof. Equivariance under a cocharacter contracting U shows theta(B)
is contained in mathfrak b and its toral projection at tu equals h(t).
Also equivariance with conjugating element g itself gives
Ad(g)theta(g)=theta(g). The operator M is lower triangular in increasing
root height, with M_ii=1. Since Ad(t)E_alpha_i=r_i E_alpha_i, the i-th
component of this fixed-vector equation is

\[
 (r_i-1)F_i+r_i\left(v_i+\sum_{j<i}M_{ij}F_j\right)=0.
\]

On the dense open subset of T where every r_i!=1 this system has a
unique solution by triangular substitution. The actual regular map
theta already supplies a solution, so the recurrence equals it there.
The resulting rational identities imply equality everywhere after
cancellation, including torus root hypertori. The construction of M uses
only [E_alpha,h]=-alpha(h) E_alpha and the constants in
[E_alpha,E_beta]=N_alpha,beta E_(alpha+beta). The exponentials terminate
because these operators strictly raise root height. This proves the
statement in arbitrary finite root type, with consistent root data.

The recurrence is a general reconstruction argument, not a new existence
theorem for an arbitrary supplied torus map. For arbitrary input h,
uncancelled poles can occur. Checking cancellation gives a regular map
on B; by itself it does not certify an extension to an equivariant map
on all of G, nor that the input comes from a representation.

The root-factor consequence remains

\[
 \frac{\partial F_i}{\partial u_i}=
 \frac{r_i\alpha_i(h(t))}{r_i-1},\qquad
 \det D_U\theta=\prod_i\frac{r_i\alpha_i(h(t))}{r_i-1}.
\]

The factors are interpreted after cancellation. Their zero loci determine
the singular vertical derivative; mixed lower-root terms determine the
actual rank there. All of this concerns t held fixed. It does not
classify the full B/G differential or the target branch/non-finite locus.

## Connection with the torus formula and normalization

Use torus coordinates z1,z2 as in the
[defining-coordinate note](BOREL_COORDINATES.md). Put D_i=z_i d/dz_i.
For a representation character chi, the logarithmic coordinate basis
H_i of mathfrak t has trace matrix

\[
 B_{ij}=D_iD_j\chi(1),\qquad
 h(t)=\sum_i (B^{-1}D\chi(t))_i H_i.
\]

This is the weight/trace-matrix formula in these coordinates. It is the
gradient of the character, not the logarithm of the character. The inverse
trace matrix is essential: unnormalized character gradients for different
representations cannot be compared as if their scalars were equal.
The normalization gives d theta_rho at 1 equal to identity. Repeating a
representation or adding trivial summands leaves theta_rho unchanged.
For a direct sum, trace matrices and character gradients both add before
the matrix is inverted.

For SL3, H1=diag(1,0,-1), H2=diag(0,1,-1), so if h has first two
diagonal entries h1,h2 its simple-root values are h1-h2 and h1+2h2.
For Sp4, in the basis (e1,e2,f2,f1), H1=diag(1,0,0,-1) and
H2=diag(0,1,-1,0), so its simple-root values are h1-h2 and 2h2.
These are epsilon/logarithmic coordinates, not simple-coroot coefficients.

## What is implemented

The [implementation](src/sage_borel/borel_reconstruction.py) has two levels:

* `reconstruct_root_coordinates` takes height-ordered roots, bracket
  constants, root characters, simple-root evaluations of h, and root
  parameters in a characteristic-zero rational function field. It checks
  data shapes, but assumes the supplied constants are consistent root
  data. It returns rational expressions before a regularity certificate.
* `reconstruct_borel_map` supplies the SL3 and Sp4 root-coordinate backends,
  derives h from a character or accepts h directly, and checks that the
  final denominators are torus monomials. It exposes positive-root and
  toral coefficients, full matrices, root factors, vertical Jacobians,
  exact rational point evaluation and triangular inverses off the
  root-factor zero locus.

Named representations are defining, dual, symmetric_square and adjoint
for both groups, and second_fundamental for Sp4. In Sp4 the dual equals
the defining representation and Sym^2 equals the adjoint representation.
The custom-character callback accepts other characters and mixtures.
It checks Weyl invariance, Laurent regularity and trace nondegeneracy;
it does not certify nonnegative irreducible multiplicities. Virtual or
rational character combinations must not be mislabeled as actual
representations merely because these algebraic checks pass.

The supplied-torus-map callback returns the first two diagonal entries
of theta(t). It must be expressed in the documented coordinates.
Successful cancellation is a Borel-regularity check, with the extension
qualification stated in the theorem. Singular fibers for general input
are not classified automatically. The earlier defining-only `fiber()`
method remains on `BorelSpringerMap` and is not reused for other
representations.

```python
from sage_borel import borel_springer_map, reconstruct_borel_map

old = borel_springer_map("Sp4")
rebuilt = reconstruct_borel_map("Sp4")
assert rebuilt.coordinates() == old.coordinates()

from_torus = reconstruct_borel_map(
    "Sp4", torus_map=lambda x, y: ((x-1/x)/2, (y-1/y)/2))
assert from_torus.coordinates() == old.coordinates()

adj = reconstruct_borel_map("SL3", "adjoint")
assert adj.vertical_rank((1, -1), (0, 0, 0)) == 1
assert adj.vertical_rank((1, -1), (2, 3, 4)) == 2
assert adj.inverse((2, 3), adj.coordinates((2, 3), (1, 2, 3))) == (1, 2, 3)

omega2 = reconstruct_borel_map("Sp4", "second_fundamental")
leading_factors = omega2.root_factors()
assert omega2.denominators_cancelled

defining = reconstruct_borel_map("SL3")
mixture = reconstruct_borel_map("SL3", character=lambda x, y:
                              defining.character + adj.character)
assert mixture.denominators_cancelled
```

## Representation-dependent collapse: two exact examples

For SL3's adjoint representation, at
t=diag(1,-1,-1), in the same root order as before, the recurrence gives

\[
 \theta(t)=0,\qquad (p,q,r)=(0,b/3,-ab/6).
\]

Thus the vertical rank is two if b!=0 and one if b=0. For q!=0 the
image has p=0, b=3q, a=-2r/q and c free, giving an affine line fiber.
For q=0 the image requires p=r=0, and a,c are free, giving an affine
plane. This is a complete description at this particular torus point,
not a full adjoint collapse-stratum classification. The defining SL3 map
has rank three at every point of every torus fiber.

For Sp4's five-dimensional V(omega_2), set x=z1 and y=z2. The coefficients
are

\[
 p=\frac{x+y}{2y}a,\quad
 q=\frac{y(x^2+1)}{2x}b,\quad
 r=\frac{xy+1}{2}c+\frac{ab}{2},\quad
 s=\frac{x(y^2+1)}{2y}d+\frac{x}{y}ac+\frac{x}{2y}a^2b.
\]

The cubic a^2 b term comes from the length-two root string: in our basis
[E_alpha1,E_alpha2]=E_(alpha1+alpha2) and
[E_alpha1,E_(alpha1+alpha2)]=2E_(2alpha1+alpha2).
No representation matrix was needed to obtain it.
The vertical determinant is

\[
 \boxed{\det D_U\theta_{\omega_2}
 =\frac{(x+y)(xy+1)(x^2+1)(y^2+1)}{16y}.}
\]

The old defining-map divisors x+y=0 and xy+1=0 persist, and there are
additional divisors x^2+1=0 and y^2+1=0. For example at x=i,y=2 and
u=1 the new rank is three, whereas the defining-map rank is four.
This rank comparison was checked exactly over Q(i). Point methods of
the ready API currently accept rational points; the symbolic formulas
can be specialized to Q(i) by a polynomial-ring homomorphism, as in the
tests. A complete fiber classification on intersections of all four
divisors is not asserted.

## Independent verification and known inputs

[Tests](tests/test_borel_reconstruction.py) verify the root brackets
from elementary matrices independently of the constants used by the
runtime. They then construct group and Lie-algebra representation
matrices directly: the dual, symmetric and exterior square tensor
actions, and the adjoint conjugation action. Full trace Gram solves
are compared as symbolic identities in every torus and root parameter.
There are nine such cases: four SL3 representations, three Sp4
representations, and a defining-plus-adjoint mixture in each group.
For Sp4, exterior square is V(omega_2) plus a trivial line; the trivial
differential is zero and does not affect trace projection.

The low-level recurrence is also checked against independent defining
SL2 and SL4 matrix formulas, using A1 and A3 root data built in the test.
Thus the recurrence is not confined to the two ready backends. These
finite checks corroborate the separate uniform proof; they do not prove
arbitrary-rank scope by extrapolation. Other checks include normalization,
regular extension across root hypertori, changed collapse loci,
triangular inverses, direct-sum normalization and invalid-input rejection.
The [verification report](borel-reconstruction-verification.json) binds
the proof, source code and tests. All arithmetic checks run inside Sage;
they are independent constructions, not external review or a second CAS.

Published inputs are the trace-projection construction and equivariance
in Kostant--Michor,
[The generalized Cayley map](https://arxiv.org/abs/math/0109066), and the
torus weight/trace-matrix formula in Rogers,
[An Explicit Determination of the Springer Morphism](https://arxiv.org/abs/1701.01538).
The fixed-vector triangular solve is a written elementary deduction from
these inputs and standard root-space theory. The formulas and new API
are explicit implementations; no literature-priority or novelty claim is
made. The source pages were checked on 2026-10-22; this was not an
exhaustive search for earlier coordinate algorithms.

## Next structural milestone

Connect the recurrence to the library's general character API and an
automatic Chevalley root-data adapter. The goal is one reconstruction
path across finite root types, with no representation-specific matrix
construction or hand-written bracket tables. Progress would mean exact
recovery in higher-rank A/C and a non-simply-laced exceptional example,
with the same cancellation and independent matrix checks where feasible.
This is more informative than adding a long table of unrelated formulas.
General singular-fiber classification and complete critical-image/branch
loci remain separate problems after the coordinate engine is available.
