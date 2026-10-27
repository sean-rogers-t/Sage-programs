# Adjoint Springer fibers: consolidated theorem and dependency audit

26 October 2026. GEO-91. This consolidates the completed case proof and
audits its scheme-theoretic implications. It is an internal mathematical
review, not external peer review, formal verification, or a priority claim.
The older notes remain dated records; their incomplete coverage statements
are superseded by the dependency table here and GEO-90.

## 1. The theorem, with all quantifiers

Let G be a simple, simply connected complex algebraic group, let g be its
Lie algebra, and let B be its Killing form. Define the adjoint
Springer--Cayley morphism theta:G -> g by

    B(theta(g),X) = tr(Ad(g) ad(X))   for every X in g.

Let G_reg be the open set of elements whose centralizer has dimension
r=rank G. For each proper Levi subgroup K, with Lie algebra k, put

    Z_K^o = {h in z(k) : Z_G(h)=K}.

**Theorem (completed by cases).** For every such K there exists a dense
Zariski open U_K in Z_K^o such that, for every complex point h of U_K,
the actual scheme

    F_h^reg = theta^{-1}(h) intersect G_reg

is smooth over C. Empty fibers are allowed. Equivalently, after a further
shrinking of U_K if needed, the family

    theta^{-1}(U_K) intersect G_reg -> U_K

is a smooth morphism. No reduction of its structure sheaf is understood.

The statement includes A, B, C, D, G2, F4, E6, E7 and E8. Low-rank
isomorphisms are counted in their simple types. Here B means Spin(2n+1),
C means Sp(2n), and D means Spin(2n), with the usual simple-rank bounds.
K=T is included; K=G and the zero-target fiber are not. The good open
may depend on K. Regular nonsemisimple source points are included.

Smoothness does not imply that all nonempty fibers have dimension
dim K-r: a transverse smooth factor may remain. On a contributing
incidence component Y_i the generic dimension is

    dim K-r + dim Y_i-dim z(k).

The E6 contact example has dimension 32 rather than 30. No general
nonemptiness, orbit-count, generic-degree or special-target classification
is part of this theorem. Nor is this regularity preservation of theta,
a theorem for other representations, or an integral Borel-image theorem.

## 2. Reduction to the full character, as schemes

The following is the common dependency of every case; see
[Theorem A of the critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md).

**Containment.** Equivariance gives Ad(g)theta(g)=theta(g), identically
as a morphism. Therefore

    (Ad(g)-1)h = -(Ad(g)-1)(theta(g)-h).

The stabilizer equations of h belong to the fiber ideal, so the actual
fiber lies in K=Z_G(h), not merely set-theoretically. Conversely theta(K)
lies in k. The Killing form restricted to k is nondegenerate: the zero
ad(h)-eigenspace is orthogonal to the other eigenspaces. Its restriction
to z(k) is nondegenerate as well.

**Coordinates.** Write S=[K,K], ell=rank S, and s=dim Z(K)^0, so ell+s=r.
The multiplication map S x Z(K)^0 -> K is a finite central etale cover.
The derived group S is simply connected. Its fundamental characters
x_1,...,x_ell give S//S=A^ell. Choose coordinates z_1,...,z_s on the
central torus and write Phi(x,z) for the FULL restricted adjoint character.
Put D_j=z_j partial/partial z_j and retain the actual schemes

    Y=V(Phi_x1,...,Phi_xell),
    mu:Y -> A^s,  mu_j=D_j Phi,
    Y_c=V(Phi_xi,D_j Phi-c_j).

The coordinates c_j=B(h,H_j), for the chosen central cocharacters H_j,
give an invertible complex-linear identification z(k)=A^s. Finite-index
changes of the cocharacter lattice are invertible here; they must still
be retained when counting covering degrees or points.

**Equality of ideals.** On K the Maurer--Cartan frame gives

    d chi_Ad - sum_j c_j d log z_j = B(theta-h,g^{-1}dg)

on the central cover. Thus the coefficient ideal of this differential is
exactly the fiber ideal. The logarithmic potential is a formal local
primitive, not a globally defined logarithm. On K_reg the adjoint
quotient is smooth; at a corresponding point q of Y_c,

    completed O_(F_h,g) = completed O_(Y_c,q)[[u_1,...,u_(dim K-r)]].

For a fiber point, equivariance also gives Z_G(g) contained in K and
hence Z_G(g)=Z_K(g). Regularity in K and G is therefore equivalent ON
THE FIBER. It need not be equivalent at an arbitrary point of K.

The finite central kernel acts freely on the central torus factor, so
the induced quotient-coordinate cover is etale. Smoothness descends.
This is different from later Weyl trace substitutions, which can ramify.

## 3. Why the generic-target step is valid

For a finite-type complex scheme Y and mu:Y -> A^s, it suffices to prove

    closure(mu(Sing Y)) is a proper subset of A^s,

where Sing Y is the nonsmooth locus of the ACTUAL scheme. Remove that
closed image, then apply characteristic-zero generic smoothness to the
remaining smooth source. Remove images of nondominant components as
well. There are finitely many components, so one dense target open works.
This argument requires no properness. Intersecting with Z_K^o retains a
dense open, since both are nonempty opens in the same irreducible center.

Conversely, smooth geometric fibers on a dense open imply this
nondominance after generic flatness and a further shrinking. In
characteristic zero the passage between the fiber formulation and the
smooth-family formulation in section 1 is consequently justified.

For equations g_1,...,g_m in a smooth ambient space of dimension m+s,
the locus defined by the equations and all m-rowed Jacobian minors is
a SUFFICIENT bad locus. A bound below s proves the required result:
off that locus the original scheme is smooth of dimension s. This does
not assume the equations form a regular sequence before checking rank.
With redundant equations or excess components this locus can be larger
than Sing Y, so the test is not necessary and must not replace the
quadratic argument in cases with smooth excess.

The [quadratic theorem](QUADRATIC_BRANCHING_FIBERS.md) uses a different
bound. At a nonsmooth point over a smooth compatible rank stratum of
dimension d, let t be the kernel dimension and p the additional mixed
Jacobian rank. The tangent dimension is t+s-p, while local dimension is
at least t+d. Hence d+p<s. Symmetry identifies p with the rank of the
central map along the affine solution fiber. Stratification then bounds
its bad target image by d+p. This includes changing rank and nilpotents;
it neither assumes Y reduced nor deletes the determinant-zero locus.

## 4. Complete dependency map

Every Levi is conjugate to a standard Levi. Explicit Weyl words identify
entire root subsystems, and equivariance identifies actual fiber schemes
and their centralizer-exact target opens. Matching abstract Dynkin types
alone is not used as a transfer certificate.

| Type | Proper standard subsets | Proof dependencies |
| --- | ---: | --- |
| A-D | All ranks | [Classical theorem](CLASSICAL_ADJOINT_FIBERS.md), the quadratic theorem, and its two square-substitution boundary exclusions |
| G2 | 3 | [G2 proof](G2_ADJOINT_FIBERS.md): torus and both root lengths of A1 |
| F4 | 15 | [F4 proof](F4_ADJOINT_FIBERS.md): ten quadratic characters and five cubic cases in three models |
| E6 | 63 | 49 quadratic, 3 root-length obstructions, 11 cubic cases; [assembly and explicit transfers](EXCEPTIONAL_QUARTIC_REDUCTION.md) |
| E7 | 127 | 88 quadratic, 7 root-length obstructions, all 32 cubic cases in the table below |
| E8 | 255 | 147 quadratic, 15 root-length obstructions, all 93 cubic cases in the table below |

For the E-type non-A cases the [root-length obstruction](EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md)
gives a derived fundamental-coordinate derivative identically 1. Their
relative critical ideals are unit ideals, proving actual regular-part
emptiness. This is stronger than a sampled lack of solutions.

| E7 cubic class | Subsets | Seed proof |
| --- | ---: | --- |
| A1^3 | 10 | [Binary criterion](BINARY_BRANCHING_FIBERS.md) |
| A2+A1^2 | 12 | [Binary criterion](BINARY_BRANCHING_FIBERS.md) |
| A1^4 | 2 | [Symmetric pencil, all rank strata](E7_SYMMETRIC_PENCIL.md) |
| A3+A1^2 | 3 | [Nonunit multiplier boundary](E7_MULTIPLIER_BOUNDARY.md) |
| A3+A2+A1 | 1 | [Unit multiplier and anchored product](ADJOINT_PENCIL_SCOPE.md) |
| A2^2+A1 | 3 | [Final families, section 1](E7_FINAL_FAMILIES.md) |
| A2+A1^3 | 1 | [Final families, section 2](E7_FINAL_FAMILIES.md) |

| E8 cubic class | Subsets | Seed proof |
| --- | ---: | --- |
| A1^3 | 21 | [Binary proof with full boundary](E8_A1_TRIPLE.md) |
| A2+A1^2 | 28 | [Binary proof and central root walls](E8_A2_A1_PAIR.md) |
| A1^4 | 7 | [Final cases, section 6](E8_FINAL_CASES.md) |
| A2^2+A1 | 8 | [Final cases, section 5](E8_FINAL_CASES.md) |
| A3+A1^2 | 10 | [Final cases, section 4](E8_FINAL_CASES.md) |
| A3+A2+A1 | 4 | [Final cases, section 3](E8_FINAL_CASES.md) |
| A4+A1^2 | 4 | [Original smooth surface, equation (46)](UNIFORM_ADJOINT_PROOF_TARGET.md) |
| A2^2+A1^2 | 2 | [Quasifiniteness, flatness and all six boundary exclusions](UNIFORM_ADJOINT_PROOF_TARGET.md) |
| A4+A2+A1 | 1 | [Coupled multiplier theorem](MULTIPLIER_BRANCHING_FIBERS.md) |
| A2+A1^3 | 8 | [Final degree-at-infinity proof](E8_LAST_ADJOINT_FIBER.md) |

The [fresh coverage audit](tools/audit_adjoint_theorem.py) rebuilds E6/E7/E8
roots in integer simple-root coordinates, enumerates every proper subset,
and replays all 136 cubic standard-subset transfers against the full root
subsystems. It matches them to the seed proof dependencies above, including
the original cubic-character archive. It does not rerun the much larger
Weyl-orbit census: completeness here follows from enumerating the standard
subsets themselves. The quadratic-versus-cubic root-degree check uses the
already proved minuscule branching shape, not an independent proof of that
root-theoretic lemma. G2/F4 character checks remain the linked case audits.

## 5. Boundary and certificate audit

**Unit eliminations.** Paired gradient elimination is an equality of
quotient rings. The central chain rule must include derivatives of all
coordinate changes and coefficient matrices. These extra terms vanish
modulo the eliminated relative equations; they cannot be omitted before
passing to that quotient. A merely generically nonzero determinant gives
only a chart, whose entire complement needs its own argument.

**Ramified traces.** Central torus isogenies are etale in characteristic
zero. A trace z+z^{-1} is not: on z=+/-1 the original logarithmic target
coordinate is zero. The analogous commuting-root walls force target
root hyperplanes. The final simultaneous-inversion trace surface has
four fixed points, at which both relevant central target coordinates
vanish. Only these justified target exclusions are used. No derived
root discriminant is inverted to discard regular nonsemisimple points.

There is an important distinction in the trace-surface arguments: mu need
not descend to the trace quotient. For example, if a=z+z^{-1}, then
mu_z=(z-z^{-1})f_a changes sign under inversion, although f descends.
In the final E8 case the dimension-two Jacobian-deficiency locus is
computed on the trace model. First pull it back along the FINITE map
from the original central torus (its dimension remains at most two),
and only then take its image under the original three-dimensional
target morphism. After removing the ramification image, etale pullback
makes the remaining actual incidence smooth. This is the precise
interpretation of the returning-to-the-group step; it does not require
an unproved target morphism on the trace surface itself.

**Flatness route.** For the E8 A2^2+A1^2 seed, the proof establishes
zero-dimensionality of EVERY allowed central-source fiber before using
the relative complete-intersection theorem. Generic reducedness plus
flatness then gives reducedness of the total algebra on that open.
All six complementary boundaries have nondominant target images. The
length drops in that model do not contradict flatness: the morphism
need not be finite over the central-source base. Finiteness of its
generic central-TARGET algebra is a separate dimension argument.

**Homogeneous certificates.** For integer homogeneous generators, ranks
of each graded coefficient matrix can only decrease modulo p. Thus the
characteristic-zero quotient Hilbert function is bounded ABOVE by the
finite-field one. Leading monomials of certified ideal members bound it
further above; the selected list need not be a full Groebner basis.
This direction of inequality is essential. Affine finite-field dimension
alone does not imply the characteristic-zero bound.

The E7 final membership proof additionally sandwiches Hilbert functions
between equal bounds to certify reverse ideal inclusion over Q. The E8
final case instead uses a degree budget: four dimension-four components
at infinity contribute 3+6+6+6=21, matching the homogeneous upper bound.
No dimension-four component can remain off infinity, so the original
affine Jacobian-deficiency locus has dimension at most two. Embedded
components cannot cancel top-dimensional degree; lower-dimensional ones
are still retained by the bad-locus exclusion.

The final E8 arithmetic has two exact finite-prime member replays and
independent native normal-jet and monomial-degree checks. The native
check is not a second implementation of the F4 ideal-membership engine.
These checks support the written geometric argument; they do not
mechanically verify it. This consolidation checks their saved dependency
bindings and adds a new coverage implementation; it does not claim to
have rerun every historical polynomial computation.

## 6. Controls, attribution and verdict

The Sp4 special-target triple point remains compatible with a GENERIC
theorem. The E6 contact excess remains compatible with SMOOTH fibers.
Altered-character examples and the actual reducible SL3 example retain
persistent nilpotents and rule out dropping the adjoint hypothesis.
An additional elementary audit control checks embedded nilpotents with
ideal (x^2,xy) over a free target coordinate: generic reducedness at
minimal primes alone does not suffice. The proof above avoids that error
by controlling the actual nonsmooth locus or proving full reducedness
with flatness.

Published inputs are separated from the application:

- [Kostant--Michor, sections 2.1 and 2.6](https://www.mat.univie.ac.at/~michor/kostant.pdf): the trace-projection map and centralizer compatibility.
- [Steinberg, Theorems 6.1 and 8.1](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf): fundamental-character coordinates and independence of their differentials at regular elements.
- [Generic smoothness in the target, Theorem 3.3](https://math.stanford.edu/~vakil/0506-216/216class5152.pdf): the characteristic-zero smooth-source step, allowing empty fibers without properness.
- [Stacks, relative complete intersections](https://stacks.math.columbia.edu/tag/00SW) and [cycles of coherent sheaves](https://stacks.math.columbia.edu/tag/02QV): flatness and generic-length multiplicities used in the two E8 arguments.

The incidence reductions, branching exclusions and case applications
are the repository deductions. Exact character identities, memberships,
leading terms and transfers are computational certificates. The source
check here verifies the cited foundational inputs; it is not a new
exhaustive search for a previously published version of this theorem.

**Audit verdict.** No additional gap was identified in the common
scheme reduction, the stated case coverage, or the characteristic-zero
certificate implications. The all-simple assertion remains established
by the recorded case proofs, with those finite computations as explicit
dependencies. A uniform root-theoretic proof remains unavailable. The
[audit report](adjoint-theorem-audit.json) records exactly what was checked.

**Recommended next step.** Prepare a compact package for independent
mathematical review, led by this statement and the final E8 degree proof.
Its purpose is to challenge the nonstandard arguments before pursuing
stronger claims. Success means an outside reader can reproduce the decisive
certificates and verify the scheme and boundary transfers. A new Levi
census or further sampled fibers would not strengthen this theorem.
