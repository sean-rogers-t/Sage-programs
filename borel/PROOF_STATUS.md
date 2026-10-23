# Proof and literature status of the integral image notes

The [classical adjoint fiber theorem](CLASSICAL_ADJOINT_FIBERS.md) closes
the B1/D2 cubic gaps and completes generic regular-fiber smoothness,
nonemptiness and dimension for every proper Levi in types A-D. In B/C/D
every GL block must have size at most four, with residual rank at most
two/one/three, respectively; these bounds are also sufficient generically.
The D statement uses n>=4; the D3 isomorphism is handled by type A.
Every nonempty generic regular fiber has pure dimension dim K-rank G.
A square-substitution criterion controls ramification, a quartic spectral
identity excludes dominating bad loci and excess quotient parameters,
and an explicit deformation at zero target proves dominance for every
permitted block pattern. The last step does not assume that the starting
source is G-regular: K-regularity suffices until the target is moved into
its proper centralizer stratum. Independent tensor-square matrices check
1,425 trace-form, 69 Springer trace and 1,425 full differential pairings;
symbolic checks cover all proof mechanisms without a rank sweep.
These are written deductions from standard inputs, not an established
novelty claim or an external referee report. Nonregular parts of fibers,
special targets, orbit counts and arbitrary representations remain separate.
No public API or PDF change accompanies this proof work.

The [quadratic branching theorem](QUADRATIC_BRANCHING_FIBERS.md) now proves
generic Morse-Bott behavior for every character quadratic in the semisimple
Levi fundamental characters, with no constant-rank hypothesis. For
X=V(Qx+f), the nonsmooth-locus target image has dimension below s: on a
compatible rank-k stratum of dimension d, its tangent dimension and the
mixed-derivative rank p force d+p<s. Generic smoothness then applies to
the remaining source. Contributing regular Springer fibers have dimension
dim K-rank G+ell-k-q, where q is the codimension of the component's
parameter image. For SL_N adjoint, generic regular fibers are nonempty
exactly when every target block has size one or two and there is at least
one singleton; with t double blocks their dimension is 2t. Excluded proper
Levis have no regular fiber points for any target in that stratum. This
does not classify wholly nonregular components. The adjoint is quadratic
on every type C Levi; in B and D the only cubic exceptions are residual
B1 or D2 with a nontrivial SL block. These were the remaining cases at
that checkpoint; the new classical adjoint theorem resolves them.
Eight exact scheme controls, 58 type A partitions, 24
dominance checks, 347 root restrictions and a Sage-free original-adjoint
matrix audit corroborate the written proofs. Standard inputs and deductions
are separated; literature priority is unverified and the public API is unchanged.

The [fundamental affine-branching classification](AFFINE_BRANCHING_CLASSIFICATION.md)
now determines exactly when the earlier affine generic-fiber criterion applies
for every fundamental representation of a simple complex simply connected
group and every proper standard Levi. The classical A-D formulas have
written all-rank proofs, including the D4 half-spin chirality exception.
The exceptional classification covers 3,373 pairs: 133 affine, 3,240 not
affine. Every positive non-torus case has a strict norm-gap certificate;
every negative case has an explicit Levi-dominant extremal-weight witness.
A separate Sage-free rational-arithmetic auditor checks all exceptional
certificates and coverage. Full-character restrictions independently check
2,453 classical and 1,462 exceptional pairs, with three additional native
branching comparisons. These are deductions from standard representation
theory and a complete finite exceptional proof, not a priority claim.
Nonaffinity does not prove singularity or nonreducedness of a Springer fiber.
Coefficient-zero codimensions and dominance still need to be determined
to compute dimensions and nonemptiness in newly covered affine cases;
the quadratic extension is now proved above. Cubic branching is the next gate.

The [affine Levi-branching criterion](AFFINE_BRANCHING_FIBERS.md) now proves
generic smoothness of the regular Springer fiber whenever the restricted
character is jointly affine-linear in all semisimple fundamental-character
coordinates on the central cover. If ell is the semisimple Levi rank and
a contributing coefficient-zero component has codimension q in the central
torus, the dimension is dim K-rank G+ell-q. The quotient-critical locus
is Morse-Bott. No smoothness assumption on the coefficient scheme is needed:
the Jacobian criterion bounds the target image of its nonsmooth locus below
the dimension of the target center. Generic smoothness applies to the
remaining incidence scheme. These are deductions from classical tools;
universal affinity for fundamental representations is not asserted.
Uniform consequences give generic regular fibers that are empty or smooth
of dimension m(m-1) for exterior^k of SL_N with one target eigenvalue of
multiplicity m<=min(k,N-k), and dimension 2p^2 for V_(omega_n) of Sp_(2n)
with zero block 2p<=2 floor(n/2) and distinct nonzero target pairs. This
includes generic reducedness in odd symplectic boundary strata and smaller
zero blocks. Nonemptiness in every admissible stratum, counts, exceptional
target values, nonlinear branchings and wholly nonregular components remain
unresolved here. Independent checks cover 158 A coefficient identities,
78 C ideals and 17 original-weight/Weyl-weight full differential ranks.
Singular and nonreduced coefficient controls corroborate the exclusions;
a nonlinear abstract potential shows persistent double structure can occur
without affinity, but is not claimed to be a fundamental character.

The [fundamental-spin counterexample](SPIN_FIBER_EXCESS.md) refutes the
candidate that every fundamental representation has isolated nondegenerate
quotient-critical points over generic semisimple targets in each Levi stratum.
For Spin_(2n+1), V_(omega_n), and every nonzero rank-two semisimple vector
target, K has derived group Spin_(2n-1). The corrected potential on its
central cover is (z+z^(-1))y-2^n lambda log z. Its critical ideal is
(z^2+1,2y+2^n lambda z), leaving n-2 independent quotient coordinates free.
The regular part of the full fiber is nonempty and smooth of dimension
n(2n-3), exceeding dim K-n by n-2. Spin_7 gives dimension 9 instead of 8.
The normal Hessian determinant is -4, so this does not refute generic
reducedness. A general unused-coordinate proposition explains the excess
whenever the corrected potential factors through fewer smooth quotient
coordinates. Original Clifford matrices independently check 244 trace
pairings and 9,724 differential entries in ranks 2..5, including exact
centralizers; Weyl characters are checked through rank 7. These are written
deductions from standard spin branching and quotient theory; priority is
unverified. Nonregular strata and universal generic Morse-Bott behavior
remain unresolved. Previous A/C stratum theorems are unchanged.

The [symplectic target-stratum theorem](SYMPLECTIC_TARGET_STRATA.md) now
treats every semisimple h with dim ker(h)=2m for the primitive middle
representation of Sp_(4m). If K=Sp_(2m) x product GL_(r_j), sum r_j=m,
there are exactly 2^m m!/(product r_j!) full-fiber components meeting the
regular locus. All have dimension 2m^2+sum r_j(r_j-1) and generic
multiplicity one; their regular fiber germs are smooth. There is no
exceptional nonzero target value producing a double component in these
strata. The first-block equations force a fixed squarefree reciprocal
spectrum, and a perfect mixed pairing eliminates every transverse variable.
For V_(omega_n) of Sp_(2n), an all-rank coefficient recurrence also proves
that any semisimple image of a regular element has kernel dimension at
most 2 floor(n/2). In odd rank n=2m+1, the maximal permitted kernel with
distinct nonzero target pairs leaves an explicit one-variable potential;
its root counts and multiplicities are not classified. These are written
deductions from known character/quotient inputs, not established priority.
Independent primitive matrices check 494 full Lie-algebra trace pairings
in Sp_4, Sp_6 and Sp_8, including all 14 external quotient assignments
for the even test ranks and two odd residual-equation checks. Symbolic
first-block ideal checks cover 136 rank/block pairs, with signed assignments
enumerated through m=5. Components and embedded structure wholly in the
nonregular locus remain outside the theorem; the previous Sp_4 embedded
thickening is unchanged. No generic-degree or central-quotient component
count is inferred.

The [simple-group branching note](SIMPLE_GROUP_CAYLEY_BRANCHING.md) proves
that fundamental representations of every simply connected simple group
have no regular points in their zero fibers, directly from Steinberg's
fundamental-character coordinates and trace differentiation. The statement
descends when the representation descends. Central points of an irreducible
map's zero fiber are reduced and isolated. A Levi branching criterion gives
necessary semisimple-target block bounds in A-D and a hook-Schur extension.
More substantially, for exterior^k of SL_(2k), on every target stratum
with centralizer S(GL_k x product GL_(r_j)), sum r_j=k, there are exactly
2 M full-fiber components meeting the regular locus, each of multiplicity
one, unless 1+I^2 lambda^2=0, when there are M components of multiplicity two.
Here M=k!/(product r_j!).
Here I=binomial(2k-2,k-1) and lambda is the eigenvalue repeated k times;
the remaining eigenvalues differ from lambda and have multiplicities r_j.
An exact quadratic after elimination proves the count in every rank.
Actual exterior matrices independently check 270 trace pairings in SL_4
and SL_6, in both regimes and every remaining target multiplicity partition.
Ordered-root orbits are independently enumerated through k=7.
Weyl weights check 20 B/C/D representations;
113 symbolic block derivatives and 78 hook expansions also agree.
These are written deductions from known inputs, not established priority.
The count is of positive-dimensional special-fiber components, not generic
degree. Target multiplicities strictly below k, analogous symplectic strata,
exceptional Levi classifications, and components wholly in the nonregular
locus remain open here.

The [general critical-potential theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md)
now treats semisimple targets across group types. For K=Z_G(h), the full
fiber is scheme-theoretically inside K and is locally the critical scheme
of the character minus a target-dependent sum of logarithms of K-characters.
Where K -> K//K is smooth with smooth quotient, the completed fiber is a
smooth factor times the quotient potential's Jacobian algebra. This holds
at regular fiber points for simply connected semisimple G. Isolated critical
points characterize collapsed quotient-fiber components; their Jacobian
lengths give generic multiplicities, with nonsingular Hessian equivalent
to multiplicity one. A general linear-in-invariants elimination lemma
recovers the A parity and C reducedness results. At nonregular points only
the critical-pullback statement persists; a PGL_2 example also shows why
regularity without the isogeny/smoothness qualification is insufficient.
Tensor-power scaling and defining SL_n fibers are explicitly attributed to
KM 2.8(3) and 6.2. Smooth character-level thickness and unbounded generic
multiplicity on each fixed semisimple group are elementary corollaries,
using generally reducible representations. These are written deductions
from established ingredients; novelty is unverified. Nonsemisimple targets
and a classification for irreducible characters remain outside the result.

The [symplectic component theorem](SYMPLECTIC_FIBER_COMPONENT.md) proves that
the target-centralizer sweep C=Z_G(theta(t)).(t U(H)) in the primitive
middle-exterior zero model is a full-fiber component of dimension
dim Z_G(theta(t))-2m and generic multiplicity one in every rank. Its
completed local ring is a power-series ring at regular-unipotent points.
A perfect mixed character pairing L=E E^T with determinant one eliminates
all transverse parameters. This uses known regular adjoint-quotient
smoothness and does not assume absence of target collisions. Exact
cyclotomic checks through m=40 find none, giving C=t U(H) in those cases;
the all-rank noncollision assertion is still unproved. In Sp_4 the full
fiber has an embedded infinitesimal thickening at the semisimple point,
with local ring C[[x,y,z,w]]/(wx,wy,wz,x^2+yz-w^2), w^2 nonzero and w^3=0.
The conjecture that the reduced formal germ equals C at every point is
proved for m=1 and on the regular locus in all ranks, but remains open
at higher-rank singular strata. Original primitive matrices verify 1396
differential pairings and ten symbolic pairings for the singular model.
Novelty has not been established.

The [symplectic realization theorem](SYMPLECTIC_CENTRALIZER_REALIZATION.md)
extends realization to genuine type C centralizers. In the fixed irreducible
primitive middle-exterior representation V_(omega_(2m)) of Sp_(4m), every
odd polynomial map on the nilpotent cone of sp_(2m) occurs, up to a nonzero
target scalar, as a translated centralizer map at a noncentral semisimple
point with centralizer Sp_(2m) x (C*)^m. This includes zero, giving a closed
2m^2-dimensional unipotent subvariety in a Springer fiber. The all-rank proof
inverts a palindromic/anti-reciprocal coefficient pairing and deforms a
squarefree root-of-unity spectrum. Independent actual primitive-kernel
matrices in dimensions 5 and 42 check 174 ambient trace pairings; coefficient
matrices are checked exactly through m=10. Classical character identities,
Cayley coordinates and nilpotent centralizers are inputs; novelty is
unverified. The realization argument alone does not give the component
theorem above. These are complex group-map statements, not integral
Borel-image claims.

The [middle-exterior full-fiber theorem](MIDDLE_EXTERIOR_FIBER_COMPONENT.md)
now identifies the zero model's component and generic scheme structure in
every rank m>=2. For even m, t U(H) itself is a component of dimension
m^2-m and is generically reduced. For odd m, the larger target centralizer
sweeps it into a component of dimension m^2-1 and generic multiplicity two;
the original collapsed variety is not a component. A scheme-theoretic
centralizer reduction and cyclic-block elimination leave 1-s^(2m) in
even m and (s^m-1)^2 for odd m. The completed local rings on the
regular-unipotent locus are explicitly proved. Independent exterior traces
check 5419 derivative pairings; exact curve, cyclotomic-collision and
higher-rank differential checks agree. Known character and regular
adjoint-quotient facts are inputs, and novelty is unverified. The local
fiber at nonregular unipotent points, other components and global fiber
classification remain outside this theorem.

The [irreducible-centralizer realization theorem](IRREDUCIBLE_CENTRALIZER_REALIZATION.md)
now proves that the fixed irreducible representation exterior^m of SL_(2m)
realizes every equivariant nilpotent polynomial map on an sl_m centralizer,
up to a nonzero target scalar, by varying a noncentral semisimple point.
The centralizer is exactly S(GL_m x (C*)^m). This includes all orders and
the zero map, whose ambient fiber contains a translated unipotent variety
of dimension m^2-m. The proof inverts a Pascal coefficient transform and
avoids a finite discriminant set. Full exterior matrices verify 500
ambient tangent pairings and 35 further symbolic SL6 pairings, with exact
determinants and centralizers. This is an arbitrary-rank written deduction
from established character identities, not a priority claim. The representation
has kernel mu_m and is infinitesimally faithful as required by KM; it is
faithful on the corresponding central quotient. Full fiber components and
transverse geometry remain unresolved here.

The [type A character-degeneration theorem](TYPE_A_CHARACTER_DEGENERATION.md)
now settles a uniform representation-realization family. For every exterior
mixture on SL_m the central translated map has an explicit binomial formula.
At -I for even m its order is m-1 minus the degree of the interpolation
polynomial of the nontrivial exterior multiplicities. All orders 1,...,m-1
occur in faithful representations; the tensor square of the full exterior
algebra realizes the zero map. Twisted primitive trace indices give the
general type A centralizer coefficient criterion. These are written
deductions using published character/index identities, with no priority
claim. Independent wedge matrices check 140 moments and 203 central
comparisons; full symbolic SL4 trace projection and an explicit SL2 tensor
square agree. The next gate concerns irreducible ambient representations
restricted to proper centralizers, not additional central mixture examples.

The 2026-10-23 [fiber-program assessment](CAYLEY_FIBER_ASSESSMENT.md)
supersedes the proposed F4 rank-table next step. KM 4.12 already defines
the translated centralizer map, and KM 4.13 supplies an isomorphism
criterion. The contracted-product fiber reduction and type A matrix-power
normal form are elementary deductions from known theory, with explicit
proofs here and no novelty claim. Actual Springer maps can have reducible
fixed-t fibers: standard plus dual for SL_(4k), at t=iI, reduces to the
matrix square map; an independent SL4 trace/ideal check verifies two
four-dimensional components. This refutes a universal affine-fiber
claim. Square-zero geometry is already studied by Melnikov. The remaining
candidate is representation-realization of degenerations, not the
abstract existence of nilpotent polynomial maps; its literature openness
is unverified. More group rank tables and the F4 rank API are paused.

The [G2 vertical-rank theorem](G2_BOREL_RANK.md) gives all ranks for the
minimal and adjoint Springer maps on tU, using a centralizer normal form
on root hypertori. Its general reduction is
rank dF_t at tz = dim([u,theta(tz)] + im d(theta|tU_t) at z).
This is a written consequence of equivariance and unipotent conjugation,
not a claimed new general theorem in the literature. Off the root
hypertori the rank is constant on tU and counts nonzero root factors.
The G2 table includes all three involutions, both order-three points,
six additional adjoint points, and the single-root curves. Full symbolic
trace calculations in explicit 7-dimensional weight-string matrices and
14-dimensional adjoint matrices agree with the recurrence. Exact minor
ideals check the normal-form rank strata; no floating-point inference
or novelty claim is made. The result concerns the vertical derivative,
not the full B/G differential, fiber schemes or target branch loci.

The [automatic root Borel API](ROOT_BOREL_API.md) connects the general
SpringerCharacterMap to finite irreducible Chevalley root data, retaining
fundamental-weight torus and simple-coroot toral coordinates. Sparse
recurrence evaluates on all rational torus points. A written finite-jet
lemma proves that highest-root-height plus one precision suffices at
root hypertori, including t=1; no numerical limiting argument is used.
Full symbolic coordinates are separately budgeted. Independent standard
matrix checks cover A3/B3/C3/D4; full adjoint trace checks cover all five
exceptional types, and conjugation checks cover G2/F4/E6/E7 minimal maps.
A dense E8 adjoint case with all 120 root parameters nonzero agrees with
the separate 248-dimensional trace calculation. Exceptional checks share
Sage's Chevalley constants; the independence claim concerns the different
computations, not an external review or separate root-data implementation.
The recurrence and jet lemma are deductions from known equivariance and
root-space inputs, with no novelty claim. Full singular-fiber/branch-locus
classification, affine groups and other isogeny-form adapters remain open here.

The [torus/root reconstruction theorem](BOREL_RECONSTRUCTION.md) is a
uniform characteristic-zero written deduction: equivariance forces
Ad(tu)theta(tu)=theta(tu), and root height makes this a triangular system
with prescribed toral part theta(t). For a genuine regular map its
apparent root-character poles cancel. The recurrence needs positive-root
brackets and the torus formula, without representation matrices. Ready
backends cover SL3/Sp4 with named or supplied characters and supplied
torus maps. Nine all-variable full-representation trace solves, additional
SL2/SL4 low-level matrix checks, normalization and collapse checks agree.
The SL3 adjoint map collapses at diag(1,-1,-1), with ranks two or one;
Sp4 omega_2 adds the divisors z1^2+1=0 and z2^2+1=0 to the defining
vertical critical locus. These are explicit deductions and an API,
not novelty claims, a general singular-fiber classification, or an
automatic adapter for every finite root type. Custom-character positivity
and global extension of arbitrary torus-map input are not certified.

The [Borel-coordinate calculation](BOREL_COORDINATES.md) implements the
defining SL3 and Sp4 Springer maps on B=T U. The SL3 fixed-t map is always
a polynomial isomorphism. The Sp4 vertical determinant is
z1(z1 z2+1)(z1+z2)/4; complete image/fiber formulas distinguish affine-line
collapse from affine-plane fibers at the two crossing torus points with
b=0. At a crossing with b!=0 the mixed ab term keeps rank three.
The crossing image is constructible and not closed. A general triangular
argument derives the diagonal factors from the torus map by equivariance.
These are characteristic-zero written deductions from standard
trace-projection theory, checked by independent symbolic full-Lie-algebra
Gram solves and exact fiber parametrizations within Sage. No novelty,
external review, arbitrary-representation implementation or full B/G
critical/branch-locus classification is claimed.

The [arbitrary-rank B/C mixture proof](BC_MIXTURE_DEGREES.md) treats a
defining and b adjoint summands, with nonnegative multiplicities and
a+b>0. On SO_(2n+1), b=0 gives degree 2^n; for b>0 the degree is
2^n(2^n-n), minus 2n binomial(n-1,k) precisely when
a=(4k-2n+1)b for 0<=k<=n-1. The positive ratios are 1,5,...,2n-3
in even rank and 3,7,...,2n-3 in odd rank. On Sp_(2n), b=0 gives
2^n and every b>0 gives 4^n, with no boundary corrections. These are
written characteristic-zero deductions: crosspolytope Morse points
give the B corrections, all hypersimplex subfaces are excluded, and
the C complete homogeneous quadratic excludes zeros on every face.
Spin pullback doubles the B degrees; pure adjoint descent to PSp halves
the C degree. Independent character/polytope checks through rank six,
sign counts through rank ten, 68 reduced original fibers in two models
and two additional isogeny controls corroborate the proof. Finite tests
do not establish its arbitrary-rank scope. Known inputs are KM degree
reduction, classical root-polytope geometry and toric intersection
theory. Novelty and external review remain unverified. No runtime API
changed. Explicit torus curves collapse to zero at every exceptional
B ratio, proving both torus and quotient maps are then non-finite;
outside those ratios they are finite flat. Next is the full non-finite
target locus, to go beyond this zero-fiber witness.

The [F4 minimal-adjoint mixture proof](F4_MIXTURE_DEGREES.md) gives the
complete formula for nonnegative multiplicities a,b with a+b>0: degree
96 if b=0; for b>0, degree 240 at a=b, 336 at a=5b, and 384 otherwise.
All three degrees (group, torus, quotient) agree. The long-root polytope
has normalized volume 384 in the actual F4 lattice. Its 24 facets have
primitive height two, with three Morse zeros each at a=b and one each
at a=5b. The first normal layer and its tangent gradient vanish there,
so each point contributes length two. Lower faces have no zeros.
This is a characteristic-zero written deduction from standard toric
intersection theory. Twenty reduced finite-field fibers from two
independent equation constructions corroborate it, including agreement
of multiplication characteristic polynomials at matching targets.
The reusable height-two lemma is stated with hypotheses; arbitrary
representations, global finiteness on the exceptional ratios, literature
priority and external review remain separate. No runtime API changed.
The next goal is a coefficient-discriminant criterion for families with
Morse facet zeros, using normal jets to determine degree jumps.

The [F4 minimal degree proof](F4_MINIMAL_DEGREE.md) gives degree 96 for
V(omega_4), its torus restriction and its adjoint quotient. All 24 facets
are octahedra with no singular zero of their equal-coefficient characters;
every lower face is a simplex. The actual lattice gives normalized volume
96 and no boundary correction. Independent highest-weight/explicit weights,
determinant volumes and four reduced specialized fibers in two equation
models corroborate the proof. The common multiplication characteristic
polynomial agrees across the models at each prime. The general no-boundary
criterion also gives global finite flatness on the torus; the F4 quotient
is finite flat as well. The full group map is not finite: an explicit
noncentral zero has a 16-dimensional conjugacy class. Odd-crosspolytope
coefficient hypotheses and the use of known properness/flatness results
are stated explicitly. Novelty and external review remain unverified.
Its proposed minimal/adjoint multiplicity follow-up is completed above.

The [E6 adjoint deformation proof](E6_ADJOINT_DEFORMATION.md) now
reproduces the boundary correction 1890 independently of the earlier
Segre formula: 810 branches approach open curve strata and 1080 approach
triple points; two-facet crossings receive none. The middle D5 root layer
has value -8 and zero full gradient on the facet critical locus. At a
crossing, tangent elimination gives two equations st U_i+eta V_i, whose
unit combination excludes eta nonzero. The actual root lattice, index-two
facet cover and final central factor three are kept separate. Fresh root,
facet, layer and local-algebra checks support this written deduction.
The degree 2754 is an existing repository result, not a new value. Two
local transfer principles are stated with hypotheses; arbitrary boundary
singularities, external review and priority remain separate. Its proposed
F4 minimal test is completed in the proof above.

The [general degree audit](GENERAL_CAYLEY_DEGREE.md) separates standard
toric degree criteria from the special half-spin product argument.
It proves a sufficient no-boundary criterion and an isolated-Morse
normal-layer formula, both deductions from known intersection methods.
An independent E6 weight-orbit/facet reconstruction confirms the already
recorded minimal degree 1080. Five exact generic fibers over QQ(b1,b2)
check examples showing that neither unweighted support nor facet
characters alone determine degree. The half-spin argument was reinspected
without a contradiction found; its positive-dimensional correction is
not asserted for arbitrary representations. The E6 adjoint transfer test
is now completed in the proof above. This adds no novelty
claim, external review, exceptional degree value or runtime API.

The [half-spin rank formula](HALF_SPIN_RANK_DEGREES.md) now gives a
written characteristic-zero proof of degree 2^(r-1)(r^2-3r-2)+4r
for either individual half-spin Springer map on Spin_(2r), r>=3.
The proof deforms the weight coefficients, counts branches on every
facet critical stratum using two product-gradient systems and a
generic cube polynomial, and excludes every deeper face by a
Puiseux-order argument at a general target. Thus Spin_12 has degree
536, with correction 840 from volume 1376. Four exact reduced
Spin_12 fibers, covering two coordinate models, two primes and both
chiralities, corroborate it. Local algebra and combinatorial checks
are separate from the all-rank proof. The actual Spin lattice is
retained, including the index-two representation kernel in even rank.
No novelty, external review or runtime API change is claimed. An
effective inverse-fiber presentation and the branch/non-finite loci
remain open here. The earlier degree notes below retain their own
scope; GEO-20 supersedes their higher-rank degree questions.

The [half-spin degree note](HALF_SPIN_DEGREE.md) proves generic degree
148 for either individual half-spin Springer map of Spin_10, including
its torus restriction and adjoint quotient. This is a characteristic-zero
written deduction: the actual Spin weight polytope has normalized volume
208, and ten coordinate facets each contribute six simple boundary
intersections; all remaining boundary strata have none. Independent
checks use the rational facet critical algebra, symbolic local Jacobians,
and twelve reduced finite-field fibers in two original torus models.
The finite-field counts corroborate rather than prove the complex degree.
Known inputs are the toric intersection formula and Kostant--Michor's
degree reduction. Priority and external review are unestablished.
The next goal is an excess-intersection correction in Spin_12, where
the facet critical loci become positive-dimensional, toward a rank law.
The earlier notes below retain their historical scope; their statements
that degree was unresolved are superseded for Spin_10 by this result.

The [half-spin boundary note](HALF_SPIN_BOUNDARY.md) now completes the
underlying zero fiber for either half-spin of Spin_(2r), r>=3. Besides
the four central points and the preceding chart classes, it is exactly
the Spin inverse image of the SO matrices having at least four ordinary
eigenvectors at each of +1 and -1. Jordan decomposition and the chiral
character restriction reduce the boundary condition to vanishing
Pfaffian cofactors on the two unipotent sign blocks. Standard nilpotent
centralizer dimensions then prove dim theta^(-1)(0)=r(2r-1)-12 for r>=4.
Thus Spin_10's zero fiber has dimension exactly 33; all seven boundary
Jordan strata are listed. Independent exact chiral-matrix tests include
nonsemisimple points and negative controls, with separate weight,
commutator-kernel and partition/sl_2 dimension checks. These are uniform
written deductions from published inputs, with no established priority
or external review. Boundary scheme structure, component incidences,
generic degree and maximal finite target locus remain separate questions.
The next priority is generic degree; further zero-fiber scheme analysis
is not a prerequisite for it. The preceding checkpoint below retains
the narrower scope of its own chart note.

The [individual half-spin note](HALF_SPIN_CAYLEY.md) treats S_sigma of
Spin_(2r), r>=3, with its actual trace normalization. On the full-spin
Cayley chart it proves theta_sigma=2c(X+sigma J(X)), where J is the
Pfaffian cofactor. Its scheme zero fiber in that chart is smooth:
two central points and 2(r-3) semisimple conjugacy classes of dimension
r(r-1). A separate boundary family has dimension r(2r-1)-12 for r>=4.
Thus Spin_10 has four 20-dimensional chart classes and a 33-dimensional
boundary zero family; no maximality assertion is made for the latter.
The all-rank arguments are written deductions from standard character,
Pfaffian and Cayley-chart inputs. Independent checks include 12
highest-weight/moment comparisons, 18 full chiral Clifford calculations,
42 exact zero witnesses, 12 full differential ranks and ten boundary
checks. Known low-rank identifications give degrees 4 and 32 for Spin_6
and Spin_8. The subsequent boundary note completes the underlying zero
fiber; Spin_10's generic degree and maximal finite target locus remain
unresolved. Priority and external review are not
established; the scalar full-spin finite-locus theorem is not transferred.

The [Spin geometry note](SPIN_CAYLEY_GEOMETRY.md) proves that for N>=5
the full polynomial-spin Springer map is finite precisely over the
complement of a_r(Y)=0, where a_r is the last even characteristic
coefficient. It includes the exact N=3,4 exceptions, a full normalized
Jacobian valid across the Cayley boundary, and corank-one regular
semisimple critical points. The essential finiteness obstruction is
the rank-at-most-two boundary image, combined with lying over for the
integral scalar model. A Spin_4 calculation demonstrates a false branch
factor produced by using that model outside its valid chart. Exact
checks include 48 weight/root Jacobians, 33 Clifford boundary Jacobians,
32 fiber-algebra derivative identities, eight boundary target matrices
and 42 regular critical invariant tuples. This is a uniform repository
deduction using published KM and standard algebra inputs; no priority
claim, external review, complete exceptional-fiber classification or
adjoint-quotient finite-locus theorem is asserted.

The [polynomial spin-character note](SPIN_POLYNOMIAL_CAYLEY_DEGREES.md)
extends the Spin tensor law to every actual character F(chi_S), with
S as specified below: degree 2 floor(N/2) deg(F). Character differentiation
and the KM double-cover chart give an explicit scalar fiber polynomial;
its leading term cannot cancel. It proves finite flatness over a stated
open target and descent to SO exactly for even F, with half the degree.
Within this family direct sums take the maximum degree and tensor
products add degrees. Nonzero roots of F' give collapsed character
levels over zero; S direct-sum S tensor S collapses chi_S=-1/2.
Independent checks include 48 scalar fibers, 24 weight-fiber algebras,
14 original Spin Laurent fibers, five SO fibers and six direct
collapsed-level witnesses. Actual negative-coefficient examples are
justified by removing a trivial summand from S tensor S. These are
written deductions from known inputs, not an all-representation formula,
novelty claim, external review or runtime API extension.

The [Spin tensor-power degree note](SPIN_CAYLEY_DEGREES.md) deduces
deg(theta_(S tensor m))=2m floor(N/2) for Spin_N, N>=3 and m>=1.
Here S is irreducible spin in odd dimension and the sum of both half-spin
representations in even dimension. It supplies a monic scalar fiber
presentation over an explicit open target and proves finite flatness
there. Even powers descend to SO_N with half this degree; odd powers do
not descend through the double cover. The argument extends the published
KM Spin Cayley chart and m=1 degree proof using their tensor identity;
it is not a claimed resolution of a published open problem. Independent
checks comprise 48 scalar polynomials, 16 weight-moment fiber algebras,
11 original Spin Laurent fibers and five original SO Laurent fibers.
No single-half-spin formula, global finiteness, external review,
established novelty, minimum Cayley degree or runtime API change is claimed.

The [symplectic generalized Cayley degree note](SYMPLECTIC_CAYLEY_DEGREES.md)
proves deg(theta_F)=2n deg(F) for Sp_(2n), n>=1, whenever F(det(I+g))
is an actual nontrivial character. It includes tensor powers of the
total exterior algebra, with degree 2nm, and polynomial mixtures with
negative coefficients when the resulting representation is actual.
The proof uses Cayley coordinates and the degree of a reduced rational
function on a generic radial line. It supplies an explicit scalar
fiber polynomial, proves surjectivity, and proves finite flatness of
degree 2n deg(F) over invertible target matrices. Within this character
subring, direct sums take the maximum of degrees and tensor products
add degrees. These are short deductions from standard ingredients;
KM's tensor-scaling identity and analogous Spin fiber method are
published inputs. Exact corroboration includes 80 scalar fibers, 24
independent weight-moment fiber algebras, ten original Laurent Groebner
fibers, 108 singular-target polynomial checks and nine full Sp4
exterior-matrix checks. No general representation degree law, minimum
Cayley degree, external review, novelty or global finiteness is claimed.
Current research remains confined to generalized Cayley maps; the earlier
suggestion of constructing an unrelated global Sp4 harmonic complement
is shelved in accordance with the user's scope clarification.

The [complete simple-group regularity note](REGULARITY_CLASSIFICATION.md)
now closes the type D, exceptional and arbitrary reducible A gaps in the
table, for actual nontrivial finite-dimensional rational representations
over C. It is a uniform written deduction with internal exact corroboration;
external review and publication priority remain unverified. When roots
are primitive in the weight lattice, positivity on real root walls and
Laurent factorization force every root coordinate to be a normalized
two-monomial difference. Root strings contradict this in B (rank>=3),
F4 and G2; Dynkin forks contradict it in D (rank>=4) and E6/E7/E8.
For SL_N, N>=3, positive weight support then leaves only copies of
defining or copies of dual defining, plus trivial summands, with no
mixture of the two. These known defining/dual maps preserve global
regularity. SL2 has the additional torus-positive determinant-power
characters but only its defining-copy family is globally positive.
Together with GEO-11 for Sp, this completes every simply connected
simple row; a descent argument covers the other simple isogeny forms.
No nontrivial central quotient admits a positive case. Exact corroboration
checks 1,942 roots in 32 root systems, both local endpoint patterns,
48 independent GT/Freudenthal A characters, 30 exact D/exceptional
witnesses including E8, four A-mixture witnesses and eight rank-one
controls. These computations supplement the all-rank proof; they do
not exhaust all representations. Earlier statements below leaving the
simple-group rows unfinished describe historical proof scopes. The
classification makes no claim about general reductive/product groups,
virtual characters, characteristic p, minimal localization of arbitrary
negative maps or publication novelty. No runtime API changes.

The [symplectic harmonic localization note](SYMPLECTIC_HARMONIC_LOCALIZATION.md)
proves that the determinant-family maps theta_m give multiplication
C[G]^G[D^(-1)] tensor theta_m^*Harm(g) -> C[G][D^(-1)] as a module
isomorphism for D=det(1+g), G=Sp_(2n), n>=2. This sufficiency is an
explicit corollary of the known graded Kostant separation theorem and
the classical Cayley isomorphism: a degree-d harmonic is multiplied by
the invariant unit 2^d delta^(-md). For any nonzero invariant f the
same multiplication is an isomorphism iff D divides f. Necessity follows
from explicit conjugate pairs on D=0 with equal Springer images but
different matrix entries. The full normalized Jacobian is
(D/4^n)^((nm-1)(2n+1)) times (1-nm+(m/2)tr(g)); its additional trace
factor need not be inverted for separation. Exact checks cover 48 full
matrix differentials, 12 independent exterior/tensor trace differentials,
12 regular critical corank-one points, 15 conjugacy pairs, and genuine
Sp4 harmonic bases through degree three. These corroborate a uniform
written argument, not a new general separation theorem, external proof
review or established novelty. Inverse-differential vector fields still
need the additional critical factor in general. No runtime API changes.

An additional internal audit of that note's Sections 2-3 found no gap
under the stated hypotheses. Section 7 proves the exact pointwise locus
theta_m^(-1)(g_reg)=G_reg intersect {D!=0} for n>=2. It adds 48
centralizer-subspace torus checks, 36 boundary-stratum checks, nine
actual exterior/tensor trace checks on the Sp4 boundary and three Sp2
rank controls. Section 8 records the important published distinction:
Richardson already gives some global complement for every simply
connected semisimple group, including Sp and the exceptional groups.
The regularity classification concerns only Springer-pulled harmonics.
Agricola--Goodman's SL2 differential-operator construction is known and
uses the same defining pullback; GL defining and orthogonal classical
Cayley-open separation follow by localizing Kostant's theorem. The
orthogonal classical Cayley map is not identified here with its defining
trace-projection Springer map. The earlier proposed explicit global
alternative for Sp4 is outside the user's current generalized-Cayley
scope and is not being pursued. It was not an unresolved existence
claim or an established new open problem.

The [mixed-character regularity classification](MIXED_REGULARITY.md)
now treats every nontrivial finite-dimensional rational representation of
Sp_(2n) and Spin_(2n+1) at n>=2, without a central-scalar hypothesis.
For Sp, torus regularity holds precisely for actual characters
A det(1+g)^m+B, with positive integers A,m and an allowed trivial-character
adjustment B. These arise from total exterior algebras and tensor powers.
Their normalized maps are polynomial multiples of the classical Cayley
transform. They still fail global regularity at an explicit regular element
with a size-two Jordan block at eigenvalue -1. For Spin of type B at
n>=3, positive-cone character derivatives and even-sign symmetry exclude
torus regularity even for mixed sums. Spin5=Sp4 has the torus-positive
symplectic exception, and still no globally positive representation.
Exact checks include nine symplectic tensor characters, 27 torus traces,
five nonsemisimple boundary ranks, actual exterior matrices against 31
Lie basis elements, and 18 mixed B witnesses evaluated directly in finite
algebras. These are uniform written deductions, not externally reviewed
or established novel results. The exterior/Cayley identity is elementary
from known inputs. Earlier statements below leaving mixtures unresolved
describe the scope of those earlier proofs; they are superseded here.

The [type B regularity proof](TYPE_B_REGULARITY.md) now excludes every
nontrivial irreducible Spin_(2n+1) representation at n>=2, including
spinorial highest weights. It also excludes all same-central-parity
Spin sums and every nontrivial SO_(2n+1) representation. At the regular
Spin torus lift with auxiliary coordinates (i,2,3,...,n), inversion of
the first coordinate is the central -1 translate. The central character
forces y_1=0 or y_2=...=y_n=0, giving a short-root zero. A direct
integer/half-integer character factorization independently proves those
zeros. Exact Freudenthal/Weyl-alternant comparisons cover 55 irreducibles,
with 13 exterior-monomial checks, seven spin sign-cube comparisons and
eight reducible controls. The basic spin formula and trace equivariance
have published antecedents in Kostant--Michor; priority for the uniform
classification is unverified. Mixed-parity Spin sums remain unresolved
here; the rank-one defining spin representation is a positive exception.
This is separate from the shelved integral Spin-kernel project.

The [full type C regularity proof](TYPE_C_REGULARITY.md) now excludes
every nontrivial irreducible Sp_(2n) representation at n>=2. A single
regular tuple (2,-2,3,-3,...), appending i in odd rank, is taken to its
central -I translate by a signed permutation. The central character
therefore forces y_1=+/-y_2. This argument also covers all representations
with a single central scalar and, by pullback, every nontrivial PSp
representation. Mixed central parities on Sp are outside its scope.
An independent Sp2-branching/positive-torus proof treats every nonfundamental
highest weight and the negative rank-one cases. The rank-one defining
case is positive as already known for SL2. Exact alternant/Freudenthal
character checks cover 79 irreducibles, with 73 universal Q(i) trace
witnesses, 68 polynomial-root witnesses, 30 direct symmetric-monomial
checks and eight same-parity reducible controls. The general proof is
not extrapolated from this finite scope. Trace equivariance, character
theory and the harmonic criterion are known inputs; priority and external
review remain unverified. Type B is addressed separately above; degrees and integral
images are not classified by this theorem.

The [type C fundamental regularity proof](TYPE_C_FUNDAMENTAL_REGULARITY.md)
excludes every omega_k representation of Sp_(2n) at n>=2. These are
primitive exterior kernels. A positive multiaffine expansion in inverse-pair
sums yields a bounded rational construction of a smaller regular torus
character zero; the full input uses quadratic algebraic entries and has
y_1=0. The defining case uses the existing Sp4 repeated-coordinate
witness extended by generic positive entries. The representation and
character inputs are classical; the defining map is already Part II
Lemma 5, and the harmonic consequence uses KM Theorem 3.2. Independent
contraction matrices in 19,676 weight blocks, Weyl alternants, complete
trace coordinates and exact quadratic-algebra witnesses corroborate all
36 fundamental cases through rank eight (35 negative and one rank-one
control). The all-rank proof is not extrapolation from these checks.
Arbitrary type C highest weights, novelty, other root systems, degrees
and integral images are outside this theorem's scope.

The [irreducible type A regularity proof](TYPE_A_REGULARITY.md), added
22 October, completes every nontrivial SL_n highest weight: torus and
global regularity preservation hold exactly at omega_1,omega_(n-1).
A nonnegative divided-root character expansion gives a nonconstant
algebraic witness polynomial whenever sum_i a_i>=2; the earlier exterior
argument covers the remaining fundamental weights. It includes the
symmetric-power criterion d=1. Character facts, the trace construction,
KM's standard case and harmonic criterion are published inputs.
Independent GT/Jacobi--Trudi character and derivative calculations, plus
full trace solves, corroborate 149 representations through n=6, including
50 symmetric powers and 134 algebraic-witness polynomials. These checks
do not replace the all-rank proof or certify novelty. Reducible cases,
other groups, generic degrees and integral images remain outside its scope.
The stated central-isogeny corollary also excludes every nontrivial
irreducible representation on SL_n/C when C is a nontrivial central
subgroup; neither positive endpoint descends. This is an elementary
corollary, not a classification for other root systems.

The [research ledger](RESEARCH_LEDGER.md) indexes the program's uniform
results, finite cases, conditional routes, counterexamples and open questions.
The detailed statements and source audits below retain their original scope.

Audit date: 2026-10-19. Scope: the type B/D Smith notes, the subsequent
type C joint-image argument, and their
normalization, lattice, geometric, and symmetric-function dependencies.
The exceptional section records the subsequent G2/F4 and E6/E7 calculations,
duality deduction and explicit primitive-image subalgebras.
This is a repository proof/source audit, not external peer review or a
formal proof certificate. It does not certify originality. It does not
re-audit every article or the earlier dissertation proof in this repository.

The [22 October exterior-power regularity proof](EXTERIOR_REGULARITY.md)
classifies all SL_n fundamental representations: torus regularity and
global regularity preservation hold exactly for the defining and dual
representations. A uniform rational construction obstructs every interior
power. The standard case is published in KM Section 6.2; the trace formula,
dual endpoint argument and interior cancellation give the recorded deduction.
The proof includes the nonsemisimple endpoints, not just torus tests.
Actual-weight trace solves independently corroborate 66 representations
through n=12 and 45 rational counterexamples, plus eight full-matrix
nonsemisimple endpoint checks. These finite checks do not replace the
all-rank proof or certify originality. Other highest weights, root systems,
degrees and integral images are outside this theorem's scope.

## How to read the claims

**Published input** means the cited source states the indicated result.
**Audited deduction** means the argument in these notes has been checked
against its hypotheses and published inputs in this audit. It is a
mathematical derivation being presented here, not a claim that a source
already states the exact Smith formula, or that the derivation is new.
**Conditional lemma** records its hypotheses explicitly.
**Conjecture** has no complete proof here.
**Finite verification** establishes only the exact recorded cases.

The previous wording did not make these distinctions sufficiently visible.
Absence of external review is not a mathematical counterexample; conversely,
agreement of computations or multiple assistant reviews is not a proof.
The conclusions below are based on the displayed arguments and checked
source statements. External expert review and any claim of priority remain
separate work.

## Classical equivariant API and symplectic higher Fitting ideals

The separate [22 October Cayley literature audit](CAYLEY_LITERATURE_AUDIT.md)
checks every bibliography entry in the supplied brief and the linked
notebook, with claim-specific primary-source locators. The common
trace map, quotient-degree comparison and regularity criterion are known
inputs. Agricola-Goodman's Theorem 4.1 already supplies a global gradient
basis of equivariant regular maps; their Theorem 5.3 gives an explicit
SL2 differential-operator harmonic complement. A general representation
classification preserving regularity was not located, but neither its
novelty nor the present openness of Popov's historical minimum-degree
problems is certified. The notebook's saved E8 inverse matrix has rank
seven and the file mixes E8/E7 settings. Part I's adjoint SL2 Example 3
has a proposed scalar correction from denominator two to four, checked
directly against its trace definition and KM Section 6.1. These are
specific source/reproducibility observations, not an invalidation of the
main source theorems. No external notebook code was executed.

The [classical backend](CLASSICAL_EQUIVARIANT.md) now evaluates defining
SL_m literal sources on every standard partial flag in equivariant
cohomology or K-theory. It uses published restriction formulas and
checks actual Laurent integrality. Full integral source generators
are enumerated only on G/B and P1; other parabolics support specified
invariant inputs. No complete type A image-module or augmentation-kernel
classification is claimed.

The [symplectic even/odd argument](SYMPLECTIC_K_FITTING.md) identifies
the entire defining Sp_(2n)/P1 integral cokernel, as a module over the
source coordinate z, with R[z]/(E(z^2+4)+zO(z^2+4),2O(z^2+4)).
Here product(s-s_i)=E(s^2)+sO(s^2). This gives an n-by-n presentation
in all ranks, a determinantal formula for every Fitting ideal, and
the closed last ideal (2) in odd rank and 2(e1,e3,...,e_(n-1))
in even rank. Explicit unimodular operations give global cyclic
decompositions for ranks one, two and three. The argument is
elementary algebra applied to the prior literal image; the ambient
K-theory, Fitting framework and source identification are inputs.
There is now a uniform stronger reduction: the cokernel is two copies
of an r-by-r block cokernel, r=floor(n/2), and an extra R/(2) if n is
odd. The block is multiplication by 2O modulo E for even rank, or
2E modulo O for odd rank. Monic division gives explicit integral
unimodular equivalences; parity and translation give the two copies.
This is an R-module splitting, not a source-coordinate module splitting.
It explains paired Smith factors under integer and DVR specializations.
For rank four, J=(e1,e3), F=product(s_i+s_j), the full Fitting list is
(16F^2), 8FJ, 4J^2, 2J, then R. Higher closed formulas were unresolved
at that stage; the subsequent Schur-minor result below advances them. The classical
Orlando resultant identity is a known input, not a claimed new result;
the specific image-module application has no priority determination.
No priority or external-review claim is made.

The software checks include BGG, Newton projective and line-bundle
restrictions, an explicit nonreduced-subword sum, and comparisons
against the original 2n projective inclusion at regular and singular
points through rank eight, plus rank-nine ordinary augmentation.
Independent coefficient variables give exact matrix certificates through
rank eight and all rank-four minor-ideal identities. The proof is
uniform; these finite checks are corroboration.

The subsequent [Schur-minor argument](SYMPLECTIC_K_MINORS.md) resolves a
uniform generator problem: for multiplication by q on A[u]/p with p
monic, the first j columns generate the entire ideal of j-minors.
This holds over arbitrary commutative A, by the published integral
Schur cyclicity of exterior powers (Laksov-Thorup; Gatto-Santiago).
It is an application of known Schubert algebra, with a full integral
argument supplied, rather than a claim of a new exterior-power theory.
The ideals are the content of a universal norm class in Schur coordinates.
For the defining symplectic image this gives all-rank partition-indexed
Fitting generators, plus complete explicit lists through ranks five,
six and seven. The cubic lists use a three-generator two-minor ideal
and adj(adj(H))=det(H)H to prove the required product containments.
There are 94 exact polynomial matrix identities and 72 independent
original-inclusion Fitting checks. Priority of the specific application is not established.

The subsequent [Hankel balancing proof](SYMPLECTIC_K_BALANCING.md)
settles the quartic simplification and proves the stronger arbitrary-rank
formula Fitt_i(C)=2^(n-i) I_floor(ell/2)(H) I_ceil(ell/2)(H),
ell=n-i-(n mod 2). A unimodular coefficient pairing makes every monic
multiplication matrix Hankel. Published generic Hankel product and monic
Groebner results (Le Dinh Nam, Theorem 3.25 and Corollary 3.26, c=1,
building on Conca) give rational containment; division by monic integer
generators proves the integral lift. This is not the false assertion
that generic unstructured matrix ideals balance in characteristic two.
There are 70 exact universal quartic Z identities, 68 original-inclusion
checks for the complete n=8,9 lists and polynomial-fiber tests in
characteristics two and three. The all-rank quantifier rests on the cited
theorems and the written lift, not on finite cases. Arithmetic primary
structure after moment substitution, minimal larger-block generators
remain; generic primary decompositions have not been carried through
arbitrary specialization. The ordinary multiplicity question from that
stage is now completed as follows.

The [closed ordinary Smith theorem](SYMPLECTIC_K_ORDINARY_SMITH.md)
identifies the augmented Hankel moments with Euler E750, VI--VII's
finite binomial continued fraction. S-to-J contraction yields recurrence
coefficients with only odd denominators. The published monic orthogonal
Hankel diagonalization (Miller-Stanton, Theorem 1) is therefore unimodular
over Z_(2), and gives every 2-primary factor, rather than only determinants.
The original triangular inclusion has power-of-two determinant, excluding
odd-primary factors. For n=2r+epsilon the paired exponents are
4+4epsilon+8j+s_2(r-j-1)-s_2(r+j+epsilon), 0<=j<r. There are n units
and one extra factor 2 for odd n. This completes all multiplicities for
the defining Sp_(2n)/P1 ordinary K-image, in all ranks. The written proof
is a deduction from classical inputs, not a claim of a new continued
fraction or Hankel theorem. Its priority is undetermined. Independent
integer half-block checks reach n=128, original inclusions n=64,
Euler-basis checks n=32, and explicit 2-local Gram identities n=32.
No equivariant primary decomposition or larger-flag formula follows.

## Other cohomology theories and affine groups: 22 October exploration

The [exploration report](OTHER_THEORIES_PROGRAM.md) records a new direction;
Spin parity is shelved at the user's request. The
[finite comparison](OTHER_THEORIES_FINITE.md) proves a common integral
K/Borel kernel and a strict support filtration for literal sources, then
constructs a specified oriented trace-coordinate lift. The universal MU
inclusion is triangular with the old integer determinant. This controls
prime support, while higher Fitting ideals and modular Tor kernels can
change. No gamma-filtration identification or arbitrary normalized-source
integrality is assumed. The [SL3 theorem](OTHER_THEORIES_SL3.md) gives
complete factors for every representation and parabolic from its source
scale and cubic correction; independent actual-character calculations
check 41 representations on all three proper flags.

The [type A coordinate comparison](OTHER_THEORIES_TYPE_A.md),
[symplectic height/K recurrence](OTHER_THEORIES_HEIGHT.md),
[quantum projective lifts](OTHER_THEORIES_QUANTUM.md) and
[affine A/C contents](OTHER_THEORIES_AFFINE.md) give written uniform
deductions with exact coefficient/source hypotheses. Formal group height,
the ambient quantum rings and Hopf models are known inputs. Independent
Euler, elliptic-curve and affine Chevalley calculations corroborate the
new image formulas. The [equivariant K note](OTHER_THEORIES_EQUIVARIANT_K.md)
uses deformed rational source freeness and reflection/localization
determinants for every parabolic; its local Smith patterns state their
regularity and DVR hypotheses. An integral Sp/P1 family is proved
separately, without promoting rational determinant formulas to general
integral source freeness. The [folding note](OTHER_THEORIES_FOLDING.md)
proves compatible trace projections under stable representations while
retaining exact integral source defects.

All are repository deductions and research prototypes. The cited ambient
characteristic-map, primitive K, oriented, quantum and affine theories are
published; the K-theoretic Peterson theorem is established. Exact
literature priority for these specified arithmetic/image formulations
has not been established. No new production runtime API, full exceptional
K Smith tables or central/rotation-compatible affine Springer construction
is asserted. The pointwise untwisted/twisted loop maps already follow by
base change; an explicit constant-loop example blocks naive flag descent.

## General central-parity kernels: positive families and counterexample

The [general kernel note](CENTRAL_PARITY_KERNEL.md) proves a sufficient
criterion using a torsion-free integral Borel presentation, a free actual
Levi source decomposition containing 1, and compatibility of the source
generators and relations with differentiation modulo two. This is a
sufficient criterion, not a necessary-and-sufficient root-data classification.

For the actual Cartesian weight lattice of Sp_(2n), symmetric-square
relations and ordinary/signed Artin bases prove the criterion for every
rank and every standard parabolic, in every degree. An odd-block type A
application works for SL_(2r) with arbitrary r, extending the scope of
kernel stability beyond the earlier odd-r quadratic-image theorem. The
kernel assertions alone do not extend the root-image Smith formula.
The positive proofs are elementary deductions from published integral
presentations, not novelty claims.

The naive index-two/fixed-odd-weight hypothesis fails on Spin7/B in degree
four. An explicit quartic is in the actual integral kernel, but its
derivative is outside the reduced primitive cubic kernel. A fresh
integer BGG calculation and an independent polynomial identity/monic
division proof establish the counterexample. Ambient mod-two cohomology
would miss it. No unit quadratic lift is available in this example.

The [odd spin extension](SPIN_PARITY_EXTENSION.md), audited 22 October
2026, disproves uniform failure on Spin(2n+1)/B. An Artin-basis base change
reduces the entire kernel in every rank to the auxiliary symmetric source
R[z,c2,...,cn], with c1=2z. Its primitive quartic relation proves that
degree-four failure occurs exactly at rank 3. A closed four-relation
presentation proves all-degree stability on Spin9/B. Independent rational
membership, modular Groebner and finite-flatness certificates prove complete
stability on Spin11/B and Spin13/B, and first failures on Spin15/B and
Spin17/B in polynomial degrees 11 and 16. These exhaust the stated spaces;
they do not classify arbitrary rank or assert novelty. The rational Borel
presentation, torsion-free Schubert basis and symmetric polynomial freeness
are known inputs; the derivative conclusions are deductions here.

The [recursive follow-up](SPIN_PARITY_RECURSION.md) proves an arbitrary-rank
image-lattice recurrence and exact derivative factorization criterion using
the classical Schur-Q root lattice and one-box Pieri rule. A lowering
operator satisfying [L,U]=(n(n+1)/2-2d) proves all-degree stability whenever
n=1,2 modulo 4; its preservation of the actual integral image lattice is
proved, rather than inferred from rational Lefschetz theory. An infinite
Schur-P derivation proves stability through degree n for every n. Native
Hermite/bit certificates check the recurrence through rank 12; a separate
quadratic straightening calculation independently gives stability in rank
11 and first failure at degree 34 in rank 12. The closed classification
for n=0,3 modulo 4 and a closed first-failure formula remain unresolved.
The Schur-P/Q presentation, Pieri rule and raising/lowering framework are
known; the central-parity conclusions are deductions, without a novelty
claim. No runtime API is changed by these proof and audit files.

This general question was the earlier structural priority and is now
shelved at the user's request. The existing all-degree E7 image theorem
is unchanged.

The [adjacent-rank boundary theorem](SPIN_PARITY_BOUNDARY.md) isolates the
whole obstruction without assuming that a derivative has descended in
earlier degrees. The boundary lattice is the entire kernel of
N_(n+1) -> N_n; lowering its largest part n+1 to n gives b_n. The derivative
obstruction is exactly (b_n L+N_n)/N_n, an elementary two-group, and its
dimension is rank over F2 of 2 F^(-1)b_n L. This is a proved general
criterion and an effective first-failure formula. A closed evaluation
depending only on n, particularly for n=0,3 modulo four, remains missing.
Classical Schur-P/Q identities supply the infinite derivation used in the
proof; they do not alone evaluate the boundary-image index. The new note
does not assert a literature-priority claim or an E7 all-degree structural
proof.

The same note proves a second, elementary source-parity formula. With
C=M[z^2], N=M[z], and e_d=log2[N_d:C_d], the exact degree-d obstruction
has dimension a_(d-1)-e_d-e_(d-1). Thus the entire kernel is preserved
if and only if log2[N:C]=2^(n-1). This does not evaluate that index in
closed form. New adjacent-rank certificates, audited by exact inclusion
and mod-two row rank without HNF software, agree with an independent
replay of the original derivative criterion. They establish the new
rank-15 first failure in polynomial degree 37 (Spin31/B, cohomology 74).
A separate Pfaffian calculation verifies 94 Schur-P lowering identities.

## Structural E7 multiplicities: scope of the new argument

The [quadratic-derivation note](E7_QUADRATIC_DERIVATION.md) isolates a
general kernel-stability criterion for the specified quadratic
presentation of N/2N over the reduced actual root image. The common E7
quadratic lift alone does not imply this criterion; an explicit algebra
counterexample is included.

Fresh primitive kernel columns prove the criterion through degree 10
uniformly for all 112 retained subsets with a fixed odd weight. The
integral relation expression is unique in that range, which justifies
passing to Levi invariants. The full-flag calculation separately reaches
degree 13, including saturated half-sums. Its independent audit uses
27 E6 forms rather than the primary Weyl orbit and its own exact bit
membership arithmetic. It reads no relative multiplicity table.

The same criterion gives an arbitrary-rank theorem for SL_(2r), r odd,
with at least one odd Levi block, over Z_(2): the usual integral Chern
presentation supplies full source/image and all-degree kernel stability.
Five complete Grassmannian images are independently checked by native
Schur/Pieri multiplication closure. This is a deduction from published
type A presentations, not a priority claim. It concerns 2-primary factors.

A subsequent [degree-eleven proof](E7_DEGREE11_LIFTING.md) handles the
first rank-one saturation correction uniformly. The mod-two cubic and
quintic are obtained by formal Frobenius derivatives of the common
quadratic. Explicit coefficient identities show that the half-sum's
derivative has an integral invariant lower-kernel lift. The argument
lifts only coefficients forced by an already invariant element; it does
not assume that every modular invariant lifts. A separate cubic
nonlifting obstruction shows that the half-sum coset is absent on E6,
the two D5 Levis and D4. Common-fiber descent covers the 108 good-Levi
cases and transfers the old full-flag derivative check through degree 13.
Thus the structural theorem now reaches degree 11 on all 112, and 13
on 108. Fresh 56-weight versus 27-form constructions independently
verify the rank-one kernel, coefficient identities, four obstructions
and family inventory. No relative Smith table is an input.

A uniform all-degree E7 kernel-stability and invariant-lifting proof is
still missing. Its all-degree multiplicity theorem retains the earlier
case/transfer proofs. Constructing the derivation from that theorem's
dimension count is valid as a corollary but would not independently
explain the counts. Runtime/API version remains 0.22.0.

## All-parabolic E7 API packaging

Version 0.22 also packages the existing complete E6 classification in the
[same relative API](E6_RELATIVE_IMAGE_API.md), preserving 3,9,27 separately.
The audit reconstructs all 63 cohomology/formula records, checks fresh
E6/P1 images in degrees 0,1,2,3,5,8 and rebuilds full A5 root/weight
invariants as integer reflection kernels in degree three, recovering Z/27.
This is implementation of the earlier image theorem, not a new E6 theorem.

The [G2 full-flag note](G2_FULL_FLAG.md) adds a complete finite image
calculation: Z[X,Y]/(X^2-3XY+3Y^2,Y^6), with the full literal minimal
source X=6omega1,Y=2omega2 transported through xi_(omega1)^B. The quotient
is free of rank twelve and all twelve basis images are independent;
this proves the integral kernel in every degree. A native root/coroot
Chevalley audit independently reconstructs all minimal/adjoint image
columns and all Smith factors. The earlier source theorem is shared.
All proper standard G2 parabolics are now covered for these literal
sources, and other irreducibles follow by the proved source-scaling law.
The [0.22 release record](release-0.22-verification.json) binds the new
checks. Literature priority remains unresolved.

Update, 21 October 2026. Version 0.21 packages the existing complete
relative E7 theorem in a [compressed summary API](E7_RELATIVE_IMAGE_API.md).
It includes all retained node sets, relative Smith multiplicities and the
normalized minimal/adjoint comparison from sources 12Q and 36X.
The gcd/lcm formulas for scaled image sums and intersections follow from
ordinary simultaneous Smith normal form. This adds no new image theorem
or literature-priority claim.

The [independent API check](relative-image-verification.json) reconstructs
the formulas and verifies scaled comparisons using integer matrix kernels;
fresh E7/P7 character-image overlaps cover degrees 0,1,2,3,6,10.
Previously audited full-source theorems, the A5 intermediary polynomial
and two native base image lists are shared inputs, explicitly identified.
The [0.21 release record](release-0.21-verification.json) binds these checks
to the runtime and installed wheel. Earlier version references below
describe their historical checkpoints. Absolute image cokernels and
structural generalizations are unchanged.

## First local-image benchmark and equivariant pilot

Update, 20 October 2026. The [new proof/source audit](LOCAL_IMAGE_AUDIT.md)
checks the [local image bound and projective family](LOCAL_IMAGE_BOUND.md)
and [equivariant specialization counterexample](EQUIVARIANT_IMAGE_PILOT.md).
The normalized source lattice lies in X(T), by differentiating its integral
Laurent pullbacks at 1. Primitive Schubert pullback then proves integrality
of every full literal Levi Borel image in the source theorem's scope.

Classical torsion-index divided differences supply integral Levi-invariant
lifts of t_G times every parabolic Schubert class. Hence each degree-d
full-source cokernel is killed by t_G e_rho^d, with
e_rho=exp(X/Lambda_rho). There is no extra Levi Weyl-order factor.
This is an explicit application of published characteristic-map machinery,
not a new general good-prime theorem.

The standard free equivariant Schubert basis and triangular induction
also recover the whole equivariant image at these good primes. Each
fixed-degree specialization kernel is finite with primes dividing
t_G e_rho. This does not force faithful specialization at the bad primes.

For defining SL_n on P^(n-1), the literal factors are n/gcd(n,d).
The proof uses elementary root-difference invariants, binomial generator
images and a p-adic gcd argument. Its root-characteristic interpretation
requires further literature comparison; priority is not established.
The equivariant SL4/P1 image gives a nonzero specialization kernel:
exactly Z/2 in polynomial degree five, with no defects in smaller degrees.
The integral proof and independent Sage/SymPy rank-52 lattice calculations
agree. Higher-degree defects and a general equivariant backend remain open.
Runtime sources and the 0.20.0 wheel are unchanged.

## Sharp E6 root-to-weight comparison

Update, 20 October 2026. The [comparison theorem](E6_ROOT_WEIGHT_COMPARISON.md)
and [stepwise audit](E6_ROOT_WEIGHT_AUDIT.md), together with the
[three-case completion](E6_REMAINING_CASES.md) and its
[audit](E6_REMAINING_AUDIT.md), prove exact 3-primary factors for all 63
proper standard E6 parabolics. The Cayley-plane image
has a complete local source-monomial basis. A common relative Schubert
basis transfers its factors to finer flags; at a good adjoint Levi,
Leray--Hirsch deconvolves full-flag multiplicities by the Levi length
polynomial, covering 60 cases. The three retained sets containing
{1,3,5,6} use a different proof: an A5+A1 subsystem at 3 and full
source coefficient lattices. The two base cases have 52 independently
audited degrees in total; the third follows from a common (1,z) basis.
The [complete tables](E6_REMAINING_TABLES.md) include all three cases.

The ring presentation and characteristic-map framework are published
inputs. Totaro's adjoint torsion index already implies the uniform
exponent bound 27 and the full-flag top-degree ratio 9. The exact
all-degree multiplicity formulas and their all-parabolic application
are deductions recorded here; no literature-priority claim is made.
An independent Python/SymPy implementation reconstructs the base
lattices, BGG coefficients, Levi Weyl polynomials and deconvolutions;
additional direct BGG checks test other flags. The completion independently
reconstructs the rational quotient and every source/image column, and has
20 further direct BGG degree/case checks. Sage exact arithmetic is shared
in this replay; the structural and literal-weight checks use plain
Python/SymPy. These are computer-assisted deductions, not an ambient H/N
classification. Runtime remains 0.20.0.

## E7/P7 simultaneous bases and the refinement exponent

Update, 20 October 2026. The [2-primary comparison](E7_ROOT_WEIGHT_COMPARISON.md)
and [audit](E7_ROOT_WEIGHT_AUDIT.md) identify the full local images as
N=R[h,2s,2t] and M=R[2h,h2,2s,2t], with R=Z_(2). A common 56-element
basis consists of h^i for i<=27 and h^i(2s),h^i(2t) for i<=13;
the root image doubles exactly the elements with odd h exponent.
Thus N^d/M^d is elementary 2-primary, with multiplicity polynomial
q[14]_(q2)+(q6+q10)[7]_(q2). The old finite tables already recorded
these factors; the extension is their explicit basis proof.

The candidate-ring basis and its closure follow by explicit identities
in Duan--Zhao's published integral ring. Four root-invariant lifts give
the candidate generators. The full-image identification is a finite
computer-assisted step: both complete invariant sources are freshly
certified in all 28 degrees, retaining the rational kernel before
integrality and image formation. An independent Girard/SymPy/bit-elimination
replay checks every source and image column; Sage exact arithmetic is
shared in that replay. Fresh BGG checks verify all seven generators,
and Chevalley multiplication checks all 56 basis classes. A plain
Python/SymPy check separately verifies lattice and ring identities.

An index-3 common E6 complement at 2, with omega7 fixed and a root
quadratic lifting h2, gives N_I=M_I+hM_I and hence 2N_I contained in M_I
for all 64 retained subsets I of {1,...,6}. This is a structural
adjoining-generator deduction, not an application of a non-torsion
Levi-fiber theorem. It supplies no multiplicities on those finer flags.
Totaro's general exponent bound 8 is published background; priority
of the sharper application and explicit basis statement remains
unestablished. Runtime and the installed wheel remain 0.20.0.

## E7/P7 refinements: a uniform truncation and a source obstruction

The [refinement note](E7_REFINEMENT_MULTIPLICITIES.md) proves the multiplicity
formula [q^(d-1)] (1-q6)(1-q8)H_I(q) for every retained I in {1,...,6}
and 0<=d<=9. The 56-weight and independent 27-form constructions yield
primitive quadratic/sextic/octic kernel certificates in this range.
A direct full-source BGG/Hermite/Smith calculation on the A5 refinement
checks all ten degrees. This is a uniform bounded-degree theorem, not
a complete finer-flag multiplicity formula.

The original all-degree target reduced by classical good-Levi deconvolution
and an integral D5 symmetry to three base polynomials: full flag, D5 and
D4. The subsequent A6 theorem below completes the full-flag polynomial. A separate degree-16 Chern invariant proves that six normalized
rational-degree E6 generators do not generate the full invariant source
at 2. Both coefficient completeness and the witness have a plain-integer
audit. This does not disprove an image transfer formula. Known inputs,
precise finite scope and unresolved novelty are recorded in the note.

## E7 full flag: complete A6 intermediary proof

The [new all-degree comparison](E7_FULL_FLAG_COMPARISON.md) proves
T_full=q P_(W(E7))/(1+q). A common A6 complement at 2 makes both full
invariant sources polynomial, and 576-dimensional exact image bases
exhaust every degree of E7/P2. The independent checks include fresh
126-root/Girard relations, plain-Python Buchberger and image closure with
primitive generation, and direct full-flag saturated kernels in degrees
10..13. The latter exhibit corrections missed by unsaturated rational
relation multiples. Both complete A6 source identifications and rational
kernel completeness are justified before taking images.

Classical common Artin/fiber bases transfer the result to the full flag
and 108 proper good-Levi parabolics; the separate P7 theorem adds a 109th.
This A6 argument completes 61/64 P7 refinements; the D5 completion below
adds two, leaving D4 as the remaining all-degree target. The degree-16 E6
invariant-ring obstruction is still
valid. The new proof uses an A6 intermediary. Published inputs, shared
arithmetic backends, finite scope and unresolved novelty are stated in
the proof note. Ambient H/N and arbitrary direct sums are separate.

## E7 D5 refinements: complete full-source comparison

The [D5 proof and audit](E7_D5_REFINEMENT.md) proves
T_D5=q P_(W(E7))/((1+q)P_(W(D5))) for retained sets {2,3,4,5,6} and
{1,2,3,4,5}, with 756 factors 2 across degrees 0..43. A ten-generator
2-local invariant source is proved complete by a mod-2 graph/Hilbert
certificate and Nakayama, rather than assuming polynomiality at 2.
The rational quotient of total dimension 1512 and both complete image
lattices are exhausted. Plain-Python source, fresh-weight polynomial
identity/Jacobian, Buchberger and image-generation audits use separate
arithmetic implementations. E6 opposition transfers the result to the
second D5 case. Classical inputs and unresolved priority are identified
in the proof. The method now covers 111 proper E7 parabolics and 63/64
P7 refinements at this stage; the D4 theorem below completes the family. This does not classify
ambient H/N or direct-sum sources and adds no runtime API.

## E7 D4 completion: all 64 P7 refinements

The [D4 proof and audit](E7_D4_REFINEMENT.md) completes retained
I={2,3,4,5}: the relative factors are all 2, with multiplicity polynomial
q P_W(E7)/((1+q)P_W(D4)) and 7,560 factors across degrees 0..51.
Its ten-generator full integral source is certified by an independent
mod-2 graph/Hilbert and lattice audit. The relative source quotient over
D5 has dimension 34, so a free rank-10 invariant-source transfer fails.
Instead all 34 integral module generators enter a rational quadric tower.
The reference root and weight lattices enlarge by equal indices in every
degree, checked by PARI HNF and independently reconstructed matrices with
native unit-pivot Smith elimination. The previously audited D5 rational
quotient and full image bases are shared inputs, explicitly identified.

Together with the P7/A6/D5 proofs, this proves
T_I=q P_W(E7)/((1+q)P_WI) for all 64 I contained in {1,...,6}.
The common source split and these Hilbert dimensions also imply that
N/2N is a free quadratic algebra on (1,w) over C=image(M -> N/2N),
with w^2=beta(B_i), choosing any omitted i in {2,5,7} and
w=beta(omega_i). The three common complements have odd indices 7,15,3,
and the integral quadratic lifts are checked separately. This C is not M/2M.
The combined scope is exactly the 112 proper E7 parabolics omitting at
least one of nodes 2,5,7. The other 15 retain all three and have N^1=M^1,
so the same polynomial cannot extend unchanged to them. This is an exhaustive computer-assisted proof for the stated
finite family, not a uniform theorem for arbitrary root data. Known
Borel/invariant/Nakayama inputs and the new finite deductions are
separated; priority, ambient H/N and a structural explanation of the
index equality remain unresolved. No runtime API changes.

## E7 completion: all 127 proper standard parabolics

The [complete E7 theorem](E7_ALL_PARABOLICS.md) determines the relative
full root/weight Borel-image groups in every degree on all 127 proper
standard parabolics. Every nonunit relative factor is two. The preceding
112-case polynomial q P_W/((1+q)P_WI) applies exactly to the retained
sets omitting at least one of 2,5,7; the other fifteen have the separate
formulas tabulated in the theorem. Their first gaps occur in degree two
(twelve cases), four (two cases), or eight (D6).

The [remaining-family note](E7_REMAINING_PARABOLICS.md) proves for E7/P1
(D6 Levi) that every nonunit relative root/weight image factor is 2, with
multiplicity polynomial q8(1+q)([10]_(q2)+q4[6]_(q2)+q12), totaling 34.
Full weight-source components are certified through top degree 33, root
components through 27; in degrees 28..33 an explicit root lower image
equals the full weight image. No unbounded invariant-algebra generation
claim is used. Native Python/FLINT source reconstruction, fresh 56-weight
moment/Jacobian identities, and native Buchberger/image audits verify
the finite proof. Higher root invariants change the truncated-generator
answer, so selected rational/integral generators alone are insufficient.

Four Weyl-conjugacy words are checked independently in root-mask and
root-matrix implementations. They reduce the initially excluded fifteen
diagrams to eleven types. The [D5+A1 completion](E7_D5_A1_COMPARISON.md)
now proves the complete E7/P6 polynomial, totaling 366 factors two.
Its full 25-generator root source uses a necessary degree-23 correction
and the degree-32 root-orbit norm. The correction's Borel image is
redundant, so the earlier independently audited lower image is the full
image. The fresh 56-weight orbit, exact Jacobian and 333 unisolvent
samples independently establish the base relations. The earlier D5
source theorem is an explicitly shared input.

The [all-degree source argument](E7_LOCAL_SOURCES.md) proves both this
25-generator source and the 18-generator D4+A1 root source. It compares
an explicitly witnessed leading-monomial algebra with the rational Weyl
invariant Hilbert series as a full rational-function identity. Independent
audits verify the product witnesses, nonnegative monoid expressions,
closure and monomial-ideal Hilbert recursion. It therefore proves
generation in all degrees, rather than extrapolating from a cutoff.

The [A5-family proof](E7_A5_FAMILY.md) completes twelve diagrams. An
odd-index A5+A2 intermediary has 293 factors two; common integral A2
source bases give two transfers. A uniform involution lemma supplies
the full even-block root sources, with 51- and 342-generator modules.
Separate rational reducers and mod-four row/column eliminations agree
on every image column through top degrees 53 and 57. Their complete
formulas have 4,911 and 29,982 factors and a common Gaussian-multinomial
expression with a halved-block correction. The expression is an exact
finite E7 deduction; a geometric folding explanation or an arbitrary-rank
generalization is not proved.

The [D4+A1 completion](E7_D4_A1_COMPARISON.md) uses the rank-ten quadric
tower over P6 and the full 34-generator weight-source module. A native
FLINT audit reconstructs all source images and checks full local lattices
in every degree through 50. It gives 3,618 factors two. The supplementary
finer-flag intersection bound can be strict; invariants after Borel are
not used as a full-source substitute.

The rational cohomology of E7/P1 is already known from Chaput--Perrin,
Proposition 5.6; the Borel presentation, signed-permutation invariants,
Weyl conjugacy and lattice primitivity criteria are classical inputs.
The relative Smith calculations and source-completeness deductions are
what is established here. Earlier targeted literature checks are retained
in the linked notes; the completion does not assert a new comprehensive
search or establish priority. Structural explanations beyond this finite
classification and ambient H/N remain separate. No runtime or public API
changes were needed for this completion.

## Exceptional integral images: first complete finite tables

The [Springer geometry extension](SPRINGER_GEOMETRY.md) concerns a separate
characteristic-zero problem. It records complete G2 minimal/adjoint quotient
maps, generic degrees 6 and 18, and the full left-translated critical
determinants. The degrees have Newton-polygon proofs and exact generic-fiber
checks; the full determinants agree coefficientwise with independent rational
representation-matrix constructions. A general root-factor formula follows
from equivariance. A standard tensor-Casimir identity gives the quadratic
coordinate uniformly, with E6/E7 minimal/adjoint checks by both Laurent and
tensor methods. The [E6 extension](E6_GEOMETRY.md) now supplies the full
basic-invariant list, all six minimal/adjoint quotient coordinates and both
full determinant products. Dominance bounds and exact unisolvent matrices
prove the large character formulas; independent 27- and 78-dimensional
matrices check all 72 root directions coefficientwise in both cases.
The [E7 extension](E7_GEOMETRY.md) supplies all seven coordinates through
degree 18 and both full 133-dimensional determinant products. Its dominance
bounds and exact unisolvent checks certify all coefficients; independent
56/133-dimensional representations check all Cartan entries and 126 root
directions for both maps. The 912-dimensional character is reconstructed
separately by root strings and Freudenthal. The [generic-degree deduction](EXCEPTIONAL_DEGREES.md)
now gives E6 degrees 1080/2754 and E7 degrees 7272/33664 for minimal/adjoint
maps on simply connected groups. It uses explicit toric boundary schemes
and their localized intersection contributions. Separate specialized-fiber
lengths provide finite corroboration, not the characteristic-zero proof.
These are exact deductions and computations, without a
novelty claim. They do not identify the minimal Cayley degree of a group.

The [irreducible lattice classification](IRREDUCIBLE_CHARACTER_LATTICES.md)
proves D_lambda=Q+Z*lambda by multiplicity-one simple-root paths, hence
Lambda_lambda=I_lambda*(Q+Z*lambda)^*. This identifies every exceptional
irreducible torus source: index times Q_long for G2/F4; index times P or Q
for E6/E7 according to whether the central character is trivial; index
times P=Q for E8. The index normalization is the published Dynkin
convention I_ad=2h^vee. The highest-weight and trace formulas are published
inputs; their source-lattice interpretation is the deduction recorded
here, without a priority claim. The fast constructor enumerates only
roots, not representation weights. A reducible F4 control shows that
index and central character alone do not classify direct sums. This
does not determine arbitrary integral Levi presentations or Borel images.
The [detailed classification audit](IRREDUCIBLE_LATTICE_AUDIT.md) found
no gap and separates the known ordinary weight-span result from the
multiplicity-weighted argument and integral Springer interpretation.
The exact formulation was not located in the reviewed literature;
priority remains unresolved. Its independent boundary checks test every
nonempty support mask on all five exceptional diagrams and explicit
controls for path selection, dual normalization and reducible inputs.

The [general trace-source theorem](TRACE_SOURCE_THEOREM.md) now proves
C_rho=Sym_Z(Lambda_rho) for every supported normalized representation map,
including direct sums with faithful differential. Saturation gives full
coefficient rank modulo every prime; the trace-moment factorization makes
the adapted pivot exponents independent over Q. Leading terms and the
least-denominator argument then prove the full, inhomogeneous intersection.
This resolves the previous representation-by-representation degree-one
generation question under these hypotheses. It does not assert that
Levi-invariant rings are generated in degree one. The saved survey and
Python-only weight/trace/minor audit check 59 representation cases,
354 prime-adapted leading bases and 18 bounded searches. The proof is
all-degree and all-prime; the numerical checks are finite. No originality
claim is made. The [2026-10-20 detailed proof/source audit](TRACE_SOURCE_AUDIT.md)
found no gap under the stated hypotheses and distinguishes coefficient
minors modulo p from exponent independence over Q. It did not locate the
exact integral equality in the reviewed literature; priority remains
unresolved.

The [E7 adjoint proof](E7_REPRESENTATION_COMPARISON.md) identifies
Sym_Z(36P) by a unit coefficient minor and a determinant-one row change
whose seven leading coefficients are all one and whose exponent
determinant is two. The certificate works in every characteristic.
Since P/Q has order 2, the coprime-index argument proves
J56^d intersect J133^d = 3^d J56^d and their sum = 3^(-d) J133^d,
for every standard Levi and degree. These are deductions without a
novelty claim. The divisor images 24 Z h and 36 Z h prove that the
full E7/P7 image rings are incomparable; some individual degrees have
adjoint-in-minimal containment. Specifically, containment is strict in
degrees 2, 4, 24 and 26; all other positive degrees through 27 are
incomparable. The quotient of their sum by the minimal image is
elementary abelian of rank zero, one or two, as listed in the comparison.

The complete E7/P7 comparison computes all 28 adjoint source lattices,
certifies their completeness, and freshly checks all 1,967 rational
Borel columns by the separate ordinary-ring/Chevalley implementation
with its documented BGG anchor dependencies. A Python-only audit
reconstructs the integral weight-lattice coordinates, compressed
polynomials by a closed partition formula, every source coefficient
constraint, eight direct Laurent intersections and image arithmetic.
The prior minimal tables are reused with six fresh minimal API checks.

The [E6 adjoint proof](E6_REPRESENTATION_COMPARISON.md) identifies the
full torus source as Sym_Z(24P), while the minimal source is Sym_Z(6Q).
A unit coefficient minor proves the linear lattice, and a unimodular
change with six independent unit-leading monomials proves the full source
in every characteristic. Since P/Q has order 3, a coprime-index argument
proves J27^d intersect J78^d = 4^d J27^d and their sum = 4^(-d) J78^d
for every standard Levi. These are checked deductions, with no novelty claim.
Complete E6/P1 tables show incomparability in all 16 positive degrees;
the divisor example alone proves incomparability of the full image rings.

The [E6/E7 public APIs](MINUSCULE_CHARACTER_API.md) reproduce all 45
historically audited minimal degree rows and compute all 17 adjoint rows.
All 2,731 rational image columns are checked by the separate
published-ring/Chevalley implementation, retaining its documented BGG
anchor dependencies. A standard-library audit reconstructs every adjoint
source coefficient lattice and eight direct Laurent certificates, and
checks all image/comparison arithmetic. The finite checks are distinguished
from the all-degree proof. The API is limited to E6/P1 and E7/P7; it does
not classify arbitrary representations or supply finite integral source
presentations.

The [general character-source criterion](CHARACTER_SOURCES.md) proves that
the entire integral torus source equals Sym_Z(Lambda_rho) exactly when
the pullbacks of a basis of its integral linear lattice are algebraically
independent modulo every prime. This is an elementary saturation argument,
not a claimed new theorem. The finite leading-term certificate is sufficient;
an inconclusive certificate or an unsuccessful bounded relation search is
not a contrary all-degree result.

The [G2 application](G2_CHARACTER_FORM.md) proves C7=Z[6omega1,2omega2]
and the two maximal-Levi invariant presentations over Z. It gives the
complete integral image presentations Z[X,Q]/(X^2-3Q,Q^3) and Z[Y]/(Y^6).
The source proof is all-degree; the image kernels follow from free quotient
bases with nonzero images in different degrees. The adjoint source is
Z[24omega1,8omega2], with a separate characteristic-three certificate.
A Python audit independent of Sage and the runtime library checks the
torus certificates, Laurent intersections and every image comparison.
These are checked deductions, with no claim of originality. The abstract
modular-obstruction example is not asserted to come from a representation.

The [F4 character-form proof](F4_CHARACTER_FORM.md) identifies the entire
26-dimensional torus source as Sym_Z(6 Q_long), hence the full integral
character source for every Levi by taking invariants. A unit-coefficient
leading-monomial certificate handles all primes simultaneously.
All four maximal-parabolic images are computed in every degree (74 rows).
Source lattices have integral-coefficient and gcd-one maximal-minor
certificates; every Borel column is checked by localization and independent
root/coroot Chevalley multiplication. A standard-library Python audit
reconstructs all source lattices and checks their certificates and image
arithmetic. This does not supply a finite integral invariant-ring
presentation or a classification of representations, and no novelty is claimed.

The [F4 adjoint comparison](F4_REPRESENTATION_COMPARISON.md) proves
C52=Sym_Z(18 Q_long), using separate characteristic-two and odd-prime
leading-term cases. Taking invariants and the common Borel map proves
A52,L^d=3^d A26,L^d and J52,L^d=3^d J26,L^d for every standard Levi and
every degree. This is a checked deduction, with no priority claim.
All 74 adjoint image rows are derived from this theorem and the audited
minimal tables. Fresh API checks cover 20 parabolic-degree pairs; 32
homogeneous/filtered Laurent certificates are independently reconstructed.
These finite checks support the implementation; they are not the
all-degree proof, and the two Springer morphisms are not scalar multiples.

The [G2/F4 note](EXCEPTIONAL_INTEGRAL.md) studies all maximal parabolics and
all degrees for two explicitly selected generator forms: least integer
image-denominator clearing, and individual primitive-image rescaling.
The complex generators are published Part III inputs. The G2 image
presentations follow by the displayed Chevalley calculation. F4 Smith lists
are exact finite computations, independently reconstructed by root/coroot
Chevalley multiplication; no general exceptional classification is claimed.
An additional F4 root-duality argument identifies the original selected image
algebras in the pairs P1/P4 and P2/P3, and identifies their odd-primary ambient
quotients degreewise. It depends on the displayed, independently checked
generator identities; it does not identify the ambient integral lattices or
their 2-primary quotients. The primitive-image forms have the same odd-primary
quotients as the original forms because their rescaling uses powers of 2.
Neither chosen form is identified with a representation-ring intersection,
and novelty has not been established. In particular, primitive individual
generator images do not guarantee saturation of the generated algebra.

The [E6/P1 and E7/P7 note](EXCEPTIONAL_MINUSCULE.md) adds complete tables
for both forms in all degrees. The known integral presentations are
Duan--Zhao Theorems 4 and 5 (with the E6 diagram involution). Fixed-point
interpolation checks every monomial against an ordinary-ring engine using
those presentations and Chevalley multiplication. Its starting generator
images and one E6 square anchor use the original BGG implementation; these
dependencies are explicit. The 90 degree/form rows contain 4744 columns
obtained from 2372 independently compared unscaled monomials and exact
rescaling. Integer identities in both directions prove that the primitive
images are Z[h,16s] and Z[h,81s,19683t]. These are deductions for the specified
forms inside published cohomology rings, not a general exceptional theorem
or a claim of novelty. Their quotients can involve only 2 and 3 respectively.

The subsequent [E6 character-form proof](E6_CHARACTER_FORM.md) identifies
the literal intersection with integral Levi characters, for simply
connected E6 and the normalized omega1 Springer map, as Sym_Z(6Q)^(W_L).
Here Q is the E6 root lattice. A unimodular change of trace-moment
coordinates has independent leading Laurent exponents over every prime
field, with characteristic 3 checked separately. Denominator reduction
then proves the full torus intersection, and equivariance gives the Levi
form. This establishes its grading in all degrees, not just a bounded
test. The argument is a repository deduction without a novelty claim.

For P1, all 17 cohomology degrees are computed using complete invariant
source lattices, with independent gcd-one minor certificates for their
completeness. Fixed-point and published-ring calculations agree. A separate
standard-library Python audit verifies the trace/leading-term argument,
ten low-degree Laurent intersections, all 51 comparison lattices and all
34 containments. The literal character image is contained in both earlier
chosen images; it already differs in degree 1. This resolves the E6/P1
source-form question.

The [E7/P7 character-form proof](E7_CHARACTER_FORM.md) subsequently proves
the analogous identification Sym_Z(12Q)^(W_L), with Q now the E7 root
lattice and the normalized 56-dimensional representation. Its all-primes
leading-term argument treats characteristic 2 separately. An integral
permutation basis for an A5 subgroup compresses coefficient tests without
changing the lattice. Complete invariant source lattices are certified by
gcd-one maximal minors, and Borel products are compared by the two earlier
independent methods. A standard-library audit separately checks the torus
certificate, permutation coordinates, low-degree Laurent intersections
and every embedded image comparison. The character and catalogue images
are incomparable: (216/5)psi5 is a character-form element with Borel image
outside the catalogue, while h is absent from the character image. This
is a deduction with explicit certificates, not a claim of novelty or of
a minimal polynomial presentation of the integral invariant ring.

## General type-C joint images

The [joint-image proof](TYPE_C_JOINT_PRIMITIVITY.md) proves, for every
a,b,k>=1 and the fixed normalized inputs c_j=e_j(M), q_i=e_i(T^2), that:

- the joint image is primitive in integral symplectic flag cohomology;
- the sum of the two separate images is primitive in the joint image
  and in ambient cohomology;
- their intersection consists exactly of the constants.

The proof identifies the joint image with pullback from
Sp(n)/(Sp(a) x U(b) x Sp(k)). Integral Lagrangian Schubert lifts give
Leray–Hirsch decompositions. Reordering the quaternionic blocks as M,A,T
identifies the two squared-factor subrings with standard Grassmannian
pullbacks in a complex type A two-step flag ring. They occupy subsets of
one integral Schubert basis. The same middle Lagrangian lifts then give
a simultaneous basis for the actual separate images. The proof includes
an integral kernel presentation and Gaussian rank formulas.

Borel's integral flag presentations, the Pragacz–Ratajski Lagrangian
basis and quaternionic description, and compatibility of Schubert classes
with projections are published inputs. Their application to these selected
joint Springer/Borel images is a deduction written and checked here.
The note's audit table records the normalization, integral generation,
compatible lifts, and block-order issues. Its
[finite verifier](experiments/type_c_joint_primitivity.py) constructs the
asserted bases and checks their embedded lattices. No claim of novelty,
external review, or extension to scaled characters or type B/D joint images
is made.

A subsequent [focused proof and literature audit](TYPE_C_JOINT_AUDIT.md)
found no mathematical gap under these hypotheses. It classifies the result
as a corollary of classical integral theory applied to the selected
Springer/Borel generators, rather than an established new theorem.
Borel's Proposition 30.2 independently supplies the intermediate ring.
The exact combined statement was not located in the searched literature;
no priority conclusion follows. The audit clarifies the global fiber lifts,
compact homogeneous-space maps, integral elimination and necessary
normalization, without changing the formulas. A separate formal-polynomial
check tests the integral base-change basis without importing the library.

## Orthogonal two-step joint images

The [B/D joint-image argument](ORTHOGONAL_JOINT_IMAGES.md) proves that
the sum of the normalized middle and tail images is primitive inside
their joint image, and that their intersection is the constants. It
assumes a,b>=1, k>=1 in B and k>=2 in D. The joint image need not be
primitive in ambient cohomology, but its nonzero Smith factors are powers
of 2. The note gives integral kernels and closed ranks, including the
free rank of the joint image modulo the separate-image sum.

In B, equality of the rational Cartesian kernels and torsion-free target
cohomology identify the integral source quotient with the type C one.
In D, the global Euler relation gives the additional joint relation
c_b p h_a(M^2 union T^2). A four-component decomposition into free
integral flag-ring modules proves the integral kernel and locates both
separate images exactly. These arguments do not use the earlier B/D
single-tail Smith multiplicity theorems.

These are repository deductions from classical inputs, without a novelty
claim or external certification. The [finite verification](experiments/orthogonal_joint_structure_verification.json)
checks exact integral kernels, ranks, image intersections and quotient
torsion, with selected independent Cartesian BGG columns. The remaining
problem is a closed formula for the ambient joint Smith multiplicities,
not the primitive-sum assertion.

## General type-B comparison completed in these notes

The subsequent [general proof](TYPE_B_GENERAL_SMITH.md) supplies the
formerly missing comparison in every rank. A balanced incidence map is a
Chow envelope over function fields. It imports unbalanced-incidence rows
integrally; their Jacobi--Trudi leading terms give a unit triangular
minor. A beta-number count matches the defect-locus Betti bound, forcing
all prefixes primitive with the required ranks. This is a complete
repository derivation under the stated normalization and earlier prefix
reduction, not a published or externally peer-reviewed theorem.

A [fresh full-proof audit](TYPE_B_GENERAL_AUDIT.md) rechecked the complete
dependency chain at baseline `5e71833`, including the earlier prefix
reduction. No gap or mathematical correction was found. The
[updated Shrawan PDF](../articles/orthogonal-images-summary/README.md)
now includes the general result and its integral argument, with the
remaining type-D questions still distinguished.
The fresh audit also preserves an independent pure-Python larger-case
check: 16,200,912 corner comparisons and 88 degreewise polynomial
identities, with exact scope and hashes in its
[report](experiments/type_b_general_audit_verification.json).

The polynomial expression also has an existing Demazure-character
interpretation: Chari--Schneider--Shereen--Wand, Theorem 3.3(ii)--(iii),
[primary text](https://arxiv.org/pdf/1310.5191). That complex character
identity does not itself prove the integral Springer comparison.

## Type-D progress after the general B proof

The [type-D unit theorem](TYPE_D_FROBENIUS.md) determines every normalized
unit multiplicity, including the exact Gaussian overlap for even tails.
Its mod-2 argument retains its independent scope.

The [general square comparison](TYPE_D_SQUARE_THEOREM.md) now proves
that every normalized D square-sector coefficient-row lattice equals
its B counterpart in every rank and degree. In particular the requested
smaller-block-three case has exactly the multiplicities of \(1,2,4,8\)
displayed in equation (13) of that note. The proof uses flagged orthogonal
incidence bundles, the published equivariant eta formulas with denominator
at most \(2^t\), and an integral coefficient retraction obtained by
compensated head/tail sign changes. The retraction preserves the nonzero
base Euler class. The resulting square ring is the relative quaternionic
Grassmann ring. Existing integral B-prefix primitivity then gives
\(2^tP_t\subseteq\mathcal R_D\); summing proves the missing reverse
row containment. A fresh independent adversarial audit found no material
gap in this argument. This is an unpublished repository deduction without
an originality claim or external expert certification.

The [two-row note](TYPE_D_TWO_ROW_SMITH.md) retains an independent direct
boundary-row proof and the general exponent bound. The
[structural comparison](TYPE_D_SQUARE_COMPARISON.md) retains common
kernels, duality, finite primitive ranges, and the failed broader shortcuts.

The [general even-tail gluing theorem](TYPE_D_EVEN_GLUING_THEOREM.md)
now proves the [full-D formula](TYPE_D_CANDIDATE.md) for every tail,
closing the remaining saturated-gluing condition. Moving the specified
tail to the head gives a precise square-plus-head-Euler source. Its
relative Euler extension has a primitive selected source, including
the nonzero base Euler relation. The eta denominator bound and an
even-head parity retraction give exact D/B row-lattice equality.
Primitive actual-head B prefixes and signed level-zero overlap
identities compute that B lattice. Comparing the resulting full index
with the two sector indices proves that every gluing direction is
already a unit direction.

The [fresh integrated audit](TYPE_D_EVEN_GLUING_AUDIT.md) records the
complete dependency checks, including large actual heads and the
two-fork boundary. A misleading intermediate claim that wider relative
Schur classes vanish was corrected: only the narrower source subset
needs to map basis by basis. No unresolved gap was found in the final
argument. These are internally audited repository deductions, not
published or externally certified results.

The independent [two-head](TYPE_D_TWO_HEAD_GLUING.md) and
[two-tail](TYPE_D_TWO_TAIL_GLUING.md) proofs, and the
[full exponent bound](TYPE_D_FULL_EXPONENT_BOUND.md), are retained.
The stronger unweighted prefix strictness conditions in
[TYPE_D_GLUING_FILTRATION.md](TYPE_D_GLUING_FILTRATION.md) remain
a separate conditional route; the new proof does not assume them.

The [new full-gluing check](experiments/type_d_even_gluing_verification.json)
passes 138 full components, 16 independent full row-module comparisons,
and 27 explicit primitive unit-gluing representatives. Both historical
rank-nine timeouts now complete. The separate actual-head verifier
checks 610 exact prefix modules and 35 signed overlap identities.
These are finite corroboration, separate from the arbitrary-rank proof.

The [three-root Chern verifier](experiments/verify_type_d_three_row_smith.py)
and [report](experiments/type_d_three_row_smith_verification.json) check
actual integer coefficient modules, the four closed formulas, independent
Cartesian D expansions, and a relative quadric with nonzero base Euler
class. The separate pure-Python
[projection verifier](experiments/verify_type_d_incidence_projection.py)
and [report](experiments/type_d_incidence_projection_verification.json)
check 5,440 monomials and their Euler/square-relation multiples, retaining
the base Euler symbol and detecting the false uncompensated projection.
Finite checks corroborate the written arguments; they do not prove the
geometric or equivariant inputs. Earlier unit and square reports retain
their explicitly bounded scopes and resource exclusions.

## Current mathematical status

Throughout, $q_i=e_i(x_{\rm tail}^2)$ are the fixed normalized inputs,
$p=\prod x_{\rm tail}$ is included in D, and polynomial degree $d$
means cohomological degree $2d$. Indices are in the saturation of the
selected rational image. They need not be indices in all cohomology.

| Claim | Status after audit | Reason or limitation |
| --- | --- | --- |
| Integral theta/eta bases, finite presentations, B/C Chern-preserving comparison | Published inputs | BKT sources below; the B-to-C powers of two are existing Schubert calculus. |
| Quaternionic Grassmannian integral Schur basis and isotropic-to-quaternionic fibration | Published inputs | Pragacz--Ratajski, Appendix A Sections 10.3, 10.5; bundle version also in Panin--Walter. |
| The polynomial $U_{n,\ell}$ as a count of square-zero invariant subspaces | Classical counting, not a new polynomial | Exact specialization of Birkhoff/Hall-type counting; see below. |
| Counting subspaces with prescribed defect $\dim(U+NU)-\dim U$ | Existing general counting framework; elementary specialization here | For $N^2=0$, defect $e$ means subspace profile $(\ell,e)$. Ram's general profile-count theorem applies over finite fields; it is not a Springer-Smith theorem. |
| All-rank type-B unit Smith multiplicities equal $U_{n,\ell}$; primitivity of $P_0$ | Audited deduction | [Unit note](TYPE_B_UNIT_SMITH.md): 2-quotient rows, bounded affine Schubert union, integral Chern generation, and primitive resolution pullback. This Smith interpretation is an additional argument, not the assertion of Birkhoff's theorem. |
| Type-B rank, exponent bound $0\le e\le\min(a,k)$, block duality, stable primitivity, finite sufficient range | Audited deductions | [Structural note](TYPE_B_SMITH.md). Integrality of pullbacks and source changes is essential. |
| The general Schur-product description of each $P_t$ | Audited deduction | [Filtration note](TYPE_B_FILTRATION.md). Cell sections and both primitive pullbacks preserve the entire integer row lattice. |
| Smith multiplicities are successive ranks of $P_t$ | Lattice lemma, now applicable in every rank | [General proof](TYPE_B_GENERAL_SMITH.md) establishes the required primitivity for all prefixes; ranks alone without that hypothesis would be insufficient. |
| Every type-B multiplicity in every rank | Complete repository deduction | [General proof](TYPE_B_GENERAL_SMITH.md): Chow-envelope containment, unit triangular minor, exact beta-number count, and matching Betti bound. No external expert review or novelty claim. |
| Largest possible factor $2^\ell$ in every rank | Audited deduction | [Highest-factor note](TYPE_B_HIGHEST_SMITH.md). A Pieri coordinate-lattice proof determines this without all-prefix primitivity. |
| Three highest multiplicities in every rank; all $\ell\le4$ | Independent repository deduction | [High-layer note](TYPE_B_HIGH_LAYERS.md): complete-symmetric spectators, a coordinate prefix and graph components. This gives an alternative proof of those layers. |
| Full positive type-B partition formula for arbitrary $\ell,e$ | Complete repository deduction; formerly conjectural | [Formula](TYPE_B_CLOSED_CONJECTURE.md) and [integral proof](TYPE_B_GENERAL_SMITH.md). Defect-stratum counting is combined with the integral triangular squeeze. |
| Type-D stable primitivity, sufficient finite range, only powers of two, Euler-sector B comparison, quadric table | Audited deductions | [Type-D note](TYPE_D_SMITH.md), with integral eta specialization, localized Reynolds projection, and the explicit BKT difference basis. |
| Every normalized type-D unit multiplicity | Complete repository deduction | [Unit proof](TYPE_D_FROBENIUS.md): mod-2 Grassmann embedding, exact Frobenius/Euler intersection, derivation, core-one SXP incidence and affine Chern generation. Its unit-layer proof is independent of the later full gluing theorem. |
| Type-D square common kernel/rank, block duality, improved primitive range and one-sided B bound | Repository deductions | [Square comparison](TYPE_D_SQUARE_COMPARISON.md). General square-sector equality is not inferred from rational rank or stable presentations. |
| Type-D square factors for $\min(a,k)\le2$; general square exponent bound | Complete repository deductions | [Two-row proof](TYPE_D_TWO_ROW_SMITH.md): weighted boundary-row elimination proves equality with B; an integral parity retraction of the Chern subring gives $e\le\min(a,k)$. The subsequent general incidence theorem covers all block sizes. |
| General D square-sector B comparison, including all smaller-block-three multiplicities | Complete repository deduction | [Incidence proof](TYPE_D_SQUARE_THEOREM.md): equivariant eta denominator bound, compensated-parity retraction, quaternionic model, and primitive B prefixes give exact row-lattice equality. |
| Explicit full-D Gaussian formula | Complete internally audited repository deduction in every rank | [Formula](TYPE_D_CANDIDATE.md) and [general proof](TYPE_D_EVEN_GLUING_THEOREM.md). The relative Euler extension, primitive actual-head prefixes, and integral row containment close the even-tail step. |
| Type-D two-flag gluing rule | Audited elementary lattice lemma | [Orthogonal note](ORTHOGONAL_SMITH.md). It needs sector exponents and flag-intersection data; it is not a closed partition classification. |
| Full even-tail D partition formula and unit-gluing exhaustion | Complete repository deduction | [General theorem and index argument](TYPE_D_EVEN_GLUING_THEOREM.md): all saturated gluing is at the two unit levels; no hidden positive-level correction. |
| Stronger unweighted typed-prefix strictness | Conditional route; finite checks only in general | [Filtered lattice criterion](TYPE_D_GLUING_FILTRATION.md). The full theorem does not rely on this stronger intermediate assertion. |
| Literal defining/vector SO/Sp lattice is scaled by $2^d$ in degree $d$ | Audited normalization deduction | [Normalization audit](INTEGRAL_LITERATURE.md). Uses the stated character lattice; it does not automatically apply to Spin. |
| 174 Smith lists, 430 prefix comparisons, 516 exact row-module equalities | Finite computational evidence | Exact scopes and input hashes are in the linked [candidate](experiments/type_b_closed_candidate_verification.json) and [filtration](experiments/type_b_filtration_verification.json) reports. None is a proof for arbitrary rank. |

No fatal mathematical defect was found in the retained derivations in this
audit. This is a report of the audit's conclusion, not a guarantee against
error. The repaired citations, explicit geometric maps, and dependence
information are intended to make those derivations easier to check.

Concrete repairs include explicit cell sections and the Sato reference
class in the unit proof, the geometric sections of the localized
characteristic-map citation in D, a source-rank notation collision in the
four-row lattice split, and an obsolete two-row scope sentence.
The factor-eight argument now states its narrower dependency on $P_2$
and the largest-factor deduction, rather than all prefix results.

The new [general verification report](experiments/type_b_general_smith_verification.json)
records 585,243 direct beta-number/corner comparisons in 72 parameter cases,
and 21 fresh exact integer LR module comparisons with the constructed
unit triangular rows, through smaller-block size eight. The [high-layer report](experiments/type_b_high_layers_verification.json)
records 495 full-graph degree components and 56 exact LR module equalities.
These finite checks support the written proofs and are not proofs by themselves.

## The unit polynomial is already classical

Let $n=2r+\epsilon$, $c=r+\epsilon$, $\epsilon\in\{0,1\}$.
In Birkhoff's subgroup-count formula use the partition types

$$
\lambda=(2^r,1^\epsilon),\qquad
\mu=(2^j,1^{\ell-2j}),\qquad
\lambda'=(c,r),\quad\mu'=(\ell-j,j).
$$

The standard formula

$$
\prod_i q^{\mu'_{i+1}(\lambda'_i-\mu'_i)}
{\lambda'_i-\mu'_{i+1}\brack\mu'_i-\mu'_{i+1}}_q
$$

specializes exactly to

$$
q^{j(c-\ell+j)}
{r\brack j}_q{c-j\brack\ell-2j}_q.
$$

Summing over $j$ gives $U_{n,\ell}$. Thus the counting polynomial
belongs to existing theory. Its proposed use for these particular Smith
units requires the separate integral comparison in the unit note.

The original reference is
[Garrett Birkhoff, *Subgroups of Abelian Groups*, Proc. London Math. Soc. 38 (1935), 385--401](https://doi.org/10.1112/plms/s2-38.1.385).
Its metadata was checked; the full original was not accessible in this
audit. The exact general formula was inspected in the
[official Sage documentation](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/q_analogues.html#sage.combinat.q_analogues.q_subgroups_of_abelian_group),
which cites Butler (1987), and in
[Joshua Maglione's original alternate proof, Theorem 1](https://joshmaglione.com/2023/12/17/Birkhoff.html).
These access distinctions matter: the audit does not pretend to have read
an inaccessible original.

This executable check compares the specialization with Sage's independently
implemented Birkhoff and Delsarte algorithms. It checks an attribution
identity, not the additional Springer-Smith comparison.

```python
from sage.all import PolynomialRing, ZZ
from sage.combinat.q_analogues import q_binomial, q_subgroups_of_abelian_group

R = PolynomialRing(ZZ, "q")
q = R.gen()
checked_types = 0
for n in range(2, 17):
    r, c = n // 2, (n + 1) // 2
    ambient_type = [2] * r + [1] * (n % 2)
    for ell in range(1, min(r, 6) + 1):
        for j in range(ell // 2 + 1):
            subgroup_type = [2] * j + [1] * (ell - 2 * j)
            term = (q ** (j * (c - ell + j))
                    * q_binomial(r, j, q) * q_binomial(c - j, ell - 2*j, q))
            for algorithm in ("birkhoff", "delsarte"):
                assert term == q_subgroups_of_abelian_group(
                    ambient_type, subgroup_type, q, algorithm=algorithm)
            checked_types += 1
assert checked_types > 100
```

## Sources and their exact role

* [Buch--Kresch--Tamvakis, isotropic Giambelli](https://math.umd.edu/~harryt/papers/isogiamfinal.pdf),
  Definition 1/Theorem 1, Sections 2.1, 6.1: integral bases, raising operators,
  the Chern-preserving B/C comparison, and signed Schubert conventions.
* [Buch--Kresch--Tamvakis, even orthogonal Giambelli](https://math.umd.edu/~harryt/papers/oggiam.pdf),
  Proposition 2, Theorem 4, Section 5.2: the explicit difference-basis comparison
  and integral eta basis/specialization. The integral upgrade used in the
  Euler-sector deduction is explained in the type-D note.
* [Pragacz--Ratajski](https://www.impan.pl/~pragacz/download/comp3868.pdf),
  Appendix A Sections 10.3, 10.5: topological quaternionic Schur calculus and
  the fibration. [Panin--Walter](https://arxiv.org/pdf/1011.0649),
  Theorems 11.2--11.3, supplies the corresponding bundle presentation.
* [Wildon](https://arxiv.org/pdf/1508.07030), equation (1): the
  Littlewood/SXP input for squared-root Schur coefficients.
* [Anderson](https://doi.org/10.4153/S0008439524000122),
  main theorem p.728 and ordinary specialization pp.729--730:
  integral Chern generation on the affine Grassmannian. Its application
  to the bounded locus and its Smith interpretation are deductions here.
* [Hahn--Nebe--Stanojkovski--Sturmfels](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/shuffle0702final.pdf),
  Section 4, Proposition 14: rectangular nilpotent invariant-subspace varieties
  already have a bounded affine-lattice model. This covers the even
  case directly; the unit note spells out its odd-case Iwahori bounds.
* [Griffin](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/springer-varieties-and-halllittlewood-polynomials/64AC25F936FDD32A3A955031B880E830),
  Section 2.4, equation (2.8): preserved partial flags are Steinberg varieties;
  the Young-invariant comparison displayed there is over $\mathbb Q$.
  It is not a substitute for the integral argument here. Spaltenstein
  flags annihilated on successive quotients have different conditions.
* [Ram, *Subspace profiles, q-Whittaker functions and Krylov methods*](https://www.mat.univie.ac.at/~slc/wpapers/FPSAC2025/115.pdf),
  Definition 1.1 and Theorem 5.2 (also stated as Theorem 1.3):
  an existing general finite-field formula for subspaces of a given
  operator profile. Our defect-$e$ count for $N^2=0$ is its
  $(\ell,e)$ case. This is a comparison of the counted objects; the
  symmetric-function formula has not been transformed term by term to
  the displayed Gaussian sum here. It supplies no integral Smith bridge.
* [Karpenko](https://sites.ualberta.ca/~karpenko/publ/gqf3.pdf),
  Theorem 2.1 and Remark 2.6, and
  [Karpenko--Merkurjev](https://sites.ualberta.ca/~karpenko/publ/duality14.pdf),
  Theorem 1.2: existing integral presentations and duality for the full
  tautological Chern subring. That ring is larger than the selected
  square-only subring; these results do not directly assert its Smith
  multiplicities or saturation.
* [Calmes--Petrov--Zainoulline](https://arxiv.org/pdf/0905.1341),
  Theorems 6.4/6.9 and Sections 13.3, 13.6--13.7, with
  [Totaro](https://www.math.ucla.edu/~totaro/papers/public_html/spin.pdf),
  Theorem 3.2: the localized characteristic-map input used in D.
* [Dimitrov--Roth](https://arxiv.org/pdf/1707.06840), Section 2.8, equation
  (2.8.1): the classical integral quadric ring underlying the small-head
  corollaries.

Parts I--III supply the Springer maps, chosen forms, finite generator
images, and stable injectivity discussed with locators in
[INTEGRAL_LITERATURE.md](INTEGRAL_LITERATURE.md). Their complex injectivity
statements should not be relabeled integral saturation theorems.

* [Tamvakis, double eta polynomials](https://math.umd.edu/~harryt/papers/dbleta.pdf),
  Definition 1 and equations (4)--(5), Theorem 1, Section 4.1:
  the equivariant integral Schubert basis and explicit Chern/Euler
  geometrization used to deduce the relative denominator bound. The
  flagged incidence comparison and selected Smith interpretation are
  additional arguments in the new D square note.

## Literature search limits

The search covered Springer/Borel selected images, integral theta/eta
presentations, tautological Chern subrings, quaternionic Grassmannians,
affine Grassmannian Chern generation, invariant-subspace/Jordan-quiver
Grassmannians, Hall--Littlewood/Birkhoff counting, plethysm and Frobenius
ranks, and subspace profiles/anti-invariant subspaces. Representative query phrases included “isotropic Grassmannians
Smith”, “Chern subring orthogonal Grassmannian”, “Springer Smith normal
form”, “Schur Smith plethysm”, and “Grassmannian Frobenius image
cohomology characteristic two”.

The reviewed sources supply substantial existing structure. No exact
published theorem matching the all-rank selected-image Smith
classification was located. Nor was an exact published Smith-unit
identification located in the sources read. Neither statement proves that
such a result is absent from the literature. Before describing any
deduction here as original research, a specialist literature review is
still needed.
