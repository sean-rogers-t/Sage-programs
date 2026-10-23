# sage-borel

A small SageMath library for exact Borel characteristic maps in a Schubert basis.
The historical programs in the repository root are preserved. This package is a
separate implementation incorporating the corrections found during their audit.

## Documentation

The guides explain both what each method returns and the mathematics that
makes the calculation valid. Start with the method reference for usage, or
the mathematical guide for the Borel/Springer construction.

| Question | Guide |
| --- | --- |
| When does Levi branching prove generic smoothness and predict excess dimension? | [Affine branching theorem, singular-coefficient exclusion and uniform A/C consequences](AFFINE_BRANCHING_FIBERS.md) |
| Do fundamental representations always give isolated quotient-critical points? | [No: an all-rank spin counterexample, smooth excess fibers, and the general unused-coordinate criterion](SPIN_FIBER_EXCESS.md) |
| Is there a general character criterion for fiber components and multiplicities? | [Corrected quotient potentials, exact local Jacobian algebras, linear elimination and the limits at nonregular points](GENERAL_CAYLEY_FIBER_CRITICAL.md) |
| What can we prove uniformly for simple groups and whole target strata? | [All-type fundamental zero fibers, classical target bounds, and exact middle-exterior component counts and multiplicities](SIMPLE_GROUP_CAYLEY_BRANCHING.md) |
| Does the component formula extend to primitive symplectic representations? | [All even-rank equality strata: signed component counts, multiplicity one, a stronger kernel bound, and the odd-rank residual equation](SYMPLECTIC_TARGET_STRATA.md) |
| What is the full fiber in the symplectic zero model? | [All-rank component and generic reducedness; a singular embedded thickening and the reduced-germ conjecture](SYMPLECTIC_FIBER_COMPONENT.md) |
| Does the realization theorem extend beyond type A? | [Every odd nilpotent map occurs on a symplectic centralizer in the fixed primitive middle-exterior family](SYMPLECTIC_CENTRALIZER_REALIZATION.md) |
| Are the collapsed unipotent varieties actual full-fiber components? | [Middle-exterior component theorem: parity of the block size determines enlargement and multiplicity](MIDDLE_EXTERIOR_FIBER_COMPONENT.md) |
| Can a single irreducible representation realize these degenerations? | [All type A nilpotent polynomial maps occur on proper centralizers in the fixed middle-exterior family](IRREDUCIBLE_CENTRALIZER_REALIZATION.md) |
| Which degeneration orders can actual representations realize? | [Type A character criterion, closed exterior-mixture classification and every-order realization in even dimension](TYPE_A_CHARACTER_DEGENERATION.md) |
| Can I reconstruct on B for larger classical and exceptional groups? | [Automatic Chevalley data, general character integration, exact torus-wall evaluation and independent checks through E8](ROOT_BOREL_API.md) |
| When do mixed root terms recover vertical rank? | [General centralizer rank reduction and the complete G2 minimal/adjoint rank rules](G2_BOREL_RANK.md) |
| Does the fiber program go beyond known theory? | [Literature assessment, uniform deductions, a reducible-fiber obstruction and the decision to pause more rank tables](CAYLEY_FIBER_ASSESSMENT.md) |
| Can I reconstruct Borel coordinates from a torus formula for other representations? | [Uniform root-commutator recurrence, SL3/Sp4 character API and independent representation checks](BOREL_RECONSTRUCTION.md) |
| Can I evaluate the Springer map on B, including unipotent coordinates? | [SL3/Sp4 defining formulas, exact collapse fibers, general triangular mechanism and API](BOREL_COORDINATES.md) |
| What are the defining-adjoint degrees in arbitrary type B/C rank? | [Closed formulas, binomial boundary corrections, Spin/SO and Sp/PSp distinctions](BC_MIXTURE_DEGREES.md) |
| How do F4 degrees change for minimal-adjoint direct sums? | [All-multiplicity formula, exceptional ratios a=b and a=5b, and height-two corrections](F4_MIXTURE_DEGREES.md) |
| What is the F4 minimal Springer degree, and is its torus map finite? | [Degree 96, odd-crosspolytope criterion, finite-flat corollary and independent fibers](F4_MINIMAL_DEGREE.md) |
| Does the half-spin deformation method transfer to an exceptional group? | [E6 adjoint correction 1890, stationary middle layer and all-orders crossing exclusion](E6_ADJOINT_DEFORMATION.md) |
| Which degree arguments extend beyond Spin, and what can fail? | [General criteria, normal-layer corrections, exact counterexamples and independent E6 audit](GENERAL_CAYLEY_DEGREE.md) |
| Can I compute classical equivariant Springer images? | [Type A partial flags in H/K and the integral symplectic projective K-cokernel](CLASSICAL_EQUIVARIANT.md) |
| What higher type C Fitting formulas are proved? | [Doubled presentation](SYMPLECTIC_K_FITTING.md), [Schur generators](SYMPLECTIC_K_MINORS.md), [all-rank integral Hankel balancing and complete ranks 1-9](SYMPLECTIC_K_BALANCING.md) |
| What are all ordinary type C projective K Smith factors? | [Closed binary-digit multiplicities, Hankel proof and known inputs](SYMPLECTIC_K_ORDINARY_SMITH.md) |
| What general results do we have, and which structural questions remain? | [Research ledger: results, partial results, counterexamples and open problems](RESEARCH_LEDGER.md) |
| Which Springer and Cayley questions are already covered by the literature? | [Audit of all eight brief references, the linked notebook, related work and current-status limits](CAYLEY_LITERATURE_AUDIT.md) |
| Which exterior-power Springer maps preserve regularity? | [All-rank SL_n classification, uniform rational counterexamples and independent trace checks](EXTERIOR_REGULARITY.md) |
| Which irreducible SL_n Springer maps preserve regularity? | [Complete highest-weight classification, symmetric-power corollary and independent character checks](TYPE_A_REGULARITY.md) |
| Which fundamental Sp_(2n) Springer maps preserve regularity? | [All-rank obstruction, constructive character zeros and contraction-kernel verification](TYPE_C_FUNDAMENTAL_REGULARITY.md) |
| What about every irreducible Sp_(2n) highest weight? | [Universal central--Weyl obstruction, independent branching proof, symmetric powers and central quotients](TYPE_C_REGULARITY.md) |
| What about type B, including spinorial highest weights? | [All irreducible Spin_(2n+1) and all SO_(2n+1) representations fail regularity for n>=2](TYPE_B_REGULARITY.md) |
| Can mixing central characters repair regularity? | [Complete B/C answer: symplectic torus-positive characters, nonsemisimple failure, and arbitrary-sum type B obstruction](MIXED_REGULARITY.md) |
| How much localization does symplectic harmonic separation need? | [Exactly det(1+g): minimal invariant principal localization and full Cayley determinant](SYMPLECTIC_HARMONIC_LOCALIZATION.md) |
| How do symplectic generalized Cayley degrees change under tensor powers and mixtures? | [Degree 2n deg(F) for actual characters F(det(1+g)), explicit fibers and representation-operation laws](SYMPLECTIC_CAYLEY_DEGREES.md) |
| How do Spin generalized Cayley degrees change under tensor powers? | [Degree 2m floor(N/2), an explicit fiber presentation, and parity-dependent descent to SO](SPIN_CAYLEY_DEGREES.md) |
| What happens when even and odd spin characters are mixed? | [Degree 2 floor(N/2) deg(F), precise SO descent, representation operations and collapsed character levels](SPIN_POLYNOMIAL_CAYLEY_DEGREES.md) |
| Where are polynomial Spin maps finite, and where do their inverses branch? | [Maximal finite target locus, full Jacobian, boundary rank bound and a false discriminant factor](SPIN_CAYLEY_GEOMETRY.md) |
| What changes for one half-spin representation? | [Pfaffian formula, complete zero fiber in the Cayley chart and a large boundary zero family](HALF_SPIN_CAYLEY.md) |
| What is the entire half-spin zero fiber, including nonsemisimple points? | [All-rank four-eigenvector criterion and sharp dimension r(2r-1)-12](HALF_SPIN_BOUNDARY.md) |
| What is the generic degree of either individual half-spin map of Spin_10? | [Degree 148: a toric intersection proof and independent original-weight checks](HALF_SPIN_DEGREE.md) |
| Is there a half-spin degree formula across ranks? | [Degree 2^(r-1)(r^2-3r-2)+4r on Spin_(2r), with a deformation proof and independent Spin_12 checks](HALF_SPIN_RANK_DEGREES.md) |
| Is the simple-group regularity table complete, including reducible representations and exceptional types? | [Complete table and uniform primitive-root obstruction, with rank-one and isogeny qualifications](REGULARITY_CLASSIFICATION.md) |
| What did the exploration of other cohomologies and affine groups produce? | [Research report, theorem inventory, known inputs and prioritized API path](OTHER_THEORIES_PROGRAM.md) |
| How does the actual Springer-character image change in K-theory? | [General filtered comparison and counterexamples](OTHER_THEORIES_FINITE.md), [complete SL3 formulas](OTHER_THEORIES_SL3.md), [primitive character data](OTHER_THEORIES_COTANGENT.md) |
| Can the earlier Cayley determinants control equivariant K-images? | [All-parabolic determinant theorem, local Smith patterns and an integral symplectic family](OTHER_THEORIES_EQUIVARIANT_K.md) |
| What can already be proved for quantum and affine images? | [Flat quantum projective families](OTHER_THEORIES_QUANTUM.md), [all-rank affine A/C divisor contents and evaluation/level images](OTHER_THEORIES_AFFINE.md) |
| What comes after the complete E6/E7 relative classification? | [Prioritized next milestones, proposed API and structural proof targets](../articles/research-program-proposal/README.md#next-milestones-after-the-e7-completion) |
| What mathematical map does each method compute? | [Methods and worked examples](METHODS.md) |
| Why do the Borel, BGG, Chevalley and coinvariant methods work? | [Mathematical framework and conventions](MATHEMATICS.md) |
| How do I multiply Schubert classes or lift them to polynomials? | [Schubert multiplication and inverse Borel calculations](SCHUBERT_CALCULUS.md) |
| Which integer classes do fixed generators produce? | [Integral images and Smith factors](INTEGRAL_IMAGES.md) |
| How do two integral images compare, and what do mixed products add? | [Lattice comparisons](INTEGRAL_COMPARISONS.md) and [two-step joint-image tables](TWO_STEP_INTEGRAL.md) |
| Can I get the B/D Smith factors without a matrix? | [Fast orthogonal formula API](ORTHOGONAL_API.md) |
| What is already known about integral Springer images? | [Literature and normalization audit](INTEGRAL_LITERATURE.md) |
| Which primes can obstruct a full literal Springer–Borel image? | [Local bound, all-Levi proof and defining SL_n projective factors](LOCAL_IMAGE_BOUND.md) |
| Can the E6 root-to-weight factors be determined sharply? | [60-case transfer theorem](E6_ROOT_WEIGHT_COMPARISON.md) and [three-case completion](E6_REMAINING_CASES.md) cover all 63 proper parabolics; [tables](E6_REMAINING_TABLES.md), [completion audit](E6_REMAINING_AUDIT.md) |
| What happens to the E7 root-to-weight factors at 2? | [Explicit P7 simultaneous bases and multiplicity polynomial](E7_ROOT_WEIGHT_COMPARISON.md); exponent 2 on all 64 refinements, with [proof/source audit](E7_ROOT_WEIGHT_AUDIT.md) |
| What is proved for finer E7/P7 flags? | [Full flag and A6 transfer](E7_FULL_FLAG_COMPARISON.md), [D5 proof](E7_D5_REFINEMENT.md), [D4 completion: all 64 refinements in all degrees](E7_D4_REFINEMENT.md) |
| What about the other E7 parabolics? | [Complete theorem for all 127 proper E7 parabolics](E7_ALL_PARABOLICS.md), including [D5+A1](E7_D5_A1_COMPARISON.md), [D4+A1](E7_D4_A1_COMPARISON.md) and the [uniform even-block formula](E7_A5_FAMILY.md) |
| What could explain the E7 multiplicities structurally? | [Quadratic-derivation criterion, bounded 112-case proof and arbitrary-rank type A model](E7_QUADRATIC_DERIVATION.md) |
| Does the first integral E7 half-sum respect Levi invariants? | [Degree-eleven lifting theorem for all 112 cases and cubic nonlifting obstructions](E7_DEGREE11_LIFTING.md) |
| When does central parity preserve an entire integral Borel kernel? | [Sufficient integral hypotheses, all-rank type C theorem and Spin7 counterexample](CENTRAL_PARITY_KERNEL.md) |
| Does the Spin7 obstruction persist for every odd spin group? | [Every-rank quartic theorem, stable Spin9 presentation and complete ranks 3–8](SPIN_PARITY_EXTENSION.md) |
| Is there an arbitrary-rank recurrence for the odd-spin kernels? | [Integral Pieri lattice recurrence, uniform stable family and remaining classification gap](SPIN_PARITY_RECURSION.md) |
| Can the remaining odd-spin obstruction be isolated at the rank boundary? | [Exact adjacent-rank boundary criterion and obstruction-rank formula](SPIN_PARITY_BOUNDARY.md) |
| How do I query all E7 relative Smith factors and compare representations? | [Compressed all-parabolic API, source scaling and minimal/adjoint comparisons](E7_RELATIVE_IMAGE_API.md) |
| Can I query all E6 relative factors, including 9 and 27? | [All 63 parabolics, representation comparisons and shared-prime scaling](E6_RELATIVE_IMAGE_API.md) |
| What is the complete literal G2 full-flag image? | [Integral kernel, all Smith factors, API and independent Chevalley audit](G2_FULL_FLAG.md) |
| Does an equivariant image always specialize faithfully? | [SL4/P1 counterexample and independent integer-lattice checks](EQUIVARIANT_IMAGE_PILOT.md) |
| Which claims are published, deductions here, or conjectures? | [Proof and literature status](PROOF_STATUS.md) |
| What general Smith-factor formulas can be proved in types B and D? | [Orthogonal formulas and proofs](ORTHOGONAL_SMITH.md) |
| Why is the returned kernel complete? How do backends and normal forms work? | [Complete presentations and proofs](PRESENTATIONS.md) |
| What do a degree matrix, rank and homogeneous relation mean? | [Degreewise images and relations](IMAGE_RELATIONS.md) |
| Which exceptional inputs and normalizations are available? | [Exceptional catalogue](EXCEPTIONAL.md) |
| What integral images do the exceptional generators produce? | [Complete G2/F4 tables and normalization comparisons](EXCEPTIONAL_INTEGRAL.md) |
| What happens for the minuscule E6/P1 and E7/P7 spaces? | [Complete tables, independent checks and image subalgebras](EXCEPTIONAL_MINUSCULE.md) |
| Which E6/P1 polynomials actually pull back to integral Levi characters? | [Character integral form and complete image comparisons](E6_CHARACTER_FORM.md) |
| What is the E7/P7 character form, and does it contain the catalogue image? | [E7 character form, sums and intersections](E7_CHARACTER_FORM.md) |
| How does the representation determine the integral polynomial-character source? | [General criterion and character API](CHARACTER_SOURCES.md), with [complete G2 sources and image rings](G2_CHARACTER_FORM.md) |
| Is the full torus character source always generated in degree one? | [Trace-source theorem, proof, further representations and independent audit](TRACE_SOURCE_THEOREM.md) |
| Can I classify exceptional irreducible sources without listing all weights? | [Closed lattice formula, proof and fast API](IRREDUCIBLE_CHARACTER_LATTICES.md) |
| What are the quotient map, generic degree and critical locus of a Springer morphism? | [Geometry API, complete G2 formulas, independent matrix audit and E6/E7 quadratic coordinates](SPRINGER_GEOMETRY.md) |
| What are all higher E6 quotient coordinates and the full determinant? | [Minimal/adjoint formulas, exact interpolation proof, root products and independent matrices](E6_GEOMETRY.md) |
| What are all higher E7 quotient coordinates and the full determinant? | [Minimal/adjoint formulas, certified degree-18 calculations and independent matrices](E7_GEOMETRY.md) |
| What are the generic degrees of the E6/E7 minimal and adjoint Springer maps? | [Degree API, toric boundary proof and independent checks](EXCEPTIONAL_DEGREES.md) |
| What is the literal F4 character form and its complete image? | [Proof, API and independent checks](F4_CHARACTER_FORM.md), [all four parabolic tables](F4_CHARACTER_TABLES.md) |
| How does the F4 adjoint representation change the integral image? | [All-degree scaling theorem and comparison API](F4_REPRESENTATION_COMPARISON.md) |
| How do I compute the certified E6/P1 and E7/P7 character images? | [Public character API](MINUSCULE_CHARACTER_API.md) and [E6 minimal/adjoint comparison](E6_REPRESENTATION_COMPARISON.md) |
| How does E7's adjoint representation change the integral image? | [E7 source proof, complete comparisons and API](E7_REPRESENTATION_COMPARISON.md) |
| How do the historical functions map to the new API? | [Migration guide](MIGRATION.md) |

[examples/mathematical_methods.py](examples/mathematical_methods.py) is an
executable companion covering projective space, the mixed C3 image, the C9
kernel and the type D Euler class.

The [exceptional generator catalogue](EXCEPTIONAL.md) adds three-parameter
calculations such as `exceptional_image("F4", 1, "e2")`, with named polynomial
inputs and computed Schubert images. It includes the corrected Part III inputs,
the missing E6/P6 Euler generator, and the displayed E8 dissertation cases.

The [image-and-relations calculator](IMAGE_RELATIONS.md) analyzes any chosen
homogeneous Levi-invariant generators: exact graded image matrices, ranks,
kernel relations and Hilbert-function prefixes. It includes a two-step
symplectic family and reproducible experiments on joint Levi-factor images.

The [complete-kernel calculator](PRESENTATIONS.md) adds `A.presentation()` for
any chosen generator algebra. It returns the full defining ideal as a reduced
Groebner basis, a quotient basis, the exact Hilbert series, and a finite
completeness certificate. No guessed degree cutoff is required.

Version 0.22 adds the all-parabolic E6 interface to
`exceptional_root_weight_comparison`, retaining factors 3,9,27. It also
adds `g2_character_algebra()` for G2/B and an adjoint selector for all
supported G2 spaces. The [G2 full-flag note](G2_FULL_FLAG.md) proves the
complete integral kernel and gives the absolute Smith factors in every degree.

Version 0.21 adds `exceptional_root_weight_comparison("E7", retained_nodes=...)`
for all 127 proper parabolics and G/G. It returns compressed relative Smith
factors, proof provenance, and comparisons of scaled root/weight sources,
including the normalized minimal and adjoint Springer representations.
The [guide](E7_RELATIVE_IMAGE_API.md) explains the mathematics and the
distinction between these relative groups and absolute Schubert-image factors.

Version 0.5 uses direct coinvariant normal forms by default for standard-numbered
types A/B/C/D. This avoids Schubert expansion when computing kernel relations;
`A.presentation(backend="bgg")` retains the original method for cross-checking.
The formerly slow C9/P8 tail example now completes with kernel `(q^9)` and
dimension 9; [the benchmark](benchmarks/coinvariant-results.json) records setup,
certification and a fresh certificate recheck.

Version 0.6 adds arbitrary Schubert multiplication, cup powers, intersection
numbers and inverse Borel calculations for full and partial flags in every
supported type. Polynomial lifts use exact degreewise linear algebra or BGG
divided differences; low-degree exceptional calculations avoid constructing
the full positive-root discriminant. See the
[guide](SCHUBERT_CALCULUS.md) and
[executable example](examples/schubert_calculus.py).

Version 0.7 adds `A.integral_degree(d)`: exact integer image lattices, Smith
factors, saturation indices, membership tests and integer polynomial preimages.
It checks the fixed generator images without silently clearing denominators.
The [integral guide](INTEGRAL_IMAGES.md) distinguishes normalized Springer
generators from literal integral polynomial characters; the
[literature audit](INTEGRAL_LITERATURE.md) explains the classical type C
saturation result and the scope of the orthogonal and exceptional questions.

The [orthogonal research note](ORTHOGONAL_SMITH.md) records deductions that normalized
type B/D factors are powers of two and that the stable images are saturated.
It includes a type B exponent bound and block symmetry, closed quadric
formulas, and an exact type D integral gluing rule. The
[general type B proof](TYPE_B_GENERAL_SMITH.md) now gives every normalized
Smith multiplicity in every rank. A Chow envelope, integral triangular
rows, beta-number counting and a Betti bound close the previously missing
comparison. The [formula](TYPE_B_CLOSED_CONJECTURE.md) matches every
recorded Smith list, and the [high-layer proof](TYPE_B_HIGH_LAYERS.md)
independently determines the three largest factors in every rank. The [proof/source audit](PROOF_STATUS.md)
identifies published ingredients and checked repository arguments separately.
The unit counting polynomial is classical Birkhoff/Hall theory; novelty of
the Smith interpretations and other deductions has not been established.

For type D, the [unit-layer proof](TYPE_D_FROBENIUS.md) now gives a
closed formula for every normalized factor equal to one, including the
exact Gaussian overlap for even tails. A separate
[square comparison](TYPE_D_SQUARE_COMPARISON.md) proves block symmetry,
a sufficient primitive range and a termwise B bound. The
[two-row square theorem](TYPE_D_TWO_ROW_SMITH.md) gives an independent
small-block proof and a general exponent bound. The subsequent
[general D square theorem](TYPE_D_SQUARE_THEOREM.md) now proves equality
of every normalized B/D square-sector factor in every rank and degree.
Its relative incidence proof uses the equivariant eta denominator bound,
an integral compensated-parity retraction preserving the base Euler
class, and the established primitive B prefixes. It gives explicit
multiplicity polynomials for factors 1, 2, 4 and 8 in the smaller-block
three case. The [general even-tail theorem](TYPE_D_EVEN_GLUING_THEOREM.md)
now proves the [full D formula](TYPE_D_CANDIDATE.md) in every rank:
all saturated gluing comes from the known unit overlap. Its relative
Euler extension and precise head-Euler B comparison are checked in the
[full audit](TYPE_D_EVEN_GLUING_AUDIT.md). The new verifier passes
138 full components and resolves both old rank-nine timeouts.
The independent small-block proofs and stronger conditional prefix
criterion remain separately labeled.
Version 0.8 packages the proved B/D formulas as
`orthogonal_smith(family, n, k, degree)`. It returns compressed exact
multiplicities, image rank, saturation-index exponents, and the D unit-gluing
count without constructing a Weyl group or Schubert matrix. Square and Euler
sectors and the scaled character normalization are explicit options.
See the [API and mathematics](ORTHOGONAL_API.md) and
[executable examples](examples/orthogonal_smith.py). The general matrix API
remains the method for membership, preimages, and arbitrary generators.

Version 0.9 adds exact integral image containment, sums, intersections and
quotients, preserving their Schubert coordinates. The
[comparison guide](INTEGRAL_COMPARISONS.md) explains the API and integer
linear algebra. The [two-step experiment](TWO_STEP_INTEGRAL.md) records
96 degree calculations in types C, B and D, including mixed-product
contributions and independent Cartesian checks. These bounded tables
are separate from the general single-tail B/D theorems.

Version 0.10 adds `springer_character_map` for a specified highest-weight
representation: normalized torus formulas, the exact integral linear
lattice, an all-primes leading-term certificate, and bounded searches for
modular obstructions. Certificates distinguish a proof from an inconclusive
test. `g2_character_algebra(1)` and `(2)` expose the proved full character
sources and work with the existing Borel, image and lattice APIs.
See the [mathematics and API](CHARACTER_SOURCES.md) and the
[full regression record](character-verification.json).

Version 0.11 adds `f4_character_source(node)` and the general
`IntegralCharacterSource`: exact degreewise intersections with integral
characters, retaining source elements with zero Borel image.
`source.degree(d).image()` gives the embedded Schubert lattice and supports
the existing comparison API. The [F4 proof](F4_CHARACTER_FORM.md) identifies
the full torus source as Sym_Z(6 Q_long) and gives all four complete
maximal-parabolic image tables.

Version 0.12 adds `f4_character_source(node, representation="adjoint")`.
The [representation comparison](F4_REPRESENTATION_COMPARISON.md) proves
that the adjoint torus source is Sym_Z(18 Q_long), and its degree-d Levi
sources and Borel images are 3^d times their minimal counterparts.
The normalized Springer morphisms themselves are different; the equality
compares lattices in common polynomial and cohomology coordinates.

Version 0.13 adds `e6_character_source(1)` and `e7_character_source(7)`,
with exact degreewise sources and embedded images. E7 uses integral
symmetric coordinates to keep coefficient matrices manageable.
`e6_character_source(1, representation="adjoint")` computes the newly
certified Sym_Z(24P) source. The minimal Sym_Z(6Q) and adjoint images
are incomparable in every positive degree on E6/P1. Their sums and
intersections satisfy the [proved all-Levi formulas](E6_REPRESENTATION_COMPARISON.md);
the [API guide](MINUSCULE_CHARACTER_API.md) explains verification and limits.

Version 0.14 adds `e7_character_source(7, representation="adjoint")`.
Its [full-source proof](E7_REPRESENTATION_COMPARISON.md) gives
Sym_Z(36P), compared with the minimal Sym_Z(12Q), and proves all-Levi
sum/intersection formulas. The complete E7/P7 tables and independent
audit distinguish degreewise containment from incomparability of the
whole image rings. The original minimal factory default is preserved.

Version 0.15 proves [degree-one generation for all normalized trace
representation maps](TRACE_SOURCE_THEOREM.md), including higher weights
and direct sums. `certify_trace_source()` verifies saturated coefficient
columns and their exact relation to weight exponents. The general
`IntegralCharacterSource` now uses this theorem by default; explicit
leading-basis requests retain their earlier behavior. The Levi-invariant
ring and its cohomology image still require their own calculations.

Version 0.16 adds `irreducible_character_lattice(type, highest_weight)`.
The [classification proof](IRREDUCIBLE_CHARACTER_LATTICES.md) gives
Lambda_lambda=I_lambda*(Q+Z*lambda)^*, using long roots of squared length
two. It computes the dimension, index, central-character order and full
linear lattice from root data alone. All exceptional irreducibles have
the listed one or two lattice shapes; a direct-sum counterexample explains
why the irreducible hypothesis is necessary. Full Laurent formulas remain
available through `springer_character_map`.

The subsequent [type C joint-image proof](TYPE_C_JOINT_PRIMITIVITY.md)
proves saturation of the normalized joint image for all positive block
sizes. A common integral basis also proves primitivity of the separate-image
sum and that the two separate images intersect only in the constants.
The note gives the integral kernel and closed rank formulas. A subsequent
[proof and literature audit](TYPE_C_JOINT_AUDIT.md) found no gap under the
stated hypotheses and classifies the results as corollaries of classical
integral theory in these Springer/Borel coordinates. No originality claim
is made; finite basis checks are recorded separately.

The [orthogonal joint-image proof](ORTHOGONAL_JOINT_IMAGES.md) extends
the primitive-sum and constant-intersection results to types B and D.
It gives integral kernel and rank formulas; D has an additional Euler
relation. Ambient joint saturation fails in general, but only powers of
2 can occur as nonzero Smith factors. Their closed multiplicities remain
further work. These are deductions from classical theory, with explicit
finite kernel and lattice checks, not a claim of originality.

The [extension roadmap](EXTENSIONS.md) separates proposed API work from
further mathematical questions. The full-paper project is on hold.

It computes

\[
\beta:\operatorname{Sym}(\mathfrak t^*)^{W_P}\longrightarrow H^{2*}(G/P;\mathbb Q),
\qquad \beta(\omega_i)=\sigma_{s_i}.
\]

Here `i` is an omitted node; on the full flag variety every node is allowed.

This is the characteristic-map stage of the research calculations. For a
representation-ring element written as the Springer pullback `q(Theta(t))`, pass
the invariant polynomial `q(omega)` in fundamental-weight coordinates to `expand`.
The independent exceptional Springer calculations remain in
[`../articles/calculations`](../articles/calculations/).

## Run it

Requires SageMath and Python 3.9 or later; verified with SageMath 9.5 / Python 3.10.
Run from this `borel/` directory without installing anything:

```bash
PYTHONPATH=src sage -python examples/basic.py
```

The exceptional catalogue has its own example:

```bash
PYTHONPATH=src sage -python examples/exceptional.py
```

For graded images and relations:

```bash
PYTHONPATH=src sage -python examples/image_relations.py
```

For complete kernel presentations:

```bash
PYTHONPATH=src sage -python examples/kernel_presentations.py
```

For the worked mathematical explanations:

```bash
PYTHONPATH=src sage -python examples/mathematical_methods.py
```

For Schubert products, polynomial lifts and intersection numbers:

```bash
PYTHONPATH=src sage -python examples/schubert_calculus.py
```

For integral images and normalization comparisons:

```bash
PYTHONPATH=src sage -python examples/integral_images.py
```

Or install into your Sage environment and import from any directory:

```bash
sage -pip install .
sage -python examples/basic.py
```

On this Windows machine these commands run inside the Ubuntu WSL installation.
Use `**` for exponentiation in Python files. In a Sage notebook, `^` is also
preparsed by Sage. This package requires the Sage runtime; ordinary Python alone
does not provide its root systems or exact polynomial arithmetic.

```python
from sage_borel import BorelMap, classical_borel

# Your old type C parameters: rank n=12, tail size k=10.
b = classical_borel("C", 12, 10)
answer = b.tail_elementary(10, 5)  # beta(e_5(x_3^2,...,x_12^2))
print(answer)
print(answer.to_terms())          # JSON-safe words and exact coefficients
print(b.last_stats)               # frontier sizes and operation counts

# The corrected type D example from the audit:
d = classical_borel("D", 4, 3)
assert d.tail_elementary(3, 1).to_terms() == [
    {"word": [2, 1], "coefficient": "-1"}
]
```

`classical_borel("A", n, k)` uses **GL(n) dimension**, matching the old type A
programs. `BorelMap("A3")` instead uses **Dynkin rank** and means SL(4).
The factory omits node `n-k`. Its ranges are A: `1 <= k < n`; B/C:
`0 <= k < n`; D: `2 <= k < n`. Use `BorelMap` directly for spin or other
parabolics.

## General inputs and parabolics

```python
b = BorelMap("A3", omitted_nodes=(1, 3))
w1, w2, w3 = b.weights
answer = b.expand((w1 + 2*w3)**3 + 7)
print(b.coefficient((w1 + 2*w3)**3, [1, 2, 3]))

f = BorelMap("F4")
print(f.expand(f.weights[3]**3))

e = BorelMap("E7", omitted_nodes=7)
print(e.divisor_product([7, 7, 7]))
```

Supported types are finite irreducible crystallographic A, B, C, D, E, F, G.
`omitted_nodes=None` gives the full flag variety; an integer gives a maximal
parabolic; an iterable gives an arbitrary parabolic; `()` gives G/G. Retained
nodes generate the Levi Weyl group. Inputs to a parabolic map must be invariant
under those nodes; the library checks this and rejects invalid inputs.

Weights use Sage's Dynkin numbering. Polynomial degree `d` means cohomological
degree `2*d`. A squared tail elementary polynomial of index `i` has polynomial
degree `2*i`, hence cohomological degree `4*i`. Type A coordinate weights have
sum zero; the GL(n) central character has zero characteristic-map image.

## API

| Operation | Meaning |
| --- | --- |
| `b.ring`, `b.weights`, `b.W` | Exact QQ polynomial ring, fundamental weights, Weyl group |
| `b.expand(f)` | Full Schubert expansion, including zero, constants and mixed degrees |
| `b.coefficient(f, word)` | One coefficient without generating other classes |
| `b.schubert(word)` | One Schubert basis class; `[]` gives the unit |
| `b.schubert_polynomial(word, method="auto")` | A homogeneous Levi-invariant polynomial lifting that class |
| `b.lift(expansion, method="auto")` | A polynomial representative of any Schubert expansion |
| `b.schubert_product(u, v)` | Product of two Schubert classes indexed by words |
| `b.multiply(a, c)` or `a*c` | Cup product of two Schubert expansions |
| `a**n` | Nonnegative integral cup power, including the unit for exponent zero |
| `b.integrate(a)`, `b.intersection_number(*classes)` | Top-degree integral or integral of a product |
| `b.reflect(f, i)`, `b.divided_difference(f, i)` | Simple reflection and exact BGG divided difference |
| `b.chevalley(i, word)` | Multiply a Schubert class by the divisor `omega_i` |
| `b.divisor_product(indices)` | Repeated divisor product, with empty product equal to 1 |
| `b.multiply_divisor(expansion, i)` | Multiply an existing sparse expansion by a divisor |
| `b.classical_coordinates()` | Standard coordinate weights for A/B/C/D |
| `b.tail_polynomial(k, i)` | Literal tail elementary input polynomial |
| `b.tail_elementary(k, i)` | Its Borel image, using the complementary shortcut when cheaper |
| `b.head_polynomial(r, i)`, `b.head_elementary(r, i)` | Unsquared head elementary polynomial and image |
| `b.euler_polynomial(k)`, `b.euler(k)` | Type D tail Euler product and image |
| `b.basis(d)` | Minimal right-coset representatives in degree `d`, generated on demand |
| `b.poincare_dual(word)` | Dual index `w0 * w * w0P` |
| `b.longest_representative()` | Top Schubert index, without enumerating W |
| `b.clear_cache()` | Release generated bases, words and per-class products |
| `to_type_a(w)` | Structural signed-permutation embedding for B/C/D |
| `exceptional_image(type, parabolic, generator)` | Compute a named exceptional generator image |
| `exceptional_generators(type, parabolic)` | List names, degrees and source metadata |
| `exceptional_polynomial(type, parabolic, generator)` | Inspect its fundamental-weight polynomial |
| `exceptional_borel(type, parabolic)` | Reusable map with named generator methods |
| `A = GeneratorAlgebra(b, generators)` | The graded map from named formal variables to the chosen classes |
| `A.degree(d)` | Exact degree-d image matrix and homogeneous kernel space |
| `D = A.integral_degree(d, include_ambient=False)` | Integer image lattice, Smith factors and saturation data for fixed generators |
| `D.contains(a)`, `D.preimage(a)` | Integer membership and an integer-coefficient formal preimage |
| `K = A.presentation(backend="auto")` | Full kernel ideal, quotient basis, Hilbert series and completeness certificate |
| `K.reduce(p)` | Canonical representative of a formal expression in the selected image |
| `K.verify(backend=None)` | Recheck completeness with the original or an explicitly selected backend |

`tail_polynomial` and `tail_elementary` default to unsquared coordinates in A,
squared coordinates in B/C/D; `squared=True/False` overrides this. Classical
coordinate helpers and the signed embedding require standard Dynkin numbering;
generic Borel operations also accept integer relabelings of Dynkin nodes.

An `Expansion` is an immutable mapping from Weyl elements to rational
coefficients. `answer[word]` returns zero for an absent valid word;
`answer.as_dict()` returns a copy; `answer.to_terms()` serializes reduced words
and rational numbers as strings. Addition, rational scalar multiplication,
arbitrary cup products and nonnegative integral powers are provided. To add a
constant class, use `b.schubert([])` for the unit. Polynomial lifts are chosen
representatives, not unique inverses; `b.expand(b.lift(answer)) == answer`.
Weyl elements
must have the map's parent; passing a list of Dynkin indices converts a word
explicitly. Parabolic divisor operations require omitted nodes.

## Why this is faster

The original routines repeatedly enumerate every Weyl element, reconstruct
roots, expand all subsets and accumulate duplicate term lists. At type C rank 8,
W already has 10,321,920 elements.

The replacement follows nonzero divided differences only. Braid-equivalent
paths merge once, rather than being counted as distinct contributions. Divisor
products use cached positive-root covers. Elementary polynomials use dynamic
recurrences; tail images can use
`e_i(tail^2) = (-1)^i h_i(head^2)` in the coinvariant quotient. This often makes
the input much smaller. [The mathematical conventions and algorithm](MATHEMATICS.md)
explain why these optimizations preserve the result.

Measured results, including setup times and complete-degree coefficient checks,
are in [`benchmarks/results.json`](benchmarks/results.json). The benchmark script
records actual timings on the local machine; these depend on hardware, degree,
input size and output size. Large dense expansions can still be expensive.

```bash
PYTHONPATH=src sage -python benchmarks/run.py
```

## Verification and migration

The [worked mathematical companion](examples/mathematical_methods.py) contains
exact assertions for the explanations in the guides. The
[documentation checker](tools/verify_documentation.py) runs the guides' Python
blocks in reading order and checks local links; its latest result is in
[documentation-verification.json](documentation-verification.json):

```bash
PYTHONPATH=src sage -python tools/verify_documentation.py
```

The independent test oracle uses explicit coordinate swaps/sign changes and
Schubert representatives derived from the positive-root discriminant. It does
not reuse the library's Chevalley root/coroot product formula. The exhaustive
suite checks all **1,088** divisor products in B2, C2, B3, C3 and D4, plus
classical generators, mixed degrees, parabolics, duality/intersection pairings,
historical counterexamples, exceptional-type sanity checks, and rank-10 word
handling. See [`verification.json`](verification.json) for counts and environment.

The [Schubert-calculus tests](tests/test_products.py) independently check
arbitrary products against Cartesian divided differences, Grassmannian
Littlewood–Richardson tables and the middle intersections on a type D quadric.
They also check both polynomial-lift methods, partial flags, relabeled nodes,
mixed classes and low-degree exceptional inputs. The exact coverage and latest
run are recorded in [products-verification.json](products-verification.json).

The [integral tests](tests/test_integral.py) independently check Smith factors
by determinantal divisors and classical coefficient matrices by Cartesian BGG.
The [integral verification record](integral-verification.json) and
[bounded Springer-image experiment](experiments/integral-springer-results.json)
record the exact cases, generator normalizations and integer preimage witnesses.

```bash
PYTHONPATH=src sage -python -m unittest discover -s tests -v
PYTHONPATH=src BOREL_EXHAUSTIVE=1 BOREL_TEST_REPORT=verification.json \
  sage -python -m unittest discover -s tests -v
```

[`MIGRATION.md`](MIGRATION.md) maps your old functions to the new API and records
the corrected examples. [`legacy_manifest.json`](legacy_manifest.json) records
canonical Git and local checkout SHA-256 hashes of the 35 unchanged historical
program files. `python verify_legacy.py` verifies repository contents;
`python verify_legacy.py --checkout-bytes` also verifies this machine's recorded
checkout bytes. These are recorded separately because Windows Git converts line
endings on checkout.
