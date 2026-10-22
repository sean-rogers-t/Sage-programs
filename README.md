# Sage-programs
This is a repository of functions I used to explicitly compute the Borel Morphism.

The [research ledger](borel/RESEARCH_LEDGER.md) tracks uniform results, exact
cases, partial arguments, counterexamples and open questions. The priority
is generalizable mathematics, with case calculations supporting that goal.
The [next milestones](articles/research-program-proposal/README.md#next-milestones-after-the-e7-completion)
prioritize packaging the complete relative formulas, explaining their
structure, and then extending the absolute and equivariant image results.

The new [sage-borel library](borel/README.md) packages the corrected calculations,
supports arbitrary parabolics and finite Lie types, and avoids full Weyl-group
enumeration. It includes independent tests and benchmarks through rank 20.
The original program files are preserved unchanged.

The all-parabolic [E7](borel/E7_RELATIVE_IMAGE_API.md) and
[E6](borel/E6_RELATIVE_IMAGE_API.md) summary APIs now package
the complete relative Smith data and minimal/adjoint representation
comparisons, with source normalization and proof provenance. It avoids
Schubert matrices for summary queries; absolute cohomology/image cokernels
remain a separate calculation.

The [G2 full-flag calculation](borel/G2_FULL_FLAG.md) completes the
minimal-representation images on every proper standard G2 parabolic.
It proves the integral image presentation and records all absolute Smith
factors, with an independent native root/coroot Chevalley audit. Adjoint
and other irreducible images follow by the proved source-scaling law.

The [exceptional catalogue](borel/EXCEPTIONAL.md) supplies named generator
calculations for G2, F4, E6 and E7, plus the displayed E8 dissertation cases.

The [E6/P1 and E7/P7 character APIs](borel/MINUSCULE_CHARACTER_API.md)
compute the certified integral sources and their complete Borel images.
The [E6 representation comparison](borel/E6_REPRESENTATION_COMPARISON.md)
proves the minimal and adjoint source formulas, their all-Levi image sums
and intersections, and records the independently checked E6/P1 tables.
The [sharp E6 root-to-weight comparison](borel/E6_ROOT_WEIGHT_COMPARISON.md)
and its [three-case completion](borel/E6_REMAINING_CASES.md) now determine
the relative 3-primary image factors on all 63 proper standard parabolics,
with [complete new tables](borel/E6_REMAINING_TABLES.md) and a separate
[proof/source audit](borel/E6_REMAINING_AUDIT.md).
The [E7 minimal/adjoint comparison](borel/E7_REPRESENTATION_COMPARISON.md)
now provides the corresponding full-source proof, complete E7/P7
comparison and independent audit.
The [E7/P7 root-to-weight basis theorem](borel/E7_ROOT_WEIGHT_COMPARISON.md)
explains all its relative 2-primary factors and proves exponent 2 for
all 64 refinements of P7, with a separate [audit](borel/E7_ROOT_WEIGHT_AUDIT.md).
The [E7 refinement extension](borel/E7_REFINEMENT_MULTIPLICITIES.md) proves
multiplicities through degree nine on all 64 refinements, with independent
checks and an explicit obstruction to a six-generator source shortcut.
The [A6 intermediary proof](borel/E7_FULL_FLAG_COMPARISON.md) now gives
all-degree full-flag multiplicities and relative factors on 109 proper
E7 parabolics, including 61 of those 64 refinements. The
[D5 completion](borel/E7_D5_REFINEMENT.md) adds both D5 cases, making
111 proper parabolics and 63/64 refinements. The
[D4 completion](borel/E7_D4_REFINEMENT.md) now proves the common formula
q P_W/((1+q)P_WI) for all 64 refinements in all degrees, and makes
112 proper parabolics in total. The [complete E7 theorem](borel/E7_ALL_PARABOLICS.md)
now determines the relative full root/weight image factors on **all 127
proper standard parabolics**. Every nonunit relative factor is two.
The completion includes [E7/P6](borel/E7_D5_A1_COMPARISON.md), with 366
factors, [D4+A1](borel/E7_D4_A1_COMPARISON.md), with 3,618, and a common
[even-block formula](borel/E7_A5_FAMILY.md) for twelve A5-family diagrams.
All-degree integral source proofs and separate native image audits support
the formulas. Published inputs, finite E7 deductions and unresolved
structural generalizations are distinguished.
The [general trace-source theorem](borel/TRACE_SOURCE_THEOREM.md) now
proves degree-one generation for every supported normalized representation
map, with a 59-case independent weight/lattice audit and a public certificate.
The [irreducible lattice classification](borel/IRREDUCIBLE_CHARACTER_LATTICES.md)
now gives the full linear source from the highest weight and trace index,
without enumerating representation weights, for all exceptional groups.
It includes a general proof, a fast API and separate exact-arithmetic checks.

The [Springer geometry API](borel/SPRINGER_GEOMETRY.md) computes adjoint
quotient coordinates and full critical determinants. It includes complete
G2 minimal/adjoint formulas of generic degrees 6 and 18, independent matrix
checks, and the quadratic coordinates for E6/E7 minimal and adjoint maps.
The [E6 geometry extension](borel/E6_GEOMETRY.md) now supplies all six
quotient coordinates and both full 78-dimensional determinant products,
with exact unisolvent certificates and independent representation matrices.
The [E7 geometry extension](borel/E7_GEOMETRY.md) supplies all seven
coordinates through degree 18 and both full 133-dimensional determinant
products, with exact certificates and independent 56/133-dimensional matrices.
The [generic-degree calculation](borel/EXCEPTIONAL_DEGREES.md) gives E6
degrees 1080/2754 and E7 degrees 7272/33664 for the minimal/adjoint Springer
maps on simply connected groups, with geometric proofs and independent audits.

The [graded image calculator](borel/IMAGE_RELATIONS.md) computes exact relations
among chosen generator images, with two-step symplectic flag experiments and
an invariant-theory derivation of their joint-image presentation.

The [proof and literature status](borel/PROOF_STATUS.md) separates published
inputs, deductions checked in this repository, finite evidence and conjectures.
The audit does not establish originality or replace external mathematical review.
The [orthogonal Smith-factor notes](borel/ORTHOGONAL_SMITH.md) record audited
deductions for type B/D restrictions, stable primitivity, quadric formulas and
an integral gluing rule. The [general type B proof](borel/TYPE_B_GENERAL_SMITH.md)
and [general type D proof](borel/TYPE_D_EVEN_GLUING_THEOREM.md) now give complete
normalized Smith multiplicity formulas in every rank. In D, all saturated
square/Euler gluing comes from the explicitly counted unit overlap.
The [fast formula API](borel/ORTHOGONAL_API.md) now exposes these B/D
multiplicities without Schubert matrices, including explicit sector and
normalization choices.
The [integrated D audit](borel/TYPE_D_EVEN_GLUING_AUDIT.md) records the proof
dependencies and independent finite checks. These are unpublished repository
deductions, without a novelty claim or external certification. The unit
counting polynomial itself is classical Birkhoff/Hall enumeration.

Research articles and the independently checked F4, E6 and E7 Springer-morphism
note are indexed in [articles/README.md](articles/README.md).

The completed dissertation proof, type A and D extensions, comparison with
Parts II and III, and exact calculation checks are in
[the combinatorial proof note](articles/combinatorial-proofs/README.md).
