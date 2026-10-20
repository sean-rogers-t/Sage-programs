# Sage-programs
This is a repository of functions I used to explicitly compute the Borel Morphism.

The new [sage-borel library](borel/README.md) packages the corrected calculations,
supports arbitrary parabolics and finite Lie types, and avoids full Weyl-group
enumeration. It includes independent tests and benchmarks through rank 20.
The original program files are preserved unchanged.

The [exceptional catalogue](borel/EXCEPTIONAL.md) supplies named generator
calculations for G2, F4, E6 and E7, plus the displayed E8 dissertation cases.

The [E6/P1 and E7/P7 character APIs](borel/MINUSCULE_CHARACTER_API.md)
compute the certified integral sources and their complete Borel images.
The [E6 representation comparison](borel/E6_REPRESENTATION_COMPARISON.md)
proves the minimal and adjoint source formulas, their all-Levi image sums
and intersections, and records the independently checked E6/P1 tables.
The [E7 minimal/adjoint comparison](borel/E7_REPRESENTATION_COMPARISON.md)
now provides the corresponding full-source proof, complete E7/P7
comparison and independent audit.
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
