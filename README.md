# Sage-programs
This is a repository of functions I used to explicitly compute the Borel Morphism.

The new [sage-borel library](borel/README.md) packages the corrected calculations,
supports arbitrary parabolics and finite Lie types, and avoids full Weyl-group
enumeration. It includes independent tests and benchmarks through rank 20.
The original program files are preserved unchanged.

The [exceptional catalogue](borel/EXCEPTIONAL.md) supplies named generator
calculations for G2, F4, E6 and E7, plus the displayed E8 dissertation cases.

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
The [integrated D audit](borel/TYPE_D_EVEN_GLUING_AUDIT.md) records the proof
dependencies and independent finite checks. These are unpublished repository
deductions, without a novelty claim or external certification. The unit
counting polynomial itself is classical Birkhoff/Hall enumeration.

Research articles and the independently checked F4, E6 and E7 Springer-morphism
note are indexed in [articles/README.md](articles/README.md).

The completed dissertation proof, type A and D extensions, comparison with
Parts II and III, and exact calculation checks are in
[the combinatorial proof note](articles/combinatorial-proofs/README.md).
