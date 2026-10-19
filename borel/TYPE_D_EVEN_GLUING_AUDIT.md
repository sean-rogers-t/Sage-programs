# Internal audit of the general even-tail gluing proof

**Audit verdict: PASS, 19 October 2026.** The integrated proof in
[TYPE_D_EVEN_GLUING_THEOREM.md](TYPE_D_EVEN_GLUING_THEOREM.md) was
read adversarially after its supporting lemmas were reviewed separately.
No material mathematical gap remains in the reviewed argument.
This is an internal research-session audit of an unpublished repository
deduction. It is neither external expert certification nor a claim of
originality in the literature.

## What the proof establishes

For the fixed normalized inputs \(q_i=e_i(x_{\rm tail}^2)\) and
\(p=\prod x_{\rm tail}\), the proof establishes the full Smith
multiplicity formula for every \(a\ge1\), \(k\ge2\). For even
\(k\), every saturated square/Euler gluing direction is represented
at the original unit levels, with the stated Gaussian overlap
polynomial. Odd tails follow from the previously proved separate
sector lists and their opposite degree parities.

All indices concern the saturation of the selected rational image.
The proof does not classify arbitrary mixtures of head Chern classes,
and does not assert symmetry between two full tail-image problems.
Literal defining/vector polynomial-character normalization is a
separate degree-dependent scaling already recorded in the notes.

## Reviewed dependencies and integral steps

The audit checked the following chain of arguments.

1. **Actual-head type-B prefixes.** Section 2 of
   [TYPE_D_HEAD_EULER_B_COMPARISON.md](TYPE_D_HEAD_EULER_B_COMPARISON.md)
   extends the integral envelope/triangular squeeze of
   [TYPE_B_GENERAL_SMITH.md](TYPE_B_GENERAL_SMITH.md) to \(H>K\).
   Restricting the unbalanced index to its feasible range and capping
   the source corner by \(K\) preserves the beta-number criterion.
   The candidate \(j=b-d\) meets every new feasibility bound.
   The defect rank bound and the complete-coordinate unit minor then
   prove actual-head primitivity; weighted block duality alone is not
   used to infer a prefix theorem.

2. **The relocated B source.** For positive incidence level, the
   square and head-Euler sources occupy the distinct integral
   Lagrangian Leray--Hirsch components \(1,c_t(U^*)\). The narrower
   Euler rectangle gives precisely the \((H,K-1)\) square prefix.
   At level zero, the full-column SXP sign changes by \((-1)^h\).
   The exact integer overlap identities, their signed identity minor,
   and the proved mod-two overlap dimension give matching rational
   and mod-two ranks. The resulting primitivity is only claimed over
   \(\mathbb Z_{(2)}\), which suffices: the corresponding full type-C
   source is primitive and B/C rescaling excludes odd-primary factors.
   The adapted-chain computation replaces exactly the stated
   number of zero exponents by ones.

3. **The formal relative Euler algebra.** Sections 1--6 of
   [TYPE_D_RELATIVE_EULER_EXTENSION.md](TYPE_D_RELATIVE_EULER_EXTENSION.md)
   prove polynomial freeness over explicit squared/Euler generators
   integrally, including characteristic two. Direct-summand base
   change makes \(\mathcal C\) projective over the flag base ring.
   Thus the uses of flatness in the annihilator sequences are valid.
   The normal form at \(Z=0\), the identities
   \(\operatorname{ann}_R(e)=cv_KR\) and
   \(ceR\cong R/(Z,v_K)\), and the quotients proving
   \(R\oplus ceR\) primitive were checked with both kernel
   containments. The formal ring injects into cohomology by rational
   Borel/Leray--Hirsch and its independently proved torsion freeness.

4. **The specified relative source.** Section 7 of that note identifies
   its two independent source images with primitive actual-head
   prefixes. The Euler quotient uses the smaller ambient quaternionic
   bundle, followed by a primitive flag pullback. Every branching
   partition has width at most \(K-1\), so the restriction preserves
   its entire coefficient submodule, basis by basis.

5. **Relative parity and denominators.** Section 8 retains the all-even
   and all-odd head parity patterns. Averaging uses even head sign
   changes, which preserve the nonzero base Euler relation \(ep=Z\).
   It is an integral coefficient retraction on the polynomial subring,
   with descent justified rationally and then by torsion-free
   cohomology. The equivariant eta formula has denominator at most
   \(2^t\), including at most one additional typed central half.
   Flagging the base supplies global integral Schubert lifts and
   integral Leray--Hirsch. These assertions were checked directly
   against Tamvakis's Definition 1, Theorem 1, and geometrization in
   Section 4.1 of
   [Double eta polynomials and equivariant Giambelli formulas](https://math.umd.edu/~harryt/papers/dbleta.pdf).

6. **The final lattice implications.** In the integrated proof,
   source primitivity in \(\mathcal C\), the retraction, and the
   denominator bound show every selected incidence factor divides
   \(2^t\). Equality of rational kernels and primitive prefix rows
   then give \(2^tN_t\subseteq T_t\); this is not a rank-only
   inference. Integral incidence restriction and the level-zero rows
   produce the reverse B/D row containment. Finally the exact index
   identity gives \(\dim W=I\). The independently proved unit
   intersection already has dimension \(I\), so it exhausts \(W\).

The earlier dependencies were also checked for their role here:
[TYPE_B_FILTRATION.md](TYPE_B_FILTRATION.md),
[TYPE_D_SQUARE_THEOREM.md](TYPE_D_SQUARE_THEOREM.md),
[TYPE_D_SQUARE_COMPARISON.md](TYPE_D_SQUARE_COMPARISON.md),
the integral Euler comparison in [TYPE_D_SMITH.md](TYPE_D_SMITH.md),
the unit theorem in [TYPE_D_FROBENIUS.md](TYPE_D_FROBENIUS.md), and
the lattice lemma in [ORTHOGONAL_SMITH.md](ORTHOGONAL_SMITH.md).
The published Schubert/Giambelli and quaternionic presentation inputs
remain identified in those notes; the Smith deductions are repository
arguments, not attributed as published Smith theorems.

## Corrections and boundary checks

An earlier relative-source draft said that all quaternionic Schur
classes of width \(K\) vanish after lowering the tail rank. This is
false for nontrivial ambient bundles: wider relative classes can reduce
to expressions involving base classes. The reviewed version uses only
the correct statement that the width-\(\le K-1\) coordinate submodule
maps isomorphically to the target Schur basis. Every selected Euler
branch lies in that submodule, so the corrected source proof is intact.

The audit separately checked \(t=0\), the trivial even-sign group at
\(t=1\), \(b=0\), and the relocated \(K=1\) two-fork boundary.
The point convention at \(K-1=0\) agrees with the quaternionic
quotient and the source ranks. No nonexistent outer-involution
splitting over a base with nonzero Euler class is assumed.

The stronger unweighted-prefix intersection assertions in
[TYPE_D_GLUING_FILTRATION.md](TYPE_D_GLUING_FILTRATION.md) remain a
separate conditional route. The completed proof uses the relocated
incidences and does not require an unproved strictness assertion from
that route.

## Independent calculation scope

The computations support implementation and transcription; none
replaces the preceding mathematical arguments. The reports overlap
and their counts should not be added as disjoint samples.

- [The actual-head and signed-overlap report](experiments/type_d_head_euler_b_verification.json)
  checks 106 \(H>K\) degree components and 610 exact integer prefix
  modules, plus 92 ordinary-Grassmannian level-zero components and
  35 exact signed overlap identities. Its verifier uses Schur/LR and
  ordinary plethysm, without a Borel or Weyl-group calculator.
- [The head-three filtration report](experiments/type_d_gluing_filtration_verification.json)
  checks 132 components, 528 exact module intersections, level-zero
  quotient detection, and the weighted module calculation. It also
  checks 27 explicit unit relations before row weighting.
- [The full even-gluing report](experiments/type_d_even_gluing_verification.json)
  computes 138 saturated gluing dimensions from full and sector Smith
  indices, including 132 head-three components, 16 independent full
  Cartesian crosschecks, and the two formerly timed-out components.
- [The higher-head filtration report](experiments/type_d_higher_gluing_filtration_verification.json)
  records 83 full components and 447 exact integer prefixes, with
  further full Cartesian and model crosschecks.
- [The extended Chern comparison report](experiments/type_d_extended_chern_comparison_verification.json)
  checks 78 exact B/D row-module equalities for the specified relocated
  source, with three independent Cartesian B/D comparison pairs.
  The component \(H=6,K=4,d=24\) reached its 90-second resource bound
  and is explicitly excluded. A separate counterexample records why
  the comparison must not be generalized to arbitrary Chern mixtures.

Each report records selected source and dependency hashes; the repository
commit supplies the full source snapshot. Sage versions and resource
exclusions are recorded separately. The proof's scope and unpublished
status are unchanged by the finite agreements.
