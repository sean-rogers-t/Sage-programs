# Springer–Borel maps: a structural research program

A research proposal and progress inventory prepared on 20 October 2026.
It organizes the published source work and current repository deductions
into four core aims and a six-part extension track, with proof routes,
milestone gates, fallback results, verification requirements and proposed
API additions. The extension track was added on 20 October 2026.

- [Research proposal PDF](output/pdf/Springer_Borel_research_proposal.pdf)
- [Editable LaTeX source](main.tex)
- [Artifact verification record](verification.json)
- [Research ledger: results, partial results and open questions](../../borel/RESEARCH_LEDGER.md)
- [Detailed proof and source status](../../borel/PROOF_STATUS.md)
- [Original papers and verification notes](../README.md)

## Research aims

| Aim | Goal | Existing starting point | First deliverable |
| --- | --- | --- | --- |
| A, primary | Local arithmetic of integral Springer–Borel images | Degree-one source theorem, irreducible lattice classification, B/D formulas and exceptional images | Precise good-prime comparison; then a sharp shared-prime example |
| B | Effective generic degrees from boundary geometry | G2/E6/E7 minimal/adjoint maps and degree calculations | Certified face-scheme method for a coherent family |
| C | Regularity preservation and harmonic localization | Uniform differential product and G2 obstruction witnesses | A reusable obstruction or sufficient condition |
| D, supporting | Geometric explanation of Smith filtrations | Type B primitive prefixes and incidence/defect comparison | An integral chain/incidence comparison proving primitivity |

## Track E: other theories and groups

| Project | Goal and first deliverable | Dependency and current status |
| --- | --- | --- |
| E1: equivariant cohomology | First literal-source pilot completed: SL4/P1 has a first specialization kernel Z/2 in polynomial degree five | Automatic faithful specialization is false; general kernel criteria and a generic backend remain open |
| E2: K-theory | Restrict the associated-bundle map to integral Springer polynomial characters and compare filtered images | Source conventions first; associated-graded agreement with ordinary images is a question, not an assumption |
| E3: oriented theories | Define a formal-group-law-compatible source/image construction, starting with connective K-theory | Use the E2 pilot; additive degree-one generation does not already prove this extension |
| E4: quantum cohomology | Construct one multiplicative lift recovering our map at q=0, then analyze its kernel and image | Define the source, coefficient ring and generator lifts before claiming a quantum Springer image |
| E5: affine Grassmannians | Transport an E4 image through the localized Peterson correspondence and identify information lost under localization | Depends on E4; target uses Pontryagin homology, not ordinary cup-product cohomology |
| E6: Kac–Moody flags and groups | Extend nil-Hecke/characteristic calculations in a specified degree range; separately define or obstruct a Springer analogue | Nil-Hecke foundations are published; trace/completion questions remain foundational |

The proposal includes the established background sources for these theories.
Their representation-dependent image questions remain proposed projects.
Existing isogeny/reductive questions remain in OPEN-10; the new ledger entries
are OPEN-17 (K-theory), OPEN-18 (formal group laws), OPEN-19 (affine transport)
and OPEN-20 (Kac–Moody). OPEN-14/15 retain their equivariant/quantum identities.

Activation proceeds E1, then the E2/E3 comparison, then E4 and E5. E6's
nil-Hecke calculations can be developed independently; they do not establish
an infinite-dimensional trace Springer morphism.

## Progress after the proposal baseline

The first Aim A/E1 work package is now recorded in the
[local image theorem and projective family](../../borel/LOCAL_IMAGE_BOUND.md),
[equivariant specialization counterexample](../../borel/EQUIVARIANT_IMAGE_PILOT.md)
and [proof/source audit](../../borel/LOCAL_IMAGE_AUDIT.md).
The full literal Levi images are integral, and their degree-d cokernels
are killed by t_G e_rho^d, without a Levi Weyl-order factor. The defining
SL_n/P1 factors are exactly n/gcd(n,d). The SL4 equivariant image has a
degree-five specialization kernel Z/2, independently checked with
integer lattices. Good primes recover the full equivariant image, and
each graded specialization kernel is finite with primes dividing t_G e_rho.
The sharp bad-prime benchmark is now supplied by the
[E6 transfer theorem](../../borel/E6_ROOT_WEIGHT_COMPARISON.md) and
[three-case completion](../../borel/E6_REMAINING_CASES.md): together they
determine root-to-weight image factors on all 63 proper standard E6
parabolics. The completion has 52 independently audited base degrees
and a structural P1 extension; [the audit](../../borel/E6_REMAINING_AUDIT.md)
separates published inputs from these deductions. The
[E7/P7 basis theorem](../../borel/E7_ROOT_WEIGHT_COMPARISON.md) then gives
a closed 2-primary multiplicity polynomial and a structural exponent-2
bound on all 64 refinements of P7. It uses a stable common complement
and a root quadratic lift; that exponent argument alone does not give
finer-flag multiplicities.
The [A6 intermediary proof](../../borel/E7_FULL_FLAG_COMPARISON.md) now
settles the full-flag polynomial q P_W/(1+q) and all-degree relative
factors on 109 proper E7 parabolics, including 61/64 refinements of P7.
The [D5 completion](../../borel/E7_D5_REFINEMENT.md) adds two cases, making
111 proper parabolics and 63/64 refinements. The
[D4 completion](../../borel/E7_D4_REFINEMENT.md) now closes the entire
64-refinement family, with polynomial q P_W/((1+q)P_WI) in every case.
This makes 112 proper E7 parabolics. Its 34-generator module certificate
also records a failure of a free rank-10 integral source transfer.
The [complete E7 theorem](../../borel/E7_ALL_PARABOLICS.md) now settles
the relative full root/weight image factors on all **127 proper standard
parabolics**, with every nonunit relative factor equal to two. The
[D5+A1 completion](../../borel/E7_D5_A1_COMPARISON.md) uses a full
25-generator local source and has 366 factors; the
[D4+A1 completion](../../borel/E7_D4_A1_COMPARISON.md) uses eighteen root
generators and a full 34-generator weight-source module and has 3,618.
The [A5-family proof](../../borel/E7_A5_FAMILY.md) gives a common
even-block formula for twelve diagrams, with a halved-block correction.
All-degree source Hilbert identities and separate native image audits
complete the finite classification. A structural explanation of the equal
index corrections and halved-block pattern, and extensions beyond E7,
remain targets.
General bad-prime
structure beyond the covered E6/E7 cases, ambient weight-image cokernels and a structural
criterion for the equivariant specialization kernel remain targets.

These subsequent notes update the live program. The proposal PDF remains
the source-grounded baseline described below, rather than a claim that its
initial targets are all still unresolved.

The proposal uses the notation of the Kumar papers for the maps and
reserves `P` for a parabolic, `X=X(T)` for the character lattice, and
`Q_rt` for the root lattice. It distinguishes literal character sources
from normalized selected generators and requires Borel-image integrality
to be justified separately.

## Next milestones after the E7 completion

Updated 21 October 2026. [PR 10](https://github.com/sean-rogers-t/Sage-programs/pull/10)
is merged into master as `df20de9`. The
[complete relative E7 theorem](../../borel/E7_ALL_PARABOLICS.md) covers all
127 proper standard parabolics; the E6 comparison already covers all 63.
This closes the finite relative-classification benchmark in Aim A.
The following items are proposed work, not additional proved results.
The full-paper project remains on hold.

For a simply connected G, write X=X(T), Q_rt for the root lattice, and
I for retained Levi nodes. In polynomial degree d, put

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}^d Q_{\rm rt})^{W_I}),
 \qquad N_I^d=\beta^{P_I}((\operatorname{Sym}^d X)^{W_I}).
\]

The completed comparison determines N_I/M_I. The degree-d normalized
E7 minimal and adjoint Springer images are 12^d M_I^d and 36^d N_I^d.
The absolute quotient H^(2d)(G/P_I;Z)/N_I^d remains a separate target.

| Priority | Work package | Concrete deliverable | Completion test |
| --- | --- | --- | --- |
| 1, bounded software step | Package the certified relative formulas | A compressed all-parabolic E7 API, then the analogous E6 interface, with representation scaling and proof provenance | Every retained set has a verified record; overlap with existing embedded-image APIs agrees independently; no Schubert matrix is built for a summary query |
| 2, main mathematical effort | Explain the two-primary image algebra structurally | A reusable integral lifting criterion and a mod-two image presentation recovering the E7 exponent and multiplicities by families | Prove the stated lattice/Levi hypotheses and source completeness; identify exactly where the criterion applies and give obstructions outside its hypotheses |
| 3, focused generalization | Explain the halved-block correction | An explicit algebraic filtration, exact sequence or complex accounting for the correction and its q^4 grading | Construct the comparison and prove its integral or local properties; a Weyl-polynomial identity alone is insufficient |
| 4, next arithmetic benchmark | Determine an uncompleted absolute weight-image cokernel | A Schubert-coordinate Smith calculation for E7/P1, degree by degree at 2 and 3 | Reuse all previously available minuscule tables, prove full required sources, and independently check the actual integral embedding |
| 5, after the ordinary-image mechanism | Resume equivariant and K-theoretic pilots | A bad-prime equivariant specialization criterion, followed by one precisely filtered K-theory comparison | Recover the ordinary construction with all source, filtration and specialization hypotheses stated |

The generic-degree/critical-locus project remains an independent Aim B/C
direction. F4 is a natural next geometry benchmark. It is not a test of
a nontrivial root-to-weight quotient: Q_rt=X in F4 and E8. Quantum,
affine and Kac–Moody projects retain the dependency gates in Track E.

### 1. Package the proved information

The current `e6_character_source` and `e7_character_source` factories
package embedded images only for E6/P1 and E7/P7. The new all-parabolic
relative theorem does not automatically extend those constructors.
Start with a separate summary API. A possible interface is:

```text
exceptional_root_weight_comparison("E7", retained_nodes=(1,2,3,4,5,7))
comparison.degree(d)
comparison.compare_scaled_sources(root_scale=a, weight_scale=b, degree=d)
```

These names are a design proposal and are not implemented. A degree
summary should return the cohomology rank, compressed relative factor
multiplicities, coefficient-ring/source conventions, proof family and
certificate provenance. The scalar comparison uses gcd/lcm on the
simultaneous relative Smith basis, as proved in the completion note.
It gives abstract relative groups; an embedded Schubert basis or an
integral polynomial preimage requires a separate calculation.

Make retained nodes explicit; validate node sets and degrees, handle G/G
and degrees outside cohomology, and avoid expanding thousands of equal
Smith factors by default. Ship E7 first from its complete 127-case
manifest, then normalize the E6 records to the same schema, retaining
all 3-power factors. Compare with existing minimal/adjoint image APIs
on their overlapping minuscule cases and test both containment and
incomparability. This is a software consequence of proved mathematics.

### 2. Replace the case tables with a structural mechanism

Work over R=Z_(2) for E7. Let A_X=(Sym_R X_R)^(W_I),
A_Q=(Sym_R (Q_rt tensor R))^(W_I), and K_I=A_X intersect ker(beta), with the
kernel taken in the rational Borel presentation. The exact target is

\[
                         2A_X\subseteq A_Q+K_I.
\]

This is the already-proved image exponent restated before taking Borel
images. The research task is a reusable proof from lattice and Levi
data, with explicit corrections in K_I, and a criterion applicable
beyond the finite E7 list. A source inclusion without these corrections
is a stronger statement and must be checked separately.

Next describe the actual root image
C_I=im(M_I -> N_I/2N_I). The completed exponent theorem identifies it
with M_I/2N_I, and

\[
 m_{I,d}=\operatorname{rank}_R N_I^d-\dim_{\mathbb F_2}C_I^d.
\]

A uniform presentation or filtration of C_I would explain the Smith
multiplicities. Start from the common quadratic lift in the 112-case
family and the triangular-involution source lemma for even blocks;
then account for the D6, D5+A1 and D4+A1 behaviors. Keep reductions of
actual integral invariants distinct from the full modular invariant ring,
and prove every use of invariant lifting before Borel.

Modular invariant theory supplies established tools here. For example,
[Campbell–Shank–Wehlau](https://arxiv.org/abs/0901.2811) give generators
and SAGBI bases for vector invariants of modular cyclic groups. Whether
their methods yield the needed integral lifting and Borel-image theorem
in our setting is a research question, not a consequence already checked.
Our block-coordinate involution also needs an explicit comparison with
their linear representations; invariant lifting through a specialization
must be proved.

### 3. Test and explain the folding pattern

The exact E7 A5-family formula has a Gaussian multinomial for ordinary
even blocks and a correction with halved blocks evaluated at q^4.
That correction can also be expressed using C3 Weyl polynomials. The
next goal is an actual map, filtration or complex producing this term,
and an explanation of the grading and integral primitivity. C_0(q) in
the current formula has signed coefficients; it cannot be declared the
Hilbert series of a graded vector space. A virtual or exact-sequence
interpretation must account for those signs.

Use the full integral involution presentation as the first algebraic
model. For a controlled test outside E7, investigate type A_(2r-1) at
two with r odd and even Levi blocks, starting with SL6 and then SL10.
Here X_R/Q_rt,R has order two. Begin with low-degree full-source kernels
on the small Grassmannian cases before attempting a general formula.
Derive the correct statement for that target; the E7 quotient relations
and its formula cannot simply be transplanted. A proved mechanism or a
precise counterexample to a proposed extension is a useful deliverable.

### 4. Start absolute images from a genuine gap in coverage

Inventory the existing E6/P1 and E7/P7 embedded character-image tables
first. Their adjoint images can be descaled using the known source
normalizations to recover existing weight-image information. The next
new benchmark is E7/P1: 126 Schubert classes, dimension 33, and audited
two-local full source components already available in the D6 calculation.

Reuse those certificates at two. At three, establish the full local
source separately from the actual D6 lattice and signed-permutation
coordinates. The existing two-local generator completeness does not
establish three-local completeness.

Embed the full weight image N_I^d in the integral Schubert basis and
determine H^(2d)(G/P_I;Z)/N_I^d at the allowed bad primes 2 and 3.
Check the embedding with a second BGG/Chevalley construction. The general
local bound already controls prime support; its existence does not give
the sharp graded Smith factors.

[Totaro's torsion-index calculations](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf)
are published constraints and reference checks. The torsion index is
defined by the top-degree characteristic image on the full flag; it
does not by itself specify the degreewise cokernels of a partial flag.
After a verified absolute benchmark, formulate a structural prime/exponent
question before expanding the case catalogue.

### Proof and literature gates

Before claiming a new structural result, compare its precise lattice,
normalization and source with the published characteristic-map,
torsion-index and modular-invariant literature. The two references above
are starting points checked for this roadmap, not a completed novelty
audit. Preserve the distinction between the known inputs, the existing
computer-assisted E6/E7 deductions and proposed generalizations.

The recommended next work package is **1**, followed by the mathematical
core **2**, with **3** as its focused test. Begin **4** after inventorying
the existing absolute tables. No new API or theorem is asserted by this
roadmap, and the frozen proposal PDF is unchanged.

## Status and provenance

The mathematical progress inventory is frozen to repository commit
`e9fcb7f3b7d28bea9c248abab1f01d8de28250aa`. PDF links to working notes are
pinned to that commit so the baseline remains identifiable. Its runtime
verification counts are historical baseline counts, not a claim that this
proposal reran the runtime suite. Subsequent status belongs in the ledger.

Published foundations, internally audited deductions, finite calculations
and proposed goals are separated. The document adds no theorem and makes
no novelty or external-review claim. Its bibliography is a starting map,
not a completed novelty search. A full proof-paper project remains on hold.

## Build

The source has an inline bibliography and compiles with Tectonic or another
LaTeX installation supporting the listed standard packages. From this folder:

```powershell
New-Item -ItemType Directory -Force output/pdf | Out-Null
tectonic --keep-logs --outdir output/pdf main.tex
Move-Item -Force output/pdf/main.pdf output/pdf/Springer_Borel_research_proposal.pdf
```

The PDF is also available through the desktop LaTeX project's live preview.
The delivered artifact was compiled, rendered and visually reviewed; the
verification record documents the exact checks and source/PDF hashes.
