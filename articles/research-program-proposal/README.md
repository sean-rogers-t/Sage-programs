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
