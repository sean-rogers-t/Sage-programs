# Audit of the three remaining E6 comparisons

Date: 20 October 2026. The statement under review is the
[three-case completion theorem](E6_REMAINING_CASES.md). Completion status
is controlled by its all-degree certificate reports, not by agreement of
selected examples. This is an internal mathematical and implementation
audit; literature priority and external review remain separate.

## Claims and dependencies

| Step | Reason it is valid | Scope or limitation |
|---|---|---|
| Integral objects | Use the actual E6 lattices Q inside X, full Levi-invariant symmetric sources, and their embedded characteristic images M inside N. Since 3X is contained in Q, N/M is finite 3-primary. | Determines N/M, not the ambient cokernel H/N. All degrees are half the cohomological degree. |
| Subsystem coordinates | Literal Cartan calculations give an A5+A1 root subsystem, root-basis determinant plus or minus 2 and weight-basis determinant plus or minus 1. | These are 3-local decompositions, not asserted integral decompositions at 2. |
| Weyl action | Each original retained simple reflection acts by the stated adjacent permutation of the x variables; the extra node 2 flips z. | Checked in the original six-dimensional weight lattice. |
| Complete weight invariants | Symmetric functions in each permutation block generate. Taking the sum-zero quotient commutes with invariants here because each permutation monomial module has H1=0 over the torsion-free coefficient ring. | Does not rely on averaging by a Weyl order divisible by 3. |
| Rational Borel kernel | Literal 27 weights give invariant moments of degrees 2,5,6,8,9,12. Their finite quotient proves independence. Classical E6 invariant degrees and Hilbert series prove generation. Apply the rational Borel theorem. | An application of published characteristic-map theory; an integral Borel presentation is not assumed. |
| Finite degree cutoff | The exact quotient has ranks 72 and 720, with top degrees 21 and 29. The Grassmannian extension has a monic 20-element relative basis. | Establishes exhaustiveness of the finite calculation, not an extrapolation from low degrees. |
| Complete root invariants | Every weight-invariant monomial is expanded in a complete permutation-invariant ring on root differences. A basis B satisfies CB integral over Z_(3) and CB has full column rank modulo 3. | The least-denominator argument proves completeness, including all rational kernel combinations before taking the image. |
| Embedded images | Every full source column is reduced in rational cohomology. Local membership and primitive-column tests verify the proposed image bases. | An arbitrary rational cohomology basis suffices for relative lattices; it need not be an integral Schubert basis. |
| Relative Smith groups | The inclusion matrix in these two image bases has factors 1,3,9,27. | Independently checked modulo 81; the bound is also tested on the inverse inclusion matrix. |
| Third parabolic | A common free basis (1,z), with z squared equal to 4c2, decomposes both images. The quotient is the base quotient plus its degree shift. | A projective-bundle/Leray--Hirsch deduction; no high-degree extrapolation for the third case. |

The primitivity certificate deserves emphasis. Integrality of CB alone
would only prove that some root invariants were found. Full column rank
modulo 3 excludes any additional rational coefficient vector with a
3-denominator, and is what proves that the full source was used.
Likewise, clearing denominators after taking the Borel quotient would
give a different problem. The implementation imposes source integrality
before it applies the quotient map.

## Independent checks

The [structural check](experiments/e6-remaining-structure.json) uses plain
Python/SymPy, a literal Cartan matrix, independently solved coordinate
equations, and exact reflection orbits. It checks the root and weight
lattices, the 27 weights, every retained reflection and the fiber pairing 2.
It also reconstructs the Poincare polynomials from the E6 and Levi Weyl
degrees. The generated tables must agree with all of their graded ranks,
not just with their total dimensions.
Pure-power leading terms independently bound every base variable;
enumerating every standard monomial within those bounds gives dimension
36 and top degree 20, so this cutoff is also checked without assuming it.

The [independent base construction](experiments/e6-independent-rational-base.json)
starts from individual branching weights and symmetrizes their literal
powers. The main implementation instead uses Newton recurrences directly
in elementary functions. The all-degree audits use SymPy grlex Groebner
bases and a separately written sparse normal-form reducer; the main
experiment uses Singular and a different monomial order. The coefficient
expansions are independently reconstructed using Newton sums. Every
source basis and every image column is replayed. Reported image bases
are verified by local membership and primitivity, rather than trusted
because the primary computation returned them.

The final inclusion factors are computed by independent Python list
operations modulo 81, in place of the primary integer Smith algorithm.
The audit also checks that 27 times the inverse inclusion matrix is
3-integral. Sage exact polynomial and matrix arithmetic is shared in
this replay; independence does not extend to every arithmetic backend.

The all-degree reports are [A5](experiments/e6-remaining-a5-audit.json)
and [A2+A2+A1](experiments/e6-remaining-a2a2a1-audit.json).
The [direct BGG check](experiments/e6-remaining-bgg-check.json) separately
computes Schubert columns in degrees 0--6 of both base cases and 0--5
of the third case: 20 degree/case comparisons. This uses a different
mathematical algorithm, but is not an all-degree Schubert recomputation.

During development, the independent SymPy reducer initially omitted
division by a Groebner polynomial's leading coefficient when SymPy
selected an integer coefficient domain. That produced a failed check.
The reducer now explicitly uses QQ and divides by the actual leading
coefficient; the affected audit was restarted. This was an implementation
error in the checker, not evidence for a counterexample to the theorem.

## Published inputs and what is deduced here

- The rational Borel characteristic theorem, Weyl invariant degrees,
  symmetric-polynomial identities, Grassmannian relations, and
  Hermite/Smith theory are established tools.
- [Duan--Zhao, Theorem 3](https://arxiv.org/abs/math/0511332) already
  gives the ambient integral E6/P2 ring (their SU(6) notation corresponds
  to our derived Levi A5). We do not claim a new cohomology presentation.
- [Totaro, Theorem 6.1](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf)
  gives adjoint E6 torsion index 54, while the simply connected index is
  6. The all-parabolic bound 27 and full-flag top comparison 9 follow
  from published torsion-index and divided-difference theory. Neither
  is a new result of these tables.
  Totaro's proof already uses a Levi isogenous to A5 times a torus to
  bound the 3-part. That structural use of an A5 Levi is also established;
  it does not itself supply the degreewise relative multiplicities here.
- The characteristic-map foundations can also be followed through
  [Demazure](https://doi.org/10.1007/BF01418790) and
  [Calmes--Petrov--Zainoulline](https://arxiv.org/abs/0905.1341).
- The deductions here are the complete degreewise relative Smith
  multiplicities in the three exceptional retained-node cases, together
  with the structural degree-shift relation and the assembly of the
  all-parabolic E6 comparison. The literal Springer scalings come from
  the earlier source theorem, ultimately using Rogers's Springer
  formula and Kumar's polynomial-character/Borel construction.

Targeted searches included E6 parabolic characteristic images,
root/weight lattices, A5 and A2+A2+A1 invariant sources, and adjoint
torsion indices. These checks did not locate the exact all-degree
relative factor lists. This is a limited literature check, not a proof
of originality or a claim that these questions were posed in the cited
papers. A computational deduction can be valid without having established
literature priority.

## Reproduction and data provenance

From borel/, using Sage for the two image calculations and audits:

```sh
python experiments/build_e6_independent_base.py
python experiments/e6_remaining_structure_check.py
PYTHONPATH=src sage -python experiments/e6_remaining_images.py --case A5
PYTHONPATH=src sage -python experiments/e6_remaining_images.py --case A2A2A1
python experiments/package_e6_remaining.py
sage -python experiments/audit_e6_remaining_images.py --case A5
sage -python experiments/audit_e6_remaining_images.py --case A2A2A1
PYTHONPATH=src sage -python experiments/e6_remaining_bgg_check.py
```

The large calculation was assembled from exact per-degree checkpoints;
its recorded elapsed time is for the last resumed segment, not total
research time. The independent audit replays every saved degree and
binds checked prefixes by individual row hashes. Resuming it requires
unchanged checker code and unchanged checked rows. Source basis choices
from earlier algorithm revisions are acceptable only because the full
membership and primitivity certificates are independently rechecked.

The large complete coefficient bases are retained in a deterministic
gzip archive referenced by the readable JSON summary. Its SHA256 is
over raw bytes; text hashes use CRLF-to-LF normalization. The archive
is proof data, not an installed API feature. Runtime and wheel remain
version 0.20.0. Documentation and release audits bind the new notes,
scripts and reports without rerunning unchanged runtime regressions.
