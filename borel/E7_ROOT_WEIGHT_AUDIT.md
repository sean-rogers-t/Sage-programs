# Proof and source audit: E7/P7 at 2

Date: 20 October 2026. Statement audited:
[the explicit simultaneous image-basis theorem](E7_ROOT_WEIGHT_COMPARISON.md).
This is an internal computer-assisted proof audit, not external peer
review, formal verification or a certification of originality.

## Scope and mathematical audit

| Step | Checked argument | Status and limit |
|---|---|---|
| Objects and prime | M and N are images of full root and weight invariant symmetric sources. Since X/Q=Z/2, N/M is finite 2-primary. | Localization at 2 recovers the entire integer relative group. The source is unscaled; literal Springer factors are restored separately. |
| Cohomology ring | Use Duan--Zhao's E7/P7 presentation, with published Schubert classes h,s,t of degrees 1,5,9. | Published input; no new cohomology-ring theorem. |
| Candidate N basis | Three towers h^i (i<=27), h^i(2s), h^i(2t) (i<=13). Six explicit ring identities prove spanning and closure over Z_(2). | Rational dimension 56 proves independence; these need not be ambient integral Schubert bases. |
| Candidate M basis | Multiply exactly the tower elements with odd h exponent by 2. The coefficients of all parity-changing reductions are even with odd denominators. | Proves this module is R[2h,h2,2s,2t]. It does not alone prove it is the full root image. |
| Actual source integrality | A5 permutation coordinates are actual root/weight bases of determinant -1. Symmetric functions in six permutation variables detect integral coefficients exactly. | The full E6 invariants sit inside this smaller symmetric ring; averaging by a Weyl order is unnecessary. |
| Four explicit lifts | A=2omega7, B=-psi2/144, C=(35AB2-psi5/5760)/81, D=(psi9/387072+1666AB4+243B2C)/19683. | Both implementations check coefficient integrality at 2; the auxiliary psi quotients are integer root polynomials. Their images are 2h,h2,2s,2t. |
| Full-source certification | In every degree 0--27 retain every rational invariant monomial. Verify proposed source columns are integral and primitive modulo 2 in coefficient space. | The least-denominator argument proves completeness, including all rational Borel-kernel combinations before taking the image. |
| Full-image certification | Reduce every certified source image into the proposed M/N basis, checking local coordinates and full row rank modulo 2. | Excludes missing generators. There are 1,967 rational invariant monomials, checked for both full sources. |
| Multiplicity formula | In the common basis, inclusion is diagonal with entries 1 or 2 according to h-exponent parity. | Counts give q[14]_(q2)+(q6+q10)[7]_(q2), with no 4 or 8. This is a basis proof, not inference from two ambient Smith lists. |
| Exhaustion | Published top degree is 27; independent leading-ideal checks recover dimension 56 and the cutoff. | The finite certificates cover the whole P7 image, not a conjectural degree range. |
| Refinement bound | Actual lattices split at 2 as V plus Rw and V plus R(2w), with index 3 and w fixed by W_E6. A common root quadratic maps to h2. | N=M+hM and 2h belongs to M, hence 2N is contained in M for all 64 retained subsets of K. This proves an exponent bound only on the refinements. |

The index-3 stable-complement calculation is independent of the torsion
index of the E6 fiber. The latter is not a unit at 2, so the earlier
common-Schubert-fiber basis argument cannot be applied automatically.
The sharper exponent bound here instead uses the adjoining-generator
argument and a root-image quadratic lift. These are different hypotheses.
The latter lift also follows without the catalogue: for the orthogonal
E6-root/omega7 basis, the invariant quadratic is F_E6+(2/3)omega7^2.
The even diagonal of the integer matrix 3C_E6^(-1) makes -(3/2)F_E6 an
integer root polynomial, and the classical rational Borel relation maps
it to h2. The structural checker verifies this calculation separately.

## Independent implementations and their shared inputs

The [primary experiment](experiments/e7_two_primary_basis.py) reconstructs
Newton coefficient formulas, certifies both complete sources, and uses
Singular normal forms in the published rational ring. It checks 28
degrees of each source and 392 candidate-module products. Its
[report](experiments/e7-two-primary-basis.json) binds the unchanged
historical source proposals by SHA256.

The [independent replay](experiments/audit_e7_two_primary_basis.py)
imports neither the main experiment nor sage_borel. Closed Girard
partition sums and binomial expansion of the 27 linear forms replace
the Newton construction. SymPy grlex Groebner bases and a separately
written sparse reducer replace Singular. Independent Python bit
elimination replaces the primary rank algorithm modulo 2. Every
degree of both coefficient lattices and every image column is checked,
including the explicit source lifts. The lattice-splitting determinants
are reconstructed from a Cartan matrix. The
[audit report](experiments/e7-two-primary-basis-audit.json) has 28 passed
degree records, each checking both full sources.

These two implementations share Sage exact polynomial and matrix
arithmetic for the large coefficient products. They also share the
published mathematical inputs and the proposed historical source bases.
The proposals are independently certified by coefficient membership and
primitivity, rather than presumed complete. This is independence of
construction and algorithms, not independence of every arithmetic backend.

The [Schubert bridge](experiments/e7-two-primary-schubert-check.json)
checks a different mathematical route. Fresh BGG divided differences
expand all seven actual catalogue polynomials. The runtime catalogue
contains input definitions, not stored image coefficients. Positive-root
Chevalley multiplication, starting only
from the three published Schubert anchors, constructs every proposed
basis class. All 28 basis degrees agree locally with the historical
embedded image lattices. It does not independently rebuild all source
coefficient lattices; that is supplied by the first two checks.

The [plain Python/SymPy check](experiments/e7-two-primary-ring-check.json)
uses no Sage or project imports. It verifies the six ring identities,
the Poincare series and multiplicity polynomial, both permutation-basis
determinants and reflection actions, both index-3 complements, and the
integer root quadratic supported only on the E6 simple roots.
It also identifies the quadratic with the classical invariant form;
the refinement exponent proof therefore does not rely on the finite
high-degree source computations.

## Published results and repository deductions

- [Duan--Zhao, Theorem 5](https://arxiv.org/abs/math/0511332): the ambient
  integral ring and Schubert generators. The elementary reductions in
  this note are applications of that presentation.
- [Kumar--Xie, Part III, Theorem 8.2(g)](https://kumar.math.unc.edu/papers/kumar86.pdf):
  rational Levi generators and the exceptional map calculation, with the
  normalization retained here. The earlier source theorem transports
  the integral character forms through the normalized Springer maps.
- [Totaro, Theorem 6.1](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf):
  adjoint E7 torsion index 24. The general relative exponent bound 8
  follows from this, the classical divided-difference lifts, and the
  fact that the root/weight quotient is 2-primary. That bound is known.
- [Kaji--Nakagawa, Introduction and Theorem 1.1](https://arxiv.org/abs/0709.3702)
  gives related integral Schubert/characteristic-map computations and
  Chow rings of the groups. Their quotient by the ideal of divisor
  classes is a different object from our additive relative lattice
  quotient N/M; its presentation does not immediately state (1)--(3).
- The adjoining-generator lemma is elementary polynomial reduction.
  The repository deductions here are its application yielding exponent
  2 on all refinements of P7, the explicit simultaneous bases for P7,
  and the resulting closed multiplicity polynomial. The old P7 factor
  lists already contained these numbers; the new work explains them.

Targeted searches covered E7 characteristic homomorphisms, integral
root/weight images, E7/E6 parabolic cohomology and torsion-index
comparisons. The explicit source and cohomology statements cited above
were checked; related Kaji--Nakagawa characteristic-map calculations
were inspected for the distinction of objects. These checks did not
locate the displayed simultaneous basis comparison. That is a limited
literature audit, not an originality claim or evidence that the cited
authors posed this exact question.

## Remaining work

The P7 relative multiplicity formula is complete. The structural bound
does not itself give multiplicities for its finer flags. The
[refinement extension](E7_REFINEMENT_MULTIPLICITIES.md) separately proves
a uniform formula through degree nine. The subsequent A6 theorem linked
below completes 61 cases. The [D5](E7_D5_REFINEMENT.md) and
[D4 completions](E7_D4_REFINEMENT.md) now settle the other three. The other retained
sets containing node 7 are outside the complement argument as stated.
Ambient weight-image cokernels H/N, general transfer corrections at
torsion primes, and arbitrary representation direct sums remain separate
questions. Runtime and the installed 0.20.0 wheel are unchanged.

The subsequent [A6 intermediary proof](E7_FULL_FLAG_COMPARISON.md)
settles the full flag and 61/64 refinements in all degrees. Its finite
image certificate and separate audits have their own stated scope.
