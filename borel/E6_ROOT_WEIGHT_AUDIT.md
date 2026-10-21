# Audit of the E6 root-to-weight comparison

Date: 20 October 2026. This is an internal stepwise proof/source audit,
with an independently written Python/SymPy calculation. It is not external
expert review, formal verification, or certification of literature priority.
The statement being audited is [the exact comparison and transfer theorem](E6_ROOT_WEIGHT_COMPARISON.md).

## Proof dependencies and checks

| Step | What was checked | Scope |
| --- | --- | --- |
| Integral image definition | Q and X are the actual E6 root and weight lattices; integral full-flag characteristic classes descend primitively to each parabolic | Full unscaled invariant sources, not selected catalogue generators |
| Primary reduction | 3X is contained in Q, so every N/M is finite and killed degreewise by 3^d | All standard E6 parabolics |
| Source splitting | The stated six vectors have root index 64, and replacing 3omega1 by omega1 gives weight index 64; retained D5 reflections are signed permutations on the common five coordinates | Exact over Z_(3); not an integral direct sum over Z |
| Invariant rings | Even signed permutations have the elementary-square and Euler invariant generators when 2 is invertible | Complete D5 invariant rings |
| Ambient ring | Duan--Zhao Theorem 4 supplies the two-generator integral Cayley-plane presentation; P6 is transported to P1 by diagram symmetry | Published input |
| Generator images | All 11 Schubert coefficients of the six invariant images were recomputed by Sage and by independent sparse binomial divided differences | Independent implementations of the same BGG formula; not independent geometric foundations |
| Weight image | h,3s,3s^2 generate exactly the weight-source image over Z_(3) | Uses unit factors 32 and 256; not an equality asserted over Z |
| Local bases | The 27 displayed root-image basis vectors are actual source monomials; two-sided local matrix witnesses and 193 generator-closure checks prove completeness | All Cayley-plane degrees, since degrees above 16 vanish |
| Smith factors | Independent determinantal-divisor valuations agree with Sage Smith reduction in all 17 degrees | Exact finite proof certificates, not extrapolation |
| Relative polynomial basis | A torsion-index polynomial for D5 gives a B0-basis b_w with exact partial_i b_w lowering identities; the same minimal-coset subset spans every parabolic invariant ring | Classical characteristic-map and nil-Hecke input, with no averaging or division by the Levi Weyl order |
| Common image basis | The relative lifts lie in the common five-coordinate sublattice and restrict to a fiber Schubert basis | Required for both images to have the same Leray--Hirsch decomposition |
| Transfer | Ambient Leray--Hirsch independence makes the quotient a direct sum of degree shifts, rather than a potentially non-split filtration | All 48 standard parabolics refining P1 or P6, all degrees |
| Good-Levi deconvolution | Common root-polynomial fiber lifts give bases over both invariant rings and over ambient parabolic cohomology; full-flag multiplicity polynomials are divisible by the Levi length polynomial | 60 proper standard parabolics whose retained sets do not contain {1,3,5,6} |
| Local isogeny criterion | C[:,I] has full column rank mod 3 exactly off that support; all proper simply connected derived Levis are of types A/D with torsion index prime to 3 | Supplies the good-Levi hypothesis, without assuming Weyl-group order invertible |
| Representation comparisons | Earlier source classification gives scalar root or weight lattices for irreducible E6 representations; min/max valuations in a common Smith basis handle shared primes | Does not classify arbitrary direct sums or the three exceptional proper parabolics |

The direct-sum assertion does not identify an arbitrary Schubert basis
of the total space with a tensor product. It constructs a particular
common cohomology basis from polynomial fiber lifts, which is enough
to compare the two embedded image lattices. This distinction is essential.

## Independent finite coverage

The [main calculation](experiments/e6-root-weight-comparison.json) rebuilds
the six invariant BGG images, all 17 graded image matrices, selected source
monomial bases and two-sided local witnesses, 193 module-closure checks,
and the 32 D5 parabolic length polynomials.

The [independent Python audit](experiments/e6-root-weight-independent.json)
imports neither Sage nor sage_borel nor the main computation. It rebuilds
the Cartan/lattice splitting, verifies all 27 minuscule coset labels by
the weight orbit, independently evaluates the 11 relevant BGG coefficients,
rebuilds all 17 quotient relation matrices and image lattices using SymPy,
and recomputes their Smith exponents from determinantal divisors. It
enumerates all 1,920 even signed permutations to check all 32 fiber
Poincare polynomials and their multiplicity convolutions. It independently
enumerates all 63 proper E6 Levi Weyl groups on the faithful 27-weight
orbit, checking the mod-3 restriction ranks and all 60 deconvolutions.

A [separate direct BGG calculation](experiments/e6-root-weight-refinements.json)
computes 22 actual image comparisons on four other spaces: full flags in
degrees 0--4, and the D4 refinement with omitted nodes {1,6} in degrees
0--6, and the retained A1+A1 flag {1,6} in degrees 0--4. The latter lies
outside the minuscule-refinement family. It also checks maximal P3,
whose retained Levi is A4+A1, in degrees 0--4, constructing the full
integral invariant sources as reflection-kernel lattices in root and
weight coordinates. It obtains columns from invariant polynomials on those spaces,
not from the base/fiber model. Their factors match the transfer theorem.
It is not a direct computation of all degrees for all 60 parabolics.

## Attribution and literature comparison

Published inputs:

- The root/weight lattices, D5 even-sign invariant presentation and
  ordinary characteristic-map construction are classical.
- [Duan--Zhao, Theorem 4](https://arxiv.org/abs/math/0511332) supplies the
  ambient integral ring. [Iliev--Manivel, Section 6 and Proposition 6.1](https://arxiv.org/abs/math/0306329)
  already compute D5-invariant images in Schubert coordinates. Our
  coordinate-specific formulas (5) should not be advertised as a newly
  discovered invariant-to-Schubert method.
- [Demazure](https://doi.org/10.1007/BF01418790) and
  [Calmes--Petrov--Zainoulline, Sections 5--6](https://arxiv.org/abs/0905.1341)
  provide the torsion-index/coinvariant/divided-difference foundations.
  The common-basis transfer is a deduction using these and Leray--Hirsch.
- [Totaro, Theorem 6.1](https://www.math.ucla.edu/~totaro/papers/public_html/torsione8.pdf)
  proves the adjoint E6 torsion index 54; the simply connected value is 6.
  Thus the full-flag top-degree ratio 9 is already known. The uniform
  all-parabolic bound 27 follows from the adjoint torsion index and
  divided differences; it is not a new torsion-index theorem.
- The interpretation as Springer images uses Rogers's normalized formula,
  Kumar's polynomial-character/Borel map, and the earlier repository
  integral source theorems. Representation scalings are 6^d and 24^d.

Work recorded in this package: a complete local module-basis proof of
the Cayley-plane relative factors, their three multiplicity polynomials,
their transfer to 48 minuscule refinements, and good-Levi deconvolution
to 60 proper standard parabolics. The base factors agree with
the earlier finite representation-comparison table; they are not newly
computed numbers. The extension here is the explicit common-basis proof
and its consequence for the larger parabolic family and shared-prime
representation comparisons.

Targeted searches covered E6 root/weight characteristic images,
Cayley-plane cohomology, adjoint torsion indices, and the published
invariant-to-Schubert calculations. The cited primary sources were read,
including Totaro's Section 6. His E6 proof uses an A5 Levi and a
Merkurjev lower bound, not the displayed all-degree comparison. No exact
prior statement of our multiplicity polynomials and 60-parabolic comparison
was located in these checks. This is a limited literature comparison,
not evidence establishing mathematical priority.

## Completion and remaining work

Sharp multiplicities for the retained sets {1,3,5,6}, {1,2,3,5,6} and
{1,3,4,5,6} are proved in the [three-case completion](E6_REMAINING_CASES.md)
and its [independent audit](E6_REMAINING_AUDIT.md). Their adjoint-form
Levis have a 3-primary isogeny obstruction, so the common-fiber-basis
argument audited above still does not apply to them. The new proof
uses full coefficient-lattice certificates and a separate P1 reduction.
Together the two proofs cover all 63 proper standard parabolics. The
trivial G/G quotient is already zero.
Ambient weight-image cokernels H/N and equivariant specialization kernels
also remain distinct problems.
