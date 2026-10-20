# Proof and literature status of the integral image notes

Audit date: 2026-10-19. Scope: the type B/D Smith notes, the subsequent
type C joint-image argument, and their
normalization, lattice, geometric, and symmetric-function dependencies.
This is a repository proof/source audit, not external peer review or a
formal proof certificate. It does not certify originality. It does not
re-audit every article or the earlier dissertation proof in this repository.

## How to read the claims

**Published input** means the cited source states the indicated result.
**Audited deduction** means the argument in these notes has been checked
against its hypotheses and published inputs in this audit. It is a
mathematical derivation being presented here, not a claim that a source
already states the exact Smith formula, or that the derivation is new.
**Conditional lemma** records its hypotheses explicitly.
**Conjecture** has no complete proof here.
**Finite verification** establishes only the exact recorded cases.

The previous wording did not make these distinctions sufficiently visible.
Absence of external review is not a mathematical counterexample; conversely,
agreement of computations or multiple assistant reviews is not a proof.
The conclusions below are based on the displayed arguments and checked
source statements. External expert review and any claim of priority remain
separate work.

## General type-C joint images

The [joint-image proof](TYPE_C_JOINT_PRIMITIVITY.md) proves, for every
a,b,k>=1 and the fixed normalized inputs c_j=e_j(M), q_i=e_i(T^2), that:

- the joint image is primitive in integral symplectic flag cohomology;
- the sum of the two separate images is primitive in the joint image
  and in ambient cohomology;
- their intersection consists exactly of the constants.

The proof identifies the joint image with pullback from
Sp(n)/(Sp(a) x U(b) x Sp(k)). Integral Lagrangian Schubert lifts give
Leray–Hirsch decompositions. Reordering the quaternionic blocks as M,A,T
identifies the two squared-factor subrings with standard Grassmannian
pullbacks in a complex type A two-step flag ring. They occupy subsets of
one integral Schubert basis. The same middle Lagrangian lifts then give
a simultaneous basis for the actual separate images. The proof includes
an integral kernel presentation and Gaussian rank formulas.

Borel's integral flag presentations, the Pragacz–Ratajski Lagrangian
basis and quaternionic description, and compatibility of Schubert classes
with projections are published inputs. Their application to these selected
joint Springer/Borel images is a deduction written and checked here.
The note's audit table records the normalization, integral generation,
compatible lifts, and block-order issues. Its
[finite verifier](experiments/type_c_joint_primitivity.py) constructs the
asserted bases and checks their embedded lattices. No claim of novelty,
external review, or extension to scaled characters or type B/D joint images
is made.

A subsequent [focused proof and literature audit](TYPE_C_JOINT_AUDIT.md)
found no mathematical gap under these hypotheses. It classifies the result
as a corollary of classical integral theory applied to the selected
Springer/Borel generators, rather than an established new theorem.
Borel's Proposition 30.2 independently supplies the intermediate ring.
The exact combined statement was not located in the searched literature;
no priority conclusion follows. The audit clarifies the global fiber lifts,
compact homogeneous-space maps, integral elimination and necessary
normalization, without changing the formulas. A separate formal-polynomial
check tests the integral base-change basis without importing the library.

## General type-B comparison completed in these notes

The subsequent [general proof](TYPE_B_GENERAL_SMITH.md) supplies the
formerly missing comparison in every rank. A balanced incidence map is a
Chow envelope over function fields. It imports unbalanced-incidence rows
integrally; their Jacobi--Trudi leading terms give a unit triangular
minor. A beta-number count matches the defect-locus Betti bound, forcing
all prefixes primitive with the required ranks. This is a complete
repository derivation under the stated normalization and earlier prefix
reduction, not a published or externally peer-reviewed theorem.

A [fresh full-proof audit](TYPE_B_GENERAL_AUDIT.md) rechecked the complete
dependency chain at baseline `5e71833`, including the earlier prefix
reduction. No gap or mathematical correction was found. The
[updated Shrawan PDF](../articles/orthogonal-images-summary/README.md)
now includes the general result and its integral argument, with the
remaining type-D questions still distinguished.
The fresh audit also preserves an independent pure-Python larger-case
check: 16,200,912 corner comparisons and 88 degreewise polynomial
identities, with exact scope and hashes in its
[report](experiments/type_b_general_audit_verification.json).

The polynomial expression also has an existing Demazure-character
interpretation: Chari--Schneider--Shereen--Wand, Theorem 3.3(ii)--(iii),
[primary text](https://arxiv.org/pdf/1310.5191). That complex character
identity does not itself prove the integral Springer comparison.

## Type-D progress after the general B proof

The [type-D unit theorem](TYPE_D_FROBENIUS.md) determines every normalized
unit multiplicity, including the exact Gaussian overlap for even tails.
Its mod-2 argument retains its independent scope.

The [general square comparison](TYPE_D_SQUARE_THEOREM.md) now proves
that every normalized D square-sector coefficient-row lattice equals
its B counterpart in every rank and degree. In particular the requested
smaller-block-three case has exactly the multiplicities of \(1,2,4,8\)
displayed in equation (13) of that note. The proof uses flagged orthogonal
incidence bundles, the published equivariant eta formulas with denominator
at most \(2^t\), and an integral coefficient retraction obtained by
compensated head/tail sign changes. The retraction preserves the nonzero
base Euler class. The resulting square ring is the relative quaternionic
Grassmann ring. Existing integral B-prefix primitivity then gives
\(2^tP_t\subseteq\mathcal R_D\); summing proves the missing reverse
row containment. A fresh independent adversarial audit found no material
gap in this argument. This is an unpublished repository deduction without
an originality claim or external expert certification.

The [two-row note](TYPE_D_TWO_ROW_SMITH.md) retains an independent direct
boundary-row proof and the general exponent bound. The
[structural comparison](TYPE_D_SQUARE_COMPARISON.md) retains common
kernels, duality, finite primitive ranges, and the failed broader shortcuts.

The [general even-tail gluing theorem](TYPE_D_EVEN_GLUING_THEOREM.md)
now proves the [full-D formula](TYPE_D_CANDIDATE.md) for every tail,
closing the remaining saturated-gluing condition. Moving the specified
tail to the head gives a precise square-plus-head-Euler source. Its
relative Euler extension has a primitive selected source, including
the nonzero base Euler relation. The eta denominator bound and an
even-head parity retraction give exact D/B row-lattice equality.
Primitive actual-head B prefixes and signed level-zero overlap
identities compute that B lattice. Comparing the resulting full index
with the two sector indices proves that every gluing direction is
already a unit direction.

The [fresh integrated audit](TYPE_D_EVEN_GLUING_AUDIT.md) records the
complete dependency checks, including large actual heads and the
two-fork boundary. A misleading intermediate claim that wider relative
Schur classes vanish was corrected: only the narrower source subset
needs to map basis by basis. No unresolved gap was found in the final
argument. These are internally audited repository deductions, not
published or externally certified results.

The independent [two-head](TYPE_D_TWO_HEAD_GLUING.md) and
[two-tail](TYPE_D_TWO_TAIL_GLUING.md) proofs, and the
[full exponent bound](TYPE_D_FULL_EXPONENT_BOUND.md), are retained.
The stronger unweighted prefix strictness conditions in
[TYPE_D_GLUING_FILTRATION.md](TYPE_D_GLUING_FILTRATION.md) remain
a separate conditional route; the new proof does not assume them.

The [new full-gluing check](experiments/type_d_even_gluing_verification.json)
passes 138 full components, 16 independent full row-module comparisons,
and 27 explicit primitive unit-gluing representatives. Both historical
rank-nine timeouts now complete. The separate actual-head verifier
checks 610 exact prefix modules and 35 signed overlap identities.
These are finite corroboration, separate from the arbitrary-rank proof.

The [three-root Chern verifier](experiments/verify_type_d_three_row_smith.py)
and [report](experiments/type_d_three_row_smith_verification.json) check
actual integer coefficient modules, the four closed formulas, independent
Cartesian D expansions, and a relative quadric with nonzero base Euler
class. The separate pure-Python
[projection verifier](experiments/verify_type_d_incidence_projection.py)
and [report](experiments/type_d_incidence_projection_verification.json)
check 5,440 monomials and their Euler/square-relation multiples, retaining
the base Euler symbol and detecting the false uncompensated projection.
Finite checks corroborate the written arguments; they do not prove the
geometric or equivariant inputs. Earlier unit and square reports retain
their explicitly bounded scopes and resource exclusions.

## Current mathematical status

Throughout, $q_i=e_i(x_{\rm tail}^2)$ are the fixed normalized inputs,
$p=\prod x_{\rm tail}$ is included in D, and polynomial degree $d$
means cohomological degree $2d$. Indices are in the saturation of the
selected rational image. They need not be indices in all cohomology.

| Claim | Status after audit | Reason or limitation |
| --- | --- | --- |
| Integral theta/eta bases, finite presentations, B/C Chern-preserving comparison | Published inputs | BKT sources below; the B-to-C powers of two are existing Schubert calculus. |
| Quaternionic Grassmannian integral Schur basis and isotropic-to-quaternionic fibration | Published inputs | Pragacz--Ratajski, Appendix A Sections 10.3, 10.5; bundle version also in Panin--Walter. |
| The polynomial $U_{n,\ell}$ as a count of square-zero invariant subspaces | Classical counting, not a new polynomial | Exact specialization of Birkhoff/Hall-type counting; see below. |
| Counting subspaces with prescribed defect $\dim(U+NU)-\dim U$ | Existing general counting framework; elementary specialization here | For $N^2=0$, defect $e$ means subspace profile $(\ell,e)$. Ram's general profile-count theorem applies over finite fields; it is not a Springer-Smith theorem. |
| All-rank type-B unit Smith multiplicities equal $U_{n,\ell}$; primitivity of $P_0$ | Audited deduction | [Unit note](TYPE_B_UNIT_SMITH.md): 2-quotient rows, bounded affine Schubert union, integral Chern generation, and primitive resolution pullback. This Smith interpretation is an additional argument, not the assertion of Birkhoff's theorem. |
| Type-B rank, exponent bound $0\le e\le\min(a,k)$, block duality, stable primitivity, finite sufficient range | Audited deductions | [Structural note](TYPE_B_SMITH.md). Integrality of pullbacks and source changes is essential. |
| The general Schur-product description of each $P_t$ | Audited deduction | [Filtration note](TYPE_B_FILTRATION.md). Cell sections and both primitive pullbacks preserve the entire integer row lattice. |
| Smith multiplicities are successive ranks of $P_t$ | Lattice lemma, now applicable in every rank | [General proof](TYPE_B_GENERAL_SMITH.md) establishes the required primitivity for all prefixes; ranks alone without that hypothesis would be insufficient. |
| Every type-B multiplicity in every rank | Complete repository deduction | [General proof](TYPE_B_GENERAL_SMITH.md): Chow-envelope containment, unit triangular minor, exact beta-number count, and matching Betti bound. No external expert review or novelty claim. |
| Largest possible factor $2^\ell$ in every rank | Audited deduction | [Highest-factor note](TYPE_B_HIGHEST_SMITH.md). A Pieri coordinate-lattice proof determines this without all-prefix primitivity. |
| Three highest multiplicities in every rank; all $\ell\le4$ | Independent repository deduction | [High-layer note](TYPE_B_HIGH_LAYERS.md): complete-symmetric spectators, a coordinate prefix and graph components. This gives an alternative proof of those layers. |
| Full positive type-B partition formula for arbitrary $\ell,e$ | Complete repository deduction; formerly conjectural | [Formula](TYPE_B_CLOSED_CONJECTURE.md) and [integral proof](TYPE_B_GENERAL_SMITH.md). Defect-stratum counting is combined with the integral triangular squeeze. |
| Type-D stable primitivity, sufficient finite range, only powers of two, Euler-sector B comparison, quadric table | Audited deductions | [Type-D note](TYPE_D_SMITH.md), with integral eta specialization, localized Reynolds projection, and the explicit BKT difference basis. |
| Every normalized type-D unit multiplicity | Complete repository deduction | [Unit proof](TYPE_D_FROBENIUS.md): mod-2 Grassmann embedding, exact Frobenius/Euler intersection, derivation, core-one SXP incidence and affine Chern generation. Its unit-layer proof is independent of the later full gluing theorem. |
| Type-D square common kernel/rank, block duality, improved primitive range and one-sided B bound | Repository deductions | [Square comparison](TYPE_D_SQUARE_COMPARISON.md). General square-sector equality is not inferred from rational rank or stable presentations. |
| Type-D square factors for $\min(a,k)\le2$; general square exponent bound | Complete repository deductions | [Two-row proof](TYPE_D_TWO_ROW_SMITH.md): weighted boundary-row elimination proves equality with B; an integral parity retraction of the Chern subring gives $e\le\min(a,k)$. The subsequent general incidence theorem covers all block sizes. |
| General D square-sector B comparison, including all smaller-block-three multiplicities | Complete repository deduction | [Incidence proof](TYPE_D_SQUARE_THEOREM.md): equivariant eta denominator bound, compensated-parity retraction, quaternionic model, and primitive B prefixes give exact row-lattice equality. |
| Explicit full-D Gaussian formula | Complete internally audited repository deduction in every rank | [Formula](TYPE_D_CANDIDATE.md) and [general proof](TYPE_D_EVEN_GLUING_THEOREM.md). The relative Euler extension, primitive actual-head prefixes, and integral row containment close the even-tail step. |
| Type-D two-flag gluing rule | Audited elementary lattice lemma | [Orthogonal note](ORTHOGONAL_SMITH.md). It needs sector exponents and flag-intersection data; it is not a closed partition classification. |
| Full even-tail D partition formula and unit-gluing exhaustion | Complete repository deduction | [General theorem and index argument](TYPE_D_EVEN_GLUING_THEOREM.md): all saturated gluing is at the two unit levels; no hidden positive-level correction. |
| Stronger unweighted typed-prefix strictness | Conditional route; finite checks only in general | [Filtered lattice criterion](TYPE_D_GLUING_FILTRATION.md). The full theorem does not rely on this stronger intermediate assertion. |
| Literal defining/vector SO/Sp lattice is scaled by $2^d$ in degree $d$ | Audited normalization deduction | [Normalization audit](INTEGRAL_LITERATURE.md). Uses the stated character lattice; it does not automatically apply to Spin. |
| 174 Smith lists, 430 prefix comparisons, 516 exact row-module equalities | Finite computational evidence | Exact scopes and input hashes are in the linked [candidate](experiments/type_b_closed_candidate_verification.json) and [filtration](experiments/type_b_filtration_verification.json) reports. None is a proof for arbitrary rank. |

No fatal mathematical defect was found in the retained derivations in this
audit. This is a report of the audit's conclusion, not a guarantee against
error. The repaired citations, explicit geometric maps, and dependence
information are intended to make those derivations easier to check.

Concrete repairs include explicit cell sections and the Sato reference
class in the unit proof, the geometric sections of the localized
characteristic-map citation in D, a source-rank notation collision in the
four-row lattice split, and an obsolete two-row scope sentence.
The factor-eight argument now states its narrower dependency on $P_2$
and the largest-factor deduction, rather than all prefix results.

The new [general verification report](experiments/type_b_general_smith_verification.json)
records 585,243 direct beta-number/corner comparisons in 72 parameter cases,
and 21 fresh exact integer LR module comparisons with the constructed
unit triangular rows, through smaller-block size eight. The [high-layer report](experiments/type_b_high_layers_verification.json)
records 495 full-graph degree components and 56 exact LR module equalities.
These finite checks support the written proofs and are not proofs by themselves.

## The unit polynomial is already classical

Let $n=2r+\epsilon$, $c=r+\epsilon$, $\epsilon\in\{0,1\}$.
In Birkhoff's subgroup-count formula use the partition types

$$
\lambda=(2^r,1^\epsilon),\qquad
\mu=(2^j,1^{\ell-2j}),\qquad
\lambda'=(c,r),\quad\mu'=(\ell-j,j).
$$

The standard formula

$$
\prod_i q^{\mu'_{i+1}(\lambda'_i-\mu'_i)}
{\lambda'_i-\mu'_{i+1}\brack\mu'_i-\mu'_{i+1}}_q
$$

specializes exactly to

$$
q^{j(c-\ell+j)}
{r\brack j}_q{c-j\brack\ell-2j}_q.
$$

Summing over $j$ gives $U_{n,\ell}$. Thus the counting polynomial
belongs to existing theory. Its proposed use for these particular Smith
units requires the separate integral comparison in the unit note.

The original reference is
[Garrett Birkhoff, *Subgroups of Abelian Groups*, Proc. London Math. Soc. 38 (1935), 385--401](https://doi.org/10.1112/plms/s2-38.1.385).
Its metadata was checked; the full original was not accessible in this
audit. The exact general formula was inspected in the
[official Sage documentation](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/q_analogues.html#sage.combinat.q_analogues.q_subgroups_of_abelian_group),
which cites Butler (1987), and in
[Joshua Maglione's original alternate proof, Theorem 1](https://joshmaglione.com/2023/12/17/Birkhoff.html).
These access distinctions matter: the audit does not pretend to have read
an inaccessible original.

This executable check compares the specialization with Sage's independently
implemented Birkhoff and Delsarte algorithms. It checks an attribution
identity, not the additional Springer-Smith comparison.

```python
from sage.all import PolynomialRing, ZZ
from sage.combinat.q_analogues import q_binomial, q_subgroups_of_abelian_group

R = PolynomialRing(ZZ, "q")
q = R.gen()
checked_types = 0
for n in range(2, 17):
    r, c = n // 2, (n + 1) // 2
    ambient_type = [2] * r + [1] * (n % 2)
    for ell in range(1, min(r, 6) + 1):
        for j in range(ell // 2 + 1):
            subgroup_type = [2] * j + [1] * (ell - 2 * j)
            term = (q ** (j * (c - ell + j))
                    * q_binomial(r, j, q) * q_binomial(c - j, ell - 2*j, q))
            for algorithm in ("birkhoff", "delsarte"):
                assert term == q_subgroups_of_abelian_group(
                    ambient_type, subgroup_type, q, algorithm=algorithm)
            checked_types += 1
assert checked_types > 100
```

## Sources and their exact role

* [Buch--Kresch--Tamvakis, isotropic Giambelli](https://math.umd.edu/~harryt/papers/isogiamfinal.pdf),
  Definition 1/Theorem 1, Sections 2.1, 6.1: integral bases, raising operators,
  the Chern-preserving B/C comparison, and signed Schubert conventions.
* [Buch--Kresch--Tamvakis, even orthogonal Giambelli](https://math.umd.edu/~harryt/papers/oggiam.pdf),
  Proposition 2, Theorem 4, Section 5.2: the explicit difference-basis comparison
  and integral eta basis/specialization. The integral upgrade used in the
  Euler-sector deduction is explained in the type-D note.
* [Pragacz--Ratajski](https://www.impan.pl/~pragacz/download/comp3868.pdf),
  Appendix A Sections 10.3, 10.5: topological quaternionic Schur calculus and
  the fibration. [Panin--Walter](https://arxiv.org/pdf/1011.0649),
  Theorems 11.2--11.3, supplies the corresponding bundle presentation.
* [Wildon](https://arxiv.org/pdf/1508.07030), equation (1): the
  Littlewood/SXP input for squared-root Schur coefficients.
* [Anderson](https://doi.org/10.4153/S0008439524000122),
  main theorem p.728 and ordinary specialization pp.729--730:
  integral Chern generation on the affine Grassmannian. Its application
  to the bounded locus and its Smith interpretation are deductions here.
* [Hahn--Nebe--Stanojkovski--Sturmfels](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/shuffle0702final.pdf),
  Section 4, Proposition 14: rectangular nilpotent invariant-subspace varieties
  already have a bounded affine-lattice model. This covers the even
  case directly; the unit note spells out its odd-case Iwahori bounds.
* [Griffin](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/springer-varieties-and-halllittlewood-polynomials/64AC25F936FDD32A3A955031B880E830),
  Section 2.4, equation (2.8): preserved partial flags are Steinberg varieties;
  the Young-invariant comparison displayed there is over $\mathbb Q$.
  It is not a substitute for the integral argument here. Spaltenstein
  flags annihilated on successive quotients have different conditions.
* [Ram, *Subspace profiles, q-Whittaker functions and Krylov methods*](https://www.mat.univie.ac.at/~slc/wpapers/FPSAC2025/115.pdf),
  Definition 1.1 and Theorem 5.2 (also stated as Theorem 1.3):
  an existing general finite-field formula for subspaces of a given
  operator profile. Our defect-$e$ count for $N^2=0$ is its
  $(\ell,e)$ case. This is a comparison of the counted objects; the
  symmetric-function formula has not been transformed term by term to
  the displayed Gaussian sum here. It supplies no integral Smith bridge.
* [Karpenko](https://sites.ualberta.ca/~karpenko/publ/gqf3.pdf),
  Theorem 2.1 and Remark 2.6, and
  [Karpenko--Merkurjev](https://sites.ualberta.ca/~karpenko/publ/duality14.pdf),
  Theorem 1.2: existing integral presentations and duality for the full
  tautological Chern subring. That ring is larger than the selected
  square-only subring; these results do not directly assert its Smith
  multiplicities or saturation.
* [Calmes--Petrov--Zainoulline](https://arxiv.org/pdf/0905.1341),
  Theorems 6.4/6.9 and Sections 13.3, 13.6--13.7, with
  [Totaro](https://www.math.ucla.edu/~totaro/papers/public_html/spin.pdf),
  Theorem 3.2: the localized characteristic-map input used in D.
* [Dimitrov--Roth](https://arxiv.org/pdf/1707.06840), Section 2.8, equation
  (2.8.1): the classical integral quadric ring underlying the small-head
  corollaries.

Parts I--III supply the Springer maps, chosen forms, finite generator
images, and stable injectivity discussed with locators in
[INTEGRAL_LITERATURE.md](INTEGRAL_LITERATURE.md). Their complex injectivity
statements should not be relabeled integral saturation theorems.

* [Tamvakis, double eta polynomials](https://math.umd.edu/~harryt/papers/dbleta.pdf),
  Definition 1 and equations (4)--(5), Theorem 1, Section 4.1:
  the equivariant integral Schubert basis and explicit Chern/Euler
  geometrization used to deduce the relative denominator bound. The
  flagged incidence comparison and selected Smith interpretation are
  additional arguments in the new D square note.

## Literature search limits

The search covered Springer/Borel selected images, integral theta/eta
presentations, tautological Chern subrings, quaternionic Grassmannians,
affine Grassmannian Chern generation, invariant-subspace/Jordan-quiver
Grassmannians, Hall--Littlewood/Birkhoff counting, plethysm and Frobenius
ranks, and subspace profiles/anti-invariant subspaces. Representative query phrases included “isotropic Grassmannians
Smith”, “Chern subring orthogonal Grassmannian”, “Springer Smith normal
form”, “Schur Smith plethysm”, and “Grassmannian Frobenius image
cohomology characteristic two”.

The reviewed sources supply substantial existing structure. No exact
published theorem matching the all-rank selected-image Smith
classification was located. Nor was an exact published Smith-unit
identification located in the sources read. Neither statement proves that
such a result is absent from the literature. Before describing any
deduction here as original research, a specialist literature review is
still needed.
