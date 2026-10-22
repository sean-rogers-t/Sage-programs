# Proof and literature status of the integral image notes

The [research ledger](RESEARCH_LEDGER.md) indexes the program's uniform
results, finite cases, conditional routes, counterexamples and open questions.
The detailed statements and source audits below retain their original scope.

Audit date: 2026-10-19. Scope: the type B/D Smith notes, the subsequent
type C joint-image argument, and their
normalization, lattice, geometric, and symmetric-function dependencies.
The exceptional section records the subsequent G2/F4 and E6/E7 calculations,
duality deduction and explicit primitive-image subalgebras.
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

## Classical equivariant API and symplectic higher Fitting ideals

The [classical backend](CLASSICAL_EQUIVARIANT.md) now evaluates defining
SL_m literal sources on every standard partial flag in equivariant
cohomology or K-theory. It uses published restriction formulas and
checks actual Laurent integrality. Full integral source generators
are enumerated only on G/B and P1; other parabolics support specified
invariant inputs. No complete type A image-module or augmentation-kernel
classification is claimed.

The [symplectic even/odd argument](SYMPLECTIC_K_FITTING.md) identifies
the entire defining Sp_(2n)/P1 integral cokernel, as a module over the
source coordinate z, with R[z]/(E(z^2+4)+zO(z^2+4),2O(z^2+4)).
Here product(s-s_i)=E(s^2)+sO(s^2). This gives an n-by-n presentation
in all ranks, a determinantal formula for every Fitting ideal, and
the closed last ideal (2) in odd rank and 2(e1,e3,...,e_(n-1))
in even rank. Explicit unimodular operations give global cyclic
decompositions for ranks one, two and three. The argument is
elementary algebra applied to the prior literal image; the ambient
K-theory, Fitting framework and source identification are inputs.
There is now a uniform stronger reduction: the cokernel is two copies
of an r-by-r block cokernel, r=floor(n/2), and an extra R/(2) if n is
odd. The block is multiplication by 2O modulo E for even rank, or
2E modulo O for odd rank. Monic division gives explicit integral
unimodular equivalences; parity and translation give the two copies.
This is an R-module splitting, not a source-coordinate module splitting.
It explains paired Smith factors under integer and DVR specializations.
For rank four, J=(e1,e3), F=product(s_i+s_j), the full Fitting list is
(16F^2), 8FJ, 4J^2, 2J, then R. Higher closed formulas were unresolved
at that stage; the subsequent Schur-minor result below advances them. The classical
Orlando resultant identity is a known input, not a claimed new result;
the specific image-module application has no priority determination.
No priority or external-review claim is made.

The software checks include BGG, Newton projective and line-bundle
restrictions, an explicit nonreduced-subword sum, and comparisons
against the original 2n projective inclusion at regular and singular
points through rank eight, plus rank-nine ordinary augmentation.
Independent coefficient variables give exact matrix certificates through
rank eight and all rank-four minor-ideal identities. The proof is
uniform; these finite checks are corroboration.

The subsequent [Schur-minor argument](SYMPLECTIC_K_MINORS.md) resolves a
uniform generator problem: for multiplication by q on A[u]/p with p
monic, the first j columns generate the entire ideal of j-minors.
This holds over arbitrary commutative A, by the published integral
Schur cyclicity of exterior powers (Laksov-Thorup; Gatto-Santiago).
It is an application of known Schubert algebra, with a full integral
argument supplied, rather than a claim of a new exterior-power theory.
The ideals are the content of a universal norm class in Schur coordinates.
For the defining symplectic image this gives all-rank partition-indexed
Fitting generators, plus complete explicit lists through ranks five,
six and seven. The cubic lists use a three-generator two-minor ideal
and adj(adj(H))=det(H)H to prove the required product containments.
There are 94 exact polynomial matrix identities and 72 independent
original-inclusion Fitting checks. Further simplification at block size
four (n=8,9), primary arithmetic structure and full ordinary Smith
multiplicities remain. Priority of the specific application is not established.

## Other cohomology theories and affine groups: 22 October exploration

The [exploration report](OTHER_THEORIES_PROGRAM.md) records a new direction;
Spin parity is shelved at the user's request. The
[finite comparison](OTHER_THEORIES_FINITE.md) proves a common integral
K/Borel kernel and a strict support filtration for literal sources, then
constructs a specified oriented trace-coordinate lift. The universal MU
inclusion is triangular with the old integer determinant. This controls
prime support, while higher Fitting ideals and modular Tor kernels can
change. No gamma-filtration identification or arbitrary normalized-source
integrality is assumed. The [SL3 theorem](OTHER_THEORIES_SL3.md) gives
complete factors for every representation and parabolic from its source
scale and cubic correction; independent actual-character calculations
check 41 representations on all three proper flags.

The [type A coordinate comparison](OTHER_THEORIES_TYPE_A.md),
[symplectic height/K recurrence](OTHER_THEORIES_HEIGHT.md),
[quantum projective lifts](OTHER_THEORIES_QUANTUM.md) and
[affine A/C contents](OTHER_THEORIES_AFFINE.md) give written uniform
deductions with exact coefficient/source hypotheses. Formal group height,
the ambient quantum rings and Hopf models are known inputs. Independent
Euler, elliptic-curve and affine Chevalley calculations corroborate the
new image formulas. The [equivariant K note](OTHER_THEORIES_EQUIVARIANT_K.md)
uses deformed rational source freeness and reflection/localization
determinants for every parabolic; its local Smith patterns state their
regularity and DVR hypotheses. An integral Sp/P1 family is proved
separately, without promoting rational determinant formulas to general
integral source freeness. The [folding note](OTHER_THEORIES_FOLDING.md)
proves compatible trace projections under stable representations while
retaining exact integral source defects.

All are repository deductions and research prototypes. The cited ambient
characteristic-map, primitive K, oriented, quantum and affine theories are
published; the K-theoretic Peterson theorem is established. Exact
literature priority for these specified arithmetic/image formulations
has not been established. No new production runtime API, full exceptional
K Smith tables or central/rotation-compatible affine Springer construction
is asserted. The pointwise untwisted/twisted loop maps already follow by
base change; an explicit constant-loop example blocks naive flag descent.

## General central-parity kernels: positive families and counterexample

The [general kernel note](CENTRAL_PARITY_KERNEL.md) proves a sufficient
criterion using a torsion-free integral Borel presentation, a free actual
Levi source decomposition containing 1, and compatibility of the source
generators and relations with differentiation modulo two. This is a
sufficient criterion, not a necessary-and-sufficient root-data classification.

For the actual Cartesian weight lattice of Sp_(2n), symmetric-square
relations and ordinary/signed Artin bases prove the criterion for every
rank and every standard parabolic, in every degree. An odd-block type A
application works for SL_(2r) with arbitrary r, extending the scope of
kernel stability beyond the earlier odd-r quadratic-image theorem. The
kernel assertions alone do not extend the root-image Smith formula.
The positive proofs are elementary deductions from published integral
presentations, not novelty claims.

The naive index-two/fixed-odd-weight hypothesis fails on Spin7/B in degree
four. An explicit quartic is in the actual integral kernel, but its
derivative is outside the reduced primitive cubic kernel. A fresh
integer BGG calculation and an independent polynomial identity/monic
division proof establish the counterexample. Ambient mod-two cohomology
would miss it. No unit quadratic lift is available in this example.

The [odd spin extension](SPIN_PARITY_EXTENSION.md), audited 22 October
2026, disproves uniform failure on Spin(2n+1)/B. An Artin-basis base change
reduces the entire kernel in every rank to the auxiliary symmetric source
R[z,c2,...,cn], with c1=2z. Its primitive quartic relation proves that
degree-four failure occurs exactly at rank 3. A closed four-relation
presentation proves all-degree stability on Spin9/B. Independent rational
membership, modular Groebner and finite-flatness certificates prove complete
stability on Spin11/B and Spin13/B, and first failures on Spin15/B and
Spin17/B in polynomial degrees 11 and 16. These exhaust the stated spaces;
they do not classify arbitrary rank or assert novelty. The rational Borel
presentation, torsion-free Schubert basis and symmetric polynomial freeness
are known inputs; the derivative conclusions are deductions here.

The [recursive follow-up](SPIN_PARITY_RECURSION.md) proves an arbitrary-rank
image-lattice recurrence and exact derivative factorization criterion using
the classical Schur-Q root lattice and one-box Pieri rule. A lowering
operator satisfying [L,U]=(n(n+1)/2-2d) proves all-degree stability whenever
n=1,2 modulo 4; its preservation of the actual integral image lattice is
proved, rather than inferred from rational Lefschetz theory. An infinite
Schur-P derivation proves stability through degree n for every n. Native
Hermite/bit certificates check the recurrence through rank 12; a separate
quadratic straightening calculation independently gives stability in rank
11 and first failure at degree 34 in rank 12. The closed classification
for n=0,3 modulo 4 and a closed first-failure formula remain unresolved.
The Schur-P/Q presentation, Pieri rule and raising/lowering framework are
known; the central-parity conclusions are deductions, without a novelty
claim. No runtime API is changed by these proof and audit files.

This general question was the earlier structural priority and is now
shelved at the user's request. The existing all-degree E7 image theorem
is unchanged.

The [adjacent-rank boundary theorem](SPIN_PARITY_BOUNDARY.md) isolates the
whole obstruction without assuming that a derivative has descended in
earlier degrees. The boundary lattice is the entire kernel of
N_(n+1) -> N_n; lowering its largest part n+1 to n gives b_n. The derivative
obstruction is exactly (b_n L+N_n)/N_n, an elementary two-group, and its
dimension is rank over F2 of 2 F^(-1)b_n L. This is a proved general
criterion and an effective first-failure formula. A closed evaluation
depending only on n, particularly for n=0,3 modulo four, remains missing.
Classical Schur-P/Q identities supply the infinite derivation used in the
proof; they do not alone evaluate the boundary-image index. The new note
does not assert a literature-priority claim or an E7 all-degree structural
proof.

The same note proves a second, elementary source-parity formula. With
C=M[z^2], N=M[z], and e_d=log2[N_d:C_d], the exact degree-d obstruction
has dimension a_(d-1)-e_d-e_(d-1). Thus the entire kernel is preserved
if and only if log2[N:C]=2^(n-1). This does not evaluate that index in
closed form. New adjacent-rank certificates, audited by exact inclusion
and mod-two row rank without HNF software, agree with an independent
replay of the original derivative criterion. They establish the new
rank-15 first failure in polynomial degree 37 (Spin31/B, cohomology 74).
A separate Pfaffian calculation verifies 94 Schur-P lowering identities.

## Structural E7 multiplicities: scope of the new argument

The [quadratic-derivation note](E7_QUADRATIC_DERIVATION.md) isolates a
general kernel-stability criterion for the specified quadratic
presentation of N/2N over the reduced actual root image. The common E7
quadratic lift alone does not imply this criterion; an explicit algebra
counterexample is included.

Fresh primitive kernel columns prove the criterion through degree 10
uniformly for all 112 retained subsets with a fixed odd weight. The
integral relation expression is unique in that range, which justifies
passing to Levi invariants. The full-flag calculation separately reaches
degree 13, including saturated half-sums. Its independent audit uses
27 E6 forms rather than the primary Weyl orbit and its own exact bit
membership arithmetic. It reads no relative multiplicity table.

The same criterion gives an arbitrary-rank theorem for SL_(2r), r odd,
with at least one odd Levi block, over Z_(2): the usual integral Chern
presentation supplies full source/image and all-degree kernel stability.
Five complete Grassmannian images are independently checked by native
Schur/Pieri multiplication closure. This is a deduction from published
type A presentations, not a priority claim. It concerns 2-primary factors.

A subsequent [degree-eleven proof](E7_DEGREE11_LIFTING.md) handles the
first rank-one saturation correction uniformly. The mod-two cubic and
quintic are obtained by formal Frobenius derivatives of the common
quadratic. Explicit coefficient identities show that the half-sum's
derivative has an integral invariant lower-kernel lift. The argument
lifts only coefficients forced by an already invariant element; it does
not assume that every modular invariant lifts. A separate cubic
nonlifting obstruction shows that the half-sum coset is absent on E6,
the two D5 Levis and D4. Common-fiber descent covers the 108 good-Levi
cases and transfers the old full-flag derivative check through degree 13.
Thus the structural theorem now reaches degree 11 on all 112, and 13
on 108. Fresh 56-weight versus 27-form constructions independently
verify the rank-one kernel, coefficient identities, four obstructions
and family inventory. No relative Smith table is an input.

A uniform all-degree E7 kernel-stability and invariant-lifting proof is
still missing. Its all-degree multiplicity theorem retains the earlier
case/transfer proofs. Constructing the derivation from that theorem's
dimension count is valid as a corollary but would not independently
explain the counts. Runtime/API version remains 0.22.0.

## All-parabolic E7 API packaging

Version 0.22 also packages the existing complete E6 classification in the
[same relative API](E6_RELATIVE_IMAGE_API.md), preserving 3,9,27 separately.
The audit reconstructs all 63 cohomology/formula records, checks fresh
E6/P1 images in degrees 0,1,2,3,5,8 and rebuilds full A5 root/weight
invariants as integer reflection kernels in degree three, recovering Z/27.
This is implementation of the earlier image theorem, not a new E6 theorem.

The [G2 full-flag note](G2_FULL_FLAG.md) adds a complete finite image
calculation: Z[X,Y]/(X^2-3XY+3Y^2,Y^6), with the full literal minimal
source X=6omega1,Y=2omega2 transported through xi_(omega1)^B. The quotient
is free of rank twelve and all twelve basis images are independent;
this proves the integral kernel in every degree. A native root/coroot
Chevalley audit independently reconstructs all minimal/adjoint image
columns and all Smith factors. The earlier source theorem is shared.
All proper standard G2 parabolics are now covered for these literal
sources, and other irreducibles follow by the proved source-scaling law.
The [0.22 release record](release-0.22-verification.json) binds the new
checks. Literature priority remains unresolved.

Update, 21 October 2026. Version 0.21 packages the existing complete
relative E7 theorem in a [compressed summary API](E7_RELATIVE_IMAGE_API.md).
It includes all retained node sets, relative Smith multiplicities and the
normalized minimal/adjoint comparison from sources 12Q and 36X.
The gcd/lcm formulas for scaled image sums and intersections follow from
ordinary simultaneous Smith normal form. This adds no new image theorem
or literature-priority claim.

The [independent API check](relative-image-verification.json) reconstructs
the formulas and verifies scaled comparisons using integer matrix kernels;
fresh E7/P7 character-image overlaps cover degrees 0,1,2,3,6,10.
Previously audited full-source theorems, the A5 intermediary polynomial
and two native base image lists are shared inputs, explicitly identified.
The [0.21 release record](release-0.21-verification.json) binds these checks
to the runtime and installed wheel. Earlier version references below
describe their historical checkpoints. Absolute image cokernels and
structural generalizations are unchanged.

## First local-image benchmark and equivariant pilot

Update, 20 October 2026. The [new proof/source audit](LOCAL_IMAGE_AUDIT.md)
checks the [local image bound and projective family](LOCAL_IMAGE_BOUND.md)
and [equivariant specialization counterexample](EQUIVARIANT_IMAGE_PILOT.md).
The normalized source lattice lies in X(T), by differentiating its integral
Laurent pullbacks at 1. Primitive Schubert pullback then proves integrality
of every full literal Levi Borel image in the source theorem's scope.

Classical torsion-index divided differences supply integral Levi-invariant
lifts of t_G times every parabolic Schubert class. Hence each degree-d
full-source cokernel is killed by t_G e_rho^d, with
e_rho=exp(X/Lambda_rho). There is no extra Levi Weyl-order factor.
This is an explicit application of published characteristic-map machinery,
not a new general good-prime theorem.

The standard free equivariant Schubert basis and triangular induction
also recover the whole equivariant image at these good primes. Each
fixed-degree specialization kernel is finite with primes dividing
t_G e_rho. This does not force faithful specialization at the bad primes.

For defining SL_n on P^(n-1), the literal factors are n/gcd(n,d).
The proof uses elementary root-difference invariants, binomial generator
images and a p-adic gcd argument. Its root-characteristic interpretation
requires further literature comparison; priority is not established.
The equivariant SL4/P1 image gives a nonzero specialization kernel:
exactly Z/2 in polynomial degree five, with no defects in smaller degrees.
The integral proof and independent Sage/SymPy rank-52 lattice calculations
agree. Higher-degree defects and a general equivariant backend remain open.
Runtime sources and the 0.20.0 wheel are unchanged.

## Sharp E6 root-to-weight comparison

Update, 20 October 2026. The [comparison theorem](E6_ROOT_WEIGHT_COMPARISON.md)
and [stepwise audit](E6_ROOT_WEIGHT_AUDIT.md), together with the
[three-case completion](E6_REMAINING_CASES.md) and its
[audit](E6_REMAINING_AUDIT.md), prove exact 3-primary factors for all 63
proper standard E6 parabolics. The Cayley-plane image
has a complete local source-monomial basis. A common relative Schubert
basis transfers its factors to finer flags; at a good adjoint Levi,
Leray--Hirsch deconvolves full-flag multiplicities by the Levi length
polynomial, covering 60 cases. The three retained sets containing
{1,3,5,6} use a different proof: an A5+A1 subsystem at 3 and full
source coefficient lattices. The two base cases have 52 independently
audited degrees in total; the third follows from a common (1,z) basis.
The [complete tables](E6_REMAINING_TABLES.md) include all three cases.

The ring presentation and characteristic-map framework are published
inputs. Totaro's adjoint torsion index already implies the uniform
exponent bound 27 and the full-flag top-degree ratio 9. The exact
all-degree multiplicity formulas and their all-parabolic application
are deductions recorded here; no literature-priority claim is made.
An independent Python/SymPy implementation reconstructs the base
lattices, BGG coefficients, Levi Weyl polynomials and deconvolutions;
additional direct BGG checks test other flags. The completion independently
reconstructs the rational quotient and every source/image column, and has
20 further direct BGG degree/case checks. Sage exact arithmetic is shared
in this replay; the structural and literal-weight checks use plain
Python/SymPy. These are computer-assisted deductions, not an ambient H/N
classification. Runtime remains 0.20.0.

## E7/P7 simultaneous bases and the refinement exponent

Update, 20 October 2026. The [2-primary comparison](E7_ROOT_WEIGHT_COMPARISON.md)
and [audit](E7_ROOT_WEIGHT_AUDIT.md) identify the full local images as
N=R[h,2s,2t] and M=R[2h,h2,2s,2t], with R=Z_(2). A common 56-element
basis consists of h^i for i<=27 and h^i(2s),h^i(2t) for i<=13;
the root image doubles exactly the elements with odd h exponent.
Thus N^d/M^d is elementary 2-primary, with multiplicity polynomial
q[14]_(q2)+(q6+q10)[7]_(q2). The old finite tables already recorded
these factors; the extension is their explicit basis proof.

The candidate-ring basis and its closure follow by explicit identities
in Duan--Zhao's published integral ring. Four root-invariant lifts give
the candidate generators. The full-image identification is a finite
computer-assisted step: both complete invariant sources are freshly
certified in all 28 degrees, retaining the rational kernel before
integrality and image formation. An independent Girard/SymPy/bit-elimination
replay checks every source and image column; Sage exact arithmetic is
shared in that replay. Fresh BGG checks verify all seven generators,
and Chevalley multiplication checks all 56 basis classes. A plain
Python/SymPy check separately verifies lattice and ring identities.

An index-3 common E6 complement at 2, with omega7 fixed and a root
quadratic lifting h2, gives N_I=M_I+hM_I and hence 2N_I contained in M_I
for all 64 retained subsets I of {1,...,6}. This is a structural
adjoining-generator deduction, not an application of a non-torsion
Levi-fiber theorem. It supplies no multiplicities on those finer flags.
Totaro's general exponent bound 8 is published background; priority
of the sharper application and explicit basis statement remains
unestablished. Runtime and the installed wheel remain 0.20.0.

## E7/P7 refinements: a uniform truncation and a source obstruction

The [refinement note](E7_REFINEMENT_MULTIPLICITIES.md) proves the multiplicity
formula [q^(d-1)] (1-q6)(1-q8)H_I(q) for every retained I in {1,...,6}
and 0<=d<=9. The 56-weight and independent 27-form constructions yield
primitive quadratic/sextic/octic kernel certificates in this range.
A direct full-source BGG/Hermite/Smith calculation on the A5 refinement
checks all ten degrees. This is a uniform bounded-degree theorem, not
a complete finer-flag multiplicity formula.

The original all-degree target reduced by classical good-Levi deconvolution
and an integral D5 symmetry to three base polynomials: full flag, D5 and
D4. The subsequent A6 theorem below completes the full-flag polynomial. A separate degree-16 Chern invariant proves that six normalized
rational-degree E6 generators do not generate the full invariant source
at 2. Both coefficient completeness and the witness have a plain-integer
audit. This does not disprove an image transfer formula. Known inputs,
precise finite scope and unresolved novelty are recorded in the note.

## E7 full flag: complete A6 intermediary proof

The [new all-degree comparison](E7_FULL_FLAG_COMPARISON.md) proves
T_full=q P_(W(E7))/(1+q). A common A6 complement at 2 makes both full
invariant sources polynomial, and 576-dimensional exact image bases
exhaust every degree of E7/P2. The independent checks include fresh
126-root/Girard relations, plain-Python Buchberger and image closure with
primitive generation, and direct full-flag saturated kernels in degrees
10..13. The latter exhibit corrections missed by unsaturated rational
relation multiples. Both complete A6 source identifications and rational
kernel completeness are justified before taking images.

Classical common Artin/fiber bases transfer the result to the full flag
and 108 proper good-Levi parabolics; the separate P7 theorem adds a 109th.
This A6 argument completes 61/64 P7 refinements; the D5 completion below
adds two, leaving D4 as the remaining all-degree target. The degree-16 E6
invariant-ring obstruction is still
valid. The new proof uses an A6 intermediary. Published inputs, shared
arithmetic backends, finite scope and unresolved novelty are stated in
the proof note. Ambient H/N and arbitrary direct sums are separate.

## E7 D5 refinements: complete full-source comparison

The [D5 proof and audit](E7_D5_REFINEMENT.md) proves
T_D5=q P_(W(E7))/((1+q)P_(W(D5))) for retained sets {2,3,4,5,6} and
{1,2,3,4,5}, with 756 factors 2 across degrees 0..43. A ten-generator
2-local invariant source is proved complete by a mod-2 graph/Hilbert
certificate and Nakayama, rather than assuming polynomiality at 2.
The rational quotient of total dimension 1512 and both complete image
lattices are exhausted. Plain-Python source, fresh-weight polynomial
identity/Jacobian, Buchberger and image-generation audits use separate
arithmetic implementations. E6 opposition transfers the result to the
second D5 case. Classical inputs and unresolved priority are identified
in the proof. The method now covers 111 proper E7 parabolics and 63/64
P7 refinements at this stage; the D4 theorem below completes the family. This does not classify
ambient H/N or direct-sum sources and adds no runtime API.

## E7 D4 completion: all 64 P7 refinements

The [D4 proof and audit](E7_D4_REFINEMENT.md) completes retained
I={2,3,4,5}: the relative factors are all 2, with multiplicity polynomial
q P_W(E7)/((1+q)P_W(D4)) and 7,560 factors across degrees 0..51.
Its ten-generator full integral source is certified by an independent
mod-2 graph/Hilbert and lattice audit. The relative source quotient over
D5 has dimension 34, so a free rank-10 invariant-source transfer fails.
Instead all 34 integral module generators enter a rational quadric tower.
The reference root and weight lattices enlarge by equal indices in every
degree, checked by PARI HNF and independently reconstructed matrices with
native unit-pivot Smith elimination. The previously audited D5 rational
quotient and full image bases are shared inputs, explicitly identified.

Together with the P7/A6/D5 proofs, this proves
T_I=q P_W(E7)/((1+q)P_WI) for all 64 I contained in {1,...,6}.
The common source split and these Hilbert dimensions also imply that
N/2N is a free quadratic algebra on (1,w) over C=image(M -> N/2N),
with w^2=beta(B_i), choosing any omitted i in {2,5,7} and
w=beta(omega_i). The three common complements have odd indices 7,15,3,
and the integral quadratic lifts are checked separately. This C is not M/2M.
The combined scope is exactly the 112 proper E7 parabolics omitting at
least one of nodes 2,5,7. The other 15 retain all three and have N^1=M^1,
so the same polynomial cannot extend unchanged to them. This is an exhaustive computer-assisted proof for the stated
finite family, not a uniform theorem for arbitrary root data. Known
Borel/invariant/Nakayama inputs and the new finite deductions are
separated; priority, ambient H/N and a structural explanation of the
index equality remain unresolved. No runtime API changes.

## E7 completion: all 127 proper standard parabolics

The [complete E7 theorem](E7_ALL_PARABOLICS.md) determines the relative
full root/weight Borel-image groups in every degree on all 127 proper
standard parabolics. Every nonunit relative factor is two. The preceding
112-case polynomial q P_W/((1+q)P_WI) applies exactly to the retained
sets omitting at least one of 2,5,7; the other fifteen have the separate
formulas tabulated in the theorem. Their first gaps occur in degree two
(twelve cases), four (two cases), or eight (D6).

The [remaining-family note](E7_REMAINING_PARABOLICS.md) proves for E7/P1
(D6 Levi) that every nonunit relative root/weight image factor is 2, with
multiplicity polynomial q8(1+q)([10]_(q2)+q4[6]_(q2)+q12), totaling 34.
Full weight-source components are certified through top degree 33, root
components through 27; in degrees 28..33 an explicit root lower image
equals the full weight image. No unbounded invariant-algebra generation
claim is used. Native Python/FLINT source reconstruction, fresh 56-weight
moment/Jacobian identities, and native Buchberger/image audits verify
the finite proof. Higher root invariants change the truncated-generator
answer, so selected rational/integral generators alone are insufficient.

Four Weyl-conjugacy words are checked independently in root-mask and
root-matrix implementations. They reduce the initially excluded fifteen
diagrams to eleven types. The [D5+A1 completion](E7_D5_A1_COMPARISON.md)
now proves the complete E7/P6 polynomial, totaling 366 factors two.
Its full 25-generator root source uses a necessary degree-23 correction
and the degree-32 root-orbit norm. The correction's Borel image is
redundant, so the earlier independently audited lower image is the full
image. The fresh 56-weight orbit, exact Jacobian and 333 unisolvent
samples independently establish the base relations. The earlier D5
source theorem is an explicitly shared input.

The [all-degree source argument](E7_LOCAL_SOURCES.md) proves both this
25-generator source and the 18-generator D4+A1 root source. It compares
an explicitly witnessed leading-monomial algebra with the rational Weyl
invariant Hilbert series as a full rational-function identity. Independent
audits verify the product witnesses, nonnegative monoid expressions,
closure and monomial-ideal Hilbert recursion. It therefore proves
generation in all degrees, rather than extrapolating from a cutoff.

The [A5-family proof](E7_A5_FAMILY.md) completes twelve diagrams. An
odd-index A5+A2 intermediary has 293 factors two; common integral A2
source bases give two transfers. A uniform involution lemma supplies
the full even-block root sources, with 51- and 342-generator modules.
Separate rational reducers and mod-four row/column eliminations agree
on every image column through top degrees 53 and 57. Their complete
formulas have 4,911 and 29,982 factors and a common Gaussian-multinomial
expression with a halved-block correction. The expression is an exact
finite E7 deduction; a geometric folding explanation or an arbitrary-rank
generalization is not proved.

The [D4+A1 completion](E7_D4_A1_COMPARISON.md) uses the rank-ten quadric
tower over P6 and the full 34-generator weight-source module. A native
FLINT audit reconstructs all source images and checks full local lattices
in every degree through 50. It gives 3,618 factors two. The supplementary
finer-flag intersection bound can be strict; invariants after Borel are
not used as a full-source substitute.

The rational cohomology of E7/P1 is already known from Chaput--Perrin,
Proposition 5.6; the Borel presentation, signed-permutation invariants,
Weyl conjugacy and lattice primitivity criteria are classical inputs.
The relative Smith calculations and source-completeness deductions are
what is established here. Earlier targeted literature checks are retained
in the linked notes; the completion does not assert a new comprehensive
search or establish priority. Structural explanations beyond this finite
classification and ambient H/N remain separate. No runtime or public API
changes were needed for this completion.

## Exceptional integral images: first complete finite tables

The [Springer geometry extension](SPRINGER_GEOMETRY.md) concerns a separate
characteristic-zero problem. It records complete G2 minimal/adjoint quotient
maps, generic degrees 6 and 18, and the full left-translated critical
determinants. The degrees have Newton-polygon proofs and exact generic-fiber
checks; the full determinants agree coefficientwise with independent rational
representation-matrix constructions. A general root-factor formula follows
from equivariance. A standard tensor-Casimir identity gives the quadratic
coordinate uniformly, with E6/E7 minimal/adjoint checks by both Laurent and
tensor methods. The [E6 extension](E6_GEOMETRY.md) now supplies the full
basic-invariant list, all six minimal/adjoint quotient coordinates and both
full determinant products. Dominance bounds and exact unisolvent matrices
prove the large character formulas; independent 27- and 78-dimensional
matrices check all 72 root directions coefficientwise in both cases.
The [E7 extension](E7_GEOMETRY.md) supplies all seven coordinates through
degree 18 and both full 133-dimensional determinant products. Its dominance
bounds and exact unisolvent checks certify all coefficients; independent
56/133-dimensional representations check all Cartan entries and 126 root
directions for both maps. The 912-dimensional character is reconstructed
separately by root strings and Freudenthal. The [generic-degree deduction](EXCEPTIONAL_DEGREES.md)
now gives E6 degrees 1080/2754 and E7 degrees 7272/33664 for minimal/adjoint
maps on simply connected groups. It uses explicit toric boundary schemes
and their localized intersection contributions. Separate specialized-fiber
lengths provide finite corroboration, not the characteristic-zero proof.
These are exact deductions and computations, without a
novelty claim. They do not identify the minimal Cayley degree of a group.

The [irreducible lattice classification](IRREDUCIBLE_CHARACTER_LATTICES.md)
proves D_lambda=Q+Z*lambda by multiplicity-one simple-root paths, hence
Lambda_lambda=I_lambda*(Q+Z*lambda)^*. This identifies every exceptional
irreducible torus source: index times Q_long for G2/F4; index times P or Q
for E6/E7 according to whether the central character is trivial; index
times P=Q for E8. The index normalization is the published Dynkin
convention I_ad=2h^vee. The highest-weight and trace formulas are published
inputs; their source-lattice interpretation is the deduction recorded
here, without a priority claim. The fast constructor enumerates only
roots, not representation weights. A reducible F4 control shows that
index and central character alone do not classify direct sums. This
does not determine arbitrary integral Levi presentations or Borel images.
The [detailed classification audit](IRREDUCIBLE_LATTICE_AUDIT.md) found
no gap and separates the known ordinary weight-span result from the
multiplicity-weighted argument and integral Springer interpretation.
The exact formulation was not located in the reviewed literature;
priority remains unresolved. Its independent boundary checks test every
nonempty support mask on all five exceptional diagrams and explicit
controls for path selection, dual normalization and reducible inputs.

The [general trace-source theorem](TRACE_SOURCE_THEOREM.md) now proves
C_rho=Sym_Z(Lambda_rho) for every supported normalized representation map,
including direct sums with faithful differential. Saturation gives full
coefficient rank modulo every prime; the trace-moment factorization makes
the adapted pivot exponents independent over Q. Leading terms and the
least-denominator argument then prove the full, inhomogeneous intersection.
This resolves the previous representation-by-representation degree-one
generation question under these hypotheses. It does not assert that
Levi-invariant rings are generated in degree one. The saved survey and
Python-only weight/trace/minor audit check 59 representation cases,
354 prime-adapted leading bases and 18 bounded searches. The proof is
all-degree and all-prime; the numerical checks are finite. No originality
claim is made. The [2026-10-20 detailed proof/source audit](TRACE_SOURCE_AUDIT.md)
found no gap under the stated hypotheses and distinguishes coefficient
minors modulo p from exponent independence over Q. It did not locate the
exact integral equality in the reviewed literature; priority remains
unresolved.

The [E7 adjoint proof](E7_REPRESENTATION_COMPARISON.md) identifies
Sym_Z(36P) by a unit coefficient minor and a determinant-one row change
whose seven leading coefficients are all one and whose exponent
determinant is two. The certificate works in every characteristic.
Since P/Q has order 2, the coprime-index argument proves
J56^d intersect J133^d = 3^d J56^d and their sum = 3^(-d) J133^d,
for every standard Levi and degree. These are deductions without a
novelty claim. The divisor images 24 Z h and 36 Z h prove that the
full E7/P7 image rings are incomparable; some individual degrees have
adjoint-in-minimal containment. Specifically, containment is strict in
degrees 2, 4, 24 and 26; all other positive degrees through 27 are
incomparable. The quotient of their sum by the minimal image is
elementary abelian of rank zero, one or two, as listed in the comparison.

The complete E7/P7 comparison computes all 28 adjoint source lattices,
certifies their completeness, and freshly checks all 1,967 rational
Borel columns by the separate ordinary-ring/Chevalley implementation
with its documented BGG anchor dependencies. A Python-only audit
reconstructs the integral weight-lattice coordinates, compressed
polynomials by a closed partition formula, every source coefficient
constraint, eight direct Laurent intersections and image arithmetic.
The prior minimal tables are reused with six fresh minimal API checks.

The [E6 adjoint proof](E6_REPRESENTATION_COMPARISON.md) identifies the
full torus source as Sym_Z(24P), while the minimal source is Sym_Z(6Q).
A unit coefficient minor proves the linear lattice, and a unimodular
change with six independent unit-leading monomials proves the full source
in every characteristic. Since P/Q has order 3, a coprime-index argument
proves J27^d intersect J78^d = 4^d J27^d and their sum = 4^(-d) J78^d
for every standard Levi. These are checked deductions, with no novelty claim.
Complete E6/P1 tables show incomparability in all 16 positive degrees;
the divisor example alone proves incomparability of the full image rings.

The [E6/E7 public APIs](MINUSCULE_CHARACTER_API.md) reproduce all 45
historically audited minimal degree rows and compute all 17 adjoint rows.
All 2,731 rational image columns are checked by the separate
published-ring/Chevalley implementation, retaining its documented BGG
anchor dependencies. A standard-library audit reconstructs every adjoint
source coefficient lattice and eight direct Laurent certificates, and
checks all image/comparison arithmetic. The finite checks are distinguished
from the all-degree proof. The API is limited to E6/P1 and E7/P7; it does
not classify arbitrary representations or supply finite integral source
presentations.

The [general character-source criterion](CHARACTER_SOURCES.md) proves that
the entire integral torus source equals Sym_Z(Lambda_rho) exactly when
the pullbacks of a basis of its integral linear lattice are algebraically
independent modulo every prime. This is an elementary saturation argument,
not a claimed new theorem. The finite leading-term certificate is sufficient;
an inconclusive certificate or an unsuccessful bounded relation search is
not a contrary all-degree result.

The [G2 application](G2_CHARACTER_FORM.md) proves C7=Z[6omega1,2omega2]
and the two maximal-Levi invariant presentations over Z. It gives the
complete integral image presentations Z[X,Q]/(X^2-3Q,Q^3) and Z[Y]/(Y^6).
The source proof is all-degree; the image kernels follow from free quotient
bases with nonzero images in different degrees. The adjoint source is
Z[24omega1,8omega2], with a separate characteristic-three certificate.
A Python audit independent of Sage and the runtime library checks the
torus certificates, Laurent intersections and every image comparison.
These are checked deductions, with no claim of originality. The abstract
modular-obstruction example is not asserted to come from a representation.

The [F4 character-form proof](F4_CHARACTER_FORM.md) identifies the entire
26-dimensional torus source as Sym_Z(6 Q_long), hence the full integral
character source for every Levi by taking invariants. A unit-coefficient
leading-monomial certificate handles all primes simultaneously.
All four maximal-parabolic images are computed in every degree (74 rows).
Source lattices have integral-coefficient and gcd-one maximal-minor
certificates; every Borel column is checked by localization and independent
root/coroot Chevalley multiplication. A standard-library Python audit
reconstructs all source lattices and checks their certificates and image
arithmetic. This does not supply a finite integral invariant-ring
presentation or a classification of representations, and no novelty is claimed.

The [F4 adjoint comparison](F4_REPRESENTATION_COMPARISON.md) proves
C52=Sym_Z(18 Q_long), using separate characteristic-two and odd-prime
leading-term cases. Taking invariants and the common Borel map proves
A52,L^d=3^d A26,L^d and J52,L^d=3^d J26,L^d for every standard Levi and
every degree. This is a checked deduction, with no priority claim.
All 74 adjoint image rows are derived from this theorem and the audited
minimal tables. Fresh API checks cover 20 parabolic-degree pairs; 32
homogeneous/filtered Laurent certificates are independently reconstructed.
These finite checks support the implementation; they are not the
all-degree proof, and the two Springer morphisms are not scalar multiples.

The [G2/F4 note](EXCEPTIONAL_INTEGRAL.md) studies all maximal parabolics and
all degrees for two explicitly selected generator forms: least integer
image-denominator clearing, and individual primitive-image rescaling.
The complex generators are published Part III inputs. The G2 image
presentations follow by the displayed Chevalley calculation. F4 Smith lists
are exact finite computations, independently reconstructed by root/coroot
Chevalley multiplication; no general exceptional classification is claimed.
An additional F4 root-duality argument identifies the original selected image
algebras in the pairs P1/P4 and P2/P3, and identifies their odd-primary ambient
quotients degreewise. It depends on the displayed, independently checked
generator identities; it does not identify the ambient integral lattices or
their 2-primary quotients. The primitive-image forms have the same odd-primary
quotients as the original forms because their rescaling uses powers of 2.
Neither chosen form is identified with a representation-ring intersection,
and novelty has not been established. In particular, primitive individual
generator images do not guarantee saturation of the generated algebra.

The [E6/P1 and E7/P7 note](EXCEPTIONAL_MINUSCULE.md) adds complete tables
for both forms in all degrees. The known integral presentations are
Duan--Zhao Theorems 4 and 5 (with the E6 diagram involution). Fixed-point
interpolation checks every monomial against an ordinary-ring engine using
those presentations and Chevalley multiplication. Its starting generator
images and one E6 square anchor use the original BGG implementation; these
dependencies are explicit. The 90 degree/form rows contain 4744 columns
obtained from 2372 independently compared unscaled monomials and exact
rescaling. Integer identities in both directions prove that the primitive
images are Z[h,16s] and Z[h,81s,19683t]. These are deductions for the specified
forms inside published cohomology rings, not a general exceptional theorem
or a claim of novelty. Their quotients can involve only 2 and 3 respectively.

The subsequent [E6 character-form proof](E6_CHARACTER_FORM.md) identifies
the literal intersection with integral Levi characters, for simply
connected E6 and the normalized omega1 Springer map, as Sym_Z(6Q)^(W_L).
Here Q is the E6 root lattice. A unimodular change of trace-moment
coordinates has independent leading Laurent exponents over every prime
field, with characteristic 3 checked separately. Denominator reduction
then proves the full torus intersection, and equivariance gives the Levi
form. This establishes its grading in all degrees, not just a bounded
test. The argument is a repository deduction without a novelty claim.

For P1, all 17 cohomology degrees are computed using complete invariant
source lattices, with independent gcd-one minor certificates for their
completeness. Fixed-point and published-ring calculations agree. A separate
standard-library Python audit verifies the trace/leading-term argument,
ten low-degree Laurent intersections, all 51 comparison lattices and all
34 containments. The literal character image is contained in both earlier
chosen images; it already differs in degree 1. This resolves the E6/P1
source-form question.

The [E7/P7 character-form proof](E7_CHARACTER_FORM.md) subsequently proves
the analogous identification Sym_Z(12Q)^(W_L), with Q now the E7 root
lattice and the normalized 56-dimensional representation. Its all-primes
leading-term argument treats characteristic 2 separately. An integral
permutation basis for an A5 subgroup compresses coefficient tests without
changing the lattice. Complete invariant source lattices are certified by
gcd-one maximal minors, and Borel products are compared by the two earlier
independent methods. A standard-library audit separately checks the torus
certificate, permutation coordinates, low-degree Laurent intersections
and every embedded image comparison. The character and catalogue images
are incomparable: (216/5)psi5 is a character-form element with Borel image
outside the catalogue, while h is absent from the character image. This
is a deduction with explicit certificates, not a claim of novelty or of
a minimal polynomial presentation of the integral invariant ring.

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

## Orthogonal two-step joint images

The [B/D joint-image argument](ORTHOGONAL_JOINT_IMAGES.md) proves that
the sum of the normalized middle and tail images is primitive inside
their joint image, and that their intersection is the constants. It
assumes a,b>=1, k>=1 in B and k>=2 in D. The joint image need not be
primitive in ambient cohomology, but its nonzero Smith factors are powers
of 2. The note gives integral kernels and closed ranks, including the
free rank of the joint image modulo the separate-image sum.

In B, equality of the rational Cartesian kernels and torsion-free target
cohomology identify the integral source quotient with the type C one.
In D, the global Euler relation gives the additional joint relation
c_b p h_a(M^2 union T^2). A four-component decomposition into free
integral flag-ring modules proves the integral kernel and locates both
separate images exactly. These arguments do not use the earlier B/D
single-tail Smith multiplicity theorems.

These are repository deductions from classical inputs, without a novelty
claim or external certification. The [finite verification](experiments/orthogonal_joint_structure_verification.json)
checks exact integral kernels, ranks, image intersections and quotient
torsion, with selected independent Cartesian BGG columns. The remaining
problem is a closed formula for the ambient joint Smith multiplicities,
not the primitive-sum assertion.

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
