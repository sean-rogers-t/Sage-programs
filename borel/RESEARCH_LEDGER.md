# Springer and Borel research ledger

Updated 2026-10-20. The research priority is **generalizable statements**:
uniform formulas, structural explanations, sharp hypotheses and counterexamples.
Individual group calculations serve as tests of those statements and remain
valuable exact results with their own stated scope.

This is the working index of results, partial results and questions for this
program. An **open question** means unresolved here, not necessarily open in
the literature. Literature priority is separate from mathematical validity.
The supporting proofs and [proof/source audit](PROOF_STATUS.md) determine the
precise hypotheses; this ledger does not replace them or re-audit them.

The [research proposal and progress inventory](../articles/research-program-proposal/README.md)
organizes these questions into four core aims and extension track E, with
proof routes and milestone gates. Its priority is the local arithmetic of
integral images, followed by effective generic degrees and regularity criteria;
an integral geometric explanation of Smith filtrations is a supporting aim.
Track E adds equivariant, K-theoretic, formal-group-law, quantum, affine and
Kac–Moody projects. The first local-image benchmark and an equivariant
specialization counterexample are now recorded in SRC-08, IMG-10/11 and
LIMIT-06. General extension backends remain unimplemented. The expanded
ledger has 78 entries, including 20 open questions. The proposal itself
adds no theorem; subsequent progress is tracked here and in its README.

## Status and maintenance

| Status | Meaning |
| --- | --- |
| Published input | An established result used by the program; the linked note identifies its source. |
| Repository theorem | A written deduction with the proof/audit recorded here, under its stated hypotheses. This does not mean external review, formal verification or established novelty. |
| Exact finite result | A completed exact calculation or identity for the specified groups, representations, spaces and degrees. Exhausting one space does not establish an arbitrary-rank theorem. |
| Partial or conditional | A proved part, a conditional implication, or limited evidence; the missing step is identified. |
| Counterexample | A verified obstruction to the specified stronger statement. |
| Implemented method | A computational capability with an explanation and verification; not automatically a new mathematical theorem. |
| Open question | A question to investigate. It is not a conjecture unless an explicit proposed statement is supplied. |

Keep the identifiers below stable. When a question is resolved, retain the
question and link it to the result that resolves it. Every substantial change
should record its scope, proof or counterexample, computational evidence and
literature status. Do not promote examples to a theorem or an unlocated citation
to a novelty claim. A proposed new calculation should name the general question
it tests. Add further entries as the program develops.

## Maps and integral objects

For simply connected simple complex G, write T for a maximal torus,
P=X(T) for the weight lattice, Q for the root lattice, and W for the Weyl
group. A standard Levi L has Weyl group W_L and parabolic P_L. Let
S_Q=Sym_Q(P tensor Q). The normalized trace Springer morphism theta_rho
has differential equal to the identity at 1. Its torus restriction gives
an injective pullback theta_rho^*:S_Q -> Q[P].

The full integral torus source, linear lattice and Levi source are

\[
\mathcal C_\rho=\{f\in S_{\mathbb Q}:\theta_\rho^*f\in\mathbb Z[P]\},
\qquad \Lambda_\rho=(\mathcal C_\rho)_1,
\qquad \mathcal A_{\rho,L}=\mathcal C_\rho^{W_L}.
\]

The polynomial-character map in the Kumar papers is the composite
beta^{P_L} after (theta_rho^*)^{-1}, where that inverse is defined on the
image of the injective pullback. It is not a birational inverse of theta_rho.
The Borel image in polynomial coordinates is
J_{rho,L}=beta^{P_L}(A_{rho,L}). Integrality in cohomology must be justified
before treating this image as an integer lattice. A chosen normalized
generator algebra need not equal the full character source A_{rho,L}.

The source grading is on the polynomial side, not Laurent total degree.
Polynomial degree d means cohomology H^(2d). For an integral image I_d,
its nonzero Smith factors describe Sat(I_d)/I_d, where
Sat(I_d)=(I_d tensor Q) intersect H^(2d)(G/P_L;Z). This is an index in the
selected rational image, and becomes an index in all cohomology only if
that rational image has full rank.

## Published foundations and the earlier papers

| ID | Statement | Status and relationship to this program | Reference |
| --- | --- | --- | --- |
| BASE-01 | Rational Borel presentations, BGG divided differences, Schubert bases and parabolic pullbacks. | Published input. These supply the ambient rings and evaluation methods. | [Framework](MATHEMATICS.md) |
| BASE-02 | Trace projection and Rogers's weight-moment formula explicitly determine theta_rho on T. | Published input. The general coordinate formula is from Rogers's paper, not a new result of these extensions. | [Formula and attribution](SPRINGER_GEOMETRY.md) |
| BASE-03 | Polynomial characters transported by inverse Springer pullback map to flag cohomology. | Published input from Kumar Part I; Parts II/III supply the classical and exceptional generator calculations used here. | [Connection to the papers](MATHEMATICS.md), [articles](../articles/README.md) |
| BASE-04 | The generalized Cayley map has the same generic degree as its map T/W -> t/W; its full differential determinant enters harmonic localization. | Published input from Kostant--Michor. Our explicit coordinates and determinant factors make these constructions calculable. | [Geometry and published locators](SPRINGER_GEOMETRY.md) |
| BASE-05 | Integral theta/eta Schubert calculus, quaternionic flag descriptions and integral Leray--Hirsch. | Published inputs to the Smith and joint-image arguments; the ambient theories are not new results here. | [Source audit](PROOF_STATUS.md), [joint-image audit](TYPE_C_JOINT_AUDIT.md) |

## General source and representation results

The two broadest source results are

\[
\mathcal C_\rho=\operatorname{Sym}_{\mathbb Z}(\Lambda_\rho),
\qquad
\Lambda_\lambda=I_\lambda(Q+\mathbb Z\lambda)^*
\quad\text{for irreducible }V(\lambda),
\]

where the dual is for the invariant inner product with long roots of
squared length two, and I_lambda is the trace index with I_ad=2h^vee.
The multiplicity-weighted lattice is D_lambda=sum_mu Z(m_mu mu), using
the weights mu and their multiplicities m_mu in V(lambda).

| ID | Statement | Scope and status | Proof or evidence |
| --- | --- | --- | --- |
| SRC-01 | The full integral torus source is generated in degree one: C_rho=Sym_Z(Lambda_rho). | Repository theorem for the normalized trace construction on simply connected simple complex groups, with faithful differential, including direct sums and inhomogeneous polynomials. | [Theorem](TRACE_SOURCE_THEOREM.md), [audit](TRACE_SOURCE_AUDIT.md) |
| SRC-02 | For every nontrivial irreducible, the multiplicity-weighted moment lattice is D_lambda=Q+Z lambda and Lambda_lambda=I_lambda(Q+Z lambda)^*. | Repository theorem in every finite simple type. The exceptional classification is a corollary, including E8; explicit weights are unnecessary for the lattice formula. | [Classification](IRREDUCIBLE_CHARACTER_LATTICES.md), [audit](IRREDUCIBLE_LATTICE_AUDIT.md) |
| SRC-03 | The literal integral Levi source is Sym_Z(Lambda_rho)^(W_L). | Corollary of SRC-01 and equivariance for every standard Levi in its scope. This invariant ring need not be generated by linear invariants. SRC-08 now proves integrality of its Borel image in this scope. | [Source theorem and scope](TRACE_SOURCE_AUDIT.md), [integrality proof](LOCAL_IMAGE_BOUND.md) |
| SRC-04 | A coprime lattice-exponent/scaling comparison gives degreewise source and Borel-image sums and intersections. | Repository lattice lemma for certified sources Sym_Z(a Lambda), Sym_Z(ac Lambda'), Lambda subset Lambda', with c coprime to the exponent of Lambda'/Lambda. | [General lemma and E7 application](E7_REPRESENTATION_COMPARISON.md), [E6 application](E6_REPRESENTATION_COMPARISON.md) |
| SRC-05 | Repeating the same representation leaves the normalized Springer map and its full source unchanged. | General consequence of multiplying both trace matrix and trace moment by the same integer. It explains the reducible counterexample LIMIT-03. | [Direct-sum limitation](IRREDUCIBLE_CHARACTER_LATTICES.md) |
| SRC-06 | The specified literal defining/vector SO/Sp form scales each degree-d normalized image lattice by 2^d. | Repository normalization deduction for the character lattices specified in the audit. It does not apply automatically to Spin or other representations. | [Normalization proof](INTEGRAL_LITERATURE.md) |
| SRC-07 | Minimal/adjoint sources admit all-Levi comparison formulas in E6 and E7; F4 admits uniform 3^d image scaling. | Repository theorems for these representation pairs, in every degree. The E6/E7 formulas illustrate SRC-04, rather than classifying every pair. | [F4](F4_REPRESENTATION_COMPARISON.md), [E6](E6_REPRESENTATION_COMPARISON.md), [E7](E7_REPRESENTATION_COMPARISON.md) |
| SRC-08 | Lambda_rho is contained in the actual character lattice X, and the full literal Levi source has integral Borel image. | Uniform repository deduction in SRC-01's normalized simply connected simple complex scope: differentiate an integral Laurent pullback at 1, then use primitive Schubert pullback. Not a statement about arbitrary rationally normalized selected generators. | [Proof](LOCAL_IMAGE_BOUND.md), [audit](LOCAL_IMAGE_AUDIT.md) |

## General integral-image results and classical families

The orthogonal single-tail results use exactly q_i=e_i(x_tail^2), with
p=product(x_tail) and p^2=q_k in D. Type B has a,k>=1; the main type D
theorems have a>=1,k>=2. Changes of integral source are explicit, not implicit.

| ID | Statement | Scope and status | Proof or evidence |
| --- | --- | --- | --- |
| IMG-01 | Combinatorial generator formulas and stable specialization arguments work in types A/C/D. | Earlier repository alternative proofs/refinements of published main results. Includes the completed dissertation argument, type A sign convention and type D Euler class; not a claim that Parts II/III lacked proofs. | [Classical proof note and comparison](../articles/combinatorial-proofs/README.md) |
| IMG-02 | The normalized type C tail comparison supplies a primitive integral Schur lattice. | Established integral geometry applied in the repository proofs. This supplies the primitive reference lattice needed to interpret B row weights. | [Orthogonal structure](ORTHOGONAL_SMITH.md), [B filtration](TYPE_B_FILTRATION.md) |
| IMG-03 | Every nonzero normalized type B Smith factor and its multiplicity has a closed positive partition/Gaussian formula. | Repository theorem for all ranks and allowed blocks, including primitivity of every required prefix. The finite tables are corroboration. | [General proof](TYPE_B_GENERAL_SMITH.md), [formula](TYPE_B_CLOSED_CONJECTURE.md), [audit](TYPE_B_GENERAL_AUDIT.md) |
| IMG-04 | Every type D square-sector Smith factor agrees with the corresponding type B factor. | Repository theorem for all allowed blocks and degrees. Earlier two-row and three-row arguments remain independent special-case proofs. | [General comparison](TYPE_D_SQUARE_THEOREM.md), [two-row proof](TYPE_D_TWO_ROW_SMITH.md) |
| IMG-05 | The full normalized type D formula follows from square/Euler sectors and unit-level gluing; there is no further positive-level gluing. | Repository theorem for all allowed blocks, with the even-tail step included. | [General theorem](TYPE_D_EVEN_GLUING_THEOREM.md), [audit](TYPE_D_EVEN_GLUING_AUDIT.md) |
| IMG-06 | Stable primitivity, explicit sufficient degree ranges, block duality and exponent bounds hold for the specified B/D images. | Uniform repository deductions. These concern the selected image, not surjectivity onto all cohomology. | [Overview and bounds](ORTHOGONAL_SMITH.md), [D exponent bound](TYPE_D_FULL_EXPONENT_BOUND.md) |
| IMG-07 | Type B Smith prefixes have an integral incidence/defect interpretation for a square-zero operator. | Repository theorem linking Schur-product rows, primitive prefixes and defect-stratum ranks. Classical Birkhoff/Hall/profile-counting polynomials are published inputs, not new enumerative polynomials. | [Defect comparison](TYPE_B_DEFECT_COMPARISON.md), [integral completion](TYPE_B_GENERAL_SMITH.md) |
| IMG-08 | On normalized type C two-step flags, the joint image is saturated, the separate-image sum is primitive, and the separate intersection is Z. | Repository corollary of classical integral geometry for every a,b,k>=1 and every degree. Integral kernels and graded ranks are also established. | [Theorem](TYPE_C_JOINT_PRIMITIVITY.md), [audit](TYPE_C_JOINT_AUDIT.md) |
| IMG-09 | On normalized B/D two-step flags, the separate-image sum is primitive inside the joint image, the intersection is Z, and ambient torsion is 2-primary. | Repository theorem for a,b>=1, k>=1 in B and k>=2 in D. Ambient joint saturation is false in general; complete ambient multiplicities remain OPEN-07. | [Orthogonal joint theorem](ORTHOGONAL_JOINT_IMAGES.md) |
| IMG-10 | Every full-source degree-d Borel-image cokernel is killed by t_G e_rho^d, where e_rho=exp(X/Lambda_rho); good primes also recover the full equivariant image. | Explicit corollary of classical torsion-index/BGG machinery and the source inclusion, without a Levi Weyl-order factor. Equivariant triangular induction proves good-prime surjectivity; each graded specialization kernel is finite with primes dividing t_G e_rho. This is not a sharp bad-prime Smith/kernel classification. | [Proof and attribution](LOCAL_IMAGE_BOUND.md), [audit](LOCAL_IMAGE_AUDIT.md) |
| IMG-11 | For defining SL_n on SL_n/P_1=P^(n-1), the full literal degree-d image is (n/gcd(n,d)) Z h^d for 0<=d<n. | Written all-n projective-family deduction: Lambda=Q, Levi invariants are elementary root differences, their images are binomial(n,j) h^j, and a valuation argument computes the monomial gcd. This root-characteristic formulation needs further literature comparison; no novelty claim. | [Family and proof](LOCAL_IMAGE_BOUND.md), [independent checks](experiments/local-image-independent.json) |
| IMG-12 | E6 root-to-weight image quotients have exact 3-primary multiplicities on all 63 proper standard parabolics. | A Cayley-plane module-basis proof and good-Levi deconvolution cover 60 cases. The remaining A5 and A2+A2+A1 cases have complete coefficient-lattice proofs, independently replayed in all 52 base degrees; A2+A2 follows by a common (1,z) basis and multiplication by 1+q. The bound 27 and full-flag top ratio 9 follow from published torsion-index results. These are computer-assisted repository deductions, with priority unresolved; ambient H/N is separate. Shared-prime scalar irreducible comparisons follow by min/max valuations. | [60-case theorem](E6_ROOT_WEIGHT_COMPARISON.md), [completion](E6_REMAINING_CASES.md), [tables](E6_REMAINING_TABLES.md), [completion audit](E6_REMAINING_AUDIT.md) |
| IMG-13 | E7/P7 has explicit simultaneous root/weight image bases at 2; the quotient is elementary 2 with polynomial q[14]_(q2)+(q6+q10)[7]_(q2). | N has towers h^i (i<=27), h^i(2s), h^i(2t) (i<=13); M doubles exactly the odd h exponents. Candidate-ring identities are deduced from the published cohomology presentation; both complete sources and all images are independently replayed in every degree. An index-3 stable common complement and a root quadratic lift give exponent 2 on all 64 retained subsets of the E6 Levi, without giving their multiplicities. General bound 8 is known; priority of the sharper deductions is unresolved. | [Basis theorem and refinement bound](E7_ROOT_WEIGHT_COMPARISON.md), [proof/source audit](E7_ROOT_WEIGHT_AUDIT.md), [independent replay](experiments/e7-two-primary-basis-audit.json) |
| IMG-14 | All 64 refinements of E7/P7 have proved relative multiplicities through degree 9: m_(I,d)=[q^(d-1)] (1-q6)(1-q8)/((1-q)^(6-|I|) product_j(1-q^(d_j))). | Uniform computer-assisted deduction in degrees 0..9 only, for the full root/weight invariant sources. Independent 56-weight and 27-form relations, primitive kernel columns, and full-source A5 BGG checks agree. IMG-15 now settles the full flag and 61 of these refinements in all degrees; D5/D4 higher factors remain. Published Borel/Levi inputs are separated from these deductions; novelty unresolved. | [Refinement theorem and audit](E7_REFINEMENT_MULTIPLICITIES.md), [independent 64-case table](experiments/e7-refinement-independent.json) |
| IMG-15 | E7 full-flag relative factors are all 2, with polynomial q P_W/(1+q). Good-Levi transfer settles 108 proper parabolics; P7 adds a 109th. | A full-source A6 intermediary gives an exhaustive 576-dimensional image certificate, fresh adjoint-ideal reconstruction and plain-Python Buchberger/image audits. Classical common Artin fiber bases transfer the result. This settles 61/64 P7 refinements; D5 pair and D4 remain. Direct saturated full-flag kernels in degrees 10..13 independently corroborate the formula and exhibit integral half-sum corrections. All irreducible scalar Q/X comparisons follow on this scope; ambient H/N, arbitrary direct sums and novelty remain separate. | [All-degree theorem and audit](E7_FULL_FLAG_COMPARISON.md), [transfer scope and tables](experiments/e7-full-flag-transfer.json) |

## General Springer and Cayley geometry

| ID | Statement | Scope and status | Proof or evidence |
| --- | --- | --- | --- |
| GEO-01 | A normalized trace map has identity differential at 1, hence is dominant and generically finite. | General consequence of the published construction with faithful differential in characteristic zero. It does not determine the degree. | [Setup and normalization](SPRINGER_GEOMETRY.md) |
| GEO-02 | The logarithmic torus derivative is J_T=M_rho^(-1) sum_mu m_mu mu mu^t e^mu. | Uniform exact differentiation of Rogers's formula. | [General differential argument](SPRINGER_GEOMETRY.md) |
| GEO-03 | The full determinant is det(J_T) times the product over positive roots of e^alpha (alpha(theta)/(e^alpha-1))^2. | Uniform repository derivation from equivariance. The quotients extend regularly as Laurent polynomials; the root directions cannot be omitted. | [Full differential formula](SPRINGER_GEOMETRY.md) |
| GEO-04 | The quadratic invariant pullback has a tensor-Casimir character formula. | Uniform deduction using the standard tensor-Casimir identity and the specified trace normalization; not a new Casimir identity. | [Formula and attribution](SPRINGER_GEOMETRY.md) |
| GEO-05 | Adjoint-quotient degree equals the full Springer degree; a central isogeny through which the map factors contributes its finite cover degree. | Published quotient-degree comparison plus ordinary multiplicativity of function-field degrees. Used explicitly in the E6/E7 calculations. | [Degree setup](EXCEPTIONAL_DEGREES.md) |

## Exact exceptional cases and checks

These entries have finite scope even when every degree of a fixed space has
been computed. They are useful tests and examples for SRC, IMG and GEO, not
substitutes for a uniform theorem.

| ID | Completed result | Scope and status | Record |
| --- | --- | --- | --- |
| CASE-01 | Independently reconstructed minimal-representation Springer coordinates for F4/E6/E7; proposed E6 denominator correction from 6 to 18. | Exact checks of the specified formulas. The correction note is a proposed correction, not a published erratum. F4/E7 coordinate checks agree with the printed formulas. | [Article verification and note for Shrawan](../articles/README.md) |
| CASE-02 | Complete selected integral-image tables for every maximal parabolic of G2/F4. | Exact finite results for two specified generator forms, with independent ordinary-ring checks; not a general exceptional Smith classification. | [Scope](EXCEPTIONAL_INTEGRAL.md), [tables](EXCEPTIONAL_INTEGRAL_TABLES.md) |
| CASE-03 | G2 minimal character-source and both maximal-parabolic image presentations; G2 adjoint torus form; F4 minimal/adjoint source forms and maximal-parabolic images. | Source identifications are written deductions; the specified complete image calculations are finite cases. An adjoint torus identification alone is not a table of all its Levi images. | [G2](G2_CHARACTER_FORM.md), [F4](F4_CHARACTER_FORM.md), [F4 comparison](F4_REPRESENTATION_COMPARISON.md) |
| CASE-04 | Complete E6/P1 and E7/P7 selected-image tables and primitive-generator subalgebras Z[h,16s] and Z[h,81s,19683t]. | Exact descriptions inside the published integral cohomology rings, with inherited relations. The algebras are not free polynomial rings. | [Structure and verification](EXCEPTIONAL_MINUSCULE.md) |
| CASE-05 | Complete E6/P1 and E7/P7 literal-character images, including minimal/adjoint embedded comparisons. | Exact finite image calculations supporting the all-Levi source/comparison statements; this is not coverage of every exceptional parabolic. | [E6 form](E6_CHARACTER_FORM.md), [E7 form](E7_CHARACTER_FORM.md), [API](MINUSCULE_CHARACTER_API.md) |
| CASE-06 | Complete G2 minimal/adjoint quotient maps, full determinants, degrees 6 and 18, and regular-to-nonregular witnesses. | Exact calculations and written degree arguments. The witnesses establish the localization obstruction for these two maps using the published criterion. | [G2 geometry](SPRINGER_GEOMETRY.md) |
| CASE-07 | Complete E6/E7 minimal/adjoint quotient coordinates and full determinant products. | Exact character formulas with support/unisolvence certificates and independent representation matrices. The determinants are explicit products rather than single expanded character polynomials. | [E6](E6_GEOMETRY.md), [E7](E7_GEOMETRY.md) |
| CASE-08 | Generic degrees E6 1080/2754 and E7 7272/33664, in minimal/adjoint order on simply connected groups. | Repository geometric deductions with explicit boundary contributions. Seven exact specialized-fiber reports corroborate them: two primes each for the first three cases, one for E7 adjoint. Finite-field agreement is not the characteristic-zero proof. | [Proof and coverage](EXCEPTIONAL_DEGREES.md), [release audit](release-0.20-verification.json) |

## Counterexamples and limits that a general theorem must respect

| ID | Stronger statement ruled out or restricted | Verified obstruction | Record |
| --- | --- | --- | --- |
| LIMIT-01 | Primitive individual generators force a saturated generated algebra. | E6/P1 and E7/P7 primitive-generator algebras still have nontrivial saturation quotients. | [Explicit subalgebras](EXCEPTIONAL_MINUSCULE.md) |
| LIMIT-02 | Changing from minimal to adjoint always gives containment or scalar rescaling. | E6/P1 and E7/P7 minimal/adjoint images are incomparable already in divisor degree. | [E6](E6_REPRESENTATION_COMPARISON.md), [E7](E7_REPRESENTATION_COMPARISON.md) |
| LIMIT-03 | Trace index and central character classify the source for reducible representations. | Three F4 minimal modules and the adjoint both have index 18, but their linear lattices are 6Q_long and 18Q_long. Irreducibility in SRC-02 matters. | [Counterexample](IRREDUCIBLE_CHARACTER_LATTICES.md) |
| LIMIT-04 | Type C joint saturation extends unchanged to B/D. | B3 with blocks (1,1,1) has degree-three factors (1,2); D4 with (1,1,2) has degree-four factors (1,1,2,2,2). | [Joint theorem and examples](ORTHOGONAL_JOINT_IMAGES.md) |
| LIMIT-05 | Equal rational kernels, abstract image rings or Smith lists identify embedded integer images. | B/C joint rings can be isomorphic with different ambient Smith factors. Embedded comparison needs actual Schubert lattices. | [Orthogonal joint images](ORTHOGONAL_JOINT_IMAGES.md), [comparison methods](INTEGRAL_COMPARISONS.md) |
| LIMIT-06 | An integral equivariant Springer image always specializes faithfully to its ordinary image. | False for the defining SL4 map on P3 with the full literal source Sym_Z(Q)^(S3). For J_T=S[4h,2h^2,4h^3+2ah], the degree-five specialization kernel is exactly Z/2, witnessed by nu=2a^2h. Polynomial degree five means cohomology ten. Lower degrees specialize faithfully; higher defects are not classified. | [Integral proof](EQUIVARIANT_IMAGE_PILOT.md), [audit](LOCAL_IMAGE_AUDIT.md), [independent rank-52 lattice check](experiments/local-image-independent.json) |
| LIMIT-07 | Six rational E6 basic invariants can be normalized to generate the whole 2-local invariant ring in the common E7 complement. | False: complete source components in degrees 2,5,6,8,9,12 force the displayed six-generator subalgebra, but an integral degree-16 Chern invariant has odd value at a common mod-2 zero of those generators. Independently checked with integer multinomial expansions. This obstructs a source-basis shortcut, not the proposed Borel-image multiplicities. | [Explicit obstruction](E7_REFINEMENT_MULTIPLICITIES.md#5-a-concrete-obstruction-to-a-polynomial-source-basis-shortcut), [independent audit](experiments/e7-refinement-invariant-obstruction-audit.json) |

## Partial results and conditional routes

| ID | Available result | Missing step or restriction | Record |
| --- | --- | --- | --- |
| PART-01 | A stronger unweighted typed-prefix criterion would imply type D gluing. | The lattice implication is proved, but its stronger strictness hypotheses are not proved for arbitrary parameters. IMG-05 is already proved by a different route and does not require them. | [Conditional criterion](TYPE_D_GLUING_FILTRATION.md) |
| PART-02 | F4 root duality compares the selected image algebras at paired maximal parabolics and their odd-primary quotients. | It does not identify ambient integer lattices or their 2-primary quotients. | [Scope in the audit](PROOF_STATUS.md) |
| PART-03 | Rational joint-image kernels and integral C/B/D two-step structural theorems are available. | Closed ambient B/D joint Smith multiplicities and a general many-block theory are not supplied by the single-tail formulas. | [Rational derivation](experiments/TWO_STEP_DERIVATION.md), [joint theorem](ORTHOGONAL_JOINT_IMAGES.md) |
| PART-04 | Exact square-completion and quadratic-module presentations preserve the E7 adjoint specialized fiber, including nilpotents. | The supplementary module/shifted and second-prime runs did not yield completed counts before finalization. Only the completed original-chart count is numerical evidence for CASE-08's E7 adjoint check. | [Methods and completed-count scope](EXCEPTIONAL_DEGREES.md), [attempt history](experiments/exceptional-degree-bounded-attempts.json) |

## Implemented general methods

| ID | Capability | Scope and mathematical role | Documentation |
| --- | --- | --- | --- |
| API-01 | Exact Borel/BGG evaluation for arbitrary standard parabolics in supported finite crystallographic types. | A single omitted node is not a restriction of the general map API. Classical coinvariant normal forms avoid many expensive Schubert expansions. | [Framework](MATHEMATICS.md), [methods](METHODS.md) |
| API-02 | Complete rational kernels, quotient bases and Hilbert series for chosen homogeneous generators. | Finite-border completeness certificates distinguish a full kernel from a guessed bounded-degree list. | [Presentations](PRESENTATIONS.md) |
| API-03 | Arbitrary Schubert multiplication, powers, intersections and inverse Borel lifts. | General exact calculations; inverse Borel returns a polynomial representative of a cohomology class. | [Schubert calculus](SCHUBERT_CALCULUS.md) |
| API-04 | Integral images, nonzero Smith factors, saturation, membership and integer preimages. | Computes the embedded image for explicitly integral generator images, preserving source normalization. | [Integral images](INTEGRAL_IMAGES.md) |
| API-05 | Image containment, equality, sums, intersections and quotient groups; compressed all-rank B/D formulas. | Embedded comparisons retain directions, while the closed orthogonal calculator avoids constructing enormous Schubert matrices. | [Comparisons](INTEGRAL_COMPARISONS.md), [orthogonal API](ORTHOGONAL_API.md) |
| API-06 | Representation sources, highest-weight-only irreducible lattices, quotient coordinates and determinant evaluation. | The general constructors and the packaged exceptional cases have different computational costs and guards. A packaged degree certificate is bound to its actual map. | [Lattice API](IRREDUCIBLE_CHARACTER_LATTICES.md), [geometry API](SPRINGER_GEOMETRY.md), [degree API](EXCEPTIONAL_DEGREES.md) |

## Open questions organized by generalization

These are questions for the program. A fresh literature check is required
before presenting any as an unsolved problem in the wider literature.
Current work centers on OPEN-02/09, with OPEN-14 as a bounded exploratory
project and OPEN-05 as the next sharp comparison target. OPEN-01 remains
the principal separate geometric direction. Row order is not a priority
ranking, and partial progress below does not close an entire question.

| ID | Question | What would constitute progress | Existing starting point |
| --- | --- | --- | --- |
| OPEN-01 | Can generic Springer degrees be computed uniformly from representation/root data, including face degeneracies? | A theorem or certified algorithm specifying the Newton lattice, all contributing face schemes and boundary multiplicities. A raw volume alone is insufficient. F4/E8 would test the mechanism. | GEO-05, CASE-08; [boundary calculations](EXCEPTIONAL_DEGREES.md) |
| OPEN-02 | Which structural conditions control the primes and exponents in integral Springer/Borel images? | First benchmark complete: IMG-10 proves full-source good-prime surjectivity and an exponent bound; IMG-11 gives sharp projective factors. IMG-12 gives the complete E6 relative comparison; IMG-13 gives an E7/P7 basis formula and exponent 2 on 64 refinements; IMG-15 gives all-degree relative factors on 109 proper E7 parabolics. The remaining target is sharper bad-prime structure beyond these cases and separate criteria for selected generator algebras. | SRC-01/02/08, IMG-03/05/09/10/11/12/13, CASE-02/04 |
| OPEN-03 | Can Hiller-type chain models explain the Smith filtration structurally? | Specify the exact chain model and source statement; construct an integral comparison to the coefficient-row prefixes and prove the necessary primitivity. Matching Betti numbers alone will not determine Smith factors. | IMG-07; [primitive B prefixes](TYPE_B_GENERAL_SMITH.md) |
| OPEN-04 | What replaces the irreducible closed source-lattice formula for arbitrary direct sums? | A formula from the summands' highest-weight data that handles overlapping weight spaces, beyond the existing exact combined coefficient-lattice calculation. | SRC-01/02/05, LIMIT-03 |
| OPEN-05 | Can representation-image sums and intersections be classified when the scaling and lattice exponent are not coprime? | IMG-12 resolves shared-prime scalar root/weight comparisons for irreducible E6 sources on all 63 proper parabolics, with G/G trivial. IMG-13 proves the E7/P7 basis formula and exponent 2 for all 64 refinements; IMG-14 gives multiplicities through degree 9 on all refinements. IMG-15 settles all degrees on 109 proper E7 parabolics, including 61/64 P7 refinements. The other 18 proper E7 parabolics, other root data and arbitrary direct sums remain. | SRC-04/07, IMG-12/13, LIMIT-02/05 |
| OPEN-06 | Can the full critical divisor and its irreducible factors be described uniformly across representations? | A structural analysis of the Cartan determinant and root factors beyond the universal product identity; identify representation-theoretic causes of components and multiplicities. | GEO-02/03, CASE-06/07 |
| OPEN-07 | What are the complete ambient B/D joint-image Smith multiplicities on two-step flags? | Closed formulas with an integral proof accounting for mixed generators. The single-tail B/D theorem does not settle this embedding. | IMG-09, PART-03, LIMIT-04 |
| OPEN-08 | Which joint-image results extend to arbitrary numbers of Levi blocks? | Precisely stated generators and an integral compatible-basis, saturation or torsion theorem for many-step flags. | IMG-08/09; [extension notes](EXTENSIONS.md) |
| OPEN-09 | Can full integral Levi invariant rings and their Borel images be computed uniformly from Lambda_rho? | SRC-08 justifies integral Borel images, and IMG-10 bounds their cokernels. IMG-11 solves the defining projective family; IMG-12 gives sharp relative E6 root/weight factors on all 63 proper parabolics, without classifying the ambient weight-image cokernel. IMG-13 supplies explicit local E7/P7 image bases and a refinement exponent bound. IMG-15 completes the relative full-flag comparison through an A6 intermediary and covers 109 proper E7 parabolics; ambient weight-image cokernels remain separate. | SRC-01/03/08, IMG-10/11/12/13, API-02/04 |
| OPEN-10 | Which parts of the source classification extend to other isogeny forms or reductive groups? | A theorem using the actual character lattice and all trace-nondegeneracy hypotheses, including central torus directions where relevant. | SRC-01/02, GEO-05; [scope audit](TRACE_SOURCE_AUDIT.md) |
| OPEN-11 | What are the target branch loci and nonproperness loci of the explicit Springer quotient maps? | Compute the images/closures of critical loci and account for fibers escaping to the boundary. A source determinant alone does not give these target loci. | GEO-03, CASE-07/08 |
| OPEN-12 | When can harmonic localization be removed for these representation maps? | Apply the exact Kostant--Michor criterion, seeking uniform conditions or regular-to-nonregular witnesses. The two G2 cases already have an obstruction. | BASE-04, CASE-06; [localization discussion](SPRINGER_GEOMETRY.md) |
| OPEN-13 | Which uniform exceptional image patterns survive beyond the currently checked parabolics? | IMG-13 proves exponent 2 on all E7/P7 refinements; IMG-14 supplies multiplicities through degree 9. IMG-15 proves the full-flag formula and settles 61/64 P7 refinements in all degrees; the remaining independent targets are D5 and D4. It covers 109 proper E7 parabolics. LIMIT-07 blocks the six-generator E6 source shortcut, while the A6 intermediary supplies a different common basis. Further predictions from SRC-02 and OPEN-02/09, tested on additional parabolics or F4/E8 representations, followed by proof or a counterexample. | CASE-02/04/05, IMG-13; [catalogue](EXCEPTIONAL.md) |
| OPEN-14 | What is the correct equivariant integral-image invariant? | First Track E1 pilot complete: IMG-10 proves good-prime recovery and prime support of graded kernels; LIMIT-06 gives an exact bad-prime specialization defect. The next target is a criterion or structural description at the allowed bad primes. A general backend is not implemented. | [Good-prime theorem](LOCAL_IMAGE_BOUND.md), [counterexample](EQUIVARIANT_IMAGE_PILOT.md), [equivariant scope](EXTENSIONS.md#equivariant-calculations) |
| OPEN-15 | Is there a useful quantum extension of the selected Springer-image construction? | Track E4: define the source, coefficient ring and multiplicative lift; prove q=0 recovery and then analyze the deformed kernel and image. Quantum products or a Schubert-basis identification alone do not provide this lift. | [Quantum scope](EXTENSIONS.md#quantum-products-and-affine-transport), [program](../articles/research-program-proposal/README.md) |
| OPEN-16 | How do these representation-map degrees compare with minimum Cayley degrees? | Identify applicable published minima/bounds and clarify what varying rho can achieve. CASE-08 gives degrees of particular maps, not minima over all equivariant rational maps. | [Degree scope](EXCEPTIONAL_DEGREES.md) |
| OPEN-17 | How does the K-theoretic image of integral Springer polynomial characters compare with our cohomological image? | Track E2: restrict the associated-bundle map R(L) -> K^0(G/P), specify induced filtrations, and prove or obstruct an associated-graded comparison. The Springer grading and augmentation/gamma filtration must not be identified without proof. | [K-theory scope](EXTENSIONS.md#other-cohomology-theories-and-groups), [program](../articles/research-program-proposal/README.md) |
| OPEN-18 | Can representation-dependent sources and images be defined compatibly with a formal group law? | Track E3: start with connective K-theory, specify formal group algebra, coefficients and completion, and prove a specialization comparison or an obstruction. Existing oriented characteristic maps do not automatically extend the additive source theorem. | [Oriented-theory scope](EXTENSIONS.md#other-cohomology-theories-and-groups), SRC-01/03 |
| OPEN-19 | What affine image statements follow from a quantum Springer-image lift under the Peterson correspondence? | Track E5: after OPEN-15, transport an image to the appropriate localized quotient of affine-Grassmannian Pontryagin homology, specify the inverted classes, and determine information lost under localization. Unlocalized integral recovery needs a separate argument. | [Quantum/affine route](EXTENSIONS.md#quantum-products-and-affine-transport), [Lam–Shimozono](https://arxiv.org/abs/0705.1386) |
| OPEN-20 | Which characteristic-image constructions extend to affine or general Kac–Moody flags, and is there a suitable Springer analogue? | Track E6: implement a specified degree-bounded nil-Hecke setting using published foundations; separately define the group, Cartan lattice and completion for any candidate Springer map and prove it well defined. Ordinary finite traces on an infinite-dimensional faithful representation are unavailable. | [Kac–Moody scope](EXTENSIONS.md#other-cohomology-theories-and-groups), [Kostant–Kumar](https://doi.org/10.1073/pnas.83.6.1543) |

## Resolved questions and earlier routes

| Earlier question | Current resolution | Remaining distinction |
| --- | --- | --- |
| Is the full torus character source always generated in degree one? | SRC-01. | This does not imply degree-one generation after Levi invariants or saturation after the Borel map. |
| Can irreducible exceptional source lattices be classified uniformly? | SRC-02, in fact in every finite simple type. | Direct sums require OPEN-04. |
| Is the general normalized type B partition formula only conjectural? | IMG-03 and its integral prefix proof close the gap. | The historical file name TYPE_B_CLOSED_CONJECTURE.md is retained; it no longer denotes a current conjecture. |
| Does the full normalized type D formula still have an even-tail gluing gap? | IMG-05 closes it in the stated range. | PART-01's stronger conditional route remains a separate question. |
| Are type C two-step joint saturation and primitive sums only experimental? | IMG-08 supplies the arbitrary-rank argument. | B/D ambient joint saturation is false; OPEN-07 asks for its factors. |
| Is the E7 adjoint degree missing a full exact fiber check? | CASE-08 includes a completed prime-32003 count, independently recounted. | The characteristic-zero proof is geometric; a second-prime E7 adjoint count is not recorded. |
| Is full literal-source Borel integrality still an extra unproved hypothesis in the normalized source-theorem scope? | SRC-08 proves it. | Arbitrary selected rationally normalized sources still require separate checks. |
| Must a conservative good-prime theorem exclude primes dividing the Levi Weyl-group order? | IMG-10 uses integral divided-difference lifts and does not require that exclusion. | Bad-prime exponents are not classified by this bound. |
| Does every equivariant image specialize faithfully to its ordinary image? | LIMIT-06 supplies an SL4/P1 counterexample, with first kernel Z/2 in degree five. | Conditions eliminating or controlling the kernel remain OPEN-14. |

## How to choose the next calculation

Start with a precise general question above and the obstruction that prevents
a proof. Choose a case that discriminates between competing mechanisms, and
retain its exact source lattice and normalization. Record the result as a
finite test until a separate argument establishes the general statement.
F4, E8, additional representations and more parabolics are useful test cases
for this purpose, not automatic substitutes for a structural result.

The full-paper project remains on hold. This ledger, the supporting proof
notes and the API can develop independently of writing a new paper.
