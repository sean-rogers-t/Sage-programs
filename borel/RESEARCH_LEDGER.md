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
Kac–Moody projects. These are proposed extensions, not completed theorems
or implemented backends. The expanded ledger has 69 entries, including
20 open questions. The proposal adds no theorem.

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
| SRC-03 | The literal integral Levi source is Sym_Z(Lambda_rho)^(W_L). | Corollary of SRC-01 and equivariance for every standard Levi in its scope. This invariant ring need not be generated by linear invariants. | [Source theorem and scope](TRACE_SOURCE_AUDIT.md) |
| SRC-04 | A coprime lattice-exponent/scaling comparison gives degreewise source and Borel-image sums and intersections. | Repository lattice lemma for certified sources Sym_Z(a Lambda), Sym_Z(ac Lambda'), Lambda subset Lambda', with c coprime to the exponent of Lambda'/Lambda. | [General lemma and E7 application](E7_REPRESENTATION_COMPARISON.md), [E6 application](E6_REPRESENTATION_COMPARISON.md) |
| SRC-05 | Repeating the same representation leaves the normalized Springer map and its full source unchanged. | General consequence of multiplying both trace matrix and trace moment by the same integer. It explains the reducible counterexample LIMIT-03. | [Direct-sum limitation](IRREDUCIBLE_CHARACTER_LATTICES.md) |
| SRC-06 | The specified literal defining/vector SO/Sp form scales each degree-d normalized image lattice by 2^d. | Repository normalization deduction for the character lattices specified in the audit. It does not apply automatically to Spin or other representations. | [Normalization proof](INTEGRAL_LITERATURE.md) |
| SRC-07 | Minimal/adjoint sources admit all-Levi comparison formulas in E6 and E7; F4 admits uniform 3^d image scaling. | Repository theorems for these representation pairs, in every degree. The E6/E7 formulas illustrate SRC-04, rather than classifying every pair. | [F4](F4_REPRESENTATION_COMPARISON.md), [E6](E6_REPRESENTATION_COMPARISON.md), [E7](E7_REPRESENTATION_COMPARISON.md) |

## General classical integral-image results

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
before presenting any as an unsolved problem in the wider literature. The
first five are the current structural priorities; no proof is asserted here.

| ID | Question | What would constitute progress | Existing starting point |
| --- | --- | --- | --- |
| OPEN-01 | Can generic Springer degrees be computed uniformly from representation/root data, including face degeneracies? | A theorem or certified algorithm specifying the Newton lattice, all contributing face schemes and boundary multiplicities. A raw volume alone is insufficient. F4/E8 would test the mechanism. | GEO-05, CASE-08; [boundary calculations](EXCEPTIONAL_DEGREES.md) |
| OPEN-02 | Which structural conditions control the primes and exponents in integral Springer/Borel images? | A proved prime-support or saturation criterion with exact source and parabolic hypotheses, plus counterexamples to overbroad versions. | SRC-01/02, IMG-03/05/09, CASE-02/04 |
| OPEN-03 | Can Hiller-type chain models explain the Smith filtration structurally? | Specify the exact chain model and source statement; construct an integral comparison to the coefficient-row prefixes and prove the necessary primitivity. Matching Betti numbers alone will not determine Smith factors. | IMG-07; [primitive B prefixes](TYPE_B_GENERAL_SMITH.md) |
| OPEN-04 | What replaces the irreducible closed source-lattice formula for arbitrary direct sums? | A formula from the summands' highest-weight data that handles overlapping weight spaces, beyond the existing exact combined coefficient-lattice calculation. | SRC-01/02/05, LIMIT-03 |
| OPEN-05 | Can representation-image sums and intersections be classified when the scaling and lattice exponent are not coprime? | A p-primary replacement for SRC-04, including the effect of the Borel kernel; compare embedded images rather than Smith lists alone. | SRC-04/07, LIMIT-02/05 |
| OPEN-06 | Can the full critical divisor and its irreducible factors be described uniformly across representations? | A structural analysis of the Cartan determinant and root factors beyond the universal product identity; identify representation-theoretic causes of components and multiplicities. | GEO-02/03, CASE-06/07 |
| OPEN-07 | What are the complete ambient B/D joint-image Smith multiplicities on two-step flags? | Closed formulas with an integral proof accounting for mixed generators. The single-tail B/D theorem does not settle this embedding. | IMG-09, PART-03, LIMIT-04 |
| OPEN-08 | Which joint-image results extend to arbitrary numbers of Levi blocks? | Precisely stated generators and an integral compatible-basis, saturation or torsion theorem for many-step flags. | IMG-08/09; [extension notes](EXTENSIONS.md) |
| OPEN-09 | Can full integral Levi invariant rings and their Borel images be computed uniformly from Lambda_rho? | Certified generators/relations and image descriptions; distinguish the already solved torus degree-one theorem from higher-degree invariant generators. | SRC-01/03, API-02/04 |
| OPEN-10 | Which parts of the source classification extend to other isogeny forms or reductive groups? | A theorem using the actual character lattice and all trace-nondegeneracy hypotheses, including central torus directions where relevant. | SRC-01/02, GEO-05; [scope audit](TRACE_SOURCE_AUDIT.md) |
| OPEN-11 | What are the target branch loci and nonproperness loci of the explicit Springer quotient maps? | Compute the images/closures of critical loci and account for fibers escaping to the boundary. A source determinant alone does not give these target loci. | GEO-03, CASE-07/08 |
| OPEN-12 | When can harmonic localization be removed for these representation maps? | Apply the exact Kostant--Michor criterion, seeking uniform conditions or regular-to-nonregular witnesses. The two G2 cases already have an obstruction. | BASE-04, CASE-06; [localization discussion](SPRINGER_GEOMETRY.md) |
| OPEN-13 | Which uniform exceptional image patterns survive beyond the currently checked parabolics? | Predictions from SRC-02 and OPEN-02/09, tested on additional parabolics or F4/E8 representations, followed by proof or a counterexample. | CASE-02/04/05; [catalogue](EXCEPTIONAL.md) |
| OPEN-14 | What is the correct equivariant integral-image invariant? | Track E1: define the integral image over the torus polynomial ring, give a module presentation/Fitting or finite graded description, and prove the ordinary-specialization comparison including any base-change defect. Ordinary Smith form does not classify modules over this coefficient ring. | [Equivariant scope](EXTENSIONS.md#equivariant-calculations), [program](../articles/research-program-proposal/README.md) |
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

## How to choose the next calculation

Start with a precise general question above and the obstruction that prevents
a proof. Choose a case that discriminates between competing mechanisms, and
retain its exact source lattice and normalization. Record the result as a
finite test until a separate argument establishes the general statement.
F4, E8, additional representations and more parabolics are useful test cases
for this purpose, not automatic substitutes for a structural result.

The full-paper project remains on hold. This ledger, the supporting proof
notes and the API can develop independently of writing a new paper.
