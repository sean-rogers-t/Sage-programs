# Extending the Springer/Borel program: exploration and next steps

22 October 2026. The Spin parity problem is shelved. This report records
the requested exploration of other cohomology theories and affine groups.
The linked notes contain statements, proofs, sources and exact checks.
They are working research deductions, not external peer review or claims
of originality. The full-paper project remains on hold. The original
exploration changed no production API. The subsequent
[classical equivariant backend](CLASSICAL_EQUIVARIANT.md) now implements
the first H/K maps and a [symplectic cokernel reduction](SYMPLECTIC_K_FITTING.md);
it does not make every proposed extension a completed backend.

## 1. What we are extending

In the notation of the Kumar papers, the literal integral source is

\[
 A_{\rho,L}=\operatorname{Sym}_{\mathbb Z}(\Lambda_\rho)^{W_L},
 \qquad \theta_\rho^*A_{\rho,L}\subset R(L),\qquad
 \xi_\rho^P=\beta^P:A_{\rho,L}\longrightarrow H^{2*}(G/P;\mathbb Z).
\]

The inverse Springer pullback has already transported the polynomial
characters to these polynomial coordinates. We now keep that source
explicit while changing the target characteristic theory. For K-theory
the canonical extension is simply

\[
                 f\longmapsto\kappa_P(\theta_\rho^*f)\in K^0(G/P).
\]

It uses the actual virtual characters. Ordinary integrality of a
rationally normalized generator is insufficient. For oriented and
quantum theories the linked notes specify the lift, rather than assuming
that every familiar ambient Schubert theory automatically supplies one.
The affine construction uses evaluation and a separately chosen level;
it does not require an infinite-dimensional trace projection.

Three representation-data levels are now visible:

| Target | Data retained by the specified map |
| --- | --- |
| Ordinary flag cohomology | The linear source lattice Lambda_rho |
| Nonequivariant K-image | Lambda_rho and the normalized character class tau_rho in (m_G/m_G^2) tensor Q |
| Equivariant K-image | The full character/torus Springer map; the earlier quotient coordinates and root Cayley quotients reappear |

The normalized primitive class tau is known K-theoretic representation
data. The contribution being explored is its effect on these specified
embedded image lattices and their specializations.

## 2. Results obtained in the exploration

| Direction | Written result | What remains |
| --- | --- | --- |
| Finite K-theory | Same integral source kernel as ordinary Borel, strict support filtration, associated graded equal to the old image, and unchanged total index | Ambient Smith extensions between degrees; they can change |
| Complete rank-two benchmark | Every SL3 representation and every parabolic has explicit K Smith factors from its source homothety and one cubic-index gcd | Higher-rank classification using the full primitive vector |
| Universal oriented lift | For even Z-torsion-free coefficients, the specified image is a filtered copy of the old one; the universal zeroth Fitting ideal is the old integer index | Higher Fitting ideals, and modular specialization kernels described by Tor_1 |
| Defining type A | For every partial flag, a change of Chern coordinates carries the complete old image to the specified oriented image | Other representations do not obey this automatic embedded comparison |
| Symplectic projective flags | A complete height-dependent modular image; an arbitrary-rank K-cokernel recurrence and closed largest factor/multiplicity | Closed multiplicities of all intermediate K factors; multigenerator height problems |
| Equivariant K-theory | A deformed source family and a determinant-ideal formula for every finite group/parabolic, directly involving the old Cayley root factors | General integral source freeness and global higher Fitting ideals |
| Horizontal local factors | Exact DVR Smith patterns at a single reflection wall and at scalar collapse; the latter uses ordinary Betti degrees | Multiple walls, larger stabilizers and arithmetic coefficient primes |
| A complete global equivariant family | For tensor powers of defining SL_m, every parabolic cokernel over Q[X] splits into cyclic modules with character powers indexed by Schubert lengths, giving every higher Fitting ideal | General representations and actual integral source factors remain separate |
| Integral equivariant family | Defining Sp_(2r)/P1 has a complete source kernel and determinant ideal 2^r product_(i<j)(s_i+s_j)^2 | Its global higher Fitting ideals |
| Quantum | Complete flat lifts for the defining SL and Sp projective families, with faithful q=0 recovery | Nontrivial partial flags with source relations and nonunique higher lifts |
| Affine | Closed all-rank A/C primitive divisor-power contents, with an evaluation/level convolution of the old finite Smith factors | More loop generators, other finite root types and additional affine parabolics |
| Folding | Compatible trace maps for stable representations, with exact integral source defects in A/C, E6/F4 and D4/G2 examples | A broader restriction/source comparison; central/rotation and flag-image questions remain |
| Plain/twisted loops and jets | Apply the loop functor to finite theta; stable representations preserve twisted loops; finite jet degrees stay the old degree and their critical determinant is its m-th power at the constant jet | These formal consequences do not give right-Iwahori descent or central/derivation directions; jet flag cohomology alone adds no new topology |

Precise statements and proofs:
[finite comparison](OTHER_THEORIES_FINITE.md),
[complete SL3 formulas](OTHER_THEORIES_SL3.md),
[defining type A](OTHER_THEORIES_TYPE_A.md),
[height and symplectic K-arithmetic](OTHER_THEORIES_HEIGHT.md),
[primitive character data](OTHER_THEORIES_COTANGENT.md),
[equivariant K/Cayley connection](OTHER_THEORIES_EQUIVARIANT_K.md),
[quantum lifts](OTHER_THEORIES_QUANTUM.md),
[affine images](OTHER_THEORIES_AFFINE.md),
[folding](OTHER_THEORIES_FOLDING.md).

## 3. Three concrete formulas

The strongest connection back to the earlier exceptional geometry is
the equivariant K determinant. For root length a, define

\[
 J_{\rho,a}=\frac{\prod_{\alpha>0,\,\mathrm{length}=a}\alpha(\theta_\rho(t))}
 {e^{\rho_a}\prod_{\alpha>0,\,\mathrm{length}=a}(1-e^{-\alpha})}.
\]

It is a regular invariant Laurent polynomial with value one at 1.
For n=|W/W_L|, N_a positive roots of length a and N_(L,a) in the Levi,
the rational determinant ideal is

\[
 \operatorname{Fitt}_0(C_{\rho,T}^P)=
       \left(\prod_a J_{\rho,a}^{\frac n2(1-N_{L,a}/N_a)}\right).
\]

The source family uses the already-computed invariant values
f_i(theta_rho(t)). On E6/P1 its determinant is J_rho^6; on E7/P7 it is
J_rho^12. On the full flag, the root part of the full Cayley determinant
is J_rho^2. The Cartan Jacobian is separate. This supplies a direct
connection between the image program and the earlier Kostant--Michor
calculations, rather than only a thematic resemblance.

For defining Sp_(2r)/P1=CP^(2r-1), a characteristic-two orientation of
height s has complete image

\[
                  R[z]/z^{\lceil2r/2^s\rceil}.
\]

In integral K-theory, the largest factor is instead
2^(4(r-1)-s_2(r-1)) for r>=2, occurring twice; r=1 has factor two.
Here s_2 counts binary ones. Thus the same source detects formal-group
behavior and integral extension data that ordinary cohomology misses.

For the affine Grassmannian of SL_n, the content c_(n,d) of the d-th
power of its primitive divisor satisfies

\[
       v_p(c_{n,d})=v_p(\lfloor d/p^a\rfloor!),\quad p^a<n\le p^{a+1}.
\]

The symplectic analogue has the same odd-prime formula with n replaced
by 2n, while its two-primary valuation is v_2(d!). Every prime eventually
appears. This is arithmetic inside established affine cohomology, not a
definition of an affine Springer map.

## 4. What is published, deduced, checked and still open

**Published foundations.** Kumar's polynomial-character construction,
Rogers's explicit trace formula, characteristic K-maps and localization,
Hodgkin/Brylinski--Zhang primitive data, oriented flag presentations,
formal-group height, quantum flag rings, Bott/Lam affine Hopf models
and both cohomological and K-theoretic Peterson correspondences. The
proof notes give primary-source links and theorem locators. We did not
discover these ambient theories, and the K-theoretic Peterson theorem
is no longer merely conjectural.

**Deductions written here.** The specified image comparisons, source
counterexamples, SL3 classification, symplectic exponent and height
applications, determinant/local Smith formulas and affine content
formulas. Many use short formal consequences of established results.
Their mathematical arguments and hypotheses are given; their literature
priority has not been established. A search failing to locate the exact
statement is not evidence of originality, nor are these asserted to be
questions posed by the cited authors.

**Independent exact corroboration.** Six complete full-flag K-image
matrices agree between Laurent Demazure and native Todd/BGG methods.
The SL3 classification is checked on 41 actual representations and all
three proper flags. The symplectic K recurrence is checked against
projective Euler characteristics through projective dimension 63.
Height checks use both p-typical logarithms and a separate supersingular
elliptic curve. Affine formulas have 6,231 type A and 4,824 type C
arithmetic checks plus independent A1/C2/C3 Chevalley calculations.
Equivariant determinants have independent rank-two tests and exceptional
reflection counts; the latter are not full exceptional K-image tables.
Folding checks compare actual weights, Laurent coordinates and integer
source lattices. Source hashes bind these reports to their scripts.
The [combined verifier](other-theories-verification.json) also checks
the all-parameter integral block reductions behind the SL3 theorem.

**Still open here.** A general embedded K/CK Smith/Fitting classification,
the symplectic intermediate factors, larger chromatic sources, global
equivariant higher Fitting ideals, geometric quantum lifts beyond the
benchmarks, unlocalized Peterson image recovery, other affine Hopf
contents and a central/rotation-compatible Kac--Moody Springer analogue.
The plain and twisted pointwise loop maps already exist by base change;
their missing affine-flag source/descent is stated explicitly.

## 5. Recommended program and API path

For a self-contained introduction to the geometry, notation, integral
results and these new targets, see the
[classical expository paper](../articles/classical-equivariant-exposition/README.md).

### Current focus: equivariant images in the classical A-D families

The first implementation is available: all-parabolic type A input
evaluation in H/K with literal integrality checks, complete source
generators on full flags and P1, and the type C projective image/cokernel.
The latter now has a proved all-rank doubled half-size presentation,
an extra R/(2) summand in odd rank, a top Fitting parity formula,
complete closed Fitting ideals through rank nine and global cyclic
decompositions in ranks one to three.
Integral Schur cyclicity now gives uniform partition-indexed generators
for every C ideal using initial-column minors. Integral Hankel balancing
now reduces each entire convolution to a single product of equal or
adjacent half-block ideals, in every rank and including at two.
The [ordinary Hankel theorem](SYMPLECTIC_K_ORDINARY_SMITH.md) now gives
every Smith multiplicity in this symplectic projective family by a closed
binary-digit formula. It applies Euler's finite binomial continued fraction
and known orthogonal Hankel diagonalization over Z_(2); the application
has no priority claim. Primary arithmetic after equivariant moment
substitution and full type A image modules remain targets.

The user's subsequent choice is to prioritize uniform classical-series
statements. This replaces the exceptional-first sequencing below; the
exceptional, quantum and affine directions remain in the inventory.
The Spin central-parity classification remains shelved.

Fix the defining/vector representation first, with its normalized trace
Springer map and literal source A_(rho,L)=Sym_Z(Lambda_rho)^(W_L).
Use SL_m in type A, Sp_(2n) in type C, and Spin_(2n+1), Spin_(2n) in
types B and D. The latter vector representations factor through SO and
have faithful differential, rather than being faithful representations
of Spin. Keep the full Spin character lattice in the coefficient ring;
changing to SO is an explicit isogeny comparison, not a silent convention.

| Family | Established starting point here | Next target, not yet asserted |
| --- | --- | --- |
| A_(m-1) | Defining ordinary/oriented all-parabolic comparison; rational equivariant K isomorphism for the defining representation and a complete tensor-power family | Actual integral equivariant source image, higher Fitting ideals and specialization for Grassmannians, then general partial flags |
| C_n | Complete integral projective source/kernel; doubled floor(n/2) cokernel blocks with an extra R/(2) for odd n; every Fitting ideal is one balanced product of Schur/initial-column ideals; complete closed lists through rank nine; all ordinary K Smith multiplicities have a closed binary-digit formula | Primary structure after equivariant moment substitution, then larger isotropic Grassmannians |
| B_n | General rational determinant theorem and earlier ordinary Borel-image calculations | Integral vector-source image on the odd quadric Spin_(2n+1)/P1, then isotropic Grassmannians |
| D_n | General rational determinant theorem and earlier ordinary Borel-image calculations | Integral vector-source image on the even quadric Spin_(2n)/P1, retaining both middle Schubert classes, then isotropic Grassmannians |

For equivariant ordinary cohomology, study the coefficient-span image
J_(rho,T)^H=S_Z beta_T^P(A_(rho,L)), with
S_Z=Sym_Z(X(T)), whenever the specified integral characteristic lift
has been established. For equivariant K-theory the unambiguous literal
character construction is

\[
 J_{\rho,T}^{K,P}
 =R(T)\,\kappa_T^P(\theta_\rho^*A_{\rho,L})
 \ \subset K_T^0(G/P).
\]

The common goal is a uniform presentation of this embedded image and
its cokernel, including higher Fitting ideals and the kernel of
augmentation to the ordinary image. The rational deformed source and
zeroth determinant are already proved; recomputing those alone is not
the next milestone. A multivariable coefficient ring does not in general
admit Smith normal form. Integer Smith factors belong to the specialized
finite lattice calculation; DVR Smith factors require a specified prime.

Implementation order: (1) package the defining type A equivariant map
with its integral source and independent restrictions/Schubert checks;
(2) add the already-proved symplectic projective model and attack its
higher Fitting ideals; (3) construct B/D quadric targets with the full
Schubert or spinor-bundle information; (4) extend to isotropic
Grassmannians and arbitrary partial flags. Type D maximal isotropic
components must be specified separately. Tensor and exterior powers,
adjoint and spin representations are later comparison inputs, not
assumed to have the defining representation's image.

Completion of the first stage means a documented classical equivariant
API returning source lattice, coefficient ring, image presentation,
cokernel/Fitting data and augmentation map, with independent exact
checks. Neither a production backend nor a uniform B/D equivariant
integral-image theorem is claimed by this roadmap update.

### Broader exploration inventory and earlier sequencing

1. **Package primitive character data and a K/CK image backend.** Start
   with the completed SL3 formulas and G2 source/image cases. Return
   source lattice, normalization, actual K matrix, support filtration,
   total factors and associated-graded comparison separately. For CK,
   return Fitting ideals rather than a misleading Smith list over Z[beta].
   Completion means a documented API with independent checks, including
   the known G2 noncontainment and same-lattice SL3 counterexample.
2. **Use the old exceptional geometry to study equivariant K-images.**
   Start with G2 and E6/P1, then E7/P7. The zeroth determinant and its
   exponents are already settled above; the new target is higher Fitting
   data and local stabilizer strata, with the actual integral source
   retained. Recomputing only J^6 or J^12 would repeat completed work.
3. **Package the affine A/C content and evaluation/level formulas.** This
   is a feasible extension that reuses all old finite parabolic Smith
   tables. Keep cup cohomology separate from Peterson Pontryagin homology.
   Then choose a higher loop generator or another root type to test a
   structural question.
4. **Quantum partial flags after the two projective baselines.** Choose
   a small Grassmannian or minuscule exceptional target where higher
   generator degrees meet quantum degrees. Specify the lift first, then
   test its kernel, flatness and integral localization losses.
5. **Chromatic multigenerator sources and twisted/affine group maps.**
   The first is a concrete medium-term problem using the universal MU
   embedding. The central/rotation extension remains a foundational
   longer-term direction: the pointwise loop map is already available,
   while an affine-flag source and topology/completion need work.
   Existing generalized Kac--Moody Schubert theories are starting inputs,
   not missing theories for us to invent.

The current classical sequence above takes precedence over this earlier
broader ordering. The Spin classification and the full-paper project
remain paused in the research sense; their files and unresolved questions
are preserved.
