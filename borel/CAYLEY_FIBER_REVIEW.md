# Review: what the Springer fiber program is trying to prove

Subsequent progress: the [binary critical-equation criterion](BINARY_BRANCHING_FIBERS.md)
now covers two further E7 Weyl families (22 standard embeddings), raising
coverage to 123/127. This review and its eight-control audit retain the
101/127 checkpoint at which the recommendation below was made.

Date: 2026-10-24. This is a proof and research-direction review of the current
character/critical-scheme strand. It does not re-audit the older integral-image
program or claim an externally reviewed classification. No new exceptional
family is promoted here, and no literature priority is established.

## 1. The mathematical destination

Let G be simple, complex and simply connected, rho an actual nontrivial
finite-dimensional representation, and normalize the Springer map by

    B_rho(theta_rho(g),X)=tr(rho(g) d rho(X)),
    B_rho(Y,X)=tr(d rho(Y) d rho(X)).

For a proper Levi K, write r=rank G and

    Z_K^circ={h in z(k): Z_G(h)=K},
    F_(rho,h)^reg=theta_rho^(-1)(h) intersect G_reg.

Regular means that the source centralizer has dimension r. The target need
not be regular: if K is larger than a maximal torus, it is not. These are
fibers of a map G -> g, not the varieties of Borels in the Springer resolution.

The aim is to determine from the full branching character chi_rho|K when
the generic F_(rho,h)^reg is empty, smooth, or has excess dimension. A later
refinement would count regular K-orbits when the transverse scheme is finite
and describe exceptional targets and multiplicities. Neither refinement is
implicit in a generic smoothness proof.

The connection to the earlier papers is at the map itself: the torus
weight formula supplies theta_rho, and the character supplies its trace
gradient. The Borel characteristic map to flag cohomology is not used in
the argument below. There is no deduction here about Smith factors or an
integral Borel image. The [published-map inventory](SPRINGER_GEOMETRY.md)
records the normalization and attribution.

## 2. A uniform criterion, and which part is already standard

On the finite central cover S x T -> K, where S=[K,K], choose fundamental
characters x_1,...,x_ell of S and coordinates z_1,...,z_s of T. Here
ell+s=r. Write phi(x,z)=chi_rho|K, D_j=z_j partial/partial z_j, and retain
the actual scheme and target map

    X=V(phi_x1,...,phi_xell) in A^ell x T,
    mu:X -> A^s,  mu=(D_1 phi,...,D_s phi).

The target coordinates are c_j=B_rho(h,H_j) in the chosen central
cocharacters. This is an invertible linear identification of the central
target space. The corrected potential is locally

    W_c=phi-sum_j c_j log z_j.

Its differential is algebraic even though its primitive uses local
logarithms. The quotient critical scheme is exactly mu^(-1)(c), including
nilpotents. At a regular group point, the
[critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md) gives

    completed O_(F,g) = completed O_(mu^(-1)(c),q)
                         [[u_1,...,u_(dim K-r)]].

The finite central action must still be divided out. It acts freely on
the central torus factor, so the quotient cover is etale. Smoothness and
local dimensions descend, but naive counts of points need not.

**Proposition (standard generic-smoothness reformulation).** For any
finite-type complex scheme X and morphism mu:X -> A^s, the following
conditions are equivalent:

1. The geometric fibers over some dense open subset of A^s are smooth,
   allowing empty fibers.
2. The closure of mu(Sing X) is a proper subset of A^s.

Here Sing X is the locus where X itself is not smooth over C, with its
underlying support. X is not replaced by its reduction. If X is empty
the equivalence is immediate. If s=0, it says precisely that X is smooth.

Proof of 2 implies 1: remove the closure of the bad image. The remaining
source is smooth. Apply characteristic-zero generic smoothness to each
dominant component and remove images of nondominant components. Finitely
many components give one dense target open.

Proof of 1 implies 2: shrink that target open using generic flatness.
The restricted morphism is flat and of finite presentation with smooth
geometric fibers, hence smooth. Its base is smooth over C, so its source
is smooth over C. This open meets no image of Sing X, proving 2.

For a contributing irreducible component X_i, after shrinking the target
open further the corresponding local group-fiber dimension is

    dim K-r + (dim X_i-s).

Nonemptiness requires at least one dominant component. Finiteness of the
transverse fiber requires dim X_i=s for every contributing component.
These are additional claims, not consequences of smoothness alone.

The proposition follows from existing generic flatness and smoothness
theory; it is not a new Springer theorem. Its value is to identify the
precise representation-theoretic task: prove that singular incidence
components cannot dominate the central target space.

Two cautions matter. Sing X is not the relative critical locus of mu on
X_sm; that second locus can still give exceptional target values. Also,
maximal minors of a chosen list of equations need not describe Sing X
when the equations are redundant or there are smooth excess components.
Full Jacobian rank is a sufficient certificate on its open set; its
complement may be larger than the true nonsmooth locus.

## 3. Proof dependency audit

The written arguments below were read for their scheme, target and
boundary steps. The last column states the limits of this review rather
than treating an old successful report as a new proof.

| Result | Essential argument and verdict | Remaining qualification |
| --- | --- | --- |
| [Critical presentation, Theorem A](GENERAL_CAYLEY_FIBER_CRITICAL.md) | Equivariance places the fiber scheme in K. The Maurer--Cartan trace identity identifies its ideal with dW. Smooth quotient coordinates give the formal product. These steps withstand this review. | The product requires a smooth quotient map; it is not asserted at nonregular points or arbitrary isogeny forms. |
| [Affine branching](AFFINE_BRANCHING_FIBERS.md) | At a nonsmooth coefficient point of local dimension d and derivative rank p, d+p<s. The entire affine-variable target image therefore has dimension below s. | Target-image exclusion, not source reducedness, is the key. Dominance and dimensions must be checked separately. |
| [Quadratic branching](QUADRATIC_BRANCHING_FIBERS.md) | Constant-rank strata have full affine solution fibers of dimension ell-k. Symmetry gives rank(B^t restricted to ker Q)=rank[Q,B]-k. At a nonsmooth incidence point this forces d+p<s. | The proof includes rank changes. It cannot be extended to all cubic potentials just by retaining the same formula. |
| [Unit multiplier elimination](MULTIPLIER_BRANCHING_FIBERS.md) | Solving the paired gradients gives an isomorphism of quotient rings. The chain rule retains derivatives of the central coefficient matrix. | A determinant nonzero in C(z) is insufficient; global elimination needs a Laurent unit. |
| [Nonunit multiplier](E7_MULTIPLIER_BOUNDARY.md) | The entire multiplier boundary is excluded by its target image; translated-subtorus zero divisors have fixed annihilators. The remaining scalar gradient has smooth source. | The divisor and boundary hypotheses are substantive. E7's application does not prove them for arbitrary characters. |
| [Spectral pencil](SPECTRAL_BRANCHING_FIBERS.md) | Etale eigensplitting gives smooth nonzero-level quadrics. The eta=0 collision target factors through a smaller-dimensional spectral base. | Missing boundary control, persistent collisions or Jordan blocks invalidate this route. |
| [Symmetric pencil](E7_SYMMETRIC_PENCIL.md) | Congruence, the Schur complement and the nonzero kernel coordinate give the exact corank-one scheme. The three dimension bounds cover the discarded strata. | A dominant ordinary branch proves nonemptiness, not absence of bad branches; the bounds supply that second fact. |
| Classical and E6 coverage | The classical spectral/trace argument and E6 dependency/conjugacy assembly were reviewed against their stated conclusions. | This review does not rerun every historical root enumeration or independently reconstruct every E6 seed proof. Their linked exact audits remain dependencies. |

No failure of the main critical presentation or the seven criterion
arguments above was found. This is an internal mathematical review,
not formal proof verification or external peer review.

The foundational Sage checker was rerun, including the nonregular SL2,
ramified PGL2 and Sp4 group-versus-quotient controls. The native E7 audit
was also rerun: it reconstructs the character from roots, checks the
collision polynomials, lower-rank identities and target normalization.
The [review controls](tools/audit_cayley_fiber_review.py) independently
check the distinctions in section 2 using exact polynomial ideals.

## 4. Current conclusions and actual limitations

| Scope | Current conclusion | What this does not settle |
| --- | --- | --- |
| A-D, adjoint, proper Levis | Generic regular fibers are smooth or empty; nonempty ones have dimension dim K-r. Block criteria classify nonemptiness in the recorded rank ranges. | Whole fibers, special targets and orbit counts. See [classical theorem](CLASSICAL_ADJOINT_FIBERS.md). |
| E6, adjoint, proper Levis | Generic smoothness or emptiness in every proper Levi stratum. | A complete nonemptiness/dimension table is not supplied by the smoothness theorem alone. See [E6 assembly](EXCEPTIONAL_QUARTIC_REDUCTION.md). |
| E6 contact Levi, derived type A5 | Nonzero central targets have smooth regular fibers of dimension 32, exceeding dim K-r=30. | A universal no-excess theorem is false. See [contact proof](EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md). |
| E7, adjoint | 101/127 proper standard Levis covered; 26 embeddings in four Weyl families remain. | Neither all E7 nor all exceptional groups is completed. |
| G2, F4, E8 | Earlier selected results remain in the ledger. This review certifies no new all-Levi classification for them. | Degree, regularity and integral-image results elsewhere do not fill gaps in this fiber theorem. |

The four unresolved E7 families are A1^3 (10 embeddings), A2+A1^2 (12),
A2^2+A1 (3), and A2+A1^3 (1). The counts enumerate standard Levi
subgroups, not conjugacy orbits inside a Springer fiber or map degrees.

Existing counterexamples delimit the goal. The same cubic monomial
support can produce smooth fibers or persistent nilpotents; see
[the explicit reductive representation](SPECTRAL_BRANCHING_FIBERS.md#3-a-sharp-obstruction-the-graph-alone-is-insufficient).
That example is not a simple-group adjoint counterexample. The actual
E6 contact example refutes no-excess but does not refute smoothness.
The nonregular and PGL2 examples refute dropping the quotient hypotheses.
These distinctions rule out several tempting overstatements.

## 5. The companion certificate: valid route, limited payoff

For the E7 A2+A1^3 model let I be the ideal generated by the three
relative derivatives of psi and the Cayley cubic equation h, in the
six variables x,y,w,A,B,C. Let J contain I and the 4 by 4 minors of
that four-equation Jacobian. The locus V(J) contains every point where
full rank fails. A bound dim V(J)<=1 suffices: outside it the incidence
scheme is smooth of dimension two, and a one-dimensional locus cannot
dominate the two-dimensional target.

The stored six polynomials have leading monomials

    x^5, y^5, w^2 B^4, w^3 C^3, A^6, B^5 C.

If each lies in the original 19-generator J, then in(J) contains the
displayed monomial ideal. Its radical is (x,y,A,wB,wC,BC), of dimension
one. The inequality dim V(J)<=1 follows in the correct direction.
No assertion that the six polynomials form a Groebner basis is needed.

The archive contains 242 computed basis members. The six selected
members have degrees five or six and 37-59 terms each. Their size does
not bound the size of ideal-membership witnesses. The native audit
checks the original generators, archive hash and leading supports, but
does NOT establish the missing memberships. Modular dimension agreement
is corroboration and does not supply the missing characteristic-zero
identities under the repository's independent-certificate standard.

The finite torus map to the Cayley cubic is etale away from four
inversion-fixed points. At those points all trace-coordinate differentials
vanish and the central target is zero. Thus the proposed descent argument
correctly excludes them for generic targets. A completed certificate,
together with the already checked dominance witness, would prove the
companion's generic regular-fiber dimension 12 and raise coverage to
102/127. It would not prove the other three families or a uniform theorem.

Assessment: the proof route is mathematically adequate, but a further
unbounded Groebner/lift run is not the next research milestone. A future
certificate attempt should retain explicit rational combinations of the
original generators (or a checkable proof DAG), verified separately.
Smaller members with the same dimension-forcing support would also suffice.
This review leaves PART-07 conditional and preserves its false promotion flags.

## 6. Literature attribution checked in this review

| Source inspected | Relevant published content | Relation to this program |
| --- | --- | --- |
| [Kostant--Michor](https://www.mat.univie.ac.at/~michor/kostant.pdf), sections 2.1, 2.4, 2.6, 3.2-3.3 and 4.12 | Trace projection, equivariance/centralizers, subgroup restriction, harmonic localization, degree reduction and the nilpotent-part map. | Supplies the map and background. The strata studied here are not a new definition of their map or a solution of a verified named open problem. |
| [Steinberg](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf), Theorems 1.2, 1.4, 1.5 and section 6 | Regular classes, quotient coordinates and independence of fundamental-character differentials in the simply connected case. | Supplies regular lifts and smooth quotient directions; the representation-dependent critical equations are additional input. |
| [Stacks: generic flatness](https://stacks.math.columbia.edu/tag/0529), [smooth morphisms, Lemma 29.35.3](https://stacks.math.columbia.edu/tag/01V4), and [Vakil, generic smoothness](https://math.stanford.edu/~vakil/216blog/FOAGjul2724public.pdf) | Generic flatness, the flat-plus-smooth-fibers criterion and characteristic-zero generic smoothness. | Section 2's equivalence is a standard consequence, not an originality claim. |
| [Tulczyjew--Urbanski](https://arxiv.org/pdf/1006.4183), section 9, Definition 1, Theorem 1 and Proposition 4 | Regular generating families, Hessian rank, clean intersection and their immersed Lagrangian images. | Closely related established geometric language. It does not supply the missing algebraic scheme-level branching/boundary verification. |

For the last connection, on X one has d phi=sum mu_j dlog z_j. Thus
(x,z) -> (z,mu) is the natural generating-family map to T* T. This
explains why mixed Hessians occur. It does not make X smooth or eliminate
nilpotents. The algebraic target-image argument must still be proved.

Targeted searches on 2026-10-24 included generalized Cayley adjoint
fibers, Cayley Levi smoothness, character critical loci and regular
generating families. No exact published all-simple-adjoint/all-Levi
classification was identified. That is a limited search result, not a
claim that the conjecture is a recognized open problem or that our
specialized results are original. Affine Springer fibers, compact real
height functions and Higgs-bundle Cayley correspondences address
different objects and were not used as substitutes.

## 7. Recommendation and next theorem target

Keep the following as a **conjecture of this program**, not a theorem:

    For every simple complex simply connected G and every proper Levi K,
    the adjoint Springer map has smooth, possibly empty regular-source
    fibers over a dense open subset of Z_K^circ.

It includes A, B, C, D, G2, F4 and E6-E8. It excludes the zero-target
stratum K=G, does not assert nonemptiness, and deliberately allows excess
dimension. By section 2 its equivalent technical target is

    closure(mu(Sing X_(Ad,K))) != A^s.

**Next milestone:** seek a root/central-charge argument proving this
nondominance for the simply laced adjoint all-A Levi class, after exact
unit-pair elimination. The existing norm theorem and quartic reduction
are the starting inputs. For E6-E8, the residual models have at most
four affine variables; that is a finite reduction, not the missing proof.

The route is to identify the possible rank-defect strata from those
branching constraints and find target relations on each stratum, retaining
the central derivatives. The E6 excess example must remain allowed; the
persistent-collision control must not be incorrectly admitted. Use the
four remaining E7 families as tests of a stated lemma, not four unrelated
invitations to compute.

A successful next step is either (i) one proved sufficient branching
criterion covering at least two unresolved Weyl families, with exact
independent checks of its applications, or (ii) an actual adjoint
counterexample identifying the failed hypothesis. If neither is obtained,
record the specific obstruction; do not extend the table automatically.
An abstract polynomial counterexample restricts a proposed lemma but
does not by itself refute the group conjecture.

The six-membership completion is a secondary certification task. E8
table expansion, a new public API and a full paper remain deferred.
The immediate purpose is now explicit: explain generic smoothness from
representation structure, before undertaking orbit counts or special-fiber
singularity classifications.
