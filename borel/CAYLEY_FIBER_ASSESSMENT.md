# Fiber program: theorem and literature assessment

Audit date: 2026-10-23. New group-by-group rank tables are paused.

## Decision

Do not treat centralizer reduction, generic affine-space fibers, or the
existence of polynomial nilpotent coordinates as a new research objective.
They are published constructions or elementary consequences of known
theory. Preserve the Borel API and the audited G2 ranks as useful explicit
work. Do not start the proposed F4 rank table or rank API merely to extend
coverage.

A more precise candidate project is to determine **which degenerations of
the translated centralizer map actually occur for representation-induced
Springer maps**, and what their fibers are. In type A centralizers, the
local map reduces to matrix powers; the remaining representation-dependent
question is which powers can occur and where. This is a proposed direction,
not a verified open problem of Kostant--Michor or a claim of novelty.

The rank reduction does not by itself solve fibers. Even an invertible
translated map leaves an incidence problem when the target centralizer
is larger than the source centralizer.

## 1. Source audit and corrected attribution

| Source and material inspected | What it establishes for this program |
| --- | --- |
| [Kostant--Michor, arXiv v3](https://arxiv.org/pdf/math/0109066), Sections 2.4--2.6, 4.6, 4.11--4.13 | Equivariance and centralizer compatibility; the unipotent isomorphism; preservation of the semisimple part; the translated centralizer map and its isomorphism under an invertible full differential. |
| [Sobaje, arXiv:1708.04153](https://arxiv.org/pdf/1708.04153), Example 2.4.2, Lemmas 2.4.3 and 3.1.1, Proposition 3.1.2 and its following discussion | Polynomial descriptions in type A; parameterization by regular-nilpotent centralizer data; homogeneous equivariant nilpotent maps. Section 3 explicitly includes characteristic zero. |
| [McNinch--Testerman](https://gmcninch.math.tufts.edu/assets/manuscripts/nilpotent-centralizers-and-springer-isomorphisms.pdf), introduction and stated centralizer/parameterization results | Relevant established structure of nilpotent centralizers and Springer isomorphisms. This paper studies isomorphisms, not arbitrary degenerate translated maps. |
| [Goodwin, Relative Springer isomorphisms](https://doi.org/10.1016/j.jalgebra.2005.04.014), publisher abstract | Existing work on induced B-equivariant isomorphisms on unipotent quotients. Full text was not retrieved; no uninspected theorem is used in the proofs below. |
| [Melnikov, arXiv:math/0312290](https://arxiv.org/pdf/math/0312290), introduction and square-zero setup | The variety of square-zero strictly upper-triangular matrices and its Borel-orbit geometry are established subjects. The counterexample below lands in this known geometry. |
| [Goddard's 2017 thesis](https://etheses.bham.ac.uk/7476/1/Goddard17PhD.pdf), Theorem 5.2.1 | Additional confirmation that the type A polynomial-isomorphism description is well known. |

The important correction is precise: KM Section 4.12 already defines

    Phi_t(z) = theta(tz) - theta(t)

on the full unipotent variety of Z_G(t), with values in its nilpotent
cone. Theorem 4.13 proves it is an isomorphism when theta(1)=0 and the
full differential d theta at t is invertible. Our use of U_t below is
the restriction to a maximal unipotent subgroup of that centralizer.
For semisimple G, theta(1)=0. A vanishing vertical factor and a singular
full differential are not interchangeable hypotheses.

Our earlier literature audit already listed KM 4.13. The oversight was
not making it an explicit baseline when proposing further centralizer
work. The G2 calculations are not retracted; their interpretation as
progress toward a new general theorem needed this qualification.

Searches covered translated Cayley maps, relative Springer isomorphisms,
equivariant nilpotent maps and matrix-power fiber geometry. No result
establishes exhaustive literature coverage or 2026 openness of the
remaining representation-realization question. No claim here is presented
as an explicitly posed open problem in KM.

## 2. Precise setup

Work over C. Let G be connected semisimple, rho a finite-dimensional
rational representation with nondegenerate trace form B_rho, and theta
its normalized trace-projection Springer map. Fix B=TU and t in T.
The actual character lattice of G is retained; no change of isogeny form
is implicit. Put

    h=theta(t),  H=Z_G(t)^0,  K=Z_G(h),
    A=Z_U(t),    C=Z_U(h),    psi(z)=theta(tz)-h  (z in A).

Here A is a subgroup of C. The symbol C in this display denotes a
unipotent group, not the coefficient field. KM's construction has a
larger domain, the full unipotent variety of H; its restriction is psi.
The latter maps A into a=Lie(A). Write c=Lie(C).

We study the morphism theta:tU -> h+u. This is not the Borel
characteristic homomorphism into cohomology, and its fibers are not the
fibers of the Springer resolution.

## 3. What the general fiber reduction actually proves

There are variety isomorphisms

    tU  = U x^A A,       [v,z] |-> vtzv^{-1},
    h+u = U x^C c,       [v,n] |-> Ad(v)(h+n).

Use the equivalence (v,z) ~ (va^{-1},aza^{-1}) in the first contracted
product, and the analogous action in the second. The Springer map is
the induced morphism

    [v,z] |-> [v,psi(z)].                                      (1)

**Proof and scheme scope.** Choose ordered root subgroup coordinates for
A and an ordered product of the complementary positive-root subgroups
as a section of U/A. Multiplication gives polynomial coordinates by
root height. In the first displayed normal-form map, a noncentral root
coordinate has linear coefficient e^{-alpha}(t)-1, up to the harmless
choice of left/right trivialization; a central root has coefficient 1
from z. Every other term uses lower heights. All these diagonal
coefficients are nonzero, so successive polynomial elimination gives
a polynomial inverse. For the Lie-algebra normal form the corresponding
coefficients are nonzero multiples of alpha(h) on the complementary
roots and 1 on c. This proves both isomorphisms without relying only on
a bijection of complex points. Equivariance proves (1).

Consequently, after conjugating a target to h+n with n in c, its fiber
is the incidence quotient

    { (v,z) in C x A : Ad(v)psi(z)=n } / A.                    (2)

This equality respects fiber schemes: it follows by base change of the
normal-form isomorphisms, with the quotient understood through the
polynomial section of C/A. It is a description, not a classification
of the incidence scheme on the right.

Two immediate consequences are worth recording as deductions, not
standalone novelty claims:

* If A=1 (t is root-regular), every nonempty fiber is a coset of C,
  hence an affine space of dimension dim C. The image is the U-orbit
  of h.
* If A=C and psi is an isomorphism, theta:tU -> h+u is an isomorphism.

There are two separate sources of complicated fibers: degeneration of
psi, and enlargement A subset C. Even when psi is an isomorphism,
(2) asks which C-conjugates of n lie in a; there is no general affine-
space conclusion from the reduction alone.

## 4. The invertible case is largely already understood

On h_der=Lie([H,H]), the bilinear form

    B_t(X,Y)=tr(rho(t) d rho(X) d rho(Y))

is H-invariant. Relative to B_rho it is a scalar c_j on each simple
ideal s_j. It gives the derivative at 1 of the translated map on that
ideal. Equivalently c_j is the regular root factor q_alpha(t) for any
root in that simple component of H.

In height-ordered coordinates on A, T-equivariance makes psi triangular:
the alpha coordinate is c_j z_alpha plus lower-height monomials.
It follows that

    psi:A -> a is an isomorphism iff every c_j is nonzero.      (3)

Sufficiency is triangular inversion; necessity follows by differentiating
an isomorphism at 1. The empty set of simple factors gives the trivial
map between points. This tests only the centralizer's unipotent part,
so it is weaker than invertibility of the full differential in KM 4.13.
It is nevertheless an elementary consequence of known structure, not a
substantial new classification by itself.

A zero c_j does not imply that the whole component maps to zero:
higher-degree terms can remain. Equation (3) also does not make the
full tU map invertible when C is larger than A.

## 5. Uniform type A normal form: what survives after degeneration

**Deduction.** Suppose h_der is a product of simple ideals sl_(m_j).
Identify their nilpotent cones by these chosen matrix realizations.
Compose the full translated KM map with exp on the nilpotent cone.
Its component on each factor has the form

    f_j(X_j)=sum_(k=1)^(m_j-1) c_(j,k) X_j^k.                 (4)

If it is not zero, let r_j be the first index with c_(j,r_j) != 0.
There is an equivariant polynomial automorphism Y_j=Q_j(X_j), preserving
the strictly upper-triangular subspace, such that

    f_j(X_j)=c_(j,r_j) Y_j^(r_j).                            (5)

Zero factors remain zero. This is a statement about the translated
centralizer map, not the removal of the incidence condition (2).

**Proof of (4).** First consider one factor. Choose a regular nilpotent
Jordan matrix J. Equivariance forces f(J) to commute with the centralizer
of J, and in particular with J. Its traceless matrix centralizer is
spanned by J,...,J^(m-1). Thus f(J)=sum c_k J^k. Both sides extend as
equivariant polynomial maps and agree on the dense regular orbit, so
they agree everywhere. For a product of factors use a tuple of regular
nilpotents. The product of their orbits is dense, and independent
conjugation in each factor shows that no dependence on the other
nilpotent variables remains. The central component is zero because the
translated map is nilpotent-valued.

**Proof of (5).** Write p(z)=c_r z^r u(z), u(0)=1, modulo z^m. Define
Q(z)=z u(z)^(1/r) modulo z^m using the finite binomial expansion. It has
linear coefficient 1 and a composition inverse modulo z^m. Substitution
in nilpotent matrices proves (5) and gives a polynomial inverse for Q.
All identities hold on the reduced nilpotent cone and on the affine
space of strictly upper-triangular matrices; their induced fiber schemes
are therefore isomorphic.

The polynomial description in (4) is standard type A centralizer theory,
consistent with the cited Springer-isomorphism literature. Extracting
the unit in (5) is elementary formal algebra. We did not locate this
exact translated-map formulation, but that is not evidence of novelty.

The actual representation enters through the c_(j,k), not through the
abstract existence of these normal forms. For example, choose Y_(j,k)
in s_j B_rho-dual to J_j^k among J_j,...,J_j^(m_j-1). Homogeneity of
(4) gives the coefficient extraction formula

    c_(j,k) = tr(rho(t) d rho(J_j)^k d rho(Y_(j,k))) / k!.     (6)

Indeed the degree-k part of theta(t exp(s J_j))-h is
s^k c_(j,k) J_j^k. Equation (6) follows by the defining trace pairing.
For k=1 this is the twisted trace scalar in (3). Higher coefficients
are extra representation data; the diagonal factors alone do not
determine them.

## 6. A controlled obstruction to the affine-fiber hope

This is one uniform counterexample family testing the proposed theorem,
not a new collection of rank examples. Let m be divisible by four,
G=SL_m, rho=V direct-sum V*, and t=i I_m. This is a central element of
the actual group. Its normalized Springer map is

    theta(g) = (g-g^{-1})_0 / 2,

where the subscript means traceless part; B_rho(X,Y)=2tr(XY).
For u=I+N in U, theta(t)=0 and

    theta(tu) = (i/2) N^2 (I+N)^{-1} = (i/2) Y^2,
    Y=N(I+N)^(-1/2).                                        (7)

The inverse coordinate change is

    N = Y^2/2 + Y sqrt(I+Y^2/4).

All series terminate in nilpotent matrices. Thus the restriction is
equivalent to the square map in every rank m divisible by four, even
though its derivative at u=1 is zero. For irreducible rho a central t
acts by a nonzero scalar and this particular cancellation cannot occur;
the direct sum and its different central characters are essential.

For m=4 write

    Y = [0 a d f; 0 0 b e; 0 0 0 c; 0 0 0 0].

The zero fiber has ideal

    I=(ab,bc,ae+dc) = (a,c) intersection (b,ae+dc).            (8)

Both factors on the right are prime. The first component is A^4; the
second is the irreducible quadric cone ae+dc=0 times the f-line, also
of dimension four. Their intersection is A^3. The fiber is reduced
and reducible, hence not an affine space. No conclusion about radicality
of the square-zero ideal in arbitrary rank is inferred from this check.

[The symbolic verifier](tools/verify_cayley_fiber_assessment.py) checks
the full 8-dimensional trace representation for SL4, the polynomial
coordinate change and its inverse, and the ideal decomposition. The
all-rank assertion (7) has the separate formal proof above. Square-zero
triangular varieties themselves belong to the established Melnikov
literature; this is an explicit occurrence in our Springer-map family.

## 7. What remains worth asking, and the stopping rule

The strongest precise candidate in type A centralizers is:

> Determine the possible orders r_j (including the zero-map case) for
> translated Springer maps induced by actual representations, and give
> uniform character criteria for the strata where those orders occur.

Equation (6) is an extraction procedure, not such a classification.
The positivity and integrality of representation multiplicities and the
central character values constrain which cancellations are realizable.
Abstract classification of equivariant polynomial maps does not answer
that representation-realization question. Applying (5) would then
transfer known matrix-power geometry instead of recomputing it.

For centralizers of other types, one must use the established homogeneous
nilpotent covariants in place of matrix powers. A claim that the first
nonzero degree alone classifies their fibers has not been proved here
and is not adopted as a conjecture without testing the possible moduli.

**Go/no-go decision:** stop expanding group rank tables and postpone
the F4 rank API. Retain (1)--(6) as explanatory deductions. A further
research step should first audit the higher trace-index/covariant
literature for the representation-realization question, then attempt
one arbitrary-rank classification with a named representation family.
Success means explicit necessary and sufficient conditions for r_j or
a uniform obstruction, not more isolated values. If those conditions
are already known, record the application and stop this branch of the
program unless it produces new fiber geometry. No claim of a new
flagship theorem is justified by this assessment.
