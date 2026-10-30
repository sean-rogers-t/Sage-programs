# The active target: a theorem for all simple groups

24 October 2026. This supersedes the recommendation to complete the E7
table. The user explicitly requested a direct attempt at a general theorem
and rejected further automatic case-by-case Levi calculations. No new
uniform smoothness theorem is claimed in this note. The 25 October
Schur-pencil continuation near the end now proves quasifiniteness on the
allowed source open for the second authorized E8 test candidate. Together
with the earlier boundary exclusions, it settles generic-target smoothness
for that candidate. This does not settle the all-simple conjecture or
restart the Levi census; earlier failed-attempt status paragraphs below
record the state at those stages.

## Exact conjecture

Let G be any simple simply connected complex algebraic group, and let
theta=theta_Ad be normalized by

    B(theta(g),X)=tr(Ad(g) ad(X)),
    B(X,Y)=tr(ad(X) ad(Y)).

For every proper Levi K, put

    Z_K^circ={h in z(k): Z_G(h)=K}.

Conjecture: there is a dense open U_K in Z_K^circ such that the actual
scheme theta^(-1)(h) intersect G_reg is smooth for every h in U_K.

Empty fibers and positive transverse excess are allowed. This includes
every finite simple type. It does not concern arbitrary representations,
the zero target with K=G, all special targets, nonregular source points,
or the Springer resolution. It is our program conjecture, not a question
being attributed to Kostant--Michor.

The [ledger](RESEARCH_LEDGER.md) and [proof status](PROOF_STATUS.md) retain
the existing results. Their coverage counts are unchanged. Completing a
finite table is not the active objective.

## A direct intrinsic attack, without individual branching expansions

Write k=s direct-sum z, where s=[k,k] and z=z(k). The ambient Killing
form B restricts nondegenerately to k and makes this sum orthogonal.
Let p_s and p_z denote its projections. Since theta(K) lies in k, define
the actual incidence scheme on the K-regular locus by

    Y={g in K_reg: p_s(theta(g))=0},
    nu:Y -> z,  nu(g)=p_z(theta(g)).

For h in Z_K^circ, nu^(-1)(h) is exactly the regular-source Springer
fiber. Indeed a solution lies in K by equivariance, and
Z_G(g) is contained in Z_G(theta(g))=K. Thus regularity in K and in G
coincide at such a solution. This is the existing centralizer reduction,
not a new conclusion about smoothness.

Define A_g:k -> k by left-translating the differential of theta|K:

    B(A_g X,Y)=tr(Ad(g) ad(X) ad(Y)).                    (1)

This identity follows by differentiating the trace definition along
g exp(tX). At every underlying point of the incidence scheme Y,

    B(A_g X,Y)-B(A_g Y,X)
      =tr(Ad(g) ad([X,Y]))
      =B(theta(g),[X,Y])=0.                             (2)

Thus A_g is symmetric with respect to B on k. Equivariance also gives

    (Ad(g)^(-1)-1)k subset ker A_g.                     (3)

These are precisely the conjugacy directions in the left tangent frame.
For g regular in K they have dimension dim K-rank G. Consequently A_g
induces a symmetric form on the rank(G)-dimensional vector-space quotient
k / (Ad(g)^(-1)-1)k. This quotient, rather than the whole k, contains the
transverse differential. No nonsingularity of this form is asserted.

The incidence tangent equations and target differential are

    T_g Y = ker(p_s A_g),
    d nu_g = p_z A_g restricted to T_g Y.                (4)

Equations (1)-(4) work for every K and every simple G at once. Their
formal trace/symmetry part also holds for other nondegenerate trace
representations, so it cannot by itself establish an adjoint-specific
theorem. These are elementary deductions from the known trace construction
and equivariance, not a claim of originality.

## Where this proof attempt stops

The [existing generic-smoothness equivalence](CAYLEY_FIBER_REVIEW.md#2-a-uniform-criterion-and-which-part-is-already-standard)
reduces the conjecture to

    closure(nu(Sing Y)) != z.                           (5)

On the regular quotient, this is the same condition as
closure(mu(Sing X)) != A^s for the character incidence scheme. A sufficient
way to prove (5) is to show that on the generic point of every smooth
stratum of the reduced support of Sing Y, d nu has rank less than dim z.
Characteristic zero then bounds the image dimension. The original Y must
retain its nilpotents; taking Y_red would lose the obstruction.

The missing step is to deduce that rank bound from the actual adjoint
trace identities and root structure. Symmetry (2), the forced conjugacy
kernel (3), and a rank calculation at one good point do not prove it.
A good point proves only dominance or smoothness of one component, not
exclusion of a different bad component.

A small algebraic control explains the problem without computing another
group example. On A^2_(x,y) times Gm_z, take

    phi=x^3+(z-1)y.

Its relative critical incidence is

    X=V(3x^2,z-1),   mu=z phi_z=zy.

Hence every target c has quotient fiber C[x]/(x^2), with y=c and z=1.
The Hessian is symmetric and the gradient/target construction is exact,
but nonreduced structure persists over every target. This is not an
adjoint character or a counterexample to the group conjecture. It shows
why a general gradient argument needs a further representation-theoretic
input. Compact-real trace-function results likewise do not automatically
prove a complex scheme statement.

## A failed uniform shortcut, with an actual adjoint obstruction

The continuation on 24 October tested the stronger possible lemma

    Y is smooth at every g for which nu(g) belongs to Z_K^circ.

Equivalently, one might hope that singular incidence points necessarily
map to targets whose centralizer is larger than K. This lemma is FALSE.
The following exact counterexample tests that uniform proof route; it
does not add a new Levi to the coverage table or refute the generic
conjecture above.

Take G=Sp4(C), rho=Ad, and the Siegel Levi K=GL2. On its finite central
cover SL2 times Gm write x=tr(s), g=diag(zs,z^(-1)s^(-t)). Put

    H=diag(1,1,-1,-1), S=z^2+z^(-2), A=S+1,
    D=z partial_z, delta(z)=D A=2(z^2-z^(-2)).

The trace form is B(X,Y)=6 tr(XY), so B(H,H)=24. From
Ad(Sp4)=Sym^2(C^4), the restricted character and central target map are

    phi=A x^2-S=A(x^2-1)+1,
    X=V(A x),  mu=D phi=delta(z)(x^2-1).                  (6)

Choose z0 with z0^2 a primitive cube root of unity. Then

    A(z0)=0, delta0=delta(z0), delta0^2=-12,
    D delta(z0)=-4.

At x=0,z=z0 the incidence has completed local equation a*x=0, where
a=A(z) is a valid local coordinate. It is a reduced node, hence singular.
Nevertheless its target is

    c0=-delta0 != 0,
    h=(c0/24)H, Z_G(h)=K, B(h,h)=-1/2.                  (7)

Thus the bad image is not confined to larger-centralizer targets, even
for an actual adjoint character and a regular semisimple source point.

### Actual fiber scheme: a triple point, not merely a degenerate Hessian

Near this point write d=delta(z). The exact equations are

    a*x=0,
    d^2=4(a^2-2a-3),
    d(x^2-1)=-delta0.                                   (8)

Since d=delta0 and a=x=0 at the point, the relevant local branches are
unambiguous. Squaring the last equation, using the second, and reducing
modulo a*x gives

    a^2-2a+6x^2-3x^4=0.

Multiplication by x gives 3x^3(2-x^2)=0, hence x^3=0 locally. The
remaining equations give exactly a=3x^2 and d=delta0(1+x^2). Conversely
these substitutions satisfy all of (8). Therefore

    completed quotient fiber ring = C[[x]]/(x^3).       (9)

This retains the nonzero class x^2. It is not a reduced point or a
double point. The regular quotient adds dim K-rank G=2 smooth directions,
so at the corresponding regular source g0 the actual completed fiber is

    O_hat_(theta^(-1)(h),g0) = C[[u1,u2,x]]/(x^3).       (10)

For an explicit source take s=diag(i,-i). All eight root values of
g0=diag(iz0,-iz0,(iz0)^(-1),(-iz0)^(-1)) differ from one, so g0 is
regular semisimple in G. The defining-matrix adjoint Springer formula

    theta(g)=[tr(g)(g-g^(-1))+(g^2-g^(-2))]/12

gives (7) directly. The transverse logarithmic Hessian at the critical
point is diag(0,4). Independently, constructing Ad(g0) and ad(X) on an
explicit ten-element basis of sp4 gives rank d theta_(g0)=7; hence the
fiber tangent dimension is 3, consistent with (10) and its reduced
dimension 2.

The [exact audit](tools/audit_uniform_adjoint_obstruction.py) checks the
local ideal by saturation and elimination, reconstructs the character
from Sym^2, and compares the defining-matrix differential with the
ten-dimensional adjoint trace definition in all 100 matrix entries.
See the [audit report](uniform-adjoint-obstruction-audit.json). These are
independent algebraic constructions within this audit, not an independent
human proof review. The local argument (8)-(10) is the proof.

### Consequence for the general attempt

The two target values in this example satisfy c^2=-12, or equivalently
h=lambda H with lambda^2=-1/48. They are a finite proper subset of the
one-dimensional central target space. Thus this obstruction is entirely
compatible with the proved classical generic theorem and with our
all-simple-adjoint conjecture. No earlier generic theorem is retracted.

It rules out an actual proposed universal lemma, not just a generic
gradient argument: even adjoint root structure does not force every bad
incidence point to have a larger-centralizer target. A successful proof
of (5) must allow such nonzero exceptional targets and bound their image
dimension. The main conjecture and that dimension bound remain unproved.

A primary-literature check also located Garibaldi--Guralnick--Rains,
[Invariant derivations and trace bounds](https://arxiv.org/abs/2312.03101),
especially Lemma 3.3 and the invariant-derivation method. They provide
tools for critical points of invariant characters, including compact
trace bounds. They do not by themselves establish smoothness of our
parameterized complex critical schemes or eliminate their nilpotents.
Kostant--Michor 4.13 assumes an invertible differential at the semisimple
source; that hypothesis cannot be used at g0, whose differential has
rank seven. No claim of literature priority is made for this calculation.

## A uniform partial theorem from logarithmic differentials

The next continuation proves [central rigidity](CENTRAL_RIGIDITY_FIBERS.md),
recorded as GEO-73. For every nondegenerate trace representation, the
representation character is constant on each reduced irreducible component
of every semisimple-target Springer fiber. Its image in K/[K,K] lies in
a coset of the subtorus specified by the integer relations among the
central target coordinates. Outside finitely many rational homogeneous
target hyperplanes, that abelianization image is finite. For a
one-dimensional center it is finite at every nonzero target.

The proof uses exact-versus-logarithmic pole orders, then constructibility
and algebraic Baire. An independent isotropic-image argument proves the
generic finite-image assertion. Both apply to full fibers, including
nonregular source points, but only assert constancy on reductions.
The finite central image has an Artinian scheme structure whose
idempotents split the actual affine fiber; no nilpotent central directions
are discarded. This is a deduction from standard tools, not a novelty claim.

This result does not prove (5). The persistent double-structure control
above already has constant character and constant central position.
The unresolved task is still the adjoint-specific exclusion of dominant
singular incidence components, now with a proved restriction on their
generic central images. No generic smoothness or new Levi coverage is
claimed, and no case table is resumed.

## Continuation: generic targets do not justify scalar central substitution

The next attempt tested whether the finite central scheme image must
be reduced at generic targets. It does not follow from the trace identity,
even with phi exactly constant on each actual fiber. Section 5 of
[central rigidity](CENTRAL_RIGIDITY_FIBERS.md) constructs

    f=x^5+y^5+x^3 y^3,    phi=z(1+f),    mu=D_z phi.

On the origin-localized relative critical algebra A, dim A=16 and
0!=f satisfies f^2=0. For every c!=0 the actual fiber is A with
z=c(1-f), so its central image is a double point. Setting z=c loses
one dimension of the algebra. Saturation and a separately constructed
16-dimensional cyclic module verify this exactly.

This is an abstract polynomial potential, not an adjoint character.
It disproves another proposed inference from the general trace/gradient
identities, not the conjecture. In particular, even proving that the
character is scheme-theoretically constant would not by itself solve
the central-nilpotent issue. No adjoint-specific exclusion of this
mechanism has been established. A quartic three-variable version,
f_4=x^4+y^4+w^4+xyw, has the same obstruction (local lengths 11 and
10). Thus the existing quartic/four-variable bound alone also cannot
complete the proof; the actual root-derived coefficient constraints
must be used. Condition (5) is still the missing theorem.

## A uniform root-coefficient theorem and exact central-coordinate elimination

The next continuation establishes GEO-74 below. It uses the actual adjoint
roots, unlike the abstract quartic controls. It does not prove (5).

### Statement and normalization

Let G be simple simply connected of simply laced type A, D or E, and K
a proper Levi whose derived group S is a product of A-type groups. Put
Z=Z(K)^0 and work on the finite central multiplication cover S times Z.
Use the ACTUAL connected torus Z, not an arbitrary finite-index list of
central cocharacters. Normalize roots to squared length two. Write x for
the fundamental characters of S and e^gamma for a character of Z. The
trace form and Springer map retain the normalization above.

The existing norm theorem gives a restricted character of affine degree
at most three in x: intrinsic terms are x^2-1 on A1 and x_1 x_m-1 on
A_m for m>=2; external terms use one fundamental character on each of
at most three distinct factors. Let chi_3 be its homogeneous cubic part.

**Theorem (cubic coefficients, centralizer roots, and scheme elimination).**

1. Every occurring cubic monomial involving a factor larger than A1 has
   coefficient exactly ONE Laurent character e^gamma, with coefficient
   one. Its dual monomial has the inverse charge. In particular these
   coefficients never vanish on Z.
2. A cubic monomial M_i involving three A1 factors has coefficient exactly

       a_i=e^(gamma_i)+e^(-gamma_i),   ||gamma_i||^2=1/2.

   The vector delta_i=2 gamma_i is an ambient root orthogonal to S.
   Together with the three A1 roots and the roots of this constituent,
   it gives a D4 root subsystem. Moreover +/-delta_i is an isolated A1
   factor of the root system of the connected centralizer of S.
3. Distinct such gamma_i are orthogonal. Their triples of A1 factors
   intersect in zero or two factors, never one. In particular their
   number p is at most dim Z. The actual integral pairings are

       <gamma_i,delta_j^vee>=1 if i=j, and 0 otherwise.

   Consequently Z splits as (Gm)^p times Z', with independent primitive
   coordinates q_i=e^(gamma_i), and the corresponding cocharacters are
   delta_i^vee. No additional isogeny is required for this splitting of Z.
4. Set u_i=q_i+q_i^(-1). The ENTIRE restricted character can be written
   Phi(x,u,w), polynomial in x,u and Laurent in coordinates w of Z'.
   Its dependence on each u_i is exactly of the form

       Phi=u_i^2+u_i A_i(x,u_others,w)+B_i(x,u_others,w).       (11)

   Thus the square coefficient is the constant 1. Its cubic part in x is

       Phi_3=sum_i u_i M_i(x)+sum_m w^(eta_m) M_m(x),         (12)

   where the second sum is over mixed-factor cubics from part 1, including
   dual terms. These mixed coefficients are independent of every q_i.
   The divisors a_i=0 on Z form a simple normal-crossing union.
5. For h in Z_K^circ, write c_i=B_Ad(h,delta_i^vee); every c_i is nonzero.
   On the entire preimage of theta^(-1)(h) in S times Z, q_i^2-1 is a
   UNIT, including on its nonreduced scheme structure. Thus q -> u is an
   etale change of central coordinates along these fibers. On the regular
   quotient the original critical equations are scheme-theoretically
   equivalent to

       Phi_(x_j)=0,
       D_(w_l) Phi=c'_l,
       (u_i^2-4) (Phi_(u_i))^2=c_i^2,                       (13)

   with inverse reconstruction

       q_i=(u_i+c_i/Phi_(u_i))/2,
       q_i^(-1)=(u_i-c_i/Phi_(u_i))/2.                      (14)

   The equations themselves make Phi_(u_i) and u_i^2-4 units, so the
   displayed division is legitimate. These are equalities of actual
   schemes, not of reductions. Retain the induced finite central-kernel
   action when descending from the cover. No smoothness of (13) is asserted.

The regular-source hypothesis is needed for the critical-quotient
presentation in part 5, not for the root statements or the unit condition.
G2 and F4 are outside this simply laced theorem; no new Levi smoothness
coverage is claimed in any type.

### Proof of coefficient sparsity

For an external constituent, let mu be its dominant S-highest weight and
gamma its orthogonal central projection. Its highest weight

    alpha=mu+gamma

is an ambient root. Its multiplicity is one because the ambient root
space has dimension one. For a cubic, mu is a sum of three fundamental
weights on distinct A-factors, so m=||mu||^2>=3/2, with equality precisely
for three A1 factors. All charges for this same mu have norm squared 2-m.

If two distinct charges gamma,eta occurred, the distinct roots
alpha=mu+gamma and beta=mu+eta would satisfy

    (alpha,beta)=m+(gamma,eta)>=m-(2-m)=2m-2.

Distinct nonopposite simply laced roots have inner product at most one.
For m>3/2 this is impossible. For m=3/2 equality forces eta=-gamma;
there can be no third charge. The A1 triple is self-dual, so inversion
of all roots guarantees that BOTH charges occur. This proves parts 1-2's
coefficient claims, without assuming independent arbitrary coefficients.

### The D4 subsystem and centralizer factors

For an A1 triple, alpha=mu+gamma and beta=mu-gamma have inner product one.
The root-string rule implies alpha-beta=2gamma=delta is a root.
It is orthogonal to every Levi root. Let beta_1,beta_2,beta_3 be the
three A1 simple roots, so mu=(beta_1+beta_2+beta_3)/2. Reflections in the
four mutually orthogonal roots beta_1,beta_2,beta_3,delta give the roots

    +/-beta_1, +/-beta_2, +/-beta_3, +/-delta,
    (epsilon_1 beta_1+epsilon_2 beta_2+epsilon_3 beta_3
                         +epsilon_4 delta)/2,  epsilon_i in {+1,-1}.

These 24 roots form a D4 subsystem. This is a consequence of standard
root strings, not a new classification of D4 or a new group example.

If beta is any root orthogonal to S, root integrality gives

    (alpha,beta)=(delta,beta)/2 in Z.

For beta different from +/-delta, (delta,beta) can only be -1,0,1.
It must therefore be zero. Hence delta is an isolated A1 factor in the
centralizer root system. This argument does not assume that the
centralizer of S is itself a Levi.

For two distinct A1 triples, let k be the number of shared factors.
Their S-highest weights have pairing k/2. Their central charges have
norm squared 1/2. Integrality of the ambient root pairing forces their
central pairing to be zero for k=0 or k=2. For k=1 it would have to be
plus or minus 1/2, so the two charges agree up to sign. Choose signs to
make them equal. The two highest roots then have inner product one,
and their difference is a root with zero central projection. Every such
root belongs to the Levi. But this difference has nonzero components on
four distinct A1 factors, and cannot be a root of that Levi. Contradiction.

Thus the charges are orthogonal. Since delta_i=2gamma_i has norm squared
two, its coroot gives precisely the stated Kronecker pairings. Both
gamma_i in X^*(Z) and delta_i^vee in X_*(Z) are integral: the first is the
restriction of the constituent's ambient root, and the second is an
ambient coroot orthogonal to the Levi. The pairings split the character
and cocharacter lattices. Explicitly,

    z' = z product_i delta_i^vee(q_i(z))^(-1)

lies in the common kernel Z' of the q_i; this gives the torus product
and proves connectedness of Z'.

### The entire character and the actual fiber equations

The root SL2 for delta_i centralizes S. Its Weyl reflection fixes S and
Z' and sends q_i to q_i^(-1), fixing the other q_j. Hence chi is separately
invariant under these inversions, and is a polynomial in the u_i over
the remaining character ring. For any ambient root alpha, its q_i
exponent is <alpha,delta_i^vee> in {-2,-1,0,1,2}. The exponents +/-2
occur only for alpha=+/-delta_i themselves. Those roots contribute
q_i^2+q_i^(-2)=u_i^2-2 with coefficient one, and no other variables.
All remaining contributions are at most linear in u_i. This proves (11).
The uniqueness of each mixed cubic coefficient and inversion invariance
force its q_i exponent to vanish. This proves (12).

Each a_i=0 is q_i^2=-1. Its logarithmic derivative is q_i-q_i^(-1), which
is nonzero there. The q_i are independent coordinates, proving the
normal-crossing assertion. This concerns the cubic-coefficient boundary,
not the singular locus of the complete critical scheme.

For the fiber assertion, delta_i is outside the Levi, so delta_i(h)!=0
when Z_G(h)=K. The invariant form makes c_i a nonzero scalar multiple
of delta_i(h). If q_i^2=1 at a point (s,z) above the fiber, the root
subgroup U_(delta_i) would commute with both s and z, hence with g=sz.
Equivariance would then make it centralize theta(g)=h, contradicting
delta_i(h)!=0. Thus q_i^2-1 has no zero on the preimage fiber and is a
unit in its coordinate ring. This argument is valid for nonreduced
fibers as well; it does not replace them by their reductions.

The chain rule gives the exact central target equation

    (q_i-q_i^(-1)) Phi_(u_i)=c_i.

Squaring and using (q_i-q_i^(-1))^2=u_i^2-4 yields (13). Conversely,
because c_i!=0, (13) makes both factors units. Formula (14) gives
q_i q_i^(-1)=1, q_i+q_i^(-1)=u_i, and recovers the original UNSQUARED
target equation. Thus no sign choices or extraneous sheets have been
introduced. Affine and w-derivatives are unchanged by this independent
central-coordinate substitution. All steps are equivariant with the
original finite central-kernel action, which must be retained for descent.

### What is proved, checked, and still missing

This proves a uniform constraint and an exact change of equations for
the actual adjoint problem. It excludes arbitrary cubic coefficient
functions, multiple vanishing of an individual cubic coefficient,
nontransverse intersections between these coefficient divisors, and
ramification of q -> u along the specified Springer fibers. It does NOT
exclude nilpotents or singularities produced by the coupled affine
critical equations (13). The main conjecture (5) remains open.

The inputs are the existing adjoint norm/branching theorem, standard
simply laced root strings and root integrality, the trace identity and
equivariance. Kostant's [Levi central-weight-space theorem](https://arxiv.org/pdf/0711.2809),
Theorem 0.1, gives the published irreducible/multiplicity-free branching
framework. The coefficient, centralizer and elimination assertions here
are written deductions; no literature-priority claim is made.

The [direct-root audit](tools/audit_cubic_root_coefficients.py) reconstructs
cubic constituents by a root-raising test, independently of the stored
fundamental-character polynomial inversion. Across the existing 136
cubic E6/E7/E8 records it checks 230 mixed cubic monomials, 78 A1-triple
binomial pairs and their D4 subsystems, and 48 pairs of orthogonal charges.
It checks the integral coroot pairings and the exact coordinate identities.
See the [audit report](cubic-root-coefficients-audit.json). These are finite
corroborations of the uniform proof, not new coverage or independent human
peer review. No additional Levi enumeration or fiber table is performed.

The next proof step is to use the constrained form (11)-(13) to bound the
image of its singular incidence locus. In particular, any reduction must
retain the equations (u_i^2-4)(Phi_(u_i))^2=c_i^2, rather than treating the
u_i as freely fixed scalar coefficients. A proof must control all coupled
components and their nilpotents, not merely the coefficient-zero boundary.

## All nonlinear coefficients: centralizer characters and their boundaries

The next continuation proves GEO-75. It extends the root argument to
EVERY external quadratic and cubic coefficient, not just the cubic part.
It still does not prove (5). The hypotheses and trace normalization are
those of GEO-74: simply laced simple simply connected G, rho=Ad, a proper
all-A Levi K, and the actual connected central torus Z=Z(K)^0.

### Uniform coefficient theorem

Let M(x) be an occurring external monomial of affine degree at least two.
Its factors belong to distinct simple factors of S=[K,K]. Let mu be its
S-highest weight, m=||mu||^2, and let Gamma be the central projections of
the ambient roots with that highest weight. Its coefficient is

    a_M(z)=sum_(gamma in Gamma) e^gamma(z).                 (15)

Each summand occurs once, and 1<=m<2. These are the external coefficients;
the intrinsic terms x^2-1 and x_1 x_m-1 are not included in (15).

**Theorem.** Exactly the following two descriptions apply.

1. If m>1, Gamma is the vertex set of a regular simplex. For p=|Gamma|,

       (gamma_i,gamma_j)=1-m for i!=j,
       p(m-1)<=1,
       ||(sum gamma_i)/p||^2=1-m+1/p.                     (16)

   If p>=2, the roots gamma_i-gamma_j form an isolated A_(p-1)
   component of the root system of H=Z_G(S)^0. Every other root of H
   is orthogonal to every gamma_i. The coefficient is the character of
   its defining p-dimensional module, with a common central character
   twist. More precisely, after an appropriate finite central cover,

       a_M=d(v_1+...+v_p),  product v_i=1,

   where d is a Laurent unit and the v_i are the standard weights of
   that A_(p-1) factor. If p=1, the coefficient is a Laurent unit.
   The zero divisor of a_M, whenever nonempty, is smooth on Z itself.

2. If m=1, M is a product of two A1 characters. For some p>=1,

       Gamma={+/-gamma_1,...,+/-gamma_p},
       (gamma_i,gamma_j)=1 if i=j, and 0 otherwise,
       a_M=sum_i(q_i+q_i^(-1)),  q_i=e^(gamma_i).          (17)

   For p>=2 the roots +/-gamma_i+/-gamma_j, i!=j, form an isolated
   D_p root subsystem of H; D_2 means A1+A1 and D_3 means A3.
   Every other root of H is orthogonal to every gamma_i. Thus (17)
   is the defining vector character of that factor on its torus.
   For p=1 it is a pair of inverse torus characters; no D_1 root
   system or additional SL2 factor is asserted.

The words "isolated" mean that all other centralizer roots are
orthogonal to the displayed subsystem. D_2 is the union of its two
isolated A1 components. The theorem concerns root systems and characters;
it does not identify H with a direct product of simply connected groups.
The q_i are independent characters, but their integral span need not be
primitive in X^*(Z). Only GEO-74's special cubic situation supplied the
stronger integral splitting. Any covers here are finite etale torus or
central covers, and their finite kernel actions must be retained.

In particular, every nonlinear external coefficient is a standard
type-A character with a central twist, a standard type-D vector
character, a single unit, or an inverse-character pair. This is a
uniform description in arbitrary simply laced rank. It is not a
classification of the linear or constant terms, or of simultaneous
couplings between different monomials.

### Proof by root pairings, without a branching table

For gamma in Gamma the corresponding highest weight alpha=mu+gamma
is an ambient root. Its root space is one dimensional, which proves
the multiplicity-one assertion. All charges have squared norm 2-m.
The existing norm theorem gives at most three fundamental weights on
distinct A-factors. Each has squared norm at least 1/2; hence m>=1.
The charges are nonzero for these external monomials, so m<2.

For distinct charges gamma,eta, the distinct nonopposite roots
alpha=mu+gamma and alpha'=mu+eta satisfy

    (alpha,alpha')=m+(gamma,eta)>=2m-2.                    (18)

If m>1, the integer on the left is positive and at most one, so it
equals one. Thus (gamma,eta)=1-m and alpha-alpha'=gamma-eta is a
root orthogonal to S. The differences form the A_(p-1) root system:
relative to their barycenter the p weights have Gram matrix
I-(1/p)11^t. The barycenter is orthogonal to all the differences,
and its norm gives (16), including the bound on p.

The roots of H are exactly the ambient roots orthogonal to S. Reflection
in such a root beta fixes mu and permutes Gamma: it fixes every S-simple
root, so preserves the highest-root raising test. Moreover

    <gamma,beta^vee>=<alpha,beta^vee> in {-1,0,1}.

The values +/-2 would require alpha=+/-beta, contrary to mu!=0.
If this reflection moves gamma to gamma', then beta is, up to sign,
gamma-gamma'. Consequently every H-root acting nontrivially on Gamma
is already in the displayed A root system. This proves isolation and
the defining-character description. After a finite central cover the
common barycenter gives the stated central twist. No integrality of
the barycenter in X^*(Z) was assumed.

If m=1, the norm lower bound forces exactly two A1 fundamental weights.
Their sum mu is self-dual. Taking negatives of ambient roots and then
the longest Weyl element of S therefore proves Gamma=-Gamma. For two
charges which are neither equal nor opposite, their scalar product is
strictly between -1 and 1. Integrality of 1+(gamma,eta) forces it to
be one, so the charges are orthogonal. Selecting one charge from each
opposite pair gives (17). For distinct axes all +/-gamma_i+/-gamma_j
are root differences of this kind. They form D_p. The same reflection
argument proves that no other H-root acts nontrivially on these axes;
a reflection cannot send gamma to -gamma because 2gamma has squared
norm four and is not an ambient root. This proves part 2.

For part 1's boundary assertion, differentiation along the actual
centralizer coroot gamma_i-gamma_j gives

    D_(gamma_i-gamma_j)^vee a_M=e^(gamma_i)-e^(gamma_j).

If all these derivatives vanished, all p summands would have the same
nonzero value, and their sum could not vanish in characteristic zero.
Thus a_M=0 is smooth, directly on the actual torus Z.

### What survives on actual Springer fibers

There is a useful uniform localization, valid also outside the simply
laced setting. Let R_H be the centralizer roots. For every h in Z_K^circ,
on the entire preimage of theta^(-1)(h) in S times Z,

    e^beta(z)-1 is a UNIT for every beta in R_H.           (19)

Indeed, if e^beta(z)=1, the root subgroup U_beta commutes with s and z.
It would commute with g=sz and, by equivariance, fix theta(g)=h. But
beta is outside the Levi root system, so beta(h)!=0. This contradiction
proves the absence of complex zeros. The algebraic Nullstellensatz then
gives a unit in the actual fiber ring, including any nilpotents. The
argument does not require a regular source point.

Write Z_H^reg for the open defined by all the units in (19). All the
central coordinates of these fibers lie there. On this open each
individual coefficient divisor in (15) is smooth, with just one possible
kind of exception: a transverse crossing of two smooth branches arising
from the D_2 vector character in (17).

To prove this, the independent characters q_i in (17) give smooth torus
coordinates after an etale cover. Their logarithmic derivatives are
q_i-q_i^(-1). A singular point of a_M=0 must therefore have

    q_i=epsilon_i in {+1,-1},  sum_i epsilon_i=0.          (20)

This requires p even. The logarithmic Hessian is diag(2 epsilon_i),
so transversely these are nondegenerate quadratic hypersurface
singularities. If p>=4, two signs agree, and the centralizer root
gamma_i-gamma_j has value q_i/q_j=1, forbidden by (19). If p=2,
the balanced signs give the remaining possibility. The exact factorization

    q_1+q_1^(-1)+q_2+q_2^(-1)
       =(q_1+q_2)(q_1 q_2+1)/(q_1 q_2)                  (21)

has transverse factors at (1,-1) and (-1,1). The two D_2 root values
q_1/q_2 and q_1 q_2 are both -1 there, so (19) does not exclude this
crossing. Extra central directions are smooth parameters throughout.

This proves reducedness and the stated local shape of each INDIVIDUAL
nonlinear-coefficient divisor on the relevant central open. It does not
prove that distinct coefficient divisors meet transversely, that their
union is a normal-crossing divisor, or that the affine critical incidence
is reduced. In particular (19)-(21) do not justify specializing the
central coordinates to scalars in a nonreduced fiber.

### Simultaneous coefficient vanishing is reduced

The isolated-factor conclusion gives more than individual boundary
control: the SAME centralizer character coordinates can be used for
all the coefficients at once.

**Corollary (simultaneous coefficient ideals).** For any collection C
of the occurring external quadratic and cubic coefficients, the scheme

    V(a_M : a_M in C) inside Z_H^reg                         (22)

is reduced, possibly empty. Etale locally it is a union of coordinate
subspaces over a smooth torus coset. This includes all simultaneous
vanishings; no genericity assumption on the chosen central point is made.
It is not an assertion about the critical ideal (Phi_(x_j)).

Proof. Take the standard finite central cover H_sc times E -> H, where
H_sc is the simply connected cover of [H,H] and E=Z(H)^0. Its maximal
torus gives a finite etale cover of Z. Fix one simple system of H and
let xi_j be the corresponding fundamental characters of H_sc.

Every nontrivial multiplicity module above is minuscule on the factor
which acts: all its weights lie in the displayed single Weyl orbit and
have multiplicity one, and their pairings with every coroot lie in
{-1,0,1}. On an irreducible root component, its highest weight is a
single fundamental weight. Indeed two nonzero fundamental coefficients,
or a coefficient exceeding one, would pair at least two with the
highest coroot, whose simple-coroot coefficients are positive.
The only reducible acting root system in the classification is D_2;
its highest weight uses one fundamental weight from each of its two
A1 components. This argument also handles D_3=A3 and D4 triality in
the SAME fixed simple system; a separate choice of vector node for
each coefficient is not being made.

Consequently, simultaneously on this cover, each coefficient is a unit
times one of

    1,  xi_j,  xi_j xi_k (j!=k),  or psi^2+1,              (23)

where each psi is a character of E alone. The last possibility comes
from the inverse-character pairs which are trivial on [H,H]; multiply
psi+psi^(-1) by the unit psi. The unit twists do not affect their ideals.
The products xi_j xi_k come only from the D_2 case, and have no repeated
factor. Distinct coefficients may use the same xi_j or the same product.

The map from the covered torus to A^(rank H_sc) times E, given by the
fundamental characters and E-coordinates, is etale on the preimage of
Z_H^reg. This is the standard simply connected torus quotient: the
fundamental-character Jacobian is a nonzero unit times the Weyl
denominator, so its zero set is exactly the root loci e^beta=1.
The fundamental-character presentation is a published input; see
[Steinberg, section 6](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf).

For any selection from (23), the xi-generators define a squarefree
monomial ideal, hence a reduced union of coordinate subspaces. The
selected binomial equations psi_i^2=-1 define a fiber of a homomorphism
from E to a torus. Such a fiber is empty or a translate of its kernel;
over C the kernel is smooth, though it can be disconnected. Thus the
binomial ideal defines a reduced disjoint union of smooth torus cosets.
These two sets of variables are independent. Their product is reduced,
with the stated coordinate-subspace description. A selected unit gives
the empty scheme. Etale pullback preserves reducedness, and faithful
flatness of the central cover descends it to Z_H^reg. This proves (22).

The coordinate-subspace union can be singular. For example (xi_1 xi_2,
xi_2 xi_3) is radical but has components of different dimensions. Nor
does (22) imply that all coefficient divisors form a normal-crossing
union: distinct central binomial constraints can be dependent. These
limitations are compatible with the assertion of reducedness.

### Audit and the remaining proof obligation

The extended [direct-root audit](tools/audit_cubic_root_coefficients.py)
reconstructs all external degree-two and degree-three coefficients in
the existing 136 cubic-character records: 728 quadratic and 308 cubic
coefficient weight sets. It compares them with the
stored character expansions, checks the complete centralizer root
systems and the simplex/cross-polytope identities, and separately checks
the boundary derivative, Hessian and factorization formulas. It fixes a
single centralizer simple system per record and checks all highest-weight
labels in it: 376 unit coefficients, 514 single fundamental-character
coefficients, 88 products of two distinct A1 fundamental characters, and
58 central binomials. Coordinate-prime decompositions independently
check the resulting squarefree monomial ideals. The finite
checks corroborate the root proof; they are not a new coverage census.
The [report](cubic-root-coefficients-audit.json) continues to mark the
generic smoothness conjecture as unproved.

Published inputs are the root-string/reflection rules, standard
defining characters, Kostant's Levi branching framework, Steinberg's
fundamental-character quotient, and elementary squarefree monomial and
characteristic-zero torus-fiber facts.
The classification and fiber-localization consequences are written
deductions here, with no priority claim. A targeted search for a direct
adjoint-character/Levi critical-fiber theorem did not identify one; this
was not a comprehensive literature audit.

The precise next obligation is still (5). This step rules out arbitrary
nonlinear coefficient functions and nonreduced simultaneous coefficient
vanishing on the relevant central open. A successful proof must now
control the AFFINE CRITICAL EQUATIONS over these reduced coefficient
strata, including the linear terms and the full central target equations.
Derivatives couple the x variables and different coefficients: a radical
coefficient ideal does not imply a radical critical ideal. Neither (22)
nor the localization establishes that the singular incidence has
nondominant target image.

## Rank-bound attempt: the nonlinear coefficients alone do not suffice

The continuation tested whether the common coordinates in (23), together
with the intrinsic quadratic terms, force (5). They do not. The following
control keeps EVERY nonlinear coefficient and the constant term of an
actual adjoint restriction unchanged. Only its linear term is changed.
It identifies a missing adjoint identity, not a counterexample to the
adjoint conjecture. This reuses an existing formula; no new Levi case or
coverage claim is involved.

### An exact control with unchanged nonlinear coefficients

Use the already proved E6 formula for K with derived group A2+A1+A2.
On its degree-six central cover, write (a,b), x, (u,v) for the fundamental
characters of SL3, SL2, SL3, respectively, and z for the central coordinate.
Set

    S=z^3+z^(-3), A=z^(-2)+zx, B=z^2+z^(-1)x,
    Phi_k=ab+uv+x^2-2+k Sx+A av+B bu.                       (24)

Phi_1 is exactly chi_Ad|K. All external nonlinear coefficients are the
Laurent units z, z^(-1), z^2, z^(-2). The intrinsic quadratic terms and
constant term are independent of k. Thus changing k from 1 to 2 does
not change any nonlinear-coefficient vanishing scheme, even before
localization. Here the derived-Levi centralizer has no roots.

Phi_2 is itself an honest representation character of the reductive K.
The module W with character Sx is the sum of the two SL2 defining
constituents of central charges +3 and -3 in Ad(G)|K. Hence it descends
to K, and

    R=Ad(G)|K direct_sum W,   dim R=82,   chi_R=Phi_2.

Its trace form is nondegenerate. This follows, for example, by restricting
to a compact form: a faithful differential has negative-definite trace
form there. R contains Ad(G)|K with faithful differential. The central
trace norm changes from 144 to 180, since the four added weights each
have central charge +/-3. Accordingly the central target coordinate is
c=B_R(h,H)=180 t for h=tH in the original primitive cocharacter direction.
R is NOT an adjoint representation of a simple ambient group, and no
extension of R from K to G is asserted.

The exact completion a'=a+B u, b'=b+A v uses AB=1+Sx+x^2 and gives

    Phi_2=a'b'+x^2+2Sx-x(x+S)uv-2.

Eliminate a',b' by their derivative equations and put X=x+S. The
remaining potential is

    F=-S^2-2+X^2(1-uv)+S Xuv.                              (25)

Near X=u=v=0 and S!=0, the relative critical ideal is generated by

    2X(1-uv)+Suv,  X(X-S)u,  X(X-S)v.

Since X-S is a unit, these equations are EXACTLY equivalent to

    2X+Suv=0,  u^2v=0,  uv^2=0.                           (26)

For the forward inclusion, Xu=Xv=0; multiplying the first derivative
by u and v then gives u^2v=uv^2=0, and Xuv=0 removes its remaining term.
Conversely (26) implies every original derivative and makes X-S a unit.
Thus the local relative critical algebra is C[u,v]/(u^2v,uv^2), with
the indicated base parameter z. The class uv is nonzero and (uv)^2=0.
This is an embedded nilpotent, not just a reduced crossing.

The full central equation must be retained. Differentiating (25) at
fixed X,u,v and reducing modulo (26) gives exactly

    mu=D_z Phi_2=-2S D_z S=-6(z^6-z^(-6)).                 (27)

The z-dependent coordinate changes contribute only multiples of the
affine critical equations, so (27) also holds for the ORIGINAL central
derivative in the nonreduced quotient. It is a nonconstant dominant
map to the one-dimensional target. Its logarithmic derivative is
-36(z^6+z^(-6)); its critical target values satisfy c^2+144=0.

For c!=0 and c^2+144!=0, every solution of (27) has S!=0 and is an
etale central point. At u=v=0 the two SL3 character pairs are (0,0),
represented by regular semisimple matrices with eigenvalues 1,omega,
omega^2. The SL2 trace is -S, with S^2!=4 because c!=0. Hence this is
a K-regular semisimple source point. The completed regular group fiber
there is

    C[[r_1,...,r_14,u,v]]/(u^2v,uv^2).                    (28)

The 14 smooth directions are dim K-rank K=20-6. The finite central
cover is etale and does not remove the nilpotent. Therefore generic
central fibers for this honest reductive representation are nonreduced.
Its singular incidence DOES dominate the target, despite identical
nonlinear coefficients and constant term to the adjoint model.

### The precise identity missing from the proposed rank argument

In (24) the unchanged bilinear pencil is

    M(x,z)=[[1,A],[B,1]],  det M=-x(x+S).

Its two eigenvalue branches are lambda=0 and lambda=-S. The affine
quadratic part g_k=x^2+kSx-2 has derivative 2x+kS. For k=1 its values
on those branches are S and -S, so neither vanishes identically. For
k=2 the second value is identically zero. The quadratic critical point
has aligned with a singular pencil branch, producing (26). The genuine
adjoint normalization prevents THIS resonance; the nonlinear-coefficient
classification alone cannot distinguish the two models.

There is a useful general arithmetic identity explaining this particular
nonresonance. Let R0=Z[Lambda] be an integral Laurent character ring,
f in R0, and A0,B0 square m by m matrices over R0. If f and det B0 are
both nonzero modulo 2, then

    det(A0-(f/2)B0)!=0 in Frac(R0).                        (29)

Indeed det(2A0-fB0) reduces modulo 2 to (-f)^m det B0, which is nonzero
in the domain F_2[Lambda]. Thus no root lambda of det(A0+xB0) can obey
2lambda+f=0. In the simply laced all-A setup used above, an external linear
A1 coefficient in the raw adjoint character, when nonzero, is a sum of
DISTINCT central Laurent characters, each
with coefficient one, by ambient root multiplicity one. It therefore
satisfies the oddness condition on f. In (24), det B0=-1 and f=S;
doubling f destroys precisely this condition.

Identity (29) is a restricted nonresonance lemma, NOT the rank bound
(5). The full character need not have a bilinear-pencil presentation,
and nonlinear elimination need not preserve these integral hypotheses.
It does not control intersections of critical branches or all central
boundary strata. No additional adjoint smoothness cases follow here.

The [extended audit](cubic-root-coefficients-audit.json) matches (24)
with the archived adjoint character, independently reconstructs the
four added root weights and the trace norms, checks both localized-ideal
inclusions in (26), and replays ALL original affine and central
derivatives in the nonreduced algebra. This is a proof-method control,
not another exceptional-group classification or a disproof of the
uniform adjoint conjecture. It refines the older
[persistent-collision control](SPECTRAL_BRANCHING_FIBERS.md#3-a-sharp-obstruction-the-graph-alone-is-insufficient)
by preserving the actual adjoint nonlinear coefficients and constant term.

The rank attempt therefore stops at a precise requirement: use the
FULL adjoint weight relations, including the multiplicity-one linear
coefficients, to exclude dominant critical-branch resonance. Reduced
coefficient strata cannot supply this exclusion on their own. The
next attempt should address that nonresonance problem directly, rather
than extend the coefficient-boundary classification or a Levi table.

## 25 October: invariant differential identities do not close the gap

This continuation tested a uniform differential-operator argument using
the full root character. It did NOT prove (5), produce an adjoint
counterexample, or increase the smoothness coverage. Two distinctions
prevent the proposed argument from being a proof.

### The root-length equation does not enforce the required multiplicities

For a simply laced root system R of rank r, normalize root lengths to 2.
The flat logarithmic torus Laplacian Delta is characterized by

    Delta(e^lambda)=(lambda,lambda)e^lambda.

Consequently the actual adjoint character satisfies

    Delta(Phi)=2(Phi-r),  Phi=r+sum_(alpha in R) e^alpha.   (30)

This elementary identity involves the flat torus operator, NOT the
radial group Casimir. It also holds for r+sum m_alpha e^alpha with
arbitrary constant root multiplicities m_alpha. In particular the
82-dimensional reductive-Levi character (24), k=2, still satisfies
(30) with r=6: its four added weights are ambient E6 roots and its
zero-weight multiplicity remains 6. Its nonreduced generic fibers
therefore rule out deriving (5) from (30), the trace identities, and
the previously established nonlinear-coefficient data together.

The altered character is invariant under W_K, but not under W_G. This
is exactly where the full adjoint symmetry has been lost. Conversely,
for an irreducible simply laced root system, W_G is transitive on R.
Thus a W_G-invariant Laurent polynomial supported on {0} union R has
the form c0+c1 sum e^alpha. Fixing the adjoint root multiplicity leaves
no deformation of its linear coefficients within this support. A generic
coefficient theorem cannot simply be specialized to the adjoint point;
preserving full ambient symmetry only supplies the scalar family, whose
nonzero members have the same smoothness question after target scaling.
This is an elementary restriction on a possible deformation proof, not
a proof that the conjecture is false or that all deformation methods fail.

### Invariant derivations do not replace the actual quotient Jacobian ideal

The current source is Garibaldi--Guralnick--Rains,
[Invariant derivations and trace bounds, v4, sections 4-5](https://arxiv.org/pdf/2312.03101v4).
Their invariant-derivation matrix has the form M=J^t A J, where J is
the torus-to-character-quotient Jacobian and A is the nondegenerate
torus form (equation 5.6 and Proposition 5.4). Its use for torus critical
points does not identify those equations with the ordinary partial
derivatives on the smooth character quotient at the discriminant.

For example, on SL2 let x=t+t^(-1). The adjoint character is x^2-1.
The invariant torus derivation becomes (x^2-4) partial_x, so its
critical equation on the quotient is

    2x(x^2-4)=0,

whereas the ordinary quotient critical equation is 2x=0. At x=2
these disagree. The semisimple representative is the identity, which
is critical on SL2 but is not a regular source. A regular unipotent
representative u=[[1,1],[0,1]] has the SAME quotient coordinate x=2,
but, with F=[[0,0],[1,0]],

    d chi_Ad at u applied to uF = 2 tr(u) tr(uF) = 4.

Hence it is not critical. This is a direct check of why the regular
group quotient, not the torus critical locus alone, is essential in
our conjecture. There is no contradiction with the cited paper's
semisimple-point method. On the torus-regular open, M is invertible
and the two systems are equivalent; that equivalence does not cover
all regular nonsemisimple group sources in the conjecture.

The existing audit now checks (30) root by root for the true and altered
characters, supplies an explicit ambient-reflection failure for the
altered multiplicities, and checks the SL2 calculation both in torus
coordinates and by matrix differentiation. These are checks of this
failed proof route, not new smoothness cases. The remaining gap is still
to use FULL ambient Weyl/root relations to bound the central image of
the singular locus of the ACTUAL relative critical scheme. No such
bound follows from the differential identities tested here.

## 25 October: what bracket preservation actually says about the differential

The next attempt uses U=Ad(g) as an automorphism of the WHOLE simple
Lie algebra, rather than only its torus eigenvalues. It gives the exact
operator formulas below, valid in every simple type and also at
nonsemisimple g. They are elementary Lie-algebra deductions from the
trace construction, not a new generic-smoothness theorem or priority claim.

Extend A_g from k to the full left-translated differential mathcal A_g
on g. Let e_i and e^i be Killing-dual bases. Define the conjugation
Casimir on End(g) by

    C_End(V)=sum_i [ad(e_i),[ad(e^i),V]].

The normalization B(X,Y)=tr(ad(X)ad(Y)) gives

    sum_i ad(e_i)ad(e^i)=Id_g.

Indeed this Casimir commutes with the irreducible adjoint action; its
trace is dim g, so its scalar is 1. The full differential is

    mathcal A_g=sum_i ad(e_i) U ad(e^i)
               =U-(1/2) C_End(U).                         (31)

To check the first equality, invariance of B and the dual-basis trace
formula give, for ANY endomorphism V,

    B(sum_i [e_i,V[e^i,X]],Y)=tr(V ad(X) ad(Y)).

Set V=U and use the trace differential (1). Expanding the double
commutator and using the Casimir identity proves the second equality.
In particular, (31) is a useful exact expression for the differential,
but the contraction identity itself does not require V to preserve
brackets. It is not an extra equation on the actual adjoint differential.

Writing * for the B-adjoint and h=theta(g), one also obtains

    mathcal A_g-mathcal A_g^*=ad(h),
    mathcal A_g U=U mathcal A_g,
    mathcal A_g (Id-U^(-1))=ad(h).                         (32)

The first identity is the trace commutator calculation (2) on all of g.
For the last one, differentiate conjugation by exp(tX). Its left tangent
at g is U^(-1)X-X, and equivariance gives the output [X,h]. Thus
mathcal A_g(U^(-1)-Id)=-ad(h). The middle identity follows by conjugating
(31) by U and changing the dual bases. Bracket preservation enters through
U ad(X) U^(-1)=ad(UX); the resulting last identity is precisely the full
equivariance identity, not an additional rank condition.

For a semisimple target with centralizer K, decompose
g=k direct-sum m, where m is the sum of the nonzero ad(h)-eigenspaces.
Both U and mathcal A_g preserve this decomposition. On m, (32) shows
that BOTH mathcal A_g and Id-U^(-1) are invertible. On k its right side
is zero, and it reduces to the already known conjugacy kernel (3).
Consequently this argument controls the complementary directions but
does not bound the remaining rank-r transverse differential inside k.

### Linearizing the automorphism equations adds no missing tangent condition

Write a first-order variation as U_epsilon=U(Id+epsilon D), epsilon^2=0.
Expanding U_epsilon[X,Y]=[U_epsilon X,U_epsilon Y] says exactly

    D[X,Y]=[DX,Y]+[X,DY].

Every such derivation is inner. Here is a short proof with the same
normalization: choose Z by B(Z,X)=tr(D ad(X)). The derivation property
gives [D,ad(X)]=ad(DX), so

    B(DX,Y)=tr([D,ad(X)]ad(Y))
           =tr(D ad([X,Y]))=B([Z,X],Y).

Nondegeneracy gives D=ad(Z). Thus these are exactly the variations
already parametrized by g exp(epsilon Z). Bracket preservation is a
substantive global constraint on U, but its linearization cannot be
counted again as an extra condition on tangent vectors already in G.
This does not rule out a future nonlinear use of the automorphism
equations; it identifies why the proposed first-order shortcut stops.

### The differential kernel need not be a Lie subalgebra

There is also a concrete obstruction to trying to close the argument
by bracket-closing ker mathcal A_g. In the SAME Sp4/GL2 family used
above, take

    g=diag(2i,-2i,-i/2,i/2),
    h=theta(g)=(-5/16) diag(1,1,-1,-1).

The source is regular semisimple and Z_G(h)=GL2. With matrix units,
put E=E_12-E_43, F=E_21-E_34 and
H_alpha=[E,F]=diag(1,-1,-1,1). The exact differential satisfies

    mathcal A_g E=mathcal A_g F=0,
    mathcal A_g H_alpha=-(7/4) H_alpha != 0.               (33)

Hence the kernel is not closed under the Lie bracket. This concerns
left-translated tangent vectors at ONE point; it does not contradict
integrability of tangent vector fields on a smooth fiber, whose left
coordinates vary with the point. It is not a counterexample to the
generic-fiber conjecture and supplies no new Levi coverage.

The [extended existing audit](tools/audit_uniform_adjoint_obstruction.py)
checks (31)-(32) against both the defining-matrix and adjoint-trace
differentials, and checks (33) by two matrix constructions. The missing
result remains (5). No target polynomial has been obtained from this
bracket argument, and the uniform conjecture remains unproved.

## 25 October: higher-order obstruction and the corrected quartic

This continues the uniform attempt, without adding a Levi case. It gives
an exact local test, but DOES NOT prove that its failure locus has
nondominant central image. The formal splitting step is standard:
[Greuel--Lossen--Shustin, Theorem I.2.47](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/greuel.pdf).
The formulas below are direct expansions of that step, not a novelty claim.

### Exact local criterion, including the central variables

Fix g in F_h intersect G_reg, with K=Z_G(h). Use the smooth local
quotient Q=K//K and the logarithm-corrected potential W_h from
[the critical-presentation theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md).
Its variables include BOTH the derived-group characters and the central
torus parameters. Holding the latter fixed would give the wrong fiber.
Translate the critical point to zero and subtract the constant. Write

    W(v,y)=(1/2) H_N(y,y)+(1/6) T_3(v+y,v+y,v+y)
           +(1/24) T_4(v+y,v+y,v+y,v+y)+O(5),

where R=ker Hess(W)(0), v belongs to R, y belongs to a complement N,
and H_N is invertible. Here T_j denotes the jth symmetric derivative
in the chosen coordinates. The formal implicit-function theorem gives
a unique y=y(v), of order at least two, solving partial_y W=0. Put

    F(v)=W(v,y(v)),       b_v=T_3(v,v,-)|_N.

Then

    y_2(v)=-(1/2) H_N^(-1)b_v,
    F_3(v)=(1/6) T_3(v,v,v),
    F_4(v)=(1/24) T_4(v,v,v,v)
                    -(1/8) b_v(H_N^(-1)b_v).            (34)

For the quartic, the quadratic term contributes
+(1/8)b_v(H_N^(-1)b_v) and the mixed cubic term contributes
-(1/4)b_v(H_N^(-1)b_v). Their sum proves the correction in (34).
Terms involving y_3 cancel or have higher degree. In particular the
raw fourth derivative restricted to R is NOT the effective quartic.
The displayed expression uses the chosen splitting; when F_3 vanishes,
F_4 is the leading obstruction, invariant up to invertible linear change
of the residual coordinates. Do not assert such invariance of F_4 alone
when F_3 is nonzero.

There is an exact isomorphism of completed local rings

    completed O_(F_h,g)
      = (C[[v]]/(partial F/partial v_i))
           [[u_1,...,u_(dim K-r)]].                     (35)

Indeed partial_y W generates the ideal (y-y(v)), and the chain rule
identifies the remaining equations with partial_v F. Consequently:

    F_h is smooth at g  if and only if  F is constant.  (36)

To see the reverse implication as well as the forward one, F has order
at least three. Thus its Jacobian ideal J lies in (v)^2. If J is nonzero,
C[[v]]/J has embedding dimension dim R but Krull dimension strictly
smaller than dim R, since C[[v]] is a domain. It is not regular. In
characteristic zero J=0 means that F is constant. Completion and the
smooth orbit variables preserve this conclusion.

Any first nonzero F_d, d>=3, therefore certifies a nonsmooth actual
fiber. This need not mean an isolated nilpotent thickening: a critical
scheme can also be reduced and singular. Vanishing of F_3 and F_4 is
only necessary, not sufficient. For example y^2/2+v^6 has both zero,
but its residual Jacobian algebra is C[[v]]/(v^5). Even for a polynomial
W, elimination can produce an infinite residual series.

### Direct relation to the adjoint trace formula

There is a coordinate description valid without choosing quotient
coordinates. On the full formal neighborhood g exp(X) in K,

    W_h(X)=tr(U exp(ad X))-tr(U)-B(h,X),   U=Ad(g).

The logarithmic correction is linear here because h annihilates [k,k].
The defining trace equation cancels the linear term. Its higher diagonal
derivatives are exactly

    T_m(X,...,X)=tr(U (ad X)^m),   m>=2.                 (37)

Polarization gives all mixed derivatives. One may apply (34) on this
FULL neighborhood, retaining orbit directions in the Hessian radical.
One may alternatively use (34) on the smooth quotient as in (35).
These are distinct coordinate computations: higher derivatives of the
quotient potential are not obtained by simply restricting (37) to an
arbitrary linear complement to a conjugacy orbit. This caution matters
particularly at regular nonsemisimple sources.

### Checks against the existing obstruction and excess fibers

In the Sp4/GL2 triple point above, put z=z_0 exp(y), delta=2i sqrt(3),
c=-delta and W=A(z)(x^2-1)+1-cy, up to its constant. Then

    H_yy=4,  T_3(x,x,y)=2delta,  T_4(x,x,x,x)=0,
    y_2=-(delta/4)x^2,  F_3=0,  F_4=(3/2)x^4.          (38)

Thus partial_x F=6x^3+O(x^5) and the transverse algebra is exactly
C[[x]]/(x^3), agreeing with the earlier elimination of the ORIGINAL
fiber equations. The raw quartic would miss this obstruction entirely.
The independent expansion uses D_y W=0, not (34), to recover y_2
and then substitutes it into W.

In the known E6 contact A5 excess fiber, the character is

    Phi=x_1 x_5+z^2+z^(-2)+(z+z^(-1))x_3.

For fixed nonzero central target c, take z_0=i, x_1=x_5=0,
x_3=c/(2i). In variables (x_1,x_5,x_3,log(z/z_0)) the Hessian is
the direct sum of [[0,1],[1,0]] and [[0,2i],[2i,-8]]. It is invertible.
The remaining variables x_2,x_4 do not occur in W AT ALL. The eliminated
solution is independent of them, so F is constant to EVERY order.
This recovers smooth excess of dimension two, not merely a vanishing
quartic. The other central lift gives the same local conclusion.

For comparison, the altered 82-dimensional Levi representation above
has, after completing a,b and solving for X, the exact relative potential

    -S^2-2-S^2 u^2 v^2/(4(1-uv)).

At a generic central point z_0, the Hessian of -S^2-2-c log(z/z_0)
is nonzero. The further central elimination changes the potential only
in degree at least eight in u,v. Its effective quartic is therefore
-S(z_0)^2 u^2 v^2/4, nonzero when S(z_0)!=0. This correctly detects
the already proved persistent nonsmooth structure. It remains a control
representation of a reductive Levi, NOT an adjoint counterexample.

### What this attempt does and does not establish

The existing audit checks (34) by a symbolic substitution with two
kernel and two complementary variables, (38) against the original
Sp4 local ideal, and the exact absence of all residual E6 excess terms.
The control also shows why finite-order vanishing is insufficient.

None of (34)-(37) forces F to be constant for the adjoint character.
Indeed actual adjoint exceptional targets already have nonconstant F.
To finish the generic conjecture one must prove that EVERY incidence
component carrying a nonconstant residual potential has nondominant
central image. The higher-order calculation supplies no such image bound.
In particular a nonzero cubic or quartic at a special point cannot be
promoted to a generic counterexample, and checking arbitrarily many
vanishing jets at selected points cannot prove the desired theorem.

This closes the proposed local higher-order diagnostic, not the uniform
proof. Another jet calculation is not the next goal. The remaining
structural task is to rule out a dominant bad component using the full
adjoint character, or exhibit one. No new result ID or coverage is added.

## 25 October: the actual generic-target algebra and the surviving obstruction

The next attempt works over the function field of the CENTRAL TARGET,
not over the function field of the source central torus. These operations
are different. This section supplies a precise algebraic reduction and
tests two shortcuts; it does not establish the adjoint-specific theorem.

Put S=[K,K], ell=rank S, s=dim Z(K)^0. On the finite central etale cover
S times Z(K)^0 -> K, use fundamental characters x_1,...,x_ell and central
torus coordinates z_1,...,z_s. Let Phi be the FULL restricted adjoint
character, and set D_j=z_j partial_z_j. With independent target variables
c_1,...,c_s, define

    L=C(c_1,...,c_s),
    A_K=L[x_1,...,x_ell,z_1^(+-1),...,z_s^(+-1)] /
          (Phi_x1,...,Phi_xell,D_1 Phi-c_1,...,D_s Phi-c_s). (39)

The coordinates c_j are paired with the chosen central cocharacters by
the ambient Killing form; this is an invertible linear change from h.
The generic h automatically has Z_G(h)=K: each root outside K restricts
to a nonzero linear form on z(k). Every quotient point has a regular
representative of K, and on the fiber regularity in K equals regularity
in G. Thus all points of Spec A_K matter, including quotient points of
regular NONSEMISIMPLE representatives. We do not discard a torus
discriminant or make a point-set replacement of the ideal.

The conjecture for this K is equivalent to

    A_K is smooth over L, including A_K=0.               (40)

Proof: the smooth regular quotient and the finite central etale cover
preserve and reflect smoothness. If the generic fiber is smooth, generic
flatness and openness of the smooth locus, followed by constructibility
of the image of its complement, give a dense target open over which the
whole morphism is smooth. Conversely such an open contains the generic
point. A constructible image avoiding the generic point of an irreducible
base is not dense, which is the needed spreading step; properness is not
required. Over the perfect field L, a finite-type algebra is smooth iff
all its local rings are regular. These are standard facts, not a new
Springer theorem; see [generic flatness](https://stacks.math.columbia.edu/tag/0529)
and [smooth algebras over fields](https://stacks.math.columbia.edu/tag/00TQ).

### What generic localization legitimately removes

Before base change to L the ring is a C[c]-algebra. This localization
inverts nonzero polynomials in the TARGET c. It removes components with
nondominant target image. It does not invert arbitrary polynomials in z
or x, remove a discriminant component of the source, or take a radical.
After substituting c=D Phi on the incidence scheme, it can of course
invert such a target polynomial expressed as a source function; the
authorization is its origin in C[c], not its appearance as a denominator.

The known E6 contact model makes the distinction exact. Over C(c) its
original equations reduce to

    (x_1,x_5,z^2+1,2z x_3-c),
    A_K = C(c)[z,x_2,x_4]/(z^2+1).                      (41)

This is smooth of dimension two and nonempty. Working instead over
C(z,c) makes z^2+1 a unit and incorrectly returns the zero ring. A
generic target can be supplied entirely by a special central position
of the SOURCE. This is the existing excess example, not new coverage.

### The finite central algebra can be split, but not silently reduced

By the central-rigidity theorem, B=image(L[z_1^(+-1),...,z_s^(+-1)] -> A_K)
is finite-dimensional over L. This follows either by base change of the
generic central-finiteness result or by applying its algebraicity argument
over an algebraic closure of L and descending the resulting relations.
If A_K is nonzero, its finite central algebra B has a nilpotent radical N
and B/N is a finite product of separable extensions of L. There is a unique
L-algebra section B/N -> B of the quotient map: B/N is formally etale,
so its identity lifts uniquely across the nilpotent ideal. See
[formally etale maps](https://stacks.math.columbia.edu/tag/00UP).

Consequently the central residue fields and their idempotents can be
lifted inside the ACTUAL ring. This gives A_K=product_i A_i, with A_i
an algebra over a finite separable field E_i/L. In each factor write

    z_j = zeta_ij (1+epsilon_ij),

where zeta_ij belongs to E_i^* and epsilon_ij is nilpotent in A_i.
This does not set epsilon_ij to zero. Smoothness of A_K is equivalent to
smoothness of every A_i over E_i and would imply that all these epsilon
vanish. The converse implication from their vanishing is FALSE.

In particular, separability of the finite central image is necessary
but is not enough to prove smoothness of the derived-character directions.
For a finite B, reducedness can be tested by squarefreeness of the minimal
polynomials over L of all its central generators. This test says nothing
by itself about the rest of A_K.

### A dominant obstruction can hide at component intersections

The previously audited branch of the altered 82-dimensional Levi
representation gives a decisive control for the proposed generic-point
argument. On this branch the central equation over L=C(c) is

    P(z)=6z^12+c z^6-6=0.                               (42)

It is squarefree over L, and its roots have S=z^3+z^(-3)!=0. After a
finite separable extension selecting one root and completion at the
previously identified singular point, the transverse algebra is

    R=E[[u,v]]/(u^2 v,u v^2).                           (43)

The central variables here have no residual nilpotent thickening: the
original central equation was replayed in the nonreduced ring, and
P'(z) is a unit. Nevertheless (43) is not smooth. More precisely,

    (u^2 v,u v^2) = (u) intersect (v) intersect (u^2,v^2),
    sqrt(0)_R=(uv),  (uv)^2=0,  Ann_R(uv)=(u,v).        (44)

Both minimal components are smooth and R is reduced at their generic
points: invert u and the ideal becomes (v); invert v and it becomes (u).
The nonzero nilpotent uv is supported at their intersection, an embedded
associated point. At that point the tangent dimension is two while the
Krull dimension is one. Even R_red is the singular crossing E[[u,v]]/(uv).
The smooth orbit factors from the earlier group model do not change this.

Thus ALL of the following can hold simultaneously over a generic target:
a reduced finite central image, a generically reduced local fiber scheme,
and smooth generic points of each of its minimal components, yet a
nonsmooth actual fiber. No statement about other branches of that
representation is needed. There is no target polynomial that removes
this intersection: the
calculation already takes place over C(c). This remains a representation
of a reductive Levi, not an adjoint representation of a simple group, and
therefore does not disprove the conjecture.

### Result of the attempt

The exact audit uses original equations over Q(c), whose certificates
remain valid over C(c). It verifies (41), the error from source-coefficient
localization, squarefreeness in (42), and both the primary decomposition
and localizations in (44). It also checks the SAME Sp4 model over Q(c)
by a unit Jacobian ideal and independently by separating its two branches.
The special-target triple point disappears legitimately under target
localization. No additional Levi cases were computed.

This leaves an exact unresolved assertion: for the FULL ADJOINT Phi,
every local ring of (39), including points on intersections of dominant
components, must be regular. The root/coefficient lemmas presently proved
do not exclude the embedded or crossing obstruction in (43). Working over
the generic target does not supply that missing argument. We have neither
a proof nor an actual adjoint counterexample, and add no theorem ID.

The next proof target should include intersections and embedded structure
from the outset. A root-theoretic description of a smooth reduced component
alone, or of the finite central algebra alone, cannot finish (40).

## 25 October: reducedness attempt and the missing depth argument

The weaker target is now that A_K in (39) is REDUCED for every simple
simply connected G and proper Levi K, with Phi the full adjoint character.
It is not asserted as a theorem. It would exclude nilpotents but would
still permit reduced singular intersections, so it is weaker than (40).

There are two genuinely different requirements. By the standard
[R_0 plus S_1 criterion](https://stacks.math.columbia.edu/tag/031R), a
Noetherian ring A is reduced exactly when:

* R_0: A_p is a field at each minimal prime p (multiplicity one at each
  component's generic point);
* S_1: A has no embedded associated primes (equivalently, the required
  depth is at least one at every positive-dimensional localization).

This criterion is a published input, not a new result. A proof for the
adjoint character must supply BOTH properties; neither has been proved
uniformly by the current root arguments.

### Why the proposed complete-intersection argument does not close this

A local complete intersection is Cohen--Macaulay and hence satisfies S_1.
It need not satisfy R_0. The already recorded abstract control
Phi=x^3+(z-1)y has generic fiber

    C(c)[x,y,z^(+-1)]/(x^2,z-1,zy-c) = C(c)[x]/(x^2).

This is a complete intersection with no embedded primes, but its minimal
component has multiplicity two. Conversely the existing altered-Levi
branch (43) satisfies R_0 and fails S_1. These controls separate the two
obligations; neither is an actual adjoint counterexample.

Nor can one deduce a regular sequence merely by counting the displayed
gradient equations. The E6 excess model has six ambient quotient variables
and six displayed equations, but its generic fiber has dimension two.
Two derived-character derivatives are identically zero, and its ideal
has the four independent generators in (41). It IS a smooth complete
intersection after removing redundancies. It is not a complete intersection
of the six displayed equations in the expected codimension six. A general
argument must establish the necessary syzygies and the ACTUAL codimension,
including every excess branch, before using the depth conclusion.

The smooth coefficient-boundary theorem GEO-75 does not give those
syzygies or a resolution of the full coupled gradient ideal. The altered
Levi branch has exactly the same nonlinear coefficient data. Thus applying
Cohen--Macaulayness or unmixedness to that data would be an unjustified step.

### An exact target-polynomial certificate for the weaker goal

There is a useful way to formulate the missing root-dependent calculation
without recomputing all target fibers. Let

    P=C[x_1,...,x_ell,z_1^(+-1),...,z_s^(+-1)],
    J=(Phi_x1,...,Phi_xell),  R=P/J,
    N=sqrt(J)/J subset R,    mu_j=D_j Phi in R.

Regard R as a C[c_1,...,c_s]-algebra by c_j -> mu_j. Then (39) is exactly
the localization of R at all nonzero target polynomials. Nilradicals
commute with localization, so

    A_K reduced
      iff N tensor_(C[c]) C(c)=0
      iff some 0!=p in C[c] satisfies p(mu) N=0
      iff p(mu) sqrt(J) subset J for some 0!=p in C[c]. (45)

For the middle implication, take a finite set of R-module generators of
N; each is killed by some nonzero target polynomial after the localization
vanishes, and the product kills N. The converse is immediate. Finiteness
as a C[c]-module is NOT assumed or needed. Empty generic fibers are allowed.
Characteristic zero makes the target function field perfect, so reducedness
here also gives geometric reducedness of the generic fiber.

In particular, a radical relative incidence ideal J would suffice, even
if that incidence were singular. It is stronger than necessary: nilpotents
confined to special target values may be killed by p. For example the
abstract potential Phi=(z-1)^2 x+z has J=((z-1)^2) and

    mu=2z(z-1)x+z,   (mu-1)(z-1)=0 in R.

Here p(c)=c-1 kills the nilradical and the generic fiber is empty. By
contrast, on the branch (43), no nonzero target polynomial kills uv,
since the branch already exists over C(c). One must not replace sqrt(J)
by J or quotient by N to manufacture the desired certificate.

The known Sp4 relative incidence ideal (A(z)x) is radical because A is
squarefree, despite its nonreduced SPECIAL fiber (9). This illustrates
the distinction: generic target localization does not create nilpotents,
whereas specialization can. The E6 contact relative incidence ideal
(x_1,x_5,z^2+1) is also radical. These are consistency checks of (45),
not additional adjoint cases or a proof for the remaining root systems.

### Outcome and stopping point of this method

The exact audit checks the two independent R_0/S_1 failures, the vertical
annihilator p=c-1, and the actual Sp4/E6 incidence calculations. It does
not certify either R_0 or S_1 for all adjoint characters.

The attempted depth proof stops at the absence of an adjoint-specific
resolution or regular-sequence argument for J after target localization;
generic multiplicity one is also not supplied by the trace identities.
Equivalently, no root-dependent polynomial certificate (45) has been
constructed uniformly. Thus the weaker reducedness assertion, as well
as the original smoothness conjecture, remains unproved. The criterion
(45) is elementary commutative algebra, not a new Springer theorem or
novelty claim. No coverage or ledger result ID is added.

Further formal reformulations of reducedness would not resolve this gap.
A useful next attempt must construct (45) from the FULL root character,
or find a surviving nilpotent for an actual adjoint character. A proof
about nonlinear coefficients alone has already been ruled out.

## 25 October: bounded search for an actual adjoint nilpotent

This search tests two existing E8 adjoint characters for a counterexample
to the UNIFORM conjecture. It is not a resumption of the exceptional
coverage census. The full root characters and the archived unit
eliminations are reconstructed independently by
[the search audit](tools/audit_adjoint_nilpotent_search.py). No changed
linear multiplicities, non-adjoint Levi representation, or special-target
singularity is admitted as a counterexample.

### A candidate ruled out with its entire central boundary retained

Take the archived E8 Levi on nodes [1,2,3,4,6,8], with derived type
A4+A1+A1. After the previously audited endpoint-pair elimination, let
the four affine coordinates be a_old,b_old,x,y and the central torus
coordinates z,w. Make the unit coordinate changes

    a=a_old*z, b=b_old/z, r=z^5*w^2.

The map (z,w)->(r,w) is finite etale of degree 5. The residual character is

    Phi = a*(x+w*y+r^-1) + b*(x+w^-1*y+r) + H,

    H = -w^2-w^-2-3
        -x^2*((w*r+(w*r)^-1)*y+y^2)
        -x*(w^2*r+(w^2*r)^-1+(r+r^-1)*y^2+2*(w+w^-1)*y)
        -(r/w+w/r)*y.

This is an equality with the actual archived adjoint character after
the earlier elimination, not a family with independent coefficients.
The affine equations for a,x eliminate them globally:

    x=-w*y-r^-1,  a=-b-H_x(-w*y-r^-1,y,w,r).

All substitutions are over the Laurent polynomial ring. In particular,
w^2-1 has NOT been inverted. The remaining potential and relative ideal are

    Psi=b*F+h,  F=A*y+r-r^-1,  A=w^-1-w,
    h=H(-w*y-r^-1,y,w,r),
    J=(F, b*A+h_y) in C[b,y,w^(+-1),r^(+-1)].             (46)

The audit verifies every affine derivative and both logarithmic central
derivatives. Before imposing affine equations, the original central
derivatives transform as

    D_z(old Phi)=5*D_r(Phi)+a*Phi_a-b*Phi_b,
    D_w(old Phi)=D_w(Phi)+2*D_r(Phi).

Thus on the critical scheme the target coordinates undergo an invertible
linear change. The subsequent elimination also preserves these derivatives
by the chain rule. The target function field is therefore unchanged up to
an invertible linear reparametrization.

**Exact conclusion for this candidate.** The entire relative critical
surface Spec(C[b,y,w^(+-1),r^(+-1)]/J) is smooth and integral.

Proof. On A!=0 the two equations solve uniquely for y and b. This open
is isomorphic to the central torus open A!=0. On A=0 they give
w=epsilon, r=eta with epsilon,eta in {1,-1}, and

    h_y = -2*(2*y+epsilon*eta)*(y^2+epsilon*eta*y-1).

Each of these four cubics is squarefree; the audit supplies exact Bezout
identities with its derivative. The Jacobian minor in columns (r,y) of
(F,b*A+h_y), restricted to this boundary, is twice the derivative of
that cubic. It is nonzero at every boundary solution. Hence the whole
surface is smooth of dimension two. The boundary is a finite disjoint
union of b-lines, of dimension one, so it contains no irreducible
component. The irreducible open A!=0 is consequently dense in the whole
surface, proving integrality. Equivalently, F is primitive linear in y,
the two equations are a regular sequence, and the dimension-one boundary
cannot contribute a component or an embedded prime.

The finite etale pullback to the original central cover remains smooth
and reduced (it need not remain integral). Localization at the nonzero
TARGET polynomials cannot create nilpotents. This candidate therefore
cannot supply a nilpotent counterexample, including on the central
boundary that a source-field calculation would discard. Generic
smoothness in characteristic zero also gives smooth-or-empty generic
target fibers here. This is a bounded candidate exclusion, not an
all-E8 or all-simple-group theorem.

The algebraic inputs are standard: complete intersections are
[Cohen--Macaulay](https://stacks.math.columbia.edu/tag/00S8), such rings
have [no embedded associated primes](https://stacks.math.columbia.edu/tag/00N2),
and [R_0 plus S_1 gives reducedness](https://stacks.math.columbia.edu/tag/031R).
In this example the explicit Jacobian argument is stronger and avoids
needing an unproved codimension assertion. No literature-priority claim
is made for this candidate calculation.

### A second candidate: source-generic screening is only a partial test

The second selected character is E8 on [1,2,3,5,6,8], of derived type
A2+A1+A2+A1. Its four-variable quartic residual is screened by
[this reproducible script](tools/search_adjoint_nilpotents.py). With
original residual coordinates x,u_old,v_old,y and central z,w, use

    u_old=u/w, v_old=v*w, t=z^3, r=z^6*w^3.

The central map is finite etale of degree 9. Unlike (46), this screen
then works over Q(t,r), explicitly inverting source polynomials. It
can detect a nilpotent on that open, but a negative answer there CANNOT
exclude a nilpotent component supported on a central boundary. Such a
component could still dominate the target through its free affine
coordinates. Any positive candidate would also require the actual
central target equations and a dominance certificate.

The initial full radical computation stopped without a certificate.
It is not counted as evidence of radicality. The smaller tests and their
precise scope are recorded in the JSON reports generated by the script.

The replacement check DOES certify that this source-generic relative
critical algebra is **finite etale of rank 38 over Q(t,r)**, hence also
over C(t,r). It uses the following exact argument, rather than inferring
a generic result from a single numerical fiber:

1. The generic monic Groebner basis has 28 elements and 38 standard
   monomials. Every coefficient is regular at (t,r)=(2,3).
2. Specializing the entire basis gives exactly the independently computed
   reduced basis of the four derivatives at (2,3), with the same leading
   monomials. The monic basis defines a finite free model of the GENERIC
   algebra over a parameter localization containing (2,3). This statement
   does not assume that the original incidence family equals that model
   at every other central point.
3. At (2,3), the ideal of the four derivatives together with their Hessian
   determinant is the unit ideal over Q. The determinant therefore has
   nonzero multiplication norm in the finite free model; over Q(t,r) it
   is a unit. The zero-dimensional Jacobian criterion gives finite
   etaleness of the generic algebra.

The Sage calculation checks all three assertions. An independent SymPy
replay reconstructs the character transformation, recomputes the 28-element
specialized basis modulo the prime 65521, counts its 38 standard monomials,
recomputes the Hessian determinant, and obtains the unit Jacobian ideal.
The rational basis coefficients are defined at this prime and the leading
monomials agree. This is an independent modular corroboration of the exact
rational specialization certificate. A separate rational radical
calculation at (2,3) also agrees.

Here 38 counts geometric points of this RELATIVE critical algebra with
generic CENTRAL SOURCE coordinates. It is not the degree of the Springer
map and not a computation of the generic CENTRAL TARGET algebra (39).
No central boundary has been certified for this second model. Consequently,
any nilpotent counterexample in it must be supported on a proper central
source locus; deciding whether such a locus supports a nilpotent that
survives target localization remains the unperformed part of this search.

**Search outcome:** no actual adjoint nilpotent counterexample has been
certified. The uniform reducedness and smoothness assertions remain open.
This does not repair the missing full-character argument identified
above, and does not justify expanding the search into a case census.

## 25 October continuation: boundary exclusions and the remaining gap

Continue ONLY with the second E8 candidate above. Write

    R=C[x,u,v,y,t^(+-1),r^(+-1)]/J,  J=(Phi_x,Phi_u,Phi_v,Phi_y),
    mu=(t*Phi_t,r*Phi_r),  L=C(c,d),
    A=R tensor_(C[c,d]) L,  c->mu_1, d->mu_2.

These are the full adjoint potential and the same finite central cover as
in the source-generic calculation. In particular, t and r are SOURCE
coordinates. No equation in them is discarded merely because c,d are
generic. This continuation has NOT ruled out the candidate completely.

### Six entire source boundaries cannot dominate the target

Exact calculations over Q (or Q(i), with both conjugate embeddings)
give dimension one for each of the following entire closed relative
critical schemes:

    Spec(R/(t-r)),        Spec(R/(t+r)),
    Spec(R/(t*r-1)),      Spec(R/(t*r+1)),
    Spec(R/(t^2+1)),      Spec(R/(r^2+1)).                 (47)

These are global calculations along the curves, including their special
points and all affine critical coordinates. They are not checks at
generic points of those curves. The direct exact Sage computations are
recorded in [the six-boundary report](adjoint-simple-boundaries.json).

Consequently, in the GENERIC TARGET algebra A the element

    b=(t^2-r^2)*((t*r)^2-1)*(t^2+1)*(r^2+1)               (48)

is a unit. Thus no nilpotent supported on these six boundaries can
survive target localization. This conclusion includes nilpotent
thickenings; it does not require the boundary schemes to be reduced.

Proof of the target conclusion. The closure of the image of a finite-type
scheme of dimension at most one in A^2 has dimension at most one, by
the standard [dimension formula](https://stacks.math.columbia.edu/tag/02JT).
For each boundary f=0 there is therefore a nonzero target polynomial
p_f(c,d) whose pullback lies in sqrt(J+(f)). Some power of that pullback
lies in J+(f). After tensoring with C(c,d) it is a unit, so A/(f)=0
and f is a unit in A. Taking the product gives (48). This is a standard
dimension/Nullstellensatz argument applied to exact incidence ideals,
not a new general dimension theorem or an explicit formula for p_f.

For independent verification, the
[certificate generator](tools/verify_adjoint_boundary_exclusions.py)
constructs polynomials in the ORIGINAL five-generator ideals, including
the torus equation z*q-1 after parametrizing a boundary. For every
coordinate two-plane, one witness has leading monomial supported in
that plane. Hence the leading-monomial ideal has dimension at most one.
This proves an upper dimension bound without assuming the selected
witnesses themselves form a Groebner basis. Each witness comes with
its explicit ideal-membership coefficients. The
[Sage-free replay](tools/audit_adjoint_boundary_exclusions.py) reconstructs
the specialized original gradients, checks every membership identity,
and checks the monomial support bound. Its report states exactly which
boundaries have completed this independent replay; the direct dimension
calculations and independent-certificate status must not be conflated.
At that stopping point only t=r and t=-r had independent replays.
The subsequent three-path continuation below completes all six:
60 original-ideal memberships on four boundaries, and two exact
original-potential involutions. See [the completed audit](adjoint-boundary-exclusions-audit.json).

### The original-ideal boundary is larger than the basis denominators

The earlier source-generic reduced basis has nine nonunit denominator
factors, three beyond the six in (47). That list alone is NOT a certified
boundary for the original ideal: expressing its basis elements in terms
of the original generators can introduce additional denominators.

A tracked standard-basis calculation now verifies both ideal inclusions
after localization. The 28 monic basis elements have explicit membership
identities in the four original affine gradients, and those original
gradients reduce to zero by the basis. The certificate is saved in
[this archive](adjoint-central-membership-forward.json.gz); its denominator
list and scope are in [the membership report](adjoint-central-boundary-lift-forward.json).
The calculation is reproducible with
[the exploration script](tools/explore_adjoint_central_boundary.py).

Let D be the squarefree product of the NONUNIT factors in that report.
There are sixteen, including the six already excluded. Over
C[t^(+-1),r^(+-1),D^-1] the original relative algebra equals its monic
standard-basis model and is finite free of rank 38. Its generic algebra
is reduced by the preceding finite-etale certificate. Torsion-freeness
therefore makes R[D^-1] reduced: a nilpotent maps to zero in the generic
algebra, and the map into that algebra is injective. Target localization
preserves this conclusion.

Thus any surviving nilpotent must be supported on one of the remaining
TEN source curves in that specific certificate. The count is a statement
about the chosen elimination certificate, not a geometric classification
of bad fibers: several factors can be artifacts of the basis computation.
Their intersections are included and must not be silently deleted.

This is the precise remaining gap. Whole-incidence and nonlinear-boundary
Groebner attempts did not finish within their bounded runs and produced no
certificate. No conclusion is inferred from those incomplete calculations.
Neither the candidate's full generic-target reducedness nor the uniform
adjoint conjecture is proved or disproved by this continuation.

The proposed next proof target was flatness of this relative incidence algebra
over the central open b!=0. Together with its reduced generic algebra,
that would imply reducedness there, and (48) would then settle the
candidate's generic-target reducedness at once. Flatness on that larger
open is a PROPOSED target, not an established consequence of the rank-38
calculation. It may fail, in which case its failure must be analyzed in
the actual target algebra rather than removed by source localization.

## 25 October flatness attempt: a genuine non-finiteness obstruction

The central open (48) has not been proved flat. The following calculation
does rule out the stronger FINITE model on that entire open, and explains
why a rank-38 basis cannot extend across every allowed source parameter.
All statements here concern the second full E8 adjoint model, on the same
finite central cover and with the original four affine derivatives.

### Exact fibers and what their lengths prove

At (t,r)=(1,2),(-1,2),(2,1),(2,-1),(2,3), the critical algebra has
dimension zero and length 38 over Q. At r=2, let t be any root of either
irreducible polynomial

    h_+(T)=T^4+3*T^3+(17/4)*T^2-3*T+1,
    h_-(T)=T^4-3*T^3+(17/4)*T^2+3*T+1.                (49)

For EACH such geometric central point, the critical algebra has dimension
zero and length 37. All these points satisfy b!=0. These are exact
number-field calculations with the original derivative ideals, not
numerical roots or counts of distinct solutions. No claim that the
specialized algebras are reduced is needed here.

The corresponding two central curves, before setting r=2, are

    (t-t^-1 + r-r^-1)^2+4=0,
    (t-t^-1 - r+r^-1)^2+4=0.                           (50)

Only the stated points on these curves have been checked: length 37 on
their entire generic points is NOT established by the specialization.
A third shared denominator curve gives length 38 at r=2 and t a root of

    T^8-(25/18)*T^6+(4937/144)*T^4-(25/18)*T^2+1.

The [exact screen](adjoint-central-flatness-curve-screen.json) saves the
number-field bases. The [independent audit](adjoint-central-flatness-audit.json)
reconstructs the original potential in SymPy, recomputes the specialized
bases at good finite-field images, compares the entire bases, and counts
standard monomials. It obtains 37,37,38, agreeing with the rational
number-field computations. This is independent MODULAR corroboration;
the characteristic-zero length statements come from the exact Sage
number-field bases. The [test script](tools/test_adjoint_central_flatness.py)
and [Sage-free auditor](tools/audit_adjoint_central_flatness.py) retain
the distinction.

**Non-finiteness deduction.** Spec(R[b^-1]) is not finite over ANY central
Zariski neighborhood containing one of the points (49). Indeed, if it
were finite there, its algebra would be a finite module M over the
corresponding integral local base. The generic rank is 38, whereas
dim(M/mM)=37. Nakayama would give 37 local generators, which after
passing to the fraction field would generate a 38-dimensional vector
space. This is impossible. This is an elementary finite-module argument,
not a new general finiteness theorem. Thus the loss of length is real;
it cannot be repaired by choosing a different finite free basis there.

The conclusion is NOT non-flatness. For example,
C[s] x C[s,s^-1] is a flat, finitely presented C[s]-algebra with generic
fiber length two and special fiber length one at s=0. Finite presentation
and finiteness as a module must not be confused. Nor do these source
lengths count Springer fibers over a target.

### The valid flatness lemma, and its unverified hypothesis

Here is the standard algebraic input in its useful precise form. Let
S=C[t^(+-1),r^(+-1),b^-1], P=S[x_1,...,x_m], and R=P/(f_1,...,f_m).
Assume R tensor_S Frac(S) is reduced. Then at every point of a
zero-dimensional closed fiber of Spec(R)->Spec(S), R is flat over S
and its local ring is reduced.

Proof. At a closed point p over s=(t_0,r_0), the m equations together
with t-t_0,r-r_0 form a system of parameters in the regular ambient
local ring of dimension m+2. A system of parameters in a regular local
ring is a regular sequence. Consequently R_p is Cohen--Macaulay of
dimension two, and the two base parameters form an R_p-regular
sequence. Miracle flatness gives flatness over S_s. Flatness over the
domain S_s makes the map to the generic algebra injective, so generic
reducedness implies reducedness of R_p. Checking closed points suffices
for this finite-type setting. More globally, if EVERY nonempty fiber
has dimension zero, the displayed presentation is a relative global
complete intersection, hence flat by the standard
[relative complete-intersection theorem](https://stacks.math.columbia.edu/tag/00SW).

In particular, the total relative scheme is flat and reduced at the
points above the explicitly checked fibers. This does not assert those
fibers themselves are reduced, or that all remaining fibers are zero
dimensional. A complete-intersection argument requires the fiber
dimension hypothesis; four displayed equations alone do not supply it.

### What the second presentation does and does not cover

A tracked basis in reverse variable order also verifies both original
ideal inclusions on its own denominator open. See the
[reverse certificate report](adjoint-central-boundary-lift-reverse.json)
and [membership archive](adjoint-central-membership-reverse.json.gz).
After removing b and torus units, the forward and reverse denominator
lists share three curve factors: the two in (50), up to Laurent units,
and the degree-twelve factor in the earlier report. Away from these
curves the two certified opens cover all but a finite set of central
points. They do NOT cover that finite set automatically.

An exact two-variable saturation calculation retains the uncovered
intersections: 22 nonempty pairwise intersection schemes, with overlaps
possible. This is a count of presentations of finite intersection
schemes, NOT a count of distinct points. The equations are recorded in
[the intersection report](adjoint-central-flatness-cover-intersections.json).
This reduces the geometric gap to three entire curves plus that finite
set; it does not remove any of them from the actual target algebra.
The reverse membership check and these saturation computations are exact
Sage calculations, without a completed independent identity replay.

Whole-family standard-basis and comprehensive-basis attempts, a third
tracked ordering, and a generic-curve calculation did not complete in
their bounded runs. They supplied no additional certificate. In
particular, the samples (49) cannot be promoted to generic or global
statements on (50).

**Outcome at that stage and next gate.** No non-flat fiber or adjoint nilpotent has
been exhibited; flatness on b!=0, generic-target reducedness, and the
uniform smoothness conjecture all remain open. The finite rank-38
extension route is ruled out at the specified points. To advance this
flatness route, one must prove that EVERY remaining allowed central
fiber is zero dimensional (including isolated exceptions), or exhibit
a positive-dimensional fiber and then analyze its target image. More
isolated length computations do not meet that gate. This is a proof-route
audit of the existing candidate, not new exceptional coverage or a claim
of a new general theorem.

## 25 October structural completion: every allowed critical fiber is finite

**Theorem, with exact scope.** Let G be complex simply connected E8,
rho=Ad, and K the proper standard Levi with retained nodes
{1,2,3,5,6,8}, of derived type A2+A1+A2+A1. Use the finite central
cover, unit eliminations, and full residual adjoint potential Phi(x,u,v,y;t,r)
specified in the preceding search and
[exact potential record](adjoint-nilpotent-source-screen.json).
For EVERY (t,r) in (C^*)^2 satisfying

    b_0=(t^2-r^2)*((t*r)^2-1)*(t^2+1)*(r^2+1) != 0,

the critical algebra C[x,u,v,y]/(Phi_x,Phi_u,Phi_v,Phi_y) is finite
dimensional, possibly zero. Equivalently, every nonempty critical fiber
on this open has dimension zero. This includes all exceptional curves
and their intersections; none is removed by an additional localization.

Consequently the entire relative critical algebra is flat and reduced
over this open. Combining this with the earlier six whole-boundary
exclusions proves that the generic CENTRAL TARGET algebra A in (39)
is finite etale or zero. In particular, (40) holds for this K. The
generic-target degree, nonemptiness, all-special-target smoothness, and
the uniform all-simple conjecture are not asserted here.

### A squarefree plane-pencil argument

Make the invertible affine-variable change u=a/t, v=t*b, and call the
resulting potential f(x,a,b,y). In this section b is an AFFINE variable;
the allowed central factor is denoted b_0 to avoid ambiguity. Set

    T=t+t^-1,  R=r+r^-1,  z=(a,b)^tr,
    f=(1/2) z^tr H z + L^tr z + Q,
    D=det(H),  N=Q*D-(1/2)L^tr adj(H)L.

Here H and L depend on x,y and the central parameters; Q is independent
of a,b. Explicitly,

    H_11=-2*(1+x/t)*(1+y/r),
    H_22=-2*(1+t*x)*(1+r*y),
    H_12=H_21=-(x^2+y^2+T*x+R*y+1).

Direct expansion of the FULL potential gives the following polynomial
identities, verified independently in characteristic zero:

    deg(N)=7,
    N_7=x^2*y^2*(x-y)*(x+y)*(T*x-R*y),
    deg_x(N)=deg_y(N)=6,  [x^6]N=[y^6]N=-1,
    deg_x(D)=deg_y(D)=4,                               (51)
    N-f*D=-(1/2)(H*z+L)^tr adj(H)(H*z+L).              (52)

The coefficients [x^6] and [y^6] here are constant polynomials in the
OTHER variable. On b_0!=0, T and R are nonzero and T!=R,-R, since

    T^2-R^2=(t^2-r^2)*((t*r)^2-1)/(t^2*r^2).

Thus the five lines x,y,x-y,x+y,T*x-R*y in N_7 are distinct, and only
x and y occur with multiplicity two.

**Elementary pencil lemma.** Every plane polynomial F_lambda=N-lambda*D
is squarefree, for every lambda in C. Indeed, if an irreducible h occurs
twice, its top homogeneous part squared divides N_7. Up to scalar, that
top part can only be x, y, or xy. The first two possibilities give
h=x+c or h=y+c; these cannot divide a polynomial whose respective
leading coefficient in the OTHER variable is the constant -1. In the
third possibility h=xy+alpha*x+beta*y+gamma, so its leading coefficient
in x is y+alpha. This too cannot divide the constant leading coefficient
of F_lambda in x. All three possibilities are impossible.

Suppose now that an irreducible positive-dimensional critical component
exists. In characteristic zero, f is constant, say lambda, on its
reduction: all its differential components vanish there. At a critical
point g=H*z+L=0, and (52) and its x,y derivatives give

    N-lambda*D=0,  N_x-lambda*D_x=0,  N_y-lambda*D_y=0.

These are polynomial identities on the critical locus; they do NOT
require D!=0. The image in the (x,y)-plane lies in the singular locus of
the squarefree plane curve F_lambda=0. That singular locus is finite.
An irreducible component must therefore have CONSTANT x and y.

This argument is useful beyond this example: a quadratic-block potential
with squarefree Schur pencil cannot have a critical component moving in
the remaining plane variables. It does not by itself exclude components
inside a singular quadratic block; the following step is essential.

### Keeping all singular quadratic blocks

At fixed x,y, the equations f_a=f_b=0 are H*z=-L. If H is invertible
there is at most one solution. If H has rank one, its solution set is
empty or a line z=z_0+s*w. A positive-dimensional critical subset would
have to be that whole line. The quadratic coefficients of f_x and f_y
along it then vanish:

    w^tr H_x w=w^tr H_y w=0.

For rank-one symmetric H, adj(H)=c*w*w^tr for some c!=0. Therefore
such a line requires D=D_x=D_y=0. Rank-zero H also has D=D_x=D_y=0.
It suffices to analyze ALL singular points of D in the plane.

Put M=x^2+T*x+1 and V_0=y^2+R*y+1. The determinant identities are

    D=4*M*V_0-(M+V_0-1)^2,
    D_x=2*(V_0-M+1)*(2*x+T),
    D_y=2*(M-V_0+1)*(2*y+R).                          (53)

They classify its singular points completely:

1. The centered point x=-T/2, y=-R/2, when D vanishes there.
2. T=2*epsilon, x=-epsilon, y=0 or -R, epsilon=+1 or -1.
3. R=2*epsilon, y=-epsilon, x=0 or -T.

For example, if 2*x+T=0 but 2*y+R!=0, (53) gives V_0=M+1 and
D=4*M=0. Hence M=0, V_0=1, T^2=4, and y=0 or -R. The other
case is symmetric. If neither derivative of M,V_0 vanishes, the two
remaining factors in (53) cannot both vanish because their sum is 2.
This proves exhaustiveness, including special parameter intersections.

### The centered point: a compact exclusion identity

At a centered singular point on b_0!=0, t^2 and r^2 are both different
from 1. If t^2=1, then D=0 at the center forces R=0, which is excluded;
the other case is symmetric. Thus H_11 is a unit on this chart, and
H has rank one. Solve f_a=0 exactly for

    a=-(H_12*b+L_1)/H_11.

Let ell_x and ell_y be the coefficients of b in f_x and f_y after
this substitution and after setting x=-T/2,y=-R/2. Define

    p=t^2, q=r^2, S=p+q, V=p*q,
    C=16*t^4*r^4*D_at_center,
    E=(V-1)*S-2*V*(V+1),
    B_x=8*t^2*r^3*(p-1)*ell_x,
    B_y=8*t^3*r^2*(q-1)*ell_y.

Expansion of the same full potential gives

    C=-(V-1)^2*S^2-8*V*(V+1)*S+4*V*(V+1)^2,
    B_x-B_y=-3*(p-q)*E,                              (54)
    (V-1)*C+4*V*(V+1)^4
       =E*(-(V-1)^2*S-2*V*(V+1)*(V+3)).               (55)

A critical line would require ell_x=ell_y=0. Because p-q is a unit,
(54) gives E=0. Also C=0. Identity (55) and V!=0 force V=-1;
then E=0 gives S=0. Thus {p,q}={1,-1}, contradicting b_0!=0.
No centered critical line exists. The identity uses the actual linear
terms of the adjoint character, not just its quadratic coefficient matrix.

The saved certificate also contains an explicit combination of the
original centered constraints equal to r^4. Independent exact replay
checks that identity as an additional algebraic check; the compact
identities (54)-(55) already give the written argument.

### Rank-zero points: two uniform quartics

The other points in (53) have H=0 and L=0. They have t=epsilon or
r=epsilon. Let z denote the remaining nonzero central coordinate and
epsilon=+1 or -1. One of f_x,f_y is

    +(z+z^-1)*(a*b-1)  or  -(z+z^-1)*(a*b-1).

Since z+z^-1!=0 on b_0!=0, a*b=1. Substituting b=a^-1 into the
other equation and multiplying by a^2 gives, respectively,

    P_0(a)=-epsilon*(a^4+1)-(z-z^-1)*(a^3-a),
    P_1(a)=epsilon*(z^-2*a^4+z^2)
              -(z+z^-3)*a^3-(z^3+z^-1)*a.             (56)

These correspond to the zero and negative-trace choices in (53).
Their leading coefficients are -epsilon and epsilon*z^-2, both units.
They therefore have finitely many roots for every allowed parameter.
This covers both signs and both central axes with the same two formulas.

Every putative positive-dimensional critical component has now been
excluded. This proves the theorem about ALL central fibers on b_0!=0.

### Flatness, reducedness, and the actual target consequence

Let S_0=C[t^(+-1),r^(+-1),b_0^-1] and let R_0 be the original relative
critical algebra over S_0. Its presentation has four equations in four
affine variables, and every nonempty fiber now has dimension zero.
The standard [relative complete-intersection theorem](https://stacks.math.columbia.edu/tag/00SW)
therefore gives flatness. The earlier source-generic finite-etale
certificate of rank 38 gives a reduced generic algebra. Flatness makes
R_0 inject into that algebra, proving R_0 reduced. The relative scheme
is Cohen--Macaulay and pure of dimension two. It need not be finite over
S_0; the previously proved length drops remain consistent with flatness.

For the actual central target, retain mu=(t*Phi_t,r*Phi_r), with the
already justified coordinate change modulo affine gradients. The earlier
six whole-boundary exclusions prove b_0 is a unit after localization at
the central TARGET function field L=C(c,d). Hence

    A=R_0 tensor_(C[c,d]) L.

This is a localization of a reduced ring, hence reduced. Moreover every
component of R_0 has dimension two. A component surviving target
localization dominates the two-dimensional target, so its generic
fiber has dimension zero, by the standard dimension formula. Thus A
is a reduced finite-dimensional L-algebra, or zero. In characteristic
zero it is a product of finite separable field extensions, hence
[finite etale](https://stacks.math.columbia.edu/tag/00U3).
The regular-quotient transfer in (39)-(40) now gives generic smoothness
of the actual regular-source Springer fibers for this K, including
regular nonsemisimple representatives.

**Verification and remaining limits.** The
[Sage verifier](tools/verify_adjoint_quasifinite.py) and
[exact certificate](adjoint-quasifinite-verification.json) reconstruct
the full potential and check (51)-(56) and the centered membership.
The [independent SymPy audit](tools/audit_adjoint_quasifinite.py) and
[audit report](adjoint-quasifinite-audit.json) reconstruct everything
over Q(t,r), without modular specialization or Sage. They check the
compact centered identity, original-constraint membership, and all
rank-zero specializations. The argument above supplies the geometry;
the programs verify its explicit algebraic hypotheses and identities.

The generic-target corollary also uses earlier input: all six boundary
dimensions have exact Sage checks. The subsequent three-path continuation
completes the independent audit of all six, using original-ideal
memberships and original-potential involutions as detailed below.
The zero-dimensionality theorem on b_0!=0 itself does not depend on those
boundary calculations. The standard algebraic implications are published
inputs; the full-potential identities and their application are deductions
here. No literature-priority claim is made.

What remains open is the uniform all-simple/all-proper-Levi conjecture,
the generic-target degree and nonemptiness for this candidate, and the
behavior at all special targets. The next structural question is whether
adjoint root/branching data force squarefree Schur pencils and analogous
linear-coefficient exclusion identities in a family of Levi models. Those
are the precise features used here; a new list of isolated fiber counts
would not supply that explanation.

## 25 October: implementing the three structural proof paths

This continuation keeps the conjecture at the start of this note unchanged.
It proves a root-subsystem classification for individual cubic interactions
and an exclusion theorem for points smooth in their relative critical
fiber. Neither is a proof of the uniform adjoint conjecture. The first is
a deduction from the earlier norm bound and standard Dynkin theory; the
second extends the argument of the existing quadratic theorem. No priority
claim or additional Levi census is made.

### Completion of the six-boundary audit

All six source-boundary dimension exclusions used in the E8 generic-target
corollary now have independent characteristic-zero certificates. Four
boundaries, t=r, t=-r, tr=1 and t^2+1=0, have fifteen explicit
original-ideal memberships each, with leading supports proving dimension
at most one. The Q(i) identities hold under both complex embeddings,
so the last case retains both t=i and t=-i.

The other two follow through the following involutions of the FULL
original potential Phi(x,u,v,y;t,r):

    sigma:(x,u,v,y;t,r) -> (-x,u,v,y;-t,r),
    tau:  (x,u,v,y;t,r) -> (y,t*u/r,r*v/t,x;r,t).

Both preserve Phi and square to the identity. Their affine Jacobians
have determinant -1, and their coefficients and inverses are Laurent
regular on the entire central torus. The gradient chain rule therefore
gives isomorphisms of the ACTUAL relative critical schemes. Sigma maps
tr=-1 to tr=1; tau maps r^2+1=0 to t^2+1=0. No source parameter, affine
divisor, nilpotent, or intersection is discarded. A scheme of dimension
at most one cannot dominate the two-dimensional central target, whether
or not these isomorphisms are used to identify its target coordinates.

The [Sage certificate report](adjoint-boundary-exclusions-verification.json)
and [independent native audit](adjoint-boundary-exclusions-audit.json)
verify the 60 memberships, the leading-support bounds, both involutions,
and the full gradient chain rule. Thus the previously stated independent
audit gap in the E8 corollary is closed. This does not calculate its
generic degree or nonemptiness, or prove the uniform conjecture.

### Path 1: the root subsystem supporting each cubic interaction

**Theorem (cubic carriers).** Let G be simple, simply connected and simply
laced, and let K be a proper Levi with derived factors of type A. Every
external cubic monomial in the FULL adjoint character couples exactly
three different factors. Up to permutation and opposition in a factor,
their types and fundamental weights are precisely among

| Coupled factors | Fundamental weights | Generated root subsystem |
| --- | --- | --- |
| A1, A1, A_m, m>=1 | endpoint weights on all three | D_(m+3) |
| A1, A2, A_m, 2<=m<=4 | endpoint weights on all three | E_(m+4) |

This is a necessary classification of interactions that actually occur,
not an assertion that each row embeds in every G. To construct the
subsystem, choose a root alpha highest for the Levi constituent producing
the cubic term. Adjoin -alpha to the simple roots of the three coupled
factors. These roots are a simple system for the displayed subsystem.
Every root generated by them is an ambient root and occurs with
multiplicity ONE in the adjoint character.

Other Levi factors are orthogonal to this subsystem. Contributions from
ambient roots OUTSIDE it may still couple its character variables to
other variables and parameters. Thus the statement does not decompose
the full potential into independent summands on disjoint variables.

**Proof.** The earlier adjoint norm theorem gives three fundamental
weights on distinct A factors and a nonzero central projection gamma,
with

    ||omega_i||^2+||omega_j||^2+||omega_k||^2 < 2.

For A_m, ||omega_i||^2=i(m+1-i)/(m+1). At least one factor is A1,
since three factors of rank at least two have sum at least 2. With one
A1 fixed, an interior fundamental weight on either remaining factor
has squared norm at least 1, while the third factor contributes at least
1/2; this contradicts the strict inequality. Hence all three weights
are endpoints. Put the factor ranks in increasing order, 1<=b<=c.
The norm inequality is

    1/(b+1)+1/(c+1) > 1/2.

If b=1 there is no further bound on c. If b>=2, then b=2 and c<=4.
This proves the two rows and works in arbitrary rank.

The Levi-highest root alpha pairs with each simple coroot as zero,
except for a one at the chosen endpoint of each coupled factor. Thus
-alpha has inner product -1 with those three endpoints and zero with
the other simple roots. It is independent of the Levi span because
gamma!=0. The resulting positive definite Gram matrix is the star with
arm lengths (1,1,m), or (1,2,m). Its Dynkin type is D_(m+3), or
E_(m+4), respectively. Its determinant can also be read without diagram
identification by a Schur complement:

    4(m+1)*(2-1-m/(m+1)) = 4,
    6(m+1)*(2-1/2-2/3-m/(m+1)) = 5-m.

Reflections in these roots preserve the ambient root system, so their
orbits give the claimed subsystem. The adjoint representation has every
ambient root with multiplicity one. Standard finite-root-system theory,
not a computation of fiber schemes, completes the proof. The classification
input can be compared with [Fomin--Reading, lectures 1--2](https://arxiv.org/abs/math/0505518).

**What this does and does not supply for Schur elimination.** It determines
the complete root-string contribution associated to one cubic term,
including the lower-degree terms of that subsystem. It explains why
retaining only the highest cubic monomial loses necessary data. But two
such subsystems can overlap; their roots must be counted once in the
ambient character, not once per subsystem. The Hessian, its adjugate,
and the resulting Schur numerator involve all these terms together.
Neither squarefreeness of N-lambda*D nor the rank-one exclusion identity
in (55) follows from the star diagram argument. That is the exact
unproved implication in Path 1. No root-defined extension of the E8
quasifiniteness theorem is claimed here.

In fact, demanding squarefree pencils uniformly would be too strong.
The established E6 contact character is

    Phi=a*b+(z+z^-1)*x+z^2+z^-2,

with two further free fundamental characters w1,w2. Taking the quadratic
block (a,b) gives D=-1 and

    N=-((z+z^-1)*x+z^2+z^-2).

At z=i or -i, Phi is -2 on the relative critical fiber, and
N-(-2)D is the ZERO polynomial, including if a free coordinate is used
as a second plane variable. These source strata support the entire
nonempty generic target fiber. This is an actual adjoint control, not
a defect removable by declaring z generic. It rules out a universal
squarefree-pencil lemma of this form; it does not rule out its use
on suitably justified finite-critical-fiber charts.

### Paths 2 and 3: smooth relative critical fibers cannot cause a dominant bad image

Here is a statement valid in every affine dimension and for arbitrary
Laurent-polynomial potentials, not only adjoint characters. It uses the
ACTUAL relative critical scheme, including embedded or nonreduced structure.

Let Phi(x,z) be regular on A^ell times T, T=(C*)^s, and set

    X=V(Phi_x1,...,Phi_xell), pi:X->T,
    mu_j=D_j Phi, D_j=z_j partial_z_j.

Define the constructible subset

    M={p in X : the actual fiber X_(pi(p)) is smooth at p}.

Smoothness here is local; the rest of that same fiber need not be smooth.
It does NOT mean that pi is smooth or flat at p.

**Theorem (smooth-relative-fiber exclusion).**

    dim closure(mu(Sing(X) intersect M)) < s.             (57)

Consequently, after deleting a proper closed target subset, every
nonsmooth point of an actual target fiber lies in

    E=Sing(X) intersect (X minus M).                      (58)

In particular, if every relative critical fiber is smooth at all its
points, then generic target fibers are smooth, possibly empty and
possibly of positive or varying dimension. More generally it suffices
to prove that mu(E) is not dense. This is a reduction of the remaining
obstruction, not a claim that E has nondense image for every adjoint
character.

**Proof.** Make the tangent and local-dimension comparisons at closed
complex points; constructibility then gives the stated image bounds.
Stratify the central torus into finitely many smooth locally
closed strata S such that X_S->S is flat. This is standard
[generic-flatness stratification](https://stacks.math.columbia.edu/tag/0H3Y).
At a point of M over S its fiber is smooth, so X_S->S is smooth locally
there, by the [flatness and smooth-fiber criterion](https://stacks.math.columbia.edu/tag/01V8).
Refine the strata in X further so that the following ranks and relative
dimension t are constant. Write d=dim S and

    H=(Phi_xi_xj), B=(D_j Phi_xi),
    k=rank H, p=rank[H | B]-k.

Smoothness of the actual relative fiber gives t=ell-k. The ambient
incidence tangent dimension is

    dim T_p X=ell+s-rank[H | B]=t+s-p.

The locally closed smooth family X_S has dimension d+t at this point,
so dim_p X>=d+t. If X is singular at this point, its tangent dimension
is strictly larger than its local Krull dimension. Therefore

    d+p < s.                                             (59)

On a fixed central fiber the differential of mu is B^tr on ker H.
Symmetry of H gives the exact linear-algebra identity

    rank(B^tr restricted to ker H)=rank[H | B]-rank H=p.   (60)

Any smooth stratum of the reduced SUPPORT of the bad set in that fiber
has target differential rank at most p. Characteristic zero bounds its
image dimension by p. Allowing the central stratum contributes at most
d, so each bad stratum has target image dimension at most d+p<s.
Equivalently first bound the image under (pi,mu), then project to the
target. Taking the finite union proves (57). Only sets used to bound image
dimensions have been reduced; no fiber ideal has been radicalized.

Finally, generic smoothness on the smooth locus X_sm excludes its
critical values and nondominant components after shrinking the target.
This needs no properness. Together with (57), any remaining nonsmooth
target-fiber point must lie in (58). The case s=0 follows directly:
X is its own relative fiber and Sing(X) intersect M is empty.

**Independent linear-algebra check.** The dual of the map in (60) is
the map from the central tangent space to (ker H)^*. For a symmetric
H, (ker H)^* identifies with the quotient of the affine cotangent
space by im H. Its rank is exactly the rank added by the columns of B
modulo im H. This recovers (60) without assuming a nonsingular Hessian
block, a regular sequence, or a constant-rank family on the whole torus.

**Exact excess check.** For the E6 character above, the relative fibers
are empty when z+z^-1!=0. At z=+/-i they are A^3 with coordinates
x,w1,w2: a=b=0, H has rank two, and B adds rank one. The actual target
equation is

    (z-z^-1)*x+2(z^2-z^-2)=c.

Thus over L=C(c), its entire generic algebra is

    L[z,w1,w2]/(z^2+1), x=-c*z/2, a=b=0.                 (61)

The quadratic z^2+1 is squarefree and z is a unit in this algebra.
Its two geometric sheets are smooth of dimension two. All six original
equations (five affine derivatives and the target derivative) give this
same ideal; their number has not been used to infer their height.
The rank argument handles the free directions on their actual source
stratum, instead of inverting its vanishing coefficient.

There is no general requirement that smooth excess be an affine-space
factor. For example Phi=(uv-1)^2+z has critical scheme the disjoint union
of the origin and uv=1, times the central torus. Generic target fibers
retain a point and Gm, both smooth. This is an abstract control, not a
new representation example or a claim of global product structure for G.

The theorem also permits a nonreduced TOTAL incidence with smooth
relative fibers. For Phi=(z-1)^2*x+z, the incidence ideal is ((z-1)^2).
Its fiber at z=1 is a smooth affine line, but the total incidence is
nonreduced. Here mu-1=(z-1)*(1+2*z*x); its square, but not itself,
vanishes in the incidence ring. The bad target support is c=1 and
generic target localization kills the entire thickening. This checks
that (57) does not silently replace the incidence by its reduction.

### Path 3 attempt: higher invariant flat differential equations do not distinguish the obstruction

The earlier failure of the flat Laplacian extends to the ENTIRE algebra
of Weyl-invariant constant-coefficient logarithmic differential operators.
This closes a possible next route; it does not concern the radial group
Casimir or operators with nonconstant coefficients.

Indeed, for a simply laced irreducible root system R and arbitrary root
multiplicities m_alpha, put F=sum m_alpha*e^alpha. For a polynomial p on
the weight space, its flat operator satisfies

    p(partial_log) e^alpha = p(alpha)e^alpha.

If p is W-invariant, transitivity of W on R makes p(alpha) a constant
lambda_p on the entire root set. Hence EVERY choice of m_alpha satisfies

    p(partial_log) F=lambda_p F.                         (63)

Equivalently, for Phi=r+F,

    (p(partial_log)-lambda_p)Phi=(p(0)-lambda_p)r.

For a nonsimply laced root system the same argument applies separately
to its long-root and short-root orbits. Such operators can distinguish
the two orbits, but cannot enforce equal multiplicities within either.

The true adjoint character and the altered 82-dimensional E6-Levi
character above have the same zero-weight multiplicity and differ only
in multiplicities of four ambient roots. Both therefore satisfy (63)
for every invariant p, while the altered character has a dominant
nonreduced target fiber. Thus even ALL these flat differential identities,
together with the already shared nonlinear coefficient data, cannot
prove the desired bound. Finite root-moment checks corroborate the
diagonal-operator calculation; the all-operator assertion follows from
the displayed proof, not from checking finitely many degrees.

Full ambient Weyl invariance of the CHARACTER itself is a different and
stronger condition. It forces equal root multiplicities and excludes the
altered control. It has not yet yielded an identity on the singular
incidence scheme. In particular, passing from the Laplacian to higher
flat invariant operators would not supply the missing input in (62).

### What remains after the three paths

The remaining adjoint-specific task is to show

    closure(mu(E)) != A^s, with E as in (58).             (62)

This focuses the problem on incidence points which are ALSO singular
in their actual relative critical fiber. It removes smooth excess from
the unresolved obstruction. It does not remove nonreduced isolated
relative fibers or singular positive-dimensional relative fibers.

Even reducedness of every relative fiber is insufficient. The abstract
potential Phi=p*q*r+z has relative gradient ideal (qr,pr,pq), the radical
ideal of three coordinate axes. The incidence and every relative fiber
are reduced, but their triple intersection is singular and maps by
mu=z onto every nonzero target. Thus (62) needs more than reducedness.
The altered honest Levi character remains a different control: its
local ideal (u^2*v,u*v^2) has the nonzero nilpotent uv, and is not in M
along the crossing. Neither control is an adjoint counterexample.

The known Sp4 special-target obstruction also survives the theorem:
its relative affine critical fibers are smooth, but its singular
incidence points have target c satisfying c^2+12=0. This is precisely
a proper bad image allowed by (57); all-special-target smoothness has
not been asserted.

**Status and next proof gate.** The root-carrier theorem supplies a
uniform description of individual cubic interactions. The exclusion
theorem supplies a uniform treatment of smooth relative excess. The
missing step is a full-adjoint identity controlling (62) when several
carriers interact, or an actual adjoint counterexample. No such identity
has been proved by this continuation. The next attempt should analyze
the gradient ideal at the generic point of a potentially dominant
component of E, retaining the ambient root multiplicities. A valid
advance is a target annihilator or rank defect there; another finite
case or a further reformulation is not a substitute.

The [native exact audit](adjoint-proof-paths-audit.json) reconstructs
carrier root subsystems from the existing character records, checks the
arbitrary-rank determinant identities, and verifies the control ideals
in characteristic zero. The separate Sage verifier checks the same
carrier statement directly from root-lattice coordinates and replays
the control ideals; its [report is recorded separately](adjoint-proof-paths-verification.json).
These finite checks corroborate the written proofs;
they are not an independent human review or evidence of novelty.

## 25 October: overlapping subsystems and a coefficient-alignment mechanism

**Outcome.** The overlap argument does NOT prove (62) uniformly. It does
explain one genuine full-coefficient cancellation and proves the
arbitrary-dimensional normal-form theorem below (GEO-78). We have not
shown that the unresolved adjoint germs admit that form. In particular,
this is not a new E6/E7/E8 coverage claim or a root-theoretic family
classification. Two proposed shortcuts fail for the reasons below.

### Root-string compatibility gives syzygies, not extra critical equations

Normalize the inner product of a simply laced root system by
(alpha,alpha)=2. Write E_alpha=e^alpha in the torus character ring and
F=sum_(beta in R) E_beta, so chi_Ad=rank(G)+F. Root strings give

    partial_alpha F=(1-E_(-alpha))*Q_alpha,
    Q_alpha=2*(E_alpha+1)+sum_((beta,alpha)=1) E_beta.    (64)

Here partial_alpha is the logarithmic derivative in the coroot
direction. Indeed the terms beta and beta-alpha pair the inner-product
1 and -1 roots; the two roots +/-alpha contribute
2*(E_alpha-E_(-alpha)). This proves (64) in the actual group algebra,
without evaluating at points or taking radicals. It uses the correct
adjoint root multiplicities. The root-string input is classical; see
[Fomin--Reading, Lectures 1--2](https://arxiv.org/abs/math/0505518).

If (alpha,beta)=-1, so alpha+beta is a root, the overlap identity is

    (1-E_(-alpha-beta))*Q_(alpha+beta)
      =(1-E_(-alpha))*Q_alpha+(1-E_(-beta))*Q_beta.       (65)

This is exactly partial_(alpha+beta)=partial_alpha+partial_beta.
Consequently, wherever the relevant root denominators are units,
adding the divided root equations for the rank-two overlaps does not
enlarge the ideal generated by the original directional derivatives.
The same statement holds for prescribed target derivatives, because
the target pairing is linear in the coroot. Differentiating (65) gives
linear combinations of the existing Jacobian equations modulo that
ideal; it supplies no independent row-rank bound there.

This is a limitation of THIS use of overlap identities, not a claim
that full Weyl invariance is useless. Nor may (64) replace the ordinary
fundamental-character gradient on the Levi discriminant. There the
torus quotient is ramified and regular nonsemisimple sources must still
be retained. Formula (65) cannot justify deleting that locus.

### A full critical point need not be critical for its overlapping subsystems

There is a small exact countercontrol to this proposed gluing step
inside the actual SL4 adjoint map, with a proper target centralizer.
Let r satisfy 10r^2+12r+5=0 and put

    g=diag(1,2,r,1/(2r)),
    y=g*tr(g^-1)-g^-1*tr(g),  theta_Ad(g)=y/8.           (66)

All four eigenvalues of g are distinct; y_1=y_2=-27/10, while y_3,y_4
are distinct from this value and from one another. Thus g is regular
semisimple and its target has precisely the proper Levi centralizer
S(GL2 x GL1 x GL1). This is an allowed centralizer-exact target.

For a subset I of the diagonal indices let chi_I be its root-subsystem
adjoint character, with constant |I|-1. Counting each root once gives

    chi_1234=chi_123+chi_124-chi_12+chi_34-1.            (67)

In the retained direction H=diag(1,-1,0,0), the full derivative is zero.
The two A2-subsystem derivatives are respectively

    -3/5+(3/2)r,  -3/5+(3/2)/(2r),

both nonzero; their sum is -3, exactly the A1 overlap derivative.
The orthogonal A1 term has zero derivative in H. Thus cancellation
between the correctly counted overlapping subsystems occurs at an
actual incidence point. Individual subsystem critical-locus theorems
cannot be applied to this point by restriction alone.

This tests a proposed lemma, not the main conjecture: the example
asserts neither a singular fiber nor new type-A coverage. It prevents
an invalid inference that would already fail in a known positive family.

### A coefficient identity that really does control the bad image

Return to the EXISTING E6 control, rather than introducing another Levi
test. With S=z^3+z^-3, A=z^-2+zx, B=z^2+x/z, its full character is

    Phi_k=ab+uv+x^2+k*S*x+Aav+Bbu-2,

where k=1 is the ambient adjoint and k=2 is the altered honest Levi
character. The exact invertible change a=a'-Bu, b=b'-Av gives

    Phi_k=a'b'+(1-uv)*x*(x+S)+(k-1)*S*x-2.             (68)

The important identity is q=-det(M)-2 for the scalar term
q=x^2+S*x-2 and M=[[1,A],[B,1]], since det(M)=-x(x+S).
It depends on the linear coefficient, not just the unchanged matrix M.

Here the root overlap DOES distinguish the coefficients. The four
roots of the added S*x module have multiplicity k; the other ambient
roots have multiplicity one. Reflection in the omitted simple node 4
sends one of these four roots to a root outside that set. Ambient Weyl
invariance forces k=1. Thus this reflection restores (68)'s
factorization; the flat invariant differential equations did not.
The native weight and independent Sage root-lattice checks both verify
the reflection across the two sets. This explains the already known
E6 coefficient correction without assuming a new root-family theorem.

The change depends on z, but its central derivative changes by a
combination of affine derivatives. It therefore preserves the actual
target equations on the incidence scheme, not just their zero sets.
Eliminating the nonsingular a'b' block and centering x gives

    f=(1-uv)*(X^2-delta(z))+Q(z),
    delta=S^2/4, Q=-2.                                 (69)

The following general deduction explains why this form works.

### Theorem: the singular image for an aligned quadratic coupling

Let Z=(C*)^s, s>=1, m>=1, and delta,Q be Laurent polynomials on Z,
with delta NOT the zero polynomial. Define

    w=1-sum_(i=1)^m u_i*v_i, P=x^2-delta,
    f=w*P+Q,
    X=V(f_x,f_u1,...,f_um,f_v1,...,f_vm),
    mu=(D_1 f,...,D_s f), D_j=z_j partial_zj.

All these are actual schemes. Let B=V(delta) in Z and write DQ and
Ddelta for the vectors of logarithmic derivatives. On underlying
constructible sets (without replacing the incidence ideal by its radical),

    mu(Sing X) subset DQ(B) union (DQ-Ddelta)(B),        (70)
    dim closure(mu(Sing X)) <= s-1.

Hence the actual generic target fiber of mu is smooth, possibly empty
and possibly of positive dimension. Adding free affine variables and
an independent nonsingular quadratic block preserves this assertion.
No squarefreeness of delta or regular-sequence hypothesis is required.

**Proof, including the boundary.** Up to nonzero constants the exact
incidence ideal is

    J=(xw, P*u_i, P*v_i),                              (71)
    mu=DQ-w*Ddelta.

Where delta!=0, either all u_i,v_i vanish, in which case x=0 and
the local ideal is (x,u_i,v_i), or some u_i or v_i is a unit. In the
latter case P=0, x is a unit and w=0; both equations eliminate smooth
coordinates. Thus X is smooth away from B.

Consider a point over B. If all u_i,v_i vanish, w is a unit and x=0;
its target is DQ-Ddelta, giving the second set in (70). We need not
declare this point singular to obtain the inclusion.

Otherwise some u_i or v_i is a unit. In that chart (71) is EXACTLY
(P,xw). Also w is a valid smooth affine coordinate, since its
differential in the u,v variables is nonzero. If x!=0 this chart is
smooth as above. At x=0, if w=0 its target is DQ, in the first set.
If w!=0, the local ideal is (x,delta). This is smooth whenever
d delta!=0, and otherwise its target is again DQ because Ddelta=0.
This accounts for every singular point, including thickened components
when delta has repeated factors. The ideal equalities used no radicals.

Because delta is a nonzero Laurent polynomial on an irreducible torus,
B has dimension at most s-1. Each displayed map B->A^s therefore has
image closure of dimension at most s-1, proving (70). Remove the
closures of these two target images, then apply characteristic-zero generic
smoothness to the remaining smooth source components. Nondominant
components can also be removed by their target-image closures.
This proves the generic assertion. The final input is the classical
[generic-smoothness theorem in the target, Vakil, Section 3.3](https://math.stanford.edu/~vakil/0506-216/216class5152.pdf).
The chart and image-bound argument above is our written deduction;
no novelty claim has been checked.

**Separate geometric check on the nonzero-vector chart.** The two
equations P=x^2-delta and xw have height two: P is a nonzero divisor,
and xw is a nonzero divisor modulo P, since delta is nonzero and P
has no common factor with x or w. This proves the regular-sequence
assertion on this chart rather than assuming it. Its Jacobian in
x,w and logarithmic central coordinates has rows

    (2x, 0, -Ddelta),  (w, x, 0).

The two-by-two minors include 2x^2 and w*D_j delta. Thus in the
actual Jacobian singular subscheme, delta=0 and mu-DQ=-w*Ddelta=0.
The target map on that subscheme factors through DQ on V(delta),
including its nilpotents. On supports this independently recovers
the first set of (70); the zero-vector locus gives the second set.
This calculation also confirms that no additional singular chart was
lost in the preceding case argument.

The delta!=0 hypothesis matters. If delta is identically zero and
Q=z for s=m=1, then near u=v=1 the incidence ideal is (x^2,xw).
The class of x is a nonzero nilpotent at w=0 and mu=z. This bad
component dominates the target. Thus the theorem is not a disguised
claim about arbitrary coupled potentials.

For (69), the support of B=V(S^2) is V(S), and both Ddelta and DQ
vanish on that support. The singular-incidence image is therefore
supported at {0}. This excludes its entire scheme over the generic
target without asserting that its special target image is reduced.
It recovers a known E6 control by the coefficient identity, but does not by itself
assert that every nonzero fiber is smooth: critical values of mu on
X_sm must also be excluded, as in the theorem's generic assertion.
The altered k=2 character has the extra S*x term in (68) and is not
covered. Its embedded nilpotent remains a valid negative control.
The separate E6 contact excess and Sp4 special-target checks remain
unchanged; neither is discarded by a source localization here.

### Exact remaining obstruction after the overlap attempt

We have NOT derived (68), or a substitute for (70), for general
overlapping D/E carriers. Correct multiplicities and all the first-order
root-string compatibilities give (64)-(65), but those compatibilities
are syzygies of the existing equations. Nor does a full critical point
lie in the critical locus of each constituent subsystem, by (66)-(67).
The missing implication is a constraint on the ENTIRE coupled scalar
term and Hessian block, strong enough to put the singular target image
inside lower-dimensional central images, while preserving the actual
ideal on singular-matrix strata. We have no such implication for the
remaining adjoint restrictions and no actual adjoint counterexample.

This continuation establishes a mechanism and its precise limits, not
the uniform conjecture. A useful next proof test is to derive a
coordinate-independent replacement for q=-det(M)+constant from the
full root character on the residual Hessian kernel. It must account
for all cross terms and imply a target-image bound without separately
critical subsystem assumptions. If it cannot be derived, the overlap
route should be reported as stalled, not expanded into more Levi tests.

The [native audit](adjoint-proof-paths-audit.json) and
[Sage verification](adjoint-proof-paths-verification.json) now check
(64)-(65) in independent root coordinates, (66) by root sums versus
matrix traces over an exact number field, (68)-(69) symbolically, and
the original chart ideals for m=1,2 with separate Groebner engines.
The arbitrary-m geometry is established by the written proof, not by
those two finite checks. No new API, PDF, census or coverage is added.

## 25 October: the scalar-block identity and its actual limits

**Outcome.** A determinant-only extension of the preceding E6 identity
fails in the already studied E8 quadratic block. There IS a uniform
differential identity from the full adjoint root character, obtained
from the classical radial Casimir identity. It transports through Schur
elimination as described below. This does not prove a new singular-image
bound, and we have not derived an algebraic coefficient relation that
would put every residual adjoint potential in GEO-78's normal form.
The distinction is essential: a general identity is not itself a proof
of generic smoothness. No new Levi calculation or coverage is claimed.

### The determinant alone does not determine the residual scalar

For a quadratic block write

    Phi(q,y,z)=(1/2)q^tr H(y,z)q+L(y,z)^tr q+Q(y,z),
    D=det H, N=QD-(1/2)L^tr adj(H)L,
    S=N/D=Q-(1/2)L^tr H^(-1)L.                         (72)

The completed scalar S, unlike the uncompleted Q, is invariant under
invertible affine changes in q over the residual base. The determinant
changes by det(A)^2 under q=Aq'+b. Thus even a proposed determinant
relation must specify its block and basis convention.

The existing E8 model already rules out a relation S=F(D,z), even an
algebraic dependence between S and D, for its recorded block q=(a,b)
over C(t,r). Reconstruct H,L,Q from the full stored potential rather
than trusting a precomputed numerator. At

    t=2, r=3, x=1, y=0,

the exact calculation gives

    D=-9/4, N=27/4,
    (S_x D_y-S_y D_x)=515/6,
    (Q_x D_y-Q_y D_x)=-525/2.                            (73)

All denominators here are units. In characteristic zero the nonzero
two-form dS wedge dD proves algebraic independence over C(t,r).
Hence neither the original Q nor the completed S can be recovered
from D and the central parameters alone. This is stronger than merely
failing to fit S=-D+constant in this presentation.

There is also a genuine pole of S along the determinant divisor. At
t=2,r=3,x=-1/2,y=0, one has

    D=0, D_x=6, N=12.                                  (74)

The quadratic system Hq+L=0 is inconsistent at this point: H has rank
one and [H|L] has rank two. This is NOT a singular-fiber counterexample.
It shows why a globally regular completion cannot be assumed. The pole
of S is unchanged by regular invertible affine changes of this block.
Neither (73) nor (74) excludes a different block, mixing residual and
quadratic coordinates, or an identity involving other matrix invariants.
They invalidate the proposed determinant-only implication, not the
uniform adjoint conjecture or the earlier completed E8 proof.

### The identity supplied by the full root character

Let Delta be the roots of any simple complex G, r its rank, and
normalize the auxiliary root inner product so long roots have square
length two. This fixes the following differential operator; it does
not change the Killing-form normalization of the Springer map.
For e^lambda in the torus group algebra define

    partial_alpha(e^lambda)=(alpha,lambda)e^lambda,
    Lap(e^lambda)=(lambda,lambda)e^lambda,
    V=sum_(alpha>0) [(1+e^(-alpha))/(1-e^(-alpha))]
                      *partial_alpha,
    Phi=r+sum_(beta in Delta)e^beta,
    Psi=sum_(beta in Delta)(beta,beta)e^beta.

The uniform identity is

    V(Phi)=2*h^vee*Phi-Psi.                             (75)

For simply laced G it becomes

    V(Phi)=kappa*Phi+2r, kappa=2(h^vee-1).              (76)

For B,C,F4,G2 the root-length term Psi must be retained; replacing
it by 2(Phi-r) would be incorrect.

**Published input.** This is the radial Casimir character equation,
not a new representation-theoretic theorem. The Weyl denominator
A_rho satisfies Lap(A_rho)=(rho,rho)A_rho, and the character formula
gives Lap(A_rho*chi_lambda)=(lambda+rho,lambda+rho)A_rho*chi_lambda.
The product rule yields (Lap+V)chi_lambda=(lambda,lambda+2rho)chi_lambda.
For the adjoint highest root this eigenvalue is 2h^vee, and
Lap(Phi)=Psi. This proves (75). See
[Samelson, Section 3.7, Lemma E and the following character-formula proof](https://pi.math.cornell.edu/~hatcher/Other/Samelson-LieAlg.pdf),
printed pages 121--122. We are rewriting a classical input for the
present critical-scheme problem, without a novelty claim.

There is a direct root-string check, independent of character inversion.
Pair beta with s_alpha(beta), put n=2(beta,alpha)/(alpha,alpha)>0,
and use

    (e^beta-e^(beta-n*alpha))/(1-e^(-alpha))
        =sum_(j=0)^(n-1)e^(beta-j*alpha).               (77)

Apply (77) to partial_alpha Phi and multiply by 1+e^(-alpha).
Summing the resulting Laurent polynomials gives (75) coefficientwise.
In the simply laced case each root has coefficient 2(h^vee-1),
while the constant coefficient is 2r*h^vee. These are the same
coefficients as the right side of (76).

Unlike the flat invariant eigen-equations discussed earlier, (75)
DOES distinguish the altered E6-Levi control. At the exact torus point
e^(alpha_i)=2 for every simple root, its defect in
(Lap+V)Phi-24Phi is

    51958967157800573 / 5821210127040 != 0.             (78)

The actual E6 adjoint defect is zero. All roots have nonzero height
at this point, so every displayed radial denominator is invertible.
The altered character is only Levi-invariant; it is not being asserted
to satisfy an ambient character theorem. Formula (78) establishes that
this necessary full-character constraint is stronger than the previous
flat-operator tests; it does not establish its sufficiency for smoothness.

### Descent to the Levi quotient without losing regular nonsemisimple points

Let x_i be derived-Levi fundamental characters and z_j the central
coordinates on the finite central cover. Put

    E_ext=product_(alpha in Delta minus Delta_K)(1-e^alpha).

It is W_K-invariant. The derivation V preserves the invariant function
field. Its apparent poles from roots in Delta_K cancel on EVERY
W_K-invariant Laurent polynomial: pair its weights under s_alpha and
apply (77), with the weight's integral coroot pairing in place of n.
All other denominators divide powers of E_ext. Consequently V descends
to a regular derivation on the quotient localized at E_ext, and has
an expression

    V=sum_i a_i(x,z)*partial_xi + sum_j b_j(x,z)*D_j.    (79)

This localization is justified on the actual centralizer-exact target
fibers. At any geometric quotient point of such a fiber, lift its
semisimple class to t in T. The ordinary character critical equations
imply theta_Ad(t)=h. If e^alpha(t)=1 for an external root, reflection
in alpha fixes t, so equivariance forces alpha(h)=0. This contradicts
Z_G(h)=K. Thus E_ext has no zero on the fiber and is a unit in its
actual coordinate ring, including any nilpotents. This is a unit
argument, not radicalization of the fiber ideal.

No internal Levi discriminant has been inverted. In particular the
quotient points represented by regular nonsemisimple elements remain
present. We use ordinary x_i-derivatives for the incidence throughout;
torus derivatives have not been substituted for them at ramification.
This descent is an elementary consequence of the root-string
divisibility and the established centralizer reduction.

A direct rank-one check is V(x)=x for SL2, x=t+t^-1. Thus
V=x*partial_x on its character quotient, regular even at x=+/-2.
The apparent root denominator has cancelled; those quotient points,
including their regular unipotent representatives, are not removed.

### What survives quadratic elimination

On D!=0 the original equations Phi_q=0 are exactly q=q_*=-H^(-1)L.
For any residual derivative partial_i, including a logarithmic central
derivative, the chain rule gives the exact envelope identity

    partial_i S=partial_i Q+q_*^tr partial_i L
                          +(1/2)q_*^tr(partial_i H)q_*.(80)

Project V to the residual coordinates by evaluating its coefficients
at q_*. Call the resulting derivation V_bar. On the common open where
D and E_ext(q_*) are units, (75) becomes

    V_bar(S)=2h^vee*S-Psi(q_*,y,z).                    (81)

For simply laced G, this reads

    V_bar(S)=kappa*S+2r,
    D*V_bar(N)-N*V_bar(D)=kappa*N*D+2r*D^2.            (82)

These formulas are a general differential relation between the scalar
Schur potential and its quadratic data. The coefficients of V_bar
involve the full root system and the substitution -H^(-1)L; they are
not determined by det H alone. The classical character identity is
the representation-specific input; (80) and its transport are ordinary
linear algebra and the chain rule. If the indicated common open is
empty, (81) supplies no information there.

Equations (81)-(82) do not extend the elimination q=q_* across D=0.
Clearing their denominators gives identities, but can introduce powers
of D and remove exactly the information needed on the singular block.
The actual equations Hq+L=0, with all residual derivatives, must remain
on that divisor. The universal bordered-determinant identity for N is
also valid for every quadratic potential; it is not an additional
adjoint constraint.

### Exact gap after using the radial identity on the full incidence

Work directly with the full incidence X=V(Phi_x), rather than deleting
the D=0 locus. On the justified E_ext-open, (75) and (79) give, in
the ACTUAL coordinate ring of X,

    b^tr mu=2h^vee*Phi-Psi.                             (83)

In simply laced type, differentiating this identity in Omega_X and
using dPhi=sum mu_j*dlog z_j gives

    b^tr dmu=mu^tr*(kappa*dlog z-db).                   (84)

Thus b is NOT automatically an annihilator of the target differential.
The right side of (84) need not vanish, and b itself could vanish.
For nonsimply laced G the corresponding right side is
2h^vee*mu^tr*dlog z-dPsi-mu^tr db. No cancellation of these terms
on every remaining singular stratum has been proved.

This is the missing implication in the current attempt. To obtain
nondominance on a smooth stratum C of the reduced support of the
remaining bad locus, one needs a nonzero covector annihilating dmu|C,
or a nonzero target polynomial on its image. Equations (75)--(84)
have not supplied either uniformly. Nilpotents remain in the original
incidence; reducing C only for an image-dimension argument cannot
be used to infer that the generic fiber is reduced.

**Next proof gate.** The determinant-only proposal should not be
pursued as a uniform lemma. Investigate (83) on the rank-deficient
incidence Hq+L=0, keeping q rather than inverting H. A successful next
step must control the right side of (84), or produce another genuine
annihilator on the singular stratum. Merely writing another character
identity, or showing one good component, would not close the gap.

The [native audit](adjoint-radial-schur-audit.json) checks (75)
coefficientwise from root strings in representatives of all nine finite
simple families, the altered-control defect, the two E8 obstructions,
and (80) with independent coefficient variables. The
[Sage verification](adjoint-radial-schur-verification.json) instead uses
ambient root-lattice inner products and exact rational torus evaluation,
then reconstructs the E8 block from its original potential and solves
the linear equations directly. The all-type proof is classical and
written above; these finite checks are corroboration. Both reports
explicitly leave the uniform singular-image bound unproved.

## 25 October: applying the radial equation on the singular quadratic block

**Outcome.** The equation can be applied without inverting det(H). Its
first derivatives yield (86)--(88) below, but these do not establish a
rank defect of the central map. An exact abstract cubic control satisfies
the radial equation, has a logarithmic-gradient radial coefficient, and
has a dominant singular image supported on det(H)=0. It is NOT an adjoint
character or a counterexample to the uniform conjecture. It disproves the
implication from those abstract differential hypotheses alone (LIMIT-16).
The remaining problem must use the compatibility of the full adjoint
character with its own ambient root discriminant, beyond the one scalar
equation. No new group coverage is claimed.

### The unreduced block equations and their derivative identity

Use the same finite central cover, characteristic-zero coefficient field,
and Killing-normalized adjoint map as above. On the justified E_ext-open,
write x=(q,y), with q the quadratic variables, and put

    F=Phi_q=Hq+L, E=Phi_y, mu_j=D_j Phi,
    X=Spec O/(F,E),
    V=a_q^tr partial_q+a_y^tr partial_y+b^tr D_z.       (85)

The letter E in (85) is the vector of residual derivatives, not the
external discriminant E_ext. All displayed identities below are identities
in the actual quotient ring; nilpotents have not been removed. In simply
laced type, differentiate V(Phi)=kappa Phi+c, c=2 rank(G), with respect
to q BEFORE taking any quotient. Restriction to F=0 gives

    H a_q+F_y a_y+F_logz b
          +(partial_q a_y)^tr E+(partial_q b)^tr mu=0. (86)

The matrices F_y and F_logz are the derivatives of F in the indicated
directions; in components they are H_i q+L_i. The convention for
partial_q a_y is that its columns are q-derivatives of a_y. For every
kernel vector eta of H, (86) therefore gives

    eta^tr F_y a_y+eta^tr F_logz b
                +E^tr partial_eta a_y
                +mu^tr partial_eta b=0.               (87)

There is no inverse of H. The identity also holds for kernel sections
over a nonreduced base. Contracting a derivative means differentiating
the coefficient functions and then pairing with eta; eta itself is not
differentiated. On X the E term vanishes, but the F_y and F_logz terms
need not. A kernel vector of H alone is therefore not an annihilator of
the full critical-equation Jacobian.

For a simultaneous version set A=Phi_xx and B=(D_j Phi_x)_j. On X,

    A a+B b+(partial_x b)^tr mu=0,
    mu^tr partial_xi b=0
       whenever xi^tr A=0 and xi^tr B=0.               (88)

**Proof.** The unrestricted derivative of the scalar radial identity
contains (partial_x a)^tr Phi_x and kappa Phi_x as well as the three
displayed terms. Both extra terms vanish modulo the actual critical
ideal. Pair with xi to get the second equation. This is the chain rule,
not an independent smoothness theorem. In nonsimply laced type the
right sides in (86) and in the first equation of (88) are respectively
-Psi_q and -Psi_x; the last equation becomes
mu^tr partial_xi b=-partial_xi Psi. The root-length term is not discarded.

Nontrivial relations among the displayed generators can also exist at
smooth positive-excess points. Thus a nonzero xi in (88) does NOT by
itself diagnose a singularity. No height or complete-intersection
assertion has been made. Conversely, reducing the equations to a
determinantal rank stratum before differentiating can lose their normal
derivatives. These are separate pitfalls, both relevant on det(H)=0.

### The root expression for the remaining coefficient

The coefficient b in (88) can be written directly from the external
roots, without Schur elimination. Choose the derived and central
coordinates on the cover in the orthogonal decomposition for the
auxiliary root inner product used in (75). Let C be the constant Gram
matrix of the central logarithmic characters. Then

    b=C D_z log(E_ext),                                (89)

where logarithmic differentiation means D_z(E_ext)/E_ext; no analytic
logarithm is chosen. Indeed, for a root pair,

    d log[(1-e^alpha)(1-e^(-alpha))]
       =[(e^alpha+1)/(e^alpha-1)] d alpha.

Summing over positive roots gives V=grad log(E_all). The internal Levi
factor has zero central derivatives, giving (89). The apparent internal
root poles have already cancelled on the Levi quotient as established
above. E_ext is still a unit on the actual centralizer-exact fibers.

Consequently (88) constrains mu by derivatives in the derived variables
of C D_z log(E_ext). These covectors depend on the source point and on
xi. Orthogonality of a target value to such a varying covector is not
an annihilation of the target differential, nor a polynomial equation
in target coordinates alone. The actual relation remains

    b^tr dmu=mu^tr(kappa dlog z-db)                     (84)

on the simply laced incidence. Neither (87), (88), nor (89) forces the
right-hand side to vanish. The exact compatibility of Phi and E_ext as
root expressions from the SAME ambient simple group has not yet been
used to prove a further rank condition.

### A control disproving sufficiency of the scalar radial argument

Here is a written algebraic obstruction valid for every nonzero constant
kappa and every constant c. It is an abstract Laurent potential, with
no assertion that it is a character. On

    U=Spec C[z,z^(-1),(z^2-1)^(-1),(z^2+1)^(-1)],
    s=z+z^(-1), u=z-z^(-1),

take derived variables q,y and define

    f=s(1+y q^2)-c/kappa,
    V=b D_z, b=kappa*s/u.                              (90)

The functions s,u are units on U. Since D_z s=u and D_z u=s,

    V(f)=kappa f+c.

The potential and the vector field are invariant under z -> z^(-1).
When kappa is a positive even integer (as in the simply laced adjoint
constants), b=D_z log(u^kappa), and u^kappa is an inversion-invariant
Laurent unit on U. Thus even this logarithmic-gradient property and
rank-one inversion symmetry do not repair the abstract implication.
This function is NOT being identified with the discriminant of an
ambient group giving f as its adjoint character. That compatibility is
precisely absent. The shift -c/kappa is likewise not asserted to be a
character constant term.

The actual critical ideal and central map are

    I=(f_q,f_y)=(q y,q^2),
    I=(q) intersect (y,q^2),
    mu=D_z f=u(1+y q^2)=u in O_X.                      (91)

The nonzero class of q squares to zero and is annihilated by y. It is
supported on C=V(q,y), an embedded component of the incidence. Away from
C, y is a unit and the scheme is the smooth q=0 locus. At C the local
dimension is two and the tangent dimension is three, so its points are
singular in the actual incidence. At fixed z the local dimension is one
and tangent dimension two, so they also lie in the remaining locus of
GEO-77: singular in their actual fixed-central critical fibers.

The quadratic block is H=2s y, L=0; it vanishes on C. Nevertheless

    mu|C=u,
    D_z(mu|C)=s != 0.                                 (92)

Its image contains the dense target open h(h^2+4)!=0: solve
z^2-hz-1=0. Its two roots are nonzero, have u=h!=0, and satisfy
s^2=h^2+4!=0. In particular the allowed localization has not erased
the persistent bad component. Over C(h), the actual fiber algebra is

    C(h)[q,y,z]/(q^2,qy,z^2-hz-1),                     (93)

with nonzero nilpotent q. The omitted source units are automatically
invertible in (93). Replacing I by its radical (q) would incorrectly
erase all the singularities under discussion. Smooth schemes over C
are reduced; see [Stacks, Section 33.25](https://stacks.math.columbia.edu/tag/04QM).
Here the stronger failure is visible directly from the local dimensions.

On C the differentiated radial right-hand side is explicitly

    mu*(kappa dlog z-db)
       =kappa*s^2/u * dlog z != 0.                    (94)

The q,y derivatives of b are zero, so the conormal constraint (88) holds
vacuously. Both dmu and b are nonzero. This proves that the hoped-for
cancellation is not a formal consequence of the radial identity, its
first derivatives, or a logarithmic-gradient coefficient.

Nor can one substitute y=0 before forming the critical scheme: f|y=0
is independent of q, but the normal equation f_y|y=0=s q^2 remains
essential. If a block of positive rank is desired, add sum p_i^2/2 and
add (kappa/2) sum p_i partial_pi to V. The critical ideal becomes
(p_i,q^2,qy), and H=diag(1,...,1,2s y) has the same remaining obstruction.
This extension is elementary and is not additional group coverage.

### Independent check and stopping point

The [native audit](adjoint-singular-radial-audit.json) checks the
unrestricted derivative identities symbolically, the root-pair
logarithmic derivative, the radial control, its monomial primary
intersection, and the nonzero right-hand side (94). The independent
[Sage/Singular verification](adjoint-singular-radial-verification.json)
imports none of those audit functions: it uses polynomial jets for
the derivative identity and exact ideal operations for (91)--(93).
It verifies that the generic fiber is one-dimensional and nonreduced,
that its bad support has length two over C(h), and that all declared
source localizations preserve that fiber. Both calculations are in
characteristic zero; the geometric proof above is a separate argument.

The three actual program controls remain unchanged: special Sp4 bad
targets require genericity; the E6 example requires allowing smooth
positive excess; and the altered Levi character requires retaining
persistent nilpotents. The new abstract control satisfies the scalar
radial equation, whereas that earlier altered character does not. This
does not make the new control an ambient adjoint representation.

**What remains unproved.** No target annihilating polynomial, uniform
rank defect, or bound on the image of the remaining ACTUAL ADJOINT
singular locus has been obtained. Equations (86)--(89) are exact
deductions, not a proof of that bound. The scalar-radial-only route
should stop here: adding more derivatives of an identity already
satisfied by (90) cannot establish the false abstract implication.

**Next substantive gate.** Use the full root-character/discriminant
compatibility to control the critical-equation conormal module over the
generic TARGET, retaining normal equations and embedded components.
A successful advance must show that this compatibility excludes a
dominant nonsmooth component, for example by producing a nonzero target
polynomial annihilating its support. Studying only the reduced support,
the rank of H, or another scalar radial identity would repeat the failed
argument. This is an unresolved research task, not an asserted theorem
or a promise that this conormal approach will succeed.

## 25 October: the full character and its root discriminant

**Outcome.** There is a uniform scheme-valid relation between the root
discriminant and the ambient character coordinates: the external-root
discriminant is, up to a nonzero constant in logarithmic coordinates,
the square of the character-quotient Jacobian. The quotient map is thus
etale on the entire allowed open, including its internal Levi
discriminant. This is a corollary of Steinberg's Jacobian formula, not
a new smoothness theorem. The full root characteristic polynomial also
couples the adjoint character and discriminant by Newton identities.
Neither statement proves nondominance of the singular critical scheme.
A proposed way to turn the spectral identities into additional critical
equations fails already for the actual SL3 adjoint character over its
generic Levi-central parameter. This is a test of an implication, not a
new Levi coverage calculation or an adjoint counterexample (LIMIT-17).

### A uniform discriminant-Jacobian formula

Let G be simple simply connected over C and K a standard Levi. Write
r=rank(G), ell=rank([K,K]), and s=r-ell. Let

    Q_G=T/W, Q_K=T/W_K, pi:Q_K -> Q_G.

Work on the finite central cover of Q_K with derived fundamental
characters x_1,...,x_ell and central torus coordinates z_1,...,z_s.
Let c_1,...,c_r be the fundamental characters of G pulled back to this
cover. Set

    J=(partial_x c | D_z c),
    E_ext=product_(alpha in Delta minus Delta_K)(1-e^alpha).

**Uniform corollary of the character Jacobian formula.** There is a
nonzero constant gamma, determined by the character-coordinate volume
forms, such that

    det(J)^2=gamma E_ext.                              (95)

Using ordinary dz rather than dlog(z) changes gamma to a Laurent unit.
In particular pi and its covered version are etale wherever E_ext is
invertible. This assertion concerns local etaleness; it does not assert
that restriction to this source open is a finite etale cover of the
whole ambient quotient.

**Proof.** For a root system R put

    A_R=e^(rho_R) product_(alpha in R_+)(1-e^(-alpha)).

Steinberg's character Jacobian formula gives the wedge of fundamental
character differentials as A_R times an invariant torus volume form.
Apply it to G and to the semisimple derived group of K, adjoining the
central logarithmic differentials in the latter case. On the torus
cover the chain rule yields

    det(J)=gamma_0 A_Delta/A_(Delta_K),

where gamma_0 != 0 is a constant lattice-volume factor. This equality
initially computed away from A_(Delta_K)=0 is an equality of rational
functions; the external product on the right is regular on the cover.
The identity A_R^2=(-1)^(|R_+|) product_(alpha in R)(1-e^alpha)
proves (95) globally. Both covered quotients are smooth of dimension r,
so the invertible Jacobian gives etaleness. The input is
[Steinberg, Theorem 6.1 and Lemma 8.2](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf).
This proof is a standard consequence of that input; no novelty is claimed.

We already proved E_ext is a unit in each actual centralizer-exact
target fiber. Equation (95) therefore supplies character coordinates
etale locally at every such quotient point, including those represented
by regular nonsemisimple group elements. It does not invert A_(Delta_K),
take radicals, or assume the critical equations are a regular sequence.

To state exactly what this buys, let P(c) be the polynomial representing
the full adjoint character in the ambient fundamental characters. Then

    Phi=P(c(x,z)),
    (Phi_x,mu)^tr=J^tr grad_c P.                       (96)

Thus the actual generic-target fiber is described, after the justified
localization, by

    J_x^tr grad P=0, J_z^tr grad P=h                   (97)

over C(h_1,...,h_s). The entries of J and their derivatives vary with
the source point. Invertibility of J is a statement about first-order
ambient coordinates, not reducedness of the critical ideal (97).
When differentiating (96), the terms involving derivatives of J remain.
No bound on those terms on the nonsmooth locus has been proved.

This distinction is substantive even when the potential is itself a
local coordinate. For the earlier abstract control (90), the coordinate
map (q,y,z) -> (q,y,f) has logarithmic Jacobian

    D_z f=u(1+yq^2).

It is invertible on an open containing its entire actual incidence:
u is a declared unit and 1+yq^2=1 in O_X. That etale chart coexists
with the persistent nilpotent q. This is only an obstruction to an
etaleness-only inference; the control is not a full root-character model.

### The spectral coupling from the same roots

Introduce a formal variable v, distinct from every source coordinate:

    P_R(v,t)=product_(alpha in R)(1-v e^alpha(t)),
    P_Delta=P_(Delta_K) P_ext,
    P_ext(1,t)=E_ext(t),
    det(1-v Ad(t))=(1-v)^r P_Delta(v,t).               (98)

These are polynomial identities over the character ring, so they
remain valid in any quotient ring with nilpotents. Put
Phi_m(t)=chi_Ad(t^m)=r+sum_alpha e^(m alpha)(t). Then

    -v partial_v log(P_Delta)=sum_(m>=1)(Phi_m-r)v^m,
    -partial_x log(P_Delta)
             =sum_(m>=1)(partial_x Phi_m)v^m/m.       (99)

The logarithms in (99) are formal series at v=0, whose constant term
is one. Evaluation at v=1 is only justified for the corresponding
rational logarithmic derivative when its denominator is a unit; one
cannot substitute v=1 term by term into an infinite formal sum. For
P_ext that denominator is E_ext, already a justified unit. The internal
factor need not be a unit and must not be inverted on the whole fiber.

The identities (98)--(99) make the common-root origin precise. They
are elementary characteristic-polynomial/Newton identities, not an
additional root-theoretic smoothness result. Only the m=1 derived
gradient is set to zero in our incidence. There is no established
reason that the other coefficients of (99) vanish there.

### Failure of spectral criticality propagation in the actual adjoint model

Test the tentative lemma that partial_x Phi_1=0 forces partial_x Phi_m=0
after the external-discriminant localization. This lemma is FALSE, even
on the generic critical incidence of an actual simple adjoint character.

For SL3 take the Levi with derived group SL2 and write

    t=diag(z a,z/a,z^(-2)), x=a+a^(-1),
    S=z^3+z^(-3), R=z^3-z^(-3).

The full adjoint root sum and external-root polynomial give

    Phi=x^2+Sx,
    P_ext(v)=(1-v z^3 x+v^2 z^6)
                   *(1-v z^(-3) x+v^2 z^(-6)),
    E_ext=(x-S)^2.                                   (100)

The ambient fundamental characters are
c_1=zx+z^(-2), c_2=x/z+z^2, so Phi=c_1 c_2-1 and

    det partial_(x,log z)(c_1,c_2)=2(S-x),
    det(J)^2=4 E_ext.                                (101)

For the second Adams character,

    Phi_2=(x^2-2)^2+(S^2-2)(x^2-2).

On the ACTUAL critical incidence 2x+S=0,

    partial_x Phi_2=-(3/2) S R^2,
    mu=D_z Phi=-(3/2) S R,
    E_ext=(9/4)S^2.                                  (102)

All three are nonzero generically. Over C(z), adding partial_x Phi_2=0
to the original critical equation produces the unit ideal. Thus it
deletes the whole generic incidence, not just some irrelevant branch.
The map mu in (102) is nonconstant, so this obstruction also persists
over a dense set of target values. For SL3 and this Levi a nonzero
central target is centralizer-exact. No false special-target assumption
is involved.

For a rational check choose z=2, x=-65/16. Then

    mu=-12285/128,
    E_ext=38025/256,
    partial_x Phi_2=-773955/1024.

The internal Levi discriminant x^2-4 is nonzero at this point as well,
so a regular semisimple source realizes the obstruction. The claim being
refuted is spectral propagation, not smoothness: this critical incidence
is smooth. We have not added a Levi to any coverage count.

### Independent checks, exact gap, and next gate

The [native audit](adjoint-character-discriminant-audit.json) constructs
the SL3 root products directly and uses Newton recurrence for the
spectral coefficients. The independent
[Sage verification](adjoint-character-discriminant-verification.json)
constructs conjugation on eight explicit traceless matrix basis elements
from an SL2 companion block inside SL3. It takes the full and external
matrix determinants and computes power traces and the quotient ideals
over Q(z). No native audit functions or diagonalizing-cover derivative
formulas are imported. Both checks are exact in characteristic zero.
The all-type proof of (95) is the written Steinberg corollary, not an
extrapolation from these SL3 checks.

The existing Sp4 exceptional singular targets, smooth E6 excess, and
altered-character nilpotent control remain unchanged. Nothing here
replaces the actual incidence by its reduced support or its torus
gradient at the internal discriminant.

**Remaining gap.** The character and discriminant are now coupled by
(95)--(99), but these identities do not establish a nonzero target
annihilator for the nonsmooth locus. Etaleness of the ambient coordinate
map does not establish regularity of (97), and imposing criticality of
higher Adams characters is invalid by (102). No uniform singular-image
nondominance theorem has been proved or disproved.

The next necessary scheme-theoretic target is generic-TARGET reducedness
of the actual ideal (97), using its full root-derived coefficients.
This would exclude persistent nilpotent obstructions, including embedded
components. It would still leave the possibility of reduced singular
fibers, which must be addressed separately before claiming smoothness.
Neither generic-source etaleness, reducedness of individual coefficient
divisors, nor a regular-sequence assumption settles this target. A useful
advance must prove a property of that ideal without adding spectral
critical equations or deleting internal-discriminant points. This is a
specific necessary research target, not a newly proved reduction in
complexity or an assertion that the proof method will succeed.

## 25 October: an explicit nilradical certificate in every type A Levi

**Outcome and scope.** The all-simple generic-target reducedness claim
is still open in this program. The R_0/S_1 criterion and the target-torsion
criterion (45) were already recorded above; repeating them would not be
progress. This attempt constructs the annihilator in (45) directly from
the FULL root character for every Levi of SL_n, in every rank. It also
describes the entire nilradical of the relative incidence on its finite
central cover. Type A generic smoothness, hence generic reducedness,
already followed from the quadratic theorem. The contribution here is
the explicit ideal and nilpotent-module calculation, not new smoothness
coverage or a literature-priority claim (GEO-79).

### Statement, with the relative incidence distinguished from its fibers

Let G=SL_n(C), K a proper Levi, and work on the usual finite central
etale cover with derived fundamental characters x and central torus Z.
Use the full adjoint character Phi and the Killing normalization of
the Springer map. Set

    P=C[Z][x], J=(Phi_x), R=P/J,
    N=sqrt(0)_R, mu_j=D_j Phi in R.                   (103)

The D_j form a basis of logarithmic central derivatives; fixed invertible
linear changes identify the mu_j with the central target coordinates.

**Type A nilradical theorem.** For every such K,

    N^2=0, and mu_j N=0 for EVERY j.                  (104)

If the Levi partition has any block of size at least three, or has a
block of size one, or is the torus, then N=0. The only nonzero nilradical
occurs for partitions (2,2,...,2), with at least two blocks. In that case
the explicit presentation below gives N=(L)/(L D), Ann(N)=D, and the
central target functions themselves lie in N. In particular their pairwise
products vanish. These assertions are on the covered incidence; the
annihilator assertion and reducedness descend by faithful flatness.

It follows directly, by localization at nonzero TARGET polynomials, that
the actual generic-target algebra A_K of (39) is reduced, possibly zero.
For a proper Levi the central dimension is positive; any target coordinate
is a nonzero polynomial and becomes a unit in the function field, killing
N by (104). Characteristic zero then gives geometric reducedness.

This does NOT assert that every special nonzero target fiber is reduced.
Localization in the target and specialization at a target are different:
a reduced incidence can have a nonreduced special fiber. Nor does (104)
by itself imply smoothness of the remaining reduced fibers.

### The full root character and the easy partitions

For V=C^n, End(V)=C plus sl_n, so the actual character is

    Phi=A B-1, A=tr_V(g), B=tr_V(g^(-1)).              (105)

Equivalently, on the torus the n^2 terms in A B consist of n constant
terms and all n(n-1) roots, each with multiplicity one. This identity
uses the full adjoint character, including its constant term and all
linear coefficients. It is the standard representation identity, not
a new character formula.

If a Levi block has size d>=3, its first and last fundamental characters
are independent coordinates. Their two critical equations are units
times B and A, respectively. All other critical equations belong to
(A,B), hence J=(A,B). Eliminate those two coordinates, whose coefficients
are Laurent units. The result is a polynomial algebra over the central
torus ring, so R is reduced. Also mu_j=(D_j A)B+A(D_j B)=0 in R.
This does not require any source-coefficient localization. If all blocks
have size one, J=0 and R is already a torus ring.

It remains to treat m blocks of size two and k singleton blocks, where
n=2m+k. Write z_i for the scalar on the ith two-dimensional block,
w_l for each singleton, and x_i for the SL2 trace. Set

    A=sum_i z_i x_i+alpha, B=sum_i z_i^(-1)x_i+beta,
    alpha=sum_l w_l, beta=sum_l w_l^(-1).              (106)

For k>=1 the central torus is
product_i z_i^2 product_l w_l=1. This is connected: one w_l can be
eliminated by a Laurent monomial. For k=0 use its identity component
product_i z_i=1; the resulting central cover still surjects onto K.
Thus the base ring B_0=C[Z] is a Laurent polynomial UFD in either case.
The finite central cover has not removed internal Levi discriminant
points or regular nonsemisimple representatives.

### An exact presentation retaining every coefficient boundary

Put a_i=z_i x_i, c_i=z_i^(-2), and u=A. Retain a_1,...,a_(m-1),
and eliminate a_m=u-alpha-sum_(i<m) a_i. The critical equations,
after multiplying by known Laurent units, are

    F_i=B+c_i A.

For i<m define d_i=c_i-c_m and put

    D=(d_1,...,d_(m-1)) subset B_0[a_1,...,a_(m-1)],
    L=sum_(i<m) d_i a_i+beta-c_m alpha.

Then F_m=2c_m u+L and F_i-F_m=d_i u. Because 2c_m is a unit, this
gives the isomorphism of ACTUAL rings

    R = B_0[a_1,...,a_(m-1)]/(L D).                  (107)

The substitution is u=-L/(2c_m). No d_i, L, or singleton expression
has been inverted. For m=1, only the unit-coefficient equation F_m is
present and R=B_0 is reduced. Below assume m>=2.

For ANY central derivative, computed in the original variables,

    mu=(D_z A)B+A(D_z B)
        =A*(D_z B-c_m D_z A) modulo J.                (108)

The last equality uses F_m=0. There is no change of derivative convention
hidden in (108). After the elimination, every mu belongs to (L) in (107).

### Radicality with singleton blocks

The ideal D is radical of codimension m-1. To see this without a
generic-source argument, set z_i=epsilon_i z_m on D, epsilon_i=+/-1.
For k>=1 each sign choice gives the torus
z_m^(2m) product_l w_l=1, obtained by solving for one singleton. These
components are disjoint, since z_m is a unit. Thus D is reduced and has
exactly the asserted codimension.

On each component of D the scalar g=beta-c_m alpha is NOT identically
zero. Indeed the determinant-compatible curve

    z_m=t^k, z_i=epsilon_i t^k, w_l=t^(-2m)

gives g=k*(t^(2m)-t^(-2k-2m)), a nonzero Laurent polynomial. Hence L
is a nonzerodivisor modulo D: the quotient by D is a product of torus
domains with polynomial variables, and L reduces to this nonzero g
on each factor.

L is also irreducible, hence squarefree, in B_0[a]. For m>=3 the d_i
have no common irreducible factor, because D has codimension at least
two. For m=2 the principal d_1 is squarefree and g vanishes identically
on none of its components, so gcd(d_1,g)=1. In either case L is primitive
over the UFD B_0, and is linear and irreducible over its fraction field.
Gauss's lemma therefore applies. Consequently

    L D = (L) intersect D                             (109)

is an intersection of radical ideals. Equation (109) follows from
(D:L)=D; it is not an assumption of transverse intersection. All
intersections and coefficient zeros remain present. Thus R is reduced
whenever k>=1.

### The actual nilpotents when there are no singleton blocks

Now k=0 and alpha=beta=0, so L=sum d_i a_i belongs to D. On D the
equations z_i=epsilon_i z_m together with product z_i=1 give
z_m^m product_(i<m) epsilon_i=1. All roots are simple in characteristic
zero, so D is radical of codimension m-1, with a finite reduced central
zero scheme. No extension to the source function field is allowed here.

For m>=3, the same primitive-linear argument makes L irreducible. For
m=2, L=d_1 a_1 and d_1 is squarefree; hence L is squarefree in that case
as well. Since L belongs to D,

    sqrt(L D)=(L), N=(L)/(L D),
    N^2=0, Ann(N)=D.                                 (110)

The radical equality follows from V(D) subset V(L) and squarefreeness
of L. For the annihilator, multiplication by L in the ambient polynomial
domain identifies (L)/(L D) with that domain modulo D. This describes
the whole nilradical, including its support at embedded components;
it is not merely a multiplicity-one assertion at minimal primes.

By (108) every mu_j belongs to N. Thus mu_j N=0 and mu_i mu_j=0.
The generic nonzero target fiber in this case is actually empty: an
invertible target coordinate cannot equal a nilpotent. Together with
the preceding partitions this proves (104) for every Levi of SL_n.
Faithfully flat descent from the central cover retains the annihilator
statement. A finite etale base change of a reduced C-algebra is reduced,
so the nonzero nilradical on the all-two-block cover is not an artifact
that could arise by pulling back a reduced incidence.

### Independent checks and relation to the earlier results

The all-rank proof above uses exact ideal decompositions and unit
elimination. It supplies both multiplicity and embedded-component
control directly; it does not assume either half of the standard
[R_0/S_1 criterion](https://stacks.math.columbia.edu/tag/031R).

The [native audit](type-a-incidence-nilradical-audit.json) checks the
full type A root-exponent character identity, the arbitrary-coefficient
elimination at six block counts, the central-derivative factorization,
and determinant-compatible nonvanishing curves. The independent
[Sage/Singular audit](type-a-incidence-nilradical-verification.json)
checks ideal intersections, radicals and ideal quotients in four
multiplier dimensions. It also constructs the fifteen-dimensional
adjoint matrix for SL4 from defining matrices and tests the ORIGINAL
Laurent gradient ideal, with z^(-1) represented by r and zr=1. It does
not replace the central torus by its function field when checking
radicality.

The SL4 control makes the source/target distinction visible. With Levi
(2,2), S=z^2+z^(-2), T=z^2-z^(-2), and traces x,y,

    Phi=x^2+y^2+Sxy-1,
    J=(2x+Sy,2y+Sx),
    y=-Sx/2, R=C[z^(+-1),x]/(T^2 x),
    N=(Tx), mu=-TSx^2, mu*N=0, mu^2=0.               (111)

Both Tx and mu are nonzero in R. The original root-derived incidence
is therefore genuinely nonreduced, while the generic TARGET fiber is
empty and reduced. This rules out trying to prove global radicality
of every actual adjoint incidence as a shortcut to the desired theorem.
For SL5 with Levi (2,2,1), the independent defining-matrix calculation
instead verifies the radical product/intersection in (109), retaining
the entire coefficient-zero locus. These are checks of the general
proof, not new entries in a Levi census.

The prior [quadratic theorem](QUADRATIC_BRANCHING_FIBERS.md#1-a-quadratic-theorem-without-a-constant-rank-hypothesis)
already proves generic regular-fiber smoothness for type A. We do not
claim (104) as new generic-smoothness coverage, a solution of the
exceptional problem, or an original result in the literature. The
explicit source nilradical and its target-coordinate annihilator are
the additional deductions recorded here. No exhaustive literature
priority audit has been performed for this module description.

### What still prevents an all-simple theorem

The decisive adjoint-specific input to (107) is the exact two-trace
factorization (105). The present exceptional root-character arguments
do not furnish such a factorization or an analogue of the radical
intersection (109) and the annihilator (110). Their nonlinear coefficient
data alone cannot do so, by the earlier altered-Levi control. None of
this constructs a target-polynomial annihilator for the remaining
exceptional nilradicals. The previously checked Sp4 special-target
singularity, smooth E6 positive excess, and persistent altered-character
nilpotents remain valid controls.

Thus the all-simple generic-target reducedness assertion remains
UNPROVED. The next useful task is to obtain a comparable ideal-module
description from the exceptional full character, allowing a general
nonzero p(mu) in (45). Requiring every mu_j to annihilate that nilradical
would be an unjustified strengthening beyond type A. A successful
argument must retain all embedded structure and establish its target
annihilator; another verification of the already-known type A generic
fiber or another R_0/S_1 reformulation would not advance this gap.

## 25 October: target-annihilator attempt and an exceptional radial obstruction

**Outcome: no exceptional target-polynomial annihilator has been derived.**
This attempt tests two proposed ways to obtain (45) from the full root
character: descend the root-discriminant radial field to the critical
algebra, and convert a source annihilator by a norm over the finite
central image. The first fails even for the ACTUAL E6 adjoint character
on a target-dominating regular-semisimple branch. The second needs a
module annihilator not supplied by central finiteness. These are precise
limits of the attempted proof, not a proof of generic reducedness, an
adjoint counterexample, or new Levi coverage. LIMIT-18 records the first
obstruction; no new positive theorem ID is assigned.

### Why radial descent was a substantive additional assumption

On the justified external-discriminant open write the root-derived field
as V=a.partial_x+b.D_z and the full critical ideal as J=(F), F=Phi_x.
In simply laced type V(Phi)=kappa Phi+2 rank(G). Differentiation gives

    V(F_i)=kappa F_i-sum_j (partial_i a_j)F_j
                      -sum_j (partial_i b_j)mu_j.       (112)

Thus V induces a derivation of P/J only if its final term vanishes
modulo J for every i. Without that condition, V cannot be treated as
an operator on the nilradical N of P/J. A scalar eigenidentity for Phi
does not provide it. Nor can one silently replace J by a V-stable
enlargement: that changes the fiber scheme.

One possible repair was to multiply V by a nonzero target polynomial
p(mu). The following actual exceptional calculation shows that even
this repair cannot hold uniformly. This is a test of a general method,
using an existing E6 character, not a new Levi census entry.

### An actual E6 branch on which no target multiple repairs descent

Use the full character (24) with k=1 for the Levi A2+A1+A2 in E6:

    Phi=ab+uv+x^2-2+Sx+(z^(-2)+zx)av+(z^2+z^(-1)x)bu,
    S=z^3+z^(-3), R=z^3-z^(-3), R^2=S^2-4.            (113)

The locus

    a=b=u=v=0, x=-S/2                                 (114)

satisfies ALL five original derived critical equations, with

    mu=D_z Phi=-3SR/2.                                (115)

This is a nonconstant rational function of z, so the branch dominates
the one-dimensional target. With the Killing normalization the actual
target is h=(mu/144)H. Any fixed nonzero scalar normalization leaves
dominance and the no-polynomial conclusion below unchanged.

We compute the root discriminant on the entire trace-zero slice before
substituting x=-S/2. The two SL3 eigenvalue sets are each {1,omega,omega^2}.
Write t,t^(-1) for the SL2 eigenvalues, so x=t+t^(-1). The full list of
external root weights from (113) is

    z^(+/-3)t^(+/-1);
    z^(-2)p_i/r_j, z*t^(+/-1)p_i/r_j;
    z^2 r_j/p_i, z^(-1)t^(+/-1)r_j/p_i.               (116)

There are 4+9+18+9+18=58 external roots; the 14 internal roots give
all 72 E6 roots. Every cube root occurs three times among p_i/r_j.
Taking the products of 1 minus these weights gives exactly

    E_ext=-R^6(S-x)^2(S-x^3+3x)^6.                    (117)

On (114),

    E_ext=-9*S^8*(S^2-4)^9/2^20,                      (118)

which is nonzero generically. The internal discriminant is also nonzero
generically: both SL3 eigenvalue sets are distinct and x^2!=4. Thus this
is an allowed regular-semisimple source branch, not a point removed by
external-root localization. For generic z, mu!=0 and the target has
exactly this proper Levi centralizer.

Use the root inner product with squared root length two, as in (75).
The SL2 and central cocharacter norms on this cover are 2 and 6. Since
E_all equals a nonzero constant times (4-x^2)E_ext on this slice,

    V(x)=(x^2-4)/2 * partial_x log[(4-x^2)E_ext],
    V(log z)=(1/6) D_z log E_ext.                     (119)

These formulas do not divide by an internal root on the quotient; the
apparent internal factor cancels in V(x). They can also be evaluated
directly on the regular torus. For the FULL generator

    F_x=2x+S+z av+z^(-1)bu,

its differential on the trace-zero slice has no components in a,b,u,v.
Consequently (117)--(119), followed by x=-S/2, give

    V(F_x)|_(114)=2*(S^2+320)/(3S) != 0.              (120)

There is no failure of the radial character identity: on this same
branch V(Phi)=22 Phi+12. The failure is preservation of the critical
ideal, a different assertion.

For every nonzero p in C[c], p(mu(z)) is a nonzero rational function,
because (115) is dominant. The product of it with (120) is nonzero in
C(z). Hence

    p(mu)V(J) is NOT contained in J for ANY 0!=p.      (121)

Indeed membership in J would vanish on the actual branch (114), which
the displayed product does not. This also excludes descent after
localization at ANY nonzero target polynomial. Thus the proposed
target-multiple repair of this radial operator fails for the full E6
root character and its own root discriminant together.

This does NOT say p(mu)N cannot vanish. The current E6 example has an
already-controlled generic fiber. Equation (121) concerns one candidate
operator on P/J, not the existence of a different nilradical certificate.
It does not exclude corrected vector fields, additional operator
identities, or a direct ideal-module proof of (45).

### Why a central norm does not supply the missing certificate

Over L=C(c), the existing central-rigidity theorem makes the image B of
the central Laurent algebra in A_K finite dimensional. It does NOT make
A_K finite dimensional. The norm construction requires a finite locally
free algebra: it is the determinant of multiplication by an element
([Stacks, norm construction](https://stacks.math.columbia.edu/tag/0BD2)).
For e in B already known to annihilate the nilradical of A_K,
Cayley--Hamilton does imply that Norm_(B/L)(e) annihilates it. But this
norm is nonzero exactly when e is a unit in B. Proving that fact for
such an e would already kill the nilradical. Central finiteness supplies
neither the required e nor the nonvanishing of its norm. A function of
derived character variables need not even belong to B.

Two existing controls prevent unjustified substitutions in that step:

* The actual E6 contact fiber is L[z,w_1,w_2]/(z^2+1). Its central
  algebra is finite and reduced, but the whole fiber has dimension two.
  The source polynomial z^2+1 has norm ZERO in that central algebra.
  Inverting it erases the entire generic target fiber. A norm of the
  whole fiber algebra as though it were finite is not available.
* The altered-Levi control has, on its existing localized branch,
  algebra B[u,v]/(u^2v,uv^2), with
  B=L[z]/(z^12+(c/6)z^6-1). This B is finite etale over L: its
  polynomial has nonzero generic discriminant and nonzero constant
  term. Nevertheless uv is a nonzero square-zero element. Its
  annihilator is (u,v), whose intersection with B is zero. No nonzero
  central element, hence no nonzero target scalar, kills uv there.
  This is the previously identified honest Levi representation, NOT
  an adjoint character of a simple ambient group.

Thus even a reduced finite central image is insufficient; nilpotents
along the remaining character directions must be controlled. Conversely
the E6 positive-excess example must remain in the proof rather than be
deleted by a source norm. The Sp4 control remains consistent: its
relative incidence is reduced, while special nonzero target fibers can
be nonreduced. A relative nilradical theorem would not rule that out.

### Independent checks, status, and a stopping rule

The [native check](adjoint-annihilator-routes-audit.json) reconstructs
the full E6 character from a separate Cartan-root enumeration, proves
(117)--(120) symbolically over characteristic zero, and checks the norm
controls. The [independent Sage check](adjoint-annihilator-routes-verification.json)
starts from Sage's ambient E6 root lattice. At z=2 it evaluates all 72
roots over an exact number field, reconstructs the root metric and
radial field, and obtains

    mu=-12285/128, V(F_x)=1647/52, E_ext!=0.

It separately verifies the original critical branch and target
dominance, and checks the retained nilpotent and positive-excess ideals
with Singular. No native audit functions are imported. The numerical
field point is an exact corroboration; the all-polynomials conclusion
(121) follows from the written dominant-branch argument.

Published inputs here are the radial/Casimir formula and finite-algebra
norm construction. The explicit failure (121) is a repository deduction,
with no literature-priority claim. There is no new positive theorem,
no new exceptional annihilator, and no increase in group coverage.

The missing step is still an actual root-derived relation in
(J:sqrt(J)) coming from a nonzero target polynomial, including the
embedded associated components. Neither a scalar radial identity nor
the central norm supplies it. Stop trying to descend this uncorrected
radial field by target scaling. A further proof needs a different
ideal/module identity using the full linear and cubic couplings. A
candidate must be tested by exact ideal membership BEFORE using it
as an operator on N; another scalar spectral identity is not progress
on that requirement.

## 26 October: a direct generator-relation annihilator for aligned couplings

**Outcome.** We now have an explicit whole-ideal decomposition, the
entire nilradical, and its EXACT annihilator for the aligned quadratic
family of GEO-78, with arbitrary numbers of pairs and central variables.
This gives a direct target-polynomial certificate for the full E6
character (113), without using the failed radial operator. It does NOT
give a normal form or annihilator for every exceptional restriction.
Generic smoothness for this family was already proved in GEO-78; the
additional deduction is the exact nilpotent-module and membership
description, not new group coverage or a novelty claim (GEO-80).

### Statement: exact nilradical and target annihilator

Work over C. Let B=C[z_1^(+/-1),...,z_s^(+/-1)], s>=1, m>=1, and
let delta,Q belong to B with delta!=0. Put

    A=B[x,u_1,...,u_m,v_1,...,v_m],
    M=(u_1,...,u_m,v_1,...,v_m), w=1-sum_i u_i v_i,
    P=x^2-delta, f=wP+Q,
    J=(xw, P*u_i, P*v_i), R=A/J,
    mu_j=D_j f=D_j Q-w D_j delta, D_j=z_j partial_zj.  (122)

The displayed J is the full derived critical ideal, up to the unit 2
and harmless signs in its generators. Write delta as a unit times
product_l h_l^(e_l), where the h_l are distinct irreducible Laurent
polynomials. Set d=product_l h_l and a=delta/d. If delta is a unit,
take d=1 and a=delta. Thus a retains multiplicities e_l-1, and may
itself be a unit. The use of the letter a here is unrelated to the
E6 fundamental character a in (113).

**Theorem.** With N=sqrt(0)_R, all the following are exact:

    J=(x,M) intersect (P,w) intersect (x,delta),        (123)
    sqrt(J)=J+(d*w*u_i,d*w*v_i),
    (J:sqrt(J))=(x,a).                                (124)

In particular N is generated by the classes of d*w*u_i and d*w*v_i,
and Ann_R(N)=(x,a)/J. If any e_l>=2, the nilpotence index of N is
exactly max_l e_l. Otherwise N=0.

For T=C[c_1,...,c_s] acting by c_j -> mu_j, its exact annihilator is

    Ann_T(N)={p(c): p(DQ) belongs to (a) in B}
      =(a, c_1-D_1Q,...,c_s-D_sQ) intersect T.         (125)

This ideal contains a NONZERO target polynomial. More explicitly, if
a is not a unit, let p_0 be any nonzero polynomial vanishing on the
closure of DQ(V(a)_red) in A^s, and set r=max_l(e_l-1). Then

    p_0(mu)^r N=0.                                    (126)

If Q is constant and N!=0, the exact target annihilator is simply
(c_1,...,c_s); EVERY mu_j annihilates N, regardless of the multiplicities
of delta. If N=0, Ann_T(N)=T and p=1 suffices.

The theorem concerns the actual relative critical algebra. It neither
reduces that algebra in advance nor claims that every special target
fiber is reduced. It allows positive-dimensional fibers. Free affine
variables and an independent quadratic block eliminated with unit
coefficients preserve these assertions. For applications on a finite
central etale cover, the annihilation assertion descends faithfully.

### Proof of the whole-ideal decomposition

In any polynomial ring over any coefficient ring, w=1-sum_i u_i v_i
is a nonzerodivisor: it is invertible in the formal power-series ring
in the u,v variables, into which the polynomial ring injects. This
also holds after taking a quotient of the coefficient ring.

Modulo P, the intersection of (w) and (x) is (xw). Indeed w is a
nonzerodivisor modulo (P,x), which is a polynomial ring over B/(delta).
It follows that

    (P,w) intersect (x,delta)=(P,xw).

The generator xw already belongs to (x,M). Modulo (x,M), P maps to
-delta in the domain B, and delta!=0. Hence ((x,M):P)=(x,M). Thus

    (x,M) intersect (P,xw)=(P*x, P*M, xw)=J.

The final equality follows from the explicit generator relation

    P*x=P*(xw)+sum_i x*v_i*(P*u_i).                   (127)

Every equality holds in A, including w=0, delta=0 and their
intersections. No division by those functions has occurred. This is
why a check on a dense source open is unnecessary here.

### Proof of the entire nilradical and colon ideal

The first factor (x,M) in (123) is prime. The second factor (P,w) is
radical. To see the latter, P is squarefree in the Laurent polynomial
UFD B[x]: it is monic and gcd(P,2x)=1 over Frac(B), since delta!=0.
Thus B[x]/(P) is reduced. On w=0 some u_i or v_i is invertible; these
principal opens cover because sum_i u_i v_i=1. In each open solve for
its paired variable. The resulting ring is a polynomial localization
over B[x]/(P), hence reduced. Therefore A/(P,w) is reduced everywhere.

The exact intersection (123) embeds R into the product of its three
factor quotients. Its nilradical injects into

    A/(x,delta)=(B/(delta))[u,v],

whose nilradical is generated by d. The image of N must also belong
to M and to (w), from the other two factors. In this quotient,

    (d) intersect M intersect (w)=d*w*M.               (128)

For precision, (d) intersect M=dM by looking at polynomial coefficients.
The element w is a nonzerodivisor modulo (d), and w is 1 modulo M.
Thus if wg lies in dM, then g lies in both (d) and M. This proves
(128). Conversely each d*w*u_i and d*w*v_i lifts to the intersection
of the radicals of all three factors in (123). Therefore (128) is
the entire image of N and proves the radical formula in (124).

Multiplication by w*u_1 is injective on (B/(delta))[u,v]. The
annihilator of d in that ring is (a), because delta=da in the domain
B. Hence the annihilator of d*w*M is exactly (a). Pulling it back
to R gives (x,a)/J, proving the colon formula. This controls the whole
nilpotent module; it makes no assumption of unmixedness and does not
discard embedded associated structure. In particular nothing has
been inferred solely from multiplicities at minimal primes.

Likewise the kth power of d*w*u_1 is zero in this injected image
exactly when delta divides d^k in B. This proves the asserted exact
nilpotence index when N is nonzero.

### Direct membership certificate in the original generators

For h any u_i or v_i, put n_h=d*w*h, F_0=xw and F_h=P*h. Then

    a*n_h=x*h*F_0-w*F_h.                              (129)

This is an explicit polynomial identity, not only a vanishing claim.
Every D_j delta is divisible by a, including when the unit factor in
delta is nonconstant. Consequently (mu_j-D_jQ)N=0 by (129). Reducing
modulo the exact colon (x,a) proves (125).

If a is a nonunit, its reduced zero locus in the s-dimensional central
torus has dimension at most s-1. The image of DQ has closure of at
most that dimension, so a nonzero p_0 as in (126) exists. The Laurent
Nullstellensatz gives rad(a) dividing p_0(DQ). Raising it to r makes
it divisible by a, proving (126) by (129). This is target localization,
not inversion of a source polynomial. The norm route of the previous
attempt is not used.

The multiplicity is material. For s=m=1, delta=(z-1)^2(z+1)^3 and
Q=z, (125) gives the exact target ideal

    Ann_T(N)=((c-1)(c+1)^2).

The support polynomial c^2-1 alone does NOT annihilate N; its square
does. The checks below verify both the successful membership and the
failed squarefree one. This distinguishes scheme annihilation from
merely finding a polynomial vanishing on the support.

### Apply the certificate to the FULL E6 adjoint character

For (113), put A_0=z^(-2)+zx and B_0=z^2+z^(-1)x. The global
triangular changes

    a'=a+B_0u, b'=b+A_0v, X=x+S/2

give the exact full-character identity

    Phi=a'b'+(1-uv)(X^2-S^2/4)-2.                    (130)

This uses A_0B_0=1+Sx+x^2 and the coefficient ONE of Sx in the
actual root character. Both audits reconstruct that full character
from the E6 roots. The earlier omitted-root reflection argument
explains the alignment; it is not assumed for other embeddings.
The equations a'=b'=0 eliminate a nonsingular pair globally. Since
S=z^(-3)(z^6+1) has simple zeros on Gm, delta=S^2/4 has exactly double
zeros. The theorem gives, with w=1-uv,

    N=(S*w*u,S*w*v), N^2=0,
    Ann_R(N)=(X,S), Ann_(C[c])(N)=(c), mu*N=0.         (131)

These are ideals/classes in the full critical quotient, with the
eliminated pair understood. The nilradical is NONZERO. On it the
central target acts by zero; the relative critical scheme itself
has not been replaced by its reduction.

In the new coordinates mu=-3SRw/2 modulo the full gradient ideal.
For h=u or v, the certificate can be printed explicitly:

    mu*(S*w*h)
      =-3*R*w*X*h*(2Xw)+6*R*w^2*((X^2-S^2/4)*h).    (132)

Both factors in parentheses are critical generators, up to signs.
Because the changes depend on z, the original central derivative
also acquires multiples of the original affine gradient generators.
The chain-rule identities and membership for that ORIGINAL derivative
are checked separately; (132) is not a change of the Springer map.

The target c here pairs h with the original central cocharacter by
the ambient Killing form (B(H,H)=144). Invertible scalar normalization
does not change (131). Finite central-cover descent and the existing
smooth regular-quotient passage preserve generic reducedness without
removing regular nonsemisimple representatives.

### Controls, independent verification, and exact limits

The assumption delta!=0 cannot be dropped. If delta=0, Q=z and m=s=1,
then J=(xw,x^2u,x^2v), x^2 belongs to J but x does not, and
(J:x)=(x,w). This embedded nilpotent is supported on w=0, and mu=z
is dominant there. No nonzero target polynomial annihilates it. The
step ((x,M):P)=(x,M) fails exactly when delta=0. The independent
ideal audit checks this control as well as the theorem's intersections.

The altered-Levi k=2 character is excluded for a different reason:
after (130) it has an extra term S*X-S^2/2, including a linear term
in X. Its critical generators do not have (122). Its previously
verified embedded nilpotent over a generic target remains intact.
The Sp4 special-target triple point is also consistent: (131) and
(125) concern relative nilpotents, not all special-fiber multiplicities.
The existing E6 contact algebra L[z,w_1,w_2]/(z^2+1) is retained;
neither that source polynomial nor any delta or w is inverted here.
The theorem itself also allows positive excess, for example delta=1,
Q=sum z_j on the w=0 branch, with free quadric directions.

The [native audit](aligned-critical-annihilator-audit.json) checks
(127), (129), (132), the original central-derivative chain rule,
an entire ideal intersection and the sharp mixed-multiplicity target
certificate with SymPy. It separately reconstructs the E6 character
from Cartan roots. The [Sage/Singular audit](aligned-critical-annihilator-verification.json)
uses Sage's ambient roots and checks exact radicals, ideal quotients
and intersections, with one to three pairs and several multiplicities.
For E6 it checks the ORIGINAL Laurent gradient ideal, representing
z^(-1) by r with zr=1, and computes its radical and annihilator colon.
It imports no native audit functions. The arbitrary-m,s proof is the
written ideal argument above, not extrapolation from these checks.

Standard inputs are UFD factorization, the Nullstellensatz and elementary
facts about [reduced rings and schemes](https://stacks.math.columbia.edu/tag/01IZ).
The explicit formulas (123)--(132) are repository deductions; literature
priority has not been established. GEO-78 already implied existence
of a target annihilator by generic smoothness. We are now identifying
that annihilator and the complete nilpotent module directly.

**Remaining gap.** No argument places all exceptional full critical
ideals in (122), or proves an analogous intersection and colon formula
for their overlapping cubic couplings. The earlier E8 determinant-only
failure and altered-Levi control prevent assuming such a reduction.
Thus the all-exceptional and all-simple generic-target reducedness
claims remain unproved. A useful next advance would be an exact
intersection or module presentation for coupled blocks that controls
their ENTIRE colon (J:sqrt(J)); it must not treat each root carrier
as separately critical. This is a concrete successful mechanism to
test against that gap, not a license to resume a Levi census.

## Working rule for the next attempt

Continue with (1)-(5), looking for an adjoint/root identity that excludes
a dominant nonsmooth incidence component. Existing quadratic, norm, and
unit-elimination theorems are available as lemmas, but their special-case
applications are not the endpoint. Any proposed extension must accommodate
the known E6 excess and must not falsely admit the persistent double
structure above.

Report a proof, an actual adjoint counterexample, or the precise remaining
obstruction. Do not call another abstract reformulation a proof of the
conjecture, and do not automatically replace an unsuccessful uniform
attempt with completion of the four E7 subsets.

Primary input: [Kostant--Michor, sections 2 and 4](https://www.mat.univie.ac.at/~michor/kostant.pdf).
A limited search on 24 October for adjoint generalized-Cayley fiber
smoothness did not supply a general theorem. This is not an exhaustive
literature audit or evidence of novelty.
