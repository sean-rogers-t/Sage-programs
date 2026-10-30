# An integral affine extension through evaluation and level

Exploration, 22 October 2026. This note constructs a specified ordinary
cohomology image on an affine flag, and computes its integral arithmetic
in types A/C. It uses the finite Springer source through evaluation.
Section 7 separately records the elementary pointwise loop extension;
the central/rotation Kac--Moody problem is different. The
ambient loop-space and Schubert frameworks are published; the displayed
image formulas are deductions from them, with no literature-priority claim.

## 1. The map and the extra parameter

Let G=SL_n(C), K=SU(n), and identify the affine flag topologically with
LK/T. Evaluation at 1 gives e:LK/T -> K/T. The homeomorphism

\[
 LK/T\cong\Omega K\times K/T,\qquad
 g(t)T\longmapsto (g(t)g(1)^{-1},g(1)T)
\]

identifies e with the second projection. Write eta for the generator of
H^2(Omega K;Z), pulled back to the affine flag. It is the level-one
Schubert divisor epsilon_0; the finite divisor classes pull back as
e^*omega_i=epsilon_i-epsilon_0 in type A. The evaluation/divisor formula
is a published input: [Mare--Mihalcea, Section 4, Proposition 4.2](https://arxiv.org/pdf/1409.3587).

For a finite-dimensional representation rho with faithful differential,
let Lambda_rho be its literal integral linear Springer source, and let
J_rho=beta^B(Sym_Z Lambda_rho) in H^(2*)(G/B;Z). Fix an integer k>=1.
The **evaluation-and-level image** is

\[
 J_{\rho,k}^{\mathrm{aff}}
    =\operatorname{im}\left[
       \operatorname{Sym}_{\mathbb Z}(\Lambda_\rho)[z]
       \longrightarrow H^{2*}(LK/T;\mathbb Z)\right],\quad
 f\longmapsto e^*\beta^B(f),\quad z\longmapsto k\eta.               \tag{1}
\]

The level k is extra data; it is not inferred from the finite normalized
Springer map. For example, repeating a representation leaves that finite
map and its source lattice unchanged but multiplies its trace index.
Any proposal to choose the trace index as level must say so explicitly.
The affine null-root character has zero first Chern class here.

The image (1) is a ring under ordinary cup product. This construction is
different from the Peterson correspondence, whose affine side is
Pontryagin **homology** after specified localization and quotient.

## 2. Content of divisor powers on the affine Grassmannian

Put u=eta on Gr_SL_n, homotopy equivalent to Omega SU(n). The published
integral Hopf model identifies its homology with Z[h_1,...,h_(n-1)] and
its cohomology with the symmetric-function quotient having basis m_lambda,
lambda_1<n. The Hopf pairing is the Hall pairing. See
[Lam, Sections 3--4](https://websites.umich.edu/~tfylam/FPSAC2006/FPSAC2006.pdf)
(extended abstract of the 2008 JAMS paper). Thus u=m_(1), and the content
of its d-th cup power in the integral cohomology lattice is

\[
 c_{n,d}=\gcd_{\lambda\vdash d,\ \lambda_1<n}
                      \frac{d!}{\lambda_1!\cdots\lambda_{\ell}!},
 \qquad c_{n,0}=1.                                                \tag{2}
\]

Content means that u^d=c_(n,d)*v_d for a primitive integral class v_d.
It does not assert that v_d itself is a Schubert class.

**Theorem.** For a prime p set a_p=max{a>=0:p^a<n}. Then

\[
 \boxed{v_p(c_{n,d})=
           v_p\!\left(\left\lfloor\frac{d}{p^{a_p}}\right\rfloor!\right).}
                                                                    \tag{3}
\]

In particular c_(2,d)=d!, c_(3,d)=d!/2^floor(d/2), and the first d with
c_(n,d)>1 is the smallest prime power at least n. Every prime eventually
divides a divisor-power content in a fixed affine type A rank.

### Proof

Fix p and let P=p^(a_p), the largest allowed p-power part. Minimizing
the valuation of a multinomial in (2) means maximizing
sum_i v_p(lambda_i!). Replace each part by the powers occurring in its
base-p expansion, repeated with their digit multiplicities. Legendre's
formula shows that this preserves the sum of factorial valuations, and
every new part is at most P. Repeatedly carry p equal powers below P
into the next power. Each carry increases that sum of factorial
valuations by one. At termination there are q=floor(d/P) parts P and
the base-p digits of r=d-qP below P. This is an allowed partition and
has the largest possible valuation, since any allowed partition can be
improved to this one. Its value is

\[
 \sum_i v_p(\lambda_i!)=\frac{d-q-s_p(r)}{p-1}.
\]

Subtract from v_p(d!)=(d-s_p(d))/(p-1), and use
s_p(d)=s_p(q)+s_p(r). The answer is (q-s_p(q))/(p-1)=v_p(q!), proving
(3). The first positive valuation occurs at d=pP=p^(a_p+1). Taking
the smallest of these prime powers proves the threshold statement. QED.

This calculation explains why a finite-flag bound involving only a fixed
finite set of torsion primes cannot simply be copied to all affine degrees.

## 3. Complete degree-wise image factors on the affine flag

Let b_i=rank H^(2i)(G/B;Z), N=dim_C G/B, and let
s_(rho,i,1),...,s_(rho,i,b_i) be the full finite Springer-image Smith
factors in degree i, including units. Let p_(n-1)(j) denote the number
of partitions of j with all parts at most n-1, with p_(n-1)(0)=1.

**Corollary.** In degree d, the image (1) has rank
sum_(i=0)^min(d,N) b_i. Its ambient quotient is isomorphic to a free
abelian group of rank

\[
 \sum_{i=0}^{\min(d,N)} b_i\bigl(p_{n-1}(d-i)-1\bigr)
\]

plus the finite group

\[
 \boxed{\bigoplus_{i=0}^{\min(d,N)}\
        \bigoplus_{j=1}^{b_i}
     \mathbb Z/\bigl(s_{\rho,i,j}\,k^{d-i}c_{n,d-i}\bigr).}         \tag{4}
\]

The displayed cyclic factors need not already be ordered invariant
factors; prime-wise sorting converts them to Smith order.

**Proof.** Integral Kunneth and the product decomposition split degree d
into finite degree i times loop degree d-i. In that summand, the image is
J_rho^i tensor Z*(ku)^(d-i). Divide u^(d-i) by its content to obtain a
primitive vector, extend it to an integral basis, and apply the finite
Smith basis. This gives (4) and the stated free complement. QED.

The source kernel in (1) is just the finite Borel kernel extended to
the polynomial variable z: the nonzero u^j are rationally independent
in their different loop degrees, and the integral target is torsion-free.
For n=2 and the full weight source, the complete degree-d image factors
are ((d-1)!,d!). This is also the full divisor image, because
e^*omega_1=epsilon_1-epsilon_0 and eta=epsilon_0.

## 4. Independent verification and remaining questions

The [arithmetic pilot](experiments/other_theories_affine_pilot.py) compares
(3) with a direct gcd recurrence
c_d=gcd_(1<=j<=min(n-1,d)) binomial(d,j)c_(d-j).
It checks 6,231 pairs n=2,...,32 and d=0,...,200. Separately it computes
the affine A1 divisor image through degree 100 from real coroots and
the infinite-dihedral Chevalley formula. That method does not use the
symmetric-function Hopf presentation and gives ((d-1)!,d!) in every
checked degree. [Report](experiments/other-theories-affine-pilot.json).

The theorem is the argument above, not an extrapolation from those tests.
Its arithmetic priority has not been established; the rank-one factorial
phenomenon is already the classical divided-power phenomenon for Omega S3.

Substantial next questions are other finite root types, additional
independent loop generators, actual affine-parabolic character sources,
and compatibility with a chosen representation-induced central extension.
The resulting maps should first be called characteristic-image maps.
A central/rotation-compatible affine Springer construction needs
additional work beyond the pointwise loop map described below.

## 5. A closed symplectic analogue

There is a type C formula as well. On Gr_Sp_(2n), normalize u to the
primitive degree-two Schubert class and let c_(C,n,d) be the content of
u^d. The published integral Hopf model is
H_*=Z[P_1,P_3,...,P_(2n-1)], with
Delta P_j=P_j tensor 1+1 tensor P_j+2 sum_(0<a<j) P_a tensor P_(j-a).
These precise statements are
[Lam--Schilling--Shimozono, Theorems 1.3 and 5.5](https://arxiv.org/abs/0710.2720).
They concern the established ambient Hopf algebra, rather than a new
Springer construction.

Iterating the coproduct and pairing with u gives
<u^j,P_j>=2^(j-1). Therefore, in the polynomial monomial homology basis,

\[
 c_{C,n,d}=\gcd_{\lambda\vdash d,\ \lambda_i\text{ odd}<2n}
       {d!\,2^{d-\ell(\lambda)}\over\prod_i\lambda_i!}.             \tag{5}
\]

**Theorem.** Its valuations are

\[
 \boxed{v_2(c_{C,n,d})=v_2(d!),\qquad
 v_p(c_{C,n,d})=v_p\!\left(\left\lfloor d/p^{a_p}\right\rfloor!\right)
 \quad(p\text{ odd}),\quad p^{a_p}<2n\le p^{a_p+1}.}              \tag{6}
\]

For odd p the type A carrying proof applies unchanged: every p-power
part used in that proof is odd and allowed, and the factors two are
units. At two, Legendre's formula makes the valuation of a term in (5)
equal to v_2(d!)+sum_i(s_2(lambda_i)-1). Each summand is nonnegative;
the all-ones partition attains zero. This proves (6) in every rank and
degree. In particular u^2 has content two in every symplectic rank.

The odd-primary contents agree with type A of size 2n. Their quotient
is precisely

\[
 {c_{C,n,d}\over c_{2n,d}}
       =2^{v_2(d!)-v_2(\lfloor d/2^a\rfloor!)},\quad 2^a<2n\le2^{a+1}.
\]

This is an arithmetic comparison of the specified primitive divisor
powers. It does not assert an unproved identification of arbitrary
affine Springer sources under a group embedding.

The direct weighted-multinomial recurrence checks (6) in another 4,824
rank/degree cases, n=1,...,24 and d<=200. An
[independent affine Chevalley audit](experiments/audit_other_theories_affine_c.py)
uses real coroots and Grassmannian minimal cosets for affine C2 through
degree 20 and C3 through degree 12; it does not use the Hopf coproduct
or factorial formula. [Audit report](experiments/other-theories-affine-c-audit.json).
The closed formula is the deduction above, and its priority remains
unestablished.

## 6. Partial affine flags and other finite groups

The evaluation-and-level construction extends to LK/K_L, the topological
model of an affine partial flag whose retained nodes lie in the standard
finite diagram. It splits as Omega K times K/K_L. Replace the full-flag
finite factors in (4) by the old factors for G/P_L; the same proof gives
the complete quotient. The type C formula (6) supplies the corresponding
loop content for the symplectic cases. Parabolics containing the affine
node require a different cut or an additional construction, and are not
implicitly covered with the same source choices.

For every finite simply connected G the product argument still works
with its actual primitive divisor content c_(G,d) and its loop Betti
numbers. Affine Schubert cells supply torsion-free even lattices. What
is missing outside the computed A/C families is the closed content
formula, not the elementary Kunneth convolution. Additional loop
generators, beyond one divisor, remain a more substantial image problem.

## 7. An actual pointwise loop map already exists

It would be too strong to say that no loop Springer map can be defined
without an infinite-dimensional trace. Any finite Springer morphism is
a morphism of C-schemes, so for every commutative C-algebra A it already
gives theta_A:G(A) -> g tensor A. Apply this to A=C[z,z^(-1)], C[[z]]
or C((z)). The resulting map is

\[
          L\theta_\rho:\gamma(z)\longmapsto\theta_\rho(\gamma(z)).
\]

It is conjugation-equivariant, has identity differential at 1, and is
an ind-morphism for algebraic loops. It is just the loop functor applied
to the existing finite map, not a new theorem about loop spaces.
It also commutes with loop rotation; extending to the rotation/derivation
direction as part of a fully conjugation-equivariant enlarged group map
is a different requirement.
For the standard loop/arc functor framework see
[Zhu, Section 1](https://arxiv.org/abs/1603.05593).

There is an equally concrete trace interpretation: extend the finite
trace-orthogonal projection End(V) -> g by scalars to A. The module
V tensor A still has finite rank over A, so its A-valued matrix trace
is defined. On Laurent loops, taking the constant Laurent coefficient
of that trace gives the usual nondegenerate loop-algebra pairing.
No ordinary C-linear trace of an infinite matrix is being used.

For a diagram-stable representation, the finite sigma-equivariance
proved in the [folding note](OTHER_THEORIES_FOLDING.md) shows that this
also preserves twisted loops gamma(zeta*z)=sigma(gamma(z)). Thus raw
untwisted and twisted loop maps are available. What is not supplied is
a compatible map on the central extension plus rotation, or the
appropriate polynomial-character source on an affine flag.

Nor does pointwise theta automatically descend through right Iwahori
cosets. This already fails for finite SL2: with
g=[[1,0],[1,1]] and b=diag(s,s^(-1)) in B, the lower-left entry of
theta(gb)-theta(g) is s-1. It lies outside Lie(B) for s!=1. The same
constant-loop example prevents descent LG/I -> Lg/Lie(I). This is why
the characteristic-map/source problem remains substantive even though
the pointwise loop morphism itself is immediate.

## 8. Finite jet approximations retain the old degrees

For A_m=C[z]/z^m, m>=1, the same base change gives a finite-dimensional
map theta_(rho,m):J_mG -> J_m g. Here J_m means jets of length m
(some conventions call this J_(m-1)). The jet group is generally
nonreductive for m>1. Its map still has identity differential at 1.

Two useful identities follow formally, with no new exceptional counting:

\[
 \deg\theta_{\rho,m}=\deg\theta_\rho,\qquad
 \Psi_{\rho,m}(\gamma)=\Psi_\rho(\gamma(0))^m.                    \tag{7}
\]

For the degree, restrict the original dominant generically finite map
to a dense target open where it is finite etale. An etale preimage of
the constant jet lifts uniquely along nilpotent thickenings, so every
one of the original generic preimages has exactly one length-m jet
lift. The number stays the old degree, rather than its m-th power.

For the determinant, use left translation on J_mG and identify its
tangent space with g tensor A_m. The differential is an A_m-linear
matrix whose constant term is the old left-translated differential at
gamma(0). As a C-linear matrix it is block triangular with that same
constant block repeated m times. Its determinant is therefore the
second formula in (7).

These are standard jet/etale functoriality deductions applied to our
existing finite morphisms, not new generic-degree computations or a
harmonic-separation theorem for a nonreductive jet group. They provide
a controlled finite approximation to the loop problem. For compatible
P the jet quotient has contractible unipotent fibers over G/P, so this
alone does not produce new ordinary flag cohomology.
