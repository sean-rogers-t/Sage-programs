# Generic Springer fibers under polynomial composition of characters

27 October 2026. A written deduction from the character critical-scheme
presentation, central rigidity, and the completed adjoint theorem. The
character chain rule and tensor-power formula are published inputs, not
new constructions. Literature priority for the assembled fiber theorem is
not established. This note concerns actual schemes over C.

## 1. The general closure theorem

Let Phi be regular on A^ell x T, T=(C*)^s, with s>0. Set

    Y=V(Phi_x1,...,Phi_xell),  mu_j=D_j Phi,  D_j=z_j partial_z_j.

Assume the image of the nonsmooth locus of the ACTUAL scheme Y under mu
has nondense closure in A^s. Equivalently, its actual generic-target
fibers are smooth after shrinking the target open: for this equivalence,
generic flatness plus smooth geometric fibers gives a smooth morphism,
while a smooth source over C has smooth generic fibers in characteristic
zero. Let f in C[t] be a
nonconstant polynomial and put

    Psi=f(Phi),
    Y_f=V(Psi_x1,...,Psi_xell),  nu_j=D_j Psi=f'(Phi)mu_j.

**Theorem.** The actual generic-target fibers of nu:Y_f -> A^s are
smooth, allowing emptiness and excess. In particular, polynomial
composition preserves the generic-smoothness property in this setting.

Proof. Write a=f'(Phi). Over a nonzero target c, some c_j is nonzero,
and the equation a mu_j=c_j makes a a unit IN THE ACTUAL FIBER RING.
Consequently it is legitimate to work on a!=0. There

    (Psi_xi)=(Phi_xi),   Y_f=Y,

as schemes. Thus the possible bad source locus is contained in Sing Y
on this open. It remains to check its image under nu, not assume that
an arbitrary variable rescaling preserves nondominance.

For any irreducible subvariety C of Y_red, the morphism

    (mu,Phi): C -> A^s x A^1

has image of the SAME dimension as mu(C). Indeed on a fixed mu fiber
the relative critical identity gives

    d Phi=sum c_j dlog z_j.

On every reduced irreducible fiber component, a rational exact
differential equal to a constant logarithmic differential forces Phi
to be constant. The proof uses a normal projective model: a pole of
order m in Phi would give a pole of order m+1 in dPhi, while logarithmic
differentials have only simple poles. There are no poles, so Phi is
constant. Equivalently, over the geometric generic point of mu(C), Phi
is algebraic over the target function field. This is the existing
[central-rigidity argument](CENTRAL_RIGIDITY_FIBERS.md#2-proof-from-the-trace-identity).

Since nu=f'(Phi)mu factors through (mu,Phi), it follows that

    dim closure(nu(C)) <= dim closure(mu(C)).

In fact there is equality whenever C is not contained in a=0. Apply
the SAME residue argument to Psi=f(Phi) and nu. The value f(Phi) is
algebraic over C(nu), so Phi is algebraic over C(nu) as well, because
f is nonconstant. On a!=0 the two finite extensions coincide:

    C(mu,Phi)=C(nu,Phi).

Both target images therefore have the same dimension. This argument
works for every irreducible subvariety C of Y_red, including a component
of the nonsmooth locus; it does not require C to be clean or Lagrangian.

Apply this to every reduced irreducible component of Sing Y. Their
nu-images are nondense by the original hypothesis. Remove these images
and the zero target, then apply characteristic-zero generic smoothness
to the remaining smooth source. This proves the actual-fiber assertion;
it does not replace either Y or its fibers by their reductions.

The use of C_red in the dimension bound only computes the image of the
underlying nonsmooth locus. Every nilpotent and embedded component on
the retained source is controlled by the actual smoothness criterion.

## 2. Exactly which nonempty components survive

For a clean component Y_i of Y, use the critical-value variety V_i and
the horizontal stabilizer criterion in
[the nonemptiness note](ADJOINT_NONEMPTINESS_CRITICAL_VALUES.md).
If Y_i is not contained in f'(Phi)=0, its open part survives as a clean
component for Psi. Its new critical-value variety is the image of V_i
under the finite map

    T x A^1 -> T x A^1,   (z,t) -> (z,f(t)).

The connected horizontal stabilizer is unchanged. One inclusion follows
by taking the image. For the converse, a connected torus preserving the
image preserves its finite inverse image and therefore each of its
finitely many irreducible components. V_i is such a component, because
the map is finite and the dimensions agree. Hence it preserves V_i.

The linear-image theorem therefore gives the SAME central target-image
subspace for the surviving component. In particular, it dominates for
Psi if and only if it dominated for Phi, and its generic excess

    dim Y_i-s

is unchanged. No new generic nonempty components arise from f'(Phi)=0:
their target image is zero and is removed.

Thus a formerly dominating component is lost precisely when Phi is
constant on it, say Phi=b, and f'(b)=0. Nonconstant critical-value
components always survive. Generic nonemptiness holds if at least one
dominating component survives. The dimensions of the surviving generic
components are unchanged; their numbers, degrees and special fibers
need not be unchanged.

The dominance and dimension portions of this conclusion also follow
independently from the equality of function fields up to finite extension
in section 1. The stabilizer argument supplies the stronger equality of
the actual linear target-image subspaces on clean components. These are
two separate checks, not two computations using the same elimination.

The condition must be applied component by component. Neither global
linear independence of the nu_j on a reducible union nor the stabilizer
of that entire union is a substitute.

For a semisimple product, the rule is particularly simple. The adjoint
character and critical equations split as sums and products, respectively.
A product component dominates exactly when every factor component does.
Its total character is constant exactly when each factor character is
constant, and then its value is the sum of those constants. If any factor
has a nonconstant value, no nonzero f' can erase that product component.
The zero-fiber calculation for a full simple-factor Levi contributes
constant -1 in type A and an empty regular fiber in every other type.
Finite central descent does not change these conclusions. Thus the
semisimple extension also retains an explicit componentwise interpretation;
it is not an assumption that a polynomial in a sum of characters splits
as a product of Springer maps.

There is also a constructive consequence within actual representations.
Given any finite set S of nonpositive rational numbers, write its distinct
elements as -p_j/q_j, p_j>=0 and q_j>0 integers, and let k=|S|. Define

    f(t)=(k+1)! integral_0^t prod_(j=1)^k (q_j u+p_j) du.

Every coefficient of f is a nonnegative integer: the integrand has
nonnegative integer coefficients, and (k+1)! clears each integration
denominator. Its derivative vanishes EXACTLY at S and is positive at
every positive representation dimension. Thus f(chi_Ad) is an actual
nondegenerate representation character that removes precisely the
dominating components with constant critical value in S, while preserving
the other dominating components. This includes simultaneous choices
across different Levis; components sharing a value cannot be separated
by this operation. Nonconstant-value components always survive.
Conversely a nonconstant polynomial with nonnegative integer coefficients
has positive derivative at every positive real number, so such operations
cannot remove a component whose constant value is positive real. No
rationality assertion about all exceptional critical values is being made.

## 3. Actual representations and normalization

Let G be connected complex semisimple, rho a rational representation
with nondegenerate trace form, d=dim rho, and chi=chi_rho. Suppose
f(chi) is an actual representation character chi_sigma and f'(d)!=0.
Then

    B_sigma=f'(d) B_rho,
    theta_sigma = f'(chi)/f'(d) theta_rho.

The first equality follows by differentiating twice at the identity:
d chi(e)=0 because the Lie algebra is semisimple, so the f'' term
vanishes. Trace differentiation proves the second equality. In
particular, every polynomial f with nonnegative integer coefficients
and positive degree is actual, via direct sums of tensor powers, and
f'(d)>0. Trivial summands are allowed.

The tensor-power specialization is already in Kostant--Michor,
[Theorem 2.8(3)](https://www.mat.univie.ac.at/~michor/kostant.pdf).
The polynomial chain rule was also already recorded in
[the repository's general critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md#6-general-multiplicity-corollaries-and-published-baselines).
The present deduction is preservation of the generic regular-fiber
property, with the exact nonemptiness qualification in section 2.

**Adjoint consequence.** For every connected complex semisimple G,
every proper Levi K, and every actual character f(chi_Ad) with
f'(dim G)!=0, the corresponding normalized Springer map has smooth
actual regular-source fibers over a dense open of Z_K^o. This follows
from the completed adjoint theorem, its central-isogeny/product extension,
and section 1. It includes all tensor powers and all direct sums of
tensor powers of the adjoint, in every simple type and rank.

Scaling a nonzero central target by a nonzero scalar preserves its
centralizer, so passing between the trace normalizations does not change
the centralizer stratum. The group-fiber presentation and regular lifts
are the same scheme-theoretic reduction as for the original theorem.

The restriction to proper K is essential for this extension. At K=G
the only target is zero, and the newly created character-level components
cannot be discarded. Tensor powers can give them arbitrary multiplicity,
as already explained by the published chain rule and the earlier
zero-fiber multiplicity corollary. Nor does the theorem extend to arbitrary
representations: the actual reducible SL3 control in
[LIMIT-20](WEYL_JET_OBSTRUCTION.md) still has generic nonreduced fibers.
Its starting character does not satisfy the required generic-smoothness
hypothesis.

More precisely, the proof gives a useful obstruction to repairing such
examples by a character polynomial. On every dominating component of
Sing Y with NONCONSTANT Phi, all nonconstant polynomial compositions
retain a dominating nonsmooth component. A composition can remove a
dominating bad component only if Phi is constant on it, at a root of
f'. For the actual SL3 control in LIMIT-20, Phi=S^4 on its dominating
nonreduced component, with S=z^3+z^-3. It is nonconstant. Hence the
failure persists for EVERY nonconstant polynomial composition that is
an actual character with nondegenerate trace form, not only its tensor
powers. This is a deduction from that existing counterexample and the
dimension-equality lemma; it is not a new irreducible-representation
counterexample.

## 4. A genuine change of nonemptiness: the E6 contact stratum

For the E6 adjoint contact Levi, use the existing full character

    Phi=a b+z^2+z^-2+(z+z^-1)u

with two further free derived coordinates. Put F=z+z^-1 and

    q=Phi+2=a b+F u+F^2.

The old relative critical ideal is exactly (a,b,F), so Phi=-2 in its
actual quotient ring. It gives the known smooth nonempty generic regular
group fibers of dimension 32, with transverse excess two.

For a pure adjoint tensor power f(t)=t^m, f'(-2)!=0. Its nonzero-target
regular fibers in this stratum are therefore the old fibers with the
constant target rescaling m(-2)^(m-1) in raw character coordinates.
After normalization the scale is (-2/78)^(m-1). Nonemptiness and the
excess two survive.

In contrast, take the ACTUAL representation

    sigma=(Ad direct-sum 1^2)^(tensor m),  m>=2.

Its character is q^m, its dimension is 80^m, its trace form is
m*80^(m-1) B_Ad, and

    theta_sigma=(q/80)^(m-1) theta_Ad.

Every nonzero centralizer-exact target in the contact stratum has EMPTY
regular-source fiber for sigma. This is an actual scheme assertion.
For a direct ideal certificate set sigma0=m q^(m-1), let c be the raw
central target, and write R=D Phi. The relative derivatives satisfy

    a*(q^m)_a+(u+F)*(q^m)_u=m q^m.

Consequently q^m belongs to the relative critical ideal. In the full
fiber quotient c=sigma0 R, so

    c^m=m^m q^(m(m-1)) R^m=0.

Since m(m-1)>=m, this uses the original ideal, with no radical. A
nonzero numerical c is a unit, so the entire quotient is zero. This
independently proves emptiness and demonstrates why unconditional
preservation of nonemptiness under tensor operations would be false.

The example does not weaken the smoothness theorem: an empty fiber is
smooth. It identifies exactly how the critical VALUE, in addition to
the critical locus, controls these generalized Cayley maps.

## 5. An all-rank type A consequence

For G=SL_N, let K be a proper Levi with defining-space block sizes b_i.
For EVERY actual character f(chi_Ad) with f nonconstant and f'(N^2-1)!=0,
the generic regular-source fiber over the centralizer-exact K stratum is
nonempty precisely when

    every b_i is at most two, and at least one b_i is one.

When nonempty, its pure dimension is twice the number of two-blocks.
In particular this holds for every nonconstant polynomial with nonnegative
integer coefficients. Thus the entire existing type A nonemptiness and
dimension classification survives these representation operations, even
though the E6 contact classification does not. Generic degrees may change.

Proof. Empty generic strata stay empty by the dimension comparison in
section 1. Suppose there are q two-blocks and t>=1 singletons, N=2q+t.
The case q=0 is the maximal torus: chi_Ad is nonconstant and its logarithmic
gradient is dominant, so section 2 gives the assertion.

For q>=1, on the usual finite central cover write the two-block coordinates
as z_i, their SL2 characters as x_i, and the singleton coordinates as v_j.
The determinant constraint is prod z_i^2 prod v_j=1. Put

    U=sum z_i x_i + A,     A=sum v_j,
    V=sum z_i^-1 x_i + B,  B=sum v_j^-1,
    Phi=UV-1.

The derived critical equations are z_i V+z_i^-1 U=0. The locus U=V=0
has central target zero, so it cannot contain the generic point of a
dominating component. Elsewhere U and V are both nonzero, and these
equations force all z_i^2 to be equal. On each sign branch write
z_i=epsilon_i z and X=sum epsilon_i x_i. Then

    Phi=X^2+(zB+z^-1 A)X+AB-1,
    X=-(zB+z^-1 A)/2,
    psi=-1-(zB-z^-1 A)^2/4.

These are actual scheme eliminations: multiply the i-th equation by z_i
and subtract the first. The result is (z_i^2-z_1^2)V. Since V is a unit
on this open, the central equations are exactly z_i^2=z_1^2. Their sign
branches are reduced and disjoint, because 2z_1 is a unit. The remaining
equation is 2zX+z^2 B+A=0 and has a unit coefficient.

The central branch S is the irreducible torus z^(2q) prod v_j=1.
The displayed single affine equation for X leaves q-1 free coordinates.
This describes the whole open part U,V!=0 of the actual critical scheme,
not just a chosen reduced subvariety. Each branch is smooth and
irreducible of dimension q+t-1. Its closure is therefore a clean
irreducible component, and no smaller hidden or embedded component
occurs on this open.

The critical value is NONCONSTANT on every branch. Restrict to the
one-parameter subtorus z=u^t and v_j=u^(-2q). It satisfies the determinant
constraint and gives

    psi=-1-t^2*(u^N-u^-N)^2/4,

which is nonconstant. Thus no component with dominant central image can
be erased by any nonzero polynomial f'. Section 2 preserves its excess.
The existing type A theorem gives pure excess zero and group-fiber
dimension 2q, completing the proof. Equivalently, on these branches
there are q-1 affine directions and q-1 central compatibility conditions;
the difference is zero. The argument applies after central-isogeny
descent as well, since every polynomial in chi_Ad descends.

Only the representation-extension statement is a new deduction in this
repository. The original adjoint classification and character factorization
are inputs. The symbolic completed-square identity and its one-parameter
restriction are independently checked in native arithmetic and Sage;
neither a finite rank table nor a generic numerical fiber is used as
the proof of the all-rank assertion.

## 6. Four-dimensional blocks: a uniform B/C/D contrast

Fix q>=1 and take the following simple simply connected groups and Levis:

| Type | G | Defining-matrix image of K | Constant b_q |
| --- | --- | --- | --- |
| B_(4q) | Spin_(8q+1) | GL_4^q times SO_1 | -q/(2q-1) |
| C_(4q) | Sp_(8q) | GL_4^q | 0 |
| D_(4q) | Spin_(8q) | GL_4^q | 0 |

For Spin, K denotes the Levi above this matrix subgroup, with its finite
central cover retained. For every h in Z_K^o the adjoint character is
identically b_q on the ACTUAL regular-source fiber. Consequently, for
every actual character f(chi_Ad) with f'(dim G)!=0:

* If f'(b_q)!=0, its regular fiber over h is exactly the adjoint regular
  fiber over h/a, where a=f'(b_q)/f'(dim G).
* If f'(b_q)=0, its regular fiber over h is empty.

The original adjoint generic fibers are nonempty, smooth and pure of
dimension 12q by the existing classical theorem. Thus the first alternative
preserves this generic dimension and nonemptiness; the second erases the
whole centralizer-exact stratum, including its possible nongeneric targets.
The statements do not concern the nonregular source.

Here is a direct scheme-level proof, independent of the stabilizer theorem.
Write T for the defining trace, S2 for the defining trace of g^2, and M_i
for a GL_4 block. The published trace construction and standard exterior/
symmetric-square identities give the matrix formulas recorded in
[the classical proof, section 3](CLASSICAL_ADJOINT_FIBERS.md#3-the-adjoint-matrix-identity-and-a-spectral-lemma).
On a regular fiber each GL_4 block is cyclic and satisfies

    p_(tau,L_i)(M_i)=0,
    p_(tau,L)(u)=u^4-tau*u^3+L*u^2+tau*u-1,
    tau=T in B/D, and tau=-T in C.

A cyclic four-dimensional block annihilated by this monic quartic has
characteristic polynomial exactly p. This holds over the actual fiber
ring: cover the regular locus by opens where a cyclic-vector matrix is
invertible, and compare the polynomial coefficients in its cyclic basis.
No distinct-root assumption or reduction is used. The Vieta identities are

    tr M_i = tr M_i^-1 = tau,
    tr M_i^2 = tau^2-2L_i,
    tr M_i^-2 = tau^2+2L_i.

For B the extra SO_1 eigenvalue is 1. Hence

    T=2q*T+1,  S2=2q*T^2+1,
    chi_Ad=(T^2-S2)/2=-q/(2q-1).

For D there is no residual block: T=2q*T forces T=0 and S2=0,
so chi_Ad=0. For C, T=-2q*T forces T=0 and S2=0, and
chi_Ad=(T^2+S2)/2=0. All divisions here are by nonzero integers over C.
These identities therefore hold in the actual fiber rings, not just on
their complex points. In particular the values are independent of q
target parameters L_i and of collisions among the quartic roots.

For the composed map a(g)theta_Ad(g), every nonzero target makes a(g) a
unit locally on the actual fiber, as in section 1. The same matrix
argument applies with the variable scalar target h/a(g) and gives
chi_Ad=b_q there. Thus a(g)=f'(b_q)/f'(dim G) is constant. If it is zero,
the unit requirement is impossible; otherwise the equations reduce
exactly to the rescaled adjoint equations. Central-isogeny descent retains
these identities and the regular source open.

This supplies unbounded-rank ACTUAL representation families with erased
strata. In C and D take Ad^(tensor m), m>=2. In B take

    ((Ad)^(direct-sum (2q-1)) direct-sum 1^q)^(tensor m), m>=2.

Its character is ((2q-1)chi_Ad+q)^m, with nonnegative integer coefficients
and derivative zero at b_q. The trace form remains nondegenerate.
The type B constants also show that constant critical values need not be
integers: the fraction -q/(2q-1) is reduced, and its denominator is
unbounded as the rank grows. This rules out a finite rank-independent
list of integer critical values as an explanation of the phenomenon.
These are counterexamples to preservation of generic nonemptiness under
tensor operations, not to generic smoothness: an empty fiber is smooth.
They show that the E6 erasure is part of a general representation-dependent
phenomenon. No new generic degree or irreducible-representation claim follows.

The quartic identity, cyclic coefficient comparison and trace relations
are checked independently using a symbolic companion matrix and reciprocal
Newton sums. The proof is in arbitrary q; finite computations corroborate
the symbolic identities and do not establish its rank range.
As a separate check of normalization and the actual map, the native audit
constructs the complete adjoint matrices of dimensions 36, 36 and 28
for rational regular semisimple elements in B4, C4 and D4. It verifies
tr(Ad(g) ad X)=tr(ad(h) ad X) on every Lie-algebra basis vector and obtains
character values -1,0,0 and target scalars 3/14,-3/20,1/4, respectively.
This check uses the defining trace projection, not the tensor-square
character or classical matrix formula as an assumed output.

## 7. Scope and next mathematical question

The input from Kostant--Michor is the trace construction and character
chain rule. The other standard inputs are residue calculus, generic
smoothness, finite morphisms and the regular adjoint quotient. The new
repository statements are the closure theorem and its componentwise
survival rule. They are written deductions with exact controls, not a
claim to settle a named published open problem or establish priority.

For the adjoint nonemptiness program, this identifies the next useful
invariant: the constant critical values on target-dominating components.
A root-theoretic description would determine which tensor-polynomial
representations retain each nonempty stratum, while the stabilizer
criterion handles the nonconstant components. A complete classification
of those values for all exceptional Levis is not established here.
