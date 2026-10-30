# Generic nonemptiness and excess from critical values and torus symmetries

27 October 2026. Written deductions with two proof routes for the central
linear-image lemma and independently implemented exact controls. The aim
is a structural criterion, not another Levi census. The completed adjoint
smoothness theorem is an input, not a claim reproved here. This is an
internal audit, not external peer review or a literature-priority claim.

## 1. The construction and the linear-image theorem

Let T=(C*)^s, with invariant derivatives D_j=z_j partial_z_j, and let
Phi(x,z) be a regular function on A^ell x T. Keep the actual scheme

    Y=V(Phi_x1,...,Phi_xell),   mu=(D_1 Phi,...,D_s Phi).

Call an irreducible component Y_i clean if it meets the smooth locus of
the ACTUAL scheme Y. Equivalently, its generic local ring has length one;
over C its reduced support is then smooth on a dense open. This condition
is not the assertion that the reduced support is generically smooth, which
would be automatic and would lose the required information.
Here "clean" is only an abbreviation for this stated scheme condition;
it is not an appeal to another use of that word in microlocal geometry.

For a clean component define its critical-value variety

    V_i = closure{(z,Phi(x,z)) : (x,z) in Y_i} in T x A^1.

Let H_i be the identity component of the subgroup of T preserving V_i
under (z,t) -> (h z,t). Thus the scalar critical value is fixed, not rescaled.

**Linear-image theorem.** The target image closure of a clean component is

    closure(mu(Y_i)) = (Lie H_i)^perp in Lie(T)^*.

In particular it is a rational LINEAR subspace, not an arbitrary nonlinear
subvariety. The component dominates the full central target precisely when
H_i is trivial. In that case its generic transverse fiber dimension is

    e_i=dim Y_i-s.

Sections 2-4 prove this criterion from algebraic symplectic geometry.
It does not assert that all components of an arbitrary critical scheme
are clean. The actual-character counterexamples and the SL4 control in
section 7 make that distinction necessary.

## 2. Exact Lagrangian image of a clean critical component

Write Q=Phi_xx, B=(D_j Phi_xi), and let k=rank Q at a general smooth point
of Y_i. Put p=rank[Q B]-k. The tangent space of Y is the kernel of [Q B].
Consequently

    dim Y_i=ell+s-k-p,
    rank(dz|Y_i)=s-p.

On the vertical tangent ker Q, the differential of mu is B^t|ker Q. By
symmetry of Q, this map has rank p. It follows that

    gamma:Y_i -> T x Lie(T)^*,   gamma=(z,mu)

has rank (s-p)+p=s at the general point. Its image closure Lambda_i has
dimension s. This calculation does not assume the original equations form
a regular sequence; it uses actual smoothness and their full Jacobian.

Let lambda=sum c_j dlog z_j be the tautological one-form on T^*T, identified
with T x Lie(T)^*. The relative critical equations give

    gamma^*lambda=d(Phi|Y_i),   gamma^*(d lambda)=0.

Generic separability makes pullback of rational differential forms injective
on the image. Thus Lambda_i is Lagrangian. Moreover lambda is exact as a
RATIONAL form on Lambda_i. Here is the descent detail, which should not be
omitted when gamma has positive-dimensional fibers.

Set F=C(Y_i), E=C(Lambda_i). Relative to E, dPhi=0, so Phi is algebraic
over E in characteristic zero. In the finite extension E(Phi)/E, trace of
differentials gives

    d(Tr(Phi)/[E(Phi):E])=lambda.

The difference between Phi and this rational primitive has zero absolute
differential, and is therefore in C. Hence Phi itself lies in E. Denote
the descended rational function by psi; then lambda=d psi on Lambda_i.

The map gamma has generic fiber dimension

    e_i=dim Y_i-s=ell-k-p >= 0.

This is the excess that can remain after the moment map on Lambda_i is
generically finite. It is distinct from a positive-dimensional fiber of
the moment projection of Lambda_i itself.

## 3. Algebraic Lagrangians in a torus cotangent bundle

**Lemma.** Let Lambda be an irreducible closed algebraic Lagrangian in T^*T,
and let D be the closure of its image under the moment projection to
Lie(T)^*. Then D is an affine translate of a rational linear subspace.
If lambda|Lambda has a rational primitive, that translate passes through
zero. In this exact case

    D=(Lie Stab_T(Lambda)^0)^perp.

Proof. Put d=dim D. At a general smooth point (z,c) of Lambda above a smooth
point c of D, symplectic orthogonality gives

    dlog_z(ker d(moment))=(T_c D)^perp.

Both sides have dimension s-d, so this is equality. Thus every general
irreducible fiber component in T has constant logarithmic tangent space.
Such a component is a translate of an algebraic subtorus: locally its
logarithm is an affine linear space; every Laurent polynomial vanishing on
that local leaf vanishes on its entire analytic exponential, and hence on
its Zariski closure, an algebraic subtorus translate. Dimension forces
equality. Taking the closed fiber, rather than an arbitrary open chart,
justifies this statement about its closure.

There are only countably many subtori of T. The tangent Gauss map of the
irreducible variety D has constructible image, consisting of rational
subspaces. Over C a positive-dimensional constructible set is uncountable.
Therefore the tangent space of D is constant on a dense open. Integrating
its annihilating linear functions shows that D is an affine linear space.
Let H be the subtorus whose Lie algebra is its orthogonal direction.
General moment fibers are unions of H-cosets, so H preserves Lambda.
Conversely every connected horizontal stabilizer has tangent in the
moment-fiber kernel. Thus H=Stab_T(Lambda)^0.

If lambda=d psi, restrict to a general H-coset in a moment fiber. This
gives d(psi|H)=c(dlog h). An invariant one-form on a torus cannot have a
nonconstant rational primitive: on each one-parameter subgroup a nonzero
multiple of dt/t has a nonzero residue, whereas a rational differential
df has residue zero. Hence c annihilates Lie H. The affine translate is
therefore exactly (Lie H)^perp. This proves the lemma.

No positivity of coefficients, polynomial degree restriction, or adjoint
root identity is used in this lemma. Its role in the present problem is a
deduction from the existing character presentation, not a replacement for
the difficult adjoint smoothness proof.

### 3.1 Independent algebraic proof of the exact case

The exact case needed here also has a proof without analytic exponentials
or the countability argument. This is useful for checking the previous
proof independently.

Put F=C(Lambda), K=C(D), d=dim D, and pass to an irreducible component of
the geometric generic moment fiber over an algebraic closure k of K.
Its function field L has transcendence degree s-d over k and is generated
by the restrictions of z_1,...,z_s. Take a smooth projective model of L.
Let M be the integer matrix whose rows are the orders of these rational
functions along the finitely many prime divisors in their divisors. Put
m=rank M. The kernel of M consists exactly of the integer vectors a for
which z^a is in k*: a rational function with no zeros or poles on a normal
projective variety is constant. Consequently the multiplicative group
generated by the z_j modulo k* has rank m. Choosing a basis up to finite
index shows that trdeg_k L <= m. Thus

    s-d <= m.

On this geometric generic fiber, the coordinates c_j of the generic
target are constants. Taking residues in

    d psi = sum c_j dlog z_j

gives M c=0. Exact rational differentials have zero residue, also at a
prime divisor in a higher-dimensional smooth model. These are integer
linear equations on the generic target, so D is contained in ker M,
and

    d <= s-m.

Both inequalities are equalities. Since D is closed and irreducible of
dimension d, it is exactly the rational linear space ker M.

The equations z^a=constant for a in ker_Z M place the geometric generic
moment fiber in a coset of the subtorus H with Lie algebra the row space
of M. Both have dimension m=s-d, so the closed fiber component is that
coset. Its conjugates over K have the same H, which is defined by an
integer character lattice. They cover the geometric generic fiber, so
it follows by descent and closure that H preserves Lambda. Conversely every
connected horizontal stabilizer lies in the moment-fiber tangent kernel,
so it is contained in H. This again proves D=(Lie H)^perp and identifies
H as the full connected stabilizer.

This proof only needs dimension s and rational exactness of lambda on
Lambda; the Lagrangian property is a consequence of exactness and supplies
the first proof. It explains directly why algebraicity and a rational
primitive, rather than a locally chosen logarithm, are essential.

## 4. Recover the stabilizer from critical values

Let S_i be the closure of the central projection of Y_i. Since dPhi
vanishes along the fibers over S_i, the critical-value variety V_i is
generically finite over S_i. On a dense smooth open its projection has
injective differential. Define the affine conormal space over that open by

    C(V_i)={(z,t,c): c(dlog z(v))=dt(v) for every v in T_(z,t)V_i}.

It is an affine bundle of rank s-dim S_i, hence irreducible of dimension s.
The map (z,Phi,mu) from Y_i lands in it, and has image of dimension s by
section 2. Its image is therefore dense. Projection to (z,c) has closure
Lambda_i. The scalar coordinate t is its rational primitive psi, so this
projection is birational at the generic point.

If a connected subtorus preserves V_i with t fixed, it preserves C(V_i)
and hence Lambda_i. Conversely section 3 gives mu_v=0 for every vector
v tangent to Stab_T(Lambda_i)^0. Since d psi=lambda, psi is invariant under
this connected torus; consequently it preserves V_i. Thus

    Stab_T(V_i)^0=Stab_T(Lambda_i)^0.

Combining this equality with the lemma proves the theorem of section 1.
This concerns the Zariski closure of the image. It does not assert
surjectivity onto that subspace, or linearity of the exceptional target
locus, discriminant, or image of the nonsmooth part of Y.

## 5. Consequence for the adjoint Springer map

Let G be complex simple and simply connected, K a proper Levi, and use the
full restricted adjoint character Phi on the finite central cover
[K,K] x Z(K)^0. Its fundamental-character quotient is A^ell x T. Identify
central targets by c_j=B(h,H_j), using the Killing form and the chosen
central cocharacters.

The completed adjoint theorem excludes the image of the ACTUAL nonsmooth
locus of Y from a dense target open. Therefore generic nonemptiness is
equivalent to dominance of at least one clean component. Sections 1-4 give:

**Adjoint nonemptiness/excess criterion.** Generic centralizer-exact targets
have regular preimages if and only if at least one clean critical component
has a critical-value variety with finite horizontal torus stabilizer.
On each contributing component the generic group-fiber dimension is

    dim K-rank G + dim Y_i-dim Z(K)^0.

Different contributing components can have different excess. No uniform
zero-excess claim is made. Smoothness of the fibers and the exclusion of
nonsmooth components are inherited from the completed adjoint theorem.
The critical-value criterion by itself cannot remove nonreduced components.

Finite central covers, actual group fibers, regular nonsemisimple points
and central-isogeny descent are handled by the existing scheme reduction.
Connected stabilizer dimension and dominance are unchanged by the finite
central isogenies: the connected inverse image of a stabilizing torus
fixes each of the finitely many lifted components, and a stabilizer
upstairs has stabilizing image downstairs. Centralizer-exactness is imposed only after obtaining
a dense target open.

This is a structural necessary-and-sufficient criterion, not yet a closed
classification directly in terms of Dynkin diagrams. The remaining work
is to determine the relevant critical components and their stabilizers
from root-derived coefficients without repeated target-fiber calculations.

## 6. Finite computations of the symmetry obstruction

When Y_i dominates T, the irreducible critical-value variety has a primitive
irreducible equation P(z,t)=0 in C[z^+/-1,t]. Let E be the set of its CENTRAL
exponents, ignoring the t degree. Then

    Lie H_i={v: <m-m',v>=0 for all m,m' in E}.

Indeed a translation preserving the principal prime ideal sends P to a
unit multiple. A Laurent monomial unit cannot shift its finite support
nontrivially, so the multiple is a scalar and all exponents have the same
weight. Therefore dominance is equivalent to affine span(E)=Q^s.

For a rational critical value N/D, first cancel common factors and use
P=D t-N. The test uses the UNION of the supports of N and D and their
differences. Using only the numerator, or failing to cancel common factors,
can give an incorrect answer. For a Laurent value f this reduces to the
linear span of its nonzero exponents, since D=1 supplies exponent zero.

For a nonprincipal V_i one may compute a reduced Groebner basis of its
ideal, with inverse variables for the torus, then impose equality of
weights within each basis polynomial. Inverse variables have opposite
weights and t has weight zero. The uniqueness of the monic reduced basis
justifies this homogeneous-weight test. Arbitrary unprocessed generators
give only a sufficient test and must not be substituted for a reduced basis.

These are finite algebraic criteria once the clean component and its
critical-value ideal have been established. They are not claimed to make
component decomposition inexpensive or to settle every exceptional case.

There is an equivalent test that does not require constructing V_i.
On a clean component set

    M_i={a in Z^s : sum a_j mu_j=0 in C(Y_i)}.

The theorem says the full ideal of algebraic relations among the mu_j is
generated by the linear forms sum a_j c_j for a in M_i. In particular:

    mu_1,...,mu_s algebraically independent
      <=> they are linearly independent over C
      <=> no nonzero integer linear combination vanishes on Y_i.

Thus a nonlinear elimination into the target is unnecessary once a clean
prime component has been identified. Reduce the finitely many mu_j modulo
its prime ideal; collect coefficients of their normal forms and test the
rank of the resulting coefficient matrix. Full rank proves dominance.
If the rank is d<s, its kernel is the complexification of the rational
lattice M_i, and the image closure has dimension exactly d. This is a
consequence of the proved linear-image theorem, not a generic property
of polynomial maps. A prime component, its generic actual smoothness,
and exact arithmetic are essential. It does not assert that the image
of the entire critical scheme, or its nonsmooth locus, is linear.

### 6.1 A complete test on a quadratic rank stratum

For a quadratic relative character write

    Phi(q,z)=q^t H(z) q/2+L(z)^t q+Q(z),   H=H^t.

Let S be an irreducible locally closed central stratum on which H has
constant rank k and Hq=-L is compatible. Assume its solution bundle has
closure Y_i equal to an irreducible component of the full relative
critical scheme. Put b=codim_T S and a=ell-k. The solution bundle is
an affine bundle of rank a over S. In particular its reduced support
has dimension s-b+a, without a regular-sequence assumption.

At a general solution (q,z), project the mixed derivative

    B_q(v)=(d_v H)q+d_v L

to coker H. The component is clean if and only if this projected map
has rank b at the general point. Indeed the full critical Jacobian is
[H B_q], with rank k+rank(B_q mod im H), whereas the codimension of
the component is k+b. Equality of Jacobian rank and codimension at a
general smooth point of the reduced support is exactly smoothness of
the actual scheme. The rank cannot exceed b there, since the reduced
solution bundle is contained in the actual zero scheme.

For a clean component its critical value is the rational function

    psi_S=Q+L^t q_0/2

on S, where q_0 is any rational solution of Hq_0=-L. This is independent
of q_0: L is orthogonal to ker H. Thus V_i is the closure of its graph,
and the criterion becomes

    generic nonemptiness <=> finite horizontal stabilizer of graph(psi_S),
    excess = a-b.

If H is invertible, this is Q-L^t H^(-1)L/2 on a dense open of T,
and section 6 gives the numerator/denominator support test. When H is
singular, the same theorem applies on every compatible stratum with its
own clean-component test. No determinant-zero source stratum is deleted.
For a quadratic potential, the already established
[quadratic generic-smoothness theorem](QUADRATIC_BRANCHING_FIBERS.md)
supplies the exclusion of nonsmooth incidence images, even without an
adjoint hypothesis.

This criterion recovers the E6 contact excess: H has rank two among five
derived variables, compatibility imposes one smooth equation z+z^-1=0,
and the mixed map has rank one. Hence a=3, b=1 and e=2. It also explains
the SL4 failure below: the relevant coefficient vanishes to order two,
so the mixed map has rank zero instead of the required one.

### 6.2 Multiplier families and finite organization of the fibers

Suppose F_1,...,F_p cut out a smooth irreducible complete intersection S in T, and

    Phi=sum_(j=1)^p u_j F_j(z)+Q(z)

has e additional free derived coordinates. Then the actual relative
critical scheme is S x A^p x A^e and is smooth. Its critical-value
variety is graph(Q|S). Consequently generic targets have nonempty fibers
if and only if this graph has finite horizontal stabilizer. Their
transverse fibers are finite disjoint unions of copies of A^e. Indeed
the affine-conormal map from S x A^p has dimension s and is generically
finite etale when dominant; the extra coordinates remain free. This is
a structural positive-excess theorem, including the E6 contact control.
For a disconnected smooth S apply the criterion separately to each
irreducible component. A finite stabilizer of the UNION is insufficient:
for Phi=u(z-1)(w-1), the two critical-value components are {z=1,t=0}
and {w=1,t=0}. Their union has finite horizontal stabilizer, but their
target images are the two coordinate axes, so a generic target is empty.

When Q|S is constant, this affine conormal is the ordinary conormal of S.
For smooth CLOSED S in T, its generic moment degree is the published
Gaussian degree

    gdeg(S)=(-1)^(dim S) chi(S,C).

This is [Franecki--Kapranov, Corollary 1.5 and Proposition 2.2](https://arxiv.org/pdf/math/9909088).
Thus in this multiplier family the generic transverse fiber has exactly
gdeg(S) copies of A^e. For S={z+w=1} in (C*)^2, chi(S)=-1, agreeing
with the explicit one-solution conormal calculation. For the E6 contact
cover S consists of two points, giving two affine planes before the
known central identification. This count is a direct application of the
published theorem. It does not extend unchanged to nonconstant Q, singular
S, or actual group-fiber component counts without checking the central
cover and the remaining smooth factors.

More generally, for a dominating clean component in the adjoint setting,
the factorization

    Y_i -> Lambda_i -> A^s

becomes a smooth surjection followed by a finite etale surjection after
shrinking the target to a dense open. For the second arrow, generic
finiteness and affineness make it finite after a target localization;
characteristic zero then removes its branch locus. The image of the
first arrow contains a dense open of Lambda_i. Its complement has a
proper closed image in A^s under the finite second arrow, so another
shrinking makes the first arrow surjective. Apply generic flatness to
Y_i -> A^s and shrink once more. The completed adjoint theorem supplies
smooth geometric fibers there; flatness and finite presentation then
make Y_i -> A^s smooth. An etale intermediate base makes Y_i -> Lambda_i
smooth as well. Smooth fibers alone were not used to assert flatness.

The finite cover therefore records the finitely many possible central
source positions and critical values; the smooth first map records excess
dimension. It need not be a product, and its fibers need not be connected
or affine. The degree of Lambda_i -> A^s is NOT generally the degree of
the full Springer map or the number of group-fiber components. For K=T,
however, Y=T and Lambda is the graph of the torus differential; its degree
is exactly the torus Springer degree, hence the full generic degree by
Kostant--Michor. No new numerical degree is computed by this statement.

### 6.3 A root-only sufficient test: a Morse zero of the torus character

There is also a useful sufficient test that avoids critical-component
decomposition. Let T_G be a maximal torus containing the central torus of
K. Suppose t in T_G satisfies

    theta_Ad(t)=0,
    alpha(t)!=1 for every root alpha of K,
    H_t=sum_(alpha in R(G)) alpha(t) alpha tensor alpha is nondegenerate.

Then generic centralizer-exact targets for K have nonempty regular fibers,
and at least one contributing component has zero transverse excess.

The second condition says t is regular semisimple IN K, not necessarily
in G. Thus T_G -> K//K is etale at t. The first condition makes the full
quotient character critical there. At a critical point the Hessian
changes by congruence under an etale coordinate change, so the full
Hessian of Phi in derived and logarithmic central coordinates is
nondegenerate. Its first ell rows make Y smooth of dimension s. Its
remaining s rows give an isomorphism on the tangent kernel of those
first rows. Hence mu:Y -> A^s is etale at this point, even if the
derived Hessian block alone is singular. This proves dominance of a
zero-excess component. Intersect its dense target image with Z_K^o and
use the regular quotient lifts; over this open K-regularity becomes
G-regularity. The zero target used to certify dominance is not itself
being asserted centralizer-exact.

The test is entirely in root data. Two useful finite-order reductions are:

* If Ad(t)^2=1, opposite roots cancel in the first derivative, so
  theta_Ad(t)=0 automatically. With F=Z_G(t)^0,

      H_t=2 B_F|_(Lie T_G)-B_G|_(Lie T_G).

  Here B_F is its intrinsic Killing form, zero on its center. The formula
  follows by separating the +1 and -1 root spaces, so it also works when
  F is reductive rather than semisimple.
* If Ad(t)^3=1 and F is semisimple of full rank, equivariance puts
  theta_Ad(t) in z(Lie F)=0. Opposite nontrivial root phases add to -1,
  giving H_t=(3 B_F-B_G)/2 on the torus. Nondegeneracy remains a required
  condition, not an automatic consequence of finite order.

For a simply laced root subsystem factor F_j, the first test reduces to
2 h^vee(F_j)!=h^vee(G), and the second to
3 h^vee(F_j)!=h^vee(G), since the root-length normalizations agree.
Central directions of F in the involution case contribute -B_G.
In unequal root lengths use the bilinear forms themselves, including the
embedding indices, rather than these simplified dual-Coxeter formulas.

This is a sufficient existence test, not an equivalence. A smooth excess
family such as the E6 contact example cannot be proved by an etale
zero-excess witness. The critical-value criterion remains the broader
statement and does not assume a Morse zero exists.

### 6.4 The type A classification recovered by the Morse test

The existing all-rank type A answer has a short structural explanation.
Let G=SL_N and K have defining-space block sizes b_1,...,b_m. For a proper
Levi, generic regular-source fibers are nonempty exactly when every
b_i<=2 and at least one block is a singleton. Their dimension is twice
the number of two-blocks; this classification is already proved in
[the classical note](CLASSICAL_ADJOINT_FIBERS.md).

For sufficiency put one + and one - eigenvalue in each two-block and put
all singleton eigenvalues in the + group. Multiply every eigenvalue by
a common scalar to obtain determinant one. If the total multiplicities
are p and q, then p-q is the number of singletons, so p!=q. The element
is K-regular semisimple and Ad(t)^2=1. Its centralizer is
S(GL_p x GL_q). On the sl_p and sl_q torus summands, 2B_F-B_G has
respective factors proportional to p-q and q-p; on the remaining
one-dimensional center it is -B_G. It is therefore nondegenerate.
The Morse test proves nonemptiness without solving the central target
equations. The known quadratic theorem gives smoothness of all generic
regular-source fibers and the existing dimension calculation gives the
stated pure dimension.

For necessity write the character as Phi=UV-1, where U and V are the
defining and dual traces. A block of size at least three has independent
first and last fundamental characters, so its derived derivatives force
U=V=0. All central derivatives then vanish. If all blocks have size two,
write U=sum z_i x_i and V=sum z_i^-1 x_i. The derived equations are
z_i V+z_i^-1 U=0. Unless U=V=0, they force every z_i^2 to have a common
value a, whence V=U/a and the same equations give 2U=0, a contradiction.
Again all central targets are zero. These cases therefore cannot meet
the centralizer-exact target open for a proper Levi. No reduction of a
nonzero-target fiber is needed: having no complex points makes its
finite-type scheme empty.

The balanced involution p=q>=2 has torus Hessian kernel of dimension
N-2, matching the failure of this Morse mechanism in the all-two-block
case. This is an explanation and independent proof of an existing
classification, not an additional rank table or a new priority claim.
In the integral torus basis diag(e_i-e_N), the exact determinant is

    det H_t=(-1)^q 2^(N-1) N^2 (p-q)^(N-2),

with the exponent-zero convention when N=2. Splitting the trace form
as above gives this formula: det B_SL_N=2^(N-1)N^N, and the relative
eigenvalues are (p-q)/N, (q-p)/N and -1 on the respective summands.
Direct root sums check the formula independently for every p+q=N,
2<=N<=12. The type A proof itself has no rank bound.

## 7. Required controls and root-character applications

The E6 contact family has, after a unit pair elimination,

    Phi=z^2+z^-2+(z+z^-1)u,

with two further free derived coordinates. Its two clean components have
z=+/-i and critical value -2. Each critical-value variety is a point, so
its horizontal stabilizer is trivial. The dimension formula gives excess
two, exactly the known nonempty smooth 32-dimensional group fiber. The
two sheets are identified by the actual central cover, as already proved.

For SL4 with Levi S(GL2 x GL2), use central parameter z and derived
characters x,y. The full character is

    Phi=x^2+y^2+(z^2+z^-2)xy-1.

With u=x+y, v=x-y this is

    Phi=(z+z^-1)^2 u^2/4-(z-z^-1)^2 v^2/4-1.

The component supported on z=1,u=0, with v generically nonzero, has
local critical equation (z-1)^2=0. It is not clean. Its reduced
critical-value variety is a point and has trivial stabilizer, but its
central target is supported at zero. Thus dropping the clean-component
hypothesis would falsely predict generic nonemptiness even for an actual
adjoint character. The clean full-torus branch x=y=0 has value -1 and
stabilizer T, correctly giving only target zero. The analogous walls
must be kept rather than silently reduced.

The Sp4 special-target singularity remains compatible with this generic
criterion. The altered-character persistent-nilpotent controls show why
the adjoint nonsmooth-image theorem is a substantive input.

The final E8 A2+A1^3 model has an exact dominance witness
(x,y,w,a,U,V)=(0,0,0,0,1,1), with a=s+s^-1 and s=i. The relative Hessian
is the 3-by-3 matrix with diagonal 2 and off-diagonal -4, of determinant
-216. Its central logarithmic Schur complement is
diag(-40,[[-16,-8],[-8,-16]]), of determinant -7680.

This has also been checked in the FULL, uneliminated character, rebuilt
independently from all 240 E8 roots. For retained nodes 12357, use the
archived variables (x1,x2,x3,x5,x7,a0,z0,z1), replace a0=s+s^-1, and
take (x1,x2,x3,x5,x7,s,z0,z1)=(-3,0,-3,0,0,i,1,1). All five derived
and three logarithmic central first derivatives vanish. The derived
Hessian determinant is 216, the central Schur determinant is -69120,
and the full Hessian determinant is -14929920. The two central monomial
coordinates U=z0^3/z1 and V=z0^6/z1 have exponent determinant three,
explaining the factor nine between the two Schur determinants.
Native SymPy and Sage over Q(i) independently obtain the original
character determinants.

The full Hessian makes the relative central map etale at this point,
hence dominant. Its target is zero, but only dominance is certified
there; the generic targets are then intersected with Z_K^o. The prior
all-boundary proof gives smoothness and pure dimension twelve for every
nonempty generic regular fiber. They are now proved NONEMPTY as well.
All eight certified Weyl-equivalent standard subsets inherit this result.
No generic map degree or number of orbit components is asserted.

### 7.1 A root-theoretic application in every simply laced rank

Let G have an irreducible simply laced root system of rank r>=2, and let
K be a Levi with derived group A1. Then generic centralizer-exact targets
have nonempty smooth regular-source fibers of pure dimension two.

Here is an all-rank argument, rather than an inference from checked ranks.
Normalize roots to have squared norm two and let alpha be the retained
root. On the central cover write x for the SL2 fundamental character.
Every other root has alpha-coroot pairing 0,+1 or -1. The full character
therefore has the form

    Phi=x^2+(r-2)+sum_(beta orthogonal alpha) z^beta + x L(z),
    L=sum_(<beta,alpha^vee>=1) z^(beta-alpha/2).

The projected weights are characters of the actual central cover; an
additional finite cover may be used to write coordinates, without changing
dominance. Every weight mu=beta-alpha/2 in L has squared norm 3/2.
Eliminating x=-L/2 is global and scheme-theoretic, since Phi_xx=2.
The critical value is the Laurent polynomial

    psi=(r-2)+sum_(beta orthogonal alpha) z^beta-L^2/4.

For each mu the coefficient of z^(2mu) is -1/4. Indeed 2mu has squared
norm six and cannot be an orthogonal root. If mu'+mu''=2mu, equality in
the strict convexity inequality for the common norm forces
mu'=mu''=mu. Thus these terms cannot cancel.

These mu span alpha-perp. To see this, let E be the span of all roots
not orthogonal to alpha, including alpha. Reflection in a root belonging
to E preserves E. Reflection in a root orthogonal to alpha permutes the
roots not orthogonal to alpha and also preserves E. Hence E is a nonzero
W-invariant subspace, and irreducibility gives E equal to the whole root
space. Projection onto alpha-perp proves the claim; pairing -1 gives
the negatives of the weights from pairing +1.

Consequently the support of psi spans all r-1 central directions. The
Laurent support criterion, already covered by Fassarella's theorem, makes
its logarithmic gradient dominant. Y is the graph x=-L/2 over T, so its
excess is zero. The group orbit contribution is dim SL2-rank SL2=2.
Generic smoothness follows directly in this graph situation, or from the
completed adjoint theorem. This proves the assertion.

The proof uses equal root lengths in the noncancellation step. It does
not silently extend to short-root A1 Levis in other types. The finite
root checks in the accompanying audit corroborate the identity and span
claims for A2-A8, D4-D8 and E6-E8; they are not the all-rank proof.

### 7.2 Orthogonal A1 Levi factors in E7 and E8, without a Levi census

For G of type E7 or E8 and every Levi whose derived type is A1^k, generic
centralizer-exact targets have nonempty smooth regular fibers, with at
least one zero-excess component of dimension 2k.

The proof applies the root-only Morse test uniformly. Any involution of
the adjoint torus has H_t=2B_F-B_G. Since the ambient system is simply
laced, this form could be singular only if a simple root subsystem
factor F_j had h^vee(F_j)=h^vee(G)/2. In E7 the latter is nine: the only
simply laced possibility would be A8, whose rank exceeds seven. In E8
it is fifteen: the only possibility would be A14, whose rank exceeds
eight. D-type dual Coxeter numbers are even, and the exceptional values
are 12,18,30. Therefore EVERY adjoint-torus involution has nondegenerate
torus Hessian in E7 and E8. This argument does not classify involution
conjugacy classes or enumerate Levis.

Conjugate K to a standard Levi. Assign value -1 to every retained simple
root and +1 to the others; these assignments define a point of the
adjoint torus and lift through the central isogeny to G. Since the
retained roots are pairwise orthogonal, the roots of K are just their
positive and negative simple roots. The lift is therefore K-regular
semisimple, and all three Morse-test hypotheses hold.

For the seven-subset E8 A1^4 family, the previously completed boundary
proof already gives pure dimension eight if nonempty. The present
root-theoretic argument proves nonemptiness, so these generic regular
fibers are smooth, nonempty and pure of dimension eight. This also
recovers the previous E7 A1^4 and E8 A1^3 nonemptiness results by a
common mechanism. It does not establish that every other exceptional
Levi admits an involution witness, nor that every component always has
zero excess. E6 demonstrates the limit: h^vee(E6)/2=6 is attained by
A5, so its analogous universal involution claim is false.

## 8. Published inputs and the scope of the audit

Fassarella's [Logarithmic Hesse's Problem](https://arxiv.org/abs/1102.1659)
already proves the Laurent-polynomial logarithmic Hessian rank statement.
The corresponding support test in section 6 is therefore not claimed as
new. Franecki--Kapranov's
[Gauss map and noncompact Riemann--Roch paper](https://arxiv.org/abs/math/9909088)
provides the Gaussian-degree framework for Lagrangians in torus cotangent
bundles. Their terminology helps distinguish finite generic fibers from
generic nonemptiness. Generic smoothness, finite central descent and the
regular adjoint quotient remain the standard geometric inputs.

The exact Lagrangian reduction, stabilizer criterion and root-character
applications are written deductions here. The proof of the exact torus
lemma has been checked by two routes: symplectic orthogonality and an
independent valuation/residue argument on a projective model. The latter
uses only torus coordinates and rational exactness; no analogous assertion
about arbitrary symplectic manifolds is needed. No claim that this lemma
or the assembled criterion is unpublished is made. The literature search
checked the above primary sources and Kostant--Michor, and did not establish
priority for the assembled theorem.

The [native exact audit](critical-value-symmetry-audit.json) and
[Sage verification](critical-value-symmetry-verification.json) have separate
implementations. They check rational support/Hessian identities, conormal
inverses, the actual SL4 nonreduced saturation, the E6 excess and composition
controls, root-derived Hessians, and the original uneliminated E8 witness.
Their source bindings make the checked versions explicit. The native E8
character is rebuilt from all 240 roots; Sage differentiates the original
character over Q(i). Finite checks corroborate the all-rank proofs and do
not replace them. The earlier E8 homogeneous membership certificate still
has its documented one-engine dependency; it is not rerun or independently
reimplemented by these new checks.

The remaining mathematical problem is a root-theoretic classification of
the contributing critical components, stabilizers and excess dimensions.
The theorem reduces that problem to symmetry and component geometry; it
does not yet solve it for every exceptional Levi. Nor does it prove
special-target smoothness, smoothness outside the regular source, generic
map degrees, irreducibility, or global harmonic separation.

## 9. Zero targets and the extension to semisimple groups

The proper-Levi restriction can now be removed after a short separate
zero-target calculation. This uses the existing
[regular zero-fiber theorem](SIMPLE_GROUP_CAYLEY_BRANCHING.md#1-scope-and-the-general-zero-fiber-theorem),
which identifies the actual regular zero fiber with the smooth pullback
of the critical scheme of the adjoint character in fundamental characters.

For simple simply connected G, these character polynomials are:

| Simple type | Adjoint character | Critical scheme |
| --- | --- | --- |
| A1 | x1^2-1 | x1=0, reduced |
| A_r, r>=2 | x1*x_r-1 | x1=x_r=0, an affine space of dimension r-2 |
| C_r, r>=2 | x1^2-x2-1 | Empty, since the x2 derivative is -1 |
| B_r (r>=3), D_r (r>=4), G2, F4, E6, E7, E8 | A fundamental character | Empty, since its derivative is 1 |

Low-rank isomorphisms handle the duplicated types. The A identity follows
from End(V)=1+Ad. In C, Ad=Sym^2 V and wedge^2 V=1+V(omega_2), giving
the displayed polynomial. In the final row the highest root is a
fundamental weight; this is standard root data. Independent Weyl-character
checks are supporting controls, not the all-rank argument.

Thus the regular zero fiber is ALWAYS smooth. It is nonempty exactly in
type A. Its dimension is two for A1 and dim G-2=(r+1)^2-3 for A_r, r>=2.
These assertions include actual schemes and regular nonsemisimple points.
They do not say the entire zero fiber is smooth or empty in other types.

**Semisimple extension.** Let G be any connected complex SEMISIMPLE
algebraic group, let theta_Ad use its nondegenerate Killing form, and let
K be ANY Levi subgroup, including K=G. There is a dense open U_K in
Z_K^o on which every actual regular-source fiber is smooth, allowing
emptiness and excess.

Proof. A finite central etale cover of G is a product of simple simply
connected groups. The Killing form and the adjoint Springer map split
factor by factor. The inverse image of K is a product of Levis K_j.
For every proper K_j use the completed simple-group theorem. A factor
with K_j=G_j contributes the zero target and the smooth regular zero
fiber just calculated. Their product is smooth over the product of the
good target opens, and finite central etale descent gives the assertion.
Regularity is factorwise and is preserved by the central isogeny.

Generic nonemptiness requires nonemptiness in every proper-Levi factor
and type A for every full simple-factor Levi. In the latter factors the
dimensions above are added. The proof does not require a new exceptional
case calculation. It also closes the full-factor qualification left in
the preceding [central-isogeny audit](E8_DEGREE_TRANSFER_AUDIT.md).

No extension to reductive groups with positive-dimensional center is
asserted with this normalization: their Killing form is degenerate, so
the definition would first require an additional choice of trace form.
