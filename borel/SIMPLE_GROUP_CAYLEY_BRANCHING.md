# Simple-group Cayley fibers: character coordinates and Levi obstructions

Date: 2026-10-23. Written characteristic-zero deductions, with exact finite
corroboration. These are consequences of established character and quotient
theory, not a claim of new literature priority. No integral Borel-image,
generic-degree, or regularity-preservation classification is asserted here.

## 1. Scope and the general zero-fiber theorem

Let G be a connected simply connected complex simple group of rank r.
Let rho be a nontrivial actual finite-dimensional representation, chi its
character, and theta its Springer/generalized Cayley map, normalized by

    B_rho(theta(g),X) = tr(rho(g) d rho(X)),
    B_rho(X,Y) = tr(d rho(X) d rho(Y)).

This trace form is nondegenerate. Write x_i=chi_(omega_i) for the r
fundamental characters. They give the quotient map

    pi:G -> G//G = A^r,        chi=P(x_1,...,x_r).

An element is regular if its centralizer has dimension r; this does not
mean it is semisimple. All assertions below concern complex group maps,
not parabolic cohomology or the minimum degree among rational Cayley maps.

**Theorem 1 (the regular part of the zero fiber).** As schemes on G_reg,

    theta^(-1)(0) intersect G_reg = pi^(-1)(Crit(P)) intersect G_reg.   (1)

Every quotient point has a regular lift. Consequently the zero fiber has
a regular point if and only if the polynomial P has a critical point in
A^r. At such a point its completed local ring is

    Jac(P at pi(g)) [[u_1,...,u_(dim G-r)]].                          (2)

In particular, if rho is fundamental, its zero fiber has no regular
elements. The same conclusion holds for a nonzero direct sum of fundamental
representations, with arbitrary nonnegative multiplicities and trivial
summands: P is then a nonconstant affine linear polynomial.
The underlying zero fiber has dimension at most dim G-3.

For any nontrivial irreducible rho, every central element z is nevertheless
a reduced isolated point of theta^(-1)(0). Thus the fundamental zero-fiber
statement is an exclusion of regular inputs, not an empty-fiber assertion.

**Proof.** Trace differentiation gives

    d chi_g(gX)=B_rho(theta(g),X).

The Maurer--Cartan frame and B_rho identify the fiber ideal with the ideal
of d chi. On G_reg, Steinberg's fundamental-character differentials are
independent. The chain rule therefore identifies this ideal with the
pullback of (P_1,...,P_r), giving (1); smooth local coordinates give (2).
Steinberg's regular cross-section gives a regular lift of every quotient
point. For an affine linear P the gradient is a nonzero constant vector.
The dimension bound is Steinberg's codimension-three theorem for the
nonregular locus; it is a bound on the underlying scheme's dimension.

For irreducible rho, Schur's lemma gives rho(z)=lambda I, lambda nonzero.
Since tr d rho(X)=0 on the simple Lie algebra, theta(z)=0. Differentiating
the defining trace identity at z, in left-translated coordinates, gives

    B_rho(d theta_z(zX),Y)=lambda B_rho(X,Y).

Thus d theta_z is invertible. The fiber is reduced and zero-dimensional
in a neighborhood of z. This last argument does not require rho fundamental.

These conclusions descend through any central isogeny on which rho is an
actual representation. A central isogeny is etale in characteristic zero,
preserves regularity, and intertwines the normalized trace maps via its
Lie algebra isomorphism. Fundamental representations that do not descend
are not representations of the quotient group and are not included there.

## 2. A branching criterion for arbitrary semisimple targets

Fix a semisimple h in Lie(G), and put K=Z_G(h), Q=K//K. The earlier
[critical-potential theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md) proves that
the entire scheme theta^(-1)(h) lies in K. At its regular points it is
controlled by the quotient potential

    W_h = chi|_Q - sum c_j log eta_j.                                (3)

Here eta_j are a basis of characters of K and
B_rho(h,X)=sum c_j d eta_j(X). Logs are local formal series; their
differentials are regular. In particular the target correction cannot
be dropped. Pass if necessary to the finite etale central cover
Z(K)^0 times the simply connected simple factors of [K,K]. Its quotient
coordinates are central torus coordinates and fundamental characters
of these factors. Character logs depend only on central coordinates.

**Theorem 2 (branching and elimination).** Suppose that in these local
quotient coordinates, x is a derived-factor coordinate and

    chi|_Q = A(y) x + B(y).

If A is nonzero at a given quotient point, no regular point above it
belongs to theta^(-1)(h). If A is an everywhere nonzero constant, the
whole fiber has no regular points.

More generally suppose that x=(x_1,...,x_p) are derived-factor
coordinates, y contains the other r-p coordinates, and

    chi|_Q = sum x_i F_i(y) + R(y).

At a regular fiber point, F(y)=0 necessarily. If dF has rank p there,
the completed transverse algebra of the fiber is the Jacobian algebra
of the restriction to F=0 of

    R(y) - sum c_j log eta_j(y).                                    (4)

There are r-2p residual variables. In particular r=2p gives a reduced
smooth fiber germ of dimension dim K-r; r=2p+1 leaves a single-variable
critical equation. If that critical point is isolated, its derivative
order is the generic multiplicity of the corresponding quotient-fiber
component. If the derivative vanishes identically, the critical point
is not isolated and no finite multiplicity claim follows.

**Proof.** The derivative of (3) in x is A or F, since all eta_j are
independent of x. This proves the exclusion. For the second assertion
use v=F(y) as part of the y coordinates. The equations in x impose
v=0, the equations in v solve uniquely for x, and the remaining equations
are the derivatives of (4) on v=0. This is the exact elimination of
Lemma C of the general note. Smoothness of the quotient at regular fiber
points is the other input. Translating x to its value at the point adds
a multiple of F to R, which restricts to zero and has no effect.

This is a sufficient branching criterion, not a classification of all
irreducible characters. It applies without a type restriction, including
exceptional groups. An actual decomposition must still establish its
linearity and rank hypotheses; they are not automatic for fundamental rho.

## 3. Uniform restrictions for classical fundamental representations

The following are necessary conditions for the fiber over a semisimple h
to contain a regular group element. They do not say the whole fiber is
empty when the condition fails, or that the conditions are sufficient.

**Theorem 3 (target block bounds).**

| Group and irreducible representation | Necessary target condition |
| --- | --- |
| SL_N, rho=exterior^k, 1<=k<N | Every eigenvalue multiplicity of h is at most min(k,N-k). |
| Sp_(2n), rho=V_(omega_k), 1<=k<=n | If dim ker(h)=2m, then m<k. |
| Spin_(2n+1), rho=exterior^k of the vector representation, 1<=k<=n-1 | If dim ker(h)=2m+1, then m<=k. |
| Spin_(2n), rho=exterior^k of the vector representation, 1<=k<=n-2, n>=4 | If dim ker(h)=2m, then m<=k+1. |

The B and D representations are pulled back from SO. The same bounds
therefore hold on SO and on every central quotient carrying the specified
representation. Half-spin and spin representations are not included in
the exterior ranges of the last two rows. The all-type fundamental
zero-fiber conclusion of Theorem 1 still applies to them on Spin.

**Proof in type A.** An eigenvalue block of size m gives a GL_m factor
of K, subject to the product-one determinant constraint among blocks.
Let a_j be its elementary characteristic coefficients. For j<m these
are independent derived quotient coordinates when determinants are fixed.
The exterior character is the coefficient of z^k in the product of all
block characteristic polynomials. If m>k, it has the form

    chi_k = a_k + expression independent of a_k.

Its derivative in a_k is 1; the determinant constraint and all target
logs are independent of a_k. Theorem 2 excludes regular fiber points.
For the dual representation, B_(rho*)=B_rho and

    theta_(rho*)(g)=-theta_rho(g^(-1)).                              (5)

Since exterior^(N-k) is dual to exterior^k, applying the same argument
to -h gives m<=N-k as well. This proves the first row.

**Proof in type C.** The zero eigenspace gives an Sp_(2m) factor. Write
A(z) for its reciprocal characteristic polynomial and E(z) for the
remaining paired blocks, with E(0)=1. If m>=k, the coefficient a_k of
A(z) is an independent quotient coordinate, and

    chi_(omega_k)=Coeff_(z^k) ((1-z^2) A(z) E(z)) = a_k + other terms. (6)

Indeed a_k first occurs at degree k; its reciprocal occurrence is at
degree 2m-k>=k, with only one occurrence if m=k. All other coefficients
of degrees <=k are independent of a_k. The primitive correction
-[z^(k-2)]A E does not involve a_k. Equivalently a_k differs from the
k-th fundamental character of Sp_(2m) by lower coefficients. The target
logs involve no Sp coordinate. Theorem 2 proves m<k.

**Proof in types B and D.** On the zero block the vector characteristic
polynomial A(z) has degree 2m+1 in B and 2m in D. Its coefficient a_k is
the k-th exterior character of that vector block. In B it is a fundamental
quotient coordinate when k<=m-1; in D when k<=m-2. On the central cover
these coordinates are independent of the spin/half-spin coordinates.
The character of the ambient exterior representation is [z^k]A(z)E(z),
whose derivative in a_k is 1, with no earlier coefficient involving a_k.
Thus m>=k+1 in B and m>=k+2 in D are excluded. The range endpoints are
essential: the next exterior characters involve spin coordinates
nonlinearly, so the same proof does not extend by one more step.

In particular, for the defining SL_N map, a semisimple image of a regular
element must be regular semisimple. For the defining Sp_(2n) map, a
semisimple image of a regular element must be invertible. These statements
do not assert preservation of regularity in Sp; nonzero eigenvalue
collisions can still occur.

## 4. A whole target stratum for middle exterior powers

**Theorem 4 (uniform component and multiplicity formula).** Let k>=2,
G=SL_(2k), rho=exterior^k, and I=binomial(2k-2,k-1), so
B_rho(X,Y)=I tr(XY). Let

    h=diag(lambda I_k, h_1,...,h_k),
    k lambda + sum h_j=0,

where all h_j differ from lambda. They may repeat: let r_1,...,r_s
be the multiplicities of their distinct values, so sum r_j=k. Thus
K=Z_G(h)=S(GL_k x product_j GL_(r_j)). Set M=k!/(product_j r_j!).
The irreducible components of the full fiber theta^(-1)(h) that meet
G_reg have dimension k^2-k+sum_j r_j(r_j-1). Their number and generic
scheme multiplicity are exactly

| Target condition | Number of components meeting G_reg | Multiplicity of each |
| --- | --- | --- |
| 1+I^2 lambda^2 != 0 | 2 M | 1 |
| 1+I^2 lambda^2 = 0 | M | 2 |

On the regular locus their completed local rings are a formal smooth factor
of the stated dimension times C, respectively C[s]/(s^2). Thus the
multiplicity-weighted number is always 2 M. This is a count of positive-dimensional special
fiber components, not the generic degree of the Springer map. Components
wholly contained in G minus G_reg are not classified by this theorem.

**Proof.** First let the remaining target eigenvalues be distinct, so
g in K has the form (A,b_1,...,b_k). Let a_j=e_j(A), D=product b_i,
and b_j^el=e_j(b_1,...,b_k). The determinant constraint gives a_k=D^(-1).
The quotient has coordinates a_1,...,a_(k-1) and the nonzero labelled
eigenvalues b_i. Its corrected character potential is

    W = sum_(j=1)^(k-1) a_j b_(k-j)^el + D^(-1)+D
        - sum_i I(h_i-lambda) log b_i.                              (8)

The a derivatives impose b_1^el=...=b_(k-1)^el=0. Consequently

    product_i (u-b_i)=u^k+(-1)^k D.

Since D is nonzero, the b_i are automatically distinct. The full map
from ordered b_i to elementary coefficients is etale there (Vandermonde
Jacobian). Thus the first k-1 elementary coefficients have independent
differentials. Theorem 2 eliminates the a_j and these k-1 coefficients.
On each remaining local curve, b_i=s beta_i and D=s^k D_0. With
C=sum I(h_i-lambda)=-2k I lambda, the residual potential is

    w(s)=D_0^(-1)s^(-k)+D_0 s^k-C log s + constant.

Its logarithmic derivative vanishes precisely when

    D-D^(-1)+2I lambda=0,
    or equivalently D^2+2I lambda D-1=0.                            (9)

Both roots of (9) are nonzero. There are two simple roots unless
1+I^2 lambda^2=0, when there is one double root. D is a local parameter
on the curve (dD/ds=kD/s is nonzero), so these are exactly the transverse
Jacobian lengths 1 and 2. For each D, the polynomial in u has k distinct
roots; assigning them to the k labelled target eigenspaces gives k!
quotient points. For each assignment the other critical equations solve
uniquely for the k-1 coordinates a_j, by the same full-rank mixed pairing.

Now allow the h_i to repeat. On the full K quotient, use the characteristic
coefficients of its remaining GL_(r_j) blocks. Before ordering roots, the
a derivatives still impose that the characteristic polynomial of their
combined k-dimensional block is u^k+(-1)^k D. Its roots are necessarily
distinct, so every possible critical point lies in the open set where
ordering these roots is an etale coordinate cover. Formula (8) applies
on that cover. Its logarithmic terms are invariant under permutations
within each target block, because the coefficients I(h_i-lambda) agree
there. The product of symmetric groups acts freely on the ordered-root
points; its order is product r_j!. Consequently the count k! for each
D becomes M, with the local lengths unchanged. There are no other
critical points hidden outside the ordered-root chart.

Every such quotient point has a regular K lift, and the critical equations
imply theta(g)=h there. At a fiber point Z_G(g) is contained in K, hence
regularity in K and G agrees. Each K quotient fiber is irreducible of
dimension dim K-r=k^2-k+sum r_j(r_j-1) and has a dense regular orbit. The isolated
critical point therefore makes its closure a full-fiber component, with
the computed multiplicity, by Corollary B of the general note. Conversely
every regular fiber point satisfies these equations, so the list is
exhaustive among components meeting G_reg.

This removes the special root-of-unity target choice from the earlier
middle-exterior calculation and includes all collisions among the remaining
target eigenvalues. In its odd-k example, r_j are pairs and one singleton,
so the dimension is k^2-1, agreeing with that earlier proof. If an h_i
equals lambda, that eigenvalue has multiplicity greater than k and
Theorem 3 excludes every regular point in the fiber. The exceptional condition
here concerns the repeated target eigenvalue lambda, not the source
eigenvalues or a representation-mixture parameter.

## 5. An additional irreducible type A family

**Corollary 5 (hook Schur representations).** Let lambda=(a,1^b), a>=1,
b>=0, d=a+b<N, and rho=S_lambda(C^N) on SL_N. If a semisimple target h
has an eigenvalue block larger than d, its fiber has no regular points.
In particular the zero fiber has no regular points. The same holds for
the dual representation (with target sign reversed).

**Proof.** In the elementary basis the homogeneous Schur function has

    s_(a,1^b) = (-1)^(a-1) e_d + Q(e_1,...,e_(d-1)).                (10)

For completeness the hook Jacobi--Trudi expansion is
sum_(j=0)^b (-1)^j h_(a+j)e_(b-j). The only term involving e_d is
(-1)^b h_d, and the generating identity H(z)E(-z)=1 makes its
coefficient (-1)^b(-1)^(d-1)=(-1)^(a-1). In a block of size m>d,
e_d of the full alphabet contains a_d with coefficient 1, and lower
e_j do not contain it. Thus the derivative is the nonzero constant
(-1)^(a-1), and Theorem 2 applies. The strict inequality d<N ensures
this is a free quotient coordinate when h=0; at d=N it is a fixed
determinant and the argument fails. No arbitrary-Schur classification
is claimed. This includes Sym^d for d<N, beyond fundamental modules.

## 6. Inputs, corroboration, and the next question

Published inputs are the trace construction and differentiation in
[Kostant--Michor, Section 2](https://arxiv.org/pdf/math/0109066),
fundamental-character quotient coordinates, regular differential
independence and codimension three in
[Steinberg, Theorems 1.3/1.5 and Section 6](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf),
classical fundamental/exterior constructions (including primitive
symplectic powers) in
[Fulton--Harris, Chapters 15--20](https://www.cimat.mx/~gil/docencia/2008/representaciones/fulton-harris.pdf),
also presented in [Etingof, Lecture 31](https://ocw.mit.edu/courses/18-755-lie-groups-and-lie-algebras-ii-spring-2024/mit18_755_s24_lec05.pdf),
and the symmetric-function identities in
[Macdonald, Chapter I, Sections 2--3](https://math.berkeley.edu/~corteel/MATH249/macdonald.pdf).
The proofs above apply these ingredients to the corrected character
potential; novelty has not been established. In particular Theorem 1 is
a short direct corollary of Steinberg and trace differentiation.

[Exact checks](simple-group-cayley-branching-verification.json) compare
actual Weyl-character weight multiplicities with exterior/primitive
formulas, verify symbolic block-coordinate derivatives and hook
coefficients, and test the necessary endpoint caveats in rank one.
Original exterior-power matrices in SL_4 and SL_6 check both the generic
and double-component cases of Theorem 4 for every partition of the remaining
target multiplicities in these ranks, against every Lie algebra basis
direction, including the normalization I. Symbolic residual polynomials
check the quadratic and its multiplicities in further ranks.
These finite tests corroborate the written proofs; they do not prove
all-rank statements or classify exceptional Levi branchings.

Theorem 4 settles all semisimple-target strata having an eigenvalue of
multiplicity k for exterior^k of SL_(2k). The next goal is a
necessary-and-sufficient criterion on the analogous symplectic boundary
strata, or for target multiplicities strictly below k in this type A family.
There the constant
derivative disappears: its replacement is an explicit coefficient
function on the other blocks. Determining whether that function can
vanish, and whether its derivative has full rank, will distinguish
actual collapsed components from excluded strata and compute their
transverse multiplicities by Theorem 2. The same branching test applies to
other simple types. This is also a precise test to
apply to exceptional Levi branchings, rather than extending a table of
individual group elements. The nonregular part of a fiber remains a
separate problem governed by the critical scheme of the pulled-back
potential, not (1) outside G_reg.
