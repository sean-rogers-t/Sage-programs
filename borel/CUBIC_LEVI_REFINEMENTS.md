# Cubic Levi refinements: product transfer and a Fricke boundary argument

Date: 2026-10-23. Written characteristic-zero deductions with independent
exact checks. No literature-priority claim. This closes the two mechanism
tests left in GEO-53/54, not every exceptional Levi.

## 1. Statements, conventions, and the common character

Let G be complex simply connected E6, rho=Ad, and use

    B_Ad(theta(g),X)=tr(Ad(g) ad(X)).                         (1)

For a Levi K and a semisimple target h with Z_G(h)=K, write
F_h^reg=theta^(-1)(h) intersect G_reg, with its actual scheme structure.
Regular means centralizer dimension six. Nothing here describes all
nonregular points of the full fiber or changes the integral Borel work.

**Theorem 1.** There is a dense open set of central targets with
centralizer K for which F_h^reg is nonempty and smooth of pure dimension:

| Retained Bourbaki nodes | Derived K | dim K | dim F_h^reg |
| --- | --- | --- | --- |
| {2,3,5,6} | A1+A1+A2 | 16 | 10 |
| {2,3,5} | A1+A1+A1 | 12 | 6 |

In both cases the quotient-critical scheme is finite and reduced over
generic targets, with zero excess dimension. No generic orbit count is
asserted for these two cases. The theorem follows from two reusable
criteria below; it is not an inference from sampled fibers.

Both refinements come from the SAME polynomial

    chi=sum_(i=1)^3 a_i b_i + a_1 a_2 a_3 + b_1 b_2 b_3 -3. (2)

Refining one SL3 pair means substituting

    a_i=z_i x_i+z_i^(-2), b_i=z_i^(-1)x_i+z_i^2.             (3)

These are the ordinary defining and dual SL3 characters restricted to
SL2 times C*. Refining one, two, or three pairs gives respectively the
preceding A2+A1+A2 case and the two cases of Theorem 1. Equation (2)
is a useful character identity, not an assertion that a maximal-rank
SL3^3 subgroup is itself a Levi. The full central-charge records verify
its pullback to the Levis under consideration.

Use z_1=p,z_2=z,z_3=q, and central cocharacters

    H_p=(-2,0,-1,0,0,0),
    H_z=(2,3,4,6,4,2),
    H_q=(0,0,0,0,-1,-2)                                    (4)

in the simple-coroot basis. The two-refinement case uses H_p,H_z;
the three-refinement case uses all three. They satisfy

    B_Ad(H_i,H_j)=144 delta_(ij).                            (5)

Thus for h=sum t_i H_i the logarithmic target parameters are c_i=144t_i.
The finite central multiplication covers have degrees 12 and 24. They
are etale, so smoothness and dimension descend. We do not count the
cover points as distinct orbits downstairs.

## 2. Exact transfer for products of quadratic factors

Let independent torus blocks z_i have logarithmic derivatives D_i,
let S_i=S_i(z_i) be Laurent polynomials, and put

    P_i=x_i(x_i+S_i),
    chi=kappa product_i P_i + sum_j u_j v_j + constant,       (6)

where kappa!=0. Additional unused affine coordinates may be retained.
The potential is corrected by -sum c_(i,l) log z_(i,l).

**Theorem 2 (product transfer).** On the target open where each block
c_i is nonzero, the critical scheme of (6) is exactly the critical
scheme on the torus for

    R(z)=kappa*(-1/4)^m product_i S_i(z_i)^2,                 (7)

with x_i=-S_i/2 and u_j=v_j=0. Any unused coordinates remain free.
This is a scheme isomorphism, without taking radicals.

Proof. Eliminate u,v. A nonzero coordinate of the i-th target block
satisfies c_(i,l)=kappa*x_i*D_(i,l)S_i*product_(j!=i)P_j.
In that localized quotient ring, x_i and every P_j with j!=i are units.
For m>=2 all P_j are units; for m=1 the derivative equation directly
sets 2x_1+S_1=0. In either case the x derivatives give 2x_i+S_i=0.
Substitution gives P_i=-S_i^2/4 and precisely (7). Conversely (7)
satisfies all original equations, including the same unit conditions.

The torus map D R has smooth source and equal source/target dimension.
Consequently its generic fibers are smooth and finite, possibly empty.
They are nonempty if D R is dominant. One sufficient dominance test is
that S=product S_i has nonnegative real coefficients and its exponent
vectors span the entire real torus space. Its logarithmic Hessian
argument is the same as for the square-discriminant suspension: S>0
and Hess(S^2) is positive definite at positive torus points. Requiring
the PRODUCT support to span is essential; separate monomials S_i
would not give a dominant map in more than one variable.

### Apply the transfer to A1+A1+A2

In the [preceding E6 character](BRANCHING_NORM_FIBERS.md#5-three-mechanism-tests-and-an-e6-theorem),
replace the first SL3 pair by a=px+p^(-2), b=p^(-1)x+p^2.
Keep u,v for the second SL3 and y for SL2. Define

    T=p^3+p^(-3), S=z^3+z^(-3),
    A=z^(-2)+zy, B=z^2+z^(-1)y.

The invertible affine change u'=u+Aa, v'=v+Bb gives

    chi=u'v'-(x^2+Tx)(y^2+Sy)-2.                            (8)

Indeed AB=1+y^2+Sy and ab=1+x^2+Tx. This identity includes all central
variables; the target derivatives must be taken after the entire
coordinate change, not by dropping its z-dependence.

For c_p c_z!=0, Theorem 2 reduces the exact quotient-critical scheme to

    D_p R=c_p, D_z R=c_z,
    R=-(p^3+p^(-3))^2*(z^3+z^(-3))^2/16.                    (9)

The product support spans R^2 with positive coefficients. Thus the
generic quotient fibers are nonempty finite and reduced. The regular
quotient lift has relative dimension dim K-6=10, proving the first
part of Theorem 1. The nonzero-block restriction is a dense target
open, not a statement about every target.

## 3. A finite Fricke critical algebra and a boundary criterion

The familiar Fricke polynomial has the form

    f=X^2+Y^2+Z^2+XYZ+A X+B Y+C Z.                            (10)

Its relative gradient algebra over Q[A,B,C] is free of rank five,
with basis 1,X,Y,Z,Z^2. A monic Groebner basis in graded reverse
lexicographic order X>Y>Z is

    X^2-Z^2+(A X-C Z)/2,
    XY+2Z+C,
    Y^2-Z^2+(B Y-C Z)/2,
    XZ+2Y+B,
    YZ+2X+A,
    Z^3+(C/2)Z^2+B X+A Y-4Z-2C+AB/2.                       (11)

These are consequences of the three gradient equations; conversely
those equations occur in (11). The fifteen S-polynomials reduce to
zero over Q[A,B,C]. Thus freeness holds at EVERY parameter, not merely
over its fraction field. This is an explicit elementary algebra
calculation for a published cubic family, not a claim that Fricke
surfaces or their critical points are new.

**Theorem 3 (Fricke transfer with controlled boundary).** Let chi be a
Laurent family in three affine variables over a torus T of dimension s.
Suppose on U={d!=0} an invertible affine change and multiplication by a
unit put its relative critical equations into the form grad f=0.
Assume one fiber over U consists of five reduced points. Suppose also
that the image under mu=D_z chi of the relative critical locus over
T minus U is nondense in A^s. Then the generic fibers of mu on the
ENTIRE original relative critical scheme are smooth and finite,
possibly empty. A point with nonsingular relative Hessian and full
target differential proves nonemptiness generically.

Proof. Over U, (11) gives a finite flat algebra of rank five. The reduced
fiber makes its discriminant nonzero. The algebra is generically
reduced; flatness over the integral smooth base makes it torsion-free,
so it is reduced everywhere. All its components have dimension s and
dominate U. Its nonsmooth locus has dimension at most s-1, hence
nondense target image. Exclude that image and apply characteristic-zero
generic smoothness on the smooth locus, also excluding nondominant
components. Finally exclude the stipulated boundary image. This
argument never replaces the original critical scheme by its radical.

A practical way to establish the boundary hypothesis is useful here.
If a boundary stratum is a translate of a proper subtorus, its relative
critical scheme is an affine bundle, and chi is constant on that
bundle, then mu lies in the FIXED annihilator of that subtorus's Lie
algebra. Its target image is proper. Remaining strata of total critical
dimension <s also have nondense image. Constancy alone on an arbitrary
curved base would not suffice: its conormal directions can vary.

Multiplication by a base unit preserves the RELATIVE gradient algebra,
but it changes the target derivatives. Throughout Theorem 3, mu means
the derivatives of the original chi. Rank five is the length over
fixed central torus parameters, NOT a Springer generic degree, generic
orbit count, or length over fixed target c.

## 4. Three refinements and the vanishing cubic coefficient

Substitute (3) in all three pairs of (2). Write

    Z=z_1 z_2 z_3, D=Z+Z^(-1), E=Z^(-1)-Z,
    delta_i=(z_i^3-z_i^(-3))/2,
    X_i=x_i+(z_i^3+z_i^(-3))/2.

Direct expansion gives

    chi=sum X_i^2 + D X_1 X_2 X_3
        +E sum_i delta_i X_j X_k
        +D sum_i delta_j delta_k X_i
        +E delta_1 delta_2 delta_3-sum delta_i^2-3,            (12)

where {i,j,k}={1,2,3} in each sum. On D!=0 translate
X_i=Y_i-E delta_i/D and set U_i=D Y_i. The relative critical equations
become those of (10), with

    A_i=4 delta_j delta_k-2E delta_i.                         (13)

The whole potential is D^(-2)f(U)+a central Laurent-rational term;
we use (13) only for the relative critical algebra.

At z_1=z_2=z_3=1, the five relative critical points have X=(0,0,0)
or X_i in {+1,-1} with product -1. Their Hessian determinants in X
are respectively 8 and -32. Thus the reduced-fiber condition holds.

Now keep, rather than discard, D=0. Each of its two components is the
translated subtorus Z=+/-i. On either component, expansion of the
Laurent expressions gives

    1+sum delta_i^2-E delta_1 delta_2 delta_3=0,
    chi=(1/2) X^t Q X-2,
    Q_ii=2, Q_ij=E delta_k (k distinct from i,j), det Q=0.   (14)

The principal two by two minors are 4(1+delta_k^2).
Where rank Q=2 the relative critical locus is the kernel line bundle
over the two-dimensional boundary. Its critical value is constantly
-2. The boundary criterion puts its target image in the line
c_1=c_2=c_3, the normal direction to Z=constant.

Rank Q<=1 requires delta_i^2=-1 for all i, equivalently z_i^3=+/-i.
This is a finite base locus. Q cannot have rank zero because its
diagonal entries are two; its critical fibers here have dimension two.
Their total dimension is therefore two, strictly below target dimension
three. These strata also have nondense target image. This handles
every point over D=0, including possible nonreduced boundary structure.

Theorem 3 proves generic smoothness and finiteness. For nonemptiness,
at z_i=1 and x_i=-1 the full Hessian in affine x and logarithmic z
coordinates is

    diag(2,2,2,-18,-18,-18).                                 (15)

The relative Hessian is invertible and the target map on its local
critical branch has invertible differential. Hence mu is dominant on
that component. The witness has c=0 and need not be G-regular; it is
used ONLY to prove dominance on the smooth quotient. Generic central
targets with Z_G(h)=K then have regular lifts by the usual quotient
theorem and equivariance. Their dimension is dim K-6=6. This proves
the second part of Theorem 1.

## 5. Independent checks, published input, and the next question

The [Sage verifier](tools/verify_cubic_levi_refinements.py) checks the
master identity against both full central-charge records, the two
transfers, the boundary identities, and all five relative Hessians.
It independently constructs regular semisimple torus points over
F_1009 and evaluates the original adjoint root trace and differential.
The ranks are 68 and 72, giving kernel dimensions 10 and 6. These
finite-field points corroborate the characteristic-zero proofs; they
are not proofs by specialization or generic-degree counts.

The [native auditor](tools/audit_cubic_levi_refinements.py), without
Sage, verifies the universal monic Groebner basis including all fifteen
S-pairs, the exact transfer identities, and the direct root evaluations
using separately enumerated Cartan reflections. It verifies the changed
central lattice and target normalizations. See the
[verification](cubic-levi-refinements-verification.json) and
[audit](cubic-levi-refinements-audit.json) reports.

The underlying trace construction is [Kostant--Michor's](https://www.mat.univie.ac.at/~michor/kostant.pdf).
Regular quotient geometry, generic smoothness and Groebner specialization
are standard inputs. The Fricke family is classical; for a recent primary
description see [Rebelo--Roeder, introduction](https://arxiv.org/pdf/2307.10962).
Their cubic-surface and dynamical questions are not the target-map
question treated here. Our use of its relative gradient algebra, the
constant-value boundary argument, and the resulting E6 fiber deductions
are written here without a priority claim. Targeted searches on
2026-10-23 located the Fricke and generalized-Cayley sources; they do
not establish novelty of these deductions.

The three original cubic mechanism tests are now resolved. This does
NOT prove all E6/E7/E8 adjoint Levi cases: identical abstract Levi types
can have different embeddings and central coefficients, and other cubic
couplings remain possible. The next structural question is whether
the product/Fricke reductions are forced by a recognizable branching
configuration, especially when an A1 factor couples to more than two
other factors. A useful next result is an embedding-sensitive criterion
or an explicit obstruction, rather than an unchecked propagation by
Levi type. Whole fibers, nonsemisimple targets, PDFs and public APIs
remain outside this step.
