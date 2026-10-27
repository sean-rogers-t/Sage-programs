# Full Weyl symmetry does not force the required critical-jet compatibility

26 October 2026. Continuation of the [shared-coupling first-jet gate](SHARED_SQUARE_COUPLINGS.md#7-attribution-remaining-root-problem-and-next-gate).
This is a counterexample to a proposed PROOF INPUT, not to the adjoint
conjecture. The representation below is actual and belongs to a simple
simply connected group, but is reducible and is not the adjoint.

## 1. Statement

Let G=SL_3(C), H=diag(1,1,-2), and K=S(GL_2 x GL_1). There is an actual
640-dimensional rational representation rho of G with nondegenerate trace
form such that, for generic lambda!=0, the regular-source fiber

    theta_rho^(-1)(lambda H) intersect G_reg

is nonreduced at regular SEMISIMPLE points. At the points exhibited below,
its completed local ring is

    C[[o_1,o_2,x]]/(x^2).                                  (1)

Here theta_rho is normalized by B_rho(theta_rho(g),X)=tr(rho(g)d rho(X)).
The same assertion holds for every tensor power rho^(tensor n), n>=1,
with its own trace normalization. Thus this is an infinite representation
family with the obstruction, not just an abstract potential or a
representation of a reductive Levi.

Genericity is WITHIN the one-dimensional proper-Levi-central target
stratum {lambda H: lambda!=0}. There is no claim of nonreduced fibers
over generic points of the whole Lie algebra. There is also no claim
about irreducible representations. The all-simple ADJOINT conjecture and
the classical ADJOINT theorem remain unchanged.

## 2. An actual representation, with full ambient Weyl symmetry

Write V_(p,q) for the irreducible SL_3 representation with highest weight
p*omega_1+q*omega_2. Let a=chi_(1,0), b=chi_(0,1), U=ab-1=chi_Ad, and

    V=a^3+b^3-3ab+5
      =chi_(3,0)+chi_(0,3)+chi_(1,1)+4.

Define the following polynomial character:

    P=V^2-6U^2.                                            (2)

Although written as a difference, it has the POSITIVE irreducible expansion

    P=13 chi_(0,0)+4 chi_(1,1)
       +6(chi_(3,0)+chi_(0,3))+3 chi_(2,2)
       +3(chi_(4,1)+chi_(1,4))
       +chi_(6,0)+chi_(0,6)+2 chi_(3,3).                    (3)

Take rho to be this direct sum. The dimension formula
(p+1)(q+1)(p+q+2)/2 gives dim rho=640. Every coefficient in (3) is
a nonnegative integer; there is no virtual or complex-coefficient
representation in the assertion. The character is a polynomial in
the two ambient fundamental characters, hence is fully S_3-invariant.

One elementary way to verify (3) is to expand V^2 using these familiar
type-A tensor-product identities (each also checked independently):

    V_(3,0)^2 = V_(6,0)+V_(4,1)+V_(2,2)+V_(0,3),
    V_(3,0)*V_(0,3) = V_(3,3)+V_(2,2)+V_(1,1)+1,
    V_(3,0)*V_(1,1) = V_(4,1)+V_(2,2)+V_(3,0)+V_(1,1),
    V_(1,1)^2 = V_(2,2)+V_(3,0)+V_(0,3)+2V_(1,1)+1,

together with their duals. Standard highest-weight and tensor-product
theory are the inputs; (2)-(3) are an explicit finite character identity.

## 3. The exact derived and central jets

Use the usual finite central cover SL_2 x Gm -> K. On its torus write

    g=diag(zp,z/p,z^(-2)), x=p+p^(-1),
    S=z^3+z^(-3), R=z^3-z^(-3), D=z partial_z.

Then a=zx+z^(-2), b=z^(-1)x+z^2, so

    U=x^2+Sx,
    V=S^2+3x^2+Sx^3.

The full restricted character f=P(a,b) is EXACTLY

    f=S^4+2S(S^2-6)x^3+3x^4+6Sx^5+S^2x^6.                (4)

Consequently

    f_x=6x^2 A,
    A=S(S^2-6)+2x+5Sx^2+S^2x^3,
    mu=Df=M(z) modulo (x^3),
    M(z)=12S^3R,
    DM=144S^2(S^2-3).                                    (5)

Thus both the linear and quadratic derived-character coefficients
vanish IDENTICALLY as functions of z, while the central target map on
x=0 is nonconstant. The cubic coefficient is nonzero generically.
This is precisely the bad compatibility that symmetry alone was hoped
to exclude. It occurs for the entire actual character (3).

This is also visible directly in trace moments. Put p=i, so x=0.
For D_0=diag(1,-1,0), direct sums over the representation weights give

    tr(rho(g)d rho(D_0))=0,
    tr(rho(g)d rho(D_0)^2)=0,
    tr(rho(g)d rho(D_0)^3)=-96 i S(S^2-6),
    tr(rho(g)d rho(H))=12S^3R.                             (6)

These are derived and central jets of the SAME fully Weyl-invariant
character; no subsystem criticality or omitted-root alteration is used.

## 4. Actual generic-target ideal and trace normalization

Work first with the raw target c=B_rho(h,H). Over L=C(c), retain the
actual derived gradient f_x and central equation Df-c. Invert A only
on its explicit chart. By (5) the resulting ideal is EXACTLY

    (f_x,Df-c)=(x^2,M(z)-c).                              (7)

All generic roots of M(z)=c lie in this chart: S(S^2-6) is nonzero at
them. Indeed the Laurent equation, multiplied by the unit z^12, is

    12(z^24+2z^18-2z^6-1)-c z^12=0.                       (8)

Its constant term is -12, and over C(c) it is coprime both to its
derivative and to (z^6+1)*((z^6+1)^2-6z^6). Thus it has 24 distinct
nonzero geometric roots, all with A(0,z)!=0. The localized quotient
has length 48 and reduced length 24 on this finite central cover.
These are lengths of the displayed transverse CHART, not the degree
of the full Springer map, the whole fiber, or counts of group orbits.

For explicit complex targets it suffices to exclude

    c=0, c^2=62208, c^2=-3888.                            (9)

The first two avoid zeros of A(0,z); the last avoids zeros of DM.
At each selected z_0 the formal implicit-function theorem eliminates
z-z_0 from (7), leaving C[[x]]/(x^2).

The trace form is nondegenerate and the normalization is explicit:

    B_rho(H,H)=10368, B_rho(D_0,D_0)=3456,
    B_rho(H,D_0)=0.

Thus B_rho=1728 tr_def on sl_3, and h=lambda H corresponds to
c=10368 lambda. Scalar normalization does not remove the nilpotent.

At p=i, g=diag(iz,-iz,z^(-2)) is regular semisimple provided S!=0:
an equality with the third eigenvalue occurs exactly when z^6=-1.
For c satisfying (9), all the chosen points are therefore regular
semisimple. At lambda!=0, Z_G(lambda H)=K exactly.

The [regular critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md)
identifies (7) with the transverse fiber scheme. Alternatively use the
regular-semisimple conjugation chart in K: the internal torus variable
p has d(p+p^(-1))/dlog p=p-p^(-1)=2i at this point, so x is a valid
local quotient coordinate. The K-orbit has dimension dim K-rank G=2.
Equivariance puts the whole fiber scheme inside K, and the quotient
map is smooth there. Adding these two smooth orbit directions gives
(1). This passage retains the actual nilpotent class x.

In particular no nonzero polynomial in the central target can annihilate
the entire nilradical of this critical incidence on the displayed open:
such a polynomial becomes a unit over C(c), where x remains nonzero.
This is a persistent nonreduced component; an embedded component is not
needed for the obstruction in this particular example.

## 5. An infinite tensor-power family

Let rho_n=rho^(tensor n), n>=1, and d=640. Since g is semisimple as a
Lie algebra, tr(d rho(X))=0. The tensor trace formulas give

    chi_(rho_n)=f^n,
    B_(rho_n)=n d^(n-1) B_rho,
    theta_(rho_n)(g)=(f(g)/d)^(n-1) theta_rho(g).            (10)

On the same x=0 locus, f=S^4 is nonzero generically. Therefore
partial_x(f^n) is x^2 times a unit on an open containing its generic
target points. Modulo x^2 its raw target is

    M_n=12n S^(4n-1) R,
    D M_n=144n S^(4n-2)*(n S^2-(4n-1)).                    (11)

This is nonconstant for every n. Only finitely many target values
come from the excluded zeros or critical points. The same formal
argument proves (1) for generic targets in the specified stratum
for every n, with c_n=10368 n 640^(n-1) lambda. Equivalently, at the
exhibited source point the normalized target coefficient is

    lambda_n=S^(4n-1)R/(864*640^(n-1)).                    (12)

No large-rank computation or decomposition of all tensor powers is
used: actual representations stay actual under tensor products, and
(10)-(11) prove the all-n assertion. These representations are still
reducible and are not adjoint representations.

## 6. What this says about the ADJOINT first-jet problem

This disproves the implication

    full ambient Weyl symmetry + actual simple-group character
       => generic reducedness on every proper-Levi-central stratum.

It therefore cannot be the missing proof input by itself. It also
strengthens the earlier altered-Levi controls: this obstruction occurs
on the simple group itself and respects its entire Weyl group.

It does NOT survive unchanged in the actual SL_3 adjoint character.
There U=x^2+Sx, and U_x at x=0 equals S. A regular semisimple source
on this locus has S!=0, so it is not even in the adjoint critical
incidence. The added higher-weight constituents in (3) create the
first- and second-jet cancellation. The adjoint's precise root support
prevents this particular construction.

More generally, in simply laced type, a Weyl-invariant Laurent sum
supported only at zero and the roots has the form a_0+m sum_alpha e^alpha,
because the roots form one Weyl orbit. If m!=0, its trace-normalized map
is theta_Ad; the constant term has no differential. This elementary
observation was already recorded in the earlier proof audit. It supplies
no NEW relation on the singular critical ideal. In nonsimply laced type
there are two root orbits, and their multiplicities must also be fixed
to recover the adjoint character.

**Exact remaining gap.** A relation using the multiplicative root
characters, their fixed adjoint multiplicities and the full derived
gradient ideal must still force the singular-incidence image to be
nondense. No such relation for arbitrary overlapping exceptional
couplings has been proved here. A symmetry-only first-jet argument
would incorrectly apply to (3), so that route must stop at this point.
The next meaningful gate is an explicit root-support-dependent syzygy
or target-polynomial certificate, retaining embedded components.
Neither another arbitrary-character example nor a Levi census closes it.

## 7. A uniform identity that really does retain adjoint root support

There is an elementary higher-jet identity stronger than mere Weyl
invariance. Let Phi=rank(G)+sum_(beta in Delta) e^beta be the ACTUAL
adjoint character, alpha any LONG root, and partial_alpha the torus
logarithmic derivative in direction alpha^vee. Then in the torus ring

    (partial_alpha^3-partial_alpha) Phi
        =6(e^alpha-e^(-alpha)),
    (partial_alpha^4-partial_alpha^2) Phi
        =12(e^alpha+e^(-alpha)).                          (13)

Proof: the integers <beta,alpha^vee> lie in {-2,-1,0,1,2}, and the
values +2 and -2 occur only for beta=alpha and -alpha. This follows
from the root-length bound and equality in Cauchy--Schwarz, since
alpha is long. Apply the two displayed interpolation polynomials to
each weight monomial. The zero-weight constant is killed. This proof
works in every finite simple type, not just the simply laced ones.

In particular the fourth-order operators

    E_alpha^+ = (partial_alpha^4-partial_alpha^2)Phi/24
                  +(partial_alpha^3-partial_alpha)Phi/12,
    E_alpha^- = (partial_alpha^4-partial_alpha^2)Phi/24
                  -(partial_alpha^3-partial_alpha)Phi/12

recover e^alpha and e^(-alpha) EXACTLY. Hence

    E_alpha^+ E_alpha^-=1,
    E_alpha^+ E_beta^+=E_(alpha+beta)^+

whenever alpha,beta,alpha+beta are long roots. These nonlinear jet
relations retain both the root support and its multiplicity one.
They fail for arbitrary Weyl-invariant representation characters: at
the SL_3 source above, e^alpha=-1, while (6) makes the first left side
of (13) nonzero generically and its purported right side zero.

These are elementary consequences of known root theory, not a novelty
claim or a new smoothness theorem. They do NOT automatically descend to
the relative critical quotient: differentiation of an equation in an
ideal need not remain in that ideal. Nor can the torus critical ideal
replace the ordinary fundamental-character gradient on an internal
Weyl discriminant. No target-rank defect or target-annihilator syzygy
has yet been deduced from (13). This is the concrete adjoint-specific
input to test next against the full gradient generators and Hessian
minors, with those two limitations retained.

## 8. Independent checks, attribution, and scope

The [native audit](weyl-jet-obstruction-audit.json) enumerates the weights
of the positive sum (3) using GL_3 interlacing patterns for partitions
(p+q,q,0). It verifies the character identity, both simple Weyl
reflections, all jets in (5)-(6), trace normalization, and coprimality
over QQ(c). It does not import Sage or a Sage character table.

The [independent Sage/Singular verification](weyl-jet-obstruction-verification.json)
constructs V^2-6 Ad^2 in Sage's Weyl character ring, checks its positive
decomposition and weights, and checks (7) from the original gradients.
It represents z^(-1) and A^(-1) by explicit inverse variables, retaining
the nilpotents. The ideal, its radical and the annihilator of x are
computed over QQ(c); lengths 48 and 24 are independently checked.
There is no generic-source-only substitution or modular specialization.

Both audits also check the monomial projectors in (13) on A2 and all
five exceptional root systems. The native code uses independently
specified Cartan matrices and rational root reflections; Sage uses
ambient root vectors and inner products. The arbitrary-rank argument
is the root-length proof, not extrapolation from these finite checks.

Published inputs are highest-weight/tensor-product theory, the trace
construction and smooth regular-semisimple quotient geometry; the
critical-presentation deduction is documented in the linked repository
note. Sage's [representation-theory tutorial](https://doc.sagemath.org/pdf/en/thematic_tutorials/thematic_tutorials.pdf)
describes its character and tensor-product computations. The explicit
positive character and fiber obstruction are deductions here. No
exhaustive literature search, novelty claim, or external peer review
is asserted. The finite identities and independent algorithms support
the written proof; they do not resolve the uniform adjoint conjecture.
