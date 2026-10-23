# Symplectic centralizers: realization of every odd nilpotent map

Date: 2026-10-23. Written arbitrary-rank theorem in characteristic zero,
with independent exact finite checks. Novelty is not established.

## 1. Statement and comparison with type A

Fix m>=1 and G=Sp_(4m)(C), of rank 2m. Let V=C^(4m) be its defining
representation and let rho=V_(omega_(2m)) be its last fundamental
representation. This is the primitive middle exterior power

    ker(contraction: exterior^(2m) V -> exterior^(2m-2) V),

not the whole exterior power. It is irreducible and has character
e_(2m)-e_(2m-2). Its trace index relative to the defining trace is

    I=binomial(4m-2,2m-1)-binomial(4m-2,2m-3)>0.

Out-of-range binomials are zero. The normalized Springer map satisfies

    I tr(theta_rho(g) X)=tr(rho(g) d rho(X)).

We use the simply connected group and full weight lattice. This
representation has kernel {+1,-1} and faithful differential; it is
faithful on PSp_(4m). No integral cohomology source or Borel
characteristic homomorphism is involved.

**Theorem.** Given any odd polynomial

    P(Y)=p_1 Y+p_3 Y^3+...+p_(2m-1) Y^(2m-1),

including P=0, there are a noncentral semisimple t in G and c in C*,
with centralizer

    H=Z_G(t)=Sp_(2m) x (C*)^m,

such that, for every nilpotent Y in sp_(2m),

    u(Y)=(I+Y)(I-Y)^(-1),
    theta_rho(t diag(u(Y),I_(2m)))-theta_rho(t)
        = c diag(P(Y),0).                                 (1)

The blocks are an orthogonal direct sum of symplectic spaces. In a
standard paired basis this is conjugate to the usual symplectic torus
notation. The Cayley coordinate u(Y) and its inverse
Y=(u-I)(u+I)^(-1) are polynomial on their nilpotent/unipotent varieties,
because the inverse series terminate. They serve as coordinates here;
they are not being identified with the normalized theta_rho itself.

Every Sp_(2m)-equivariant algebraic map from its nilpotent cone to itself
is an odd polynomial of the displayed form. Thus the fixed irreducible
ambient representation realizes all these maps, up to a nonzero target
scalar, by varying t. Every order 1,3,...,2m-1 and the zero map occurs.
The entire translated centralizer map is covered, since the torus factors
of H have no nontrivial unipotent elements.

Compared with [type A realization](IRREDUCIBLE_CENTRALIZER_REALIZATION.md),
the restriction to odd powers is essential: it comes from the symplectic
Lie algebra, not from a missing part of the calculation. This is a type C
realization theorem. The [type A full-fiber component and multiplicity
theorem](MIDDLE_EXTERIOR_FIBER_COMPONENT.md) does not automatically extend.

## 2. Character calculation on the centralizer

Choose t=I_(2m) direct-sum B, where B has distinct nonzero eigenvalues
in reciprocal pairs, none equal to +1 or -1. Write

    Q(z)=det(I+zB)=sum_(h=0)^(2m) q_h z^h,
    q_h=q_(2m-h),       q_0=q_(2m)=1.                      (2)

These conditions give exactly H=Sp_(2m) x (C*)^m. For g in G, the
fundamental-character generating expression is

    chi_rho(g)=[z^(2m)] (1-z^2) det(I+zg).                (3)

Let M^* be symplectic adjoint, so g^*=g^(-1) and Y^*=-Y. The trace
projection to the symplectic Lie algebra is M -> (M-M^*)/2. Differentiate
(3) using the defining trace pairing, and apply this projection. On the
first block, with g=u(Y), the relevant matrix function is

    f_z(Y)=z(1+Y)/((1+z)+(z-1)Y).

For j>=1 its Y^j coefficient is

    2z/(1+z)^2 * ((1-z)/(1+z))^(j-1).

Taking (f_z(Y)-f_z(-Y))/2 keeps exactly the odd powers. Also
det(I+z u(Y))=(1+z)^(2m). It follows that the translated first block is
sum_(j odd, 1<=j<=2m-1) b_j Y^j, where

    b_j = 2/I * [z^(2m-1)] Q(z) S_j(z),
    S_j(z)=(1-z)^j(1+z)^(2m-j).                           (4)

The external block's gradient is unchanged: its dependence on the
unipotent block is only through the constant determinant (1+z)^(2m).
At t the first projected block is zero. Hence (4) is an equality for
the full translated ambient map, not merely its first-block projection.
All rational expressions above become polynomial after clearing the
determinant, or after the finite nilpotent expansion.

## 3. Why the coefficient map is invertible

Let A be the space of degree-at-most-2m polynomials satisfying

    z^(2m) S(1/z)=-S(z).

It has dimension m. The m polynomials S_1,S_3,...,S_(2m-1) form a basis:
they satisfy this symmetry, and independence follows by dividing a
relation by (1+z)^(2m) and using the distinct powers of
(1-z)/(1+z).

For S=sum s_l z^l in A, use (2) to rewrite the functional in (4):

    [z^(2m-1)] Q S
       = sum_(l=0)^(m-1) s_l (q_(l+1)-q_(l-1)),           (5)
    q_(-1)=0.

The independent spectral parameters are q_1,...,q_m; q_0=1 is fixed.
Their variations map to

    delta q_1, delta q_2, delta q_3-delta q_1,
    ..., delta q_m-delta q_(m-2).

This linear map is triangular with diagonal one. Thus (4) is an affine
isomorphism from the m-dimensional palindromic spectral-coefficient
space to the m odd-map coefficients. This proves more than a rank
test on particular examples.

Its unique zero is

    Q_0(z)=1+z^2+...+z^(2m)
          =(1-z^(2m+2))/(1-z^2).                         (6)

Indeed q_(l+1)=q_(l-1), starting with q_0=1 and q_(-1)=0, gives
alternating coefficients 1,0,1,0,...,1.

For prescribed P, solve the invertible linear system for a palindromic
R with zero constant and top coefficients, so that Q=Q_0+R gives b_j=p_j.
Then Q_c=Q_0+cR gives b_j=c p_j.

## 4. The spectrum is symplectic and semisimple

The polynomial Q_0 is squarefree and Q_0(1)=Q_0(-1)=m+1. Its roots are
the (2m+2)th roots of unity other than +1 and -1. Consequently

    disc(Q_c) Q_c(1) Q_c(-1)

is a polynomial in c which is nonzero at zero. Choose c!=0 outside
its finite zero set. Q_c is then squarefree, reciprocal, has constant
and leading coefficients one, and has no roots +1 or -1. Writing
Q_c=product(1+beta_a z), the beta_a are distinct nonzero reciprocal
pairs avoiding +1 and -1. They are the eigenvalues of a regular
semisimple symplectic matrix B. Thus t=I_(2m) direct-sum B belongs
to Sp_(4m), has exactly the stated centralizer, and (1) follows.

For P=0, take Q_0 itself. This gives a finite-order semisimple point
whose translated map is zero on the whole unipotent variety U(H).
The full ambient Springer fiber over theta(t) therefore contains the
closed subvariety t U(H), of dimension

    dim Sp_(2m)-rank Sp_(2m)=2m^2.                         (7)

This containment alone does not establish component status, local
multiplicity, or the dimension of the whole fiber. The subsequent
[component theorem](SYMPLECTIC_FIBER_COMPONENT.md) proves that its
target-centralizer sweep is a generically reduced component; it also
finds singular nonreduced structure already in rank one.

## 5. Why odd polynomials exhaust the equivariant maps

A regular nilpotent J in sp_(2m), acting on its defining 2m-dimensional
space, is a single Jordan block. If F is equivariant then exp(sJ) fixes
F(J), so F(J) commutes with J. Thus F(J) is a polynomial in J of degree
at most 2m-1. The condition F(J)^*=-F(J), together with J^*=-J, removes
all even powers. Equivariance and density of the regular nilpotent orbit
extend this polynomial identity to the whole nilpotent cone.
This is the established regular-centralizer argument, included to make
the realization theorem's scope precise.

If the first surviving power is r, then r is odd. Write
P(Y)=p_r Y^r v(Y^2), v(0)=1. The finite odd-polynomial change
Z=Y v(Y^2)^(1/r) has an odd-polynomial inverse, preserves the symplectic
Lie algebra and nilpotence, and gives P(Y)=p_r Z^r. Hence these orders
also give the matrix-power normal forms. No even order can occur for
an equivariant map on this symplectic nilpotent cone.

## 6. Published inputs and verification

The primitive-exterior construction, character subtraction and trace
index are classical; see Fulton--Harris, Theorem 17.5, and the source
audit in [our fundamental type C note](TYPE_C_FUNDAMENTAL_REGULARITY.md).
The prior audit inspected the textbook directly. A new web fetch in
this turn was blocked by the PDF size limit; no new full-text inspection
of that textbook is claimed. The primitive representation must not be
replaced by the full middle exterior power.

[Kostant--Michor, Sections 2 and 4](https://arxiv.org/pdf/math/0109066)
provide the trace-projection construction and translated centralizer
maps. The Cayley coordinates and regular-nilpotent polynomial argument
are standard inputs. The invertible pairing (4)-(5), the good zero
spectrum (6), and the resulting fixed-irreducible-family realization
are deductions written here. Targeted searches did not identify this
particular formulation, but novelty has not been established.

The [independent verifier](tools/verify_symplectic_centralizer_realization.py)
checks the coefficient isomorphism in symbolic families, then uses
actual primitive representation matrices in small ranks. It constructs
the contraction kernel directly and restricts compound group matrices
and signed-wedge Lie actions to it. Thus the main trace check does not
merely subtract exterior characters or reuse the coefficient formula.
It tests the entire ambient symplectic Lie algebra, including directions
outside the centralizer, and checks the trace index and symplectic forms.
The calculations share Sage arithmetic and earlier wedge helpers;
they are independent corroboration, not an external audit or proof of
all ranks by finite testing.

## 7. What this advances, and the next question

The realization framework now extends to a genuinely symplectic
centralizer in a fixed irreducible type C family. It permits precisely
the odd powers allowed by the Lie algebra, with unbounded orders as
m grows, and includes total collapse. This is not just a type A Levi
embedded in a symplectic group.

The [follow-up proof](SYMPLECTIC_FIBER_COMPONENT.md) now settles the
component problem for Q_0 and its generic multiplicity. It formulates
a conjecture about reduced support at the remaining singular points,
allowing the embedded infinitesimal structure found in Sp_4. The
realization argument alone does not imply these stronger statements.
