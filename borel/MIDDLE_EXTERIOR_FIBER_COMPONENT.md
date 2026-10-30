# Full Springer-fiber components in the middle-exterior zero model

Date: 2026-10-23. Arbitrary-rank characteristic-zero proof, with independent
finite exact matrix checks. Novelty is not established.

## 1. Setup and theorem

Let m>=2, G=SL_(2m)(C), rho=exterior^m C^(2m),
D=binomial(2m-2,m-1), and alpha=i. The normalized trace Springer map
uses B_rho(X,Y)=D tr(XY). The representation is irreducible,
infinitesimally faithful, and has central kernel mu_m; it becomes faithful
on SL_(2m)/mu_m. We work on the simply connected group, with the full
weight lattice. There is no integral cohomology or Borel characteristic
map assertion here.

Let gamma_1,...,gamma_m be all roots of x^m+1, beta_j=alpha gamma_j, and

    t=diag(alpha I_m,beta_1,...,beta_m),
    H=Z_G(t)=S(GL_m x (C*)^m),
    h=theta(t),     K=Z_G(h),     F=theta^(-1)(h).

F denotes the scheme-theoretic fiber. The [realization theorem](IRREDUCIBLE_CENTRALIZER_REALIZATION.md)
proved that theta(tu)=h for every u in the unipotent variety U(H).
Write C=K . (t U(H)), where the dot means conjugation.

**Theorem (component and generic scheme structure).** C is closed and
irreducible and is an irreducible component of the underlying reduced
fiber F_red. The following distinctions hold.

| | m even | m odd |
| --- | --- | --- |
| K | H | S(GL_m x GL_2^((m-1)/2) x GL_1) |
| C | t U(H) | A strictly larger conjugation sweep of t U(H) |
| Dimension d of C | m^2-m | m^2-1 |
| Generic multiplicity in F | 1 | 2 |
| Tangent dimension along the regular-unipotent part of C | m^2-m | m^2 |
| Rank of d theta there | 3m^2+m-1 | 3m^2-1 |

More precisely, at every point g of C whose first block has one Jordan
block, the completed local ring of F is

    C[[z_1,...,z_d]]                         if m is even,
    C[[z_1,...,z_d,epsilon]]/(epsilon^2)      if m is odd.   (1)

For odd m the original t U(H) is therefore not a component. The
larger component has one nonreduced transverse direction at its generic
point. This is a full ambient-fiber statement, not merely a statement
about the map restricted to H. It does not classify other components
of F, intersections at nonregular points, or the generic degree of theta.

## 2. The target centralizer and the parity split

For any g, let e_j(g) be the elementary characteristic coefficients and
put

    F_m(g)=sum_(k=1)^m (-1)^(k-1) e_(m-k)(g) g^k,
    theta(g)=(F_m(g)-tr(F_m(g))/(2m) I)/D.                 (2)

This is the determinant-character gradient from the previous proof.
In particular g commutes with theta(g).

The external block B_0=diag(beta_j) has characteristic polynomial
x^m+alpha^m, so e_1(B_0)=...=e_(m-1)(B_0)=0. The first block of
F_m(t) is alpha^m I_m. Its difference from the jth external entry is

    K_j=(alpha-beta_j)^m=alpha^m(1-gamma_j)^m.             (3)

Thus the corresponding eigenvalue difference for h is K_j/D, which
is nonzero. There is no enlargement mixing the first block with the
external eigenspaces.

All gamma_j have absolute value one. Equality of (1-gamma_j)^m and
(1-gamma_k)^m forces equality of their absolute values, hence
gamma_k is either gamma_j or gamma_j^(-1). For inverse roots,

    (1-gamma^(-1))^m = (-1)^(m+1)(1-gamma)^m.              (4)

For even m distinct external target eigenvalues remain distinct. For
odd m they pair by inversion, with the single self-inverse root -1
left over. This proves the descriptions of K in the table.

The group t U(H) consists of the first block alpha times an arbitrary
unipotent matrix, with the external diagonal fixed. For odd m, each
GL_2 factor in K sweeps its two distinct beta eigenvalues through their
closed semisimple conjugacy class of dimension two. Consequently C is
the product of a translated GL_m unipotent variety and (m-1)/2 such
classes. This proves closedness, irreducibility and the dimension
formula before the component assertion. Equivariance gives C subset F.

## 3. Why the ambient fiber can be computed inside K

The polynomial identity [g,theta(g)]=0 implies

    [g,h]=[g,h-theta(g)].                                  (5)

Hence the linear equations defining the matrix centralizer of h belong
to the ideal of the fiber. This containment is scheme-theoretic, including
infinitesimal points: F is a closed subscheme of K, not just a set of
commuting C-points.

Fix g_0 in C with a regular first block. In K write g=diag(A,B).
The first block A has size m and is cyclic in an open neighborhood of
g_0. The external B is diagonal when m is even and has the indicated
2-by-2 blocks when m is odd. Near B_0 these blocks are regular semisimple.
After an etale coordinate change their eigenvalues can be labeled
b_1,...,b_m near beta_1,...,beta_m. Their conjugation coordinates contribute
zero free variables for even m and m-1 free variables for odd m.
The target is scalar on each such block, so these conjugation coordinates
do not occur in the fiber equations.

Let a_l=e_l(A). For cyclic A the characteristic-coefficient map is smooth
of relative dimension m^2-m. This can be seen directly from a companion
matrix section: every cyclic A is conjugate to its companion matrix,
whose characteristic coefficients vary independently. Thus formal local
coordinates consist of a_1,...,a_m and m^2-m orbit variables. This is
the standard regular adjoint-quotient fact, not a new smoothness theorem.

The first block of theta(g) must be scalar. Since I,A,...,A^(m-1) are
linearly independent on a cyclic chart, reducing (2) by the characteristic
polynomial of A makes this condition equivalent, as an ideal, to

    e_d(A direct-sum B)=a_d,       1<=d<=m-1.

The convolution e_d(A+B)=sum_l a_(d-l)e_l(B), with a_0=1, is triangular
with leading coefficient one in e_d(B). These equations are therefore
equivalent to

    e_1(B)=...=e_(m-1)(B)=0.                              (6)

One can check the ideal equivalence on an invertible minor of the columns
I,A,...,A^(m-1); it does not follow merely from pointwise cyclicity.
After (6), the first block of F_m(g) is exactly a_m I_m.

Since the beta_j are distinct, (6) is locally the one-dimensional smooth
family

    b_j=s beta_j,       s near 1.

Indeed its characteristic polynomial is x^m+constant, and its simple
roots have a common scaling on the chosen etale branch. The group
determinant constraint then gives

    a_m=alpha^m s^(-m).                                  (7)

All transformations so far are changes of coordinates or ideal
equivalences on a neighborhood of g_0; no radical has been taken.

## 4. Exact elimination and the transverse equation

Under (6), the first-block/external differences in F_m(g) are
det(A-b_j I). Since theta(g) and h are both traceless, equality of
these m differences, together with the scalar block conditions, is
equivalent to theta(g)=h. Thus the remaining equations are

    det(A-s beta_j I)=(alpha-beta_j)^m,     j=1,...,m.      (8)

Expand

    det(A-s beta I)=sum_(k=0)^m (-1)^k a_(m-k) s^k beta^k.

Use beta_j^m=-alpha^m and the invertible Vandermonde evaluation on the
m distinct beta_j. Comparing coefficients beta^k, k=1,...,m-1, yields

    a_l=binomial(m,l) alpha^l s^(l-m),     1<=l<=m-1.      (9)

The constant coefficient, after (7), is

    s^(-m)+eta s^m=1+eta,          eta=(-1)^(m+1).

Since s is a unit, its ideal is equivalently generated by

    R_m(s)=1+eta s^(2m)-(1+eta)s^m
          = 1-s^(2m)          for even m,
          = (s^m-1)^2         for odd m.                   (10)

Equations (6)-(9) eliminate all transverse coordinates except s; the
orbit variables remain free. In characteristic zero, at s=1, (10)
is a unit times s-1 for even m and a unit times (s-1)^2 for odd m.
This proves the completed local rings (1).

Their reduced local rings have exactly the dimension already computed
for C, and are irreducible. Since C is closed and irreducible, it is
an irreducible component of F_red. The regular-first-block locus is
dense in C; the scheme multiplicities there are respectively one and
two. The tangent dimensions and differential ranks follow from (1).
In particular the multiplicity conclusion for odd m uses the full ideal
elimination, not just a calculation of tangent rank.

## 5. A transverse curve that independently tests the obstruction

Let A(s) be the companion matrix with coefficients (7),(9), B(s)=s B_0,
and g(s)=diag(A(s),B(s)). It lies in G and g(1) belongs to the regular
part of C. Direct substitution into the whole map gives

    theta(g(s))-h = alpha^m/(2D) s^(-m) R_m(s)
                    diag(I_m,-I_m).                      (11)

This gives a nonzero first-order obstruction when m is even, and an
exactly second-order obstruction when m is odd. In the latter case the
extra tangent direction cannot be continued to a reduced one-parameter
curve in the fiber. Equation (11) alone would not establish the full
local scheme; that conclusion requires Section 4.

The [independent verifier](tools/verify_middle_exterior_fiber_component.py)
compares the full derivative obtained from (2) against the original
trace construction in actual irreducible exterior matrices. It checks
all tangent pairings for m=2,3,4; polynomial differential ranks and target
centralizers for m=2,...,6; full transverse curve identities for m=2,...,8;
and exact cyclotomic target collisions through m=16. Curve checks for
m=2,3 also use compound exterior matrices directly. All computations are
exact and share Sage arithmetic; this is independent calculation, not
an external proof audit or a separate-CAS verification.

## 6. Relation to published work and remaining scope

[Kostant--Michor, Sections 2 and 4](https://arxiv.org/pdf/math/0109066)
provide the equivariant trace map, centralizer behavior, and translated
unipotent maps. The determinant identity used in (2) is classical and was
audited in the [preceding proof](IRREDUCIBLE_CENTRALIZER_REALIZATION.md).
Smoothness of the adjoint quotient at cyclic matrices and the regular
nilpotent orbit dimension are standard; the companion-matrix argument
above suffices in this setting. For the broader established context see
[Lehn--Namikawa--Sorger, introduction](https://doi.org/10.1112/S0010437X11005550).

The parity-dependent component and local multiplicity formulas are
deductions proved here from those inputs. Targeted literature searches
did not identify this particular formulation; no exhaustive review,
priority claim, or resolution of a named published open problem is asserted.

We have identified a component and its generic scheme structure, not
the whole fiber. At nonregular unipotent points the cyclic-coordinate
argument no longer applies: other components may meet C, and further
embedded or nonreduced structure is not ruled out. Nor does this
special fiber determine the generic degree or minimum Cayley degree.

**Next step:** determine the local fiber at the nonregular points of C,
starting with the semisimple point t. The aim is to decide whether the
generic simple/double structure extends over the nilpotent singularities
or whether extra components or embedded structure occur. A structural
description there would complete this component's geometry; another
list of generic-rank checks would not.
