# Adjoint pencil scope and elimination of linear character constraints

Date: 2026-10-23. Written characteristic-zero arguments with independent
exact checks. This note both limits the preceding spectral strategy and
proves a further actual adjoint-fiber result. Novelty remains unverified.

## 1. The original-coordinate pencil is more restrictive than expected

Let G be simple complex simply connected of simply laced type, and let K
be a proper standard Levi whose derived factors are all of type A. Choose
the fundamental character x of an A1 factor. Assume at least one other
derived factor is present. Keep the original fundamental characters and
the full central-charge coefficients on [K,K] times T -> K.

Call the adjoint restriction an original-coordinate pair pencil if it is

    chi_Ad|K=x^2+b(z)x+q(z)+u^t(A(z)+x B(z))v,             (1)

after only partitioning/permuting the remaining active characters into
equally sized vectors u,v. Absent variables may be free. No translations,
rational eliminations, or changes mixing fundamental characters are
allowed in this definition.

**Theorem 1 (scope classification).** Under these hypotheses, (1) occurs
if and only if G is E6 and K is the Levi obtained by omitting Bourbaki
node 4. Its derived type is A2+A1+A2, and x is the node-2 character.
Consequently the noncollision and boundary hypotheses for this precise
adjoint pencil class hold by the already established E6 argument. This
does NOT classify pencils obtainable after arbitrary affine changes,
other representations, or Levis with D/E factors.

This is a limitation of the proposed method, not a new broad adjoint
smoothness theorem. The arbitrary-size spectral theorem remains valid,
but its raw-coordinate adjoint applications do not form an expanding
exceptional family.

### Proof from omitted-node incidence and the highest root

The [adjoint norm argument](BRANCHING_NORM_FIBERS.md#3-the-adjoint-cubic-obstruction-has-a-restricted-shape)
shows that every nonzero central-charge term is a positive coefficient
times a product of fundamental characters on distinct A-factors. Thus
these external terms cannot cancel an intrinsic square from another A1
factor. Equation (1) therefore permits exactly one A1 factor.

For an omitted simple root alpha_j, consider the adjoint central-charge
space of alpha_j. Its lowest Levi weight is the restriction of alpha_j:
minus one fundamental weight on each retained component adjacent to j.
No omitted node has two neighbors in a single retained component, since
the Dynkin diagram is a tree. Highest-lowest duality and the standard
irreducibility theorem for Levi central-root spaces give a product of
one fundamental character for each such component. In particular, if j
touches exactly ONE component other than the distinguished A1, this
charge supplies either a lone leaf variable or x times a leaf variable.
Both are forbidden by (1). Distinct central charges cannot cancel them.

Here is the resulting tree argument, valid across A and D in every rank.
Take a shortest path from the distinguished A1 component to any other
retained component. Its interior nodes are omitted. The last omitted
node on this path touches that other component. To avoid exactly one
noncentral neighbor, it must touch a second retained noncentral component.
It also has its edge back toward the A1. Thus it is a degree-three node.
The connected simply laced finite diagrams have at most one such node.
The two arms away from the A1 each begin with a retained component of
size at least two. Type A has no branch node. Type D has arm lengths
1,1,n-3, so two such arms are impossible.

The E6,E7,E8 arm lengths are 1,2,2; 1,2,3; 1,2,4. Hence the A1 lies
on the length-one arm, the branch node is omitted, and both longer arms
start retained. Each longer arm must be entirely retained: an omitted
tail touches only one noncentral component, while two retained components
of size at least two separated by a gap would require arm length at least
five. Thus the only candidates are the three maximal Levis obtained by
omitting node 4, of types A1+A2+A_(r-3), r=6,7,8.

The ambient highest root now distinguishes them. Its coefficients and
fundamental highest-weight label are:

| Type | Highest root in simple roots | Fundamental label | Coefficient at node 4 |
| --- | --- | --- | --- |
| E6 | (1,2,2,3,2,1) | omega_2 | 3 |
| E7 | (2,2,3,4,3,2,1) | omega_1 | 4 |
| E8 | (2,3,4,6,5,4,3,2) | omega_8 | 6 |

This vector is also highest for the Levi. In E7 and E8 its nonzero
central-charge constituent is a SINGLE fundamental character of a
noncentral factor. The norm argument makes that character monomial
literal, with no lower-character cancellation. This violates (1).
In E6 the same term is the distinguished A1 character, which is allowed;
the [verified E6 formula](SPECTRAL_BRANCHING_FIBERS.md#4-recover-the-e6-suspension-and-check-its-boundary)
has (1) and completes the converse. This is a proof from diagram structure,
not a conclusion drawn from a rank sweep.

## 2. Linear character terms can impose useful constraints

Excluding linear terms was therefore too restrictive. The following
elementary scheme lemma allows us to use them without discarding equations.

**Lemma 2 (multiplier elimination).** Over a Laurent coefficient ring,
let a,t be affine variables, y the remaining variables, and suppose

    chi=t*(a-f(y,z))+G(a,y,z),                              (2)

where G is independent of t. The entire critical scheme of
chi-sum c_j log z_j is isomorphic to that of

    G(f(y,z),y,z)-sum c_j log z_j,
    a=f(y,z), t=-G_a(f(y,z),y,z).                          (3)

The same holds when the first term is t*lambda(z)*(a-f), with lambda a
Laurent unit, after rescaling t. Additional unused coordinates are retained.

Proof. The t derivative sets a=f and the a derivative determines t.
The remaining y derivatives become G_y+G_a f_y; the central logarithmic
derivatives become D_z G+G_a D_z f. These are exactly the derivatives of
the substituted potential in (3), in the quotient ring itself. Thus this
is an ideal elimination, including multiplicities. Substituting a=f but
forgetting the equation for t or the chain-rule terms would be incorrect.

**Lemma 3 (anchored product transfer).** Let delta be a nonzero constant
and S a Laurent polynomial on a positive-dimensional torus. For

    chi=sum a_i b_i+(x^2-delta)*(y^2+S(z)y)+constant,       (4)

the critical scheme over ANY nonzero target vector c reduces exactly to

    x=0, y=-S/2, a_i=b_i=0,
    D_j(delta*S^2/4)=c_j for every j.                     (5)

Proof. After eliminating a_i,b_i, a nonzero target coordinate satisfies
(x^2-delta)*y*D_j S=c_j. Thus x^2-delta and y are units in the localized
quotient. The y derivative implies 2y+S=0, so y^2+Sy=-y^2 is a unit.
The x derivative forces x=0. Substituting gives (5), and the converse is
immediate. This uses no radicals and does not remove a possible rank-drop
component without proof. It asserts nothing about c=0.

Generic smoothness and nonemptiness now reduce to the explicit torus
gradient map in (5). This extends the earlier product argument to a
factor with a nonzero constant shift and no separate torus target block.
It is an arbitrary-Laurent-family lemma; the application below is one
actual simple-group case.

## 3. An E7 adjoint theorem from the forced linear term

Let G be simply connected E7, rho=Ad, K the Levi omitting node 4, and
normalize theta by B_Ad(theta(g),X)=tr(Ad(g) ad(X)). In simple-coroot
coordinates put

    H=(4,6,8,12,9,6,3), h=tH, c=B_Ad(h,H)=432t.          (6)

Here dim K=27, rank G=7. For t!=0 the target centralizer is exactly K.
On the cover SL3 times SL2 times SL4 times C*, write (a,b) for the
SL3 defining and dual characters, x for the SL2 defining character,
and (w,v,u) for the three SL4 fundamental characters IN THAT ORDER.
The central cocharacter parameter is z.

**Theorem 4.** For c!=0 the quotient-critical scheme on the actual Levi
quotient is

    Spec C[W,W^(-1)]/(3W^2-cW-3).                        (7)

For c!=0 and c^2+36!=0, theta_Ad^(-1)(h) intersect G_reg is the disjoint
union of two smooth regular K-orbits, each of dimension 20. In particular
the generic regular fiber is nonempty, smooth and has zero excess.
At c^2+36=0 the quotient scheme is a double point. This describes the
regular-source part, not all nonregular points of the whole group fiber,
and the number two is NOT the generic degree of E7 -> e7.

### Character identity and exact reduction

The positive charge pieces of the adjoint restriction are:

| Charge | Character | Dimension |
| --- | --- | --- |
| 1 | a x u | 24 |
| 2 | b v | 18 |
| 3 | x w | 8 |
| 4 | a | 3 |

Their negative pieces are dual. The charge-zero part is
ab+uw+x^2-2. The dimensions sum to 27+2*(24+18+8+3)=133. These are
the usual Levi central-root modules: their highest roots are dominant
on the retained nodes, and the displayed fundamental-module dimensions
exhaust the adjoint dimension. Both verification routes independently
reconstruct every root weight, retaining the central charges. Thus

    chi=ab+uw+x^2-2+x*(z*a*u+z^(-1)*b*w)
        +v*(z^2*b+z^(-2)*a)+x*(z^3*w+z^(-3)*u)
        +z^4*a+z^(-4)*b.                                 (8)

The middle SL4 character v is a multiplier. Its equation gives
a=-z^4 b, and the a equation gives

    v=-z^2*(b+x*z*u+z^4).                                 (9)

Lemma 2 eliminates these two variables. In the remaining potential set

    u'=u+x*(z^(-1)*b+z^3),
    w'=w+x*(-z^5*b+z^(-3)),
    y=z^2*b, S=z^6-z^(-6).

These are invertible changes over C[z,z^(-1)], and exact expansion yields

    chi_reduced=u'w'+(x^2-1)*(y^2+Sy)-2.                 (10)

Lemma 3 proves, for every c!=0,

    x=u'=w'=0, y=-S/2, 3*(z^12-z^(-12))=c.              (11)

In the original coordinates this is

    x=u=w=0, a=z^2*S/2, b=-S/(2z^2),
    v=-(z^6+z^(-6))/2.                                   (12)

All transformations are defined on the entire central torus. In
particular there is no unproved spectral-discriminant boundary left
over. The cancellation is justified by the nonzero TARGET equation,
not by declaring S or a source determinant generically nonzero.

### Central descent and group fibers

The six retained simple-coroot columns together with H have determinant
12, and H is primitive. The multiplication cover has kernel of order 12,
projecting injectively to the central z-coordinate, hence onto mu12.
After (11)-(12) there is precisely one set of affine coordinates per z,
and the kernel action preserves the critical scheme. Its invariant
coordinate is W=z^12. Descending (11) gives (7) scheme-theoretically.
The quadratic discriminant is c^2+36, and its constant term is nonzero.

For each quotient point there is one regular K-conjugacy class of lifts.
The smooth regular quotient has relative dimension dim K-rank G=20.
Equivariance forces Z_G(g) to lie in K when theta(g)=h, so a regular
K-lift is G-regular. The regular-quotient critical presentation and finite
etale central descent therefore prove all assertions of Theorem 4.
At the double-point targets, the actual regular fiber retains the
corresponding nonreduced structure. No assertion about t=0 follows.

## 4. Checks, known inputs, and what remains

The [Sage verifier](tools/verify_adjoint_pencil_scope.py) checks the
exceptional diagram candidates, the highest-root obstructions, all
charge records, the exact coordinate identities, the original central
derivative and the universal nonzero-target ideal in Lemma 3. A direct
adjoint-root evaluation over F_1033, independent of the character
elimination, gives differential rank 113 and kernel dimension 20 at a
regular semisimple source. This finite-field point corroborates the
characteristic-zero proof; it does not establish genericity by sampling.

The [native audit](tools/audit_adjoint_pencil_scope.py), without Sage,
uses Cartan reflections and subset-weight models to reconstruct the
characters, checks the scalar equation and central mu12 action, and
recounts the root differential. Reports:
[verification](adjoint-pencil-scope-verification.json),
[audit](adjoint-pencil-scope-audit.json).

Published inputs include [Kostant's Levi-root-space theorem](https://arxiv.org/pdf/0711.2809),
the ADE diagram and highest-root classification, and the trace
construction of [Kostant--Michor](https://www.mat.univie.ac.at/~michor/kostant.pdf).
The root-data convention can be checked against the
[Sage reference](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/root_system.html).
Regular quotient geometry and the chain-rule eliminations are standard
tools. The scope classification and fiber deductions are written here
without a priority claim. Targeted source checks on 2026-10-23 did not
establish novelty or identify these as a named published open problem.

The next structural goal is to treat the linear character terms as a
system of constraints before asking for a spectral normal form. The
highest-root argument forces such terms in E7 and E8, so this is a
root-theoretically motivated change of method. A useful next theorem
would give an elimination criterion for several coupled multipliers,
with the actual E8 omission-of-node-4 character as the discriminating
test. E8, arbitrary affine reparametrizations, other Levi embeddings and
the all-simple-adjoint smoothness conjecture remain unresolved here.
