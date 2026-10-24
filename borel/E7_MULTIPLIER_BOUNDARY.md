# E7 Weyl families and a multiplier with a controlled boundary

Characteristic zero; 2026-10-23. Written deductions, with independent exact
checks. No literature-priority claim. This concerns the adjoint Springer
map of the simple complex simply connected group E7, not arbitrary
representations or its full-map generic degree.

## 1. Precise progress

Normalize the actual map by

    B_Ad(theta_Ad(g),X)=tr(Ad(g) ad(X)),
    B_Ad(Y,X)=tr(ad(Y) ad(X)).

For a Levi K and h in z(k) with Z_G(h)=K, write

    F_h^reg=theta_Ad^(-1)(h) intersect G_reg.

**Theorem 1.** If K has derived type A3+A1+A1, there is a nonempty
Zariski-open subset of its central targets on which F_h^reg is nonempty,
smooth, and pure of dimension 16. There are exactly three standard
embeddings of this type in E7, and they are Weyl conjugate. On a finite
central cover, the generic critical quotient is the logarithmic-gradient
fiber of the single rational torus function (6) below.

This completes one additional Weyl family, not all E7 cases. The complete
grouping of the 32 cubic all-A standard Levis is:

| Derived type | Representative retained nodes | Standard embeddings | Full root-subsystem orbit size | Status here |
| --- | --- | --- | --- | --- |
| A1^3 | 125 | 10 | 3780 | Open |
| A2+A1^2 | 1235 | 12 | 15120 | Open |
| A2^2+A1 | 12356 | 3 | 10080 | Open |
| A1^4 | 1257 | 2 | 3780 | Open |
| A2+A1^3 | 12357 | 1 | 5040 | Open |
| A3+A1^2 | 12457 | 3 | 7560 | Theorem 1 |
| A3+A2+A1 | 123567 | 1 | 5040 | Earlier node-4 theorem |

Node strings use Bourbaki numbering. Distinct rows have different root
system types; each row is one Weyl class, as established by explicit words
on the complete root subsystems. Orbit sizes count all conjugate root
subsystems, not standard Levis or Springer orbits.

Together with the 88 quadratic cases and seven non-A emptiness cases,
generic smoothness or emptiness is now proved for 99 of the 127 proper
standard E7 Levis. The remaining 28 form the five open rows above.

Every cubic class has a representative inside the retained-node subsystem
123567. For A3+A1^2 use 12567; the others in the table already lie there.
Consequently all seven are Levi refinements of one known maximal Levi.
This is a route for comparing their characters, not a theorem that
generic smoothness descends to refinements.

## 2. A general transfer with a nonunit multiplier

Let T=(C*)^s, s>=2. Let A,B be Laurent polynomials, A not identically zero,
Q(x,z) a polynomial in x with Laurent coefficients, and C a constant. Set

    psi=(1-y^2)Q(x,z)+v(A(z)x+B(z))+C.                       (1)

Use the actual relative-critical scheme X=V(psi_x,psi_y,psi_v), with
mu_j=D_j psi and D_j=z_j partial/partial z_j. On U={A!=0} put

    x_0=-B/A,  K(z)=Q(x_0,z),
    v_0=-(1-y^2)Q_x(x_0,z)/A.

**Theorem 2.** Suppose K is not identically zero, every irreducible
component of {K=0} in U is an open subset of a translated codimension-one
subtorus of T, and the target image of X intersect {A=0} is nondense.
Then the generic fibers of mu are smooth and finite, possibly empty.
More precisely, outside a proper closed target subset they are exactly

    x=x_0, v=v_0, y=0, D K=c, A K!=0.                       (2)

They are nonempty if D K is dominant. An invertible logarithmic Hessian
of K at one point of {A K!=0} is a sufficient dominance certificate.
A sufficient boundary condition is dim(X intersect {A=0})<s; this uses
the underlying topological dimension, not reducedness of its equations.

**Proof.** Invert A only temporarily. The equation psi_v=0 solves for x,
and psi_x=0 solves for v, over the original coefficient ring localized at
A. Substitution leaves the exact scheme equation -2Ky=0. By the chain
rule the original central derivatives become (1-y^2)D K. Terms arising
from x_0 and v_0 are multiplied by the vanished relative derivatives.
This verifies scheme-level elimination and retains the target map.

On {K=0}, K vanishes identically along the corresponding translated
subtorus, so D K annihilates its fixed tangent Lie algebra. Hence mu lies
in the same fixed one-dimensional annihilator. This includes y=0,
intersections of components, and multiple components of the zero divisor.
There are finitely many such lines, all proper since s>=2. Exclude their
closures and the assumed nondense boundary image. Then K and A are units
on the target fibers, and -2Ky=0 gives y=0 as an equality of schemes.

The surviving problem is the morphism D K from the smooth open torus
{A K!=0} of dimension s to A^s. If nondominant its generic fibers are
empty. If dominant, characteristic-zero generic smoothness gives smooth
zero-dimensional fibers on a dense target open. Each such finite-type
zero-dimensional fiber is finite. Dominance also ensures nonemptiness
after shrinking the target open. The logarithmic Hessian is the
differential of D K. This proves the assertions.

This criterion allows a nonunit determinant, but only after its entire
boundary has been controlled. It does not assert that all quartics or all
linear multipliers have generically smooth fibers. The theorem concerns
generic fibers; it does not claim that D K is finite over the entire target.

## 3. The E7 character and exact reduction

Use retained nodes 12457. On the finite cover, x and y are the A1
characters at nodes 1 and 7; a,v,b are the A3 characters along 2,4,5.
For omitted nodes 3 and 6 the primitive central cocharacters are

    H_r=(3,4,6,8,6,4,2), H_s=(2,3,4,6,5,4,2)

in simple-coroot coordinates. Put R=r+r^(-1) and
S=r^2 s^2+r^(-2)s^(-2). The complete restricted adjoint character is

    chi=ab+a F+b G+H,
    F=y(r^(-1)s^(-1)x+s(1+r^2)),
    G=y(rsx+s^(-1)(1+r^(-2))),
    H=x^2+y^2+RSx+v(Rx+S)+R^2-3.                            (3)

This is checked against all 126 roots and seven zero weights; it is not
a selected part of the character. Complete the product ab+aF+bG.
The equations in a,b give a=-G, b=-F globally, without localization.
The exact residual potential is

    psi=(1-y^2)(x^2+RSx+R^2-1)+v(Rx+S)-2.                  (4)

Thus Theorem 2 applies with A=R, B=S. On R!=0,

    x=-S/R,
    v=-(1-y^2)(2x+RS)/R,
    K=(R^2-1)(R^2-S^2)/R^2.                                (5)

For the geometry make the finite etale torus change t=r^2s^2, so

    R=r+r^(-1), S=t+t^(-1),
    K(r,t)=((r+r^(-1))^2-1)
           *((r+r^(-1))^2-(t+t^(-1))^2)/(r+r^(-1))^2.       (6)

The map (r,s)->(r,t) has degree two. It is used to simplify the equations,
not to identify the quotient by the entire finite central kernel. No
Springer orbit count follows from this extra cover.

The zero divisor of K on R!=0 consists of translated subtori:
R=+/-1 gives finitely many constant r, and S=+/-R gives
t=+/-r or t=+/-r^(-1). Their target images lie in the fixed normal
lines described in Theorem 2, so none contributes to a generic target.

We must also retain R=0. The v-equation then forces S=0. In (r,t)
coordinates these are the four base points r^2=t^2=-1. Over any one,
the relative equations in x,y reduce to

    x(1-y^2)=0, y(x^2-1)=0,                                 (7)

with v free. Their support consists of (0,0) and the four points
(+/-1,+/-1). The entire boundary therefore has dimension one, whereas
the central target dimension is two. Its target-image closure is proper.
Any scheme structure there is excluded together with its support; no
radical replaces the actual target-fiber ideal in the argument.

Finally, at r=2,t=3, K=-49/12 and the logarithmic Hessian is

    [ 3889/225   -256/75 ]
    [ -256/75  -2296/75 ],  determinant=-9125752/16875 !=0.  (8)

This proves dominance and hence nonemptiness as well as smoothness.
It is a quotient dominance witness, not an assertion about a particular
regular source point at its target. The regular-quotient theorem supplies
regular lifts for generic targets with centralizer exactly K. The smooth
orbit factor has dimension dim K-7=12+2+2=16. Finite central etale descent
and Weyl conjugacy prove Theorem 1 for all three standard embeddings.

Normalization is retained: for h=u H_r+w H_s the original target is

    (c_r,c_s)^t = [[216,144],[144,144]] (u,w)^t.              (9)

If d=(D_r K,D_t K) uses (r,t), then c_r=d_r+2d_t and c_s=2d_t.
The trace matrix in (9) is obtained independently by summing the root
pairings; no informal rescaling of theta is used.

## 4. Evidence, attribution, and the next gate

The [Sage verifier](tools/verify_e7_multiplier_boundary.py) enumerates the
full Weyl subsystem orbits, supplies all 32 conjugacy words, reconstructs
the seed character, and checks both eliminations and every original
central derivative. The [native audit](tools/audit_e7_multiplier_boundary.py)
uses independent Cartan-reflection roots, subset character models and
SymPy algebra. It verifies the census, every conjugacy word, the boundary
ideal, logarithmic Hessian, trace matrix, and exact residual identity.
See the [verification](e7-multiplier-boundary-verification.json) and
[audit](e7-multiplier-boundary-audit.json).

The trace-projection map and its centralizer restriction are published
inputs from [Kostant--Michor, sections 2.4 and 2.6](https://www.mat.univie.ac.at/~michor/kostant.pdf).
Regular conjugacy quotients, generic smoothness, the chain rule and image
dimension bounds are standard. The multiplier criterion and this E7
application are deductions written here. Rechecking that primary source
does not establish their novelty or a literature-wide classification.

The next useful gate is to analyze the A2+A1^3 representative 12357 as a
refinement of 123567, and determine whether a comparable boundary-control
argument handles its three-variable residual. A successful criterion
should then be tested on the A1^4 refinement. The purpose is to resolve
two related families by one mechanism. Similar degree-four bounds alone
would not count as a solution. Five cubic Weyl families remain open;
whole fibers, special targets, global regularity preservation, generic
degrees and all E8 Levis are separate questions.
