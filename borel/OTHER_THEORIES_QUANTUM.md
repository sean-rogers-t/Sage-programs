# A quantum projective-space benchmark and the affine transport route

Exploration, 22 October 2026. We specify a quantum lift of the literal
Springer source, prove its complete image for an arbitrary-rank family,
and separate this from the existing Peterson correspondences. The ambient
quantum rings and correspondences are published mathematics. The image
formula below is an elementary deduction; no priority claim is made.

## 1. Why a quantum lift needs a source

Let G be finite, simply connected and simple, and let rho have faithful
differential. On G/B the literal source is Sym_Z Lambda_rho, a polynomial
ring on a lattice basis. Assign each linear element its ordinary divisor
class inside the Schubert module underlying QH^*(G/B;Z), then use quantum
multiplication to extend. This gives a canonical Z[q]-algebra map from
Sym_Z(Lambda_rho)[q] and recovers the ordinary Springer/Borel map at q=0.
It does not assert that the ordinary kernel remains a kernel, or that
the image specializes injectively. On a partial flag, the Levi source
can need higher generators and have relations, so those lifts require
separate specification and verification.

For example on SL2/P1, z=2omega_1 maps to 2h and z^2 maps to 4q. The old
relation z^2=0 no longer holds. Thus a ring map from the old ordinary
image with unchanged relations is already impossible in this first case.

## 2. The specified projective-space lift

Take G=SL_m(C), rho the defining module, and P=P1, so G/P=CP^(m-1).
The literal linear source is Q. A permutation basis for Q is
r_j=epsilon_1-epsilon_j, 2<=j<=m. The full Levi source is the polynomial
ring A=Z[E_1,...,E_(m-1)], where E_j=e_j(r_2,...,r_m).
Its ordinary images and complete ordinary lattice are proved in
[the projective calculation](LOCAL_IMAGE_BOUND.md):

\[
 \beta(E_j)=\binom mj h^j,\qquad
 J_H^r=g_r\mathbb Z h^r,\quad
 g_r=\frac{m}{\gcd(m,r)}\quad(0\le r<m),\quad g_0=1.              \tag{1}
\]

The published quantum presentation is

\[
       QH^*(CP^{m-1};\mathbb Z)=\mathbb Z[q,h]/(h^m-q),
       \qquad \deg h=1,\quad\deg q=m.                             \tag{2}
\]

For background and the projective-space example, see
[Fulton--Pandharipande, Notes on stable maps and quantum cohomology, Section 8](https://people.math.ethz.ch/~rahul/sm.pdf).
Our lift is E_j -> binomial(m,j) h^j. It is the unique graded lift of
these chosen generator images: j<m, so there is no positive-q term in
their degrees. Since A is polynomial, this assignment defines the map
without an additional relation-lifting assumption.

**Theorem.** The complete image, as an embedded Z[q]-module, is

\[
             \boxed{J_q=\bigoplus_{r=0}^{m-1}g_r\mathbb Z[q]h^r.} \tag{3}
\]

It is free over Z[q], specializes faithfully to J_H at q=0, and has
ambient quotient direct-sum_(r=0)^(m-1) (Z/g_r)[q]. For the basis
y_r=g_r h^r, y_0=1, multiplication is

\[
 y_a y_b=\frac{g_a g_b}{g_{(a+b)\bmod m}}
               q^{\lfloor(a+b)/m\rfloor}y_{(a+b)\bmod m}.         \tag{4}
\]

Thus (4), for 1<=a,b<m and with y_0=1, is also a complete presentation
using these image generators: it reduces every monomial to their free
Z[q]-span, and that span embeds independently in (2).

### Proof

The ordinary image computation supplies an integral polynomial in the
E_j of each degree r<m with image g_r h^r. No quantum reduction occurs
in those degrees, so the module on the right of (3) is contained in J_q.
Conversely g_(a+b mod m) divides g_a*g_b. To see this prime by prime,
put e=v_p(m); the exponent in g_r is max(0,e-v_p(r)), with the
convention v_p(0)=infinity. The inequality
v_p(a+b)>=min(v_p(a),v_p(b)) proves the required divisibility, including
reduction modulo m. Therefore the proposed module is a subalgebra.
It contains all generator images in (1) and hence contains J_q. This
proves equality. Independence follows from the basis 1,h,...,h^(m-1)
of (2); the other assertions follow immediately. QED.

There is a useful intrinsic description. The Euler derivation
delta(h)=h, delta(q)=mq preserves h^m-q, and

\[
                   J_q=\{f:\delta f\in m\,QH^*(CP^{m-1};\mathbb Z)\}.
                                                                    \tag{5}
\]

Indeed a normal-form term a_r(q)h^r has derivative
(m q a_r'(q)+r a_r(q))h^r. Divisibility by m is exactly divisibility
of a_r by g_r. Formula (5) is a congruence description of this one
specified image, not a claimed derivation criterion for every flag.
Inverting q does not remove the integer torsion in its ambient quotient.

## 3. What Peterson transport supplies

Lam--Shimozono's integral equivariant Peterson theorem identifies a
specified localized quotient of affine-Grassmannian Pontryagin homology
with localized quantum cohomology of G/P:
[Theorem 10.16](https://arxiv.org/pdf/0705.1386).
Applying that published isomorphism to the specified image (3) gives a
corresponding affine homology subalgebra after localization. This is a
consequence of an existing theorem; it does not construct a new
Springer morphism on an affine group. Recovering a preferred unlocalized
homology lattice requires another argument.

The K-theoretic Peterson correspondence is also established; it should
not be labelled an outstanding conjecture based on older paper titles.
See [Kato, Annals 2025](https://doi.org/10.4007/annals.2025.202.3.2)
and [Chow--Leung's alternative proof, including the parabolic formulation](https://arxiv.org/abs/2201.12951).
A prospective quantum K-image still needs its own specified source map
and a check of its q=0 specialization.

For affine flags themselves there are already rigorously constructed
quantum products on particular divisor/evaluation subrings, including
the affine quantum deformation of finite-flag cohomology in
[Mare--Mihalcea](https://arxiv.org/abs/1409.3587).
That result should not be described as a construction of every possible
full quantum cohomology ring of an ind-variety.

## 4. Checks and the actual next gap

The [independent integer recurrence](experiments/other_theories_quantum_pilot.py)
enumerates the gcd of coefficients of all actual source monomials in
each degree, then adds powers of q. It checks m=2,...,80 through degree
8m against (3), and separately checks all multiplication constants in
(4). [Report](experiments/other-theories-quantum-pilot.json).
These corroborate the proof and do not replace it.

The next substantial test is an integral Grassmannian or exceptional
partial flag whose source has relations and generators in degrees at
least the smallest quantum parameter degree. There the graded lift may
not be unique, and image flatness does not follow from (3). The questions
are which lift is geometrically justified, whether specialization has a
kernel, and whether its torsion can be controlled before localization.

## 5. The defining symplectic source on the same projective space

For G=Sp_(2r), P=P1 and the defining representation, use the literal
source Z[z,q_1,...,q_(r-1)] with z=2epsilon_1. Write Q for the quantum
parameter to distinguish it from these source variables. The generator
degrees 1,2,...,2r-2 are all smaller than deg Q=2r. The unique graded
lift therefore sends z to 2h and q_j to (-4h^2)^j in
Z[Q,h]/(h^(2r)-Q). Its entire image and complete source kernel are

\[
 J_Q=\mathbb Z[Q,z]/(z^{2r}-2^{2r}Q),\qquad
 \ker=(q_j-(-z^2)^j,\ z^{2r}-2^{2r}Q).
\]

Indeed reduction by these relations leaves the independent basis
1,z,...,z^(2r-1), whose target images are 1,2h,...,2^(2r-1)h^(2r-1).
This proves completeness and faithful Q=0 specialization. The image is
flat over Z[Q], with its old diagonal ambient cokernel. It is different
from the defining SL_(2r) image (3) on the same projective space, because
the group and literal Springer source are different. No compatibility
between those two trace projections is implicitly assumed.
