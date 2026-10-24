# Degrees of symplectic polynomial-determinant generalized Cayley maps

Date: 22 October 2026. Status: uniform written deduction with exact finite
corroboration. External review and publication priority are not established.
The proof uses the classical Cayley chart, character differentiation and
the elementary degree of a rational function in one variable. This is a
tractable representation family, not a general degree algorithm.

## 1. Representations, normalization and theorem

Work over C. Let G=Sp_(2n), n>=1, acting on V=C^(2n), and put

\[
 E=\Lambda^\bullet V,\qquad D(g)=\det(I+g)=\chi_E(g),\qquad K=4^n.
\]

Let F in Z[d] have degree q>=1, and assume F(D) is the character of an
actual nontrivial rational representation rho_F. Nonnegative coefficients
of F suffice, but are not required: for example (D-1)^q is the character
of (E minus one trivial summand) to the q-th tensor power. The representation
is allowed to be reducible and to have mixed central characters.
No virtual representation with negative irreducible multiplicities is used.
There is no parabolic or Borel quotient in this calculation.

Use the normalized trace form and map

\[
 B_F(A,H)=\operatorname{tr}(d\rho_F(A)d\rho_F(H)),\qquad
 B_F(\theta_F(g),H)=\operatorname{tr}(\rho_F(g)d\rho_F(H)).
\]

**Theorem.** The full generalized Cayley map theta_F:G -> g is surjective,
and its generic degree is

\[
 \boxed{\deg(\theta_F)=2nq=2n\deg F.}                       \tag{1}
\]

Over the open subset of invertible matrices in g=sp_(2n), it is finite
flat of degree 2nq. The torus and adjoint-quotient maps have the same
generic degree. In particular, for rho_m=E^(tensor m),

\[
 \boxed{\deg(\theta_{\rho_m})=2nm\quad(n,m\ge1).}          \tag{2}
\]

The degree in (1) is that of this specified representation-induced map.
It is not the minimum Cayley degree of G, and does not establish a
minimum over all representations. The finite-flat assertion is restricted
to the stated open target; global finiteness is not asserted.

## 2. The normalized map from the character

Since G is semisimple, dD at identity is zero. Differentiating F(D)
twice there gives B_F=F'(K) B_E. Both actual nontrivial trace forms are
positive multiples of a fixed invariant form, so F'(K)>0.
Differentiation at an arbitrary g gives

\[
 \theta_F(g)=\frac{F'(D(g))}{F'(K)}\theta_E(g).
\]

The previously checked [exterior trace formula](SYMPLECTIC_HARMONIC_LOCALIZATION.md)
therefore yields the polynomial identity

\[
 \theta_F(g)=\frac{2F'(D(g))}{K F'(K)}
                    (g-I)\operatorname{adj}(g+I).            \tag{3}
\]

On G_D, introduce classical Cayley coordinates

\[
 X=(g-I)(g+I)^{-1},\quad
 g=(I+X)(I-X)^{-1},\quad
 \delta(X)=\det(I-X),\quad D(g)=K/\delta(X).
\]

This is an equivariant isomorphism G_D -> g_delta. Equation (3) becomes

\[
 \theta_F(g)=\frac{2F'(K/\delta(X))}{F'(K)\delta(X)}X.       \tag{4}
\]

For F=d^m, it reduces to 2 delta^(-m) X, with all normalization factors
cancelled. Copies of a representation and added trivial summands do not
change its normalized map.

## 3. A radial rational-function proof of the generic degree

We use an elementary observation. For a rational map on a vector space
of the form X -> u(X)X, the ratios of its coordinates are unchanged.
On a chart write X=xZ, with one linear coordinate of Z equal to one.
The other coordinates of Z generate a purely transcendental field L.
The source function field is L(x); the target subfield is L(r(x)), where
r(x)=x u(xZ). If r=A/B is nonconstant with coprime polynomials A,B in
L[x], then

\[
 [L(x):L(r(x))]=\max(\deg A,\deg B).
\]

For completeness, A(T)-yB(T) is irreducible over L(y): it is primitive
and linear in y, and a factor independent of y would divide both A and B.
Its degree in T is max(deg A,deg B). This proves the observation without
assuming irreducibility of a numerical fiber specialization.

Apply it to (4). Write F(d)=sum_(k=0)^q a_k d^k, let ell>=1 be the
smallest k with a_k!=0, and put p(x)=delta(xZ). On a generic direction,
det Z!=0, so p has degree 2n and p(0)=1. Define

\[
 H(p)=\sum_{k=1}^q k a_k K^{k-1}p^{q-k}.
\]

Then

\[
 r(x)=\frac{2x H(p(x))}{F'(K)p(x)^q}.
\]

The numerator and denominator are coprime: x is coprime to p because
p(0)=1, and H(0)=q a_q K^(q-1)!=0. Their degrees are respectively
1+2n(q-ell) and 2nq. Since 2n ell>1, the latter is larger.
Thus the function-field degree is exactly 2nq. Cayley coordinates are
birational coordinates on the group, so this is the full group degree.

For the torus the same proof works in its diagonal Cayley coordinates.
For the quotient, apply Kostant--Michor's Corollary 3.3, or cancel the
two Weyl field degrees in the torus/quotient diagram. There is no extra
division by the Weyl group order.

## 4. An explicit fiber polynomial, surjectivity and a finite-flat open set

Given Y in sp_(2n), put p_Y(s)=det(I-sY) and

\[
 P_Y(s)=F'(K)p_Y(s)^q
       -2s\sum_{k=1}^q k a_k K^{k-1}p_Y(s)^{q-k}.           \tag{5}
\]

Every root s of P_Y has s!=0, since P_Y(0)=F'(K). It also has p_Y(s)!=0:
at a zero of p_Y the second term is -2s q a_q K^(q-1), which is nonzero.
Consequently

\[
 g=(I+sY)(I-sY)^{-1}                                      \tag{6}
\]

is a well-defined symplectic matrix, and (4)-(5) give theta_F(g)=Y.
When Y!=0, every preimage in G_D arises this way: the scalar in (4)
cannot vanish, so X=sY for a unique s.

The eigenvalues of Y occur in opposite pairs. If p_Y is nonconstant
of degree 2r, the same leading-term argument gives deg P_Y=2rq>=2.
If p_Y=1 (Y is nilpotent), then P_Y=F'(K)(1-2s).
Thus P_Y is always nonconstant and has a complex root, proving
surjectivity. This existence argument also applies to Y=0. For general
F, it does not assert that (5) describes the entire zero fiber: zeros
of F'(D) can give additional zero-image loci inside G_D.

For invertible Y, P_Y has degree 2nq with leading coefficient
F'(K)(det Y)^q, a unit on that target open set. There are no preimages
on D=0: equation (3) gives an image of rank at most one there, whereas
an invertible 2n-by-2n matrix has rank 2n.
Scheme-theoretically, the adjugate determinant identity makes det(theta_F)
divisible by D^(2n-1); thus D is a unit after inverting det(theta_F).

In fact the preimage of this target open set has coordinate ring

\[
 \mathbb C[\mathfrak g][(\det Y)^{-1}][s]/(P_Y(s)).          \tag{7}
\]

To check this scheme-theoretically, p_Y is a unit modulo P_Y. Indeed,
H(p)=H(0)+p H_1(p) with H(0) a nonzero constant, so the ideal (P_Y,p_Y)
contains s; since p_Y(0)=1 it contains 1. Equation (5) then makes s and
the scalar u in (4) inverse units. Formula (6) defines an inverse map
to the map taking g to (theta_F(g),1/u(g)). These formulas are regular
on the stated opens and are mutual inverses over the coordinate rings.
After dividing P_Y by its unit leading coefficient, (7) is free with
basis 1,s,...,s^(2nq-1). This proves finite flatness and the fiber length
including multiplicities. In characteristic zero a generic fiber is reduced.

For F=d^m the polynomial is particularly simple, up to a constant unit:

\[
 \boxed{P_Y(s)=\det(I-sY)^m-2s.}                          \tag{8}
\]

This is also an explicit inversion procedure for a generic fiber.

## 5. Degree laws for representation operations in this family

Let rho_F and rho_J be two actual nontrivial representations of this
form, with q=deg F and r=deg J. Their characters under direct sum and
tensor product are (F+J)(D) and (FJ)(D). The leading coefficient of
each polynomial is positive: in torus coordinates, the weight
(q,...,q) occurs in F(D) with multiplicity a_q, since the corresponding
coefficient in D^q is one and lower powers cannot contain that weight.
Thus leading terms cannot cancel in F+J. Formula (1) implies

\[
 \begin{aligned}
 \deg\theta_{\rho_F\oplus\rho_J}
     &=\max(\deg\theta_{\rho_F},\deg\theta_{\rho_J}),\\
 \deg\theta_{\rho_F\otimes\rho_J}
     &=\deg\theta_{\rho_F}+\deg\theta_{\rho_J},\\
 \deg\theta_{\rho_F^{\otimes m}}&=m\deg\theta_{\rho_F}.
 \end{aligned}                                           \tag{9}
\]

These laws concern representation operations, not composition of maps.
They are proved only for characters polynomial in D. Most Sp
representations at n>=2 are outside this subring. At n=1, D=2+chi_V
generates the character ring, and the degree statement recovers the
elementary SL2 rule: twice the largest highest weight in an actual sum.

The special torus-regularity-preserving characters A D^m+B are a
subfamily. Formula (1) also covers polynomial mixtures which fail torus
regularity; degree and regularity preservation are different properties.

## 6. Independent computation, published inputs and limitations

The [exact verifier](tools/verify_symplectic_cayley_degree.py) and
[report](symplectic-cayley-degrees-verification.json) use three routes:
univariate fiber arithmetic, exterior-weight/tensor-convolution moment
equations evaluated in the entire fiber algebra, and small multivariable
Groebner fibers formed directly from those Laurent moments. The last
route saturates only the torus coordinates; it does not use Cayley
coordinates or assume that the boundary contributes no solutions.
The report records 80 scalar cases (ranks 1-8), 24 independent weight
fiber calculations, ten original Laurent Groebner fibers, 108
singular-target polynomial checks and nine full exterior-matrix trace
checks on Sp4 fibers. The matrix targets include an invertible
nonsemisimple element, a regular nilpotent and a singular semisimple element.
Counts are finite corroboration of Sections 3-4, not an exhaustive
proof over all ranks or representations.

[Kostant--Michor](https://arxiv.org/pdf/math/0109066) supplies the trace
construction, tensor-power identity (Theorem 2.8), and quotient-degree
comparison (Corollary 3.3). Their Section 7.15 already computes a Spin
degree by reducing a scalar-Cayley fiber to one determinant equation.
Our radial argument is closely related to that published method; the
SL2 degree check follows from their Section 6.1 formula. The classical
Cayley chart and one-variable rational degree are standard inputs.

The conclusions recorded here are the polynomial-determinant family
formula, its explicit fiber presentation, surjectivity, finite flatness
over invertible targets, and the restricted operation laws (9).
A targeted search on 22 October 2026 for generalized Cayley symplectic
degrees and exterior/tensor character degrees found no matching explicit
statement. That limited search does not establish novelty or present
openness in the literature. No claim about arbitrary representations,
minimum Cayley degree, integral Borel images or global finiteness follows.

The next degree question is whether the scalar-Cayley reduction gives
similarly uniform tensor-power laws for Kostant--Michor's Spin family,
with its double cover and central parity treated explicitly. The broad
exceptional boundary-correction problem remains open in this program.
