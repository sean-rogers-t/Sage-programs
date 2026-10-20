# Springer maps on adjoint quotients and their critical determinants

Date: 2026-10-20. These are exact characteristic-zero calculations for the
normalized trace-projection Springer morphism. They retain invariant functions
that disappear under the Borel map. The formulas and proofs below are repository
deductions from published inputs; no originality or external-review claim is made.

## Maps, coordinates, and published inputs

Let G be a simply connected simple complex group, T a maximal torus, W its
Weyl group, and P its weight lattice. For a nontrivial representation rho,
write theta_rho:G -> g for trace projection, identifying d rho(g) with g.
The normalization is d theta_rho at 1 = identity. Fundamental weights omega_i
give coordinates w_i on t; simple coroots give the dual basis on t itself.
The torus coordinates t_i=e^{omega_i} are multiplicative, not the w_i.

For weight multiplicities m_mu, Sean's explicit formula is

\[
 M_\rho=\sum_\mu m_\mu\mu\mu^t,\qquad
 \theta_\rho(t)=M_\rho^{-1}\sum_\mu m_\mu\mu e^\mu(t).
\]

This is Theorem 1 of [Rogers's preprint](https://arxiv.org/pdf/1701.01538v1).
The trace projection is the generalized Cayley map of
[Kostant--Michor](https://arxiv.org/pdf/math/0109066), Sections 2.1 and 2.4.
Their notation Q(X) means rational functions on the complex variety X, not
a choice of integral form over the rational numbers.

For homogeneous basic invariants f_1,...,f_r of Q[t]^W, their pullbacks
describe the map T/W -> t/W. On the simply connected group, fundamental
characters chi_i=character of V(omega_i) are polynomial coordinates on T/W.
The API checks Weyl invariance, fundamental degrees and a nonzero Jacobian
for supplied f_i. Algebraic independence and matching Hilbert series then
prove that they generate Q[t]^W. Characteristic-zero base change gives C.

Kostant--Michor Corollary 3.3 equates the generic degree on adjoint quotients
with the generic degree of G -> g. The torus map has this same degree:
for regular semisimple h, equivariance forces a preimage g of h to commute
with h, hence to lie in its maximal torus. Equivalently, the two Weyl-group
field extensions in the restriction diagram have equal order and cancel.

The Borel map beta^P kills the positive-degree G-invariants. Thus these
quotient calculations study additional geometry of theta, while the
[integral-source classification](IRREDUCIBLE_CHARACTER_LATTICES.md) determines
which polynomial inputs come from integral characters. Neither computation
alone determines the other.

## A general full-differential formula

In logarithmic torus coordinates put

\[
 J_T(t)=\left(t_j\frac{\partial\theta_i}{\partial t_j}\right)_{ij}
       =M_\rho^{-1}\sum_\mu m_\mu\mu\mu^t e^\mu(t).
\]

This is an r by r matrix. Kostant--Michor's Cayley determinant is the
dim(g) by dim(g) determinant using left-translated tangent vectors.

For a root alpha set a=e^alpha(t), h=theta(t). Equivariance gives

\[
 d\theta_t([X_\alpha,t])=[X_\alpha,h].
\]

Under left translation, [X_alpha,t] becomes (a^{-1}-1)X_alpha and the
right side is -alpha(h)X_alpha. The eigenvalue of the differential on that
root direction is therefore a alpha(h)/(a-1). Opposite roots contribute

\[
 a\left(\frac{\alpha(\theta(t))}{a-1}\right)^2.
\]

Consequently the full determinant restricted to T is

\[
 \boxed{\Psi_\rho|_T=\det J_T\prod_{\alpha>0}
 e^\alpha\left(\frac{\alpha(\theta)}{e^\alpha-1}\right)^2.}
\]

Although written rationally, each root-direction quotient is a Laurent
polynomial: it is a matrix coefficient of the regular left-translated
differential. The displayed formula first holds on the dense locus where
the denominators are nonzero and extends regularly. The implementation
checks every division exactly after removing Laurent monomial shifts.
This also handles nonprimitive roots in P, for example alpha=2 omega in A1.
The normalization gives Psi(1)=1. Its zero set is the critical locus in G;
it is not a computed branch locus in the target. Invariant functions depend
only on the semisimple part, so their restriction to T determines this locus.

This is a direct derivation from equivariance, not a new general theorem
attributed to Kostant--Michor. Their Theorem 3.2 uses this full determinant
to localize the harmonic-separation identity. The calculations here do not
prove that localization can be removed.

## Complete G2 quotient maps

Use Sage/Bourbaki numbering: omega_1 gives the 7-dimensional representation
and omega_2 the 14-dimensional adjoint. On t choose

\[
 q=3w_1^2-3w_1w_2+w_2^2,\qquad
 z=w_1^2(w_1-w_2)^2(2w_1-w_2)^2.
\]

Their degrees are 2 and 6. They are Weyl invariant and their Jacobian is
nonzero, so they are basic generators. On T/W put x=chi_1 and y=chi_2.
At the identity (x,y)=(7,14). Define

\[
 D=x^2+2x-4y-7,\qquad A=2x^2+13x+y+13.
\]

For the minimal representation the entire quotient map is

\[
 \boxed{\theta_7^*q=\frac{x^2-4x-y-7}{12},\qquad
 \theta_7^*z=\frac{DA^2}{46656}.}
\]

The three different determinants are

\[
 \det J_T=\frac{y-2}{12},\qquad
 \det\frac{\partial(\theta_7^*q,\theta_7^*z)}{\partial(x,y)}
 =-\frac{(y-2)(x+1)A}{10368},
\]
\[
 \boxed{\Psi_7=\frac{(y-2)(x+1)^2A^2}{35831808}.}
\]

Thus the reduced critical locus is the union of the inverse images in G
of y=2, x=-1, and A=0 in G//G. The latter two factors occur squared in
the determinant. Psi_7(7,14)=1 checks the scalar normalization.

For the adjoint representation define

\[
 H=-9x^3-5x^2+26xy+21x+30y+25,
\]
\[
 L=x^3+3x^2-4xy+2y^2-5x+4y+1,
\]
\[
 S=9x^3-27x^2y-11x^2-18xy+18y^2+23x+17y-13.
\]

Then

\[
 \boxed{\theta_{14}^*q=
 \frac{-3x^3-5x^2+6xy+3y^2+11x+2y-7}{192},\qquad
 \theta_{14}^*z=\frac{DS^2}{191102976}.}
\]
\[
 \det J_T=-H/192,\quad
 \det\frac{\partial(\theta_{14}^*q,\theta_{14}^*z)}{\partial(x,y)}
 =-\frac{HLS}{679477248},\quad
 \boxed{\Psi_{14}=-\frac{HL^2S^2}{9618527719784448}.}
\]

Again Psi_14(7,14)=1. Both cases are recorded, expanded and factored, in
[the exact calculation report](experiments/springer-geometry-results.json).

## Generic degrees: proofs and independent computations

The two generic degrees are respectively **6 and 18**.

For the minimal quotient set theta_7^*q=a and theta_7^*z=b with a,b
algebraically independent. Eliminate y=x^2-4x-7-12a. The remaining equation is

\[
 (-x^2+6x+7+16a)((x+1)(x+2)-4a)^2-1728b=0.
\]

It has degree six in x. Over Q(a), a nonconstant polynomial F(x)-1728b
is irreducible in Q(a)[x,b]: its quotient ring is Q(a)[x], an integral domain.
It is primitive as a polynomial in x over Q(a)[b] because its leading
coefficient is a nonzero constant. Gauss's lemma gives irreducibility over
Q(a,b). Thus the function-field extension has degree exactly six.

There is a second proof covering both representations. Multiply theta(t)-c
by M to obtain the moment equations. The minimal nonzero weights are the
six vertices of a hexagon of normalized area 6 in the weight lattice. For
the adjoint, the twelve roots have a long-root hexagonal convex hull of
normalized area 18, with the short roots strictly inside it. Zero lies
strictly inside both polygons. Every edge contains exactly two support
weights, and these two vectors have nonzero determinant (absolute values
1 and 3 respectively). The two vector face equations on an edge are
therefore an invertible linear system in two nonzero Laurent monomials;
they have no torus solution. At a vertex the same assertion is immediate.
The target constants affect no proper face. Bernstein's nondegeneracy
criterion therefore gives the normalized areas as the generic fiber lengths.
The common-support formulation permits individual coordinate equations to
have smaller supports; see [Mondal, Theorems VII.5, VII.7 and Corollary
VII.10](https://arxiv.org/pdf/1806.05346). Dominance and generic finiteness
also follow from d theta_1=id. No genericity of the representation's fixed
weight coefficients is being assumed.

The audit separately computes Groebner bases over Q(a,b) for the quotient
fibers and obtains dimensions 6 and 18. It checks every edge of each Newton
polygon and independently constructs the full differential as follows:

1. Build the 7-dimensional module from its simple-root sl2 strings.
2. Verify all Chevalley and Serre relations with rational matrices.
3. Generate a 14-dimensional Lie basis by matrix commutators; construct its
   adjoint representation by solving for commutator coordinates.
4. For each representation, form B_ij=tr(X_j X_i) and
   A_ij(t)=tr(rho(t) X_j X_i). The differential matrix is B^{-1}A.
5. Compare its full determinant coefficientwise with the root-factor formula,
   and its Cartan block with J_T. Both comparisons pass in both cases.

See [the independent audit](experiments/audit_springer_geometry.py) and
[its report](experiments/springer-geometry-audit.json). The representation
matrix construction does not use WeylCharacterRing or the root-factor formula.
The generic-degree proof is specific to these polygons; no automatic
weight-polytope-volume formula is claimed in higher rank.

The degree here belongs to a **specified trace map**. The group invariant
Cay(G2), which minimizes over all dominant equivariant rational maps, is
already known to equal 2; see [Lemire--Popov--Reichstein, Example
10.30](https://www.math.uni-bielefeld.de/LAG/man/150.pdf).
[Their later paper](https://arxiv.org/pdf/math/0608473) also discusses a
degree-at-most-six construction for G2. We make no priority claim for
the degrees or the explicit factorizations here.

## Exact obstruction to removing harmonic localization in G2

Kostant--Michor Theorem 3.2 says that their separation identity without
localization holds precisely when theta maps regular orbits to regular
orbits. Both G2 trace maps here fail that condition.

For the minimal representation take (t_1,t_2)=(-1,2). None of the twelve
root characters equals 1, so this is regular semisimple. Its image in
simple-coroot coordinates is (-1/2,-3/4), annihilated by the root
alpha_2=(-3,2), and is consequently nonregular.

For the adjoint representation take (t_1,t_2)=(-1,eta), where
3 eta^2+2 eta+3=0. Again every root character differs from 1: the only
excluded eta values along t_1=-1 are 1,-1,i,-i, none of which satisfies
this quadratic. Along this line the alpha_1=(2,-1) root quotient is

\[
 \frac{\alpha_1(\theta_{14})}{e^{\alpha_1}-1}
 =-\frac{(t_2+1)(3t_2^2+2t_2+3)}{24t_2}.
\]

It vanishes at eta, giving a nonregular Lie image of a regular semisimple
group element. The audit checks these witnesses over QQ and the exact
quadratic number field, including all twelve root characters. This is a
consequence of the published criterion and the explicit formulas, not a
claim about a different choice of harmonic complement or a different map.
Use exact Sage rational/number-field arguments when evaluating Laurent
polynomials at points; Python integer powers with negative exponents can
otherwise introduce floating-point values.

## A uniform quadratic-coordinate formula and E6/E7 checks

Normalize the invariant bilinear form so long roots have squared length two.
Let q(h)=(h,h)/2, I_lambda the trace index, and
c_lambda=(lambda,lambda+2 rho_W) the Casimir eigenvalue. If
V_lambda tensor V_lambda = sum N_nu V_nu, then

\[
 \boxed{\theta_\lambda^*q=
 \frac{1}{4I_\lambda^2}\sum_\nu N_\nu(c_\nu-2c_\lambda)\chi_\nu.}
\]

To prove it, take an orthonormal Lie basis X_i. The tensor-product Casimir
is C tensor 1 + 1 tensor C + 2 sum rho(X_i) tensor rho(X_i).
Trace this identity against rho(t) tensor rho(t). The last term gives twice
the squared invariant norm of the character gradient. Root directions
contribute zero to that gradient on T, since tr(rho(t)rho(X_alpha))=0.
The trace-projection equation identifies the gradient with I_lambda theta.
Dividing by 4 I_lambda^2 proves the formula. This is the standard Casimir
tensor identity applied to this map, not a new Casimir identity.

For simply laced types q=(1/2) w^t C w, with C the Cartan matrix. In
fundamental characters (Sage/Bourbaki numbering), the four results are:

| Type and representation | Pullback of q |
| --- | --- |
| E6, V(omega_1), dimension 27 | chi_1^2/54 - chi_3/36 - 5 chi_6/36 |
| E6, adjoint V(omega_2), dimension 78 | chi_2^2/576 - chi_1 chi_6/144 - chi_2/192 - chi_4/576 - 1/64 |
| E7, V(omega_7), dimension 56 | chi_7^2/192 - chi_1/24 - chi_6/144 - 5/48 |
| E7, adjoint V(omega_1), dimension 133 | chi_1^2/1296 - 5 chi_1/648 - chi_3/1296 - 5 chi_6/1296 - 19/1296 |

Each is checked both by explicit Laurent substitution/reconstruction and
by the tensor-Casimir computation, using different formulas. These checks
share Sage's representation data; they are not independent representation
matrix constructions. See [the experiment](experiments/exceptional_quadratic_geometry.py)
and [the report](experiments/exceptional-quadratic-geometry.json).
Version 0.17 computed only the quadratic coordinate for E6/E7. The subsequent
[E6 extension](E6_GEOMETRY.md) now gives all six coordinates and both full
determinants as explicit products for the minimal and adjoint maps, with
independent matrix certificates. The [E7 extension](E7_GEOMETRY.md) supplies
all seven coordinates through degree 18 and both full determinant products,
with separate 56/133-dimensional matrix constructions and exact certificates.
The [generic-degree proof and API](EXCEPTIONAL_DEGREES.md) give 1080 and 2754
for E6, and 7272 and 33664 for E7, in minimal/adjoint order. The proof subtracts
the boundary Segre contributions from Newton volumes and restores the center
factor for adjoint representations. There is no E8 geometry result here.

## API and reproducibility

```python
from sage_borel import (springer_character_map, SpringerGeometry,
                       g2_adjoint_quotient, quadratic_character_pullback)
s = springer_character_map("G2", (1, 0))
geometry = SpringerGeometry(s)
quotient = g2_adjoint_quotient()
assert quotient.degrees == (2, 6)
assert quotient.generic_degree() == 6
assert g2_adjoint_quotient(representation="adjoint").generic_degree() == 18
assert geometry.cayley_determinant()(1, 1) == 1
x, y = geometry.character_ring.gens()
A = 2*x*x+13*x+y+13
assert geometry.character_coordinates(geometry.cayley_determinant()) == (
    (y-2)*(x+1)**2*A*A/35831808)
quadratic = quadratic_character_pullback("E6", (1,0,0,0,0,0))
assert quadratic  # rational Weyl-character-ring element in irreducible characters
```

`SpringerGeometry(s)` accepts the existing map, including direct sums.
`torus_jacobian()` returns the logarithmic r by r matrix.
`root_factors()` returns positive roots and their exact Laurent quotients.
`cayley_determinant()` returns the full determinant restricted to T.
`character_coordinates(f)` converts any Weyl-invariant Laurent input to a
polynomial in fundamental characters, constructing only the characters it
actually needs. `adjoint_quotient(invariants)` checks a supplied basic list
and returns both torus and character-coordinate pullbacks.
`AdjointQuotient.character_jacobian()` differentiates that quotient map.
`generic_degree()` uses the bound geometric proof for packaged E6/E7 maps;
otherwise it computes an exact generic fiber over a rational function field,
with a default rank guard of two. `algorithm='elimination'` forces that
calculation, including for the packaged maps. The determinant defaults to rank
at most four. Explicit guards also bound intermediate Laurent supports,
conversion steps, and the dimensions of needed fundamental characters.
They do not impose wall-clock or peak-memory limits on individual Sage calls.
`quadratic_character_pullback(type, highest_weight)` uses the uniform
tensor-Casimir formula and guards the input dimension before tensor decomposition.

Run inside Sage from borel/:

```sh
PYTHONPATH=src sage -python experiments/springer_geometry.py
PYTHONPATH=src sage -python experiments/audit_springer_geometry.py
PYTHONPATH=src sage -python experiments/exceptional_quadratic_geometry.py
PYTHONPATH=src sage -python -m unittest tests.test_geometry -v
```

The regression tests include the standard SL2 and SL3 maps, whose generic
degrees 2 and 3 are computed in Kostant--Michor Section 6. None of the
historical programs is changed.
