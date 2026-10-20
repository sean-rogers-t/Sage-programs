# Generic degrees of the E6 and E7 Springer maps

Date: 2026-10-20. We use simply connected complex groups, Bourbaki/Sage
numbering, and the normalized trace-projection Springer morphisms from
[the E6](E6_GEOMETRY.md) and [E7](E7_GEOMETRY.md) calculations. The results
below are deductions using toric intersection theory and the explicit maps.
They are not claims of originality or external peer review. Exact finite
checks of the proof's root, lattice, scheme and intersection inputs are in
[the geometry audit](experiments/exceptional-degree-geometry.json).

| group and representation | lattice for the toric calculation | normalized volume | boundary contribution | central cover | degree on the simply connected group |
| --- | --- | ---: | ---: | ---: | ---: |
| E6, V(omega_1), dimension 27 | P | 1080 | 0 | 1 | **1080** |
| E6, adjoint V(omega_2), dimension 78 | Q | 2808 | 1890 | 3 | **2754** |
| E7, V(omega_7), dimension 56 | P | 9792 | 2520 | 1 | **7272** |
| E7, adjoint V(omega_1), dimension 133 | Q | 21312 | 4480 | 2 | **33664** |

The last column is `(volume - boundary contribution) * central cover`.
In particular a Newton volume by itself gives the wrong answer in three
of these cases.

## API

```python
from sage_borel import (exceptional_springer_degree,
                       e6_adjoint_quotient, e7_adjoint_quotient)
assert e6_adjoint_quotient().generic_degree() == 1080
assert e6_adjoint_quotient(representation='adjoint').generic_degree() == 2754
assert e7_adjoint_quotient().generic_degree() == 7272
assert e7_adjoint_quotient(representation='adjoint').generic_degree() == 33664
c = exceptional_springer_degree('E6', representation='adjoint')
assert (c.normalized_volume, c.boundary_contribution, c.central_cover) == (2808,1890,3)
assert c.torus_degree == 918 and c.degree == 2754
```

These helpers use the proof below. `generic_degree(algorithm='elimination',
max_rank=7)` instead attempts an exact rational-function-field elimination;
that computation can be very expensive. A certificate is bound to the
packaged invariant degrees and character polynomials; modifying those
coordinates requires a fresh calculation. No optional fiber solver is
needed for the installed API.

## Degree and relation to the Springer/Borel maps

For a dominant generically finite Springer map theta_rho:G -> g, its degree
means the characteristic-zero field extension

\[
 \deg\theta_\rho=[\mathbb C(G):\theta_\rho^*\mathbb C(\mathfrak g)].
\]

[Kostant--Michor, Corollary 3.3](https://arxiv.org/pdf/math/0109066)
identifies this with the degree on adjoint quotients. The square of torus
and Weyl quotient maps also gives

\[
 \deg\theta_\rho
 =\deg(T\longrightarrow\mathfrak t)
 =\deg(T/W\longrightarrow\mathfrak t/W).
\]

Both horizontal Weyl quotient maps have degree |W|, which cancels. Thus one
does not divide the torus degree by |W|. This computes the degree of these
particular regular Springer maps, rather than the minimum over all rational
equivariant maps. The calculation concerns theta before the Borel quotient;
positive-degree G-invariants vanish after taking the Borel map.

Write P for the weight lattice, Q for the root lattice, and
T_M=Hom(M,C^*) for a character lattice M. For the adjoint representation,
the torus map factors through T_P -> T_Q. This is an isogeny of degree
[P:Q]=3 for E6 and 2 for E7. The center factors in the table restore the
degree for the simply connected group.

Let chi_rho be the representation character. In fundamental-weight
coordinates the explicit Springer formula is

\[
 \theta_\rho(t)=M_\rho^{-1}D\chi_\rho(t),\qquad
 D\chi_\rho=\sum_\mu m_\mu\mu e^\mu.
\]

The constant invertible matrix M_rho does not change a generic degree.
For the adjoint map the zero weights do not contribute, so its calculation
uses the sum over roots with coefficient one. The minimal characters also
have all weight multiplicities one.

The fundamental characters chi_i identify T_P/W with affine r-space:
this is [Steinberg, Theorem 6.1 and Lemma 6.3](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf).
In particular their tuple separates Weyl orbits of torus points. This
published input is used both for quotient fibers and for the critical orbit
calculation below. Normalization gives d theta at the identity equal to the
identity, so the maps in question are dominant and generically finite.

## The toric intersection calculation

Let Delta be the convex hull of the nonzero source weights, in the lattice
M indicated in the table. The origin is interior. On the normal projective
toric variety X_Delta let L be the ample toric line bundle of Delta and
H=c1(L). For a generic moment value b, the equations

\[
 D_i\chi_\rho-b_i=0\quad(1\leq i\leq r)
\]

are r sections of L. Their total intersection number is

\[
 \int_{X_\Delta}H^r=r!\operatorname{vol}_M(\Delta).
\]

This is the standard normalized-volume formula for toric intersection
numbers; see [Fulton--Sturmfels](https://arxiv.org/pdf/alg-geom/9403002).
An explicit determinant triangulation in the audit independently checks
the polyhedral volume routine.

On a proper exposed face F, the constant target terms disappear. If the
exposing primitive normal has positive height h, its linear combination
of the moment sections is h*chi_F. Thus common boundary zeros are exactly
the singular zeros chi_F=D chi_F=0 of the face character. All other face
strata have no common zero.
Multiplication of a face character by a nonzero monomial preserves this
singular-zero scheme: (m*f,D(m*f))=(f,Df). Thus the affine weight
configurations below may be translated before choosing their lattice covers.

If B is the boundary zero scheme for generic b, its contribution is

\[
 \int\{c(L^{\oplus r})\cap s(B,X_\Delta)\}_0.
\]

Subtract this contribution from the normalized volume to obtain the torus
degree. This is the localized top-Chern/Segre-class formula, not an
assumption that all boundary solutions are isolated. See the intersection
construction in [Aluffi--Marcolli, Section 2.5](https://www.its.caltech.edu/~matilde/AluffiMarcolliMathemAmplitudes.pdf).
All boundary zero schemes below lie in the smooth locus of X_Delta.
One can resolve the remaining toric singularities without changing the zero
scheme; the pulled-back line bundle has the same top intersection number.
In characteristic zero the generic torus fiber is reduced, since the
Springer map is dominant and generically etale.

## Facets, lattices and volumes

The complete facet descriptions used below are as follows. Coweights in
the table are primitive normals in the dual of the indicated lattice.

| polytope | facet-normal orbits | facets and their shapes | heights |
| --- | --- | --- | --- |
| E6 minimal in P | roots; W(3 omega_6^vee) | 72 five-simplices; 27 five-cross-polytopes | 1; 2 |
| E6 roots in Q | W(omega_1^vee); W(omega_6^vee) | 54 five-dimensional demicubes | 1 |
| E7 minimal in P | roots; W(2 omega_2^vee) | 126 six-cross-polytopes; 576 six-simplices | 1; 3 |
| E7 roots in Q | W(omega_2^vee); W(omega_7^vee) | 576 six-simplices; 56 E6 minuscule polytopes | 2; 1 |

The audit generates these normal orbits by integral reflections from
hardcoded Cartan matrices, verifies that their half-spaces have exactly
the desired vertices, and checks the dimension of every proposed facet.
The cone volumes over one facet in each orbit give

\[
\begin{aligned}
 V_{6,\min}&=72\cdot3+27\cdot32=1080,\\
 V_{6,\mathrm{ad},Q}&=54\cdot52=2808,\\
 V_{7,\min}&=576\cdot3+126\cdot64=9792,\\
 V_{7,\mathrm{ad},Q}&=576\cdot2+56\cdot360=21312.
\end{aligned}
\]

Simplex faces have no singular zero: the vertices are linearly independent,
so their nonzero monomials cannot have zero vector moment.

For a d-cross-polytope, let L_F be the lattice generated by vertex
differences and M_F its saturation in the facet lattice. The standard
cross character on a cover is sum_i(x_i+x_i^{-1}). Its singular zeros have
x_i=+/-1 and equally many signs of each kind. The cover from Z^d to the
root lattice D_d has degree two, identifying simultaneous sign reversal.
Consequently the number on T_{M_F} is

\[
 [M_F:L_F]\binom d{d/2}/2
\]

for even d, and zero for odd d. The Hessian is nonsingular: on the cover
its logarithmic Hessian is diagonal with entries +/-2. Isogenies are etale
over C. Smith forms in the audit determine all lattice indices used here.

At a height-one facet with a singular zero having nonsingular Hessian,
take a normal coordinate s and tangent coordinates u. The sections have
forms f_F(u)-b_n*s+O(s^2) and Df_F(u)+O(s). The tangent equations solve
u=O(s); the first then has leading term -b_n*s. For generic b_n this is
a reduced boundary point of intersection multiplicity one.

## E6 minimal and E7 minimal

Every proper face of the E6 minuscule polytope is a simplex or lies in
a five-cross-polytope. Odd cross dimension excludes singular zeros, so
there is no boundary contribution. Therefore the E6 minimal degree is 1080.

For E7, only the interiors of the 126 six-cross-polytope facets contribute.
Their difference lattice has index two in its saturation. Each has
2*binomial(6,3)/2=20 ordinary singular zeros. All facet heights are one,
so each contributes one. Lower faces are simplices and have no zero.
Thus the E7 minimal degree is 9792-126*20=7272.

## E7 adjoint: the E6 critical orbit

The only possible boundary zeros of the E7 root polytope lie inside its
56 E6 minuscule facets. The other facets are simplices. The facet at
omega_7^vee has roots with seventh simple-root coefficient one. Restriction
to the first six Dynkin nodes identifies its 27 weights with V_E6(omega_1).
The facet difference lattice is the full E6 root lattice.

We need the singular zeros of the E6 character chi_1, not the entire
fiber of the E6 Springer map. These satisfy theta_27=0 and chi_1=0.
Using the already certified six E6 quotient coordinates, form the ideal

\[
 I=(F_{27,2},F_{27,5},F_{27,6},F_{27,8},F_{27,9},F_{27,12},\chi_1).
\]

Exact rational reduction proves that its radical is

\[
 (\chi_1,\chi_2+3,\chi_3,\chi_4-9,\chi_5,\chi_6).
\]

For reproducibility, the audit verifies inclusion in this maximal ideal
and nilpotence powers (1,7,6,7,5,6) for its six generators modulo I.
Thus this is an exact support certificate, not an inference from a few
finite-field points. The target invariants are basic, so their simultaneous
zero on t forces theta_27=0.

Let zeta be a primitive cube root. The torus point with fundamental-weight
coordinates t_i=e^{omega_i}(t)=zeta^{v_i}, v=(0,0,0,0,1,2), has representation-character
tuple (0,-3,0,9,0,0). Its Weyl orbit has 240 points, verified by integral
reflections modulo three. Independently, its 18 centralizing roots split
into three orthogonal A2 systems, so its Weyl stabilizer has order 6^3=216
and its orbit has size 51840/216=240. Here the connectedness of semisimple
centralizers in simply connected groups is the standard input recorded in
[Steinberg, Remark 2.10](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf).
Fundamental characters separate Weyl orbits;
the unique character tuple therefore proves that these are **all** singular
zeros. The Cartan determinant from E6_GEOMETRY.md is 1/64 at this tuple,
so all 240 points have nonsingular Hessian. The center acts freely on the
torus and preserves the set, leaving 240/3=80 points in the root torus.
Every proper face of the E6 facet is already nondegenerate by the preceding
E6 minimal argument.

Each E7 adjoint facet has height one in Q and contributes 80 reduced points.
The root-torus degree is 21312-56*80=16832. Restoring the central double
cover gives degree **33664** on simply connected E7.

## E6 adjoint: the curve contribution

Here boundary zeros are not isolated. Each of the 54 facets is a D5
half-spin configuration with 16 weights. On its double cover the face
character, up to a monomial and an invertible coordinate change, is

\[
 f(x)=\frac{\prod_{i=1}^5(1+x_i)+\prod_{i=1}^5(1-x_i)}2.
\]

The cover involution is x -> -x. The audit explicitly identifies an actual
E6 facet with this half-spin configuration and checks its difference lattice.

The exact ideal (f,partial f/partial x_i) is the intersection of the 30
line ideals obtained by fixing two coordinates to +1 and two others to -1.
For a direct support proof, away from x_i=+/-1 the two products are nonzero;
f=0 makes them opposite, and no derivative can then vanish. With zeros in
both products, all derivatives vanish precisely when each product has at
least two zero factors. Equality of the ideals, including their reduced
scheme structures, is separately checked over Q.

The involution pairs these 30 curves, giving 15 per facet. The 20 points
with three signs of one kind and two of the other are also paired, giving
10 triple points per facet. Each belongs to three curves. The transverse
Hessian along a curve has rank four away from its two parameters +/-1;
at these parameters its rank drops to two. The half-spin polynomial near
a triple point is contact equivalent to

\[
 v_1v_2+u_1u_2u_3.
\]

Indeed its two product terms are a unit times v_1v_2 and a unit times
u_1u_2u_3; dividing by the second unit and rescaling v_2 gives the stated
form. With a generic nonzero normal target, elimination of the normal
coordinate and the two Morse variables gives the boundary zero ideal

\[
 (z_1,z_2,z_3,u_1u_2,u_1u_3,u_2u_3).
\]

Adding a multiple of f to its derivatives does not change this local ideal:
the weighted Euler identity expresses f in its Jacobian ideal with
coefficients vanishing at the point. This checks the full section scheme,
not only the reduced critical set of a face character.

The polytope has 270 codimension-two four-cross-polytope faces, each in
exactly two facets. Their difference lattice is saturated, so each has
binomial(4,2)/2=3 critical points. At each point one curve from each facet
meets. The normal cone is unimodular and the four tangent Hessian directions
are nonsingular. In normal coordinates s,t, solving those tangent equations
leaves normal equations with leading terms -b_1*s*t and -b_2*s*t.
The exact root-section Schur-complement check is recorded in the audit.
The two axes are the two incident facet curves, so both normal equations
are divisible by s*t to all orders. For generic b their remaining factor
is a unit. The boundary zero ideal is therefore

\[
 (z_1,z_2,z_3,z_4,st).
\]

All lower faces are simplices, hence have no zero. This proves that the
complete boundary scheme B is reduced, with

\[
 N_C=54\cdot15=810\text{ rational curves},\qquad
 N_3=54\cdot10=540\text{ triple points},\qquad
 N_2=270\cdot3=810\text{ double points}.
\]

Each curve closes to a P1: its free coordinate has two toric endpoints.
The root monomials along it have exponents zero or one in that coordinate,
so H.C=1. The toric variety is smooth along B; all its facet heights in Q
are one, so L=-K_X. Thus deg c1(N_C/X)=H.C-2=-1.

For clarity, the Segre corrections at the crossing points must also be
included. Two intersecting coordinate lines in P2 have Segre class
2H-4H^2, compared with 2H-2H^2 for the sum of their individual classes:
the point correction is -2. Three concurrent noncoplanar coordinate lines
in P3 have Segre class 3H^2-10H^3, compared with 3H^2-6H^3: the correction
is -4. The latter follows from the projective degrees (1,2,1,0) of
(xy,xz,yz) and

\[
 s(B,\mathbb P^n)=1-\sum_i\frac{g_iH^i}{(1+2H)^{i+1}}.
\]

This is [Aluffi, Proposition 3.1](https://www.math.fsu.edu/~aluffi/archive/paper173.pdf).
Normal cones are local; adding the common regular normal directions does
not change these point corrections relative to the separate branches.
More intrinsically, c(TX) cap s(B,X) is Fulton's intrinsic class of B
([Aluffi, Section 3.9](https://www.math.fsu.edu/~aluffi/archive/paper173.pdf)).
Its one-dimensional part is the sum of the branch classes. Therefore the
difference from the separate-branch Fulton classes is supported at the
singular points, and multiplication by the invertible class c(TX) does not
change its zero-dimensional coefficient. The two local models give exactly
the corrections computed above, independent of their ambient regular normal
directions. This justifies computing on the normalized components and gluing.
Since B has dimension one, only these point terms and the first normal Chern
classes enter the zero-dimensional Segre calculation.

A separate smooth curve contributes 6H.C-deg c1(N_C/X)=7. The complete
boundary contribution is consequently

\[
 \boxed{7N_C-4N_3-2N_2=7\cdot810-4\cdot540-2\cdot810=1890.}
\]

The root-torus degree is 2808-1890=918. The central triple cover gives
degree **2754** on simply connected E6.

## Independent checks and reproduction

[exceptional_degree_geometry.py](experiments/exceptional_degree_geometry.py)
checks all facet half-spaces and vertices, determinant volumes, Smith
indices, smooth normal cones, the actual D5 and E6 facet characters,
the half-spin Jacobian ideal, curve/triple incidence, local Hessians and
Segre coefficients, and the exact E6 critical-orbit support certificate.
Its derivation of the degree table uses these geometric inputs, not any
specialized-fiber length.

[exceptional_degree_fibers.py](experiments/exceptional_degree_fibers.py)
provides a separate exact finite-field fiber calculation from the full
quotient polynomials. Singular and msolve use different Groebner-basis
implementations. The reports label these as specialized fibers: agreement
at primes is corroboration, rather than a proof that a characteristic-zero
generic degree equals a specialized length. Optional elimination of one
character uses a constant nonzero coefficient in the quadratic coordinate,
so it is a global isomorphism of the fiber schemes.
The optional [msolve executable](https://msolve.proj.lip6.fr/binaries/index.html)
is version 0.11.0, used with its exact linear-algebra option `-l 2`;
it is not a package dependency or included binary.

```sh
PYTHONPATH=src sage -python experiments/exceptional_degree_geometry.py
PYTHONPATH=src sage -python experiments/exceptional_degree_fibers.py E6 minimal 32003 20261020 weighted
PYTHONPATH=src sage -python experiments/exceptional_degree_fibers.py E6 adjoint 32003 20261020 drl
# Optional msolve executable, obtained separately from its official project:
MSOLVE_BINARY=/path/to/msolve PYTHONPATH=src sage -python experiments/exceptional_degree_fibers.py E7 minimal 65521 20261021 msolve eliminate
```

The published inputs are the degree comparison, toric intersection formula,
fundamental-character quotient theorem, localized top-Chern formula and
Segre-class machinery. The exceptional
numbers and their boundary calculations are deductions recorded here.
This note does not settle priority or the minimum Cayley degrees of E6/E7.
