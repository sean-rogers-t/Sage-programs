# General degree criteria and an independent E6 audit

2026-10-22. The goal is to separate the general toric degree method from
the special product-character argument for half-spin representations.
The E6 minimal degree 1080 was already established in
[the exceptional degree note](EXCEPTIONAL_DEGREES.md). Reobtaining it here
is an independent structural check, not a new exceptional degree result.
No minimum Cayley degree, novelty, external review or new API is claimed.

## 1. Precise setting and relation to the Springer morphism

Work over C. Let M be the actual rank-r character lattice of a torus T,
N=Hom(M,Z), and let

\[
 \chi=\sum_{a\in A}c_a e^a,\qquad
 \Gamma_\chi=D\chi=\sum_{a\in A}c_a a e^a\in M\otimes\mathbb C.
\]

All displayed support coefficients are nonzero. Assume that the Newton
polytope Delta=conv(A) is full dimensional, that 0 is in its interior,
and that Gamma_chi is dominant. The latter holds if its Jacobian is
invertible at one point. The constant coefficient of chi is irrelevant
to Gamma, though the constant target monomial remains in the equations.
Normalized volume means r! times volume in which a fundamental
parallelepiped of M has volume one.

For a representation rho of a connected complex semisimple group G,
take M=X*(T), and c_a equal to weight multiplicities. Assume rho is
infinitesimally faithful. Its trace form is nondegenerate and

\[
 \theta_\rho|_T=M_\rho^{-1}\Gamma_\chi,
 \qquad M_\rho=\sum_a c_a aa^{\mathsf T}.
\]

This is an invertible linear change of target, so it preserves degree.
Here dominance follows from d(theta)_1=Id. The same degree applies to
G -> g and T/W -> t/W. For a regular semisimple target, equivariance
forces a preimage into its torus centralizer; the two Weyl quotient
degrees cancel. This agrees with
[Kostant--Michor, Corollary 3.3](https://arxiv.org/pdf/math/0109066).
It is the degree of this specified map, not the minimum over all
dominant equivariant rational maps. No parabolic or cohomology ring
enters this calculation.

Changing M matters. If the weights generate M0 of finite index in M,
the map factors through T_M -> T_M0 of degree [M:M0]. One must either
calculate volume in M throughout or restore that cover degree at the
end, never both.

## 2. A general theorem with a checkable hypothesis

**Theorem A (no boundary zeros).** In the setting above, suppose every
proper face F of Delta has no solution of

\[
 \chi_F=0,\qquad d\chi_F=0
\]

on its face torus, where chi_F is the face character after any monomial
translation into the face difference lattice. Then

\[
 \deg\Gamma_\chi=\operatorname{Vol}_M(\Delta).
\]

**Proof.** Let X_Delta be the normal projective toric variety and L its
polytope line bundle. The r equations D_i chi-b_i are sections of L.
On the orbit of a proper face the target term vanishes, since 0 is
interior. Pairing the remaining vector equation with an exposing
normal gives h chi_F=0, with h nonzero. The tangent equations then
give d chi_F=0. Conversely these equations annihilate the initial
vector equation. Thus the hypothesis excludes every boundary zero.
For general b, dominance and characteristic zero give a finite reduced
torus fiber. Its length equals the top intersection c1(L)^r, which is
the normalized volume. This intersection is supported in the smooth
open torus, so singularities of X_Delta cause no extra local terms. QED.

This is a direct application of standard toric intersection theory,
not a claimed new general theorem. The volume input is
Cox--Little--Schenck, *Toric Varieties*, Theorem 13.4.3
([authors' book page](https://dacox.people.amherst.edu/toric.html)).
We assert a sufficient criterion here; no converse is needed.

The same hypothesis gives global finite flatness of the torus gradient,
as proved in the [F4 follow-up](F4_MINIMAL_DEGREE.md#5-what-generalizes).
This also applies to the E6 minimal example. Finiteness of the full group
map is a separate question; the F4 note gives an explicit counterexample.

The test uses the face polynomial with its actual coefficients, not
just the shape of the polytope. It is unchanged by a monomial unit or
by passing to a finite torus cover in characteristic zero.

## 3. Isolated boundary zeros: the normal-layer correction

**Theorem B (isolated Morse face zeros).** Retain the setting of
Theorem A. Suppose all proper-face singular zeros lie in facet orbits,
are finite in number there, and have nonsingular tangent Hessians.
For each such point p on a facet of primitive height h, use a local
normal character s and tangent torus coordinates u. Write

\[
 \chi=s^{-h}F(s,u),\qquad
 F=f_0(u)+s f_1(u)+s^2f_2(u)+\cdots.
\]

Thus f_0(p)=df_0(p)=0. In coordinates dual to (s,u), the cleared
fiber equations are

\[
 D_uF-b_u s^h=0,\qquad
 sF_s-hF-b_s s^h=0.                         \tag{1}
\]

The tangent Hessian and the formal implicit function theorem give a
unique solution u=u(s;b_u) near p to the first equations. Define

\[
 m_p=\operatorname{ord}_{s=0}
 \left[(sF_s-hF)(s,u(s;b_u))-b_s s^h\right]
\]

at a general target. Then

\[
 1\le m_p\le h,\qquad
 \deg\Gamma_\chi=\operatorname{Vol}_M(\Delta)-\sum_p m_p. \tag{2}
\]

In particular h=1 implies m_p=1. More generally, if all intermediate
layers f_1,...,f_(h-1) vanish identically, then m_p=h.

**Proof.** A normal toric variety is smooth along every facet orbit:
locally it is a torus times A1. In the completed local ring the tangent
equations eliminate r-1 regular parameters. The remaining equation
has length its order in s. The tangent solution is independent of b_s,
so the coefficient of s^h in the remaining equation is a fixed
expression in b_u minus b_s. It is nonzero for general b_s. This proves
the bound and isolatedness. If intermediate layers vanish, u-p is
O(s^h), f_0(u) is O(s^(2h)), and the direct term s^h f_h cancels in
sF_s-hF. The first surviving generic coefficient is -b_s s^h.
There are no other boundary points by hypothesis. The isolated zero
scheme is in the smooth locus, and its total intersection length is
c1(L)^r; subtracting the boundary lengths proves (2). QED.

All orders refer to the actual lattice M. A cover used to simplify a
face character must not silently change point counts or normal orders.
The theorem is an elementary local intersection deduction. It does not
handle positive-dimensional critical loci such as the half-spin facets.

## 4. Why both coefficients and intermediate layers matter

The following are actual positive-multiplicity characters of T=(C*)^2,
with M=Z^2. Their trace moment matrices are invertible. Thus their
logarithmic gradients differ from their generalized Cayley maps only
by an invertible constant matrix.

First take

\[
 \chi_c=(x+x^{-1})(y+y^{-1}+c),\qquad c=1,2,3.
\]

All three have the same support and square polytope of normalized volume
8. For c=1,3 all faces have nonsingular zero sets and the degree is 8.
For c=2 the two x-facets have the double zero y=-1, each of height one,
so their Morse corrections are one each and the degree is 6.
Consequently the polytope, even with its complete unweighted support,
does not determine the degree.

For the normal-layer distinction, set

\[
 \psi_a=(x^2+x^{-2})(y+y^{-1}+2)+a(x+x^{-1}),\qquad a=0,1.
\]

Both have the same rectangle of normalized volume 16 and exactly the
same facet characters. Each x-facet has height two and the Morse zero
y=-1. With s=x^-1 at one facet,

\[
 F=(1+s^4)(y+y^{-1}+2)+a(s+s^3).
\]

For a=0 the missing intermediate layer gives m_p=2. For a=1 the normal
equation has leading term -s, giving m_p=1. The other x-facet is
identical by inversion, and all other faces are nondegenerate. Hence

| Character | Normalized volume | Total boundary correction | Generic degree |
| --- | ---: | ---: | ---: |
| chi_1 | 8 | 0 | 8 |
| chi_2 | 8 | 2 | 6 |
| chi_3 | 8 | 0 | 8 |
| psi_0 | 16 | 4 | 12 |
| psi_1 | 16 | 2 | 14 |

These are counterexamples to overly coarse degree rules, not new
simple-group classification results. Independently of the toric
argument, the verifier computes all five generic fibers exactly over
QQ(b1,b2), clearing monomials and imposing zxy=1. Thus the degrees
above are checked in characteristic zero with algebraically independent
targets, rather than inferred from finite-field samples.

## 5. E6 minimal: transfer succeeds without a boundary correction

For simply connected E6, rho=V(omega_1), use M=P(E6), with fundamental
weights as the integral basis and all 27 weights of multiplicity one.
The fresh verifier constructs the orbit using Sage's root-system
implementation, not the earlier custom reflection code, and enumerates
the actual facets rather than starting with the earlier facet table.

| Facet type | Number | Primitive height | Normalized cone volume |
| --- | ---: | ---: | ---: |
| Five-simplex | 72 | 1 | 3 |
| Five-crosspolytope | 27 | 2 | 32 |

The simplex character has no singular zero, because its affine
exponents are independent. On a finite torus cover, a crosspolytope
character, up to a monomial, is sum_(i=1)^5(z_i+z_i^-1). At a critical
point every z_i is +/-1. Its value is twice a sum of five signs,
which cannot vanish. Every proper subface of either facet type is a
simplex. This exhausts all proper faces, not only facets.

Theorem A therefore gives

\[
 \deg\theta_{V(\omega_1)}=72\cdot3+27\cdot32=1080.
\]

Every facet cone is also triangulated independently by determinants;
the total agrees with Sage's polytope volume. The weights generate
P(E6), so there is no missing central cover. These checks do not import
the packaged degree constant or the earlier quotient equations.

## 6. Audit of what the half-spin argument actually generalizes

The [half-spin proof](HALF_SPIN_RANK_DEGREES.md) uses the same general
compactification and conservation of number, but its additional
local claims require special hypotheses. Reinspection found no
contradiction in that written argument. This is an internal audit,
not an independent human review or a formal proof certificate.

| Step in that proof | Required input | Transfer status |
| --- | --- | --- |
| Compactification and total intersection | Actual character lattice, interior origin, dominance | General; Theorem A and the boundary intersection setup |
| Face critical strata | Face character equal, up to units, to two products on disjoint variable sets | Special to the stated product-character geometry; not automatic for a representation |
| Transverse multiplicity (k-1)(l-1) | After elimination the perturbed product-gradient equations have eta times units on their right sides | First-jet freedom is proved there from parity supports; not implied by polytope volume |
| f! free-coordinate solutions | Restricted perturbation has independent full-cube coefficients and the target ratios are general | Checked there; other strata need their own critical-point counts |
| No deeper-face escapes | Full section, including adjacent layers, lies in the two-product ideal; two fixed-normal target equations are independent | Essential special hypothesis; it cannot be replaced by counting facet zeros |
| Cover and branch counts | Actual Spin lattice and local finite covers tracked separately | Must be redone for a different group or central quotient |

The key valuation implication is valid because (A,B) is contained in
m_trans J_trans(F); hence the terms of F and its parameter derivatives
have strictly greater order than the transverse gradient. This is an
identity for the full local product expression, not merely its initial
face. It supports the deeper-face argument precisely where that full
expression is available. In using contact changes of coordinates, the
free equations must still be read as restrictions of the original
logarithmic system: multiplying a function by a nonconstant unit does
not by itself preserve its nonzero-target gradient equations.

The E6 minimal test uses only the general part. It neither tests nor
proves a universal formula for positive-dimensional boundary corrections.
The E6 adjoint boundary is a more demanding existing benchmark: its
five-dimensional half-spin facets meet along deeper faces, and the
existing proof obtains correction 1890 using Segre classes. A new
deformation derivation would have to reproduce that number with the
actual Q(E6) face lattices and all normal layers. One cannot simply
copy the Spin facet contribution 70: lattice quotients and normal
equations must be checked first.

## 7. Verification, status and next step

```bash
PYTHONPATH=src sage -python tools/verify_general_cayley_degree.py
```

The [script](tools/verify_general_cayley_degree.py) and
[report](general-cayley-degree-verification.json) distinguish independent
finite E6 geometry from generic characteristic-zero fiber calculations.

Published ingredients are the generalized Cayley construction and degree
comparison, toric intersection theory, and standard implicit-function
and intersection-multiplicity arguments. Theorems A/B are deductions
from these ingredients, presented with explicit usable hypotheses;
no priority claim is made. The E6 number is a recheck of an existing
repository deduction. The two elementary character families establish
the stated limitations by proof and exact generic elimination. There
has been no exhaustive literature search for these formulations.

**Follow-up completed:** the [E6 adjoint deformation proof](E6_ADJOINT_DEFORMATION.md)
now obtains 1890 from open-curve and triple strata and excludes all crossing
limits. It verifies middle-layer stationarity and the actual facet cover.
The subsequent [F4 minimal test](F4_MINIMAL_DEGREE.md) now gives degree 96
and the finite-flat torus corollary. The original transfer goal
below records the question this follow-up resolves.

**Original transfer goal:** determine whether the two-product deformation
argument reproduces the existing E6 adjoint correction 1890 without
using its Segre-class computation. This is a transfer test for a
general method, not a search for another Spin rank. Success means
accounting for every facet and deeper-face contribution in Q(E6),
obtaining root-torus degree 918 and restoring the central factor three
to give 2754. Failure should identify the missing local hypothesis.
Only after this test should the method be promoted as a general
positive-dimensional boundary algorithm.
