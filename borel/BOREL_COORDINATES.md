# Springer maps on Borel subgroups: SL3 and Sp4

## Scope and normalization

We work over C (all displayed identities are defined over Q), with the
**defining representations** of SL3 and Sp4, their diagonal maximal torus
T, and upper triangular Borel subgroup B=T U. We write g=tu, in that order.
The coordinates below are epsilon characters of T, not the
fundamental-weight coordinates of the exceptional-character API.

The Springer map is trace projection: theta(g) is the unique element of
the Lie algebra satisfying

\[
 \operatorname{tr}(\theta(g)X)=\operatorname{tr}(gX)
 \quad\text{for every }X\in\mathfrak g.
\]

Thus theta(g)=g-tr(g)I/3 for SL3 and theta(g)=(g-g^{-1})/2 for Sp4.
It has d theta at 1 equal to the identity. These are the standard
generalized Cayley maps, with no additional scalar rescaling.

This note determines theta on B and the **fixed-t** polynomial map
U -> mathfrak u obtained by subtracting theta(t). Its derivative is called
the vertical derivative. This is different from the full derivative of
B -> mathfrak b, or G -> mathfrak g. No full target branch or non-finite
locus is asserted here. The cohomological Borel homomorphism in `BorelMap`
is a different map from this restriction to a Borel subgroup.

## A triangular mechanism that generalizes

Let G be connected reductive in characteristic zero, with a representation
whose trace form is nondegenerate on mathfrak g, and let theta be its trace
projection map. Fix positive roots and root vectors E_alpha, and order
the positive roots by increasing height. Use root coordinates

\[
 u=\prod_{\alpha>0}x_\alpha(u_\alpha),\qquad
 \theta(tu)=\theta(t)+\sum_{\alpha>0}F_\alpha(t,u)E_\alpha.
\]

Then

\[
 F_\alpha(t,u)=q_\alpha(t)u_\alpha+
 P_\alpha(t,(u_\beta)_{\operatorname{ht}\beta<\operatorname{ht}\alpha}),
 \qquad
 q_\alpha(t)=\frac{e^\alpha(t)\alpha(\theta(t))}{e^\alpha(t)-1}.
\]

The quotient denotes its regular Laurent-polynomial extension across
e^alpha=1, not evaluation of an uncancelled rational expression.
In particular q_alpha(1)=1. Consequently

\[
 \det D_U\theta(tu)=\prod_{\alpha>0}q_\alpha(t).
\]

When every q_alpha(t) is nonzero, triangular back substitution makes the
fixed-t map a polynomial isomorphism. If a factor vanishes, the vertical
Jacobian is singular everywhere on that torus fiber. Counting vanishing
factors alone does **not** give the rank: lower-root mixed monomials matter.

Proof. A cocharacter contracting U to 1 shows theta(B) is contained in
mathfrak b: equivariance and the existence of the contraction limit forbid
negative-root components. T-equivariance makes F_alpha a polynomial of
T-weight alpha. Its only term involving u_alpha is q_alpha u_alpha;
every other contributing positive root has smaller height. Differentiating
theta(exp(s E_alpha)t exp(-s E_alpha)) at s=0 gives
q_alpha(e^{-alpha}(t)-1)=-alpha(theta(t)), proving the displayed formula
on a dense open set. The coefficient q_alpha was regular to begin with,
so this identity extends everywhere. Finally d theta at 1 is identity
because trace projection restricts to identity on mathfrak g.

This elementary consequence of equivariance explains what the two
implementations test. It is not claimed as a new literature theorem or
as a complete arbitrary-representation algorithm.

## SL3: no unipotent collapse

Put t=diag(z1,z2,(z1 z2)^{-1}) and order the root vectors as
E12,E23,E13, corresponding to alpha1, alpha2, alpha1+alpha2. Then

\[
 u=(I+aE_{12})(I+bE_{23})(I+cE_{13})
 =\begin{pmatrix}1&a&ab+c\\0&1&b\\0&0&1\end{pmatrix}.
\]

If m=(z1+z2+(z1 z2)^{-1})/3, the toral diagonal is
(z1-m,z2-m,(z1 z2)^{-1}-m). The positive Lie-root coefficients are

\[
 (p,q,r)=(z_1a,\ z_2b,\ z_1(c+ab)).
\]

The vertical determinant is z1^2 z2, a unit on T. For every t,

\[
 a=p/z_1,\qquad b=q/z_2,\qquad c=r/z_1-ab
\]

is a polynomial inverse. This includes torus points with repeated
eigenvalues. No unipotent direction collapses anywhere in B for this map.

## Sp4: formulas and complete fixed-t fibers

Use the vector basis (e1,e2,f2,f1), with symplectic form

\[
 J=\begin{pmatrix}0&0&0&1\\0&0&1&0\\0&-1&0&0\\-1&0&0&0\end{pmatrix},
 \qquad t=\operatorname{diag}(z_1,z_2,z_2^{-1},z_1^{-1}).
\]

For simple roots alpha1=e1-e2 and alpha2=2e2, the ordered root vectors are
E12-E34, E23, E13+E24, E14. They correspond to alpha1, alpha2,
alpha1+alpha2, 2 alpha1+alpha2. Each squares to zero, so x_alpha(v)=I+vE_alpha.
Their product is

\[
 u=\begin{pmatrix}
 1&a&ab+c&d+ac\\0&1&b&c\\0&0&1&-a\\0&0&0&1
 \end{pmatrix}.
\]

Set A=(z1+z2^{-1})/2 and C=(z1+z2)/2. The four Lie-root coefficients are

\[
 \boxed{p=Aa,\qquad q=z_2b,\qquad
 r=Cc+\tfrac12z_1ab,\qquad s=z_1(d+ac).}
\]

The full matrix is

\[
 \theta(tu)=\begin{pmatrix}
 h_1&p&r&s\\0&h_2&q&r\\0&0&-h_2&-p\\0&0&0&-h_1
 \end{pmatrix},\qquad h_i=(z_i-z_i^{-1})/2.
\]

The vertical Jacobian and determinant are

\[
 D_U\theta=\begin{pmatrix}
 A&0&0&0\\0&z_2&0&0\\z_1b/2&z_1a/2&C&0\\z_1c&0&z_1a&z_1
 \end{pmatrix},\qquad
 \det D_U\theta=\frac{z_1}{4}(z_1z_2+1)(z_1+z_2).
\]

The b and d pivots are always nonzero. Eliminating them leaves the block
[[A,0],[z1 b/2,C]] for the remaining rank calculation. This proves the
following classification over C, at every point, not only generically.
Targets in the table mean positive-root coordinates with diagonal theta(t).

| Torus stratum | Source condition | Vertical rank | Exact image condition | Nonempty fiber |
| --- | --- | --- | --- | --- |
| A C != 0 | Any u | 4 | Every (p,q,r,s) | One point |
| A=0, C!=0 | Any u | 3 | p=0 | Affine line |
| A!=0, C=0 | Any u | 3 | 2 A z2 r = z1 p q | Affine line |
| A=C=0 | b!=0 | 3 | p=0, q!=0; r,s arbitrary | Affine line |
| A=C=0 | b=0 | 2 | p=q=r=0; s arbitrary | Affine plane |

The divisors are z1 z2=-1 and z1+z2=0. Their intersection is exactly
(z1,z2)=(1,-1),(-1,1). At either crossing, the exact image is
{p=0,q!=0} union {p=q=r=0}, with s unrestricted. It is **not closed**;
its closure is p=0. In particular a target with p=q=0 and r!=0 has no
preimage, despite belonging to the image closure.

Here is an explicit proof of every fiber assertion. Always set b=q/z2
and d=s/z1-ac. If A C!=0, set a=p/A and c=(r-z1 ab/2)/C.
If A=0,C!=0, impose p=0 and choose a freely, then solve c by that same
formula. If A!=0,C=0, a=p/A is fixed, the displayed image relation is
necessary and sufficient, and c is free. At A=C=0 and b!=0, impose p=0,
set a=2r/(z1 b), and choose c freely. At A=C=0 and b=0, impose p=r=0
and choose both a and c freely. These parametrizations and their inverse
coordinate projections are regular, so the fibers are the stated affine
spaces, with their reduced scheme structures.

The infinitesimal collapsed directions can also be read exactly. On
A=0,C!=0 they are spanned in (a,b,c,d) coordinates by
(1,0,-z1 b/(2C),-c+a z1 b/(2C)). On C=0,A!=0, and at the crossings with
b!=0, they are spanned by (0,0,1,-a). At the crossings with b=0 there
are two independent directions (1,0,0,-c) and (0,0,1,-a).
Thus the mixed term ab recovers a direction at the crossing when b!=0.
It would be incorrect to report rank two on the entire crossing fiber.

At t=1 all root factors are one: the restriction U -> mathfrak u is a
polynomial isomorphism in both groups, consistent with the established
unipotent Springer isomorphism. This does not force the other tU fibers
to be isomorphisms, as the Sp4 calculation demonstrates.

## API

[Implementation](src/sage_borel/borel_coordinates.py) exports
`BorelSpringerMap` and `borel_springer_map`. Supported group strings are
exactly `SL3` and `Sp4`, with defining representations. No arguments to
formula/matrix methods returns symbolic rational functions over Q.
Point evaluation, inverses and fibers currently accept rational data.
`fiber()` gives `None` for an empty fiber, otherwise its affine-space
dimension, free polynomial parameters and source coordinates.

```python
from sage_borel import borel_springer_map

sl = borel_springer_map("SL3")
assert sl.inverse((2, 3), sl.coordinates((2, 3), (1, 2, 3))) == (1, 2, 3)

sp = borel_springer_map("Sp4")
formulas = sp.coordinates()                 # positive Lie-root coefficients
source_matrix = sp.group_matrix()          # tu
springer_matrix = sp.lie_matrix()          # theta(tu)
vertical_jacobian = sp.vertical_jacobian()  # t held fixed
assert sp.vertical_rank((1, -1), (2, 3, 5, 7)) == 3
assert sp.vertical_rank((1, -1), (2, 0, 5, 7)) == 2
assert sp.fiber((1, -1), (0, 0, 0, 7))["dimension"] == 2
assert sp.fiber((1, -1), (0, 0, 1, 7)) is None
```

## Independent checks, published inputs, and next goal

[Tests](tests/test_borel_coordinates.py) form the ordered root products
independently of the API matrices. They compute trace projection by
solving the full trace Gram system: an 8-element basis of sl3 and a
10-element basis of sp4 constructed from the condition JX symmetric.
This does not use either closed formula to define the projection. The
result is compared both to the API and to g-tr(g)I/3 or (g-g^{-1})/2.
The tests are identities in all torus/root variables over Q, not sampled
floating-point comparisons. They also check symbolic stratum ranks,
the general root-factor identity, fiber parametrizations identically in
their free parameters, empty fibers, and removable factors at t=1.
[Verification report](borel-coordinates-verification.json) binds the
implementation, tests and this proof note. These are independent
constructions within Sage, not an external review or separate CAS audit.

The trace-projection construction, equivariance and unipotent isomorphism
are published inputs: Kostant--Michor,
[The generalized Cayley map](https://arxiv.org/abs/math/0109066), especially
the construction in Section 2 and Theorem 4.6. Rogers,
[An explicit determination of the Springer morphism](https://arxiv.org/abs/1701.01538),
provides the torus weight/trace-matrix formula. The triangular argument
and this coordinate/fiber classification are written deductions from
these inputs; no novelty claim or comprehensive priority search is made.

The next useful step is a root-data implementation of the triangular
coefficients for more representations. Its goal is to determine when
vanishing root factors cause additional rank loss and when mixed root
terms recover directions, without constructing large representation
matrices. A successful first milestone is reproducing these two complete
fiber classifications from torus formulas and root commutators, then
testing a non-defining representation. An unrestricted subgroup-fiber
classification, full B/G differential analysis and target branch locus
remain separate tasks.
