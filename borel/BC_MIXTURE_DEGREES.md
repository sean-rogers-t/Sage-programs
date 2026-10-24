# Arbitrary-rank defining-adjoint Springer degrees in types B and C

2026-10-22. These are characteristic-zero written deductions, corroborated
by independent finite calculations. The degree is for the specified
trace-projection morphism, not the minimum Cayley degree of a group.
No novelty or external-review claim is made.

## 1. The two theorems

Work over C. For type B let G=SO_(2n+1), with defining module V of
dimension 2n+1. For type C let G=Sp_(2n), with defining module V of
dimension 2n. In both cases n>=1, and Ad denotes the adjoint module
of dimension n(2n+1). For integers a,b>=0 with a+b>0, put

\[
 \rho_{a,b}=V^{\oplus a}\oplus\operatorname{Ad}^{\oplus b}.
\]

Let theta_(a,b):G -> g be its Springer/generalized Cayley morphism,
with the trace-projection normalization. Denote its generic degree by
d_B(n;a,b) or d_C(n;a,b). The same degree applies to its torus restriction
T -> t and induced adjoint-quotient map T/W -> t/W.

**Type B theorem.** If b=0, then d_B(n;a,0)=2^n. For b>0,

\[
 \boxed{
 d_B(n;a,b)=2^n(2^n-n)-
 2n\sum_{k=0}^{n-1}\binom{n-1}{k}
   \mathbf1_{\{a=(4k-2n+1)b\}}.}
 \tag{1}
\]

At most one summand is nonzero. For n>=2 the exceptional positive ratios
are exactly

\[
 \frac ab=
 \begin{cases}
 1,5,9,\ldots,2n-3,&n\text{ even},\\
 3,7,11,\ldots,2n-3,&n\text{ odd}.
 \end{cases}                                             \tag{2}
\]

At a ratio r in (2), k=(r+2n-1)/4 and the boundary correction is
2n binomial(n-1,k). There are floor(n/2) exceptional ratios.
The pure adjoint case a=0 is never exceptional. At n=1 the defining
and adjoint SO_3 modules coincide, and the formula gives degree 2
for every mixture, with no exceptional ratios.

**Type C theorem.** There are no exceptional ratios for b>0:

\[
 \boxed{
 d_C(n;a,b)=
 \begin{cases}2^n,&b=0,\\4^n,&b>0.\end{cases}}
 \tag{3}
\]

The boundary correction is zero in both cases of (3). The torus and
adjoint-quotient maps are finite flat of the indicated rank. In type B
the same finite-flat conclusion holds when b=0 or the ratio is not
exceptional. At every exceptional type B ratio, both maps are NOT
finite: Section 4 exhibits a nonconstant torus curve collapsed to zero.
No finiteness of the full group maps is asserted.

The statements cover all ranks and all nonnegative multiplicities by
the proofs below. Small-rank tests are corroboration, not the basis for
extrapolating (1)--(3). There is no parabolic or cohomology coefficient
ring in these theorems.

## 2. Lattices, central quotients and normalization

For SO_(2n+1) and Sp_(2n), use the actual character lattice M=Z^n,
with z_i=e^(e_i). Put

\[
 x_i=z_i+z_i^{-1},\qquad w_i=z_i-z_i^{-1},\qquad X=\sum_i x_i.
\]

These satisfy w_i^2=x_i^2-4. Conversely z_i=(x_i+w_i)/2, with inverse
(x_i-w_i)/2. Thus these equations present the SAME torus, without
a cover and without excluding the points w_i=0.

For type B, the defining character is 1+sum x_i and the adjoint
character is n+sum x_i+sum_(i<j)x_i x_j. For type C, the defining
character is sum x_i and the adjoint character is

\[
 n+\sum_i(z_i^2+z_i^{-2})+\sum_{i<j}x_i x_j
 =\sum_i x_i^2+\sum_{i<j}x_i x_j-n.
\]

Consequently, after dropping constants, the logarithmic gradients are

\[
 \Gamma_{B,i}=w_i\bigl[a+b+b(X-x_i)\bigr],\qquad
 \Gamma_{C,i}=w_i\bigl[a+b(X+x_i)\bigr].                \tag{4}
\]

The epsilon-coordinate trace moment matrices are respectively

\[
 M_B=\bigl(2a+(4n-2)b\bigr)I_n,\qquad
 M_C=\bigl(2a+4(n+1)b\bigr)I_n.                       \tag{5}
\]

For example, in type C the roots +/-2e_i contribute 8I_n and the
roots +/-e_i +/-e_j contribute 4(n-1)I_n to the adjoint trace form.
The normalized maps are Gamma_B or Gamma_C divided by the scalar
in (5). These scalars are positive for the stated multiplicities.
The identity differential is the identity, so the maps are dominant.

Equivalently, on the entire group,

\[
 \theta^B_{a,b}=
 \frac{a\theta_V+(2n-1)b\theta_{\rm Ad}}{a+(2n-1)b},
 \qquad
 \theta^C_{a,b}=
 \frac{a\theta_V+2(n+1)b\theta_{\rm Ad}}{a+2(n+1)b}.
 \tag{6}
\]

Thus these are mixtures of the Springer maps already under study,
weighted by their trace forms. Simultaneously multiplying a,b by a
positive integer, or adding trivial summands, leaves the map unchanged.

**Spin cover.** Pulling the type B representation back along
Spin_(2n+1) -> SO_(2n+1) doubles all degrees in (1), including b=0.
Its character lattice is Z^n+Z*(1/2,...,1/2), while every weight of
these defining/adjoint mixtures belongs to Z^n. The map factors through
the central cover of degree two. The degree for the Spin group, its
torus and its own adjoint quotient is therefore 2*d_B. This statement
concerns the pulled-back vector representation, not the spin module.

**Adjoint type C group.** If a=0, the representation descends to PSp_(2n).
Its degree there is 4^n/2=2^(2n-1), since Q(C_n) consists of the integer
vectors with even coordinate sum and has index two in Z^n. If a>0,
the defining summand prevents descent, regardless of the parity of a.
The original Sp degree (3) already includes its central cover factor.

For a useful rank-two control, Spin_5 is isomorphic to Sp_4, and both
pure adjoint formulas give degree 16; their adjoint-form groups SO_5 and
PSp_4 both give degree 8. Their defining modules are different and
should not be identified in mixtures.

## 3. Type B polytope and volume

Suppose n>=2 and b>0. Write c=a+b. The nonzero character is

\[
 c\sum_i(z_i+z_i^{-1})+
 b\sum_{i<j}(z_i+z_i^{-1})(z_j+z_j^{-1}).               \tag{7}
\]

Its Newton polytope is

\[
 \Delta_B=\operatorname{conv}\{\pm e_i\pm e_j:i<j\}
 =\{u:|u_i|\le1,\ \sum_i|u_i|\le2\}.                 \tag{8}
\]

The short weights +/-e_i are centers of its coordinate facets. For
n>=3 the facets have two types:

- The 2n coordinate facets u_i=+/-1, translates of (n-1)-dimensional
  crosspolytopes, with their short-weight center included in the support.
  Their primitive lattice height in Z^n is one.
- The 2^n sign facets sum epsilon_i u_i=2, copies of the hypersimplex
  Delta(2,n), with no additional weights. Their primitive height is two.

For n=2 the polytope is a square with only the four coordinate facets;
the would-be sign inequalities touch vertices and are redundant as
facet inequalities. This low-rank distinction does not change (1).

This root-polytope description is established combinatorics; see
Ardila--Beck--Hosten--Pfeifle--Seashore, Proposition 18 and Sections 4--5,
[Root polytopes and growth series of root lattices](https://arxiv.org/abs/0809.5123).
Their D_n root lattice has index two in the Z^n lattice used here.
The volume must be adjusted accordingly.

Here is a direct volume derivation in Z^n. In one orthant, (8) is the
simplex sum u_i<=2, u_i>=0, truncated by u_i<=1. The original simplex
has volume 2^n/n!. Each of the n removed parts has volume 1/n!;
intersections of two removed parts have measure zero. The 2^n orthants
therefore give normalized volume

\[
 \operatorname{Vol}_{\mathbb Z^n}(\Delta_B)
 =2^n(2^n-n).                                         \tag{9}
\]

Normalized volume means n! times ordinary lattice volume. Formula (9)
also gives 2 at n=1, though that case uses the interval of short weights.

## 4. Type B: all boundary zeros and their multiplicities

On a coordinate facet, translate by its short-weight center. In its
actual tangent lattice Z^(n-1), the face character is

\[
 f=c+b\sum_{j=1}^{n-1}(t_j+t_j^{-1}).                  \tag{10}
\]

The critical equations force t_j=+/-1. If exactly k signs are negative,
then f=0 precisely when

\[
 a+b+2b(n-1-2k)=0
 \quad\Longleftrightarrow\quad a=(4k-2n+1)b.           \tag{11}
\]

There are binomial(n-1,k) such points on each of the 2n coordinate
facets. The logarithmic Hessian is diagonal with entries 2b*t_j and
is invertible. Each point is an isolated Morse zero. Proper faces of
a crosspolytope are simplices, with no short-weight centers in them,
and their characters have no singular zeros.

It remains essential to check the sign facets and THEIR lower faces.
After changing signs, a sign-facet character is b*e_2(t_1,...,t_n),
where e_2=sum_(i<j)t_i t_j. If its zero were singular on the face torus,
Euler's homogeneous identity implies that all ambient logarithmic
derivatives vanish as well. Since each t_i is nonzero, this gives

\[
 T-t_i=0\quad\text{for every }i,\qquad T=\sum_i t_i.
\]

Thus t_i=T and (n-1)T=0, impossible on the torus over C.
This argument is legitimate on the face torus: one may set one
homogeneous coordinate to 1, or pull back to the homogeneous coordinate
torus and use Euler's identity at a zero.

A face of Delta(2,n) fixes some coordinates to 0 or 1 in its polytope
description. Its polynomial is therefore one of:

- e_2 in a smaller subset of at least two variables;
- a monomial times a nonempty sum of variables;
- a single nonzero monomial.

The same gradient argument excludes the first; the latter two cannot
have singular zeros. Hence (10)--(11) list ALL boundary zeros. In
particular there are no corrections on hypersimplex facets or their
intersections, even when those lower faces are not simplices.

Finally, every coordinate facet has primitive height one. In a primitive
normal coordinate s, write the character as s^-1 F(s,t). At a Morse
face zero p, the tangent fiber equations solve t=t(s) uniquely. The
remaining cleared normal equation is

\[
 (sF_s-F)(s,t(s))-\eta_s s=0.
\]

Its linear coefficient is nonzero at a general target, so its local
intersection length is one. There is no intermediate normal layer at
height one. This is the isolated-Morse criterion proved in the
[general degree note](GENERAL_CAYLEY_DEGREE.md#3-isolated-boundary-zeros-the-normal-layer-correction).

The total compactified intersection has length (9). Subtracting one for
each of the 2n binomial(n-1,k) boundary points proves (1) for b>0.
For b=0 the coordinates of the normalized map are
(z_i-z_i^-1)/2; each generic coordinate equation is quadratic and has
two roots. Their product gives 2^n. At n=1, defining and adjoint agree
and give the same quadratic map for every mixture. This proves all
rank and parameter cases of the type B theorem.

**An actual non-finite fiber at every exceptional ratio.** Fix k minus
signs epsilon_2,...,epsilon_n satisfying (11), and let

\[
 z_1=q\in\mathbb C^*,\qquad z_j=\epsilon_j\quad(j>1).
\]

For j>1, w_j=0. The remaining component is
Gamma_(B,1)=(q-q^-1)[a+b+2b sum_(j>1) epsilon_j]=0 by (11).
Thus the entire curve lies in the zero fiber. Its image in T/W is
nonconstant, since the invariant sum x_i=q+q^-1+2 sum epsilon_j
varies. The quotient map also has a positive-dimensional zero fiber.
This proves non-finiteness at the exceptional ratios and completes the
finite-versus-non-finite classification of the torus and quotient maps
within this family. It does not determine their maximal finite target
open subsets; that is the next question below.

## 5. Type C: why no ratio causes a correction

For b>0 the Newton polytope is the doubled crosspolytope

\[
 \Delta_C=\operatorname{conv}\{\pm2e_i\}
 =\{u:\sum_i|u_i|\le2\}.
\]

The defining weights +/-e_i lie strictly inside it, so varying a does
not affect ANY boundary face character. A sign facet, with variables
t_i=z_i^(epsilon_i), has character

\[
 b\,h_2(t)=b\left(\sum_i t_i^2+\sum_{i<j}t_i t_j\right).
 \tag{12}
\]

Here h_2 is the complete homogeneous quadratic. Every smaller face has
the same form in a nonempty subset of m variables. At a hypothetical
singular zero, homogeneity again makes all ambient logarithmic
derivatives vanish. They yield

\[
 t_i(T+t_i)=0,
 \quad T=\sum_i t_i.
\]

Every t_i is nonzero, so t_i=-T and (m+1)T=0. This contradicts the
torus condition in characteristic zero. Thus EVERY proper face has
nonsingular zero set, for every a>=0 and b>0.

The volume of Delta_C in Z^n is 2^n*2^n=4^n: each of its 2^n orthants
is a simplex of normalized volume 2^n. The no-boundary criterion gives
degree 4^n. For b=0 the same coordinatewise quadratic calculation as
type B gives 2^n. This proves (3), including n=1.

The no-boundary criterion also gives the stated finite flatness of
the torus maps: the graph is closed in its projective toric
compactification over the target, hence proper, and the map is affine,
hence finite. Miracle flatness applies because both tori and affine
targets are smooth. For the SO and Sp quotients, the invariant ring
is C[x_1,...,x_n]^(S_n), a polynomial ring in the elementary symmetric
functions of x_i. The finite quotient map is again flat. This uses
standard algebraic geometry, as detailed in
[the F4 note](F4_MINIMAL_DEGREE.md#5-what-generalizes).

## 6. Examples and the connection to F4

These are SO degrees in the B rows; multiply by two for Spin.

| Family | Degree with b=0 | Usual degree with b>0 | Exceptional a/b : degree |
|---|---:|---:|---|
| B2 | 4 | 8 | 1 : 4 |
| B3 | 8 | 40 | 3 : 34 |
| B4 | 16 | 192 | 1 : 168; 5 : 184 |
| B5 | 32 | 864 | 3 : 824; 7 : 854 |
| B6 | 64 | 3712 | 1 : 3592; 5 : 3652; 9 : 3700 |
| C_n | 2^n | 4^n | None |

For B4 the exceptional ratios 1 and 5 coincide with those of our F4
family, because both have centered three-dimensional crosspolytope
face characters (10). Their corrections differ: B4 on SO has eight
height-one coordinate facets, whereas F4 has twenty-four height-two
facets. On the Spin cover the eight B4 facets have height two and no
intermediate layer; the total correction doubles. F4 also has a
nonempty intermediate half-root layer, whose vanishing was separately
checked in the [F4 family proof](F4_MIXTURE_DEGREES.md).

This distinction explains why the facet critical values alone predict
the ratios but cannot determine the correction without lattice and
normal-layer information.

## 7. Independent calculations and precise scope

The [verifier](tools/verify_bc_mixture_degrees.py) checks explicit weight
lists against independent highest-weight character constructions for
both types in ranks 2 through 6. It verifies the trace forms, root
support lattices, actual polytopes, facet types and normalized volumes.
Direct sign enumeration through rank 10 agrees with the binomial
corrections. In those ranks, original weight sums and trace coordinates
also both verify the zero-curve identities at every exceptional ratio.
Exact gradient matrices J-I and J+I have determinants
(-1)^(m-1)(m-1) and m+1, respectively, confirming the e_2/h_2 exclusion
mechanism through 12 variables. These bounded calculations corroborate
the symbolic proofs for arbitrary rank.

Original fibers are constructed in two ways:

1. Enumerate root and defining weights in the actual Z^n lattice,
   differentiate their Laurent monomials, clear denominators and
   saturate by the product of the torus coordinates.
2. Use the independent equations (4), together with
   w_i^2-x_i^2+4=0. These polynomial equations already present the torus;
   no division by w_i or by a possibly vanishing mixture coefficient
   is performed.

The same normalized epsilon target is used in both constructions,
with coordinates (1,3,7,13) truncated to the rank. Multiplication by
sum i*z_i, expressed in the second model as sum i*(x_i+w_i)/2, is
compared using its full characteristic polynomial. Squarefreeness
certifies geometric reducedness of the sampled fibers.

The checks cover seventeen rank/multiplicity cases over two primes
(32003 and 32009), each in both constructions: 68 fibers. They include
all B2--B4 exceptional ratios, pure adjoints, nonexceptional mixtures,
minimal-only controls and rank-one cases. Two additional checks verify
degree 68 on the Spin_7 cover of the B3 ratio-three case and degree 8
on PSp_4 for the adjoint representation.

See the [exact report](bc-mixture-degrees-verification.json). The
computations share Sage/Singular arithmetic; independence is in the
equation constructions, character construction and geometric argument.
Finite-field lengths are not used to infer characteristic-zero generic
degree, and no arbitrary-rank theorem is inferred from finite ranks.

```bash
PYTHONPATH=src sage -python tools/verify_bc_mixture_degrees.py
```

## 8. Published inputs, deductions and next goal

The representation-theoretic trace construction and the equality of
group, torus and quotient generic degrees are published inputs:
[Kostant--Michor, Corollary 3.3](https://arxiv.org/pdf/math/0109066).
The root-polytope geometry and volumes are established combinatorics,
with an elementary volume derivation included above; see the cited
[root-polytope paper](https://arxiv.org/abs/0809.5123).
The top-intersection calculation and isolated Morse local lengths use
standard toric intersection theory and the implicit function theorem.

Our deductions here are the complete mixture formulas (1)--(3), the
rank-dependent exceptional ratios, their binomial corrections and the
explicit isogeny comparisons. The defining-only degree is an immediate
product of quadratic coordinate equations, not a claimed new result.
The type C result is a direct no-boundary application once (12) is
identified. Type B requires the additional critical-value count and
the exclusion of every hypersimplex subface. These claims are not
described as questions explicitly posed by Kostant--Michor.

Targeted searches on 2026-10-22 for defining/adjoint generalized Cayley
degrees and root-polytope volume formulas found the cited standard
inputs, but did not locate a published statement of these mixture
formulas. This is a limited search, not a priority certificate.
Neither theorem classifies arbitrary classical representations or
virtual coefficients, and positive-dimensional face singularities
still require other methods.

**Recommended next goal:** determine the non-finite target locus at
the exceptional type B ratios using the explicit trace-coordinate
equations (4). The degree formula now counts a general fiber, but does
not say where solutions can escape to infinity or special fibers gain
dimension. Success would be equations for that target locus, with a
proof distinguishing actual escaping branches from fixed base points
of a compactification. This advances the geometry of these same maps
rather than extending a degree table by more examples.
