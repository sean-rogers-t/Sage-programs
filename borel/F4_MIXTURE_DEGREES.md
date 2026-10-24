# Degrees of all minimal-adjoint mixtures in F4

2026-10-22. A characteristic-zero proof for a two-parameter family, with
independent character, lattice, volume and original-fiber calculations.
These are deductions recorded here from published representation theory
and toric intersection theory. No priority or external-review claim is made.

## 1. The theorem and its precise scope

Let G=F4(C), and use the numbering in the Kumar-paper notes: V(omega_4)
is the 26-dimensional minimal representation and V(omega_1) is the
52-dimensional adjoint representation. For integers a,b >= 0, a+b > 0,
put

\[
 \rho_{a,b}=V(\omega_4)^{\oplus a}\oplus V(\omega_1)^{\oplus b}.
\]

Let theta_(a,b):G -> g be the Springer/generalized Cayley morphism
defined by trace projection for this direct-sum representation.

**Theorem.** Its generic degree is

\[
 \boxed{
 d(a,b)=
 \begin{cases}
 96,&b=0,\quad a>0,\\
 240,&b>0,\quad a=b,\\
 336,&b>0,\quad a=5b,\\
 384,&b>0,\quad a\notin\{b,5b\}.
 \end{cases}}
 \tag{1}
\]

This is also the degree of the torus restriction T -> t and of the
induced map T/W -> t/W. It is not the minimum Cayley degree of F4.
Everything is over C, with the actual lattice X*(T)=P(F4)=Q(F4).
No parabolic, integral Borel image or cohomology theory enters (1).

The proof covers all nonnegative integral a,b, rather than extrapolating
from the fiber examples. Multiplying both multiplicities by the same
positive integer leaves the normalized map itself unchanged.
Adding trivial representations also leaves it unchanged.

## 2. Character and trace normalization

Work in epsilon coordinates with orthonormal e_1,...,e_4. Write

\[
 S=\sum_{w\text{ short root}}e^w,
 \qquad L=\sum_{w\text{ long root}}e^w.
\]

The 24 short roots are the eight +/-e_i and the sixteen half-sign sums
(+/-e_1 +/-e_2 +/-e_3 +/-e_4)/2. The 24 long roots are +/-e_i +/-e_j
for i<j. The minimal character is S+2, while the adjoint character is
S+L+4. Thus, with c=a+b,

\[
 \chi_{a,b}=cS+bL+(2a+4b),\qquad
 M_{a,b}=(6a+18b)I_4,
 \tag{2}
\]

and

\[
 \theta_{a,b}|_T=
 \frac{c\sum_{w\text{ short}}w e^w+
             b\sum_{w\text{ long}}w e^w}{6a+18b}.
 \tag{3}
\]

In particular the map of the mixture is a trace-weighted combination:

\[
 \theta_{a,b}=\frac{a\theta_{\min}+3b\theta_{\mathrm{ad}}}{a+3b}.
 \tag{4}
\]

Equation (4) holds on G as well as T, by the defining trace pairing.
For degree purposes the nonzero scalar in (3) can be removed, leaving
the logarithmic gradient of cS+bL. Dominance follows from
d(theta_(a,b))_1=Id. Constants in the character do not contribute.

The lattice is M=Z^4+Z*delta, delta=(1/2,1/2,1/2,1/2), of Euclidean
covolume 1/2. F4 has trivial center. All volumes and normal orders below
are already in M; there is no further central factor.

The construction and equality of the three generic degrees are published
inputs; see [Kostant--Michor, Corollary 3.3](https://arxiv.org/pdf/math/0109066).
The family calculation below supplies the particular values in (1).

## 3. The Newton polytope and all proper faces

For b=0 the normalized map is the minimal map, whose degree 96 is proved
in [the minimal-representation note](F4_MINIMAL_DEGREE.md).
Assume b>0 from now on. Since c=a+b>0, the support consists of all
48 roots, and its convex hull Delta is the long-root 24-cell.
The short roots are exactly the centers of its 24 octahedral facets.

The supporting hyperplanes are

\[
 x_i=\pm1\quad(8\text{ facets}),\qquad
 \pm x_1\pm x_2\pm x_3\pm x_4=2\quad(16\text{ facets}).
\]

The dual lattice is N={n in Z^4 : sum n_i is even}. Consequently the
primitive normals are +/-2e_i and the sixteen sign vectors, respectively.
Every facet has primitive height TWO. Treating e_i as a primitive normal
would give an incorrect boundary correction.

The Weyl group acts transitively on these 24 facets: their normal
directions are the short-root directions. For the representative x_1=1,
the vertices are e_1 +/- e_j (j=2,3,4), with short-root center e_1.
All lower-dimensional faces are simplices; they contain no other support
weights. Their characters cannot have singular zeros, by affine
independence of their exponent vectors.

The normalized volume is 384. An explicit computation avoids a lattice
convention ambiguity: triangulate the representative octahedron along
one opposite pair. The four tetrahedra, coned from zero, each have
epsilon-coordinate determinant of absolute value 2. Dividing by the
lattice covolume 1/2 gives four normalized volumes of 4. Every facet
therefore contributes 16, so

\[
 \operatorname{Vol}_M(\Delta)=24\cdot16=384.             \tag{5}
\]

The verifier also calculates this volume from the Weyl orbit of omega_1
in the fundamental-weight basis, independently of the epsilon basis.

## 4. Exactly two exceptional ratios

Use tangent coordinates x=e^(e_2), y=e^(e_3), z=e^(e_4) on the facet
x_1=1, and translate by its center. Its character becomes

\[
 g=c+b(x+x^{-1}+y+y^{-1}+z+z^{-1}).                      \tag{6}
\]

These are coordinates on the actual face torus: M intersect {x_1=0}
is Z e_2 + Z e_3 + Z e_4. No finite-cover division is needed.
The logarithmic critical equations give x,y,z in {+1,-1}.
At a sign triple epsilon the value is

\[
 c+2b(\epsilon_1+\epsilon_2+\epsilon_3).
\]

Since c/b=(a+b)/b >= 1, a critical value is zero exactly in these cases:

| Ratio | Sign triples at the singular zero | Points per facet |
|---|---|---:|
| a=b, so c=2b | Two minus signs and one plus sign | 3 |
| a=5b, so c=6b | Three minus signs | 1 |

There are no other ratios for nonnegative multiplicities. Equivalently,
putting r=c/b, the product over the eight critical values, divided by b^8,
is

\[
 (r-6)(r-2)^3(r+2)^3(r+6).
\]

The logarithmic Hessian of g is diagonal with entries 2b*epsilon_i;
it is invertible. All these singular zeros are isolated Morse points.
Multiplication by a monomial unit preserves the singular-zero condition
and nondegeneracy of the Hessian at such a point.

## 5. Why each boundary point has length two

The normal order requires the actual character lattice and the next
weight layer. Put

\[
 s=e^{-\delta},\quad x=e^{e_2},\ y=e^{e_3},\ z=e^{e_4},
 \quad D=xyz,\quad A=(1+x)(1+y)(1+z),
\]

\[
 H=x+x^{-1}+y+y^{-1}+z+z^{-1},
 \qquad K=\sum_{u<v\text{ in }\{x,y,z\}}(u+u^{-1})(v+v^{-1}).
\]

The four characters (-delta,e_2,e_3,e_4) form a basis of M. Thus s is
a primitive local parameter for the facet divisor, not a ramified-cover
coordinate. Since e^(e_1)=s^-2/D, direct substitution gives

\[
 S=s^{-2}/D+s^{-1}A/D+H+sA+s^2D,
 \qquad L=(s^{-2}/D+s^2D)H+K.
\]

After clearing the pole of order two,

\[
 F=s^2(cS+bL)=f_0+s f_1+s^2f_2+s^3f_3+s^4f_4,
\]

\[
 f_0=(c+bH)/D,\quad f_1=cA/D,\quad f_2=cH+bK,
 \quad f_3=cA,\quad f_4=D(c+bH).                        \tag{7}
\]

At every singular zero in Section 4, at least two of x,y,z are -1.
Consequently A and all its first derivatives vanish there:

\[
 f_1(p)=0,\qquad df_1(p)=0.                            \tag{8}
\]

For an arbitrary target eta in the coordinate basis dual to (s,x,y,z),
the cleared gradient equations near p are

\[
 D_{x,y,z}F=\eta_{\mathrm{tan}}s^2,
 \qquad sF_s-2F=\eta_s s^2.                            \tag{9}
\]

The invertible tangent Hessian and formal implicit function theorem
solve the first three equations uniquely for u(s)=(x(s),y(s),z(s)).
By (8), u(s)-p=O(s^2). Therefore f_0(u(s))=O(s^4) and
f_1(u(s))=O(s^4). In the remaining equation the direct s^2 f_2 term
cancels identically from sF_s-2F. Its leading generic term is

\[
 -\eta_s s^2+O(s^3).
\]

For a general target eta_s is nonzero, so the local intersection length
is exactly two. This holds simultaneously at all the boundary points,
since there are only finitely many facet normals to avoid.

There are no lower-face zeros, by Section 3. For a general target the
torus fiber is finite and reduced by dominance in characteristic zero.
The entire compactified intersection is therefore zero dimensional and
supported in the smooth torus and facet orbits. Its length is the toric
top intersection 384. Subtracting the boundary lengths yields

\[
 \begin{aligned}
 a=b:&\quad384-24\cdot3\cdot2=240,\\
 a=5b:&\quad384-24\cdot1\cdot2=336,\\
 \text{other }a\ge0, b>0:&\quad384-0=384.
 \end{aligned}
\]

This proves (1). The intersection input is standard toric geometry;
see Cox--Little--Schenck, Theorem 13.4.3
([authors' book page](https://dacox.people.amherst.edu/toric.html)),
and the precise local-intersection argument in
[the general degree note](GENERAL_CAYLEY_DEGREE.md#3-isolated-boundary-zeros-the-normal-layer-correction).

## 6. What generalizes, and what does not

**A reusable height-two lemma.** Let a dominant torus gradient have a
full-dimensional Newton polytope containing zero in its interior.
Suppose every boundary zero lies in a facet orbit, is isolated and has
invertible tangent Hessian. At a primitive-height-two facet write the
cleared character as f_0+s f_1+s^2 f_2+... . If f_1(p)=df_1(p)=0 at
each boundary zero, its generic intersection length is two. When all
facets carrying zeros satisfy these hypotheses, the degree is normalized
volume minus twice the number of boundary zeros on the actual face tori.
The proof is (9); no exceptional-group hypothesis is used.

**Coefficient families.** A centered crosspolytope face with character
c+b sum_(i=1)^d(z_i+z_i^-1) becomes singular precisely when
c/b=-2 sum epsilon_i. Thus candidate degree-drop ratios and point counts
can be read from sign sums. One must still check the lattice, all other
faces, Hessians and normal layers. F4 illustrates why the ratios alone
do not determine the correction: the primitive height is two and the
intermediate short-root layer has a special vanishing property.

**Finite-flat corollary away from the exceptional ratios.** For b>0 and
a not in {b,5b}, the no-boundary argument makes T -> t finite flat of
rank 384; the same holds for T/W -> t/W. The minimal-only case has rank
96. This follows from the proper graph argument and miracle flatness
in [the minimal note](F4_MINIMAL_DEGREE.md#5-what-generalizes).
We make no global-finiteness assertion at the two exceptional ratios,
or for the full group map. Boundary base points alone would not settle
those separate questions.

Nothing here proves a degree formula for arbitrary F4 irreducibles,
other exceptional groups, or arbitrary complex/virtual multiplicities.
In particular, negative ratios can introduce different layer behavior
or a degenerate trace form. Positivity of actual multiplicities is a
stated hypothesis of (1).

## 7. Independent checks and their limits

The [verifier](tools/verify_f4_mixture_degree.py) checks:

- Both weight lists against independently computed highest-weight
  characters, including zero multiplicities, dimensions and trace forms.
- Facets, extra face weights, primitive heights, actual lattice charts
  and determinant volumes; a separate fundamental-weight orbit gives
  normalized volume 384 again.
- Exact rational identities (7), the critical-value polynomial, tangent
  Hessians and first normal-layer jets at the two exceptional ratios.
  Scaled pairs (2,2) and (10,2) check the same local calculations.
- Original saturated torus fibers from two independent equation
  constructions at the SAME normalized epsilon target (1,2,3,5):
  one enumerates weights in the basis (e_1,e_2,e_3,delta), while the other
  differentiates the factored product character described below.

For the latter construction set z_i=e^(e_i), u=e^delta, u^2=product z_i.
Then

\[
 S=\sum_i(z_i+z_i^{-1})+\frac{\prod_i(1+z_i)}{u},
 \qquad L=\sum_{i<j}(z_i+z_i^{-1})(z_j+z_j^{-1}).
\]

Differentiate in epsilon coordinates, divide by 6a+18b and saturate by
all torus coordinates. The relation u^2=product z_i presents the actual
torus and must not be interpreted as a final extra factor two.

The cases (a,b)=(0,1),(1,1),(5,1),(2,1),(1,0) are checked over both
F_32003 and F_32009, in both models: twenty fibers in total. Their
lengths are respectively 384,240,336,384,96. An independent integer
staircase count agrees with the Sage quotient basis. Multiplication by
u+2z_1+3z_2 has squarefree characteristic polynomial in every fiber,
certifying geometric reducedness; the entire characteristic polynomials
agree between the two models for each case and prime.

These finite-field calculations corroborate the separate proof; they
do not establish characteristic-zero generic degrees by specialization.
The two constructions share Sage/Singular arithmetic, so this is not
an independent software-engine audit. See the
[machine-readable report](f4-mixture-degree-verification.json).

```bash
PYTHONPATH=src sage -python tools/verify_f4_mixture_degree.py
```

## 8. Literature status and next goal

Published inputs are the trace-projection construction and equivariant
degree reduction of Kostant--Michor, the standard weight systems, toric
intersection theory and the formal implicit function theorem. The
specific all-multiplicity formula (1), its two exceptional ratios and
their corrections are deductions proved above. They are not described
as open questions explicitly posed by Kostant--Michor. Targeted searches
on 2026-10-22 for F4 minimal/adjoint generalized Cayley degrees did not
locate this family formula; that limited search is not a novelty audit.

**Recommended next goal:** formulate a coefficient-discriminant theorem
for representation families with isolated Morse facet zeros. The purpose
is to predict degree jumps from face critical values and primitive normal
jets before solving full fibers. Success would be a precise, reusable
criterion determining the exceptional parameter locus and correction
multiplicities, with F4 as its complete worked example and an example
showing what changes when the first intermediate layer has nonzero jet.
This advances a general method rather than merely adding another group
to a numerical table. Positive-dimensional face singularities remain a
separate problem.
