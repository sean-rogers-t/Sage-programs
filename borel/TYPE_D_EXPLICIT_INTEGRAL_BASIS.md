# A constructive full type-D integral basis

29 October 2026. This is a constructive corollary of the existing
full even-tail gluing theorem, together with the explicit square-prefix
basis. It does not reprove the difficult relative Euler comparison.
No external-review or literature-priority claim is made.

## 1. Fixed source, degrees, and earlier dependencies

Let X=OG(a,2(a+k)), a>=1, k>=2. Use the normalized vector Springer inputs
q_i=e_i(x_tail^2) and p=product(x_tail), with p^2=q_k. Polynomial degree
d is cohomological degree 2d. Put H=k and K=a. These letters describe
the original tail, moved to the head on a complete flag for the
comparison argument; they do not interchange the full Euler problems.

The integral source modulo its kernel has the basis

    s_lambda(x_tail^2),        lambda in K^H,
    p s_mu(x_tail^2),          mu in (K-1)^H.

Thus as an abelian group it is R_(K,H) plus p R_(K-1,H). Products are
computed by ordinary rectangular Schur multiplication: the square
sector acts on the Euler sector through the projection R_(K,H) ->
R_(K-1,H), and two Euler factors multiply using p^2=e_H(y). Multiplication
by e_H identifies the narrower rectangle with its shifted image in
R_(K,H). This retains the actual Euler source and its narrower kernel.

The earlier [full gluing proof](TYPE_D_EVEN_GLUING_THEOREM.md) identifies
the coefficient-row lattice in these common coordinates as

    Lambda_D = M_0 + sum_(t>=1) 2^t N_t,
    N_t = P_t(H,K) plus P_t(H,K-1).

The degree shift H is understood on the second summand. M_0 is the
ordinary Grassmannian coefficient matrix of the displayed source,
after replacing tail variables by the H ordinary Chern roots on
Gr(H,H+K). The square-prefix theorem applies to the ACTUAL head H,
even when H>K. The full gluing proof handles even H. For odd H the
sectors occupy opposite degree parities; the primitive typed-difference
map in [the Euler comparison](TYPE_D_SMITH.md#exact-euler-sector-comparison-with-type-b)
and the square row equality give the same blockwise statement.

## 2. An explicit construction using only integer triangles and F2

Apply the [simultaneous basis construction](CLASSICAL_ADAPTED_SMITH_BASIS.md)
to each square-prefix factor, in each required source weight. Let V_0
be the block-diagonal row matrix and e_i its exponent labels. Empty
sectors are omitted; R_(0,H)=Z is the point convention. Then

    N_t = span_Z {the rows of V_0 with e_i<=t}.

Write Z_0={i:e_i=0} and Z_+={i:e_i>0}. Form C=M_0 V_0^(-1).
Every column of C indexed by Z_+ is zero INTEGERS, not merely modulo
two. Indeed both projections of M_0 lie in their separate P_0.
The square projection is its defining ordinary Grassmannian matrix.
For the Euler projection, multiplication by the ordinary top Chern
class shifts every Schur partition by (1^H); it identifies the
smaller rectangle and kills terms outside it. Thus its row module
is exactly P_0(H,K-1).

Reduce C restricted to Z_0 modulo two and take reduced row echelon
form. If its rank is r, lift its r nonzero rows by choosing entries
0 and 1. Complete them by the coordinate unit rows at nonpivot
columns. This is an INTEGRAL UNIMODULAR basis of Z^(Z_0): after
ordering pivot columns first its matrix is block unit triangular.
Adjoin the unchanged coordinate rows indexed by Z_+ and call the
result U. Set V=U V_0. Give its rows the labels

    0 for the lifted echelon rows;
    1 for the nonpivot complement in Z_0;
    e_i for the unchanged positive-level rows.

Call these labels epsilon_j. This construction requires no Smith
normal form, Weyl enumeration, or arbitrary-rank search. It is a finite
formula in ordinary symmetric functions followed by binary elimination.
The source partitions and a fixed pivot order make it deterministic.

**Full basis theorem.** If f_j is the source element whose coefficient
column is column j of V^(-1), then

    b_j = beta(f_j)/2^epsilon_j

are an integral basis of the saturation of the FULL normalized D image.
The image itself has basis 2^epsilon_j b_j. The b_j may mix square and
Euler terms; the statement does not identify them with single Schubert
classes or declare their Schubert coefficients positive.

## 3. Why binary data determine the full integer answer

The earlier row formula becomes

    Lambda_D = M_0 + 2 N_0
               + direct_sum_(e_i>0) 2^e_i Z (V_0)_i.

This follows because the positive prefix sum contributes 2N_0 and
weight 2^e_i to each other basis row. In particular, the unit-level
problem is EXACTLY the full inverse image of row(C mod 2) inside N_0.
There can be no additional odd index or undiscovered higher power
of two at this step: 2N_0 is already present. The lifted echelon
rows and twice their complement therefore generate M_0+2N_0 over Z.
Hence

    Lambda_D = direct_sum_j 2^epsilon_j Z V_j.

Dualizing the row lattice proves the claimed integral saturation basis,
just as in the square case. This argument uses the full earlier
row-lattice identity; applying binary elimination to an unrelated
matrix without that identity would not prove its integral Smith form.

The number of complement rows is the previously proved unit-gluing
multiplicity. Thus this gives explicit vectors realizing the old
formula, not a new multiplicity formula. For odd H the complement
is empty. At even H it consists precisely of the old unit-level
gluing corrections; no positive-level basis vector is changed.

## 4. Products and an independent binary description

Expand products in the integral source basis f_j using p^2=e_H(y)
and the two rectangular Schur rules. If f_i f_j=sum_l C_ij^l f_l,
then

    b_i b_j = sum_l 2^(epsilon_l-epsilon_i-epsilon_j) C_ij^l b_l.

These coefficients are integers because the saturation is a subring.
This gives a constructive ring presentation by an additive basis and
its entire multiplication tensor, not yet a compact minimal list of
generators and relations.

There is also a separate implementation of the binary gluing matrix.
In characteristic two, the ordinary source map is

    (f,g) -> f^2 + e_H g^2 in R_(K,H) tensor F2,

where the g term is evaluated in the narrower rectangle before the
e_H shift. Thus the binary matrix can be built from ordinary Pieri
products and Frobenius, with no plethysm conversion. This is the
existing [Frobenius description](TYPE_D_FROBENIUS.md), used here to
independently reconstruct the chosen integral basis.

## 5. Verification and limits

[Sage verification](type-d-adapted-basis-verification.json) constructs
124 nonempty degree components across twelve block pairs. In each it
checks the FULL integer row-module equality, not just Smith counts.
It also reconstructs 28 actual D divided-difference matrices for four
block pairs, checking that the divided columns are integral, primitive,
and reconstruct every original column.

[Native replay](type-d-adapted-basis-audit.json) independently implements
permutation Jacobi--Trudi, Pieri products, Frobenius, and binary Gaussian
elimination. Integer Hermite forms check primitivity of the exported
actual D columns. This shares the actual D matrix input with the Sage
calculation; it does not claim a second independent D BGG engine.
Both implementations also compute every product of adapted classes
in seven small block pairs, using separate symmetric-function and
Pieri arithmetic, and obtain identical counts of nonzero integral
structure constants. These checks include the narrower Euler quotient
and p^2=e_H(y), rather than treating the sectors as a direct-product ring.

The additional [Chevalley audit](type-d-basis-chevalley-audit.json)
independently reconstructs all 28 actual D matrices from Cartesian roots,
signed permutations and complete-flag divisor multiplication. It imports
neither Sage nor a root library nor BGG. Exact coefficient-row multisets
agree, so this check closes the shared-BGG-input limitation of the first
arithmetic replay for these 28 components.

The arbitrary-rank proof is Sections 1--3 with the cited prior theorem,
not the number of computations. General representations, full integral
Levi-character sources, arbitrary partial flags, a geometrically positive
basis and novelty remain outside the conclusion. The original normalized
B/C/D single-tail arithmetic is now constructive at the basis level.
