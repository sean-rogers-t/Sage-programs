# The closed full type-D formula

**Status, 19 October 2026: complete internally audited repository deduction.**
The [general even-tail gluing theorem](TYPE_D_EVEN_GLUING_THEOREM.md)
closes the remaining condition in this note. Together with the earlier
square and Euler comparisons, it proves every multiplicity in the stated
scope. The historical filename is retained for existing links and
experiments. This is unpublished work, without an originality claim or
external expert review; see the [full audit](TYPE_D_EVEN_GLUING_AUDIT.md).

## 1. The fixed image problem

Use \(X_D=\operatorname{OG}(a,2(a+k))\), \(a\ge1,k\ge2\), and
the normalized inputs
\[
q_i=e_i(x_{\rm tail}^2),\qquad p=\prod x_{\rm tail},\qquad p^2=q_k.
\]
Polynomial degree \(d\) is cohomological degree \(2d\).
All factors refer to the integer input image inside its selected
rational image's saturation, with unit factors retained.
No assertion concerns the index in all ambient cohomology.

For the proved type-B multiplicities define
\[
\mathcal B_{a,k,e}(q)
=q^{e\lceil(a+k)/2\rceil}
{\lfloor(a+k)/2\rfloor\brack e}_q
U_{a+k-2e,\min(a,k)-e}(q).
\tag{1}
\]
It is zero for \(e\) outside \(0\le e\le\min(a,k)\).
The variable records polynomial degree \(2m\), as in the
[general B proof](TYPE_B_GENERAL_SMITH.md). At the point boundary set
\(\mathcal B_{0,k,0}=1\), \(\mathcal B_{0,k,e}=0\) for \(e>0\).
Write \(\mathcal U_{a,k}=\mathcal B_{a,k,0}\).

Put \(\mathcal I_{a,k}=0\) for odd \(k\). For \(k=2h\), put
\[
u=\lceil a/2\rceil,\quad v=\lfloor a/2\rfloor,\qquad
\mathcal I_{a,2h}(q)=
q^{h(v+1)}{h+u-1\brack h}_q.
\tag{2}
\]
The independent [unit theorem](TYPE_D_FROBENIUS.md) proves that this
counts the overlap between the two original image reductions modulo 2.
The subsequent gluing theorem proves that these classes exhaust the
gluing between their saturated lattices as well.

## 2. The full answer

Let \(M_D(a,k;e,d)\) count the factor \(2^e\) in polynomial
degree \(d\). The proved formula is
\[
\boxed{
\begin{aligned}
\sum_d M_D(a,k;e,d)z^d\ =\ &
\mathcal B_{a,k,e}(z^2)
+z^k\mathcal B_{a-1,k,e}(z^2)\\
&+(\delta_{e,1}-\delta_{e,0})\mathcal I_{a,k}(z^2).
\end{aligned}}
\tag{3}
\]

**Equation (3) is proved for every \(a\ge1,k\ge2,e,d\).**
Its \(e=0\) case was independently proved earlier.
Equivalently, the complete D unit polynomial is
\[
\boxed{
\mathcal U_{a,k}(z^2)
+z^k\mathcal U_{a-1,k}(z^2)
-\mathcal I_{a,k}(z^2).
}\tag{4}
\]
The entire square-sector Smith list now equals the B list in every rank,
by [TYPE_D_SQUARE_THEOREM.md](TYPE_D_SQUARE_THEOREM.md).
The entire Euler-sector Smith list, including its higher powers, is
already proved to equal the shifted smaller B list.

For odd \(k\), the two sectors have opposite degree parities and
there is no same-degree gluing. For even \(k\), formula (3) shows
that each overlap replaces a pair of unit factors by \(1,2\).
The relative incidence comparison supplies the integral proof that
every overlap has this form.

## 3. The even-tail integral step

Let \(L_\pm\) be the square and Euler image lattices in a degree,
\(S_\pm\) their respective saturations, and
\[
S=((L_+\oplus L_-)\otimes\mathbb Q)\cap H^{2d}(X_D;\mathbb Z),
\qquad W=S/(S_+\oplus S_-).
\]
The existing [two-flag gluing lemma](ORTHOGONAL_SMITH.md#lattice-lemma)
embeds \(W\) in both \(S_\pm/2S_\pm\).
Let \(W_0^\pm\) be the pullback of the reductions of \(L_\pm\)
in these two spaces.

The square comparison is now a proved theorem: every normalized D
square-sector list equals its B list, by the relative incidence argument
in [TYPE_D_SQUARE_THEOREM.md](TYPE_D_SQUARE_THEOREM.md). The Euler-sector
list is the shifted smaller B list by its separate integral comparison.
The condition needed to complete (3) when \(k\) is even was:

**Exhaustion by unit gluing:** \(W=W_0^+\cap W_0^-\).
This is now the theorem in
[TYPE_D_EVEN_GLUING_THEOREM.md](TYPE_D_EVEN_GLUING_THEOREM.md).
Its proof moves the specified tail to the head and compares the precise
square-plus-head-Euler source with type B. A relative Euler extension
has a primitive selected source; an integral parity retraction and the
equivariant eta denominator bound yield the reverse row containment.
The selected B calculation uses primitive actual-head prefixes and
signed level-zero overlap identities. The resulting full index,
compared with the two independently known sector indices, proves that
\(\dim W\) equals the coefficient of (2). The unit theorem already
supplies that many independent unit-gluing directions.

Consequently every gluing pair has exponents \((0,0)\), and
there are exactly the overlaps counted by (2). The integral gluing lemma
replaces each such pair by \((0,1)\), giving the correction in (3).
This proves the formula for even tails.

For odd \(k\), square inputs have even polynomial degree and Euler
multiples have odd polynomial degree. The two sectors therefore never
occur in the same homogeneous component. The two now-proved sector
comparisons give (3) with no gluing hypothesis, for every head size and
every degree. This is a complete repository deduction for odd tails.

The mod-2 unit theorem alone did not establish the square comparison;
the new relative integral argument supplies it. The unit count alone did not identify
positive-level gluing; the separate integral comparison is necessary.

## 4. Checks and concrete examples

For \(a=k=2\), degree \(4\), the square list is \((1,2)\)
and the Euler list is \((1)\). The proved overlap has dimension one.
The full list is \((1,2,2)\); this agrees with (3).

For \(a=5,k=4\), degree \(16\), the proved unit overlap has
dimension two. It has explicit Schur-square representatives indexed
by \(\mu=(5,3),(4,4)\) in the ordinary Grassmannian model.
The original fundamental-weight BGG calculation at this degree timed
out. The new independent Cartesian calculation now completes it, with
full saturated gluing dimension two. Degree 14, the other old rank-nine
timeout, also completes, with gluing dimension one.

The [pure-Python candidate verifier](experiments/verify_type_d_candidate.py)
imports no Sage, Borel engine, or prior verifier. Its
[report](experiments/type_d_candidate_verification.json) compares (3)
with 70 historical and 14 fresh exact full Smith lists, and separately
checks the corresponding square/Euler lists and gluing dimensions.
All 84 comparisons pass; the two historical rank-nine timeouts remain
explicitly excluded from that historical data set. The source and data
byte hashes are recorded.

The [new full-gluing verifier](experiments/verify_type_d_even_gluing.py)
and [report](experiments/type_d_even_gluing_verification.json) compute
132 complete head-three components and six larger Cartesian components,
with 16 independent full row-module cross-checks, 27 explicit primitive
unit-gluing representatives, and both old timeouts resolved. They compute
the saturated gluing dimension from three actual Smith indices,
separately from the original mod-two overlap. All 138 comparisons pass.

The [ordinary-Grassmann experiment](experiments/type_d_frobenius_search.py)
and [its report](experiments/type_d_frobenius_search_results.json)
corroborate the proved unit formula using 189 degree components,
96 explicit intersection identities, and the exact BGG data.
The separate
[upper-bound verifier](experiments/verify_type_d_unit_upper_bound.py)
and [its report](experiments/type_d_unit_upper_bound_verification.json)
test the proved derivation/core-one argument on 147 components, with
40 exact integer row-module/primitivity checks and 41 signed SXP rows.
Historical uses of the word candidate in experiment names do not change
the scope of their finite evidence or replace the subsequent general proof.

From the repository root:

```sh
python borel/experiments/verify_type_d_candidate.py
```

For the literal defining/vector character form, multiply every factor
in degree \(d\) by \(2^d\). All formulas above use the fixed normalized
form from [TYPE_D_SMITH.md](TYPE_D_SMITH.md).
