# A closed full type-D formula: odd tails proved and even-tail gluing open

**Status, 19 October 2026:** the unit multiplicities in this note are
proved in [TYPE_D_FROBENIUS.md](TYPE_D_FROBENIUS.md). The
[general square comparison](TYPE_D_SQUARE_THEOREM.md) now determines
every square-sector factor. Together with the Euler comparison this
proves the full formula for every odd tail. For even tails, the higher
multiplicities remain a **conjecture** subject to the single gluing
hypothesis below. It matches 84 exact recorded full Smith lists;
two timed-out calculations are excluded. This is unpublished repository
work, without an originality claim or external expert review.

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
The [unit theorem](TYPE_D_FROBENIUS.md) proves that this counts exactly
the overlap between the two **original image** reductions modulo 2,
not yet all gluing between their saturated lattices.

## 2. The proposed full answer

Let \(M_D(a,k;e,d)\) count the factor \(2^e\) in polynomial
degree \(d\). The candidate is
\[
\boxed{
\begin{aligned}
\sum_d M_D(a,k;e,d)z^d\ \stackrel{?}{=}\ &
\mathcal B_{a,k,e}(z^2)
+z^k\mathcal B_{a-1,k,e}(z^2)\\
&+(\delta_{e,1}-\delta_{e,0})\mathcal I_{a,k}(z^2).
\end{aligned}}
\tag{3}
\]

**What is proved in every rank in (3):** all exponents for odd \(k\),
and the \(e=0\) assertion for every \(a,k,d\).
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
there is no same-degree gluing. For even \(k\), formula (3) predicts
that each overlap replaces a pair of unit factors by \(1,2\).
This interpretation is the conditional argument below, not a proof
that every integral overlap has this form.

## 3. The remaining even-tail integral step

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
Thus the following is the sole remaining sufficient hypothesis for (3)
when \(k\) is even:

**Exhaustion by unit gluing:** \(W=W_0^+\cap W_0^-\).
The unit theorem proves the dimension of this intersection is the
relevant coefficient of (2). Equality with the full \(W\) is still
unproved: saturation can introduce additional mod-2 classes.

Under this hypothesis every gluing pair has exponents \((0,0)\), and
there are exactly the overlaps counted by (2). The integral gluing lemma
replaces each such pair by \((0,1)\), giving the correction in (3).
This proves the conditional statement for even tails.

For odd \(k\), square inputs have even polynomial degree and Euler
multiples have odd polynomial degree. The two sectors therefore never
occur in the same homogeneous component. The two now-proved sector
comparisons give (3) with no gluing hypothesis, for every head size and
every degree. This is a complete repository deduction for odd tails.

The mod-2 unit theorem alone did not establish the square comparison;
the new relative integral argument supplies it. Nor does a unit count
identify all positive-level gluing for even tails.

## 4. Checks and concrete examples

For \(a=k=2\), degree \(4\), the square list is \((1,2)\)
and the Euler list is \((1)\). The proved overlap has dimension one.
The full list is \((1,2,2)\); this agrees with (3).

For \(a=5,k=4\), degree \(16\), the proved unit overlap has
dimension two. It has explicit Schur-square representatives indexed
by \(\mu=(5,3),(4,4)\) in the ordinary Grassmannian model.
The direct full integer BGG calculation at this degree timed out;
it is not counted as a full Smith verification.

The [pure-Python candidate verifier](experiments/verify_type_d_candidate.py)
imports no Sage, Borel engine, or prior verifier. Its
[report](experiments/type_d_candidate_verification.json) compares (3)
with 70 historical and 14 fresh exact full Smith lists, and separately
checks the corresponding square/Euler lists and gluing dimensions.
All 84 comparisons pass; two fresh rank-nine timeouts are explicitly
excluded. The source and data byte hashes are recorded.

The [ordinary-Grassmann experiment](experiments/type_d_frobenius_search.py)
and [its report](experiments/type_d_frobenius_search_results.json)
corroborate the proved unit formula using 189 degree components,
96 explicit intersection identities, and the exact BGG data.
The separate
[upper-bound verifier](experiments/verify_type_d_unit_upper_bound.py)
and [its report](experiments/type_d_unit_upper_bound_verification.json)
test the proved derivation/core-one argument on 147 components, with
40 exact integer row-module/primitivity checks and 41 signed SXP rows.
The ordinary-Grassmann experiment's historical use of the word candidate does not change the scope of
the subsequent written proof or promote the full Smith conjecture.

From the repository root:

```sh
python borel/experiments/verify_type_d_candidate.py
```

For the literal defining/vector character form, multiply every factor
in degree \(d\) by \(2^d\). All formulas above use the fixed normalized
form from [TYPE_D_SMITH.md](TYPE_D_SMITH.md).
