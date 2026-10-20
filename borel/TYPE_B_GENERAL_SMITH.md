# A general type-B Smith multiplicity formula

**Status: complete repository derivation, 19 October 2026.** This note
supplies the integral comparison missing from the earlier conjecture.
Its arguments have been checked within this research session, including
independent reviews and exact finite calculations. It is unpublished,
has not received external expert review, and makes no originality claim.
The published inputs and their relationship to Parts II and III remain
as recorded in [PROOF_STATUS.md](PROOF_STATUS.md).

The proof uses the already audited integral Schur-product description of
the symplectic prefixes. It does not infer saturation from a rational
character identity, or assume the missing Chern-generation statement.
Instead a Chow envelope supplies additional integral rows; a triangular
minor and a matching geometric rank bound determine the entire lattice.

## 1. Statement and conventions

Put \(n=a+k\), \(\ell=\min(a,k)\), \(A=\max(a,k)\),
\(r=\lfloor n/2\rfloor\), \(c=\lceil n/2\rceil\),
and \(\epsilon=c-r\in\{0,1\}\).
Use the normalized inputs \(q_i=e_i(x_{\rm tail}^2)\) on
\(\operatorname{OG}(a,2n+1)\), as in [TYPE_B_SMITH.md](TYPE_B_SMITH.md).
All indices are in the saturation of the selected rational image.
Write \(G(v,u)={v\brack u}_q\), zero outside \(0\le u\le v\), and
\[
U_{n,L}(q)=\sum_{j=0}^{\lfloor L/2\rfloor}
q^{j(\lceil n/2\rceil-L+j)}
G(\lfloor n/2\rfloor,j)G(\lceil n/2\rceil-j,L-2j),
\qquad U_{0,0}=1.
\]

**Deduction.** If \(M_{n,\ell,e}(m)\) is the multiplicity of the nonzero
Smith factor \(2^e\) in polynomial degree \(2m\), then
\[
\boxed{\quad
\sum_m M_{n,\ell,e}(m)q^m
=q^{ec}G(r,e)\,U_{n-2e,\ell-e}(q),
\qquad 0\le e\le\ell.\quad} \tag{1}
\]
The parameters inside \(U\) are recomputed at its arguments.
Polynomial degree \(2m\) is cohomological degree \(4m\).
The original image is zero in odd polynomial degrees.
For the literal defining/vector character form, multiply every factor
in polynomial degree \(d\) by \(2^d\), as in the
[normalization audit](INTEGRAL_LITERATURE.md).

The polynomial expression itself has existing counting and representation
theoretic interpretations; Section 8 distinguishes those from (1).
The assertion here concerns the additional integral Springer-image
interpretation.
The partition statistic below counts multiplicities; it does not assert
that each original Schur generator is itself an individual Smith vector.

## 2. The precise lattice problem

By integral block duality put the smaller block at the head.
In degree \(2m\) use source basis
\(s_\lambda(y_1^2,\ldots,y_\ell^2)\),
\(\lambda\subseteq A^\ell\), \(|\lambda|=m\).
The C matrix has full primitive row lattice \(M_m\), and the B matrix
has rows \(2^{\nu(w)}C_w\). Identify \(M_m\) with the degree-\(m\)
part of the ordinary Grassmann ring \(R_{A,\ell}\), by the Schur basis.
For \(P_t=\operatorname{span}_{\mathbb Z}\{C_w:\nu(w)\le t\}\), the
[prefix theorem](TYPE_B_FILTRATION.md) gives
\[
P_t=\operatorname{span}_{\mathbb Z}
\{s_\mu s_\alpha s_\beta\text{ in }R_{A,\ell}\}, \tag{2}
\]
where, for \(b=\ell-t\), \(p=\lceil b/2\rceil\), \(s=\lfloor b/2\rfloor\),
\[
\mu\subseteq A^t,\quad
\alpha\subseteq(c-p)^p,\quad
\beta\subseteq(r-s)^s.
\]
All products have total degree \(m\).

It suffices to prove every \(P_t\) primitive, with
\[
\sum_m\operatorname{rank}P_t(m)q^m
=\sum_{e=0}^t q^{ec}G(r,e)U_{n-2e,\ell-e}. \tag{3}
\]
Indeed a chain of primitive subgroups of a free abelian group admits
an adapted integral basis. The weighted row lattice is
\(\sum_t2^tP_t\), and its exponent-\(t\) multiplicity is
\(\operatorname{rank}P_t-\operatorname{rank}P_{t-1}\).

## 3. A geometric upper bound and an integral envelope

Choose a maximal-rank square-zero \(N:V\to V\) on \(\mathbb C^n\),
with \(I=\operatorname{im}N\), \(K=\ker N\), of dimensions \(r,c\).
Define
\[
W_t=\{U\in\operatorname{Gr}(\ell,V):
\operatorname{rank}(N:U\to V/U)\le t\}.
\]
The defect-stratum calculation in the
[candidate note](TYPE_B_CLOSED_CONJECTURE.md) gives free even cohomology
and the Betti polynomial on the right of (3).
Briefly, for defect \(e\), put \(D=NU\), \(B=U\cap K\),
and \(j=\dim(D\cap B)\). Choose
\[
D\in\operatorname{Gr}(e+j,r),\quad
C\in\operatorname{Gr}(j,D),\quad
(B/C)\bmod D\in\operatorname{Gr}(\ell-e-2j,K/D).
\]
The two affine lift spaces have total dimension
\(ec+j(c-\ell+j)\).
Grassmann Schubert filtrations and these affine bundles give the
compact-support polynomial
\[
\sum_j q^{ec+j(c-\ell+j)}
G(r,e+j)G(e+j,j)G(c-e-j,\ell-e-2j)
=q^{ec}G(r,e)U_{n-2e,\ell-e}.
\]
Closed rank filtrations make the free even ranks add. Since \(W_t\)
is projective, their sum computes its ordinary Betti numbers.

For \(b=\ell-t\), consider
\[
Y_{\rm bal}=\{(F,S,U):F\subset I,\ \dim F=s;\
F\subset S\subset N^{-1}F,\ \dim S=b;\
S\subset U,\ \dim U=\ell\}.
\tag{4}
\]
Over \(\operatorname{Gr}(s,r)\), the bundle
\(E=N^{-1}F/F\) has rank \(c\) and total Chern class \(1\):
\(0\to K/F\to E\to F\to0\).
Then choose \(S/F\in\operatorname{Gr}(p,E)\), followed by
\(U/S\in\operatorname{Gr}(t,V/S)\).
Its Schur branching rows are exactly (2).
The cohomology restriction
\(H^*(\operatorname{Gr}(\ell,n))\to H^*(W_t)\to H^*(Y_{\rm bal})\)
has these coefficient rows. Therefore, degreewise,
\[
\operatorname{rank}P_t(m)\le
\operatorname{rank}H^{2m}(W_t;\mathbb Z). \tag{5}
\]
This inequality does not require ambient Chern generation.

The proper map \(f:Y_{\rm bal}\to W_t\) is also a **Chow envelope**
over every extension field of \(\mathbb C\). Here is a short proof.
For \(U\) of defect at most \(t\), the invariant space
\(U\cap N^{-1}U\) has dimension at least \(b\).
A square-zero Jordan basis chooses an invariant \(b\)-plane \(S\)
inside it, over the same field.
If \(j=\operatorname{rank}N|_S\le\lfloor b/2\rfloor=s\), then
\[
\dim(S\cap I)\ge\max(j,b-j-\epsilon)\ge s.
\]
For the last inequality, if \(j=s\) use the first term; if \(j<s\)
use \(b-j-\epsilon\ge b-s\ge s\).
Choose \(NS\subset F\subset S\cap I\) of dimension \(s\).
This supplies a rational lift at the generic point of every irreducible
subvariety \(Z\subset W_t\). Closing its graph gives a subvariety of
\(Y_{\rm bal}\) birational onto \(Z\). Thus \(f_*\) surjects on
integral cycle classes, and hence on integral Chow groups.
The [defect comparison note](TYPE_B_DEFECT_COMPARISON.md) records
the stronger primitive cohomology pullback, but it is not needed for (5).

## 4. Extra integral rows from unbalanced incidence bundles

For every \(0\le j\le b\), define \(Y_j\) by the same incidence (4)
with \(\dim F=b-j\), \(\dim S/F=j\).
All these dimensions are permitted since \(b\le\ell\le r\).
Write \(g_j:Y_j\to W_t\) for the forgetful map and \(i:W_t\to
\operatorname{Gr}(\ell,V)\) for the ambient inclusion.
Its three integral Grassmann-bundle bases give coefficient rows
\[
s_\mu s_\alpha s_\beta,\quad
\mu\subseteq A^t,\quad
\alpha\subseteq(c-j)^j,\quad
\beta\subseteq(r-b+j)^{b-j}. \tag{6}
\]
The middle Grassmann bundle has total ambient Chern class \(1\), so
out-of-width middle Schur classes vanish. The top fiber has width \(A\);
its terms already fit since the source partition has width at most \(A\).
There are no omitted base-dependent fiber reductions.

**Every row in (6) belongs to \(P_t\) integrally.**
To see this in degree \(2m\), pair ambient Schur classes with an integral
cycle \(D\in\operatorname{CH}_m(Y_j)\) dual to a cohomology basis element.
Such dual cycles exist because an iterated Grassmann bundle has an
integral Schur Chow basis and a perfect integral intersection pairing.
Lift the cycle \((g_j)_*D\) on \(W_t\) through the Chow envelope \(f\).
The projection formula then expresses its pairings, and hence the
coefficient row (6), as an integer combination of rows on \(Y_{\rm bal}\).
Smoothness of \(W_t\) is unnecessary: ambient Chern classes act on its
cycles and satisfy the proper projection formula.

## 5. A unit triangular minor

The products \(h_\gamma\), \(\gamma\subseteq A^\ell\), form an integral
basis of \(R_{A,\ell}\). Kostka/Jacobi--Trudi triangularity and
\(h_i=0\) for \(i>A\) also allow \(s_\mu\) in (6) to be replaced by
\(h_\mu\), with \(\mu\subseteq A^t\).
For
\[
0\le j\le\left\lfloor\frac{b+\epsilon}{2}\right\rfloor,\qquad
\Delta_j=((c-j)^j,(r-b+j)^{b-j}),
\]
let \(\mathcal K_t\) be the union of the partition shapes
\[
\gamma\subseteq(A^t,\Delta_j). \tag{7}
\]
For \(b=0\), take \(\mathcal K_\ell=A^\ell\).
The range of \(j\) makes the displayed corner weakly decreasing;
all its parts lie between \(0\) and \(A\).

For each \(\gamma\in\mathcal K_t\), choose one permissible \(j\);
let \(\mu\) be its largest \(t\) parts, \(\alpha\) the next \(j\),
and \(\beta\) its remaining \(b-j\).
The integral row \(v_\gamma=h_\mu s_\alpha s_\beta\) of Section 4
has coefficient \(1\) on \(h_\gamma\). Every other surviving complete
coordinate has strictly larger energy
\(\sum_i\gamma_i^2\).
Indeed a nonidentity permutation in Jacobi--Trudi replaces
\(\alpha_i\) by \(\alpha_i-i+\sigma(i)\); its energy difference is
\[
2\sum_i\alpha_i(\sigma(i)-i)+\sum_i(\sigma(i)-i)^2>0.
\]
The first sum is nonnegative by rearrangement and the second is positive.
The same holds for \(\beta\), and multiplication adds energies.
Terms with a negative index or an index greater than \(A\) vanish.

Order all degree-\(m\) coordinates by energy, breaking ties arbitrarily.
Project the selected rows to their columns in \(\mathcal K_t\).
The resulting square matrix is triangular with diagonal \(1\).
Consequently their span \(Q_t(m)\) is a split, primitive sublattice of
\(M_m\), of rank \(\#\{\gamma\in\mathcal K_t:|\gamma|=m\}\).
It is not being asserted that \(P_t\) is a coordinate lattice.

## 6. Counting the pivot partitions

For \(\gamma\subseteq A^\ell\), its beta numbers
\(\beta_i=\gamma_i+\ell-i\) are an \(\ell\)-element subset of
\(\{0,\ldots,n-1\}\), in decreasing order. Put
\[
d=\#\{\beta_i<r\},\qquad
e(\gamma)=\#\{\beta_i\ge c+d\}. \tag{8}
\]
Then
\[
\boxed{\ \gamma\in\mathcal K_t\quad\Longleftrightarrow\quad
e(\gamma)\le t.\ } \tag{9}
\]

For clarity, if \(b=0\) this is automatic. Otherwise write \(J=\ell-d\).
If \(J\le t\), all bottom \(b\) beads are below \(r\), so the corner
\(j=0\) works. Since \(c+d\ge r\), also \(e(\gamma)\le J\le t\).
If \(J>t\), the bottom \(b\) beads have \(j_0=J-t=b-d\) entries
at least \(r\). The corner inequalities are
\[
J\le t+j,\qquad
\beta_{t+1}\le c+b-j-1.
\]
When \(j=0\) the latter is redundant; an empty second block likewise
imposes no condition. Since \(j\ge j_0\), any permissible corner
implies \(\beta_{t+1}<c+d\), equivalent to \(e(\gamma)\le t\).
Conversely this inequality places the \(j_0\) middle beads in
\([r,c+d-1]\), an interval of length \(d+\epsilon=b-j_0+\epsilon\).
Hence \(2j_0\le b+\epsilon\), and the corner \(j_0\) is permissible
and satisfies both inequalities. This proves (9), including empty
interval and \(b=1\) cases.

Count partitions of exact statistic \(e\). If there are \(d\) beads
in \([0,r-1]\), there are \(j=\ell-d-e\) in \([r,c+d-1]\), and \(e\)
in \([c+d,n-1]\). The weighted-subset formula and
\(|\gamma|=\sum_i\beta_i-\ell(\ell-1)/2\) give
\[
\sum_{\gamma:e(\gamma)=e}q^{|\gamma|}
=q^{ec}G(r,e)
\sum_{j\ge0}q^{j(R-L+j)}
G(R,L-j)G(L-j+\epsilon,j), \tag{10}
\]
where \(R=r-e\), \(L=\ell-e\).
In deriving (10), the three choices contribute
\(G(r,d)G(d+\epsilon,j)G(r-d,e)\);
their minimum-degree shift is \(ec+j(r-\ell+j)\).
Use \(G(r,d)G(r-d,e)=G(r,e)G(r-e,d)\).

For \(\epsilon=0\), the flag identity changes the inner product to
\(G(R,j)G(R-j,L-2j)\), giving \(U_{2R,L}\).
For \(\epsilon=1\), use
\[
G(L-j+1,j)=G(L-j,j-1)+q^jG(L-j,j).
\]
Shift \(j=s+1\) in the first term, apply the same flag identity,
and combine with the second using
\[
G(R+1-s,L-2s)=G(R-s,L-2s)
+q^{R+1-L+s}G(R-s,L-2s-1).
\]
This gives \(U_{2R+1,L}\). Thus (10) is exactly the right side of (1),
and summing over \(e\le t\) proves that the pivot count in (7)
equals the Betti polynomial in (3).

## 7. The integral squeeze closes the gap

The selected rows give
\[
Q_t(m)\subseteq P_t(m)\subseteq M_m.
\]
Section 5 makes \(Q_t(m)\) primitive. Section 6 and the geometric
upper bound (5) force
\[
\operatorname{rank}Q_t(m)
=\operatorname{rank}P_t(m)
=\operatorname{rank}H^{2m}(W_t;\mathbb Z).
\]
A primitive sublattice contains every integral point of its rational
span. Therefore \(P_t(m)=Q_t(m)\), proving every prefix primitive
with the exact ranks (3). The adapted-basis Smith lemma in Section 2
now gives (1).

As a consequence, the ambient restriction
\(H^*(\operatorname{Gr}(\ell,n);\mathbb Z)\to H^*(W_t;\mathbb Z)\)
is surjective: its composite with the primitive incidence pullback
has a primitive image of the full Betti rank.
Here transposition preserves Smith invariants, so primitive coefficient
rows give a primitive image of the restriction matrix as well.
This is a consequence of the argument, not a hypothesis used in it.

## 8. Published structure, verification, and limitations

The integral B/C comparison, quaternionic bases, and Schur-product
reduction retain their [audited dependencies](PROOF_STATUS.md).
The polynomial \(U\) is classical Birkhoff/Hall counting.
There is a stronger published representation-theoretic interpretation:
Chari--Schneider--Shereen--Wand,
[*Modules with Demazure Flags and Character Formulae*](https://arxiv.org/pdf/1310.5191),
SIGMA 10 (2014), 032, Theorem 3.3(ii)--(iii), gives
\[
[W_{\rm loc}(n\omega):D(2,(n-2e)\omega)]_q
=q^{ec}G(r,e).
\]
The weight \(n-2\ell\) subspace in the remaining Demazure character has
polynomial \(U_{n-2e,\ell-e}\), by their alternating character identity.
Thus the character expression predates these notes.
That source works over \(\mathbb C\); it does not by itself determine
the integral Springer-image lattice. Sections 3--7 supply that comparison.

The [general verifier](experiments/verify_type_b_general_smith.py)
and [its report](experiments/type_b_general_smith_verification.json)
check the beta-number criterion, weighted counts, and actual integral
LR modules with the constructed triangular rows.
The [high-layer note](TYPE_B_HIGH_LAYERS.md) gives an independent
elementary argument for the three highest layers.
Finite checks support implementation and transcription; the proof
of the unrestricted statement is the integral squeeze above.

This note completes the type-B argument under the stated normalization
and earlier audited reduction. It does not settle the separate general
type-D partition classification or establish novelty of these deductions.
