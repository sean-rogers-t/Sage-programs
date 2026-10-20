# The general even-tail gluing theorem

**Status: internally audited complete repository deduction,
19 October 2026.** The [dependency audit](TYPE_D_EVEN_GLUING_AUDIT.md)
records the independent reviews and one corrected relative-bundle
wording issue. This is an unpublished repository deduction from
the published Schubert calculus inputs identified below. It makes no
originality claim and has not received external expert review.

## 1. Statement and normalization

Put \(n=a+k\), with \(a\ge1\) and \(k\ge2\), and let
\(X_D=\operatorname{OG}(a,2n)\). Use the fixed normalized Springer inputs
\[
q_i=e_i(x_{a+1}^2,\ldots,x_n^2),\qquad
p=x_{a+1}\cdots x_n,\qquad p^2=q_k.
\]
The square image is \(L_+=\mathbb Z[q_1,\ldots,q_k]\); the Euler
image is \(L_-=p\mathbb Z[q_1,\ldots,q_k]\), with their cohomological
kernels understood. All statements below are homogeneous. Polynomial
degree \(d\) means cohomological degree \(2d\).

Let \(S_\pm\) be the saturations of \(L_\pm\) in integral cohomology,
and let \(S\) be the saturation of their sum. The rational sectors are
independent. The integral gluing group is
\[
W=S/(S_+\oplus S_-).
\]
The [two-flag lattice lemma](ORTHOGONAL_SMITH.md#lattice-lemma)
embeds this elementary 2-group in each \(S_\pm/2S_\pm\).
Let \(W_0^\pm\) denote the inverse images in \(W\) of the reductions
of \(L_\pm\).

**Theorem.** For every even tail \(k=2h\),
\[
                         W=W_0^+\cap W_0^-.
\tag{1}
\]
Thus all gluing occurs between unit levels of the two sector lattices.
Its degree-\(2m\) dimension is the coefficient of \(q^m\) in
\[
I_{a,2h}(q)=q^{h(\lfloor a/2\rfloor+1)}
               {h+\lceil a/2\rceil-1\brack h}_q.
\tag{2}
\]
There is no additional gluing at positive Smith levels.

Here is the resulting complete Smith formula. Write \(G(v,u)={v\brack u}_q\),
zero outside \(0\le u\le v\), and
\[
U_{N,L}(q)=\sum_{j=0}^{\lfloor L/2\rfloor}
 q^{j(\lceil N/2\rceil-L+j)}
 G(\lfloor N/2\rfloor,j)G(\lceil N/2\rceil-j,L-2j).
\]
Set \(U_{0,0}=1\), and put
\[
B_{a,k,e}(q)=q^{e\lceil(a+k)/2\rceil}
 G(\lfloor(a+k)/2\rfloor,e)
 U_{a+k-2e,\min(a,k)-e}(q).
\tag{3}
\]
This is zero for \(e\notin[0,\min(a,k)]\); at \(a=0\) use
\(B_{0,k,0}=1\) and \(B_{0,k,e}=0\) for \(e>0\).
Set \(I_{a,k}=0\) when \(k\) is odd.

If \(M_D(a,k;e,d)\) counts the nonzero factor \(2^e\) in degree \(d\),
then for all \(a\ge1,k\ge2\),
\[
\boxed{\ \sum_d M_D(a,k;e,d)z^d
 =B_{a,k,e}(z^2)+z^kB_{a-1,k,e}(z^2)
   +(\delta_{e,1}-\delta_{e,0})I_{a,k}(z^2).\ }
\tag{4}
\]
Units are retained. These are factors in the saturation of the selected
rational Springer image, not an index in all cohomology. For the
literal defining/vector polynomial-character form, multiply each
degree-\(d\) factor by \(2^d\), as in [TYPE_D_SMITH.md](TYPE_D_SMITH.md).

## 2. Relocate the specified tail, keeping its source algebra

Put \(H=k\), \(K=a\), so \(n=H+K\). Pull back to the complete flag
and permute the \(k\) specified tail roots to the head. Complete-flag
pullback is a primitive integral injection, and a Weyl permutation is
an integral cohomology automorphism. Both operations preserve the
coefficient-row lattice of any fixed list of source polynomials.
The same list descends to
\[
X'_D=\operatorname{OG}(H,2n).
\]
In common integral source coordinates it is
\[
\mathcal F:\quad
s_\lambda(x_1^2,\ldots,x_H^2),\quad \lambda\subseteq K^H;
\qquad
e_Hs_\mu(x_1^2,\ldots,x_H^2),\quad
\mu\subseteq(K-1)^H,
\quad e_H=\prod_{i=1}^H x_i.
\tag{5}
\]
These rectangle kernels and independent rational sectors follow from
the square and Euler comparisons in
[TYPE_D_SMITH.md](TYPE_D_SMITH.md) and
[TYPE_D_SQUARE_COMPARISON.md](TYPE_D_SQUARE_COMPARISON.md).
In particular (5) is an integral basis of the original source modulo
its kernel. The narrower Euler rectangle is essential.

Let \(\mathcal R_D,\mathcal R_B\) be the coefficient-row lattices of
this same list on \(X'_D\) and
\(X'_B=\operatorname{OG}(H,2n+1)\). Integral restriction along
\(X'_D\hookrightarrow X'_B\) gives
\[
                         \mathcal R_D\subseteq\mathcal R_B.
\tag{6}
\]
We will prove equality. This comparison is specific to (5); it is
false for arbitrary mixtures of head Chern polynomials.

Moving the specified tail to the head is not a symmetry between two
full tail-image problems. In particular the original Euler generator
is kept fixed throughout this change of coordinates.

## 3. The type-B lattice and its prefixes

The complete integral calculation is in
[TYPE_D_HEAD_EULER_B_COMPARISON.md](TYPE_D_HEAD_EULER_B_COMPARISON.md).
For clarity, the inputs needed here are as follows.

Let \(P_t(H,K)\) denote the unweighted symplectic square prefix with
at most \(t\) negative head entries. Each prefix is primitive for the
actual head size \(H\), including \(H>K\). Section 2 of that note
extends the envelope and unit-triangular proof to this case; a block
exchange of weighted Smith lists alone would not establish it.
The point boundary \(K=0\) has only its degree-zero class.

For \(t>0\), the full unweighted prefix of (5) is exactly
\[
N_t=P_t(H,K)\oplus P_t(H,K-1),
\tag{7}
\]
with Euler degrees shifted by \(H\). On the symplectic incidence
bundle the two sources occupy the distinct integral Lagrangian
Leray--Hirsch components \(1,c_t(U^*)\).

Let \(M_0\) be the actual combined level-zero row lattice. It is the
row lattice of (5) on the ordinary Grassmannian
\(\operatorname{Gr}(H,n)\). Exact B/C row rescaling gives
\[
                  \mathcal R_B=M_0+\sum_{t=1}^H2^tN_t.
\tag{8}
\]
This equality is integral and does not require primitivity of \(M_0\).

For its Smith list, the level-zero sum is primitive over
\(\mathbb Z_{(2)}\) and has rank equal to the sum of its two separate
prefix ranks minus \(I_{K,H}\). The crucial integral overlap identities
are, with \(H=2h\) and \(v=\lfloor K/2\rfloor\),
\[
\psi_2(s_\theta)=(-1)^h e_H
                    \psi_2(s_{\theta-(1^h)})
 \quad\text{in }H^*(\operatorname{Gr}(H,n);\mathbb Z),
\quad
\theta\subseteq K^h,\quad\theta_h\ge v+1.
\tag{9}
\]
The signed SXP proof and the rational/mod-two rank squeeze are given
in Section 4 of the companion note. An adapted basis for the nested
primitive chain over \(\mathbb Z_{(2)}\) computes (8). Its list is the
union of \(B(H,K)\) and the shifted \(B(H,K-1)\), with exactly
\(I_{K,H}\) zero exponents replaced by ones. The full symplectic
source is primitive, so B/C rescaling also excludes odd-primary
factors. This proves (4) for the selected B source (5).

## 4. Relative orthogonal incidences and the integral Euler algebra

Fix \(t>0\), \(b=H-t\), \(N=K+t\). Over
\(\operatorname{Gr}(b,n)\), let \(S\) be the tautological bundle and
\(W=\mathbb C^n/S\). Retain its full ordinary flag bundle
\(F=\operatorname{Fl}(W)\), and put \(B=H^*(F;\mathbb Z)\).
The flagged orthogonal incidence is
\[
Y_t=\operatorname{OG}(t,W\oplus W^*)\longrightarrow F,
\qquad j_t:Y_t\longrightarrow X'_D.
\]
The orthogonal form is hyperbolic. Write \(U\) for the residual
isotropic \(t\)-plane and use the notation
\[
c=c_b(S^*),\quad e=c_t(U^*),\quad
p=e(U^\perp/U),\quad Z=\pm c_N(W).
\]
Choose compatible Euler signs, so
\[
               e^2=u_t,\qquad p^2=v_K,\qquad ep=Z,
\tag{10}
\]
where \(u_i,v_j\) are the squared head and tail characteristic classes.
In general \(Z\ne0\); it must not be discarded.

Let \(R\) be the relative quaternionic Grassmann ring over \(B\):
\[
R=B[u_1,\ldots,u_t,v_1,\ldots,v_K]/
 \left[\left(\sum u_i s^i\right)\left(\sum v_j s^j\right)
                       =\prod_{\ell=1}^N(1+w_\ell^2s)\right].
\]
The [relative Euler extension lemma](TYPE_D_RELATIVE_EULER_EXTENSION.md)
proves the following integral assertions:

1. \(\mathcal C=R[e,p]/(e^2-u_t,p^2-v_K,ep-Z)\) injects into
   \(H^*(Y_t;\mathbb Z)\). It is projective over \(B\) and torsion free.
2. \(\mathcal D=R\oplus ceR\) is primitive in \(\mathcal C\), and
   \(ceR\cong R/(Z,v_K)\).
3. The first summand is the quaternionic incidence ring for \((H,K)\).
   The second identifies with the corresponding ring over
   \(\operatorname{Gr}(b,n-1)\), with residual tail rank \(K-1\).

The proof does not take invariants of an integral quotient. It proves
polynomial freeness over the explicit type-D invariant generators,
uses direct summands under base change to rule out hidden torsion,
and computes the additive quotients and annihilators. In particular
\(\mathcal C/(R+eR)\cong R/(u_t,Z)\) and
\((R+eR)/(R+ceR)\cong R/(c)\) are free abelian. This is the
additional integral step required for the Euler source.

The restrictions of (5) lie in \(\mathcal D\), because \(e_H\)
restricts to \(ce\). Schur branching identifies their coefficient
rows in the two summands with the two primitive prefixes in (7).
For the second summand, restriction to the quotient by \((Z,v_K)\)
maps the submodule with quaternionic Schur indices of width at most
\(K-1\) isomorphically, basis by basis, onto the target.
Every source Euler term lies in this submodule: its partition has
width at most \(K-1\), and every branching subpartition has that bound.
No assertion that wider relative Schur classes vanish is needed.
Thus these coefficient rows agree integrally with
\(P_t(H,K-1)\), not merely rationally.

Consequently the source image \(L_t\) of (5) is primitive in
\(\mathcal D\), hence in \(\mathcal C\), in its rational span; its
row lattice calculated there is \(N_t\).

## 5. Denominators and a parity retraction give the reverse containment

Let \(A_t\) be the subring of \(H^*(Y_t;\mathbb Z)\) generated by
base classes, the Chern classes of \(U^*\), and \(p\).
The relative eta denominator lemma in
[TYPE_D_SQUARE_THEOREM.md, Section 4](TYPE_D_SQUARE_THEOREM.md#4-the-relative-eta-denominator-bound)
gives
\[
                         2^tH^*(Y_t;\mathbb Z)\subseteq A_t.
\tag{11}
\]
Its published input is Tamvakis's
[double eta Giambelli theorem](https://math.umd.edu/~harryt/papers/dbleta.pdf),
Definition 1, Theorem 1, and Section 4.1. A finite fiber partition has
length at most \(t\); its denominator is at most \(2^t\), including
the distinguished typed factor. Retaining \(F\) supplies the torus
reduction and integral equivariant Schubert Leray--Hirsch classes.

There is an integral additive retraction
\[
                            \rho:A_t\longrightarrow\mathcal C.
\tag{12}
\]
Reduce any representative to \(F_0(z)+pF_1(z)\). In both symmetric
polynomials retain precisely monomials whose \(t\) head exponents
are all even or all odd. They become square polynomials and \(e\)
times square polynomials; use \(ep=Z\) for the last product.
This coefficient deletion is integral. Rationally it is averaging
over the even head sign changes. Those Weyl changes preserve the
relative squared characteristic relations and the global Euler
relation, with base fixed. Thus the operation preserves the rational
characteristic kernel. Torsion-free cohomology and injective
complete-flag pullback make it independent of the representative
integrally. It fixes \(\mathcal C\). For \(t=1\) the sign group is
trivial and the same statement holds.

Let \(h\) be any integral cohomology class in the rational span of
\(L_t\). By (11), \(2^th\in A_t\); by (12), which fixes that span,
\(2^th\in\mathcal C\). Primitivity of \(L_t\) in \(\mathcal C\)
then gives \(2^th\in L_t\).
Every nonzero Smith factor of the source map on \(Y_t\) therefore
divides \(2^t\).

Let \(T_t\) be its coefficient-row lattice in integral cohomology.
Its rational kernel agrees with the map into \(\mathcal C\), whose
primitive row lattice is \(N_t\). Hence
\(\operatorname{Sat}(T_t)=N_t\). Transposing the Smith form gives
\[
                         2^tN_t\subseteq T_t\subseteq\mathcal R_D.
\tag{13}
\]
The second inclusion is integral restriction along \(j_t\).
At \(t=0\) the incidence is the ordinary Grassmannian, and its rows
give \(M_0\subseteq\mathcal R_D\) directly. Summing (13) and using
(8) proves \(\mathcal R_B\subseteq\mathcal R_D\). Together with (6),
\[
                              \mathcal R_D=\mathcal R_B.
\tag{14}
\]
This proves the full formula (4) for even tails. For odd tails the
two original sectors occupy opposite degree parities; the earlier
square and Euler theorems already prove (4).

When \(K=1\), the relocated D space omits both fork nodes. The
incidence proof still applies, with \(K-1=0\) interpreted as the
point quaternionic Grassmannian. It makes no assertion that this
boundary is a single-node maximal parabolic.

## 6. The full formula forces exhaustion by unit gluing

The independently proved square-sector list is \(B(a,k)\); the Euler
list is the shifted \(B(a-1,k)\). In a fixed degree, let \(v_+\),
\(v_-\), and \(v\) be the sums of their base-two Smith exponents
and of the full list, respectively. The exact index identity is
\[
[S:L_+\oplus L_-]
   =[S:S_+\oplus S_-]\,[S_+:L_+]\,[S_-:L_-].
\]
Consequently
\[
                         \dim_{\mathbb F_2}W=v-v_+-v_-.
\tag{15}
\]
Formula (4) makes the right hand side exactly the coefficient \(I_m\)
of (2). The separate [unit theorem](TYPE_D_FROBENIUS.md) already
proves \(\dim(W_0^+\cap W_0^-)=I_m\). This is a subspace of \(W\);
their equal dimensions prove (1).

Thus the proof identifies all saturated gluing, rather than inferring
it from the original mod-two overlap alone.

## 7. Relation to the papers and independent checks

The Springer/Borel map and its polynomial generators come from the
papers already documented in [PROOF_STATUS.md](PROOF_STATUS.md).
The additional assertion here describes the integral divisibility of
that fixed image inside its rational saturation. It does not change
the rational Springer formulas or establish a new formula for the
Springer morphism itself.

Published ingredients include B/C Schubert rescaling, ordinary and
quaternionic Grassmannian characteristic presentations, and typed eta
Giambelli formulas. The general B prefix theorem, relative integral
comparison, and the full Smith multiplicity interpretation are
repository deductions. Their novelty in the literature is not claimed.

The independent two-head and two-tail proofs remain in
[TYPE_D_TWO_HEAD_GLUING.md](TYPE_D_TWO_HEAD_GLUING.md) and
[TYPE_D_TWO_TAIL_GLUING.md](TYPE_D_TWO_TAIL_GLUING.md).
The [full exponent bound](TYPE_D_FULL_EXPONENT_BOUND.md) is another
independent consequence of the parity method.
The stronger prefix-intersection formulation in
[TYPE_D_GLUING_FILTRATION.md](TYPE_D_GLUING_FILTRATION.md) is a
separate conditional reduction; (14) does not silently promote every
one of its intermediate assertions to a theorem.

Exact finite checks are recorded in the
[full-gluing report](experiments/type_d_even_gluing_verification.json),
[extended head-Chern comparison](experiments/type_d_extended_chern_comparison_verification.json),
[actual-head prefix report](experiments/type_d_head_euler_b_verification.json),
and [higher-prefix report](experiments/type_d_higher_gluing_filtration_verification.json). They corroborate the argument and include the
two formerly timed-out rank-nine components. Finite calculations
are not used as a substitute for (10)--(14).
