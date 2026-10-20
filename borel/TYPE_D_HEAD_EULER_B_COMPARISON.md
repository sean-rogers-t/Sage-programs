# The relocated type-B head-Euler image

**Status: repository deduction, 19 October 2026.**
This note computes the type-B lattice that occurs after moving a
specified even type-D tail to the head. The signed overlap identity
and the use of coefficients localized at two are explicit below.
The separate comparison from this lattice back to type D requires its
relative incidence argument; it is not inferred from this calculation.

## 1. The selected algebra

Let \(H=2h\ge2\), \(K\ge1\), and \(n=H+K\). On
\(\operatorname{OG}(H,2n+1)\), use head roots \(x_1,\ldots,x_H\)
and put \(e_H=x_1\cdots x_H\). Select the source classes
\[
s_\lambda(x^2),\quad \lambda\subseteq K^H,
\qquad
e_Hs_\mu(x^2),\quad \mu\subseteq(K-1)^H.
\tag{1}
\]
Their degrees are \(2|\lambda|\) and \(H+2|\mu|\).
This is a specified sublattice of the head Euler image: the narrower
rectangle on the second family matters.

For comparison with the original full D problem, take
\(H=k\), \(K=a\). A permutation of complete-flag roots sends the
original specified tail Euler class to this head Euler class. The
square and Euler source rectangles in (1) then have the required
rational ranks. This does not exchange the two full tail-image problems.

Let \(P_t(H,K)\) denote the actual-head symplectic square prefix:
the coefficient rows of \(s_\lambda(x^2)\) with at most \(t\)
negative entries. Prefixes for the point source \(K=0\) mean the
rank-one degree-zero lattice and zero lattices otherwise.

## 2. Primitivity for the actual head

The integral triangular squeeze in
[TYPE_B_GENERAL_SMITH.md](TYPE_B_GENERAL_SMITH.md) extends to arbitrary
\(H,K\), without first moving the smaller block to the head. Here
are the changes needed when \(H>K\).

Write
\[
r=\lfloor n/2\rfloor,\quad c=\lceil n/2\rceil,\quad
\epsilon=c-r,\quad b=H-t.
\]
The balanced Chow envelope remains defined: its dimensions
\(s=\lfloor b/2\rfloor\) and \(p=\lceil b/2\rceil\)
satisfy \(s\le r\), \(p\le c\). For an invariant \(b\)-plane
\(S\), the same argument gives
\(\dim(S\cap I)\ge\max(j,b-j-\epsilon)\ge s\), where
\(j=\operatorname{rank}N|_S\). Thus its integral Chow-envelope
argument works for every actual head.

The unbalanced incidence with factors
\[
s_\eta s_\alpha s_\beta,\qquad
\eta\subseteq K^t,\quad
\alpha\subseteq(c-j)^j,\quad
\beta\subseteq(r-b+j)^{b-j}
\tag{2}
\]
is used only for feasible indices
\[
\max(0,b-r)\le j\le\min\left(b,c,
                         \left\lfloor\frac{b+\epsilon}{2}\right\rfloor\right).
\tag{3}
\]
The envelope projection formula still puts every such row in
\(P_t(H,K)\). A coordinate partition \(\gamma\subseteq K^H\)
is selected when its last \(b\) parts fit
\[
                   ((c-j)^j,(r-b+j)^{b-j})
\tag{4}
\]
for an index (3). When writing this as a single corner inside the source
rectangle, cap the entries of (4) by \(K\). This changes no test on
\(\gamma\).

The unit triangular minor is unchanged: choose the first \(t\)
parts of \(\gamma\) for \(\eta\), the next \(j\) for
\(\alpha\), and the rest for \(\beta\), and use
\(h_\eta s_\alpha s_\beta\). Its leading complete coordinate has
coefficient one and all other coordinates have greater energy.

The beta-number criterion also remains unchanged. Put
\(\beta_i=\gamma_i+H-i\),
\(d=\#\{\beta_i<r\}\), \(J=H-d\), and
\[
                      e(\gamma)=\#\{\beta_i\ge c+d\}.
\]
Then (4) holds for a feasible \(j\) exactly when
\(e(\gamma)\le t\). To check the only new feasibility issues:
if \(J\le t\), then \(b\le d\le r\), so \(j=0\) is
feasible. If \(J>t\), the usual choice
\(j_0=J-t=b-d\) satisfies \(j_0\ge b-r\) and
\(j_0\le c\). The same middle-bead interval gives
\(2j_0\le b+\epsilon\). Hence \(j_0\) satisfies (3).

The weighted beta-number count and the defect-stratum upper bound use
the same formulas as the general B proof. They remain valid for actual
head \(H\): an out-of-range Grassmann choice contributes zero. The
invariant-subspace polynomial \(U_{N,L}\) is symmetric under
\(L\leftrightarrow N-L\), by annihilator duality with the transpose
square-zero operator. Thus for \(H>K\),
\(U_{n-2e,H-e}=U_{n-2e,K-e}\), giving the same known B counts.
The primitive triangular sublattice and this matching rank upper bound
again squeeze the entire prefix. Consequently every actual-head
\(P_t(H,K)\) is primitive. This conclusion concerns prefix modules;
it is not obtained by assuming that block exchange preserves the
individual signed Schubert rows.

## 3. Every positive prefix splits into the two sectors

Use the symplectic incidence resolution for the \(\nu\le t\)
Schubert union. For \(b=H-t\), it is a bundle over
\(\operatorname{Gr}(b,n)\). Forgetting the residual complex
Lagrangian gives its quaternionic Grassmann bundle \(Q_t\), with
fiber of the forgetful map \(\operatorname{LG}(t,2t)\).
Both pullbacks used here are primitive integral injections, as in
[TYPE_B_FILTRATION.md](TYPE_B_FILTRATION.md).

The square source restricts to a class pulled back from \(Q_t\).
The Euler source restricts to
\[
                    c_b(S^*)\,c_t(U^*)\,s_\mu(S^{*2}+U^2).
\tag{5}
\]
For \(t>0\), the classes \(1\) and \(c_t(U^*)\) belong to an
integral Leray--Hirsch basis for the Lagrangian fiber: the second is
its primitive special Schubert class. The remaining integral Schubert
classes are also universal Chern polynomials, so this basis is global.
Thus the two families in (1) occupy independent integral coefficient
components, and their prefix row lattice is a direct sum.

For the Euler component, multiplying by \(c_b(S^*)\) shifts each
ordinary base Schur partition by the full \(b\)-column. Equivalently,
before that shift the width bound is one smaller. Its SXP two-quotient
therefore uses \(n-1\), not \(n\). Moreover every residual
\(s_\eta(U^2)\) has width at most \(K-1\), because its source
partition \(\mu\) does. The resulting coefficient products are
exactly the prefix products for \((H,K-1)\). Consequently the full
unweighted prefix for every \(t>0\) is
\[
                  P_t(H,K)\oplus P_t(H,K-1),
\tag{6}
\]
with the Euler source degree shifted by \(H\).
The smaller width in (1) is what makes (6) hold; it would be false
for the entire wider head-Euler source.

## 4. The level-zero signed overlap

At \(t=0\), the resolution is the ordinary Grassmannian
\(G=\operatorname{Gr}(H,n)\). Its two source images are
\[
\psi_2(s_\lambda),\qquad e_H\psi_2(s_\mu),
\tag{7}
\]
in its integral Schur ring, where \(\psi_2 f(x)=f(x^2)\).
Each image separately is primitive. For the second, multiplication
by \(e_H\) is the primitive full-column injection from
\(\operatorname{Gr}(H,n-1)\); use the primitive square prefix there.

Put \(u=\lceil K/2\rceil\), \(v=\lfloor K/2\rfloor\).
For every \(\theta\subseteq K^h\) with
\(\theta_h\ge v+1\), there is the exact integral identity
\[
          \psi_2(s_\theta)
             =(-1)^h e_H\psi_2(s_{\theta-(1^h)})
                       \quad\text{in }H^*(G;\mathbb Z).
\tag{8}
\]
This strengthens the corresponding identity modulo two used in
[TYPE_D_FROBENIUS.md](TYPE_D_FROBENIUS.md).

To verify the signs, the SXP coefficient on a target \(\gamma\)
with empty two-core is
\(\varepsilon_2(\gamma)c_{\alpha\beta}^{\theta}\).
The truncation bounds are \(\alpha\subseteq u^h\) and
\(\beta\subseteq v^h\). If \(\gamma\) has no full \(H\)-column,
then \(\alpha_h=0\); the LR bound
\(\theta_h\le\alpha_h+\beta_1\le v\) rules it out.
Otherwise subtraction of that column changes the two-quotient to
\((\beta,\alpha-(1^h))\), and determinant twisting gives
\[
c_{\alpha\beta}^{\theta}
       =c_{\beta,\alpha-(1^h)}^{\theta-(1^h)}.
\]
The two-sign of an empty-core partition is
\((-1)^{\gamma_2+\gamma_4+\cdots+\gamma_H}\): horizontal dominoes
contribute an even number of boxes to even rows and vertical dominoes
contribute one. Subtracting the full column therefore changes that
sign by \((-1)^h\), proving (8). Empty-core support is preserved by
this full even column shift; all other coefficients vanish on both
sides.

These intersection classes have an identity minor up to row signs:
choose the row with two-quotient
\(\alpha=\theta-(v^h)\), \(\beta=(v^h)\). The determinant twist
then singles out exactly \(\theta\). They are independent over
\(\mathbb Z\), and their number in half-degree \(m\) is
\[
I_m=[q^m]\;q^{h(v+1)}{h+u-1\brack h}_q.
\tag{9}
\]

The proved unit theorem says that the intersection of the two
reductions of (7) modulo two has dimension exactly \(I_m\).
Their separate images are primitive, so their rational and mod-two
ranks agree. The independent rational intersection (8) has dimension
at least \(I_m\), while matrix rank over \(\mathbb Q\) is at least
rank over \(\mathbb F_2\), giving the reverse bound. Thus the
rational intersection also has dimension exactly \(I_m\), and the
combined matrix in (7) has identical rational and mod-two ranks.

It follows that this combined image is primitive **over
\(\mathbb Z_{(2)}\)** in its rational span. No assertion excluding
odd-primary torsion in this particular prefix sum is needed here.

## 5. The selected B Smith list

Let \(N_t\) denote (6), including the separately defined direct sum
\(N_0=P_0(H,K)\oplus P_0(H,K-1)\). Let \(M_0\) be the actual
combined level-zero row lattice. By Section 4 it is primitive over
\(\mathbb Z_{(2)}\), of rank \(\operatorname{rank}N_0-I_m\).
For \(t>0\) the actual unweighted prefixes are \(N_t\), which are
primitive by Sections 2--3. These form the nested chain
\[
                       M_0\subset N_1\subset\cdots\subset N_H.
\]
Its weighted row lattice is
\(M_0+\sum_{t\ge1}2^tN_t\), by the exact B/C row rescaling.
An adapted basis over the PID \(\mathbb Z_{(2)}\) therefore gives
the union of the two sector Smith lists, with exactly \(I_m\)
zero exponents replaced by ones.

The two sector lists are those of \(B(H,K)\) and the degree-shifted
\(B(H,K-1)\), equivalently \(B(K,H)\) and the shifted
\(B(K-1,H)\). The full selected B lattice has no odd-primary
factors: after inverting two, B/C row rescaling is invertible, while
the full C source is primitive. The latter follows from the integral
Lagrangian Leray--Hirsch components \(1,e_H\), with the narrower
Euler rectangle a primitive Schur subset of the quaternionic base.
Thus the two-local calculation determines every integer Smith factor
of the selected B image (1).

This proves the desired full multiplicity expression for the relocated
B source. Transferring that answer to the original D source requires
the separate integral D/B comparison for (1).

## 6. Independent finite checks

The Schur/LR verifier
[verify_type_d_head_euler_b.py](experiments/verify_type_d_head_euler_b.py)
and its
[report](experiments/type_d_head_euler_b_verification.json)
check 610 actual-head prefix modules in 106 degree components with
\(H>K\), including their unit triangular minors and exact integer
balanced-module equality. They also check 92 level-zero components
by ordinary symmetric-function plethysm, including 35 exact signed
identities (8), their rational and mod-two overlaps, and two-local
primitivity. The verifier imports no Borel or Weyl-group calculator.
These finite checks corroborate the written arguments; they are not
their all-rank justification.
