# Three highest type-B Smith layers in every rank

**Status: repository deductions, 19 October 2026.** These arguments extend
the earlier audited notes using their integral Schur-product reduction.
They are not published Smith theorems or claims of originality. The
[general formula](TYPE_B_GENERAL_SMITH.md) is proved separately; this
note gives an independent elementary argument for the highest layers.
See [PROOF_STATUS.md](PROOF_STATUS.md) for the published inputs.

Put \(A=\max(a,k)\), \(\ell=\min(a,k)\), \(n=A+\ell\),
\(r=\lfloor n/2\rfloor\), \(c=\lceil n/2\rceil\), and
\(\epsilon=c-r\in\{0,1\}\). The variable \(q\) records half polynomial
degree, equivalently one quarter of the original cohomological degree.
Write \(G(v,u)={v\brack u}_q\), \([v]=1+\cdots+q^{v-1}\), and
\[
U_{n,b}=\sum_{j=0}^{\lfloor b/2\rfloor}
q^{j(\lceil n/2\rceil-b+j)}
G(\lfloor n/2\rfloor,j)G(\lceil n/2\rceil-j,b-2j).
\]
The parameters in \(U\) are recomputed at its arguments.

## Conclusions and exact scope

Let \(E_e=\sum_m M_e(m)q^m\), where \(M_e(m)\) counts the nonzero
Smith factors \(2^e\) in polynomial degree \(2m\) of the normalized image.
The following identities hold without assumptions on lower prefixes:
\[
\begin{aligned}
E_\ell&=q^{\ell c}G(r,\ell),\\
E_{\ell-1}&=q^{(\ell-1)c}G(r,\ell-1)[c-\ell+1]
&&(\ell\ge2),\\
E_{\ell-2}&=q^{(\ell-2)c}G(r,\ell-2)U_{n-2\ell+4,2}
&&(\ell\ge3).
\end{aligned} \tag{1}
\]
The first line is the [earlier largest-factor deduction](TYPE_B_HIGHEST_SMITH.md).
The next two lines are proved here. They agree with the general candidate.
The [unit deduction](TYPE_B_UNIT_SMITH.md) supplies \(E_0=U_{n,\ell}\).
Consequently **all factors are determined when \(\ell\le4\)**.
For \(\ell=4\), the only remaining factor count after (1) and \(E_0\)
is \(E_1\); the total-rank identity
\(\sum_e E^{\rm geom}_e=G(n,\ell)\) identifies it with the candidate.
The layers \(1\le e\le\ell-3\) for \(\ell\ge5\) are outside this
elementary argument and are covered by the separate general proof.

## 1. An integral change to complete-symmetric spectators

Use the [Schur-product prefix theorem](TYPE_B_FILTRATION.md):
\[
P_t=\operatorname{span}_{\mathbb Z}
\{s_\mu s_\alpha s_\beta\}\ \subset R_{A,\ell},
\quad
R_{A,\ell}=H^{2*}(\operatorname{Gr}(\ell,n);\mathbb Z),
\]
where \(\mu\subseteq A^t\), and, for \(b=\ell-t\),
\[
\alpha\subseteq(c-\lceil b/2\rceil)^{\lceil b/2\rceil},
\qquad
\beta\subseteq(r-\lfloor b/2\rfloor)^{\lfloor b/2\rfloor}.
\]
All spans are taken degreewise.

In this ring \(h_j=0\) for every \(j>A\). The products
\(h_\gamma=\prod_i h_{\gamma_i}\), \(\gamma\subseteq A^\ell\),
are an integral basis: their Schur expansions are Kostka unitriangular;
terms outside the rectangle vanish. The same argument restricted to
partitions in \(A^t\) shows
\[
\operatorname{span}_{\mathbb Z}\{s_\mu:\mu\subseteq A^t\}
=\operatorname{span}_{\mathbb Z}\{h_\mu:\mu\subseteq A^t\}. \tag{2}
\]
Indeed Jacobi--Trudi preserves length at most \(t\), and a term with
an index greater than \(A\) vanishes. Restricting the triangular transition
matrix therefore preserves its unit diagonal. Thus the unrestricted
factor \(s_\mu\) can be replaced by \(t\) complete-symmetric spectator
indices in \([0,A]\), padded by zeros.

## 2. A high-layer splitting lemma

Let \(M\) be the source coefficient lattice in one degree, and let \(L\)
be the weighted B row lattice. The full C row lattice is \(M\), so
\(2^\ell M\subseteq L\). For any primitive prefix \(P_t\),
\[
2^tP_t\subseteq L,\qquad L\subseteq P_t+2^{t+1}M.
\]
Choose an integral complement \(M=P_t\oplus Q\). The \(P_t\)-component
of a row of weight greater than \(t\) lies in
\(2^{t+1}P_t\subseteq L\) and may be subtracted. Hence \(L\) splits
between these two summands. Its factors in \(P_t\) divide \(2^t\);
its factors in \(Q\) are divisible by \(2^{t+1}\) and have full rank.
Therefore
\[
\#\{\text{exponents greater than }t\}
=\operatorname{rank}M-\operatorname{rank}P_t. \tag{3}
\]
This needs primitivity of this prefix alone.

## 3. The next-largest factor: a coordinate lattice

For \(b=2\), equation (2) makes the rows
\(h_\mu h_u h_v\), with \(\mu\subseteq A^{\ell-2}\),
\(u\le c-1\), \(v\le r-1\). Sorting the indices gives the exact equality
\[
P_{\ell-2}=
\bigoplus_{\gamma\subseteq(A^{\ell-2},c-1,r-1)}
\mathbb Z h_\gamma. \tag{4}
\]
Conversely every displayed coordinate occurs by taking the largest
\(\ell-2\) parts as \(\mu\). This proves primitivity in every rank.

If \(r=c\), the complement consists of partitions with
\(\gamma_{\ell-1}\ge c\). Those with \(\gamma_\ell<c\) contribute
\(q^{(\ell-1)c}G(r-1,\ell-1)[c]\); those with
\(\gamma_\ell\ge c\) contribute \(E_\ell\).
The flag identity
\([r]G(r-1,\ell-1)=G(r,\ell-1)[r-\ell+1]\)
gives the second line of (1).

If \(r=c-1\), split the complement into
\(\gamma_{\ell-1}\ge c,\ \gamma_\ell\le c-2\), and
\(\gamma_\ell\ge c-1\). Its polynomial is
\[
q^{(\ell-1)c}G(r-1,\ell-1)[r]
+q^{\ell(c-1)}G(r+1,\ell).
\]
Subtract \(E_\ell\), use
\(G(r+1,\ell)=G(r,\ell-1)+q^\ell G(r,\ell)\),
and then \([c-\ell]+q^{c-\ell}=[c-\ell+1]\).
This gives the same second line of (1), by (3).

## 4. The third-largest factor: a graph lattice

For \(b=3\), put \(H=c-1\). The rows become
\[
h_\mu(h_uh_v-h_{u+1}h_{v-1})h_w,
\quad \mu\subseteq A^{\ell-3},\
0\le v\le u\le H-1,\quad 0\le w\le r-1. \tag{5}
\]
A \(v=0\) row is a single coordinate, called an anchor.
Every other row is a difference of two complete-symmetric coordinates;
both terms survive since \(u+1\le H\le A\).
On a connected component these differences span the sum-zero lattice;
if the component has an anchor, they span its full coordinate lattice.
A spanning tree proves both assertions integrally. Thus
**\(P_{\ell-3}\) is primitive in every rank**.

For its rank, call a part exceptional if it exceeds \(H\).
Its sorted multiset is fixed along every edge and has size
\(j\le\ell-3\). Fix that multiset and write \(s=\ell-j\ge3\).
The remaining small parts have the following bounds:

| Parity of \(n\) | Small coordinates | Highest anchored small degree |
| --- | --- | --- |
| Even | \(H^s\) | \((s-1)H-1\) |
| Odd | \((H^{s-2},(H-1)^2)\) | \((s-1)H-2\) |

Here the odd bound records that \(w\le H-1\), and at most one \(H\)
can come from the active Jacobi--Trudi pair.

Every fixed-degree small sector is connected. If largest and smallest
small entries \(a,b\) have gap at least two, replace them by
\((a-1,b+1)\). This is the reverse edge in (5), with
\(u=a-1,\ v=b+1\). In the odd case, when \(a=H\), at least two
small entries are at most \(H-1\); after selecting the minimum, another
can serve as \(w\). All remaining entries are spectators.
The sum of squares decreases until the unique balanced partition is
reached. It stays within the displayed shape: the maximal small sum
is \(sH\) or \(sH-2\), respectively.

An anchor has one zero, one entry at most \(H-1\), one at most
\(r-1\), and \(s-3\) other entries at most \(H\).
Sums of these integer intervals attain every degree up to the bound
in the table. Above it, a nonempty connected sector loses exactly one
rank. Coordinates absent from every row contribute zero, as required.

The polynomial for exceptional multisets of size \(j\) is
\[
Z_j=q^{jc}G(r-\ell+j,j).
\]
The rank polynomial of the small sector is
\[
K_s=
\begin{cases}
G(c+s-1,s)-q^{(s-1)(c-1)}[c],&\epsilon=0,\\
G(c+s-1,s)-q^{(s-1)(c-1)-1}([c-1]+q[c]),&\epsilon=1.
\end{cases} \tag{6}
\]
For odd \(n\), subtract from \(H^s\) the coordinates with \(s-1\)
or \(s\) entries equal to \(H\), then subtract the unanchored degrees.
This gives exactly the second line. Therefore the closed rank is
\[
F_{\ell-3}=\sum_{j=0}^{\ell-3}Z_jK_{\ell-j}. \tag{7}
\]

## 5. Gaussian simplification, rather than numerical fitting

Partitioning the entire rectangle by its exceptional multiset gives
\[
G(n,\ell)=\sum_{j=0}^{\ell}
Z_jG(c+\ell-j-1,\ell-j).
\]
Also \(\sum_{j=0}^L q^jG(d+j,j)=G(d+L+1,L)\), by Gaussian Pascal.
Consequently the rank deficit \(D=G(n,\ell)-F_{\ell-3}\) is
\[
\begin{aligned}
D={}&q^{(\ell-2)c}G(r-2,\ell-2)G(c+1,2)\\
&+q^{(\ell-1)c}G(r-1,\ell-1)[c]+E_\ell\\
&+q^{(\ell-1)(c-1)-\epsilon}
  \begin{cases}[c],&\epsilon=0,\\[c-1]+q[c],&\epsilon=1\end{cases}
  G(r-2,\ell-3).
\end{aligned} \tag{8}
\]
By (3), \(D=E_{\ell-2}+E_{\ell-1}+E_\ell\).
For completeness the simplification to (1) uses the following identities.
Put \(k=\ell-2\) and \(h=c-k\).

For even \(n\), the second term of (8) equals \(E_{\ell-1}\).
The remaining identity is
\[
G(c-2,k)G(c+1,2)+q^{h-1}[c]G(c-2,k-1)
=G(c,k)G(h+1,2). \tag{9}
\]
After cancelling Gaussian factorials, it is
\([c-1][h+1]=[c+1][h-1]+q^{h-1}[k][2]\).
Multiply by \((1-q)^2\) to check this directly.
Since \(U_{2h,2}=G(h,2)+q^{h-1}[h]=G(h+1,2)\), this proves (1).

For odd \(n\), first use
\[
[c]G(c-2,\ell-1)-[c-\ell+1]G(c-1,\ell-1)
=-q^{c-\ell}G(c-2,\ell-2).
\]
After subtracting the two highest layers from (8) and dividing by
\(q^{kc}\), the result is
\[
G(c-3,k)G(c+1,2)
+q^{h-2}([c-1]+q[c])G(c-3,k-1)
-q^{2c-k-2}G(c-2,k).
\]
Apply (9) with \(c-1\) in place of \(c\), then
\(G(c+1,2)-G(c,2)=q^{c-1}[c]\) and
\(G(c-2,k)=G(c-3,k-1)+q^kG(c-3,k)\).
It becomes
\[
G(c-1,k)G(h,2)+q^{h-1}[c-1]G(c-2,k)
=G(c-1,k)\bigl(G(h,2)+q^{h-1}[h-1]\bigr).
\]
The bracket is \(U_{2h-1,2}\). This proves the third line of (1)
for both parities, without primitivity assumptions on any smaller prefix.

## Verification and next question

The [high-layer verifier](experiments/verify_type_b_high_layers.py)
checks exact graph components, closed polynomials, and independent
integer LR row modules. Its [report](experiments/type_b_high_layers_verification.json)
records the finite scopes and source hashes. These checks corroborate
the proof and its transcription; they do not establish arbitrary-rank
statements by themselves.

For \(b=4\), both remaining Schur factors can have two rows. A row is
then a four-term product of two Jacobi--Trudi differences.
The graph argument above does not prove primitivity in that case.
The [general proof](TYPE_B_GENERAL_SMITH.md) instead uses unbalanced
incidence bundles, a Chow envelope, and a unit triangular minor to
handle every prefix.
