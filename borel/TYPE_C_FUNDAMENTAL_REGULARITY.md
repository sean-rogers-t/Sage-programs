# Regularity obstruction for every type C fundamental Springer map

Date: 22 October 2026. Status: a written arbitrary-rank deduction with
independent exact finite corroboration. The representation construction,
character identities and trace formula are published inputs. Novelty and
external review of this regularity classification remain unverified.
No PDF or LaTeX article is changed.

## 1. Statement, group and normalization

Let G=Sp_(2n)(C) and V_(omega_k) its fundamental representation,
1<=k<=n. For k<=n this is the kernel of symplectic contraction on
exterior^k(C^(2n)); it is not the full exterior power when k>=2.
All these nontrivial representations have faithful differential; a finite
central kernel is allowed. Let theta_k:G -> sp_(2n) be the trace Springer
map, normalized by

\[
 \operatorname{tr}(d\rho_k(\theta_k(g))d\rho_k(X))
 =\operatorname{tr}(\rho_k(g)d\rho_k(X)).
\]

Write

\[
 t(z)=\operatorname{diag}(z_1,\ldots,z_n,z_1^{-1},\ldots,z_n^{-1}),
 \quad
 \theta_k(t(z))=\operatorname{diag}(y_1,\ldots,y_n,-y_1,\ldots,-y_n).
\]

Regular means centralizer dimension n. Thus torus regularity is
z_i!=1,-1 and z_i!=z_j,z_j^(-1) for i!=j. Image regularity is
y_i!=0 and y_i!=y_j,-y_j for i!=j. Equivalently, with
x_i=z_i+z_i^(-1), input regularity is x_i!=2,-2 and x_i!=x_j.

**Theorem.** For every n>=2 and every k=1,...,n, theta_k fails to
preserve torus regularity. In particular it fails global regularity
preservation. Each failure has an explicit algebraic regular semisimple
torus witness. At n=1 the sole fundamental representation is the defining
SL2 representation and does preserve global regularity.

The theorem classifies all fundamental highest weights in type C.
It does not classify other irreducible highest weights or reducible
representations. There are no remaining nonsemisimple positive cases
to check at n>=2: every case already fails on T_reg.

## 2. Established inputs and coordinate factorization

The primitive exterior-power construction and character identity are
classical; see [Fulton--Harris, Representation Theory, Theorem 17.5](https://www.cimat.mx/~gil/docencia/2008/representaciones/fulton-harris.pdf),
printed page 260 (page 138 in the scanned two-page PDF, visually checked).
Contraction is surjective in this range by the usual exterior Lefschetz
identity: if L wedges with sum_i e_i wedge f_i, then C=L^* and
[C,L]=(n-d)id on exterior degree d. For k>=2 put d=k-2<n. Then L is injective,
since ||Lv||^2=||Cv||^2+(n-d)||v||^2; hence C onto degree d is
surjective. This explains the character subtraction below. These are
standard linear-algebra/representation inputs.
The normalized trace construction is from
[Kostant--Michor, Section 2.1](https://arxiv.org/pdf/math/0109066), and
the general explicit weight formula is from
[Rogers, Theorem 1](https://arxiv.org/pdf/1701.01538v1).
These inputs are not new results of this note.

Let E_z(u)=product_i(1+z_i u)(1+z_i^(-1)u). With e_j the elementary
symmetric functions in the 2n-letter alphabet, put e_j=0 for j<0.
Then

\[
 \chi_k^{(n)}=e_k-e_{k-2}
 =[u^k]\,(1-u^2)\prod_{i=1}^n(1+x_i u+u^2).
\]

Use chi_0^(r)=1, including rank zero. Relative to the defining trace
form on sp_(2n), the index of V_(omega_k) is

\[
 I_{n,k}=\binom{2n-2}{k-1}-\binom{2n-2}{k-3}>0.
\]

The out-of-range binomials are zero. For H_i with diagonal 1 at i,
-1 at n+i and zero elsewhere, the trace Gram matrix is 2I_(n,k) times
the identity. Counting k-subsets containing exactly one of this pair
gives this formula for exterior powers, and subtracting the k-2 trace
gives the primitive representation index. Hence

\[
 \boxed{\quad
 y_i=\frac{z_i-z_i^{-1}}{2I_{n,k}}\,
 \chi_{k-1}^{(n-1)}(z_1,\ldots,\widehat z_i,\ldots,z_n).
 \quad}                                                       \tag{1}
\]

Indeed differentiating (1+x_i u+u^2) in log(z_i) gives
(z_i-z_i^(-1))u, so differentiation of the character generating
function gives the lower-rank character in (1). This is a specialization
and elementary differentiation of known character/weight formulas,
not a new construction of the Springer map.

For completeness, the short-root coordinates can also be factored.
For i!=j put s=k-1 and

\[
 R(u)=(1-u^2)\prod_{a\ne i,j}(1+x_a u+u^2),\qquad
 A=[u^s]\,(1+u^2)R(u),\quad B=[u^{s-1}]R(u).
\]

With eta_i=z_i-z_i^(-1), (1) gives 2I y_i=eta_i(A+x_j B)
and 2I y_j=eta_j(A+x_i B). Direct simplification yields

\[
 2I(y_i-y_j)=(z_i-z_j)
       \left[(1+(z_i z_j)^{-1})A+\frac{2(z_i+z_j)}{z_i z_j}B\right],
\]
\[
 2I(y_i+y_j)=(z_i z_j-1)
       \left[\frac{z_i+z_j}{z_i z_j}A
             +\frac{2(z_i z_j+1)}{z_i z_j}B\right].            \tag{2}
\]

The proof below needs only the long-root coordinate 2y_i for k>=2
and the first short-root coordinate for k=1. Neither (1) nor (2)
identifies the vanishing locus of the torus Jacobian.

## 3. A constructive zero lemma for smaller fundamental characters

**Lemma.** For every r>=1 and 1<=s<=r, there is a regular torus
element of Sp_(2r) at which chi_s^(r)=0. Its inverse-pair sums x_i
can be chosen rational, distinct, different from 2,-2, and also all
different from 5/2. The z_i themselves need only be quadratic algebraic.

**Proof.** First expand the generating function in the x variables:

\[
 \chi_s^{(r)}(x)=
 \sum_{j=0}^{\lfloor s/2\rfloor}
 \left[\binom{r-s+2j}{j}-\binom{r-s+2j}{j-1}\right]
 e_{s-2j}(x_1,\ldots,x_r).                                  \tag{3}
\]

To derive it, choose a=s-2j factors contributing x_i u and j of the
remaining r-a factors contributing u^2, and subtract the coefficient
in degree s-2. Since r-s+2j>=2j, every bracket in (3) is strictly
positive. Thus this polynomial is multiaffine with nonnegative
coefficients, contains e_s with coefficient one, and all its monomials
have total degree congruent to s modulo two.

Set x_j=tj for 2<=j<=r and write

\[
 \chi_s^{(r)}(x_1,2t,3t,\ldots,rt)=A(t)x_1+B(t).
\]

For t>0, A(t)>0: its leading contribution is
t^(s-1)e_(s-1)(2,3,...,r), which is positive also at s=1.
Here the empty e_0 is one. Also B(t)>=0. The powers in A(t) have
parity s-1, those in B(t) have parity s. Therefore B(t)-2A(t)
is a nonzero polynomial of degree at most s: its nonzero A terms
cannot cancel the opposite-parity B terms.

Among the s+1 integers t=2,3,...,s+2, choose the first with
B(t)-2A(t)!=0; such an integer must exist. Put

\[
 x_1=-B(t)/A(t),\qquad x_j=tj\quad(2\le j\le r).
\]

Then x_1<=0 and x_1!=-2, while every other x_j is distinct and at
least four. Consequently all x_i are distinct, none equals 2,-2
or 5/2, and chi_s^(r)(x)=0. Choose any root z_i of
z_i^2-x_i z_i+1=0 for each i. The roots are nonzero, never 1,-1,
and no two chosen entries can be equal or inverse, because their
pair sums x_i are distinct. These are the required regular torus
coordinates. The construction includes r=s=1: x_1=0 and z_1=+i
or -i. This proves the lemma with at most s+1 rational trials.

## 4. Failure for k>=2

Apply the lemma with r=n-1 and s=k-1. Let w_1,...,w_(n-1) be any
of its chosen torus coordinates and take

\[
 z=(2,w_1,\ldots,w_{n-1}).
\]

The pair sum for the added entry is 2+1/2=5/2, excluded from all the
other pair sums by the lemma. Thus t(z) is regular in Sp_(2n).
Equation (1) gives

\[
 y_1=\frac{2-1/2}{2I_{n,k}}\chi_{k-1}^{(n-1)}(w)=0.
\]

The long root 2epsilon_1 therefore vanishes on the Lie algebra image,
which is nonregular. This proves every k>=2 in every n>=2, with an
explicit bounded recipe using rational pair sums and quadratic roots.

For example n=2,k=2 gives z=(2,i). Its four eigenvalues
2,i,1/2,-i are distinct, but y=(0,5i/4), so the image has a zero
long-root coordinate. The input is algebraic, not rational; the pair
sums (5/2,0) are rational.

## 5. Failure for the defining representation and the rank-one control

The defining map is already explicit in
[Kumar--Rogers, Part II, Lemma 5](https://kumar.math.unc.edu/papers/kumar81.pdf),
printed page 347. The coordinate formula used here is that same trace
map, after putting the inverse diagonal entries in paired order.
For k=1, I=1 and (1) becomes y_i=(z_i-z_i^(-1))/2. For n>=2 take

\[
 z=(2,-1/2,3,4,\ldots,n).
\]

The first two entries are neither equal nor inverse; the remaining
positive entries are distinct, at least three, and introduce no equality,
inverse equality or entry 1,-1. This is a regular torus element, but
y_1=y_2=3/4. The short root epsilon_1-epsilon_2 vanishes on the
image. This extends the already checked Sp4 defining witness uniformly;
it is not a new Sp4 example.

For n=1 there is no such pair of coordinates. Sp2=SL2 and the defining
Springer formula preserves the full centralizer, as in KM Section 6.2.
This proves the theorem, with the rank-one exception explicitly included.

## 6. Consequences and attribution boundary

KM Theorem 3.2 already identifies global regularity preservation with
the multiplication isomorphism

\[
 \mathbb C[G]^G\otimes\theta_k^*\operatorname{Harm}(\mathfrak g)
 \longrightarrow\mathbb C[G].
\]

Thus for n>=2 none of these fundamental representations supplies that
particular harmonic complement without localization. This does not rule
out different harmonic complements, change the polynomial-character or
Borel maps in the Kumar papers, or compute any generic degree or Smith
factor. The statement does not classify arbitrary type C highest weights.

For even k, V_(omega_k) descends to PSp_(2n); the same failure transfers
along the central isogeny, since the trace maps and Lie centralizer
dimensions correspond. Odd k does not descend. This is a consequence
of the stated fundamental cases, not a classification of all PSp
irreducibles.

The [verifier](tools/verify_type_c_fundamental_regularity.py) constructs
the actual symplectic contraction matrices in exterior monomial bases,
computes their rational ranks separately in every weight block, and
reconstructs trace projection from the resulting kernel multiplicities.
This does not assume the primitive character subtraction to construct
the multiplicities. Numerical rational character values are additionally
checked by the type C Weyl alternant determinant. The long- and short-root
formulas and trace indices are checked against these kernel weights.

The algebraic witnesses are verified exactly in the nonzero free algebra
Q[w_1,...,w_r]/(w_j^2-x_j w_j+1), using w_j^(-1)=x_j-w_j.
Actual kernel-weight traces vanish there; no floating-point roots or
assumed generic numerical fibres are used. The finite scope, further
symbolic ballot-expansion checks, bounded-search trials and source hashes
are recorded in [the report](type-c-fundamental-regularity-verification.json).
These corroborate Sections 2--5; they are not an extrapolation proof,
formal proof certificate or external peer review.

The character construction is classical, the coordinate differentiation
is elementary, and the fundamental-character zero lemma supplies the
written all-rank obstruction. Searches on 22 October 2026 combined
Springer/Cayley regularity with symplectic/fundamental character zero
terms. They recovered the published Part II defining formula, but no
matching all-fundamental classification was located. The Fulton--Harris
kernel theorem was visually inspected in the scanned source; the book
has no extractable text layer. This limited search is not an exhaustive
literature audit. No claim of original publication
priority is made merely because an exact matching classification was
not located. The [broader literature audit](CAYLEY_LITERATURE_AUDIT.md)
retains its scope and limitations.

The next question at this stage was nonfundamental irreducible type C
highest weights. The subsequent [all-highest-weight proof](TYPE_C_REGULARITY.md)
settles it negatively: a universal central--Weyl witness works for every
irreducible at n>=2, and a separate rank-one branching argument handles
all nonfundamental weights. The positive multiaffine expansion (3) above
remains specific to fundamentals; its coefficient positivity is not
silently extended to arbitrary symplectic Schur characters.
