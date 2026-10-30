# Regularity classification for all irreducible type C Springer maps

Subsequent result: the [arbitrary-representation classification](MIXED_REGULARITY.md)
resolves the mixed-central-character question left outside this note.
Mixed sums can preserve torus regularity, precisely for characters
A det(1+g)^m+B, but every nontrivial representation fails global regularity
at rank at least two. The present note retains its narrower theorem and proofs.

Date: 22 October 2026. Status: a uniform written deduction, with an
internal argument audit and independent exact character/weight checks.
The trace construction, highest-weight theory, complete reducibility and
SL2 characters are published inputs. Publication priority and external
review of this classification remain unverified. The coefficient field
is C; finite central kernels are allowed. No integral image, generic
degree, positive-characteristic or mixed-central-parity classification is
asserted here.

## 1. Statement

Let G=Sp_(2n)(C), and let V_lambda be a nontrivial irreducible rational
representation. Write its dominant highest weight in orthonormal
coordinates as the integer partition

\[
 p=(p_1\ge p_2\ge\cdots\ge p_n\ge0),\qquad L=p_1,\quad D=\sum_i p_i.
\]

In fundamental-weight notation lambda=sum_i a_i omega_i, this means
p_j=sum_(i=j)^n a_i and L=sum_i a_i. Its differential is faithful.
Let theta_lambda:G -> sp_(2n) be the normalized trace Springer map.

**Theorem.** If n>=2, no nontrivial irreducible V_lambda preserves torus
regularity. Every such map has an algebraic regular semisimple torus
element with nonregular image. Hence none preserves global regularity.
If n=1, the map preserves torus and global regularity precisely for the
defining representation, p=(1).

Consequently every Sym^d(C^(2n)), d>=1, fails when n>=2. At n=1 only
d=1 is positive. The theorem is not merely a computation of symmetric
powers: the proof covers every dominant highest weight.

In fact, the n>=2 obstruction holds for every nontrivial finite-dimensional
representation on which -I acts by a single scalar, even if it is reducible.
This condition includes all irreducibles and direct sums of the same
central parity. Mixed central parities are not classified here.

## 1a. A universal central--Weyl witness

The following elementary observation gives a short proof at every n>=2.
Let c be central, rho(c)=s id, t in T, and w in W with w(t)=ct. The
trace identity gives theta_rho(ct)=s theta_rho(t), and equivariance gives

\[
 w\theta_\rho(t)=s\theta_\rho(t).
\]

Thus if t is regular but the s-eigenspace of w on t contains no regular
Lie algebra element, regularity preservation fails. This is an elementary
consequence of known trace equivariance, not a new equivariance theorem.

For Sp_(2n) let c=-I and r=floor(n/2). Take the torus coordinates

\[
 z=(2,-2,3,-3,\ldots,r+1,-(r+1))
 \quad\text{if } n=2r,
\]
\[
 z=(2,-2,3,-3,\ldots,r+1,-(r+1),i)
 \quad\text{if } n=2r+1.
\]

All paired entries have absolute value greater than one and are distinct;
the optional i has absolute value one. None is 1,-1 and no two chosen
coordinates are equal or inverse. Thus the input t(z) is regular.
Let w swap each displayed pair and, in odd rank, invert the final
coordinate. This is an element of the signed-permutation Weyl group.
It sends every coordinate to its negative, including i^(-1)=-i, so
w(t)=ct. On Lie algebra coordinates it swaps the same pairs and
negates the final coordinate in odd rank.

For an irreducible representation s=(-1)^D. More generally s=+1 or -1
whenever -I acts by one scalar. Since n>=2, there is a first pair, and
the eigenvector relation gives

\[
 \boxed{\qquad y_1=s y_2.\qquad}
\]

For s=+1, the short root epsilon_1-epsilon_2 vanishes; for s=-1,
epsilon_1+epsilon_2 vanishes. Either way the image is nonregular.
This proves the rank-at-least-two theorem with a single input independent
of the highest weight. Its entries are rational in even rank and lie in
Q(i) in odd rank. Sections 2--5 retain a different highest-weight proof
for the nonfundamental cases and also settle the rank-one statement.

## 2. Character coordinates and elementary root bounds

Use torus and Lie coordinates

\[
 t=\operatorname{diag}(z_1,\ldots,z_n,z_1^{-1},\ldots,z_n^{-1}),
 \quad \theta_\lambda(t)=\operatorname{diag}(y_1,\ldots,y_n,-y_1,\ldots,-y_n),
 \quad x_i=z_i+z_i^{-1}.
\]

Input regularity is x_i!=2,-2 and x_i!=x_j. An output with any y_i=0
is nonregular because the long root 2epsilon_i vanishes.

The weight character is invariant under each z_i -> z_i^(-1), so it
is a polynomial F_p(x_1,...,x_n). To see this without a determinant
formula, use the monic polynomials z^a+z^(-a) in x=z+z^(-1), generated
by the recurrence T_(a+1)=x T_a-T_(a-1). The invariants of all independent
inversions are C[x_1,...,x_n]. Permutation invariance also makes F symmetric.

Every weight mu satisfies |mu_i|<=L. Indeed H_1 with standard diagonal
entries 1,-1 on the first symplectic pair is dominant: every positive
root has nonnegative value on it. The highest-weight ordering gives
mu_1<=p_1, while Weyl sign/permutation symmetry gives both bounds for
every coordinate. Therefore F_p has degree at most L separately in
each x_i, and degree exactly L in x_1.

Let I_p>0 be the trace index relative to the defining trace form;
the Cartan trace Gram matrix in the pair basis is 2I_p times identity.
The published trace/weight construction gives

\[
 \boxed{\quad y_i=\frac{z_i-z_i^{-1}}{2I_p}\,
                   \frac{\partial F_p}{\partial x_i}(x).\quad}       \tag{1}
\]

This is elementary differentiation of the character, using the
[Kostant--Michor trace definition](https://arxiv.org/pdf/math/0109066)
and [Rogers's explicit weight formula](https://arxiv.org/pdf/1701.01538v1).
It is not a new construction of theta. Positivity of I_p follows from
2I_p=sum_mu m_mu mu_1^2>0.

## 3. Rank-one restriction supplies the needed positivity

Restrict V_lambda to the commuting block subgroup
Sp2 x Sp_(2n-2). For n=1 the second factor is trivial. Complete
reducibility gives an actual representation decomposition

\[
 V_\lambda=\bigoplus_{a=0}^{L}\operatorname{Sym}^a(\mathbb C^2)\otimes M_a,
 \qquad M_L\ne0.                                             \tag{2}
\]

The last assertion follows because the highest weight p has H_1 value
L, and there are no larger H_1 weights. Let P_a(x)=U_a(x/2) denote the
SL2 character, with P_0=1, P_1=x, P_(a+1)=x P_a-P_(a-1).
The character identity from (2) is

\[
 F_p(x_1,x')=\sum_{a=0}^L P_a(x_1) A_a(x'),\qquad A_a=\chi_{M_a}.
\]

Although A_a can have negative coefficients when expressed in x', its
value is strictly positive whenever M_a!=0 and all other z_j are
positive real: it is then a sum of positive weight monomials. For every
a>=1, P_a'(x)>0 on x>=2. For x=z+z^(-1)>2, this follows by pairing
the weights of Sym^a(C^2):

\[
 P_a'(x)=\sum_{\substack{b>0\\b\text{ a weight}}}
                  b\frac{z^b-z^{-b}}{z-z^{-1}}>0.
\]

At x=2 its limit is binomial(a+2,3)>0. Thus, for every choice of
x_j>2, j>=2,

\[
 Q_p(x_1,x'):=\partial_{x_1}F_p(x_1,x')>0\quad\text{when }x_1\ge2,
 \qquad \deg_{x_1}Q_p=L-1.                                  \tag{3}
\]

The leading coefficient is L A_L(x')>0. These facts use actual
multiplicity modules, not an unsupported positivity claim about every
coefficient of a symplectic Schur polynomial.

## 4. The central character prevents a boundary-only obstruction

The central element -I acts on V_lambda by (-1)^D. Replacing every
z_i by -z_i replaces every x_i by -x_i, so

\[
 F_p(-x)=(-1)^D F_p(x),\qquad
 Q_p(-x)=(-1)^{D-1}Q_p(x).
\]

In particular

\[
 Q_p(-2,x')=(-1)^{D-1}Q_p(2,-x').                            \tag{4}
\]

Equation (3) shows Q_p(2,x') is not the zero polynomial. Hence (4)
shows Q_p(-2,x') is not the zero polynomial either. This nonvanishing
as a polynomial is the step needed to exclude x_1=-2 as a forced root.
It does not assert Q_p(-2,x')!=0 at every positive tuple.

There is a bounded constructive choice. For j=2,...,n let z_j range
over the L+1 disjoint positive integers

\[
 S_j=\{2+(j-2)(L+1),\ldots,1+(j-1)(L+1)\}.
\]

Their pair sums x_j=z_j+z_j^(-1) are distinct, greater than two, and
give L+1 distinct values for each variable. The nonzero polynomial
Q_p(-2,x') has degree at most L in each variable, so it cannot vanish
on the entire product grid (successive univariate interpolation proves
this). Choose a tuple with Q_p(-2,x')!=0. At most (L+1)^(n-1)
exact rational trials are needed. At n=1 the empty grid is valid and
(4) directly gives Q_p(-2)!=0.

## 5. Every L>=2 now has a regular-torus witness

Fix the tuple from Section 4 and put q(X)=Q_p(X,x'). It has degree
L-1>=1 and rational coefficients. Moreover

\[
 q(-2)\ne0,\qquad q(2)>0,\qquad q(x_j)>0\quad(j\ge2).
\]

Choose any complex algebraic root xi of q. Then xi differs from
2,-2 and every other pair sum x_j. Choose z_1 as either root of
z_1^2-xi z_1+1=0. All z_i are nonzero, none is 1,-1, and no two
are equal or inverse. Thus t(z) is regular, but (1) gives y_1=0.
This proves failure for every highest weight with L>=2, at all ranks.

If L=1, the nonzero partition p=(1^k,0^(n-k)) is fundamental. For
n>=2, the [fundamental-type-C proof](TYPE_C_FUNDAMENTAL_REGULARITY.md)
already supplies exact all-rank witnesses for every k. For n=1 the
sole L=1 case is the defining SL2 representation, whose Springer map
preserves its entire centralizer. This completes all cases of the theorem.

## 6. Symmetric powers and an explicit example

Sym^d(C^(2n)) is the irreducible highest weight d omega_1. Its character
is h_d in the 2n eigenvalue alphabet and its trace index is

\[
 I_{n,d}=\binom{2n+d}{2n+1}.
\]

One can obtain the same obstruction polynomial directly by taking
the indicated coefficient in

\[
 q(X)=[u^{d-1}]\,(1-Xu+u^2)^{-2}
                      \prod_{j=2}^n(1-x_j u+u^2)^{-1}.       \tag{5}
\]

This is the elementary character derivative, not a different Springer
map. The general proof supplies a valid tuple and root for every d>=2;
the d=1 case at n>=2 is the previously checked defining obstruction.

For Sym^2(C^4), the adjoint representation of Sp4, take z_2=2 and
xi=-5/4. Here F=x_1^2+x_2^2+x_1x_2-2, I=6, and q(X)=2X+5/2.
With z_1=(-5+i sqrt(39))/8,

\[
 \theta_{2\omega_1}\bigl(\operatorname{diag}(z_1,2,z_1^{-1},1/2)\bigr)
       =\operatorname{diag}(0,15/32,0,-15/32).
\]

The input is regular and the zero output coordinate is a long-root
obstruction. The verifier checks the example directly from weights.

## 7. Consequences, known inputs and verification boundary

All negative witnesses transfer to the adjoint central quotient PSp_(2n)
for representations that descend. In particular for n>=2 no nontrivial
irreducible representation of either isogeny form preserves regularity.
More strongly, every representation of PSp_(2n) pulls back with the
single central scalar +1, so every nontrivial finite-dimensional PSp
representation is negative, including reducible ones. This holds also
at n=1: use z_1=i and w inversion, so w(t)=-t while w(y)=-y; the
scalar +1 forces y=0. The transfer uses equality
of Lie centralizer dimensions and naturality of trace projection under
a central isogeny, not a new descent theorem.

By the published KM Theorem 3.2, at n>=2 no irreducible type C trace
map gives the unlocalized multiplication isomorphism using the particular
space theta_lambda^* Harm(g). This does not exclude other harmonic
complements or invalidate the polynomial-character/Borel construction.
It says nothing about the minimum degree over all rational Cayley maps.

The same universal input is also a critical point of the full map.
For a root alpha vanishing on theta(t), differentiate equivariance:
d theta_t([E_alpha,t])=[E_alpha,theta(t)]=0. The input vector is
nonzero because t is regular. Thus the full differential has a kernel.
This is a consequence of the standard root-direction differential
identity, not a claim that the torus Jacobian must be singular there.

The central--Weyl proof also applies to representations with a single
central parity; it does not prove the result for a direct sum mixing the
two central characters, because theta_rho(ct) is then not a scalar
multiple of theta_rho(t). No assertion about such mixed sums is made.

The mathematical inputs above are standard highest-weight theory,
rank-one characters, complete reducibility, the central character,
the Rogers trace formula and the previously written fundamental proof.
Fulton--Harris Chapter 17 supplies the classical symplectic representation
setup. Section 1a gives the universal central--Weyl deduction; Sections
3--5 give a separate branching proof for the nonfundamental cases.
Targeted searches
on 22 October 2026 for Cayley/Springer regularity and symplectic symmetric
powers or irreducibles did not locate an exact matching classification;
many search results concern unrelated finite-group characters or graph
Cayley maps. This is not a certificate of novelty or a complete literature
audit. The [earlier audit](CAYLEY_LITERATURE_AUDIT.md) retains its qualifications.

The [exact checker](tools/verify_type_c_regularity.py) evaluates the
universal witness directly from actual weights in Q(i), and compares the
Chebyshev/Weyl alternant slice with Freudenthal weight multiplicities,
reconstructs the actual rank-one strings, checks full Cartan trace
coordinates, central parity, bounded-grid witnesses and their root
coordinates in an exact quadratic algebra over Q[X]/(q). Direct symmetric
monomial enumeration adds a third character check for symmetric powers.
The [report](type-c-regularity-verification.json) records the finite scope
and source hashes. The arbitrary-rank proof is Sections 2--5, not an
extrapolation from these finite checks or an external peer review.

Next direction: type B and type D. Their character lattices and spin
representations need separate treatment; the polynomial pair-sum and
central-parity argument above does not automatically classify them.
