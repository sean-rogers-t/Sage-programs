# Arbitrary representations: torus regularity versus global regularity

Date: 22 October 2026. Status: written deduction with internal proof audit
and exact independent corroboration. Neither external review nor publication
novelty is established. Work is over C, for finite-dimensional rational
representations with nonzero differential. No assumption of irreducibility
or of one central character is made.

## 1. Answer to the mixed-character question

For a representation rho of a simple group G, write theta_rho for its
trace-projection Springer morphism, normalized by

\[
 \operatorname{tr}(d\rho(\theta_\rho(g))d\rho(H))
 =\operatorname{tr}(\rho(g)d\rho(H)).
\]

We distinguish the following properties:

* **Torus regularity:** theta_rho(T_reg) is contained in t_reg.
* **Global regularity:** theta_rho(G_reg) is contained in g_reg.

The first tests regular semisimple inputs only. It need not imply the
second: regular nonsemisimple elements must also be tested. Regular means
that the centralizer has the rank of G, not that the element is semisimple.

**Theorem C.** Let G=Sp_(2n), n>=2, and let chi_rho(t)=F(x), where
x_i=z_i+z_i^(-1) in the defining torus coordinates. Put

\[
 P(x)=\prod_{i=1}^n(x_i+2)=\det(1+t).
\]

Then rho preserves torus regularity if and only if

\[
 \boxed{\quad F(x)=A P(x)^m+B,\quad m\ge1,\quad A\in\mathbb Z_{>0}.
 \quad}                                                     \tag{1}
\]

Here B is an integer such that the displayed expression is an actual
representation character; it may be negative because P^m has trivial
summands which can be removed. Every such character occurs, for example
P^m is the character of (exterior^bullet C^(2n))^(tensor m).
Copies and trivial summands do not change the normalized map.

**Nevertheless no nontrivial representation of Sp_(2n), n>=2, preserves
global regularity.** For the torus-positive characters (1), an explicit
regular nonsemisimple input proves the failure. Thus mixed central
characters do repair torus regularity, but cannot repair global regularity.

**Theorem B.** If G=Spin_(2n+1), n>=3, no nontrivial representation
preserves even torus regularity, including arbitrary mixed-character sums.
At n=2, Spin_5 is Sp_4: torus-positive representations exist as in
Theorem C, but none preserves global regularity. Consequently no
nontrivial representation of Spin_(2n+1), n>=2, preserves global regularity.

These statements strengthen the earlier
[type C](TYPE_C_REGULARITY.md) and [type B](TYPE_B_REGULARITY.md)
single-central-character obstructions. Those notes remain proofs of their
stated subcases; their mixed-character exclusions are superseded here.
This does not classify type D, exceptional groups, rational Cayley maps,
minimum Cayley degree or integral Springer/Borel images.

## 2. Two elementary ingredients

Let T_k(v) be the monic polynomial defined by
T_k(u+u^(-1))=u^k+u^(-k) for k>0, with T_0=1 for an exponent-zero
coordinate. Inversion-invariant characters are positive linear combinations
of products of these polynomials: group their actual weight monomials
by independent sign orbits. The coefficients are weight multiplicities.

For v>=2, T_k(v)>0 and T_k'(v)>0 for k>0. For u>1,

\[
 T_k'(u+u^{-1})=k\frac{u^k-u^{-k}}{u-u^{-1}}>0,
 \qquad T_k'(2)=k^2>0.
\]

Consequently every nontrivial B/C character polynomial used below has
strictly positive partial derivatives on the closed cone where all its
variables are at least 2. There is a weight nonzero in each coordinate,
by faithfulness of the differential and Weyl symmetry. This also proves
strict positivity when other coordinates or the differentiated coordinate
equal 2. No positivity of all monomial coefficients in the new variables
is being assumed.

The second ingredient is the elementary hypersurface-unit lemma:
if a nonzero polynomial Q has no zero off the union of finitely many
linear hyperplanes, then every irreducible factor of Q defines one of
those hyperplanes. Over C this follows from unique factorization and
the Nullstellensatz: an irreducible zero hypersurface contained in a
finite union must equal one of its members. In particular Q is a constant
times a product of the defining linear forms. Finite checks are not
needed for this all-degree lemma.

## 3. Necessity in type C

Write chi=F(x), a symmetric polynomial. The full Cartan trace Gram matrix
is S id, S=sum_mu m_mu mu_1^2>0, so the known trace formula gives

\[
 y_i=S^{-1}(z_i-z_i^{-1})F_i(x),\qquad F_i=\partial_iF.       \tag{2}
\]

The regular torus in x-space is exactly

\[
 U_C=\{x:x_i\ne2,-2,\ x_i\ne x_j\ (i<j)\}.
\]

Every point of U_C has torus lifts, obtained from z_i^2-x_i z_i+1=0.
If theta preserves torus regularity, its long-root coordinates y_i
never vanish there. By (2), F_i has no zero on U_C. The hypersurface-unit
lemma allows only factors x_j-2, x_j+2 and x_j-x_k. Positivity on
x_j>=2 excludes x_j-2 and x_j-x_k: each has a zero in that cone.
Only x_j+2 factors remain. Permutation symmetry therefore gives integers
a,b>=0 and c>0, independent of i, with

\[
 F_i=c(x_i+2)^a\prod_{j\ne i}(x_j+2)^b.                    \tag{3}
\]

If b>0, equality of mixed partial derivatives for i!=j forces a=b-1:
after cancelling common nonzero monomials, it says
((x_i+2)/(x_j+2))^(a-b+1)=1 identically. Integrating (3) gives
F=(c/b)P^b+B, the required form.

If b=0, integration instead gives

\[
 F=\frac{c}{a+1}\sum_i(x_i+2)^{a+1}+B.
\]

This alternative cannot preserve short-root regularity. Set
h(X)=(X-2)(X+2)^(2a+1). Then y_i^2=(c/S)^2 h(x_i).
Choose x_1>2 and consider h(X)-h(x_1). Its degree is 2a+2>=2;
x_1 is a simple root because h'(x_1)>0. It has another root x_2,
different from x_1 and from +/-2 because h(+/-2)=0!=h(x_1).
Choose other x_j generically distinct, excluding +/-2 and x_1,x_2.
The torus lifts have y_1^2=y_2^2!=0. Inverting the second lift changes
only the sign of y_2, so a lift has y_1=y_2. Its input is regular but
its output is not. This excludes the b=0 alternative.

Finally A=c/b is a positive integer: the weight (m,...,m) in P^m
is a highest weight m omega_n of multiplicity one, since all coordinate
exponents are at most m. Its multiplicity in F is A. B is then integral
by uniqueness of the irreducible-character expansion. More precisely,
B may be any integer at least -A times the trivial multiplicity of P^m.

## 4. Sufficiency and the exact symplectic map

For F=A P^m+B, differentiating at identity x=(2,...,2) gives

\[
 S=2A m4^{nm-1},\qquad
 y_i=\frac{2P(x)^m}{4^{nm}}\frac{z_i-1}{z_i+1}.             \tag{4}
\]

The trace index relative to the defining trace form is S/2. Formula (4)
does not depend on A or B. On U_C its common scalar is nonzero and
z_i!=+/-1. The differences and sums of the rational factors are

\[
 \frac{z_i-1}{z_i+1}-\frac{z_j-1}{z_j+1}
 =\frac{2(z_i-z_j)}{(z_i+1)(z_j+1)},
\]
\[
 \frac{z_i-1}{z_i+1}+\frac{z_j-1}{z_j+1}
 =\frac{2(z_i z_j-1)}{(z_i+1)(z_j+1)}.
\]

All output roots are nonzero at every regular torus input. This proves
sufficiency and completes the torus classification.

On the open set det(1+g)!=0, the full map is

\[
 \theta_m(g)=\frac{2}{4^{nm}}\det(1+g)^m(g-1)(g+1)^{-1}.
                                                               \tag{5}
\]

It extends polynomially to all of G as

\[
 \boxed{\quad\theta_m(g)=\frac{2}{4^{nm}}\det(1+g)^{m-1}
                 (g-1)\operatorname{adj}(g+1).\quad}       \tag{6}
\]

The formula follows either by differentiating chi(g)=A det(1+g)^m+B
in the trace identity, or from (4), equivariance and density of regular
semisimple conjugacy classes. It is a polynomial multiple of the known
classical Cayley transform, with normalization d theta_e=id.

## 5. Why the globally regular property still fails

Use a symplectic direct sum of n planes. Let

\[
 g_0=J_2(-1)\oplus\bigoplus_{j=2}^n
                \operatorname{diag}(r_j,r_j^{-1}),\qquad
 J_2(-1)=\begin{pmatrix}-1&1\\0&-1\end{pmatrix},
 \quad r_j=j.
\]

This is a regular element: its semisimple centralizer is
Sp_2 times an (n-1)-dimensional torus, and its unipotent part in Sp_2
is regular. Thus its centralizer has dimension 1+(n-1)=n.
It is not semisimple, and det(1+g_0)=0.

For m>=2, (6) gives theta_m(g_0)=0. For m=1, put
D=prod_(j=2)^n(1+r_j)(1+r_j^(-1)). Formula (6) gives

\[
 \theta_1(g_0)=\frac{4D}{4^n}
     \begin{pmatrix}0&1\\0&0\end{pmatrix}\oplus0_{2n-2}.
\]

This rank-one nilpotent has Jordan type (2,1,...,1); its centralizer
has dimension 2n^2-n>n when n>=2. Alternatively its centralizer already
contains Sp_(2n-2) and the one-dimensional centralizer in Sp_2, whose
combined dimension exceeds n. Hence it is nonregular. This proves
global failure for every representation, including the torus-positive
ones, not just for the exterior-algebra example.

## 6. Type B: arbitrary central mixtures at rank at least three

Use the Spin torus half-cover a_i with z_i=a_i^2, as in the
[type B note](TYPE_B_REGULARITY.md), and now set b_i=a_i+a_i^(-1).
All weights have integer doubled coordinates. Inversion invariance
therefore expresses every character as a polynomial F(b). The half-cover
kernel consists of even sign changes of a. Hence F is invariant under
even sign changes of b, as well as permutations. In other words, this
polynomial has the signed-permutation symmetry D_n. This is a coordinate
symmetry of the character, not a claim that the group has type D.

Let S=sum_mu m_mu mu_1^2>0. Differentiation in a gives

\[
 y_i=\frac{a_i-a_i^{-1}}{2S}F_i(b).
\]

Regularity of the B torus is equivalent to

\[
 U_B=\{b:b_i\ne2,-2,\ b_i\ne b_j,-b_j\ (i<j)\}.
\]

Indeed z_i=1 iff b_i^2=4, and z_i=z_j or z_i z_j=1 iff
b_i^2=b_j^2. Every point of U_B lifts to a regular Spin torus input.
Since the roots epsilon_i occur, preservation requires F_i to be
nonzero on U_B. The hypersurface-unit lemma allows only factors
b_j+/-2 and b_j+/-b_k. The positive-cone argument excludes b_j-2
and b_j-b_k. Even sign changes preserve the zero divisor of F_i,
because F_i(epsilon b)=epsilon_i F_i(b). They send b_j+2 to b_j-2
by flipping j and another index. For n>=3, they also send b_j+b_k
to b_j-b_k by flipping k and an index outside {j,k}. Thus the remaining
allowed factors are excluded too.

Every F_i must therefore be a nonzero constant. But flipping i and
another index sends F_i to -F_i, which a nonzero constant cannot satisfy.
This contradiction proves Theorem B for arbitrary representations,
without assuming a scalar action of the center. In fact the argument
shows gcd(F_i, product_j(b_j^2-4) product_(j<k)(b_j^2-b_k^2))=1:
each derivative is nonconstant and its zeros occur away from the
torus boundary divisor.

The index outside {j,k} is essential. At n=2 the argument does not
apply, and Spin_5=Sp_4 supplies the exact exception. If V is the
five-dimensional vector representation and S_spin the four-dimensional
spin representation, the actual representation

\[
 V\oplus2S_{\rm spin}\oplus3\mathbb C
\]

has character (b_1+b_2)^2, corresponding to P for Sp_4, and preserves
torus regularity. Removing its three trivial summands gives the same
map. More generally the torus-positive Spin_5 characters are precisely
A(b_1+b_2)^(2m)+B, with the same actual-character condition as (1).
The symplectic isomorphism then proves global failure in rank two.

## 6a. The exact good locus for the symplectic family

The formulas give more than one bad input. For every m>=1 and n>=2,

\[
 \theta_m^{-1}(\mathfrak g_{\rm reg})
 =G_{\rm reg}\cap\{g:\det(1+g)\ne0\}.                     \tag{7}
\]

On the indicated open set put C(g)=(g-1)(g+1)^(-1). The rational
inverse is g=(1+C(g))(1-C(g))^(-1), so g and C(g) have exactly the
same group centralizer. Multiplication by the nonzero scalar
2 det(1+g)^m/4^(nm) does not change that centralizer. Thus theta_m
preserves and reflects regularity throughout this open set, including
regular nonsemisimple inputs.

On its complement, theta_m is zero for m>=2. For m=1 the adjugate
formula has matrix rank at most one. A rank-one element of sp_(2n)
is nilpotent because its trace is zero; all nonzero such elements have
Jordan type (2,1,...,1) and are nonregular for n>=2. Zero is nonregular
as well. This proves (7). This is an immediate consequence of the
explicit formula and the known classical Cayley inverse, not a separate
novel construction or a claim about the whole Jacobian-zero locus.

## 7. Audit, finite verification and literature status

The [exact verifier](tools/verify_mixed_regularity.py) and
[report](mixed-regularity-verification.json) use separate constructions:
Freudenthal characters versus sign-orbit/Chebyshev polynomials;
direct exterior subsets and tensor convolution versus binomial weights;
actual induced exterior matrices for the nonsemisimple boundary input;
and exact centralizer equations. Mixed type B examples use polynomial
factor coprimality and direct trace evaluation in a finite multiquadratic
algebra, not numerical root searches. These are finite corroborations
of the arguments above, not proofs by extrapolation.

Known inputs are the trace/character construction, the symplectic
classical Cayley transform, elementary polynomial factorization, the
root/weight lattices and standard regular-centralizer theory.
[Kostant--Michor, Theorem 3.2 and Sections 4 and 7](https://arxiv.org/pdf/math/0109066)
already give the harmonic criterion, Jordan-decomposition machinery,
orthogonal failure examples and the spin/Cayley relation. Under their
criterion our global negative result means their particular pulled-back
harmonic complement cannot give the unlocalized factorization for these
representations. It does not rule out a different harmonic complement.

In particular KM Section 7.8 explicitly recalls
tr(exterior^bullet a)=det(1+a) in the orthogonal setting. The same
elementary exterior-character identity used above is a known input,
not a new character formula. Our classification argument is the
positivity/factorization/integrability necessity proof, with the separate
nonsemisimple test needed for global regularity.

The targeted literature check on 22 October 2026 used searches for
Springer-morphism regularity, reducible generalized Cayley representations,
and exterior-algebra/determinant symplectic Cayley maps. The nearby
[Lu--Ye paper, version 2, Introduction and Problem 1.2](https://arxiv.org/pdf/2411.02071)
classifies a different property: applicability of the classical matrix
Cayley transform in a chosen representation, chiefly irreducibly. It
does not state the character criterion proved here. No matching arbitrary-
representation classification was located in this limited search.
That is not an exhaustive citation search or evidence sufficient to
claim novelty. The elementary exterior/Cayley identity alone should
not be advertised as a new construction.

The main advance within this investigation is resolving the previously
unclassified central mixtures and separating the torus question from
the full-group question. Rank one, type D and exceptional groups remain
outside these two theorems. No PDF or LaTeX manuscript is updated here.
