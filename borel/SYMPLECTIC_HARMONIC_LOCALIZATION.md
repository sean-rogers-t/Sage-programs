# Sharp harmonic localization for the symplectic determinant family

Date: 22 October 2026. Status: written deduction with exact corroboration,
not externally reviewed or established novel. The sufficiency argument is
a direct consequence of the published Kostant separation theorem and the
classical Cayley isomorphism. The minimal-support argument and explicit
Jacobian calculation are recorded separately below. Everything is over C.

## 1. Setup and statement

Let G=Sp_(2n)(C), n>=2, g=sp_(2n), N=dim g=n(2n+1), and put

\[
 D(g)=\det(I+g),\qquad
 \theta_m(g)=\frac{2}{4^{nm}}D(g)^{m-1}(g-I)\operatorname{adj}(g+I),
 \qquad m\ge1.
\]

These are exactly the normalized Springer maps for the torus-positive
characters A D^m+B classified in [the mixed-representation note](MIXED_REGULARITY.md).
Write A_G=C[G], R_G=A_G^G. Let H=Harm(g) be the usual graded Kostant
harmonic **vector space** in C[g], and H_m=theta_m^*H in C[G].
The coefficient ring R_G is the full conjugation-invariant coordinate
ring, as in Kostant--Michor; it is not restricted to theta_m^*C[g]^G.
All tensor products below are over C. They express a module decomposition,
not a tensor product of algebras: H and H_m are not being called subalgebras.

**Theorem.** Multiplication induces an R_G[D^(-1)]-linear, G-equivariant
isomorphism

\[
 \boxed{\quad R_G[D^{-1}]\otimes H_m
          \ \xrightarrow{\ \sim\ }\ A_G[D^{-1}].\quad}     \tag{1}
\]

Moreover, for every nonzero invariant f in R_G, multiplication
R_G[f^(-1)] tensor H_m -> A_G[f^(-1)] is an isomorphism **if and only
if D divides f in R_G**. Thus D=0 is the minimal closed support that
must be removed among these invariant principal localizations. Powers
of D give the same localization; no distinguished denominator exponent
or classification of arbitrary nonprincipal open sets is asserted.

The normalized full Cayley determinant, computed using the left-invariant
input frame gH and the same fixed Lie algebra output basis, is

\[
 \boxed{\quad
 \Psi_m(g)=\left(\frac{D(g)}{4^n}\right)^{(nm-1)(2n+1)}
      \left(1-nm+\frac m2\operatorname{tr}g\right).
 \quad}                                                    \tag{2}
\]

Consequently (1) holds across an additional critical divisor that is
removed by the full-Jacobian localization in Kostant--Michor. This is
about multiplication/separation. Their inverse-differential vector
fields generally still have poles on that additional divisor.

## 2. A direct proof of sufficiency using known separation

The classical Cayley map and its inverse are

\[
 X=C(g)=(g-I)(g+I)^{-1},\qquad
 g=(I+X)(I-X)^{-1}.
\]

They give an equivariant isomorphism

\[
 G_D\simeq \mathfrak g_\delta,\qquad
 \delta(X)=\det(I-X),\qquad D(g)=\frac{4^n}{\delta(X)}.
\]

Here det(I+X)=det(I-X) for symplectic X, so the inverse indeed lies
in G. In these coordinates the Springer map is simply

\[
 \theta_m(g)=2\delta(X)^{-m}X.                              \tag{3}
\]

Kostant's published separation theorem gives a graded harmonic space H
such that multiplication C[g]^G tensor H -> C[g] is an isomorphism.
Since delta is invariant, localizing gives

\[
 C[\mathfrak g]^G[\delta^{-1}]\otimes H
       \simeq C[\mathfrak g][\delta^{-1}].
\]

For a homogeneous h in H of degree d, (3) gives the exact equality

\[
 (\theta_m^*h)(g)=2^d\delta(X)^{-md}h(X).                   \tag{4}
\]

The multiplier is an invariant unit. Choose a homogeneous vector-space
basis of H and rescale each basis member by its corresponding unit.
This changes neither the span nor the linear independence over the
localized invariant ring. Transporting the multiplication isomorphism
through C proves (1), both surjectivity and injectivity.

Although H is infinite-dimensional, there is no infinite-sum operation:
polynomials and tensor elements involve finitely many basis members.
Homogeneity is essential; H is chosen with its standard grading.
Localization at an invariant commutes with invariants here because the
denominator is fixed and the coordinate rings are domains. This identifies
R_G[D^(-1)] with C[g]^G[delta^(-1)] under the Cayley isomorphism.

This proof does not need theta_m to be locally invertible everywhere
on G_D. It is a formal use of known separation and homogeneous rescaling,
not a new proof of Kostant's theorem.

## 3. Necessity: a conjugacy orbit that the proposed generators cannot separate

Use torus pair sums x_i=z_i+z_i^(-1). The usual restriction theorem gives

\[
 R_G=\mathbb C[e_1(x),\ldots,e_n(x)],\qquad
 D=\sum_{k=0}^n2^{n-k}e_k(x),\quad e_0=1.
\]

In particular D is prime in this polynomial ring: it is monic linear
in e_n. If D does not divide f, the restriction of f to the quotient
divisor D=0 is not zero.

A dense family in that divisor is represented by the regular elements

\[
 g_0=J_2(-1)\oplus\bigoplus_{j=2}^n\operatorname{diag}(r_j,r_j^{-1}),
 \qquad J_2(-1)=\begin{pmatrix}-1&1\\0&-1\end{pmatrix},
\]

where the remaining pair sums are generic, distinct and different from
2,-2. They parameterize the multisets with one pair sum equal to -2;
these are dense in D=0 on the quotient. Thus some such regular g_0
satisfies f(g_0)!=0.

The previous explicit formula gives theta_m(g_0)=0 for m>=2, and for
m=1 a nonzero nilpotent supported only on the first symplectic plane.
Let k act identically on that plane, interchange the two eigendirections
in the second plane by a symplectic matrix, and fix all other planes.
Then g_1=kg_0k^(-1) differs from g_0, but

\[
 \theta_m(g_1)=\theta_m(g_0),\qquad f(g_1)=f(g_0)\ne0.
\]

Every invariant coefficient takes the same value on g_0,g_1, and every
member of H_m takes the same value because it factors through theta_m.
Therefore every element in the image of the proposed multiplication
map on G_f takes equal values on these two points. Matrix coordinate
functions distinguish them. Hence the multiplication map is not surjective.

Conversely if D divides f, then D is a unit in R_G[f^(-1)], and further
localization of (1) proves the desired isomorphism. This proves the iff
statement. The assumption n>=2 is used in the second symplectic plane;
the globally positive SL2 defining case is not a counterexample to the
theorem with its stated rank restriction.

## 4. The full Jacobian: the extra factor is explicit

Set F(X)=2 delta(X)^(-m)X. Differentiation gives

\[
 dF_X(Y)=2\delta^{-m}
       \left(Y-m\delta^{-1}d\delta_X(Y)X\right).
\]

The rank-one determinant lemma therefore gives, writing E for the Euler
vector field X -> X on g,

\[
 \det dF_X=2^N\delta^{-mN}\left(1-m\frac{E\delta}{\delta}\right).
                                                               \tag{5}
\]

For the group-to-algebra Cayley coordinate change, differentiate on a
left-invariant tangent vector gH:

\[
 dC_g(gH)=\tfrac12(I+X)H(I-X),\qquad
 \det(dC_g)=2^{-N}\delta^{2n+1}.                             \tag{6}
\]

For completeness, identify sp_(2n) with Sym^2(C^(2n)) using the
symplectic form. The transformation H -> (I+X)H(I-X) then becomes
S -> (I+X)S(I+X)^t. Its determinant is
det(I+X)^(2n+1)=delta^(2n+1); the factor 1/2 contributes 2^(-N).
Thus (6) uses the actual full Lie algebra differential, not merely the
torus derivative.

Finally

\[
 \frac{E\delta}{\delta}=-\operatorname{tr}(X(I-X)^{-1})
       =n-\tfrac12\operatorname{tr}g.
\]

Multiplying (5) and (6) yields (2) on G_D, hence everywhere on G by
density, since both sides are regular functions. The exponent is
mN-(2n+1)=(nm-1)(2n+1)>=0. The normalization is Psi_m(I)=1.

This determinant calculation also agrees with the tensor-power scaling
identity already published in Kostant--Michor Theorem 2.8(3). Neither
tensor-power scaling nor the rank-one determinant lemma is a new input.

## 5. Harmonic separation across ramification

Write Q_m(g)=1-nm+(m/2)tr(g). In R_G the factors D and Q_m are
distinct primes when n>=2: D depends linearly on e_n and Q_m depends
linearly on e_1. On G_D the critical locus is precisely Q_m=0.
It is nonempty even among regular semisimple inputs. Explicitly set

\[
 x_i=2-\frac{4i}{mn(n+1)},\qquad 1\le i\le n,
\]

and choose roots z_i of z_i^2-x_i z_i+1=0. These pair sums are distinct
and lie strictly between 0 and 2, so the torus input is regular and
D=prod_i(x_i+2)!=0. Their sum is 2n-2/m, so Q_m=0. The image remains
regular by the centralizer-preserving Cayley formula. In Cayley coordinates
(5) is an identity minus a rank-one operator with one-dimensional kernel
spanned by X at these points, so the full differential has corank one.

Thus criticality of the map is strictly stronger than failure of regularity
preservation. Harmonic multiplication remains an isomorphism on G_D,
including Q_m=0. By contrast pulling back a constant vector field through
the inverse differential generally requires inverting Q_m as well;
the statement about separation does not extend that inverse across ramification.

## 6. Exact corroboration and what follows from published work

The [verifier](tools/verify_symplectic_harmonic_localization.py) and
[report](symplectic-harmonic-localization-verification.json) check the
Jacobian using direct full matrix differentials and, independently,
actual induced exterior-representation matrices and trace projection.
They include nonsemisimple inputs, regular critical inputs with singular
differential but regular image, explicit conjugacy pairs on D=0, and
low-degree genuine symplectic harmonic bases. They corroborate the
uniform argument; finite checks do not prove the full decomposition.

Known inputs include Kostant's separation theorem and the classical
Cayley isomorphism. [Kostant--Michor Section 3.2](https://arxiv.org/pdf/math/0109066)
already prove harmonic separation after inverting the full Jacobian, and
their proof explains why regularity preservation is the relevant issue.
Their Theorem 2.8 gives direct-sum and tensor-product scaling of the maps;
Section 7.8 explicitly recalls the exterior-character determinant identity.
Equation (1) here follows directly from those standard ingredients and
homogeneous rescaling. It should be described as an explicit corollary,
not a new general separation theorem.

The additional conclusions recorded here are the minimal invariant
principal localization for this particular family, the exact normalized
Jacobian factorization, and the separation/ramification comparison.
A targeted search on 22 October 2026 for symplectic exterior Cayley
Jacobian formulas and harmonic determinant localization located no matching
explicit statement; this limited check does not establish novelty or an
open problem in the literature. No claim is made about other harmonic
complements, positive characteristic, integral cohomology or minimum
Cayley degree. No PDF, LaTeX manuscript or runtime library API is changed.

## 7. Adversarial audit of the torus/global distinction

The following distinctions were rechecked directly on 22 October 2026.
No gap was found in Sections 2-3 under their stated hypotheses. This is
an internal audit, not external validation of the full regularity classification.

On the defining torus, the roots are z_i^2 and z_i/z_j, z_i z_j, with
their inverses. Thus a regular torus element has z_i != -1 for every i,
and lies in G_D. This explains why the determinant family is torus-positive.
It does not imply that every regular group element lies in G_D:
the single Jordan block at -1 in Section 3 is regular, with Lie
centralizer dimension n, even though its semisimple part is not regular.

There is an exact pointwise statement, stronger than a generic assertion:

\[
 \theta_m^{-1}(\mathfrak g_{\mathrm{reg}})
       =G_{\mathrm{reg}}\cap G_D,\qquad n\ge2.                \tag{7}
\]

On G_D, Cayley inversion recovers g from X. Since theta_m(g) is a nonzero
scalar times X, their centralizers are equal. On D=0, the adjugate has
rank at most one; theta_1 has rank at most one and theta_m=0 for m>=2.
A nonzero rank-one element in sp_(2n) is square-zero: writing it as
v ell, the symplectic condition forces ell proportional to the pairing
with v, so ell(v)=0. Its adjoint orbit has dimension 2n (the minimal
nilpotent orbit), and its centralizer has dimension 2n^2-n>n for n>=2.
Zero also has centralizer dimension greater than n. This proves (7).
In particular every boundary image is nonregular, not merely the
generic boundary images used to prove necessity of localization.

The minimal-support proof does not replace the full invariant ring by
the image of Lie-algebra invariants. Its coefficient ring is R_G.
The primeness of D is needed in that ring, where its monic linear
expression in e_n proves it. The dense quotient boundary family then
supplies a conjugate pair with f!=0 for every invariant f not divisible
by D. Equal theta images obstruct surjectivity independently of any
choice of harmonic basis. Sufficiency uses invariant-unit rescaling
of a homogeneous Kostant basis, and needs neither an inverse of d theta
nor regularity of every point of G_D.

The strengthened exact audit adds 48 regular-torus checks of equality
of centralizer subspaces, 36 boundary-stratum checks (a Jordan block,
a semisimple -1 plane, and two -1 Jordan blocks), nine independent
actual exterior/tensor trace projections on the Sp4 boundary, and
three Sp2 controls. The Sp2 m=1 control stays regular on that boundary;
this verifies why n>=2 is essential. These checks supplement, rather
than prove, the all-rank arguments above.

## 8. Other separation results and the scope of the obstruction

Here separation means a multiplication isomorphism over the full
invariant coordinate ring. Existence of a G-stable complement, an
explicit complement, and a differential-operator description of that
complement are different requirements.

* **Abstract global complements for all simply connected semisimple
  groups are known.** Richardson's result, recalled with precise
  hypotheses in [Agricola--Goodman, Theorem 5.1](https://arxiv.org/pdf/math/0207161),
  applies to conjugation on G. Thus Sp and the exceptional groups have
  some H with C[G]^G tensor H -> C[G] an isomorphism. The negative
  regularity classification rules out H=theta_rho^*Harm(g) in its scope;
  it does not rule out these other complements.
* **An explicit differential-operator realization for SL2 is published.**
  [Agricola--Goodman, Theorem 5.3 and Remark 5.1](https://arxiv.org/pdf/math/0207161)
  give such a harmonic space and identify it with the defining Springer
  pullback. This is not an additional simple-group positive case.
* **Every admissible Springer map has Jacobian-localized separation**
  under [Kostant--Michor Theorem 3.2](https://arxiv.org/pdf/math/0109066).
  The Sp result above improves the required localization for this
  particular family; it does not supply the first localized existence theorem.
* **GL_N has an elementary global example outside the simple-group table.**
  For the defining representation theta(g)=g in gl_N. Localize the
  usual graded Kostant decomposition of C[gl_N] at the invariant det:

  \[
    \mathbb C[GL_N]^{{GL_N}}\otimes\operatorname{Harm}(\mathfrak{gl}_N)
       \simeq\mathbb C[GL_N].
  \]

  Harmonics here are restricted to GL_N. This is a direct known-input
  consequence, not a new classification of reductive representations.
* **Classical Cayley open sets give more explicit local complements.**
  For SO_N as well as Sp_(2n), C(g)=(g-I)(g+I)^(-1) identifies G_D
  with g_delta, delta=det(I-X). Transporting the graded Kostant basis
  gives C[G]^G[D^-1] tensor C^*Harm(g) -> C[G][D^-1]. For SO_N this
  rational classical Cayley pullback has not been identified here with
  the defining trace-projection Springer pullback. These two maps must
  not be interchanged. Products of established examples likewise give
  separation by tensoring their complements; no claim about arbitrary
  representations of a product follows.

More generally, if an equivariant isomorphism C:U -> g_delta is given
and theta=u C for an invariant unit u in C[U]^G, then the same rescaling
proof works: a homogeneous harmonic of degree d becomes u^d C^*h.
This transfer observation explains the examples rather than asserting
a new existence theorem. A scalar-Cayley formula after a finite covering,
such as the Spin construction in KM Section 7.15, additionally requires
checking the covering's coordinate ring; an isomorphism on SO cannot
simply be treated as an isomorphism on Spin.

These alternative complements are background, not the current research
scope. At the user's request the program remains tied to generalized
Cayley maps and their fixed spaces theta_rho^*Harm(g). The next degree
step is now proved in [the symplectic degree note](SYMPLECTIC_CAYLEY_DEGREES.md):
every actual character F(D) has generic degree 2n deg(F), with an
explicit fiber equation. Sharp localization for other generalized
Cayley maps and degree laws for the Spin family remain relevant next
questions. No independent global-complement construction is being pursued.
