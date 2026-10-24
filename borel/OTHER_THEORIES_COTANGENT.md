# Character differentials, primitive K-classes and Springer image data

Exploration, 22 October 2026. This note explains a connection to established
representation-ring and K-theory constructions. It supplies finite data
for the representation-dependent image problem; it does not classify
Springer morphisms as maps of varieties or determine their generic degrees.

## 1. The same character differential

Fix a finite simply connected simple group G and one invariant form
normalization, with long roots of length squared two. Write I_rho for
the trace index and chi_rho for its character. On the torus,

\[
 d\chi_\rho=\sum_\mu m_\mu e^\mu\otimes\mu
        \in R(T)\otimes X.
\]

The explicit Springer formula says exactly that theta_rho^*(ell) is
the contraction of this logarithmic differential with ell under the
invariant form, divided by I_rho. This is the published trace/weight
formula from Rogers's paper, expressed as a differential.

The differential itself is an established K-theoretic object:
[Brylinski--Zhang, Theorem 3.2](https://arxiv.org/abs/dg-ga/9710035)
identifies the exterior algebra of representation-ring differentials
with conjugation-equivariant K-theory of G in the simply connected case.
This is a connection to **odd equivariant K-theory of G**, not an
identification with K^0(G/P). The latter enters through the separate
associated-bundle map.

## 2. What survives in a flag K-image

Let m_G be the augmentation ideal of R(G). The associated-bundle map
sends every global G-character to its dimension. Consequently the
coefficients of d chi_rho reduce at augmentation. The surviving datum is

\[
        \tau_\rho=I_\rho^{-1}[\chi_\rho-\dim\rho]
                    \in(m_G/m_G^2)\otimes\mathbb Q.                \tag{1}
\]

The [finite comparison](OTHER_THEORIES_FINITE.md) proves that the pair
(Lambda_rho,tau_rho) determines the embedded literal K-image for every
parabolic. This is sufficient comparison data; it is not claimed that
distinct pairs always have distinct image subrings.

Hodgkin's theorem identifies m_G/m_G^2 with primitive odd K-classes of
the compact group. This precise formulation and an explicit algorithm
for its coordinates appear in
[Davis, From representation theory to homotopy groups, Section 2, equation (2.1)](https://www.lehigh.edu/~dmd1/e8.pdf).
Thus (1) can equally be viewed as a normalized primitive K-class. The
primitive groups, their Adams operations and the differentiated
representation-ring algorithm are published results, not new invariants
introduced by this project.

If chi_rho=P(chi_1,...,chi_r) in fundamental characters and d_i their
dimensions, its coordinates before division by I_rho are
(partial P/partial chi_i)(d_1,...,d_r). An effective recursive way to
obtain them is to differentiate

\[
 [V_\mu][V_{\omega_i}]=[V_{\mu+\omega_i}]
                           +\sum_{\nu<\mu+\omega_i}m_\nu[V_\nu].
\]

The leading highest weight occurs once; every other weight is lower in
dominance order, so the recursion terminates. This is the algorithm
already used by Davis. A trace-index functional provides an independent
consistency check: I_rho=sum_i a_i I_(omega_i).

### Uniform representation operations

Dimensions d and trace indices give the exact identities

\[
 \tau_{\rho\oplus\sigma}
      ={I_\rho\tau_\rho+I_\sigma\tau_\sigma\over I_\rho+I_\sigma},
 \qquad
 \tau_{\rho\otimes\sigma}
      ={d_\sigma I_\rho\tau_\rho+d_\rho I_\sigma\tau_\sigma
                          \over d_\sigma I_\rho+d_\rho I_\sigma}.   \tag{2}
\]

In particular every positive tensor power rho^(tensor k) has the same
normalized tau as rho. This does not say their integral source lattices
are equal. It says that their rational K-coordinate maps are equal;
on a common literal integral subsource they agree integrally.
Their actual Springer morphisms on G can still differ, because global
characters become dimensions only after taking the associated K-class.

The same sufficiency extends to the specified trace-coordinate lift in
any complex-oriented theory, by passage through MU. Rationally its
coordinate is a linear combination of directional derivatives of the
homogeneous moments of chi_rho(exp x)/I_rho, with coefficients from the
formal exponential of the chosen orientation. A difference in m_G^2
expands into products of positive W-invariant series. Differentiation
leaves a positive invariant factor, which vanishes in the rational
Borel quotient. Hence equal tau gives equal universal MU coordinates
rationally, then integrally by torsion-freeness, and under every
orientation specialization. This assertion is for the specified lift
on a common integral source, not for every possible generalized Springer
construction or an enlarged modular invariant source.

## 3. A G2 check against the published primitive calculation

Use fundamental representations of dimensions 7 and 14, with indices
2 and 8. The representation-ring Adams operation psi_R^2 on
indecomposables, with columns in that order, is

\[
             \begin{pmatrix}12&-208\\-2&56\end{pmatrix}.            \tag{3}
\]

It follows directly from the exterior-square identities in Davis's
Section 2. Its degree-two and degree-six left eigenvectors can be
chosen (2,8) and (1,-26), with eigenvalues 4 and 64. Thus the ratio of
the normalized degree-six indices is (-26/8)/(1/2)=-13/2.

The [independent Todd/BGG audit](experiments/audit_other_theories_k.py)
computes the degree-five Springer coordinate corrections from actual
weight moments and obtains exactly this ratio. In its fundamental-weight
and length-five Schubert bases, the minimal matrix is diag(1/60,-1/60)
and the adjoint matrix is diag(-13/120,13/120).
This cross-check is independent of the Laurent Demazure computation of
the complete K-images. Note the operation convention: Davis also uses
K^1 rather than K^(-1), which changes an Adams factor by k. Formula
(3) is explicitly the representation-ring operation, so no grading
conventions are silently identified.

## 4. Exceptional data and an API path

The [Sage research prototype](experiments/other_theories_cotangent_pilot.py)
implements the published recursion and checks trace indices and tensor
squares. With fundamental reduced-character classes e_i in Sage/Bourbaki
numbering, examples are

| G and highest weight | Dimension | I | Class before division by I |
| --- | ---: | ---: | --- |
| G2, 2omega_1 | 27 | 18 | 13e_1-e_2 |
| F4, 2omega_4 | 324 | 162 | 51e_4-e_1-e_3 |
| E6, omega_1+omega_6 | 650 | 300 | 27e_1-e_2+27e_6 |
| E7, 2omega_7 | 1463 | 660 | 112e_7-e_1-e_6 |
| E8, 2omega_8 | 27000 | 13500 | 495e_8-e_1-e_7 |

These are representation calculations, not new exceptional flag-image
Smith tables. [Report](experiments/other-theories-cotangent-pilot.json).
For example 27 tensor 27* in E6 decomposes as 1+78+650, so its line
above can be verified immediately by differentiating that tensor identity.
The more general prototype check uses the degree-two trace functional;
it is not an independent all-higher-index audit in every exceptional type.

A practical API extension would return Lambda_rho, I_rho, the integral
indecomposable vector and tau_rho with their normalization, then a separate
K-image object with its support filtration, Smith/Fitting data and
comparison map. Keeping those stages separate makes clear which data
are already computed and which embedded arithmetic remains to be done.
