# A prime-local bound for literal Springer–Borel images

Status, 20 October 2026: written repository deduction from the normalized
trace source theorem and classical integral BGG/characteristic-map theory.
The good-prime assertion is a consequence of the established machinery,
not a claim of a new characteristic-map theorem. The projective-family
application below is elementary; its exact literature priority is not
established. See the [stepwise audit](LOCAL_IMAGE_AUDIT.md).

This completes the first conservative local benchmark in Aim A and
removes the unnecessary exclusion of primes dividing the Levi Weyl-group
order. It does not give every bad-prime Smith multiplicity.

## 1. Scope and notation

Let G be simply connected, simple and complex, with maximal torus T,
character lattice X, Weyl group W and positive-root count N.
Let theta_rho be normalized trace projection for a representation with
faithful differential, including reducible representations.
Use the [degree-one source theorem](TRACE_SOURCE_THEOREM.md):

\[
\mathcal C_\rho=\operatorname{Sym}_{\mathbb Z}(\Lambda_\rho),\qquad
\Lambda_\rho=\{\ell\in X\otimes\mathbb Q:
                         \theta_\rho^*\ell\in\mathbb Z[X]\}.
\]

For a standard parabolic P with Levi Weyl group W_L, put

\[
A_{\rho,L}=\operatorname{Sym}_{\mathbb Z}(\Lambda_\rho)^{W_L},
\quad J_{\rho,L}=\beta^P(A_{\rho,L}).
\]

These are the **full literal character sources**, not arbitrary selected
normalized generator algebras. Polynomial degree d maps to H^(2d).

Define the ordinary torsion index t_G for the specified weight lattice by

\[
\beta^B(\operatorname{Sym}_{\mathbb Z}^{N}X)
             =t_G\mathbb Z[\mathrm{pt}].
\]

Equivalently choose a homogeneous integral F of degree N with
partial_(w0)(F)=t_G. Existence is the definition of the positive generator.
For the additive divided differences and the torsion-index construction,
see [Demazure](https://doi.org/10.1007/BF01418790) and
[Calmès–Petrov–Zainoulline, §§5.1–5.2 and Theorem 6.4](https://bcalmes.perso.math.cnrs.fr/pro/recherche/publis/cobordismflag/cobordismflag.pdf).
The argument below is in ordinary integral cohomology; it does not
discard coefficient hypotheses for arbitrary oriented theories.

## 2. The linear source lies in the actual character lattice

**Lemma.** Lambda_rho is contained in X.

Write theta_rho^*(ell)=sum_mu a_mu e^mu with every a_mu integral.
Its differential at the identity, in logarithmic torus coordinates, is
sum_mu a_mu mu, an element of X. Since d theta_rho at 1 is the identity,
that differential is ell itself. Therefore ell belongs to X.

This direct argument does not require irreducibility or an index table.
The full-lattice statement from the source theorem now gives finite index.
Let

\[
e_\rho=\exp(X/\Lambda_\rho),\qquad
e_\rho X\subseteq\Lambda_\rho\subseteq X.
\]

In particular the comparison index delta in the proposal can be replaced
by this exponent for the present normalized sources; no rational
superlattice of X is needed here.

## 3. Integrality after taking Levi invariants

The full-flag characteristic map sends Sym_Z(X) to H^(2*)(G/B,Z).
The pullback pi^*:H^*(G/P,Z) -> H^*(G/B,Z) is primitive: in the integral
Schubert bases it sends each class indexed by a minimal right coset
representative to the corresponding full-flag class. Thus an integral
full-flag class lying in its rational image has an integral preimage.

For f in Sym_Z(X)^(W_L), its rational Borel image beta^P(f) pulls back
to the integral beta^B(f). Primitivity proves beta^P(f) is integral.
Since Lambda_rho is contained in X, every J_(rho,L) above is consequently
integral. This settles integrality in the stated **full-source** scope;
it is not a claim about an arbitrary rationally normalized generator list.

## 4. A torsion-index multiple of every parabolic Schubert class

Let w be a minimal right coset representative in W/W_L and put

\[
F_w=\partial_{w^{-1}w_0}F.
\]

Integral divided differences preserve Sym_Z(X). For i retained in W_L,
w has no right descent at i, so s_i is a left descent of w^(-1)w0.
The nil-Hecke relations give partial_i(F_w)=0, hence s_i(F_w)=F_w.
Therefore F_w is an **integral Levi-invariant polynomial**.

If v has length equal to length(w), the BGG coefficient formula and the
nil-Coxeter multiplication rule give

\[
\partial_v F_w=
\begin{cases}
t_G,&v=w,\\
0,&v\ne w.
\end{cases}
\]

Indeed the composite has total operator length N, and the only Weyl
element of that length is w0; a nonreduced composite is zero.
Consequently

\[
\beta^P(F_w)=t_G\sigma_w.
\]

This supplies parabolic invariant lifts directly. It uses no averaging
by |W_L|, so primes dividing that order are not additional exclusions.

## 5. General bound and good-prime theorem

**Theorem.** In every polynomial degree d and for every standard Levi,

\[
\boxed{
t_G e_\rho^d H^{2d}(G/P,\mathbb Z)
       \ \subseteq\ J_{\rho,L}^d\
       \subseteq H^{2d}(G/P,\mathbb Z).
}
\tag{1}
\]

The second inclusion is Section 3. For the first, F_w is homogeneous of
degree d=length(w), and
e_rho^d Sym_Z^d(X) is contained in Sym_Z^d(Lambda_rho).
Its invariance is unchanged by scalar multiplication. Applying beta^P
to e_rho^d F_w supplies t_G e_rho^d sigma_w for each Schubert basis vector.

Thus the full-source image has full rational rank. Its cokernel is finite
and killed by t_G e_rho^d; each nonzero Smith factor divides this integer.
In particular every occurring prime divides t_G e_rho. If
p does not divide t_G e_rho, then

\[
\boxed{
J_{\rho,L}\otimes\mathbb Z_{(p)}
       =H^{2*}(G/P,\mathbb Z_{(p)}).
}
\tag{2}
\]

There is no |W_L| factor in (1) or (2). Degree zero is exactly Z and has
no cokernel, regardless of the harmless t_G factor in the bound.

The formula need not be sharp at a bad prime. Its purpose is to prove
a uniform restriction before computing matrices. Further progress must
explain which allowed primes and exponents actually occur.

## 6. A complete projective family for the defining SL_n representation

Let G=SL_n, n>=2, rho its defining representation, P=P_1, and
G/P=projective (n-1)-space. Keep G simply connected: the representation
has not been replaced by a representation of PGL_n.

Use additive diagonal weights x_1,...,x_n, sum x_i=0, with x_1=omega_1,
and let Q be the root lattice. The standard trace formula gives

\[
\theta_\rho^*(x_i)=z_i-\frac1n\sum_j z_j,\qquad
\prod_j z_j=1.
\]

The preceding integrality lemma first puts ell in X.
Represent ell by integral coordinates (c_1,...,c_n) modulo the all-ones
vector. Its Laurent coefficients are c_i-(sum c_j)/n. They are integral
exactly when sum c_j is divisible by n, which is the criterion ell in Q.
Therefore Lambda_rho=Q and e_rho=n. The full polynomial source is Sym_Z(Q).

The differences r_j=x_1-x_j, 2<=j<=n, are an integral basis of Q.
The Levi group S_(n-1) permutes them, so

\[
A_{\rho,L}=\mathbb Z[E_1,\ldots,E_{n-1}],\qquad
E_j=e_j(r_2,\ldots,r_n).
\]

The integral presentation H^*(P^(n-1),Z)=Z[h]/(h^n), with beta(x_1)=h,
and the shift identity for elementary symmetric polynomials give

\[
\beta(E_j)=\binom n j h^j,\qquad 1\le j<n.
\tag{3}
\]

One way to see (3) is to compare coefficients in
prod_i(u+x_1-x_i). In ordinary cohomology the positive symmetric
functions e_k(x_1,...,x_n) vanish; the remaining coefficient is
binomial(n,j) h^j. These are classical characteristic-map calculations.

**Proposition.** The complete literal image is

\[
\boxed{
J_{\rho,L}^d=\frac{n}{\gcd(n,d)}\,\mathbb Z h^d
                 \quad(0\le d<n).
}
\tag{4}
\]

For d=0 the coefficient is one; for d>=n the target is zero.
Every degree has a single nonzero factor. This is an additional explicit
application, not a universal formula for all type A parabolics.

### Proof of the gcd assertion

Let g_d be the gcd of all products of binomial(n,j) with positive indices
whose sum is d. These are exactly the coefficients of the degree-d source
monomials in (3). Fix a prime p and put a=v_p(n).
The identity j*binomial(n,j)=n*binomial(n-1,j-1) gives

\[
v_p\binom n j\ge\max(0,a-v_p(j)).
\]

In a partition of d, some part has p-valuation no larger than v_p(d).
Thus every product has valuation at least max(0,a-v_p(d)).
If a=0, the product binomial(n,1)^d attains valuation zero.
If a>0, write d=q*p^a+r with 0<=r<p^a.
For 1<=r<p^a the low a base-p digits of n-1 are all p-1, so
binomial(n-1,r-1) is a p-unit, and
v_p binomial(n,r)=a-v_p(r).
Also binomial(n,p^a) is a p-unit whenever that index is used.
The partition consisting of q copies of p^a and, when r>0, the part r
therefore attains the lower bound. All indices used are below n because
d<n. This proves v_p(g_d)=max(0,v_p(n)-v_p(d)), hence (4).

The unit claim also follows directly from
binomial(n-1,r-1)=prod_(i=1)^(r-1)((n-i)/i):
for i<p^a, v_p(n-i)=v_p(i). Thus this proof needs no additional
assumption about the higher base-p digits of n.

For example SL_4/P_1 has factors (1,4,2,4), substantially sharper than
the degreewise bound 4^d. Primes coming from the source are genuine:
arbitrarily large prime n already occurs in divisor degree, although
the simply connected type A torsion index is one.

## 6a. Equivariant good primes and the support of specialization defects

Let S=Sym_Z(X) be the equivariant coefficient ring, using a separate copy
of the Cartan variables, and define

\[
J^T_{\rho,L}=S\cdot\beta_T^P(A_{\rho,L})
                  \subset H_T^{2*}(G/P,\mathbb Z).
\]

This is an S-subalgebra. It is integral: the full-flag characteristic
classes are integral, and equivariant Schubert pullback to G/B is
primitive over S, by the same basis argument as in Section 3.
The standard equivariant Schubert basis is homogeneous and free over S,
and forgetful specialization sends it to the ordinary basis.
These are published equivariant Schubert-theory inputs; see
[Calmès–Zainoulline–Zhong](https://arxiv.org/abs/1409.7111).

**Corollary.** If p does not divide t_G e_rho, then

\[
J^T_{\rho,L}\otimes\mathbb Z_{(p)}
             =H_T^{2*}(G/P,\mathbb Z_{(p)}).
\tag{5}
\]

Indeed put R=Z_(p). Since e_rho is a unit, Lambda_rho tensor R equals
X tensor R. Localization commutes with taking invariants here: invariants
are a kernel of the finitely many reflection-difference maps, and
localization is flat. Thus F_w/t_G is a legitimate localized Levi source.
Its equivariant image specializes to sigma_w, so homogeneity gives

\[
\beta_T^P(F_w/t_G)=\sigma_w^T+
     \sum_{\ell(v)<\ell(w)} c_{wv}(t)\sigma_v^T.
\]

The coefficients have positive parameter degree. Induction on length
puts each lower class, and then sigma_w^T, in the localized image.
This proves (5) without averaging by |W_L|.

Let I=S_+ and K be the kernel of the natural surjection
J^T_(rho,L)/I J^T_(rho,L) -> J_(rho,L). In each fixed total polynomial
degree, K_d is a finite abelian group and its primes divide t_G e_rho.
Finiteness follows from finite generation of the graded ambient groups
and vanishing of the rational kernel. The prime assertion follows by
localizing (5), where specialization is the usual ambient isomorphism.
No sharp bad-prime exponent for K is asserted.

Thus faithful specialization is proved at good primes, while
[the SL4 example](EQUIVARIANT_IMAGE_PILOT.md) shows a genuine integral
failure at an allowed bad prime.

## 7. Exact checks, literature status and next step

The [Sage verifier](experiments/local_image_equivariant_pilot.py) computes
top torsion-index polynomials in seven small types and checks all their
parabolic Schubert lifts. It checks (3) by BGG on SL_n/P_1 for 2<=n<=8
and the coefficient gcd through n=80. The
[independent SymPy/stdlib verifier](experiments/verify_local_image_independent.py)
rechecks the coefficient gcd and the separate equivariant counterexample.
See the [Sage report](experiments/local-image-equivariant-pilot.json) and
[independent report](experiments/local-image-independent.json).
The calculations corroborate these proofs; they do not replace them.

The good-prime framework and torsion-index lifts are published theory
applied here. The projective image also agrees with the root-characteristic
problem for the adjoint root datum, so it should be compared with that
literature before claiming priority. Targeted searches on 20 October 2026
did not establish a precise prior citation for (4) in this Springer
formulation. That is not evidence of originality.

The main remaining arithmetic target is a sharp bad-prime comparison
outside this rank-one cohomology family. The
[equivariant pilot](EQUIVARIANT_IMAGE_PILOT.md) already shows why faithful
specialization cannot be assumed as a general method for that target.

The existing lattice API supplies e_rho without enumerating representation
weights. For example:

```python
from sage_borel import irreducible_character_lattice
minimal = irreducible_character_lattice("E6", (1, 0, 0, 0, 0, 0))
adjoint = irreducible_character_lattice("E6", (0, 1, 0, 0, 0, 0))
assert minimal.linear_basis_matrix.elementary_divisors()[-1] == 18
assert adjoint.linear_basis_matrix.elementary_divisors()[-1] == 24
```
