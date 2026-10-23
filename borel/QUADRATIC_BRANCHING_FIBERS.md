# Quadratic Levi branching and generic Springer fibers

Date: 2026-10-23. Characteristic-zero written deductions with exact independent
checks. Standard inputs are the Jacobian criterion, dimension of images,
generic smoothness, and classical character identities. Literature priority
is unverified. The type A adjoint representation is generally not fundamental.

## 1. A quadratic theorem without a constant-rank hypothesis

Let T=(C*)^s, x in A^ell, and

    chi(x,z)=a(z)+f(z)^t x + (1/2) x^t Q(z)x,                    (1)

where a,f,Q have Laurent-polynomial coefficients and Q is symmetric. Write
D_j=z_j partial/partial z_j. Keep the actual scheme

    X=V(Qx+f) in A^ell x T,
    mu:X -> A^s,  mu_j=D_j chi.                                 (2)

For c in A^s the fiber mu^(-1)(c) is exactly the critical scheme of

    W_c=chi-sum_j c_j log z_j.                                  (3)

Only the logarithmic differential is used; no global choice of logarithms
is needed. The following assertion includes ell=0 and Q=0.

**Theorem 1.** There is a dense open U in A^s such that every fiber over U
is smooth, possibly empty. Its local dimension on a contributing component
X_0 is dim X_0-s. The critical scheme is Morse-Bott: the Hessian kernel is
exactly its tangent space. No constant-rank condition on Q, reducedness
condition on X, or smoothness condition on rank-drop loci is required.

If X_0 has generic Q-rank k and its image closure in T has codimension q,
then

    dim mu^(-1)(c) = ell-k-q,                                   (4)
    rank Hess(W_c) = s+k+q.

Only components whose map to A^s is dominant contribute. In particular,
a negative value of ell-k-q cannot correspond to a contributing component.

### 1.1 The key exclusion argument

Partition the compatible parameter locus

    {z: rank Q(z)=rank[Q(z) | f(z)]}

into finitely many smooth locally closed strata Y with constant Q-rank k.
On each such stratum X_Y -> Y is an affine-space bundle of relative
dimension t=ell-k. Indeed, invert a k by k minor, solve for k variables,
and use compatibility for the remaining rows. Put d=dim Y.

At (x,z) in X_Y define the ell by s mixed-derivative matrix

    B_ij = D_j (Qx+f)_i,
    p = rank[Q | B]-k.

The tangent-space dimension of the **ambient incidence scheme X** is

    dim T_(x,z) X = ell+s-rank[Q | B] = t+s-p.                  (5)

The local Krull dimension of X is at least t+d, because X_Y is a locally
closed affine bundle through that point. If X is nonsmooth there, the
Jacobian criterion gives

    t+s-p > dim_(x,z) X >= t+d, hence d+p<s.                    (6)

On the affine fiber X_z, the differential of mu is B^t restricted to
ker Q. Symmetry of Q identifies its rank with p:

    rank(d(mu|X_z)) = rank(B^t|ker Q) = rank[Q | B]-k.           (7)

Stratify the nonsmooth locus further by p. In characteristic zero the
image of each of its fibers over Y has dimension at most p: take smooth
strata of the reduced underlying fiber and use the differential-rank
bound. Consequently its total target image has dimension at most d+p<s.
There are finitely many strata. Their image closures form a proper closed
subset of A^s.

This excludes all nonsmooth incidence points from generic target fibers,
including points of nonreduced components and intersections of components.
The scheme X has not been replaced by its radical. Reduction is used only
to bound images of its bad locus; the actual fiber ideals remain (2)-(3).

### 1.2 Finish the proof and compute dimensions

Apply characteristic-zero generic smoothness to mu on X_sm, and exclude
the images of its nondominant components. After removing the bad target
set above, all remaining fibers of the original scheme X are smooth.
Their dimensions are dim X_0-s. At a critical point the tangent equations
are the Hessian of W_c, so smoothness identifies its kernel with the
tangent space of the fiber. This is the Morse-Bott assertion.

For the generic parameter point of the image of X_0, the compatible
linear equations give an affine space of dimension ell-k. An irreducible
component of X cannot have only a proper part of that affine space as
its generic fiber: the full affine bundle over a suitable dense stratum
would have a closure containing that component. Thus

    dim X_0 = (s-q)+(ell-k),

proving (4). This completes the proof. If s=0, X is itself an affine
space or empty, giving the same conclusion directly.

### 1.3 Constant-rank elimination and why it is not enough alone

On a constant-rank open set, choose a nondegenerate k by k symmetric
block A in Q after a fixed linear change if needed. Split x=(y,u) and
write Q=[[A,C],[C^t,D]], f=(b,e). Completing the square gives

    y' = y+A^(-1)(Cu+b),
    chi = (1/2)y'^t A y' + a-(1/2)b^t A^(-1)b
                  + (e-C^t A^(-1)b)^t u,

because rank k forces the Schur complement D-C^t A^(-1)C=0.
After eliminating y', this is the previous affine theorem, with ell-k
remaining variables. Theorem 1 additionally controls rank-changing loci
using the full derivatives B; restricting chi to such a locus before
differentiating would lose equations.

For example chi=(z-1)x^2/2 has generic critical scheme

    z=1, x^2=2c,

so the rank-drop locus supplies the entire nonempty generic fiber.
Discarding det Q=0 would give an incorrect empty answer.

The degree restriction is substantive for arbitrary potentials:
chi=x^3+z gives critical ideal (x^2,z-c) for every c!=0, which is
nonreduced. This is an abstract control, not a representation character
or a counterexample for fundamental Springer maps.

## 2. Consequence for the Springer map

Let G be simple complex simply connected and rho an actual nontrivial
representation. Normalize its map by

    B_rho(theta_rho(g),X)=tr(rho(g)d rho(X)),
    B_rho(Y,X)=tr(d rho(Y)d rho(X)).

Fix a proper Levi K, put S=[K,K], ell=rank S, s=dim Z(K)^0 and
r=rank G=ell+s. Work on S x Z(K)^0 -> K with fundamental characters
x of S and torus coordinates z. Suppose chi_rho|K has degree at most
two in x. There is a dense open set of central targets h with Z_G(h)=K
for which

    theta_rho^(-1)(h) intersect G_reg

is smooth, possibly empty. For a contributing incidence component with
generic rank k and parameter-image codimension q its dimension is

    dim K-r+ell-k-q.                                           (8)

This follows from Theorem 1 and the previously proved
[regular-quotient critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md).
The target coordinates c are related to h by the nondegenerate trace
form. Smoothness and dimensions descend under the finite central etale
cover. The statement recovers the [affine theorem](AFFINE_BRANCHING_FIBERS.md)
at Q=0. Wholly nonregular source components and special targets are not
classified by this assertion.

## 3. Type A adjoint branching in every rank

Let G=SL_N and rho=Ad on sl_N. For N>=3 its highest weight is
omega_1+omega_(N-1), so it is not one of the fundamental modules in
the preceding affinity classification. The trace form and map are

    B_Ad(X,Y)=2N tr(XY),
    theta_Ad(g)=[tr(g^(-1))g-tr(g)g^(-1)]/(2N).                 (9)

To verify (9), use End(C^N)=C plus sl_N, Ad(g)=g tensor g^(-t),
and ad(X)=X tensor 1-1 tensor X^t. Multiplying and taking traces gives
tr(Ad(g)ad(X))=tr(gX)tr(g^(-1))-tr(g)tr(g^(-1)X).
The scalar summand contributes zero. The displayed numerator is traceless.
This fixes the normalization, rather than just proportionality of maps.

For K=S(product GL_(b_i)), sum b_i=N, write g_i=z_i s_i on the central
cover, with det s_i=1 and product z_i^(b_i)=1. Use the identity component
of this torus when the determinant equation is disconnected. Let x_(i,j)
be the j-th fundamental character of SL_(b_i). Set

    A=tr(g)=a_0+u^t x,       a_0=sum_(b_i=1) z_i,
    B=tr(g^(-1))=b_0+v^t x, b_0=sum_(b_i=1) z_i^(-1),

where u has coefficient z_i at x_(i,1), and v has coefficient z_i^(-1)
at x_(i,b_i-1). For b_i=2 these are the same variable. Then

    chi_Ad=AB-1,
    Q=uv^t+vu^t, f=b_0 u+a_0 v, a=a_0 b_0-1.                  (10)

Thus the branching is quadratic for every Levi and rank Q<=2.
More precisely:

- If some b_i>=3, u and v are independent everywhere and rank Q=2.
- If every nontrivial block has size two, let t be their number.
  For t=1 the rank is one. For t>=2 it is two except on the locus
  z_i^2=z_j^2 for all two-blocks, where it is one.
- If t=0 there are no semisimple variables and Q is empty.

In the two-block case rank never falls to zero, because Q has diagonal
entries 2. The rank-one locus has codimension t-1 when size-one blocks
are present. These are statements in arbitrary rank, not interpolations.

## 4. Complete generic regular-fiber theorem for SL_N adjoint

**Theorem 2.** Let K be a proper Levi with block sizes b_i. For a generic
semisimple target h with centralizer K, the regular part of theta_Ad^(-1)(h)
is nonempty if and only if every b_i is at most two and at least one b_i
equals one. In that case it is smooth of pure dimension

    2t,  where t=number of two-blocks.                          (11)

For the excluded block patterns the regular part is empty for **every**
target with centralizer K, not merely generically. No assertion is made
that the entire fiber is empty: nonregular source points can remain.

### 4.1 Exclude larger blocks and the all-two-block pattern

The semisimple equations are u B+v A=0. If a block has size >=3,
independence of u,v forces A=B=0. All central derivatives of AB-1 then
vanish, so h=0. This contradicts Z_G(h)=K proper.

If there are only two-blocks, on the rank-two locus the same argument
applies. On the rank-one locus write z_i=epsilon_i w with epsilon_i=+/-1.
There are no external scalar contributions, so A=w U and B=w^(-1)U,
where U=sum epsilon_i x_i. The semisimple equation forces U=0, again
giving A=B=0 and target zero. This proves the exclusion.

An independent matrix explanation of the larger-block bound comes from
(9): on a target eigenvalue block, g satisfies a scalar polynomial of
degree at most two unless A=B=0. A regular GL_b block is cyclic and has
minimal polynomial of degree b, so b<=2 is necessary for nonzero targets.

### 4.2 Parametrize the contributing locus

Suppose there are t>=1 two-blocks and p>=1 singleton blocks, N=2t+p.
Every contributing point lies on the rank-one locus, with common
two-block determinant w^2. Write their traces as w s_i and the singleton
entries as w v_j. Retain the finite determinant cover

    w^N product_j v_j=1.

The semisimple equation and the traces give

    2 sum_i s_i + sum_j(v_j+v_j^(-1))=0,
    L=(1/2)sum_j(v_j-v_j^(-1)),
    A=w L, B=-L/w, H=-L/(2N).

Equation (9) now gives target eigenvalues

    lambda_i=H s_i       on each two-block,
    nu_j=H(v_j+v_j^(-1)) on each singleton.                     (12)

Their trace condition 2 sum lambda_i+sum nu_j=0 is automatic.
These formulas retain all ambient derivatives, including directions
transverse to the equal-determinant locus.

### 4.3 Prove dominance, rather than infer it from dimensions

It suffices to show that the map v -> (nu_1,...,nu_p) is dominant;
the t-1 independent s_i then give the other target coordinates when H!=0.
At v_1=...=v_p=2,

    H=-3p/(8N),
    (D_k nu_j)=-(9p I_p+25 J_p)/(16N),                          (13)

where J_p is the all-ones matrix. Its eigenvalues are -9p/(16N),
with multiplicity p-1, and -34p/(16N), with multiplicity one.
They are nonzero. The determinant cover is finite etale, so it does
not change dominance. Thus the target image contains a dense open set;
intersect it with the open condition that the block target eigenvalues
are distinct. The equality of singleton eigenvalues at the test point
does not cause a problem: it proves dominance before taking that open.

Regular lifts exist: choose a cyclic GL_2 matrix for each trace and
determinant. It is K-regular, and for a lift lying over a target with
centralizer K, equivariance gives Z_G(g)=Z_K(g), hence G-regularity.

Here ell=t, k=1 and q=t-1, so (8) gives dim K-r=2t. All contributing
components have that dimension. If t=0, K=T; d theta at the identity
is the identity by the trace normalization, so the torus map is dominant
and generically etale. This gives nonempty generic regular fibers of
dimension zero. Theorem 2 follows.

## 5. Transfer to classical adjoint representations

The general quadratic theorem does not depend on type A. For types B,C,D
write the defining module on a Levi central cover as

    V = direct_sum_i (E_i plus E_i^*) plus W,

with central twists understood. The adjoint modules are exterior^2 V
in orthogonal types and Sym^2 V in type C. Use the standard direct-sum
identities for these powers. Products of defining characters of different
Levi factors have degree at most two whenever each defining character
is fundamental or constant. Internal SL_b adjoint/symmetric-square
characters are quadratic, and exterior squares are fundamental or trivial.

Let m be the residual classical rank and t the number of GL blocks of
size at least two. As in the affinity classification, absorb Spin_2 into
the central torus in type D, so m=0 or m>=2 there.

**Theorem 3 (coverage).** The adjoint character is of degree at most two
in all semisimple fundamental-character coordinates precisely as follows:

| Type and group | Quadratic coverage |
| --- | --- |
| C_n, Sp_(2n) | Every proper Levi |
| B_n, Spin_(2n+1) | Every proper Levi except m=1 and t>0 |
| D_n, Spin_(2n) | Every proper Levi except m=2 and t>0 |

Thus generic regular Springer fibers are smooth or empty for every pair
covered by this table. Formula (8) applies, but their ranks, compatibility
codimensions, nonemptiness and exceptional target sets are not classified
here. The adjoint representations are fundamental for B_n at n>=3 and
D_n at n>=4, so those rows give applications to the original fundamental
representation question. C adjoint and B_2 adjoint are not fundamental.

Proof of the qualifications. In residual C_m, the vector character is
fundamental and the adjoint character is x_1^2-x_2-1 at m>=2, or x_1^2-1
at m=1. In residual B_m with m>=2 the vector is fundamental; the adjoint
is fundamental for m>=3 and quadratic for B_2. In B_1=SL_2 the vector
character is y^2-1, not y. A nontrivial GL block therefore creates the
nonzero cubic term (z_i x_(i,1)+z_i^(-1)x_(i,b_i-1)) y^2.
No other term in exterior^2 V can cancel this multidegree.

In residual D_m with m>=3 the vector is fundamental (also D_3=A_3);
the adjoint is fundamental at m>=4 and quadratic at m=3. At m=2 the
group is SL_2 x SL_2, its vector character is yz, and its adjoint is
y^2+z^2-2. A nontrivial GL block creates a nonzero cubic x_i yz term.
When t=0 these cross terms have central coefficients and degree two.
When m=0 the orthogonal residual module is zero in D and a trivial
line in B. This proves both sufficiency and necessity in the table.

The cubic exclusions are exclusions from the theorem's hypothesis, not
claims of singular fibers. They provide a sharply defined next family
for testing whether the smoothness conclusion extends beyond quadratics.

## 6. Verification and literature status

The [verification program](tools/verify_quadratic_branching_fibers.py)
checks scheme ideals in exact arithmetic; all type A Levi block partitions
for N=2,...,8; symbolic dominance determinants; and classical branching
coverage against characters reconstructed independently from roots.
The [matrix audit](tools/audit_type_a_adjoint_fibers.py) constructs Ad and
ad on traceless matrix units directly and checks the trace normalization,
Springer coordinates, and full differentials in SL_3, SL_4 and SL_5.
It does not call Sage's character or quotient code. Reports are
[quadratic checks](quadratic-branching-fibers-verification.json) and
[direct matrix checks](type-a-adjoint-fibers-audit.json).

The completed checks comprise eight exact scheme controls, 58 proper
type A block partitions through N=8, 24 exact dominance determinants,
and 347 classical root-restriction polynomials through rank six. The last
set contains 37 cubic exclusions, all agreeing with Theorem 3. The
root-restriction check forgets central characters; the written tensor
argument proves the central-coefficient statement. The direct matrix
audit checks 865 trace-form pairings, 142 Springer trace pairings and
2,882 full differential pairings at eight points. One source is regular
but nonsemisimple. A separate SL_4 regular zero-fiber control confirms
that the target-zero case belongs to K=G, outside the proper-Levi theorem.

Proof-audit points: (6) uses the dimension of the entire incidence scheme,
not just a chosen component; (7) uses symmetry of Q; image bounds use
characteristic zero; generic smoothness is applied only after excluding
the nonsmooth source locus. Formula (8) still adds the smooth group-orbit
directions, and the adjoint type A result does not turn an empty regular
part into a claim about emptiness of the whole fiber. Neither test-point
counts nor successful specialized Jacobians replace the all-rank arguments.

Known inputs and checks performed:

- [Kostant--Michor, The generalized Cayley map](https://www.mat.univie.ac.at/~michor/kostant.pdf)
  provides the trace-projection construction and general framework. Its
  SL_N example in section 6.2 is the standard representation, not the
  adjoint-family theorem above. Section 6.1 also treats SL_2 representations.
- The [Stacks Project treatment of smooth morphisms](https://stacks.math.columbia.edu/tag/01V4)
  and the classical characteristic-zero generic smoothness theorem supply
  the smoothness tools. The mixed-derivative dimension bound above is a
  written deduction, not a replacement for these standard results.
- [Etingof's classical fundamental-representation notes](https://ocw.mit.edu/courses/18-755-lie-groups-and-lie-algebras-ii-spring-2024/mit18_755_s24_lec05.pdf)
  provide the exterior and symplectic fundamental models. Adjoint tensor
  descriptions, low-rank isomorphisms and square completion are classical.

Targeted searches on 2026-10-23 included quadratic critical schemes and
generic smoothness, quadratic generating families and Morse-Bott behavior,
and adjoint generalized Cayley fibers for SL_N. They did not identify an
exact published version of the statements here. This is not an exhaustive
priority audit and does not establish originality. Work on real trace
Morse-Bott functions and on Higgs-bundle Cayley correspondences is not
automatically a result about these complex trace-projection fibers.

The output is proof notes and research checks, not a new public API,
a classification of entire fibers, or a computation of their component
counts or generic map degrees. The next focused question is the actual
cubic B_1 and D_2 residual-Levi families in Theorem 3; determine whether
their semisimple critical schemes have a structure extending the exclusion
argument, or supply a genuine representation-theoretic obstruction.
