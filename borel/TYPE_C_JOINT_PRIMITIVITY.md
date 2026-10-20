# Integral joint images and primitive sums in type C

Status: a complete geometric deduction under the normalization below, using
classical integral cohomology and Leray–Hirsch. The argument and its hypotheses
are audited at the end of this note. This is an unpublished repository proof;
no priority or originality claim is made. The computations corroborate the
argument but are not used to infer its validity in arbitrary rank.

## Statement and conventions

Let a,b,k be positive integers, n=a+b+k, and let

$$X=\operatorname{Sp}(n)/(U(a)\times U(b)\times\operatorname{Sp}(k)).$$

Here Sp(r) denotes the compact quaternionic unitary group. Thus X is the
compact realization of the complex symplectic flag variety of isotropic
subspaces U contained in V of dimensions a and a+b in C^(2n). Its omitted
nodes in the library are (a,a+b). Split the Cartesian coordinates into
head A, middle M, and tail T, with sizes a,b,k.

Use precisely the normalized polynomial generators

$$c_j=e_j(M),\quad 1\le j\le b,\qquad
q_i=e_i(T^2),\quad 1\le i\le k.$$

The library calls q_i `e{i}`. The letter q avoids confusing these classes
with the elementary symmetric function operation. Their polynomial degrees
are j and 2i; degree d means cohomology H^(2d). All cohomology and all
polynomial spans in this note have integer coefficients.

Inside H^(2*)(X;Z), define

$$L=\beta(\mathbb Z[c_1,\ldots,c_b]),\qquad
M_0=\beta(\mathbb Z[q_1,\ldots,q_k]),\qquad
J=\beta(\mathbb Z[c_1,\ldots,c_b,q_1,\ldots,q_k]).$$

For a subgroup N of a free abelian group H, **primitive** (or **saturated**)
means H/N is torsion-free, equivalently N=(N tensor Q) intersect H.

**Theorem.** For every a,b,k>=1 and every d>=0:

1. J_d is primitive in H^(2d)(X;Z). All its nonzero Smith factors equal 1.
2. L_d+(M_0)_d is primitive in J_d, and also in H^(2d)(X;Z).
3. L intersect M_0=Z times the identity: the degreewise intersection is zero
   for d>0 and is Z for d=0.

There is a single integral basis of J in which both L and M_0 are spans of
subsets. This stronger assertion proves the second and third statements.
It also explains why mixed products can contribute new free directions
without creating torsion in J/(L+M_0).

## Published inputs

The following standard results are used with integral coefficients.

- [Borel, Sur la cohomologie des espaces fibres principaux et des espaces homogenes](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Borel-Sur.pdf),
  Proposition 31.1 and the following numbered case (31.2), journal p.202,
  give the integral block-symmetric presentations for complex and
  quaternionic partial flags. The latter is case (31.2), not Proposition
  31.2, which concerns a different homogeneous space.
- [Pragacz–Ratajski](https://arxiv.org/pdf/alg-geom/9602019), Theorem 2.1(i),
  gives the Lagrangian Schubert representatives Q-tilde_lambda in the Chern
  classes of the dual tautological bundle, indexed by strict partitions with
  parts at most r. Appendix A, §§10.1–10.3, identifies quaternionic flag
  cohomology with the corresponding complex flag ring with degrees doubled,
  and describes quaternionic Grassmannian generators. Their §§10.4–10.5
  discuss the associated Lagrangian fibrations.
- [Brion, Lectures on the geometry of flag varieties](https://arxiv.org/pdf/math/0410240),
  §1.2, Example 1.2.3(6), and §1.3, supplies the integral Schubert basis and
  its compatibility with parabolic projections. Using opposite Schubert
  classes, pullback preserves the minimal coset label and codimension.
- The integral Leray–Hirsch theorem: if global classes restrict to an integral
  basis of a fiber, multiplication by those classes identifies total-space
  cohomology with a direct sum of shifted copies of base cohomology. We use
  actual global Chern-polynomial lifts and include the unit among them.

The polynomial-to-Chern-class interpretation agrees with the characteristic
map in [Kumar–Rogers, Part II](https://kumar.math.unc.edu/papers/kumar81.pdf),
Proposition 11 and the normalized Z-form in the proof of Theorem 16. The
present theorem concerns joint products from two selected Levi factors;
it is not quoted as a theorem stated in that paper.

## 1. An integral quaternionic base

Enlarge both unitary blocks to quaternionic ones and put

$$Q=\operatorname{Sp}(n)/(\operatorname{Sp}(a)\times
\operatorname{Sp}(b)\times\operatorname{Sp}(k)).$$

This space parametrizes ordered mutually orthogonal quaternionic blocks
(A,M,T). Write p_i^A,p_j^M,p_l^T for their symplectic characteristic
classes, normalized to restrict to elementary functions of squared roots.
In particular p_l^T=q_l. Its integral ring is

$$R=H^*(Q;\mathbb Z)
 =\frac{\mathbb Z[p_1^A,\ldots,p_a^A,p_1^M,\ldots,p_b^M,
                    p_1^T,\ldots,p_k^T]}
 {(\operatorname{coeff}_{u^r}(P_A(u)P_M(u)P_T(u)-1):1\le r\le n)},\tag{1}$$

where P_A(u)=1+sum p_i^A u^i, and similarly for M,T. Each p_i has
cohomological degree 4i. These classes are signed even Chern classes of
the underlying complex symplectic bundles; the signs are units.

This is an integral presentation, not an invariant-ring calculation over Q.
One may also justify it by projecting Q to HGr(k,n): the fiber is
HGr(b,n-k). Schur polynomials in the middle characteristic classes lift
the fiber's integral Schubert basis. Leray–Hirsch gives generation and
free rank n!/(a!b!k!). The Whitney relation gives the surjection from (1);
the identical complex partial-flag presentation is free of the same rank,
so this surjection is an isomorphism. In particular R is free abelian.

Reorder the blocks as **M,A,T**. The degree-doubling identification sends
R to the integral cohomology of

$$F=\operatorname{Fl}(b,n-k;n;\mathbb C).$$

It identifies the block characteristic classes and therefore the two
subrings

$$R_M=\mathbb Z[p_1^M,\ldots,p_b^M]\subset R,\qquad
R_T=\mathbb Z[p_1^T,\ldots,p_k^T]\subset R$$

with the two standard Grassmannian pullback images in F. For the second
projection, the rank-k quotient bundle generates the same ring as the
rank-(n-k) tautological bundle, by their total-Chern inverse relation.
This observation is essential: in the original A,M,T order the middle
block is a quotient of successive flag steps, not the first tautological
subspace for a standard forgetting map.

Let

$$D=\{b,n-k\},\qquad
W^D=\{w\in S_n:\operatorname{Des}_R(w)\subset D\}.$$

Transfer the opposite Schubert basis of F to a basis {sigma_w:w in W^D}
of R. Its polynomial degree in the type C grading is 2 length(w). Then

$$R_M=\bigoplus_{\operatorname{Des}_R(w)\subset\{b\}}\mathbb Z\sigma_w,
\qquad
R_T=\bigoplus_{\operatorname{Des}_R(w)\subset\{n-k\}}\mathbb Z\sigma_w.\tag{2}$$

Both are spans of subsets of the same integral basis. Their sum is primitive
in R. Since a>0 gives b<n-k, the only common index has no descents, hence
is the identity permutation. Consequently R_M intersect R_T=Z.

## 2. Restore the unsquared middle block

Set

$$Y=\operatorname{Sp}(n)/(\operatorname{Sp}(a)\times U(b)\times
                         \operatorname{Sp}(k)).$$

The projection rho:Y->Q has fiber Sp(b)/U(b)=LG(b,2b). Let E_M be the
rank-b complex tautological middle bundle on Y, with the dual convention
chosen so c_j(E_M^*) corresponds to c_j=e_j(M). For each strict partition
lambda with parts at most b, let

$$u_\lambda=\widetilde Q_\lambda(c_1(E_M^*),\ldots,c_b(E_M^*)),
\qquad u_\varnothing=1.$$

These are global integral classes restricting to the integral Schubert
basis of each Lagrangian fiber. Integral Leray–Hirsch gives

$$H^*(Y;\mathbb Z)=\bigoplus_{\lambda\text{ strict},\,\lambda_1\le b}
                         u_\lambda\,\rho^*R.\tag{3}$$

Equation (3) is a direct-sum statement for graded abelian groups/modules;
it is not a claim that the cohomology ring is a tensor product of the
cohomology rings of base and fiber.

The squared middle classes p_i^M are integer polynomials in c_1,...,c_b.
Indeed, for C(v)=1+sum c_j v^j,

$$C(v)C(-v)=\sum_{i=0}^b(-1)^ip_i^M v^{2i}.\tag{4}$$

The relation P_AP_MP_T=1 expresses every p_i^A integrally in p^M,p^T.
Thus rho^*R is contained in Z[c_1,...,c_b,q_1,...,q_k]. Equation (3)
then proves that these generators generate **all** of H^*(Y;Z).

Now consider the projection

$$g:Y\longrightarrow Z_M:=\operatorname{Sp}(n)/(U(b)\times
                                      \operatorname{Sp}(n-b)), $$

which remembers only the middle complex isotropic b-plane. It fits into
the commuting square with rho:Y->Q, the projection Q->HGr(b,n), and
Z_M->HGr(b,n). The two vertical fibers in this square are LG(b,2b),
and their Chern-polynomial lifts u_lambda agree under g^*.

Applying Leray–Hirsch to Z_M->HGr(b,n) shows that H^*(Z_M;Z) is generated
by its middle Chern classes: the base classes p^M are polynomials (4),
and the fiber lifts u_lambda are also such polynomials. Therefore its
pullback is exactly the middle-only algebra, with the **same** lifts as (3):

$$L_Y:=\mathbb Z[c_1,\ldots,c_b]
      =\bigoplus_\lambda u_\lambda\rho^*R_M\ \subset H^*(Y).\tag{5}$$

The pullback is injective: the corresponding base pullback R_M->R is
injective by (2), and (3) and its analogue for Z_M use identical fiber
bases. The tail-only algebra is

$$M_Y:=\mathbb Z[q_1,\ldots,q_k]=u_\varnothing\rho^*R_T.\tag{6}$$

Consequently the classes

$$\mathcal B=\{u_\lambda\rho^*\sigma_w:
                 \lambda\text{ strict},\lambda_1\le b,\ w\in W^D\}\tag{7}$$

form an integral basis of H^*(Y), with both L_Y and M_Y spanned by subsets:

- L_Y uses every lambda and the indices with descents contained in {b}.
- M_Y uses lambda=empty and descents contained in {n-k}.

Their sum is therefore primitive, and their intersection is Z. More
explicitly there is a graded abelian-group decomposition

$$\frac{H^*(Y)}{L_Y+M_Y}\cong
\frac{R}{R_M+R_T}\ \oplus\!
\bigoplus_{\lambda\ne\varnothing}u_\lambda\frac{R}{R_M}.\tag{8}$$

Every summand on the right is free by (2). No assertion that sums of
arbitrary primitive sublattices are primitive is being used.

## 3. Restore the head and identify the actual image

The projection pi:X->Y has fiber LG(a,2a). The rank-a head bundle on X
supplies global integral Q-tilde representatives v_mu of its Schubert
basis, including v_empty=1. Leray–Hirsch gives

$$H^*(X;\mathbb Z)=\bigoplus_\mu v_\mu\,\pi^*H^*(Y;\mathbb Z).\tag{9}$$

Since H^*(Y) is free by (3), the unit summand pi^*H^*(Y) is primitive
in H^*(X). The characteristic-class identification and the generation
proved after (4) show that this summand is precisely J. This proves
the first claim of the theorem. Equations (5)–(8) identify L and M_0
inside J and prove the remaining claims. Finally, a subgroup primitive
in J is primitive in H^*(X) because both successive quotients are free.
All decompositions preserve degree, so every assertion holds separately
in H^(2d). This completes the proof.

## Integral kernel as a corollary

The earlier [rational joint presentation](experiments/TWO_STEP_DERIVATION.md)
also holds over Z for these normalized inputs:

$$J\cong\frac{\mathbb Z[c_1,\ldots,c_b,q_1,\ldots,q_k]}
 {(h_{a+1}(M^2\cup T^2),\ldots,h_n(M^2\cup T^2))}.\tag{10}$$

To check that this adds no hidden torsion or extra relations, put
D_b=Z[c_1,...,c_b] and A_b=Z[p_1^M,...,p_b^M], with the map (4).
The integral Q-tilde identities give

$$D_b=\bigoplus_{\lambda\text{ strict},\lambda_1\le b}A_bu_\lambda.$$

Indeed Q-tilde_I for all partitions with parts at most b is an integral
basis of D_b; removing a pair of equal parts i multiplies it by p_i^M.
Unique pair-removal leaves a strict partition. These are Proposition 4.7
and the factorization in Propositions 4.2–4.3 of Pragacz–Ratajski.

The natural map R tensor_(A_b) D_b -> H^*(Y) takes this free module basis
to the Leray–Hirsch basis (3), so it is an integral ring isomorphism.
Eliminate p_i^A from its presentation (1): the coefficient equations give
p_i^A=(-1)^i h_i(M^2 union T^2), followed by the relations in (10).
The elimination is triangular with unit coefficients. This proves (10)
integrally, without averaging over a Weyl group or dividing by its order.

## Closed rank formulas

Let [n choose r]_q be the Gaussian binomial and set

$$F_b(t)=\prod_{j=1}^b(1+t^j),\qquad
Q_{a,b,k}(t)={n\brack b}_{t^2}{n-b\brack k}_{t^2}.$$

The basis (7) gives the following polynomials, with t recording polynomial
degree, not cohomological degree:

$$P_J(t)=F_b(t)Q_{a,b,k}(t),\qquad
P_L(t)=F_b(t){n\brack b}_{t^2},\qquad
P_{M_0}(t)={n\brack k}_{t^2}.\tag{11}$$

In each degree, the Smith list of J consists of [t^d]P_J ones. The
number of unit factors for the separate-image sum is the coefficient
of t^d in P_L+P_(M_0)-1, and

$$J_d/(L_d+(M_0)_d)\cong\mathbb Z^{r_d},\qquad
r_d=\operatorname{coeff}_{t^d}(P_J-P_L-P_{M_0}+1).\tag{12}$$

The constant corrections account for the shared identity. Formula (11)
is identical to the earlier rational Hilbert series. The new integral
argument establishes that the ranks count primitive lattice bases.

For a=b=k=1 the quotient polynomial is t^3+t^5+t^6+t^7. In particular,
in degree 3 the new direction is represented by beta(c_1 q_1).

## Audit, limitations and verification

| Potential gap | Resolution in the proof |
| --- | --- |
| Rational equality mistaken for integral equality | Global integral fiber bases and integral quaternionic cohomology are used throughout. |
| Joint generators might miss some integral classes of Y | Equations (3) and (4), plus elimination of p^A, prove integral generation. |
| Middle algebra might be smaller than the asserted summands | The commuting square with Z_M identifies the same Chern-polynomial fiber lifts and proves (5). |
| Separate primitive sublattices might have a nonprimitive sum | Equations (2) and (7) give a common integral basis, a stronger property. |
| Middle block is not a standard forgetting map in the original ordering | Reorder the quaternionic blocks as M,A,T before identifying the two Grassmannian images. |
| Two basis subsets might share positive-degree indices | a>0 makes the two permitted descent positions distinct; their common permutation is the identity. |
| Leray–Hirsch mistaken for a ring tensor product | Only its module decomposition is used; (10) has a separate base-change argument retaining characteristic-class relations. |
| Scaling changes the lattice | The theorem uses c_j=e_j(M), q_i=e_i(T^2). Scaled character forms require separate analysis and need not be saturated. |
| Orthogonal examples might contradict the proof | The integral Lagrangian Chern-polynomial basis is specific to type C. No B/D joint primitivity theorem is asserted here. |

The [verification script](experiments/type_c_joint_primitivity.py) checks
the rank formulas against the previous type C tables, tests additional
block sizes, and constructs the basis (7) from ordinary type A Schubert
polynomials in the reordered squared roots and integral Q-tilde Pfaffians.
It checks the actual type C Schubert lattices, both basis subsets, all
containments, and the Smith factors. C3 basis columns also have an
independent Cartesian BGG check. Exact ranges and source hashes are in
[the verification report](experiments/type_c_joint_primitivity_verification.json).

The published ingredients and the combination established here are separated
above. A systematic search for this exact joint-image formulation, or a
claim that it is new, is outside the scope of this proof check.
