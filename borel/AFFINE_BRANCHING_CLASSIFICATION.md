# Which fundamental representations have affine Levi branching?

Date: 2026-10-23. Written all-rank classical proofs and a complete finite,
certificate-assisted exceptional classification. All groups and characters
are over C. The proof uses standard highest-weight theory, classical
exterior/spin models and weight-polytope convexity. No novelty or exhaustive
literature-priority claim is made.

This determines the exact scope, for fundamental representations, of the
[affine-branching generic-fiber theorem](AFFINE_BRANCHING_FIBERS.md).
It does not classify singular Springer fibers: failure of its sufficient
hypothesis is not a proof that the fiber is singular or nonreduced.

## 1. Definition and general criteria

Let G be a simple complex simply connected algebraic group, with maximal
torus T, simple roots indexed by 1,...,r, and fundamental weights omega_k.
Let V=V(omega_k). For a proper node subset I, let K=K_I be the standard
Levi and S=[K,K]. The empty subset gives K=T. Work on the finite central
cover S x Z(K)^0 -> K. The group S is simply connected; its fundamental
characters x_i, i in I, freely generate C[S]^S as a polynomial algebra.
Choose character coordinates z on the central torus. Say that (V,K) is
**jointly affine** if

    chi_V|_K = a(z) + sum_(i in I) x_i f_i(z),                    (1)

with Laurent-polynomial coefficients. This excludes products between
different simple factors as well as powers within a factor. The definition
uses fundamental characters, not arbitrary polynomial quotient coordinates.
There is no change to an integral Borel source or to a cohomology map here.

**Character criterion.** For any finite-dimensional actual representation V,
(1) holds if and only if every irreducible S-constituent is either trivial,
or fundamental on exactly one simple factor and trivial on all others.
Equivalently, every constituent highest weight has nonnegative Dynkin
coordinates whose sum is at most one.

Proof. Decompose V into S-irreducibles with central-torus multiplicity
spaces. Irreducible S-characters are linearly independent over the central
Laurent ring. Thus (1) holds exactly for the stated constituents. Forgetting
the central grading does not alter the test: multiplicities are nonnegative,
so forbidden constituents cannot cancel. This also explains why a tensor
product of fundamentals on two factors is forbidden. QED.

Two criteria make the exceptional calculation a finite proof with small
certificates. They apply to any irreducible V(lambda), not only fundamentals.
All norms below use a positive definite Weyl-invariant form on the real
root space, with the Levi inheriting that same form.

**Extremal obstruction.** Suppose mu=w lambda for w in the Weyl group,
mu_i=<mu,alpha_i^vee> >= 0 for i in I, and sum_(i in I) mu_i > 1.
Then branching is not jointly affine.

Proof. The extremal weight mu occurs with multiplicity one. Every weight
of V(lambda) lies in conv(W lambda), hence has norm at most ||lambda||.
For i in I, ||mu+alpha_i||^2 > ||mu||^2 because (mu,alpha_i)>=0.
Consequently no raising operator for S can move its weight vector: it is
an S-highest vector. Its highest weight violates the character criterion.
QED. A certificate consists just of a Weyl word and its resulting weight;
one word may obstruct many node subsets.

**Norm-gap criterion.** Let varpi_i be the fundamental weights of S,
realized in its root span, and for nonempty I define

    M_I = min_(i,j in I, i<=j) ||varpi_i+varpi_j||^2.

If ||lambda||^2 < M_I, then branching is jointly affine.

Proof. Restriction of a weight to S is its orthogonal projection to the
Levi root span. Thus every constituent highest weight nu has squared norm
at most ||lambda||^2. Fundamental-weight Gram entries are nonnegative:
positive within an irreducible component and zero between components.
If nu is neither zero nor a single fundamental weight, its nonnegative
integer expansion contains some varpi_i+varpi_j, allowing i=j. The remaining
terms can only increase squared norm. Therefore ||nu||^2>=M_I, a
contradiction. QED.

The norm gap is a sufficient test, not a proposed necessary condition for
arbitrary representations. Remarkably, for the exceptional fundamental
representations the norm-gap and extremal tests together decide every pair.
The empty subset is always affine, independently of a norm test.

## 2. Complete classical classification

For A_(N-1), write G=SL_N and the Levi block sizes as b_1,...,b_s, with
sum b_i=N. Its derived group is product SL_(b_i). Put

    t = number of blocks b_i >= 2.

For B_n, C_n and D_n, use Levi block data (b_1,...,b_s;m), where
sum b_i+m=n: the derived group has factors SL_(b_i), together with the
residual group Spin_(2m+1), Sp_(2m), or Spin_(2m), respectively.
These descriptions are on the central cover; determinant half-powers
are torus characters there, not an assertion about descent to SO.
For D, absorb a residual Spin_2 torus into a size-one GL block, so m=0
or m>=2. Spin_4 is A_1 x A_1 and must be treated as two simple factors.
The two possible fork embeddings may interchange half-spin chirality.

**Theorem (all ranks).** The following conditions are necessary and sufficient
for joint affinity for every proper standard Levi. The ordinary rank-one
case SL_2 is immediate. Use B_n,C_n at n>=2 and D_n at n>=4.

| Group and highest weight | Exact condition |
| --- | --- |
| SL_N, omega_k, 1<=k<N | k=1 or k=N-1 or t<=1 |
| Sp_(2n), omega_k, 1<=k<=n | k=1 or t=0 |
| Spin_(2n+1), omega_k, 1<=k<n | (k=1 or t=0) and (m=0 or m>k) |
| Spin_(2n+1), spin omega_n | t+1_(m>0)<=1 |
| Spin_(2n), omega_k, 1<=k<=n-2 | (k=1 or t=0) and (m=0 or m>=k+2) |
| Spin_(2n), either half-spin, m>=2 | t=0 |
| Spin_(2n), either half-spin, m=0 | t<=1, or the D_4 exception described below |

### 2.1 Exterior powers in type A

The standard exterior decomposition is

    exterior^k(direct_sum E_i)
      = direct_sum_(sum a_i=k) tensor_i exterior^(a_i) E_i,
        0<=a_i<=b_i.

The i-th factor is derived-trivial for a_i=0 or b_i, and fundamental
otherwise. Joint affinity is therefore equivalent to the impossibility
of two partial blocks. If k=1 or N-1, two partial blocks are impossible.
If t<=1 they are also impossible. Conversely, choose two blocks of size
at least two. Requiring both allocations to be interior gives all total
sums from 2 through N-2 when the other allocations vary freely. Thus
every 2<=k<=N-2 has a forbidden constituent. This proves the first row.

### 2.2 Primitive symplectic exterior powers

The defining module restricts as a direct sum of E_i, E_i^*, and the
residual symplectic defining module, so k=1 is affine.
For k>=2, if b_i>=2, take an extremal weight of omega_k with one +1
and one -1 in that GL block and k-2 other occupied positions. The
Weyl orbit consists of signed permutations of (1^k,0^(n-k)). Sort within
each Levi factor to make it dominant. The chosen SL block still has
largest coordinate 1 and smallest -1, so its Dynkin-coordinate sum is
two. The extremal obstruction proves necessity of t=0.

If t=0 the only possible derived factor is Sp_(2m). Every weight of the
primitive exterior module has coordinates in {-1,0,1}. A dominant weight
for Sp_(2m) is therefore (1^j,0^(m-j)), which is fundamental or zero.
This proves sufficiency, including m=0.

### 2.3 Odd orthogonal exterior powers

For k<n, the fundamental module is exterior^k of the vector module.
The same extremal GL-block argument excludes t>0 when k>=2. If
0<m<=k, choose an extremal weight having all m residual coordinates
equal to 1 and place the remaining k-m occupied coordinates outside.
It is Levi-dominant after sorting. The residual highest weight is
(1,...,1)=2 omega_m, including the B_1 case. This is forbidden.

Conversely, with t=0 and m>k, every residual dominant weight is
(1^j,0^(m-j)) with j<=k<m, hence fundamental or zero. If m=0 the
derived group is trivial. For k=1 arbitrary SL blocks are allowed:
the defining-module direct sum gives only one-factor fundamentals
unless the residual block is B_1, which was already excluded.

### 2.4 Even orthogonal exterior powers

For k<=n-2 the same GL obstruction applies. There are at least two
zero positions, so the even-sign requirement of the D Weyl group does
not restrict the signed extremal weights used here. If 2<=m<=k+1,
choose exactly m-1 residual coordinates equal to 1, the last zero,
and k-m+1 occupied coordinates outside. The resulting residual highest
weight is omega_(m-1)+omega_m. This is forbidden; for m=2 it is a
tensor product of the two A_1 fundamentals, not a single fundamental.

If m>=k+2 and t=0, each dominant residual weight has j<=k<=m-2
ones and otherwise zero, so it is fundamental or zero. Signs can be
removed because at least two residual positions are zero. The vector
case with arbitrary SL blocks again follows from the direct sum.
This proves the exterior D row, including the m=0 alternative.

### 2.5 Spin, half-spin and the D_4 exception

On the central cover the odd spin module restricts as the tensor product
of the GL exterior algebras and, when m>0, the residual odd spin module,
with central determinant twists. Each interior GL exterior degree is
a fundamental constituent of its SL factor; the residual spin is always
fundamental. Hence two or more nontrivial derived factors produce a
forbidden tensor constituent, and at most one gives only allowed ones.
This is exactly t+1_(m>0)<=1.

For even half-spin with m>=2, exterior-degree parity chooses one of the
two residual half-spin modules. Either is fundamental on one residual
simple factor (also for D_2). Any nontrivial SL block produces a forbidden
tensor constituent with that factor. With no such block all constituents
are allowed, proving the t=0 condition.

For m=0 the two half-spin modules are the even and odd parts of the
tensor product of GL exterior algebras, again with central twists.
At most one nontrivial SL block is harmless. With two or more such
blocks, a block of size >=3 permits both parities while retaining two
interior allocations. A size-one block also permits a parity change.
If all blocks have size two and there are at least three, choose either
two or three interior allocations to realize both parities. The sole
exception is exactly two blocks of size two: two interior allocations
have even parity, whereas an odd allocation has at most one active factor.

Thus for D_4 with GL_2 x GL_2 Levi and no residual D factor, the odd
half-spin is affine and the even half-spin is not. In standard numbering,
I={1,3} allows omega_3 and excludes omega_4; I={1,4} interchanges them.
The different Levi I={3,4} has residual D_2 and no nontrivial GL block;
both half-spins are affine there. This completes the proof in all ranks.

## 3. Complete exceptional classification

Use Bourbaki/Sage node numbering. In E_n the chain is
1--3--4--5--6--7--8, truncated at n, with node 2 attached to node 4.
In F_4 the chain is 1--2=>3--4, with 1,2 long and 3,4 short.
In G_2 node 1 is short and node 2 long. A connected subset means its
induced Dynkin diagram is connected. All subsets in the table are proper;
the empty subset is always allowed and is omitted from the middle column.

**Theorem (finite exceptional classification).** The following table gives
all and only the affine pairs. Counts include the empty subset and count
standard node subsets before identifying conjugate Levi subgroups.

| Group and fundamental node | Allowed nonempty I | Affine pairs / all proper Levis |
| --- | --- | --- |
| G_2, 1 (dimension 7) | {2} | 2/3 |
| G_2, 2 (adjoint 14) | none | 1/3 |
| F_4, 1 (adjoint 52) | {1,2,3} | 2/15 |
| F_4, 2 or 3 | none | 1/15 each |
| F_4, 4 (dimension 26) | {1}, {2}, {1,2}, {2,3}, {1,2,3}, {2,3,4} | 7/15 |
| E_6, 1 or 6 (dimension 27) | every connected I | 25/63 each |
| E_6, 2 (adjoint 78) | connected D_4 or D_5 | 4/63 |
| E_6, 3, 4 or 5 | none | 1/63 each |
| E_7, 1 (adjoint 133) | connected D_4, D_5, D_6 or E_6 | 6/127 |
| E_7, 2 (dimension 912) | {1,2,3,4,5,6} | 2/127 |
| E_7, 3, 4, 5 or 6 | none | 1/127 each |
| E_7, 7 (dimension 56) | every connected I | 34/127 |
| E_8, 1 (dimension 3875) | {1,2,3,4,5,6,7} | 2/255 |
| E_8, 2, 3, 4, 5, 6 or 7 | none | 1/255 each |
| E_8, 8 (adjoint 248) | connected D_4, D_5, D_6, D_7, E_6 or E_7 | 8/255 |

Proof with certificates. For each allowed nonempty subset the verification
file records the ambient highest-weight norm, the inherited Levi
fundamental-weight Gram matrix, and a strict norm gap. For every excluded
subset it supplies a Weyl word producing a Levi-dominant extremal weight
with coordinate sum greater than one. The two proved criteria establish
the answers. The independent auditor reconstructs the Cartan matrices
from the diagrams above and verifies that every proper subset is covered
exactly once for each fundamental node. Thus this is exhaustive within
the stated finite scope; it does not rely on extrapolation from a sample.

There are 3,373 pairs: 133 affine and 3,240 nonaffine. The 27 empty-subset
cases are automatic, and the other 106 positive cases have strict norm
certificates. It is unnecessary to construct huge characters. For example,
start with E_8 omega_4, of dimension 6,899,079,264, and apply simple
reflections successively at 4,5,6,7,8. The resulting ambient Dynkin labels are

    (0,1,1,0,0,0,0,-1).

For the E_7 Levi on nodes 1,...,7 this is a highest-weight obstruction
omega_2+omega_3. The five-reflection certificate replaces a billions-of-
dimensions character computation for this pair.

The minuscule positive cases have another familiar explanation. A
minuscule weight restricts with pairing in {-1,0,1} against every Levi
coroot. On a simple Levi factor a nonzero dominant constituent must
therefore be a minuscule fundamental weight: pairing with the highest
coroot excludes two positive Dynkin coefficients or a coefficient above
one. For E_6's 27 and E_7's 56, the certificates also show that every
disconnected nonempty proper subset fails. That last assertion is specific
to these pairs, not a general minuscule theorem.

## 4. Independent checks and reproducibility

The [Sage verifier](tools/verify_affine_branching_classification.py) and
[canonical certificate file](affine-branching-classification-verification.json)
provide the following checks:

- 2,453 classical full-weight restriction comparisons through rank six,
  independent of the all-rank block formulas. Restriction aggregates
  weights and subtracts irreducible characters, rejecting any forbidden
  highest weight. The finite comparisons support, rather than replace,
  the written all-rank proof.
- All 3,373 exceptional pairs decided by the two certificate criteria.
- A separate full-weight restriction check of 1,462 exceptional pairs:
  all G_2, F_4 and E_6 fundamentals; E_7 nodes 1,2,6,7; E_8 nodes 1,8.
- Three checks against Sage's built-in Levi branching: F_4 omega_1 to B_3,
  E_7 omega_2 to E_6, and E_8 omega_1 to E_7.

The [Sage-free auditor](tools/audit_affine_branching_certificates.py) uses
only Python rational arithmetic. It independently specifies Cartan
matrices, replays every extremal reflection word, recomputes every norm
gap, checks coverage of every proper subset, and recomputes Weyl dimensions
and the dimension sums of positive restrictions. Its
[audit report](affine-branching-classification-audit.json) binds the exact
sources by SHA256. Dimension sums alone do not establish branching
identities; the norm certificates independently prove affinity, and full
weight restrictions check the listed decompositions.

Run from borel under Sage, then ordinary Python:

```sh
PYTHONPATH=src:tools sage -python tools/verify_affine_branching_classification.py
python tools/audit_affine_branching_certificates.py
```

The native branching identities, omitting central characters, are

    F4 omega1 | B3 = 1 + 2 V(omega1) + V(omega2) + 2 V(omega3),
    E7 omega2 | E6 = V(omega1) + 2 V(omega2) + V(omega3)
                      + V(omega5) + V(omega6),
    E8 omega1 | E7 = 1 + 3 V(omega1) + 2 V(omega2)
                      + V(omega6) + 2 V(omega7).

These are checks of classical branching identities, not newly constructed
representations or proposed new branching laws.

## 5. Published inputs, deductions and the next research gate

The fundamental-character polynomial presentation comes from standard
simply connected quotient theory; see
[Steinberg, Regular elements of semisimple algebraic groups](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf).
The classical exterior, primitive exterior and spin models are standard;
see [Etingof's lecture on classical fundamental representations](https://ocw.mit.edu/courses/18-755-lie-groups-and-lie-algebras-ii-spring-2024/mit18_755_s24_lec05.pdf).
The [Sage branching documentation](https://doc.sagemath.org/html/en/thematic_tutorials/lie/branching_rules.html)
describes the independent native branching implementation. Weight-polytope
convexity and highest-weight theory are classical inputs.

The character criterion and the two certificates are elementary deductions
from those inputs. The classical classification is a written deduction;
the exceptional table is a complete finite proof using those criteria and
audited rational certificates. We have not established whether precisely
this classification has already been tabulated. None of these validity
checks is a novelty certificate.

For a pair in the table, the earlier affine theorem applies to the
Springer map defined by the representation's trace form. For a generic
semisimple target h with centralizer K, its regular fiber is empty or
smooth. A contributing coefficient-zero component of codimension q gives

    dim(theta_rho^(-1)(h) intersect G_reg)
      = dim K - rank G + rank S - q.

The classification determines applicability, not q or dominance of the
coefficient incidence map. Actual central characters must still be retained
to compute those data. It also does not address wholly nonregular fiber
components, special target parameters, or integral Borel-image lattices.

**Next structural gate:** extend the smoothness criterion to restricted
characters with quadratic dependence on the semisimple quotient variables.
The first target is a theorem using the rank and degeneration of the
quadratic coefficient matrix over the central torus, with a precise test
for generic Morse-Bott behavior or a representation-theoretic counterexample.
That would enlarge the geometric theorem beyond the now-classified affine
cases, rather than merely add more examples to this completed table.
