# Generic smoothness and excess dimension from affine Levi branching

Date: 2026-10-23. Characteristic-zero written deductions, independently
corroborated by exact calculations. The inputs are classical generic
smoothness, the Jacobian criterion, fundamental-character coordinates and
the trace definition of the Springer map. No literature priority is claimed.

## 1. Group, representation and exact hypothesis

Let G be a simple complex simply connected algebraic group, rho an actual
finite-dimensional representation with nondegenerate trace form, and

    B(X,Y)=tr(d rho(X) d rho(Y)),
    B(theta_rho(g),X)=tr(rho(g) d rho(X)).

Simple here allows a finite center. Type B and D simply connected groups
are Spin; the SO defining representations can be pulled back to them.
The present statement concerns characteristic-zero group fibers, with no
parabolic or integral Borel source being changed.

Fix a proper Levi subgroup K of G, write S=[K,K], and put

    ell=rank S,  s=dim Z(K)^0,  r=rank G=ell+s.

The central multiplication cover S x Z(K)^0 -> K is finite etale, S is
simply connected, and its adjoint quotient has coordinates

    (x,z) in A^ell x (C*)^s.

Here x_i are the fundamental characters of all simple factors of S, and
z_j are torus characters. Retain this cover; its finite kernel may identify
different quotient-critical pieces downstairs.

**Affine branching hypothesis.** The restricted character has the form

    chi_rho(x,z)=a(z)+sum_(i=1)^ell x_i f_i(z),                     (1)

where a and f_i are Laurent polynomials in z. This means *jointly*
affine-linear in x: products x_i x_j and powers x_i^2 are excluded.
The hypothesis is not automatic for fundamental representations.

For h in the center of Lie(K), the trace form gives independent linear
target parameters c_1,...,c_s, so that the corrected local potential is

    W_c(x,z)=a(z)+sum x_i f_i(z)-sum c_j log z_j.                   (2)

The logarithms are local primitives; all differentials used below are
global algebraic forms. Independence of the c_j follows from nondegeneracy
of B on the center of Lie(K). Requiring Z_G(h)=K removes root hyperplanes.

Define the **coefficient zero scheme**, with its actual scheme structure,

    Z = V(f_1,...,f_ell) subset T=(C*)^s.                          (3)

## 2. The generic-fiber theorem

**Theorem.** Under (1), there is a Zariski-dense open set of targets
h in the center of Lie(K), with Z_G(h)=K, such that:

1. F_h intersect G_reg is smooth, possibly empty, where
   F_h=theta_rho^(-1)(h) and regular means centralizer dimension r.
2. Every quotient-critical point lies over a smooth point of the scheme Z.
   If its component of Z has codimension q in T, the local dimension is

       dim(F_h intersect G_reg)=dim K-r+ell-q.                    (4)

3. The quotient-critical scheme has local dimension ell-q. The Hessian
   of W_c has rank s+q, so its kernel is exactly the tangent space of the
   critical locus. Equivalently, it is nondegenerate in directions normal
   to that locus (Morse-Bott). The regular fiber is reduced there.

Components of different codimensions can give different local dimensions.
Only components whose target map is dominant contribute generically.
The theorem does not assert nonemptiness or a component count.

Neither reducedness nor smoothness of the entire coefficient scheme Z is
assumed. Its bad locus cannot contribute over generic target parameters.

### 2.1 Exact incidence equations

Write D_j=z_j partial/partial z_j. The critical equations of (2) are

    f_i(z)=0,
    mu_j(x,z):=D_j a(z)+sum_i x_i D_j f_i(z)=c_j.                  (5)

Thus the scheme Crit(W_c) is exactly the fiber over c of

    mu:A^ell x Z -> A^s.                                         (6)

No radical has been taken. At a fixed z the possible target values form
an affine subspace of dimension p=rank(D_j f_i)(z).

### 2.2 Why singular coefficient points disappear generically

At a closed point z of Z let d=dim_z Z and p=rank(D_j f_i)(z).
The embedding dimension of Z at z is s-p. By the Jacobian/regular-local-ring
criterion, if Z is not smooth there, then

    s-p>d.                                                       (7)

Stratify the nonsmooth locus into finitely many locally closed pieces on
which d and p are constant. Such a stratum has dimension at most d.
Its image under (6), allowing arbitrary x, therefore has dimension at most

    d+p<s.                                                       (8)

Indeed the x images at each fixed z are affine spaces of dimension p.
The closure of the union of these finitely many images is a proper closed
subset of the target A^s. This excludes singular points, embedded structure
and generically nonreduced coefficient components from generic target fibers.
It does not discard them for special targets.

### 2.3 Generic smoothness on the remaining incidence scheme

The remaining source A^ell x Z_sm is smooth. Characteristic-zero generic
smoothness in the target gives a dense open set over which (6) is smooth;
images of nondominant components can be excluded. On a component of Z_sm
of codimension q, the dimension of its source in (6) is ell+s-q, hence
its nonempty smooth fibers have dimension ell-q.

Their defining ideal is the gradient ideal of W_c by (5). Therefore their
tangent spaces are the kernels of its Hessian. Smoothness and the computed
dimension imply Hessian rank (ell+s)-(ell-q)=s+q. Symmetry then gives
nondegeneracy on the quotient by the tangent space of the critical locus.

The [general critical-potential theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md)
adjoins dim K-r smooth directions at every regular group-fiber point.
This gives (4). Conversely every quotient-critical point has a regular K
lift, and the trace identity places that lift in F_h. Equivariance then
gives Z_G(g)=Z_K(g), so that lift is G-regular. Finite etale descent preserves
the local smoothness and dimensions. Finally intersect the target open set
with Z_G(h)=K. This proves the theorem.

### 2.4 Local elimination, a second way to read the result

At a smooth point of Z of codimension q, choose q independent f_i as
normal coordinates v_1,...,v_q. Locally the remaining f_i are combinations
of these. An invertible change of the x variables makes (2) equal to

    sum_(i=1)^q x'_i v_i + a(z)-sum c_j log z_j,

independent of ell-q other x variables. Eliminating the q normal coordinates
and their q multipliers leaves the critical equations of

    a|_Z - sum c_j log(z_j|_Z),                                  (9)

together with ell-q free coordinates. For generic c the residual critical
points are isolated and nondegenerate. The incidence proof above supplies
this genericity conclusion and shows why singular coefficient points do
not supply hidden generic components.

## 3. Uniform applications

### 3.1 The spin counterexample is the positive-excess case

For G=Spin_(2n+1), rho=V_(omega_n), and the rank-two semisimple target
stratum, S=Spin_(2n-1), ell=n-1 and s=1. Only the spin coordinate occurs:

    chi=(z+z^(-1))y.

Thus Z is the two reduced points z=+i,-i on the central cover, with
codimension q=1. Formula (4) gives excess ell-q=n-2. The stronger
[explicit spin proof](SPIN_FIBER_EXCESS.md) establishes smoothness for every
nonzero target, not just generic targets, and explains the deck involution.

### 3.2 Type A: every exterior power with one repeated target eigenvalue

Let G=SL_N, rho=exterior^k(C^N), and h semisimple with one eigenvalue of
multiplicity m and q=N-m other, mutually distinct eigenvalues, all different
from the repeated one. For m=1 this simply means a regular semisimple target.
Then K=S(GL_m x (C*)^q), ell=m-1 and s=q.

If m>min(k,N-k), there are no regular fiber points, by the earlier block
bound. Suppose m<=min(k,N-k). Write b_1,...,b_q for external eigenvalues
on the group quotient, e_j their elementary symmetric functions, and
a_i for the first-block characteristic coefficients. The determinant
constraint sets a_m=(product b_j)^(-1). Exterior-power branching gives

    chi=sum_(i=0)^m a_i e_(k-i)(b),  a_0=1.

On the central cover, a_i is a central unit times the i-th fundamental
character of SL_m, for 1<=i<m. Thus the affine hypothesis holds and

    Z: e_(k-m+1)=...=e_(k-1)=0.                                  (10)

These are m-1 independent coefficient coordinates, with indices between
1 and q-1. The symmetric quotient of the external torus is finite, with
coordinate ring C[e_1,...,e_q,e_q^(-1)]. Every component of the inverse
image of (10) has codimension m-1: the coefficient subspace gives the
dimension upper bound, and m-1 equations give the lower bound. Passing
through the additional finite etale central cover changes no dimensions.

**Corollary A.** For generic targets in this stratum, the regular fiber is
empty or smooth of pure dimension

    dim K-r = m(m-1).                                            (11)

There is no excess dimension and every regular germ is reduced. This covers
all N,k,m in the stated range, not just SL_(2k) with m=k. The earlier
middle-exterior theorem gives stronger counts and exceptional target values
in that particular family; those do not follow from generic smoothness.

### 3.3 Type C: every zero-block size for the last fundamental module

Let G=Sp_(2n), rho=V_(omega_n), and let a semisimple target have a zero
eigenspace of dimension 2p and q=n-p distinct nonzero eigenvalue pairs.
Then K=Sp_(2p) x (C*)^q, ell=p and s=q. The case p=0 is the ordinary
generic regular-semisimple-target case. Assume p>=1.

Write A(t) for the reciprocal characteristic polynomial of the first block
and Q(t)=sum b_j t^j for the reciprocal polynomial of the external block.
The primitive character is Coeff_(t^n)((1-t^2)A(t)Q(t)). The coefficients
of A are affine-linear combinations of the fundamental Sp_(2p) characters,
so the affine hypothesis holds. The previous
[all-rank coefficient identity](SYMPLECTIC_TARGET_STRATA.md) gives

    Z: b_j-b_(j-2)=0,  q-p+1<=j<=q,                              (12)

when p<=q; for p>q the coefficient ideal is the unit ideal. Here b_0=1
and b_j=0 for j<0. The p conditions in (12) are independent affine-linear
conditions on b_1,...,b_q. The external torus quotient by signed permutations
is finite onto A^q. As in type A, the dimension bound from this finite map
and the bound from p equations imply every component of Z has codimension p.

**Corollary C.** If p>floor(n/2), no regular fiber points occur. If
p<=floor(n/2), then for generic targets in this stratum the regular fiber
is empty or smooth of pure dimension

    dim K-r = 2p^2.                                              (13)

This proves generic reducedness on the regular locus in every permitted
zero-block size, including the previously unclassified odd-rank boundary
and smaller zero blocks. It does not determine whether every admissible
stratum has a nonempty generic regular fiber, nor count its components or
classify its exceptional target values. Repeated nonzero target pairs add
semisimple Levi factors; joint affinity then requires a separate check.
The earlier even-rank equality theorem remains stronger on its own strata.

## 4. Known inputs, checks, limits and the next gate

Generic smoothness in the target is classical; see
[Vakil, Classes 51-52, Theorem 3.3](https://math.stanford.edu/~vakil/0506-216/216class5152.pdf).
The incidence construction is the usual Lagrange-multiplier/critical-point
construction. When a|_Z is constant and Z is smooth and closed in a torus,
the residual problem (9) is the logarithmic master-function problem in
[Huh, The maximum likelihood degree of a very affine variety, Theorem 1](https://web.math.princeton.edu/~huh/MLDegree.pdf).
That theorem gives generic nondegeneracy and an Euler-characteristic count
in its setting. No Euler-characteristic count is inferred here for arbitrary
nonconstant a, singular Z, or a smooth locus not closed in the original torus.

The contributions written here are the application to Springer fibers,
the exclusion of nonsmooth coefficient points by (8), the dimension rule,
and the A/C branching deductions. These follow from known tools; this note
does not establish that their particular formulation is unpublished.

[Independent verification](affine-branching-fibers-verification.json) checks
the type A coefficient identities through N=12 and the type C ideals through
n=13. Original exterior weights and separately generated Weyl characters
check torus trace moments and Hessians, plus all root-direction ranks, at
regular semisimple source points. The target centralizer strata are checked
explicitly. These computations do not use the coefficient equations to
compute the images or differentials. Symbolic controls include singular
and nonreduced coefficient schemes and a nonlinear potential with persistent
double structure. The all-rank/generic conclusions rest on the written proof.

Affinity is essential to this proof. For example W=x^3+z-c log z has
critical ideal (x^2,z-c) for every c!=0, so it stays nonreduced. This is
an abstract polynomial control, not a demonstrated fundamental character
of a simple group. Universal generic smoothness for all fundamental
representations therefore remains unresolved here. Components wholly in
G minus G_reg and integral images are also outside the theorem.

**Next gate:** determine which fundamental-representation/Levi pairs satisfy
joint affinity, and isolate the first unavoidable nonlinear branching terms.
The goal is to extend the criterion with a precise hypothesis or find a
genuine representation-theoretic counterexample to generic smoothness.
Counting the A/C residual critical points is a separate, now well-defined
problem on coefficient varieties; repeated finite checks alone would not
settle either goal.
