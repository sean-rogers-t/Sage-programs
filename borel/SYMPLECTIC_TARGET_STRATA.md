# Primitive symplectic Springer maps: target strata and reduced components

Date: 2026-10-23. Characteristic-zero proofs and separate exact checks.
These are deductions from classical characters, trace projection and regular
adjoint quotients. Literature novelty is not established. All component
counts below refer to Sp, not automatically to its central quotient.

## 1. Setup and an all-rank obstruction

Let G=Sp_(2n)(C), rho=V_(omega_n), the last fundamental representation.
It is the primitive middle exterior power, the kernel of contraction from
exterior^n of the defining space to exterior^(n-2). Its character and trace
normalization are

    chi(g)=Coeff_(z^n) ((1-z^2) det(1+zg)),
    I_n=binomial(2n-2,n-1)-binomial(2n-2,n-3),
    I_n tr(theta(g) X)=tr(rho(g) d rho(X)).                         (1)

Negative lower binomial entries give zero. Use the simply connected group,
its actual representation and full weight lattice. These are complex
group-map assertions, not integral Borel-image or cohomology assertions.
A regular group element has centralizer dimension n; it need not be
semisimple. A fiber always means the full scheme-theoretic fiber.

Fix a semisimple target h with dim ker(h)=2p. Set q=n-p. Choosing one
representative lambda_j from each distinct nonzero eigenvalue pair,
with multiplicity r_j, gives

    K=Z_G(h)=Sp_(2p) x product_j GL_(r_j),     sum r_j=q.           (2)

The lambda_j are nonzero and distinct up to sign. Their choice of signs
does not affect any count. The general critical-potential theorem places
the entire fiber in K. On its regular locus the fiber is the critical
scheme of the restricted character minus the target logarithms on K//K.

**Theorem 1.** If theta^(-1)(h) has a regular point, then

    p <= floor(n/2).                                              (3)

More precisely, let A be its Sp_(2p) block and let Q(z)=det(1+zB) be the
reciprocal polynomial of its entire remaining symplectic block. Put

    Q(z)=sum_(j=0)^(2q) b_j z^j,  b_j=b_(2q-j), b_0=1,
    C(z)=(1-z^2)Q(z)=sum c_j z^j.

For p<=q, the first-block critical equations are exactly the linear
coefficient constraints

    c_(q-p+1)=...=c_q=0,
    equivalently b_j-b_(j-2)=0 for q-p+1<=j<=q,                    (4)

where b_(-1)=0. These p equations have independent linear parts in
b_1,...,b_q. If p>q, the same first-block equations generate the unit
ideal. In particular there are no regular fiber points in that range.

**Proof.** For p>0 write the first-block characteristic polynomial as

    A_a(z)=1+z^(2p)+sum_(i=1)^p a_i v_i(z),
    v_i=z^i+z^(2p-i) (i<p),   v_p=z^p.

Its coefficients a_i are quotient coordinates on Sp_(2p). At a regular
point of K they have independent differentials. There are no target logs
in these coordinates, because h is zero on this factor. The derivatives
of the restricted character with respect to a_i are

    F_p=c_q,
    F_(p-j)=c_(q+j)+c_(q-j)=c_(q-j)-c_(q+2-j),  1<=j<p.           (5)

The second equality follows from anti-reciprocity
c_(2q+2-l)=-c_l; coefficients outside the polynomial range are zero.
The middle coefficient c_(q+1) is therefore zero. Starting with F_p and
F_(p-1), successive equations impose c_q,c_(q-1),...,c_(q-p+1)=0.
This is a triangular equivalence of ideals, not only pointwise vanishing.
When p>q this list includes c_0=1. When p<=q, the listed c_j=b_j-b_(j-2)
have independent linear parts, by their triangular dependence on b_j.
For p=0 the constraints are empty and (3) is automatic.

This strengthens the earlier elementary bound p<n for this representation.
It is not a sufficiency theorem for every p<q. Even-rank equality is settled
below; smaller kernels can leave further critical equations.

## 2. All even-rank equality strata

**Theorem 2 (component count and multiplicity).** Let n=2m, m>=1, so
G=Sp_(4m) and rho=V_(omega_(2m)). Let h be any semisimple target with
dim ker(h)=2m. Write

    K=Sp_(2m) x product_j GL_(r_j),   sum r_j=m,
    N=2^m m!/(product_j r_j!),
    d=2m^2+sum_j r_j(r_j-1).                                      (6)

Then theta^(-1)(h) has exactly N irreducible components that meet G_reg.
Each has dimension d and generic scheme multiplicity one. At every
regular point of the fiber its completed local ring is C[[u_1,...,u_d]].
There is no exceptional value of the nonzero target eigenvalues giving a
double component within this target stratum.

These are full-fiber components, but components wholly in the nonregular
locus and embedded structure there are not classified. The statement does
not say the full fiber is reduced everywhere. The known Sp_4 singular
thickening remains compatible with this theorem.

**Proof: fixed external spectrum.** Equation (4), with p=q=m, gives
c_1=...=c_m=0. Reciprocity and c_0=1 force

    Q(z)=Q_0(z)=1+z^2+...+z^(2m).                                 (7)

Its 2m eigenvalues are the (2m+2)-th roots of unity other than +1,-1.
They are distinct, in m reciprocal pairs, with no pair fixed by inversion.
Thus every quotient critical point lies in the open set where these
external eigenvalues can be ordered in etale local coordinates. This
excludes hidden critical points at the external spectral discriminant.

**Proof: exact elimination.** Write the two reciprocal polynomials of
degree 2m in coefficient coordinates a=(a_1,...,a_m), b=(b_1,...,b_m).
The restricted character is bilinear. Its mixed matrix is

    L_ij=Coeff_(z^(2m)) ((1-z^2)v_i(z)v_j(z))= (E E^T)_ij,
    E_ii=1, E_(i,i+2)=-1, all other entries zero.                 (8)

Direct coefficient comparison proves (8) in every rank; det L=1. Let
b^0 be the coefficients of (7). The identity
(1-z^2)Q_0=1-z^(2m+2) gives chi(a,b^0)=1 for every a. Consequently

    chi_a=L(b-b^0),     chi_b=L^T a+d_0,                          (9)

where d_0 is a constant vector. On an ordered eigenvalue chart choose one
eigenvalue beta_i from each reciprocal pair. The beta_i assigned to a
target block with eigenvalue lambda_j all have target log coefficient
2 I_(2m) lambda_j: the defining trace on a reciprocal pair is 2 lambda_j
d log beta_i. Let c be this length-m vector and define

    J_(i,l)=beta_i partial b_l/partial beta_i.

The matrix J is invertible at (7), since the external spectrum is regular
semisimple in Sp_(2m). The remaining critical equations are precisely

    J (L^T a+d_0)=c.                                             (10)

For every assignment of the external eigenvalues, (10) has a unique
solution a. No condition on c beyond (2) has been used. Near this point
the full critical ideal is (b-b^0,a-a^0): first impose b=b^0 with (9),
then solve (10). Thus its local algebra has length one. Equivalently the
perfect mixed pairing eliminates all 2m transverse variables. There is
no residual one-variable equation of the type responsible for double
components in the type A middle-exterior family.

**Proof: count and full components.** To specify the external quotient
point, assign the m reciprocal pairs to the m positive target slots and
choose an orientation in each pair. There are 2^m m! assignments. Within
a GL_(r_j) block, order is immaterial, giving the divisor product r_j!.
The action is free because all source eigenvalues are distinct. Hence N
is exact; all the quotient critical points are accounted for.

Every quotient point has a regular K lift, and the critical equations
force theta(g)=h there. At a fiber point equivariance implies Z_G(g)
is contained in K, so regularity in K and G is equivalent. A K quotient
fiber is irreducible of dimension dim K-2m=d and has a dense regular
orbit. The isolated reduced critical point makes its closure a full-fiber
component, of multiplicity one, by the general critical-potential theorem.
Conversely every regular fiber point lies over one of the points listed.
This completes the count and the local-ring assertion.

This theorem includes arbitrary collisions among the nonzero target
eigenvalues. When all are distinct up to sign, N=2^m m! and d=2m^2;
when all belong to one nonzero pair of multiplicity m, N=2^m and
d=3m^2-m. The previous root-of-unity/unipotent example identifies one
component in this general family; its singular behavior was studied
separately. The present argument counts every component meeting G_reg.

## 3. The odd-rank boundary has one residual variable

The preceding all-rank recurrence also gives a precise next problem.
Let n=2m+1, p=m, q=m+1, m>=1. The largest kernel permitted by (3)
has dimension 2m. Equation (4) now imposes c_2=...=c_(m+1)=0 and leaves
one spectral parameter t=b_1. The external polynomial must be

    Q_t(z)=sum_(j=0)^(m+1) z^(2j) + t sum_(j=0)^m z^(2j+1).       (11)

Consider targets with K=Sp_(2m) x (C*)^(m+1), so their nonzero
eigenvalues lambda_i are distinct up to sign. Every point in a fiber
has external eigenvalues beta_i distinct up to inversion and avoiding
plus or minus one. Indeed, otherwise an external symplectic root
subgroup would centralize g without centralizing h, contradicting
Z_G(g) subset Z_G(h). Thus all fiber points lie in the external regular
spectral chart. Locally choose source eigenvalues beta_i(t), beta_i(t)^(-1)
so that Q_t(z)=product_i (1+beta_i(t)z)(1+beta_i(t)^(-1)z).
They are roots of Q_t(-x), not of Q_t(x); the sign matters in odd rank.

**Proposition 3 (exact odd-rank reduction).** At every regular point of
such a fiber the completed transverse algebra is the one-variable
Jacobian algebra of

    w_h(t)=t-2 I_(2m+1) sum_(i=1)^(m+1) lambda_i log beta_i(t).    (12)

In particular the critical equation is

    1=2 I_(2m+1) sum_i lambda_i beta_i'(t)/beta_i(t).              (13)

The roots beta_i are chosen on an etale local chart; (12) is a formal
potential there, not a globally single-valued function of t. The parameter
t is a character coefficient, not a group element. On this chart simple
zeros of w_h' give multiplicity-one components; an isolated zero of order
mu gives multiplicity mu. An identically zero derivative would give a
nonisolated critical locus and is not assigned a finite multiplicity.
No all-rank root count or exceptional-target classification is proved here.

**Proof.** The p=m first-block equations have rank m in the external
coefficient coordinates, by (5), and the external torus-to-coefficient
map is etale on the specified chart. The general linear-elimination
lemma therefore removes 2m of the total 2m+1 quotient variables. To
find the residual character, set the first-block variables a_i to zero,
so A_base=1+z^(2m). On (11),

    (1-z^2)Q_t=1+t z-t z^(2m+3)-z^(2m+4).

The coefficient of z^(2m+1) after multiplying by A_base is t. The target
logs involve only the external beta_i, with coefficients 2 I lambda_i.
This proves (12), and differentiation proves (13). The component and
multiplicity conclusion follows only when the critical point is isolated.

## 4. What is known, what was deduced, and how it was checked

The primitive character is classical; see
[Etingof, Lecture 31.1](https://ocw.mit.edu/courses/18-755-lie-groups-and-lie-algebras-ii-spring-2024/mit18_755_s24_lec05.pdf).
The trace construction and equivariance are published in
[Kostant--Michor, Section 2](https://arxiv.org/pdf/math/0109066).
Regular quotient theory is from
[Steinberg, Theorem 1.5 and Section 6](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf).
The repository's [general critical-potential proof](GENERAL_CAYLEY_FIBER_CRITICAL.md)
and [perfect symplectic pairing](SYMPLECTIC_FIBER_COMPONENT.md) provide
the local geometry. These are the inputs to the all-target component
count, strengthened kernel bound, and odd-rank reduction written here.
No exhaustive literature search or priority claim is made for these deductions.

[Independent exact checks](symplectic-target-strata-verification.json)
include original primitive contraction-kernel matrices, all Lie-algebra
trace pairings in the tested groups, character comparisons from Weyl
weights, symbolic first-block ideal reductions and signed root-assignment
enumeration. The all-rank assertions rely on the written proofs, not on
extrapolating these checks. No runtime API or PDF is changed.

The next mathematical goal is to classify the critical points of (12)
uniformly, including the exceptional target parameters where its critical
equation acquires a multiple zero. This determines whether higher
multiplicities occur in the odd-rank symplectic boundary family and yields
the missing counterpart of the type A quadratic. A useful next result
would be a root-count or multiplicity bound from (11)-(13), together with
the exact hypotheses needed to include targets with repeated nonzero
eigenvalue pairs. Nonregular fiber structure remains a separate problem.
