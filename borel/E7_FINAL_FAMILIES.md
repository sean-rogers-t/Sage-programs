# Completing the remaining E7 adjoint Levi families

26 October 2026. Characteristic zero. These are case proofs for the actual
adjoint Springer--Cayley map, not a uniform theorem for arbitrary characters.
The user has explicitly resumed case-by-case completion. Published geometric
inputs, written deductions and exact computational checks are distinguished.

## 1. The A2+A2+A1 family

**Theorem.** Let G be the complex simply connected group E7 and let K be a
Levi subgroup of derived type A2+A2+A1 in the Weyl class of retained nodes
12356 (Bourbaki numbering). Normalize theta by

    B_Ad(theta(g),X)=tr(Ad(g)ad(X)).

There is a dense open subset U of z(k), consisting of h with Z_G(h)=K,
such that the actual scheme theta^(-1)(h) intersect G_reg is nonempty,
smooth and pure of dimension 14 for every h in U. This covers the three
standard embeddings in this Weyl class.

### Exact reduction

Use the finite central cover and fundamental characters (a,x,b,y,z),
where (a,b) and (y,z) belong to the two A2 factors and x to A1. With the
central cocharacters recorded in the root certificate, the central torus
coordinates are p,q. Eliminate the first A2 endpoint pair by the global
unit-paired reduction in [the quartic note](EXCEPTIONAL_QUARTIC_REDUCTION.md).
Put r=p^3 and t=q^2. This change is a finite etale isogeny of tori, and
the residual potential depends only on r,t. It is

    f=-2+(1/2)v^t H v+ell^t v,       v=(y,z)^t,
    a0=x+r, b0=x+r^(-1), P=a0*b0,
    H=[[-2*a0/(rt), -P],[-P, -2*rt*b0]],
    ell_1=-t*(x^2+rx+r^2-1)-(rx+1)/(r^2 t^2),
    ell_2=-rt^2*(x+r)-(x^2+x/r+r^(-2)-1)/t.             (1)

This is the full character after elimination, including its constant and
linear terms. The original target equations are retained: in these
coordinates they are mu=(r*f_r,t*f_t). In the p,q basis they are
(3*mu_r,2*mu_t). The trace pairing gives an invertible linear conversion
to the coordinates of h. No derivative of an eliminated variable is
discarded: its contribution is a multiple of a vanished critical equation.

We study X=V(f_x,f_y,f_z) in A^3 times (C*)^2, with its actual scheme
structure. Set

    D=det H=P*(4-P),
    k=r*x^2+(r^2+1)*x-3*r,
    D=-(x+r)*(rx+1)*k/r^2.                              (2)

### The invertible-block chart

On D!=0 solve v=-H^(-1)ell, giving the exact residual potential

    R=-2-ell^t adj(H)ell/(2D).

The remaining equation is R_x=0. Its polynomial numerator can be defined
without a long coefficient table by

    J=r^6*t^3*D^2*R_x in Q[x,r,t].                       (3)

Its x-degree is nine and its leading coefficient is 2*r^6*t^3, a Laurent
unit. At (r,t)=(2,1),

    J=2*(4x^2+10x+13)
         *(16x^7+120x^6+184x^5-160x^4-308x^3+140x^2+261x-10).

This degree-nine polynomial is coprime to its derivative over Q. Therefore
the discriminant of J in x is nonzero; J is squarefree over C(r,t).
The unit leading coefficient excludes any nonunit Laurent content, so J
is squarefree in C[r^+/-1,t^+/-1,x]. Consequently X on D!=0 is a reduced
hypersurface of pure dimension two. Its nonsmooth locus has dimension at
most one and hence cannot dominate the two-dimensional central target.
This argument proves reducedness over C, not only over a number field.

### The entire singular-block boundary

We now retain all of D=0; it cannot simply be inverted globally.

On x=-r with r^2!=1, H has rank one and consistency H v+ell=0 requires

    r^2*t^3=1.

The same is true on x=-r^(-1) with r^2!=1. Each consistent base has
dimension one. The determinant has a simple x-zero here. After eliminating
the nonzero quadratic direction, the remaining x-equation is a quadratic
in the kernel coordinate with nonzero leading coefficient: that coefficient
is a unit multiple of D_x. Thus these critical boundary loci have dimension
at most one, including multiple roots of the resulting quadratic.

On k=0 the matrix has rank one everywhere, since P=4. Modulo k, the two
components of adj(H)ell are

    2*(r^2+1)*(rx-1)*(r^2*t^3+1)/(r^2*t),
    2*(x-r)*(r^2+1)*(r^2*t^3+1)/(r^3*t^2).             (4)

The consistent base is contained in the union of

    k=0, r^2+1=0;
    k=0, r^2*t^3+1=0;
    r^2=1, x=r.

Each has dimension at most one. Where D_x!=0 the same kernel-quadratic
argument bounds the critical locus by one. Where D_x=0, k_x=0 and

    x=-(r+r^(-1))/2,       r^4+14r^2+1=0.

These r-values satisfy r^2!=1 and r^2!=-1; (4) then forces r^2*t^3=-1.
The consistent base is finite. Even if the last equation vanishes
identically in the kernel coordinate, the critical locus has dimension
at most one. This includes every collision of determinant factors.

Finally H has rank zero exactly when r=epsilon and x=-epsilon,
epsilon=+1 or -1. Here ell=0 for every t, and the full potential restricted
to this locus is the constant -2. Directly, t*f_t=0 there for all y,z.
Thus the entire rank-zero critical locus maps into the target line mu_t=0,
regardless of its dimension or scheme structure. No assumption of a
regular sequence or reduced boundary was used.

### Smooth generic targets and nonemptiness

Exclude the closures of the boundary images above and of the nonsmooth
locus on D!=0. Each is a proper closed subset of the target. Over the
complement X is smooth of dimension two. Characteristic-zero generic
smoothness, applied to the finitely many components and removing images
of nondominant ones, now gives smooth zero-dimensional generic fibers.
Embedded or nilpotent structure on the excluded boundary cannot survive
this target restriction; on D!=0 reducedness was proved for the actual
hypersurface, not just its support.

For nonemptiness, at

    (x,y,z,r,t)=(-4,-3,-3,1,1)

the relative gradient vanishes, det H=-45, and the affine Hessian has
determinant -810. The Schur complement giving the central logarithmic
Hessian along the local critical branch is

    [[336/5,36],[36,54]],       determinant=11664/5.

Thus the target map is dominant on this smooth branch. The witness itself
need not be centralizer-exact: dominance supplies a dense target open
which can be intersected with that condition. The
[regular-quotient theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md) transfers the
actual quotient-fiber assertion to regular elements of G. Its smooth orbit
factor has dimension 6+6+2=14. Finite central descent and the previously
certified [Weyl conjugacies](E7_MULTIPLIER_BOUNDARY.md) complete the theorem.

## 2. The A2+A1+A1+A1 family and all-Levi E7 completion

**Theorem.** For the E7 Levi with retained nodes 12357, of derived type
A2+A1+A1+A1, generic centralizer-exact central targets have nonempty,
smooth actual regular-source fibers of pure dimension 12. Together with
section 1 and the previous proofs, this establishes generic smoothness
(allowing emptiness) for EVERY proper Levi of E7.

The exact paired-A2 elimination and root reconstruction are in
[the symmetric-pencil note](E7_SYMMETRIC_PENCIL.md#3-status-of-the-companion-a2a13-refinement).
Its residual potential, including the central-only terms, is

    psi=x^2+y^2+w^2-x^2*y^2-x^2*w^2-y^2*w^2
        -xyw*(Ax+By+Cw)-ABxy-ACxw-BCyw-ABC-3,
    h=A^2+B^2+C^2-ABC-4=0,                              (7)

where A=U+U^(-1), B=V+V^(-1), C=UV+(UV)^(-1), with U=q,V=p^3*q.
Let I in Q[x,y,w,A,B,C] be generated by psi_x,psi_y,psi_w,h and ALL
four-by-four minors of the Jacobian of these four equations in all six
variables. The earlier exploratory calculation selected six polynomials
with leading monomials

    x^5, y^5, w^2*B^4, w^3*C^3, A^6, B^5*C.            (8)

Their memberships in the ORIGINAL ideal I are now certified, as follows.

### A homogeneous Hilbert-function membership certificate

Here is the exact algebraic verification principle, including the direction
of inclusion that was formerly missing. Homogenize the original nineteen
generators with a new variable H, obtaining the homogeneous ideal I^h
in Q[x,y,w,A,B,C,H]. This is the ideal of the homogenized generators;
we do NOT assume that it is saturated or is the projective closure.

The [certificate archive](e7-homogeneous-membership-certificate.json.gz)
contains 2,347 monic homogeneous polynomials G. The exact checks establish:

1. G is a Groebner basis over Q, by rational Buchberger verification.
2. Every original homogenized generator reduces to zero by G. Hence
   I^h is contained in J=(G) over Q.
3. At a recorded prime p dividing none of the relevant denominators,
   the leading monomials of G survive, and the ideal of the ORIGINAL
   homogenized generators over F_p equals (G mod p). This reverse
   inclusion is checked by direct finite-field standard-basis arithmetic.

For each degree d, finite coefficient-matrix ranks give

    dim_Q(Q[variables]/I^h)_d
       <= dim_Fp(F_p[variables]/I^h_p)_d.

The inclusion in item 2 gives the opposite bound by the Hilbert function
of J. Monicity, the rational Groebner identities and denominator control
make G mod p a Groebner basis with the same leading monomials. Therefore
J and J_p have identical Hilbert functions. Item 3 sandwiches the Hilbert
function of I^h between two equal functions. Thus I^h=J over Q in every
degree. This is a proof of reverse inclusion over Q, not a probabilistic
inference from finite-field dimensions. Homogeneity is essential.

Dehomogenizing H=1 now gives the ORIGINAL affine ideal I. All six selected
polynomials reduce exactly to zero by the dehomogenized G. Consequently
they belong to I. The generator uses prime 32003 and slimgb; the
[separate verifier](tools/verify_e7_companion_membership.py) uses a different
admissible prime and the std algorithm, and independently reconstructs the
original Jacobian-minor ideal. Both perform exact rational Buchberger
verification. No long rational multiplier expression is assumed or needed.

The initial ideal of I contains (8). Their radical monomial support is

    (x,y,A,wB,wC,BC),

whose dimension is one: at most one of w,B,C can remain free. Hence
dim V(I)<=1, also after extending scalars to C. This supplies precisely
the dimension bound that the earlier note left conditional.

### Geometric finish, including the actual scheme

Away from V(I), the four equations psi_x,psi_y,psi_w,h have Jacobian
rank four in A^6. Their actual incidence is therefore smooth of dimension
two. The torus-to-Cayley-surface map is finite and etale away from its four
fixed points U,V in {+1,-1}. On this etale open the same statement holds
for the full torus critical scheme. The inverse image of V(I) has dimension
at most one, so its target image cannot dominate the two-dimensional target.
At the four fixed points all three trace differentials dA,dB,dC vanish;
the full central target is zero there and is excluded by centralizer-
exactness. These exclusions retain all nilpotents on every surviving
chart: no radical replacement of the fiber equations has been made.

Generic smoothness on the surviving smooth two-dimensional incidence gives
smooth finite generic quotient fibers. For nonemptiness, the actual torus
branch x=y=w=0 has relative Hessian determinant -216 at U=V=1, and the
central logarithmic Hessian of -ABC-3 is

    [[-16,-8],[-8,-16]],        determinant=192.

This proves dominance before removing the proper bad images. The witness
is allowed to have target zero; it proves dominance, not centralizer-
exactness at that specific point. The regular quotient and finite central
descent then give the group-fiber assertion, with orbit dimension
6+2+2+2=12. The existing 123 cases, section 1's three cases, and this one
case exhaust the 127 proper standard subsets of E7.

## 3. Status and verification

The trace-gradient construction, regular quotient geometry, characteristic-
zero generic smoothness and the reduced-hypersurface criterion are standard
inputs. The reduction and boundary proof above are repository deductions,
with no assertion of literature priority. The
[Sage verifier](tools/verify_e7_final_families.py) reconstructs the full
adjoint character from root data; the
[independent audit](tools/audit_e7_final_families.py) uses native Cartan-root
reconstruction and separate SymPy identities. Neither checks only a
specialized fiber: the boundary identities hold over the full Laurent ring.

Section 1 adds three standard subsets to the previous 123, and section 2
adds the last one: E7 coverage is now 127/127. The earlier exploration
archive deliberately retains its historical false promotion flag; the new
homogeneous membership certificate and verification supersede that pending
status. All-proper-Levi E7 generic smoothness is now proved by cases.
The all-simple-group conjecture is NOT complete: E8 still has unresolved
coverage. No literature-wide novelty assertion is made.
