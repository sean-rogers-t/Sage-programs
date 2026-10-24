# Generic degree of the individual half-spin Springer maps of Spin_10

Updated 2026-10-22. This is a characteristic-zero written deduction with
independent exact checks. It is not an externally reviewed proof or a
claim of literature priority. It follows the
[all-rank zero-fiber classification](HALF_SPIN_BOUNDARY.md).

## The statement and its purpose

Let G=Spin_10(C), let S_+ and S_- be its two irreducible half-spin
representations, and let theta_sigma:G -> so_10 be the Springer map
defined by orthogonal trace projection for S_sigma. Thus

    B_sigma(theta_sigma(g),Y) = tr(S_sigma(g) dS_sigma(Y)),
    B_sigma(X,Y) = tr(dS_sigma(X) dS_sigma(Y)).

**Theorem.** Each of theta_+ and theta_- has generic degree **148**.
Their restrictions T -> t and the induced maps T/W -> t/W have the
same degree. Equivalently, over a general complex target, the full
Springer map has 148 distinct inverse images.

The proof is the exact intersection calculation

    148 = 208 - 10*6.

Here 208 is the normalized volume of the half-spin weight polytope in
the actual character lattice of Spin_10. Each of ten coordinate facets
contributes six reduced boundary intersections. Every other boundary
stratum contributes zero.

This answers a different question from the preceding zero-fiber work:
the zero fiber has dimension 33, whereas a general fiber has 148 points.
It measures the number of generic inverse branches. It does not assert
that the map is finite everywhere, classify its branch locus, or give
the minimum degree among all equivariant rational maps.

## 1. Normalization and the actual Spin lattice

Use orthonormal weight coordinates epsilon_1,...,epsilon_5. The weights
of S_+ are mu=(e_1,...,e_5)/2 with e_i in {1,-1} and product(e_i)=1.
All have multiplicity one. Their second moment matrix is 4I, so on T

    theta_+(t) = (1/4) sum_mu mu e^mu(t).

Consequently the logarithmic gradient of chi=sum_mu e^mu differs from
theta_+ by an invertible constant transformation of its target. Such a
transformation does not change the degree.

The character lattice is

    M = P(D5) = Z^5 + Z delta,    delta=(1/2,...,1/2).

A basis is epsilon_1,...,epsilon_4,delta. In this basis the exponent of
the weight e/2 is

    ((e_1-e_5)/2,...,(e_4-e_5)/2,e_5).

These sixteen exponents generate M. In particular, we are not counting
on an auxiliary cover of the Spin torus or replacing M by the root
lattice. The covolume of M in the epsilon coordinates is 1/2.

The other half-spin is the dual representation. The trace definition
gives theta_-(g)=-theta_+(g^(-1)); hence the two degrees agree. We prove
the result for S_+ and check both chiralities computationally.

## 2. The compactified intersection number is 208

Let P be the convex hull of the sixteen weights. It is the
five-dimensional demihypercube. Its facets are:

* Ten coordinate facets mu_i=+1/2 or -1/2, each a four-dimensional
  crosspolytope with eight vertices.
* Sixteen simplex facets, each with five vertices, cutting off the
  odd-sign corners of the cube [-1/2,1/2]^5.

The corner simplices have disjoint interiors and Euclidean volume 1/5!.
Thus the normalized lattice volume is

    Vol_M(P) = 5! * (1-16/5!) / (1/2) = 208.

Every proper face of dimension at most three is a simplex: it is a
proper face of one of the crosspolytope facets or of a simplex facet.
The origin lies in the interior of P.

Let X_P be the normal projective toric variety of the normal fan, using
the lattice M, and let L_P be the ample Cartier line bundle of P. The
five logarithmic-gradient equations, with a constant target subtracted,
are sections of L_P. The standard toric intersection theorem gives

    integral_X_P c_1(L_P)^5 = Vol_M(P) = 208.

We use the actual common zero scheme of these five sections. Their
special coefficients have boundary zeros, so a generic-coefficient
Bernstein count cannot simply be substituted for the torus degree.

## 3. Only coordinate facet orbits can contain boundary zeros

On the torus orbit of a proper face F, the initial gradient equations
are, up to a common invertible character and a constant change of rows,

    sum_(mu in F) mu e^mu = 0.

The constant target monomial disappears because 0 is interior to P.
If F is a simplex, its vertices are affinely independent and lie in a
supporting hyperplane not containing 0. They are therefore linearly
independent. The displayed vector equation would force every monomial
coefficient to vanish, which is impossible on a torus orbit.

It follows that no simplex facet orbit and no orbit of codimension at
least two contains a common zero. This argument rules out zeros even
on singular strata of X_P; we will only need local intersection
calculations at its smooth codimension-one orbits.

## 4. Six critical points on each coordinate facet

The Weyl group is transitive on the ten coordinate facets. Consider
mu_5=1/2. Write z_i=e^(epsilon_i), i=1,...,4, and b=e^delta. Put

    D = z_1 z_2 z_3 z_4,
    E = 1 + sum_(i<j) z_i z_j + D,
    O = sum_i z_i + sum_(i<j<k) z_i z_j z_k.

The full character, including the opposite facet, is exactly

    chi = b E/D + b^(-1) O.

The tangent character lattice of this facet is Z^4, with basis
epsilon_1,...,epsilon_4. A transverse parameter is t=b^(-1).
The primitive normal is -2 epsilon_5^*: it pairs t with 1. The
constant target monomial therefore vanishes to order one, and the
opposite-facet terms to order two. This is where using the actual
Spin character lattice is essential.

At t=0 the five equations reduce to

    E=0,    z_i partial_i E=0  (i=1,...,4).

They have exactly the six points with z_i=+1 or -1 and two of each
sign. Here is a direct exhaustion, over C. Write s=sum_i z_i. From
the derivative equations,

    z_i(s-z_i)+D=0.

Every z_i is a root of x^2-sx-D. If all four coordinates equal a,
the derivative equation gives a^2=-3, but then E=-8. If there are two
distinct roots a,b, occurring p and 4-p times, Vieta's formula gives

    (p-1)a+(3-p)b=0,    ab=-D.

The cases p=1 or 3 force a root to be zero, impossible on the torus.
For p=2 we get b=-a, and ab=-D gives a^2=1. These six sign points
indeed satisfy all the equations.

At each point, the logarithmic Hessian of E has diagonal zero and
off-diagonal entries 1+z_i z_j. Grouping the signs gives two blocks
[[0,2],[2,0]], so its determinant is 16. The six critical points are
nondegenerate. The independent saturated critical algebra over Q has
length six as well.

## 5. Each point has compactified intersection multiplicity one

Let y_1,...,y_4,y_b be the target coordinates for logarithmic
derivatives with respect to z_1,...,z_4,b. Divide the sections by b
and clear the unit D. On the open neighborhood with z_i invertible,
the equations are

    f_i = z_i partial_i E - E + t^2 D z_i partial_i O - y_i tD,
    f_b = E - t^2 D O - y_b tD.

These are local equations on X_P, not equations on a ramified cover:
the z_i and t form the indicated lattice coordinates, with t allowed
to vanish. At any of the six points with t=0, the Jacobian with
respect to z_1,...,z_4,t has determinant

    -16 y_b.

One can see this from the invertible four-by-four Hessian block and
the final row: its four z-derivatives vanish and its t-derivative is
-y_b. Thus the point is a simple common zero whenever y_b is nonzero.
The same argument on the other nine facets imposes only finitely many
nonzero linear conditions on the target. A general target meets all
of them. Total boundary length is therefore 10*6=60.

## 6. Completing the degree argument

The normalized Springer map has derivative equal to the identity at
the identity. Its torus restriction is dominant and generically
finite. In characteristic zero its general torus fiber is reduced.

For a general target the common zero set on X_P consists of this
finite torus fiber and the sixty boundary points just classified.
There are no other boundary zeros. All common zeros lie in the smooth
locus, and the five Cartier sections meet properly. Their intersection
number is the sum of these local lengths. Consequently

    deg(theta_+ restricted to T) + 60 = 208.

For a regular semisimple Y in t, any g with theta_+(g)=Y commutes
with Y: equivariance and conjugation by g give Ad_g(Y)=Y. Thus
g belongs to Z_G(Y)=T. Scheme-theoretically the transverse root
directions are nonsingular for such Y, or equivalently one may use
the regular semisimple parametrizations G x^(N_G(T)) T_reg and
G x^(N_G(T)) t_reg. Hence the torus and group generic degrees agree.
The degree on T/W also agrees, by the degree tower in the commutative
diagram with the two Weyl quotients, both of degree |W|. This is
consistent with Kostant--Michor's Corollary 3.3.

## 7. Independent exact checks and their limits

Run from borel with Sage:

```bash
PYTHONPATH=src sage -python tools/verify_half_spin_degree.py
```

The [verifier](tools/verify_half_spin_degree.py) writes the
[verification report](half-spin-degree-verification.json). It checks:

* Exact polytopes in the P(D_r) basis for r=3,4,5 and both
  chiralities, with volumes 4,32,208. It checks facet types and linear
  independence of every simplex boundary face at r=5, primitive
  normals, and the separate cube-minus-corners volume calculation.
* The saturated Q-algebra E=z_i partial_i E=0 has length six.
  All six points and all six local Jacobian determinants -16y_b
  are checked exactly with symbolic target variables.
* Twelve exact specialized fibers over two finite fields, using both
  the five independent P(D5) torus coordinates and the original
  Spin model b^2=z_1...z_5. Every fiber has length 148.
  A squarefree multiplication characteristic polynomial checks that
  each tested fiber is geometrically reduced.

The original Spin-coordinate equations use the normalized formula

    theta_i = (1/8) sum_(product e=sigma) e_i b prod_j z_j^((e_j-1)/2),
    b^2 = prod_j z_j.

Both models are saturated to keep every torus coordinate nonzero.
Leading monomials and the multiplication characteristic polynomials
are saved in the report. The models give separate equation-level
checks; they share Sage/Singular arithmetic, so they are not separate
computer algebra implementations. The staircase length count is
performed by the repository's native monomial algorithm.

A useful negative control is also retained. Over F_32009, the original
normalized target (1,2,3,4,5) gives length 148 but only 147 geometric
points, for each chirality. A multiplication operator has minimal
polynomial of degree 148 and a characteristic polynomial with a single
repeated factor of total gcd degree one. Thus this fiber is a cyclic
nonreduced algebra, with one double point. This is an exact certificate,
not an inference from failure to find a separating linear form. The
twelve reduced-fiber checks use different targets, including
(2,3,5,7,11) in the original normalized coordinates. Specialization and
genericity are deliberately kept separate.

The low-rank polytope numbers agree with the previously justified
Spin_6 and Spin_8 degrees. In rank four the weights span a sublattice
of index two; retaining the actual Spin lattice preserves the correct
degree 32 rather than inadvertently counting on its image quotient.

Finite-field specializations do not prove the generic degree over C.
The proof is Sections 1--6. A rational specialized Groebner experiment
was stopped without a certificate; no such certificate is claimed.
The proof has been internally audited, not externally reviewed or
formally verified. No runtime API changed.

## 8. Published inputs, deduction, and next goal

The trace-projection construction and the equality of group and
adjoint-quotient degrees are published in
[Kostant--Michor, Sections 2--3](https://arxiv.org/pdf/math/0109066).
Their Section 7 full-spin representation has degree 10 on Spin_10;
that representation is S_+ direct sum S_-, rather than either
individual sixteen-dimensional summand studied here.

The toric top-intersection formula is standard: see Cox--Little--Schenck,
*Toric Varieties*, Theorem 13.4.3
([book](https://doi.org/10.1090/gsm/124)). Half-spin weights, duality and
the D5 character lattice are standard representation theory. The
specific boundary exhaustion, local multiplicities and resulting
degree 148 are the written calculation here. Targeted searches did
not locate this value explicitly, but were not a comprehensive
priority audit; we make no novelty claim.

**Next goal: derive a boundary correction for Spin_12, organized toward
a rank-dependent half-spin degree formula.** The reason is structural:
coordinate facets of a higher-rank half-spin polytope are lower-rank
half-spin polytopes, so the preceding critical-locus classification
can control which solutions escape to the toric boundary. For Spin_12
the facet character is of D5 type. Its zero-character critical locus
has positive-dimensional pieces (already visible by fixing two torus
coordinates at each of +1 and -1 and varying the fifth). Therefore
the subtraction of sixty isolated points does not extend unchanged.
An excess-intersection or resolved-base-locus calculation is needed,
including compatibility where boundary strata meet.

Progress means an exact boundary contribution, justified in the
correct lattice and matched to an independent general-fiber count.
An all-rank law, maximal finite target loci, and boundary scheme
structures remain open here. A larger table alone would not resolve
those questions.
