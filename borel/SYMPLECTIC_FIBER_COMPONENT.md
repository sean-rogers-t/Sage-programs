# Symplectic collapsed components: generic reducedness and singular thickening

Date: 2026-10-23. Characteristic zero. Written all-rank component proof;
separate finite exact checks. Literature novelty is not established.

## 1. Setup and theorem

For m>=1 let G=Sp_(4m)(C) and rho=V_(omega_(2m)), the irreducible
primitive middle exterior power. Use the simply connected group, its full
weight lattice and the normalized trace Springer map

    I tr(theta(g) X)=tr(rho(g) d rho(X)),
    I=binomial(4m-2,2m-1)-binomial(4m-2,2m-3).

The differential is faithful; the representation has central kernel {+1,-1}.
These are group-map and complex fiber statements, not assertions about
integral cohomology, Borel images, or minimum Cayley degree.

Let B_0 be the regular semisimple element of Sp_(2m) with

    Q_0(z)=det(1+z B_0)=1+z^2+...+z^(2m).

Thus its eigenvalues beta are all (2m+2)th roots of unity except +1,-1.
On a symplectic orthogonal direct sum put

    t=I_(2m) direct-sum B_0,
    H=Z_G(t)=Sp_(2m) x (C*)^m,
    h=theta(t), K=Z_G(h), F=theta^(-1)(h),
    C=K . (t U(H)).

The dot denotes conjugation and F is the scheme-theoretic fiber.
The [realization theorem](SYMPLECTIC_CENTRALIZER_REALIZATION.md) gives
theta(tu)=h for every unipotent u in H.

**Theorem.** C is a closed irreducible component of F_red, and its generic
scheme multiplicity in F is one, for every m>=1. At every point of C whose
first symplectic block is regular unipotent,

    completed local ring of F = C[[s_1,...,s_d]],
    d=dim K-2m.                                           (1)

In particular there is no type A odd-rank double structure at these points.
The theorem allows coincidences among the external eigenvalues of h.
It does not assume an unproved all-rank absence of those coincidences.

More explicitly the first block of h is zero, and its external eigenvalue
on the beta eigenspace is

    lambda(beta)=(beta^2-1)(1-beta)^(2m)/(2I).              (2)

These numbers are nonzero and lambda(beta^(-1))=-lambda(beta).
For each unordered nonzero target pair {lambda,-lambda}, let r_lambda
be the multiplicity of lambda in the external block. Then

    K=Sp_(2m) x product_{pairs} GL_(r_lambda),
    sum_{pairs} r_lambda=m,
    d=2m^2+sum_{pairs} r_lambda(r_lambda-1).                (3)

If the external target eigenvalues are all distinct, K=H, C=t U(H),
and d=2m^2. Exact checks through m=40 find no external collisions; the proof
of (1)-(3) does not rely on extending a finite collision check to all ranks.

## 2. A perfect coefficient pairing

For A,B in Sp_(2m) write their reciprocal characteristic polynomials as

    D_a(z)=det(1+zA)=1+z^(2m)+sum_(j=1)^m a_j v_j(z),
    Q_q(z)=det(1+zB)=1+z^(2m)+sum_(j=1)^m q_j v_j(z),
    v_j=z^j+z^(2m-j) (j<m),   v_m=z^m.

The primitive character restricted to these two blocks is

    chi(A,B)=[z^(2m)] (1-z^2) D_a(z) Q_q(z).              (4)

It is bilinear in a and q. Let L be its mixed coefficient matrix:

    L_ij=[z^(2m)] (1-z^2)v_i(z)v_j(z).

If E is the m-by-m matrix with E_ii=1, E_(i,i+2)=-1 and all other
entries zero, direct coefficient extraction gives

    L=E E^T,             det L=1.                        (5)

For example the diagonal entries of L are 2 except for its last two
entries, which are 1; entries two positions off the diagonal are -1.
For m=1 the sole entry is 1. Formula (5), rather than a computed rank
table, proves invertibility in every rank.

Let q^0 be the coefficients of Q_0. Since

    (1-z^2)Q_0(z)=1-z^(2m+2),

equation (4) is identically 1 for every A when q=q^0. Consequently

    partial chi / partial a = L(q-q^0).                  (6)

In particular theta_A(A,B_0)=0 for every A in Sp_(2m), not only unipotent A.
Conversely, for fixed q=q^0 the variation of the external coefficient
gradient when a changes from a^0 to a is

    (partial chi / partial q)(a,q^0)
      -(partial chi / partial q)(a^0,q^0)=L(a-a^0),       (7)

where a^0_j=binomial(2m,j).

## 3. Target centralizer and the full fiber

Differentiate (4) on the beta, beta^(-1) pair at A=I. The coefficient
identity above cancels the factor 1-z^2 against Q_0. The two logarithmic
derivatives are obtained from

    [z^(2m)] (1+z)^(2m) z beta/(1+z beta)
       =1-(1-beta)^(2m).

Taking their difference and dividing by 2I gives

    ((1-beta^(-1))^(2m)-(1-beta)^(2m))/(2I),

which is (2) because beta^(2m+2)=1. Nonvanishing shows that the first
symplectic block cannot mix with the external block in K. The standard
centralizer of a semisimple symplectic Lie algebra element gives (3).

The polynomial identity [g,theta(g)]=0 implies that all entries of
[g,h]=[g,h-theta(g)] belong to the ideal defining F. Thus F is a closed
subscheme of K, including infinitesimal points. This is why working
inside K here controls the full ambient fiber.

Inside K, C is the product of the Sp_(2m) unipotent variety and the
closed semisimple K_external-conjugacy class of B_0. Both factors are
irreducible. The first has dimension 2m^2 and the second has dimension
sum r_lambda(r_lambda-1). This proves closedness and (3) before asserting
that C is a component.

## 4. Scheme-theoretic elimination at regular unipotent points

Fix g_0 in C with first block A_0 regular unipotent. The invariant
coefficient map A -> (a_1,...,a_m) is smooth at A_0. This is the standard
regular-element criterion for the adjoint quotient: the fundamental
characters of Sp_(2m) are a_j-a_(j-2), a triangular change of generators,
and their differentials are independent at a regular element.

Consequently the gradients of the a_j are independent on an open
neighborhood. By the nondegenerate trace pairing and (6), the first block
of theta(g)=h generates precisely the ideal (q-q^0). This is ideal
equality: choose an invertible minor of the matrix of the m gradients
to recover the coefficients L(q-q^0) from the entries of theta_A.
Conversely every entry is a combination of those coefficients.

Near B_0 the external matrix is regular semisimple. After an etale
coordinate change one can label reciprocal eigenvalue pairs b_j,b_j^(-1).
The map from these m torus coordinates to q is etale, since B_0 has
distinct eigenvalues avoiding +1,-1. Its Jacobian in logarithmic torus
coordinates is therefore invertible. The remaining external coordinates
are the K_external orbit directions.

The equations q=q^0 fix the labeled eigenvalues to beta_j. They leave
only the orbit directions in that block. The target h is scalar on
each GL_(r_lambda) block, so those directions do not occur in the fiber
equations. For these fixed eigenvalues, (7) and the invertible torus
Jacobian show that the external equations generate precisely (a-a^0).
The normalization by I and the factor two in the torus trace pairing are
nonzero constants and do not affect these ideal equalities.

Thus in etale/formal coordinates the entire fiber ideal is

    (q-q^0, a-a^0).                                      (8)

There are 2m transverse coordinates, and all remaining coordinates are
free. This proves (1). Since the regular-unipotent part is dense in the
closed irreducible C, C is a component and has generic multiplicity one.
Neither a radical nor a tangent-space approximation was used in (8).

The structural reason for reducedness is the perfect mixed pairing (5):
the first block fixes all external spectral coefficients and the external
block fixes all first-block invariants. No residual scalar equation
remains, unlike the type A middle-exterior model.

## 5. A singular point already has nonreduced structure

Generic reducedness does not extend to the whole component, even for m=1.
Here G=Sp_4, rho has dimension 5, B_0=diag(i,-i), I=2 and K=Sp_2 x C*.
Write a=tr(A), q=b+b^(-1). The character is chi=1+a q. Write

    A=(a/2)I+X,        X=[[x,y],[z,-x]].

At t=(I,B_0), the first-block fiber equations are qX=0. The external
equation is a(b-b^(-1))=4i. Choosing the formal branch b near i gives

    a=4/sqrt(4-q^2),
    x^2+yz=a^2/4-1=q^2/(4-q^2).

The formal coordinate w=q/sqrt(4-q^2) has nonzero linear coefficient.
The completed full-fiber local ring is therefore exactly

    C[[x,y,z,w]] / (wx, wy, wz, x^2+yz-w^2).              (9)

Its radical is (w,x^2+yz): its reduced germ is the usual nilpotent cone
of sl_2, equivalently the unipotent variety of Sp_2 near I. But w^2!=0
and w^3=0. The nonzero element w^2 is killed by x,y,z,w, so there is an
embedded associated point at the vertex. Away from that vertex the
ring is reduced. This is an explicit obstruction to claiming that the
component carries only the reduced unipotent scheme everywhere.

## 6. Conjecture and a concrete test

**Conjecture (no additional reduced branches along the collapsed component).**
For every m>=1 and every g in C, the closed immersion C -> F induces

    (completed local ring O_(F,g))_red
        = completed local ring O_(C,g).                  (10)

The intended equality is the natural quotient identification, not an
unspecified isomorphism. Thus no other reduced fiber component passes
through C. Infinitesimal structure supported on its singular unipotent
strata is allowed; (9) shows that it must be allowed. This conjecture
does not claim that C is the whole global fiber.

What supports (10): the theorem proves it on the entire regular-unipotent
locus in every rank, and (9) proves it also at the remaining point when
m=1. Higher-rank singular strata have not been settled. Finite generic
rank checks do not provide evidence about those missing strata.

A useful sufficient test is isolation of t in the torus fiber
(theta restricted to T)^(-1)(h). The standard Jordan-part compatibility
of the equivariant map gives theta(g)_s=theta(g_s): the translated
unipotent map in the centralizer is nilpotent and commutes with theta(g_s).
If the local semisimple spectral solution is isolated, nearby reduced
fiber points have semisimple part K-conjugate to t and hence lie in C.
This would prove the conjecture set-theoretically locally; the reduced
formal identification then follows for these complex algebraic local
rings. Scheme thickness would still require separate analysis.

An additional, simpler question is whether (2) has distinct external
values in every rank. A positive answer would simplify (3) to C=t U(H)
and d=2m^2 uniformly. No such assertion is needed for the theorem or (10).

The next substantive step is a local torus-fiber isolation argument,
starting with the m=2 singular strata. Success would remove the regularity
hypothesis in the statement about reduced support, rather than adding
another generic-rank example.

## 7. Sources, independent checks and limits

[Kostant--Michor, Sections 2 and 4](https://arxiv.org/pdf/math/0109066)
give the trace construction and translated centralizer maps.
[Steinberg, Theorem 1.5](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf)
gives independence of fundamental-character differentials exactly on
the regular locus. Both sources were consulted for this note. The
primitive-character identity and regular centralizer descriptions are
classical inputs. The pairing, component deduction, and singular local
calculation are proved above; their novelty has not been established.
Targeted searches did not locate this formulation and are not exhaustive.

The [verifier](tools/verify_symplectic_fiber_component.py) checks (5) through
m=20, target values and collisions in exact cyclotomic fields through
m=40, and full ambient differentials for m=1,...,4. In dimensions 5 and
42 it independently computes the differential from the original trace
construction on actual contraction-kernel representation matrices and
compares all 1396 tangent pairings with polynomial differentiation. It also
checks the entire m=1 centralizer chart by ten symbolic original trace
pairings before checking the radical and nonzero nilpotent socle in (9).
These calculations share Sage and earlier signed-wedge helpers; they
are independent formulas, not an external proof audit or a novelty review.
