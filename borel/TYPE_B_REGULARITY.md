# A uniform regularity obstruction in type B, including spin representations

Subsequent result: the [arbitrary-representation classification](MIXED_REGULARITY.md)
resolves the mixed-central-character question left outside this note.
No representation preserves torus regularity at rank at least three;
Spin5=Sp4 admits torus-positive mixtures, but no nontrivial representation
preserves global regularity at any rank at least two. The present note
retains its narrower central-scalar theorem and proof.

Date: 22 October 2026. Status: written repository deduction, internally
audited and corroborated by exact calculations using two character methods.
The trace construction, highest-weight theory and Spin character lattice
are known inputs. Publication priority and external review are unverified.
All groups and representations here are algebraic over C. This note concerns
regularity, not integral Borel images or minimum Cayley degrees.

## 1. The theorem and its precise scope

Let G=Spin_(2n+1)(C), n>=2, with covering map pi:G -> SO_(2n+1)(C)
and central kernel {1,c}. Let rho be a nontrivial finite-dimensional
representation on which c acts by one scalar s in {1,-1}. In particular,
rho may be any nontrivial irreducible representation, or a direct sum of
irreducibles with the same central parity. The differential is injective,
since the Lie algebra is simple, and its trace form is nondegenerate.
Define the normalized Springer map by

\[
 B_\rho(\theta_\rho(g),H)
 =\operatorname{tr}(\rho(g)d\rho(H)),\qquad
 B_\rho(X,H)=\operatorname{tr}(d\rho(X)d\rho(H)).
\]

**Theorem.** Every such rho has a regular semisimple torus element t
whose image theta_rho(t) is nonregular. One t, independent of rho,
works for every representation satisfying the stated central-scalar
hypothesis. Hence none of these maps preserves global regularity.

**Corollary.** Every nontrivial finite-dimensional representation of
SO_(2n+1)(C) fails regularity preservation, including reducible ones.
Its pullback to Spin has central scalar +1. Every nontrivial irreducible
of Spin fails, including the spin representation and all higher spinorial
highest weights. The theorem leaves reducible Spin representations
mixing the two central characters unclassified.

Rank one is different: Spin_3 is SL_2. Its two-dimensional spin
representation is the defining SL_2 representation and preserves
regularity; all other nontrivial irreducibles fail by the
[type A classification](TYPE_A_REGULARITY.md). The present rank-at-least-two
proof must not be used to rule out that exception.

## 2. Torus coordinates and the Spin lift

In orthonormal coordinates epsilon_i, the roots are

\[
 \{\pm\epsilon_i,\ \pm\epsilon_i\pm\epsilon_j:i<j\}.
\]

Write an SO torus point as diag(z_1,...,z_n,1,z_1^(-1),...,z_n^(-1)).
It is regular semisimple exactly when

\[
 z_i\ne1,\qquad z_i\ne z_j,\quad z_i z_j\ne1\quad(i<j).
\]

In particular z_i=-1 is permitted. Regular here means centralizer
dimension n; it does not require a connected centralizer in SO. This
distinction matters: the Weyl group can fix such a regular SO point.
The lift to the simply connected Spin group will not be fixed.

The Spin character lattice is

\[
 P=\mathbb Z^n\ \cup\ ((1/2,\ldots,1/2)+\mathbb Z^n).
\]

Indeed integrality against the simple coroots epsilon_i-epsilon_(i+1)
and 2epsilon_n means that consecutive coordinates differ by integers
and 2mu_n is integral. The root lattice is Z^n, so its quotient in P
has order two; the integer-coordinate characters are exactly the ones
that descend to SO.

An auxiliary tuple a in (C^*)^n specifies a Spin torus point by
mu(t(a))=prod_i a_i^(2 mu_i); its SO coordinates are z_i=a_i^2.
The kernel of a -> t(a) consists of sign tuples whose product is +1.
Changing just one a-coordinate's sign multiplies t by c. Indeed all
integer characters have ratio +1 and all half-integral characters
have ratio -1. These characters separate torus points, so this verifies
the claimed central translate on the actual Spin torus, rather than
only on its SO image.

Set

\[
 \boxed{\quad a=(i,2,3,\ldots,n),\qquad
 z=(-1,4,9,\ldots,n^2).\quad}
\]

Every root evaluates differently from 1: no coordinate equals 1,
the coordinates are distinct, and no pair has product 1. Thus t(a)
and its SO projection are regular semisimple.

Let w be the Weyl element which inverts the first coordinate and fixes
the others, the reflection in the short root epsilon_1. Since i^(-1)=-i,

\[
 w(t(a))=t(-i,2,3,\ldots,n)=ct(a).
\]

This calculation also resolves the possible sign ambiguity of the lift:
both lifts satisfy the same identity because w fixes c.

## 3. Proof by central character and equivariance

The trace identity and rho(c)=s id give theta_rho(ct)=s theta_rho(t).
Equivariance gives theta_rho(w(t))=w theta_rho(t). Therefore

\[
 w\theta_\rho(t)=s\theta_\rho(t).
\]

On Lie torus coordinates y=(y_1,...,y_n), w sends y to
(-y_1,y_2,...,y_n). If s=+1, this forces y_1=0. The short root
epsilon_1 vanishes. If s=-1, it forces y_2=...=y_n=0. Since n>=2,
the short root epsilon_2 vanishes. In both cases the output is
nonregular, proving the theorem for arbitrary highest weights and
same-parity sums. The proof does not extrapolate finite computations.

In fundamental-weight notation lambda=sum a_j omega_j, all its weights
are integer-coordinate if a_n is even and half-integer-coordinate if
a_n is odd. Hence s=(-1)^(a_n). No assumption on the other coefficients
or on the size of the representation is needed.

## 4. Independent character-factor audit and normalization

Write the weights as mu with multiplicities m_mu. Weyl invariance gives
the trace Gram matrix

\[
 \sum_\mu m_\mu\mu_i\mu_j=S_\rho\delta_{ij},\qquad
 S_\rho=\sum_\mu m_\mu\mu_1^2>0.
\]

Then the actual normalized coordinates, with epsilon_i(H_j)=delta_ij,
are

\[
 y_i=S_\rho^{-1}\sum_\mu m_\mu\mu_i
                 \prod_j a_j^{2\mu_j}.
\]

Relative to the vector trace form one may write S_rho=2I_rho, but
I_rho need not be integral in this normalization: the Spin_5 spin
representation has S=1 and I=1/2. Confusing the two indices changes
the map's scale, although not the zero obstruction.

There is also a direct character proof of the forced zeros, independent
of the central-translate calculation. Put x_i=z_i+z_i^(-1).
For integer weights, invariance under each coordinate inversion implies
chi=F(x_1,...,x_n). Thus

\[
 y_1=S_\rho^{-1}(z_1-z_1^{-1})\partial_{x_1}F=0
 \quad\text{when }z_1=-1.
\]

For half-integral weights, the inversion-invariant Laurent characters
factor as

\[
 \chi(a)=\prod_i(a_i+a_i^{-1})F(x_1,\ldots,x_n).
\]

To see this, group monomials into independent sign orbits. In one
coordinate each odd-power sum a^(2k+1)+a^(-(2k+1)) is
(a+a^(-1)) times a polynomial in a^2+a^(-2), by the elementary
Chebyshev recurrence. The multivariable assertion follows orbit by
orbit, including negative exponents. At a_1=i the first factor is zero.
Every logarithmic derivative in another coordinate retains that factor,
so y_j=0 for j>=2. This proves the same obstruction directly from the
character, including for a same-parity reducible sum.

For the basic spin representation the weights are all (sign_1/2,...,
sign_n/2), each once. Therefore

\[
 S_{\rm spin}=2^{n-2},\qquad
 y_i=2^{1-n}(a_i-a_i^{-1})\prod_{j\ne i}(a_j+a_j^{-1}).
\]

At our point y_1=2^(2-n)i prod_(j=2)^n(j+j^(-1)) is nonzero, while
y_2=...=y_n=0. For Spin_5 this is (5i/2,0), a nonzero singular
output. For Spin_3 the same expression has only one coordinate 2i,
which is regular, illustrating exactly why n>=2 is necessary.

## 5. Differential and harmonic consequences, with limits

The full group differential is singular at t. Choose any root alpha
which vanishes on theta_rho(t). Equivariance differentiated along
conjugation gives

\[
 d\theta_{\rho,t}([E_\alpha,t])
 =[E_\alpha,\theta_\rho(t)]=0.
\]

The input vector is nonzero because alpha(t)!=1. This proves singularity
of the full Cayley differential. It does not assert that the torus
Jacobian is singular, or determine the entire determinant or branch locus.

Kostant--Michor's Theorem 3.2 consequently rules out their particular
unlocalized harmonic factorization with Harm(G)=theta_rho^*Harm(g)
under these hypotheses. It does not rule out a different harmonic
complement, a localized factorization, or a different rational Cayley map.
In particular this is not a bound on the minimum Cayley degree of G.

## 6. Exact corroboration and provenance

[The verifier](tools/verify_type_b_regularity.py) compares Freudenthal
weight multiplicities with an independent type B Weyl alternant and
its exact logarithmic derivatives, both at positive generic points and
at the universal Q(i) witness. It checks the complete trace Gram matrix,
integer/half-integer central parity and a vanishing output root. Exterior
powers of the vector representation and the explicit sign-cube spin
weights supply additional independently enumerated character checks.
Same-parity reducible sums and the rank-one exception are controls.
See the [machine-readable report](type-b-regularity-verification.json)
for the exact finite range. These computations corroborate the proof;
they do not establish the arbitrary-rank claim by exhaustion.
The current run covers 55 irreducibles (34 integer-weight and 21
spinorial), with 110 character/derivative comparisons and 55 universal
witness checks. It adds 13 exterior-monomial character checks, seven
independent spin sign-cube comparisons, eight same-parity reducible
sums and eight closed spin-formula controls, including rank one.

Known inputs include the trace-projection definition and equivariance in
[Kostant--Michor, Sections 2--3](https://arxiv.org/pdf/math/0109066),
and standard type B root and weight lattices. Their Section 7 already
studies the basic spin representation and its relation to the matrix
Cayley transform. The character factor above is consistent with that
published calculation; no novelty is claimed for the basic spin formula.
The uniform obstruction for arbitrary highest weights is presented here
as a deduction. No exhaustive literature search or priority determination
for that classification has been completed. This extends the
[type C central--Weyl argument](TYPE_C_REGULARITY.md) to a different
root system and does not resume the shelved integral Spin-kernel project.
