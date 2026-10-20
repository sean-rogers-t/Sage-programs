# Full type-D images with a two-dimensional tail

**Status: internally audited repository deduction, 19 October 2026.**
This argument determines the full normalized image for every head size
when the specified tail has dimension two. It uses the integral
two-head Chern calculation and the previously proved unit formula.
It is unpublished; novelty and external expert certification are not claimed.

## 1. Statement and an explicit partition rule

Let \(A\geq1\), \(n=A+2\), and
\[
X=\operatorname{OG}(A,2n),\qquad
q_1=x_{A+1}^2+x_{A+2}^2,\qquad p=x_{A+1}x_{A+2}.
\]
The normalized source algebra is \(\mathbb Z[q_1,p]\), since
\(q_2=p^2\). Both displayed generators have polynomial degree two.
As throughout these notes, polynomial degree is half the cohomological
degree, and Smith factors refer to saturation in the selected rational
image, with units retained.

**Theorem.** In each even polynomial degree \(d\), list the partitions
\(\lambda=(r,s)\), allowing \(s=0\), satisfying
\[
r+s=d,\qquad 0\leq s\leq r\leq2A+1,
\qquad r=s>A\ \text{is forbidden}.
\tag{1}
\]
The full normalized Smith list consists of one factor
\[
                         2^{\#\{\lambda_i>A\}}
\tag{2}
\]
for each such partition. Include the empty partition in degree zero.
Odd-degree lists are empty. In particular every factor divides four.

The full formula in [TYPE_D_CANDIDATE.md](TYPE_D_CANDIDATE.md) therefore
holds for \(k=2\) and every \(A\). The entire saturated gluing is
unit gluing. Its dimension is one in each degree \(d=2m\) with
\[
                   \lfloor A/2\rfloor+1\leq m\leq A,
\tag{3}
\]
and zero in every other degree.

## 2. Moving the two specified roots to the head

Consider also
\[
                         X'=\operatorname{OG}(2,2n).
\]
Its complementary tail has dimension \(A\). Pull back both varieties
to a complete type-D flag manifold. These pullbacks are primitive:
global integral Schubert classes give a Leray--Hirsch basis over each
partial flag manifold, including the class one, so each pullback has an
integral additive left inverse. Composing a lattice map with a primitive
injection preserves its Smith factors in its rational-image saturation.

A permutation of coordinate roots is a type-D Weyl group element.
The right Weyl action on the compact complete flag manifold induces an
integral cohomology automorphism and acts by that permutation on the
characteristic classes. Move the two roots in the original specified
tail to the head of \(X'\), and call them \(x,y\). Then the generators
become
\[
                       q_1=x^2+y^2,\qquad p=xy.
\]
Put \(h=x+y\) and \(e=xy\). Since \(q_1=h^2-2e\), their integral
image is exactly
\[
                         \mathbb Z[h^2,e].
\tag{4}
\]
This is the entire even-polynomial-degree part of the head Chern
subring \(\mathbb Z[h,e]\) of \(X'\). Every monomial of even
polynomial degree has an even power of \(h\), because \(\deg h=1\)
and \(\deg e=2\).

The argument compares the original tail algebra to a different algebra
on \(X'\), namely its even head Chern algebra. It does not identify the
two full tail algebras under block exchange. The selected Euler
generator has become the head top Chern class.

When \(A=1\), the complementary \(SO(2)\) parabolic of \(X'\)
omits both fork nodes. Its ordinary isotropic Grassmannian and complete
flag pullback still have the same integral properties. Thus this
boundary does not require a single-node parabolic interpretation.

## 3. The integral Chern lattice gives the factors directly

Apply the two-root Chern model with tail parameter \(A\), as proved in
[TYPE_D_TWO_ROW_SMITH.md](TYPE_D_TWO_ROW_SMITH.md) and detailed in
[TYPE_D_TWO_HEAD_GLUING.md](TYPE_D_TWO_HEAD_GLUING.md#2-the-two-root-chern-model-and-the-quotient-by-square-classes).
In even polynomial degree, its retained unnormalized Giambelli classes
\(\Gamma_\lambda\) form an integral basis of (4). The retained
partitions are precisely (1): they are \(A\)-strict, have at most two
parts, and have first part at most \(2A+1\).

This Chern-model argument applies to every \(A\geq1\). The C and D
two-row raising thresholds differ only at total size \(2A+1\), which
is odd. The missing even symplectic classes are
\(\Gamma_{(2A+2,2j)}=e\,s_{(A,j)}\), \(0\leq j\leq A\).
Thus this use of the Chern model does not require an even new tail.

For clarity, the integral basis assertion follows from an integral
symplectic Schubert basis before quotienting. The missing classes have
first part \(2A+2\) and are exactly the coordinate basis of the even
Chern kernel. Thus quotienting removes those basis vectors and leaves
the retained classes as an integral basis. This is stronger than a
rational basis assertion.

The even orthogonal Giambelli formula gives
\[
                 \Gamma_\lambda=2^{\nu(\lambda)}Y_\lambda,
\qquad \nu(\lambda)=\#\{\lambda_i>A\}.
\tag{5}
\]
Here \(Y_\lambda\) is an individual Schubert class for a type-zero
partition, or the sum of the two typed Schubert classes otherwise.
The \(Y_\lambda\) form the integral plus-eigenlattice basis. This
eigenlattice is primitive in the ambient integral cohomology: it is
the kernel of the integral map \(\iota-1\).

The rational span of (4) is the entire even-degree plus eigenspace.
Consequently (5) is already a diagonal matrix from an integral source
basis to a primitive target basis. Sorting its diagonal entries gives
exactly (2). All comparisons in the preceding section preserve these
factors, proving the partition rule for the original full tail image.

The published input to (5) is the typed Schubert/Giambelli formula of
[Buch--Kresch--Tamvakis, Theorems 1 and 4](https://math.umd.edu/~harryt/papers/oggiam.pdf).
The selected-image Smith conclusion is the deduction made here.

## 4. Identifying the gluing with the proposed formula

Let \(q\) record \(m=d/2\). The factors four in (2) correspond to
distinct \(r,s\in\{A+1,\ldots,2A+1\}\) of the same parity.
For \(N\geq0\), the elementary Gaussian identity
\[
       \sum_{0\leq j<i\leq N-1}q^{i+j}
                         =q{N\brack2}_q
\tag{6}
\]
follows by writing the pair as the partition \((i-1,j)\) in the
\(2\)-by-\((N-2)\) rectangle. Set the right side to zero for
\(N<2\).

If \(A=2v\), the odd and even pairs give respectively
\[
                  q^{2v+2}{v+1\brack2}_q,
\qquad           q^{2v+3}{v\brack2}_q.
\tag{7}
\]
If \(A=2v+1\), their sum is
\[
             (q^{2v+3}+q^{2v+4}){v+1\brack2}_q.
\tag{8}
\]
These are exactly
\[
                    \mathcal B_{A,2,2}(q)
                         +q\mathcal B_{A-1,2,2}(q),
\tag{9}
\]
by the proved general B formula. For \(A=1,2\), terms outside the
allowed B exponent range vanish; (7)--(9) still agree.

The two sector ranks add to the full rational-image rank. Their lists
are known in every degree from the square and Euler comparisons.
The [general unit theorem](TYPE_D_FROBENIUS.md) gives the full unit
count as the sum of the two sector unit counts minus
\[
\mathcal I_{A,2}(q)
  =q^{\lfloor A/2\rfloor+1}
       {\lceil A/2\rceil\brack1}_q.
\tag{10}
\]
There are no exponents above two, by (2), and (9) shows that the full
number of factors four is the sum of the sector numbers. The rank and
unit count therefore force the number of factors two to be the sum
of the sector numbers plus (10). This proves every multiplicity in
the proposed full formula.

Finally, the index identity from the
[integral gluing lemma](ORTHOGONAL_SMITH.md#lattice-lemma) gives
\[
\dim_{\mathbb F_2}W
   =v_2([S:L])-v_2([S_+:L_+])-v_2([S_-:L_-])
   =[q^m]\mathcal I_{A,2}(q).
\]
The unit theorem already supplies that many independent unit-unit
gluing classes. They exhaust \(W\), proving (3) and the assertion
about the gluing flags.

## 5. Scope

Together with the two-head gluing theorem, this settles the full
even-tail formula whenever \(\min(a,k)\leq2\). Odd tails were
already settled by the separate sector comparisons and their opposite
degree parities. The subsequent
[general gluing theorem](TYPE_D_EVEN_GLUING_THEOREM.md) now settles
every even tail. The present argument remains an independent proof
for the two-tail family.

For the literal defining/vector polynomial-character form, multiply
every degree-\(d\) factor by \(2^d\). This note uses only the fixed
normalized generators above.

The independent ordinary-Python
[partition verifier](experiments/verify_type_d_two_tail_gluing.py) and
[report](experiments/type_d_two_tail_gluing_verification.json) check
1,680 degree components through head size 40 against the separate B
two-row rule and compare 33 recorded exact full Smith matrices. No
Sage, Borel engine, or other verifier is imported by that script.
These finite comparisons corroborate the proof, rather than establish
its all-rank claim.

From the repository root:

```sh
python borel/experiments/verify_type_d_two_tail_gluing.py
```
