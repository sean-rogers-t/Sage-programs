# Integral images for G2 and F4

This is the first complete exceptional-group experiment: every cohomological
degree for both G2 maximal parabolics and all four F4 maximal parabolics.
The objects are explicitly chosen generator Z-algebras. They are not asserted
to be the intersection of the polynomial representation algebra with an
integral representation ring. Numbering is Sage/Bourbaki throughout.

The subsequent [G2 character-form proof](G2_CHARACTER_FORM.md) now identifies
the literal seven-dimensional character source and both complete image
rings. Those images are properly contained in both forms studied here in
every positive nonzero degree. The subsequent
[F4 character-form proof](F4_CHARACTER_FORM.md) identifies its literal
source as Sym_Z(6 Q_long)^(W_L) and supplies complete comparisons for all
four maximal parabolics in [new tables](F4_CHARACTER_TABLES.md).

## Relation to the Springer papers

[Kumar--Xie, Part III](https://kumar.math.unc.edu/papers/kumar86.pdf),
Theorems 5.2 and 6.2, supplies complex-algebra generators and their Schubert
images. Our [audited catalogue](EXCEPTIONAL.md) implements their input
polynomials, retaining the normalization and correction metadata.
This experiment studies integer combinations of all products of those
generators. Knowing the individual generator images over C does not by itself
give the integral image of the algebra they generate.

The ambient integral rings are classical objects. For example,
[Duan--Zhao](https://arxiv.org/pdf/math/0511332) gives presentations for the
two end-node F4 Grassmannians. These tables are computations for the chosen
Springer inputs inside such rings, not new presentations of the ambient
cohomology and not an originality claim. We have not identified the exact
selected lattices with a published integral representation-ring construction.

The [complete readable tables](EXCEPTIONAL_INTEGRAL_TABLES.md) give every
Smith factor in all 172 degree/form rows. All nonunit factors in these
G2 cases are powers of 2; all primes occurring in these F4 cases are 2 or 3.
This is an exhaustive statement about these six spaces and two specified
forms, not a general assertion for exceptional groups or other normalizations.

## Exactly which integral forms?

Write f_i for a catalogue polynomial, beta for the Borel characteristic map,
and sigma_w for an integral Schubert class. For nonzero beta(f_i), write

\[
\beta(f_i)=\sum_w a_{iw}\sigma_w,\qquad
D_i=\operatorname{lcm}_w\operatorname{denom}(a_{iw}),\qquad
g_i=\gcd_w(D_i a_{iw}).
\]

We compare two explicit choices:

1. **Image denominators cleared:** generators D_i f_i. The multiplier D_i
   is the least positive integer making this individual image integral.
2. **Individually primitive images:** generators (D_i/g_i) f_i. Each
   nonzero image is now a primitive vector in its Schubert lattice.

Zero images, if present, are left unchanged. Primitivity of the individual
vectors does not imply that their products generate a saturated subalgebra.
Neither construction canonically determines a character lattice.

For either chosen list F_i, the degree-d lattice is

\[
L_d=\sum_{\sum_i m_i\deg F_i=d}\mathbb Z\,
\beta(F_1^{m_1}\cdots F_r^{m_r})
\subset H^{2d}(G/P;\mathbb Z).
\]

All tables have full rational rank. Their nonzero Smith factors s_j therefore
describe the complete finite quotient H^{2d}(G/P;Z)/L_d; factors 1 are retained
in the tables. Degree d means cohomological degree 2d.

For G2/P1 only, D_quadratic=4, while its primitive-image multiplier is 4/3.
All F4 catalogue images are already integral. The primitive form divides
F4/P1 e2 by 2 and e3 by 4, and divides F4/P2 e3 by 2. It makes no changes
on F4/P3 or F4/P4.

## API

```python
from sage.all import QQ
from sage_borel import exceptional_algebra, exceptional_borel

# Exact catalogue inputs are the default; no silent denominator clearing.
g2 = exceptional_algebra("G2", 1, scales={"quadratic": 4})
assert g2.integral_degree(3).smith_factors == (2,)

b = exceptional_borel("F4", 1)
original = b.generator_algebra()
primitive = b.generator_algebra(scales={"e2": QQ(1)/2, "e3": QQ(1)/4})
assert original.integral_degree(4).smith_factors == (2, 48)
assert primitive.integral_degree(4).smith_factors == (1, 48)
assert original.integral_degree(4).compare(primitive.integral_degree(4)).left_in_right
```

Scales use canonical catalogue names and positive exact rational values;
omitted entries have scale 1. The returned object is a `GeneratorAlgebra`,
so the existing degree, kernel, integral-image, membership, preimage and
lattice-comparison methods apply. Calling `integral_degree` on the default
G2/P1 algebra raises an error because its quadratic image is fractional.
The same helper accepts the catalogued E6/E7/E8 cases; the full tables here
only concern G2/F4, and the E8 catalogue remains partial.

## A complete elementary G2 calculation

Both spaces have one Schubert class sigma_d in each degree d=0,...,5.
Let h=beta(omega_i) on G2/P_i, and let q denote the catalogue quadratic.
The elementary Chevalley products give the following precise descriptions.

On P1, beta(q)=-3h^2/4, and

\[
(1,h,h^2,h^3,h^4,h^5)
=(\sigma_0,\sigma_1,\sigma_2,2\sigma_3,2\sigma_4,2\sigma_5).
\]

Thus both cleared and primitive-image forms generate exactly Z[h], whose
integral presentation is Z[h]/(h^6).

On P2 put y=beta(q)=sigma_2. Here

\[
h^2=3y,\quad hy=2\sigma_3,\quad y^2=2\sigma_4,
\quad hy^2=2\sigma_5.
\]

Consequently the image is Z[h,y]/(h^2-3y,y^3). Indeed the quotient has
the free additive basis 1,h,y,hy,y^2,hy^2, and the displayed nonzero images
lie in distinct degrees, proving that there are no additional relations.
This is an elementary deduction from the Chevalley products; the independent
calculation checks every product column in both cases.

For either parabolic and either selected form, the Smith factors in degrees
0,...,5 are therefore **1,1,1,2,2,2**. This is a finite calculation and proof,
not a new uniform family theorem.

## An explicit F4 obstruction that survives primitive rescaling

On F4/P1 in polynomial degree 4, use the Schubert basis with words
(2,3,2,1) and (4,3,2,1). The source monomials e2 and h^4 have matrix

\[
\begin{pmatrix}-10&2\\28&4\end{pmatrix}.
\]

Its entries have gcd 2 and its determinant has absolute value 96, so its
Smith factors are (2,48). After replacing e2 by e2/2, the matrix becomes

\[
\begin{pmatrix}-5&2\\14&4\end{pmatrix},
\]

with factors (1,48). Thus the quotient still contains a Z/48 summand,
including 3-torsion, even though all individual generator images in the
second form are primitive. This torsion is in the quotient by the selected
image, not in the cohomology group itself. The determinant computation is
also an elementary independent check of this row's Smith arithmetic.

## Why the F4 calculation is smaller than the input list suggests

The experiment first removes an image generator only when it obtains and
verifies an integer polynomial expression in retained generators.
For each middle F4 parabolic, write h for the common linear generator,
A for the first quadratic generator image and B for the second. Put Q=A+B.
In cohomology the exact relations are

\[
\begin{array}{c|cc}
 & A & B\\\hline
P_2 & -11h^2+9Q & 11h^2-8Q\\
P_3 & 11h^2-8Q & -11h^2+9Q.
\end{array}
\]

Since Q=A+B in the original algebra, the substitutions are integral in both
directions. The computation can use h,Q,e3 instead of two redundant linear
generators, two quadratics and a cubic. No rational change of lattice is
hidden in this optimization. The JSON records the forward and reverse
expressions used for each case.

## F4 duality explains the matching odd-primary quotients

There is a structural deduction for the original F4 catalogue forms. Write
J_i for their image algebra in H_i=H*(F4/P_i;Z). Then J_4 and J_1 are
isomorphic as graded integral algebras, as are J_3 and J_2. These isomorphisms
do not identify the ambient integral Schubert lattices. For each odd prime p,
the p-primary parts of (H_i/J_i) in each degree agree with those for 5-i.
The same odd-primary statement holds for the primitive-image forms.
The quotients in this assertion are degreewise additive groups.

Here is the argument, including the finite identities on which it depends.
Define a rational polynomial substitution

\[
T(\omega_1)=2\omega_4,\quad T(\omega_2)=2\omega_3,
\quad T(\omega_3)=\omega_2,\quad T(\omega_4)=\omega_1.
\]

Put c_1=c_2=2 and c_3=c_4=1. The F4 Cartan matrix gives
T s_i=s_{5-i}T and T(alpha_i)=c_i alpha_{5-i}. Thus

\[
\partial_{5-i}T=c_iT\partial_i.
\]

It follows by iterating divided differences that T sends sigma_w to
2^{nu(w)} sigma_{rev(w)}, where rev replaces every node i by 5-i and
nu(w) counts letters 1 or 2 in a reduced word. This count is independent
of the chosen word: odd braid relations join only nodes of the same
root length. The substitution preserves the rational global invariant
ideal and exchanges the two relevant Levi-invariant subrings. It therefore
induces a graded ring map H_i -> H_{5-i}; the displayed Schubert formula
proves that it is integral. Moreover T^2 is multiplication by 2^d in
polynomial degree d, so this ambient map is an isomorphism after inverting 2.

At the end nodes, put h=beta(y4), u=beta(e2). The exact relations
beta(e1)=-h^2 and beta(e3)=-h^2 u show that J_i=Z[h,u]. At the middle
nodes the preceding integral changes of generators show J_i=Z[h,Q,r],
where r=beta(e3). Direct substitution and the Borel relations give

| Direction | Images of the displayed generators in the target |
| --- | --- |
| P4 -> P1 | (h,u) -> (h,u) |
| P1 -> P4 | (h,u) -> (2h,16u) |
| P3 -> P2 | (h,Q,r) -> (h,Q,r) |
| P2 -> P3 | (h,Q,r) -> (2h,4Q,8r) |

These are identities of cohomology images; they are not asserted as literal
polynomial equalities before taking the global invariant quotient. They,
the generator reductions, and the root-intertwining identities are checked
in the [duality verifier](experiments/exceptional_f4_duality.py), including
independent Chevalley checks of the image identities.

The first and third maps are injective because their ambient maps are
injective, and surjective onto the selected image because they map the
displayed generating sets onto each other. This proves the integral algebra
isomorphisms. After inverting 2 the ambient maps are also isomorphisms, hence
so are the quotient maps. This proves the odd-primary assertion. Finally,
the F4 primitive-image forms only divide generators by powers of 2, so they
agree with the original forms after inverting 2 as well.

More explicitly, matching monomials in the end-node generators or in the
middle-node generators gives the original image matrix for P1 or P2 from
that for P4 or P3 by reversing Schubert words and multiplying each source
row by 2^{nu(w)}. The verifier checks this identity for every stored degree.
This explains the paired odd-primary factors without asserting that row
scaling multiplies the individual Smith factors in the same way. It is a
deduction from root duality and the checked generator identities, with no
claim that this formulation is new.

## Reproduce and audit

```bash
PYTHONPATH=src sage -python experiments/exceptional_integral_images.py
PYTHONPATH=src sage -python experiments/exceptional_chevalley_audit.py
PYTHONPATH=src sage -python experiments/exceptional_f4_duality.py
PYTHONPATH=src sage -python -m unittest tests.test_exceptional_integral -v
```

The [complete dataset](experiments/exceptional-integral-results.json) records
all degree matrices, source monomials, Schubert words, scaling choices,
integral generator reductions and Smith lists. The first script uses the
library's BGG divided differences and exact integer matrix operations.
The second reconstructs every saved column using positive-root/coroot
Chevalley multiplication, including top degree, without calling library
expansion or multiplication. Its [audit report](experiments/exceptional-chevalley-audit.json)
binds the checks to the dataset and source hashes.

The independent algorithms share the catalogue polynomial definitions and
Sage root-system data. The existing catalogue tests separately reconstruct
the inputs and compare the journal audit fixtures. Agreement is computational
evidence with these explicit shared inputs, not external proof certification.

The completed run checked 172 degree matrices and independently reconstructed
all 1,650 source columns. The duality audit additionally checked 64 root
intertwining identities, 20 Chevalley image identities and 1,556 row-scaling
entries across all 37 paired degree matrices. The library regression suite
passed all 98 tests, including four tests for the new exceptional algebra API.

## Next questions

The natural continuation is to explain the observed integral indices through
known exceptional cohomology presentations, then extend the same explicit
normalization audit to E6/P1 and E7/P7. A separate, deeper question is whether
one can specify an intrinsic integral form on the representation side and
relate its image to these coordinate-defined forms. Changing the chosen
generators merely to remove their common divisors does not settle that question.
