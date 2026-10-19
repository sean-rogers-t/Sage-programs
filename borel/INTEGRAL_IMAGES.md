# Integral images of fixed Springer/Borel generators

This guide describes `A.integral_degree(d)` for a `GeneratorAlgebra`. It computes the integer span of products of the exact chosen generators, in polynomial degree `d` (cohomological degree `2*d`). The [literature audit](INTEGRAL_LITERATURE.md) separates published results, deductions and questions not settled by the sources reviewed.

## Fix the integral form first

For the specific orthogonal tail sources, the
[closed-form API](ORTHOGONAL_API.md) computes Smith multiplicities and index
exponents without constructing this guide's coefficient matrices. Use this
general matrix API for membership, preimages, or different input generators.

Let `f_j` be the chosen homogeneous Levi-invariant Cartan polynomials, and let `a_j=beta(f_j)`. When every `a_j` has integral Schubert coefficients, the chosen map restricts to

$$
\phi_{\mathbb Z}:\mathbb Z[y_1,\ldots,y_m]\longrightarrow H^{2*}(G/P;\mathbb Z),\qquad y_j\longmapsto a_j.
$$

The grading is `deg(y_j)=deg(f_j)`. The library checks every generator image for integrality, including when the requested degree is zero or has no monomials. Rational coefficients in a Cartan polynomial are allowed if its Borel image is integral. A nonintegral Schubert image raises `ValueError`; denominators are never silently cleared.

This chosen Z-form is not automatically the intersection of a polynomial representation algebra with the ordinary integral representation ring. For the defining/vector Springer map in types B/C, Part II's normalized inputs are `e_i(x^2)`, while the literal tail polynomial-character intersection uses `4^i*e_i(x^2)`. For SO in type D it additionally uses `2^k*product(x)` for the Euler generator. These scale the degree-`d` image by `2^d`; the audit gives the source locators and proof. This concerns the stated tail factors, not arbitrary representations or full Levi lattices.

## Smith factors and what the index means

List all weighted degree-`d` monomials in the formal generators. Let `M_d` be their integer Schubert coefficient matrix, with monomials as columns. The degreewise image is its column lattice `L_d`. A Smith decomposition gives

$$
U M_d V=D,\qquad D=\operatorname{diag}(s_1,\ldots,s_r,0,\ldots),\quad 0<s_1\mid\cdots\mid s_r.
$$

Here `U,V` are unimodular integer changes of basis. Put

$$
\operatorname{sat}(L_d)=(L_d\otimes\mathbb Q)\cap H^{2d}(G/P;\mathbb Z).
$$

Then

$$
\operatorname{sat}(L_d)/L_d\cong\bigoplus_{i=1}^r\mathbb Z/s_i\mathbb Z,\qquad
[\operatorname{sat}(L_d):L_d]=\prod_i s_i.
$$

A saturated image can still be a proper rational subspace of cohomology. If the ambient degree has rank `b_d`, the full cokernel has free rank `b_d-r` and the displayed finite torsion factors. These factors measure a quotient by the selected image; the integral cohomology of the complex flag variety itself is torsion free.

Only nonzero Schubert-support rows are stored. Adding zero rows does not change Smith factors or saturation in that rational span. `include_ambient=True` requests `b_d` by generating the ambient Schubert basis in that degree; it is optional because that enumeration can cost more than computing the selected image.

## Public API

| Method or property | Meaning |
| --- | --- |
| `D=A.integral_degree(d, include_ambient=False)` | Exact integral image report in one polynomial degree |
| `D.matrix`, `D.monomials`, `D.schubert_indices` | Immutable integer image matrix and its labels |
| `D.rank`, `D.kernel_dimension` | Image rank and rank of the degreewise kernel |
| `D.smith_factors` | Positive nonzero invariant factors, including ones |
| `D.torsion_factors` | Only invariant factors greater than one |
| `D.saturation_index`, `D.is_saturated` | Index in the saturation and whether it equals one |
| `D.image_lattice_basis`, `D.saturation_lattice_basis` | Integer basis columns in the stored Schubert-support coordinates |
| `D.ambient_dimension`, `D.free_cokernel_rank` | Ambient rank and free cokernel rank, or `None` when not requested |
| `D.full_ambient` | Whether the rational image has full ambient rank; this does not imply integral surjectivity |
| `D.contains(a)` | Whether a homogeneous Schubert expansion belongs to this integer image |
| `D.preimage(a)` | One integer-coefficient formal polynomial mapping to `a`; raises if none exists |
| `D.as_dict()` | JSON-safe matrices, labels and invariant factors |
| `D.lattice()`, `D.compare(other)` | Embedded lattice and exact comparisons, sums and intersections; see [the comparison guide](INTEGRAL_COMPARISONS.md) |

Membership requires the same flag-variety context and degree; zero is allowed in every degree. A rational class with noninteger Schubert coefficients returns `False`. The preimage is in the existing rational formal polynomial ring, with every coefficient integral. It need not be unique when the map has relations.

## A normalization check on projective space

```python
from sage_borel import BorelMap, GeneratorAlgebra

b = BorelMap("A3", omitted_nodes=1)
h = b.schubert([1])
A = GeneratorAlgebra(b, {"x": 2*b.weights[0]})
x, = A.domain.gens()
D = A.integral_degree(2, include_ambient=True)
assert D.smith_factors == (4,)
assert D.saturation_index == 4 and D.full_ambient
assert not D.contains(h**2)
assert D.contains(4*h**2)
assert D.preimage(4*h**2) == x**2
assert A.evaluate(D.preimage(4*h**2)) == 4*h**2
```

Both `h` and `2*h` generate the same rational algebra. Their integer degree-`d` spans differ by `2^d`. Redundant generators can change the lattice too: the pair `2*h,3*h` generates `h` by an integer Bezout combination, so computing only rational pivot columns would give the wrong lattice.

## Your type D generators, including Euler

```python
from sage_borel import classical_borel

d = classical_borel("D", 4, 2)
A_D = GeneratorAlgebra(d, {"q": d.tail_polynomial(2, 1),
                         "u": d.euler_polynomial(2)})
D4 = A_D.integral_degree(4, include_ambient=True)
assert D4.smith_factors == (1, 2, 2)
assert D4.saturation_index == 4
assert D4.rank == 3 and D4.ambient_dimension == 4
assert D4.free_cokernel_rank == 1
assert not D4.is_saturated
```

Thus the image in cohomological degree eight has rank three, its saturation quotient is `(Z/2)^2`, and the full cokernel also has one free direction. Here the omitted second squared elementary generator is `u^2` literally, so adding it does not change the selected Z-algebra. This example is independently checked against Cartesian divided differences and determinantal divisors.

## Normalized and literal type C forms

```python
c = classical_borel("C", 3, 1)
q = c.tail_polynomial(1, 1)
normalized = GeneratorAlgebra(c, {"q": q})
literal = GeneratorAlgebra(c, {"Q": 4*q})
assert normalized.integral_degree(2).smith_factors == (1,)
assert normalized.integral_degree(4).smith_factors == (1,)
assert literal.integral_degree(2).smith_factors == (4,)
assert literal.integral_degree(4).smith_factors == (16,)
```

For the natural normalized type C tail, all nonzero Smith factors are one at finite and stable rank. This follows from the published quaternionic-Grassmannian fibration and integral cohomology, as explained with a primitivity proof in the [literature audit](INTEGRAL_LITERATURE.md). It is a classical geometric consequence, not a theorem newly established by finite computation. The literal representation-ring normalization changes the answer as shown above.

## Exceptional inputs and nonintegral normalizations

```python
from sage_borel import exceptional_borel

f = exceptional_borel("F4", 4)
A_F = GeneratorAlgebra(f, {info.name: f.generator_polynomial(info.name)
                           for info in f.generators()})
assert A_F.integral_degree(4).smith_factors == (1, 12)

g = exceptional_borel("G2", 1)
A_G = GeneratorAlgebra(g, {info.name: g.generator_polynomial(info.name)
                           for info in g.generators()})
try:
    A_G.integral_degree(0)
except ValueError as error:
    assert "nonintegral Schubert coefficient" in str(error)
else:
    raise AssertionError("G2/P1's original quadratic normalization is nonintegral")
```

The G2/P1 catalogue quadratic maps to `-3/4*sigma_(2,1)`. Multiplying it by four defines an explicitly different input form. The experiment records that alternative separately and does not identify it with the literal integral representation-ring intersection. No automatic classification of exceptional integral polynomial characters is implemented.

## Why the lattice bases and preimages are correct

Write `E` for the first `r` columns of `U^-1`. Since `V` is unimodular, multiplying `M` by `V` does not change its column lattice. Consequently `E*diag(s_i)` is a basis of the image, and `E` is a basis of its saturation. All these matrices are exact and immutable.

For a target coordinate vector `a`, membership is equivalent to: its coordinates outside stored support vanish; the last `b-r` entries of `U*a` vanish; and each of its first `r` entries is divisible by the corresponding `s_i`. Dividing those entries gives a vector `y` with `D*y=U*a`. Then `V*y` is an integer vector of formal-monomial coefficients, since `M*(V*y)=a`. This is the returned preimage. The implementation checks the exact Smith transformation and agreement with the rational rank.

## Scope, cost and reproducibility

For uniform mathematical results beyond the bounded calculations, see
[Orthogonal Springer images](ORTHOGONAL_SMITH.md). The accompanying proofs
give powers-of-two restrictions, stable saturation, complete quadric cases,
the type B weighted-minor formula and a type D integral gluing rule.
The [general type B proof](TYPE_B_GENERAL_SMITH.md) now establishes the
closed [all-rank formula](TYPE_B_CLOSED_CONJECTURE.md), including the
formerly missing integral comparison. The
[general D square theorem](TYPE_D_SQUARE_THEOREM.md) now identifies
every normalized square-sector factor with B, including all factors
for smaller block three. The
[general even-tail gluing theorem](TYPE_D_EVEN_GLUING_THEOREM.md)
now completes the full classification: the only gluing correction
replaces the known number of unit exponents by exponent one.
Its integral relative Euler argument and proof dependencies are
recorded in the [new audit](TYPE_D_EVEN_GLUING_AUDIT.md).
These are deductions in the notes, not an attribution of the exact Smith
formulas to the cited papers or a claim of originality. The
[claim-by-claim proof/source audit](PROOF_STATUS.md) records the distinctions.

The method computes one exact degree at a time. A list of bounded reports is not a proof of a uniform formula or a complete integral ring presentation. `A.presentation()` continues to compute the complete kernel over QQ; its Groebner basis must not be relabeled an integral kernel certificate. Integral image calculations currently require Schubert-coordinate matrices rather than the rational coinvariant backend, whose coordinates need not form an integral Schubert lattice.

`A.clear_cache()` clears integral reports and generator integrality checks alongside the existing caches. Smith transformations and integer preimages can grow large even when the invariant factors are short. The default avoids ambient-basis generation; it does not promise fast high-rank dense expansions.

The [independent tests](tests/test_integral.py) compare Smith factors with gcds of minors and compare small coefficient matrices with Cartesian BGG oracles. The [verification record](integral-verification.json) contains exact counts. The [bounded experiment](experiments/integral_springer_images.py) computes normalized and literal classical forms, selected exceptional forms, and explicit preimage/obstruction witnesses; its [recorded results](experiments/integral-springer-results.json) cover the stated finite degrees only.

Run inside SageMath from `borel/`:

```bash
PYTHONPATH=src sage -python examples/integral_images.py
PYTHONPATH=src sage -python -m unittest tests.test_integral -v
PYTHONPATH=src sage -python experiments/integral_springer_images.py
```
