# Complete kernels of chosen Borel generators

`GeneratorAlgebra.presentation()` computes the entire rational kernel ideal of
the chosen generator map, a monomial basis of its image, and its finite Hilbert
series. Its completeness certificate uses exact Schubert images and a finite
monomial border. No user-supplied degree cutoff is needed.

All degrees below are **polynomial degrees**, half the cohomological degrees.
The formal variables have the positive degrees of their input Cartan
polynomials. The computation applies to any `GeneratorAlgebra`, including
classical and exceptional types and arbitrary parabolics. Its cost depends on
the selected image, which can still be large.

## A first example

```python
from sage_borel import BorelMap, GeneratorAlgebra

b = BorelMap("A2", omitted_nodes=1)  # projective plane
A = GeneratorAlgebra(b, {"h": b.weights[0]})
h, = A.domain.gens()
P = A.presentation()

assert P.domain is A.domain
assert P.generator_degrees == (1,)
assert P.relations == (h**3,)
assert P.standard_monomials == (1, h, h**2)
assert P.dimension == 3
assert P.hilbert_coefficients == (1, 1, 1)
print(P.hilbert_series)              # t^2 + t + 1
assert P.reduce(h**4 + 2*h**2) == 2*h**2
assert P.kernel_ideal() == A.domain.ideal(h**3)
assert P.certificate.certified
assert P.verify().certified
```

The quotient `A.domain / P.kernel_ideal()` is the **selected image algebra**.
It need not be all of `H*(G/P; QQ)`. The basis consists of formal monomials whose
Schubert images are independent; the API does not require the complete ambient
Schubert basis to be enumerated.

## API

`A.presentation()` returns a frozen `KernelPresentation`. Repeated calls use the
cached presentation; `A.clear_cache()` clears it together with the algebra's
other computed data.

| API | Result |
| --- | --- |
| `P.domain` | The existing formal polynomial ring `A.domain` |
| `P.generator_degrees` | Positive polynomial degrees of the formal variables |
| `P.relations` | Tuple containing a reduced monic Groebner basis of the complete kernel |
| `P.standard_monomials` | Tuple of monomials giving a rational basis of the quotient and image |
| `P.dimension` | Total rational dimension of the selected image |
| `P.hilbert_coefficients` | Image dimensions in polynomial degrees zero through the last nonzero degree |
| `P.hilbert_series` | The finite Hilbert polynomial in `QQ[t]` |
| `P.certificate` | The exact completeness checks recorded for this presentation |
| `P.kernel_ideal()` | The defining ideal in `P.domain` |
| `P.reduce(p)` | The unique normal form of a formal polynomial in the standard monomial basis |
| `P.verify()` | Recheck the presentation and return a `KernelCertificate` |
| `P.verify().certified` | Whether the completeness and output-consistency checks pass |
| `P.as_dict()` | JSON-safe exact source, relation, basis and certificate data |

The Groebner basis uses the domain's existing ordinary degree reverse
lexicographic term order. The mathematical grading is supplied separately by
`generator_degrees`: ordinary total degree and the grading degree need not
agree. Inputs can have algebraic dependencies or zero cohomology image; neither
prevents the construction.

This complements the bounded methods in [IMAGE_RELATIONS.md](IMAGE_RELATIONS.md).
`A.degree(d).relations` is a vector-space basis of one homogeneous kernel
component. `P.relations` instead generates the **entire kernel ideal**. The
number of Groebner basis elements need not be the minimum number of ideal
generators.

## Two-step symplectic example

```python
from sage_borel import two_step_symplectic

F = two_step_symplectic(n=3, b=1, k=1)
A = F.joint_algebra()
c1, e1 = A.domain.gens()
P = A.presentation()

assert P.generator_degrees == (1, 2)
assert set(P.relations) == {c1**4 + c1**2*e1 + e1**2, e1**3}
assert P.dimension == 12
assert P.hilbert_coefficients == (1, 1, 2, 2, 2, 2, 1, 1)
assert P.reduce(c1**4) == -c1**2*e1 - e1**2
assert P.reduce(e1**3) == 0
assert P.verify().certified
```

Here the selected image has presentation

$$
\mathbb Q[c_1,e_1]/(c_1^4+c_1^2e_1+e_1^2,\ e_1^3),
\qquad \deg(c_1)=1,\quad \deg(e_1)=2.
$$

Its standard monomials are `c1**i * e1**j` for `0 <= i < 4` and
`0 <= j < 3`, arranged in discovery order. Thus its Hilbert series is

$$
(1+t+t^2+t^3)(1+t^2+t^4)
=1+t+2t^2+2t^3+2t^4+2t^5+t^6+t^7.
$$

The ambient two-step flag has cohomology dimension 24. This selected
middle-plus-tail image has dimension 12. The independent
[two-step derivation](experiments/TWO_STEP_DERIVATION.md) supplies a separate
mathematical presentation for the family, rather than being an assumption of
the generic algorithm.

## The exact staircase construction

Write

$$
R=\mathbb Q[y_1,\ldots,y_r],\qquad \deg(y_i)=d_i>0,
\qquad \phi:R\longrightarrow H^*(G/P;\mathbb Q).
$$

The map substitutes the chosen Cartan polynomials and applies the Borel map.
Its kernel is an ideal homogeneous for the positive weights `d_i`.

Start with the candidate monomial `1`. Process candidates in increasing
weighted degree, breaking ties by the domain's ordinary degree reverse
lexicographic order. Retain a monomial when its exact Schubert image is
linearly independent of the retained images in the same weighted degree.
Only retained monomials enqueue their products by the formal variables.

If a candidate image is dependent, exact rational linear algebra expresses it
in the earlier retained images. Record the monic relation

$$
m-\sum_{b\in B}q_b b,
\qquad \phi(m)=\sum_{b\in B}q_b\phi(b),
$$

where the sum uses retained monomials of the same weighted degree. Its leading
monomial is `m`. Candidates divisible by an already discovered leading
monomial are discarded. Continue until no candidate remains. This constructs
the standard monomial staircase and the relations at its boundary.

The construction is FGLM-style: it discovers an ideal basis from exact linear
dependencies of monomial images in a finite-dimensional algebra. It is not
literally the original change-of-order FGLM algorithm, whose input includes an
existing Groebner basis.

### Why the existing term order is preserved

The discovery order first compares positive weighted degree and then compares
ordinary degree reverse lexicographic order. It is a global monomial order:
it is well ordered and multiplication preserves comparisons. Inside one
weighted homogeneous component, it agrees with the domain's term order.

Every recorded relation is weighted homogeneous, so it has the same leading
monomial in both orders. More generally, let `f` belong to the kernel and take
the weighted homogeneous component containing its ordinary-degrevlex leading
monomial. That component also belongs to the kernel, and its leading monomial
in the discovery order is the ordinary-degrevlex leading monomial of `f`.
Consequently, a Groebner basis for the discovery order is also a Groebner basis
for the existing domain order. No new formal polynomial ring or conversion of
the public variables is required.

The tails of the recorded relations use only standard monomials. Their
leading monomials do not divide one another, and every leading coefficient is
one. The resulting basis is therefore reduced and monic.

## Why the certificate proves completeness

Let `B` be the reported standard monomials, `G` the reported relations, and
`L` the monomial ideal generated by the leading monomials of `G`. The
`KernelCertificate` records five conditions:

| Condition | Meaning |
| --- | --- |
| `relation_images_zero` | Every relation in `G` has zero exact Schubert image |
| `basis_images_independent` | The images of `B` are linearly independent in every weighted degree |
| `standard_monomials_reduced` | No monomial in `B` is divisible by a leading monomial in `G` |
| `border_closed` | For every `b` in `B` and variable `y_i`, the monomial `b*y_i` lies in `B` or in `L` |
| `contains_unit` | The monomial `1` belongs to `B` |

The implementation additionally checks `groebner_basis_reduced` (monic
homogeneous relations, minimal leading terms and tails in the standard basis)
and `summary_data_correct` (the formal domain, weights and Hilbert data match
the reported basis). `P.verify()` recomputes the geometric top degree rather
than trusting the recorded value. `certified` is true when all these checks
hold. The five mathematical conditions above imply
completeness, rather than merely checking that some proposed relations vanish.

First, `B` is precisely the set of monomials outside `L`. To see this, suppose
there were a monomial outside both `B` and `L`, and choose one of smallest
ordinary total degree. It is not the unit. Remove one variable factor. The
parent is outside `L`, because monomial ideals are closed under multiplication;
by minimality the parent belongs to `B`. Border closure now contradicts the
choice of the original monomial. The reverse inclusion follows from
`standard_monomials_reduced`.

Division by the monic relations terminates in a linear combination of `B`.
Let `J` be the ideal generated by `G`. The zero-image condition gives
`J <= ker(phi)`. For any `f` in `ker(phi)`, its remainder is a linear
combination of `B` with zero image. Independence forces that remainder to be
zero, so `f` belongs to `J`. Therefore

$$
J=\ker(\phi),\qquad R/J\cong\operatorname{im}(\phi).
$$

The same argument shows that `G` is a Groebner basis: a nonzero kernel
polynomial with a standard leading monomial would retain that term after
division and could not have zero remainder. This proof also gives the unique
normal form returned by `P.reduce(p)`.

The certificate uses the library's exact Borel evaluator. Independent
coinvariant-ring calculations in the tests check the implementation against a
different computational method; the certificate itself does not claim to be
a second implementation of the Borel map.

## Termination and cost

Let `D` be the dimension of the selected image. Each retained monomial adds an
independent image, so there are exactly `D` retained monomials. Only these
monomials enqueue variable multiples. At most `1 + r*D` distinct candidates
are therefore created, including the initial unit; pruning can reduce the
number of actual image evaluations.

There is also a geometric degree bound. Put

$$
N=\dim_{\mathbb C}(G/P)=\ell(w_0w_{0,P}),
$$

the highest nonzero polynomial degree of the cohomology ring. A retained
monomial has nonzero image and hence weighted degree at most `N`. Every
candidate is the unit or a variable multiple of a retained monomial, so no
candidate has degree greater than `N + max(generator_degrees)`. Images above
`N` are zero by the grading. These bounds prove finite termination; completion
is detected by the exhausted border, not by assuming a guessed degree is
sufficient.

This is output-sensitive in the number of selected image basis monomials,
not a promise of fast calculations in every rank or exceptional type.
Schubert supports, rational coefficients and input polynomials can grow, and
the selected image itself can approach the full cohomology ring. In
particular, supporting type E8 as an input does not make every E8 presentation
practical. The degreewise API remains useful when only a bounded portion is
needed.

An exploratory `C9/P8` calculation with the single tail generator `x9**2`
exceeded a 45-second process limit while other validation was running. No
completed certificate or timing is claimed for that case. Even a small image
can have expensive Schubert supports. Computing classical images directly in
a coinvariant quotient is a possible future optimization.

## Verification records

The [new independent suite](presentation-verification.json) checks 13 complete
elimination ideals in types A/C/D/G2, three exceptional ample-divisor kernels,
65 normal forms and 32 vanishing relations. It also rejects forged degree
bounds, incomplete bases and inconsistent Hilbert data. The complete 39-test
regression suite passed under SageMath 9.5; package and historical-file checks
are recorded in [presentation-validation.json](presentation-validation.json).

The [runnable example](examples/kernel_presentations.py) writes
[exact sample presentations](presentation-examples.json) for projective space,
Gr(2,4), the two-step C3 image, a D4 Euler example, the full G2 flag ring and an
F4 divisor image. Run it with `--output PATH` to reproduce the JSON data.

## Mathematical context

The rational Borel presentation and finite-dimensional flag cohomology are
classical. This feature computes a certified presentation of a user-selected
subalgebra of that known cohomology ring. It makes no novelty claim for the
underlying Borel theorem or for the staircase linear-algebra method.

* A. Borel, *Sur la cohomologie des espaces fibres principaux et des espaces
  homogenes de groupes de Lie compacts*, Annals of Mathematics 57 (1953),
  115-207, especially Theorem 26.1. [Original paper](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Borel-Sur.pdf).
* J.-C. Faugere, P. Gianni, D. Lazard and T. Mora, *Efficient Computation of
  Zero-dimensional Groebner Bases by Change of Ordering*, Journal of Symbolic
  Computation 16 (1993), 329-344. [Publisher and DOI](https://www.sciencedirect.com/science/article/pii/S0747717183710515).
