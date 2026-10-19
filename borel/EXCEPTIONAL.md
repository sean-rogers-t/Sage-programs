# Exceptional generator images

Version 0.2 adds named Levi-invariant generators to the Borel engine for G2,
F4, E6, E7 and the displayed E8 dissertation cases. Inputs are polynomials in
fundamental weights. Their images are computed by the library, rather than read
from a stored answer table.

## Three-parameter calculation

```python
from sage_borel import exceptional_image

answer = exceptional_image("F4", parabolic=1, generator="e2")
print(answer.to_terms())
```

The three arguments are the Cartan type, omitted Dynkin node, and generator
name. They play the role of your classical `(n,k,h)` arguments, while allowing
the different invariant generators of each exceptional Levi subgroup.

For F4/P1, `e2` is the elementary polynomial of index 2 in the squares of

\[
y_1=-\omega_1+2\omega_4,\quad
y_2=-\omega_1+2\omega_3-2\omega_4,\quad
y_3=-\omega_1+2\omega_2-2\omega_3.
\]

Its image is

\[
-10\sigma_{s_2s_3s_2s_1}+28\sigma_{s_4s_3s_2s_1}.
\]

This is polynomial degree 4 and cohomological degree 8.

## Discover and inspect inputs

```python
from sage_borel import exceptional_borel, exceptional_generators

print([g.name for g in exceptional_generators("F4", 1)])

b = exceptional_borel("F4", 1)
info = b.generator_info("e2")
f = b.generator_polynomial("e2")
answer = b.generator_image("e2")

assert dict(answer) == dict(b.expand(f))
print(f)
print(info)
```

The returned map also provides the ordinary `BorelMap` operations, including
`expand`, `coefficient`, `basis`, `poincare_dual` and `last_stats`. Reuse it when
working with several generators on one flag variety. `exceptional_polynomial`
is the corresponding standalone convenience function for inspecting an input.

Names distinguish coordinate alphabets when a Levi has several factors.
Unambiguous elementary generators use `e1`, `e2`, and so on; ambiguous names
include their coordinate alphabet, such as `e2_y1_y2` or
`e1_squared_y1_y2_y3_y4_y5`. Linear generators use
`y1`, `y4` or `omega7`, as appropriate. Euler products use `euler`; the E6
invariants in E7/P7 use `psi2`, `psi5`, `psi6`, `psi8`, `psi9` and `psi12`.
List the available names before selecting one. An ambiguous shorthand is not
silently assigned to a different Levi factor.

The [complete name and degree index](EXCEPTIONAL_GENERATORS.md) lists every
catalogued input with its source location. `info.as_dict()` returns a mutable,
JSON-compatible copy of the immutable metadata. `b.catalogue_coverage` records
whether the named generator list is complete for the covered parabolics.

## Coverage

| Type | Catalogue scope |
| --- | --- |
| G2 | Both maximal parabolics; all four Part III rows |
| F4 | All four maximal parabolics; all 18 Part III rows |
| E6 | All six maximal parabolics; all 41 Part III rows and the omitted P6 Euler image |
| E7 | All seven maximal parabolics; all 55 Part III rows, including the complete degree-12 invariant |
| E8 | Nine dissertation rows on P1 and P2 only; partial generator lists |

There are 128 catalogue entries in total. The E8 catalogue does not claim to
generate every Levi invariant ring, even on P1/P2. Arbitrary invariant
polynomials on any E8 parabolic can still be evaluated through `BorelMap`.

## Normalization and corrections

Dynkin nodes, words and fundamental weights use standard Sage/Bourbaki
numbering. The catalogue uses the normalized Springer morphisms from the
trace-projection method. In particular, the six E6 coordinates use denominator
18 in Part III Lemma 7.1, correcting its printed denominator 6.

For an invariant polynomial `q`, the character-side generator is
`q(Theta_1(t),...,Theta_r(t))`. Its image under the research morphism is computed
as `beta(q(omega_1,...,omega_r))`. A Laurent expression on the torus is therefore
not the input to this polynomial Borel API. The catalogue supplies the latter
polynomial explicitly.

E7/P7 `psi` generators follow the journal's unscaled power-sum convention.
The arXiv version's divided inputs must not be mixed with the journal's outputs.
Correction metadata also records the four numerical journal discrepancies,
the missing plus sign in the E6/P3 fifth elementary row, the supplementary
E6/P6 Euler generator, and the E8/P2 row whose repeated `e3` label should be `e4`.
These remain proposed local corrections, rather than a published erratum.

## Sources and verification

For integral questions, keep each catalogue normalization explicit. Some
complex-algebra generators have fractional Schubert images: the G2/P1
quadratic, for example, maps to `-3/4*sigma_(2,1)`. Such a form is rejected by
`GeneratorAlgebra.integral_degree`, even though its rational Borel expansion is
valid. Rescaling defines a different chosen Z-form. See the
[integral guide](INTEGRAL_IMAGES.md) and
[literature audit](INTEGRAL_LITERATURE.md); no automatic identification of
exceptional integral polynomial-character intersections is implemented.

The polynomial inputs come from the previously audited journal and dissertation
data. Source identifiers, source hashes, normalization choices and corrections
are packaged with each catalogue entry. Source publications are
[Kumar--Xie, Part III](https://kumar.math.unc.edu/papers/kumar86.pdf)
([DOI](https://doi.org/10.1016/j.jpaa.2023.107539)) and Xie's dissertation,
printed pages 98--99 (PDF pages 106--107) for the E8 Schubert supplement.

Only polynomial definitions and metadata are included in the runtime data file.
Reference Schubert coefficients live in the test fixtures. The exceptional
suite recomputes all entries, checks their Levi invariance, degree and Schubert
indices, and compares the entire expansion against the frozen prior audit.
Selected inputs are also checked by a separate root/coroot Chevalley product
calculation. The prior audit used an independently implemented BGG evaluator
and cross-checked the corrected numerical rows by Chevalley multiplication.

```bash
PYTHONPATH=src sage -python -m unittest tests.test_exceptional -v
```

The recorded execution details are in
[`exceptional-verification.json`](exceptional-verification.json).

The recorded run passed all five test methods: 128 complete expansions with
356 nonzero coefficient positions, 128 independently reconstructed input
polynomials, 661 ambient-coordinate Levi reflection checks and eight independent
Chevalley comparisons. The complete E7/P7 `psi12` input has 4,893 fundamental-weight
monomials and gives three Schubert terms; its recorded expansion took 0.91 seconds
after input construction, with maximum frontier size three. Timings depend on
the machine and input.
