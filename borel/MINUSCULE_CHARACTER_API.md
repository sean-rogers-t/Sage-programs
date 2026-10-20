# E6/P1 and E7/P7: certified integral character APIs

Version 0.13 packages the previously proved and audited minimal character
sources and the E6 adjoint comparison.
Version 0.14 also adds the E7 adjoint source.

| Factory | Representation | Full torus source |
|---|---|---|
| `e6_character_source(1)` | V(omega1), dimension 27 | Sym_Z(6Q) |
| `e6_character_source(1, representation="adjoint")` | V(omega2), dimension 78 | Sym_Z(24P) |
| `e7_character_source(7)` | V(omega7), dimension 56 | Sym_Z(12Q) |
| `e7_character_source(7, representation="adjoint")` | V(omega1), dimension 133 | Sym_Z(36P) |

P is the full weight lattice of the simply connected group, and Q its root
lattice. The source for the indicated Levi is the W_L-invariant subring of
the displayed torus source. Numbering and generator normalizations are
unchanged from the catalogue and the Kumar papers.

These helpers support only E6/P1 and E7/P7. They reject other nodes rather
than imply that every exceptional parabolic has been certified. The
selectors for E6 and E7 accept exactly `"minimal"` and `"adjoint"`.

## Using the API

```python
from sage_borel import e6_character_source, e7_character_source

minimal = e6_character_source(1)
adjoint = e6_character_source(1, representation="adjoint")
e7 = e7_character_source(7)

source_degree = minimal.degree(2)
assert len(source_degree.source_basis) == 2
assert source_degree.image().smith_factors == (108,)
assert e7.degree(2).image().smith_factors == (144,)

a = minimal.degree(1).image()
b = adjoint.degree(1).image()
comparison = a.compare(b)
assert not comparison.left_in_right
assert not comparison.right_in_left
assert comparison.sum.smith_factors == (6,)
assert comparison.intersection.smith_factors == (72,)
assert comparison.sum.quotient_by(a).torsion_factors == (3,)
assert comparison.sum.quotient_by(b).torsion_factors == (4,)
```

`degree(d)` means polynomial degree d, corresponding to H^(2d).
Its `monomials` are the rational invariant monomials; columns of
`basis_matrix` express `source_basis` in those monomials. Source elements
with zero Borel image are retained. `image()` returns an embedded integral
Schubert lattice, supporting `contains`, `compare` and `quotient_by`.
See [lattice comparisons](INTEGRAL_COMPARISONS.md).

The same comparison method accepts the earlier catalogue images:

```python
from sage_borel import exceptional_algebra

catalogue = exceptional_algebra("E7", 7).integral_degree(5)
old_comparison = e7.degree(5).image().compare(catalogue)
assert not old_comparison.left_in_right
assert not old_comparison.right_in_left
```

This reproduces the earlier E7 degree-five incomparability witness.
It does not change the catalogue's integral normalization.

The formal variables, in order, are:

- E6: y1, e1, e2, euler, e3, e4, of degrees 1, 2, 4, 5, 6, 8.
- E7: omega7, psi2, psi5, psi6, psi8, psi9, psi12, of degrees
  1, 2, 5, 6, 8, 9, 12.

These are rational polynomial generators. The returned degreewise integer
bases can involve rational combinations of their monomials. Individual
generator rescaling alone need not produce the full character source.
The helpers do not claim a finite integral generator-and-relation
presentation.

## Mathematics and computation

The [E6 minimal proof](E6_CHARACTER_FORM.md) and
[E7 proof](E7_CHARACTER_FORM.md) establish the respective full torus
sources by algebraic independence of integral Laurent pullbacks modulo
every prime. The [E6 adjoint proof](E6_REPRESENTATION_COMPARISON.md)
uses a separate unit-leading-coefficient certificate.
The [E7 adjoint proof](E7_REPRESENTATION_COMPARISON.md) separately certifies
Sym_Z(36P), with all-Levi sum/intersection formulas.

For E6 the coefficient constraints are computed in a root basis for the
minimal representation and a weight basis for the adjoint. For E7 the
integral A5 change of coordinates gives Z[E1,...,E6,u] for Q and
Z[E1,...,E6,z] for P, with u=2z+E1. Here Ej is an elementary symmetric
polynomial in six permuted coordinates. This replaces large root-variable
expansions with a smaller weighted coefficient matrix. The integral
coordinate identification is proved in the E7 note; it is not an
approximation.

Exact Hermite reduction finds the preimage of the integer coefficient
lattice. The source basis is then multiplied by 6^d, 24^d, 12^d or 36^d as
appropriate. Every source basis still includes the Borel kernel.

Images use ordinary fixed-point interpolation: evaluate invariant
polynomials at Weyl translates of a regular point, invert the triangular
Schubert restriction matrix, and retain the coefficients of degree d.
Homogeneity makes these the ordinary Borel coefficients. Products are
evaluated pointwise, avoiding dense expansion of high-degree products.
This backend uses the same normalization as the direct BGG engine.
It is not a public equivariant calculation API.

The runtime imports no experimental scripts or saved image tables.
`rational_algebra` retains the usual formal-polynomial interface; an
explicit call to `substitute` can still create a large Cartan polynomial.

## Bounds and checks

The default source limits are 500 rational monomials and 2,000,000 possible
coefficient-matrix entries. The latter counts weighted monomials in the E7
compressed model. These are allocation checks, not runtime guarantees.
For a larger requested degree, raise limits explicitly:

For example, `e7.degree(27, max_matrix_entries=10000000)` requests the
top-degree source with a larger coefficient budget.

The complete [API replay](experiments/minuscule_character_api.py) checks:

- All 17 E6 minimal and 28 E7 minimal source lattices and embedded images
  against their historical audited tables.
- All 17 E6 adjoint degrees, with coefficient integrality and gcd-one
  maximal-minor certificates proving completeness of each source lattice.
- All 2,731 rational Borel columns against the separate
  [published-ring/Chevalley implementation](experiments/exceptional_published_rings.py).
  That implementation retains its documented BGG generator-image and
  E6 square-anchor dependencies; it shares no fixed-point interpolation.

The [Python-only audit](experiments/audit_e6_adjoint_api.py) independently
reconstructs the adjoint root system, torus certificate, source coefficient
constraints in all 17 degrees, eight direct Laurent intersections and all
image/comparison arithmetic. It compares the 45 minimal source lattices
with the previously audited tables rather than claiming to rederive all
historical proofs.

Results: [full API replay](experiments/minuscule-character-api-results.json),
[adjoint torus and Laurent checks](experiments/e6-adjoint-certificate-results.json),
[independent audit](experiments/e6-adjoint-api-audit.json).
Historical reports retain their original dependency hashes; new reports
bind the exact historical tables used.

The subsequent [E7 adjoint comparison](E7_REPRESENTATION_COMPARISON.md)
adds all 28 adjoint degrees and 1,967 fresh ordinary-ring/Borel-column
comparisons. Its Python-only audit reconstructs the integral P model and
compressed polynomials with a separate closed partition formula, all
source coefficient constraints, eight direct Laurent lattices and all
embedded comparisons. The previously certified minimal image tables
are reused, with six fresh minimal API degree checks.
