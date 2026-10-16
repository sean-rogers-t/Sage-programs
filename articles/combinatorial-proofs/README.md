# Combinatorial proofs for classical Springer-induced maps

This research note completes the unfinished combinatorial argument for
Theorem 5.5 in Sean Rogers's dissertation and develops the corresponding
type A and type D arguments. It includes finite-rank generator images,
explicit Schubert expansions, the type D Euler generator, stable injectivity
proofs, source corrections, and independently computed exact checks.

- [Read the PDF](Combinatorial_proofs_types_A_C_D.pdf)
- [Editable, standalone LaTeX source](Combinatorial_proofs_types_A_C_D.tex)
- [Calculation scripts and recorded results](verification/README.txt)
- [File hashes and provenance](PROVENANCE.json)

The note was prepared for Sean Rogers on 16 October 2026. It is a research
review and alternative proof, rather than a published erratum or a claim of
new main theorems. The original article PDFs remain in [the articles index](../README.md).

## Comparison with the published papers

| Result | Published argument | Argument in this note |
| --- | --- | --- |
| Type C finite formula | Part II, Proposition 11 and Lemma 13, pp. 349-351: a tautological-bundle and total-Chern-class argument. The paper credits its shorter proof to L. Mihalcea. | Coinvariants, complete symmetric functions, exhaustive divided-difference paths, and Giambelli; gives the explicit alternating two-row Schubert expansion. |
| Type C stable injectivity | Part II, Theorem 16, pp. 353-354: an integral presentation, reduction modulo 2, and a degreewise free-kernel argument. | Specialize `c_p` to the elementary symmetric polynomial `e_p(y_1,...,y_k)`; the generator images become the algebraically independent `e_i(y_1^2,...,y_k^2)`. |
| Type D finite formula | Part III, Lemmas 3.4-3.5 and Theorem 3.6, pp. 9-11: the analogous bundle argument. | Divided differences at the fork explain the two middle classes and factors of two; Giambelli gives the typed two-row expansion, and a separate calculation gives the Euler image. |
| Type D stable injectivity | Part III, Theorem 3.8, pp. 14-15: first proves injectivity on the square-symmetric subring modulo 2, then handles the Euler generator by parity. | Specialization replaces the first stage. The final argument from `f^2 = e_k g^2` is the same parity argument already present in Part III, expressed using a prime valuation. |

The published proofs of these main classical results are already complete.
Part II already has the correct modulo-2 squares and zero relations.
Part III already has the correct fork half-sum weights and the valid
torsionfree scalar-extension statement. The dissertation errors identified
in the note should not be attributed to the published proofs.

The symmetric-function identities and Buch-Kresch-Tamvakis Giambelli formulas
are established inputs. This note supplies the completed dissertation
argument, explicit Schubert refinements, and a different stable injectivity
argument. Its type A calculation also identifies the sign needed to compare
the Springer-induced map with the positive-generator stable isomorphism.

For the type labels used in Part III, the Euler image is
`(-1)^k (tau_k - tau_k')`. The note also gives the explicit Weyl words, so its
formula can be interpreted without choosing which middle class is primed.

## Reproduce the calculations

Use SageMath, rather than ordinary Python. The recorded calculations used
SageMath 9.5 on Ubuntu/WSL. From the repository root:

```sh
sage -python articles/combinatorial-proofs/verification/tmp/proof-acd-a/check_type_a.py
sage -python articles/combinatorial-proofs/verification/tmp/proof-acd-d/check_d.py
sage -python articles/combinatorial-proofs/verification/tmp/proof-acd-check/verify.py
```

The archived directory layout is preserved so every import resolves within
`verification/`. No historical calculation routine, external source article,
or unpublished input file is required. Each script writes its JSON output
beside itself. The `tmp/` names preserve the original calculation paths;
these checked-in files are part of the reproducible evidence.

All recorded main checks passed with exact rational arithmetic:

| Check | Recorded cases |
| --- | ---: |
| Type A head and tail generators, ranks 2-8 | 168 |
| Type A Schur images in degrees 0-6, including 145 vanishing cases | 251 |
| Type C head and squared-tail formulas, ranks 2-8 | 203 |
| Type D squared-tail formulas, ranks 4-8 | 75 |
| Type D special classes and Euler images, ranks 4-6 | 84 |

These groups include overlapping auxiliary checks and are not a total of
distinct mathematical cases. The combined C and D audit compares every
Schubert coefficient in the relevant degree, including all proposed zero
coefficients. The older `theorem55_results.json` also records intentional
negative checks showing failures of a literal reading of the dissertation's
printed signs. Those failures are evidence for the corrections, not failures
of the corrected formulas. Finite calculations corroborate the uniform
proofs and do not replace them.

## Build the PDF

The `.tex` is standalone and needs no external bibliography, image, or data
file. From this directory, use a current LaTeX installation, for example:

```sh
latexmk -pdf Combinatorial_proofs_types_A_C_D.tex
# Alternatively:
tectonic Combinatorial_proofs_types_A_C_D.tex
```

The included 14-page PDF was compiled with Tectonic 0.17.0. The source also
compiled in the Codex LaTeX editor. The rendered pages were visually checked.

## Sources

- Sean Rogers, *The Springer Morphism, Polynomial Representation Rings, and
  the Cohomology Ring of Grassmannians*, UNC Chapel Hill dissertation, 2018.
  Theorem 5.5, printed pp. 70-74; type A results, pp. 54-57; stable type C
  discussion, pp. 75-76.
- Jiale Xie, *Polynomial Representation Ring versus Cohomology Ring of Flag
  Varieties for Classical and Exceptional Groups*, UNC Chapel Hill dissertation,
  2022. Proposition 4.6 and Theorem 4.8, printed pp. 35-39.
- Kumar and Rogers, [Part II](../Kumar_Rogers_Levi_subgroups_flag_varieties_II_JAlgebra_2020.pdf),
  *Journal of Algebra* 556 (2020), 340-362.
  [DOI](https://doi.org/10.1016/j.jalgebra.2020.02.029).
- Kumar and Xie, [Part III](../Kumar_Xie_Levi_subgroups_flag_varieties_III_JPAA_2024.pdf),
  *Journal of Pure and Applied Algebra* 228 (2024), 107539.
  [DOI](https://doi.org/10.1016/j.jpaa.2023.107539).
- Buch, Kresch, and Tamvakis, *A Giambelli formula for isotropic Grassmannians*,
  [author preprint](https://www.math.umd.edu/~harryt/papers/isogiamfinal.pdf),
  and *A Giambelli formula for even orthogonal Grassmannians*,
  [author preprint](https://www.math.umd.edu/~harryt/papers/oggiam.pdf).

The note's bibliography and locator table distinguish printed dissertation
pages from the page numbers in the supplied PDF copies.
