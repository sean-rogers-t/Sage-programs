# Research articles and Springer-morphism verification

This collection accompanies the Sage calculations in this repository. The four
original article PDFs are preserved byte for byte. The verification note is a
separate local review prepared on 16 October 2026; its proposed corrections are
not a published erratum and do not modify the source articles.

## Articles

| Article | Local PDF | Published citation and original source |
| --- | --- | --- |
| Sean Rogers, *An explicit determination of the Springer morphism* | [PREPRINT, arXiv:1701.01538v1](Rogers_Springer_morphism_PREPRINT_arXiv_1701.01538v1.pdf) | *Communications in Algebra* **46** (2018), 4233–4242. [DOI](https://doi.org/10.1080/00927872.2017.1407418). The included PDF is the [6 January 2017 arXiv preprint](https://arxiv.org/pdf/1701.01538v1), not the final journal version. |
| Shrawan Kumar, *Representation ring of Levi subgroups versus cohomology ring of flag varieties* | [Part I](Kumar_Levi_subgroups_flag_varieties_I_MathAnn_2016.pdf) | *Mathematische Annalen* **366** (2016), 395–415. [DOI](https://doi.org/10.1007/s00208-015-1325-6). [Author's original PDF](https://kumar.math.unc.edu/papers/kumar76.pdf). |
| Shrawan Kumar and Sean Rogers, *Representation ring of Levi subgroups versus cohomology ring of flag varieties II* | [Part II](Kumar_Rogers_Levi_subgroups_flag_varieties_II_JAlgebra_2020.pdf) | *Journal of Algebra* **556** (2020), 340–362. [DOI](https://doi.org/10.1016/j.jalgebra.2020.02.029). [Author's original PDF](https://kumar.math.unc.edu/papers/kumar81.pdf). |
| Shrawan Kumar and Jiale Xie, *Representation ring of Levi subgroups versus cohomology ring of flag varieties III* | [Part III](Kumar_Xie_Levi_subgroups_flag_varieties_III_JPAA_2024.pdf) | *Journal of Pure and Applied Algebra* **228** (2024), 107539. [DOI](https://doi.org/10.1016/j.jpaa.2023.107539). [Author's original PDF](https://kumar.math.unc.edu/papers/kumar86.pdf). |

## Verification note for Shrawan

- [PDF ready to send](Springer_morphisms_for_Shrawan.pdf)
- [Editable Word document](Springer_morphisms_for_Shrawan.docx)
- [Standalone exact Sage calculation](calculations/reproduce_springer.py)

The same note is also available in Google Drive as a
[PDF](https://drive.google.com/file/d/1sxtIN2PXOAt3zDO6B4j-mW8O97xWyVE3/view)
and an [editable Word file](https://docs.google.com/document/d/1_VLNE71kTONKKvXP7XZT1W3LkXDhQ7LM/edit).

The note gives all Springer coordinates for the minimal-Dynkin-index
representations F4/omega4, E6/omega1 and E7/omega7, the weight matrices and their
inverses, independent exact calculations, and proposed corrections to Part III.
The weight-matrix method is Sean Rogers's method from Theorems 1 and 2 of
*An explicit determination of the Springer morphism*. Here, "independent" means
that the representation weights and calculations were regenerated separately
from the published exceptional tables and Jiale's calculation code; it does not
claim a new mathematical method.

The coordinate checks agree with the printed F4 and E7 formulas. For E6 the
proposed correction is to replace all six denominators 6 in Lemma 7.1 by 18,
leaving their numerators unchanged. The note also records separate proposed
Schubert-table corrections. Consult the note for conventions, exact locations,
and supporting calculations.

Run the self-contained reproduction with SageMath (tested with SageMath 9.5):

```sh
sage -python articles/calculations/reproduce_springer.py
```

It uses exact rational and Laurent-polynomial arithmetic, builds weights by Weyl
orbits, checks multiplicities against irreducible characters, and verifies that
the resulting morphisms vanish at the identity and have derivative equal to the
identity. It imports no published coordinate formula or historical calculation.

## Combinatorial proofs in types A, C, and D

The [combinatorial proof note](combinatorial-proofs/README.md) includes the
completed type C dissertation argument, type A and D extensions, the Euler
generator, stable injectivity arguments, and a comparison with Parts II and
III. The published classical main theorems already have proofs; this note
provides alternative arguments and explicit Schubert refinements.

- [LaTeX-generated PDF](combinatorial-proofs/Combinatorial_proofs_types_A_C_D.pdf)
- [Editable LaTeX](combinatorial-proofs/Combinatorial_proofs_types_A_C_D.tex)
- [Exact Sage checks and recorded results](combinatorial-proofs/verification/README.txt)

## Integral image summary for Shrawan

The [integral Springer–Borel image note](orthogonal-images-summary/README.md)
connects the orthogonal Smith calculations to Parts II and III. It labels
published inputs, deductions in the companion notes, and conjectures
separately, with normalization conventions, proof sketches, references, and
the integral comparison problem at its baseline. It makes no originality claim.
Subsequent work supplies a [general type-B proof](../borel/TYPE_B_GENERAL_SMITH.md);
the PDF remains the labeled earlier snapshot.

- [Five-page LaTeX PDF](orthogonal-images-summary/output/pdf/Springer_integral_images_for_Shrawan.pdf)
- [Editable LaTeX source](orthogonal-images-summary/main.tex)
- [Full proof and source status](../borel/PROOF_STATUS.md)

[PROVENANCE.json](PROVENANCE.json) records source URLs, SHA-256 hashes and file
sizes. The current PDF and Word note supersede the earlier local note filenames;
only this final pair is included here.
