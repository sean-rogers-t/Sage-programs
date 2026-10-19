# Integral Springer-Borel images: note for Shrawan

This nine-page LaTeX note connects the selected integral image problem
to Kumar-Rogers, Part II, and Kumar-Xie, Part III. It now includes the
general normalized type-B Smith multiplicity formula and its integral
argument, replacing the earlier snapshot's conjectural status.
Published inputs, repository deductions, and remaining type-D questions
are labeled separately. Originality and external expert review remain
unestablished.

The linked PDF remains the audited type-B-focused snapshot below. The
subsequent [general type-D square proof](../../borel/TYPE_D_SQUARE_THEOREM.md)
now establishes all square-sector factors, including smaller block three.
The subsequent [general even-tail gluing proof](../../borel/TYPE_D_EVEN_GLUING_THEOREM.md)
also establishes the full D formula in every rank, with every saturated
gluing direction at unit level. These later deductions and their
[integrated audit](../../borel/TYPE_D_EVEN_GLUING_AUDIT.md) are not yet
incorporated into that PDF.

- [PDF prepared for Shrawan](output/pdf/Springer_integral_images_for_Shrawan.pdf)
- [Editable, standalone LaTeX](main.tex)
- [Full general type-B proof](../../borel/TYPE_B_GENERAL_SMITH.md)
- [Fresh dependency-by-dependency proof audit](../../borel/TYPE_B_GENERAL_AUDIT.md)
- [Claim-by-claim source and proof status](../../borel/PROOF_STATUS.md)
- [Articles and earlier verification notes](../README.md)

The note states the Gaussian and beta-number formulas, explains the
integral Schur-prefix reduction, constructs Chow-envelope lifts of
unbalanced incidence rows, proves a unit triangular minor, counts its
pivots, and derives every Smith multiplicity by an integral rank squeeze.
It also records the pre-existing counting and Demazure-character
interpretations. At that PDF's baseline the general D square comparison
was unresolved; its later proof is linked above. The proved B formula
also determines the D Euler-sector factors via the earlier comparison.

The fresh audit found no mathematical gap or correction in the general
proof or its earlier prefix reduction. It explicitly distinguishes the
topological quaternionic construction and its Leray-Hirsch basis from
the cited algebraic analogue. Finite checks corroborate the formulas;
they do not replace the arbitrary-rank proof.

Updated on 19 October 2026 against the audited mathematical baseline
`5e718334a0f8ed9b31166aa3bf42185a60d05c34`. Its full-proof link is immutable;
the PDF also links the current audit and status map on the research
branch. The original articles and historical Sage programs are preserved.

## Build and verification

`main.tex` requires no external bibliography, images, or data. With Tectonic
on the command line, run the following in this directory (PowerShell):

```powershell
New-Item -ItemType Directory -Force output/pdf | Out-Null
tectonic --outdir output/pdf main.tex
Move-Item -Force output/pdf/main.pdf output/pdf/Springer_integral_images_for_Shrawan.pdf
```

The delivered PDF was built with Tectonic 0.17.0 and also compiled
successfully in the desktop LaTeX editor. All nine final pages were
rendered with Poppler and visually inspected. Glyph bounds, bibliography
references, external link annotations, and the TeX log were checked.
No library runtime code or API changed for this document.
