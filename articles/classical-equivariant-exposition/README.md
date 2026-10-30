# Springer-Borel images in classical type: expository paper

A twelve-page exposition, prepared 22 October 2026, connecting the notation
and maps of Kumar, Kumar-Rogers and Kumar-Xie to the current classical
equivariant program. This is an expository research draft, not a new
publication or a replacement for the companion proofs.

Subsequent progress: the [classical equivariant API](../../borel/CLASSICAL_EQUIVARIANT.md)
now packages the first type A maps and the symplectic projective model.
The [even/odd cokernel proof](../../borel/SYMPLECTIC_K_FITTING.md)
now proves a uniform doubled half-size presentation and complete closed
rank-four Fitting formulas. The PDF remains the dated
exposition of its starting point.

- [Read the compiled PDF](output/pdf/Classical_Springer_Borel_equivariant_exposition.pdf)
- [Editable LaTeX source](main.tex)
- [PDF and source verification](verification.json)
- [Selected mathematical example checks](examples-verification.json)

The PDF contains one full-size manuscript page per PDF page, with an
explicit SinglePage viewer setting. Four-page contact sheets in the
temporary QA folder are inspection aids, not the article PDF.

## Contents

The paper explains groups, weights, root and character lattices, retained
and omitted parabolic nodes, the parameters n and k, Bruhat decomposition,
opposite Schubert classes and their degrees. It presents representation
rings and the rational Borel quotient, with an integral Grassmannian
example, then develops the trace-projection Springer morphism and the
published polynomial-character maps.

The integral section states the normalized type C saturation result, the
complete type B partition formula and the type D Euler/gluing formula.
It distinguishes these forms from literal integral character sources and
indices in the selected saturation from full ambient cokernels.

The remaining sections introduce equivariant cohomology, ordinary and
equivariant K-theory, specialization defects, the finite filtered image
comparison, the rational deformed source and determinant theorem, and
the completed rational type A tensor-power and integral symplectic
projective families. The concluding program targets classical integral
images and higher Fitting ideals.

Published foundations, internally audited repository deductions and
proposed extensions are labeled separately. Originality and external
expert review are not asserted. The bibliography points to primary
literature and companion notes pinned to commit
2bc20702f142d6124a22ab97aaa511caf227f916.

## Build and checks

With Tectonic available, run from this directory:

    New-Item -ItemType Directory -Force tmp/pdfs,output/pdf | Out-Null
    tectonic --keep-logs --outdir tmp/pdfs main.tex
    pdftoppm -r 65 -png tmp/pdfs/main.pdf tmp/pdfs/page
    python verify_pdf.py
    Copy-Item tmp/pdfs/main.pdf output/pdf/Classical_Springer_Borel_equivariant_exposition.pdf

The PDF checker uses pypdf, pdfplumber and Pillow. If changing the page
count, remove obsolete page renderings before running it. It checks
43 equation tags, resolved citations, local companion-note targets,
page geometry and the absence of overfull boxes or undefined references.
The delivered PDF was compiled with Tectonic 0.17.0 and checked with the
desktop LaTeX compiler. All twelve rendered pages were visually reviewed.

From the repository's borel directory, the selected example checks run as:

    PYTHONPATH=src sage -python ../articles/classical-equivariant-exposition/verify_examples.py

These compare actual BGG Schubert matrices with the closed B/D formulas,
reconstruct the SL4 root-difference image, and compute the Sp4 projective
K inclusion directly in Z[v]/(v^4). They check the examples, not every
general proof in the paper. The mathematical source and audit records
remain [the proof-status guide](../../borel/PROOF_STATUS.md) and
[the research ledger](../../borel/RESEARCH_LEDGER.md).
