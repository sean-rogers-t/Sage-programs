# Springer–Borel images and generalized Cayley maps

A 31-page research statement for mathematical collaborators, dated 29 October 2026.

- [Research statement PDF](output/pdf/Springer_Borel_Cayley_Research_Statement.pdf)
- [LaTeX main source](main.tex)
- [Source and artifact verification](verification.json)
- [Claim inventory](claim-inventory.json)
- [Rendered-page review](visual-review.json)

The statement covers published foundations, integral torus and Levi sources,
ordinary classical and exceptional image results, adapted bases and saturation,
equivariant and K-theory specialization, quantum/oriented/affine extensions,
and the generalized Cayley geometry program. It ends with prioritized structural
proof targets and a result-status inventory.

The mathematical baseline is proof commit
`9f6c000a4a97a5f7540424ba1f8f2881cc55a2ae`. Working-note links are pinned to
that commit. Published results, internally proved deductions, computational
dependencies, and conjectures are distinguished. This document does not claim
external review or establish novelty, and does not replace the full proofs.
The prior proposal and expositions are preserved as historical documents.

## Sources

`main.tex` includes `foundations.tex`, `integral.tex`, `theories.tex`,
`geometry.tex`, `program.tex`, and `references.tex`. Edit these sources together.
The PDF uses single-page viewing and ordinary letter pages, not an imposed
four-page sheet layout.

## Build and audit

With Tectonic available:

```powershell
New-Item -ItemType Directory -Force output/pdf | Out-Null
tectonic --keep-logs --outdir output/pdf main.tex
Copy-Item output/pdf/main.pdf output/pdf/Springer_Borel_Cayley_Research_Statement.pdf
python verify_statement.py
```

`verify_statement.py` checks local proof-note targets, claim-status assignments,
citation keys, source/PDF hashes, page count, and compiler diagnostics. These
are editorial and provenance checks, not a fresh mathematical proof audit.
Render the PDF and inspect its pages after substantive layout changes; a clean
text extraction alone does not verify visual quality. Temporary renders and
compiler intermediates are not committed.
