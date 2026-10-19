# Integral Springer–Borel images: note for Shrawan

This five-page discussion note explains how the selected integral image problem
refines the maps and stable injectivity in Kumar–Rogers, Part II, and Kumar–Xie,
Part III. Published inputs, deductions in the companion notes, and conjectures
are labeled separately. It makes no claim of originality or external review.

- [PDF prepared for Shrawan](output/pdf/Springer_integral_images_for_Shrawan.pdf)
- [Editable, standalone LaTeX](main.tex)
- [Claim-by-claim proof and source audit](../../borel/PROOF_STATUS.md)
- [Articles and earlier verification notes](../README.md)

The note covers the normalized lattices, finite and stable B/D deductions,
the type-B unit and largest-factor formulas, the classical counting inputs,
the type-D Euler sector and gluing, and the precise unresolved integral
comparison behind the general multiplicity conjecture. Unit Smith factors
are retained along with the finite cokernel. The literal SO/Sp character
normalization is discussed separately from the normalized form.

Prepared on 19 October 2026 against research baseline
`0adbc65288539d5aef59b3dfeeb905b2420ce051`. The PDF is a summary with proof
sketches; the linked companion notes contain the full arguments and exact
finite verification records. The original articles are preserved separately.

## Build and verification

`main.tex` requires no external bibliography, images, or data. With Tectonic
on the command line, run the following in this directory (PowerShell):

```powershell
New-Item -ItemType Directory -Force output/pdf | Out-Null
tectonic --outdir output/pdf main.tex
Move-Item -Force output/pdf/main.pdf output/pdf/Springer_integral_images_for_Shrawan.pdf
```

The delivered PDF was built with Tectonic 0.17.0. The source also compiled
successfully in the Codex LaTeX editor. All five pages were rendered with
Poppler and visually inspected; extracted glyph bounds, bibliography
references, and the TeX log were checked. No library code changed for this
document, and finite computations are distinguished from uniform proofs.
