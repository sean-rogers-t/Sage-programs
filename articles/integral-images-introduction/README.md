# Integral Springer-Borel images: manuscript introduction

A three-page standalone introduction in the notation of Kumar-Rogers,
Part II, and Kumar-Xie, Part III. It summarizes the published starting
point, defines the normalized integral-image problem, and states the
general type B and type D formulas from the audited companion notes.

- [Compiled introduction](output/pdf/Springer_Borel_introduction.pdf)
- [Editable LaTeX](main.tex)
- [General B proof](../../borel/TYPE_B_GENERAL_SMITH.md)
- [General D proof](../../borel/TYPE_D_EVEN_GLUING_THEOREM.md)
- [D proof audit](../../borel/TYPE_D_EVEN_GLUING_AUDIT.md)

The notation follows the original papers: `n,k`, the varieties
`OG(n-k,2n+1)` and `OG(n-k,2n)`, the B map `bar xi_{n,k}`, the D map
`xi^D_{n,k}`, generators `e_i` and `bar h_{1,k}`, and quotient bundles
`Q^B` and `Q`. The type D range is restricted to `n >= 4` to match
Part III. New multiplicity notation is defined explicitly.

The variable `z` records half the cohomological degree. In particular,
the type D Euler term uses `z^k B_{n-1,k,e}(z)`. Relative to the longer
repository notes, the parameter conversion is `a=n-k`. Unit Smith
factors are retained and all indices concern saturation in the selected
rational image.

This is an unpublished manuscript draft, not a full proof paper.
The theorem statements follow the internally audited notes at commit
`b16f582e29947a402672d160b363d8c9a2f60c11`; originality and external
expert review remain unestablished. The earlier nine-page Shrawan note
has been preserved separately.

## Build

The source is self-contained and uses an inline bibliography. With
Tectonic installed, run from this directory:

```powershell
New-Item -ItemType Directory -Force output/pdf | Out-Null
tectonic --outdir output/pdf main.tex
Move-Item -Force output/pdf/main.pdf output/pdf/Springer_Borel_introduction.pdf
```

The delivered source compiled with Tectonic 0.17.0 and the desktop
LaTeX editor. All three rendered pages were visually checked; the
TeX log, page bounds, references, and PDF links were inspected.
