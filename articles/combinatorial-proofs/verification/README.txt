Verification for Combinatorial_proofs_types_A_C_D.tex
===================================================

The scripts in this directory are byte-identical copies of the independently
written research checks reported in the note. The three main checks were
rerun from this repository layout before publication, and their resulting
JSON files are included. The older Theorem 5.5 sign audit is preserved with
its original results. No historical Springer-calculation code is imported.
The directory structure is preserved so imports work without changes.

Requirement: SageMath (the recorded runs used SageMath 9.5 on Ubuntu/WSL).
These scripts require Sage's Python, rather than ordinary CPython.

From this directory, run:

    sage -python tmp/proof-acd-a/check_type_a.py
    sage -python tmp/proof-acd-d/check_d.py
    sage -python tmp/proof-acd-check/verify.py

Each script writes its JSON results beside itself. The recorded results are:

    tmp/proof-acd-a/type_a_checks.json
    tmp/proof-acd-d/d_results.json
    tmp/proof-acd-check/verification.json

Dependencies included here:

    papers/table_audit/recompute.py
        Root-data BGG calculator, using exact rational arithmetic.
        Imported as a library by the scripts above; its unrelated standalone
        full-table audit is not required for this note.

    tmp/dissertation-proof-audit/check_theorem55.py
        Type C partition-label helper imported by verify.py, plus the earlier
        dissertation sign audit. Its original JSON results are also included.

The combined verification checks every Schubert coefficient at the indicated
degree and parabolic, not only a proposed list of nonzero coefficients.
The finite calculations corroborate the uniform proofs in the note; they
are not substitutes for those proofs.

To rebuild the PDF, compile the standalone .tex in the parent directory with a current
LaTeX installation (for example, latexmk -pdf or tectonic). No external .bib,
image, or data file is needed. The delivered PDF was generated with Tectonic
0.17.0, and the source also compiled successfully in Codex's LaTeX editor.
