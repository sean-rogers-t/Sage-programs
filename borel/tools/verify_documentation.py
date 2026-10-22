"""Run the trusted repository's Python documentation examples inside Sage.

From borel/: PYTHONPATH=src sage -python tools/verify_documentation.py
Python fences run in reading order, with a fresh namespace per document.
The checker executes the repository's Python fences. Shell snippets and
external links are outside its scope.
"""

import contextlib
import io
import json
from pathlib import Path
import re

from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[1]
DOCUMENTS = ("../README.md", "../articles/research-program-proposal/README.md",
             "README.md", "RESEARCH_LEDGER.md", "AGENTS.md",
             "LOCAL_IMAGE_BOUND.md", "EQUIVARIANT_IMAGE_PILOT.md", "LOCAL_IMAGE_AUDIT.md",
             "E6_ROOT_WEIGHT_COMPARISON.md", "E6_ROOT_WEIGHT_AUDIT.md",
             "E6_REMAINING_CASES.md", "E6_REMAINING_AUDIT.md", "E6_REMAINING_TABLES.md", "E6_RELATIVE_IMAGE_API.md",
             "E7_ROOT_WEIGHT_COMPARISON.md", "E7_ROOT_WEIGHT_AUDIT.md",
             "E7_REFINEMENT_MULTIPLICITIES.md", "E7_FULL_FLAG_COMPARISON.md", "E7_D5_REFINEMENT.md", "E7_D4_REFINEMENT.md", "E7_REMAINING_PARABOLICS.md", "E7_D5_A1_COMPARISON.md",
             "E7_LOCAL_SOURCES.md", "E7_A5_FAMILY.md", "E7_D4_A1_COMPARISON.md", "E7_ALL_PARABOLICS.md", "E7_RELATIVE_IMAGE_API.md", "E7_QUADRATIC_DERIVATION.md", "E7_DEGREE11_LIFTING.md",
             "MATHEMATICS.md", "METHODS.md", "SCHUBERT_CALCULUS.md", "PRESENTATIONS.md",
             "ORTHOGONAL_API.md", "EXTENSIONS.md", "EXCEPTIONAL_INTEGRAL.md", "EXCEPTIONAL_INTEGRAL_TABLES.md",
             "EXCEPTIONAL_MINUSCULE.md", "EXCEPTIONAL_MINUSCULE_TABLES.md", "E6_CHARACTER_FORM.md", "E6_CHARACTER_TABLES.md",
             "E7_CHARACTER_FORM.md", "E7_CHARACTER_TABLES.md",
             "CHARACTER_SOURCES.md", "G2_CHARACTER_FORM.md", "G2_FULL_FLAG.md", "F4_CHARACTER_FORM.md",
             "F4_CHARACTER_TABLES.md", "F4_REPRESENTATION_COMPARISON.md",
             "MINUSCULE_CHARACTER_API.md", "E6_REPRESENTATION_COMPARISON.md",
             "E7_REPRESENTATION_COMPARISON.md",
             "TRACE_SOURCE_THEOREM.md", "TRACE_SOURCE_AUDIT.md",
             "IRREDUCIBLE_CHARACTER_LATTICES.md", "IRREDUCIBLE_LATTICE_AUDIT.md",
             "SPRINGER_GEOMETRY.md",
             "E6_GEOMETRY.md",
             "E7_GEOMETRY.md",
             "EXCEPTIONAL_DEGREES.md",
             "INTEGRAL_COMPARISONS.md", "TWO_STEP_INTEGRAL.md",
             "TYPE_C_JOINT_PRIMITIVITY.md", "TYPE_C_JOINT_AUDIT.md", "ORTHOGONAL_JOINT_IMAGES.md",
             "IMAGE_RELATIONS.md", "INTEGRAL_IMAGES.md", "INTEGRAL_LITERATURE.md",
             "ORTHOGONAL_SMITH.md", "TYPE_B_SMITH.md", "TYPE_D_SMITH.md",
             "TYPE_D_FROBENIUS.md", "TYPE_D_SQUARE_COMPARISON.md", "TYPE_D_CANDIDATE.md",
             "TYPE_D_TWO_ROW_SMITH.md", "TYPE_D_SQUARE_THEOREM.md",
             "TYPE_D_EVEN_GLUING_THEOREM.md", "TYPE_D_EVEN_GLUING_AUDIT.md",
             "TYPE_D_RELATIVE_EULER_EXTENSION.md", "TYPE_D_HEAD_EULER_B_COMPARISON.md",
             "TYPE_D_TWO_HEAD_GLUING.md", "TYPE_D_TWO_TAIL_GLUING.md",
             "TYPE_D_FULL_EXPONENT_BOUND.md", "TYPE_D_GLUING_FILTRATION.md",
             "TYPE_B_FILTRATION.md", "TYPE_B_TWO_ROW_SMITH.md",
             "TYPE_B_THREE_ROW_UNITS.md", "TYPE_B_THREE_ROW_SMITH.md",
             "TYPE_B_HIGHEST_SMITH.md", "TYPE_B_UNIT_SMITH.md",
             "TYPE_B_CLOSED_CONJECTURE.md", "TYPE_B_GENERAL_SMITH.md",
             "TYPE_B_GENERAL_AUDIT.md",
             "TYPE_B_HIGH_LAYERS.md", "TYPE_B_DEFECT_COMPARISON.md", "PROOF_STATUS.md",
             "experiments/TWO_STEP_DERIVATION.md")


def main():
    reports = []
    links = 0
    for name in DOCUMENTS:
        path = ROOT / name
        source = path.read_text()
        namespace = {"__name__": "documentation_example"}
        blocks = list(re.finditer(r"^```python[ \t]*\n(.*?)^```[ \t]*$",
                                  source, flags=re.MULTILINE | re.DOTALL))
        for number, block in enumerate(blocks, 1):
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(block.group(1), str(path) + ":block" + str(number),
                             "exec"), namespace)
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", source):
            if "://" in target:
                continue
            filename, _, anchor = target.partition("#")
            linked = path.parent / filename if filename else path
            assert linked.exists(), (name, target)
            if anchor and linked.is_file() and linked.suffix == ".md":
                headings = re.findall(r"^#+\s+(.+)$", linked.read_text(), flags=re.MULTILINE)
                anchors = {re.sub(r"[^\w\- ]", "", title.lower()).replace(" ", "-")
                           for title in headings}
                assert anchor in anchors, (name, target)
            links += 1
        reports.append({"file": name, "python_blocks_executed": len(blocks)})
    report = {
        "result": "passed", "sage_version": SAGE_VERSION,
        "command": "PYTHONPATH=src sage -python tools/verify_documentation.py",
        "documents": reports, "local_links_and_anchors_checked": links,
        "python_blocks_executed": sum(r["python_blocks_executed"] for r in reports),
        "scope": "Python fences executed in reading order in a fresh namespace per document; local links and anchors checked. Shell snippets and external URLs were not automatically executed or crawled.",
    }
    (ROOT / "documentation-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
