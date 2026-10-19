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
DOCUMENTS = ("README.md", "MATHEMATICS.md", "METHODS.md", "SCHUBERT_CALCULUS.md", "PRESENTATIONS.md",
             "IMAGE_RELATIONS.md", "experiments/TWO_STEP_DERIVATION.md")


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
