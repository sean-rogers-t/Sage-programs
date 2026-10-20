"""Record formula-only timings; no Schubert or Weyl construction is performed."""

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.env import SAGE_VERSION
from sage_borel import __version__, orthogonal_smith


def main():
    root = Path(__file__).resolve().parents[1]
    files = ("src/sage_borel/orthogonal.py", "src/sage_borel/operations.py",
             "benchmarks/orthogonal.py")
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sage_version": SAGE_VERSION,
        "library_version": __version__,
        "scope": "formula evaluation only; Sage import time excluded; no claim of a universal time bound",
        "source_sha256": {name: sha256((root / name).read_bytes()).hexdigest()
                          for name in files},
        "cases": [],
    }
    for family, n, k, degree in (("D", 9, 4, 14), ("D", 9, 4, 16),
                                ("B", 40, 15, 200), ("D", 40, 16, 200),
                                ("B", 200, 80, 1000), ("D", 200, 80, 1000)):
        start = perf_counter()
        result = orthogonal_smith(family, n, k, degree)
        elapsed = perf_counter() - start
        report["cases"].append({**result.as_dict(), "seconds": round(elapsed, 9)})
    destination = root / "benchmarks" / "orthogonal-results.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
