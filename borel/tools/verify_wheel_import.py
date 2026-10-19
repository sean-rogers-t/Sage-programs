"""Assert that a Sage process imports the intended isolated wheel installation."""

import argparse
import json
from pathlib import Path

import sage_borel


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("expected_root", type=Path)
parser.add_argument("expected_version")
arguments = parser.parse_args()
source = Path(sage_borel.__file__).resolve()
assert arguments.expected_root.resolve() in source.parents, source
assert sage_borel.__version__ == arguments.expected_version, sage_borel.__version__
print(json.dumps({"version": sage_borel.__version__, "imported_from": str(source),
                  "isolated_import_verified": True}))
