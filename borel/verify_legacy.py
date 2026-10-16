"""Verify that the 35 historical program files have not changed.

Run with ordinary Python or Sage: python verify_legacy.py
The --record option rebuilds the manifest only after comparison with the frozen
Git baseline; it refuses to record changed historical files as preserved.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
BASELINE = "7849117c562de02d11b388e0b6e3126e86f3474c"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--checkout-bytes", action="store_true",
                        help="also require the exact recorded local checkout bytes")
    args = parser.parse_args()
    repo = HERE.parent
    manifest_path = HERE / "legacy_manifest.json"
    if args.record:
        names = subprocess.check_output(
            ["git", "ls-tree", "--name-only", "-z", BASELINE], cwd=repo).decode().split("\0")
        names = sorted(name for name in names if Path(name).suffix in (".sage", ".py", ".txt"))
        if len(names) != 35:
            raise RuntimeError("expected exactly 35 historical programs")
        entries = []
        for name in names:
            original = subprocess.check_output(["git", "show", BASELINE + ":" + name], cwd=repo)
            current = (repo / name).read_bytes()
            git_blob = subprocess.check_output(
                ["git", "rev-parse", BASELINE + ":" + name], cwd=repo).decode().strip()
            expected_checkout = subprocess.check_output(
                ["git", "cat-file", "--filters", BASELINE + ":" + name], cwd=repo)
            if current not in (original, expected_checkout):
                raise RuntimeError("historical file differs from the baseline: " + name)
            entries.append({"path": name, "git_blob": git_blob,
                            "git_blob_bytes": len(original), "git_blob_sha256": sha256(original),
                            "checkout_bytes": len(current), "checkout_sha256": sha256(current)})
        manifest = {"git_baseline": BASELINE, "files": entries,
                    "historical_file_count": len(entries), "git_contents_preserved": True,
                    "hash_semantics": "Canonical Git blob and local Windows checkout recorded separately; Git autocrlf changes checkout line endings, not repository contents."}
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    else:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for entry in manifest["files"]:
            current = (repo / entry["path"]).read_bytes()
            original = subprocess.check_output(
                ["git", "cat-file", "blob", entry["git_blob"]], cwd=repo)
            expected_checkout = subprocess.check_output(
                ["git", "cat-file", "--filters", manifest["git_baseline"] + ":" + entry["path"]], cwd=repo)
            if sha256(original) != entry["git_blob_sha256"] or current not in (original, expected_checkout):
                raise RuntimeError("historical file changed: " + entry["path"])
            if args.checkout_bytes and (len(current) != entry["checkout_bytes"] or
                                        sha256(current) != entry["checkout_sha256"]):
                raise RuntimeError("local checkout bytes changed: " + entry["path"])
    print("All {} historical program files match the frozen baseline.".format(len(manifest["files"])))


if __name__ == "__main__":
    main()
