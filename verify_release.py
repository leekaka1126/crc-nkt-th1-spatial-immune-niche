#!/usr/bin/env python
"""Verify every frozen package file against SHA256_MANIFEST.json."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    manifest = json.loads((ROOT / "SHA256_MANIFEST.json").read_text(encoding="utf-8"))
    expected = manifest["files"]
    actual = {path.relative_to(ROOT).as_posix(): path for path in ROOT.rglob("*")
              if path.is_file() and path.name != "SHA256_MANIFEST.json"
              and "__pycache__" not in path.parts and ".venv" not in path.parts}
    if set(actual) != set(expected):
        raise AssertionError(f"File set changed: missing={sorted(set(expected)-set(actual))}, "
                             f"extra={sorted(set(actual)-set(expected))}")
    changed = [name for name, path in actual.items()
               if path.stat().st_size != expected[name]["bytes"]
               or sha256(path) != expected[name]["sha256"]]
    if changed:
        raise AssertionError(f"Hash mismatch: {changed}")
    print(f"PASS: {len(actual)} files match {manifest['release']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
