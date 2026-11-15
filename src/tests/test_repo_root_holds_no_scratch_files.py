"""TEST.745: the repository root holds no scratch script. `m.py`, a 19-line measuring helper, was committed there by accident in the GEN.1309 A commit."""

import os
import subprocess

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
ALLOWED_ROOT_PYTHON = {"setup.py"}
"""set: The Python files the root may hold. A script belongs under `scripts/`, a test under `src/tests/`."""


def test_the_root_holds_only_the_python_files_it_is_meant_to():
    tracked = subprocess.run(["git", "ls-files", "--", "*.py"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    stray = sorted(path for path in tracked if "/" not in path and path not in ALLOWED_ROOT_PYTHON)
    assert not stray, f"scratch Python files in the repository root (move them under scripts/ or src/tests/, or delete them): {stray}"
