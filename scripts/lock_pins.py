#!/usr/bin/env python3
"""
Reads requirements.lock and requirements-server.lock (written by
scripts/lock-requirements.sh): `read_lock(path)` returns each pin with its
marker and hashes. Standard library only. Imported by
src/tests/test_install_python_deps.py, which checks that the locks cover
setup.py's requirements and that every pin has hashes.
"""

import re


def normalize(name):
    """PEP 503 name: lower case, runs of -_. become one -."""
    return re.sub(r"[-_.]+", "-", name).lower()


def read_lock(path):
    """
    Returns [(name, version, marker, [hashes])] for each pin in the lock,
    in file order. A pin is one logical line (continuations joined):
    `name==version [; marker] --hash=sha256:... ...`.
    """
    with open(path, encoding="utf-8") as f:
        text = f.read()
    logical = re.sub(r"\\\n", " ", text)
    pins = []
    for line in logical.splitlines():
        line = line.split(" #", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        hashes = re.findall(r"--hash=(\S+)", line)
        requirement = re.sub(r"--hash=\S+", "", line).strip()
        match = re.fullmatch(r"([A-Za-z0-9][A-Za-z0-9._-]*)==([^\s;]+)\s*(?:;\s*(.+))?", requirement)
        if not match:
            raise SystemExit(f"error: can't read this line of {path}: {line}")
        name, version, marker = match.groups()
        pins.append((normalize(name), version, (marker or "").strip(), hashes))
    return pins
