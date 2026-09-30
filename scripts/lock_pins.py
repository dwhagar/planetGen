#!/usr/bin/env python3
"""
Reads requirements.lock (written by scripts/lock-requirements.sh) for
scripts/install-python-deps.sh. Standard library only, so it runs under
the bare system Python before anything is installed.

    lock_pins.py constraints LOCK
        Prints the lock as a pip constraints file: one "name==version"
        line per pin, with its environment marker and without its hashes
        (pip resolves against these, so every version it picks is a
        locked one).

    lock_pins.py resolve LOCK OUT [PIP_FLAG ...] -- REQUIREMENT ...
        Asks this Python's pip what installing REQUIREMENTs would add
        (`pip install --dry-run --report`), holding everything it adds
        to the locked version, and writes those pins with their hashes
        from the lock to OUT, a requirements file for `pip install
        --no-deps --require-hashes`. A library already installed is left
        out of the constraints, so a copy that already satisfies what
        needs it (say, apt's) stays in use instead of being replaced by
        the locked version; one pip must replace anyway is held to the
        lock on a second pass. Exits 1 if pip picks anything the lock
        has no hashes for, 3 if this pip has no --report (older than
        22.2), and pip's own status if it can't resolve.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
from importlib import metadata


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


def constraints(lock):
    for name, version, marker, _ in read_lock(lock):
        print(f"{name}=={version}" + (f" ; {marker}" if marker else ""))


def _installed(name):
    for candidate in (name, name.replace("-", "_")):
        try:
            metadata.distribution(candidate)
            return True
        except metadata.PackageNotFoundError:
            pass
    return False


def _dry_run(pip_flags, requirements, constraint_lines):
    """The (name, version) pairs pip would install, or an exit status."""
    with tempfile.TemporaryDirectory() as tmp:
        constraints = os.path.join(tmp, "constraints.txt")
        report = os.path.join(tmp, "report.json")
        with open(constraints, "w", encoding="utf-8") as f:
            f.write("\n".join(constraint_lines) + "\n")
        help_text = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--help"],
            capture_output=True, text=True,
        ).stdout
        if "--report" not in help_text:
            return 3
        status = subprocess.run(
            [sys.executable, "-m", "pip", "install", *pip_flags, "--dry-run", "--quiet",
             "--report", report, "-c", constraints, *requirements],
        ).returncode
        if status:
            return status
        with open(report, encoding="utf-8") as f:
            return [
                (normalize(item["metadata"]["name"]), item["metadata"]["version"])
                for item in json.load(f)["install"]
            ]


def resolve(lock, out, pip_flags, requirements):
    pins = read_lock(lock)
    hashes = {(name, version): h for name, version, _, h in pins}
    held = {name for name, _, _, _ in pins if not _installed(name)}
    for _ in range(3):
        lines = [
            f"{name}=={version}" + (f" ; {marker}" if marker else "")
            for name, version, marker, _ in pins if name in held
        ]
        chosen = _dry_run(pip_flags, requirements, lines)
        if isinstance(chosen, int):
            return chosen
        loose = {name for name, _ in chosen} - held
        if not loose:
            break
        held |= loose
    missing = [f"{name}=={version}" for name, version in chosen if not hashes.get((name, version))]
    if missing:
        print("error: requirements.lock has no hashes for " + " ".join(missing)
              + "; run scripts/lock-requirements.sh.", file=sys.stderr)
        return 1
    with open(out, "w", encoding="utf-8") as f:
        for name, version in chosen:
            f.write(f"{name}=={version} " + " ".join(f"--hash={h}" for h in hashes[(name, version)]) + "\n")
    return 0


def main(argv):
    if len(argv) == 3 and argv[1] == "constraints":
        constraints(argv[2])
        return 0
    if len(argv) >= 5 and argv[1] == "resolve" and "--" in argv[4:]:
        split = argv.index("--", 4)
        return resolve(argv[2], argv[3], argv[4:split], argv[split + 1:])
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
