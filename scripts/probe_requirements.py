#!/usr/bin/env python3
"""
Checks planetGen's Python requirements as this interpreter sees them, for
scripts/install-python-deps.sh (Linux and macOS) and
scripts/deploy-common.ps1 (Windows). Standard library only.

    probe_requirements.py SPEC ...        (SPEC is name>=floor)

Prints "<state> <spec> <version or error> <directory>" per spec: state is
ok, missing (not installed), old (below its floor) or broken (installed
but the import raised). Imports each one for real, since an installed
distribution whose own dependencies are missing is not usable either.
"""
import importlib
import re
import sys
from importlib import metadata


def key(version):
    parts = []
    for piece in version.split("."):
        match = re.match(r"\d+", piece)
        if not match:
            break
        parts.append(int(match.group()))
        if match.group() != piece:
            break
    return tuple(parts)


# Distributions whose import name isn't the pip name with "-" as "_".
IMPORT_NAMES = {"scikit-image": "skimage"}


for spec in sys.argv[1:]:
    name, floor = spec.split(">=")
    dist = None
    for candidate in (name, name.replace("-", "_")):
        try:
            dist = metadata.distribution(candidate)
            break
        except metadata.PackageNotFoundError:
            pass
    if dist is None:
        print("missing", spec, "-", "-")
        continue
    where = str(dist.locate_file("")).rstrip("/") or "-"
    if key(dist.version) < key(floor):
        print("old", spec, dist.version, where)
        continue
    try:
        importlib.import_module(IMPORT_NAMES.get(name, name.replace("-", "_")))
    except Exception as exc:  # any import failure makes it unusable
        print("broken", spec, f"{type(exc).__name__}:{exc}".replace(" ", "_").replace("\n", "_"), where)
        continue
    print("ok", spec, dist.version, where)
