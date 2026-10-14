#!/usr/bin/env python3
"""
Checks planetGen's Python requirements as this interpreter sees them, for
scripts/install-python-deps.sh (Linux and macOS). Standard library only.

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


def _pymysql_logs_in_to_mysql_8():
    """pymysql can log in with caching_sha2_password and sha256_password (MySQL 8's default accounts) only with the
    `cryptography` package; without it the first connection fails ("'cryptography' package is required")."""
    from pymysql import _auth
    if not _auth._have_cryptography:
        raise ImportError("pymysql_has_no_cryptography_for_caching_sha2_password")


# What a distribution has to be able to do beyond importing.
EXTRA_CHECKS = {"pymysql": _pymysql_logs_in_to_mysql_8}


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
        if name in EXTRA_CHECKS:
            EXTRA_CHECKS[name]()
    except Exception as exc:  # any import failure makes it unusable
        print("broken", spec, f"{type(exc).__name__}:{exc}".replace(" ", "_").replace("\n", "_"), where)
        continue
    print("ok", spec, dist.version, where)
