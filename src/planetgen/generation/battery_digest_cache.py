# planetgen/generation/battery_digest_cache.py

"""
A kept copy of the generator battery digest (PERF.209, PERF.250)
================================================================

The battery digest (`generator_battery`) takes from half a minute to over a minute, and every process that stores or extends a galaxy asks for it
(a fill, a plan, `fingerprint`, the update's record), so a run on four workers worked it out four times and the next run did it again over the same
code. `running_digest` keeps the digest in a folder under a key made from everything it reads: the name, size and modification time of every file of
the `planetgen` package and of the downloaded datasets (`datasets.fetch.cache_dir`), the Python, numpy and scipy versions and the platform. A changed
file, a new download, a new library or a new `git checkout` (which touches the times) changes the key, so the digest is worked out again and a changed
generator is never given an old answer.

The folder is what `ENV` names (the test suite sets it; set but empty, none is kept); unset, it is the folder an install keeps its kept star tables in
(`band_table_cache.root()`, made by `update.sh`), so the first run after an update pays for the digest once and every later process reads it. A host
without the folder works the digest out as before. The file is readable by every user (the update may be run by one and the web jobs by another), and
copies of older code are removed after `STALE_DAYS`.
"""

import hashlib
import json
import os
import platform
import sys
import tempfile
import time

ENV = "PLANETGEN_BATTERY_DIGEST_CACHE"
"""str: The environment variable that names the folder the digest is kept in (the test suite sets it; nothing else does)."""

STALE_DAYS = 7
"""int: A kept digest of older code not touched for this many days is removed when a new one is kept."""

_PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def folder():
    """The folder the digest is kept in, or `None` to keep none: what `ENV` names (set but empty: none), else the folder of an install's kept star
    tables (`band_table_cache.root()`) where it exists."""
    if ENV in os.environ:
        return os.environ[ENV] or None
    from planetgen.generation import band_table_cache      # it imports this module
    return band_table_cache.root()


def _library_versions():
    """The versions of the libraries whose arithmetic the battery's floats come from, and the platform they run on."""
    found = {"python": sys.version, "implementation": platform.python_implementation(), "system": platform.system(),
             "machine": platform.machine(), "libc": list(platform.libc_ver())}
    for name in ("numpy", "scipy"):
        try:
            found[name] = __import__(name).__version__
        except ImportError:
            found[name] = None
    return found


def _stamp(root):
    """`(relative path, size, modification time in ns)` of every file below `root`, sorted, leaving out compiled files."""
    found = []
    for folder, names, files in os.walk(root):
        names[:] = sorted(name for name in names if name != "__pycache__")
        for name in sorted(files):
            path = os.path.join(folder, name)
            try:
                info = os.stat(path)
            except OSError:
                continue
            found.append((os.path.relpath(path, root), info.st_size, info.st_mtime_ns))
    return found


def key():
    """The hex key of the code and data the battery reads now."""
    from planetgen.datasets import fetch

    stamp = {"libraries": _library_versions(), "package": _stamp(_PACKAGE), "datasets": _stamp(fetch.cache_dir())}
    return hashlib.sha256(json.dumps(stamp, sort_keys=True).encode("utf-8")).hexdigest()


def _path(folder, found):
    return os.path.join(folder, f"battery-{found}.json")


def digest(compute, battery_digest, case_digest):
    """
    The battery digest `compute()` makes, or the one kept for the same key in `folder()` (`battery_digest` and `case_digest` are
    `generator_battery`'s classes). With no folder it is `compute()`. A folder that cannot be read or written is ignored.
    """
    where = folder()
    if not where:
        return compute()
    found = key()
    path = _path(where, found)
    try:
        with open(path, encoding="utf-8") as handle:
            kept = json.load(handle)
        return battery_digest(kept["overall"], {name: case_digest(**case) for name, case in kept["cases"].items()})
    except (OSError, ValueError, KeyError, TypeError):
        pass
    result = compute()
    try:
        os.makedirs(where, exist_ok=True)
        text = json.dumps({"overall": result.overall, "cases": {name: vars(case) for name, case in result.cases.items()}})
        descriptor, temporary = tempfile.mkstemp(dir=where, suffix=".tmp")
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.chmod(temporary, 0o644)          # whoever runs the update, the web user's jobs read it
        os.replace(temporary, path)
        _remove_stale(where, path)
    except OSError:
        pass
    return result


def _remove_stale(where, keep):
    """Removes the other kept digests of `where` that have not been touched for `STALE_DAYS`."""
    cutoff = time.time() - STALE_DAYS * 86400
    for name in os.listdir(where):
        old = os.path.join(where, name)
        try:
            if name.startswith("battery-") and name.endswith(".json") and old != keep and os.stat(old).st_mtime < cutoff:
                os.remove(old)
        except OSError:
            continue
