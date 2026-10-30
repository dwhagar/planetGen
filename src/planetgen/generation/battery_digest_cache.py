# planetgen/generation/battery_digest_cache.py

"""
A kept copy of the generator battery digest for a test run (PERF.209)
=====================================================================

The battery digest (`generator_battery`) takes about half a minute and every test process that stores a galaxy asks for it, so a run on four
workers worked it out four times and the next run did it again over the same code. When `ENV` names a folder, `running_digest` keeps the digest
there under a key made from everything it reads: the name, size and modification time of every file of the `planetgen` package and of the
downloaded datasets (`datasets.fetch.cache_dir`), and the Python and numpy versions. A changed file, a new download or a new `git checkout`
(which touches the times) changes the key, so the digest is worked out again and a changed generator is never given an old answer; the
production code never sets `ENV`, so it always works the digest out.
"""

import hashlib
import json
import os
import sys
import tempfile

ENV = "PLANETGEN_BATTERY_DIGEST_CACHE"
"""str: The environment variable that names the folder the digest is kept in (the test suite sets it; nothing else does)."""

_PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


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
    import numpy

    from planetgen.datasets import fetch

    stamp = {"python": sys.version, "numpy": numpy.__version__, "package": _stamp(_PACKAGE),
             "datasets": _stamp(fetch.cache_dir())}
    return hashlib.sha256(json.dumps(stamp, sort_keys=True).encode("utf-8")).hexdigest()


def _path(folder, found):
    return os.path.join(folder, f"battery-{found}.json")


def digest(compute, battery_digest, case_digest):
    """
    The battery digest `compute()` makes, or the one kept for the same key in the folder `ENV` names (`battery_digest` and `case_digest`
    are `generator_battery`'s classes). With `ENV` unset it is `compute()`. A folder that cannot be read or written is ignored.
    """
    folder = os.environ.get(ENV)
    if not folder:
        return compute()
    found = key()
    path = _path(folder, found)
    try:
        with open(path, encoding="utf-8") as handle:
            kept = json.load(handle)
        return battery_digest(kept["overall"], {name: case_digest(**case) for name, case in kept["cases"].items()})
    except (OSError, ValueError, KeyError, TypeError):
        pass
    result = compute()
    try:
        os.makedirs(folder, exist_ok=True)
        text = json.dumps({"overall": result.overall, "cases": {name: vars(case) for name, case in result.cases.items()}})
        descriptor, temporary = tempfile.mkstemp(dir=folder, suffix=".tmp")
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(temporary, path)
    except OSError:
        pass
    return result
