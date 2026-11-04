# planetgen/generation/band_table_cache.py

"""
Kept copies of the bright-star sampling tables for a test run (PERF.218)
========================================================================

A sampling table (`star_population._bright_table`, `_mass_band_table`) is the bound of every cell of the IMF grid at one luminosity threshold,
population, age law and composition, and a first sector fill builds hundreds of them: 717,000 distinct track evaluations, about four minutes on
a slow machine, in every process that fills a sector (PERF.188, TEST.464). The work is distinct, not repeated, so no cache inside the process
helps (it would need over a gigabyte more memory to keep what the tables share). When `ENV` names a folder, a table is kept there as a small
JSON file and read back by the next process that asks for the same table, so the workers of a test run and the next run build each table once.

A file is named by a key of the table's arguments, of every file of the `planetgen` package and of the downloaded datasets (name, size and
modification time: `battery_digest_cache.key`), of the Python and numpy versions, and of every constant of `tuning`, so a changed file, a new
download, a `git checkout` or a test that sets a constant never reads another's table. Numbers are written the way Python writes them (shortest
text that reads back as the same float), so a table read back is the table that was made. The production code never sets `ENV`.
"""

import functools
import hashlib
import json
import os
import shutil
import tempfile
import time

from planetgen import tuning
from planetgen.generation import battery_digest_cache

ENV = "PLANETGEN_BAND_TABLE_CACHE"
"""str: The environment variable that names the folder the tables are kept in (the test suite sets it; nothing else does)."""

STALE_DAYS = 7
"""int: A sibling folder of an older code key not touched for this many days is removed when a new one is made."""


@functools.lru_cache(maxsize=1)
def _code_key():
    """The key of the code and data this process reads, worked out once (its files do not change during a run)."""
    return battery_digest_cache.key()


def _tuning_key():
    """A short key of every constant of `tuning` (a test that sets one gets tables of its own)."""
    constants = sorted((name, repr(value)) for name, value in vars(tuning).items() if name.isupper())
    return hashlib.sha256(json.dumps(constants).encode("utf-8")).hexdigest()


def _folder(root):
    folder = os.path.join(root, _code_key()[:20])
    if not os.path.isdir(folder):
        os.makedirs(folder, exist_ok=True)
        cutoff = time.time() - STALE_DAYS * 86400
        for name in os.listdir(root):
            old = os.path.join(root, name)
            try:
                if old != folder and os.path.isdir(old) and os.stat(old).st_mtime < cutoff:
                    shutil.rmtree(old, ignore_errors=True)
            except OSError:
                continue
    return folder


def kept(kind, arguments, compute, encode, decode):
    """
    The table `compute()` makes for `arguments` (a tuple whose `repr` names the table), or the one kept for the same key in the folder `ENV`
    names. `encode(table)` is the JSON-able form kept and `decode(form)` the table again. With `ENV` unset it is `compute()`; a folder that
    cannot be read or written, or a kept file that does not read, is ignored.
    """
    root = os.environ.get(ENV)
    if not root:
        return compute()
    try:
        name = hashlib.sha256(f"{kind}|{repr(arguments)}|{_tuning_key()}".encode("utf-8")).hexdigest()
        path = os.path.join(_folder(root), f"{kind}-{name}.json")
    except OSError:
        return compute()
    try:
        with open(path, encoding="utf-8") as handle:
            return decode(json.load(handle))
    except (OSError, ValueError, KeyError, TypeError):
        pass
    table = compute()
    try:
        descriptor, temporary = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(encode(table), handle)
        os.replace(temporary, path)
    except OSError:
        pass
    return table
