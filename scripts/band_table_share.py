# scripts/band_table_share.py

"""
Shares the kept sampling tables between lane machines (PERF.232)
================================================================

A first sector fill builds 357 bright-star sampling tables, about 100 s of the cold run of the slow tests (`generation/band_table_cache.py`). The
kept copy is keyed by the name, size and modification time of every file, so a fresh checkout on another machine never reads another machine's
tables, though the tables are the same. This script keeps a copy in a shared folder under a key of the file *contents* (the code of the
`planetgen` package, the Python and numpy versions and the names and sizes of the downloaded datasets), so a copy made by other code is never
read. It is a convenience for the lane machines only; production reads and writes its own folder and never runs this.

    python scripts/band_table_share.py restore    # before a test run: copy the shared tables for this code in
    python scripts/band_table_share.py save       # after a run that built them: copy the folder in for this code

The shared folder is `PLANETGEN_BAND_TABLE_SHARE`, else `/mnt/project-files/cache/band-tables`; the local one is what `band_table_cache.root()`
gives (set `PLANETGEN_BAND_TABLE_CACHE` as the tests do). A shared copy for other code is left alone for `STALE_DAYS` and then removed by `save`.
"""

import argparse
import hashlib
import json
import os
import shutil
import sys
import time

import numpy

from planetgen.datasets import fetch
from planetgen.generation import band_table_cache

SHARE_ENV = "PLANETGEN_BAND_TABLE_SHARE"
"""str: The environment variable that names the shared folder."""

DEFAULT_SHARE = "/mnt/project-files/cache/band-tables"
"""str: The shared folder of the lane machines."""

STALE_DAYS = 7
"""int: A shared copy for other code not touched for this many days is removed by `save`."""

def share_root():
    """The shared folder."""
    return os.environ.get(SHARE_ENV) or DEFAULT_SHARE


def content_key(package=None, datasets=None):
    """The hex key of the contents of the `planetgen` package (sorted relative paths and file bytes), the Python and numpy versions, and the names and sizes of the downloaded datasets."""
    package = package or os.path.dirname(os.path.dirname(os.path.abspath(band_table_cache.__file__)))
    datasets = datasets if datasets is not None else fetch.cache_dir()
    digest = hashlib.sha256()
    digest.update(f"{sys.version}|{numpy.__version__}".encode("utf-8"))
    for folder, names, files in os.walk(package):
        names[:] = sorted(name for name in names if name != "__pycache__")
        for name in sorted(files):
            if name.endswith((".pyc", ".pyo")):
                continue
            path = os.path.join(folder, name)
            digest.update(os.path.relpath(path, package).encode("utf-8"))
            with open(path, "rb") as handle:
                digest.update(hashlib.sha256(handle.read()).digest())
    listing = []
    for folder, names, files in os.walk(datasets):
        names[:] = sorted(names)
        for name in sorted(files):
            path = os.path.join(folder, name)
            try:
                listing.append((os.path.relpath(path, datasets), os.stat(path).st_size))
            except OSError:
                continue
    digest.update(json.dumps(listing).encode("utf-8"))
    return digest.hexdigest()[:20]


def _local_folder():
    base = band_table_cache.root()
    if not base:
        raise SystemExit("No local table folder: set PLANETGEN_BAND_TABLE_CACHE to one.")
    os.makedirs(base, exist_ok=True)
    return band_table_cache._folder(base)


def _copy_new(source, target):
    os.makedirs(target, exist_ok=True)
    count = 0
    for name in os.listdir(source):
        if name.endswith(".json") and not os.path.exists(os.path.join(target, name)):
            shutil.copyfile(os.path.join(source, name), os.path.join(target, name) + ".tmp")
            os.replace(os.path.join(target, name) + ".tmp", os.path.join(target, name))
            count += 1
    return count


def restore(share=None):
    """Copies the shared tables for this code into the local folder, and returns how many were new (0 when there is no shared copy for this code)."""
    source = os.path.join(share or share_root(), content_key())
    if not os.path.isdir(source):
        return 0
    return _copy_new(source, _local_folder())


def save(share=None):
    """Copies the local tables into the shared folder for this code, removes shared copies of other code not touched for `STALE_DAYS`, and returns how many were new."""
    share = share or share_root()
    key = content_key()
    count = _copy_new(_local_folder(), os.path.join(share, key))
    os.utime(os.path.join(share, key))
    cutoff = time.time() - STALE_DAYS * 86400
    for name in os.listdir(share):
        old = os.path.join(share, name)
        if name != key and os.path.isdir(old) and os.stat(old).st_mtime < cutoff:
            shutil.rmtree(old, ignore_errors=True)
    return count


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("action", choices=("restore", "save"))
    arguments = parser.parse_args(argv)
    done = restore() if arguments.action == "restore" else save()
    print(f"band tables {arguments.action}d: {done} new")


if __name__ == "__main__":
    main()
