# planetgen/generation/band_table_cache.py

"""
Kept copies of the bright-star sampling tables (PERF.218, PERF.188)
===================================================================

A sampling table (`star_population._bright_table`, `_mass_band_table`) is the bound of every cell of the IMF grid at one luminosity threshold,
population, age law and composition, and a first sector fill builds hundreds of them: 717,000 distinct track evaluations, about four minutes on
a slow machine, in every process that fills a sector (PERF.188, TEST.464). The work is distinct, not repeated, so no cache inside the process
helps (it would need over a gigabyte more memory to keep what the tables share). When `ENV` names a folder, a table is kept there as a small
JSON file and read back by the next process that asks for the same table, so the workers of a run and the next run build each table once.

A file is named by a key of the table's arguments, of the contents of the files a table is made from (`INPUT_FILES`: the 28 files a build runs and the 38 they import directly, of 1,166; PERF.281), of the downloaded
datasets (names and sizes), of the Python and numpy versions, and of the `tuning` constants those files read (`_tuning_names`), so a changed input, a new download or a test that sets a constant never reads another's
table, and an update that changes none of them keeps the tables of the run before it (until PERF.281 the key was the name, size and time of every file of the package and every constant, so each update started the next
run's tables cold: 16 minutes a sky share, 35 s a tier, on a stalled run). A table is built under a lock of its own, so the workers of a run that start on the same tables build each once. Numbers are written the way
Python writes them (shortest text that reads back as the same float), so a table read back is the table that was made.

In production (PERF.188, Boss's call 11:24 AM EST on 2026-11-04) the folder is `DEFAULT_ROOT`, which `update.sh` creates at setup
(`examples/apache/create-cache-dir.sh`, owned by the web user whose workers write it): the first sector fill of a galaxy went from 267 s and
1,842 MB to 36 s and 534 MB with the same output. The tables of an older version of the code can never be read again, so a new code key removes every
other folder of the root; a folder named by `ENV` (a test run) keeps its siblings for `STALE_DAYS`. A host without the folder builds the tables as before.
"""

import ast
import contextlib
import functools
import hashlib
import json
import os
import shutil
import sys
import tempfile
import time

from planetgen import tuning
from planetgen.generation import battery_digest_cache, battery_progress, bodybeat
from planetgen.queue import unitbeat

ENV = "PLANETGEN_BAND_TABLE_CACHE"
"""str: The environment variable that names the folder the tables are kept in (the test suite sets it); set but empty, the tables are not kept. Unset, the folder is `DEFAULT_ROOT` where it exists."""

DEFAULT_ROOT = "/var/cache/planetgen/band-tables"
"""str: The folder the tables are kept in on an install (made by `update.sh`, beside the tile cache and the datasets, which hold the code key's inputs); without it nothing is kept."""

STALE_DAYS = 7
"""int: In a folder named by `ENV`, a sibling folder of an older code key not touched for this many days is removed when a new one is made."""


def root():
    """The folder the tables are kept in, or `None` to keep none: what `ENV` names, else `DEFAULT_ROOT` when that folder exists."""
    if ENV in os.environ:
        return os.environ[ENV] or None
    return DEFAULT_ROOT if os.path.isdir(DEFAULT_ROOT) else None


INPUT_FILES = (
    "datasets/axes.py",
    "datasets/baraffe_isochrones.py",
    "datasets/evolved_tracks.py",
    "datasets/fetch.py",
    "datasets/isochrones.py",
    "datasets/mamajek_sequence.py",
    "datasets/white_dwarf_models.py",
    "db/statementmark.py",
    "galaxy/age_law.py",
    "galaxy/scales.py",
    "generation/band_table_cache.py",
    "generation/battery_digest_cache.py",
    "generation/battery_progress.py",
    "generation/bodybeat.py",
    "generation/close_binary_trend.py",
    "generation/mass_ratio.py",
    "generation/member_record.py",
    "generation/multistar.py",
    "generation/phenomena/hosts.py",
    "generation/phenomena/system_types.py",
    "generation/possibility_table.py",
    "generation/primary_mass_function.py",
    "generation/star_population.py",
    "physics/accretion.py",
    "physics/age_metallicity.py",
    "physics/alpha_tracks.py",
    "physics/bright_piece_kernel.py",
    "physics/brown_dwarf.py",
    "physics/constants.py",
    "physics/degenerate_radius.py",
    "physics/elements.py",
    "physics/giant_branch.py",
    "physics/helium_core_white_dwarf.py",
    "physics/hierarchy_stability.py",
    "physics/imf.py",
    "physics/main_sequence.py",
    "physics/mass_boundaries.py",
    "physics/metallicity.py",
    "physics/ms_lifetime.py",
    "physics/neutron_star.py",
    "physics/pair_history.py",
    "physics/post_main_sequence.py",
    "physics/pre_main_sequence.py",
    "physics/remnant_mass.py",
    "physics/remnant_progenitors.py",
    "physics/roche.py",
    "physics/stability.py",
    "physics/star_orbit_alignment.py",
    "physics/stellar_abundances.py",
    "physics/stellar_composition.py",
    "physics/stellar_evolution.py",
    "physics/system_tree.py",
    "physics/unit_registry.py",
    "physics/wolf_rayet.py",
    "queue/beatjudge.py",
    "queue/runrecord.py",
    "queue/unitbeat.py",
    "queue/usage.py",
    "util/detmath.py",
    "util/draw.py",
    "util/format.py",
    "util/log.py",
    "util/request_odds.py",
    "util/ticker.py",
    "util/user_agent.py",
)
"""tuple: The files of the `planetgen` package a table's numbers can come from (PERF.281), relative to the package: every file that a build of a table ran over every population, threshold, mass range and
composition (`tests/test_band_table_cache.py` traces the builds and fails if one runs a file that is not here), and the files those import directly, which hold the constants and helpers they use without
calling anything. A change to any other file (the web, the database, a planet's physics) cannot change a table, so it keeps the tables an update found; `tuning.py` is read by the names these files use (`_tuning_names`).
(a list of what the tables read, not a physical quantity)"""

_PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@functools.lru_cache(maxsize=1)
def _tuning_names():
    """The `tuning` constants the files of `INPUT_FILES` read: every name they take of the module (`tuning.NAME`, whatever they call it, `from planetgen.tuning import NAME`), sorted. A file that reads the module
    some other way (`getattr`, `vars`) gives every constant."""
    names, whole = set(), False
    for relative in INPUT_FILES:
        if relative == "generation/band_table_cache.py":
            continue                                  # this file reads the module to make the key, not to make a table
        with open(os.path.join(_PACKAGE, relative), encoding="utf-8") as handle:
            tree = ast.parse(handle.read())
        aliases = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "planetgen" and any(alias.name == "tuning" for alias in node.names):
                aliases.update(alias.asname or alias.name for alias in node.names if alias.name == "tuning")
            elif isinstance(node, ast.ImportFrom) and node.module == "planetgen.tuning":
                names.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                aliases.update(alias.asname for alias in node.names if alias.name == "planetgen.tuning" and alias.asname)
            elif isinstance(node, ast.Assign) and isinstance(node.value, ast.Name) and node.value.id in aliases:
                aliases.update(target.id for target in node.targets if isinstance(target, ast.Name))
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id in aliases:
                names.add(node.attr)
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ("getattr", "vars") and node.args \
                    and isinstance(node.args[0], ast.Name) and node.args[0].id in aliases:
                whole = True
    if whole:
        return tuple(sorted(name for name in vars(tuning) if name.isupper()))
    return tuple(sorted(name for name in names if name.isupper() and hasattr(tuning, name)))


@functools.lru_cache(maxsize=1)
def _code_key():
    """The key of the code and data a table is made from, worked out once (its files do not change during a run): the contents of `INPUT_FILES` (not their times, so a `git checkout` of the same code keeps the tables), the
    sizes of the downloaded datasets, and the Python and numpy versions. An update that changes nothing these read leaves the key, and the tables kept under it, as they were (PERF.281; the key was the name, size and modification
    time of every file of the package, so every update started the next run's tables cold)."""
    from planetgen.datasets import fetch

    digest = hashlib.sha256()
    digest.update(json.dumps(battery_digest_cache._library_versions(), sort_keys=True).encode("utf-8"))
    for relative in INPUT_FILES:
        with open(os.path.join(_PACKAGE, relative), "rb") as handle:
            digest.update(relative.encode("utf-8") + hashlib.sha256(handle.read()).digest())
    digest.update(json.dumps([(name, size) for name, size, _mtime in battery_digest_cache._stamp(fetch.cache_dir())]).encode("utf-8"))
    return digest.hexdigest()


def _tuning_key():
    """A short key of the `tuning` constants the table files read (a test that sets one gets tables of its own)."""
    constants = [(name, repr(tuning.__dict__[name])) for name in _tuning_names()]
    return hashlib.sha256(json.dumps(constants).encode("utf-8")).hexdigest()


def _folder(base):
    folder = os.path.join(base, _code_key()[:20])
    if not os.path.isdir(folder):
        os.makedirs(folder, exist_ok=True)
        cutoff = time.time() - (STALE_DAYS * 86400 if ENV in os.environ else 0)
        for name in os.listdir(base):
            old = os.path.join(base, name)
            try:
                if old != folder and os.path.isdir(old) and os.stat(old).st_mtime < cutoff:
                    shutil.rmtree(old, ignore_errors=True)
            except OSError:
                continue
    return folder


def _read(path, decode):
    """The table kept at `path`, or `None` when none is kept or it does not read."""
    try:
        with open(path, encoding="utf-8") as handle:
            return decode(json.load(handle))
    except (OSError, ValueError, KeyError, TypeError):
        return None


_built = [0]
"""list: How many tables this process has built (not read back) so far, for the line the unit says while it builds one."""

LOCK_POLL_SECONDS = 0.2
"""float: How often a process that waits for another's table looks at its lock again; machinery, not physics, so it is defined here and has no registry row (TEST.715)."""


@contextlib.contextmanager
def _building(path):
    """
    Holds the table's own lock while one process builds it (PERF.281): the workers of a run that start together on the same tables (a run's first, or the first after an update that changed what a table reads) each
    built every table before, three at once for work that one did, and now one builds and the others read what it kept. The lock is on a file beside the table and the system lets it go when the process ends, however
    it ends. Where locks cannot be had (no `fcntl`, a folder that cannot be written) the table is built without one, as before.
    """
    handle = None
    try:
        import fcntl
        handle = open(path + ".lock", "a", encoding="utf-8")
        with unitbeat.waiting("a sampling table another process is building"):
            while True:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    bodybeat.settled()              # a unit that waits for a table is alive, as one that builds it is (PERF.273)
                    time.sleep(LOCK_POLL_SECONDS)
    except (ImportError, OSError):
        if handle is not None:
            handle.close()
            handle = None
    try:
        yield
    finally:
        if handle is not None:
            handle.close()


def kept(kind, arguments, compute, encode, decode):
    """
    The table `compute()` makes for `arguments` (a tuple whose `repr` names the table), or the one kept for the same key in the folder `ENV`
    names. `encode(table)` is the JSON-able form kept and `decode(form)` the table again. With `ENV` unset it is `compute()`; a folder that
    cannot be read or written, or a kept file that does not read, is ignored.
    """
    base = root()
    if not base:
        return compute()
    try:
        name = hashlib.sha256(f"{kind}|{repr(arguments)}|{_tuning_key()}".encode("utf-8")).hexdigest()
        path = os.path.join(_folder(base), f"{kind}-{name}.json")
    except OSError:
        return compute()
    found = _read(path, decode)
    if found is not None:
        return found
    with _building(path):
        found = _read(path, decode)               # a process that held the table's lock built it while this one waited
        if found is not None:
            return found
        _built[0] += 1
        unitbeat.phase(f"Building the {kind} sampling table {_built[0]:,} of this process (not kept yet by any other)")      # a table built cold is seconds of work with nothing else moving
        table = compute()
        battery_progress.table_built()      # OPS.184: the digest process counts the tables it works out
        try:
            descriptor, temporary = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
            with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                json.dump(encode(table), handle)
            os.replace(temporary, path)
        except OSError:
            pass
    return table
