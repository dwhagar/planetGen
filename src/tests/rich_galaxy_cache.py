# tests/rich_galaxy_cache.py

"""
One build of the rich test galaxy for every process and machine that tests the same code (TEST.788)
===================================================================================================

The `rich_galaxy` fixture builds its galaxy by running the real commands, which takes about half an hour, and each test process that uses it builds
it again: on a merge train three files (`test_db_column_round_trip`, `test_db_check_constraints`, `test_cellkey_motion`) do, on three machines, from the
same commit. With `PLANETGEN_RICH_GALAXY_CACHE` set to a folder every machine can reach, the first process to need the galaxy of a given code state builds it
and leaves a dump of its database there; the others load the dump (seconds) and read the same rows. The galaxy is the same galaxy, with the same systems:
nothing is made smaller, and a process without the variable builds as before.

The dump is named by a key over everything that decides what the build makes: every file of the package, the support modules the build runs through, the
list of the datasets it reads, the versions of the database server, Python, numpy and numba. A different commit, or a file changed for a test, is another key and builds its own;
a stale dump can only be found by a key that is the same.
"""

import contextlib
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

ENVIRONMENT_VARIABLE = "PLANETGEN_RICH_GALAXY_CACHE"
"""str: The folder that holds the dumps. Set to a folder every machine can reach, the dumps are shared by them; set but empty, no cache (TEST.807); unset, the folder of this machine (`default_folder`)."""

STALE_DUMP_DAYS = 3
"""float: A dump nobody has loaded for this long belongs to a code state that is gone (the key follows every file of the package, so each commit has its own) and is removed when a new one is kept (TEST.807)."""

LOCK_STALE_SECONDS = 3 * 3600
"""float: A lock older than this belongs to a build that died; the next process takes the build over."""

WAIT_SECONDS = 3 * 3600
"""float: The longest a process waits for another's build before it builds its own, without caching."""

POLL_SECONDS = 15.0

SUPPORT_FILES = ("rich_galaxy_support.py", "bughunt_support.py", "fuzz_support.py", "rich_galaxy_cache.py")
"""tuple: The modules of `tests/` the build runs through, in the key."""


def default_folder():
    """The folder a machine keeps the dumps in when `ENVIRONMENT_VARIABLE` is not set: `planetgen-rich-galaxy` in the machine's temporary folder."""
    return os.path.join(tempfile.gettempdir(), "planetgen-rich-galaxy")


def cache_folder(environ=None):
    """
    The folder the dumps are kept in, or `None` when the cache is switched off. TEST.807: unset, a machine keeps its own (`default_folder`), so the test files of a gate that
    each take the galaxy in a process of their own (a shard to a file, a train gate on one machine) build it once between them instead of once each, which was 38 of the 62 minutes of
    `test_cellkey_motion` on Train 51; a folder every machine reaches shares it between machines; set but empty switches the cache off.
    """
    value = (os.environ if environ is None else environ).get(ENVIRONMENT_VARIABLE)
    if value is None:
        return default_folder()
    return value.strip() or None


def prune(folder, now=time.time):
    """Removes the dumps of `folder` nobody has loaded for `STALE_DUMP_DAYS`, and the lock folders left by a build that died; what it cannot remove it leaves."""
    for name in os.listdir(folder):
        path = os.path.join(folder, name)
        limit = LOCK_STALE_SECONDS if name.endswith(".lock") else STALE_DUMP_DAYS * 86400
        with contextlib.suppress(OSError):
            if name.startswith("rich-") and now() - os.path.getmtime(path) > limit:
                if os.path.isdir(path):
                    os.rmdir(path)
                else:
                    os.remove(path)


MACHINE_VARIABLES = ("MYSQL", "REDIS", "CACHE", "DATA_DIR", "TEST_")
"""tuple: Parts of the names of `PLANETGEN_*` variables that say where a machine's servers and folders are and not what the generator makes; the others (settings the generator reads) are in the key."""


def key_of(config, package_root, tests_root, data_root=None, extra_files=(), environ=None):
    """
    The key of the galaxy a build of this code state makes: a hash of every file under `package_root` (name and bytes, `__pycache__` left out), the files of `SUPPORT_FILES` under
    `tests_root`, the `extra_files` the build reads (a system file, the site's `config.json`) where they exist, the names and sizes of the files under `data_root` (the datasets the generator reads),
    the `PLANETGEN_*` variables that are settings, and the versions of the server of `config`, Python, numpy and numba.
    """
    import numba
    import numpy
    import pymysql

    digest = hashlib.sha256()

    def add_tree(root, names_only=False):
        for folder, folders, files in os.walk(root):
            folders[:] = sorted(name for name in folders if name != "__pycache__")
            for name in sorted(files):
                if name.endswith((".pyc", ".pyo")):
                    continue
                path = os.path.join(folder, name)
                digest.update(os.path.relpath(path, root).encode() + b"\0")
                if names_only:
                    digest.update(str(os.path.getsize(path)).encode())
                else:
                    with open(path, "rb") as handle:
                        digest.update(hashlib.sha256(handle.read()).digest())

    add_tree(package_root)
    for name in SUPPORT_FILES:
        with open(os.path.join(tests_root, name), "rb") as handle:
            digest.update(name.encode() + hashlib.sha256(handle.read()).digest())
    for path in extra_files:
        if os.path.isfile(path):
            with open(path, "rb") as handle:
                digest.update(os.path.basename(path).encode() + hashlib.sha256(handle.read()).digest())
    if data_root and os.path.isdir(data_root):
        add_tree(data_root, names_only=True)
    for name, value in sorted((os.environ if environ is None else environ).items()):
        if name.startswith("PLANETGEN_") and not any(part in name for part in MACHINE_VARIABLES):
            digest.update(f"{name}={value}".encode())
    connection = pymysql.connect(host=config.host, port=config.port, user=config.user, password=config.password)
    try:
        with connection.cursor() as cursor:
            cursor.execute("SELECT VERSION()")
            server = cursor.fetchone()[0]
    finally:
        connection.close()
    digest.update(f"{server}|{sys.version}|numpy {numpy.__version__}|numba {numba.__version__}".encode())
    return digest.hexdigest()[:32]


def _client_arguments(config):
    return [f"--host={config.host}", f"--port={config.port}", f"--user={config.user}", "--protocol=TCP"]


def _client_environment(config):
    return dict(os.environ, MYSQL_PWD=config.password or "")


_NEGATIVE_ZERO = re.compile(rb"'(?:[^'\\]|\\.)*'|(?<=[(,])-0(?=[,)])", re.S)
"""re.Pattern: A quoted string (left alone) or a bare `-0` between the separators of a row of an INSERT."""


def keep_negative_zero(line):
    """
    One line of a dump with every `-0` value of an INSERT written `-0e0`. `mysqldump` writes a DOUBLE or FLOAT of minus zero as `-0`, which loading reads as the integer 0 negated, plus zero; as `-0e0` it is
    a number with an exponent, a double, and loading keeps the sign, so the loaded rows are the built ones to the last bit.
    """
    if not line.startswith(b"INSERT INTO"):
        return line
    return _NEGATIVE_ZERO.sub(lambda found: found.group(0) if found.group(0).startswith(b"'") else b"-0e0", line)


def dump(config, path):
    """Writes the database of `config` as gzipped SQL (tables, rows, triggers, routines) to `path`, through a temporary file renamed at the end so that a reader never sees half a dump."""
    temporary = f"{path}.part{os.getpid()}"
    try:
        with gzip.open(temporary, "wb", compresslevel=3) as target:
            process = subprocess.Popen(["mysqldump", *_client_arguments(config), "--single-transaction", "--routines", "--triggers", "--no-tablespaces",
                                        "--skip-comments", "--hex-blob", "--default-character-set=utf8mb4", "--max-allowed-packet=512M", config.database],
                                       stdout=subprocess.PIPE, env=_client_environment(config))
            for line in process.stdout:
                target.write(keep_negative_zero(line))
            if process.wait() != 0:
                raise RuntimeError(f"mysqldump of {config.database} failed ({process.returncode})")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def load(config, path):
    """Loads the dump at `path` into the (empty) database of `config`."""
    process = subprocess.Popen(["mysql", *_client_arguments(config), "--default-character-set=utf8mb4", "--max-allowed-packet=512M", config.database],
                               stdin=subprocess.PIPE, env=_client_environment(config))
    with gzip.open(path, "rb") as source:
        try:
            shutil.copyfileobj(source, process.stdin, 1 << 20)
        finally:
            process.stdin.close()
    if process.wait() != 0:
        raise RuntimeError(f"loading {path} into {config.database} failed ({process.returncode})")


@contextlib.contextmanager
def _lock(folder, key, sleep=time.sleep, now=time.time):
    """Takes the build of `key`: `True` when this process holds it (a folder made by `mkdir`, which one process can make), `False` when another does."""
    path = os.path.join(folder, f"rich-{key}.lock")
    try:
        os.mkdir(path)
    except FileExistsError:
        try:
            if now() - os.path.getmtime(path) > LOCK_STALE_SECONDS:
                os.rmdir(path)               # a build that died: take it over
                os.mkdir(path)
            else:
                yield False
                return
        except (FileNotFoundError, FileExistsError):
            yield False
            return
    try:
        yield True
    finally:
        with contextlib.suppress(OSError):
            os.rmdir(path)


def build_or_load(config, build, folder, key, sleep=time.sleep, now=time.time, wait_seconds=None):
    """
    Fills the empty database of `config` with the galaxy of `key`: loads the dump in `folder` when there is one; otherwise the process that takes the lock calls `build()` and leaves a dump
    for the others, and the others wait for it (up to `wait_seconds`, then they build without caching). Returns `"loaded"`, `"built and kept"` or `"built"`.
    """
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, f"rich-{key}.sql.gz")
    deadline = now() + (WAIT_SECONDS if wait_seconds is None else wait_seconds)
    while True:
        if os.path.exists(path):
            load(config, path)
            with contextlib.suppress(OSError):
                os.utime(path)          # a dump in use is not stale
            return "loaded"
        with _lock(folder, key, sleep, now) as mine:
            if mine:
                if os.path.exists(path):          # another process finished between the look and the lock
                    load(config, path)
                    return "loaded"
                build()
                dump(config, path)
                prune(folder, now)
                return "built and kept"
        if now() > deadline:
            build()
            return "built"
        sleep(POLL_SECONDS)
