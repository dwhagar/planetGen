"""planetGen: procedural galaxy, star system and planet generation.

The package layout is docs/design/architecture/library-migration.md, section 6.
"""

import os as _os


NUMBA_CACHE_DEFAULT = "/var/cache/planetgen/numba"
"""str: The folder numba keeps its compiled kernels in on an install (made by `create-cache-dir.sh`, owned by Apache's user). Since SEC.37 the code tree is root's and read-only to the service user, so numba's default (a `__pycache__` beside each source file) cannot be written and `@njit(cache=True)` raises "no locator available"."""


def _writable(folder):
    return _os.path.isdir(folder) and _os.access(folder, _os.W_OK | _os.X_OK)


def choose_numba_cache_dir(environ=None, here=None):
    """The folder numba should cache in, or `None` to leave numba's own default.

    `NUMBA_CACHE_DIR` already set wins (`None`). Else `NUMBA_CACHE_DEFAULT` when this process can write it. Else numba's default when the
    package folder is writable (a development checkout). Else a private folder per user under the temporary directory, so a read-only
    code tree never stops a kernel from loading (it only compiles again on each start).
    """
    environ = _os.environ if environ is None else environ
    here = here or _os.path.dirname(_os.path.abspath(__file__))
    if environ.get("NUMBA_CACHE_DIR"):
        return None
    if _writable(NUMBA_CACHE_DEFAULT):
        return NUMBA_CACHE_DEFAULT
    if _writable(here) and (not _os.path.isdir(_os.path.join(here, "__pycache__")) or _writable(_os.path.join(here, "__pycache__"))):
        return None
    import tempfile
    folder = _os.path.join(tempfile.gettempdir(), "planetgen-numba-%d" % _os.getuid())
    _os.makedirs(folder, mode=0o700, exist_ok=True)
    return folder


def _set_numba_cache_dir():
    """Points `NUMBA_CACHE_DIR` somewhere writable before any module imports numba, so every entry point (CLI, job runner, workers, maintenance scripts) gets it."""
    try:
        folder = choose_numba_cache_dir()
    except OSError:
        return
    if folder:
        _os.environ["NUMBA_CACHE_DIR"] = folder


_set_numba_cache_dir()
