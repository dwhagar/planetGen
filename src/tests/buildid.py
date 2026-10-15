# tests/buildid.py

"""
The id of a test run's build (TEST.150): two runs with the same id ran the same code, tests, packages, servers,
command and selection on the same day, so what passed in one says what would pass in the other (`scripts/test.py
--resume`, TEST.151).

The id is a hash of named parts, kept next to the id so a resume that finds nothing can say which part changed:

- `tree`: the git tree of HEAD (every tracked file, lock files and `pytest.ini` among them);
- `dirty`: the content of every uncommitted or untracked file that git does not ignore;
- `packages`: every installed distribution and its version;
- `python`: the interpreter's version;
- `mysql`, `redis`: the servers' versions as the tests reach them (`unknown` when they do not answer);
- `command`: the `scripts/test.py` command, its areas and the pytest arguments, and `selection`: the commit
  `--changed` compares with;
- `env`: every `PLANETGEN_*` and `PYTEST_*` variable (except the log paths and the worker's own), hashed;
- `date`: the UTC day, so a run never resumes across midnight.

Parts are hashes or short strings; none shows a value from the environment.
"""

import datetime
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys

BIG_FILE_BYTES = 50_000_000
"""int: A file bigger than this is hashed by its size and modification time, not its content."""
IGNORED_ENV = ("PLANETGEN_FAILURE_LOG", "PLANETGEN_PROGRESS_LOG", "PYTEST_XDIST_WORKER", "PYTEST_XDIST_WORKER_COUNT",
               "PYTEST_XDIST_TESTRUNUID", "PYTEST_CURRENT_TEST", "PYTEST_VERSION")


def _git(root, *args):
    result = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, errors="replace")
    return result.stdout if result.returncode == 0 else ""


def _short(text):
    return hashlib.sha256(text.encode("utf-8", "replace")).hexdigest()[:12]


def _tree(root):
    return _git(root, "rev-parse", "HEAD^{tree}").strip() or "none"


def _file_hash(path):
    try:
        size = os.path.getsize(path)
        if size > BIG_FILE_BYTES:
            return f"big:{size}:{os.stat(path).st_mtime_ns}"
        digest = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return "gone"


def _dirty(root):
    """A hash of the uncommitted changes: git's status of each path and the content now on disk."""
    status = _git(root, "status", "--porcelain=v1", "-z", "--untracked-files=all")
    entries = [entry for entry in status.split("\0") if entry]
    lines = []
    for entry in sorted(entries):
        code, path = entry[:2], entry[3:]
        lines.append(f"{code} {path} {_file_hash(os.path.join(root, path))}")
    return _short("\n".join(lines)) if lines else "clean"


def _packages():
    names = sorted(f"{(d.metadata['Name'] or '').lower()}=={d.version}" for d in importlib.metadata.distributions())
    return _short("\n".join(names))


def _mysql_version(environ):
    try:
        import pymysql
        conn = pymysql.connect(
            host=environ.get("PLANETGEN_TEST_MYSQL_HOST", environ.get("PLANETGEN_MYSQL_HOST", "127.0.0.1")),
            port=int(environ.get("PLANETGEN_TEST_MYSQL_PORT", environ.get("PLANETGEN_MYSQL_PORT", 3306))),
            user=environ.get("PLANETGEN_TEST_MYSQL_USER", environ.get("PLANETGEN_MYSQL_USER", "planetgen")),
            password=environ.get("PLANETGEN_TEST_MYSQL_PASSWORD", environ.get("PLANETGEN_MYSQL_PASSWORD", "")),
            connect_timeout=5)
        try:
            with conn.cursor() as cur:
                cur.execute("SELECT VERSION()")
                return str(cur.fetchone()[0])
        finally:
            conn.close()
    except Exception:  # no server, no driver or no right to ask: the id says so and goes on
        return "unknown"


def _redis_version(environ):
    url = environ.get("PLANETGEN_TEST_REDIS_URL")
    if not url:
        return "none"
    try:
        import redis
        client = redis.Redis.from_url(url, socket_timeout=3, socket_connect_timeout=3)
        return str(client.info("server")["redis_version"])
    except Exception:
        return "unknown"


def _environment(environ):
    items = sorted(f"{key}={value}" for key, value in environ.items()
                   if key.startswith(("PLANETGEN_", "PYTEST_")) and key not in IGNORED_ENV)
    return _short("\n".join(items))


def selection_base(root, base):
    """The commit `--changed` compares with, as a hash: `base` resolved, or where the branch left origin/main."""
    if base:
        return _git(root, "rev-parse", base).strip() or base
    return _git(root, "merge-base", "HEAD", "origin/main").strip() or "none"


def build_parts(root, command_key, pytest_args, base="", environ=None, now=None, servers=None):
    """{part name: value} for the run `scripts/test.py <command_key>` is about to make. `servers` (a callable
    returning `{"mysql": ..., "redis": ...}`) replaces asking the servers, for tests."""
    environ = os.environ if environ is None else environ
    now = now or datetime.datetime.now(datetime.timezone.utc)
    found = servers() if servers else {"mysql": _mysql_version(environ), "redis": _redis_version(environ)}
    return {
        "tree": _tree(root), "dirty": _dirty(root), "packages": _packages(), "python": sys.version.split()[0],
        "mysql": found["mysql"], "redis": found["redis"],
        "command": _short(json.dumps([command_key, list(pytest_args)])), "selection": selection_base(root, base),
        "env": _environment(environ), "date": now.strftime("%Y-%m-%d"),
    }


def build_id(parts):
    return hashlib.sha256(json.dumps(parts, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def differences(then, now):
    """The names of the parts that differ between two runs' parts, in a fixed order."""
    return [name for name in sorted(set(then) | set(now)) if then.get(name) != now.get(name)]
