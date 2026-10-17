# tests/reuse.py

"""
Resuming across builds (TEST.208, on top of TEST.148 to TEST.152).

`scripts/test.py --resume` used to continue only a run of the identical build id, and that id hashes the whole
checkout, so one merge of main after a cancelled gate started the next from zero. Here a test that passed in an
earlier run is reused (left out of this one) when nothing it depends on has changed since:

- the environment is the same: installed packages, Python, MySQL and Redis versions, the `PLANETGEN_*` and
  `PYTEST_*` variables and the UTC day (`ENVIRONMENT_PARTS`; a different environment reuses nothing);
- none of the files that differ between the earlier run's tree and this one is one the test used. The files a test
  used are the ones the dependency map records (`tests/depmap.py`, written to the progress log as a `deps` line when
  the test passes), and the same rules as `scripts/test.py changed` decide: a changed test file, a file the test
  ran or opened, a file that imports one it ran, anything under `src/planetgen/`, `src/html/` or `scripts/` for a
  test that starts a subprocess, and everything for a change to the test setup (`conftest.py`, the support
  modules, fixtures, `pytest.ini`, the lock files).

Each run's header keeps a manifest of the tracked and untracked (not git-ignored) files, `path -> hash`, so two
runs can be compared without git still holding the earlier commit.
"""

import json
import os
import subprocess

from tests import buildid, depmap, progress_log

ENVIRONMENT_PARTS = ("packages", "python", "mysql", "redis", "env", "date")
"""tuple: The parts of a build id (`tests/buildid.py`) that must be the same for an earlier pass to count."""


def _git(root, *args):
    result = subprocess.run(["git", *args], cwd=root, capture_output=True)
    return result.stdout if result.returncode == 0 else b""


def manifest(root):
    """`{repository-relative path: hash}` of every tracked file and every untracked file git does not ignore, as
    they are on disk: a file git has unmodified by its blob id, any other by its content."""
    listing = _git(root, "ls-files", "-s", "-z").decode("utf-8", "replace")
    tracked = {}
    for entry in listing.split("\0"):
        if entry:
            meta, path = entry.split("\t", 1)
            tracked[path] = meta.split()[1]
    status = _git(root, "status", "--porcelain=v1", "-z", "--untracked-files=all").decode("utf-8", "replace")
    pieces = [piece for piece in status.split("\0") if piece]
    touched, index = set(), 0
    while index < len(pieces):
        code, path = pieces[index][:2], pieces[index][3:]
        touched.add(path)
        if code[0] in "RC":     # a rename or copy names the old path in the next piece
            index += 1
            touched.add(pieces[index])
        index += 1
    result = {}
    for path, blob in tracked.items():
        if path not in touched:
            result[path] = blob
    for path in touched:
        full = os.path.join(root, path)
        if os.path.isfile(full):
            result[path] = "d:" + buildid._file_hash(full)
    return result


def changed_between(then, now):
    """The paths added, removed or different between two manifests."""
    return {path for path in set(then) | set(now) if then.get(path) != now.get(path)}


def read_deps(path):
    """`{node id: {"files": [...], "opaque": bool}}` from a progress log's `deps` lines (the last one of each test)."""
    found = {}
    try:
        with open(path, encoding="utf-8") as f:
            for text in f:
                try:
                    line = json.loads(text)
                except ValueError:
                    continue
                if line.get("kind") == "deps":
                    found[line["nodeid"]] = line
    except OSError:
        pass
    return found


def _logs(directory, skip=None):
    """`[(path, header, tests, deps)]` of the progress logs in `directory` that have a manifest, newest first."""
    found = []
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return found
    for name in names:
        path = os.path.join(directory, name)
        if not name.endswith(".jsonl") or path == skip:
            continue
        header, tests = progress_log.read(path)
        if header.get("manifest") is not None and header.get("parts"):
            found.append((path, header, tests, read_deps(path)))
    found.sort(key=lambda item: (item[1].get("started", ""), os.path.getmtime(item[0])), reverse=True)
    return found


def environment_differences(then, now):
    """The environment parts that differ between two runs' build parts."""
    return [name for name in ENVIRONMENT_PARTS if then.get(name) != now.get(name)]


def reusable(directory, parts, current, skip=None):
    """
    The passes of earlier runs that this run may count: `{node id: {"test": its progress line, "deps": its deps
    line, "build": the run's build id, "started": when that run started}}`, the newest run's pass of each test, and
    `(that, the names of the runs looked at, the environment parts that kept the others out)`.

    Args:
        directory (str): The progress logs.
        parts (dict): This run's build parts (`buildid.build_parts`).
        current (dict): This run's `manifest`.
        skip (str, optional): A log not to read (this run's own).
    """
    found, looked, blocked = {}, [], []
    for path, header, tests, deps in _logs(directory, skip):
        different = environment_differences(header["parts"], parts)
        if different:
            blocked.extend(name for name in different if name not in blocked)
            continue
        looked.append(os.path.basename(path))
        passed = [nodeid for nodeid, line in tests.items()
                  if line.get("outcome") == "passed" and nodeid in deps and nodeid not in found]
        if not passed:
            continue
        changed = changed_between(header["manifest"], current)
        names = sorted({name for nodeid in passed for name in deps[nodeid].get("files", [])})
        number = {name: i for i, name in enumerate(names)}
        data = {"format": depmap.FORMAT, "files": names, "tests": {
            nodeid: {"f": [number[name] for name in deps[nodeid].get("files", [])], "o": bool(deps[nodeid].get("opaque"))}
            for nodeid in passed}}
        must_run, _reason = depmap.select(data, changed, passed)
        for nodeid in passed:
            if nodeid not in must_run:
                found[nodeid] = {"test": tests[nodeid], "deps": deps[nodeid], "build": header.get("build"),
                                 "started": header.get("started", "")}
    return found, looked, blocked
