# tests/worker_patches.py

"""
Patches that reach the units of a run (PERF.21, TEST.74).

A run's units run in this process under the inline queue (`tests/inline_queue.py`), so
`monkeypatch.setattr` reaches them. `patch_everywhere` keeps the call sites of the tests that were written for
worker processes: it replaces `owner.name` here, by a module-level factory named as `"module:function"` and called as
`factory(real, **params)`; `params` must be JSON. Shared state (call counts, captured calls) lives in files under the
test's `tmp_path`.
"""

import importlib
import json

from planetgen.queue import work as workQueue

import fcntl

class Interrupted(RuntimeError):
    """The fault the patches below inject."""


def _resolve(ref):
    module_name, _, name = ref.partition(":")
    return getattr(importlib.import_module(module_name), name)


def patch_everywhere(monkeypatch, owner, name, factory, **params):
    """Replaces `owner.name` (owner a module) by the factory's replacement until `monkeypatch.undo()` or until the
    returned `undo()` is called (which takes back only this patch)."""
    real = getattr(owner, name)
    monkeypatch.setattr(owner, name, _resolve(factory)(real, **params))

    def undo():
        setattr(owner, name, real)

    return undo


def _lock(f):
    """Locks `f` until it is closed."""
    fcntl.flock(f, fcntl.LOCK_EX)


def bump(path):
    """Adds one to the count in `path`, across processes; returns it."""
    with open(path, "a+", encoding="utf-8") as f:
        _lock(f)
        f.seek(0)
        count = int(f.read() or 0) + 1
        f.seek(0)
        f.truncate()
        f.write(str(count))
        f.flush()
    return count


def count(path):
    try:
        with open(path, encoding="utf-8") as f:
            return int(f.read() or 0)
    except FileNotFoundError:
        return 0


def append_json(path, value):
    """Appends one JSON line to `path`, across processes."""
    with open(path, "a", encoding="utf-8") as f:
        _lock(f)
        f.write(json.dumps(value) + "\n")
        f.flush()


def read_json_lines(path):
    try:
        with open(path, encoding="utf-8") as f:
            return [json.loads(line) for line in f if line.strip()]
    except FileNotFoundError:
        return []


def workers():
    """The worker count generation runs use in this test session
    (`PLANETGEN_WORKERS`, which TEST.74's runs set to 1, 2 or 4)."""
    return workQueue.worker_count()


# --- Factories ----------------------------------------------------------------

def fail_on_call(real, counter, call_number):
    """Raises `Interrupted` on the `call_number`-th call (1-based) counted
    across every process; every other call goes through."""

    def wrapper(*args, **kwargs):
        if bump(counter) == call_number:
            raise Interrupted(f"{real.__name__} interrupted at call {call_number}")
        return real(*args, **kwargs)

    return wrapper
