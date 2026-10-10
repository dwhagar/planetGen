# tests/worker_patches.py

"""
Patches that reach the work queue's worker processes too (PERF.21,
TEST.74).

`monkeypatch.setattr` only changes the test's own process; a run with
two or more workers does its sectors and layers in spawned processes
that import everything afresh. `patch_everywhere` patches this process
and also records the patch in the environment, where every worker that
starts while it is in place picks it up (`workQueue._run_worker_hooks`
calls `install`). `monkeypatch.undo()` takes both back.

A replacement is made by a module-level factory, named as
`"module:function"` and called as `factory(real, **params)`, so a worker
can build it as well; `params` must be JSON. Shared state (call counts,
captured calls) lives in files under the test's `tmp_path`, so every
process sees the same count.
"""

import importlib
import json
import os

from planetgen.queue import work as workQueue

import fcntl

PATCHES_ENV = "PLANETGEN_TEST_WORKER_PATCHES"
HOOK = "tests.worker_patches:install"


class Interrupted(RuntimeError):
    """The fault the patches below inject."""


def _resolve(ref):
    module_name, _, name = ref.partition(":")
    return getattr(importlib.import_module(module_name), name)


def _apply(spec):
    owner = importlib.import_module(spec["owner"])
    real = getattr(owner, spec["name"])
    setattr(owner, spec["name"], _resolve(spec["factory"])(real, **spec["params"]))


def install():
    """Applies every recorded patch (a worker's start, through
    `$PLANETGEN_WORKER_HOOKS`)."""
    for spec in json.loads(os.environ.get(PATCHES_ENV) or "[]"):
        _apply(spec)


def patch_everywhere(monkeypatch, owner, name, factory, **params):
    """Replaces `owner.name` (owner a module) here and in every worker
    started from now until `monkeypatch.undo()` or until the returned
    `undo()` is called (which takes back only this patch)."""
    spec = {"owner": owner.__name__, "name": name, "factory": factory, "params": params}
    real = getattr(owner, name)
    saved = {key: os.environ.get(key) for key in (PATCHES_ENV, workQueue.WORKER_HOOKS_ENV)}
    monkeypatch.setattr(owner, name, _resolve(factory)(real, **params))
    specs = json.loads(os.environ.get(PATCHES_ENV) or "[]") + [spec]
    monkeypatch.setenv(PATCHES_ENV, json.dumps(specs))
    hooks = [hook for hook in (os.environ.get(workQueue.WORKER_HOOKS_ENV) or "").split(",") if hook]
    if HOOK not in hooks:
        monkeypatch.setenv(workQueue.WORKER_HOOKS_ENV, ",".join(hooks + [HOOK]))

    def undo():
        setattr(owner, name, real)
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

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
