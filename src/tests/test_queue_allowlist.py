# tests/test_queue_allowlist.py

"""
SEC.45: the template runs only the functions the program queues. A job in Redis names the function the template calls; with pickle gone
(SEC.39) the name is still data, so any name that is not on `planetgen.queue.allowlist` is refused (OPS.235: RQ and its `GuardedJob` are gone, `check_function` is the guard).

- every place that queues work names a function on the list (a sweep of the source, so a new kind of work fails here until it is
  added);
- every function on the lists exists;
- a job naming a function off the list is refused and runs nothing, in `check_function` and through `api_jobs.execute`.
"""

import ast
import importlib
import os

import pytest

from planetgen.queue import allowlist, api_jobs

SOURCE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "planetgen")


def _trees():
    for folder, _dirs, files in os.walk(SOURCE):
        for name in sorted(files):
            if name.endswith(".py"):
                path = os.path.join(folder, name)
                module = os.path.relpath(path, os.path.dirname(SOURCE))[:-3].replace(os.sep, ".")
                yield module, path, ast.parse(open(path, encoding="utf-8").read())


def _call_name(call):
    function = call.func
    return function.attr if isinstance(function, ast.Attribute) else getattr(function, "id", None)


def _resolve(module_name, node):
    """`module:name` of the function a call's argument names (`api_jobs.regenerate_sector`, or a name the module defines or
    imports), or `None` for a parameter the caller passes on."""
    module = importlib.import_module(module_name)
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
        owner = getattr(module, node.value.id, None)
        target = getattr(owner, node.attr, None) if owner is not None else None
    elif isinstance(node, ast.Name):
        target = getattr(module, node.id, None)
    else:
        return None
    if target is None or not callable(target):
        return None
    return f"{target.__module__}:{target.__qualname__}"


API_SITES = {"submit": 0, "run_queued": 0, "_queued_edit": 2}
"""dict: Calls that hand a function to `api_jobs.submit`, and the position of that function among their arguments."""


def _api_calls():
    found = []
    for module, path, tree in _trees():
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or _call_name(node) not in API_SITES:
                continue
            name = _call_name(node)
            if name == "submit" and not (isinstance(node.func, ast.Attribute) and ast.unparse(node.func.value) == "api_jobs"
                                         or isinstance(node.func, ast.Name) and module.endswith("api_jobs")):
                continue                  # another submit (a thread pool's)
            position = API_SITES[name]
            if len(node.args) > position:
                found.append((module, path, node.lineno, node.args[position]))
    return found


def test_every_function_the_api_queues_is_on_the_list():
    seen = set()
    for module, path, line, argument in _api_calls():
        reference = _resolve(module, argument)
        if reference is None:
            assert isinstance(argument, ast.Name) and argument.id == "function", f"{path}:{line} queues {ast.unparse(argument)}"
            continue                      # a helper's parameter: its callers are in the sweep too
        seen.add(reference)
        assert reference in allowlist.API_FUNCTIONS, f"{path}:{line} queues {reference}, which is not on API_FUNCTIONS"
    assert seen, "the sweep found no call that queues an API job"
    assert seen <= allowlist.API_FUNCTIONS


def test_every_function_on_the_lists_exists():
    for reference in allowlist.API_FUNCTIONS:
        module, _, name = reference.partition(":")
        assert callable(getattr(importlib.import_module(module), name)), reference
    for dotted in allowlist.FUNCTIONS:
        module, _, name = dotted.rpartition(".")
        assert callable(getattr(importlib.import_module(module), name)), dotted


def test_the_names_the_program_gives_its_jobs_are_on_the_list():
    from planetgen.web.lib import tileprefetch
    assert tileprefetch.JOB in allowlist.API_FUNCTIONS        # a light job of the maintenance lane, called by `api_jobs.execute`
    assert "planetgen.queue.lightjobs.run" in allowlist.FUNCTIONS        # the process a light job runs in; the API function it calls is checked by `api_jobs.execute`


# --- The guard itself ----------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["os.system", "os.getcwd", "builtins.eval", "subprocess.run", "planetgen.queue.units.run_unit_x",
                                  "planetgen.queue.units", "", "planetgen.cli.worker.main", None, 5])
def test_a_job_naming_a_function_off_the_list_is_refused_and_runs_nothing(name):
    with pytest.raises(allowlist.ForbiddenFunction):
        allowlist.check_function(name)


def test_a_job_naming_a_function_on_the_list_is_let_through():
    for name in allowlist.FUNCTIONS:
        assert allowlist.check_function(name) == name


def test_the_api_dispatcher_calls_only_the_api_functions(monkeypatch):
    ran = []
    monkeypatch.setattr("os.getcwd", lambda: ran.append("os.getcwd") or "/")
    for reference in ("os:getcwd", "os:system", "builtins:exit", "planetgen.queue.units:run_unit", "no colon", "", None, 5):
        with pytest.raises(allowlist.ForbiddenFunction):
            api_jobs.execute(reference, [])
    assert ran == []
    reference = api_jobs.function_reference(api_jobs.run_command)
    assert reference in allowlist.API_FUNCTIONS
    assert "result" in api_jobs.execute(reference, [["true"], dict(os.environ), "/", 30, False])     # a listed one still runs


def test_the_refusal_is_logged(monkeypatch):
    from planetgen.util import log
    seen = []
    monkeypatch.setattr(log, "error", lambda message, *a, **k: seen.append(str(message)))
    with pytest.raises(allowlist.ForbiddenFunction):
        allowlist.check_function("os.system")
    assert seen and "os.system" in seen[0] and "not a function this program queues" in seen[0]


# --- The worker's own environment can add functions (the test suite does) -----------------------------------------------

def test_a_function_added_in_the_workers_environment_is_allowed_and_nothing_else(monkeypatch):
    monkeypatch.setenv(allowlist.EXTRA_ENV, "math:pow, time:*, nonsense, :x, y:")
    assert allowlist.check_function("math.pow") == "math.pow"
    assert allowlist.check_function("time.sleep") == "time.sleep"
    assert allowlist.check_reference("math:pow") == "math:pow"
    for refused in ("math.sqrt", "os.system", "mathpow", "x", "nonsense"):
        with pytest.raises(allowlist.ForbiddenFunction):
            allowlist.check_function(refused)
    with pytest.raises(allowlist.ForbiddenFunction):
        allowlist.check_reference("math:sqrt")


def test_without_the_environment_variable_only_the_lists_apply(monkeypatch):
    monkeypatch.delenv(allowlist.EXTRA_ENV, raising=False)
    assert allowlist.extra_functions() == []
    with pytest.raises(allowlist.ForbiddenFunction):
        allowlist.check_function("math.pow")
