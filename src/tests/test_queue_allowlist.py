# tests/test_queue_allowlist.py

"""
SEC.45: a worker runs only the functions the program queues. A job in Redis names the function a worker calls; with pickle gone
(SEC.39) the name is still data, so the worker refuses any name that is not on `planetgen.queue.allowlist`.

- every place that queues work names a function on the list (a sweep of the source, so a new kind of work fails here until it is
  added);
- every function on the lists exists;
- a job, or a callback, naming a function off the list is refused and runs nothing, in `GuardedJob` and through
  `api_jobs.execute`;
- through a real Redis and a worker: a refused job ends in the failed registry with the reason, and nothing ran.
"""

import ast
import importlib
import os

import pytest

from planetgen.queue import allowlist, api_jobs, redisqueue, requires, units

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
    for name in (units.UNIT_FUNCTION, units.FAILURE_FUNCTION, units.OPEN_FUNCTION, units.OPEN_FAILURE_FUNCTION, requires.FAILED_FUNCTION,
                 requires.SUCCEEDED_FUNCTION):
        assert name in allowlist.FUNCTIONS
    from planetgen.web.lib import tileprefetch
    assert tileprefetch.JOB in allowlist.FUNCTIONS
    assert f"{api_jobs.execute.__module__}.{api_jobs.execute.__name__}" in allowlist.FUNCTIONS


def test_every_enqueue_in_the_source_queues_a_function_on_the_list():
    """`queue.enqueue(X, ...)`, `enqueue_many([prepare_data(X, ...)])` and `requires.enqueue(queue, X, ...)`: X is a constant
    on the list, a module constant that is, `api_jobs.execute`, or a helper's own parameter."""
    known = {"UNIT_FUNCTION": units.UNIT_FUNCTION, "OPEN_FUNCTION": units.OPEN_FUNCTION, "JOB": "planetgen.web.lib.tileprefetch.build_nodes",
             "execute": "planetgen.queue.api_jobs.execute"}
    checked = 0
    for module, path, tree in _trees():
        if module == "planetgen.queue.redisqueue":
            continue                      # the queue class itself, not a place that queues work
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or _call_name(node) not in ("enqueue", "prepare_data"):
                continue
            if not isinstance(node.func, ast.Attribute):
                continue
            owner = ast.unparse(node.func.value)
            position = 1 if owner.endswith("requires") else 0
            if len(node.args) <= position:
                continue
            argument = node.args[position]
            if isinstance(argument, ast.Name) and argument.id == "function" and module == "planetgen.queue.requires":
                continue                  # `requires.enqueue` hands its caller's function on; its callers are swept below
            value = argument.value if isinstance(argument, ast.Constant) else known.get(getattr(argument, "id", None))
            assert value in allowlist.FUNCTIONS, f"{path}:{node.lineno} queues {ast.unparse(argument)}"
            checked += 1
    assert checked >= 4         # the Generate page's job, a node prefetch, a run's units and the API's `execute`


def test_nothing_calls_requires_enqueue_with_a_function_that_is_not_on_the_list():
    for module, path, tree in _trees():
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and _call_name(node) == "enqueue" and ast.unparse(node.func).endswith("requires.enqueue"):
                argument = node.args[1]
                value = argument.value if isinstance(argument, ast.Constant) else getattr(argument, "id", None)
                assert value in allowlist.FUNCTIONS or value in ("UNIT_FUNCTION", "JOB"), f"{path}:{node.lineno}"


# --- The guard itself ----------------------------------------------------------------------------------------------------

class _Connection:
    """RQ's job class only keeps the connection until it is saved; nothing here is."""


def _job(function, **options):
    return allowlist.GuardedJob.create(function, connection=_Connection(), serializer=redisqueue.JsonSerializer, **options)


@pytest.mark.parametrize("name", ["os.system", "os.getcwd", "builtins.eval", "subprocess.run", "planetgen.queue.units.run_unit_x",
                                  "planetgen.queue.units", "", "planetgen.cli.worker.main"])
def test_a_job_naming_a_function_off_the_list_is_refused_and_runs_nothing(name):
    job = _job("os.getcwd")                    # made with a harmless name, then told another, as a writer to Redis would
    job.func_name = name
    with pytest.raises(allowlist.ForbiddenFunction):
        job.func
    with pytest.raises(allowlist.ForbiddenFunction):
        job._execute()                 # what `perform` calls once it has touched Redis


def test_a_job_naming_a_function_on_the_list_resolves_it():
    job = _job(units.UNIT_FUNCTION, args=(0, 0, 0, "x", {}, None))
    assert job.func is units.run_unit


@pytest.mark.parametrize("callback", ["success_callback", "failure_callback", "stopped_callback"])
def test_a_callback_off_the_list_is_refused(callback):
    job = _job(units.UNIT_FUNCTION, args=(0, 0, 0, "x", {}, None))
    setattr(job, f"_{callback}_name", "os.system")
    setattr(job, f"_{callback}", allowlist.rq.job.UNEVALUATED)
    with pytest.raises(allowlist.ForbiddenFunction):
        getattr(job, callback)
    setattr(job, f"_{callback}_name", requires.SUCCEEDED_FUNCTION)
    assert callable(getattr(job, callback))


def test_a_method_job_is_refused():
    job = _job(units.UNIT_FUNCTION, args=(0, 0, 0, "x", {}, None))
    job.instance = object()
    with pytest.raises(allowlist.ForbiddenFunction):
        job.func


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


def test_every_queue_the_site_makes_uses_the_guarded_job_class():
    for name in (*redisqueue.QUEUES, "planetgen-run-0123", "planetgen-web-0123"):
        assert redisqueue.queue(name, connection=object()).job_class is allowlist.GuardedJob


def test_the_worker_is_built_with_the_guarded_job_class():
    source = open(os.path.join(SOURCE, "cli", "worker.py"), encoding="utf-8").read()
    assert "job_class=allowlist.GuardedJob" in source


# --- Through a real Redis ---------------------------------------------------------------------------------------------

def test_a_worker_fails_a_job_that_names_a_function_off_the_list_and_runs_nothing(redis_server, tmp_path):
    import rq
    connection = redisqueue.connect(redis_server)
    name = f"planetgen-allowlist-test-{os.getpid()}"
    marker = tmp_path / "ran"
    queue = rq.Queue(name, connection=connection, serializer=redisqueue.JsonSerializer)
    try:
        # A writer to Redis who is not planetGen: it names a function of its own, with nothing from this program.
        queue.enqueue("pathlib.Path.touch", str(marker), job_id=f"{name}-evil", job_timeout=30)
        queue.enqueue("os.getcwd", job_id=f"{name}-evil2", job_timeout=30,
                      on_failure=rq.job.Callback("os.getcwd"))
        rq.SimpleWorker([queue], connection=connection, serializer=redisqueue.JsonSerializer,
                        job_class=allowlist.GuardedJob).work(burst=True, logging_level="CRITICAL")
        assert not marker.exists()
        for suffix in ("evil", "evil2"):
            job = rq.job.Job.fetch(f"{name}-{suffix}", connection=connection, serializer=redisqueue.JsonSerializer)
            assert job.get_status() == "failed"
            assert "not a function this program queues" in (job.latest_result().exc_string or "")
    finally:
        queue.empty()
        for suffix in ("evil", "evil2"):
            try:
                rq.job.Job.fetch(f"{name}-{suffix}", connection=connection, serializer=redisqueue.JsonSerializer).delete()
            except Exception:  # noqa: BLE001
                pass


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
