# planetgen/queue/allowlist.py

"""
Which functions a queued job may name (SEC.45).

A job in Redis names the function a worker calls (`planetgen.queue.units.run_unit`), and RQ imports and runs whatever name it
finds. Pickle is gone from every queue (SEC.39), but the name is still data in Redis, so anything that can write there could
get a worker to run any function the worker can import. The Redis ACL user and its password (`examples/redis/planetgen.conf`)
keep other processes out; this list is the second wall: a worker refuses a job whose function is not one the program queues,
logs it, and fails it without running anything.

Two lists, both fixed, both tested against every place the program queues work (`tests/test_queue_allowlist.py`):

- `FUNCTIONS`: the dotted names a job (or one of its callbacks) names. RQ imports these.
- `API_FUNCTIONS`: `module:name` references that the job `planetgen.queue.api_jobs.execute` calls for the API's long
  operations. That one function is a dispatcher, so what it may call is a list of its own.

To queue a new kind of work, add its function to the right list; a test fails until it is there.

What this does not do: it does not check the arguments of an allowed function. `api_jobs.run_command` runs the command line it is
given, so someone who can write to Redis and reaches it still runs that command; the Redis password is what keeps them out.
"""

import os

import rq

EXTRA_ENV = "PLANETGEN_QUEUE_EXTRA_FUNCTIONS"
"""str: A worker's own environment variable that adds functions to both lists: `module:name` entries, or `module:*` for every
function of a module, separated by commas. Unset in production; the test suite sets it for the plain functions its queue tests
run (`math:pow`). Only whoever starts the worker can set it, and that process could run any code anyway; a writer to Redis
cannot."""

FUNCTIONS = frozenset({
    "planetgen.queue.units.run_unit",                   # a run's unit of work (PERF.95)
    "planetgen.queue.units.unit_failed",                # its failure callback
    "planetgen.queue.requires.job_failed",              # the callbacks every job gets (PERF.247)
    "planetgen.queue.requires.job_succeeded",
    "planetgen.queue.api_jobs.execute",                 # the API's long operations (dispatches through `API_FUNCTIONS`)
    "planetgen.web.job_runner.run",                     # a Generate page job
    "planetgen.web.lib.tileprefetch.build_nodes",       # a visitor's node prefetch and rebuild
})
"""frozenset: The dotted names a queued job or callback may name."""

API_FUNCTIONS = frozenset({
    "planetgen.api.edits:change_star_job",
    "planetgen.api.edits:delete_phenomenon_job",
    "planetgen.api.edits:delete_sector_job",
    "planetgen.api.edits:edit_body_job",
    "planetgen.api.edits:regenerate_phenomenon_job",
    "planetgen.api.routes:create_system_job",
    "planetgen.api.routes:generate_uncharted_system_job",
    "planetgen.api.routes:update_system_job",
    "planetgen.queue.api_jobs:build_gravity_grid",
    "planetgen.queue.api_jobs:generate_block",
    "planetgen.queue.api_jobs:generate_neighborhood",
    "planetgen.queue.api_jobs:regenerate_sector",
    "planetgen.queue.api_jobs:run_command",
    "planetgen.queue.pagejobs:build_answer",
    "planetgen.web.lib.estimatejob:work",               # the Generate size estimate (ADM.141)
})
"""frozenset: The `module:name` references `api_jobs.execute` may call."""


class ForbiddenFunction(Exception):
    """A queued job named a function that is not on the list; nothing was run."""


def extra_functions():
    """The `(module, name)` pairs `EXTRA_ENV` adds (`name` is `*` for the whole module)."""
    found = []
    for entry in os.environ.get(EXTRA_ENV, "").split(","):
        module, _, name = entry.strip().partition(":")
        if module and name:
            found.append((module, name))
    return found


def _is_extra(module, name):
    return any(module == m and n in ("*", name) for m, n in extra_functions())


def _refuse(what, name):
    from planetgen.util import log
    message = f"refused: a queued job named {what} {name!r}, which is not a function this program queues"
    log.error(message)
    raise ForbiddenFunction(message)


def check_function(name, what="the function"):
    """`name` when a job may name it, else raises `ForbiddenFunction` (and logs it)."""
    module, _, function = name.rpartition(".") if isinstance(name, str) else ("", "", "")
    if not isinstance(name, str) or (name not in FUNCTIONS and not _is_extra(module, function)):
        _refuse(what, name)
    return name


def check_reference(reference):
    """`reference` (`module:name`) when `api_jobs.execute` may call it, else raises `ForbiddenFunction` (and logs it)."""
    module, _, function = reference.partition(":") if isinstance(reference, str) else ("", "", "")
    if not isinstance(reference, str) or (reference not in API_FUNCTIONS and not _is_extra(module, function)):
        _refuse("the API function", reference)
    return reference


class GuardedJob(rq.job.Job):
    """
    RQ's job, except that the function it runs and its callbacks must be on `FUNCTIONS`, whatever the job in Redis says. A
    refused job fails in the worker's usual way (the failed registry holds the reason) and runs nothing. Every queue
    (`redisqueue.TimedQueue`) and worker (`cli/worker`) uses this class.
    """

    @property
    def func(self):
        if getattr(self, "instance", None) is not None:
            _refuse("a method of an object", self.func_name)
        check_function(self.func_name)
        return super().func

    @property
    def success_callback(self):
        if self._success_callback_name:
            check_function(self._success_callback_name, "the success callback")
        return super().success_callback

    @property
    def failure_callback(self):
        if self._failure_callback_name:
            check_function(self._failure_callback_name, "the failure callback")
        return super().failure_callback

    @property
    def stopped_callback(self):
        if self._stopped_callback_name:
            check_function(self._stopped_callback_name, "the stopped callback")
        return super().stopped_callback
