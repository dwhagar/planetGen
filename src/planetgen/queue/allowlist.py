# planetgen/queue/allowlist.py

"""
Which functions a queued job may name (SEC.45).

A forked job (a unit, a page step, a light job) names the function the template process calls (`planetgen.queue.forkedjob.run`), and that name is data in Redis (OPS.235: RQ itself is gone, the
template imports and runs the name). Pickle is gone from every queue (SEC.39), but anything that can write there could get the template to run any function it can import. The Redis ACL user and its password (`examples/redis/planetgen.conf`)
keep other processes out; this list is the second wall: the template refuses a job whose function is not one the program queues,
logs it, and runs nothing (`check_function`).

Two lists, both fixed, both tested against every place the program queues work (`tests/test_queue_allowlist.py`):

- `FUNCTIONS`: the dotted names a job names. The template imports these.
- `API_FUNCTIONS`: `module:name` references that the job `planetgen.queue.api_jobs.execute` calls for the API's long
  operations. That one function is a dispatcher, so what it may call is a list of its own.

To queue a new kind of work, add its function to the right list; a test fails until it is there.

What this does not do: it does not check the arguments of an allowed function. `api_jobs.run_command` runs the command line it is
given, so someone who can write to Redis and reaches it still runs that command; the Redis password is what keeps them out.
"""

import os

EXTRA_ENV = "PLANETGEN_QUEUE_EXTRA_FUNCTIONS"
"""str: A worker's own environment variable that adds functions to both lists: `module:name` entries, or `module:*` for every
function of a module, separated by commas. Unset in production; the test suite sets it for the plain functions its queue tests
run (`math:pow`). Only whoever starts the worker can set it, and that process could run any code anyway; a writer to Redis
cannot."""

FUNCTIONS = frozenset({
    "planetgen.queue.forkedjob.run",                    # a unit as a job process the template forked (OPS.226)
    "planetgen.web.warmup.warm_app",                    # the warm set of the maintenance lane's template (OPS.229): the site's app, built once
    "planetgen.queue.units.warm_star_tables",           # the warm set of a template (OPS.226): the sampling tables every job forked from it shares
    "planetgen.queue.lightjobs.run",                    # a light job as a process the template forked (OPS.230): the API's long operations, dispatched through `API_FUNCTIONS`
})
"""frozenset: The dotted names a queued job or callback may name."""

API_FUNCTIONS = frozenset({
    "planetgen.api.edits:change_star_job",
    "planetgen.api.edits:delete_phenomenon_job",
    "planetgen.api.edits:delete_sector_job",
    "planetgen.api.edits:edit_body_job",
    "planetgen.api.edits:regenerate_phenomenon_job",
    "planetgen.api.routes:create_system_job",
    "planetgen.api.lookup:name_search_job",
    "planetgen.api.routes:generate_uncharted_system_job",
    "planetgen.api.routes:update_system_job",
    "planetgen.queue.api_jobs:build_gravity_grid",
    "planetgen.queue.api_jobs:generate_block",
    "planetgen.queue.api_jobs:generate_neighborhood",
    "planetgen.queue.api_jobs:regenerate_sector",
    "planetgen.queue.api_jobs:run_command",
    "planetgen.queue.pagejobs:build_answer",
    "planetgen.web.system_page:generate_in_process",    # the one-off system page (OPS.244)
    "planetgen.web.lib.estimatejob:work",               # the Generate size estimate (ADM.141)
    "planetgen.web.lib.tileprefetch:build_nodes",       # a visitor's node prefetch and rebuild (OPS.230: a light job on the maintenance lane)
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


def _checked(name, allowed, split, what):
    """`name` when it is in `allowed` or `EXTRA_ENV` adds it (`split` cuts it into module and function), else `_refuse`."""
    module, _, function = split(name) if isinstance(name, str) else ("", "", "")
    if not isinstance(name, str) or (name not in allowed and not _is_extra(module, function)):
        _refuse(what, name)
    return name


def check_function(name, what="the function"):
    """`name` when a job may name it, else raises `ForbiddenFunction` (and logs it)."""
    return _checked(name, FUNCTIONS, lambda text: text.rpartition("."), what)


def check_reference(reference):
    """`reference` (`module:name`) when `api_jobs.execute` may call it, else raises `ForbiddenFunction` (and logs it)."""
    return _checked(reference, API_FUNCTIONS, lambda text: text.partition(":"), "the API function")


