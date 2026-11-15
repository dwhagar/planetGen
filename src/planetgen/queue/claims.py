# planetgen/queue/claims.py

"""
What each kind of job claims of the host (PERF.254): the memory it is expected to hold, the cores it keeps busy and the lane it is on,
worked out from the RQ job itself (`claim_for`) for the memory gate (`planetgen.queue.memorygate`).

**No job is exempt.** Every function a queue may run (`allowlist.FUNCTIONS`) is here, either with a claim of its own
(`FUNCTION_CLAIMS`) or because it is a callback that runs inside the job it follows (`CALLBACKS`); every function the API's dispatcher may
call (`allowlist.API_FUNCTIONS`) has a claim (`API_CLAIMS`); every command with units (`commands.REGISTRY`) has the figure of its stage
(`workers.STAGE_WORKER_MEMORY_MB`, measured) or one by its nature (`UNIT_BY_NATURE_MB`). `tests/test_memory_gate.py` fails for a function,
an API reference or a command that is on its list and not here, so a new kind of job cannot slip past the gate.

The figures are first claims. The measured ones (the stage table, from `units.peak_resident_mb`) are what a worker and its horse held, shared
pages counted once, with a margin; the others are marked `unmeasured` and stand until the gate has seen a job of the kind (it keeps the
largest peak, `memorygate.learn`, and claims that with a margin). A stand-in is not a measurement: the registry of stand-ins does not list
them because they are replaced by a measurement in the first run, and `memorygate.snapshot` shows what was learned.
"""

from planetgen.queue import commands, memorygate, pause, redisqueue, workers
from planetgen.queue.memorygate import Claim
from planetgen.util import log

UNIT_FUNCTION = "planetgen.queue.units.run_unit"
OPEN_FUNCTION = "planetgen.queue.units.open_stage"
EXECUTE_FUNCTION = "planetgen.queue.api_jobs.execute"

UNIT_CPU = 1.0
"""float: The cores a unit keeps busy: a unit is one process computing."""

CUT_MB = 1500
"""int: What the cut of a stage (the opening, `wrapper.advance`) claims before a measurement says otherwise: Lane 13 measured 164 MB for the
process, 166 MB for the 573,713 targets of the backfill and 940 MB after the first sky plan (2026-11-08), so 1,500 MB covers the one
full cut seen. It claims the bulk lane to itself (`Claim.exclusive`) and is learned per command it opens (`cut:<command>`)."""

SUPERVISOR_MB = 700
"""int: A job that supervises (a run's starter inside an API job, the Generate page's job): its own process is small; the units it queues and the
processes it starts hold the memory, and the gate counts them under the worker (`memorygate.held_mb` takes every process under it)."""

UNIT_BY_NATURE_MB = {
    "plan-skeleton": 1000, "finish-phenomena": 1000, "refresh-run-blocks": 1000, "clean-nebula-field": 700, "clear-scatter": 700,
    "generate-system": workers.BULK_WORKER_MEMORY_MB,       # the largest worker seen held 1.9 GB (measured with no sampling tables built, PERF.266)
    "run-population": 2100, "run-territories": 1000, "reset-galaxy": 700,
    "check-math": 700, "check-datasets": 700, "check-revision": 700, "check-tables": 700, "check-tables-extended": 1000, "check-orphans": 1000,
    "check-ids": 700, "check-values": 1000, "check-counts": 700, "check-version-keys": 700, "check-storage": 700, "check-systems": 1000,
    "check-object-index": 1000, "warm-object-index": 1000, "warm-galaxy-map": 1000, "warm-lists": 700, "warm-search": 1000,
    "end-map-rebuild-maintenance": 700,
}
"""dict: What a unit of a command that is one unit by nature (`CommandSpec.whole`) claims, MB, where `workers.STAGE_WORKER_MEMORY_MB` has no figure.
Unmeasured: first claims by what the command reads, until the gate has learned the kind."""

FUNCTION_CLAIMS = {
    "planetgen.web.job_runner.run": (SUPERVISOR_MB, 0.5),
    "planetgen.web.lib.tileprefetch.build_nodes": (700, 1.0),
}
"""dict: The functions a queue runs that need no more than a figure of their own: `{dotted name: (MB, cores)}`. Unmeasured."""

SPECIAL_FUNCTIONS = frozenset({UNIT_FUNCTION, OPEN_FUNCTION, EXECUTE_FUNCTION})
"""frozenset: The functions whose claim comes from their arguments: a unit's command, the command a cut opens, the API function `execute` calls."""

CALLBACKS = frozenset({
    "planetgen.queue.units.unit_failed",
    "planetgen.queue.units.open_failed",
    "planetgen.queue.requires.job_failed",
    "planetgen.queue.requires.job_succeeded",
})
"""frozenset: The functions that run as a job's callbacks, in the horse of the job they follow: they are inside its claim and take no claim of their own."""

API_CLAIMS = {
    "planetgen.api.edits:change_star_job": (workers.BULK_WORKER_MEMORY_MB, 1.0),
    "planetgen.api.edits:delete_phenomenon_job": (700, 0.5),
    "planetgen.api.edits:delete_sector_job": (700, 0.5),
    "planetgen.api.edits:edit_body_job": (700, 0.5),
    "planetgen.api.edits:regenerate_phenomenon_job": (700, 0.5),
    "planetgen.api.routes:create_system_job": (workers.BULK_WORKER_MEMORY_MB, 1.0),
    "planetgen.api.routes:generate_uncharted_system_job": (workers.BULK_WORKER_MEMORY_MB, 1.0),
    "planetgen.api.routes:update_system_job": (workers.BULK_WORKER_MEMORY_MB, 1.0),
    "planetgen.queue.api_jobs:build_gravity_grid": (1000, 1.0),
    "planetgen.queue.api_jobs:generate_block": (SUPERVISOR_MB, 0.5),
    "planetgen.queue.api_jobs:generate_neighborhood": (SUPERVISOR_MB, 0.5),
    "planetgen.queue.api_jobs:regenerate_sector": (workers.SECTOR_WORKER_MEMORY_MB, 1.0),
    "planetgen.queue.api_jobs:run_command": (SUPERVISOR_MB, 0.5),
    "planetgen.queue.pagejobs:build_answer": (700, 1.0),
    "planetgen.web.lib.estimatejob:work": (SUPERVISOR_MB, 0.5),
}
"""dict: What each function the API's dispatcher (`api_jobs.execute`) may call claims, `{module:name: (MB, cores)}`. Unmeasured, except that a
sector or a star drawn in the job is the stage figure of the same work."""

DEFAULT_MB = 1000
"""int: What the worker claims for a job whose kind has no figure (a gap `tests/test_memory_gate.py` is there to close): it is let through, with this claim and an error in the log."""


class NoClaim(LookupError):
    """A job has no claim: its function, API reference or command is on no table here (the gate refuses to guess)."""


def lane_of(queue_name):
    """The lane of the queue `queue_name`: bulk (a generation run's units, and the first generation queue), maintenance (the three queues of the
    maintenance lane), interactive, long, or other."""
    if queue_name in (redisqueue.BULK_QUEUE, redisqueue.GENERATION_QUEUE):
        return memorygate.LANE_BULK
    if queue_name in redisqueue.maintenance_lane():
        return memorygate.LANE_MAINTENANCE
    if queue_name == redisqueue.INTERACTIVE_QUEUE:
        return memorygate.LANE_INTERACTIVE
    if queue_name == redisqueue.LONG_QUEUE:
        return memorygate.LANE_LONG
    return memorygate.LANE_OTHER


def unit_table_mb(kind):
    """What a unit of command `kind` claims before anything is learned, MB: the stage's measured figure, else the one by its nature.

    Raises:
        NoClaim: The command is inline or unknown, or has neither.
    """
    if kind not in commands.REGISTRY or commands.spec(kind).inline:
        raise NoClaim(f"{kind!r} is not a command with units")
    figure = workers.stage_worker_mb(kind) or UNIT_BY_NATURE_MB.get(kind)
    if figure is None:
        raise NoClaim(f"the command {kind!r} has no memory figure: add it to workers.STAGE_WORKER_MEMORY_MB or claims.UNIT_BY_NATURE_MB")
    return figure


def api_table(reference):
    """`(MB, cores)` of the API function `reference` (`module:name`).

    Raises:
        NoClaim: It is not in `API_CLAIMS`.
    """
    if reference not in API_CLAIMS:
        raise NoClaim(f"the API function {reference!r} has no claim: add it to claims.API_CLAIMS")
    return API_CLAIMS[reference]


def function_table(function):
    """`(MB, cores)` of the queued function `function` (a dotted name).

    Raises:
        NoClaim: It is not in `FUNCTION_CLAIMS`.
    """
    if function not in FUNCTION_CLAIMS:
        raise NoClaim(f"the queued function {function!r} has no claim: add it to claims.FUNCTION_CLAIMS (or CALLBACKS, if it runs inside another job)")
    return FUNCTION_CLAIMS[function]


def typed_mb(kind):
    """The typed figure of a claim kind (`unit:<command>`, `cut:<command>`, `api:<ref>`, `fn:<name>`), MB, or `None` for a kind with none (OPS.158: the page sets it beside the learned one)."""
    try:
        prefix, _, name = kind.partition(":")
        if prefix == "unit":
            return float(unit_table_mb(name))
        if prefix == "cut":
            return float(CUT_MB)
        if prefix == "api":
            return float(api_table(name)[0])
        if prefix == "fn":
            return float(function_table(name)[0])
    except NoClaim:
        pass
    return None


def default_claim(job):
    """The claim of a job of a kind with no figure: `DEFAULT_MB` and a core, on the job's lane, under the kind `fn:<function>`."""
    return Claim(f"fn:{job.func_name}", float(DEFAULT_MB), 1.0, lane_of(job.origin))


def cap_of(kind):
    """The most workers of claim kind `kind` the old static table allows (`workers.stage_worker_cap`), `None` for a kind it never limited."""
    if kind.startswith("unit:"):
        return workers.stage_worker_cap(kind.partition(":")[2])
    return None


def cut_command(conn, run_id, index):
    """The command that stage `index` of the run runs (what a cut of it opens), or `None` when the run has no such stage. Reads the step as the run
    record holds it: nothing is resolved, so a stage whose parameters wait for an earlier one is named all the same."""
    import json
    from planetgen.queue import runrecord
    try:
        declared = json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")))[index]
        return next(commands.Step.from_dict(item).kind for item in runrecord.read_steps(conn, run_id) if item["key"] == declared["key"])
    except Exception as exc:  # noqa: BLE001 -- a run that is gone or a stage that is not there is cut under the default claim
        log.swallowed("queue/claims", "cut_command", exc)
        return None


def claim_for(job, conn):
    """
    The claim of the RQ job `job` (any queue): the figure of its kind (`kind` below), raised to what the gate has learned of the kind and
    to the larger claim a retry was given (`memorygate.claim_mb`).

    Kinds: `unit:<command>` (a run's unit), `cut:<command>` (the opening of a stage that runs that command), `api:<module:name>` (an API job)
    and `fn:<dotted name>` (any other queued function).

    Raises:
        NoClaim: The job is of a kind with no figure.
    """
    function = job.func_name
    lane = lane_of(job.origin)
    if function == UNIT_FUNCTION:
        kind = job.args[3]
        mb, cpu, exclusive = unit_table_mb(kind), UNIT_CPU, False
        kind = f"unit:{kind}"
    elif function == OPEN_FUNCTION:
        command = cut_command(conn, job.args[0], job.args[1])
        mb, cpu, exclusive = CUT_MB, 1.0, True
        kind = f"cut:{command or 'unknown'}"
    elif function == EXECUTE_FUNCTION:
        reference = job.args[0]
        (mb, cpu), exclusive = api_table(reference), False
        kind = f"api:{reference}"
    else:
        (mb, cpu), exclusive = function_table(function), False
        kind = f"fn:{function}"
    solo = pause.solo_of(conn)
    token = job.args[0] if function in (UNIT_FUNCTION, OPEN_FUNCTION) else None        # OPS.156: a run holds one figure per kind
    return Claim(kind, memorygate.claim_mb(conn, kind, mb, job_id=job.id, token=token), cpu, lane, exclusive, solo=bool(solo and solo["job"] == job.id))
