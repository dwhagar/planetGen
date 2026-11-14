"""OPS.197: no sleep in command code. A job that sleeps while it waits shows no activity at all, so a wait must sit inside a blocking primitive that says so (`unitbeat.waiting`:
`db/contention.py`, the gate, the sky lock). This test lists every call of `sleep` in the package by the function that makes it. Each must be named in `ALLOWED` with the reason it is
not a silent stretch inside a command; a new one fails here, and an entry that no longer matches a call fails too (reports/job-engine-items-stage-1-2-2026-11-13.md, item P7).

The test finds `time.sleep(...)`, a bare `sleep(...)` (the clock a function takes as a parameter, `sleep=time.sleep`) and `concurrent.futures.wait(...)`, the timed wait the battery's digest loop polls with. A loop that polls without any sleep is out of its reach. Every other
waiting call (`join`, `wait`) is bounded by its own timeout and is not a poll.
"""

import ast
import pathlib

PACKAGE = pathlib.Path(__file__).resolve().parent.parent / "planetgen"

BLOCKING = "a blocking primitive that declares its wait"
RETRY = "a pause of under a second between two tries of a transaction MySQL asked to be run again (bounded by SECTOR_SAVE_ATTEMPTS)"
OUTSIDE = "runs in a process that is not a unit (the web server, the CLI, the starter that follows a run), so no beat is expected of it"

ALLOWED = {
    ("db/contention.py", "_wait_while_closed"): BLOCKING,
    ("db/contention.py", "run_with_retries"): BLOCKING,
    ("queue/memorygate.py", "_decision_lock"): BLOCKING + " (the gate's lock, shown as `gate`)",
    ("queue/memorygate.py", "_acquire"): BLOCKING + " (waiting for memory, shown as `gate`)",
    ("generation/run_plan.py", "_run_stale_unit"): RETRY,
    ("db/store.py", "_save_with_retries"): RETRY,
    ("db/migrations/helpers.py", "_with_retries"): "the migration runner, not a unit; a pause between two tries of one statement",
    ("generation/generator_battery.py", "start_running_digest.finish"): "the battery digest wait of OPS.184: it waits on the digest's own thread and shows its progress line every poll",
    ("datasets/fetch.py", "_download_from"): "the pause between two tries of one download at setup; the line `trying again (n of N)` is shown through `progress.note` first",
    ("datasets/gaia_nearby.py", "esa_async"): "the poll of ESA's own job during the setup download; each poll is a request the download's progress line reports",
    ("util/log.py", "RelayReader.finish"): "a tenth of a second for the last lines of a child to arrive",
    ("queue/unitstuck.py", "request_dump"): "the starter waits a tenth of a second for a stack dump to land",
    ("admin/updateguard.py", "wait_until_idle"): OUTSIDE + " (update.sh)",
    ("web/jobs.py", "_end_orphan"): OUTSIDE,
    ("web/generate_page.py", "generate_job_stream.events"): OUTSIDE,
    ("web/lib/singleflight.py", "wait"): OUTSIDE,
    ("queue/wrapper.py", "_await"): OUTSIDE + " (the follower of a run)",
    ("queue/supervisor.py", "SystemClock.sleep"): OUTSIDE + " (the supervisor's own clock, between two turns)",
    ("queue/supervisor.py", "Supervisor.run"): OUTSIDE + " (the supervisor's loop, between two turns)",
    ("queue/units.py", "_end_stage"): "a pause (half a second, doubling, at most four tries) between two tries to queue the next stage's opening (OPS.186), the last unit's own work",
    ("queue/units.py", "_follow"): OUTSIDE + " (the follower of a run)",
    ("queue/runline.py", "start_when_free"): OUTSIDE + " (the waiting line)",
    ("queue/stopall.py", "kill_processes"): OUTSIDE + " (stop-all)",
    ("queue/api_jobs.py", "wait"): OUTSIDE + " (the API's wait for a job)",
    ("queue/warmstart.py", "drain"): OUTSIDE + " (the warm-up starter)",
    ("cli/worker.py", "worker_class.WarmWorker._wait_while_paused"): OUTSIDE + " (a worker that holds no unit while paused)",
    ("cli/warm.py", "_wait_for_update"): OUTSIDE,
    ("cli/warm.py", "_cancel_other_warm_up"): OUTSIDE,
    ("cli/warm.py", "_wait_for_generation"): OUTSIDE,
    ("cli/warm.py", "_settle_with_other_runs"): OUTSIDE,
    ("cli/reset_jobs.py", "reset_jobs"): OUTSIDE,
}


def _sleepers():
    """`{(file, function): lines}` of every call named `sleep` in the package."""
    found = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))

        def walk(node, stack):
            for child in ast.iter_child_nodes(node):
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    walk(child, stack + [child.name])
                    continue
                if isinstance(child, ast.Call):
                    function = child.func
                    name = function.attr if isinstance(function, ast.Attribute) else getattr(function, "id", None)
                    polls = name == "wait" and isinstance(function, ast.Attribute) and getattr(function.value, "attr", "") == "futures"      # concurrent.futures.wait(..., timeout=...) in a loop
                    if name == "sleep" or polls:
                        found.setdefault((path.relative_to(PACKAGE).as_posix(), ".".join(stack)), []).append(child.lineno)
                walk(child, stack)

        walk(tree, [])
    return found


def test_every_sleep_in_the_package_is_a_listed_wait():
    unlisted = {key: lines for key, lines in _sleepers().items() if key not in ALLOWED}
    assert not unlisted, ("A sleep in command code shows no activity, so the job looks stuck while it waits. Wait inside a blocking primitive that declares it (`unitbeat.waiting`), "
                          f"or, if the function does not run in a unit, list it in ALLOWED with the reason: {unlisted}")


def test_the_list_holds_no_entry_that_no_longer_sleeps():
    gone = sorted(set(ALLOWED) - set(_sleepers()))
    assert not gone, f"These entries of ALLOWED no longer match a sleep; delete them: {gone}"
