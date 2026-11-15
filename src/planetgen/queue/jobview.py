# planetgen/queue/jobview.py

"""
The jobs of a run and the history of runs, read from the ledger (OPS.224, L10 of `docs/design/operations/job-engine.md`).

The Queue page's live part stays what the beats and the run record say (a unit's counter, its memory, the verdict on it: they change every second and are Redis's). What
the engine has *decided and recorded* is in the control database's ledger (`engine_runs`, `engine_jobs`, `engine_attempts`, OPS.215), and this module is the one place
that turns it into what a page shows:

* `history`: the ledger's runs, newest first, filtered by state, kind or galaxy, each with its jobs counted by state.
* `jobs`: one run's jobs, filtered by state, stage and command, in the order they were planned, each with its newest attempt and, while it is not running, **why it waits**.
* `job`: one job with every attempt it has had (how each ended, what it cost, the error and the stack file it left).
* `why_waits`: the words for a job that is not running (reused: the reason `memorygate.decide` gave when a worker held the job at the gate, kept in `engine_jobs.waited`
  by `attempts.note_wait`; the retry hold of `not_before`; the jobs it needs).

Every function takes the control connection and does plain reads; no page touches a Redis key to know what a job is.
"""

import json
from datetime import timezone

from planetgen.queue import ledger, planrows, runrecord
from planetgen.util.format import format_duration_seconds

DEFAULT_LIMIT, MAX_LIMIT = 50, 200
"""int: The jobs or runs a page of the lists holds by default, and at most."""

ATTEMPT_WORDS = {
    "running": "running", "done": "finished", "failed": "failed", "killed_memory": "killed for memory", "timed_out": "ran out of time",
    "stuck": "stopped as stuck", "replaced": "replaced by a later attempt", "lost": "lost: its process was gone",
}
"""dict: An attempt's outcome in words."""

NEED_WORDS = {"units": "every unit of its stage"}
"""dict: What a `needs` entry that is not a command's name stands for."""


def _epoch(value):
    """A database time as Unix seconds, or `None`."""
    return None if value is None else value.replace(tzinfo=timezone.utc).timestamp()


def _json(text):
    """The JSON a TEXT column holds, or `None` for none and for text that is not JSON."""
    if not text:
        return None
    try:
        return json.loads(text)
    except ValueError:
        return None


def engine_run(conn, run_id):
    """The `engine_runs` id the run record names, or `None` for a run that has not been cut yet."""
    found = conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)
    return int(runrecord._s(found)) if found else None


def _clamp(limit, offset):
    return max(1, min(int(limit), MAX_LIMIT)), max(0, int(offset))


def counts(control, run):
    """`{state: jobs}` of the run's jobs, every state of `ledger.JOB_STATES` present."""
    found = {state: 0 for state in ledger.JOB_STATES}
    for row in control.execute("SELECT state, COUNT(*) AS n FROM engine_jobs WHERE run_id = ? GROUP BY state", (run,)).fetchall():
        found[row["state"]] = row["n"]
    return found


def history(control, state=None, kind=None, galaxy=None, limit=DEFAULT_LIMIT, offset=0):
    """
    The ledger's runs, newest first: `{"total", "items": [{"id", "kind", "galaxy", "state", "held_reason", "created", "started", "finished", "jobs": {state: n}, "attempts"}], "facets": {"state": {value: runs},
    "kind": {...}}}` (the facets count every run, whatever the filters).
    `attempts` is how many attempts the run's jobs have had in all (more than its jobs when some were retried).
    """
    limit, offset = _clamp(limit, offset)
    where, args = [], []
    for column, value in (("state", state), ("kind", kind), ("galaxy", galaxy)):
        if value:
            where.append(f"{column} = ?")
            args.append(value)
    clause = f" WHERE {' AND '.join(where)}" if where else ""
    total = control.execute(f"SELECT COUNT(*) AS n FROM engine_runs{clause}", tuple(args)).fetchone()["n"]
    rows = control.execute(f"SELECT id, kind, galaxy, state, held_reason, created_at, started_at, finished_at FROM engine_runs{clause} ORDER BY id DESC LIMIT ? OFFSET ?",
                           (*args, limit, offset)).fetchall()
    items = []
    for row in rows:
        by_state = counts(control, row["id"])
        spent = control.execute("SELECT COALESCE(SUM(attempt), 0) AS n FROM engine_jobs WHERE run_id = ?", (row["id"],)).fetchone()["n"]
        items.append({"id": row["id"], "kind": row["kind"], "galaxy": row["galaxy"], "state": row["state"], "held_reason": row["held_reason"],
                      "created": _epoch(row["created_at"]), "started": _epoch(row["started_at"]), "finished": _epoch(row["finished_at"]),
                      "jobs": by_state, "attempts": int(spent)})
    facets = {column: {row["value"]: row["n"] for row in control.execute(f"SELECT {column} AS value, COUNT(*) AS n FROM engine_runs WHERE {column} IS NOT NULL GROUP BY {column}").fetchall()}
              for column in ("state", "kind")}
    return {"total": total, "items": items, "facets": facets}


def _needs_text(control, run, stage, needs):
    """What a waiting job's `needs` still lack, in words, or `None` when everything it needs is done."""
    lacking = []
    for need in needs or []:
        if str(need).startswith(planrows.AFTER):
            before = str(need)[len(planrows.AFTER):]
            row = control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, before, planrows.RELEASE)).fetchone()
            if row is None or row["state"] != "done":
                lacking.append(f"stage `{before}` to end")
        elif need == "units":
            left = control.execute("SELECT COUNT(*) AS n FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key <> '' AND state <> 'done'", (run, stage)).fetchone()["n"]
            if left:
                total = control.execute("SELECT COUNT(*) AS n FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key <> ''", (run, stage)).fetchone()["n"]
                lacking.append(f"{NEED_WORDS['units']} ({total - left:,} of {total:,} done)")
        else:
            row = control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, stage, need)).fetchone()
            if row is None or row["state"] != "done":
                lacking.append(f"`{need}`")
    return "waits for " + " and ".join(lacking) if lacking else None


def why_waits(control, row, now):
    """
    Why the job is not running, in words, or `""` for a job that is running or over (`row` as `jobs` reads it). In this order: a retry hold (`not_before` in the future), the
    last reason a worker held it (`waited`: the memory gate's own words), a job it still needs, and otherwise that it is queued and no worker has taken it.
    """
    if row["state"] in ("running", "done", "failed", "canceled"):
        return ""
    held = _epoch(row["not_before"])
    if held is not None and held > now:
        last = control.execute("SELECT outcome FROM engine_attempts WHERE job_id = ? AND attempt = ?", (row["id"], row["attempt"])).fetchone()
        after = f" after its attempt {row['attempt']} {ATTEMPT_WORDS.get(last['outcome'], last['outcome'])}" if last else ""
        return f"held back for a retry for another {format_duration_seconds(held - now)}{after}"
    if row["waited"]:
        return row["waited"]
    needed = _needs_text(control, row["run_id"], row["stage"], _json(row["needs"]))
    return needed or "queued: no worker has taken it yet"


def _job_view(control, row, now):
    return {"id": row["id"], "stage": row["stage"], "command": row["command"], "item": row["item_key"], "state": row["state"], "attempts": row["attempt"],
            "class": row["job_class"], "claim_mb": row["claim_mb"], "waits": why_waits(control, row, now),
            "last": None if row["attempt"] == 0 or row.get("outcome") is None else
                    {"outcome": row["outcome"], "words": ATTEMPT_WORDS.get(row["outcome"], row["outcome"]), "error_kind": row.get("error_kind"),
                     "peak_mb": row.get("peak_mb"), "cpu_seconds": row.get("cpu_seconds"), "started": _epoch(row.get("started_at")), "ended": _epoch(row.get("ended_at"))},
            "created": _epoch(row["created_at"]), "updated": _epoch(row["updated_at"])}


def jobs(control, run, now, state=None, stage=None, command=None, limit=DEFAULT_LIMIT, offset=0):
    """
    One run's jobs in the order they were planned: `{"engine_run", "counts": {state: n}, "total", "items": [{"id", "stage", "command", "item", "state", "attempts",
    "class", "claim_mb", "waits", "last", "created", "updated"}]}`. `last` is the newest attempt (`outcome`, `words`, `error_kind`, `peak_mb`, `cpu_seconds`, `started`, `ended`);
    `total` counts the jobs that pass the filters.
    """
    limit, offset = _clamp(limit, offset)
    where, args = ["j.run_id = ?"], [run]
    for column, value in (("state", state), ("stage", stage), ("command", command)):
        if value:
            where.append(f"j.{column} = ?")
            args.append(value)
    clause = " AND ".join(where)
    total = control.execute(f"SELECT COUNT(*) AS n FROM engine_jobs j WHERE {clause}", tuple(args)).fetchone()["n"]
    rows = control.execute(
        "SELECT j.id, j.run_id, j.stage, j.command, j.item_key, j.needs, j.job_class, j.claim_mb, j.state, j.attempt, j.not_before, j.waited, j.created_at, j.updated_at,"
        " a.outcome, a.error_kind, a.peak_mb, a.cpu_seconds, a.started_at, a.ended_at"
        f" FROM engine_jobs j LEFT JOIN engine_attempts a ON a.job_id = j.id AND a.attempt = j.attempt WHERE {clause} ORDER BY j.id LIMIT ? OFFSET ?",
        (*args, limit, offset)).fetchall()
    return {"engine_run": run, "counts": counts(control, run), "total": total, "items": [_job_view(control, row, now) for row in rows]}


def job(control, run, job_id, now):
    """
    One job of the run with every attempt it has had, oldest first, or `None` when the run has no such job: the job as `jobs` gives it, `arguments` (what it was told to do) and
    `attempt_rows`: `{"attempt", "pid", "outcome", "words", "started", "ended", "peak_mb", "cpu_seconds", "error_kind", "error_text", "stack_file"}`.
    """
    row = control.execute(
        "SELECT j.id, j.run_id, j.stage, j.command, j.item_key, j.needs, j.arguments, j.job_class, j.claim_mb, j.state, j.attempt, j.not_before, j.waited, j.created_at, j.updated_at,"
        " a.outcome, a.error_kind, a.peak_mb, a.cpu_seconds, a.started_at, a.ended_at"
        " FROM engine_jobs j LEFT JOIN engine_attempts a ON a.job_id = j.id AND a.attempt = j.attempt WHERE j.run_id = ? AND j.id = ?", (run, job_id)).fetchone()
    if row is None:
        return None
    attempts = control.execute("SELECT attempt, pid, outcome, started_at, ended_at, peak_mb, cpu_seconds, error_kind, error_text, stack_file"
                               " FROM engine_attempts WHERE job_id = ? ORDER BY attempt", (job_id,)).fetchall()
    return {**_job_view(control, row, now), "arguments": _json(row["arguments"]),
            "attempt_rows": [{"attempt": item["attempt"], "pid": item["pid"], "outcome": item["outcome"], "words": ATTEMPT_WORDS.get(item["outcome"], item["outcome"]),
                              "started": _epoch(item["started_at"]), "ended": _epoch(item["ended_at"]), "peak_mb": item["peak_mb"], "cpu_seconds": item["cpu_seconds"],
                              "error_kind": item["error_kind"], "error_text": item["error_text"], "stack_file": item["stack_file"]} for item in attempts]}
