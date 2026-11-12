# planetgen/api/workqueue.py

"""
Admin endpoints behind the work queue page (ADM.10, `web/queue_page.py`), reading and controlling the runs in Redis
(PERF.99, `planetgen/queue/runrecord.py`): the queue's state (workers active, the server's load, the run that owns the
machine, the whole-queue pause), the runs with their stages, and the controls on them (cancel, clear a run that stopped
answering, delete a finished one). Every route needs a logged-in admin past the forced credential change; every change
is written to the audit and activity logs (`audit`).

Rerunning a run (Retry) starts a Generate page job, so it lives with the page (`web/queue_page.py`), not here.
"""

import re
import time

from flask import Blueprint, g, jsonify, request

from planetgen.queue import jobbeat, load as systemLoad, pause as queuePause, profileclaim, redisqueue, requires, runcontrol, runrecord, runview, stopall, unitwatch, warmphases, work as workQueue, workers as queueWorkers

from .authz import audit, require_admin
from .common import ApiError, require_json_body
from .routes import _paginate
from .schemas import QueueAction, RunAction, TuningAction, WorkerCount, parse_body

bp = Blueprint("workqueue", __name__, url_prefix="/api/admin/work")

EVENTS_SHOWN = 300
"""int: How many of a run's last events its detail carries."""


def _connection():
    """A Redis client that answers, or `None`."""
    try:
        connection = redisqueue.connect()
        connection.ping()
        return connection
    except Exception:  # noqa: BLE001 -- any connection failure
        return None


def _run_id(value):
    if not isinstance(value, str) or not runrecord.RUN_ID_RE.match(value):
        raise ApiError("No such run.", 404)
    return value


def _status(connection, now):
    """The queue's state for the page: the pause, the run that owns the machine and whether it has gone quiet, the
    workers and the runs live now."""
    paused = queuePause.state(connection)
    holder = runrecord.active_run(connection)
    held = runrecord.read_run(connection, holder) if holder else None
    live = [run_id for kind in runrecord.KINDS for run_id in [runrecord.active_run(connection, kind)] if run_id]
    heard = None if held is None else max(now - (held["updated"] or held["started"] or now), 0.0)
    return {
        "paused": paused["paused"], "paused_by": paused["by"], "paused_at": paused["at"], "paused_reason": paused["reason"],
        "solo_unit": paused["solo"]["unit"] if paused["solo"] else None,
        "holder": holder, "holder_title": held["title"] if held else None, "heartbeat_age_s": heard,
        "stale": runrecord.is_stale(held, now) if held else False,
        "workers_active": sum(queueWorkers.live_worker_count(redisqueue.queue(name, connection))
                              for name in (redisqueue.BULK_QUEUE, redisqueue.MAINTENANCE_QUEUE)),
        "runs_active": len(live),
    }


@bp.route("")
@require_admin(fresh=True)
def work():
    """
    `GET /api/admin/work[?limit=&offset=]` -- the queue and one page of runs, newest first: `{"available", "status":
    {...`_status`, "load": systemLoad.load_average()}, "items": [`runview.summary`], "total", "limit", "offset"}`.
    `available` is false (and the rest empty) when no Redis server answers.
    """
    limit, offset = _paginate(request.args)
    load = systemLoad.load_average()
    connection = _connection()
    if connection is None:
        return jsonify({"available": False, "status": {"load": load}, "stopped": stopall.state(), "items": [],
                        "total": 0, "limit": limit, "offset": offset})
    now = time.time()
    status = _status(connection, now)
    status["load"] = load
    ids = runrecord.recent_runs(connection, limit=10_000)
    items = []
    for run_id in ids[offset:offset + limit]:
        run = runrecord.read_run(connection, run_id)
        if run is not None:
            items.append(runview.summary(run, now))
    return jsonify({"available": True, "status": status, "stopped": stopall.state(), "items": items,
                    "total": len(ids), "limit": limit, "offset": offset})


DATABASE_RE = re.compile(r"^[A-Za-z0-9_$]{1,64}$")
"""re.Pattern: A database name as `GET /api/admin/work/warm` takes it."""


@bp.route("/warm")
@require_admin()
def warm():
    """
    `GET /api/admin/work/warm?database=NAME` -- what the cache warm-up has done on that database, phase by phase
    (`warmphases.read`, PERF.141): `{"available", "run", "phases": [{"name", "text", "families", "state",
    "state_text", "line", "built", "cached", "left_out", "seconds", "started_at", "finished_at"}], "auto":
    {"enabled", "wanted", "holder"}}`.
    """
    database = request.args.get("database") or ""
    if not DATABASE_RE.match(database):
        raise ApiError("`database` must name a database.", 400)
    from planetgen.web import warmup    # the phases are the web layer's (`planetgen warm` builds them)
    return jsonify(warmphases.read(database, warmup.PHASES))


@bp.route("/<run_id>")
@require_admin(fresh=True)
def run_detail(run_id):
    """
    `GET /api/admin/work/<id>` -- one run with its stages and the last events of its log:
    `{"run": `runview.detail`, "events": [{"id", "t", "level", "kind", "message"}], "log": bool}`; an event also holds the
    `stage` (from 0), `stages`, `unit` and `seconds` it was logged with, when it has them, which is what
    `runrecord.event_record` makes a log record of, so a page draws the same line the download holds (UX.123). `log` is
    whether the log can be downloaded (`/log`).
    """
    run_id = _run_id(run_id)
    connection = _connection()
    run = runrecord.read_run(connection, run_id) if connection is not None else None
    if run is None:
        raise ApiError("No such run.", 404)
    events = runrecord.read_events(connection, run_id, count=100_000)[-EVENTS_SHOWN:]
    return jsonify({"run": runview.detail(run, time.time()), "log": True,
                    "events": [{"id": entry_id, "t": float(fields.get("t", 0)), "level": fields.get("level", "info"),
                                "kind": fields.get("kind", "log"), "message": fields.get("message", ""),
                                **{name: fields[name] for name in runrecord.EVENT_FIELD_NAMES if name in fields}}
                               for entry_id, fields in events]})


@bp.route("/jobs", methods=["GET"])
@require_admin()
def running_jobs():
    """
    `GET /api/admin/work/jobs` -- the jobs that are not units of a run and are running, from their beats (OPS.142): `{"jobs": [{"kind" ("api", "web" or "tile-prefetch"), "job", "label",
    "state" ("working" or "dead"), "reason", "phase", "running_seconds", "beat_age_seconds"}]}`, the longest running first. A job that died stays listed, as `dead`, for an hour.
    """
    connection = _connection()
    if connection is None:
        return jsonify({"jobs": [], "available": False})
    return jsonify({"jobs": jobbeat.states(connection), "available": True})


@bp.route("/<run_id>/units")
@require_admin()
def run_units(run_id):
    """
    `GET /api/admin/work/<id>/units` -- what each running unit of the run is doing (OPS.138): `{"units": [{"stage", "unit", "state", "severity", "reason", "actionable", "agreed",
    "counter", "rate", "projected_finish_seconds", "still_seconds", "over_limit", "phase", "waiting_on", "done", "estimate", "restarts", "resident_mb", "beat_age_seconds",
    "running_seconds", "progress"}]}`, as the process that follows the run last judged them (`unitwatch`, the judge of `beatjudge`). `progress` is the unit's own bar (PERF.147, `taskprogress.as_dict`: `show` once the unit is expected or has been running over 15 s, `percent`,
    `eta_seconds`, `source` stored, blended or live). `state` is `working`, `slow`, `waiting`, `quiet`, `looping` or `dead`; only the stage being followed has units, and a run that is over has none.
    """
    run_id = _run_id(run_id)
    connection = _connection()
    run = runrecord.read_run(connection, run_id) if connection is not None else None
    if run is None:
        raise ApiError("No such run.", 404)
    units = [{"stage": index, "unit": unit, **shown} for index in range(len(run["stages"])) for unit, shown in sorted(unitwatch.read(connection, run_id, index).items())]
    return jsonify({"units": units})


@bp.route("/<run_id>/log")
@require_admin(fresh=True)
def run_log(run_id):
    """`GET /api/admin/work/<id>/log` -- the run's whole log as `{"text": ...}`, one line to an event: from Redis while its
    record lasts, else the file kept in the jobs directory (the last `runcontrol.KEEP_LOGS`)."""
    run_id = _run_id(run_id)
    connection = _connection()
    if connection is not None and runrecord.read_run(connection, run_id) is not None:
        text = runrecord.log_text(connection, run_id)
    else:
        path = runcontrol.log_path(run_id)
        try:
            with open(path, "r", encoding="utf-8") as handle:
                text = handle.read()
        except (OSError, TypeError):
            raise ApiError("No such run.", 404) from None
    return jsonify({"text": text})


@bp.route("/<run_id>/control", methods=["POST"])
@require_admin(fresh=True)
def control(run_id):
    """
    `POST /api/admin/work/<id>/control` `{"action": "cancel"|"clear"}`. Cancel stops a live run now (its waiting units are
    taken off the queue, its started ones sent the stop command, the machine is freed); clear ends a run that has stopped
    answering (`runcontrol.is_clearable`) and frees the machine. Returns `{"ok": bool, "canceled", "stopped"}`: false
    when the run is over or, for clear, still answering.
    """
    action = parse_body(RunAction, require_json_body()).action
    run_id = _run_id(run_id)
    connection = _connection()
    if connection is None:
        raise ApiError("Redis is not answering.", 503)
    run = runrecord.read_run(connection, run_id)
    if run is None:
        raise ApiError("No such run.", 404)
    who = g.admin_user["username"]
    stopped = (runcontrol.cancel_run if action == "cancel" else runcontrol.clear_run)(connection, run_id, who=who)
    if stopped is not None:
        audit(f"work.{action}", target=f"run:{run_id}", detail=run["title"][:200])
    return jsonify({"ok": stopped is not None, **(stopped or {})})


@bp.route("/<run_id>/delete", methods=["POST"])
@require_admin(fresh=True)
def delete(run_id):
    """`POST /api/admin/work/<id>/delete` deletes a finished run's record. `{"ok": bool}`: false for a live run."""
    run_id = _run_id(run_id)
    connection = _connection()
    if connection is None:
        raise ApiError("Redis is not answering.", 503)
    ok = runrecord.delete_run(connection, run_id)
    if ok:
        audit("work.delete", target=f"run:{run_id}")
    return jsonify({"ok": ok})


@bp.route("/clear-finished", methods=["POST"])
@require_admin(fresh=True)
def clear_finished():
    """`POST /api/admin/work/clear-finished` deletes every finished run's record (ADM.41). `{"cleared": <count>}`;
    live runs stay."""
    connection = _connection()
    if connection is None:
        raise ApiError("Redis is not answering.", 503)
    cleared = runcontrol.clear_finished(connection)
    if cleared:
        audit("work.clear", detail=f"{cleared} finished runs")
    return jsonify({"cleared": cleared})


@bp.route("/queue", methods=["POST"])
@require_admin(fresh=True)
def queue():
    """
    `POST /api/admin/work/queue` `{"action": "pause"|"resume"|"stop-all"|"resume-all"}` pauses
    the whole queue (the workers of the generation and maintenance lanes take no job
    until it's resumed; running tasks finish; PERF.152) or resumes it. `{"ok": bool}`: false when it
    already was.

    `stop-all` (PERF.186, `stopall.stop_all`) stops everything: the hold, every run and job canceled, every queue
    emptied, every planetgen worker and command process on the host ended, and nothing starts until `resume-all`
    lifts the hold. It works with no Redis answering. `{"ok": true, "summary": {"runs", "queued", "job", "killed",
    "survivors", "already"}}`; `resume-all` returns `{"ok": whether it was stopped}`.
    """
    action = parse_body(QueueAction, require_json_body()).action
    if action == "stop-all":
        summary = stopall.stop_all(g.admin_user["username"], connection=_connection())
        audit("work.queue.stop-all", target="queue",
              detail=f"{summary['runs']} runs, {summary['queued']} queued, {summary['killed']} processes")
        return jsonify({"ok": True, "summary": summary})
    if action == "resume-all":
        ok = stopall.resume_all(_connection())
        if ok:
            audit("work.queue.resume-all", target="queue")
        return jsonify({"ok": ok})
    connection = _redis()
    if action == "pause":
        ok = queuePause.pause(g.admin_user["username"], connection)
    else:
        ok = queuePause.resume(connection)
    if ok:
        audit(f"work.queue.{action}", target="queue")
    return jsonify({"ok": ok})


def _workers_state(connection):
    queue = redisqueue.queue(redisqueue.BULK_QUEUE, connection)
    return {"wanted": queueWorkers.wanted(connection), "live": queueWorkers.live_worker_count(queue),
            "ludicrous": queueWorkers.ludicrous_count(), "cpus": workQueue.cpu_count(),
            "host_cap": queueWorkers.host_worker_cap(), "host_running": queueWorkers.host_worker_processes(),
            "host_refused": int(connection.get(queueWorkers.REFUSED_KEY) or 0)}


def _redis():
    try:
        connection = redisqueue.connect()
        connection.ping()
    except Exception as exc:  # noqa: BLE001 -- any connection failure
        raise ApiError(f"Redis is not answering: {exc}", 503) from exc
    return connection


@bp.route("/workers", methods=["GET"])
@require_admin()
def workers_state():
    """`GET /api/admin/work/workers` -- `wanted` (the bulk worker count), `live` (bulk workers running now),
    `ludicrous` (the count "Ludicrous Speed" sets) and `cpus` (ADM.15)."""
    return jsonify(_workers_state(_redis()))


@bp.route("/held", methods=["GET"])
@require_admin()
def held_jobs():
    """`GET /api/admin/work/held` -- the jobs held for a job they require (PERF.247, PERF.248): `{"total", "items": [{"id", "queue",
    "label", "waits_for", "run", "stage"}]}`, the first 100; a job that requires nothing is never here."""
    return jsonify(requires.held_overview(_redis()))


@bp.route("/workers", methods=["POST"])
@require_admin(fresh=True)
def set_workers():
    """
    `POST /api/admin/work/workers` `{"count": n}` or `{"ludicrous": true}` sets the bulk worker count
    (ADM.15): running work picks it up as workers start and exit (lowering it asks the extra workers to finish
    their unit and exit). Returns the new state; written to the audit log.
    """
    body = parse_body(WorkerCount, require_json_body())
    connection = _redis()
    count = queueWorkers.ludicrous_count() if body.ludicrous else body.count
    queueWorkers.set_wanted(connection, count)
    audit("work.workers.set", target="workers", detail=f"count={count}" + (" ludicrous" if body.ludicrous else ""))
    return jsonify(_workers_state(connection))


@bp.route("/tuning", methods=["GET"])
@require_admin()
def tuning():
    """`GET /api/admin/work/tuning` -- what the queue has learned per job kind (OPS.158): `self_tuning` (whether the job profile is read), `kinds` (each `kind`, `typed_mb`, `figure_mb`
    the claim from the profile, `pinned_mb`, `samples`, `last_sample`, `changed` the newest change made to it) and `changes` (the newest 20 changes the engine or an admin made)."""
    return jsonify(profileclaim.overview(_redis()))


@bp.route("/tuning", methods=["POST"])
@require_admin(fresh=True)
def set_tuning():
    """
    `POST /api/admin/work/tuning` `{"action": "pin", "kind": "unit:x", "mb": 1500}` fixes the claim of a job kind, `"unpin"` lets the profile decide again and `"clear"` forgets everything
    learned about the kind (OPS.158). Logged with the other changes of the engine and written to the audit log. Returns the new overview.
    """
    body = parse_body(TuningAction, require_json_body())
    connection = _redis()
    who = g.admin_user["username"]
    if body.action == "pin":
        profileclaim.pin(connection, body.kind, body.mb, who)
    elif body.action == "unpin":
        if not profileclaim.unpin(connection, body.kind, who):
            raise ApiError("That kind has no pin.", 404)
    else:
        profileclaim.clear(connection, body.kind, who)
    audit(f"work.tuning.{body.action}", target=body.kind, detail=f"mb={body.mb}" if body.mb is not None else "")
    return jsonify(profileclaim.overview(connection))
