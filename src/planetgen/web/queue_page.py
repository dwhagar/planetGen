# planetgen/web/queue_page.py

"""
The admin work queue page (ADM.10, PERF.99): `/admin/queue` lists the runs in Redis (`planetgen/queue/runrecord.py`)
newest first, under the number of workers active and the server's load ("x / x / x" over 1, 5 and 15 minutes), the run
that owns the machine and the whole-queue pause; `/admin/queue/<id>` shows one run, its stages with their counters and
the last lines of its log. Every control (cancel, clear a run that stopped answering, retry, delete, pause the whole
queue) first shows a confirmation page (`/admin/queue/confirm/<action>`), then posts to `/admin/queue/action`, and is
written to the activity log.

Reads and controls go through the API (`api/workqueue.py`); Retry starts a Generate page job (`web/jobs.py`), like the
Generate page's own forms.
"""

import time

from flask import Response, abort, current_app, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import Column, Result, Table, time_text
from planetgen.queue import progress_rate, runline, runrecord, warmphases
from planetgen.admin import activity_log
from planetgen.util import log, logrecord
from planetgen.util.format import counted

from . import bp, jobs, tables
from .admin_pages import came_from, format_bytes, safe_next, _api_message, _cookie_header, _flash, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb, db_name

KIND_LABELS = {"generation": "Generation", "maintenance": "Maintenance", "check": "Database check"}

STATE_LABELS = {
    "running": "Running",
    "cancelling": "Cancelling",
    "done": "Done",
    "failed": "Failed",
    "cancelled": "Cancelled",
}

STAGE_STATE_LABELS = {"waiting": "Waiting", "running": "Running", "done": "Done"}

ACTION_TEXT = {
    "cancel": ("Cancel", "Its running units are stopped and nothing more is started; the units that finished keep "
                         "what they saved. The machine is free for the next run at once."),
    "clear": ("Clear", "It has stopped answering: its process is gone or it has not been heard from for too long. It is "
                       "ended as failed and the machine is freed so the next run can start."),
    "retry": ("Retry", None),
    "delete": ("Delete", "The run's record and log are deleted. Nothing it generated is touched."),
    "fill-habitability": ("Regenerate habitability",
                          "Re-scores the habitability of every planet and moon from the conditions already stored: the "
                          "PHI-6 index, the six domains with their tiers, the survival equipment levels, the note and "
                          "the double-star flux swing. It does not re-draw a planet, so after a change that alters how "
                          "planets are drawn a regeneration is needed instead. It runs as a queue run you can follow "
                          "here, stop and start again; the result is the same however often it runs."),
    "fill-map-luminosity": ("Fill map luminosity",
                            "Fills the map luminosity of the scattered phenomena that do not have it yet, from their "
                            "stored masses, so the Galaxy Map's glow includes them. It only fills rows still empty. It "
                            "runs as a queue run you can follow here, stop and start again."),
    "clear-finished": ("Clear finished runs", "Every finished run's record is deleted, leaving the ones still running. "
                                              "Nothing they generated is touched."),
    "queue-pause": ("Pause the queue", "No run starts or takes another unit until you resume the queue. "
                                       "Running units finish first."),
    "queue-resume": ("Resume the queue", "Runs waiting for the queue start again, one at a time."),
    "queue-stop-all": ("Stop everything", "Cancels every run and the Generate page's job, empties every queue, and ends "
                                          "every planetgen worker and command process on this server with the processes "
                                          "they started. Nothing starts again, not even from the map, until you press "
                                          "Resume. What the runs already saved stays; the runs must be started again."),
    "queue-resume-all": ("Resume", "Lifts the stop: work can start again. Nothing is started by this."),
}


QUEUE_WIDE_ACTIONS = ("clear-finished",)
"""tuple[str]: The controls that act on the whole queue, not one run."""

FILL_COMMANDS = {"fill-habitability": "rescore-habitability", "fill-map-luminosity": "fill-map-luminosity"}
"""dict: The Queue page's start buttons (ADM.129) and the `planetgen` command each starts as a job, with the rules of
any generation run: refused while another run is live, one run lock, the audit entry, Stop and Retry."""


def _is_site_wide(action):
    """Whether `action` acts on no particular run (so its confirmation page takes no `?run=`)."""
    return action.startswith("queue-") or action in QUEUE_WIDE_ACTIONS or action in FILL_COMMANDS

RETRYABLE = ("plan", "scatter", "clear-scatter", "fill-map-luminosity", "rescore-habitability", "clean-nebula-field",
             "galaxy", "sectors", "link", "backfill", "settle", "sector", "population", "territories")
"""tuple[str]: The `planetgen` commands whose run Retry starts again (the command line the run recorded)."""


def _duration(seconds):
    """`"1 h 02 m"`, `"3 m 05 s"`, `"12.4 s"`, or `""`."""
    if seconds is None:
        return ""
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.1f} s"
    whole = int(seconds)
    hours, rest = divmod(whole, 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return f"{hours} h {minutes:02d} m"
    return f"{minutes} m {secs:02d} s"


def _eta_text(eta, live):
    """The time left as a range around the estimate, "estimating" while a running run has none (PERF.33, never dashes),
    `""` for a run that isn't running."""
    if not live:
        return ""
    if eta is None or eta.get("seconds") is None:
        return "estimating"
    if eta.get("stalled"):
        return "waiting on progress"
    low, high = progress_rate.eta_range(eta["seconds"])
    if eta.get("rough") and eta.get("low") is not None:       # PERF.184: a guess from the host, not from timings
        low, high = eta["low"], eta["high"]
    low_text, high_text = _duration(low), _duration(high)
    text = high_text if low_text == high_text else f"{low_text} to {high_text}"
    return f"{text} (not timed yet)" if eta.get("rough") else text


def _actions(run):
    """The controls an admin can use on `run` now."""
    if run["live"]:
        return ["cancel", "clear"] if run["stale"] else ["cancel"]
    actions = []
    if run["state"] in ("failed", "cancelled") and _retry_plan(run) is not None:
        actions.append("retry")
    actions.append("delete")
    return actions


def _progress_text(run):
    stage = run["stage"]
    if stage is None:
        return ""
    where = f"stage {stage['index'] + 1} of {stage['count']}: {stage['label']}" if run["live"] else \
        f"{stage['count']} stage{'s' if stage['count'] != 1 else ''}"
    units = f", {run['units_done']} of {run['units']} {counted('unit', run['units'])}" if run["units"] else ""
    return where + units


def run_view(run):
    """A run (the API's `runview.summary`) as the pages show it."""
    return {
        "id": run["id"],
        "title": run["title"],
        "kind": KIND_LABELS.get(run["kind"], run["kind"]),
        "who": run["who"],
        "state": run["state"],
        "state_label": STATE_LABELS.get(run["state"], run["state"]),
        "error": run["error"],
        "started_at": _iso(run["started"]),
        "finished_at": _iso(run["finished"]),
        "duration": _duration(run["seconds"]),
        "progress": _progress_text(run),
        "percent": run["percent"],
        "eta": _eta_text(run["eta"], run["live"]),
        "stale": run["stale"],
        "live": run["live"],
        "job_url": url_for("web.generate_job", job_id=run["job"]) if run["job"] and jobs.JOB_ID_RE.match(run["job"]) else None,
        "url": url_for("web.admin_queue_run", run_id=run["id"]),
        "log_url": url_for("web.admin_queue_log", run_id=run["id"]),
        "actions": [(action, ACTION_TEXT[action][0]) for action in _actions(run)],
    }


def _retry_plan(run):
    """
    What Retry on `run` runs again, as `{"kind", "title", "steps", "database", "describe"}`, or `None`: for a run a
    Generate page job started, that job's steps from the one that didn't finish; for a run started anywhere else, the
    same `planetgen` command line. Either way the whole command is retried, not the one run: generation skips what is
    already there, so only the missing work is done again.
    """
    if run["kind"] != "generation":
        return None
    if run.get("job") and jobs.JOB_ID_RE.match(run["job"]):
        try:
            found = jobs.remaining_steps(run["job"])
        except OSError:
            found = None
        if found is not None:
            job, steps = found
            return {"kind": job.get("kind") or "galaxy", "title": f"Retry: {job.get('title') or run['title']}",
                    "steps": steps, "database": job.get("database"),
                    "describe": ", then ".join(step["label"] for step in steps)}
    argv = [str(arg) for arg in run.get("argv") or []]
    if argv and argv[0] in RETRYABLE:
        label = " ".join(["planetgen", *argv])
        return {"kind": argv[0], "title": f"Retry: {label}"[:200],
                "steps": [{"label": label, "argv": [jobs.python_executable(), *jobs.GENERATE_COMMAND, *argv]}],
                "database": run.get("database"), "describe": label}
    return None


def _stopped_view(stopped):
    """The stop-all hold for the page: `stopped`, who and when (ISO)."""
    stopped = stopped or {}
    return {"stopped": bool(stopped.get("stopped")), "by": stopped.get("by"), "at": _iso(stopped.get("at"))}


def _status_view(status):
    load = status.get("load") or {}
    return {
        "workers_active": status.get("workers_active", 0),
        "runs_active": status.get("runs_active", 0),
        "load": load.get("text") or "unknown",
        "load_label": "Load average, 1 / 5 / 15 minutes",
        "paused": status.get("paused", False),
        "paused_by": status.get("paused_by"),
        "paused_at": _iso(status.get("paused_at")),
        "holder": status.get("holder"),
        "holder_title": status.get("holder_title"),
        "heartbeat": _duration(status.get("heartbeat_age_s")),
        "stale": status.get("stale", False),
        "holder_url": url_for("web.admin_queue_run", run_id=status["holder"]) if status.get("holder") else None,
    }


def _load(cookie_header, run_id):
    """`(detail, run)` from the API -- the whole answer and the run in it -- or `(None, None)`."""
    try:
        body = apiclient.admin_work_run(cookie_header, run_id)
    except apiclient.NotFoundError:
        return None, None
    return body, body["run"]


def _job_cells(run):
    """One run as the Runs table's cells."""
    view = run_view(run)
    links = []
    for action, label in view["actions"]:
        links += [" "] if links else []
        links.append({"text": label, "href": url_for("web.admin_queue_confirm", action=action, run=view["id"])})
    if view["job_url"]:
        links += [" "] if links else []
        links.append({"text": "Job", "href": view["job_url"]})
    started = time_text(view["started_at"]) or "not yet"
    stale = " (not answering)" if view["stale"] else ""
    return [
        {"text": view["title"], "href": view["url"]},
        {"text": view["state_label"] + stale},
        {"text": started},
        {"text": view["duration"] or ""},
        {"text": view["eta"]},
        {"text": view["progress"]},
        {"text": "", "parts": links} if links else {"text": ""},
    ]


def _jobs_load(state, limit, offset, want_facets):
    """The runs, newest first, `limit` from `offset`."""
    _identity, bounce = _require_admin()
    if bounce is not None:
        abort(403)
    body = apiclient.admin_work(_cookie_header(), limit=limit, offset=offset)
    return Result([_job_cells(run) for run in body["items"]], body["total"], None)


JOBS_TABLE = tables.register(Table(
    "queue-jobs", "Runs",
    [Column("job", "Run", sortable=False), Column("status", "Status", sortable=False), Column("started", "Started", sortable=False),
     Column("duration", "Duration", sortable=False), Column("left", "Time left", sortable=False),
     Column("progress", "Progress", sortable=False),
     Column("actions", "Actions", sortable=False)],
    _jobs_load, prefix="jobs_", noun=("run", "runs"),
))


def _warm_view(warm):
    """The cache warm-up panel (PERF.141): its phases, the run they belong to and what the automatic start waits for, or `None`."""
    if warm is None or not warm["available"]:
        return None
    run = warm["run"] or {}
    auto = warm["auto"]
    return {
        "run_state": run.get("state"),
        "run_text": warmphases.STATE_TEXT.get(run.get("state"), "No run recorded"),
        "run_started": time_text(_iso(run.get("started_at"))) or "",
        "run_finished": time_text(_iso(run.get("finished_at"))) or "",
        "run_reason": run.get("reason") or "",
        "run_error": run.get("error") or "",
        "phases": [{**phase, "started": time_text(_iso(phase["started_at"])), "finished": time_text(_iso(phase["finished_at"])),
                    "number": number} for number, phase in enumerate(warm["phases"], 1)],
        "auto_enabled": auto["enabled"], "auto_wanted": auto["wanted"], "auto_running": auto["holder"] is not None,
        "index_size": format_bytes(warm["index"]["size_bytes"]) if warm["index"]["size_bytes"] else "",
        "index_kinds": [{**row, "built": time_text(_iso(row["built_at"])),
                         "rows_text": f"{row['rows']:,}" if row["rows"] is not None else ""}
                        for row in warm["index"]["kinds"]],
    }


def _iso(seconds):
    """A time given in seconds since the epoch as the ISO text `time_text` reads, or `None`."""
    if not seconds:
        return None
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(float(seconds)))


@bp.route("/admin/queue")
def admin_queue():
    """The queue's state and the runs, newest first (a data table)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    error = None
    jobs_table = None
    workers = warm = None
    try:
        try:
            workers = apiclient.admin_work_workers(cookie_header)
        except apiclient.ApiError:
            workers = None                      # Redis down: the page says so below; the rest still loads
        try:
            warm = apiclient.admin_work_warm(cookie_header, db_name())
        except apiclient.ApiError:
            warm = None
        body = apiclient.admin_work(cookie_header, limit=1, offset=0)
        if body["available"]:
            jobs_table = tables.render(JOBS_TABLE, request.path, anchor="jobs")
    except apiclient.ApiError as exc:
        body, error = {"available": False, "status": {}, "items": [], "total": 0}, _api_message(exc)
    flashed = _take_flash()
    return _render(
        "admin_queue.html", flashed=True, title="Work Queue", section="admin_queue",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Queue")],
        description="The generation runs on this server, and the controls to cancel, clear, retry or delete them.",
        admin=admin, available=body["available"], queue=_status_view(body["status"]), stopped=_stopped_view(body.get("stopped")), jobs_table=jobs_table, workers=workers,
        warm=_warm_view(warm),
        error=error or flashed.get("error"), message=flashed.get("message"),
    )


def _event_view(event):
    """The event as the logging layer's display line (`logrecord.display_line`, UX.123): the run's download holds the same
    record as a file line, so the page and the file show one event the same way."""
    return {"line": logrecord.display_line(runrecord.event_record(event))}


def _stage_rows(run):
    return [{
        "number": number,
        "label": stage["label"],
        "state": stage["state"],
        "state_label": STAGE_STATE_LABELS[stage["state"]],
        "units": f"{stage['done']} of {stage['units']} {counted('unit', stage['units'])}" if stage["units"] else "",
        "started": _iso(stage["started"]),
        "duration": _duration(stage["seconds"]),
        "longest": _duration(stage["longest_unit"]) if stage.get("longest_unit") else "",
        "peak": f"{stage['peak_mb']:.0f} MB" if stage.get("peak_mb") else "",
    } for number, stage in enumerate(run["stages"], 1)]


def _run_back_url():
    """Where Cancel returns to: the page the run's page was opened from, not a confirmation step or the page itself."""
    back = came_from(request.path)
    return "" if back.startswith("/admin/queue/confirm/") else back


@bp.route("/admin/queue/<run_id>")
def admin_queue_run(run_id):
    """One run: its facts, its stages and the last lines of its log."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    if not runrecord.RUN_ID_RE.match(run_id):
        abort(404, description="No such run.")
    body, run = _load(_cookie_header(), run_id)
    if run is None:
        abort(404, description="No such run.")
    flashed = _take_flash()
    return _render(
        "admin_queue_run.html", flashed=True, title=run["title"][:80] or "Run", section="admin_queue",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Queue", "admin_queue"), crumb("Run")],
        admin=admin, run=run_view(run), stages=_stage_rows(run), events=[_event_view(event) for event in body["events"]],
        back_url=_run_back_url(),
        error=flashed.get("error"), message=flashed.get("message"),
    )


@bp.route("/admin/queue/<run_id>/log")
def admin_queue_log(run_id):
    """A run's whole log as a plain-text download."""
    _admin, bounce = _require_admin()
    if bounce is not None:
        abort(403)   # a file download, not a page: a plain 403 like the other json_only views
    if not runrecord.RUN_ID_RE.match(run_id):
        abort(404, description="No such run.")
    try:
        text = apiclient.admin_work_run_log(_cookie_header(), run_id)
    except apiclient.NotFoundError:
        abort(404, description="No such run.")
    response = Response(text, mimetype="text/plain")
    response.headers["Content-Disposition"] = f'attachment; filename="run-{run_id}.log"'
    response.headers["Cache-Control"] = "no-store"
    return response


admin_queue_log.json_only = True  # not a page: tests/test_web_a11y.py skips it


def _confirm_text(action, run):
    """`(heading, explanation)` for the confirmation page."""
    label, explanation = ACTION_TEXT[action]
    if _is_site_wide(action):
        return f"{label}?", explanation
    if action == "retry":
        plan = _retry_plan(run)
        return (f"Retry “{run['title']}”?",
                f"Starts a Generate page job that runs {plan['describe']} on {db_name()}. Work that is "
                f"already done is skipped.")
    return f"{label} “{run['title']}”?", explanation


@bp.route("/admin/queue/confirm/<action>")
def admin_queue_confirm(action):
    """The confirmation page every control goes through first (`?run=<id>`)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    run_id = request.args.get("run") or ""
    if action not in ACTION_TEXT:
        abort(404, description="No such action.")
    run = None
    if not _is_site_wide(action):
        if not runrecord.RUN_ID_RE.match(run_id):
            abort(404, description="No such run.")
        _body, run = _load(_cookie_header(), run_id)
        if run is None:
            abort(404, description="That run is no longer there.")
        if action == "retry" and _retry_plan(run) is None:
            return _flash(_see_other(url_for("web.admin_queue_run", run_id=run_id)),
                          error="That run can't be retried from here.")
    heading, explanation = _confirm_text(action, run)
    back = url_for("web.admin_queue_run", run_id=run_id) if run else url_for("web.admin_queue")
    return _render(
        "admin_queue_confirm.html", title="Confirm", section="admin_queue",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Queue", "admin_queue"), crumb("Confirm")],
        admin=admin, heading=heading, explanation=explanation, action=action, run_id=run_id,
        button=ACTION_TEXT[action][0], back=back, next_url=safe_next(request.args.get("next"), default=""),
    )


def _start_retry(admin, plan, retry_of):
    """Starts `plan` (`_retry_plan`) as a Generate page job; behind a running job it waits in the list and starts by
    itself (PERF.154). Returns `(message, error)`."""
    database = db_name()
    if plan["database"] and plan["database"] != database:
        return None, (f"That run wrote {plan['database']}, but this site shows {database}. "
                      f"Run it again from the command line instead.")
    env = jobs.mysql_env(current_app.config["MYSQL_CONFIG"], database)
    try:
        submitted = jobs.submit(plan["kind"], plan["title"], plan["steps"], env=env, admin=admin.get("username"),
                                database=database)
    except runline.LineFull as exc:
        return None, str(exc)
    except OSError as exc:
        log.exception(f"Could not start {plan['title']}: {exc}")
        return None, f"The job could not be started: {exc}"
    if submitted.job_id is None:
        activity_log.event("GEN", "job.queue", user=admin.get("username"), kind=plan["kind"], db=database,
                           title=plan["title"], **({"retry_of": retry_of} if retry_of else {"started_from": "queue page"}))
        return f"{submitted.running} is running, so {plan['title']} waits its turn (place {submitted.position}) and starts by itself.", None
    activity_log.event("GEN", "job.start", user=admin.get("username"), job=submitted.job_id, kind=plan["kind"], db=database,
                       title=plan["title"], **({"retry_of": retry_of} if retry_of else {"started_from": "queue page"}))
    return f"Started {plan['title']}.", None


def _start_fill(admin, cookie_header, action):
    """
    Starts the `planetgen` command behind a start button (`FILL_COMMANDS`, ADM.129) as a Generate page job on this
    site's database, as a retry does (`_start_retry`). Refused while everything is stopped; behind a live run it waits its
    turn in the list (PERF.154). Returns `(message, error)`.
    """
    label = ACTION_TEXT[action][0]
    body = apiclient.admin_work(cookie_header, limit=1, offset=0)
    if (body.get("stopped") or {}).get("stopped"):
        return None, "Everything is stopped. Press Resume first; nothing starts until then."
    command = FILL_COMMANDS[action]
    plan = {"kind": command, "title": label, "database": db_name(),
            "steps": [{"label": f"planetgen {command}",
                       "argv": [jobs.python_executable(), *jobs.GENERATE_COMMAND, command]}]}
    message, error = _start_retry(admin, plan, None)
    return message, error


@bp.route("/admin/queue/action", methods=["POST"])
def admin_queue_action():
    """Runs a confirmed control, then back to the run (or the list)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    action = request.form.get("action") or ""
    run_id = request.form.get("run") or ""
    target = url_for("web.admin_queue")
    message = error = None
    try:
        if action == "queue-pause":
            ok = apiclient.admin_work_queue(cookie_header, "pause")
            message = "The queue is paused." if ok else "The queue was already paused."
        elif action == "queue-stop-all":
            done = apiclient.admin_work_stop_all(cookie_header)
            message = (f"Everything is stopped: {done['runs']} run(s) cancelled"
                       f"{' and the Generate job' if done['job'] else ''}, {done['queued']} queued job(s) removed, "
                       f"{done['killed']} process(es) ended. Nothing starts until you press Resume.")
            if done["survivors"]:
                error = (f"{len(done['survivors'])} process(es) belong to another user and could not be ended from here; "
                         f"run 'sudo planetgen queue stop-all' on the server (pids {', '.join(map(str, done['survivors']))}).")
        elif action == "queue-resume-all":
            ok = apiclient.admin_work_queue(cookie_header, "resume-all")
            message = "Resumed: work can start again." if ok else "Nothing was stopped."
        elif action == "queue-resume":
            ok = apiclient.admin_work_queue(cookie_header, "resume")
            message = "The queue is running again." if ok else "The queue wasn't paused."
        elif action == "clear-finished":
            cleared = apiclient.admin_work_clear_finished(cookie_header)
            message = f"Cleared {cleared} finished {'run' if cleared == 1 else 'runs'}." if cleared \
                else "There were no finished runs to clear."
        elif action in ("workers-set", "workers-ludicrous"):
            try:
                count = int(request.form.get("count") or "") if action == "workers-set" else None
            except ValueError:
                count = 0
            if action == "workers-set" and count < 1:
                error = "The worker count is a whole number of at least 1."
            else:
                state = apiclient.admin_work_workers(cookie_header, count=count, ludicrous=action == "workers-ludicrous")
                message = f"The bulk queue now asks for {state['wanted']} {'worker' if state['wanted'] == 1 else 'workers'}."
        elif action in FILL_COMMANDS:
            message, error = _start_fill(admin, cookie_header, action)
            if message:
                target = url_for("web.generate", _anchor="current-job")
        elif action in ("cancel", "clear", "retry", "delete"):
            _body, run = _load(cookie_header, run_id) if runrecord.RUN_ID_RE.match(run_id) else (None, None)
            if run is None:
                return _flash(_see_other(target), error="That run is no longer there.")
            target = url_for("web.admin_queue_run", run_id=run_id)
            if action == "delete":
                if apiclient.admin_work_delete(cookie_header, run_id):
                    return _flash(_see_other(url_for("web.admin_queue")), message="Deleted the run's record.")
                error = "Only a finished run can be deleted."
            elif action == "retry":
                plan = _retry_plan(run)
                if plan is None:
                    error = "That run can't be retried from here."
                else:
                    message, error = _start_retry(admin, plan, run_id)
                    if message:
                        target = url_for("web.generate", _anchor="current-job")
            else:
                ok = apiclient.admin_work_control(cookie_header, run_id, action)
                if ok:
                    message = f"{'Cancelled' if action == 'cancel' else 'Cleared'} “{run['title']}”."
                    if action == "cancel":     # back to the page the run's page was opened from, as the Generate job page does
                        target = safe_next(request.form.get("next"), default=target)
                elif action == "clear":
                    error = f"“{run['title']}” is still answering, or has finished. Cancel it instead."
                else:
                    message = f"“{run['title']}” has finished already."
        else:
            error = "Unknown action."
    except apiclient.ApiError as exc:
        error = _api_message(exc)
    # The API wrote each change to the audit and activity logs; a retry is logged as the job it starts (`_start_retry`).
    return _flash(_see_other(target), **({"error": error} if error else {"message": message}))
