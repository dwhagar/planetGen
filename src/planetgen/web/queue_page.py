# planetgen/web/queue_page.py

"""
The admin work queue page (ADM.10): `/admin/queue` lists the job trees
(ADM.12) newest first, under the number of workers active and the
server's load ("x / x / x" over 1, 5 and 15 minutes); `/admin/queue/<id>` shows one tree, expandable from the master
job down to single tasks. Every control (pause, resume, cancel, retry,
delete, pause the whole queue, clear a stale lease) first shows a
confirmation page (`/admin/queue/confirm/<action>`), then posts to
`/admin/queue/action`, and is written to the activity log.

Reads and controls go through the API (`api/workqueue.py`); Retry starts
a Generate page job (`web/jobs.py`), like the Generate page's own forms.
"""

import re
import time

from flask import abort, current_app, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import Column, Result, Table, time_text
from planetgen.queue import progress_rate, work as workQueue
from planetgen.admin import activity_log
from planetgen.util import log

from . import bp, jobs, tables
from .admin_pages import _api_message, _cookie_header, _flash, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb, db_name

KIND_LABELS = {
    "web-job": "Generate page job",
    "step": "Step",
    "queue": "Work queue",
    "skeleton": "Skeleton",
    "bright-stars": "Bright stars",
    "phenomena": "Phenomena",
    "population": "Population",
    "plan": "Plan",
    "galaxy": "Galaxy",
    "sector": "Sector",
    "system": "System",
    "phenomenon": "Phenomenon",
}

STATUS_LABELS = {
    "waiting": "Waiting",
    "running": "Running",
    "paused": "Paused",
    "done": "Done",
    "failed": "Failed",
    "cancelled": "Cancelled",
    "interrupted": "Interrupted",
}

SECTOR_KEY_RE = re.compile(r"^(-?\d+),(-?\d+),(-?\d+)$")
"""re.Pattern: A sector task's key, `ring,layer,slot`."""

ACTION_TEXT = {
    "pause": ("Pause", "Its running tasks finish, then it waits in standby without holding up other jobs, "
                       "until you resume it."),
    "resume": ("Resume", "It takes its place in line again and carries on."),
    "cancel": ("Cancel", "Its running tasks finish and nothing more is started. Cancelling a whole run, "
                         "a step or a Generate page job ends that run."),
    "retry": ("Retry", None),
    "delete": ("Delete", "The job's record and timings are deleted. Nothing it generated is touched."),
    "clear-finished": ("Clear finished jobs", "Every finished job's record and timings are deleted, "
                                              "leaving the ones still running. Nothing they generated is touched."),
    "queue-pause": ("Pause the queue", "No job starts or takes another task until you resume the queue. "
                                       "Running tasks finish first."),
    "queue-resume": ("Resume the queue", "Jobs waiting for the queue start again, one at a time."),
    "clear-lease": ("Clear the stale lease", "The run holding it has stopped answering. Its unfinished "
                                             "work is marked cancelled so the next job can start."),
}


QUEUE_WIDE_ACTIONS = ("clear-lease", "clear-finished")
"""tuple[str]: The controls that act on the whole queue, not one job."""


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


def _actions(node, root):
    """The controls an admin can use on `node` now."""
    actions = []
    if node["live"]:
        if node["control"] == "pause" or node["state"] == "paused":
            actions.append("resume")
        elif node["control"] != "cancel":
            actions.append("pause")
        if node["control"] != "cancel":
            actions.append("cancel")
    elif node["status"] in ("failed", "cancelled", "interrupted") and _retry_plan(root, node) is not None:
        actions.append("retry")
    if node is root and not node["live"]:
        actions.append("delete")
    return actions


def _node_view(node, root):
    totals = node["totals"]
    seconds = node["seconds"] if node["seconds"] is not None else (node["age_seconds"] if node["live"] else None)
    view = {
        "id": node["id"],
        "kind": KIND_LABELS.get(node["kind"], node["kind"]),
        "title": node["title"],
        "status": node["status"],
        "status_label": STATUS_LABELS.get(node["status"], node["status"]),
        "asked": {"pause": "pause asked", "cancel": "cancel asked"}.get(node["control"]) if node["live"] else None,
        "started_at": node["started_at"] or node["created_at"],
        "finished_at": node["finished_at"],
        "duration": _duration(seconds),
        "work": _duration(totals["work_seconds"]) if totals["work_seconds"] else "",
        "tasks": totals["tasks"],
        "queued": totals["queued"],
        "running": totals["running"],
        "done": totals["done"],
        "failed": totals["failed"],
        "cancelled": totals["cancelled"],
        "eta": _eta_text(totals["eta_seconds"], node["live"]),
        "workers": node["workers"],
        "holder": node["holder"],
        "web_job_url": url_for("web.generate_job", job_id=node["web_job_id"])
        if node["web_job_id"] and node["kind"] == "web-job" and jobs.JOB_ID_RE.match(node["web_job_id"]) else None,
        "url": url_for("web.admin_queue_tree", node_id=root["id"], _anchor=f"node-{node['id']}"),
        "actions": [(action, ACTION_TEXT[action][0]) for action in _actions(node, root)],
        "live": node["live"],
        "children": [_node_view(child, root) for child in node["children"]],
        "task_rows": [{
            "id": task["id"],
            "key": task["task_key"],
            "kind": task["kind"],
            "status": task["state"],
            "status_label": STATUS_LABELS.get(task["state"], task["state"]),
            "started_at": task["started_at"],
            "seconds": _duration(task["seconds"]),
            "error": task["error"],
            "retry": task["state"] in ("failed", "cancelled") and _task_retry_argv(task) is not None
            and not root["live"],
        } for task in node["tasks"]],
    }
    return view


def _eta_text(seconds, live):
    """The time left as a range around the estimate, "estimating" while a running job has none (PERF.33, never
    dashes), `""` for a job that isn't running."""
    if seconds is None:
        return "estimating" if live else ""
    low, high = progress_rate.eta_range(seconds)
    low_text, high_text = _duration(low), _duration(high)
    return high_text if low_text == high_text else f"{low_text} to {high_text}"


def _remaining_seconds(root):
    """
    Seconds a running job has left, or `None` when it can't be told (ADM.39):
    the job tree's own estimate from its finished tasks' pace, else the
    estimate a Generate page job's runner publishes (PERF.7), counted down
    from when it was written. A job that isn't running has none.
    """
    if not root["live"]:
        return None
    eta = root["totals"]["eta_seconds"]
    if eta is not None:
        return eta
    web_job_id = root.get("web_job_id")
    if root["kind"] != "web-job" or not web_job_id or not jobs.JOB_ID_RE.match(web_job_id):
        return None
    job = jobs.get_job(web_job_id)
    progress = (job or {}).get("progress") or {}
    if job is None or job.get("finished_at") or progress.get("eta_s") is None:
        return None
    updated_at = progress.get("updated_at") or time.time()
    return max(0.0, float(progress["eta_s"]) - max(0.0, time.time() - float(updated_at)))


def _status_view(status):
    load = status.get("load") or {}
    return {
        "workers_active": status.get("workers_active", 0),
        "runs_active": status.get("runs_active", 0),
        "load": load.get("text") or "unknown",
        "load_label": "Load average, 1 / 5 / 15 minutes",
        "paused": status.get("paused", False),
        "paused_by": status.get("paused_by"),
        "paused_at": status.get("paused_at"),
        "holder": status.get("holder"),
        "job_id": status.get("job_id"),
        "heartbeat": _duration(status.get("heartbeat_age_s")),
        "stale": status.get("stale", False),
    }


def _task_retry_argv(task):
    """`planetgen` arguments that fill one failed sector again, or
    `None` for a task that can't be rerun alone."""
    match = SECTOR_KEY_RE.match(task.get("task_key") or "")
    if task.get("kind") != "sector" or not match:
        return None
    ring, layer, slot = match.groups()
    return ["galaxy", "--ring", ring, "--layer", layer, "--slot", slot]


def _retry_plan(root, node=None):
    """
    What Retry on `node` runs again, as `{"kind", "title", "steps",
    "database", "describe"}`, or `None`: for a Generate page job, its
    steps from the one that didn't finish; for a run started anywhere
    else, the same `planetgen` command line. Either way the whole job
    is retried, not one node of it: generation skips what is already
    there, so only the missing work is done again.
    """
    if root["kind"] == "web-job" and root["web_job_id"]:
        try:
            found = jobs.remaining_steps(root["web_job_id"])
        except OSError:
            found = None
        if found is None:
            return None
        job, steps = found
        return {"kind": job.get("kind") or "galaxy", "title": f"Retry: {job.get('title') or root['title']}",
                "steps": steps, "database": job.get("database"),
                "describe": ", then ".join(step["label"] for step in steps)}
    if root["argv"] and root["kind"] in ("plan", "galaxy", "sector", "population"):
        argv = [str(arg) for arg in root["argv"]]
        label = " ".join(["planetgen", *argv])
        return {"kind": root["kind"], "title": f"Retry: {label}"[:200],
                "steps": [{"label": label, "argv": [jobs.python_executable(), *jobs.GENERATE_COMMAND, *argv]}],
                "database": root["database_name"], "describe": label}
    return None


def _find(tree, node_id):
    if tree["id"] == node_id:
        return tree
    for child in tree["children"]:
        found = _find(child, node_id)
        if found is not None:
            return found
    return None


def _load(cookie_header, node_id):
    """`(tree, node)` from the API, or `(None, None)`."""
    try:
        body = apiclient.admin_work_tree(cookie_header, node_id)
    except apiclient.NotFoundError:
        return None, None
    tree = body["tree"]
    return tree, _find(tree, node_id)


def _job_cells(root):
    """One job (a tree's root) as the Jobs table's cells."""
    node = _node_view(root, root)
    links = []
    for action, label in node["actions"]:
        links += [" "] if links else []
        links.append({"text": label, "href": url_for("web.admin_queue_confirm", action=action, node=node["id"])})
    if node["web_job_url"]:
        links += [" "] if links else []
        links.append({"text": "Log", "href": node["web_job_url"]})
    started = time_text(node["started_at"]) or "not yet"
    asked = f" ({node['asked']})" if node["asked"] else ""
    return [
        {"text": node["title"], "href": node["url"]},
        {"text": node["status_label"] + asked},
        {"text": started},
        {"text": node["duration"] or ""},
        {"text": _eta_text(_remaining_seconds(root), root["live"])},
        {"text": f"{node['done']} of {node['tasks']}" if node["tasks"] else ""},
        {"text": "", "parts": links} if links else {"text": ""},
    ]


def _jobs_load(state, limit, offset, want_facets):
    """The job trees, newest first, `limit` from `offset`."""
    _identity, bounce = _require_admin()
    if bounce is not None:
        abort(403)
    body = apiclient.admin_work(_cookie_header(), limit=limit, offset=offset)
    rows = []
    for root in body["items"]:
        root.setdefault("children", [])
        root.setdefault("tasks", [])
        rows.append(_job_cells(root))
    return Result(rows, body["total"], None)


JOBS_TABLE = tables.register(Table(
    "queue-jobs", "Jobs",
    [Column("job", "Job", sortable=False), Column("status", "Status", sortable=False), Column("started", "Started", sortable=False),
     Column("duration", "Duration", sortable=False), Column("left", "Time left", sortable=False),
     Column("progress", "Progress", sortable=False),
     Column("actions", "Actions", sortable=False)],
    _jobs_load, prefix="jobs_", noun=("job", "jobs"),
))


@bp.route("/admin/queue")
def admin_queue():
    """The queue's state and the job trees, newest first (a data table)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    error = None
    jobs_table = None
    try:
        body = apiclient.admin_work(cookie_header, limit=1, offset=0)
        if body["available"]:
            jobs_table = tables.render(JOBS_TABLE, request.path, anchor="jobs")
    except apiclient.ApiError as exc:
        body, error = {"available": False, "status": {}, "items": [], "total": 0}, _api_message(exc)
    flashed = _take_flash()
    return _render(
        "admin_queue.html", flashed=True, title="Work Queue", section="admin_queue",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Queue")],
        description="The generation jobs on this server, and the controls to pause, resume, cancel or retry them.",
        admin=admin, available=body["available"], queue=_status_view(body["status"]), jobs_table=jobs_table,
        error=error or flashed.get("error"), message=flashed.get("message"),
    )


@bp.route("/admin/queue/<node_id>")
def admin_queue_tree(node_id):
    """One job tree, every node with its timing and controls."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    tree, _node = _load(_cookie_header(), node_id)
    if tree is None:
        abort(404, description="No such job.")
    flashed = _take_flash()
    return _render(
        "admin_queue_tree.html", flashed=True, title=tree["title"][:80] or "Job", section="admin_queue",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Queue", "admin_queue"), crumb("Job")],
        admin=admin, root=_node_view(tree, tree), error=flashed.get("error"), message=flashed.get("message"),
    )


def _confirm_text(action, tree, node, task):
    """`(heading, explanation)` for the confirmation page."""
    label, explanation = ACTION_TEXT[action]
    if action.startswith("queue-") or action in QUEUE_WIDE_ACTIONS:
        return f"{label}?", explanation
    if task is not None:
        argv = _task_retry_argv(task)
        return (f"Retry sector {task['task_key']}?",
                f"Starts a Generate page job that runs planetgen {' '.join(argv[:1] + argv[1:])} on "
                f"{db_name()}.")
    if action == "retry":
        plan = _retry_plan(tree, node)
        return (f"Retry “{tree['title']}”?",
                f"Starts a Generate page job that runs {plan['describe']} on {db_name()}. Work that is "
                f"already done is skipped.")
    return f"{label} “{node['title']}” and everything under it?", explanation


@bp.route("/admin/queue/confirm/<action>")
def admin_queue_confirm(action):
    """The confirmation page every control goes through first
    (`?node=<id>`, and `&task=<id>` to retry one sector)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    node_id = request.args.get("node") or ""
    if action not in ACTION_TEXT:
        abort(404, description="No such action.")
    tree = node = task = None
    if not (action.startswith("queue-") or action in QUEUE_WIDE_ACTIONS):
        if not workQueue.NODE_ID_RE.match(node_id):
            abort(404, description="No such job.")
        tree, node = _load(_cookie_header(), node_id)
        if node is None:
            abort(404, description="That job is no longer there.")
        task_id = request.args.get("task")
        if task_id:
            task = next((row for row in node["tasks"] if str(row["id"]) == task_id), None)
            if task is None or _task_retry_argv(task) is None:
                return _flash(_see_other(url_for("web.admin_queue_tree", node_id=tree["id"])),
                              error="That task can't be retried on its own.")
        elif action == "retry" and _retry_plan(tree, node) is None:
            return _flash(_see_other(url_for("web.admin_queue_tree", node_id=tree["id"])),
                          error="That job can't be retried from here.")
    heading, explanation = _confirm_text(action, tree, node, task)
    back = url_for("web.admin_queue_tree", node_id=tree["id"]) if tree else url_for("web.admin_queue")
    return _render(
        "admin_queue_confirm.html", title="Confirm", section="admin_queue",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Queue", "admin_queue"), crumb("Confirm")],
        admin=admin, heading=heading, explanation=explanation, action=action, node_id=node_id,
        task_id=task["id"] if task else None, button=ACTION_TEXT[action][0], back=back,
    )


def _start_retry(admin, plan, retry_of):
    """Starts `plan` (`_retry_plan`, or a task's) as a Generate page job.
    Returns `(message, error)`."""
    database = db_name()
    if plan["database"] and plan["database"] != database:
        return None, (f"That job wrote {plan['database']}, but this site shows {database}. "
                      f"Run it again from the command line instead.")
    env = jobs.mysql_env(current_app.config["MYSQL_CONFIG"], database)
    try:
        job_id = jobs.start_job(plan["kind"], plan["title"], plan["steps"], env=env, admin=admin.get("username"),
                                database=database)
    except jobs.JobBusy as exc:
        running = exc.job["title"] if exc.job else "Another job"
        return None, f"{running} is still running. Wait for it to finish, or cancel it first."
    except OSError as exc:
        log.exception(f"Could not start a retry of {retry_of}: {exc}")
        return None, f"The job could not be started: {exc}"
    activity_log.event("GEN", "job.start", user=admin.get("username"), job=job_id, kind=plan["kind"], db=database,
                      title=plan["title"], retry_of=retry_of)
    return f"Started {plan['title']}.", None


@bp.route("/admin/queue/action", methods=["POST"])
def admin_queue_action():
    """Runs a confirmed control, then back to the tree (or the list)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    action = request.form.get("action") or ""
    node_id = request.form.get("node") or ""
    target = url_for("web.admin_queue")
    message = error = None
    try:
        if action == "queue-pause":
            ok = apiclient.admin_work_queue(cookie_header, "pause")
            message = "The queue is paused." if ok else "The queue was already paused."
        elif action == "queue-resume":
            ok = apiclient.admin_work_queue(cookie_header, "resume")
            message = "The queue is running again." if ok else "The queue wasn't paused."
        elif action == "clear-finished":
            cleared = apiclient.admin_work_clear_finished(cookie_header)
            message = f"Cleared {cleared} finished {'job' if cleared == 1 else 'jobs'}." if cleared \
                else "There were no finished jobs to clear."
        elif action == "clear-lease":
            holder = apiclient.admin_work_clear_lease(cookie_header)
            message = f"Cleared the lease held by {holder}." if holder else "The lease wasn't stale."
        elif action in ("pause", "resume", "cancel", "retry", "delete"):
            tree, node = _load(cookie_header, node_id)
            if node is None:
                return _flash(_see_other(target), error="That job is no longer there.")
            target = url_for("web.admin_queue_tree", node_id=tree["id"])
            if action == "delete":
                if apiclient.admin_work_delete(cookie_header, tree["id"]):
                    return _flash(_see_other(url_for("web.admin_queue")), message="Deleted the job's record.")
                error = "Only a finished job can be deleted."
            elif action == "retry":
                task_id = request.form.get("task")
                if task_id:
                    task = next((row for row in node["tasks"] if str(row["id"]) == task_id), None)
                    argv = _task_retry_argv(task) if task else None
                    if argv is None:
                        error = "That task can't be retried on its own."
                    else:
                        label = " ".join(["planetgen", *argv])
                        plan = {"kind": "galaxy", "title": f"Retry: {label}", "database": tree["database_name"],
                                "steps": [{"label": label, "argv": [jobs.python_executable(),
                                                                    *jobs.GENERATE_COMMAND, *argv]}]}
                        message, error = _start_retry(admin, plan, f"{node_id}/{task_id}")
                else:
                    plan = _retry_plan(tree, node)
                    if plan is None:
                        error = "That job can't be retried from here."
                    else:
                        message, error = _start_retry(admin, plan, node_id)
                if message:
                    target = url_for("web.generate", _anchor="current-job")
            else:
                ok = apiclient.admin_work_control(cookie_header, node_id, action)
                if action == "cancel" and node["kind"] in ("web-job", "step") and node["web_job_id"]:
                    # A step that isn't a planetgen run (a reset) only
                    # stops through the Generate page's own cancel.
                    try:
                        ok = jobs.cancel_job(node["web_job_id"]) or ok
                    except OSError as exc:
                        log.error(f"Could not cancel job {node['web_job_id']}: {exc}")
                past = {"pause": "pause", "resume": "resume", "cancel": "cancel"}[action]
                message = (f"Asked “{node['title']}” to {past}." if ok
                           else f"“{node['title']}” has finished or already does that.")
        else:
            error = "Unknown action."
    except apiclient.ApiError as exc:
        error = _api_message(exc)
    # The API wrote each change to the audit and activity logs; a retry
    # is logged as the job it starts (`_start_retry`).
    return _flash(_see_other(target), **({"error": error} if error else {"message": message}))
