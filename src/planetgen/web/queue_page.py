# planetgen/web/queue_page.py

"""
The admin work queue page (ADM.10): `/admin/queue` lists the job trees
(ADM.12) newest first, under the number of workers active and the
server's load ("x / x / x" over 1, 5 and 15 minutes; CPU percent on
Windows); `/admin/queue/<id>` shows one tree, expandable from the master
job down to single tasks. Every control (pause, resume, cancel, retry,
delete, pause the whole queue, clear a stale lease) first shows a
confirmation page (`/admin/queue/confirm/<action>`), then posts to
`/admin/queue/action`, and is written to the activity log.

Reads and controls go through the API (`api/workqueue.py`); Retry starts
a Generate page job (`web/jobs.py`), like the Generate page's own forms.
"""

import re

from flask import abort, current_app, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.pagination import fetch_page, parse_page
from planetgen.queue import work as workQueue
from planetgen.admin import activity_log
from planetgen.util import log

from . import bp, jobs
from .admin_pages import _api_message, _cookie_header, _flash, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb, db_name, pager

KIND_LABELS = {
    "web-job": "Generate page job",
    "step": "Step",
    "queue": "Work queue",
    "skeleton": "Skeleton",
    "bright-stars": "Bright stars",
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
    "queue-pause": ("Pause the queue", "No job starts or takes another task until you resume the queue. "
                                       "Running tasks finish first."),
    "queue-resume": ("Resume the queue", "Jobs waiting for the queue start again, one at a time."),
    "clear-lease": ("Clear the stale lease", "The run holding it has stopped answering. Its unfinished "
                                             "work is marked cancelled so the next job can start."),
}


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
        "eta": _duration(totals["eta_seconds"]),
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


def _status_view(status):
    load = status.get("load") or {}
    return {
        "workers_active": status.get("workers_active", 0),
        "runs_active": status.get("runs_active", 0),
        "load": load.get("text") or "unknown",
        "load_label": "CPU in use, 1 / 5 / 15 minutes" if load.get("kind") == "cpu"
        else "Load average, 1 / 5 / 15 minutes",
        "paused": status.get("paused", False),
        "paused_by": status.get("paused_by"),
        "paused_at": status.get("paused_at"),
        "holder": status.get("holder"),
        "job_id": status.get("job_id"),
        "heartbeat": _duration(status.get("heartbeat_age_s")),
        "stale": status.get("stale", False),
    }


def _task_retry_argv(task):
    """`generate.py` arguments that fill one failed sector again, or
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
    else, the same `generate.py` command line. Either way the whole job
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
        label = " ".join(["generate.py", *argv])
        return {"kind": root["kind"], "title": f"Retry: {label}"[:200],
                "steps": [{"label": label, "argv": [jobs.python_executable(), jobs.GENERATE_SCRIPT, *argv]}],
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


@bp.route("/admin/queue")
def admin_queue():
    """The queue's state and the job trees, newest first
    (`?page=N`)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    error = None
    try:
        body, page = fetch_page(
            lambda limit, offset: apiclient.admin_work(cookie_header, limit=limit, offset=offset),
            parse_page(request.args.get("page")),
        )
    except apiclient.ApiError as exc:
        body, page, error = {"available": False, "status": {}, "items": [], "total": 0}, 1, _api_message(exc)
    flashed = _take_flash()
    rows = []
    for root in body["items"]:
        root.setdefault("children", [])
        root.setdefault("tasks", [])
        rows.append(_node_view(root, root))
    return _render(
        "admin_queue.html", flashed=True, title="Work Queue", section="admin_queue",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Queue")],
        description="The generation jobs on this server, and the controls to pause, resume, cancel or retry them.",
        admin=admin, available=body["available"], queue=_status_view(body["status"]), rows=rows,
        pager=pager("page", page, body["total"], label="Job pages"), error=error or flashed.get("error"),
        message=flashed.get("message"),
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
    if action.startswith("queue-") or action == "clear-lease":
        return f"{label}?", explanation
    if task is not None:
        argv = _task_retry_argv(task)
        return (f"Retry sector {task['task_key']}?",
                f"Starts a Generate page job that runs generate.py {' '.join(argv[:1] + argv[1:])} on "
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
    if not (action.startswith("queue-") or action == "clear-lease"):
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
                        label = " ".join(["generate.py", *argv])
                        plan = {"kind": "galaxy", "title": f"Retry: {label}", "database": tree["database_name"],
                                "steps": [{"label": label, "argv": [jobs.python_executable(),
                                                                    jobs.GENERATE_SCRIPT, *argv]}]}
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
                    # A step that isn't a generate.py run (a reset) only
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
