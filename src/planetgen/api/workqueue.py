# planetgen/api/workqueue.py

"""
Admin endpoints behind the work queue page (ADM.10, `web/queue_page.py`):
the queue's state (workers active, the server's load, the lease, the
whole-queue pause), the job trees (ADM.12, `planetgen/queue/work.py`)
and the controls on them. Every route needs a logged-in admin past the
forced credential change; every change is written to the audit and
activity logs (`audit`).

Rerunning a job (Retry) starts a Generate page job, so it lives with the
page (`web/queue_page.py`), not here.
"""

from flask import Blueprint, g, jsonify, request

from planetgen.queue import load as systemLoad, work as workQueue

from .authz import audit, require_admin
from .common import ApiError, get_control_db, require_json_body
from .routes import _paginate
from .schemas import QueueAction, WorkAction, parse_body

bp = Blueprint("workqueue", __name__, url_prefix="/api/admin/work")


def _available(conn):
    """Whether the control schema has the job tree (v7)."""
    try:
        conn.execute("SELECT control FROM work_jobs LIMIT 1").fetchall()
        conn.execute("SELECT paused FROM work_lease LIMIT 1").fetchall()
        conn.rollback()
        return True
    except Exception:  # noqa: BLE001 -- update.sh not run yet
        conn.rollback()
        return False


def _node_id(value):
    if not isinstance(value, str) or not workQueue.NODE_ID_RE.match(value):
        raise ApiError("No such job.", 404)
    return value


def _summary(node):
    """A root as the list shows it: its row and totals, no children."""
    return {key: value for key, value in node.items() if key not in ("children", "tasks")}


@bp.route("")
@require_admin(fresh=True)
def work():
    """
    `GET /api/admin/work[?limit=&offset=]` -- the queue and one page of
    job trees, newest first: `{"available", "status": {...queue_status,
    "load": systemLoad.load_average()}, "items": [roots with totals],
    "total", "limit", "offset"}`. `available` is false (and the rest
    empty) before `update.sh` brings the control schema to v7.
    """
    limit, offset = _paginate(request.args)
    load = systemLoad.load_average()
    conn = get_control_db()
    if not _available(conn):
        return jsonify({"available": False, "status": {"load": load}, "items": [], "total": 0,
                        "limit": limit, "offset": offset})
    status = workQueue.queue_status(conn)
    status["load"] = load
    roots, total = workQueue.list_roots(conn, limit=limit, offset=offset)
    items = []
    for root in roots:
        tree = workQueue.load_tree(conn, root["id"], max_tasks=0)
        items.append(_summary(tree if tree is not None else root))
    return jsonify({"available": True, "status": status, "items": items, "total": total,
                    "limit": limit, "offset": offset})


@bp.route("/<node_id>")
@require_admin(fresh=True)
def tree(node_id):
    """`GET /api/admin/work/<id>` -- the whole tree holding that node
    (`workQueue.load_tree`), as `{"tree": ..., "node_id": id}`."""
    conn = get_control_db()
    node = workQueue.get_node(conn, _node_id(node_id))
    if node is None:
        raise ApiError("No such job.", 404)
    return jsonify({"tree": workQueue.load_tree(conn, node["root_id"]), "node_id": node_id})


@bp.route("/<node_id>/control", methods=["POST"])
@require_admin(fresh=True)
def control(node_id):
    """
    `POST /api/admin/work/<id>/control` `{"action": "pause"|"resume"|
    "cancel"}` asks a live node and everything under it to pause (finish
    the running tasks, then stand by without blocking other jobs), resume
    or cancel. Returns `{"ok": bool}`: false when the node has finished
    or already does that.
    """
    action = parse_body(WorkAction, require_json_body()).action
    conn = get_control_db()
    node = workQueue.get_node(conn, _node_id(node_id))
    if node is None:
        raise ApiError("No such job.", 404)
    ok = workQueue.request_control(conn, node_id, action)
    if ok:
        audit(f"work.{action}", target=f"job:{node_id}", detail=node["title"][:200])
    return jsonify({"ok": ok})


@bp.route("/<node_id>/delete", methods=["POST"])
@require_admin(fresh=True)
def delete(node_id):
    """`POST /api/admin/work/<id>/delete` deletes a finished job tree
    (its root's id) with its tasks. `{"ok": bool}`: false for a live
    tree or a node that isn't a root."""
    conn = get_control_db()
    ok = workQueue.delete_tree(conn, _node_id(node_id))
    if ok:
        audit("work.delete", target=f"job:{node_id}")
    return jsonify({"ok": ok})


@bp.route("/clear-finished", methods=["POST"])
@require_admin(fresh=True)
def clear_finished():
    """`POST /api/admin/work/clear-finished` deletes every finished job
    tree with its tasks (ADM.41). `{"cleared": <count>}`; live jobs stay."""
    conn = get_control_db()
    cleared = workQueue.clear_finished(conn)
    if cleared:
        audit("work.clear", detail=f"{cleared} finished jobs")
    return jsonify({"cleared": cleared})


@bp.route("/queue", methods=["POST"])
@require_admin(fresh=True)
def queue():
    """
    `POST /api/admin/work/queue` `{"action": "pause"|"resume"}` pauses
    the whole queue (no run takes the lease or a task until it's resumed;
    running tasks finish) or resumes it. `{"ok": bool}`: false when it
    already was.
    """
    action = parse_body(QueueAction, require_json_body()).action
    conn = get_control_db()
    if action == "pause":
        ok = workQueue.pause_queue(conn, g.admin_user["username"])
    else:
        ok = workQueue.resume_queue(conn)
    if ok:
        audit(f"work.queue.{action}", target="queue")
    return jsonify({"ok": ok})


@bp.route("/lease/clear", methods=["POST"])
@require_admin(fresh=True)
def clear_lease():
    """`POST /api/admin/work/lease/clear` frees a lease whose holder
    stopped refreshing it. `{"cleared": holder or null}`."""
    holder = workQueue.clear_stale_lease(get_control_db())
    if holder:
        audit("work.lease.clear", target="lease", detail=holder[:200])
    return jsonify({"cleared": holder})
