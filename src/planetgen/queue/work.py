# planetgen/queue/work.py

"""
The job tree of a run and the machine's worker count.

Every job is a tree (ADM.12, control schema v7). `job_node` opens a node
under the one this process has open (a `planetgen` run is the root; its
phases, such as the bright stars, are nodes under it). A process started
by the Generate page's job runner hangs its root under the runner's step
node (`PARENT_ENV_VAR`). Each node keeps its state, start, end and
duration; `load_tree` reads a whole tree back with every parent's
totals, timings and ETA added up from its children. The nodes are what
the admin queue page and the API show, pause, resume and cancel.

The units of work themselves are RQ jobs on Redis (`planetgen.queue.units`,
`planetgen.queue.wrapper`); this file keeps `worker_count` (at most about
80% of the machine's cores), `task_seed` (a task's `random` seed from the
run's seed and its own key, so what it generates doesn't depend on how many
workers there are or which finished first) and the control database's
`work_lease`, which says which run the machine is working for: a run holds it
while its workers run, refreshes it every `HEARTBEAT_SECONDS`, and frees it
when it finishes; a lease not refreshed for `STALE_SECONDS` belongs to a run
that died, and the next run takes it.
"""

import contextlib
import hashlib
import ipaddress
import math
import os
import re
import secrets
import socket
import threading
import time

from planetgen.admin import updateguard
from planetgen.queue import pause, progress_rate
from planetgen.util import log

CPU_SHARE = 0.8
"""float: The share of the machine's cores the pool may use."""

WORKERS_ENV_VAR = "PLANETGEN_WORKERS"
"""str: Environment variable overriding the worker count (`--workers`
on the command line wins over it); 0 or "auto" means `CPU_SHARE` of the
cores."""

WORKER_NICENESS = 10
"""int: How much `os.nice` lowers a worker's priority on Linux and
macOS."""

CHUNKS_PER_WORKER = 8
"""int: Most items one chunk of a bulk stage holds is `1 / (workers * CHUNKS_PER_WORKER)` of them. Every
task a pool runs costs about 60 database statements and a fork, which is more than a bright-star shell costs to
draw, so thousands of tiny shells go out in few chunks (PERF.79)."""

HEARTBEAT_SECONDS = 5
STALE_SECONDS = 30
KEEP_DAYS = 7
ARGV_LENGTH = 1024
"""int: The longest command-line argument a job node keeps."""

CANCELLED_EXIT_CODE = 130
"""int: The exit status of a run cancelled from the admin queue page
(`Cancelled`); the job runner reads it as "cancelled"."""

PARENT_ENV_VAR = "PLANETGEN_WORK_PARENT"
"""str: Environment variable naming the job tree node a process's root
node goes under (the job runner sets it to each step's node)."""

LIVE_STATES = ("waiting", "running", "paused")
"""tuple: Node states of a run that hasn't finished (yet)."""

NODE_ID_RE = re.compile(r"^[0-9]{8}-[0-9]{6}-[0-9a-f]{8}$")
"""re.Pattern: A node id (`new_node_id`)."""


def cpu_count():
    """The cores this process may run on (its affinity mask where the OS
    has one), at least 1."""
    if hasattr(os, "process_cpu_count"):
        count = os.process_cpu_count()
    elif hasattr(os, "sched_getaffinity"):
        count = len(os.sched_getaffinity(0))
    else:
        count = os.cpu_count()
    return max(1, count or 1)


def _is_local_host(host):
    if not host:
        return False
    if host.lower() in ("localhost", socket.gethostname().lower()):
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def worker_count(requested=None, mysql_host=None):
    """
    How many worker processes a run uses.

    Args:
        requested (int | str | None): `--workers`, else the
            `PLANETGEN_WORKERS` environment variable; a positive number
            is used as given, and 0, "auto" or nothing means `CPU_SHARE`
            of `cpu_count()`.
        mysql_host (str, optional): The database host. When it's this
            machine, one fewer worker, so MySQL keeps a core of the 80%.

    Returns:
        int: At least 1.
    """
    if requested is None:
        requested = os.environ.get(WORKERS_ENV_VAR)
    if isinstance(requested, str):
        requested = requested.strip().lower()
        requested = 0 if requested in ("", "auto") else int(requested)
    if requested:
        return max(1, int(requested))
    count = max(1, math.floor(CPU_SHARE * cpu_count()))
    if count > 1 and _is_local_host(mysql_host):
        count -= 1
    return count


def lower_priority(niceness=WORKER_NICENESS):
    """Lowers this process's scheduling priority (by `WORKER_NICENESS`, or by `niceness`) so everything else on
    the machine comes first. Never raises."""
    try:
        os.nice(niceness)
    except Exception:  # noqa: BLE001 -- a normal-priority worker is still a worker
        pass




def task_seed(run_seed, key):
    """The `random` seed for one task: the run's seed and the task's key,
    hashed, so it doesn't depend on worker count or finish order."""
    digest = hashlib.sha256(f"{run_seed}:{key}".encode("utf-8")).digest()
    return int.from_bytes(digest[:16], "big")


class Cancelled(SystemExit):
    """Raised in a run when an admin cancelled a node above its work
    queue (ADM.10): the tasks already running finish, nothing more is
    handed out, and the run ends as cancelled with
    `CANCELLED_EXIT_CODE`."""

    def __init__(self, title):
        super().__init__(CANCELLED_EXIT_CODE)
        self.title = title

    def __str__(self):
        return f"{self.title} was cancelled from the admin queue page."











STATS_KIND = {"sector": "sector", "mass-stars": "scatter-mass", "bright-stars": "scatter-bright",
              "phenomena": "phenomena"}
"""dict: Work queue task kind -> the generation stats kind its speeds are recorded under."""


class _Task:
    __slots__ = ("kind", "key", "fn", "payload", "weight", "on_done", "row_id", "future", "seed", "expected")

    def __init__(self, kind, key, fn, payload, weight, on_done, expected=None):
        self.kind, self.key, self.fn, self.payload = kind, key, fn, payload
        self.weight, self.on_done, self.expected = weight, on_done, expected
        self.row_id = self.future = self.seed = None


class _NoStore:
    """Bookkeeping when there is no control database: nothing kept."""

    available = False

    def create_node(self, node, state):
        pass

    def node_chain(self, node_id):
        return []

    def beat_nodes(self, ids):
        pass

    def finish_node(self, node_id, state):
        pass

    def controls(self, chain):
        return None

    def set_state(self, job_id, state):
        pass

    def heartbeat(self, job_id, holder):
        pass

    def add_tasks(self, job_id, tasks):
        pass

    def start_tasks(self, tasks):
        pass

    def finish_task(self, job_id, task, state, seconds=None, error=None):
        pass


class _ControlStore:
    """`work_jobs`/`work_tasks`/`work_lease` in the control database."""

    def __init__(self, config):
        from planetgen.db import store as _db

        self._db = _db
        self.config = config
        conn = self._connect()
        try:
            # v7's columns too: before `update.sh` brings the control
            # schema up to date, runs go on without the rows and lease.
            conn.execute("SELECT 1 FROM work_lease LIMIT 1").fetchall()
            conn.execute("SELECT parent_id, control FROM work_jobs LIMIT 1").fetchall()
            conn.execute("SELECT holder FROM work_lease LIMIT 1").fetchall()
        finally:
            conn.close()

    def _connect(self):
        return self._db.get_control_connection(self.config)

    available = True

    def create_node(self, node, state):
        """Inserts `node`'s row (after deleting finished trees older than
        `KEEP_DAYS` when it's a root of this process)."""
        conn = self._connect()
        try:
            with conn:
                if node.parent_id is None:
                    conn.execute(
                        "DELETE FROM work_jobs WHERE parent_id IS NULL AND created_at < NOW(6) - INTERVAL ? DAY"
                        " AND (state NOT IN ('waiting', 'running', 'paused')"
                        " OR heartbeat_at < NOW(6) - INTERVAL ? SECOND)",
                        (KEEP_DAYS, STALE_SECONDS),
                    )
                conn.execute(
                    "INSERT INTO work_jobs (id, parent_id, root_id, kind, title, holder, state, workers,"
                    " web_job_id, database_name, created_at, started_at, heartbeat_at)"
                    " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NOW(6), "
                    + ("NOW(6)" if state == "running" else "NULL") + ", NOW(6))",
                    (node.id, node.parent_id, node.root_id, node.kind[:32], node.title[:255], node.holder, state,
                     node.workers, node.web_job_id, node.database),
                )
                for position, argument in enumerate(node.argv or ()):
                    conn.execute("INSERT INTO work_job_args (job_id, position, value) VALUES (?, ?, ?)",
                                 (node.id, position, argument))
        finally:
            conn.close()

    def node_chain(self, node_id):
        """`[(id, root_id)]` from `node_id` up to its root, or `[]` when
        there's no such node."""
        chain = []
        conn = self._connect()
        try:
            while node_id is not None and len(chain) < 32:
                row = conn.execute("SELECT id, parent_id, root_id FROM work_jobs WHERE id = ?",
                                   (node_id,)).fetchone()
                if row is None:
                    break
                chain.append((row["id"], row["root_id"] or row["id"]))
                node_id = row["parent_id"]
            conn.rollback()
        finally:
            conn.close()
        return chain

    def beat_nodes(self, ids):
        if not ids:
            return
        conn = self._connect()
        try:
            with conn:
                conn.execute(
                    f"UPDATE work_jobs SET heartbeat_at = NOW(6) WHERE id IN ({', '.join('?' * len(ids))})",
                    list(ids),
                )
        finally:
            conn.close()

    def finish_node(self, node_id, state):
        conn = self._connect()
        try:
            with conn:
                conn.execute(
                    "UPDATE work_jobs SET state = ?, finished_at = NOW(6), heartbeat_at = NOW(6),"
                    " seconds = TIMESTAMPDIFF(MICROSECOND, COALESCE(started_at, created_at), NOW(6)) / 1e6"
                    " WHERE id = ?",
                    (state, node_id),
                )
        finally:
            conn.close()

    def controls(self, chain):
        """
        What an admin asked of the nodes in `chain` (ADM.10): `{"cancel":
        [ids asked to cancel], "pause": bool}`, `pause` also when the
        whole queue is paused.
        """
        conn = self._connect()
        try:
            rows = conn.execute(
                f"SELECT id, control FROM work_jobs WHERE id IN ({', '.join('?' * len(chain))})"
                " AND control IS NOT NULL",
                list(chain),
            ).fetchall()
            conn.rollback()
        finally:
            conn.close()
        return {
            "cancel": [row["id"] for row in rows if row["control"] == "cancel"],
            "pause": any(row["control"] == "pause" for row in rows) or pause.is_paused(),
        }

    def set_state(self, job_id, state):
        conn = self._connect()
        try:
            with conn:
                conn.execute("UPDATE work_jobs SET state = ?, heartbeat_at = NOW(6) WHERE id = ?", (state, job_id))
        finally:
            conn.close()

    def heartbeat(self, job_id, holder):
        conn = self._connect()
        try:
            with conn:
                conn.execute("UPDATE work_jobs SET heartbeat_at = NOW(6) WHERE id = ?", (job_id,))
                conn.execute(
                    "UPDATE work_lease SET heartbeat_at = NOW(6) WHERE id = 1 AND holder = ?", (holder,)
                )
        finally:
            conn.close()

    def add_tasks(self, job_id, tasks):
        conn = self._connect()
        try:
            with conn:
                for task in tasks:
                    cur = conn.execute(
                        "INSERT INTO work_tasks (job_id, kind, task_key, weight, state, created_at)"
                        " VALUES (?, ?, ?, ?, 'queued', NOW(6))",
                        (job_id, task.kind, str(task.key)[:255], float(task.weight)),
                    )
                    task.row_id = cur.lastrowid
                conn.execute(
                    "UPDATE work_jobs SET tasks_queued = tasks_queued + ? WHERE id = ?", (len(tasks), job_id)
                )
        finally:
            conn.close()

    def start_tasks(self, tasks):
        ids = [task.row_id for task in tasks if task.row_id is not None]
        if not ids:
            return
        conn = self._connect()
        try:
            with conn:
                conn.execute(
                    f"UPDATE work_tasks SET state = 'running', started_at = NOW(6)"
                    f" WHERE id IN ({', '.join('?' * len(ids))})",
                    ids,
                )
        finally:
            conn.close()

    def finish_task(self, job_id, task, state, seconds=None, error=None):
        if task.row_id is None:
            return
        conn = self._connect()
        try:
            with conn:
                conn.execute(
                    "UPDATE work_tasks SET state = ?, finished_at = NOW(6), seconds = ?, error = ?"
                    " WHERE id = ?",
                    (state, seconds, error[:4000] if error else None, task.row_id),
                )
                column = {"done": "tasks_done", "failed": "tasks_failed"}.get(state)
                if column:
                    conn.execute(
                        f"UPDATE work_jobs SET {column} = {column} + 1, heartbeat_at = NOW(6) WHERE id = ?",
                        (job_id,),
                    )
        finally:
            conn.close()


def _open_store(control_config):
    if control_config is None:
        return _NoStore()
    try:
        return _ControlStore(control_config)
    except Exception as exc:  # noqa: BLE001 -- no control database: run without the lease
        log.debug(f"Work queue: running without the control database's lease and task rows ({exc}).")
        return _NoStore()


# ---------------------------------------------------------------------
# The job tree (ADM.12)
# ---------------------------------------------------------------------

def new_node_id():
    """A node id: the time it was made and 8 random hex digits."""
    return time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(4)


def end_state(exc_type):
    """The state a node ends in when its block raised `exc_type`
    (`None`: it didn't): Ctrl+C, SIGTERM and `SystemExit` cancel it,
    anything else fails it."""
    if exc_type is None:
        return "done"
    return "cancelled" if issubclass(exc_type, (KeyboardInterrupt, SystemExit)) else "failed"


class Node:
    """
    One node of a job tree, open in this process (`open_node`). Its row
    is written through `store`, which is `_NoStore` (nothing written)
    without a control database.

    Attributes:
        id (str): The row's id (`new_node_id`).
        parent_id (str | None): The node above it, `None` for a root.
        root_id (str): The tree's root (itself for a root).
        chain (list[str]): Its own id and every ancestor's, nearest
            first.
    """

    def __init__(self, kind, title, store, parent_id=None, root_id=None, chain=(), holder=None, workers=0,
                 web_job_id=None, database=None, argv=None):
        self.id = new_node_id()
        self.kind, self.title, self.store = kind, title, store
        self.parent_id = parent_id
        self.root_id = root_id or self.id
        self.chain = [self.id, *chain]
        self.holder = holder or _process_holder()
        self.workers = int(workers)
        self.web_job_id, self.database = web_job_id, database
        self.argv = None if argv is None else [str(arg)[:ARGV_LENGTH] for arg in argv]
        self.closed = False


_open_nodes = []
"""list[Node]: This process's open nodes, innermost last."""

_beat_lock = threading.Lock()
_beat_thread = None


def _process_holder():
    return f"{socket.gethostname()}:{os.getpid()}"


def current_node():
    """The innermost node open in this process, or `None`."""
    return _open_nodes[-1] if _open_nodes else None


def _beat_open_nodes():
    """Keeps every open node's heartbeat fresh, so a reader can tell a
    live run from one that died (`STALE_SECONDS`)."""
    global _beat_thread
    while True:
        time.sleep(HEARTBEAT_SECONDS)
        with _beat_lock:
            nodes = [node for node in _open_nodes if node.store.available]
            if not nodes:
                _beat_thread = None
                return
        by_store = {}
        for node in nodes:
            by_store.setdefault(id(node.store), (node.store, []))[1].append(node.id)
        for store, ids in by_store.values():
            try:
                store.beat_nodes(ids)
            except Exception as exc:  # noqa: BLE001 -- retried next beat
                log.debug(f"Job tree heartbeat failed: {exc}")


def _start_beat():
    global _beat_thread
    with _beat_lock:
        if _beat_thread is None:
            _beat_thread = threading.Thread(target=_beat_open_nodes, name="job-tree-heartbeat", daemon=True)
            _beat_thread.start()


def open_node(kind, title, control_config=None, state="running", workers=0, holder=None, argv=None,
              database=None, web_job_id=None, parent_id=None):
    """
    Opens a job tree node under this process's innermost open node, or,
    for this process's first, under `parent_id` (default: the
    `PARENT_ENV_VAR` node) or as a root. Close it with `close_node`
    (or use `job_node`). Never raises for the bookkeeping: without a
    control database, or when the row can't be written, the node exists
    here only.

    Args:
        kind (str): What sort of node ("plan", "bright-stars", "queue",
            "web-job", "step" ...), the key time measurements group by.
        title (str): What it is doing, for the admin page.
        control_config (MySQLConfig, optional): Used for this process's
            first node; inner nodes use their parent's store.
        state (str): "running", or "waiting" for a queue that waits for
            the lease before it starts.
        workers (int): Its pool's worker processes (0: no pool).
        argv (list[str], optional): The command line that would run it
            again (a root's), without `--mysql-*` options.
        database (str, optional): The galaxy database it writes.
        web_job_id (str, optional): The Generate page job it belongs to.
        parent_id (str, optional): See above.

    Returns:
        Node: The open node.
    """
    parent = current_node()
    if parent is not None and parent.store.available:
        # A run whose node was cancelled from the admin page (ADM.10)
        # starts nothing more.
        try:
            controls = parent.store.controls(parent.chain) or {}
        except Exception as exc:  # noqa: BLE001 -- carry on
            log.debug(f"Job tree: can't read the admin's requests ({exc}).")
            controls = {}
        if controls.get("cancel"):
            raise Cancelled(title)
    if parent is not None:
        node = Node(kind, title, parent.store, parent.id, parent.root_id, parent.chain, holder, workers,
                    web_job_id or parent.web_job_id, database or parent.database, argv)
    else:
        updateguard.refuse_while_updating("A run")      # OPS.53: this process's first node is the run's start
        store = _open_store(control_config)
        parent_id = parent_id or os.environ.get(PARENT_ENV_VAR) or None
        chain = []
        if parent_id and NODE_ID_RE.match(parent_id) and store.available:
            try:
                chain = store.node_chain(parent_id)
            except Exception as exc:  # noqa: BLE001 -- start a tree of its own
                log.debug(f"Job tree: can't read node {parent_id} ({exc}).")
        if chain:
            node = Node(kind, title, store, chain[0][0], chain[0][1], [node_id for node_id, _root in chain],
                        holder, workers, web_job_id, database, argv)
        else:
            node = Node(kind, title, store, holder=holder, workers=workers, web_job_id=web_job_id,
                        database=database, argv=argv)
    try:
        node.store.create_node(node, state)
    except Exception as exc:  # noqa: BLE001 -- the run goes on without its row
        log.debug(f"Job tree: can't record node {node.title!r} ({exc}).")
        node.store = _NoStore()
    with _beat_lock:
        _open_nodes.append(node)
    if node.store.available:
        _start_beat()
    return node


def close_node(node, state):
    """Ends `node` (and any node still open inside it) in `state`
    ("done", "failed" or "cancelled"), recording its end and duration."""
    with _beat_lock:
        inner = []
        if node in _open_nodes:
            index = _open_nodes.index(node)
            inner = _open_nodes[index + 1:]
            del _open_nodes[index:]
    for other in reversed(inner):
        _finish(other, state)
    _finish(node, state)


def _finish(node, state):
    if node.closed:
        return
    node.closed = True
    try:
        node.store.finish_node(node.id, state)
    except Exception as exc:  # noqa: BLE001 -- it goes stale on its own
        log.debug(f"Job tree: can't record the end of node {node.title!r} ({exc}).")


def record_task(node, kind, key, state="done", seconds=None, error=None):
    """
    Writes one finished task row under `node` -- how a node that is not a work queue shows a line of
    its log on the admin page (the cache warm-up's one line per family, PERF.104).

    Args:
        node (Node): An open node.
        kind (str): The task's kind (its column on the page).
        key (str): The line shown (cut to 255 characters).
        state (str): `"done"` or `"failed"`.
        seconds (float, optional): How long it took.
        error (str, optional): What went wrong, for a failed one.
    """
    task = _Task(kind, key, None, None, 1.0, None)
    try:
        node.store.add_tasks(node.id, [task])
        node.store.start_tasks([task])
        node.store.finish_task(node.id, task, state, seconds=seconds, error=error)
    except Exception as exc:  # noqa: BLE001 -- the log line is lost, the run is not
        log.debug(f"Job tree: can't record the task {key!r} ({exc}).")


@contextlib.contextmanager
def job_node(kind, title, control_config=None, **options):
    """`with job_node("bright-stars", "Bright stars"):` -- `open_node`
    for the block, closed in the state the block ends in (`end_state`).
    Yields the `Node`."""
    node = open_node(kind, title, control_config, **options)
    try:
        yield node
    except BaseException as exc:
        close_node(node, end_state(type(exc)))
        raise
    close_node(node, "done")


# ---------------------------------------------------------------------
# Reading job trees (ADM.10, ADM.12)
# ---------------------------------------------------------------------

_TASK_STATES = ("queued", "running", "done", "failed", "cancelled")

# Times come back as Unix times (floats), whatever the server's zone.
_NODE_COLUMNS = ("id, parent_id, root_id, kind, title, holder, state, workers, tasks_total, web_job_id,"
                 " database_name, control, UNIX_TIMESTAMP(created_at) AS created_at,"
                 " UNIX_TIMESTAMP(started_at) AS started_at, UNIX_TIMESTAMP(finished_at) AS finished_at, seconds,"
                 " UNIX_TIMESTAMP(heartbeat_at) AS heartbeat_at,"
                 " heartbeat_at < NOW(6) - INTERVAL ? SECOND AS stale,"
                 " TIMESTAMPDIFF(MICROSECOND, COALESCE(started_at, created_at), NOW(6)) / 1e6 AS age_seconds")


def _attach_argv(conn, nodes):
    """Sets each `_node_dict`'s `argv` to its command line (`work_job_args`),
    or leaves `None` for a node that has none."""
    by_id = {node["id"]: node for node in nodes}
    if not by_id:
        return
    marks = ", ".join("?" * len(by_id))
    for row in conn.execute(
        f"SELECT job_id, value FROM work_job_args WHERE job_id IN ({marks}) ORDER BY job_id, position",
        list(by_id),
    ).fetchall():
        node = by_id[row["job_id"]]
        node["argv"] = (node["argv"] or []) + [row["value"]]


def _float_or_none(value):
    return None if value is None else float(value)


def _node_dict(row):
    node = dict(row)
    for key in ("created_at", "started_at", "finished_at", "heartbeat_at", "seconds", "age_seconds"):
        node[key] = _float_or_none(node.get(key))
    node["stale"] = bool(node["stale"])
    node["live"] = node["state"] in LIVE_STATES and not node["stale"]
    # A run that died without finishing shows as interrupted.
    node["status"] = "interrupted" if node["state"] in LIVE_STATES and node["stale"] else node["state"]
    node["root_id"] = node["root_id"] or node["id"]
    node["argv"] = None
    node["children"] = []
    node["tasks"] = []
    return node


def list_roots(conn, limit=50, offset=0):
    """
    The newest job trees' root rows (`_node_dict`s without their
    children), newest first, and how many roots there are.

    Args:
        conn (Connection): A control database connection.

    Returns:
        tuple[list[dict], int]: The page of roots and the total.
    """
    total = conn.execute("SELECT COUNT(*) AS n FROM work_jobs WHERE parent_id IS NULL").fetchone()["n"]
    rows = conn.execute(
        f"SELECT {_NODE_COLUMNS} FROM work_jobs WHERE parent_id IS NULL"
        " ORDER BY created_at DESC, id DESC LIMIT ? OFFSET ?",
        (STALE_SECONDS, int(limit), int(offset)),
    ).fetchall()
    roots = [_node_dict(row) for row in rows]
    _attach_argv(conn, roots)
    return roots, int(total)


def recorded_step_seconds(conn):
    """
    How long each kind of web-job step took in earlier runs (PERF.33): the
    mean seconds of every finished "step" node, by its title.

    Returns:
        dict: `{title: mean seconds}`.
    """
    rows = conn.execute(
        "SELECT title, AVG(seconds) AS mean FROM work_jobs WHERE kind = 'step' AND state = 'done'"
        " AND seconds IS NOT NULL GROUP BY title").fetchall()
    conn.rollback()
    return {row["title"]: float(row["mean"]) for row in rows}


def load_tree(conn, root_id, max_tasks=200, recorded_task_seconds=None, unstarted_steps=None):
    """
    One whole job tree, every node with its totals added up from its
    children (`_roll_up`).

    Args:
        conn (Connection): A control database connection.
        root_id (str): The root's id (any node's `root_id`).
        max_tasks (int): The most task rows listed per queue node (the
            counts always cover them all); failed tasks come first.
        recorded_task_seconds (callable, optional): `(task_kind, workers)`
            -> the seconds one task of that kind took a worker in earlier
            runs (PERF.32), or `None`; the time left blends it with the
            run's own pace (PERF.33).

    Returns:
        dict | None: The root `_node_dict`, its `children` nested, each
            queue node's `tasks` (`id`, `kind`, `task_key`, `state`,
            `started_at`, `finished_at`, `seconds`, `error`), and on
            every node `totals` (`_roll_up`). `None` for an unknown id.
    """
    rows = conn.execute(
        f"SELECT {_NODE_COLUMNS} FROM work_jobs WHERE id = ? OR root_id = ? ORDER BY created_at, id",
        (STALE_SECONDS, root_id, root_id),
    ).fetchall()
    nodes = {row["id"]: _node_dict(row) for row in rows}
    root = nodes.get(root_id)
    if root is None:
        return None
    _attach_argv(conn, nodes.values())
    for node in nodes.values():
        parent = nodes.get(node["parent_id"])
        if parent is not None and node is not root:
            parent["children"].append(node)
    ids = list(nodes)
    marks = ", ".join("?" * len(ids))
    counts = {}
    for row in conn.execute(
        f"SELECT job_id, state, COUNT(*) AS n, COALESCE(SUM(seconds), 0) AS seconds FROM work_tasks"
        f" WHERE job_id IN ({marks}) GROUP BY job_id, state",
        ids,
    ).fetchall():
        counts.setdefault(row["job_id"], {})[row["state"]] = (int(row["n"]), float(row["seconds"]))
    kinds = {}
    for row in conn.execute(
        f"SELECT job_id, kind, COUNT(*) AS n FROM work_tasks WHERE job_id IN ({marks}) GROUP BY job_id, kind", ids,
    ).fetchall():
        if int(row["n"]) > kinds.get(row["job_id"], (None, 0))[1]:
            kinds[row["job_id"]] = (row["kind"], int(row["n"]))
    for node_id, node in nodes.items():
        node["own_counts"] = counts.get(node_id, {})
        node["recorded_task_seconds"] = None
        if recorded_task_seconds is not None and node_id in kinds:
            node["recorded_task_seconds"] = recorded_task_seconds(kinds[node_id][0], node["workers"] or 1)
        if node["own_counts"] and max_tasks:
            node["tasks"] = [dict(task) for task in conn.execute(
                "SELECT id, kind, task_key, state, UNIX_TIMESTAMP(started_at) AS started_at,"
                " UNIX_TIMESTAMP(finished_at) AS finished_at, seconds, error FROM work_tasks"
                " WHERE job_id = ? ORDER BY state = 'failed' DESC, id LIMIT ?",
                (node_id, int(max_tasks)),
            ).fetchall()]
            for task in node["tasks"]:
                for key in ("started_at", "finished_at", "seconds"):
                    task[key] = _float_or_none(task[key])
    conn.rollback()
    if unstarted_steps is not None and root["live"]:
        root["unstarted_step_seconds"] = list(unstarted_steps(root))
    _roll_up(root)
    return root


def _roll_up(node):
    """
    Fills `node["totals"]` from its own tasks and its children's totals:
    `tasks` (planned), `queued` (planned, not finished or running yet),
    `running`, `done`, `failed`, `cancelled`, `work_seconds` (time
    workers spent on finished tasks), `started_at`/`finished_at` (the
    earliest start and latest end below it) and `eta_seconds` (time left
    at the measured pace, `None` when there's nothing to measure by or
    the node isn't live).
    """
    own = node.pop("own_counts", {})
    totals = {state: own.get(state, (0, 0.0))[0] for state in _TASK_STATES}
    totals["work_seconds"] = sum(seconds for _n, seconds in own.values())
    recorded = sum(totals[state] for state in _TASK_STATES)
    planned = max(node["tasks_total"] or 0, recorded)
    totals["queued"] += planned - recorded
    totals["tasks"] = planned
    unstarted = node.pop("unstarted_step_seconds", None) or []
    eta = None
    if node["live"] and planned:
        remaining = planned - totals["done"] - totals["failed"] - totals["cancelled"]
        live_rate = None
        if totals["done"] and own.get("done", (0, 0.0))[1] > 0:
            live_rate = totals["done"] / own["done"][1]       # tasks a worker second
        recorded = node.pop("recorded_task_seconds", None)
        # PERF.33: the pace measured so far blended with the recorded one.
        rate = progress_rate.blended_rate(live_rate, totals["done"], 1.0 / recorded if recorded else None)
        if rate:
            eta = remaining / rate / max(node["workers"] or 1, 1)
    starts = [node["started_at"]] if node["started_at"] else []
    ends = [node["finished_at"]] if node["finished_at"] else []
    child_etas = []
    for child in node["children"]:
        _roll_up(child)
        sub = child["totals"]
        for key in _TASK_STATES + ("tasks", "work_seconds"):
            totals[key] += sub[key]
        if sub["started_at"]:
            starts.append(sub["started_at"])
        if sub["finished_at"]:
            ends.append(sub["finished_at"])
        if sub["eta_seconds"] is not None:
            child_etas.append(sub["eta_seconds"])
    if child_etas:
        eta = (eta or 0.0) + sum(child_etas)
    # PERF.33: each step not started yet adds what earlier runs of it took.
    known = [seconds for seconds in unstarted if seconds is not None]
    if known:
        eta = (eta or 0.0) + sum(known)
    totals["eta_partial"] = len(known) < len(unstarted)
    totals["eta_seconds"] = eta if node["live"] else None
    totals["started_at"] = min(starts) if starts else None
    totals["finished_at"] = max(ends) if ends and not node["live"] else None
    node["totals"] = totals


def timing_by_kind(conn):
    """
    Measured times per kind of node and of task, from every finished
    one still kept (`KEEP_DAYS`): what the time estimates can use.

    Returns:
        dict: `nodes` and `tasks`, each `{kind: {"count", "mean_seconds",
            "min_seconds", "max_seconds"}}`.
    """
    result = {}
    for name, sql in (
        ("nodes", "SELECT kind, COUNT(*) AS n, AVG(seconds) AS mean, MIN(seconds) AS lo, MAX(seconds) AS hi"
                  " FROM work_jobs WHERE state = 'done' AND seconds IS NOT NULL GROUP BY kind"),
        ("tasks", "SELECT kind, COUNT(*) AS n, AVG(seconds) AS mean, MIN(seconds) AS lo, MAX(seconds) AS hi"
                  " FROM work_tasks WHERE state = 'done' AND seconds IS NOT NULL GROUP BY kind"),
    ):
        result[name] = {
            row["kind"]: {"count": int(row["n"]), "mean_seconds": float(row["mean"]),
                          "min_seconds": float(row["lo"]), "max_seconds": float(row["hi"])}
            for row in conn.execute(sql).fetchall()
        }
    conn.rollback()
    return result


# ---------------------------------------------------------------------
# Managing the queue from the admin page (ADM.10)
# ---------------------------------------------------------------------

ACTIONS = ("pause", "resume", "cancel")


def queue_status(conn):
    """
    The queue's state for the admin page.

    Returns:
        dict: `paused`, `paused_by`, `paused_at`; `holder`, `job_id`,
            `heartbeat_age_s` and `stale` of the lease (holder `None`
            when free); `workers_active` (worker processes of every live,
            running work queue, 1 for a one-worker run) and `runs_active`
            (live root runs).
    """
    conn.execute("INSERT IGNORE INTO work_lease (id) VALUES (1)")
    conn.commit()
    lease = conn.execute(
        "SELECT holder, job_id, TIMESTAMPDIFF(MICROSECOND, heartbeat_at, NOW(6)) / 1e6 AS age FROM work_lease WHERE id = 1"
    ).fetchone()
    paused = pause.state()
    live = conn.execute(
        "SELECT kind, parent_id, state, workers FROM work_jobs"
        " WHERE state IN ('waiting', 'running', 'paused') AND heartbeat_at >= NOW(6) - INTERVAL ? SECOND",
        (STALE_SECONDS,),
    ).fetchall()
    conn.rollback()
    age = float(lease["age"]) if lease["age"] is not None else None
    return {
        "paused": paused["paused"],
        "paused_by": paused["by"],
        "paused_at": _float_or_none(paused["at"]),
        "holder": lease["holder"],
        "job_id": lease["job_id"],
        "heartbeat_age_s": age,
        "stale": lease["holder"] is not None and (age is None or age > STALE_SECONDS),
        "workers_active": sum(max(int(row["workers"]), 1) for row in live
                              if row["kind"] == "queue" and row["state"] == "running"),
        "runs_active": sum(1 for row in live if row["parent_id"] is None),
    }


def get_node(conn, node_id):
    """One node's row (`_node_dict`, without children), or `None`."""
    row = conn.execute(f"SELECT {_NODE_COLUMNS} FROM work_jobs WHERE id = ?", (STALE_SECONDS, node_id)).fetchone()
    conn.rollback()
    return _node_dict(row) if row else None


def request_control(conn, node_id, action):
    """
    Asks a live node and its subtree to pause, resume or cancel (ADM.10);
    the run reads it between its units. Resume clears a
    pause; a cancel stays.

    Returns:
        bool: Whether a live node took the request.
    """
    if action not in ACTIONS:
        raise ValueError(f"unknown action {action!r}")
    value = {"pause": "pause", "resume": None, "cancel": "cancel"}[action]
    with conn:
        if action == "resume":
            cur = conn.execute(
                "UPDATE work_jobs SET control = NULL WHERE id = ? AND control = 'pause'"
                " AND state IN ('waiting', 'running', 'paused')",
                (node_id,),
            )
            # A paused descendant resumes with it.
            node = conn.execute("SELECT root_id FROM work_jobs WHERE id = ?", (node_id,)).fetchone()
            if node is not None:
                ids = _subtree_ids(conn, node["root_id"] or node_id, node_id)
                if len(ids) > 1:
                    conn.execute(
                        f"UPDATE work_jobs SET control = NULL WHERE control = 'pause'"
                        f" AND id IN ({', '.join('?' * len(ids))})",
                        ids,
                    )
                    return True
        else:
            cur = conn.execute(
                "UPDATE work_jobs SET control = ? WHERE id = ? AND state IN ('waiting', 'running', 'paused')"
                " AND heartbeat_at >= NOW(6) - INTERVAL ? SECOND AND (control IS NULL OR control <> 'cancel')",
                (value, node_id, STALE_SECONDS),
            )
        return cur.rowcount > 0


def _subtree_ids(conn, root_id, node_id):
    """`node_id` and every node below it in its tree."""
    rows = conn.execute("SELECT id, parent_id FROM work_jobs WHERE id = ? OR root_id = ?",
                        (root_id, root_id)).fetchall()
    children = {}
    for row in rows:
        children.setdefault(row["parent_id"], []).append(row["id"])
    found, todo = [], [node_id]
    while todo:
        current = todo.pop()
        found.append(current)
        todo.extend(children.get(current, []))
    return found


def clear_stale_lease(conn):
    """Frees the lease when its holder stopped refreshing it
    (`STALE_SECONDS`), and marks that run's live nodes cancelled. Returns
    the holder it cleared, or `None` (free, or still live)."""
    with conn:
        row = conn.execute(
            "SELECT holder, job_id FROM work_lease WHERE id = 1 AND holder IS NOT NULL"
            " AND (heartbeat_at IS NULL OR heartbeat_at < NOW(6) - INTERVAL ? SECOND) FOR UPDATE",
            (STALE_SECONDS,),
        ).fetchone()
        if row is None:
            return None
        conn.execute("UPDATE work_lease SET holder = NULL, job_id = NULL, heartbeat_at = NULL WHERE id = 1")
        if row["job_id"]:
            conn.execute(
                "UPDATE work_tasks SET state = 'cancelled', finished_at = NOW(6)"
                " WHERE job_id = ? AND state IN ('queued', 'running')",
                (row["job_id"],),
            )
            conn.execute(
                "UPDATE work_jobs SET state = 'cancelled', finished_at = NOW(6),"
                " seconds = TIMESTAMPDIFF(MICROSECOND, COALESCE(started_at, created_at), heartbeat_at) / 1e6"
                " WHERE id = ? AND state IN ('waiting', 'running', 'paused')"
                " AND heartbeat_at < NOW(6) - INTERVAL ? SECOND",
                (row["job_id"], STALE_SECONDS),
            )
        return row["holder"]


def delete_tree(conn, root_id):
    """Deletes a finished (or dead) job tree with its tasks. Returns
    whether it did: a live tree, or a node that isn't a root, stays."""
    with conn:
        cur = conn.execute(
            "DELETE FROM work_jobs WHERE id = ? AND parent_id IS NULL"
            " AND (state NOT IN ('waiting', 'running', 'paused') OR heartbeat_at < NOW(6) - INTERVAL ? SECOND)",
            (root_id, STALE_SECONDS),
        )
        return cur.rowcount > 0


def clear_finished(conn):
    """Deletes every finished (or dead) job tree with its tasks. Returns
    how many there were; live trees stay."""
    with conn:
        cur = conn.execute(
            "DELETE FROM work_jobs WHERE parent_id IS NULL"
            " AND (state NOT IN ('waiting', 'running', 'paused') OR heartbeat_at < NOW(6) - INTERVAL ? SECOND)",
            (STALE_SECONDS,),
        )
        return cur.rowcount


