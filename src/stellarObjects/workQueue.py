# stellarObjects/workQueue.py

"""
The generation work queue (PERF.8): runs independent units of work (a
sector to fill, a layer of bright stars to scatter) on a pool of worker
processes, at most about 80% of the machine's cores and at a lower
scheduling priority than everything else on it.

A `WorkQueue` belongs to the run that queues the work -- a `generate.py`
run, started from the command line or by the Generate page's job
runner. That run is the supervisor for as long as it has work: there is
no daemon, service or scheduler to install, and nothing runs once the
work is done. The run hands tasks to the pool as fast as the pool takes
them, and each finished task comes back to it as a signal (`on_done`,
with the time the task took) it uses for its progress bars and ETA.

    with WorkQueue("Sectors (ring 12)", workers=6, control_config=cfg) as queue:
        for address in addresses:
            queue.submit("sector", key, fill_sector, payload, on_done=report)

Only one run's pool uses the machine at a time. The control database's
`work_lease` row says which (see `control_schema.sql`, v5): a run holds
it while its pool runs, refreshes it every `HEARTBEAT_SECONDS`, and frees
it when it finishes; a second run started meanwhile waits for it, so two
runs never use 160% of the CPUs between them. A lease not refreshed for
`STALE_SECONDS` belongs to a run that died, and the next run takes it.
`work_jobs`/`work_tasks` keep one row per run and per task (state,
timings, a short result), which is what the Generate page and later
runs can read back. Without the control database (not created yet, or
no grant on it) the queue still runs, just without the lease or the
rows.

One worker means no pool at all: every task runs in this process, in
order, exactly as generation always did (still recorded in the control
database, without the lease, when the run passes one).

Every job is a tree (ADM.12, control schema v7). `job_node` opens a
node under the one this process has open (a `generate.py` run is the
root; its phases, such as the bright stars, are nodes under it), and a
`WorkQueue` is a node whose leaves are its tasks. A process started by
the Generate page's job runner hangs its root under the runner's step
node (`PARENT_ENV_VAR`). Each node keeps its state, start, end and
duration; `load_tree` reads a whole tree back with every parent's
totals, timings and ETA added up from its children.

Workers are separate processes (`multiprocessing`'s spawn start method,
so Linux, macOS and Windows behave the same) because generation is pure
Python and threads would share one core. Each task seeds `random` from
the run's seed and its own key (`task_seed`), so what a task generates
doesn't depend on how many workers there are or which finished first.
"""

import concurrent.futures
import contextlib
import hashlib
import ipaddress
import json
import math
import multiprocessing
import os
import random
import re
import secrets
import signal
import socket
import threading
import time

from stellarObjects import log

CPU_SHARE = 0.8
"""float: The share of the machine's cores the pool may use."""

WORKERS_ENV_VAR = "PLANETGEN_WORKERS"
"""str: Environment variable overriding the worker count (`--workers`
on the command line wins over it); 0 or "auto" means `CPU_SHARE` of the
cores."""

WORKER_NICENESS = 10
"""int: How much `os.nice` lowers a worker's priority on Linux and
macOS. Windows workers run at `BELOW_NORMAL_PRIORITY_CLASS`."""

HEARTBEAT_SECONDS = 5
STALE_SECONDS = 30
WAIT_POLL_SECONDS = 2
KEEP_DAYS = 7

_BELOW_NORMAL_PRIORITY_CLASS = 0x4000

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


def lower_priority():
    """Lowers this process's scheduling priority so everything else on
    the machine comes first. Never raises."""
    try:
        if os.name == "nt":
            import ctypes

            kernel32 = ctypes.windll.kernel32
            kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), _BELOW_NORMAL_PRIORITY_CLASS)
        else:
            os.nice(WORKER_NICENESS)
    except Exception:  # noqa: BLE001 -- a normal-priority worker is still a worker
        pass


def task_seed(run_seed, key):
    """The `random` seed for one task: the run's seed and the task's key,
    hashed, so it doesn't depend on worker count or finish order."""
    digest = hashlib.sha256(f"{run_seed}:{key}".encode("utf-8")).digest()
    return int.from_bytes(digest[:16], "big")


def _worker_init(log_level, debug_file):
    lower_priority()
    os.environ.pop("PLANETGEN_PROGRESS_FILE", None)
    log.set_component("worker")
    try:
        log.configure(log_level, debug_file=debug_file, console=False)
    except OSError:
        log.configure(log_level, console=False)


def _run_task(fn, payload, seed):
    random.seed(seed)
    started = time.monotonic()
    result = fn(payload)
    return result, time.monotonic() - started


class _Task:
    __slots__ = ("kind", "key", "fn", "payload", "weight", "on_done", "row_id", "future", "seed")

    def __init__(self, kind, key, fn, payload, weight, on_done):
        self.kind, self.key, self.fn, self.payload = kind, key, fn, payload
        self.weight, self.on_done = weight, on_done
        self.row_id = self.future = self.seed = None


class _NoStore:
    """Bookkeeping when there is no control database: nothing kept, no
    lease, so the run never waits."""

    available = False

    def create_node(self, node, state):
        pass

    def node_chain(self, node_id):
        return []

    def beat_nodes(self, ids):
        pass

    def set_total(self, job_id, count):
        pass

    def finish_node(self, node_id, state):
        pass

    def take_lease(self, job_id, holder):
        return True, None

    def heartbeat(self, job_id, holder):
        pass

    def add_tasks(self, job_id, tasks):
        pass

    def start_tasks(self, tasks):
        pass

    def finish_task(self, job_id, task, state, seconds=None, result=None, error=None):
        pass

    def finish_job(self, job_id, holder, state):
        pass


class _ControlStore:
    """`work_jobs`/`work_tasks`/`work_lease` in the control database."""

    def __init__(self, config):
        from stellarObjects import _db

        self._db = _db
        self.config = config
        conn = self._connect()
        try:
            # v7's columns too: before `update.sh` brings the control
            # schema up to date, runs go on without the rows and lease.
            conn.execute("SELECT 1 FROM work_lease LIMIT 1").fetchall()
            conn.execute("SELECT parent_id, control FROM work_jobs LIMIT 1").fetchall()
            conn.execute("SELECT paused FROM work_lease LIMIT 1").fetchall()
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
                    " web_job_id, database_name, argv, created_at, started_at, heartbeat_at)"
                    " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NOW(6), "
                    + ("NOW(6)" if state == "running" else "NULL") + ", NOW(6))",
                    (node.id, node.parent_id, node.root_id, node.kind[:32], node.title[:255], node.holder, state,
                     node.workers, node.web_job_id, node.database, node.argv),
                )
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

    def set_total(self, job_id, count):
        conn = self._connect()
        try:
            with conn:
                conn.execute("UPDATE work_jobs SET tasks_total = COALESCE(tasks_total, 0) + ? WHERE id = ?",
                             (int(count), job_id))
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

    def take_lease(self, job_id, holder):
        """Takes the lease when it's free, already this run's, or stale.
        Returns `(taken, current_holder_row)`."""
        conn = self._connect()
        try:
            with conn:
                conn.execute("INSERT IGNORE INTO work_lease (id) VALUES (1)")
                row = conn.execute(
                    "SELECT holder, job_id, heartbeat_at < NOW(6) - INTERVAL ? SECOND AS stale"
                    " FROM work_lease WHERE id = 1 FOR UPDATE",
                    (STALE_SECONDS,),
                ).fetchone()
                if row["holder"] not in (None, holder) and not row["stale"]:
                    conn.execute("UPDATE work_jobs SET heartbeat_at = NOW(6) WHERE id = ?", (job_id,))
                    return False, dict(row)
                conn.execute(
                    "UPDATE work_lease SET holder = ?, job_id = ?, heartbeat_at = NOW(6) WHERE id = 1",
                    (holder, job_id),
                )
                # Runs that died without finishing: their unfinished tasks
                # were rolled back with them.
                dead = [r["id"] for r in conn.execute(
                    "SELECT id FROM work_jobs WHERE state IN ('waiting', 'running', 'paused') AND id <> ?"
                    " AND heartbeat_at < NOW(6) - INTERVAL ? SECOND",
                    (job_id, STALE_SECONDS),
                ).fetchall()]
                for dead_id in dead:
                    conn.execute(
                        "UPDATE work_tasks SET state = 'cancelled', finished_at = NOW(6)"
                        " WHERE job_id = ? AND state IN ('queued', 'running')",
                        (dead_id,),
                    )
                    conn.execute(
                        "UPDATE work_jobs SET state = 'cancelled', finished_at = NOW(6),"
                        " seconds = TIMESTAMPDIFF(MICROSECOND, COALESCE(started_at, created_at), heartbeat_at) / 1e6"
                        " WHERE id = ?",
                        (dead_id,),
                    )
                conn.execute(
                    "UPDATE work_jobs SET state = 'running', started_at = NOW(6), heartbeat_at = NOW(6)"
                    " WHERE id = ?",
                    (job_id,),
                )
                return True, None
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

    def finish_task(self, job_id, task, state, seconds=None, result=None, error=None):
        if task.row_id is None:
            return
        summary = None
        if result is not None:
            try:
                summary = json.dumps(result, default=str)[:4000]
            except (TypeError, ValueError):
                summary = None
        conn = self._connect()
        try:
            with conn:
                conn.execute(
                    "UPDATE work_tasks SET state = ?, finished_at = NOW(6), seconds = ?, result = ?, error = ?"
                    " WHERE id = ?",
                    (state, seconds, summary, error[:4000] if error else None, task.row_id),
                )
                column = {"done": "tasks_done", "failed": "tasks_failed"}.get(state)
                if column:
                    conn.execute(
                        f"UPDATE work_jobs SET {column} = {column} + 1, heartbeat_at = NOW(6) WHERE id = ?",
                        (job_id,),
                    )
        finally:
            conn.close()

    def finish_job(self, job_id, holder, state):
        conn = self._connect()
        try:
            with conn:
                conn.execute(
                    "UPDATE work_tasks SET state = 'cancelled', finished_at = NOW(6)"
                    " WHERE job_id = ? AND state IN ('queued', 'running')",
                    (job_id,),
                )
                conn.execute(
                    "UPDATE work_jobs SET state = ?, finished_at = NOW(6), heartbeat_at = NOW(6),"
                    " seconds = TIMESTAMPDIFF(MICROSECOND, COALESCE(started_at, created_at), NOW(6)) / 1e6"
                    " WHERE id = ?",
                    (state, job_id),
                )
                conn.execute(
                    "UPDATE work_lease SET holder = NULL, job_id = NULL, heartbeat_at = NULL"
                    " WHERE id = 1 AND holder = ?",
                    (holder,),
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
        self.argv = None if argv is None else json.dumps(list(argv))[:60000]
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
    if parent is not None:
        node = Node(kind, title, parent.store, parent.id, parent.root_id, parent.chain, holder, workers,
                    web_job_id or parent.web_job_id, database or parent.database, argv)
    else:
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

_NODE_COLUMNS = ("id, parent_id, root_id, kind, title, holder, state, workers, tasks_total, web_job_id,"
                 " database_name, argv, control, created_at, started_at, finished_at, seconds, heartbeat_at,"
                 " heartbeat_at < NOW(6) - INTERVAL ? SECOND AS stale,"
                 " TIMESTAMPDIFF(MICROSECOND, COALESCE(started_at, created_at), NOW(6)) / 1e6 AS age_seconds")


def _node_dict(row):
    node = dict(row)
    node["stale"] = bool(node["stale"])
    node["live"] = node["state"] in LIVE_STATES and not node["stale"]
    # A run that died without finishing shows as interrupted.
    node["status"] = "interrupted" if node["state"] in LIVE_STATES and node["stale"] else node["state"]
    node["root_id"] = node["root_id"] or node["id"]
    try:
        node["argv"] = json.loads(node["argv"]) if node["argv"] else None
    except ValueError:
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
    return [_node_dict(row) for row in rows], int(total)


def load_tree(conn, root_id, max_tasks=200):
    """
    One whole job tree, every node with its totals added up from its
    children (`_roll_up`).

    Args:
        conn (Connection): A control database connection.
        root_id (str): The root's id (any node's `root_id`).
        max_tasks (int): The most task rows listed per queue node (the
            counts always cover them all); failed tasks come first.

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
    for node_id, node in nodes.items():
        node["own_counts"] = counts.get(node_id, {})
        if node["own_counts"]:
            node["tasks"] = [dict(task) for task in conn.execute(
                "SELECT id, kind, task_key, state, started_at, finished_at, seconds, error FROM work_tasks"
                " WHERE job_id = ? ORDER BY state = 'failed' DESC, id LIMIT ?",
                (node_id, int(max_tasks)),
            ).fetchall()]
    conn.rollback()
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
    eta = None
    if node["live"] and planned:
        remaining = planned - totals["done"] - totals["failed"] - totals["cancelled"]
        if totals["done"]:
            per_task = own.get("done", (0, 0.0))[1] / totals["done"]
            eta = remaining * per_task / max(node["workers"] or 1, 1)
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


class WorkQueue:
    """
    One run's tasks and the worker pool that runs them; see this module's
    docstring. Use as a context manager: entering takes the lease
    (waiting for another run's pool to finish first), leaving waits for
    every submitted task, then frees the lease.

    Args:
        title (str): What the run is doing, for `work_jobs` and the
            "waiting" message.
        workers (int): Worker processes (`worker_count`); 1 runs every
            task in this process, in order, with no pool and no lease.
        control_config (MySQLConfig, optional): The control database
            holding the lease and task rows; `None` runs without them.
        run_seed (int, optional): Seeds every task (`task_seed`); drawn
            from `random` when not given.
        log_level (str): Workers' log severity (`log.NORMAL` etc.); they
            never write to the console, only the debug log/file.
        debug_file (str, optional): `--debug FILE`, appended to by the
            workers too.
        on_wait (callable, optional): Called once with the lease
            holder's row when this run has to wait for another one.
    """

    def __init__(self, title, workers=1, control_config=None, run_seed=None, log_level=log.NORMAL,
                 debug_file=None, on_wait=None):
        self.title = title
        self.workers = max(1, int(workers))
        self.control_config = control_config
        self.run_seed = random.getrandbits(128) if run_seed is None else run_seed
        self.log_level = log_level
        self.debug_file = debug_file
        self.on_wait = on_wait
        self.job_id = None
        self.node = None
        self.holder = f"{socket.gethostname()}:{os.getpid()}:{secrets.token_hex(4)}"
        self._store = _NoStore()
        self._executor = None
        self._waiting = []
        self._running = set()
        self._stop = threading.Event()
        self._beat = None
        self._failure = None
        self._old_sigterm = None
        self.submitted = 0
        self.finished = 0

    @property
    def parallel(self):
        return self.workers > 1

    def __enter__(self):
        self.node = open_node("queue", self.title, self.control_config,
                              state="waiting" if self.parallel else "running", workers=self.workers,
                              holder=self.holder)
        self._store = self.node.store
        self.job_id = self.node.id
        if not self.parallel:
            return self
        waited = False
        while True:
            taken, current = self._store.take_lease(self.job_id, self.holder)
            if taken:
                break
            if not waited:
                waited = True
                log.normal(f"Waiting for another generation run to finish ({current['holder']}, job "
                           f"{current['job_id']}): only one run's workers use the machine at a time.")
                if self.on_wait:
                    self.on_wait(current)
            time.sleep(WAIT_POLL_SECONDS)
        self._catch_sigterm()
        self._beat = threading.Thread(target=self._heartbeat, name="work-queue-heartbeat", daemon=True)
        self._beat.start()
        self._executor = concurrent.futures.ProcessPoolExecutor(
            max_workers=self.workers,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_worker_init,
            initargs=(self.log_level, self.debug_file),
        )
        return self

    def expect(self, count):
        """Says `count` more tasks are coming, so the job tree can show
        how many are still to go and an ETA before they're queued."""
        if count and self.node is not None:
            try:
                self._store.set_total(self.job_id, count)
            except Exception as exc:  # noqa: BLE001 -- bookkeeping only
                log.debug(f"Work queue: could not record the task count of job {self.job_id}: {exc}")

    def _catch_sigterm(self):
        """A SIGTERM (Cancel on the Generate page, `kill`) ends the run
        through `__exit__`, so the lease is freed and the job marked
        cancelled at once rather than after `STALE_SECONDS`."""
        self._old_sigterm = None
        if threading.current_thread() is not threading.main_thread() or not hasattr(signal, "SIGTERM"):
            return

        def _on_sigterm(signum, _frame):
            raise SystemExit(128 + signum)

        try:
            self._old_sigterm = signal.signal(signal.SIGTERM, _on_sigterm)
        except (ValueError, OSError):
            self._old_sigterm = None

    def _restore_sigterm(self):
        if self._old_sigterm is not None:
            try:
                signal.signal(signal.SIGTERM, self._old_sigterm)
            except (ValueError, OSError):
                pass
            self._old_sigterm = None

    def _heartbeat(self):
        while not self._stop.wait(HEARTBEAT_SECONDS):
            try:
                self._store.heartbeat(self.job_id, self.holder)
            except Exception as exc:  # noqa: BLE001 -- retried next beat
                log.debug(f"Work queue heartbeat failed: {exc}")

    def submit(self, kind, key, fn, payload, weight=1.0, on_done=None):
        """
        Queues one task. `fn(payload)` runs in a worker (`fn` must be a
        module-level function, and `payload` picklable); `on_done(result,
        seconds, weight)` then runs in this process, from `submit` or
        `drain`, never at the same time as another `on_done`. With one
        worker, both run right here before `submit` returns.

        Raises:
            Exception: The first task's own exception, once every task
                already running has finished (the rest are cancelled).
        """
        task = _Task(kind, key, fn, payload, weight, on_done)
        self.submitted += 1
        if not self.parallel:
            self._run_here(task)
            return
        self._raise_failure()
        task.seed = task_seed(self.run_seed, key)
        self._waiting.append(task)
        # Keep every worker busy with one more task ready behind it;
        # anything beyond that waits here, not in the pool.
        while len(self._waiting) + len(self._running) > 2 * self.workers:
            self._dispatch()
            self._collect(block=True)
        self._dispatch()
        self._collect(block=False)

    def _run_here(self, task):
        """One worker: runs `task` in this process, recording its row
        when there's a control database."""
        recorded = self._store.available and self.node is not None
        if recorded:
            self._store.add_tasks(self.job_id, [task])
            self._store.start_tasks([task])
        started = time.monotonic()
        try:
            result = task.fn(task.payload)
        except BaseException as exc:
            if recorded:
                state = "failed" if isinstance(exc, Exception) else "cancelled"
                self._store.finish_task(self.job_id, task, state, error=f"{type(exc).__name__}: {exc}")
            raise
        seconds = time.monotonic() - started
        self.finished += 1
        if recorded:
            self._store.finish_task(self.job_id, task, "done", seconds=round(seconds, 3), result=result)
        if task.on_done:
            task.on_done(result, seconds, task.weight)

    def _dispatch(self):
        room = 2 * self.workers - len(self._running)
        batch, self._waiting = self._waiting[:max(room, 0)], self._waiting[max(room, 0):]
        if not batch:
            return
        self._store.add_tasks(self.job_id, batch)
        self._store.start_tasks(batch)
        for task in batch:
            task.future = self._executor.submit(_run_task, task.fn, task.payload, task.seed)
            self._running.add(task)

    def _collect(self, block):
        if not self._running:
            return
        by_future = {task.future: task for task in self._running}
        done, _pending = concurrent.futures.wait(
            by_future, timeout=None if block else 0, return_when=concurrent.futures.FIRST_COMPLETED,
        )
        for future in done:
            task = by_future[future]
            self._running.discard(task)
            self.finished += 1
            try:
                result, seconds = future.result()
            except BaseException as exc:  # noqa: BLE001 -- re-raised by _raise_failure
                self._store.finish_task(self.job_id, task, "failed", error=f"{type(exc).__name__}: {exc}")
                if self._failure is None:
                    self._failure = exc
                    self._waiting.clear()
                continue
            self._store.finish_task(self.job_id, task, "done", seconds=round(seconds, 3), result=result)
            if task.on_done:
                task.on_done(result, seconds, task.weight)

    def _raise_failure(self):
        """After a task failed: lets the tasks already running finish
        (each one's sector is its own transaction), then raises the
        first failure."""
        if self._failure is None:
            return
        failure = self._failure
        self._waiting.clear()
        while self._running:
            self._collect(block=True)
        self._failure = None
        raise failure

    def drain(self):
        """Waits for every submitted task (running each `on_done`)."""
        if not self.parallel:
            return
        while self._waiting or self._running:
            self._raise_failure()
            self._dispatch()
            self._collect(block=True)
        self._raise_failure()

    def __exit__(self, exc_type, exc, tb):
        if not self.parallel:
            self._close(end_state(exc_type))
            return False
        state = "done"
        try:
            if exc_type is None:
                self.drain()
            else:
                state = end_state(exc_type)
                self._waiting.clear()
        except BaseException as raised:
            state = end_state(type(raised))
            raise
        finally:
            if self._executor is not None:
                self._executor.shutdown(wait=True, cancel_futures=True)
            self._stop.set()
            self._restore_sigterm()
            self._close(state)
        return False

    def _close(self, state):
        """Records the job's end (and frees the lease) and closes its
        node."""
        if self.node is None:
            return
        try:
            self._store.finish_job(self.job_id, self.holder, state)
        except Exception as store_exc:  # noqa: BLE001 -- the lease goes stale on its own
            log.debug(f"Work queue: could not record the end of job {self.job_id}: {store_exc}")
        self.node.closed = True
        close_node(self.node, state)
