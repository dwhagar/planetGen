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
order, exactly as generation always did.

Workers are separate processes (`multiprocessing`'s spawn start method,
so Linux, macOS and Windows behave the same) because generation is pure
Python and threads would share one core. Each task seeds `random` from
the run's seed and its own key (`task_seed`), so what a task generates
doesn't depend on how many workers there are or which finished first.
"""

import concurrent.futures
import hashlib
import ipaddress
import json
import math
import multiprocessing
import os
import random
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

    def create_job(self, job_id, title, holder, workers):
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
            conn.execute("SELECT 1 FROM work_lease LIMIT 1").fetchall()
        finally:
            conn.close()

    def _connect(self):
        return self._db.get_control_connection(self.config)

    def create_job(self, job_id, title, holder, workers):
        conn = self._connect()
        try:
            with conn:
                conn.execute(
                    "DELETE FROM work_jobs WHERE created_at < NOW(6) - INTERVAL ? DAY", (KEEP_DAYS,)
                )
                conn.execute(
                    "INSERT INTO work_jobs (id, title, holder, state, workers, created_at, heartbeat_at)"
                    " VALUES (?, ?, ?, 'waiting', ?, NOW(6), NOW(6))",
                    (job_id, title[:255], holder, workers),
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
                    "SELECT id FROM work_jobs WHERE state IN ('waiting', 'running') AND id <> ?"
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
                        "UPDATE work_jobs SET state = 'cancelled', finished_at = NOW(6) WHERE id = ?",
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
                    "UPDATE work_jobs SET state = ?, finished_at = NOW(6), heartbeat_at = NOW(6) WHERE id = ?",
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
        if not self.parallel:
            return self
        self._store = _open_store(self.control_config)
        self.job_id = time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(4)
        self._store.create_job(self.job_id, self.title, self.holder, self.workers)
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
            started = time.monotonic()
            result = fn(payload)
            self.finished += 1
            if on_done:
                on_done(result, time.monotonic() - started, weight)
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
            return False
        state = "done"
        try:
            if exc_type is None:
                self.drain()
            else:
                state = "cancelled" if issubclass(exc_type, (KeyboardInterrupt, SystemExit)) else "failed"
                self._waiting.clear()
        except BaseException:
            state = "failed"
            raise
        finally:
            if self._executor is not None:
                self._executor.shutdown(wait=True, cancel_futures=True)
            self._stop.set()
            self._restore_sigterm()
            try:
                self._store.finish_job(self.job_id, self.holder, state)
            except Exception as store_exc:  # noqa: BLE001 -- the lease goes stale on its own
                log.debug(f"Work queue: could not record the end of job {self.job_id}: {store_exc}")
        return False
