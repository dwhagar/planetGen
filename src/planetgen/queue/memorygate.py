# planetgen/queue/memorygate.py

"""
One memory and load gate for every job of the work queue (PERF.254).

Before PERF.254 the size of a worker was a number typed into `workers.STAGE_WORKER_MEMORY_MB`, measured once on one galaxy, and
nothing read the host while a run went. A stage that grew past its figure was reported after it ended, a command with no figure
was uncapped, and the cut of a stage was outside the table altogether: the full-galaxy runs of 2026-10-31 and 2026-11-08 were
killed by the kernel that way. This is the gate that reads the host instead.

**Every job holds a claim.** A job is an RQ job on any queue and lane (a run's units, the cut of a stage, the API's and the
pages' jobs, the Generate page's, the node prefetch). Before a worker prepares a job it takes a *claim* on the gate
(`acquire`): the memory the job is expected to hold, in MB, the cores it keeps busy, and the lane it is on (`claims.claim_for`
works the three out from the job). The claim is a Redis hash entry under the worker's name, with an expiry that the worker's
heartbeat renews, so a worker killed with `kill -9` gives its claim up by itself (`holders`).

**The rule.** A job is admitted when its claim fits what the host has left:

    free = MemAvailable - reserve - kept for the reserved places - the part of every other claim not yet in memory

where a claim's part not yet in memory is `max(0, claim - what its worker and horse hold now)` (the proportional set size of the
worker and every process under it, `psutil`), so a job that has grown to its claim is counted by what the host already lost to
it, not twice. The job's own worker is credited with what it already holds (an idle worker with the star tables is 1.4 GB of the
2.1 GB a scatter claims). Cores follow the same rule: the claims' cores are held to `CPU_SHARE` of the host's, and an ordinary job
waits while the one-minute load per core is over `LOAD_LIMIT`. The places reserved for the interactive and long lanes
(`workers.reserved_places`) keep their memory against the bulk and the maintenance lanes, so a visitor's search is never queued
behind a generation run's units. A job that finds nobody else holding a claim is admitted whatever it claims: a job that does not
fit alone cannot be helped by waiting.

**It cannot stall a run.** The wait is bounded (`MAX_WAIT_SECONDS`); after it the job is admitted if fewer jobs of its kind run
than the old static table allows (`workers.stage_worker_cap`), which is the rule that stood before the gate, and otherwise it goes on
waiting for the table cap to let it in. Every decision is written to the log (an admission at debug level, a wait and a late admission as lines the
run's log shows).

**Claims are learned.** A claim starts from the table (`claims`) and rises to the largest peak seen for the kind: the worker's
heartbeat samples what the worker and its horse hold while a job runs (`refresh`), and `release` keeps the peak (`learn`). A job
killed for memory is retried once with a larger claim (`grow`), and the run stops with the job named if it dies again
(`units.horse_killed`).

Nothing here changes what a run generates: the work seed does not depend on how many jobs run at once or in what order
(`work.task_seed`).
"""

import contextlib
import dataclasses
import json
import os
import socket
import time
import uuid
from dataclasses import dataclass

import psutil

from planetgen.queue import work
from planetgen.util import log

MB = 2 ** 20

PREFIX = "planetgen:gate"
"""str: The front of the gate's Redis keys (a test gives it a prefix of its own)."""

RESERVE_MIN_MB = 1024
RESERVE_SHARE = 0.08
"""The memory the gate never hands out (MB): the larger of `RESERVE_MIN_MB` and this share of the host's. MySQL's buffer pool grows, Redis,
the web server and the kernel's own needs come out of it. Stand-in: a design edge, not a measurement."""

KEEP_CPU_PER_PLACE = 0.5
KEEP_CPU_SHARE = 0.25
"""The cores each reserved place keeps free against the bulk and maintenance lanes (OPS.130), and the most of the host's cores that all of them keep: a
generation run or a warm-up that fills the cores still leaves the site's jobs a core to start on."""

KEEP_MB_PER_PLACE = 700
"""int: What each place reserved for the interactive and long lanes keeps free against the bulk and maintenance lanes, MB (a lean
worker's figure, `workers.LEAN_WORKER_MEMORY_MB`)."""

CPU_SHARE = 1.0
INTERACTIVE_CPU_SHARE = 1.25
"""The share of the host's cores the claims may add up to: for the bulk and maintenance lanes, and for the interactive and long lanes, which may go a quarter over so
a visitor's job is not held by a generation run that fills the cores."""

LOAD_LIMIT = 1.25
INTERACTIVE_LOAD_LIMIT = 2.0
"""The one-minute load average per core over which an ordinary job waits, and an interactive or long job."""

POLL_SECONDS = 5.0
"""float: How often a waiting job reads the host again."""

MAX_WAIT_SECONDS = 300.0
"""float: How long a job waits before the table cap decides (`wait_is_over`)."""

CLAIM_SECONDS = 120.0
"""float: How long a claim lives without the worker renewing it (`refresh`, at every heartbeat)."""

GROW_FACTOR = 1.5
"""float: How much a job killed for memory claims more on its retry (`grow`)."""

LEARN_MARGIN = 1.10
"""float: A learned peak is claimed with this margin on top."""

MIN_CLAIM_MB = 150
"""int: The smallest claim."""

LOCK_SECONDS = 10
"""int: How long the decision lock lives, should its holder die inside it."""

LANE_BULK, LANE_MAINTENANCE, LANE_INTERACTIVE, LANE_LONG, LANE_OTHER = "bulk", "maintenance", "interactive", "long", "other"
LANES = (LANE_BULK, LANE_MAINTENANCE, LANE_INTERACTIVE, LANE_LONG, LANE_OTHER)
"""The lanes a claim sits in."""

SHORT_LANES = (LANE_INTERACTIVE, LANE_LONG)
"""tuple: The lanes whose places are reserved (`workers.reserved_places`)."""

DISABLE_ENV = "PLANETGEN_MEMORY_GATE"
"""str: Set to `0` to turn the gate off: every claim is admitted at once (the old behaviour: the table cap alone)."""


def _key(name):
    return f"{PREFIX}:{name}"


def enabled():
    """Whether the gate admits by the rule (not turned off by `DISABLE_ENV`)."""
    return os.environ.get(DISABLE_ENV, "1").strip().lower() not in ("0", "off", "false", "no")


@dataclass(frozen=True)
class Claim:
    """What one job asks of the host: its `kind` (the key claims are learned under), `mb` of memory, `cpu` cores, the `lane` it is on and
    whether it needs the bulk lane to itself (`exclusive`: the cut of a stage). `solo` (OPS.145) is the claim of the restarted unit that runs alone: it asks for the whole
    memory budget of the lanes that are not the site's (`solo_mb`) and waits for every other claim of those lanes to end, and none of them is admitted while it is held."""

    kind: str
    mb: float
    cpu: float = 1.0
    lane: str = LANE_OTHER
    exclusive: bool = False
    solo: bool = False


@dataclass(frozen=True)
class Host:
    """What the host has: `total_mb`, `available_mb` (the kernel's MemAvailable, or what a cgroup limit leaves, the smaller), `cores`
    and `load_per_core`, the one-minute load average over them."""

    total_mb: float
    available_mb: float
    cores: int
    load_per_core: float


@dataclass(frozen=True)
class Holder:
    """A claim that is held: who (`name`, `pid`, `host`), the `claim`, what its processes hold `held_mb` now and the most they held (`peak_mb`)."""

    name: str
    claim: Claim
    pid: int
    host: str
    held_mb: float
    peak_mb: float
    job: str = ""


@dataclass(frozen=True)
class Decision:
    """The gate's answer: `admit`, and the numbers it came from (`free_mb` is what the claim was held against; `reason` says why not)."""

    admit: bool
    reason: str
    need_mb: float
    free_mb: float
    memory: bool = False


@dataclass(frozen=True)
class Admission:
    """What `acquire` gives back: the seconds it `waited`, whether the wait `timed_out` and the table cap let it in, and the `decision`."""

    waited: float
    timed_out: bool
    decision: Decision


def _cgroup_available_mb():
    """What the cgroup the process runs in leaves of its memory limit, MB, or `None` when there is no limit (or none to read)."""
    for base in ("/sys/fs/cgroup", "/sys/fs/cgroup/memory"):
        try:
            with open(f"{base}/memory.max", encoding="ascii") as handle:
                limit = handle.read().strip()
            with open(f"{base}/memory.current", encoding="ascii") as handle:
                current = int(handle.read().strip())
            if limit.isdigit():
                return max(0.0, (int(limit) - current) / MB)
        except (OSError, ValueError):
            continue
    return None


def read_host():
    """The `Host` now: memory from the kernel's MemAvailable (less than a cgroup's limit leaves, when there is one), cores and the load
    per core. A test patches this."""
    memory = psutil.virtual_memory()
    available = memory.available / MB
    capped = _cgroup_available_mb()
    if capped is not None:
        available = min(available, capped)
    cores = work.cpu_count()
    try:
        load = psutil.getloadavg()[0] / cores
    except (OSError, AttributeError):
        load = 0.0
    return Host(memory.total / MB, available, cores, load)


def held_mb(pid):
    """What process `pid` and every process under it hold of the host's memory, MB: the proportional set size, so a page a worker shares with
    its horse is counted once; the resident size where the system gives no Pss. 0.0 for a process that is gone. A test patches this."""
    try:
        root = psutil.Process(pid)
        total = 0.0
        for process in [root, *root.children(recursive=True)]:
            try:
                info = process.memory_full_info()
                total += getattr(info, "pss", None) or info.rss
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                continue
        return total / MB
    except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
        return 0.0


def alive(pid):
    """Whether a process with this id runs here. A test patches this."""
    return psutil.pid_exists(pid)


def _decode(value):
    return value.decode("utf-8") if isinstance(value, bytes) else value


def holders(conn, now=None):
    """
    The claims held now, as `Holder`s, with what each holds in memory. A claim whose time ran out, or whose process is gone from this
    host, is dropped here (a worker killed with `kill -9` never calls `release`).
    """
    now = time.time() if now is None else now
    here = socket.gethostname()
    found = []
    for name, raw in conn.hgetall(_key("claims")).items():
        name = _decode(name)
        try:
            row = json.loads(_decode(raw))
            claim = Claim(row["kind"], float(row["mb"]), float(row["cpu"]), row["lane"], bool(row.get("exclusive")), bool(row.get("solo")))
            pid, host = int(row["pid"]), row["host"]
            expires, peak = float(row["expires"]), float(row.get("peak", 0.0))
        except (ValueError, KeyError, TypeError):
            conn.hdel(_key("claims"), name)
            continue
        if expires < now or (host == here and not alive(pid)):
            conn.hdel(_key("claims"), name)
            continue
        current = held_mb(pid) if host == here else 0.0
        found.append(Holder(name, claim, pid, host, current, max(peak, current), row.get("job", "")))
    return found


def reserve_mb(host):
    """The memory the gate keeps back of `host`'s, MB."""
    return max(RESERVE_MIN_MB, RESERVE_SHARE * host.total_mb)


def kept_mb(held, places=None):
    """
    What the places reserved for the interactive and long lanes keep free against the other lanes, MB: `places` (default
    `queue.interactive_workers` plus `queue.long_workers`) at `KEEP_MB_PER_PLACE` each, less what those lanes' claims hold already.
    """
    if places is None:
        from planetgen.util.settings import get_settings
        section = get_settings().queue
        places = section.interactive_workers + section.long_workers
    in_use = sum(holder.claim.mb for holder in held if holder.claim.lane in SHORT_LANES)
    return max(0.0, places * KEEP_MB_PER_PLACE - in_use)


def kept_cores(held, host, places=None):
    """What the places reserved for the interactive and long lanes keep free of the host's cores against the other lanes (OPS.130): `places` at
    `KEEP_CPU_PER_PLACE` each, less what those lanes' claims hold already, and never more than `KEEP_CPU_SHARE` of the host's cores."""
    if places is None:
        from planetgen.util.settings import get_settings
        section = get_settings().queue
        places = section.interactive_workers + section.long_workers
    in_use = sum(holder.claim.cpu for holder in held if holder.claim.lane in SHORT_LANES)
    return max(0.0, min(places * KEEP_CPU_PER_PLACE, KEEP_CPU_SHARE * host.cores) - in_use)


def solo_mb(host, held=(), places=None):
    """The memory a solo claim asks for, MB: the host's, less what it keeps for itself (`reserve_mb`) and for the interactive and long lanes (`kept_mb`), so the site keeps its places."""
    return max(float(MIN_CLAIM_MB), host.total_mb - reserve_mb(host) - kept_mb(held, places))


def with_solo_mb(claim, host, places=None):
    """`claim`, a solo one, with the whole budget as its memory (any other claim as it is)."""
    return dataclasses.replace(claim, mb=max(claim.mb, solo_mb(host, places=places))) if claim.solo else claim


def decide(claim, host, held, own_mb=0.0, places=None, name=None):
    """
    The gate's rule (pure: the host and the claims held are given). `held` are the other claims (`holders`); `own_mb` is what the
    worker asking already holds (credited against the claim); `name` is the asking worker, whose own old claim is never counted.
    """
    others = [holder for holder in held if holder.name != name]
    short = claim.lane in SHORT_LANES
    need = max(0.0, claim.mb - own_mb)
    if not others:
        return Decision(True, "alone: nothing else holds a claim", need, host.available_mb)
    if claim.solo:
        busy = sorted({holder.claim.kind for holder in others if holder.claim.lane not in SHORT_LANES})
        if busy:
            return Decision(False, f"a solo run waits for {', '.join(busy[:3])} to end", need, host.available_mb)
        return Decision(True, "solo: no claim of the bulk or maintenance lanes is held", need, host.available_mb)
    if not short and any(holder.claim.solo for holder in others):
        return Decision(False, "a solo run holds the host", need, host.available_mb)
    if claim.exclusive and any(holder.claim.lane == claim.lane for holder in others):
        names = sorted({holder.claim.kind for holder in others if holder.claim.lane == claim.lane})
        return Decision(False, f"waits for the {claim.lane} lane to be clear of {', '.join(names[:3])}", need, host.available_mb)
    unrealised = sum(max(0.0, holder.claim.mb - holder.held_mb) for holder in others if not (short and holder.claim.solo))      # the site's places are outside its budget
    keep = 0.0 if short else kept_mb(others, places)
    reserve = reserve_mb(host) / (2.0 if short else 1.0)
    free = host.available_mb - reserve - keep - unrealised
    if need > free:
        return Decision(False, f"needs {need:,.0f} MB, {free:,.0f} MB are free ({host.available_mb:,.0f} available, {reserve:,.0f} kept for the host, "
                               f"{keep:,.0f} for the short lanes, {unrealised:,.0f} still to be taken by the claims held)", need, free, True)
    cap = host.cores * (INTERACTIVE_CPU_SHARE if short else CPU_SHARE) - (0.0 if short else kept_cores(others, host, places))
    busy = sum(holder.claim.cpu for holder in others)
    if busy + claim.cpu > cap:
        return Decision(False, f"needs {claim.cpu:g} core(s), {busy:g} of {cap:g} are claimed", need, free)
    limit = INTERACTIVE_LOAD_LIMIT if short else LOAD_LIMIT
    if host.load_per_core > limit:
        return Decision(False, f"the load is {host.load_per_core:.2f} a core, over {limit:g}", need, free)
    return Decision(True, f"fits: needs {need:,.0f} MB of {free:,.0f} free", need, free)


@contextlib.contextmanager
def _decision_lock(conn, wait=LOCK_SECONDS, sleep=time.sleep):
    """The one lock that makes a decision and the claim it writes one step: held by the process that decides, for at most `LOCK_SECONDS`."""
    token = uuid.uuid4().hex
    deadline = time.monotonic() + wait
    while not conn.set(_key("lock"), token, nx=True, ex=LOCK_SECONDS):
        if time.monotonic() > deadline:
            break           # the lock's holder died inside it and the expiry has not come: decide without it rather than stall
        sleep(0.05)
    try:
        yield
    finally:
        if _decode(conn.get(_key("lock"))) == token:
            conn.delete(_key("lock"))


def _write(conn, name, claim, pid, job, now):
    conn.hset(_key("claims"), name, json.dumps({
        "kind": claim.kind, "mb": claim.mb, "cpu": claim.cpu, "lane": claim.lane, "exclusive": claim.exclusive, "solo": claim.solo, "pid": pid,
        "host": socket.gethostname(), "job": job, "started": now, "expires": now + CLAIM_SECONDS, "peak": 0.0}))


def wait_is_over(conn, claim, held, cap_of):
    """
    Whether the bounded wait is over for good: fewer jobs of this kind run than the old static table allows (`cap_of(kind)`, the most
    workers of the kind at once, `None` for a kind the table never limited).
    """
    cap = cap_of(claim.kind)
    return cap is None or sum(1 for holder in held if holder.claim.kind == claim.kind) < cap


def acquire(conn, claim, name, job="", pid=None, cap_of=None, on_wait=None, clock=time.monotonic, sleep=time.sleep, max_wait=MAX_WAIT_SECONDS,
            scheduler=None):
    """
    Takes `claim` for the worker `name` (the job `job`), waiting until the host has room (`decide`), or, once `max_wait` seconds have gone,
    until the old table cap lets it in (`wait_is_over`). `on_wait(decision, waited)` is called at every wait, for the heartbeat and for
    a worker to stop (it may raise). With the gate off (`DISABLE_ENV`) it admits at once.

    `scheduler` (OPS.130, `planetgen.queue.scheduler.Scheduler`) decides which waiting job goes next: it is asked in place of `decide`, under the
    same lock, and told when the job is done waiting however that ends (`done`).

    Returns:
        Admission: How long it waited, and whether the table cap admitted it.
    """
    try:
        return _acquire(conn, claim, name, job, os.getpid() if pid is None else pid, cap_of or (lambda kind: None), on_wait, clock, sleep, max_wait, scheduler)
    finally:
        if scheduler is not None:
            scheduler.done(conn, name)


def _acquire(conn, claim, name, job, pid, cap_of, on_wait, clock, sleep, max_wait, scheduler):
    started = clock()
    announced = float("-inf")
    memory_since = None
    asked = claim
    while True:
        waited = clock() - started
        with _decision_lock(conn):
            held = holders(conn)
            host = read_host()
            claim = with_solo_mb(asked, host)
            if enabled():
                if scheduler is None:
                    decision = decide(claim, host, held, own_mb=held_mb(pid), name=name)
                else:
                    decision = scheduler.decide(conn, claim, name, host, held, own_mb=held_mb(pid), waited=waited, pid=pid)
            else:
                decision = Decision(True, "the gate is off", claim.mb, host.available_mb)
            timed_out = False
            if not decision.admit and waited >= max_wait and wait_is_over(conn, claim, [h for h in held if h.name != name], cap_of):
                decision, timed_out = Decision(True, f"waited {waited:,.0f} s; the table cap lets it in ({decision.reason})", decision.need_mb, decision.free_mb), True
            if decision.admit:
                _write(conn, name, claim, pid, job, time.time())
                log.debug(f"gate admits {claim.kind} ({claim.mb:,.0f} MB, {claim.cpu:g} core(s), {claim.lane}) for {name} after {waited:.1f} s: {decision.reason}", area="queue")
                if timed_out:
                    log.warning(f"The memory gate let {claim.kind} in after {waited:,.0f} s of waiting, by the table cap: {decision.reason}")
                return Admission(waited, timed_out, decision)
        if decision.memory:
            memory_since = clock() if memory_since is None else memory_since
            _adapt_to_wait(conn, claim, name, clock() - memory_since)
        else:
            memory_since = None
        if waited - announced >= 60.0:
            announced = waited
            log.normal(f"The memory gate holds {claim.kind} ({claim.mb:,.0f} MB) for {name}: {decision.reason}.")
        if on_wait is not None:
            on_wait(decision, waited)
        sleep(POLL_SECONDS)


def _adapt_to_wait(conn, claim, name, memory_waited):
    """A bulk job has waited `memory_waited` seconds for memory: tells the adaptive control (`queue.adaptive`), which lowers the worker limit once it is long. Never raises."""
    if claim.lane != LANE_BULK or not enabled():
        return
    try:
        from planetgen.queue import adaptive
        adaptive.on_memory_wait(conn, name, claim.kind, memory_waited)
    except Exception as exc:  # noqa: BLE001 -- the control is an optimisation; a wait at the gate goes on without it
        log.debug("gate: the adaptive worker control failed: %s", exc, area="queue")


def refresh(conn, name, now=None):
    """
    Renews the claim of `name` and notes the most its processes have held (called at the worker's heartbeat while a job runs). Returns
    what they hold now, MB, or 0.0 when the worker holds no claim.
    """
    raw = conn.hget(_key("claims"), name)
    if raw is None:
        return 0.0
    now = time.time() if now is None else now
    row = json.loads(_decode(raw))
    current = held_mb(int(row["pid"]))
    row["peak"] = max(float(row.get("peak", 0.0)), current)
    row["expires"] = now + CLAIM_SECONDS
    _note_usage(row, int(row["pid"]), current, now)
    conn.hset(_key("claims"), name, json.dumps(row))
    return current


def _note_usage(row, pid, current, now):
    """OPS.155: keeps beside the peak what the claim needs for an average of the memory and the processor (cores) over its life: the sum and count of the readings, the
    processor seconds at the first and the last look and the busiest interval between two looks. Best effort."""
    try:
        from planetgen.queue import usage
        if current > 0.0:
            row["sum"], row["n"] = float(row.get("sum", 0.0)) + current, int(row.get("n", 0)) + 1
        used = usage.cpu_seconds(pid)
        if "cpu0" not in row:
            row.update(cpu0=used, cpu_last=used, t0=now, t_last=now)
        elif now - float(row["t_last"]) >= 1.0:
            row["cpu_peak"] = max(float(row.get("cpu_peak", 0.0)), max(used - float(row["cpu_last"]), 0.0) / (now - float(row["t_last"])))
            row["cpu_last"], row["t_last"] = used, now
    except Exception as exc:        # noqa: BLE001 -- an average is not worth a heartbeat
        log.swallowed("queue", "noting the usage of a gated job at its heartbeat", exc)


def release(conn, name):
    """Gives up the claim of `name` and keeps the peak its processes held as what the kind may need (`learn`). Returns the peak, MB."""
    raw = conn.hget(_key("claims"), name)
    if raw is None:
        return 0.0
    conn.hdel(_key("claims"), name)
    row = json.loads(_decode(raw))
    peak = float(row.get("peak", 0.0))
    if peak > 0.0:
        learn(conn, row["kind"], peak)
    if not str(row["kind"]).startswith("unit:"):        # OPS.155: a unit's own row is written by the unit, with its items
        from planetgen.queue import jobprofiles
        lived = float(row["t_last"]) - float(row["t0"]) if "t0" in row else 0.0
        jobprofiles.record(row["kind"], seconds=max(0.0, time.time() - float(row.get("started", time.time()))), peak_mb=peak or None,
                           avg_mb=float(row["sum"]) / int(row["n"]) if row.get("n") else None, cpu_peak=row.get("cpu_peak"),
                           cpu_avg=(float(row["cpu_last"]) - float(row["cpu0"])) / lived if lived >= 1.0 else None)
    from planetgen.queue import scheduler
    if not row.get("solo"):         # OPS.145: a solo run is not a turn of its kind, so the jobs it held back are not behind it afterwards
        scheduler.charge(conn, row["kind"], max(0.0, time.time() - float(row.get("started", time.time()))) * float(row.get("cpu", 1.0)))
    return peak


def learned(conn):
    """The largest peak seen for each kind, MB (`{kind: mb}`)."""
    return {_decode(kind): float(_decode(value)) for kind, value in conn.hgetall(_key("learned")).items()}


def learn(conn, kind, peak_mb):
    """Raises what is known of the kind's size to `peak_mb` if that is more than has been seen. Returns what is known now."""
    seen = conn.hget(_key("learned"), kind)
    top = max(float(_decode(seen)) if seen is not None else 0.0, float(peak_mb))
    if seen is None or top > float(_decode(seen)):
        conn.hset(_key("learned"), kind, f"{top:.3f}")
    return top


def forget(conn):
    """Drops everything learned (the peaks and the retry claims): the next run claims from the table again."""
    conn.delete(_key("learned"))
    for key in conn.scan_iter(match=_key("override:*")):
        conn.delete(key)


def claim_mb(conn, kind, table_mb, job_id=None, known=None, token=None):
    """
    What a job of `kind` claims, MB (from the job profile when the machine has history for it, `profileclaim`, else): the table's figure (`table_mb`), raised to the peak seen for the kind with its margin
    (`LEARN_MARGIN`), and to the larger claim a retry was given (`grow`). A kind with no table figure (`table_mb` `None`) claims what has
    been seen, or `MIN_CLAIM_MB` if nothing has: the caller gives a default for it.
    """
    known = learned(conn) if known is None else known
    peak = known.get(kind)
    figure = table_mb
    if peak:
        figure = peak * LEARN_MARGIN if table_mb is None else max(table_mb, peak * LEARN_MARGIN)
    from planetgen.queue import profileclaim
    profiled = profileclaim.claim(conn, kind, table_mb, token=token)      # OPS.156: the profile's figure, when this machine has history; None gives the figure above
    if profiled is not None:
        figure = profiled
    if job_id:
        retried = conn.get(_key(f"override:{job_id}"))
        if retried is not None:
            figure = max(figure or 0.0, float(_decode(retried)))
    return max(float(figure if figure is not None else MIN_CLAIM_MB), float(MIN_CLAIM_MB))


def grow(conn, kind, job_id, claimed_mb, peak_mb=0.0):
    """
    A job of `kind` was killed for memory while it claimed `claimed_mb`: its retry (`job_id`) claims `GROW_FACTOR` times that, or what it
    was seen to hold, whichever is more, and the kind's size is learned from it. Returns the new claim, MB.
    """
    bigger = max(claimed_mb * GROW_FACTOR, peak_mb * LEARN_MARGIN)
    conn.set(_key(f"override:{job_id}"), f"{bigger:.3f}", ex=86400)
    learn(conn, kind, bigger / LEARN_MARGIN)
    return bigger


def kill_note(conn, name):
    """What the worker `name` held when it was killed, for a message: its claim, the peak seen and what is held now, or `None`."""
    raw = conn.hget(_key("claims"), name)
    if raw is None:
        return None
    row = json.loads(_decode(raw))
    return {"kind": row["kind"], "claimed_mb": float(row["mb"]), "peak_mb": float(row.get("peak", 0.0)), "job": row.get("job", "")}


@contextlib.contextmanager
def hold(conn, claim, name=None, **options):
    """
    A claim around a stretch of work that is not an RQ job (a test gate, an import, an inline command): waits for room, holds the claim and
    gives it up, learning the peak. `name` is unique to the holder (default `kind:pid:random`).
    """
    name = name or f"{claim.kind}:{os.getpid()}:{uuid.uuid4().hex[:6]}"
    if "scheduler" not in options:
        from planetgen.queue import scheduler
        options["scheduler"] = scheduler.Scheduler()
    acquire(conn, claim, name, **options)
    try:
        yield name
    finally:
        release(conn, name)


def snapshot(conn):
    """What the gate holds, for a status page or a report: the host, the claims and what has been learned."""
    held = holders(conn)
    host = read_host()
    return {"host": host, "reserve_mb": reserve_mb(host), "claims": held, "learned": learned(conn),
            "claimed_mb": sum(holder.claim.mb for holder in held), "held_mb": sum(holder.held_mb for holder in held)}
