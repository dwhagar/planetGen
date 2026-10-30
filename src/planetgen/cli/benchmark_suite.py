# planetgen/cli/benchmark_suite.py

"""
`planetgen benchmark --suite`: a performance profile of every generation option, one repeatable command.

`planetgen benchmark` (PERF.31) times a small galaxy from the plan to the fill. The suite times each part of a galaxy
run on its own, at a size you pick, in scratch databases (nothing is written to the real one), and writes the result as
JSON for the next run to compare with:

- **the plan** (the skeleton: the shape, the clusters and satellites, the outline), and the cost of its parts;
- **every scatter kind**, each by itself on the stored plan: the massive stars, the bright stars (with the halo and
  the satellites' and clusters' stars), the phenomena, the nebula field and the dimmer band;
- **the fill**, one `planetgen galaxy --shell` run, with the time of each stage it records (the sectors, the
  neighbour link, the backfill, the sector paths) and the systems it made;
- **the population pass**, **the habitability re-score**, **the map-luminosity fill** (after emptying the column it
  fills), **the nebula block totals** and **the cache warm-up**.

Every step is a child process of the real command, so a worker's cost is in it. For each step the suite records the
wall time, the CPU the whole machine spent (from `/proc/stat`, with the database server's and Redis's own share
apart), the statements and rows the database server counted, the peak memory of the generation processes (the sum and
the largest, sampled twice a second), what the step made (stars, systems, rows) and its time per thing made, and the
stage rows the command stored. With `--repeats N` the whole pipeline runs N times in fresh databases and each step
reports its fastest, median and slowest run. The machine must be the suite's alone while it runs: CPU and database
counters are the machine's, and the Redis queue is the one the real commands use. Linux only (`/proc`).

The sizes are the galaxy's stellar mass (the `plan --mass` lever), which sets how many components, shells and stars
the steps meet: `tiny` is for tests, `small`, `medium` (a default-size galaxy) and `large` (a Milky-Way-size one).
"""

import argparse
import json
import os
import platform
import re
import signal
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import uuid

from planetgen import _version
from planetgen.cli import benchmark
from planetgen.db import store

SUITE_VERSION = 1
"""int: The layout of the suite's JSON; a change that breaks a reader of an earlier result raises it."""

DEFAULT_SEED = benchmark.DEFAULT_SEED

SIZES = {
    "tiny": {"mass_solar": 2.0e6, "sectors": 2, "floor": 1.0e4,
             "note": "for tests: a few thousand stars, seconds a step"},
    "small": {"mass_solar": 3.0e7, "sectors": 6, "floor": None,
              "note": "a dwarf galaxy: about 370 shells, 160 globular clusters"},
    "medium": {"mass_solar": 4.0e9, "sectors": 12, "floor": None,
               "note": "a default-size galaxy (the seed's own draw, 3e9 to 6e9 solar masses): about 1,085 shells, 22,000 clusters"},
    "large": {"mass_solar": 5.0e10, "sectors": 24, "floor": None,
              "note": "a Milky-Way-size galaxy"},
}
"""dict: The sizes: the galaxy's stellar mass (solar masses), the sectors the fill makes, and the bright-star floor
(solar luminosities) of the scatters (`None`: the generator's own default)."""

FILL_SCALE_LENGTHS = 2.0
"""float: The fill starts in the shell this many disc scale lengths out: in the disc, not in the crowded core."""

DEFAULT_STEP_TIMEOUT_MINUTES = 240
"""int: How long one step may run before it is stopped and recorded as not finished."""

SAMPLE_SECONDS = 0.5
"""float: How often the sampler reads the processes' memory."""

STEP_ORDER = ("plan", "scatter-mass", "scatter-luminosity", "scatter-phenomena", "scatter-nebula-field", "scatter-band",
              "fill-nebula-blocks", "fill", "population", "territories", "rescore-habitability", "fill-map-luminosity", "warm")
"""tuple: The steps, in the order they run (each needs the ones before it)."""

STEP_LABELS = {
    "plan": "Plan the galaxy (skeleton, clusters, satellites, outline)",
    "scatter-mass": "Scatter the massive stars",
    "scatter-luminosity": "Scatter the bright stars",
    "scatter-phenomena": "Scatter the phenomena",
    "scatter-nebula-field": "Draw the nebula field",
    "scatter-band": "Add the dimmer band of bright stars",
    "fill-nebula-blocks": "Make the nebula block totals",
    "fill": "Fill sectors (planetgen galaxy --shell)",
    "population": "Run the population pass",
    "territories": "Recompute the territories",
    "rescore-habitability": "Re-score the habitability",
    "fill-map-luminosity": "Fill the map luminosity",
    "warm": "Warm the caches",
}
"""dict: What each step is called in the report."""

# What each step makes, as `(noun, SQL)`: the SQL counts the rows in the scratch database, and the step made the
# difference between the count after and before.
MADE = {
    "plan": ("components", "SELECT COUNT(*) FROM galaxy_components"),
    "scatter-mass": ("stars", "SELECT COUNT(*) FROM bright_stars"),
    "scatter-luminosity": ("stars", "SELECT COUNT(*) FROM bright_stars"),
    "scatter-phenomena": ("phenomena", "SELECT COUNT(*) FROM phenomenon_scatter"),
    "scatter-nebula-field": ("clouds", "SELECT COUNT(*) FROM nebula_field"),
    "scatter-band": ("stars", "SELECT COUNT(*) FROM bright_stars"),
    "fill-nebula-blocks": ("blocks", "SELECT COUNT(*) FROM nebula_block_stats"),
    "fill": ("systems", "SELECT COUNT(*) FROM star_systems"),
    "population": ("polities", "SELECT COUNT(*) FROM polities"),
    "territories": ("owned systems", "SELECT COUNT(*) FROM system_owners"),
    "rescore-habitability": ("planets and moons", "SELECT (SELECT COUNT(*) FROM planets) + (SELECT COUNT(*) FROM moons)"),
    "fill-map-luminosity": ("phenomena", "SELECT COUNT(*) FROM phenomenon_scatter"),
    "warm": ("sectors", "SELECT COUNT(*) FROM sectors"),
}
"""dict: Per step, what it makes and how to count it. A step that changes rows without adding any (the re-score, the
map-luminosity fill) is measured against the rows it works on, so `count` is the table's size after the step."""

NO_DATABASE_ARGS = frozenset({"warm"})
"""frozenset[str]: The steps whose command takes no database options: it reads the database the settings name, which
the suite sets in the environment."""

ABSOLUTE_COUNT = frozenset({"rescore-habitability", "fill-map-luminosity", "warm", "fill-nebula-blocks", "territories"})
"""frozenset[str]: The steps whose count is the table's size after them, not the rows they added."""


def _positive_int(text):
    """`argparse` type for a count of at least one."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid int value: {text!r}") from None
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, got {text!r}")
    return value


def _suffix(text):
    """`argparse` type for the suffix of an earlier run's scratch databases: the eight hex digits the suite named them by."""
    if not re.fullmatch(r"[0-9a-f]{8}", text):
        raise argparse.ArgumentTypeError(f"must be the eight hex digits the suite printed, got {text!r}")
    return text


def _positive_mass(text):
    """`argparse` type for a mass in solar masses: a finite number above zero."""
    from planetgen.cli.generate import finite_float
    value = finite_float(text)
    if value <= 0.0:
        raise argparse.ArgumentTypeError(f"must be above zero, got {text!r}")
    return value


def add_suite_arguments(parser):
    """Adds the suite's options to the `benchmark` parser."""
    group = parser.add_argument_group("the performance suite (--suite)")
    group.add_argument("--suite", action="store_true",
                       help="Run the performance suite of every generation option instead of the small fill: "
                            "the plan, each scatter kind, the fill, the population pass, the re-score, the "
                            "map-luminosity fill, the nebula block totals and the warm-up, one step at a time.")
    group.add_argument("--size", choices=sorted(SIZES), default="small",
                       help="With --suite: the galaxy's size (its stellar mass): " + "; ".join(
                           f"{name} {size['mass_solar']:g} solar masses" for name, size in SIZES.items()) + ".")
    group.add_argument("--steps", metavar="STEP[,STEP...]", default=None,
                       help="With --suite: only these steps (default all, in order; a step needs the ones before it "
                            "to have run, so a list starts from the plan): " + ", ".join(STEP_ORDER) + ".")
    group.add_argument("--repeats", type=_positive_int, default=1, metavar="N",
                       help="With --suite: run the whole pipeline N times in fresh databases and report each step's "
                            "fastest, median and slowest run (default 1).")
    group.add_argument("--step-timeout", type=_positive_int, default=DEFAULT_STEP_TIMEOUT_MINUTES, metavar="MINUTES",
                       help=f"With --suite: stop a step that runs longer than this and record it as not finished "
                            f"(default {DEFAULT_STEP_TIMEOUT_MINUTES}).")
    group.add_argument("--mass", type=_positive_mass, default=None, metavar="M_SUN",
                       help="With --suite: the galaxy's stellar mass in solar masses, instead of the size's.")
    group.add_argument("--resume", type=_suffix, metavar="SUFFIX", default=None,
                       help="With --suite: carry on in the scratch databases of an earlier run that was stopped "
                            "(the suite prints their names, planetgen_suite_<SUFFIX>_r0 and "
                            "planetgen_suite_control_<SUFFIX>); give --steps the steps left to run.")
    group.add_argument("--fill-sectors", type=_positive_int, default=None, metavar="N",
                       help="With --suite: the sectors the fill makes, instead of the size's.")


def parse_steps(text):
    """The steps `--steps` names, in run order.

    Raises:
        SystemExit: A step that does not exist, or none.
    """
    if text is None:
        return list(STEP_ORDER)
    wanted = [part.strip() for part in str(text).split(",") if part.strip()]
    unknown = [part for part in wanted if part not in STEP_ORDER]
    if not wanted or unknown:
        raise SystemExit(f"planetgen benchmark --suite: --steps wants some of {', '.join(STEP_ORDER)}"
                         + (f"; not {', '.join(unknown)}." if unknown else "."))
    return [step for step in STEP_ORDER if step in wanted]


# --- The machine -----------------------------------------------------------------------------------------------------

def _clock_ticks():
    return os.sysconf("SC_CLK_TCK")


def machine_cpu_seconds():
    """The seconds the machine's CPUs have been busy (user, nice, system, interrupts; not idle, waiting for input or
    output, or stolen) since it started, or `None` where `/proc/stat` is not there."""
    try:
        with open("/proc/stat", encoding="ascii") as handle:
            fields = handle.readline().split()
    except OSError:
        return None
    user, nice, system, idle, iowait, irq, softirq = (int(value) for value in fields[1:8])
    return (user + nice + system + irq + softirq) / _clock_ticks()


def _processes():
    """`[(pid, command line)]` of every process that can be read."""
    found = []
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        try:
            with open(f"/proc/{name}/cmdline", "rb") as handle:
                found.append((int(name), handle.read().replace(b"\0", b" ").decode("utf-8", "replace").strip()))
        except OSError:
            continue
    return found


def _process_cpu_seconds(pid):
    """The seconds of CPU `pid` has used (user and system), or 0 when it is gone."""
    try:
        with open(f"/proc/{pid}/stat", encoding="ascii") as handle:
            fields = handle.read().rsplit(")", 1)[1].split()
    except OSError:
        return 0.0
    return (int(fields[11]) + int(fields[12])) / _clock_ticks()


def _rss_mb(pid):
    try:
        with open(f"/proc/{pid}/status", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    return 0.0


def _server_pids():
    """The pids of the database server and Redis when they run on this machine."""
    found = []
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        try:
            with open(f"/proc/{name}/comm", encoding="utf-8") as handle:
                if handle.read().strip() in ("mysqld", "mariadbd", "redis-server"):
                    found.append(int(name))
        except OSError:
            continue
    return found


def _is_generation_process(command):
    """Whether a process is the generation's: its command line names planetgen (the suite's own, which names
    benchmark, is left out by pid)."""
    return "planetgen" in command


class Sampler(threading.Thread):
    """Reads the generation processes' memory while a step runs: the most the processes held together and the most
    one of them held. A generation process is any whose command line names planetgen, other than the suite's."""

    def __init__(self, interval=SAMPLE_SECONDS):
        super().__init__(daemon=True)
        self.interval = interval
        self.stopped = threading.Event()
        self.peak_total_mb = 0.0
        self.peak_process_mb = 0.0
        self.peak_processes = 0
        self._own = os.getpid()

    def sample(self):
        total = largest = 0.0
        count = 0
        for pid, command in _processes():
            if pid == self._own or not _is_generation_process(command):
                continue
            rss = _rss_mb(pid)
            total += rss
            largest = max(largest, rss)
            count += 1
        self.peak_total_mb = max(self.peak_total_mb, total)
        self.peak_process_mb = max(self.peak_process_mb, largest)
        self.peak_processes = max(self.peak_processes, count)

    def run(self):
        while not self.stopped.wait(self.interval):
            self.sample()

    def finish(self):
        self.stopped.set()
        self.join()
        self.sample()


def _environment_has(pid, name, value):
    """Whether process `pid` was started with `name=value` in its environment; false for one that is gone or not ours to read."""
    try:
        with open(f"/proc/{pid}/environ", "rb") as handle:
            return f"{name}={value}".encode() in handle.read().split(b"\0")
    except OSError:
        return False


def stop_workers(control):
    """Ends the queue workers of this run's control database `control`, so the next step starts with cold ones (a worker
    keeps what it warmed, and a step that follows another within the workers' idle time would otherwise start warm). A
    worker of any other control database, another test's or the machine's own, is left running (TEST.302): the run's
    workers are the ones started with its `PLANETGEN_CONTROL_DATABASE`, which the steps it runs hand on."""
    own = os.getpid()
    for pid, command in _processes():
        if pid != own and "planetgen.cli.worker" in command and _environment_has(pid, store.CONTROL_DB_ENV_VAR, control):
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass


def _kill_group(process):
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(process.pid, sig)
        except OSError:
            return
        try:
            process.wait(timeout=30)
            return
        except subprocess.TimeoutExpired:
            continue


# --- A step ----------------------------------------------------------------------------------------------------------

def measure(argv, env, timeout_seconds, admin=None, counts_sql=None):
    """
    Runs `argv` as a child process and measures it.

    Returns:
        dict: `wall_seconds`, `cpu_seconds` (the machine's, less the servers' own), `server_cpu_seconds` (the database
            server's and Redis's, or `None`), `peak_rss_total_mb`, `peak_rss_process_mb`, `peak_processes`,
            `database` (the server's counters while it ran), `exit_code` (`None` when stopped), `timed_out` and
            `log_tail`.
    """
    servers = _server_pids()
    before_machine = machine_cpu_seconds()
    before_servers = sum(_process_cpu_seconds(pid) for pid in servers)
    before_status = benchmark._status(admin) if admin is not None else None
    sampler = Sampler()
    with tempfile.TemporaryFile() as log:
        started = time.perf_counter()
        process = subprocess.Popen(argv, env=env, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        sampler.start()
        timed_out = False
        try:
            process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill_group(process)
        except BaseException:
            _kill_group(process)
            raise
        finally:
            wall = time.perf_counter() - started
            sampler.finish()
        log.seek(0)
        tail = log.read()[-4000:].decode("utf-8", "replace")
    after_machine = machine_cpu_seconds()
    after_servers = sum(_process_cpu_seconds(pid) for pid in servers)
    after_status = benchmark._status(admin) if admin is not None else None
    server_cpu = (after_servers - before_servers) if servers else None
    cpu = None if before_machine is None else after_machine - before_machine - (server_cpu or 0.0)
    return {
        "wall_seconds": wall, "cpu_seconds": cpu, "server_cpu_seconds": server_cpu,
        "peak_rss_total_mb": sampler.peak_total_mb, "peak_rss_process_mb": sampler.peak_process_mb,
        "peak_processes": sampler.peak_processes,
        "database": None if before_status is None else {name: after_status[name] - before_status[name]
                                                        for name in benchmark.STATUS_VARIABLES},
        "exit_code": None if timed_out else process.returncode, "timed_out": timed_out, "log_tail": tail,
    }


def _count(admin, database, sql):
    with admin.cursor() as cur:
        cur.execute(f"USE `{database}`")
        cur.execute(sql)
        row = cur.fetchone()
    return int(row[0] or 0)


def _stage_rows(config, control, database, after_id):
    """`(rows, last id)`: the stage rows the commands stored for `database` with an id above `after_id`."""
    import pymysql
    from planetgen.db import login
    conn = login.connect(host=config.host, port=config.port, user=config.user, password=config.password,
                         database=control, autocommit=True, cursorclass=pymysql.cursors.DictCursor)
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT id, command, stage_key, label, skipped, seconds, metrics FROM generation_stage_runs "
                        "WHERE database_name = %s AND id > %s ORDER BY id", (database, after_id))
            rows = cur.fetchall()
    finally:
        conn.close()
    stages = [{"command": row["command"], "key": row["stage_key"], "label": row["label"],
               "skipped": bool(row["skipped"]), "seconds": float(row["seconds"]),
               "metrics": json.loads(row["metrics"] or "{}")} for row in rows]
    return stages, (rows[-1]["id"] if rows else after_id)


class Pipeline:
    """One run of the steps in a scratch database of its own: builds each step's command, runs it, and records what
    it made."""

    def __init__(self, config, control, database, size, steps, seed, workers, timeout_seconds, mass=None, sectors=None):
        self.config, self.control, self.database = config, control, database
        self.size = size
        self.steps = steps
        self.seed = seed
        self.workers = workers
        self.timeout = timeout_seconds
        spec = SIZES[size]
        self.mass = spec["mass_solar"] if mass is None else mass
        self.sectors = spec["sectors"] if sectors is None else sectors
        self.floor = spec["floor"]
        self.env = benchmark._child_env(config, control, workers)
        self.env["PLANETGEN_MYSQL_DATABASE"] = database
        self.env["PLANETGEN_WARM_ON_CHANGE"] = "0"      # the automatic warm-up would run beside the next step
        self.stage_id = 0
        self.admin = benchmark._admin(config)
        scratch = store.MySQLConfig(host=config.host, port=config.port, user=config.user, password=config.password,
                                    database=database)
        store.get_connection(scratch).close()           # the schema, so the counts below have tables to read
        store.close_pool(scratch)

    def _connection(self):
        return ["--mysql-host", self.config.host, "--mysql-port", str(self.config.port),
                "--mysql-user", self.config.user, "--mysql-database", self.database]

    def _floor_args(self):
        return [] if self.floor is None else ["--bright-star-min-luminosity", f"{self.floor:g}"]

    def command(self, step):
        """The `planetgen` command line of `step` (what follows `python -m planetgen.cli.generate`), without the database."""
        if step == "plan":
            return ["plan", "--seed", self.seed, "--mass", f"{self.mass:g}"]
        if step == "scatter-mass":
            return ["scatter", "--mass", *self._floor_args()]
        if step == "scatter-luminosity":
            return ["scatter", "--luminosity", *self._floor_args()]
        if step == "scatter-phenomena":
            return ["scatter", "--phenomena"]
        if step == "scatter-nebula-field":
            return ["scatter", "--nebula-field"]
        if step == "scatter-band":
            return ["scatter", "--band", "--down-to", f"{self._band_floor():g}"]
        if step == "fill-nebula-blocks":
            return ["fill-nebula-blocks"]
        if step == "fill":
            return ["galaxy", "--shell", str(self._fill_shell()), "--limit", str(self.sectors), "--yes"]
        if step == "population":
            return ["population"]
        if step == "territories":
            return ["territories"]
        if step == "rescore-habitability":
            return ["rescore-habitability"]
        if step == "fill-map-luminosity":
            return ["fill-map-luminosity"]
        if step == "warm":
            return ["warm"]
        raise ValueError(step)

    def _band_floor(self):
        """The floor one level dimmer than the scatter's: a tenth of it."""
        from planetgen import tuning
        base = self.floor if self.floor is not None else float(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
        return base / 10.0

    def _fill_shell(self):
        """The shell the fill starts in: `FILL_SCALE_LENGTHS` disc scale lengths out, inside the outline."""
        with self.admin.cursor() as cur:
            cur.execute(f"USE `{self.database}`")
            cur.execute("SELECT disk_scale_length_pc, edge_pc, outer_shell_index FROM galaxy_shape WHERE id = 1")
            scale, edge, outer = cur.fetchone()
        return max(1, min(int(round(FILL_SCALE_LENGTHS * scale / edge)), int(outer) - 1))

    def prepare(self, step):
        """What a step needs done before it runs, outside its time."""
        if step == "scatter-nebula-field":                 # an earlier step drew the field with the phenomena
            subprocess.run([sys.executable, "-m", "planetgen.cli.generate", "clear-scatter", "--nebula-field",
                            *self._connection()], env=self.env, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, check=False)
        if step == "fill-map-luminosity":
            with self.admin.cursor() as cur:
                cur.execute(f"USE `{self.database}`")
                cur.execute("UPDATE phenomenon_scatter SET map_luminosity_sol = NULL")
            self.admin.commit()

    def run_step(self, step):
        """Runs `step` and returns its record."""
        self.prepare(step)
        _skipped, self.stage_id = _stage_rows(self.config, self.control, self.database, self.stage_id)
        stop_workers(self.control)
        noun, sql = MADE[step]
        before = _count(self.admin, self.database, sql)
        words = self.command(step)
        argv = [sys.executable, "-m", "planetgen.cli.generate", *words, *([] if step in NO_DATABASE_ARGS else self._connection())]
        result = measure(argv, self.env, self.timeout, self.admin)
        stop_workers(self.control)
        after = _count(self.admin, self.database, sql)
        stages, self.stage_id = _stage_rows(self.config, self.control, self.database, self.stage_id)
        made = after if step in ABSOLUTE_COUNT else after - before
        record = {
            "step": step, "label": STEP_LABELS[step], "command": words,
            **result, "made": {"noun": noun, "count": made},
            "seconds_per": None if not made or result["timed_out"] else result["wall_seconds"] / made,
            "stages": [stage for stage in stages if not stage["skipped"]],
        }
        if record["exit_code"] not in (0, None):
            record["failed"] = True
        return record

    def run(self):
        records = []
        for step in self.steps:
            record = self.run_step(step)
            records.append(record)
            print(progress_line(record), flush=True)
            if record.get("failed") or record["timed_out"]:
                break
        self.admin.close()
        return records


def finished(runs, repeats):
    """Whether every repeat ran and none of its steps failed or stopped: only then are the scratch databases dropped."""
    return len(runs) == repeats and not any(
        record.get("failed") or record["timed_out"] or record["exit_code"] for run in runs for record in run)


def progress_line(record):
    """The line printed when a step ends, so a run that is stopped before its report keeps what finished."""
    state = "NOT FINISHED" if record["timed_out"] else "FAILED" if record.get("failed") or record["exit_code"] else "done"
    made = record["made"]
    return (f"  {state}: {record['label']} in {_seconds(record['wall_seconds'])}"
            f" ({made['count']:,} {made['noun']}, client CPU {_seconds(record['cpu_seconds'])})")


def summarize(runs):
    """Each step's fastest, median and slowest wall time over the repeats, with the first run's record."""
    summary = []
    for index, first in enumerate(runs[0]):
        times = [run[index]["wall_seconds"] for run in runs if index < len(run) and not run[index]["timed_out"]
                 and not run[index].get("failed")]
        summary.append({
            "step": first["step"], "label": first["label"], "runs": len(times),
            "wall_seconds": {"min": min(times), "median": statistics.median(times), "max": max(times)} if times else None,
            "spread": None if len(times) < 2 else (max(times) - min(times)) / statistics.median(times),
        })
    return summary


def rank(records):
    """`[(step, wall_seconds, share of the total)]`, the costliest first, from one run's records."""
    total = sum(record["wall_seconds"] for record in records) or 1.0
    ranked = sorted(records, key=lambda record: -record["wall_seconds"])
    return [{"step": record["step"], "wall_seconds": record["wall_seconds"], "share": record["wall_seconds"] / total}
            for record in ranked]


def machine_facts(admin, workers):
    with admin.cursor() as cur:
        cur.execute("SELECT VERSION()")
        version = cur.fetchone()[0]
    return {"cpus": os.cpu_count(), "workers": workers, "python": platform.python_version(),
            "platform": platform.platform(), "mysql": version, "release": _version.__version__,
            "buffer_pool_bytes": benchmark._buffer_pool_bytes(admin)}


def _drop_scratch(admin, args, runs, scratch, control, suffix):
    """Drops the scratch databases once every repeat ran and none of its steps failed or stopped, unless `--keep` asks
    to leave them (PERF.175); otherwise says how to carry on (`--resume`) or what to drop."""
    if not finished(runs, args.repeats):
        print(f"A step did not finish: the scratch databases are kept; carry on with --resume {suffix} "
              f"(drop planetgen_suite_{suffix}_r0 and {control} when done).", flush=True)
    elif args.keep:
        print(f"--keep: the scratch databases are kept (drop planetgen_suite_{suffix}_r0 and {control} when done).",
              flush=True)
    else:
        with admin.cursor() as cur:
            for database in scratch:
                cur.execute(f"DROP DATABASE IF EXISTS `{database}`")
            cur.execute(f"DROP DATABASE IF EXISTS `{control}`")


def run_suite(args):
    """`planetgen benchmark --suite`: runs, prints and (with `--output`) saves the report."""
    if args.repeats < 1:
        raise SystemExit("planetgen benchmark --suite: --repeats wants a positive whole number.")
    if not os.path.isdir("/proc"):
        raise SystemExit("planetgen benchmark --suite: needs /proc, so it runs on Linux.")
    steps = parse_steps(args.steps)
    workers = benchmark._workers(args.workers)[0]
    config = store.mysql_config_from_args(args)
    admin = benchmark._admin(config)
    resumed = args.resume is not None
    if resumed and args.repeats != 1:
        raise SystemExit("planetgen benchmark --suite: --resume carries on one run, so --repeats stays 1.")
    suffix = args.resume if resumed else uuid.uuid4().hex[:8]
    control = f"planetgen_suite_control_{suffix}"
    scratch = []
    runs = []
    try:
        if resumed:
            with admin.cursor() as cur:
                cur.execute("SELECT SCHEMA_NAME FROM information_schema.SCHEMATA WHERE SCHEMA_NAME = %s",
                            (f"planetgen_suite_{suffix}_r0",))
                if not cur.fetchall():
                    raise SystemExit(f"planetgen benchmark --suite: no scratch database planetgen_suite_{suffix}_r0 "
                                     f"to resume.")
        else:
            with admin.cursor() as cur:
                cur.execute(f"CREATE DATABASE `{control}`")
            benchmark._make_control_schema(config, control)
        print(f"Scratch databases: planetgen_suite_{suffix}_r0 and {control} (to carry on after a stop: --resume {suffix})",
              flush=True)
        facts = machine_facts(admin, workers)
        startup = None
        for repeat in range(args.repeats):
            database = f"planetgen_suite_{suffix}_r{repeat}"
            if not resumed:
                with admin.cursor() as cur:
                    cur.execute(f"CREATE DATABASE `{database}`")
            scratch.append(database)
            pipeline = Pipeline(config, control, database, args.size, steps, args.seed, workers,
                                args.step_timeout * 60, mass=args.mass, sectors=args.fill_sectors)
            runs.append(pipeline.run())
            if startup is None:
                startup = benchmark._startup_seconds(config, database, pipeline.env)
    finally:
        _drop_scratch(admin, args, runs, scratch, control, suffix)
        admin.close()
    result = {
        "suite_version": SUITE_VERSION, "size": args.size, "mass_solar": args.mass or SIZES[args.size]["mass_solar"],
        "seed": args.seed, "repeats": args.repeats, "machine": facts, "worker_startup": startup,
        "runs": runs, "summary": summarize(runs), "ranking": rank(runs[0]) if runs[0] else [],
    }
    print(format_report(result))
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
    return result


# --- The report ------------------------------------------------------------------------------------------------------

def _seconds(value):
    if value is None:
        return "-"
    if value >= 3600:
        return f"{value / 3600:.2f} h"
    if value >= 120:
        return f"{value / 60:.1f} min"
    return f"{value:.1f} s"


def _rate(record):
    if record["seconds_per"] is None:
        return ""
    return f"{record['seconds_per'] * 1000:,.2f} ms per {record['made']['noun'].rstrip('s')}"


def format_report(result):
    """The suite's result as text."""
    facts = result["machine"]
    lines = [f"Performance suite: size {result['size']} (galaxy of {result['mass_solar']:g} solar masses), seed "
             f"{result['seed']}, release {facts['release']}, {facts['cpus']} CPUs, {facts['workers']} worker(s), "
             f"MySQL {facts['mysql']}, {result['repeats']} run(s).", ""]
    lines.append("  step                                           wall      share   cpu      peak memory      made")
    first = result["runs"][0]
    total = sum(record["wall_seconds"] for record in first) or 1.0
    for record in first:
        status = " NOT FINISHED" if record["timed_out"] else (" FAILED" if record.get("failed") else "")
        cpu = "-" if record["cpu_seconds"] is None else _seconds(record["cpu_seconds"])
        lines.append(f"  {record['label'][:44]:44} {_seconds(record['wall_seconds']):>9} {100 * record['wall_seconds'] / total:6.1f}% "
                     f"{cpu:>8} {record['peak_rss_total_mb']:8.0f} MB ({record['peak_rss_process_mb']:.0f} one)  "
                     f"{record['made']['count']:,} {record['made']['noun']}  {_rate(record)}{status}")
        for stage in record["stages"]:
            if stage["key"] in ("job", "scatter"):
                continue
            notes = ", ".join(f"{key} {value}" for key, value in sorted(stage["metrics"].items())
                              if isinstance(value, (int, float)))
            lines.append(f"      stage {stage['label'][:38]:38} {_seconds(stage['seconds']):>9}  {notes}")
        counters = record["database"] or {}
        shown = ", ".join(f"{name} {value:,}" for name, value in counters.items() if value)
        if shown:
            lines.append(f"      database: {shown}")
        if record.get("failed"):
            lines.append("      log: " + record["log_tail"][-600:].replace("\n", "\n           "))
    if result["repeats"] > 1:
        lines += ["", "  run-to-run (fastest / median / slowest):"]
        for entry in result["summary"]:
            wall = entry["wall_seconds"]
            if wall is not None:
                lines.append(f"  {entry['label'][:44]:44} {_seconds(wall['min']):>9} / {_seconds(wall['median']):>9} / "
                             f"{_seconds(wall['max']):>9}  spread {100 * (entry['spread'] or 0):.0f}%")
    startup = result.get("worker_startup")
    if startup:
        lines += ["", f"  worker start-up: {startup['process_seconds']:.2f} s a process (imports {startup['import_seconds']:.2f} s, "
                      f"first connection {startup['connect_seconds']:.2f} s)"]
    lines += ["", "  ranked by wall time:"]
    for entry in result["ranking"]:
        lines.append(f"    {entry['share'] * 100:5.1f}%  {_seconds(entry['wall_seconds']):>9}  {STEP_LABELS[entry['step']]}")
    return "\n".join(lines)
