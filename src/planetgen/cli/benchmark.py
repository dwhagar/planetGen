# planetgen/cli/benchmark.py

"""
`planetgen benchmark` (PERF.31): where a generation run spends its time, from
the plan to a finished fill.

It makes a scratch database and control database, runs the real commands in
child processes with a fixed galaxy seed (`plan`, its scatters, then `galaxy` over a small
shell with the scatter on), and prints for each worker count asked:

- every stage's wall time and share of the run (the numbered stages every
  generation command records in `generation_stage_runs`, UX.89), with what
  the stage reports it did (shells, sectors, objects);
- the database's share: the statements and rows the server counted while the
  run went (`Com_*` and `Innodb_rows_*`), per second of the run;
- the server's `innodb_buffer_pool_size`, the size and row count of each of the run's largest tables and the
  worker start-up cost (a new Python process importing what a sector worker imports and opening its first
  database connection), so a result can be compared with the next (PERF.47);
- with `--profile`, the functions the fill spent most time in, by
  cumulative time, from a `cProfile` of one single-worker run.

The scratch databases are dropped at the end (`--keep` leaves them for a look).
The result also goes to `--output` as JSON, so a later run can be compared.
"""

import argparse
import cProfile  # noqa: F401 -- named so the profile run's module is importable here
import json
import os
import pstats
import subprocess
import sys
import tempfile
import time
import uuid

import pymysql

from planetgen.db import login, store
from planetgen.util import log

DEFAULT_SEED = "0000000000000000000000000000BEEF"
DEFAULT_SECTORS = 2
DEFAULT_SHELL = 12

STATUS_VARIABLES = ("Questions", "Com_select", "Com_insert", "Com_update", "Com_delete", "Com_commit",
                    "Innodb_rows_read", "Innodb_rows_inserted", "Innodb_rows_updated", "Innodb_rows_deleted",
                    "Innodb_data_writes", "Innodb_row_lock_waits")
"""tuple: The server counters read before and after a run."""

PROFILE_TOP = 25
"""int: Functions listed by `--profile`."""

TABLES_LISTED = 8
"""int: The run's largest tables the report lists (all of them go in the JSON)."""

STARTUP_RUNS = 3
"""int: How many fresh processes the worker start-up cost is measured over (the median is reported)."""

_STARTUP_PROBE = """
import json, sys, time
began = time.perf_counter()
from planetgen.generation import run_common, run_galaxy
from planetgen.db import store
imported = time.perf_counter()
store.get_connection(store.MySQLConfig(database=sys.argv[1])).close()
print(json.dumps({"import_seconds": imported - began, "connect_seconds": time.perf_counter() - imported}))
"""
"""str: What a worker does before its first task, timed in a process of its own: the imports of a sector worker and the
first database connection."""


def add_benchmark_arguments(parser):
    parser.add_argument("--sectors", type=int, default=DEFAULT_SECTORS, metavar="N",
                        help=f"Sectors the fill makes (default {DEFAULT_SECTORS}).")
    parser.add_argument("--shell", type=int, default=DEFAULT_SHELL, metavar="I",
                        help=f"The shell the fill starts in (default {DEFAULT_SHELL}); a lower one is denser.")
    parser.add_argument("--workers", default="1", metavar="N[,N...]",
                        help="Worker counts to run at, one run each (default 1).")
    parser.add_argument("--seed", default=DEFAULT_SEED, metavar="HEX", help="The galaxy seed (32 hex digits).")
    parser.add_argument("--max-shell", type=int, default=12, metavar="I",
                        help="The plan's outermost shell (default 12, a small galaxy).")
    parser.add_argument("--profile", action="store_true",
                        help="Also profile one single-worker fill and list its costliest functions.")
    parser.add_argument("--keep", action="store_true", help="Leave the scratch databases in place.")
    parser.add_argument("--output", metavar="FILE", help="Write the results to FILE as JSON.")
    store.add_mysql_connection_args(parser)
    from planetgen.cli import benchmark_suite
    benchmark_suite.add_suite_arguments(parser)


def _workers(text):
    try:
        counts = [int(part) for part in str(text).split(",") if part.strip()]
    except ValueError:
        counts = []
    if not counts or any(count < 1 for count in counts):
        raise SystemExit("planetgen benchmark: --workers wants positive whole numbers, like 1,2,4.")
    return counts


def _admin(config):
    return login.connect(host=config.host, port=config.port, user=config.user, password=config.password,
                           autocommit=True)


def _status(connection):
    with connection.cursor() as cur:
        cur.execute("SHOW GLOBAL STATUS")
        rows = {name: value for name, value in cur.fetchall()}
    return {name: int(rows.get(name, 0)) for name in STATUS_VARIABLES}


def _buffer_pool_bytes(connection):
    """The server's `innodb_buffer_pool_size`, in bytes (None where the server doesn't have the variable)."""
    with connection.cursor() as cur:
        cur.execute("SELECT @@innodb_buffer_pool_size")
        row = cur.fetchone()
    return None if row is None or row[0] is None else int(row[0])


def _table_sizes(connection, database):
    """Every table of `database` as `{"table", "rows", "data_bytes", "index_bytes"}`, the largest first. `rows` is
    the server's estimate. MySQL 8 caches these figures for a day unless told not to."""
    with connection.cursor() as cur:
        cur.execute("SET SESSION information_schema_stats_expiry = 0")
        cur.execute("SELECT table_name, table_rows, data_length, index_length FROM information_schema.tables "
                    "WHERE table_schema = %s AND table_type = 'BASE TABLE'", (database,))
        rows = cur.fetchall()
    tables = [{"table": name, "rows": int(count or 0), "data_bytes": int(data or 0), "index_bytes": int(index or 0)}
              for name, count, data, index in rows]
    tables.sort(key=lambda table: (-(table["data_bytes"] + table["index_bytes"]), table["table"]))
    return tables


def _startup_seconds(config, database, env, runs=STARTUP_RUNS):
    """The worker start-up cost: the median over `runs` fresh Python processes of the seconds each took to import
    what a sector worker imports and to open its first connection to `database`, plus the whole process's wall
    time, as `{"import_seconds", "connect_seconds", "process_seconds"}`."""
    found = []
    for _ in range(runs):
        started = time.perf_counter()
        done = subprocess.run([sys.executable, "-c", _STARTUP_PROBE, database], env=env, stdin=subprocess.DEVNULL,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        seconds = time.perf_counter() - started
        if done.returncode:
            raise SystemExit(f"planetgen benchmark: the start-up probe failed ({done.returncode}):\n{done.stdout[-2000:]}")
        found.append({**json.loads(done.stdout.strip().splitlines()[-1]), "process_seconds": seconds})
    return {name: sorted(run[name] for run in found)[len(found) // 2]
            for name in ("import_seconds", "connect_seconds", "process_seconds")}


def _child_env(config, control, workers):
    env = dict(os.environ)
    env.update({"PLANETGEN_MYSQL_HOST": config.host, "PLANETGEN_MYSQL_PORT": str(config.port),
                "PLANETGEN_MYSQL_USER": config.user, "PLANETGEN_MYSQL_PASSWORD": config.password or "",
                "PLANETGEN_CONTROL_DATABASE": control, "PLANETGEN_WORKERS": str(workers),
                "PLANETGEN_GENERATION_STATS": "1"})
    return env


def _run(argv, env, label):
    started = time.perf_counter()
    done = subprocess.run(argv, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True)
    seconds = time.perf_counter() - started
    if done.returncode:
        raise SystemExit(f"planetgen benchmark: {label} failed ({done.returncode}):\n{done.stdout[-2000:]}")
    return seconds


def _make_control_schema(config, control):
    """Builds the control schema in the scratch control database, so the stage rows have somewhere to go."""
    previous = os.environ.get(store.CONTROL_DB_ENV_VAR)
    os.environ[store.CONTROL_DB_ENV_VAR] = control
    try:
        store.get_control_connection(store.control_mysql_config(config), ensure_schema=True).close()
    finally:
        if previous is None:
            os.environ.pop(store.CONTROL_DB_ENV_VAR, None)
        else:
            os.environ[store.CONTROL_DB_ENV_VAR] = previous


def _stages(config, control, database):
    conn = login.connect(host=config.host, port=config.port, user=config.user, password=config.password,
                           database=control, autocommit=True, cursorclass=pymysql.cursors.DictCursor)
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT command, stage_n, stage_total, stage_key, label, skipped, seconds, metrics "
                        "FROM generation_stage_runs WHERE database_name = %s ORDER BY id", (database,))
            rows = cur.fetchall()
    finally:
        conn.close()
    stages = []
    for row in rows:
        stages.append({"command": row["command"], "key": row["stage_key"], "label": row["label"],
                       "skipped": bool(row["skipped"]), "seconds": float(row["seconds"]),
                       "metrics": json.loads(row["metrics"] or "{}")})
    return stages


def _profile(config, control, args, workdir):
    """One single-worker fill under cProfile, in a database of its own; returns `(seconds, [(cumulative, calls,
    function)])` for the project's own functions."""
    database = f"planetgen_bench_profile_{uuid.uuid4().hex[:8]}"
    admin = _admin(config)
    with admin.cursor() as cur:
        cur.execute(f"CREATE DATABASE `{database}`")
    try:
        env = _child_env(config, control, 1)
        base = [sys.executable, "-m", "planetgen.cli.generate"]
        conn_args = ["--mysql-host", config.host, "--mysql-port", str(config.port), "--mysql-user", config.user,
                     "--mysql-database", database]
        for command in _plan_commands(args):
            _run(base + command + conn_args, env, "the plan")
        out = os.path.join(workdir, "fill.prof")
        seconds = _run([sys.executable, "-m", "cProfile", "-o", out, "-m", "planetgen.cli.generate",
                        *_fill_args(args), *conn_args], env, "the profiled fill")
        stats = pstats.Stats(out)
    finally:
        with admin.cursor() as cur:
            cur.execute(f"DROP DATABASE IF EXISTS `{database}`")
        admin.close()
    rows = []
    for (filename, line, name), (_cc, calls, _tt, cumulative, _callers) in stats.stats.items():
        if "planetgen" in filename and "/tests/" not in filename:
            short = filename.split("planetgen/", 1)[-1]
            rows.append((cumulative, calls, f"{short}:{line} {name}"))
    rows.sort(reverse=True)
    return seconds, rows[:PROFILE_TOP]


def _plan_commands(args):
    """The commands of the plan: the skeleton, then the scatters on it, each its own command (PERF.121)."""
    from planetgen.generation import run_common
    floor = ["--bright-star-min-luminosity", "100000"]
    return [["plan", "--seed", args.seed, f"--edge={args.max_shell * run_common._edge_pc():g}pc"],
            ["scatter", "--mass", *floor], ["scatter", "--luminosity", *floor], ["scatter", "--phenomena"]]


def _fill_args(args):
    return ["galaxy", "--shell", str(args.shell), "--limit", str(args.sectors), "--yes"]


def run_benchmark(args):
    """`planetgen benchmark`: runs, prints and (with `--output`) saves the report."""
    if getattr(args, "suite", False):
        from planetgen.cli import benchmark_suite
        return benchmark_suite.run_suite(args)
    config = store.mysql_config_from_args(args)
    counts = _workers(args.workers)
    results = []
    admin = _admin(config)
    suffix = uuid.uuid4().hex[:8]
    control = f"planetgen_bench_control_{suffix}"
    scratch = []
    try:
        with admin.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{control}`")
        _make_control_schema(config, control)
        for workers in counts:
            database = f"planetgen_bench_{suffix}_w{workers}"
            with admin.cursor() as cur:
                cur.execute(f"CREATE DATABASE `{database}`")
            scratch.append(database)
            env = _child_env(config, control, workers)
            conn_args = ["--mysql-host", config.host, "--mysql-port", str(config.port), "--mysql-user", config.user,
                         "--mysql-database", database]
            base = [sys.executable, "-m", "planetgen.cli.generate"]
            before = _status(admin)
            plan_seconds = sum(_run(base + command + conn_args, env, "the plan") for command in _plan_commands(args))
            fill_seconds = _run(base + _fill_args(args) + conn_args, env, "the fill")
            after = _status(admin)
            total = plan_seconds + fill_seconds
            results.append({
                "workers": workers, "plan_seconds": plan_seconds, "fill_seconds": fill_seconds,
                "stages": _stages(config, control, database),
                "database": {name: after[name] - before[name] for name in STATUS_VARIABLES},
                "buffer_pool_bytes": _buffer_pool_bytes(admin),
                "tables": _table_sizes(admin, database),
                "worker_startup": _startup_seconds(config, database, env),
                "total_seconds": total,
            })
        profile = None
        if args.profile:
            with tempfile.TemporaryDirectory() as workdir:
                profile = _profile(config, control, args, workdir)
    finally:
        if not args.keep:
            with admin.cursor() as cur:
                for database in scratch:
                    cur.execute(f"DROP DATABASE IF EXISTS `{database}`")
                cur.execute(f"DROP DATABASE IF EXISTS `{control}`")
        admin.close()
    print(format_report(args, results, profile))
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump({"sectors": args.sectors, "shell": args.shell, "seed": args.seed, "runs": results,
                       "profile": None if profile is None else {"seconds": profile[0], "top": profile[1]}},
                      f, indent=2)


def _megabytes(count):
    return f"{count / 1048576:,.1f} MB"


def _setup_lines(run):
    """The report lines about the machine a run was made on: the buffer pool, the largest tables and the worker
    start-up cost (nothing for a result saved before they were recorded)."""
    lines = []
    if run.get("buffer_pool_bytes") is not None:
        lines.append(f"  innodb_buffer_pool_size: {_megabytes(run['buffer_pool_bytes'])}")
    tables = run.get("tables") or []
    if tables:
        total = sum(table["data_bytes"] + table["index_bytes"] for table in tables)
        lines.append(f"  tables: {len(tables)}, {_megabytes(total)} in all; the largest:")
        for table in tables[:TABLES_LISTED]:
            lines.append(f"    {table['table'][:36]:36} {table['rows']:>10,} rows  data {_megabytes(table['data_bytes']):>11}"
                         f"  indexes {_megabytes(table['index_bytes']):>11}")
    startup = run.get("worker_startup")
    if startup:
        lines.append(f"  worker start-up: {startup['process_seconds']:.2f} s a process (imports {startup['import_seconds']:.2f} s, "
                     f"first connection {startup['connect_seconds']:.2f} s), median of {STARTUP_RUNS}")
    return lines


def format_report(args, results, profile):
    """The report as text (see this module's docstring)."""
    lines = [f"Benchmark: plan to shell {args.max_shell}, then {args.sectors} sector(s) from shell {args.shell}, "
             f"seed {args.seed}."]
    for run in results:
        total = run["total_seconds"] or 1.0
        lines += ["", f"{run['workers']} worker(s): plan {run['plan_seconds']:.1f} s, fill {run['fill_seconds']:.1f} s, "
                      f"total {run['total_seconds']:.1f} s", "  stage                                         seconds   share"]
        for stage in run["stages"]:
            if stage["skipped"]:
                continue
            notes = ", ".join(f"{key} {value}" for key, value in sorted(stage["metrics"].items()))
            lines.append(f"  {stage['command'][:6]:6} {stage['label'][:38]:38} {stage['seconds']:8.2f} "
                         f"{100 * stage['seconds'] / total:6.1f}%  {notes}")
        counters = run["database"]
        per = max(run["total_seconds"], 1e-9)
        lines.append("  database: " + ", ".join(f"{name} {value:,} ({value / per:,.0f}/s)"
                                                for name, value in counters.items() if value))
        lines += _setup_lines(run)
    if profile is not None:
        seconds, rows = profile
        lines += ["", f"Profile of one single-worker fill ({seconds:.1f} s with the profiler on): "
                      f"cumulative seconds, calls, function"]
        lines += [f"  {cumulative:8.2f} {calls:10,d}  {name}" for cumulative, calls, name in rows]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Times a small galaxy from the plan to the finished fill.")
    add_benchmark_arguments(parser)
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)
    run_benchmark(args)


if __name__ == "__main__":
    main()
