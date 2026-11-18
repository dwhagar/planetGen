# planetgen/cli/smoke.py

"""
`python -m planetgen.cli.smoke`: the nightly New galaxy smoke run, with the memory of every step (PERF.183).

It makes a scratch database and control database, builds the same job the Generate page's New galaxy button builds
(`generate_page.build_job`: reset, plan, sectors around a random start of `tuning.SMOKE_RADIUS_PC`, link, the three
scatters, the bright-star backfill and the sector paths), and runs it with the job runner, which watches the memory of
MySQL, the planetgen workers, the step and the host during each step (`util/memwatch.py`). It then prints one line per
step and fails (exit 1) when any step peaked above a budget (`tuning.SMOKE_BUDGET_MB`, or `--budget GROUP=MB`), which is
how a change that makes a step hold far more memory is caught the night it merges. Exit 2: a step failed, so there is
no memory to judge. The scratch databases are dropped at the end (`--keep` leaves them).

The "workers" figure counts every planetgen worker on the host, the site's own included, as the host's memory does:
run it on a quiet host. The result also goes to `--output` as JSON. A systemd timer for it is in
`examples/maintenance/planetgen-smoke.timer`.

    python -m planetgen.cli.smoke [--budget GROUP=MB ...] [--radius-pc N] [--output FILE] [--keep]
"""

import argparse
import json
import os
import sys
import tempfile
import time
import uuid

from planetgen import tuning
from planetgen.cli import benchmark
from planetgen.db import login, store
from planetgen.queue import pagerun
from planetgen.util import log
from planetgen.util.format import format_duration_seconds, format_megabytes
from planetgen.web import generate_page, job_runner, jobs

GROUPS = ("mysql", "workers", "step", "host")
"""tuple: The groups a budget is set for."""


def add_smoke_arguments(parser):
    """Adds the smoke run's options (and the MySQL connection args)."""
    parser.add_argument("--budget", action="append", default=[], metavar="GROUP=MB",
                        help=f"Override a budget in megabytes (groups: {', '.join(GROUPS)}); repeatable.")
    parser.add_argument("--radius-pc", type=float, default=tuning.SMOKE_RADIUS_PC,
                        help=f"The radius of the neighborhood to fill (default {tuning.SMOKE_RADIUS_PC:g} pc).")
    parser.add_argument("--workers", type=int, default=2, help="Worker processes for the steps (default 2).")
    parser.add_argument("--output", metavar="FILE", help="Also write the report to FILE as JSON.")
    parser.add_argument("--keep", action="store_true", help="Leave the scratch databases and the job folder.")
    store.add_mysql_connection_args(parser)


def parse_budgets(items):
    """
    The budgets: `tuning.SMOKE_BUDGET_MB` with each `GROUP=MB` of `items` laid over it.

    Raises:
        ValueError: An item that is not `GROUP=MB` with a known group and a positive number.
    """
    budgets = dict(tuning.SMOKE_BUDGET_MB)
    for item in items:
        group, _eq, number = item.partition("=")
        try:
            value = float(number)
        except ValueError:
            value = 0.0
        if group not in GROUPS or value <= 0:
            raise ValueError(f"'{item}' is not GROUP=MB (groups: {', '.join(GROUPS)}).")
        budgets[group] = value
    return budgets


def budget_of(label, group, budgets):
    """The budget in MB of `group` in the step `label`: the step's own (`tuning.SMOKE_STEP_BUDGET_MB`) unless the
    group's budget was set on the command line (`budgets` differs from the default)."""
    if budgets[group] == tuning.SMOKE_BUDGET_MB[group]:
        return tuning.SMOKE_STEP_BUDGET_MB.get(label, {}).get(group, budgets[group])
    return budgets[group]


def over_budget(entries, budgets):
    """
    What peaked above a budget: `[(label, group, peak MB, budget MB)]` over the job's `memory.json` `entries`.
    """
    over = []
    for entry in entries:
        for group in GROUPS:
            peak = entry["peaks"].get(group, 0.0)
            budget = budget_of(entry["label"], group, budgets)
            if peak > budget:
                over.append((entry["label"], group, peak, budget))
    return over


def format_report(entries, budgets, over):
    """The report: one line per step with its peaks, then the verdict."""
    lines = ["Memory peaks of the New galaxy smoke run, in MB (budget in brackets):"]
    for entry in entries:
        peaks = entry["peaks"]
        cells = ", ".join(f"{group} {peaks.get(group, 0.0):,.0f} ({budget_of(entry['label'], group, budgets):,.0f})"
                          for group in GROUPS)
        lines.append(f"  {entry['label']} [{format_duration_seconds(entry['seconds'])}]: {cells}")
    if over:
        lines.append("")
        lines.append("OVER BUDGET:")
        for label, group, peak, budget in over:
            lines.append(f"  {label}: {group} peaked at {format_megabytes(peak)}, the budget is {format_megabytes(budget)}.")
    else:
        lines.append("")
        lines.append("Every step stayed within its budget.")
    return "\n".join(lines)


def _form(database, radius_pc):
    return {"confirm": database, "radius_pc": f"{radius_pc:g}", "neighborhoods": "1", "disregard_missing_data": "1"}


def run_smoke(args):
    """
    Runs the smoke job in scratch databases and judges its memory.

    Returns:
        int: 0 within every budget, 1 over a budget, 2 a step failed.
    """
    budgets = parse_budgets(args.budget)
    config = store.mysql_config_from_args(args)
    suffix = uuid.uuid4().hex[:8]
    control, database = f"planetgen_smoke_control_{suffix}", f"planetgen_smoke_{suffix}"
    admin = login.connect(host=config.host, port=config.port, user=config.user, password=config.password,
                          autocommit=True)
    root = tempfile.mkdtemp(prefix="planetgen-smoke-")
    try:
        with admin.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{control}`")
            cur.execute(f"CREATE DATABASE `{database}`")
        benchmark._make_control_schema(config, control)
        scratch = store.MySQLConfig(host=config.host, port=config.port, user=config.user, password=config.password,
                                    database=database)
        _kind, title, steps = generate_page.build_job("new_galaxy", _form(database, args.radius_pc), database)
        env = jobs.mysql_env(scratch, database)
        env.update({"PLANETGEN_CONTROL_DATABASE": control, "PLANETGEN_WORKERS": str(args.workers),
                    "PLANETGEN_WARM_ON_CHANGE": "0"})    # no cache warm-up of a database about to be dropped
        job_id = f"smoke-{suffix}"
        job_dir = os.path.join(root, job_id)
        os.makedirs(job_dir)
        os.environ.update(env)          # this process runs the job's steps itself: the ledger it writes to and every step see the scratch databases
        control_connection = pagerun.connect()
        try:
            pagerun.start(control_connection, {"id": job_id, "kind": "new_galaxy", "title": title, "admin": "smoke", "database": database, "created_at": time.time(),
                                               "cwd": jobs.REPO_DIR, "root": root, "draft": []}, steps)
            code = job_runner.run_to_end(job_id)
            found = pagerun.read(control_connection, job_id)
        finally:
            control_connection.close()
        entries = _entries(job_dir)
        if code != 0:
            log.error(f"The smoke job failed: {(found or {}).get('error')}\nIts log: {os.path.join(job_dir, 'output.log')}")
            args.keep = True
            return 2
        over = over_budget(entries, budgets)
        print(format_report(entries, budgets, over))
        if args.output:
            with open(args.output, "w", encoding="utf-8") as f:
                json.dump({"budgets": budgets, "steps": entries, "over": over}, f, indent=2)
        return 1 if over else 0
    finally:
        if not args.keep:
            with admin.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{database}`")
                cur.execute(f"DROP DATABASE IF EXISTS `{control}`")
            import shutil
            shutil.rmtree(root, ignore_errors=True)
        admin.close()


def _entries(job_dir):
    """The job's `memory.json` entries, `[]` when there are none."""
    try:
        with open(os.path.join(job_dir, job_runner.MEMORY_NAME), "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return []


def main(argv=None):
    parser = argparse.ArgumentParser(description="Runs a small New galaxy and judges the memory of every step.")
    add_smoke_arguments(parser)
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)
    try:
        parse_budgets(args.budget)
    except ValueError as exc:
        parser.error(str(exc))
    return run_smoke(args)


if __name__ == "__main__":
    sys.exit(main())
