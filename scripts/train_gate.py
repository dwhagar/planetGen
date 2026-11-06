#!/usr/bin/env python3
"""
The merge train's gate: the test files a train must pass, run one file at a time on one or more machines.

    python scripts/train_gate.py files [--base origin/main] [--out gate-files.txt]
    python scripts/train_gate.py shard --files gate-files.txt --of 3             print the three shards, balanced by last durations
    python scripts/train_gate.py run --files gate-files.txt --shard 2/3 --results gate-results.txt
    python scripts/train_gate.py reds --results gate-results.txt                 reds that the known-red list does not explain
    python scripts/train_gate.py learn --results gate-results.txt                keep each file's duration for the next shard split
    python scripts/train_gate.py known --add tests/test_x.py "red on main since 6.2, TEST.123"

* `files` lists the test files a train's changes reach: the tests the branch changed, those that name a module the branch changed
  (modules that nearly every test reaches are left to the fixed list below), and a fixed list of guard tests (golden, migrations,
  registry, references, wiring). It replaces the per-train script that used to live in a scratchpad.
* `shard` splits a file list over `--of N` machines so that each has about the same total time (longest first onto the lightest
  shard), from `scripts/data/gate_durations.json`; a file it has no time for counts as the median.
* `run` is resumable: a file already named in the results file is skipped, so a restart (the container is reclaimed about every
  three hours) only needs the same command again. Each line is `PASS file seconds` or `FAIL file seconds rc=N last-line`. Every
  shard may write its own results file; `reds` and `learn` read any number of them.
* `reds` prints only the failures that `scripts/data/known_reds.json` does not list (a red that is red on a clean main is known and
  is not the train's), and says which known reds passed (their entry is stale).
* `known` records a clean-main red with the reason; run the gate once from a clean `origin/main` worktree and `reds` from that
  results file with `--record` to make the list.

Run the services first with `scripts/lane_up.sh`. The gate never starts a second pytest on a machine: a second process shares the
machine's Redis and fails the queue tests (RunStopped).
"""

import argparse
import json
import os
import re
import statistics
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
DURATIONS = os.path.join(ROOT, "scripts", "data", "gate_durations.json")
KNOWN = os.path.join(ROOT, "scripts", "data", "known_reds.json")

BIG_MODULES = frozenset((
    "__init__", "tuning", "store", "models", "system", "planet", "planets", "query", "sheet", "build", "units", "draw", "config",
    "constants", "fetch", "settings", "edits", "generate", "catalogue", "components", "equipment", "mathcheck"))
"""frozenset: Modules that almost every test reaches: a change to one does not pull in every file that mentions it."""

GUARD_TESTS = (
    "test_wiring_guard", "test_burn_down_register", "test_references", "test_references_stand_ins", "test_stand_in_values",
    "test_citations_outside_the_registry", "test_american_spelling", "test_db_old_schemas", "test_db_schema_current",
    "test_db_migration_rerun", "test_db_columns_and_indexes", "test_db_column_round_trip", "test_db_check_constraints",
    "test_property_sheet", "test_object_uids", "test_galaxy_dials", "test_galaxy_gen", "test_gen_stored_trees",
    "test_golden_seed", "test_golden_update", "test_kind_tables", "test_fingerprint", "test_storage_checksum",
    "test_plain_words", "test_precision", "test_column_units", "test_api_link_ids", "test_no_own_log_formatting")
"""tuple: Tests that guard what every train changes: golden files, migrations, registries, references, units, wording."""

SLOW_CENSUS = ("test_prevalence", "test_fuzz_system_generation", "test_fuzz_planet_inputs", "test_generation_mix",
               "test_gen_kind_factors", "test_web_browser", "test_js_unit", "test_systems.py")
"""tuple: Name parts of the files the gate leaves to the census run and the end-of-version bug hunt."""

DEFAULT_LIMIT_S = 3400
"""int: Seconds one file may take before the gate gives up on it."""


def _git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def _load(path, default):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return default


def _save(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=1, sort_keys=True)
        handle.write("\n")


def gate_files(base="origin/main"):
    """The test files a train's changes reach, as paths under `src/` (`tests/test_x.py`), sorted."""
    names = _git("diff", "--name-only", base, "HEAD").split()
    chosen = {"tests/" + os.path.basename(n) for n in names if n.startswith("src/tests/test_") and n.endswith(".py")}
    modules = [os.path.basename(n)[:-3] for n in names
               if n.startswith("src/planetgen/") and n.endswith(".py") and "/migrations/" not in n]
    patterns = []
    for module in modules:
        if module in BIG_MODULES:
            continue
        stem = re.escape(module)
        patterns.append(re.compile(rf"(from\s+planetgen[\w.]*\s+import\s+[^\n]*\b{stem}\b|planetgen\.[\w.]*\b{stem}\b|import\s+{stem}\b)"))
    tests_dir = os.path.join(SRC, "tests")
    for name in sorted(os.listdir(tests_dir)):
        if not (name.startswith("test_") and name.endswith(".py")):
            continue
        with open(os.path.join(tests_dir, name), encoding="utf-8", errors="ignore") as handle:
            text = handle.read()
        if any(pattern.search(text) for pattern in patterns):
            chosen.add("tests/" + name)
    chosen.update(f"tests/{name}.py" for name in GUARD_TESTS)
    chosen = {path for path in chosen if os.path.exists(os.path.join(SRC, path))}
    return sorted(path for path in chosen if not any(part in path for part in SLOW_CENSUS))


def shards(files, count, durations=None):
    """`count` lists of files with about the same total time: the longest file goes onto the lightest list each time."""
    durations = _load(DURATIONS, {}) if durations is None else durations
    known = [v for v in durations.values() if v > 0]
    typical = statistics.median(known) if known else 60.0
    weight = {f: float(durations.get(f, typical)) for f in files}
    lists = [[] for _ in range(count)]
    totals = [0.0] * count
    for f in sorted(files, key=lambda x: -weight[x]):
        lightest = totals.index(min(totals))
        lists[lightest].append(f)
        totals[lightest] += weight[f]
    return [sorted(lst, key=lambda x: -weight[x]) for lst in lists], totals


def read_results(paths):
    """`{file: (status, seconds, detail)}` from any number of results files; a later line for a file replaces an earlier one."""
    found = {}
    for path in paths:
        try:
            lines = open(path, encoding="utf-8").read().splitlines()
        except OSError:
            continue
        for line in lines:
            match = re.match(r"^(PASS|FAIL) (\S+)(?: (\d+(?:\.\d+)?))?(?: (.*))?$", line)
            if match:
                found[match.group(2)] = (match.group(1), float(match.group(3) or 0.0), match.group(4) or "")
    return found


def run(files, results, shard=(1, 1), workers=3, limit=DEFAULT_LIMIT_S, logdir=None, say=print):
    """Runs the files of `shard` (k, n) one at a time, appending a line per file to `results`; files already there are skipped."""
    durations = _load(DURATIONS, {})
    mine = shards(files, shard[1], durations)[0][shard[0] - 1] if shard[1] > 1 else list(files)
    logdir = logdir or os.path.join(os.path.dirname(os.path.abspath(results)), "gate-logs")
    os.makedirs(logdir, exist_ok=True)
    done = read_results([results])
    env = dict(os.environ, PYTHONPATH=SRC)
    env.setdefault("PLANETGEN_DATA_DIR", "/tmp/claude/data_base")
    for path in mine:
        if path in done:
            continue
        log = os.path.join(logdir, path.replace("/", "_") + ".log")
        started = time.monotonic()
        with open(log, "w", encoding="utf-8") as handle:
            try:
                code = subprocess.run(
                    [sys.executable, "-m", "pytest", path, "-n", str(workers), "--dist", "worksteal", "--tb=short", "-rfE",
                     "-o", f"timeout={min(limit, 1500)}", "-p", "no:cacheprovider", "-q"],
                    cwd=SRC, env=env, stdout=handle, stderr=subprocess.STDOUT, timeout=limit).returncode
            except subprocess.TimeoutExpired:
                code = 124
                handle.write(f"\nthe gate gave up after {limit} s\n")
        seconds = time.monotonic() - started
        with open(log, encoding="utf-8", errors="replace") as handle:
            last = (handle.read().strip().splitlines() or [""])[-1][:100]
        line = f"PASS {path} {seconds:.0f}" if code in (0, 5) else f"FAIL {path} {seconds:.0f} rc={code} {last}"
        with open(results, "a", encoding="utf-8") as handle:
            handle.write(line + "\n")
        say(line)
    with open(results, "a", encoding="utf-8") as handle:
        handle.write(f"done {shard[0]}/{shard[1]}\n")


def reds(paths, record=False):
    """`(new, known, stale)`: failures the known-red list does not explain, failures it does, and its entries that passed."""
    results = read_results(paths)
    known = _load(KNOWN, {})
    failed = sorted(f for f, (status, _s, _d) in results.items() if status == "FAIL")
    new = [f for f in failed if f not in known]
    explained = [f for f in failed if f in known]
    stale = sorted(f for f in known if f in results and results[f][0] == "PASS")
    if record:
        for f in new:
            known[f] = results[f][2] or "red on a clean main"
        _save(KNOWN, known)
    return new, explained, stale


def learn(paths):
    """Keeps each file's seconds: the new time, or the mean of the old and new one, so one odd run does not move a shard."""
    durations = _load(DURATIONS, {})
    for f, (_status, seconds, _detail) in read_results(paths).items():
        if seconds > 0:
            old = durations.get(f)
            durations[f] = round(seconds if old is None else (old + seconds) / 2.0, 1)
    _save(DURATIONS, durations)
    return len(durations)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("files")
    p.add_argument("--base", default="origin/main")
    p.add_argument("--out")
    p = sub.add_parser("shard")
    p.add_argument("--files", required=True)
    p.add_argument("--of", type=int, required=True)
    p = sub.add_parser("run")
    p.add_argument("--files", required=True)
    p.add_argument("--results", required=True)
    p.add_argument("--shard", default="1/1")
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--limit", type=int, default=DEFAULT_LIMIT_S)
    p = sub.add_parser("reds")
    p.add_argument("--results", nargs="+", required=True)
    p.add_argument("--record", action="store_true", help="add the new reds to the known list (a run from a clean main)")
    p = sub.add_parser("learn")
    p.add_argument("--results", nargs="+", required=True)
    p = sub.add_parser("known")
    p.add_argument("--add", nargs=2, metavar=("FILE", "WHY"))
    p.add_argument("--remove", metavar="FILE")
    args = parser.parse_args(argv)

    if args.command == "files":
        files = gate_files(args.base)
        text = "\n".join(files) + "\n"
        if args.out:
            with open(args.out, "w", encoding="utf-8") as handle:
                handle.write(text)
        else:
            sys.stdout.write(text)
        print(f"{len(files)} files", file=sys.stderr)
    elif args.command == "shard":
        files = open(args.files, encoding="utf-8").read().split()
        lists, totals = shards(files, args.of)
        for index, (lst, total) in enumerate(zip(lists, totals), 1):
            print(f"# shard {index}/{args.of}: {len(lst)} files, about {total / 60:.0f} min")
            print("\n".join(lst))
    elif args.command == "run":
        files = open(args.files, encoding="utf-8").read().split()
        k, n = (int(x) for x in args.shard.split("/"))
        run(files, args.results, (k, n), args.workers, args.limit)
    elif args.command == "reds":
        new, explained, stale = reds(args.results, args.record)
        print(f"new reds: {len(new)}")
        for f in new:
            print("  " + f)
        print(f"known reds that failed again: {len(explained)}")
        if stale:
            print("known reds that passed (take them off the list): " + ", ".join(stale))
        return 1 if new else 0
    elif args.command == "learn":
        print(f"{learn(args.results)} durations kept")
    elif args.command == "known":
        known = _load(KNOWN, {})
        if args.add:
            known[args.add[0]] = args.add[1]
        if args.remove:
            known.pop(args.remove, None)
        _save(KNOWN, known)
        for f, why in sorted(known.items()):
            print(f"{f}: {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
