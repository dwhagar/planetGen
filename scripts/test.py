#!/usr/bin/env python3
"""
Picks and runs a slice of the test suite (TEST.139, docs/design/test-suite-rebuild.md).

    scripts/test.py quick                 not database, browser or heavy tests (the fast loop)
    scripts/test.py area web db           the tests of those areas (src/tests/areas.py)
    scripts/test.py browser               the headless-browser tests alone
    scripts/test.py changed [--base REF]  only the tests the changes since REF can affect (a fast loop, not a gate)
    scripts/test.py budget                the big-galaxy query budget (opt-in, builds 2,000,000 systems once; not a gate)
    scripts/test.py full                  everything: the merge gate; also refreshes the map `changed` reads
    scripts/test.py full --check          ... and lists failing tests `changed` would not have selected
    scripts/test.py full -- -x -k login   anything after `--` goes to pytest as it is

Every command except `quick` and `browser` (`changed` and `full` included) runs in two phases: the tests that
need no browser on every core, then the browser tests on `--browser-workers`
workers (default 2). Software WebGL needs the CPU the other workers would
have, and running both at once is what timed the browser tests out (TEST.126).

Failures: every failing test's full report (traceback, worker, machine load) is kept in `.test-failures.jsonl`,
and when the run ends the failed tests are run again alone, on one worker, each with the verdict "passes alone"
(it only failed under load or beside other tests) or "fails alone" (a real failure). `--no-rerun` skips that.

The first line printed is the start time and how long this command took last
time on this checkout; the last line is `PASS` or `FAIL` with the duration and the counts, or `INTERRUPTED` when the
run was stopped (Ctrl-C, or a stop from outside) -- a run whose output ends with none of the three was killed
without a chance to say so. The duration is stored in `.test-runs.json` (not tracked by git) for the next run.

Needs the same setup as plain pytest (docs/testing.md): the `test` extra, a
MySQL or MariaDB server and Redis through the PLANETGEN_TEST_* variables.
"""

import argparse
import datetime
import json
import os
import re
import signal
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_FILE = os.path.join(ROOT, ".test-runs.json")
FAILURES_FILE = os.path.join(ROOT, ".test-failures.jsonl")
PROGRESS_DIR = os.path.join(ROOT, ".test-progress")
KEEP_PROGRESS_SECONDS = 24 * 3600
AREAS = ("web", "api", "db", "generation", "physics", "names", "wiki", "ops", "other")


def build_phases(command, areas, workers, browser_workers):
    """The pytest runs for a command: a list of (label, marker expression, worker count)."""
    if command == "quick":
        return [("quick", "not db and not browser and not heavy", workers)]
    if command == "browser":
        return [("browser", "browser", browser_workers)]
    if command == "budget":
        return [("budget", "budget", 1)]
    selected = ""
    if command == "area":
        selected = " or ".join("area_" + area for area in areas)
    parts = []
    if selected:
        parts = [("non-browser", f"({selected}) and not browser", workers)]
        if "web" in areas:  # the browser tests are all `test_web_*` files
            parts.append(("browser", f"({selected}) and browser", browser_workers))
    else:
        parts = [("non-browser", "not browser", workers), ("browser", "browser", browser_workers)]
    return parts


def depmap_options(command, base, check):
    """The dependency-map options for a command: `changed` selects, `full` refreshes the map (and, with `--check`,
    lists the failures `changed` would have missed). Others leave the map alone."""
    if command == "changed":
        return ["--changed", base]
    if command == "full":
        return ["--record-deps"] + (["--deps-check", base] if check else [])
    return []


def pytest_command(marker, workers, extra):
    command = [sys.executable, "-m", "pytest", "-m", marker, "-p", "no:cacheprovider", "--dist", "worksteal"]
    if workers and workers > 1:
        command += ["-n", str(workers)]
    return command + list(extra)


SUMMARY_COUNT = re.compile(r"(\d+) (passed|failed|errors?|skipped|xfailed|xpassed|deselected)")


def summary_counts(line):
    """{word: count} from a pytest summary line such as `=== 3 failed, 15419 passed, 22 skipped in 1081.44s ===`,
    or {} when the line is not one."""
    if not (line.startswith("=") and re.search(r" in [\d.]+s", line)):
        return {}
    counts = {}
    for number, word in SUMMARY_COUNT.findall(line):
        counts["errors" if word.startswith("error") else word] = int(number)
    return counts


def run_phase(command, environment):
    """Runs one pytest command, passing its output through line by line. Returns (exit code, counts from its
    summary line). If this process is interrupted the child is stopped first and the interrupt goes on."""
    counts = {}
    child = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, errors="replace", bufsize=1)
    try:
        for line in child.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            counts = summary_counts(line.strip()) or counts
        return child.wait(), counts
    except BaseException:
        child.terminate()
        try:
            child.wait(timeout=30)
        except subprocess.TimeoutExpired:
            child.kill()
        raise


def describe_counts(counts):
    order = ("passed", "failed", "errors", "skipped", "xfailed", "xpassed", "deselected")
    return ", ".join(f"{counts[word]} {word}" for word in order if counts.get(word)) or "no tests ran"


def start_progress_log(command_key, extra, base, directory=None, now=None):
    """Starts the progress log of this run (TEST.150): works out the build id, empties the log of that id if an
    earlier run left one, writes its header line, removes logs older than a day, and returns (id, path, parts).
    `src/tests/progress_log.py` appends a line for every test that finishes."""
    sys.path.insert(0, os.path.join(ROOT, "src"))
    try:
        from tests import buildid
    finally:
        sys.path.pop(0)
    directory = directory or PROGRESS_DIR
    now = now or datetime.datetime.now(datetime.timezone.utc)
    parts = buildid.build_parts(ROOT, command_key, extra, base, now=now)
    identifier = buildid.build_id(parts)
    os.makedirs(directory, exist_ok=True)
    cutoff = time.time() - KEEP_PROGRESS_SECONDS
    for name in os.listdir(directory):
        old = os.path.join(directory, name)
        try:
            if name.endswith(".jsonl") and os.path.getmtime(old) < cutoff:
                os.remove(old)
        except OSError:
            pass
    path = os.path.join(directory, identifier + ".jsonl")
    header = {"kind": "header", "build": identifier, "command": command_key, "parts": parts,
              "started": now.strftime("%Y-%m-%dT%H:%M:%SZ")}
    with open(path, "w", encoding="utf-8") as f:
        f.write(json.dumps(header) + "\n")
    return identifier, path, parts


def load_runs():
    try:
        with open(RUNS_FILE, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def save_run(key, seconds, ok):
    runs = load_runs()
    if not ok and runs.get(key, {}).get("ok"):
        return  # a run that stopped early says nothing about how long the command usually takes
    runs[key] = {"seconds": round(seconds), "ok": ok,
                 "at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")}
    try:
        with open(RUNS_FILE, "w", encoding="utf-8") as f:
            json.dump(runs, f, indent=1)
    except OSError:
        pass


def failed_tests(path=FAILURES_FILE):
    """The distinct node ids in the failure log, in the order they first failed, and each one's first report."""
    found = {}
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                try:
                    record = json.loads(line)
                except ValueError:
                    continue
                found.setdefault(record["nodeid"], record)
    except OSError:
        pass
    return found


def rerun_failures(environment):
    """Runs each failed test again alone (one process, no other tests, no xdist) and prints the verdicts."""
    failed = failed_tests()
    if not failed:
        return
    print(f"--- failures: {len(failed)} test(s) failed; running each again alone", flush=True)
    verdicts = {}
    for nodeid, record in failed.items():
        command = [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-q", "-x", nodeid]
        result = subprocess.run(command, cwd=ROOT, env=dict(environment, PLANETGEN_FAILURE_LOG=os.devnull, PLANETGEN_PROGRESS_LOG=""),
                                capture_output=True, text=True)
        verdicts[nodeid] = "passes alone" if result.returncode == 0 else "FAILS ALONE"
    print("--- verdicts (full reports in .test-failures.jsonl)")
    for nodeid, record in failed.items():
        load = f", load {record['load']}" if record.get("load") is not None else ""
        print(f"  {verdicts[nodeid]:13} {nodeid}  [{record['when']} {record['seconds']}s on {record['worker']}{load}]")


def minutes(seconds):
    return f"{seconds // 60} min {seconds % 60} s"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("quick", "area", "browser", "budget", "changed", "full"))
    parser.add_argument("areas", nargs="*", help="for `area`: " + ", ".join(AREAS))
    parser.add_argument("-n", "--workers", type=int, default=os.cpu_count() or 1,
                        help="workers for the non-browser tests (default: every core)")
    parser.add_argument("--no-rerun", action="store_true", help="do not run the failed tests again alone")
    parser.add_argument("--browser-workers", type=int, default=2, help="workers for the browser tests (default 2)")
    parser.add_argument("--base", default="", help="for `changed` and `--check`: the commit to compare with "
                        "(default: where this branch left origin/main)")
    parser.add_argument("--check", action="store_true",
                        help="for `full`: list the failing tests `changed` would not have selected")
    argv = sys.argv[1:] if argv is None else list(argv)
    extra = []
    if "--" in argv:  # everything after `--` is pytest's, whatever it looks like
        split = argv.index("--")
        argv, extra = argv[:split], argv[split + 1:]
    args = parser.parse_args(argv)
    if args.command == "area":
        unknown = [a for a in args.areas if a not in AREAS]
        if not args.areas or unknown:
            parser.error("`area` needs one or more of: " + ", ".join(AREAS))
    elif args.areas:
        parser.error(f"`{args.command}` takes no areas")
    if args.check and args.command != "full":
        parser.error("`--check` goes with `full`")
    extra = depmap_options(args.command, args.base, args.check) + extra

    key = args.command + (" " + " ".join(args.areas) if args.areas else "")
    previous = load_runs().get(key)
    started = time.time()
    identifier, progress_path, _ = start_progress_log(key, extra, args.base)
    print(f"scripts/test.py {key}: started {datetime.datetime.now(datetime.timezone.utc):%H:%M}Z"
          + (f", last time {minutes(previous['seconds'])}" if previous else ", no earlier run recorded")
          + f", build {identifier}", flush=True)

    status = 0
    totals = {}
    try:
        os.remove(FAILURES_FILE)
    except OSError:
        pass
    environment = dict(os.environ, PLANETGEN_FAILURE_LOG=FAILURES_FILE, PLANETGEN_PROGRESS_LOG=progress_path)
    if args.command == "budget":
        environment["PLANETGEN_BUDGET"] = "1"
    phase = "before the first phase"

    def interrupted(signum, frame):
        raise KeyboardInterrupt(signal.Signals(signum).name)

    previous_handler = signal.signal(signal.SIGTERM, interrupted)  # a stop from outside ends like Ctrl-C: with the INTERRUPTED line
    try:
        for label, marker, workers in build_phases(args.command, args.areas, args.workers, args.browser_workers):
            phase = label
            print(f"--- {label}: -m {marker!r}, {workers} worker(s)", flush=True)
            code, counts = run_phase(pytest_command(marker, workers, extra), environment)
            for word, number in counts.items():
                totals[word] = totals.get(word, 0) + number
            if code not in (0, 5):  # 5: nothing collected in this phase
                status = status or code
        if status and not args.no_rerun:
            rerun_failures(environment)
    except KeyboardInterrupt as stop:
        seconds = round(time.time() - started)
        print(f"scripts/test.py {key}: INTERRUPTED ({stop.args[0] if stop.args else 'Ctrl-C'}) in the {phase} phase "
              f"after {minutes(seconds)}; not finished, {describe_counts(totals)} so far", flush=True)
        return 130
    finally:
        signal.signal(signal.SIGTERM, previous_handler)
    seconds = time.time() - started
    save_run(key, seconds, status == 0)
    print(f"scripts/test.py {key}: {'PASS' if status == 0 else 'FAIL'} in {minutes(round(seconds))}: "
          f"{describe_counts(totals)}", flush=True)
    return status


if __name__ == "__main__":
    sys.exit(main())
