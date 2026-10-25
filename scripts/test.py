#!/usr/bin/env python3
"""
Picks and runs a slice of the test suite (TEST.139, docs/design/operations/test-suite-rebuild.md).

    scripts/test.py quick                 not database, browser or heavy tests (the fast loop)
    scripts/test.py area web db           the tests of those areas (src/tests/areas.py)
    scripts/test.py browser               the headless-browser tests alone
    scripts/test.py changed [--base REF]  only the tests the changes since REF can affect (a fast loop, not a gate)
    scripts/test.py budget                the big-galaxy query budget (opt-in, builds 2,000,000 systems once; not a gate)
    scripts/test.py env                   only the check below, and the repairs it makes; runs no test
    scripts/test.py full                  everything: the end-of-version run (TEST.261, SEC.33, GEN.301), not the gate of an
                                          ordinary merge (that is `changed`, with `--resume` after an interruption); it
                                          also refreshes the map `changed` reads
    scripts/test.py full --check          ... and lists failing tests `changed` would not have selected
    scripts/test.py full -- -x -k login   anything after `--` goes to pytest as it is
    scripts/test.py <command> --global-docs  also the freshness tests of the generated global documents (TODO's and the stamp's)
    scripts/test.py <command> --resume    continue an earlier run: skip the tests that passed in it (same build) or in
                                          an earlier run whose environment is the same and whose files that test
                                          used have not changed since (TEST.208)

Every command except `quick` and `browser` (`changed` and `full` included) runs in three phases (TEST.126, TEST.344):
the tests that need no browser and no real queue worker on every core, then the `real_queue` tests (they start real
worker processes and read the shared queues by name, so two at once kill each other's workers) on one worker, then the
browser tests on `--browser-workers` workers (default 2). Software WebGL needs the CPU the other workers would have,
and running both at once is what timed the browser tests out.

The two parallel phases are also capped by what the machine has free when the phase starts (`cap_workers`): at most one
worker per `--worker-mb` MB of available memory (browser workers: `--browser-worker-mb`) and, for the first phase, no
more than the cores the current load leaves idle, never fewer than one. A lane machine that is also running other
lanes' gates then gets fewer workers instead of timeouts and OOM kills; the phase line says when it capped. `--no-worker-cap` turns the cap off.

Containers (TEST.376): before every run, `src/tests/lane_services.py` loads the `PLANETGEN_TEST_*` variables from
/etc/profile.d/planetgen-test.sh when they are missing, starts MySQL and Redis when they do not answer, purges the
MySQL binlogs when the disk is nearly full and checks `argon2-cffi` and (for the browser tests) `playwright` 1.56.0;
what it cannot repair stops the run before it starts (exit 2, `--no-preflight` skips the check). While a phase runs a
watchdog asks MySQL and Redis every 15 seconds; one that is gone is restarted at once, the phase's pytest is stopped and
started again without the tests that already passed (the progress log, as `--resume` does), up to 5 times a phase, so a
container kill costs the tests that were in flight and not the run. The closing line says how many times it happened.

Failures: every failing test's full report (traceback, worker, machine load) is kept in `.test-failures.jsonl`,
and when the run ends the failed tests are run again alone, on one worker, each with the verdict "passes alone"
(it only failed under load or beside other tests) or "fails alone" (a real failure). `--no-rerun` skips that.

Resuming (TEST.148, TEST.208): every run keeps `.test-progress/<build id>.jsonl`, a line per finished test and, for
a test that passed, the files it used. After a run was killed or stopped, the same command with `--resume` skips the
tests that passed in the same build, and the tests that passed in an earlier run (any command) of the same
environment (packages, Python, MySQL and Redis versions, PLANETGEN_* variables, day) when none of the files that
changed since is one the test used (the rules of `changed`), so a merge of main after a cancelled gate reruns only
what it touched. A different environment reuses nothing and every test runs. It says RESUMING first and RESUMED on
its last line, with how many tests were reused and from which runs. Never automatic, and not part of the merge gate.

The first line printed is the start time and how long this command took last
time on this checkout; the last line is `PASS` or `FAIL` with the duration and the counts, or `INTERRUPTED` when the
run was stopped (Ctrl-C, or a stop from outside) -- a run whose output ends with none of the three was killed
without a chance to say so. The duration is stored in `.test-runs.json` (not tracked by git) for the next run.

Needs the same setup as plain pytest (docs/reference/testing.md): the `test` extra, a
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
AREAS = ("web", "api", "db", "generation", "physics", "names", "ops", "other")


WORKER_MB = 1000
"""int: What one non-browser test worker may need (MB): its Python process, its share of the MySQL server and its
database. An estimate (a test process peaked at 260 to 460 MB), tunable with `--worker-mb`."""

BROWSER_WORKER_MB = 2500
"""int: What one browser test worker may need (MB): its Python process plus Chromium with software WebGL. An estimate,
tunable with `--browser-worker-mb`."""


def available_memory_mb():
    """The memory the kernel says is available (MB), or None where /proc/meminfo is not readable."""
    try:
        with open("/proc/meminfo") as meminfo:
            for line in meminfo:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) // 1024
    except (OSError, ValueError):
        pass
    return None


def cap_workers(requested, available_mb, load, cores, worker_mb):
    """
    The worker count for a phase: `requested`, but no more than the memory available allows (`worker_mb` each) and no
    more than the cores the 1-minute `load` leaves idle, and never fewer than one. A reading that is `None` limits
    nothing.

    Returns:
        tuple: `(workers, why)`: the count and a short reason it was cut, or `""` when it was not.
    """
    requested = max(requested, 1)
    limits = []
    if available_mb is not None:
        limits.append((max(available_mb // worker_mb, 1), f"{available_mb} MB free at {worker_mb} MB a worker"))
    if load is not None and cores:
        limits.append((max(cores - round(load), 1), f"load {load:.1f} on {cores} cores"))
    workers, why = min(limits, key=lambda limit: limit[0], default=(requested, ""))
    return (workers, why) if workers < requested else (requested, "")


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
        parts = [("non-browser", f"({selected}) and not browser and not real_queue", workers),
                 ("real-queue", f"({selected}) and real_queue and not browser", 1)]
        if "web" in areas:  # the browser tests are all `test_web_*` files
            parts.append(("browser", f"({selected}) and browser", browser_workers))
    else:
        parts = [("non-browser", "not browser and not real_queue", workers), ("real-queue", "real_queue and not browser", 1),
                 ("browser", "browser", browser_workers)]
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


def run_phase(command, environment, on_child=None):
    """Runs one pytest command, passing its output through line by line. Returns (exit code, counts from its
    summary line). If this process is interrupted the child is stopped first and the interrupt goes on."""
    counts = {}
    child = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, errors="replace", bufsize=1)
    if on_child:
        on_child(child)
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


def _services():
    sys.path.insert(0, os.path.join(ROOT, "src"))
    try:
        from tests import lane_services
    finally:
        sys.path.pop(0)
    return lane_services


def reap_strays():
    """
    TEST.393: ends the queue workers the run left behind (tests start real ones in sessions of their own, and a test that
    was stopped never tells them to stop; each keeps about 300 MB until it is idle) and returns the sentence that says
    how many and how old, or "" when there were none.
    """
    try:
        sys.path.insert(0, os.path.join(ROOT, "src"))
        try:
            from tests import stray_workers
        finally:
            sys.path.pop(0)
        found = stray_workers.reap_all()
    except Exception as exc:  # noqa: BLE001 -- the report must never fail a run
        return f"no stray worker check ({exc})"
    if not found:
        return ""
    oldest = max(worker["seconds"] for worker in found)
    return (f"{len(found)} queue worker process(es) a test started and did not stop (the oldest had run "
            f"{minutes(oldest)}); a test that starts real workers has to stop them")


MAX_RESTARTS = 5
"""int: How many times one phase is started again after a service was restarted under it."""


class ServicesLost(RuntimeError):
    """A service could not be brought back, or kept dying; the run stops with the reason."""


def passed_nodeids(progress_path):
    """The ids of the tests the progress log records as passed (so far)."""
    return {nodeid for nodeid, line in _progress_modules()[1].read(progress_path)[1].items()
            if line.get("outcome") == "passed"}


def skipping(extra, passed, path):
    """`extra` (pytest options) with `--skip-passed` naming a list of the tests it already named plus `passed`, written
    to `path`."""
    extra = list(extra)
    names = set(passed)
    if "--skip-passed" in extra:
        at = extra.index("--skip-passed")
        try:
            with open(extra[at + 1], encoding="utf-8") as earlier:
                names |= set(json.load(earlier))
        except (OSError, ValueError, IndexError):
            pass
        del extra[at:at + 2]
    with open(path, "w", encoding="utf-8") as skip:
        json.dump(sorted(names), skip)
    return extra + ["--skip-passed", path]


def run_watched(label, marker, workers, extra, environment, progress_path, say=print):
    """
    Runs one phase under the services watchdog. A MySQL or Redis that went away is restarted by the watchdog, which
    stops the phase; it is then started again without the tests that passed so far, at most `MAX_RESTARTS` times.

    Returns:
        tuple: `(exit code, counts, restarts)`, the counts adding the tests that passed in the stopped attempts.

    Raises:
        ServicesLost: A service could not be restarted, or the phase kept being cut short.
    """
    services = _services()
    totals, restarts = {}, []
    while True:
        child = []
        watchdog = services.Watchdog(environment, on_trip=lambda found, started=child: started and started[0].terminate())
        before = len(passed_nodeids(progress_path))
        watchdog.start()
        try:
            code, counts = run_phase(pytest_command(marker, workers, extra), environment, on_child=child.append)
        finally:
            watchdog.stop()
        if not watchdog.tripped.is_set() or counts:     # a phase that printed its summary was over before the stop
            restarts += watchdog.events
            for word, number in counts.items():
                totals[word] = totals.get(word, 0) + number
            return code, totals, restarts
        restarts += watchdog.events
        done = len(passed_nodeids(progress_path)) - before
        totals["passed"] = totals.get("passed", 0) + done
        if watchdog.failed:
            raise ServicesLost(f"{watchdog.what}; the {label} phase cannot go on")
        if len(restarts) > MAX_RESTARTS:
            raise ServicesLost(f"the services were restarted {len(restarts)} times in the {label} phase "
                               f"({'; '.join(restarts)}); something keeps killing them")
        say(f"--- {label}: {watchdog.what}; the phase is stopped and started again without the {done} tests that "
            f"passed in it (restart {len(restarts)} of at most {MAX_RESTARTS})")
        extra = skipping(extra, passed_nodeids(progress_path), progress_path + ".restart.skip")


def _progress_modules():
    sys.path.insert(0, os.path.join(ROOT, "src"))
    try:
        from tests import buildid, progress_log, redis_cleanup, reuse
    finally:
        sys.path.pop(0)
    return buildid, progress_log, redis_cleanup, reuse


def check_without_changes(base):
    """The refusal for `full --check` when nothing differs from `base`, else "". The check lists the failures the
    changed-only selection would have missed, so with no changes it would list every failure (TEST.294: a check on a
    checkout that is already main compares main with itself)."""
    sys.path.insert(0, os.path.join(ROOT, "src"))
    try:
        from tests import depmap
    finally:
        sys.path.pop(0)
    changed = depmap.changed_files(base)
    if changed is None or changed:
        return ""
    return ("no file differs from " + (base or "where this branch left origin/main") + ", so every failure would "
            "read as one the selection misses. Pass --base <the commit before the change>, e.g. the parent of the "
            "merge you are checking.")


def clear_test_redis():
    """Removes what earlier runs left on the test Redis (TEST.152, src/tests/redis_cleanup.py); the sentence to print, or ""."""
    return _progress_modules()[2].clean()[1]


def explain_no_resume(buildid, progress_log, directory, command_key, parts):
    """Why there is nothing to resume from: that no earlier run of `command_key` left a log, or the newest one and the
    parts of the build that differ from it."""
    newest = None
    for name in os.listdir(directory):
        if not name.endswith(".jsonl"):
            continue
        path = os.path.join(directory, name)
        header, _ = progress_log.read(path)
        if header.get("command") == command_key and (newest is None or os.path.getmtime(path) > newest[0]):
            newest = (os.path.getmtime(path), header)
    if newest is None:
        return f"no earlier run of `{command_key}` has left a progress log in {os.path.relpath(directory, ROOT)}"
    header = newest[1]
    changed = buildid.differences(header.get("parts", {}), parts)
    return (f"the last `{command_key}` run (build {header.get('build')}, started {header.get('started')}) was another "
            f"build: {', '.join(changed) or 'no part differs, but it recorded no passed test'} changed since")


class Progress:
    """The progress log of one run (TEST.150) and, when `resume` is asked for, what earlier runs had passed that this
    one may skip (TEST.151, TEST.208). `problem` says why a resume cannot be made; `reused` lists, per earlier run,
    `(build id, start time, tests)`; `note` says what kept tests from being reused."""

    def __init__(self, identifier, path, parts, skip_file=None, passed=0, earlier_started="", problem=None,
                 reused=(), note=""):
        self.identifier, self.path, self.parts = identifier, path, parts
        self.skip_file, self.passed, self.earlier_started, self.problem = skip_file, passed, earlier_started, problem
        self.reused, self.note = list(reused), note


def resume_sources(progress):
    """Where a resumed run's skipped tests passed: `N in the run started T`, and for another build, which one."""
    return "; ".join(f"{count} in the run started {when}" + ("" if build == progress.identifier else f" of build {build}")
                     for build, when, count in progress.reused)


def _write_lines(path, lines, mode="a"):
    with open(path, mode, encoding="utf-8") as f:
        for line in lines:
            f.write(json.dumps(line) + "\n")


def start_progress_log(command_key, extra, base, directory=None, now=None, resume=False):
    """Starts the progress log of this run: works out the build id, removes logs older than a day, and writes
    the header of the log of that id (emptying one an earlier run of the same build left, unless `resume`).
    `src/tests/progress_log.py` appends a line for every test that finishes.

    With `resume`, the node ids to leave out are written to a skip list for pytest (`--skip-passed`): the tests that
    passed in the log of the same build, and the tests that passed in an earlier run of the same environment whose
    files they used have not changed since (`src/tests/reuse.py`, TEST.208). Those are carried into this run's log, so
    a second resume counts them too. When no earlier run of this command left a log, the result's `problem` says so
    and nothing is written; when there are logs but nothing can be reused, the run goes on with every test and `note`
    says what kept them out."""
    buildid, progress_log, _, reuse = _progress_modules()
    directory = directory or PROGRESS_DIR
    now = now or datetime.datetime.now(datetime.timezone.utc)
    parts = buildid.build_parts(ROOT, command_key, extra, base, now=now)
    identifier = buildid.build_id(parts)
    os.makedirs(directory, exist_ok=True)
    cutoff = time.time() - KEEP_PROGRESS_SECONDS
    for name in os.listdir(directory):
        old = os.path.join(directory, name)
        try:
            if name.endswith((".jsonl", ".skip")) and os.path.getmtime(old) < cutoff:
                os.remove(old)
        except OSError:
            pass
    path = os.path.join(directory, identifier + ".jsonl")
    started = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    current = reuse.manifest(ROOT)
    header = {"kind": "header", "build": identifier, "command": command_key, "parts": parts, "started": started,
              "manifest": current}
    if not resume:
        _write_lines(path, [header], "w")
        return Progress(identifier, path, parts)
    earlier, tests = progress_log.read(path)
    same = {nodeid for nodeid, line in tests.items() if line.get("outcome") == "passed"}
    across, _looked, blocked = reuse.reusable(directory, parts, current, skip=path)
    across = {nodeid: found for nodeid, found in across.items() if nodeid not in same}
    if not earlier and not across:
        if not any(progress_log.read(os.path.join(directory, name))[0].get("command") == command_key
                   for name in os.listdir(directory) if name.endswith(".jsonl") and name != identifier + ".jsonl"):
            return Progress(identifier, path, parts, problem=explain_no_resume(buildid, progress_log, directory,
                                                                               command_key, parts))
        why = (f"the environment changed since the earlier runs ({', '.join(blocked)})" if blocked
               else "no earlier pass is left whose files are unchanged")
        _write_lines(path, [header], "w")
        return Progress(identifier, path, parts, note=why)
    if earlier and not same and not across:
        return Progress(identifier, path, parts, problem=f"the log of build {identifier} records no test that passed")
    wanted = sorted(same | set(across))
    skip_file = os.path.join(directory, identifier + ".skip")
    with open(skip_file, "w", encoding="utf-8") as f:
        json.dump(wanted, f)
    carried = []
    for nodeid, found in sorted(across.items()):
        carried += [dict(found["test"], reused_from=found["build"]), dict(found["deps"], reused_from=found["build"])]
    sources = {}
    if same:
        sources[(identifier, earlier.get("started", ""))] = len(same)
    for found in across.values():
        key = (found["build"], found["started"])
        sources[key] = sources.get(key, 0) + 1
    if earlier:
        _write_lines(path, [{"kind": "resume", "at": started, "skipped": len(wanted), "reused": len(across)}] + carried)
    else:
        _write_lines(path, [header, {"kind": "resume", "at": started, "skipped": len(wanted), "reused": len(across)}]
                     + carried, "w")
    reused = [(build, when, count) for (build, when), count in sorted(sources.items(), key=lambda item: item[0][1])]
    return Progress(identifier, path, parts, skip_file, len(wanted), reused[0][1] if reused else "", reused=reused)


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


TIMEOUT_TEXT = "from pytest-timeout"
"""str: What pytest-timeout writes when it stops a test (`Failed: Timeout (>900.0s) from pytest-timeout.`)."""


def timed_out_tests(path=FAILURES_FILE):
    """The node ids in the failure log that pytest-timeout stopped (TEST.153), in the order they failed."""
    return [nodeid for nodeid, record in failed_tests(path).items() if TIMEOUT_TEXT in record.get("text", "")]


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
    parser.add_argument("command", choices=("quick", "area", "browser", "budget", "changed", "full", "env"))
    parser.add_argument("areas", nargs="*", help="for `area`: " + ", ".join(AREAS))
    parser.add_argument("-n", "--workers", type=int, default=os.cpu_count() or 1,
                        help="workers for the non-browser tests (default: every core)")
    parser.add_argument("--no-rerun", action="store_true", help="do not run the failed tests again alone")
    parser.add_argument("--browser-workers", type=int, default=2, help="workers for the browser tests (default 2)")
    parser.add_argument("--worker-mb", type=int, default=WORKER_MB,
                        help=f"memory one non-browser worker may need, for the worker cap (default {WORKER_MB})")
    parser.add_argument("--browser-worker-mb", type=int, default=BROWSER_WORKER_MB,
                        help=f"memory one browser worker may need, for the worker cap (default {BROWSER_WORKER_MB})")
    parser.add_argument("--no-worker-cap", action="store_true",
                        help="run the phases on exactly the workers asked for, whatever memory and load say")
    parser.add_argument("--base", default="", help="for `changed` and `--check`: the commit to compare with "
                        "(default: where this branch left origin/main)")
    parser.add_argument("--resume", action="store_true",
                        help="continue the last run of this same build: skip the tests that passed in it (never "
                        "automatic; the last line says the run was resumed)")
    parser.add_argument("--no-preflight", action="store_true",
                        help="do not check (and repair) MySQL, Redis, the test variables and the packages first, and "
                        "do not watch the servers during the run (TEST.376)")
    parser.add_argument("--check", action="store_true",
                        help="for `full`: list the failing tests `changed` would not have selected")
    parser.add_argument("--global-docs", action="store_true",
                        help="also run the tests that compare the generated global documents with their sources "
                        "(PLANETGEN_GLOBAL_DOCS=1; the TODO thread's and the release stamp's, DOC.69)")
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
    if args.check:
        refusal = check_without_changes(args.base)
        if refusal:
            print(f"scripts/test.py full --check: {refusal} Nothing was run.", flush=True)
            return 2
    if args.command == "env" or not args.no_preflight:
        problems = _services().preflight(os.environ, browser=args.command != "quick",
                                         say=lambda text: print(f"scripts/test.py {args.command}: {text}", flush=True))
        if problems:
            for problem in problems:
                print(f"scripts/test.py {args.command}: PREFLIGHT FAILED: {problem}", flush=True)
            print(f"scripts/test.py {args.command}: Nothing was run (--no-preflight runs anyway).", flush=True)
            return 2
        if args.command == "env":
            print("scripts/test.py env: MySQL, Redis, the test variables and the packages are in order.", flush=True)
            return 0
    if args.global_docs:
        os.environ["PLANETGEN_GLOBAL_DOCS"] = "1"
    extra = depmap_options(args.command, args.base, args.check) + extra

    key = args.command + (" " + " ".join(args.areas) if args.areas else "")
    previous = load_runs().get(key)
    started = time.time()
    progress = start_progress_log(key, extra, args.base, resume=args.resume)
    identifier, progress_path = progress.identifier, progress.path
    if progress.problem:
        print(f"scripts/test.py {key}: cannot resume build {identifier}: {progress.problem}. Nothing was run; start the "
              "run again without --resume.", flush=True)
        return 2
    print(f"scripts/test.py {key}: started {datetime.datetime.now(datetime.timezone.utc):%H:%M}Z"
          + (f", last time {minutes(previous['seconds'])}" if previous else ", no earlier run recorded")
          + f", build {identifier}", flush=True)
    if args.resume and progress.skip_file:
        print(f"scripts/test.py {key}: RESUMING build {identifier}: the {progress.passed} tests that passed earlier are "
              f"skipped ({resume_sources(progress)}), the rest run now", flush=True)
        extra_for_pytest = extra + ["--skip-passed", progress.skip_file]
    elif args.resume:
        print(f"scripts/test.py {key}: RESUMING build {identifier}: nothing from the earlier runs can be reused "
              f"({progress.note}); every test runs", flush=True)
        extra_for_pytest = extra
    else:
        extra_for_pytest = extra
    said = clear_test_redis()
    if said:
        print(f"scripts/test.py {key}: {said}", flush=True)

    status = 0
    totals = {}
    restarts = []
    try:
        os.remove(FAILURES_FILE)
    except OSError:
        pass
    environment = dict(os.environ, PLANETGEN_FAILURE_LOG=FAILURES_FILE, PLANETGEN_PROGRESS_LOG=progress_path)
    if args.command == "budget":
        environment["PLANETGEN_BUDGET"] = "1"
    phase = "before the first phase"
    first_phase = True

    def interrupted(signum, frame):
        raise KeyboardInterrupt(signal.Signals(signum).name)

    previous_handler = signal.signal(signal.SIGTERM, interrupted)  # a stop from outside ends like Ctrl-C: with the INTERRUPTED line
    try:
        for label, marker, workers in build_phases(args.command, args.areas, args.workers, args.browser_workers):
            phase = label
            capped = ""
            if workers > 1 and not args.no_worker_cap:
                # the load at the start is other people's work; later it is this run's own previous phase still fading
                load = os.getloadavg()[0] if first_phase and hasattr(os, "getloadavg") else None
                per_worker = args.browser_worker_mb if label == "browser" else args.worker_mb
                workers, why = cap_workers(workers, available_memory_mb(), load, os.cpu_count() or 1, per_worker)
                capped = f" (capped: {why})" if why else ""
            first_phase = False
            print(f"--- {label}: -m {marker!r}, {workers} worker(s){capped}", flush=True)
            if args.no_preflight:
                code, counts = run_phase(pytest_command(marker, workers, extra_for_pytest), environment)
            else:
                refused = _services().preflight(environment, browser=label == "browser", say=lambda text: print(
                    f"scripts/test.py {key}: {text}", flush=True))   # a server that fell over between two phases
                if refused:
                    raise ServicesLost("; ".join(refused))
                code, counts, found = run_watched(label, marker, workers, extra_for_pytest, environment, progress_path)
                restarts += found
            for word, number in counts.items():
                totals[word] = totals.get(word, 0) + number
            if code not in (0, 5):  # 5: nothing collected in this phase
                status = status or code
        if status and not args.no_rerun:
            if not args.no_preflight:
                _services().preflight(environment, browser=False, say=lambda text: print(
                    f"scripts/test.py {key}: {text}", flush=True))
            rerun_failures(environment)
    except ServicesLost as lost:
        seconds = round(time.time() - started)
        print(f"scripts/test.py {key}: FAIL in {minutes(seconds)}: {lost}. Stopped in the {phase} phase; build "
              f"{identifier} keeps what passed, and the same command with --resume continues it once the "
              "services are back.", flush=True)
        return 3
    except KeyboardInterrupt as stop:
        seconds = round(time.time() - started)
        _, logged = _progress_modules()[1].read(progress_path)
        done = sum(1 for line in logged.values() if line.get("outcome") == "passed")
        print(f"scripts/test.py {key}: INTERRUPTED ({stop.args[0] if stop.args else 'Ctrl-C'}) in the {phase} phase "
              f"after {minutes(seconds)}; not finished, {done} tests had passed (build {identifier}; the same command "
              "with --resume continues it)", flush=True)
        return 130
    finally:
        signal.signal(signal.SIGTERM, previous_handler)
        left = reap_strays()
        if left:
            print(f"scripts/test.py {key}: ended {left}", flush=True)
    seconds = time.time() - started
    if not args.resume:  # a resumed run is short because it skipped tests; it says nothing about the usual time
        save_run(key, seconds, status == 0)
    timed_out = timed_out_tests()
    if timed_out:
        print(f"--- TIMED OUT: {len(timed_out)} test(s) ran past the per-test limit and were stopped:", flush=True)
        for nodeid in timed_out:
            print(f"  {nodeid}", flush=True)
    resumed = (f" (RESUMED build {identifier}: {progress.passed} more passed earlier, {resume_sources(progress)})"
               if args.resume and progress.skip_file else
               f" (RESUMED build {identifier}: 0 more passed earlier, {progress.note})" if args.resume else "")
    print(f"scripts/test.py {key}: {'PASS' if status == 0 else 'FAIL'} in {minutes(round(seconds))}: "
          f"{describe_counts(totals)}" + (f" ({len(timed_out)} TIMED OUT)" if timed_out else "")
          + (f" ({len(restarts)} service restart(s): {'; '.join(restarts)})" if restarts else "") + resumed, flush=True)
    return status


if __name__ == "__main__":
    sys.exit(main())
