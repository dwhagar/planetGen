# tests/depmap.py

"""
Which tests a change touches (TEST.139, docs/design/test-suite-rebuild.md).

A pytest plugin (loaded by `conftest.py`) with two jobs:

**Recording** (`--record-deps`, also switched on by `--changed`). While each
test runs, it notes the repository files the test used:

- Python files, through `sys.monitoring` (Python 3.12 and later): the file of
  every function that starts. A function is reported once per test, then
  switched off until the next test begins, so the cost is not measurable.
- Any other file the test opens (templates, static scripts, `schema.sql`,
  fixtures, docs), through an audit hook (`sys.addaudithook`, any Python).
- A test that starts a subprocess or forks is "opaque": what the child runs is
  not seen, so it depends on everything under `src/planetgen/`.

The result is `.pytest-deps.json` at the repository root (ignored by git): a
table of file names and, for every test that passed, the numbers of its files.
A run updates the entries of the tests it ran and keeps the rest. A test that
failed or skipped is not recorded (a short run would look like a small
dependency), so it counts as new to the map and is always selected.

**Selecting** (`--changed [BASE]`). The files that differ from BASE (default:
where this branch left `origin/main`, plus uncommitted and untracked files)
choose the tests to run:

- every test in a changed test file;
- every test whose recorded files include a changed file, or include a file
  that imports a changed file (one import hop, because `sys.monitoring` sees
  functions that ran, not a constant read from a module; test files count, so a
  test that only imports a constant is still caught);
- every opaque test, when anything under `src/planetgen/`, `src/html/` or
  `scripts/` changed;
- every test new to the map;
- every test, when the test setup itself changed (`conftest.py`, the support
  modules, `pytest.ini`, `setup.py`, the lock files, `tests/fixtures/`), or
  when there is no map yet (that run also records one).

`--deps-check [BASE]` is for a full run: it selects nothing, but at the end it
lists the tests that failed although `--changed` would not have run them. That
list is the proof that the selection can be trusted; it should stay empty.

Every shortcut errs towards running more. The full suite stays the merge gate.
"""

import ast
import collections
import json
import os
import subprocess
import sys
import threading

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAP_PATH = os.path.join(ROOT, ".pytest-deps.json")
PROOF_PATH = os.path.join(ROOT, ".pytest-deps-proof.json")
FORMAT = 2

SELECT_EVERYTHING = (
    "src/tests/conftest.py", "src/tests/depmap.py", "src/tests/areas.py", "pytest.ini", "setup.py",
    "requirements.lock", "requirements-server.lock", "config.schema.json",
)
SELECT_EVERYTHING_PREFIXES = ("src/tests/fixtures/",)
SELECT_EVERYTHING_SUFFIX = ("_support.py",)
OPAQUE_PREFIXES = ("src/planetgen/", "src/html/", "scripts/")
OPAQUE_EVENTS = frozenset({"subprocess.Popen", "os.fork", "os.forkpty", "os.posix_spawn", "os.spawn", "os.system",
                           "os.exec"})
IGNORED_PARTS = ("/__pycache__/", "/.git/", "/.pytest_cache/", "/.hypothesis/", "/node_modules/", "/.test-runs.json")

_recording = False
_STATE = {}  # the plugin's state for the hooks that get no config
_current = None  # the set of files the running test has used, or None between tests
_opaque = False
_lock = threading.Lock()
_numbers = {}  # absolute file name -> repository-relative name, or None for a file outside the repository


def _relative(path):
    name = _numbers.get(path, 0)
    if name == 0:
        name = None
        if path.startswith(ROOT + os.sep):
            name = path[len(ROOT) + 1:].replace(os.sep, "/")
            if any(part in "/" + name for part in IGNORED_PARTS) or name.endswith((".pyc", ".json-journal")):
                name = None
        _numbers[path] = name
    return name


def _on_function_start(code, _offset, _disable=getattr(getattr(sys, "monitoring", None), "DISABLE", None)):
    used = _current
    if used is not None:
        name = _relative(code.co_filename)
        if name is not None:
            used.add(name)
    return _disable


def _audit(event, args):
    used = _current
    if used is None:
        return
    if event == "open":
        path = args[0]
        if isinstance(path, str) and not path.startswith(("/proc", "/dev", "/sys")):
            name = _relative(path if os.path.isabs(path) else os.path.abspath(path))
            if name is not None:
                used.add(name)
    elif event in OPAQUE_EVENTS:
        global _opaque
        _opaque = True


_audit_installed = False


def _start_recording():
    """Turn on sys.monitoring and the audit hook; False when this Python cannot (before 3.12)."""
    global _recording, _audit_installed
    monitoring = getattr(sys, "monitoring", None)
    if monitoring is None:
        return False
    for tool in (monitoring.PROFILER_ID, monitoring.OPTIMIZER_ID, 4, 5):
        try:
            monitoring.use_tool_id(tool, "planetgen-depmap")
        except ValueError:
            continue
        monitoring.register_callback(tool, monitoring.events.PY_START, _on_function_start)
        monitoring.set_events(tool, monitoring.events.PY_START)
        break
    else:
        return False
    if not _audit_installed:
        sys.addaudithook(_audit)
        _audit_installed = True
    _recording = True
    return True


def pytest_addoption(parser):
    group = parser.getgroup("depmap", "which tests a change touches (TEST.139)")
    group.addoption("--record-deps", action="store_true", default=False,
                    help="record the files every test uses into .pytest-deps.json")
    group.addoption("--changed", nargs="?", const="", default=None, metavar="BASE",
                    help="run only the tests the changes since BASE can affect (default BASE: where this branch "
                         "left origin/main); with no map yet, runs everything and records one")
    group.addoption("--deps-check", nargs="?", const="", default=None, metavar="BASE",
                    help="run everything, then list the failing tests that --changed would not have selected")
    group.addoption("--deps-file", default=MAP_PATH, help="where the map lives (default .pytest-deps.json)")


# ---------------------------------------------------------------- the map ----

def load_map(path):
    """The stored map as {"files": [names], "tests": {nodeid: {"f": [numbers], "o": bool, "s": seconds}}} or None."""
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or data.get("format") != FORMAT:
        return None
    return data


def merge_map(old, results, file_exists):
    """Fold this run's `results` ({nodeid: (names, opaque, seconds)}) into the stored map.
    Entries of tests that were not run stay, unless their test file is gone (`file_exists(name)` is false)."""
    names = list(old["files"]) if old else []
    number = {name: i for i, name in enumerate(names)}
    tests = dict(old["tests"]) if old else {}
    for nodeid in [n for n in tests if not file_exists(n.split("::", 1)[0])]:
        del tests[nodeid]
    for nodeid, (used, opaque, seconds) in results.items():
        numbers = []
        for name in sorted(used):
            if name not in number:
                number[name] = len(names)
                names.append(name)
            numbers.append(number[name])
        tests[nodeid] = {"f": numbers, "o": bool(opaque), "s": round(seconds, 2)}
    return {"format": FORMAT, "files": names, "tests": tests}


def write_map(path, data):
    temporary = f"{path}.{os.getpid()}.tmp"
    with open(temporary, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(",", ":"))
    os.replace(temporary, path)


# ------------------------------------------------------------ the selection ----

def _git(*args):
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.splitlines() if result.returncode == 0 else None


def changed_files(base=""):
    """The repository-relative names that differ from `base` (default: the merge base with origin/main or main),
    including uncommitted and untracked files. None when git cannot say."""
    candidates = [base] if base else ["origin/main", "main"]
    start = None
    for candidate in candidates:
        found = _git("merge-base", candidate, "HEAD")
        if found:
            start = found[0]
            break
    if start is None:
        return None
    differing = _git("diff", "--name-only", start)
    untracked = _git("ls-files", "--others", "--exclude-standard")
    if differing is None or untracked is None:
        return None
    return set(differing) | set(untracked)


def importers_of(changed):
    """For each changed Python file under src/, the files under src/ (tests included) whose import statements name it
    (one hop). A test that only reads a constant from a module starts none of its functions, so the import is the
    only trace it leaves."""
    names = {}
    packages = {}  # a changed __init__.py runs whenever anything inside its package is imported
    for path in changed:
        if path.startswith("src/") and path.endswith(".py"):
            module = path[len("src/"):-3].replace("/", ".")
            if module.endswith(".__init__"):
                module = module[:-len(".__init__")]
                packages[module + "."] = path
            names[module] = path
    found = collections.defaultdict(set)
    if not names:
        return found
    for directory, _dirs, files in os.walk(os.path.join(ROOT, "src")):
        for file_name in files:
            if not file_name.endswith(".py"):
                continue
            full = os.path.join(directory, file_name)
            try:
                with open(full, encoding="utf-8") as f:
                    tree = ast.parse(f.read())
            except (OSError, SyntaxError, ValueError):
                continue
            relative = os.path.relpath(full, ROOT).replace(os.sep, "/")
            package = os.path.relpath(directory, os.path.join(ROOT, "src")).replace(os.sep, ".")
            for node in ast.walk(tree):
                imported = []
                if isinstance(node, ast.Import):
                    imported = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom):
                    if node.level:
                        base_parts = package.split(".")
                        base_parts = base_parts[:len(base_parts) - (node.level - 1)]
                        module = ".".join(base_parts + ([node.module] if node.module else []))
                    else:
                        module = node.module or ""
                    imported = [module] + [module + "." + alias.name for alias in node.names]
                for module in imported:
                    if module in names:
                        found[names[module]].add(relative)
                    for prefix, path in packages.items():
                        if module.startswith(prefix):
                            found[path].add(relative)
    return found


def select(data, changed, all_nodeids):
    """The node ids to run, as (selected_set, reason), given the stored map and the changed names."""
    if data is None:
        return set(all_nodeids), "no map yet"
    for name in changed:
        if (name in SELECT_EVERYTHING or name.startswith(SELECT_EVERYTHING_PREFIXES)
                or (name.startswith("src/tests/") and name.endswith(SELECT_EVERYTHING_SUFFIX))):
            return set(all_nodeids), f"{name} changed"
    importers = importers_of(changed)
    wanted = set(changed)
    for name in changed:
        wanted |= importers.get(name, set())
    number = {name: i for i, name in enumerate(data["files"])}
    wanted_numbers = {number[name] for name in wanted if name in number}
    opaque_change = any(name.startswith(OPAQUE_PREFIXES) for name in changed)
    changed_test_files = {name for name in changed if name.startswith("src/tests/") and name.endswith(".py")}
    chosen = set()
    for nodeid in all_nodeids:
        entry = data["tests"].get(nodeid)
        file_part = nodeid.split("::", 1)[0]
        if (entry is None or file_part in changed_test_files
                or (entry["o"] and opaque_change) or wanted_numbers.intersection(entry["f"])):
            chosen.add(nodeid)
    return chosen, "selected"


# --------------------------------------------------------------- the plugin ----

class _State:
    def __init__(self):
        self.results = {}  # nodeid -> (files, opaque, seconds) of the tests that passed
        self.bad = set()  # nodeids that failed or skipped in any phase
        self.times = collections.defaultdict(float)
        self.would_select = None  # --deps-check: what --changed would have run
        self.note = ""
        self.workers = {}  # controller under xdist: what the workers sent back
        self.worker_notes = []


def _is_worker(config):
    return hasattr(config, "workerinput")


def pytest_configure(config):
    state = config._depmap = _State()
    _STATE["state"] = state
    config._depmap_recording = False
    if _is_worker(config):
        wanted = bool(config.workerinput.get("depmap_record"))
    else:
        wanted = bool(config.getoption("--record-deps") or config.getoption("--changed") is not None
                      or config.getoption("--deps-check") is not None)
    if wanted:
        config._depmap_recording = _start_recording()
        if not config._depmap_recording and not _is_worker(config):
            state.note = "depmap: recording needs Python 3.12 or later, so the map is not updated"


def pytest_configure_node(node):
    """xdist: tell each worker whether to record."""
    node.workerinput["depmap_record"] = bool(node.config._depmap_recording)


def pytest_collection_modifyitems(config, items):
    """Applies --changed (deselects) or --deps-check (notes). Under xdist every worker collects, so every worker
    makes the same choice from the same map and git state."""
    state = config._depmap
    changed_option = config.getoption("--changed")
    check_option = config.getoption("--deps-check")
    if changed_option is None and check_option is None:
        return
    base = changed_option if changed_option is not None else check_option
    data = load_map(config.getoption("--deps-file"))
    changed = changed_files(base)
    nodeids = [item.nodeid for item in items]
    if changed is None:
        chosen, reason = set(nodeids), "git could not list the changes"
    else:
        chosen, reason = select(data, changed, nodeids)
    if check_option is not None:
        state.would_select = chosen
        state.note = f"depmap: --changed would run {len(chosen)} of {len(items)} tests ({reason})"
        return
    skipped = [item for item in items if item.nodeid not in chosen]
    state.note = f"depmap: --changed selected {len(chosen)} of {len(items)} tests ({reason}"
    if changed is not None:
        state.note += f"; {len(changed)} changed files"
    state.note += ")"
    if data is not None and skipped:
        saved = sum(data["tests"].get(item.nodeid, {}).get("s", 0) for item in skipped)
        state.note += f"; about {round(saved)} s of recorded test time skipped"
    if skipped:
        config.hook.pytest_deselected(items=skipped)
        items[:] = [item for item in items if item.nodeid in chosen]
    if not _is_worker(config):
        reporter = config.pluginmanager.get_plugin("terminalreporter")
        if reporter is not None:
            reporter.write_line(state.note)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    global _current, _opaque
    config = item.config
    if not config._depmap_recording:
        yield
        return
    _current, _opaque = set(), False
    sys.monitoring.restart_events()
    try:
        yield
    finally:
        used, opaque = _current, _opaque
        _current = None
        state = config._depmap
        if item.nodeid not in state.bad:
            state.results[item.nodeid] = (used, opaque, state.times.get(item.nodeid, 0.0))


def pytest_runtest_logreport(report):
    """Adds up each phase's time, and notes any test that did not pass in every phase (a failure or a skip).
    (A `pytest_runtest_makereport` hook of ours never runs: pytest's own is first-result-wins and comes first.)"""
    state = _STATE.get("state")
    if state is None:
        return
    state.times[report.nodeid] += report.duration
    if report.outcome != "passed":
        state.bad.add(report.nodeid)


def pytest_sessionfinish(session, exitstatus):
    config = session.config
    state = config._depmap
    if _is_worker(config):
        output = {"note": state.note}
        if config._depmap_recording:
            output["results"] = {k: [sorted(v[0]), v[1], v[2]] for k, v in state.results.items()}
        if state.would_select is not None:
            output["would_select"] = sorted(state.would_select)
        config.workeroutput["depmap"] = output
        return
    if not config._depmap_recording:
        return
    collected = {**state.results, **state.workers}
    if config.option.collectonly or not collected:
        return
    path = config.getoption("--deps-file")
    write_map(path, merge_map(load_map(path), collected, lambda name: os.path.exists(os.path.join(ROOT, name))))
    state.note = f"depmap: recorded {len(collected)} tests into {os.path.relpath(path, ROOT)}"


def pytest_testnodedown(node, error):
    output = getattr(node, "workeroutput", {}).get("depmap")
    if not output:
        return
    state = node.config._depmap
    for nodeid, (used, opaque, seconds) in output.get("results", {}).items():
        state.workers[nodeid] = (set(used), opaque, seconds)
    if output.get("note") and output["note"] not in state.worker_notes:
        state.worker_notes.append(output["note"])
    if "would_select" in output and state.would_select is None:
        state.would_select = set(output["would_select"])


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    state = config._depmap
    if _is_worker(config):
        return
    notes = [state.note] if state.note else []
    notes += [n for n in state.worker_notes if n not in notes]
    if config.getoption("--deps-check") is not None and state.would_select is not None:
        stats = terminalreporter.stats
        failed = {report.nodeid for key in ("failed", "error") for report in stats.get(key, [])}
        missed = sorted(failed - state.would_select)
        terminalreporter.write_sep("-", "depmap check")
        for note in notes:
            terminalreporter.write_line(note)
        terminalreporter.write_line(f"{len(failed)} tests failed; {len(missed)} of them would not have been selected")
        for nodeid in missed:
            terminalreporter.write_line(f"  MISSED: {nodeid}")
        _record_proof(len(state.would_select), len(failed), missed)
    elif notes:
        terminalreporter.write_sep("-", "depmap")
        for note in notes:
            terminalreporter.write_line(note)


def _record_proof(selected, failed, missed):
    """Appends this run to `.pytest-deps-proof.json`, the running record of whether the selection ever missed a failure."""
    try:
        with open(PROOF_PATH, encoding="utf-8") as f:
            proof = json.load(f)
    except (OSError, ValueError):
        proof = {"runs": []}
    proof["runs"].append({"selected": selected, "failed": failed, "missed": missed})
    try:
        with open(PROOF_PATH, "w", encoding="utf-8") as f:
            json.dump(proof, f, indent=1)
    except OSError:
        pass
