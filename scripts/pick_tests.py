#!/usr/bin/env python3
"""
Picks the tests a PR needs, runs them against a time limit, and re-runs only the files that failed (TEST.722).

    python scripts/pick_tests.py                       list the picked test files for this branch against origin/main
    python scripts/pick_tests.py --run                 run them (default limit 20 minutes), in priority order
    python scripts/pick_tests.py --run --minutes 45    another limit
    python scripts/pick_tests.py --failed              run again only the files that failed in the last run
    python scripts/pick_tests.py --failed-from LOG     the same, from a pytest log (needs `-rf` or the default summary)

The order is the order of trust, and the limit cuts from the bottom:

1. the test files the PR edits or adds;
2. the test files whose name holds the name of a changed module (`test_*<module>*.py`), which is what the gate used by hand;
3. the test files that import a changed module (`import planetgen.a.b`, `from planetgen.a import b`, `from planetgen.a.b import x`),
   fewest imports of the module first, because a file that names only that module is the likelier home of its tests.

`--run` starts a batch only before the limit has passed and prints what it left undone ("left to the end-of-version run"), so a PR
that touches a module everything imports costs the stated time, not the day. A run writes `.pick-last.json` (the files that
failed and the files left undone); `--failed` reads it, so a fix is proved by its failed files, not by the whole set again.
The files listed in `scripts/data/pregate_set.txt` are not picked here: `scripts/pregate.py` runs them.

Reused: `pregate.changed_files` for the diff, `pregate._lines` for the pre-gate set; nothing else is re-implemented.
"""

import argparse
import json
import os
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import pregate  # noqa: E402  (the shared diff and set-file readers)

TESTS = ROOT / "src" / "tests"
LAST_FILE = ROOT / ".pick-last.json"
"""pathlib.Path: What the last `--run` left: the failed files and the files not reached."""
DEFAULT_MINUTES = 20
BATCH = 6
"""int: Test files per pytest process. Small enough that the limit is checked often, large enough that start-up is not the cost."""


def module_names(path):
    """`('a.b.c', 'c')` for the file `a/b/c.py` under the package folder: the dotted name and the base name; `__init__` names its package."""
    parts = pathlib.PurePosixPath(path).with_suffix("").parts
    if parts[:2] != ("src", "planetgen"):
        return None
    dotted = list(parts[1:])
    if dotted[-1] == "__init__":
        dotted.pop()
    if len(dotted) < 2:
        return None
    return ".".join(dotted), dotted[-1]


def imports_module(text, dotted):
    """How many import statements in `text` name the module `dotted` (`import a.b.c`, `from a.b import c`, `from a.b.c import x`)."""
    package, _, leaf = dotted.rpartition(".")
    pattern = (r"^\s*(?:from\s+%(pkg)s\s+import\s+[^\n]*\b%(leaf)s\b|from\s+%(full)s(?:\.\w+)*\s+import\b|import\s+%(full)s\b)"
               % {"pkg": re.escape(package), "leaf": re.escape(leaf), "full": re.escape(dotted)})
    return len(re.findall(pattern, text, flags=re.MULTILINE))


def pick(changed, tests_dir=TESTS, skip=()):
    """The picked test files as `[(repo-relative path, reason)]`, in priority order, each file once, `skip` left out."""
    out, seen = [], set(skip)

    def add(name, reason):
        if name not in seen:
            seen.add(name)
            out.append((name, reason))

    edited = sorted(name for name in changed if re.fullmatch(r"src/tests/test_[^/]+\.py", name))
    for name in edited:
        if (ROOT / name).exists() or tests_dir != TESTS:
            add(name, "edited by the PR")
    modules = [names for names in (module_names(name) for name in changed if name.endswith(".py")) if names]
    files = sorted(tests_dir.glob("test_*.py"))
    prefix = "src/tests/" if tests_dir == TESTS else str(tests_dir) + "/"
    for dotted, leaf in modules:
        for path in files:
            if leaf != "__init__" and leaf in path.stem[len("test_"):]:
                add(prefix + path.name, "name holds " + leaf)
    ranked = []
    for dotted, _leaf in modules:
        for path in files:
            try:
                count = imports_module(path.read_text(encoding="utf-8"), dotted)
            except (OSError, UnicodeDecodeError):
                continue
            if count:
                ranked.append((count, path.name, dotted))
    for count, filename, dotted in sorted(ranked):
        add(prefix + filename, "imports " + dotted)
    return out


def failed_files(log_text):
    """The test files with a `FAILED` or `ERROR` line in a pytest log (`FAILED <test file>::<name> - ...`), sorted."""
    return sorted(set(re.findall(r"^(?:FAILED|ERROR)\s+(src/tests/test_[^\s:]+\.py)", log_text, flags=re.MULTILINE)))


def run_batches(files, deadline, workers="auto", runner=subprocess.run, clock=time.monotonic):
    """Runs `files` in batches until `deadline` (monotonic seconds). Returns `(failed files, files not reached)`."""
    failed, left = [], list(files)
    while left and clock() < deadline:
        batch, left = left[:BATCH], left[BATCH:]
        result = runner([sys.executable, "-m", "pytest", *batch, "-n", str(workers), "-q", "-rf", "-p", "no:cacheprovider", "--no-header"],
                        cwd=ROOT, capture_output=True, text=True)
        print((result.stdout.splitlines() or [""])[-1], flush=True)
        if result.returncode not in (0, 5):
            failed += failed_files(result.stdout) or [name for name in batch]
    return sorted(set(failed)), left


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--base", default="origin/main", help="the branch to compare with (default origin/main)")
    parser.add_argument("--run", action="store_true", help="run the picked files against the limit")
    parser.add_argument("--minutes", type=float, default=DEFAULT_MINUTES, help="the limit for --run (default %d)" % DEFAULT_MINUTES)
    parser.add_argument("--workers", default="auto", help="pytest-xdist workers per batch (default auto)")
    parser.add_argument("--failed", action="store_true", help="run again only the files the last run failed")
    parser.add_argument("--failed-from", metavar="LOG", help="run again only the files that failed in this pytest log")
    args = parser.parse_args(argv)

    if args.failed_from or args.failed:
        if args.failed_from:
            files = failed_files(pathlib.Path(args.failed_from).read_text(encoding="utf-8", errors="replace"))
        else:
            try:
                files = json.loads(LAST_FILE.read_text(encoding="utf-8")).get("failed", [])
            except (OSError, ValueError):
                files = []
        if not files:
            print("no failed files recorded: nothing to run again")
            return 0
        print("running again only the failed files:\n  " + "\n  ".join(files), flush=True)
        failed, left = run_batches(files, float("inf"), args.workers)
        LAST_FILE.write_text(json.dumps({"failed": failed, "left": left}, indent=1), encoding="utf-8")
        print("still failing: " + (", ".join(failed) or "none"))
        return 1 if failed else 0

    skip = set(pregate._lines(pregate.SET_FILE))
    picked = pick(pregate.changed_files(args.base), skip=skip)
    if not args.run:
        print("picked %d test files (the order is the order they run in)" % len(picked))
        for name, reason in picked:
            print("  %-60s %s" % (name, reason))
        return 0
    deadline = time.monotonic() + args.minutes * 60
    failed, left = run_batches([name for name, _ in picked], deadline, args.workers)
    LAST_FILE.write_text(json.dumps({"failed": failed, "left": left}, indent=1), encoding="utf-8")
    if left:
        print("limit of %.0f minutes reached; left to the end-of-version run (%d files):\n  %s" % (args.minutes, len(left), "\n  ".join(left)))
    print("failed files: " + (", ".join(failed) or "none"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
