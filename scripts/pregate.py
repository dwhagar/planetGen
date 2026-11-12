#!/usr/bin/env python3
"""
The pre-gate quick lane (TEST.720): the checks that turned most PRs red this week, in minutes, before the long gate run.

    python scripts/pregate.py                     the lane against origin/main
    python scripts/pregate.py --base REF          against another base
    python scripts/pregate.py --no-tests          only the static checks (seconds)
    python scripts/pregate.py --workers 2         fewer test workers on a machine that is busy
    python scripts/pregate.py --list              the checks and the test files the lane runs

The author runs it before pushing and the gate runs it before the long run. It needs no database. Steps, each printed as `ok`
or `FAIL`, and the exit status is 1 when any failed:

1. merges: the PR head merges with the base without a conflict (`git merge-tree`); when every conflicted file is a generated
   one, it names `scripts/regenerate_generated.py --resolve` as the repair.
2. release-note: `bump_version.py --check-pr` (the note exists and is well formed).
3. duplicate-definitions: no module-level `def` or `class` defined twice in a changed file (the second silently replaces the first;
   OPS.142 shadowed `_duration` of the Queue page this way). pyflakes' F811 is not used: it flags every pytest fixture import.
4. ruff: undefined names and syntax errors (`F821`, `F822`, `F823`, `E9`) in the changed Python files.
5. generated: `scripts/regenerate_generated.py --check`, where the script exists.
6. tests: the files of `scripts/data/pregate_set.txt` (the reds that kept recurring: spelling, citations, settings, own log
   formatting, constants, docs links and index, wiring guard, architecture map) plus the test files the PR edits that need no
   database, on every free core (with `--train` also the train fast set). Tests that need MySQL, a browser or a real queue worker are the gate's long run.

The Train Manager's train pre-gate shares the lane: `--train` adds `scripts/data/train_fast_set.txt` (slower: it holds database tests,
which are skipped here for want of a database).
"""

import argparse
import ast
import collections
import os
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SET_FILE = ROOT / "scripts" / "data" / "pregate_set.txt"
TRAIN_SET_FILE = ROOT / "scripts" / "data" / "train_fast_set.txt"
GENERATED = ("docs/generated/", "docs/REFERENCES.", "docs/reference/references-data.json", "docs/references/",
             "docs/design/research/stand-in-values.md")
"""tuple: Path prefixes of files a script rebuilds, so a conflict in them is repaired by regenerating, not by hand."""


def _run(args, **kwargs):
    return subprocess.run(args, cwd=ROOT, capture_output=True, text=True, **kwargs)


def _lines(path):
    """The test files a set file lists (blank lines and `#` comments skipped), repo-relative under `src/`."""
    if not path.exists():
        return []
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line if line.startswith("src/") else "src/" + line)
    return out


def changed_files(base):
    """Files changed on this branch since it left `base`, as repo-relative paths (deleted files left out)."""
    result = _run(["git", "diff", "--name-only", "--diff-filter=d", f"{base}...HEAD"])
    return [line for line in result.stdout.splitlines() if line]


def merges_cleanly(base):
    """`(ok, conflicted files)`: whether HEAD merges with `base` without a conflict."""
    result = _run(["git", "merge-tree", "--write-tree", "--name-only", base, "HEAD"])
    if result.returncode == 0:
        return True, []
    names = []
    for line in result.stdout.splitlines()[1:]:
        if not line.strip():
            break
        names.append(line.strip())
    return False, names


def duplicate_definitions(files):
    """`['path: name at lines 91 and 332']` for each module-level def or class defined more than once in a changed file."""
    found = []
    for name in files:
        if not name.endswith(".py"):
            continue
        try:
            tree = ast.parse((ROOT / name).read_text(encoding="utf-8"))
        except (OSError, SyntaxError, UnicodeDecodeError):
            continue
        seen = collections.defaultdict(list)
        for node in tree.body:
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue
            decorated = any(isinstance(d, ast.Attribute) or getattr(d, "id", "") in ("overload", "register") or getattr(getattr(d, "func", None), "id", "") == "register"
                            for d in node.decorator_list)
            if not decorated:
                seen[node.name].append(node.lineno)
        for key, lines in seen.items():
            if len(lines) > 1:
                found.append("%s: %s at lines %s" % (name, key, " and ".join(str(n) for n in lines)))
    return found


def step(label, ok, detail=""):
    print("%-22s %s%s" % (label, "ok" if ok else "FAIL", ("\n" + detail.rstrip()) if detail and not ok else ""))
    return ok


def test_files(changed, train=False):
    """The pytest files this lane runs: the fixed set, with `train` the train fast set, and the PR's own edited test files."""
    wanted = _lines(SET_FILE) + (_lines(TRAIN_SET_FILE) if train else [])
    wanted += [name for name in changed if name.startswith("src/tests/test_") and name.endswith(".py")]
    out, seen = [], set()
    for name in wanted:
        if name not in seen and (ROOT / name).exists():
            seen.add(name)
            out.append(name)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--base", default="origin/main", help="the branch to compare with (default origin/main)")
    parser.add_argument("--no-tests", action="store_true", help="only the static checks")
    parser.add_argument("--workers", default="auto", help="pytest-xdist workers for the tests step (default auto)")
    parser.add_argument("--train", action="store_true", help="also the Train Manager's train fast set (slower: it holds database tests)")
    parser.add_argument("--list", action="store_true", help="print the checks and the test files, run nothing")
    args = parser.parse_args(argv)
    changed = changed_files(args.base)
    tests = test_files(changed, args.train)
    if args.list:
        print("changed files: %d\n%s" % (len(changed), "\n".join(tests)))
        return 0
    results = []

    clean, conflicted = merges_cleanly(args.base)
    hint = ""
    if not clean:
        generated = [name for name in conflicted if name.startswith(GENERATED)]
        hint = "conflicts in: %s" % ", ".join(conflicted)
        if conflicted and len(generated) == len(conflicted):
            hint += "\nall generated: python scripts/regenerate_generated.py --resolve"
    results.append(step("merges with " + args.base, clean, hint))

    note = _run([sys.executable, "scripts/bump_version.py", "--check-pr", args.base])
    results.append(step("release-note", note.returncode == 0, (note.stdout + note.stderr)[-600:]))

    dups = duplicate_definitions(changed)
    results.append(step("duplicate-definitions", not dups, "\n".join(dups)))

    py = [name for name in changed if name.endswith(".py")]
    if py:
        flags = ["check", "--select", "F821,F822,F823,E9", "--no-cache", "--output-format", "concise", *py]
        ruff = _run([sys.executable, "-m", "ruff", *flags])
        if "No module named" in ruff.stderr:
            try:
                ruff = _run(["ruff", *flags])
            except FileNotFoundError:
                ruff = None
        if ruff is None:
            print("%-22s skipped (ruff is not installed: pip install ruff)" % "ruff (names, syntax)")
        else:
            results.append(step("ruff (names, syntax)", ruff.returncode == 0, ruff.stdout[-800:] + ruff.stderr[-300:]))

    regenerate = ROOT / "scripts" / "regenerate_generated.py"
    if regenerate.exists():
        gen = _run([sys.executable, str(regenerate), "--check"])
        results.append(step("generated files", gen.returncode == 0, (gen.stdout + gen.stderr)[-600:]))

    if not args.no_tests and tests:
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
        env.pop("PLANETGEN_REQUIRE_MYSQL", None)
        run = subprocess.run([sys.executable, "-m", "pytest", *tests, "-n", args.workers, "-q", "-p", "no:cacheprovider", "--no-header", "-x", "--maxfail=5"],
                             cwd=ROOT, env=env, capture_output=True, text=True)
        tail = "\n".join(run.stdout.splitlines()[-25:])
        results.append(step("tests (%d files)" % len(tests), run.returncode == 0, tail))
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
