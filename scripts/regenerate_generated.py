#!/usr/bin/env python3
"""
Regenerates, checks, or repairs after a merge the generated files a lane's pull request carries (TEST.721, Lane 03's gate-time report).

    python scripts/regenerate_generated.py --list               the generated files and the builder of each
    python scripts/regenerate_generated.py                      rebuild every one (what to run after editing a source)
    python scripts/regenerate_generated.py --check              fail and name the ones that are out of date (the pre-gate quick lane runs this)
    python scripts/regenerate_generated.py --only stand-in-values reference-pages     only those builders
    python scripts/regenerate_generated.py --resolve            in a merge with conflicts: rebuild the conflicted generated files and `git add` them

Why: a generated file (the stand-in list, the references pages, the data tables built from a script) is a copy of its sources, so two pull requests that each change a source
change the copy too, and the second to merge conflicts on a file nobody can merge by hand. Nine such conflicts in one week each cost one to two hours of round trip. A conflict in a
generated file has one right resolution, rebuild it from the merged sources; `--resolve` does exactly that and nothing else:

1. Resolve and `git add` every conflict in a source file first (for example `docs/design/research/gap-filling-plan.md`).
2. Run `python scripts/regenerate_generated.py --resolve`. It refuses (exit 2) while any conflict outside the generated set is left, because a builder reads the merged sources.
   For each conflicted generated path it takes our side (so no conflict marker is read), rebuilds the builder that owns the path, checks it, and `git add`s what the builder wrote.
3. `git commit` to finish the merge.

The global documents (`docs/generated`, `docs/site`) are not in this set: only the TODO thread rebuilds them (`scripts/build_global_docs.py`, DOC.69). Each builder here is the
existing script of the file; this one only knows which file belongs to which script, so a new generated file is one line in `GENERATED`.
"""

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@dataclass(frozen=True)
class Generated:
    """One generated file (or folder, `path` ending in a slash) and the script that builds it; every script takes `--check` and writes with no argument."""

    name: str
    script: str
    paths: tuple


GENERATED = (
    Generated("stand-in-values", "build_stand_in_values.py", ("docs/design/research/stand-in-values.md",)),
    Generated("reference-pages", "build_reference_pages.py", ("docs/references/",)),
    Generated("api-compatibility", "build_api_compatibility.py", ("docs/reference/api-compatibility.md",)),
    Generated("log-vocabulary", "build_log_vocabulary.py", ("scripts/logrecord.sh",)),
    Generated("life-class-rules", "build_life_class_rules.py", ("docs/design/habitability/life-class-rules.md",)),
    Generated("gas-table", "build_gas_table.py", ("src/planetgen/physics/data/gases.json",)),
    Generated("materials-table", "build_materials_table.py", ("src/planetgen/physics/data/materials.json",)),
    Generated("atmosphere-recipes", "build_atmosphere_recipes.py", ("src/planetgen/physics/data/atmosphere_recipes.json",)),
    Generated("fudge-table", "build_fudge_table.py", ("src/planetgen/physics/data/fudge_table.json",)),
    Generated("unit-registry", "build_unit_registry.py", ("src/html/static/unitregistry.js",)),
    Generated("unit-ladders", "build_unit_ladders.py", ("src/html/static/unitladders.js",)),
)
"""tuple: The generated files a lane's pull request may carry. Not here: the global documents (TODO only, DOC.69), the data built from downloads, and `docs/history` (a snapshot)."""


def owner_of(path):
    """The `Generated` entry that owns `path` (a repository-relative path with forward slashes), or None."""
    for entry in GENERATED:
        for owned in entry.paths:
            if path == owned or (owned.endswith("/") and path.startswith(owned)):
                return entry
    return None


def _run(entry, *flags):
    environment = dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src") + os.pathsep + os.environ.get("PYTHONPATH", ""))
    return subprocess.run([sys.executable, os.path.join(ROOT, "scripts", entry.script), *flags], cwd=ROOT, env=environment, capture_output=True, text=True)


def _git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)


def unmerged_paths():
    """The paths with an unresolved merge conflict, repository-relative."""
    result = _git("diff", "--name-only", "--diff-filter=U")
    return [line for line in result.stdout.splitlines() if line]


def command_regenerate(entries):
    failed = []
    for entry in entries:
        result = _run(entry)
        print(f"{'ok  ' if result.returncode == 0 else 'FAIL'} {entry.name}: {entry.script}")
        if result.returncode != 0:
            failed.append(entry.name)
            print((result.stdout + result.stderr).strip()[-600:])
    return 1 if failed else 0


def command_check(entries):
    stale = []
    for entry in entries:
        result = _run(entry, "--check")
        if result.returncode != 0:
            stale.append(entry)
            print(f"out of date: {entry.name} ({', '.join(entry.paths)}): run python scripts/regenerate_generated.py --only {entry.name}")
    print("generated files are up to date" if not stale else f"{len(stale)} of {len(entries)} generated file sets are out of date")
    return 1 if stale else 0


def command_resolve(entries):
    conflicts = unmerged_paths()
    if not conflicts:
        print("no merge conflict")
        return 0
    wanted = {entry.name: entry for entry in entries}
    generated = [(path, owner_of(path)) for path in conflicts]
    others = [path for path, owner in generated if owner is None or owner.name not in wanted]
    if others:
        print("resolve these conflicts first (they are not generated files), `git add` them, and run --resolve again:")
        for path in others:
            print(f"  {path}")
        return 2
    builders = sorted({owner.name: owner for _, owner in generated}.values(), key=lambda entry: entry.name)
    for path, _ in generated:
        _git("checkout", "--ours", "--", path)          # no conflict marker is left for the builder to read
    status = 0
    for entry in builders:
        result = _run(entry)
        if result.returncode != 0:
            print(f"FAIL {entry.name}: {(result.stdout + result.stderr).strip()[-600:]}")
            status = 1
            continue
        if _run(entry, "--check").returncode != 0:
            print(f"FAIL {entry.name}: the rebuilt file still fails --check")
            status = 1
            continue
        paths = [path for path in entry.paths if os.path.exists(os.path.join(ROOT, path))]
        _git("add", "--", *paths)
        print(f"rebuilt and added: {entry.name} ({', '.join(paths)})")
    if not status:
        print("done: `git commit` finishes the merge")
    return status


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--list", action="store_true", help="list the generated files and their builders")
    parser.add_argument("--check", action="store_true", help="fail naming the generated files that are out of date")
    parser.add_argument("--resolve", action="store_true", help="in a merge with conflicts, rebuild and `git add` the conflicted generated files")
    parser.add_argument("--only", nargs="+", metavar="NAME", help="only these builders (names from --list)")
    args = parser.parse_args(argv)
    entries = list(GENERATED)
    if args.only:
        unknown = sorted(set(args.only) - {entry.name for entry in GENERATED})
        if unknown:
            parser.error(f"unknown builder(s): {', '.join(unknown)}; see --list")
        entries = [entry for entry in GENERATED if entry.name in args.only]
    if args.list:
        for entry in entries:
            print(f"{entry.name:20} {entry.script:34} {', '.join(entry.paths)}")
        return 0
    if args.resolve:
        return command_resolve(entries)
    if args.check:
        return command_check(entries)
    return command_regenerate(entries)


if __name__ == "__main__":
    sys.exit(main())
