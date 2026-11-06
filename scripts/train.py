#!/usr/bin/env python3
"""
The merge-train runner's steps as one command, the script rule 10 of the merge-train runbook asked for.

    python scripts/train.py status                      the next free numbers on origin/main: migration, epoch, section, citation
    python scripts/train.py merge heads.txt [--label "Train 13"]
                                                        merge each head (one `branch` per line, stacking order) into this branch
    python scripts/train.py collisions heads.txt        which migration, epoch, stand-in section and citation numbers two heads (or a head
                                                        and main) both claim, before any merge
    python scripts/train.py numbers                     renumber this branch's migrations and epochs after main's (allocate_numbers)
    python scripts/train.py rebuild                     rebuild the generated files (references, pages, stand-in values, spelling)
    python scripts/train.py golden --reason "..." [--affects plan,scatter,fill,population]
                                                        regenerate the golden files and the epoch's reason (about ten minutes)
    python scripts/train.py checks                      the guards a train must pass before its PR: wiring, numbers, release note
    python scripts/train.py gate [--of N]               print the gate's file list, split over N machines when asked

A head that conflicts in a file the train rebuilds (references, stand-in values) takes this branch's side and the file is rebuilt by
`rebuild`; a conflict anywhere else stops the merge with the files named, and the runner resolves it by hand (`git status`) and runs
the command again with the heads that are left. Heads already merged are skipped.
"""

import argparse
import json
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GENERATED = re.compile(r"^(docs/design/research/stand-in-values\.md|docs/references/.*\.json|docs/REFERENCES\.(md|html)|"
                       r"docs/reference/references-data\.json)$")
"""re.Pattern: The files `rebuild` writes; a merge conflict in one takes this branch's side."""


def sh(*args, check=True, capture=False, cwd=ROOT):
    done = subprocess.run(list(args), cwd=cwd, text=True, capture_output=capture)
    if check and done.returncode != 0:
        if capture:
            sys.stderr.write(done.stdout + done.stderr)
        raise SystemExit(f"train: `{' '.join(args)}` failed ({done.returncode})")
    return done.stdout if capture else done.returncode


def main_numbers(ref="origin/main"):
    """The next free numbers on `ref`: `{migration, epoch, section, citation}`."""
    tree = sh("git", "ls-tree", "--name-only", ref, "src/planetgen/db/migrations/versions/", capture=True).split()
    migrations = [int(m.group(1)) for n in tree if (m := re.search(r"/(\d{4})_[^/]+\.py$", n))]
    epochs = json.loads(sh("git", "show", f"{ref}:src/tests/golden/epochs.json", capture=True))["epochs"]
    plan = sh("git", "show", f"{ref}:docs/design/research/gap-filling-plan.md", capture=True)
    sections = [int(m.group(1)) for m in re.finditer(r"^#{2,4} 6\.(\d+)\b", plan, re.M)]
    refs = json.loads(sh("git", "show", f"{ref}:docs/reference/references-data.json", capture=True))
    last = epochs[-1]["epoch"] if isinstance(epochs[-1], dict) and "epoch" in epochs[-1] else len(epochs)
    return {"migration": f"{max(migrations) + 1:04d}", "epoch": last + 1, "section": f"6.{max(sections) + 1}", "citation": f"R{refs['id_next']}"}


def status(args):
    sh("git", "fetch", "-q", "origin", "main")
    for key, value in main_numbers().items():
        print(f"next {key}: {value}")


def merge(args):
    sh("git", "fetch", "-q", "origin")
    heads = [h.strip() for h in open(args.heads, encoding="utf-8") if h.strip() and not h.startswith("#")]
    for head in heads:
        ref = f"origin/{head}"
        if sh("git", "merge-base", "--is-ancestor", ref, "HEAD", check=False) == 0:
            print(f"== {head}: already merged")
            continue
        print(f"== {head}")
        sh("git", "merge", "--no-ff", "-m", f"{args.label}: merge {head}", ref, check=False, capture=True)
        unmerged = sh("git", "diff", "--name-only", "--diff-filter=U", capture=True).split()
        if not unmerged:
            continue
        other = [f for f in unmerged if not GENERATED.match(f)]
        if other:
            print("STOP: conflicts to resolve by hand, then run `git commit` and this command again:\n  " + "\n  ".join(other))
            raise SystemExit(1)
        for f in unmerged:
            sh("git", "checkout", "--ours", "--", f)
            sh("git", "add", f)
        sh("git", "commit", "-qm", f"{args.label}: merge {head} (generated files rebuilt later)")


def _claims(head):
    """`{kind: set}` of the numbers `origin/<head>` adds beyond its merge base with main: migrations, sections, citations, epochs."""
    ref = f"origin/{head}"
    base = sh("git", "merge-base", "origin/main", ref, capture=True).strip()

    def show(rev, path):
        done = subprocess.run(["git", "show", f"{rev}:{path}"], cwd=ROOT, text=True, capture_output=True)
        return done.stdout if done.returncode == 0 else ""

    names = sh("git", "diff", "--name-only", "--diff-filter=A", base, ref, "--", "src/planetgen/db/migrations/versions/", capture=True).split()
    plan = "docs/design/research/gap-filling-plan.md"
    section = re.compile(r"^#{2,4} (6\.\d+)\b", re.M)
    cites = "docs/reference/references-data.json"

    def ids(rev):
        try:
            return set(json.loads(show(rev, cites)).get("ids", {}))
        except ValueError:
            return set()

    def epochs(rev):
        try:
            return {e["epoch"] for e in json.loads(show(rev, "src/tests/golden/epochs.json"))["epochs"]}
        except (ValueError, KeyError):
            return set()

    return {"migration": {m.group(1) for n in names if (m := re.search(r"/(\d{4})_", n))},
            "section": set(section.findall(show(ref, plan))) - set(section.findall(show(base, plan))),
            "citation": ids(ref) - ids(base), "epoch": epochs(ref) - epochs(base)}


def _claims_of(ref):
    """Every migration, section, citation and epoch number that `ref` holds."""
    tree = sh("git", "ls-tree", "--name-only", ref, "src/planetgen/db/migrations/versions/", capture=True).split()
    plan = sh("git", "show", f"{ref}:docs/design/research/gap-filling-plan.md", capture=True)
    refs = json.loads(sh("git", "show", f"{ref}:docs/reference/references-data.json", capture=True))
    epochs = json.loads(sh("git", "show", f"{ref}:src/tests/golden/epochs.json", capture=True))["epochs"]
    return {"migration": {m.group(1) for n in tree if (m := re.search(r"/(\d{4})_", n))},
            "section": set(re.findall(r"^#{2,4} (6\.\d+)\b", plan, re.M)), "citation": set(refs.get("ids", {})),
            "epoch": {e["epoch"] for e in epochs}}


def collisions(args):
    heads = [h.strip() for h in open(args.heads, encoding="utf-8") if h.strip() and not h.startswith("#")]
    sh("git", "fetch", "-q", "origin")
    claims = {head: _claims(head) for head in heads}
    taken = {"migration": set(), "section": set(), "citation": set(), "epoch": set()}
    clash = 0
    for kind in taken:
        owners = {}
        for head, found in claims.items():
            for number in found[kind]:
                owners.setdefault(number, []).append(head)
        for number, who in sorted(owners.items()):
            if len(who) > 1:
                clash += 1
                print(f"CLASH {kind} {number}: {', '.join(who)}")
    on_main = _claims_of("origin/main")
    for head, found in claims.items():
        for kind, numbers_ in found.items():
            for number in sorted(numbers_ & on_main[kind]):
                clash += 1
                print(f"TAKEN ON MAIN {kind} {number}: {head} (renumber it)")
    print(f"{clash} clash(es)" if clash else "no number is claimed twice or already taken on main")
    print("numbers on main now:", ", ".join(f"{k} {v}" for k, v in main_numbers().items()))


def numbers(args):
    sh(sys.executable, "scripts/allocate_numbers.py")


def rebuild(args):
    for command in (("build_references.py", "render"), ("build_reference_pages.py",), ("build_stand_in_values.py",),
                    ("american_spelling_sweep.py", "apply")):
        script = command[0]
        if not os.path.exists(os.path.join(ROOT, "scripts", script)):
            continue
        sh(sys.executable, f"scripts/{script}", *command[1:])
    print("rebuilt; commit what changed")


def golden(args):
    sh(sys.executable, "-m", "tests.golden_update", "--reason", args.reason, "--affects", args.affects, cwd=os.path.join(ROOT, "src"))


def checks(args):
    sh("git", "fetch", "-q", "origin", "main")
    out = sh(sys.executable, "scripts/wiring_guard.py", "--check", capture=True)
    problems = [line for line in out.splitlines() if line.startswith("PROBLEM")]
    print(f"wiring guard: {len(problems)} problem(s)")
    for line in problems:
        print("  " + line)
    sh(sys.executable, "scripts/allocate_numbers.py", "--base", "origin/main", "--check")
    sh(sys.executable, "scripts/bump_version.py", "--check-pr", "origin/main")
    if problems:
        raise SystemExit(1)


def gate(args):
    out = sh(sys.executable, "scripts/train_gate.py", "files", capture=True)
    files = out.split()
    if args.of and args.of > 1:
        open(os.path.join(ROOT, "gate-files.tmp"), "w", encoding="utf-8").write("\n".join(files) + "\n")
        print(sh(sys.executable, "scripts/train_gate.py", "shard", "--files", "gate-files.tmp", "--of", str(args.of), capture=True))
        os.remove(os.path.join(ROOT, "gate-files.tmp"))
    else:
        print("\n".join(files))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status").set_defaults(run=status)
    p = sub.add_parser("merge")
    p.add_argument("heads")
    p.add_argument("--label", default="Train")
    p.set_defaults(run=merge)
    p = sub.add_parser("collisions")
    p.add_argument("heads")
    p.set_defaults(run=collisions)
    sub.add_parser("numbers").set_defaults(run=numbers)
    sub.add_parser("rebuild").set_defaults(run=rebuild)
    p = sub.add_parser("golden")
    p.add_argument("--reason", required=True)
    p.add_argument("--affects", default="plan,scatter,fill,population")
    p.set_defaults(run=golden)
    sub.add_parser("checks").set_defaults(run=checks)
    p = sub.add_parser("gate")
    p.add_argument("--of", type=int, default=0)
    p.set_defaults(run=gate)
    args = parser.parse_args(argv)
    args.run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
