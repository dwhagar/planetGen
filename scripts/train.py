#!/usr/bin/env python3
"""
The merge-train runner's steps as one command, the script rule 10 of the merge-train runbook asked for.

    python scripts/train.py status                      the next free numbers on origin/main: migration, epoch, section, citation
    python scripts/train.py merge heads.txt [--label "Train 13"]
                                                        merge each head (one `branch` per line, stacking order) into this branch
    python scripts/train.py collisions heads.txt        which migration, epoch, stand-in section and citation numbers two heads (or a head
                                                        and main) both claim, before any merge
    python scripts/train.py renumber HEAD OLD NEW       in the lines a head added, change a stand-in section (6.193 to 6.196) or a citation
                                                        (R5136 to R5140) to a free number, after its merge
    python scripts/train.py numbers                     renumber this branch's migrations and epochs after main's (allocate_numbers)
    python scripts/train.py rebuild                     rebuild the generated files (references, pages, stand-in values, spelling)
    python scripts/train.py affects                     the stages the train's release notes say they alter (their `generation-affects:` lines)
    python scripts/train.py golden --reason "..." [--affects plan,scatter,fill,population|all]
                                                        regenerate the golden files and the epoch's reason (about ten minutes); with no
                                                        --affects, the stages are the union of the notes' `generation-affects:` lines
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bump_version  # noqa: E402 -- the note rules (`generation-affects:`) live there
GENERATED = re.compile(r"^(docs/design/research/stand-in-values\.md|docs/references/.*\.json|docs/REFERENCES\.(md|html)|"
                       r"docs/reference/references-data\.json|src/tests/golden/[^/]+)$")
"""re.Pattern: The files `rebuild` and `golden` write; a merge conflict in one takes this branch's side."""


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


def renumber(args):
    """Changes `OLD` to `NEW` in the lines that `origin/<head>` added (against its merge base), wherever those lines are now."""
    ref = f"origin/{args.head}"
    base = sh("git", "merge-base", "origin/main", ref, capture=True).strip()
    section = re.fullmatch(r"6\.\d+", args.old) and re.fullmatch(r"6\.\d+", args.new)
    citation = re.fullmatch(r"R\d+", args.old) and re.fullmatch(r"R\d+", args.new)
    if not (section or citation):
        raise SystemExit("train: renumber takes two stand-in sections (6.193 6.196) or two citations (R5136 R5140)")
    token = re.compile(rf"(?<![\w.]){re.escape(args.old)}(?![\d])")
    names = sh("git", "diff", "--name-only", "--diff-filter=AM", base, ref, capture=True).split()
    changed = 0
    for name in names:
        if GENERATED.match(name) or name.endswith((".gz", ".svg", ".png", ".pyc")):
            continue
        path = os.path.join(ROOT, name)
        if not os.path.exists(path):
            continue
        diff = sh("git", "diff", "-U0", base, ref, "--", name, capture=True)
        added = {line[1:] for line in diff.splitlines() if line.startswith("+") and not line.startswith("+++") and token.search(line)}
        if not added:
            continue
        lines = open(path, encoding="utf-8").read().split("\n")
        for index, line in enumerate(lines):
            if line in added:
                lines[index] = token.sub(args.new, line)
                changed += 1
        open(path, "w", encoding="utf-8").write("\n".join(lines))
    print(f"{args.old} -> {args.new}: {changed} line(s) changed in the files {args.head} touched; check `git diff`, then `rebuild`")


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


def train_notes(ref="origin/main"):
    """`{file name: text}` of the release notes this branch adds to (or changes in) `changes/` beyond `ref`: the train's items."""
    names = sh("git", "diff", "--name-only", "--diff-filter=AM", ref, "HEAD", "--", "changes/", capture=True).split()
    notes = {}
    for name in names:
        base = os.path.basename(name)
        path = os.path.join(ROOT, name)
        if base.endswith(".md") and base != "README.md" and os.path.exists(path):
            notes[base] = open(path, encoding="utf-8").read()
    return notes


def stages_of_notes(notes):
    """
    The stages a train's epoch alters, from its notes' `generation-affects:` lines (the rule in `changes/README.md`).

    A note with no such line says nothing, so it counts as every stage (what the train did for all of them before): it is named in
    `undeclared` so the lane or TODO can add the line. `generation-affects: none` alters nothing and adds nothing.

    Returns:
        tuple: `(stages, undeclared, problems)`: the stages in epoch order, the notes without a line, and the notes whose line is wrong.
    """
    stages, undeclared, problems = set(), [], []
    for name, text in sorted(notes.items()):
        try:
            declared = bump_version.parse_generation_affects(text)
        except ValueError as exc:
            problems.append(f"{name}: {exc}")
            continue
        if declared is None:
            undeclared.append(name)
        else:
            stages.update(declared)
    if undeclared:
        stages.update(bump_version.GENERATION_STAGES)
    return tuple(stage for stage in bump_version.GENERATION_STAGES if stage in stages), undeclared, problems


def affects(args):
    stages, undeclared, problems = stages_of_notes(train_notes())
    for problem in problems:
        print("PROBLEM " + problem)
    for name in undeclared:
        print(f"UNDECLARED {name}: no `generation-affects:` line, so it counts as every stage")
    print("stages: " + (",".join(stages) if stages else "none declared"))
    if problems:
        raise SystemExit(1)


def golden(args):
    if args.affects in (None, ""):
        stages, undeclared, problems = stages_of_notes(train_notes())
        if problems:
            raise SystemExit("train: fix these notes first:\n  " + "\n  ".join(problems))
        if not stages:
            raise SystemExit("train: every note says `generation-affects: none` (or there is no note), yet golden is being rewritten; "
                             "give --affects by hand")
        for name in undeclared:
            print(f"train: {name} has no `generation-affects:` line, so every stage counts")
        chosen = ",".join(stages)
        print(f"train: the epoch alters {chosen}")
    else:
        chosen = ",".join(bump_version.GENERATION_STAGES) if args.affects == "all" else args.affects
    sh(sys.executable, "-m", "tests.golden_update", "--reason", args.reason, "--affects", chosen, cwd=os.path.join(ROOT, "src"))


def checks(args):
    sh("git", "fetch", "-q", "origin", "main")
    out = sh(sys.executable, "scripts/wiring_guard.py", "--check", capture=True)
    problems = [line for line in out.splitlines() if line.startswith("PROBLEM")]
    print(f"wiring guard: {len(problems)} problem(s)")
    stages, undeclared, bad = stages_of_notes(train_notes())
    print(f"notes: stages {','.join(stages) or 'none'}; {len(undeclared)} note(s) without a `generation-affects:` line, {len(bad)} wrong")
    problems += [f"PROBLEM note {line}" for line in bad]
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
    p = sub.add_parser("renumber")
    p.add_argument("head")
    p.add_argument("old")
    p.add_argument("new")
    p.set_defaults(run=renumber)
    sub.add_parser("numbers").set_defaults(run=numbers)
    sub.add_parser("rebuild").set_defaults(run=rebuild)
    p = sub.add_parser("golden")
    p.add_argument("--reason", required=True)
    p.add_argument("--affects", default=None, help="stages, comma separated, or `all`; left out, the union of the notes' `generation-affects:` lines")
    p.set_defaults(run=golden)
    sub.add_parser("affects").set_defaults(run=affects)
    sub.add_parser("checks").set_defaults(run=checks)
    p = sub.add_parser("gate")
    p.add_argument("--of", type=int, default=0)
    p.set_defaults(run=gate)
    args = parser.parse_args(argv)
    args.run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
