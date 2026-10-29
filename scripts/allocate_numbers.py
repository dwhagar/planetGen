#!/usr/bin/env python3
"""
Gives a branch's generator epochs and schema migrations the next free numbers after main's (OPS.71, OPS.56).

    python scripts/allocate_numbers.py              renumber this branch's epochs and migrations to follow origin/main
    python scripts/allocate_numbers.py --check      only say whether they already do (exit 1 and the reason when not)
    python scripts/allocate_numbers.py --dry-run    say what would change, change nothing

Epoch and migration numbers go to whichever change merges first, so a branch that waited on a long gate finds main has taken its
numbers. This does the renumbering that used to be done by hand, in one pass, on a tree that holds the branch's changes (committed,
merged with main or not, even mid-merge or mid-rebase with `epochs.json` in conflict):

* **Migrations.** The files in `src/planetgen/db/migrations/versions/` that main does not have are the branch's own, in number order.
  They are renamed to follow main's newest (`NNNN_name.py`), their `revision` and `down_revision` rewritten (the first one's
  `down_revision` is main's head), and the `vNNN` and `NNNN` mentions inside them changed to match.
* **Epochs.** The entries of `src/tests/golden/epochs.json` whose reason main does not have are the branch's own. The file becomes
  main's with those entries after it, numbered on from main's last. When the file is in conflict the two sides are read from the index
  (stages 2 and 3), so no conflict has to be resolved by hand first.
* **Old-schema fixture and README line** (`src/tests/fixtures/old_schemas/`). A migration is tested from the schema before it, so the
  first migration needs `schema_v<main's head>.sql.gz`: `schema.sql` as main ships it, without its comment lines and blank lines.
  The fixture the branch made for its old number is replaced by that one; the fixtures of any later own migrations move with them.
* **Mentions in the branch's other files.** In the lines the branch added (against main), `v<old>`, the four-digit `<old>` and
  `epoch <old>` are changed to the new numbers, everywhere but the golden files. Each change is printed, so a line that merely
  talks about a number of main's can be put back.

It fails (exit 2) and changes nothing when a number it would give is already taken. What it cannot do is redo the golden files: they
carry the epoch, so after an epoch moves take main's golden files and run `python -m tests.golden_update` (about a minute) as usual;
it prints that reminder when an epoch moved. Run it after merging or rebasing onto main, then run the targeted tests (migrations,
db schema, golden) and merge.
"""

import argparse
import gzip
import json
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VERSIONS = "src/planetgen/db/migrations/versions"
EPOCHS = "src/tests/golden/epochs.json"
GOLDEN = "src/tests/golden/"
FIXTURES = "src/tests/fixtures/old_schemas"
SCHEMA = "src/planetgen/db/schema.sql"
MIGRATION = re.compile(r"^(\d{4})_(.+\.py)$")
FIXTURE = re.compile(r"^schema_v(\d+)\.sql\.gz$")
CONFLICT = re.compile(r"^(<<<<<<<|=======|>>>>>>>)", re.M)


class Refused(Exception):
    """The numbers cannot be given; the message says why."""


def _git(root, *args, check=True):
    done = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)
    if check and done.returncode != 0:
        raise Refused(f"git {' '.join(args)} failed: {done.stderr.strip()}")
    return done.stdout


def _show(root, ref, path):
    """The text of `path` at `ref`, or None when it is not there."""
    done = subprocess.run(["git", "show", f"{ref}:{path}"], cwd=root, capture_output=True)
    return done.stdout.decode("utf-8") if done.returncode == 0 else None


def _format_epochs(root, data):
    sys.path.insert(0, os.path.join(root, "src"))
    try:
        from planetgen.util import generator_epoch
    finally:
        sys.path.pop(0)
    return generator_epoch.format_epochs(data)


def _migrations(names):
    return {int(m.group(1)): name for name in names if (m := MIGRATION.match(name))}


def _base_migrations(root, base):
    listing = _git(root, "ls-tree", "--name-only", base, f"{VERSIONS}/").split()
    return _migrations(os.path.basename(name) for name in listing)


def _epoch_sides(root):
    """The `epochs.json` texts to read the branch's entries from: the working file, or both sides of a conflict."""
    path = os.path.join(root, EPOCHS)
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    if not CONFLICT.search(text):
        return [json.loads(text)]
    sides = [_show(root, f":{stage}", EPOCHS) for stage in (2, 3)]
    if not all(sides):
        raise Refused(f"{EPOCHS} has conflict markers and the index holds no two sides of it; resolve it, or take main's file")
    return [json.loads(side) for side in sides]


def plan(root, base):
    """
    What to renumber, as a dict: `base_head`, `migrations` (old number to new, in order), `names`, `epochs` (own entries, in order),
    `epoch_map`, `base_epochs` (main's data) and `shift` (whether anything moves).

    Raises:
        Refused: When main cannot be read or the branch's numbers cannot be told apart.
    """
    base_data_text = _show(root, base, EPOCHS)
    if base_data_text is None:
        raise Refused(f"{base} holds no {EPOCHS}; fetch it first (git fetch origin main)")
    base_data = json.loads(base_data_text)
    base_reasons = {entry["reason"] for entry in base_data["epochs"]}
    base_files = _base_migrations(root, base)
    if not base_files:
        raise Refused(f"{base} holds no migrations under {VERSIONS}")
    base_head = max(base_files)

    here = sorted((int(m.group(1)), m.group(0)) for name in os.listdir(os.path.join(root, VERSIONS)) if (m := MIGRATION.match(name)))
    own_pairs = [(number, name) for number, name in here if name not in base_files.values()]
    own = [number for number, _ in own_pairs]
    if len(set(own)) != len(own):
        raise Refused(f"the branch has two migrations numbered alike: {', '.join(name for _, name in own_pairs)}")
    migrations = {old: base_head + 1 + index for index, old in enumerate(own)}

    seen = {}
    for data in _epoch_sides(root):
        for entry in data["epochs"]:
            if entry["reason"] not in base_reasons and entry["reason"] not in seen:
                seen[entry["reason"]] = entry
    epochs = sorted(seen.values(), key=lambda entry: entry["epoch"])
    first = len(base_data["epochs"]) + 1
    epoch_map = {entry["epoch"]: first + index for index, entry in enumerate(epochs)}
    shift = any(old != new for old, new in migrations.items()) or any(old != new for old, new in epoch_map.items())
    return {"base_head": base_head, "migrations": migrations, "names": dict(own_pairs), "epochs": epochs,
            "epoch_map": epoch_map, "base_data": base_data, "shift": shift, "base_files": base_files}


def check(root, base):
    """Whether the numbers already follow main's: `(ok, message)`."""
    found = plan(root, base)
    if not found["shift"]:
        return True, (f"ok: {len(found['migrations'])} migration(s) from {found['base_head'] + 1} and "
                      f"{len(found['epochs'])} epoch(s) from {len(found['base_data']['epochs']) + 1} already follow {base}")
    moves = [f"migration {old:04d} -> {new:04d}" for old, new in found["migrations"].items() if old != new]
    moves += [f"epoch {old} -> {new}" for old, new in found["epoch_map"].items() if old != new]
    return False, f"not following {base}: " + ", ".join(moves)


def _replace_numbers(text, migrations, epochs):
    """`text` with `v<old>`, `<old:04d>` and `epoch <old>` changed to the new numbers, in one pass (so 174 -> 175 and 175 -> 176 do not chain)."""
    patterns = []
    if migrations:
        patterns.append(r"(?P<v>\bv(?P<vn>%s)\b)" % "|".join(str(n) for n in migrations))
        patterns.append(r"(?P<z>(?<![0-9A-Za-z.])(?P<zn>%s)(?![0-9A-Za-z.]))" % "|".join(f"{n:04d}" for n in migrations))
    if epochs:
        patterns.append(r"(?P<e>\b(?P<ep>[Ee]poch) (?P<en>%s)\b)" % "|".join(str(n) for n in epochs))
    if not patterns:
        return text

    def swap(match):
        found = match.groupdict()
        if found.get("v"):
            return f"v{migrations[int(found['vn'])]}"
        if found.get("z"):
            return f"{migrations[int(found['zn'])]:04d}"
        return f"{found['ep']} {epochs[int(found['en'])]}"
    return re.sub("|".join(patterns), swap, text)


def _added_lines(root, base, path):
    """The line numbers (1-based) of `path` that differ from `base`: the branch's own."""
    out = _git(root, "diff", "-U0", "--no-color", base, "--", path, check=False)
    lines = set()
    for match in re.finditer(r"^@@ -\S+ \+(\d+)(?:,(\d+))? @@", out, re.M):
        start, count = int(match.group(1)), int(match.group(2) if match.group(2) is not None else 1)
        lines.update(range(start, start + count))
    return lines


def _own_files(root, base):
    """Tracked or new text files that differ from `base`, other than the ones handled by name."""
    names = _git(root, "diff", "--name-only", base).split()
    names += _git(root, "ls-files", "--others", "--exclude-standard").split()
    skip = (VERSIONS + "/", GOLDEN, FIXTURES + "/", "changes/")
    return sorted({n for n in names if not n.startswith(skip) and os.path.isfile(os.path.join(root, n))})


def apply(root, base, dry_run=False, say=print):
    """
    Renumbers the branch's migrations, epochs, fixture and mentions to follow `base`.

    Returns:
        bool: Whether anything moved.

    Raises:
        Refused: A number it would give is already taken, or main cannot be read.
    """
    found = plan(root, base)
    migrations, epoch_map = found["migrations"], found["epoch_map"]
    if not found["shift"]:
        say(check(root, base)[1])
        return False
    for new in migrations.values():
        if new in found["base_files"] and found["base_files"][new] not in found["names"].values():
            if new <= found["base_head"]:
                raise Refused(f"migration {new:04d} is already main's ({found['base_files'][new]})")
    moved_migrations = {old: new for old, new in migrations.items() if old != new}
    moved_epochs = {old: new for old, new in epoch_map.items() if old != new}
    todo = []

    # migrations: rename, then rewrite revision ids and mentions
    folder = os.path.join(root, VERSIONS)
    previous = found["base_head"]
    for old, new in migrations.items():
        old_name = found["names"][old]
        new_name = f"{new:04d}_{MIGRATION.match(old_name).group(2)}"
        if new != old and os.path.exists(os.path.join(folder, new_name)) and old_name != new_name:
            raise Refused(f"{VERSIONS}/{new_name} already exists")
        with open(os.path.join(folder, old_name), encoding="utf-8") as handle:
            text = handle.read()
        text = _replace_numbers(text, moved_migrations, moved_epochs)
        text = re.sub(r'^(revision\s*=\s*)["\']\d+["\']', lambda m: f'{m.group(1)}"{new:04d}"', text, flags=re.M)
        text = re.sub(r'^(down_revision\s*=\s*)["\']\d+["\']', lambda m: f'{m.group(1)}"{previous:04d}"', text, flags=re.M)
        previous = new
        todo.append(("migration", old_name, new_name, text))

    # fixtures: the first own migration is tested from main's head schema; later ones keep theirs, moved
    fixtures_dir = os.path.join(root, FIXTURES)
    base_fixtures = set(_git(root, "ls-tree", "--name-only", base, f"{FIXTURES}/").split())
    base_fixture_names = {os.path.basename(name) for name in base_fixtures}
    own_fixtures = {int(FIXTURE.match(n).group(1)): n for n in os.listdir(fixtures_dir)
                    if FIXTURE.match(n) and n not in base_fixture_names}
    fixture_ops = []
    if migrations:
        first_old = min(migrations)
        needed = found["base_head"]
        if f"schema_v{needed}.sql.gz" in base_fixture_names:
            raise Refused(f"main already holds schema_v{needed}.sql.gz, so its head is not {needed}")
        schema = _show(root, base, SCHEMA)
        if schema is None:
            raise Refused(f"{base} holds no {SCHEMA}")
        for number, name in sorted(own_fixtures.items()):
            if number < first_old:                       # the old number's own fixture, made from the schema main had then
                fixture_ops.append(("drop", name, None, None))
            elif number + 1 in migrations:               # between two own migrations: moves with the one it precedes
                target = f"schema_v{migrations[number + 1] - 1}.sql.gz"
                if target != name:
                    fixture_ops.append(("move", name, target, None))
        body = "".join(line + "\n" for line in schema.splitlines() if line.strip() and not line.lstrip().startswith("--"))
        fixture_ops.append(("make", None, f"schema_v{needed}.sql.gz", body))

    # the other files: only the lines the branch added
    edits = []
    for path in _own_files(root, base):
        full = os.path.join(root, path)
        try:
            with open(full, encoding="utf-8") as handle:
                lines = handle.read().split("\n")
        except UnicodeDecodeError:
            continue
        mine = _added_lines(root, base, path) if _show(root, base, path) is not None else set(range(1, len(lines) + 1))
        changed = False
        for index in sorted(mine):
            if index - 1 >= len(lines):
                continue
            after = _replace_numbers(lines[index - 1], moved_migrations, moved_epochs)
            if after != lines[index - 1]:
                edits.append(f"{path}:{index}: {lines[index - 1].strip()[:90]}  ->  {after.strip()[:90]}")
                lines[index - 1] = after
                changed = True
        if changed:
            todo.append(("file", path, path, "\n".join(lines)))

    # epochs.json: main's entries, then the branch's, numbered on
    data = found["base_data"]
    entries = [dict(entry) for entry in found["epochs"]]
    for entry in entries:
        entry["epoch"] = epoch_map[entry["epoch"]]
    renewed = {**data, "epochs": data["epochs"] + entries}
    epochs_text = _format_epochs(root, renewed)

    for kind, old, new, _text in todo:
        if kind == "migration":
            say(f"migration {old} -> {new}")
    for old, new in moved_epochs.items():
        say(f"epoch {old} -> {new}")
    for kind, old, new, _body in fixture_ops:
        say(f"fixture {kind} {old or ''} {new or ''}".rstrip())
    for line in edits:
        say("edit " + line)
    if dry_run:
        return True

    for kind, old, new, text in todo:
        if kind == "migration":
            if old != new:
                os.remove(os.path.join(folder, old))
            with open(os.path.join(folder, new), "w", encoding="utf-8") as handle:
                handle.write(text)
        else:
            with open(os.path.join(root, new), "w", encoding="utf-8") as handle:
                handle.write(text)
    makes = [(new, body) for kind, _old, new, body in fixture_ops if kind == "make"]
    staging = []
    for kind, old, new, _body in fixture_ops:
        if kind == "drop":
            os.remove(os.path.join(fixtures_dir, old))
        elif kind == "move":                             # through temporary names, so a shift by one does not overwrite itself
            temp = os.path.join(fixtures_dir, old + ".moving")
            os.rename(os.path.join(fixtures_dir, old), temp)
            staging.append((temp, os.path.join(fixtures_dir, new)))
    for temp, final in staging:
        os.rename(temp, final)
    for new, body in makes:
        with open(os.path.join(fixtures_dir, new), "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as handle:
            handle.write(body.encode("utf-8"))
    _readme(root, base, migrations, found["base_head"], say)
    with open(os.path.join(root, EPOCHS), "w", encoding="utf-8") as handle:
        handle.write(epochs_text)
    if moved_epochs:
        say("the golden files carry the epoch: take main's golden files and run `python -m tests.golden_update` (about a minute) "
            "so their digests follow the new number")
    return True


def _readme(root, base, migrations, base_head, say):
    """Moves the branch's own lines of the fixtures README to the new versions, the first one's to main's head."""
    path = os.path.join(root, FIXTURES, "README.md")
    base_text = _show(root, base, f"{FIXTURES}/README.md") or ""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    own = [i for i, line in enumerate(lines) if re.match(r"- v\d+:", line) and line not in base_text.split("\n")]
    if not own or not migrations:
        return
    head = _git(root, "rev-parse", "--short", base).strip()
    first_old = min(migrations)
    for index in own:
        number = int(re.match(r"- v(\d+):", lines[index]).group(1))
        if number == first_old - 1:
            rest = re.sub(r"^- v\d+: \S+", "", lines[index])
            lines[index] = f"- v{base_head}: {head}{rest}"
        elif number + 1 in migrations:
            lines[index] = re.sub(r"^- v\d+", f"- v{migrations[number + 1] - 1}", lines[index])
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    say("fixtures README: own lines moved")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Give this branch's epochs and migrations the next free numbers after main's.")
    parser.add_argument("--base", default="origin/main", help="main's ref (default origin/main)")
    parser.add_argument("--check", action="store_true", help="only report whether the numbers already follow main's")
    parser.add_argument("--dry-run", action="store_true", help="say what would change and change nothing")
    parser.add_argument("--root", default=ROOT, help="the repository root (default: this one)")
    args = parser.parse_args(argv)
    try:
        if args.check:
            ok, message = check(args.root, args.base)
            print(message)
            return 0 if ok else 1
        apply(args.root, args.base, dry_run=args.dry_run)
    except Refused as error:
        print(f"refused: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
