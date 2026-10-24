# tests/golden_update.py

"""
Bumps the generator epoch after the golden files were changed on purpose (OPS.56):

    python -m tests.golden_update --reason "what changed in the generated output and why" --affects fill

A PR that changes what the generator produces for a seed changes the golden files under `src/tests/golden/`, and
`scripts/bump_version.py --check-pr` then requires three things together: those files, the epoch bumped by one, and
a line `generation-output: changed` in the PR's `changes/` note. `--affects` names the stages the change altered
(`plan`, `scatter`, `fill`, `population`, comma separated, the deepest ones: GEN.764 redoes exactly those in a galaxy made
under the epoch before; `plan` means such a galaxy needs a full regeneration); it is required for a new epoch. This tool does the first two: it rewrites the golden
files (TEST.77: `tests/test_golden_seed.py` run with `PLANETGEN_WRITE_GOLDEN=1`, about ten minutes, it needs the test
database) and appends the next epoch to `src/tests/golden/epochs.json`, with the reason and no release yet (stamping
the release names it). After a merge conflict in the golden files, take main's files and run it again: it writes the
files afresh and leaves the epoch the change already holds.

It refuses under CI (a golden file is never to be regenerated there, or a regression would rewrite its own proof),
when the golden files exist and none differs from main after the rewrite (nothing to bump for) and for a second bump
in the same change (the epoch already went up; the new golden changes ride on it). `--no-write` skips the rewrite for
golden files that were already written.
"""

import argparse
import json
import os
import subprocess
import sys

from planetgen.util import generator_epoch

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GOLDEN_DIR = "src/tests/golden"
EPOCHS_FILE = f"{GOLDEN_DIR}/epochs.json"


class Refused(Exception):
    """The bump is not allowed; the message says why."""


def _git(root, *args):
    return subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, check=True).stdout


def _changed_golden_files(root):
    """The golden files other than `epochs.json` that differ from where this branch left `origin/main` (committed
    or not), or are new."""
    try:
        base = _git(root, "merge-base", "origin/main", "HEAD").strip()
    except (OSError, subprocess.CalledProcessError):
        base = "HEAD"
    paths = _git(root, "diff", "--name-only", base, "--", GOLDEN_DIR).split()
    paths += _git(root, "ls-files", "--others", "--exclude-standard", "--", GOLDEN_DIR).split()
    return [path for path in paths if os.path.basename(path) != "epochs.json"]


def _golden_files_exist(root):
    """Whether any golden file other than `epochs.json` is checked in."""
    names = _git(root, "ls-files", "--", GOLDEN_DIR).split()
    return any(os.path.basename(name) != "epochs.json" for name in names)


_SAY_AFFECTS = ("say which stages the change altered: --affects {stages} (any of them, comma separated; `plan` means a "
                "galaxy made under the epoch before needs a full regeneration)")


def check_affects_given(affects, root=ROOT):
    """Refuses at once, before the long rewrite of the golden files, when a new epoch is due and `affects` is not given."""
    with open(os.path.join(root, EPOCHS_FILE), encoding="utf-8") as handle:
        data = json.load(handle)
    if data.get("stages") is not None and affects is None and not data["epochs"][-1].get("pending"):
        raise Refused(_SAY_AFFECTS.format(stages=",".join(data["stages"])))


def bump(reason, root=ROOT, environ=None, affects=None):
    """
    Appends the next epoch.

    Returns:
        int: The new epoch, or the current one when this change already bumped it.

    Raises:
        Refused: Under CI, with no reason, or with no changed golden file.
    """
    environ = os.environ if environ is None else environ
    if environ.get("CI"):
        raise Refused("refusing under CI: the golden files are changed on a developer's machine, in the PR")
    if not reason or not reason.strip():
        raise Refused('give the reason for the new epoch: --reason "..."')
    if not _changed_golden_files(root) and _golden_files_exist(root):
        raise Refused(f"no golden file under {GOLDEN_DIR}/ differs from main, so there is nothing to bump "
                      "the epoch for; regenerate the golden files first")
    path = os.path.join(root, EPOCHS_FILE)
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    epochs = data["epochs"]
    stages = data.get("stages")
    if affects is not None and stages is not None and (not affects or any(s not in stages for s in affects)):
        raise Refused(f"--affects names one or more of {', '.join(stages)}, got {','.join(affects) or 'nothing'}")
    if epochs[-1].get("pending"):
        if affects is not None:
            epochs[-1]["affects"] = list(affects)
            _write_epochs(path, data)
        return epochs[-1]["epoch"]
    if stages is not None and affects is None:
        raise Refused(_SAY_AFFECTS.format(stages=",".join(stages)))
    entry = {"epoch": epochs[-1]["epoch"] + 1, "first_release": None, "pending": True, "reason": reason.strip()}
    if affects is not None:
        entry["affects"] = list(affects)
    epochs.append(entry)
    _write_epochs(path, data)
    return epochs[-1]["epoch"]


def _write_epochs(path, data):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(generator_epoch.format_epochs(data))


def write_golden(root=ROOT, environ=None):
    """
    Rewrites the golden files: runs the golden test with `PLANETGEN_WRITE_GOLDEN=1`, which writes them instead of
    comparing.

    Raises:
        Refused: Under CI, or when the golden test fails to run.
    """
    environ = dict(os.environ if environ is None else environ)
    if environ.get("CI"):
        raise Refused("refusing under CI: the golden files are changed on a developer's machine, in the PR")
    environ["PLANETGEN_WRITE_GOLDEN"] = "1"
    environ["PYTHONPATH"] = os.pathsep.join(filter(None, [os.path.join(root, "src"), environ.get("PYTHONPATH")]))
    done = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                           "src/tests/test_golden_seed.py"], cwd=root, env=environ)
    if done.returncode != 0:
        raise Refused("the golden test did not run to the end, so the golden files are not rewritten")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reason", required=True, help="what changed in the generated output, and why")
    parser.add_argument("--affects", help="the stages the change altered, comma separated: plan, scatter, fill, "
                                          "population (required for a new epoch)")
    parser.add_argument("--no-write", action="store_true", help="the golden files were already rewritten")
    args = parser.parse_args(argv)
    affects = [part.strip() for part in args.affects.split(",") if part.strip()] if args.affects else None
    try:
        check_affects_given(affects)
        if not args.no_write:
            write_golden()
        epoch = bump(args.reason, affects=affects)
    except Refused as exc:
        print(f"golden_update: {exc}", file=sys.stderr)
        return 1
    print(f"golden_update: generator epoch {epoch}. Add the line 'generation-output: changed' to this PR's note "
          "under changes/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
