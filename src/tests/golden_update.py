# tests/golden_update.py

"""
Bumps the generator epoch after the golden files were changed on purpose (OPS.56):

    python -m tests.golden_update --reason "what changed in the generated output and why"

A PR that changes what the generator produces for a seed changes the golden files under `src/tests/golden/`, and
`scripts/bump_version.py --check-pr` then requires three things together: those files, the epoch bumped by one, and
a line `generation-output: changed` in the PR's `changes/` note. This tool does the second: it appends the next epoch
to `src/tests/golden/epochs.json`, with the reason and no release yet (stamping the release names it).

It refuses under CI (a golden file is never to be regenerated there, or a regression would rewrite its own proof),
when the golden files have not changed since main (nothing to bump for) and for a second bump in the same change (the epoch
already went up; the new golden changes ride on it).
"""

import argparse
import json
import os
import subprocess
import sys

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


def bump(reason, root=ROOT, environ=None):
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
    if not _changed_golden_files(root):
        raise Refused(f"no golden file under {GOLDEN_DIR}/ differs from main, so there is nothing to bump "
                      "the epoch for; regenerate the golden files first")
    path = os.path.join(root, EPOCHS_FILE)
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    epochs = data["epochs"]
    if epochs[-1].get("pending"):
        return epochs[-1]["epoch"]
    epochs.append({"epoch": epochs[-1]["epoch"] + 1, "first_release": None, "pending": True,
                   "reason": reason.strip()})
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")
    return epochs[-1]["epoch"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reason", required=True, help="what changed in the generated output, and why")
    args = parser.parse_args(argv)
    try:
        epoch = bump(args.reason)
    except Refused as exc:
        print(f"golden_update: {exc}", file=sys.stderr)
        return 1
    print(f"golden_update: generator epoch {epoch}. Add the line 'generation-output: changed' to this PR's note "
          "under changes/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
