# planetgen/util/generator_epoch.py

"""
The generator epoch (OPS.56; design in `docs/design/generation/generation-determinism.md` section 5).

An integer that goes up by one in a change that alters what the generator produces for a seed, and only then. Two
checkouts with the same epoch (and the same lock) build the same galaxy; a stored galaxy made under an older epoch
says which release to check out to rebuild it (`first_release`).

The record is `src/tests/golden/epochs.json`, beside the golden files whose changes require a bump
(`scripts/bump_version.py --check-pr` holds the rule, `python -m tests.golden_update --reason "..."` makes the bump).
It is read from the checkout, so a deployment must ship the `src/tests/golden` folder (a checkout installed with
`pip install -e` does).
"""

import functools
import json
import re
from pathlib import Path

EPOCHS_FILE = Path(__file__).resolve().parents[2] / "tests" / "golden" / "epochs.json"
"""Path: Where the epochs are recorded."""


@functools.lru_cache(maxsize=1)
def epoch_history():
    """
    Every epoch, oldest first.

    Returns:
        tuple: One dict per epoch with `epoch` (1, 2, 3, ...), `first_release` (the release that first generated
        under it; None for epoch 1, which is every release from the start, and for an epoch whose release is not
        stamped yet) and `reason`.

    Raises:
        FileNotFoundError: The checkout has no `src/tests/golden/epochs.json`.
        ValueError: The epochs are not numbered 1, 2, 3, ... in order.
    """
    entries = json.loads(EPOCHS_FILE.read_text(encoding="utf-8"))["epochs"]
    if [entry["epoch"] for entry in entries] != list(range(1, len(entries) + 1)):
        raise ValueError(f"{EPOCHS_FILE} must list the epochs 1, 2, 3, ... in order")
    return tuple(entries)


def generator_epoch():
    """The current generator epoch, an integer."""
    return epoch_history()[-1]["epoch"]


def first_release(epoch):
    """
    The first release that generated under `epoch`: None for epoch 1 (since the start) and for an epoch not yet
    stamped with its release.

    Raises:
        ValueError: There is no such epoch.
    """
    history = epoch_history()
    if not 1 <= epoch <= len(history):
        raise ValueError(f"no generator epoch {epoch}; the epochs are 1 to {len(history)}")
    return history[epoch - 1].get("first_release")


def stages():
    """The stages of generation an epoch can change, in the order a galaxy is made: `plan`, `scatter`, `fill`,
    `population` (the `stages` list of `epochs.json`)."""
    return tuple(json.loads(EPOCHS_FILE.read_text(encoding="utf-8"))["stages"])


def affects(epoch):
    """
    The stages the change that made `epoch` altered (GEN.764): a subset of `stages()`, empty for epoch 1.

    Raises:
        ValueError: There is no such epoch, or the entry names no stage list.
    """
    history = epoch_history()
    if not 1 <= epoch <= len(history):
        raise ValueError(f"no generator epoch {epoch}; the epochs are 1 to {len(history)}")
    entry = history[epoch - 1]
    if "affects" not in entry:
        raise ValueError(f"epoch {epoch} in {EPOCHS_FILE} has no `affects` list")
    return tuple(entry["affects"])


def stages_between(recorded, current=None):
    """
    The stages a result made under epoch `recorded` has to redo to match epoch `current` (default the running one):
    the union of the `affects` of the epochs after `recorded`, up to and including `current`, in the order of
    `stages()`. Empty when `recorded` is `current`.

    Raises:
        ValueError: `recorded` is after `current`, or an epoch is unknown.
    """
    current = generator_epoch() if current is None else current
    if recorded > current:
        raise ValueError(f"a result made under epoch {recorded} is newer than epoch {current}")
    touched = {stage for epoch in range(recorded + 1, current + 1) for stage in affects(epoch)}
    return tuple(stage for stage in stages() if stage in touched)


def needs_full_regeneration(recorded, current=None):
    """Whether a galaxy made under epoch `recorded` cannot be brought to `current` in part: an epoch since changed
    the plan (GEN.764)."""
    return "plan" in stages_between(recorded, current)


_SHORT_LIST = re.compile(r'("(?:affects|stages)": )\[\n((?:\s+"[a-z]+",?\n)+)\s+\]')


def format_epochs(data):
    """
    The text of `epochs.json` for `data`: indented JSON with the short `affects` and `stages` lists on one line each,
    so every writer (the epoch bump, the release stamp) leaves the file the same way.
    """
    def collapse(match):
        return match.group(1) + "[" + ", ".join(re.findall(r'"[a-z]+"', match.group(2))) + "]"
    return _SHORT_LIST.sub(collapse, json.dumps(data, indent=2)) + "\n"
