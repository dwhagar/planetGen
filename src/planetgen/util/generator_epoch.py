# planetgen/util/generator_epoch.py

"""
The generator epoch (OPS.56; design in `docs/design/generation-determinism.md` section 5).

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
