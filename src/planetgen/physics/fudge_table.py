# planetgen/physics/fudge_table.py

"""
The fudge table as a derived file: where it lives, whether it is current, and the refusal when it is not (GEN.1317)
===================================================================================================================

The class fit reads the width of every class edge from the fudge table (`class_fit`, GEN.1313, GEN.1314). The table is a function of three inputs only: the stored stand-in by
end-result matrix (`physics/data/fudge_matrix.json`), the class table (`planet_class_fit.json`) and the stand-in registry rows of the gap plan. Boss's rule (6:28 PM EST, 2026-11-06):
update checks the table against the generator epoch and builds it again when the epoch has changed.

**The file.** `scripts/build_fudge_table.py --update` writes `fudge_table.json` into `PLANETGEN_DATA_DIR` beside the downloaded datasets, with a `header`:

* `epoch`: the generator epoch the table was built at (`util.generator_epoch`),
* `registry_hash`: a hash of the registry rows (every table row of the gap plan that starts `| ID |`),
* `matrix_hash` and `class_hash`: the hashes of the matrix and the class table files,
* `matrix_epoch`: the epoch the matrix was measured at (`null` for the matrix kept in the repo, which was measured before the stamp existed).

The committed `physics/data/fudge_table.json` is the fixture: the same document without a header, which the tests, the goldens and every machine without a built table read. It is rebuilt by
`scripts/build_fudge_table.py` and a test fails when it differs, as for `stand-in-values.md`.

**The check.** `state` compares the header of the file with `current_inputs`: `current` (nothing to do, one hash of three small files), `stale` (a header field differs) or `missing`. A table is
built in a fraction of a second from the stored matrix (the *fast mode*): update never waits for it. The matrix itself is the sensitivity run (about 1.5 hours on four cores, `matrix_age_epochs` says how
old it is) and is refreshed on purpose, not by update.

**The refusal.** Bulk generation and a Generate-page run refuse to start with a missing or stale table and say to run update.sh (`require_current`, called by `run_galaxy.require_datasets`), the way
a missing dataset is refused (GEN.590). Reading, the class pages and the tests are not blocked: `class_fit` takes the widths from `widths_file`, the table when it is current and the fixture otherwise.
Stand-ins are the rows FT6 and FT7 of section 6.261 of the gap plan.
"""

import hashlib
import json
import os
import re

from planetgen.datasets import fetch
from planetgen.physics import elements
from planetgen.util import generator_epoch

NAME = "fudge_table.json"
"""str: The file's name, in the fixture folder and in `PLANETGEN_DATA_DIR`."""

FIXTURE = os.path.join(elements.DATA_DIRECTORY, NAME)
"""str: The committed fixture (no header)."""

MATRIX = os.path.join(elements.DATA_DIRECTORY, "fudge_matrix.json")
CLASSES = os.path.join(elements.DATA_DIRECTORY, "planet_class_fit.json")
PLAN = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "docs", "design", "research", "gap-filling-plan.md"))
"""str: The stand-in registry, read for its rows (`registry_hash`)."""

MATRIX_MAX_AGE_EPOCHS = 50
"""int: A stored matrix measured more than this many epochs ago is called old by `matrix_age_epochs` and update says to refresh it. Source tier T9, tweak level L3 (a round number)."""

_ROW = re.compile(r"^\| ([A-Z]+[0-9]+) \|")


class FudgeTableStale(RuntimeError):
    """The fudge table in the data folder is missing or was built for another epoch or other inputs; the message says what to run."""


def data_file():
    """The path of the built table: `NAME` in the folder the datasets are kept in."""
    return os.path.join(fetch.cache_dir(), NAME)


def _hash_file(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()[:16]


def registry_hash(plan=PLAN):
    """A hash of the registry rows of the gap plan (every table row that starts `| ID |`, in file order), so any edit of a stand-in's tier, tweak, range or text makes the table stale."""
    digest = hashlib.sha256()
    with open(plan, encoding="utf-8") as handle:
        for line in handle:
            if _ROW.match(line):
                digest.update(line.rstrip("\n").encode("utf-8"))
    return digest.hexdigest()[:16]


def current_inputs():
    """The header a table built now would carry (without `matrix_epoch`, which is the matrix's own)."""
    return {"epoch": generator_epoch.generator_epoch(), "registry_hash": registry_hash(), "matrix_hash": _hash_file(MATRIX), "class_hash": _hash_file(CLASSES)}


def header_of(path):
    """The `header` of the table at `path`, or `None` for a file without one, a file that is not there or a file that cannot be read."""
    try:
        with open(path, encoding="utf-8") as handle:
            header = json.load(handle).get("header")
    except (OSError, ValueError):
        return None
    return header if isinstance(header, dict) else None


def state(path=None):
    """
    `(state, reason)` of the built table: `"current"`, `"stale"` or `"missing"`, and a sentence saying which header field differs for a stale one.
    """
    path = path or data_file()
    if not os.path.exists(path):
        return "missing", f"{path} is not there"
    header = header_of(path)
    if header is None:
        return "stale", f"{path} has no header"
    wanted = current_inputs()
    differs = [key for key in wanted if header.get(key) != wanted[key]]
    if differs:
        return "stale", "built for another " + ", ".join(differs) + " (it holds " + ", ".join(f"{key} {header.get(key)}" for key in differs) + ")"
    return "current", ""


def require_current():
    """Raises `FudgeTableStale` when the built table is missing or stale (bulk generation refuses to start, GEN.1317)."""
    found, reason = state()
    if found != "current":
        raise FudgeTableStale(f"the fudge table is {found} ({reason}), so nothing was generated; update.sh builds it "
                              f"('python3 scripts/build_fudge_table.py --update' does it by hand, in under a second)")


def widths_file():
    """The file the class fit reads its widths from: the built table when it is current, the fixture otherwise."""
    return data_file() if state()[0] == "current" else FIXTURE


def matrix_age_epochs(header=None):
    """How many epochs ago the stored matrix was measured (`None` when the table does not say)."""
    header = header if header is not None else header_of(data_file())
    measured = (header or {}).get("matrix_epoch")
    return None if measured is None else generator_epoch.generator_epoch() - int(measured)
