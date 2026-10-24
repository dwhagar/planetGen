# planetgen/physics/design_flags.py

"""
The design-flag registry (GEN.525, docs/design/design-flags.md)
===============================================================

Some rows of the atmosphere recipes stand for a quantity no ensemble can define: no exoplanet has a measured nitrous oxide,
phosphine or nitrogen dioxide, and whether a world makes them depends on whether it has life. Decision D4 of
`docs/design/tunable-gaps.md` keeps each as a labelled design bracket that is off unless a named flag (the life stage, the
stellar energetic particles, lightning) allows it, and never tunes it. The data file `design_flags.json` lists those rows; this
module reads it and `check` proves the recipes honour it: the row exists, is labelled as a design part, takes its bracket
from the constant the registry names, is zero when none of its enabling flags is set, is not marked tunable and has no target
in the tuning-target table.

Use:

    flag("m3").enabled_by              # ("life_chemical", "sep_reducing")
    check()                            # () when the recipes honour the registry, else the problems
"""

import json
import os
from dataclasses import dataclass

from planetgen.physics import elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "design_flags.json")
"""str: The data file."""


@dataclass(frozen=True)
class DesignFlag:
    """
    One registered design flag.

    Attributes:
        row, table (str): The recipe row and the table it is in (`atmosphere-<letter>`).
        stands_for (str): What the row's presence stands for.
        constant (str): The design constant of the row that holds the bracket.
        enabled_by (tuple): The subject quantities (flags) of which at least one must be true for the row to be non-zero.
        why_not_tuned (str): Why no statistic can define the bracket.
        label (str): The source label (always `DESIGN`).
    """

    row: str
    table: str
    stands_for: str
    constant: str
    enabled_by: tuple
    why_not_tuned: str
    label: str


def load(path=DATA_FILE):
    """Reads and checks the registry at `path`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the design-flag registry {path}: {error}") from error
    flags = {}
    for raw in document.get("flags") or []:
        owner = f"design flag {raw.get('row')!r}"
        for key in ("row", "table", "stands_for", "constant", "why_not_tuned", "label"):
            if not isinstance(raw.get(key), str) or not raw[key]:
                raise DataError(f"{owner}: {key} is missing")
        enabled = raw.get("enabled_by")
        if not isinstance(enabled, list) or not enabled or not all(isinstance(q, str) and q for q in enabled):
            raise DataError(f"{owner}: enabled_by must name at least one flag")
        if raw["label"] != "DESIGN":
            raise DataError(f"{owner}: the label of a design flag is DESIGN")
        if raw["row"] in flags:
            raise DataError(f"{owner}: listed twice")
        flags[raw["row"]] = DesignFlag(raw["row"], raw["table"], raw["stands_for"], raw["constant"], tuple(enabled),
                                       raw["why_not_tuned"], raw["label"])
    if not flags:
        raise DataError("the design-flag registry is empty")
    return flags


FLAGS = load()
"""dict: Row id -> `DesignFlag`."""


def flag(row_id):
    """The registered design flag of recipe row `row_id` (`KeyError` if it is none)."""
    return FLAGS[row_id]


def check(recipes=None, targets=None):
    """
    The problems found when the recipes and the tuning-target table are held against the registry: a tuple of sentences, empty
    when they honour it.
    """
    if recipes is None:
        from planetgen.physics.atmosphere_recipes import RECIPES as recipes
    if targets is None:
        from planetgen.generation.tuning_harness import TARGETS as targets
    from planetgen.generation.possibility_table import Subject, VERIFIED
    problems = []
    rows_with_targets = {t.row for t in targets}
    for item in FLAGS.values():
        try:
            entry = recipes.entry(item.row)
        except KeyError:
            problems.append(f"{item.row}: not a row of the atmosphere recipes")
            continue
        if entry.state == VERIFIED:
            problems.append(f"{item.row}: a design flag cannot be verified")
        if item.constant not in entry.constants:
            problems.append(f"{item.row}: has no design constant {item.constant}")
        if entry.tunable:
            problems.append(f"{item.row}: a design flag is never tunable")
        if f"{item.table}:{item.row}" in rows_with_targets:
            problems.append(f"{item.row}: a design flag has a target in the tuning-target table")
        for quantity in item.enabled_by:
            if quantity not in entry.weight.text if entry.weight is not None else "":
                problems.append(f"{item.row}: the weight does not read the flag {quantity}")
        off = Subject("A", {q: False for q in item.enabled_by})
        try:
            if recipes.raw_weight(entry, off, gated=False) != 0:
                problems.append(f"{item.row}: not zero when no flag is set")
        except Exception as error:  # a weight that cannot be read without its flags is also a problem
            problems.append(f"{item.row}: cannot be read with no flag set ({error})")
    return tuple(problems)
