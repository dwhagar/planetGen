# planetgen/physics/age_needs.py

"""
How old a star must be to host a body (GEN.686, docs/design/habitability/condition-based-planet-classes.md)
=========================================================================================================

A planet or moon needs the time its settling takes: a young star can host a giant, an airless rock, a small barren world
or a molten one, but a world with air and a settled surface (and one with liquid water most of all) needs the star to be
old enough, and longer when the star evolves slowly. The rows of `physics/data/age_needs.json` name that time as a tier for
the body's stored conditions, first match wins, and the last row has none, so every body gets a tier. It reads the finished
body's conditions, not a class letter: the class is an output (`planet_classes`), and so is this.

Use:

    tier(planet)                              # "watered"
    minimum_age_gy(planet, ["fast", "normal"])  # 4.0
"""

import json
import os

from planetgen.physics import elements, planet_classes

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "age_needs.json")


class AgeNeedsError(ValueError):
    """A rule file that is refused or a body the rules cannot read."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    tiers, rows = document["tiers"], document["rows"]
    if not rows or rows[-1]["when"]:
        raise AgeNeedsError("the age rules do not end with a row that has no conditions")
    if any(not row["when"] for row in rows[:-1]):
        raise AgeNeedsError("an age row before the last has no conditions, so the rows after it never match")
    for row in rows:
        if row["tier"] not in tiers:
            raise AgeNeedsError(f"age row {row['name']!r} names the unknown tier {row['tier']!r}")
    return document


DOCUMENT = _load()
TIERS = DOCUMENT["tiers"]
ROWS = DOCUMENT["rows"]
SCALES = ("fast", "normal", "slow")


def tier_of_values(values):
    """The tier name for a column-named mapping of stored values (as `planet_classes.values_from_body` builds)."""
    found = planet_classes.quantities(values)
    for row in ROWS:
        if planet_classes.holds(row["when"], found, f"age row {row['name']!r}"):
            return row["tier"]
    raise AgeNeedsError("no age row matched")  # unreachable: the last row has no conditions


def tier(body):
    """The tier name for a finished planet or moon."""
    return tier_of_values(planet_classes.values_from_body(body))


def minimum_age_gy(body, scales):
    """The youngest age (Gy) a star evolving on any of `scales` can have and still host `body`: the highest of the body's
    tier ages on those scales (0.0 for none)."""
    tier_ages = TIERS[tier(body)]
    return max((tier_ages[scale] for scale in scales if scale in tier_ages), default=0.0)
