# planetgen/physics/strange_solvents.py

"""
The ammonia of a body's water and its silicon-rich crust (GEN.894, docs/design/habitability/life-classes.md 8.3)
================================================================================================================

Class 8 (life in ammonia-water) gates on `ammonia_water_mass_fraction`, Class 9 (silicon-bearing life in acid) on `silicon_rich_crust`. Both are worked
out from the body's bulk composition (`bulk_composition`, GEN.428), so nothing is stored and an edit moves them with the rest. Every number is in
`physics/data/strange_solvents.json` with its tier and tweak level.

* **Ammonia.** The share of the water's mass that is ammonia: `nh3_share_of_ice` of the body's non-water ice over its water, so a body with as
  much other ice as water holds a tenth of the water's mass as ammonia. 0 for a body with no water; `None` for a body saved before the bulk
  composition (the class then names the quantity as missing).
* **Silicon.** 1 when the rock's molar Mg/Si is at or under `silicon_rich_mg_si_max`, else 0; `None` when the body stores no ratio.
"""

import json
import os

from planetgen.physics import elements

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "strange_solvents.json")


class StrangeSolventError(ValueError):
    """A rule file that is refused."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) \
                or not entry.get("basis") or "value" not in entry:
            raise StrangeSolventError(f"{name}: every number needs a value, a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")
    return document


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]
VERSION = DOCUMENT["version"]


def number(name):
    """The value of the rule file's number `name`."""
    return NUMBERS[name]["value"]


def ammonia_water_mass_fraction(values):
    """The ammonia, as a share of the water's mass, of a body's column-named `values` (`water_mass_fraction`, `ice_mass_fraction`); 0 with no
    water, `None` when the body stores no ice share."""
    ice = values.get("ice_mass_fraction")
    water = values.get("water_mass_fraction")
    if ice is None or water is None:
        return None
    if water <= 0.0:
        return 0.0
    return number("nh3_share_of_ice") * ice / water


def silicon_rich_crust(values):
    """1.0 when the rock's molar Mg/Si (`rock_mg_si_ratio`) is at or under `silicon_rich_mg_si_max`, else 0.0; `None` when none is stored."""
    ratio = values.get("rock_mg_si_ratio")
    if ratio is None:
        return None
    return 1.0 if ratio <= number("silicon_rich_mg_si_max") else 0.0
