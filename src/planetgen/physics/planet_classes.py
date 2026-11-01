# planetgen/physics/planet_classes.py

"""
Planet and moon classes as a function of the finished body (GEN.645, GEN.647, docs/design/habitability/condition-based-planet-classes.md)
============================================================================================================================

A class says what a world is. It is an output: `classify` reads the conditions a finished planet or moon stores (its mass
and gravity, its mean surface temperature, the pressure and mixture of its air, its water, and the variation of its
starlight) and returns one letter. Generation never draws a class (GEN.646): it builds the body, then calls `classify` on it, so
every generated body's stored letter is this one; the class sweep (`src/tests/class_sweep_cli.py`) reports the shares.

The rules are the rows of `physics/data/planet_classes.json`, in order, one list for giants and one for terrestrial
bodies. The first row whose conditions all hold names the class. The last row of each list has no conditions, so every body
gets exactly one letter (`tests/test_planet_classes.py` also proves it over generated bodies). A row's conditions are
`min` (at least) and `lt` (below) on a named quantity, or `equals` and `not_equals` on a text one.

Quantities (all derived here from the stored columns, named as `db.models` names them):

* `mass_earth`, `gravity_g`, `temperature_k`, `pressure_pa` (0 for no air), `ocean_fraction`, `water_mass_fraction`,
  `hydrosphere`, `flux_variation` (no class reads it since GEN.1135; the At a glance flag does), `ocean_ph` (`None` without an ocean; no class reads it, the
  pigments of life do);
* `share_o2` ... `share_so2`: a stored gas's partial pressure over the surface pressure;
* `reducing_share`: the methane and hydrogen shares together.

Letters no row uses are the `RESERVE`: a letter is never retired without asking Boss, and free letters are kept for new
designations.

Use:

    classify(planet)               # "M"
    classify(row_dict)             # the same, from the stored columns of a planets or moons row
    describe("M")                  # the Class row: name, sentence, kind, basis
"""

import json
import os
from collections.abc import Mapping

from planetgen.physics import atmosphere, elements, unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "planet_classes.json")

TERRESTRIAL = "t"
GIANT = "g"

GASES = tuple(atmosphere.SPECIES)
"""tuple: The gases whose partial pressure a body stores (`share_<gas>` quantities)."""


class ClassError(ValueError):
    """A body that cannot be classified (a quantity the rules read is missing) or a rule file that is refused."""


def values_from_body(body):
    """A planet or moon (in memory) as the column-named dict `classify` reads."""
    values = {
        "body_type": getattr(body, "body_type", None),
        "mass_kg": getattr(body, "mass", None),
        "gravity_g": getattr(body, "gravity", None),
        "surface_temperature_k": getattr(body, "surface_temperature", None),
        "atmospheric_pressure_pa": getattr(body, "atmospheric_pressure", None),
        "hydrosphere": getattr(body, "hydrosphere", None),
        "ocean_fraction": getattr(body, "ocean_fraction", None),
        "water_mass_fraction": getattr(body, "water_mass_fraction", None),
        "flux_variation": getattr(body, "flux_variation", None),
        "ocean_ph": getattr(body, "ocean_ph", None),
        "roundness_ratio": getattr(body, "roundness_ratio", None),
    }
    values[atmosphere.GAS_FIELD] = atmosphere.gas_partials_pa(body)
    return values


def quantities(values):
    """The named quantities the rules read, from a body's stored `values` (a column-named mapping)."""
    body_type = values.get("body_type")
    if body_type not in (TERRESTRIAL, GIANT):
        raise ClassError(f"body_type {body_type!r} is neither 't' nor 'g'")
    mass_kg = values.get("mass_kg")
    if mass_kg is None:
        raise ClassError("the body has no mass")
    pressure = values.get("atmospheric_pressure_pa") or 0.0
    found = {
        "body_type": body_type,
        "mass_earth": _kg_to_mearth(mass_kg),
        "pressure_pa": pressure,
        "gravity_g": values.get("gravity_g"),
        "temperature_k": values.get("surface_temperature_k"),
        "hydrosphere": values.get("hydrosphere") or "dry",
        "ocean_fraction": values.get("ocean_fraction") or 0.0,
        "water_mass_fraction": values.get("water_mass_fraction") or 0.0,
        "flux_variation": values.get("flux_variation") or 0.0,
        "ocean_ph": values.get("ocean_ph"),
        "roundness_ratio": values.get("roundness_ratio"),
    }
    for gas in GASES:
        partial_pa = atmosphere.partial_kpa(values, gas) * 1000.0
        found[f"share_{gas}"] = min(1.0, partial_pa / pressure) if pressure > 0.0 else 0.0
    found["reducing_share"] = found["share_ch4"] + found["share_h2"]
    return found


class Class:
    """One class: a letter, its name and sentence, the body kind it covers and the conditions that select it."""

    __slots__ = ("letter", "name", "kind", "sentence", "basis", "when", "core")

    def __init__(self, row, core):
        self.letter = row["letter"]
        self.name = row["name"]
        self.kind = row["kind"]
        self.sentence = row["sentence"]
        self.basis = row["basis"]
        self.when = tuple(row["when"])
        self.core = self.letter in core
        for condition in self.when:
            keys = set(condition) - {"quantity", "min", "lt", "equals", "not_equals"}
            if keys or "quantity" not in condition or len(condition) < 2:
                raise ClassError(f"class {self.letter}: condition {condition!r} is not quantity plus min, lt or equals")

    def matches(self, found):
        return holds(self.when, found, owner=f"class {self.letter}")


def holds(conditions, found, owner="a rule", missing=ClassError):
    """Whether every condition (`quantity` with `min`, `lt` or `equals`) holds for the named quantities in `found`.
    A quantity the body does not store raises `missing` (a rule that reads it is wrong for the body), or, with
    `missing=None`, fails the condition (the body lacks what the rule asks for, e.g. an ocean's pH on a dry world)."""
    for condition in conditions:
        name = condition["quantity"]
        value = found.get(name)
        if value is None:
            if missing is None:
                return False
            raise missing(f"{owner} reads {name!r}, which the body does not store")
        if "min" in condition and value < condition["min"]:
            return False
        if "lt" in condition and value >= condition["lt"]:
            return False
        if "equals" in condition and value != condition["equals"]:
            return False
        if "not_equals" in condition and value == condition["not_equals"]:
            return False
    return True


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    core = tuple(document["core"])
    classes = tuple(Class(row, core) for row in document["rows"])
    for kind in (TERRESTRIAL, GIANT):
        rows = [c for c in classes if c.kind == kind]
        if not rows or rows[-1].when:
            raise ClassError(f"the {kind!r} rules do not end with a row that has no conditions")
        if any(not c.when for c in rows[:-1]):
            raise ClassError(f"a {kind!r} row before the last has no conditions, so the rows after it never match")
    letters = [c.letter for c in classes]
    if len({(c.letter, c.kind) for c in classes}) != len(letters):
        raise ClassError("a class letter is used by two rows of one kind")
    if any(len({c.name for c in classes if c.letter == letter}) != 1 for letter in letters):
        raise ClassError("the rows of one letter in the two kinds have different names")
    reserve = tuple(document["reserve"])
    if set(reserve) & set(letters):
        raise ClassError("a reserved letter is also a class")
    return document, classes, reserve


DOCUMENT, CLASSES, RESERVE = _load()
"""tuple: The `Class` rows in rule order; `RESERVE` the free letters kept for new designations."""

BY_LETTER = {}
for _c in CLASSES:
    BY_LETTER.setdefault(_c.letter, _c)     # GEN.1142: class T has a giant row and a terrestrial one; the first is the letter's own
del _c
LETTERS = tuple(sorted(BY_LETTER))
CORE_LETTERS = tuple(DOCUMENT["core"])
VERSION = DOCUMENT["version"]


class AsteroidType:
    """The object type of a body under a roundness ratio of 1 (GEN.1160): not a planet class, so no letter. It can be a moon."""

    __slots__ = ("letter", "name", "kind", "sentence", "basis", "when", "core")

    def __init__(self, row):
        self.letter = None
        self.name = row["name"]
        self.kind = row["kind"]
        self.sentence = row["sentence"]
        self.basis = row["basis"]
        self.when = tuple(row["when"])
        self.core = False

    def matches(self, found):
        return found["roundness_ratio"] is not None and holds(self.when, found, owner="the Asteroid type")


ASTEROID = AsteroidType(DOCUMENT["asteroid"])
"""AsteroidType: What a terrestrial body whose roundness ratio is under 1 is. `classify_values` gives it and `classify` gives `None` for it."""


def classify_values(values):
    """The `Class` of a body's stored `values` (a column-named mapping): the first row of its kind whose conditions hold. A terrestrial body under a roundness ratio of 1 is `ASTEROID`, which has no letter."""
    found = quantities(values)
    if found["body_type"] == TERRESTRIAL and ASTEROID.matches(found):
        return ASTEROID
    for cls in CLASSES:
        if cls.kind == found["body_type"] and cls.matches(found):
            return cls
    raise ClassError("no rule matched")  # unreachable: the last row of each kind has no conditions


def classify(body):
    """The class letter of a finished planet or moon (an object, or a column-named mapping of its stored values); `None` for an Asteroid."""
    values = body if isinstance(body, Mapping) else values_from_body(body)
    return classify_values(values).letter


def describe(letter):
    """The `Class` row of `letter`."""
    try:
        return BY_LETTER[letter]
    except KeyError:
        raise ClassError(f"{letter!r} is not a class letter") from None


def _number(value, percent):
    if percent:
        return f"{value * 100:.12g} percent"
    return f"{value:,.12g}"


def condition_words(condition):
    """A rule condition as plain words, e.g. `surface pressure below 10 Pa` or `oxygen share of the air at least 10 percent`."""
    name = condition["quantity"]
    if name.startswith("share_"):
        gas = name[len("share_"):]
        entry = DOCUMENT["quantities"]["share_<gas>"]
        label, unit, percent = entry["label"].replace("<gas>", atmosphere.species_name(gas)), "", True
    else:
        entry = DOCUMENT["quantities"][name]
        label, unit, percent = entry["label"], entry.get("unit", ""), entry.get("percent", False)
    suffix = f" {unit}" if unit else ""
    parts = []
    if "min" in condition:
        parts.append(f"{label} at least {_number(condition['min'], percent)}{suffix}")
    if "lt" in condition:
        parts.append(f"{label} below {_number(condition['lt'], percent)}{suffix}")
    if "equals" in condition:
        parts.append(f"{label} is {condition['equals']}")
    if "not_equals" in condition:
        parts.append(f"{label} is not {condition['not_equals']}")
    return " and ".join(parts)


def rule_words(cls):
    """The conditions of a `Class` as a list of plain phrases ([] for the row that takes the rest of its kind)."""
    return [condition_words(condition) for condition in cls.when]
