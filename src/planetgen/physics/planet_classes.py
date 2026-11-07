# planetgen/physics/planet_classes.py

"""
Planet and moon classes as a function of the finished body (GEN.645, GEN.647, GEN.1207, docs/design/habitability/condition-based-planet-classes.md)
==============================================================================================================================================

A class says what a world is. It is an output: `classify` reads the conditions a finished planet or moon stores (its mass
and gravity, its mean surface temperature, the pressure and mixture of its air, its water, and the variation of its
starlight) and returns one letter. Generation never draws a class (GEN.646): it builds the body, then calls `classify` on it, so
every generated body's stored letter is this one; the class sweep (`src/tests/class_sweep_cli.py`) reports the shares.

Since GEN.1207 the class is the **best fit**: every class has its own positive conditions in `physics/data/planet_class_fit.json`, each is graded, and
`physics/class_fit.py` picks the class with the highest percent match (no row relies on another having said no). `classify_detail` also gives the match, the
runner-up and its match, and the digit of the class code. `physics/data/planet_classes.json` holds what a class is called and says (name, sentence, basis), the
Asteroid object type (a body under a roundness ratio of 1, which has no letter) and the quantities the conditions are written in.

Quantities (all derived here from the stored columns, named as `db.models` names them):

* `mass_earth`, `gravity_g`, `temperature_k`, `pressure_pa` (0 for no air), `ocean_fraction`, `water_mass_fraction`,
  `hydrosphere`, `flux_variation` (no class reads it since GEN.1135; the At a glance flag does), `ocean_ph` (`None` without an ocean; no class reads it, the
  pigments of life do);
* `share_o2` ... `share_so2`: a stored gas's partial pressure over the surface pressure; `o2_partial_pa`: the oxygen's own partial pressure, Pa (GEN.1321: the classes judge oxygen by pressure);
* `reducing_share`: the methane and hydrogen shares together.

Letters no class uses are the `RESERVE`: a letter is never retired without asking Boss, and free letters are kept for new
designations.

Use:

    classify(planet)               # "M"
    classify(row_dict)             # the same, from the stored columns of a planets or moons row
    classify_detail(planet)        # Classification("M", "t", 0.97, "O", 0.41, 0)
    describe("M")                  # the Class row: name, sentence, kind, basis
"""

import json
import os
from collections.abc import Mapping

from planetgen.physics import atmosphere, elements, unit_registry as units
from planetgen.util.wordlist import join_words

_kg_to_mearth = units.converter("kg", "Mearth")

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "planet_classes.json")

TERRESTRIAL = "t"
GIANT = "g"

BAND_COUNT = 7
"""int: How many latitude bands a body stores (`regional_temperature`, `latitude_bands`)."""

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
        "relative_humidity": getattr(body, "relative_humidity", None),
        "groundwater": getattr(body, "groundwater", None),
        "aquifer_depth_km": getattr(body, "aquifer_depth_km", None),
        "aquifer_thickness_km": getattr(body, "aquifer_thickness_km", None),
        "ice_mass_fraction": getattr(body, "ice_mass_fraction", None),
        "surface_liquid_ammonia_wt": getattr(body, "surface_liquid_ammonia_wt", None),
    }
    for number in range(1, BAND_COUNT + 1):         # the seven latitude bands (GEN.1108), which the regional, day and night temperatures of a class read (GEN.1221)
        for field in ("mean_k", "swing_k"):
            values[f"band{number}_{field}"] = getattr(body, f"band{number}_{field}", None)
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
    found["o2_partial_pa"] = atmosphere.partial_kpa(values, "o2") * 1000.0       # GEN.1321: breathable oxygen is a pressure, not a share
    return found


class Class:
    """One class: a letter, its name and sentence and the body kind it covers. Its conditions are its row of `class_fit` (`class_fit.CLASSES`)."""

    __slots__ = ("letter", "name", "kind", "sentence", "basis", "core")

    def __init__(self, row, core):
        self.letter = row["letter"]
        self.name = row["name"]
        self.kind = row["kind"]
        self.sentence = row["sentence"]
        self.basis = row["basis"]
        self.core = self.letter in core


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
"""tuple: The `Class` rows; `RESERVE` the free letters kept for new designations."""

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


def classify_detail(body):
    """
    The `class_fit.Classification` of a finished planet or moon (an object, or a column-named mapping of its stored values): its letter and kind, the percent
    match as a fraction, the runner-up letter and match, and the digit of the class code. A terrestrial body under a roundness ratio of 1 is `ASTEROID`, which
    has no letter, and is returned as `None`.
    """
    values = body if isinstance(body, Mapping) else values_from_body(body)
    found = quantities(values)
    if found["body_type"] == TERRESTRIAL and ASTEROID.matches(found):
        return None
    from planetgen.physics import class_fit        # class_fit reads this module's quantities
    return class_fit.classify(values)


def classify_values(values):
    """The `Class` of a body's stored `values` (a column-named mapping): its best-fit class. A terrestrial body under a roundness ratio of 1 is `ASTEROID`, which has no letter."""
    detail = classify_detail(values)
    if detail is None:
        return ASTEROID
    return class_of(detail.letter, detail.kind)


def classify(body):
    """The class letter of a finished planet or moon (an object, or a column-named mapping of its stored values); `None` for an Asteroid."""
    detail = classify_detail(body)
    return None if detail is None else detail.letter


def class_of(letter, kind):
    """The `Class` row of a letter in a body kind."""
    for cls in CLASSES:
        if cls.letter == letter and cls.kind == kind:
            return cls
    raise ClassError(f"{letter!r} is not a class of kind {kind!r}")


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
    elif name.startswith("phase_"):
        entry = DOCUMENT["quantities"]["phase_<family>"]
        label, unit, percent = entry["label"].replace("<family>", name[len("phase_"):].replace("_", "-")), "", False
    else:
        entry = DOCUMENT["quantities"][name]
        label, unit, percent = entry["label"], entry.get("unit", ""), entry.get("percent", False)
    suffix = f" {unit}" if unit else ""
    parts = []
    if "min" in condition:
        parts.append(f"{label} at least {_number(condition['min'], percent)}{suffix}")
    if "lt" in condition:
        parts.append(f"{label} below {_number(condition['lt'], percent)}{suffix}")
    if "equals" in condition and name.startswith("phase_") and condition["equals"] == 1.0:
        parts.append(f"the {name[len('phase_'):].replace('_', '-')} family's liquid is a liquid at the surface temperature and pressure")
    elif "equals" in condition and name == "groundwater_present" and condition["equals"] == 1.0:
        parts.append("liquid groundwater under a dry surface")
    elif "equals" in condition:
        parts.append(f"{label} is {condition['equals']}")
    if "not_equals" in condition:
        parts.append(f"{label} is not {condition['not_equals']}")
    return join_words(parts)
