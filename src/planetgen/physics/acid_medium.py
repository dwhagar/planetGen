# planetgen/physics/acid_medium.py

"""
Where a world holds liquid sulfuric acid (GEN.888, docs/design/habitability/life-classes.md 8.3)
=================================================================================================

Life Classes 5 and 9 live in concentrated sulfuric acid, as the cloud droplets of a thick air or as pools on warm dry ground. Nothing else in the
generator made sulfuric acid, so the `sulfuric` solvent of the habitability index was a dead branch. This module reads a body's stored values
and says whether the acid is there, where, and how strong. Nothing is stored: it is worked out from the air, so an edit moves it with the rest.

The model, with every number in `physics/data/acid_medium.json` and its tier and tweak level:

* **Supply.** Sulfur dioxide in the air (`so2_min_share`), an oxidiser to make it sulfur trioxide (oxygen and carbon dioxide together, `oxidant_min_share`)
  and a trace of water to turn that into the acid (`water_min_share`). A world short of any one makes no acid.
* **Strength.** The acid takes up water until its water activity equals the air's: the water's partial pressure over its saturation pressure at the
  level's temperature (`aerosols.saturation_pa`), the water capped at saturation (a cold trap). The strength is a line in the natural log of the water
  activity through the Venus cloud levels (`reference_*`, `water_activity_slope_per_wt_percent`), held between 0 and 100 percent. The acid counts as the
  medium of life when it is `concentrated_wt_percent` or stronger.
* **Cloud.** The level of the dry adiabat of the air (`aerosols.profile`, the same profile the water and methane clouds use) that is in the liquid
  window (`liquid_low_k` to `liquid_high_k`) and between `cloud_min_pa` and `cloud_max_pa`, and whose acid is concentrated. The warmest such level is the one
  reported.
* **Pool.** Acid at the ground: the surface in the liquid window with a concentrated acid at its own temperature and pressure.

`level(values)` is the number the life class rules read as `acid_medium`: 0 no liquid acid, 1 a cloud layer, 2 a pool (a pool world usually has a cloud too; the
pool is the larger number because a class that needs the ground reads it).
"""

import json
import math
import os

import numpy as np

from planetgen.physics import aerosols, air_column, atmosphere, constants, elements

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "acid_medium.json")

NONE, CLOUD, POOL = 0.0, 1.0, 2.0
"""float: The values of `level`: no liquid acid, a cloud layer, a pool."""


class AcidMediumError(ValueError):
    """A rule file that is refused."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) \
                or not entry.get("basis") or "value" not in entry:
            raise AcidMediumError(f"{name}: every number needs a value, a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")
    return document


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]
VERSION = DOCUMENT["version"]


def number(name):
    """The value of the rule file's number `name`."""
    return NUMBERS[name]["value"]


class AcidMedium:
    """
    What `medium` found.

    Attributes:
        level (float): `NONE`, `CLOUD` or `POOL`.
        layer (str): `"cloud"`, `"pool"` or `None`.
        temperature_k, pressure_pa (float): Where the acid is (the cloud level, or the ground), `None` for none.
        weight_percent (float): The acid's strength there, `None` for none.
    """

    __slots__ = ("level", "layer", "temperature_k", "pressure_pa", "weight_percent")

    def __init__(self, level=NONE, layer=None, temperature_k=None, pressure_pa=None, weight_percent=None):
        self.level, self.layer = level, layer
        self.temperature_k, self.pressure_pa, self.weight_percent = temperature_k, pressure_pa, weight_percent

    def __repr__(self):
        return f"AcidMedium({self.layer}, {self.weight_percent}% at {self.temperature_k} K, {self.pressure_pa} Pa)"


def weight_percent(water_activity):
    """The strength of the acid, percent by weight, in equilibrium with air of `water_activity` (0 to 1): the line in ln(a) of the rule file,
    held between 0 and 100. No water at all is the strongest acid."""
    if water_activity <= 0.0:
        return 100.0
    strength = number("reference_wt_percent") - math.log(water_activity / number("reference_water_activity")) \
        / number("water_activity_slope_per_wt_percent")
    return max(0.0, min(100.0, strength))


def _water_activity(water_pa, temperature_k):
    """The water's partial pressure over its saturation pressure at `temperature_k`, at most 1 (more is a cloud of water, not a drier acid)."""
    saturation = aerosols.saturation_pa("H2O", temperature_k)
    if not math.isfinite(saturation) or saturation <= 0.0:
        return 1.0
    return min(1.0, water_pa / saturation)


def skin_k(values, surface_k):
    """The temperature the air turns isothermal at: 2^-1/4 of the airless grey temperature of the absorbed flux (`energy_flux_w_m2`, the
    reference albedo); the surface temperature's own skin factor when the body gives no flux."""
    flux = values.get("energy_flux_w_m2")
    if not flux or flux <= 0.0:
        return None
    absorbed = flux * (1.0 - number("reference_bond_albedo")) / 4.0
    return aerosols.SKIN_FACTOR * (absorbed / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25


def medium(values):
    """
    The liquid sulfuric acid a body's stored `values` (a column-named mapping with its gases: `atmosphere.GAS_FIELD`) hold, as an `AcidMedium`.
    A body with no air, no gravity or no surface temperature has none.
    """
    gas = atmosphere.gas_partials_pa(values)
    temperature = values.get("surface_temperature_k")
    gravity = (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
    total = sum(gas.values())
    if not gas or not temperature or gravity <= 0.0 or total <= 0.0:
        return AcidMedium()
    share = {formula: pressure / total for formula, pressure in gas.items()}
    if share.get("SO2", 0.0) < number("so2_min_share") or share.get("O2", 0.0) + share.get("CO2", 0.0) < number("oxidant_min_share") \
            or share.get("H2O", 0.0) < number("water_min_share"):
        return AcidMedium()
    column = air_column.from_partials(air_column.partials_pa(values), gravity)
    mixing = share.get("H2O", 0.0)
    low, high = number("liquid_low_k"), number("liquid_high_k")
    strong = number("concentrated_wt_percent")

    def acid_at(pressure_pa, level_k):
        # The water at a level is the surface mixing ratio, but never more than that level's saturation (a cold trap).
        saturation = aerosols.saturation_pa("H2O", level_k)
        water = mixing * pressure_pa
        if math.isfinite(saturation):
            water = min(water, saturation)
        return weight_percent(_water_activity(water, level_k))

    if low <= temperature <= high:
        wt = acid_at(total, temperature)
        if wt >= strong:
            return AcidMedium(POOL, "pool", temperature, total, wt)
    pressure, levels = aerosols.profile(column, temperature, skin_k=skin_k(values, temperature))
    best = None
    for p, t in zip(pressure, levels):
        if not (low <= t <= high and number("cloud_min_pa") <= p <= number("cloud_max_pa")):
            continue
        wt = acid_at(float(p), float(t))
        if wt >= strong and (best is None or t > best[1]):
            best = (float(p), float(t), wt)
    if best is None:
        return AcidMedium()
    return AcidMedium(CLOUD, "cloud", best[1], best[0], best[2])


def level(values):
    """The `acid_medium` the life class rules read: 0 no liquid acid, 1 a cloud layer, 2 a pool (see the module note)."""
    return medium(values).level
