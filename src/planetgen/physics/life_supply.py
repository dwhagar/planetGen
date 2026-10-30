# planetgen/physics/life_supply.py

"""
What a world supplies the life that eats chemistry (GEN.887, docs/design/habitability/life-classes.md 8.3)
===========================================================================================================

A chemical life class (`physics/life_classes.py`: the hydrogen eaters of Class 3, the sulfur, iron and nitrogen eaters of Class 4, the
dense-CO2 and cold-methane lives of Classes 6 and 7) lives on a fuel and an oxidant that the world has to deliver. These are the
supplies, each a share of what Earth gives (1 is Earth's), worked out from the columns a body stores, so nothing is stored and an edit
of the body moves them with the rest:

* `h2_supply`: hydrogen from the reaction of water with rock (serpentinisation), from the heat of the water's circulation, the redox
  state of the mantle and the water-rock contact;
* `radiolytic_h2_supply`: hydrogen from the radiolysis of the water in rock, from the ground dose and the contact;
* `donor_supply`: reduced gas and iron (hydrogen sulfide, ferrous iron, hydrogen) from the volcanism of the heat and the mantle's
  redox state;
* `nitrate_supply` and `perchlorate_supply`: the oxidants of Classes 4b and 4c, from the nitrogen in the air and the oxygen,
  and the chloride of the ocean;
* `hydrogen_acetylene_supply`: the haze organics and hydrogen of a methane sea's air (Class 7), from the methane share;
* `sulfur_supply`: the reduced sulfur of the acid eaters (Classes 5 and 9), from the air's sulfur dioxide share (GEN.894); it needs no water
  in contact with rock, since the acid is the medium.

Every number is in `physics/data/life_supply.json` with its tier (T1 to T9) and tweak level (L0 to L4). A body with no liquid water
in contact with rock gets 0 of the rock's supplies. Carbon dioxide's phase (Class 6) is read from the temperature and pressure by
`life_classes.co2_phase`; the acid of Classes 5 and 9 is `physics/acid_medium.py` and the ammonia of Class 8 `physics/strange_solvents.py`.

Use:

    supplies(values)        # {"h2_supply": 0.6, ...} for a column-named mapping of a body's values
"""

import json
import os

from planetgen.physics import atmosphere, elements, habitability, habitability_world

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "life_supply.json")

QUANTITIES = ("h2_supply", "radiolytic_h2_supply", "donor_supply", "nitrate_supply", "perchlorate_supply", "hydrogen_acetylene_supply", "sulfur_supply")
"""tuple: The names of the supplies, as the rule file of the classes reads them."""

REDOXES = ("reduced", "intermediate", "oxidized")
DEFAULT_REDOX = "intermediate"
"""str: The redox state of a body that does not store one (Earth's)."""


class LifeSupplyError(ValueError):
    """A supply file that is refused."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    numbers = document["numbers"]
    for name, entry in numbers.items():
        if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) or not entry.get("basis"):
            raise LifeSupplyError(f"{name}: every number needs a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")
        if "value" not in entry and not all(key in entry for key in REDOXES):
            raise LifeSupplyError(f"{name}: a number is a value or one value for each of {REDOXES}")
    return document, numbers


DOCUMENT, NUMBERS = _load()
VERSION = DOCUMENT["version"]


def number(name):
    """The value of the rule file's number `name`."""
    return NUMBERS[name]["value"]


def redox_factor(name, redox):
    """The multiplier `name` (`h2_redox` or `donor_redox`) gives a mantle of the redox state `redox` (Earth's if unknown)."""
    return NUMBERS[name][redox if redox in REDOXES else DEFAULT_REDOX]


def contact(values):
    """
    How much of the rock's chemistry the water reaches, 0 to 1: the whole of it under a surface, ice-covered or hycean ocean, less for
    an aquifer (by how deep it is) and 0 with no liquid water. A world with a very deep ocean (`deep_water_fraction`) has ice on its floor.
    """
    state = values.get("hydrosphere")
    if state in ("surface ocean", "hycean", "ice-covered ocean"):
        if (values.get("water_mass_fraction") or 0.0) >= number("deep_water_fraction"):
            return number("contact_deep_water")
        return number("contact_ocean")
    ground = habitability_world.groundwater(values)
    if ground == "reachable":
        return number("contact_aquifer_reachable")
    if ground == "deep":
        return number("contact_aquifer_deep")
    return 0.0


def heat_share(values):
    """The heat flow (its own and its tide's) as a share of Earth's, at most `heat_share_max`."""
    heat = (values.get("internal_heat_flux_w_m2") or 0.0) + (values.get("tidal_flux_w_m2") or 0.0)
    return min(number("heat_share_max"), max(0.0, heat) / number("earth_heat_flux_w_m2"))


def habitat_temperature_k(values):
    """
    The temperature, K, of the water the chemical classes live in: the surface's under an open ocean or on a dry surface, the melting point
    of ice under an ice lid, and the surface's plus the rise with depth (the heat flow over the crust's conductivity) in an aquifer. The
    classes that live below the surface read this and not the surface temperature.
    """
    surface = values.get("surface_temperature_k")
    if surface is None:
        return None
    if values.get("hydrosphere") == "ice-covered ocean":
        return max(surface, number("ice_ocean_base_k"))
    if values.get("hydrosphere") in ("surface ocean", "hycean") or habitability_world.groundwater(values) == "none":
        return surface
    heat = (values.get("internal_heat_flux_w_m2") or 0.0) + (values.get("tidal_flux_w_m2") or 0.0)
    depth_m = (values.get("aquifer_depth_km") or 0.0) * 1000.0
    return min(number("habitat_max_k"), surface + depth_m * heat / number("crust_conductivity_w_m_k"))


def supplies(values):
    """`{quantity: share of Earth's}` for every `QUANTITIES` of a rocky body's column-named `values`; 0 for a supply the body lacks the
    water or the chemistry for."""
    reach = contact(values)
    heat = heat_share(values)
    redox = values.get("mantle_redox")
    pressure_kpa = (values.get("atmospheric_pressure_pa") or 0.0) / 1000.0
    o2_kpa = atmosphere.partial_kpa(values, "o2")
    ch4_share = atmosphere.partial_kpa(values, "ch4") / pressure_kpa if pressure_kpa > 0.0 else 0.0
    so2_share = atmosphere.partial_kpa(values, "so2") / pressure_kpa if pressure_kpa > 0.0 else 0.0
    nitrogen = habitability_world.inventories(values)["N"] if values.get("surface_temperature_k") is not None else 0.0
    nitrogen_term = min(1.0, nitrogen / habitability.X_REF_EARTH)
    nitrate = number("nitrate_abiotic") * (1.0 + 2.0 * min(1.0, o2_kpa / number("nitrate_oxygen_kpa")))
    perchlorate_ocean = number("perchlorate_brine") if values.get("ocean_class") == "chloride brine" else number("perchlorate_other")
    return {
        "h2_supply": reach * min(1.0, heat) * redox_factor("h2_redox", redox),
        "radiolytic_h2_supply": reach * min(1.0, (values.get("dose_ground_msv_yr") or 0.0) / number("radiolytic_dose_ref_msv_yr")),
        "donor_supply": reach * heat * redox_factor("donor_redox", redox),
        "nitrate_supply": reach * nitrogen_term * min(1.0, nitrate),
        "perchlorate_supply": reach * number("perchlorate_base") * perchlorate_ocean,
        "hydrogen_acetylene_supply": min(1.0, ch4_share / number("haze_methane_share")),
        "sulfur_supply": min(1.0, so2_share / number("sulfur_so2_ref_share")),
    }
