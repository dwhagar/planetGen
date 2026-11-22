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
  and the chloride of the ocean (the abiotic nitrate of lightning and impacts and the perchlorate of the radiolysis of surface chloride do not reach an ocean under an ice lid: `nitrate_sealed_share` and `perchlorate_sealed_share`, GEN.1823);
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

from planetgen.physics import atmosphere, constants, elements, habitability, habitability_world, surface_pressure, water_layer

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
        elements.check_rule_entry(name, entry, LifeSupplyError)
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
    How much of the rock's chemistry the water reaches, 0 to 1: the whole of it under a surface, brine-seep, ice-covered or hycean ocean, less for
    an aquifer (by how deep it is) and 0 with no liquid water. A world with a very deep ocean (`deep_water_fraction`) has ice on its floor.
    """
    state = values.get("hydrosphere")
    if state in ("surface ocean", "brine seep", "hycean", "ice-covered ocean"):
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
    """The heat flow (its own and its tide's) as a share of the Earth twin's in the interior model (`surface_pressure.HEAT_FLUX_EARTH_W_M2`, GEN.1636), at most `heat_share_max`."""
    heat = (values.get("internal_heat_flux_w_m2") or 0.0) + (values.get("tidal_flux_w_m2") or 0.0)
    return min(number("heat_share_max"), max(0.0, heat) / surface_pressure.HEAT_FLUX_EARTH_W_M2)


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


LIQUID_OCEANS = ("surface ocean", "brine seep", "hycean", "ice-covered ocean")
"""tuple: The hydrospheres whose water is a liquid ocean (under a lid or not)."""


def habitat_pressure_pa(values, place="top"):
    """
    The pressure, Pa, of the place a chemical class lives in (GEN.1819, slice 1: a record of the habitat, no gate reads it). The same stored fields as `habitat_temperature_k`:

    * an ocean, open or under a lid: its `top` is the air's pressure plus the weight of the ice lid (`ice_shell_km` of ice at the body's gravity); its `floor` (the rock a hydrogen eater needs) is the pressure at the
      bottom of the water layer (`water_base_pressure_pa`, the compressed column of `water_layer`); where the liquid lies on high-pressure ice (`hp_ice_km`: the ocean is ice-sealed, no rock) the floor of the
      liquid is the pressure at which water at the lid's melting point freezes into that ice (`water_layer.liquidus_pressure_pa`); an ocean with no stored base takes the weight of `ocean_depth_km` of water;
    * an aquifer: the air's pressure plus the weight of the rock to its depth (`crust_density_kg_m3`), for the top and the floor alike;
    * anything else (a dry or ice-covered surface, a vapor world): the air's pressure.

    `None` for a body with no gravity where the weight is needed.
    """
    surface = values.get("atmospheric_pressure_pa") or 0.0
    gravity = (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
    if values.get("hydrosphere") in LIQUID_OCEANS:
        shell_m = (values.get("ice_shell_km") or 0.0) * 1000.0
        top = surface + shell_m * constants.ICE_DENSITY_KG_M3 * gravity
        if place == "top":
            return top
        if values.get("hp_ice_km"):
            return min(values.get("water_base_pressure_pa") or float("inf"), water_layer.liquidus_pressure_pa(number("ice_ocean_base_k")))
        base = values.get("water_base_pressure_pa")
        if base:
            return base
        depth_m = (values.get("ocean_depth_km") or 0.0) * 1000.0
        return top + depth_m * constants.LIQUID_WATER_DENSITY_KG_M3 * gravity
    if habitability_world.groundwater(values) in ("reachable", "deep") and values.get("aquifer_depth_km"):
        return surface + values["aquifer_depth_km"] * 1000.0 * number("crust_density_kg_m3") * gravity
    return surface


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
    sealed = values.get("hydrosphere") == "ice-covered ocean"
    open_share = number("nitrate_sealed_share") if sealed else 1.0       # GEN.1823: lightning and impacts reach an open sea, not one under a lid
    nitrate = number("nitrate_abiotic") * (open_share + 2.0 * min(1.0, o2_kpa / number("nitrate_oxygen_kpa")))
    perchlorate_ocean = number("perchlorate_brine") if values.get("ocean_class") == "chloride brine" else number("perchlorate_other")
    return {
        "h2_supply": reach * min(1.0, heat) * redox_factor("h2_redox", redox),
        "radiolytic_h2_supply": reach * min(1.0, (values.get("dose_ground_msv_yr") or 0.0) / number("radiolytic_dose_ref_msv_yr")),
        "donor_supply": reach * heat * redox_factor("donor_redox", redox),
        "nitrate_supply": reach * nitrogen_term * min(1.0, nitrate),
        "perchlorate_supply": reach * number("perchlorate_base") * perchlorate_ocean * (number("perchlorate_sealed_share") if sealed else 1.0),
        "hydrogen_acetylene_supply": min(1.0, ch4_share / number("haze_methane_share")),
        "sulfur_supply": min(1.0, so2_share / number("sulfur_so2_ref_share")),
    }
