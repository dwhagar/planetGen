# planetgen/physics/species_limits.py

"""
What a species can live in (POP.12)
====================================

A species' environment limits are the human visitor's table (`protection_limits.LIMIT_ROWS`) with some rows
replaced: `overrides`, a dict `{limit name: value}`, empty for a species with the human limits. The six
habitability domains are then scored against them (`habitability.phi4` with `habitability.bands(overrides)`),
and nothing else in the score changes (design note habitability-index.md section 10).

`generated` makes a species' limits from the conditions of its homeworld, so every species passes its own
homeworld (Boss's default of 2026-10-15). `assess` says whether a species could live on a world, unaided or with
the equipment its tech level can make.

This module reads no database; callers pass the `habitability.World`.
"""

import math

from planetgen import tuning
from planetgen.physics import habitability, habitability_world
from planetgen.physics.protection_limits import LIMITS, limit

# Edges of each side of each band, outermost first. A side moves as a whole, so the spacing between its edges
# (the Green and Yellow widths) is the human one.
_PRESSURE_LOW = ("pressure_low_outer_kpa", "armstrong_kpa", "mask_min_kpa", "pressure_no_gear_kpa")
_PRESSURE_HIGH = ("pressure_blue_max_kpa", "pressure_narcosis_kpa", "pressure_hyperbaric_max_kpa", "pressure_outer_kpa")
_COLD = ("cold_outer_c", "cold_suit_c", "cold_clothing_c", "comfort_low_c")
_HOT = ("comfort_high_c", "heat_garment_c", "heat_max_c", "heat_shielded_c", "hot_outer_c")
_DOSE = ("dose_blue_sv_yr", "dose_green_sv_yr", "dose_yellow_sv_yr", "dose_outer_sv_yr", "dose_single_year_sv")

# Bands in order (each must not fall as it goes inward), for `validate`.
_ORDERED = (
    ("pressure", _PRESSURE_LOW + _PRESSURE_HIGH),
    ("temperature", _COLD + _HOT),
    ("dose", _DOSE[:4]),
)


def validate(overrides):
    """
    Raises `ValueError` saying what is wrong with `overrides` (`{limit name: value}`): a name that is not a limit,
    a value that is not a finite number above zero (a temperature may be negative), or edges of a band out of order.
    Returns `overrides` as a new dict of floats.
    """
    out = {}
    for name, value in (overrides or {}).items():
        if name not in LIMITS:
            raise ValueError(f"{name!r} is not an environment limit")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{name} must be a number, got {value!r}")
        if LIMITS[name].unit != "C" and value <= 0:
            raise ValueError(f"{name} must be above zero, got {value!r}")
        out[name] = float(value)
    for label, names in _ORDERED:
        values = [(name, out.get(name, limit(name))) for name in names]
        for (first, low), (second, high) in zip(values, values[1:]):
            if low > high:
                raise ValueError(f"{label} limits are out of order: {first} ({low:g}) is above {second} ({high:g})")
    return out


def _shifted(names, factor=None, offset=None):
    """`{name: value}` of the human rows `names` multiplied by `factor` or moved by `offset`."""
    if factor is not None:
        return {name: limit(name) * factor for name in names}
    return {name: limit(name) + offset for name in names}


def generated(world):
    """
    The `overrides` that give a species the limits its homeworld (`habitability.World`) needs: each band the
    human limits already hold the homeworld inside, with a margin (`tuning.SPECIES_LIMIT_MARGIN`), is left as it is;
    one that does not is widened on the side the homeworld lies past, every edge of that side moving by the same
    factor (pressure, dose) or amount (temperature), so the Green and Yellow bands keep their human width.
    Oxygen follows the inspired O2 of the homeworld's air. Only rows that differ from the human value are returned.
    """
    out = {}
    margin = tuning.SPECIES_LIMIT_MARGIN
    # Pressure (log): scale a side by the factor that puts the homeworld inside Blue with the margin.
    pressure = world.pressure_kpa
    if pressure > 0.0:
        factor = (pressure / margin["pressure_factor"]) / limit("pressure_no_gear_kpa")
        if factor < 1.0:
            out.update(_shifted(_PRESSURE_LOW, factor=factor))
        factor = (pressure * margin["pressure_factor"]) / limit("pressure_blue_max_kpa")
        if factor > 1.0:
            out.update(_shifted(_PRESSURE_HIGH, factor=factor))
    # Temperature (linear): move a side by the number of degrees the homeworld's cycle lies past Blue.
    cold, hot = habitability.temperature_extremes_c(world)
    offset = (cold - margin["temperature_c"]) - limit("comfort_low_c")
    if offset < 0.0:
        out.update(_shifted(_COLD, offset=offset))
    offset = (hot + margin["temperature_c"]) - limit("comfort_high_c")
    if offset > 0.0:
        out.update(_shifted(_HOT, offset=offset))
    if -20.0 <= hot <= 50.0:
        wet = habitability.wet_bulb_c(hot, world.relative_humidity)
        if wet >= limit("wet_bulb_c"):
            out["wet_bulb_c"] = wet + margin["wet_bulb_c"]
    # Dose (log): all the edges scale together.
    if world.dose_sv_yr > 0.0:
        factor = (world.dose_sv_yr * margin["dose_factor"]) / limit("dose_blue_sv_yr")
        if factor > 1.0:
            out.update(_shifted(_DOSE, factor=factor))
    # Oxygen: inspired O2 of the homeworld's air, with the species' own Armstrong edge.
    here = habitability.bands(out)
    inspired = habitability.inspired_o2_kpa(world, here)
    if 0.0 < inspired < 2.0 * limit("inspired_o2_min_kpa"):
        out["inspired_o2_min_kpa"] = inspired / 2.0
    chronic = inspired * margin["o2_factor"]
    if chronic > limit("inspired_o2_chronic_kpa"):
        ratio = limit("inspired_o2_acute_kpa") / limit("inspired_o2_chronic_kpa")
        out["inspired_o2_chronic_kpa"] = chronic
        out["inspired_o2_acute_kpa"] = chronic * ratio
    return {name: round(value, 6) for name, value in out.items() if value != limit(name)}


def assess(world, overrides, tech_level=None, native=False):
    """
    Whether a species with the limits `overrides` could live on `world`, and what it would need.

    Returns:
        dict: `domains` (`{domain: (score, tier)}` for the six habitability domains and `"PHI-4"` the overall score,
            as `habitability.phi4`), `equipment` (the kit its limits ask for, `habitability.equipment`), `unaided`
            (no domain Red), `equipped` (it can live there with the kit: Water is not Red, and the species' `tech_level`
            reaches the kit's `tuning.SPECIES_EQUIPMENT_MIN_TECH`; a species with no tech level has no kit) and
            `failing` (the domains that are Red, in page order).

    `native` is true for the species' own homeworld: Suitability is not judged there, because the water and the food
    of its world are the species' own (POP.14); a Red Suitability alone no longer keeps it off its own world.
    """
    bands = habitability.bands(overrides)
    domains = habitability.phi4(world, bands)
    kit = habitability.equipment(world, bands)
    failing = [name for name in habitability.DOMAINS if domains[name][1] == "Red" and not (native and name == "Suitability")]
    tier = habitability_world.EQUIPMENT_NAMES.index(kit)
    equipped = (domains["Water"][1] != "Red" and tech_level is not None
                and tech_level >= tuning.SPECIES_EQUIPMENT_MIN_TECH[tier])
    return {"domains": domains, "equipment": kit, "unaided": not failing, "equipped": equipped, "failing": failing}


def explain(world, overrides, failing):
    """
    The `failing` domains (`assess`) of `world` in plain in-universe words, one sentence each, for the species whose
    limits are `overrides`: "the air is too thin", "it is too hot" and so on. No numbers, so a message is the same
    whatever the units.
    """
    bands = habitability.bands(overrides)
    cold, hot = habitability.temperature_extremes_c(world)
    reasons = {
        "Pressure": "the air is too thin to hold" if world.pressure_kpa < bands.pressure[4] else "the air is crushingly dense",
        "Temperature": ("it is too cold" if cold < bands.temperature[3] and hot <= bands.temperature[4] else
                        "it is too hot" if hot > bands.temperature[4] and cold >= bands.temperature[3] else
                        "it swings too far between hot and cold"),
        "Air": "the air cannot be breathed",
        "Water": "there is no water it could drink",
        "Radiation": "the radiation is far beyond what it can survive",
        "Suitability": "the surface offers too little water, food and shelter to live on",
    }
    return [reasons[name] for name in failing]
