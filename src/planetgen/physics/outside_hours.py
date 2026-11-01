# planetgen/physics/outside_hours.py

"""
The hours a person can be outside in one local day (GEN.768, Boss 2026-10-24): the equipment the survival rules list for
a world (`equipment.compute`) says what is worn, and this says for how long it can be worn in a day. The answer is the
smallest of five limits, and the one that binds is named beside it:

* **daylight**: the hours of daylight of the day (`daylight_hours`, GEN.767);
* **dose**: the day's share of the radiation allowance: the hours a year the listed suit allows (`equipment` `outside_hours`)
  over the work days of the year, or no limit when the stored hours are `None`;
* **heat** and **cold**: the time a worn garment keeps the core temperature inside its limits, a one-node heat balance of
  the body under the garment's insulation, its cooling or heating power, the ambient temperature at the day's hottest and
  coldest (`exposure_hours`);
* **air**: the breathing gear's supply against a shift (`protection_limits.supply_hours`);
* **uv**: the day's ultraviolet dose against what the worn garment is rated for: a person may take the band edge the
  garment stands for (`uv_cover_index` bare, `uv_suit_index` in sun-cover clothing, `uv_lethal_index` in a full-cover
  suit, any index in a sealed suit) for the daylight hours, and the dose is linear in the index, so past the edge the
  hours shrink in proportion.

A body that needs a vehicle has 0 hours and the binding limit `vehicle`. The numbers are the named rows of
`protection_limits` (the body, the garments and the shift), each with a source tier and a tweak level in the stand-in
registry of the gap-filling plan (section 6.8, GEN.725); nothing here is drawn, so a world's hours never change except by
its stored columns changing.

Use:

    hours, limit = outside_hours.outside_hours_day(world, levels, daylight_hours, days_per_year, temperatures)

A body's stored columns give the inputs (`for_body`): the equipment levels `equipment.compute` lists, the stored `daylight_hours`
(GEN.767) and the day's coldest and hottest temperature from the cycle of the light (`habitability.temperature_extremes_c`)
widened by half the stored day-night swing of the equator band (`band4_swing_k`). A tidally locked world has no day of its own: a person stays in the
terminator belt, where the star sits on the horizon for the whole 24-hour shift (`LOCKED_BELT_LIGHT_HOURS`) at the belt's stored
temperature (`band4_mean_k`, the middle ring).
"""

import math
from collections import namedtuple

from planetgen.physics import constants, equipment, habitability, regional_temperature
from planetgen.physics.protection_limits import GARMENTS, limit, supply_hours

LIMITS = ("daylight", "dose", "heat", "cold", "air", "uv", "vehicle")
"""tuple: The binding limits by the order of the ties (the first of equal hours wins) and by stored number."""

FIELDS = ("outside_hours_day", "outside_hours_limit")
"""tuple: What a body stores of this module (GEN.768): `outside_hours_day`, the hours of one local day a person can be outside, and
`outside_hours_limit`, the name in `LIMITS` of the limit that sets them. Both are `None` where the stored columns do not say."""

LIMIT_WORDS = {
    "daylight": "the hours of daylight", "dose": "the radiation dose allowance", "heat": "heat", "cold": "cold",
    "air": "the breathing gas it carries", "uv": "ultraviolet light", "vehicle": "a body that needs a vehicle outside",
}
"""dict: A limit and the plain words a page says it in."""

Temperatures = namedtuple("Temperatures", "cold_k hot_k")
"""The coldest and the hottest ambient temperature of the day outside, kelvin."""

LOCKED_BELT_LIGHT_HOURS = 24.0
"""float: The hours of light of a shift in the terminator belt of a tidally locked world, where the star stays near the horizon all the
time (stand-in row OH11: a person who works there has the light the whole shift; the belt's glancing light and its temperature are the
terminator's, not a daily cycle)."""

HOURS_PER_DAY = 24.0
HOURS_PER_YEAR = equipment.HOURS_PER_YEAR
SECONDS_PER_HOUR = 3600.0


def _film_w_m2_k(pressure_kpa, temperature_k):
    """The surface film coefficient of the outside of a garment: radiation to the surroundings plus convection that
    scales with the square root of the gas density relative to Earth's sea level (free convection, Nu ~ Ra^(1/4) gives
    a density power of 1/2)."""
    radiation = 4.0 * limit("garment_emissivity") * constants.STEFAN_BOLTZMANN_CONSTANT * temperature_k ** 3
    density = max(pressure_kpa, 0.0) / limit("earth_pressure_kpa") * limit("earth_temperature_k") / temperature_k
    return radiation + limit("convection_earth_w_m2_k") * math.sqrt(density)


def _resistance(garment, pressure_kpa, temperature_k):
    """The thermal resistance from skin to surroundings, m2 K/W: the garment's insulation plus the film."""
    return garment.insulation_clo * limit("clo_m2_k_w") + 1.0 / _film_w_m2_k(pressure_kpa, temperature_k)


def saturation_kpa(temperature_k):
    """The saturation vapour pressure of water at `temperature_k`, kPa (Magnus, over water above 0 C and over ice below)."""
    celsius = temperature_k - constants.CELSIUS_ZERO_K
    if celsius >= 0.0:
        return 0.61094 * math.exp(17.625 * celsius / (celsius + 243.04))
    return 0.61121 * math.exp(22.587 * celsius / (celsius + 273.86))


def _sweat_w(garment, pressure_kpa, temperature_k, relative_humidity):
    """What the evaporation of sweat can carry away, W (ISO 7933): the vapour pressure between the wet skin and the air over
    the vapour resistance of the garment and the film (the Lewis relation ties the film's to its convection), at most what
    the body can sweat; none inside a sealed garment (the sweat has nowhere to go)."""
    if garment.sealed:
        return 0.0
    skin_k = limit("skin_temperature_c") + constants.CELSIUS_ZERO_K
    wet = saturation_kpa(skin_k) - min(max(relative_humidity, 0.0), 1.0) * saturation_kpa(temperature_k)
    if wet <= 0.0:
        return 0.0
    convection = _film_w_m2_k(pressure_kpa, temperature_k) - 4.0 * limit("garment_emissivity") \
        * constants.STEFAN_BOLTZMANN_CONSTANT * temperature_k ** 3
    lewis = limit("lewis_k_kpa")
    resistance = garment.insulation_clo * limit("clo_m2_k_w") / (limit("clothing_permeability") * lewis) \
        + 1.0 / (lewis * max(convection, 1e-9))
    return min(limit("sweat_max_w"), limit("body_area_m2") * wet / resistance)


def exposure_hours(hot_garment, cold_garment, pressure_kpa, temperatures, relative_humidity):
    """
    `(heat_hours, cold_hours)` a person can work outside before the core temperature leaves its limits
    (`core_hot_limit_c`, `core_cold_limit_c`): in `hot_garment` (a `protection_limits.Garment`, the one the ladder lists at the
    day's hottest temperature) at that temperature, and in `cold_garment` at the coldest (`temperatures`); `None` for a limit
    the body never reaches. The balance is one node: the body's heat capacity against the working metabolism, the dry flow
    through the garment and the film, the evaporation of sweat and the garment's cooling or heating power. In air colder
    than the skin an unsealed garment is worn open (`vent_factor` of its insulation), as a person does who is working.
    """
    area = limit("body_area_m2")
    capacity = limit("body_mass_kg") * limit("body_specific_heat_j_kg_k")
    skin_k = limit("skin_temperature_c") + constants.CELSIUS_ZERO_K
    work = limit("work_metabolism_w")
    core = limit("core_normal_c")
    heat = cold = None
    hot_k, cold_k = temperatures.hot_k, temperatures.cold_k
    if hot_k > 0:
        worn = hot_garment if hot_k >= skin_k or hot_garment.sealed else \
            hot_garment._replace(insulation_clo=hot_garment.insulation_clo * limit("vent_factor"))
        flow = area * (hot_k - skin_k) / _resistance(worn, pressure_kpa, hot_k)
        rate = work + flow - _sweat_w(worn, pressure_kpa, hot_k, relative_humidity) - hot_garment.cooling_w
        if rate > 0:
            heat = capacity * (limit("core_hot_limit_c") - core) / rate / SECONDS_PER_HOUR
    if cold_k > 0:
        rate = area * (skin_k - cold_k) / _resistance(cold_garment, pressure_kpa, cold_k) - work - cold_garment.heating_w
        if rate > 0:
            cold = capacity * (core - limit("core_cold_limit_c")) / rate / SECONDS_PER_HOUR
    return heat, cold


def garment_of(level, kind):
    """The `Garment` worn for Body Protection `level` and `kind` (an `equipment.BODY_KINDS` name); the nearest lower
    level of the same kind, else the level's first row, when the pair has no row of its own."""
    if (level, kind) in GARMENTS:
        return GARMENTS[(level, kind)]
    same = [key for key in GARMENTS if key[1] == kind and key[0] <= level]
    if same:
        return GARMENTS[max(same)]
    return GARMENTS[max(key for key in GARMENTS if key[0] <= level and key[1] == "none")]


def uv_edge(level):
    """The UV index the garment of Body Protection `level` is rated for: the band edge it stands for, `inf` in a sealed
    suit (level 3 and up)."""
    if level >= 3:
        return math.inf
    if level == 2:
        return limit("uv_lethal_index")
    if level == 1:
        return limit("uv_suit_index")
    return limit("uv_cover_index")


def _air_hours(supply):
    return supply_hours(supply)


def work_days_per_year(rotation_period_hours, period_years, locked=False):
    """
    The work days of a year, the local solar days: the year over the solar day `1 / (1 / rotation - 1 / orbit)` of a world that
    turns forward on its own. A tidally locked world (or one with no rotation period) has no day of its own and counts 24-hour
    shifts, the same length a person's day is everywhere.
    """
    year_h = period_years * HOURS_PER_YEAR
    if locked or not rotation_period_hours or rotation_period_hours <= 0:
        return year_h / HOURS_PER_DAY
    turn = 1.0 / rotation_period_hours - 1.0 / year_h
    if turn <= 0:
        return year_h / HOURS_PER_DAY
    return year_h * turn


def outside_hours_day(world, levels, daylight_hours, days_per_year, temperatures):
    """
    `(hours, limit)`: the hours of one local day a person can be outside in the equipment `levels` lists for `world`, and which
    of `LIMITS` sets them. `levels` is the dict of stored equipment levels (`equipment.FIELDS`), `daylight_hours` the day's
    light (GEN.767), `days_per_year` the work days of the year (the local solar days; a tidally locked world counts
    24-hour shifts), `temperatures` the `Temperatures` of the day. A body that needs a vehicle (Body Protection 6) has 0 hours
    and the limit `vehicle`; a day with no daylight has 0 hours and the limit `daylight`.
    """
    level = int(levels.get("eq_body") or 0)
    if level >= 6:
        return 0.0, "vehicle"
    candidates = [("daylight", float(daylight_hours))]
    yearly = levels.get("outside_hours")
    if yearly is not None and days_per_year and days_per_year > 0:
        candidates.append(("dose", yearly / days_per_year))
    hot_garment = garment_of(*equipment.body_at(world, temperatures.hot_k - constants.CELSIUS_ZERO_K))
    cold_garment = garment_of(*equipment.body_at(world, temperatures.cold_k - constants.CELSIUS_ZERO_K))
    heat, cold = exposure_hours(hot_garment, cold_garment, world.pressure_kpa, temperatures, world.relative_humidity)
    if heat is not None:
        candidates.append(("heat", heat))
    if cold is not None:
        candidates.append(("cold", cold))
    air = _air_hours(int(levels.get("eq_supply") or 0))
    if air is not None:
        candidates.append(("air", air))
    index, edge = world.uv_index or 0.0, uv_edge(level)
    if index > edge:
        candidates.append(("uv", float(daylight_hours) * edge / index))
    name, hours = min(candidates, key=lambda item: (item[1], LIMITS.index(item[0])))
    return max(0.0, hours), name


MIDDLE_BAND = 4
"""int: The band of the equator (a locked world's terminator ring), `regional_temperature.column(MIDDLE_BAND, field)`."""


def for_body(world, levels, values):
    """
    `(hours, limit)` of a body from its stored columns, or `(None, None)` when they do not say enough (no rotation period, no
    stored `daylight_hours` on a world that turns on its own, no terminator temperature on a locked one: a row made before
    GEN.767 keeps `None`).

    Args:
        world (habitability.World): The body as `habitability_world.world_from` builds it.
        levels (dict): The equipment levels `equipment.compute(world)` lists for it.
        values (dict): Its column-named values: `rotation_period_hours`, `period_years`, `daylight_hours`, `band4_swing_k`,
            `band4_mean_k`.
    """
    rotation, period = values.get("rotation_period_hours"), values.get("period_years")
    if not rotation or not period:
        return None, None
    if regional_temperature.is_locked(rotation, period):
        belt = values.get(regional_temperature.column(MIDDLE_BAND, "mean_k"))
        if not belt:
            return None, None
        days = work_days_per_year(rotation, period, locked=True)
        return outside_hours_day(world, levels, LOCKED_BELT_LIGHT_HOURS, days, Temperatures(belt, belt))
    daylight = values.get("daylight_hours")
    if daylight is None:
        return None, None
    cold_c, hot_c = habitability.temperature_extremes_c(world)
    half_swing = (values.get(regional_temperature.column(MIDDLE_BAND, "swing_k")) or 0.0) / 2.0
    temperatures = Temperatures(cold_c + constants.CELSIUS_ZERO_K - half_swing, hot_c + constants.CELSIUS_ZERO_K + half_swing)
    return outside_hours_day(world, levels, daylight, work_days_per_year(rotation, period), temperatures)


def fields(world, levels, values):
    """`{FIELDS name: value}` of a body: `for_body`'s hours and limit, both `None` when its columns do not say enough."""
    hours, limit_name = for_body(world, levels, values)
    return {"outside_hours_day": hours, "outside_hours_limit": limit_name}
