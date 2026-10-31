# planetgen/physics/protection.py

"""
The one habitability describer (GEN.218): the protective gear a visitor needs on a world, and the rows and the
paragraph that say where each of the six PHI-6 categories falls and why (GEN.208, GEN.216), in Boss's short wording.
`rating` returns the rows (domain, colour, score, sentence) that `describe` shows as a table; `paragraph` is the
same rows as one paragraph for the text formats. Numbers go through `util/format.three_figures_or_tiny`, the
formatter the unit ladders use, and the gases are named as the gas table names them (`atmosphere.species_name`).

The gear ladder is the research notes' (Environmental.Extremes.and.Protection.Research.md,
Thermal.Respiratory.Protection.Research.md): modular gear (clothing, masks), intermediate
systems (counterpressure and chemical-barrier suits) and fully enclosed suits, one ladder
per stressor, with the limits of `protection_limits` (one table, the more permissive value
wherever the notes and the built table differ). A suit built for one stressor is not one built
for another: a heated environment suit does not cool, a cooled one does not shield from radiation,
and solar-storm doses need a shelter, not a suit.

Everything here is a pure function of a `habitability.World`.
"""

from dataclasses import replace

from planetgen.physics import atmosphere, equipment, habitability, latitude_bands, regional_temperature
from planetgen.physics.habitability import (
    ARMSTRONG_KPA, MASK_MIN_KPA, gas_inputs, inspired_o2_kpa, tier, water_inputs, wet_bulb_c)
from planetgen.physics.protection_limits import limit
from planetgen.util.format import three_figures_or_tiny


def _gas(symbol):
    """A gas in a sentence: `CO2` is "carbon dioxide" (the gas table's name)."""
    return atmosphere.species_name(symbol.lower())


def _join(items, conjunction="and"):
    """`a`, `a and b`, `a, b and c`."""
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + f" {conjunction} " + items[-1]


def _kpa(value):
    return f"{three_figures_or_tiny(value)} kPa"


def _celsius(value):
    return f"{three_figures_or_tiny(value)} C"


def _dose(sv_yr):
    return f"{three_figures_or_tiny(sv_yr)} Sv a year" if sv_yr >= 1.0 else f"{three_figures_or_tiny(sv_yr * 1000.0)} mSv a year"


def pressure_need(world):
    """`(rank, gear, why)` for the pressure: the gear class the ambient pressure needs, or `None` gear."""
    p = world.pressure_kpa
    why = f"the air pressure is {_kpa(p)}"
    if p < ARMSTRONG_KPA:
        return "suit", "a full gas-pressurised suit", why
    if p < MASK_MIN_KPA:
        return "suit", "a counterpressure or pressurised suit", why
    if p < limit("pressure_no_gear_kpa"):
        return "mask", "a breathing mask", why
    if p <= limit("pressure_narcosis_kpa"):
        return "none", None, why
    if p <= limit("pressure_hyperbaric_max_kpa"):
        return "suit", "a hyperbaric pressure suit", why
    return "suit", "a rigid armoured pressure suit", why


def temperature_need(world):
    """`(rank, gear, why)` for the temperature, with the wet bulb when it matters."""
    t = world.temperature_c
    why = f"the air temperature is {_celsius(t)}"
    if t < limit("cold_suit_c"):
        return "suit", "an environment suit with active heating", why
    if t < limit("cold_clothing_c"):
        return "garment", "an insulated suit", why
    if t < limit("comfort_low_c"):
        return "garment", "warm clothing", why
    if t > limit("heat_shielded_c"):
        return "suit", "a heat-shielded suit with a liquid cooling garment", why
    if t > limit("heat_garment_c"):
        return "suit", "an actively cooled environment suit", why
    bulb = wet_bulb_c(t, world.relative_humidity) if -20.0 <= t <= 50.0 else None
    if t > limit("comfort_high_c") or (bulb is not None and bulb >= limit("wet_bulb_c")):
        if bulb is not None and bulb >= limit("wet_bulb_c"):
            why = f"{why} and the wet bulb reaches {_celsius(bulb)}"
        return "garment", "cooling garments", why
    return "none", None, why


def air_need(world):
    """`(rank, gear, why)` for the air: the worst of its oxygen and the gases a filter can take out."""
    if not habitability.has_atmosphere(world):
        return "suit", "an environment suit with its own air", f"there is no atmosphere ({_kpa(world.pressure_kpa)})"
    need = ("none", None, "the oxygen and the other gases are within their long-term limits")
    o2 = inspired_o2_kpa(world)
    if o2 < limit("inspired_o2_min_kpa"):
        reaches = "no oxygen reaches the lungs" if o2 < 0.05 else f"only {_kpa(o2)} of oxygen reaches the lungs"
        need = ("mask", "a breathing mask with supplied oxygen", reaches)
    elif o2 > limit("inspired_o2_chronic_kpa"):
        need = ("suit", "an environment suit with its own air", f"{_kpa(o2)} of oxygen reaches the lungs, more than months of breathing allow")
    too_much = [(gas, kpa) for gas, kpa, band in gas_inputs(world) if habitability.band_score(kpa, band, log=True) < 1.0]
    if too_much:
        names = " and ".join(_gas(gas) for gas, _kpa_value in too_much)
        worst_score = min(habitability.band_score(kpa, band, log=True) for _g, kpa, band in gas_inputs(world))
        if worst_score < 0.4:
            return "suit", "an environment suit with its own air", f"{names} is far past what a filter can take out"
        if need[0] in ("none", "mask"):
            return "scrubber", "a mask with a scrubber", f"{names} is above its long-term limit"
    return need


def water_need(world):
    """`(rank, gear, why)` for the water: whether it can be drunk as it is."""
    if world.solvent != "water":
        return "none", "it has to be brought or made from ice", "there is no liquid water to support life"
    inputs = water_inputs(world)
    if not inputs:
        return "none", None, "the water is drinkable as it is"
    score, label = min((habitability.band_score(value, band), name) for name, value, band in inputs)
    if score >= 1.0:
        return "none", None, "the water is drinkable as it is"
    if score < 0.4:
        return "none", "it is undrinkable without heavy treatment", f"its {label} is far outside the safe range"
    return "none", "it needs treatment before drinking", f"its {label} is outside the safe range"


def radiation_need(world):
    """`(rank, gear, why)` for the dose: a person is shielded by shelter, not by a suit."""
    why = f"the radiation dose is a mean of {_dose(world.dose_sv_yr)}"
    uv = habitability.uv_score(world)
    if uv is not None and uv < habitability.band_score(world.dose_sv_yr, habitability.HUMAN.dose, log=True):
        why = f"the ground ultraviolet is {three_figures_or_tiny(world.uv_index)} times Earth's DNA-damaging dose"
        colour = tier(uv)
        edge = {"Blue": "uv_cover_index", "Green": "uv_suit_index", "Yellow": "uv_lethal_index"}.get(colour)
        band = f"inside the {colour} band (up to {three_figures_or_tiny(limit(edge))} times Earth's)" if edge else (
            "past the Yellow band")
        if colour in ("Blue", "Green"):
            return "none", None, f"{why}, {band}"
        if colour == "Yellow":
            return "garment", "a shielded shelter for long stays", f"{why}, {band}"
        return "hardened", "a hardened, shielded habitat", f"{why}, past the Yellow band (up to {three_figures_or_tiny(limit('uv_lethal_index'))} times Earth's)"
    colour = tier(habitability.radiation_score(world))
    if colour == "Blue":
        return "none", None, f"{why}, inside the Blue band (up to {_dose(limit('dose_blue_sv_yr'))})"
    if colour == "Green":
        return "none", None, f"{why}, inside the Green band (up to {_dose(limit('dose_green_sv_yr'))})"
    if colour == "Yellow":
        return "garment", "a shielded shelter for long stays", (
            f"{why}, inside the Yellow band (up to {_dose(limit('dose_yellow_sv_yr'))})")
    return "hardened", "a hardened, shielded habitat", (
        f"{why}, past the Yellow band (up to {_dose(limit('dose_yellow_sv_yr'))}), and no suit shields from it")


def requirements(world):
    """`{category: (rank, gear, why)}` for pressure, temperature, air, water and radiation."""
    return {"Pressure": pressure_need(world), "Temperature": temperature_need(world), "Air": air_need(world),
            "Water": water_need(world), "Radiation": radiation_need(world)}


def gear_clause(levels):
    """The survival equipment as the two lines that follow a rating's headline, from its stored levels
    (`equipment.compute`): "Required: a + b. Recommended: c."; a line with nothing in it is left off."""
    return " ".join(line for line in equipment.lines(levels) if line)


_TEMPERATURE_RANK = {"none": 0, "garment": 1, "suit": 2}

_BAND_PLACES = ("at the north pole", "in the northern middle latitudes", "in the northern tropics", "at the equator",
                "in the southern tropics", "in the southern middle latitudes", "at the south pole")
_RING_PLACES = ("at the centre of the day side", "on the day side", "at the edge of the day side", "at the terminator",
                "at the edge of the night side", "on the night side", "at the centre of the night side")


def _band_means(regions):
    """`[(place, kelvin)]` of the bands a body stores, north to south (day to night for a world that turns once an orbit), and whether it is locked."""
    locked = regions.get(regional_temperature.column(1, "state")) in latitude_bands.RING_STATES
    places = _RING_PLACES if locked else _BAND_PLACES
    means = [(places[number - 1], regions.get(regional_temperature.column(number, "mean_k"))) for number in range(1, regional_temperature.BAND_COUNT + 1)]
    return [(place, kelvin) for place, kelvin in means if kelvin], locked


def regional_sentence(world, regions):
    """
    The regional sentence of the Temperature row (GEN.656, GEN.1108): where the world is colder or hotter than its mean, as
    "The equator averages 12 C and the poles -22 C", the day and night swing when it is 5 K or more, and, for a place
    whose temperature calls for heavier gear than the mean does, "At the north pole, where it is -70 C, a heated sealed
    suit is recommended." `regions` maps the `regional_temperature.FIELDS` to their values (a stored row or a body); `""` when
    there are none.
    """
    if not regions or not any(regions.get(field) for field in regional_temperature.FIELDS):
        return ""
    places, locked = _band_means(regions)
    if not places:
        return ""
    celsius = {place: kelvin - 273.15 for place, kelvin in places}
    if locked:
        text = (f"The day side averages {_celsius(celsius[_RING_PLACES[1]])}, the terminator {_celsius(celsius[_RING_PLACES[3]])} "
                f"and the night side {_celsius(celsius[_RING_PLACES[5]])}.")
    else:
        text = (f"The equator averages {_celsius(celsius['at the equator'])} and the poles "
                f"{_celsius(min(celsius['at the north pole'], celsius['at the south pole']))}.")
    swing = regions.get(regional_temperature.column(4, "swing_k"))
    if not locked and swing and swing >= 5.0:
        text += (f" Day and night differ by about {three_figures_or_tiny(swing)} K from the coldest hour to the hottest "
                 f"(\u00b1{three_figures_or_tiny(swing / 2.0)} K)")
        span = regional_temperature.daily_range(regions.get(regional_temperature.column(4, "mean_k")), swing)
        if span:
            text += (f", so the equator runs from {_celsius(span[0] - 273.15)} to {_celsius(span[1] - 273.15)} "
                     f"over a day")
        text += "."
    mean_rank = _TEMPERATURE_RANK.get(temperature_need(world)[0], 0)
    worst = None
    for place, _kelvin in places:
        rank, gear, _ = temperature_need(replace(world, temperature_c=celsius[place], flux_variation=0.0))
        if gear and _TEMPERATURE_RANK.get(rank, 0) > mean_rank and (worst is None or _TEMPERATURE_RANK[rank] > worst[0]):
            worst = (_TEMPERATURE_RANK[rank], place, celsius[place], gear)
    if worst:
        text += f" {worst[1][0].upper() + worst[1][1:]}, where it is {_celsius(worst[2])}, {worst[3]} is recommended."
    return text


def rating(world, compact_host=None, regions=None):
    """
    The habitability rating of a scored rocky body as data (GEN.208, GEN.216, UX.104): `overall` (the
    PHI-6 value), `gear` (the survival equipment a visitor needs, as `gear_clause`) and `rows`, one per category in `DOMAINS`
    order: `domain`, `tier`, `score` and `text`, the one short sentence of the number and its
    consequence. `regions` (the stored regional temperatures, `regional_temperature.FIELDS`) adds the regional sentence to
    the Temperature row.
    """
    scores = habitability.phi4(world)
    needs = requirements(world)
    levels = equipment.compute(world)
    listed = equipment.items(levels)
    colours = {name: scores[name][1] for name in habitability.DOMAINS}
    rows = []
    for name in habitability.DOMAINS:
        value, colour = scores[name]
        sentence = _suitability_clause(world) if name == "Suitability" else _clause(
            name, world, needs[name], compact_host, [item for item in listed if name in item.rows(colours)])
        if name == "Temperature" and regions:
            sentence = " ".join(part for part in (sentence, regional_sentence(world, regions)) if part)
        rows.append({"domain": name, "tier": colour, "score": value, "text": sentence})
    return {"overall": scores["PHI-4"], "gear": gear_clause(levels), "rows": rows}


def _listed(items):
    """The `Item`s of one category as the two lines that follow its number: "Required: a + b. Recommended: c."; a line with
    nothing in it is left off (never "none")."""
    required = [item.name for item in items if item.severity == equipment.REQUIRED]
    recommended = [item.name for item in items if item.severity == equipment.RECOMMENDED]
    return " ".join(line for line in (("Required: " + " + ".join(required) + ".") if required else "",
                                      ("Recommended: " + " + ".join(recommended) + ".") if recommended else "") if line)


def _clause(name, world, need, compact_host, items):
    """
    One category's sentence: the number, then the gear it asks for. `items` are the `equipment.items` whose PHI-6 domain
    is this category, so the row says what the headline's equipment line says (a Green row may ask for a mask). Water has
    no equipment of its own: it says whether the water can be drunk.
    """
    rank, gear, why = need
    text = why[0].upper() + why[1:]
    listed = _listed(items)
    if name == "Water":
        text += f", so {gear}." if gear and world.solvent == "water" else "."
    elif listed:
        text += f". {listed}"
    elif name in ("Pressure", "Temperature") and rank == "none":
        text += ", so no gear is needed."
    else:
        text += "."
    if name == "Radiation" and compact_host:
        text += f" A planet of a {compact_host} is rated as a lethal dose."
    return text


def _suitability_clause(world):
    """
    Suitability's sentence: which of food, water and shelter are readily available. The water rung is how reachable the
    water is; whether it can be drunk is the Water row's, so reachable water that row asks to be treated is not called
    plainly available: the sentence says it is reachable and repeats the treatment (GEN.691).
    """
    rungs = {"water": habitability.water_access_score(world), "food": habitability.food_score(world),
             "shelter": habitability.shelter_score(world)}
    none = [name for name in ("food", "water", "shelter") if tier(rungs[name]) == "Red"]
    limited = [name for name in ("food", "water", "shelter") if tier(rungs[name]) == "Yellow"]
    good = [name for name in ("food", "water", "shelter") if tier(rungs[name]) in ("Blue", "Green")]
    treatment = water_need(world)[1] if "water" in good and world.solvent == "water" else None
    if treatment:
        good.remove("water")
        reachable = f"the water is reachable but {treatment}"
    if not none and not limited:
        if not treatment:
            return "Water, food and shelter are all readily available."
        text = f"Food and shelter are readily available, and {reachable}."
    else:
        if none and limited:
            text = f"There is no {_join(none, 'or')} readily available, and only limited {_join(limited)}"
        elif none:
            text = f"There is no {_join(none, 'or')} readily available"
        else:
            text = f"There is only limited {_join(limited)} readily available"
        if good:
            text += f", but there is {_join(good)} available" + (f", and {reachable}." if treatment else ".")
        else:
            text += f", but {reachable}." if treatment else "."
    ceiling = habitability.suitability_ceiling(habitability.core_scores(world))
    held = tier(habitability.suitability_score(world)) == ceiling != "Blue" and tier(min(rungs.values())) != ceiling
    if held:
        text += f" The other categories hold it to {ceiling}."
    return text
