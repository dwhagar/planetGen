# planetgen/physics/protection.py

"""
The protective gear a visitor needs on a world, and the paragraph that says where each
habitability category falls and why (GEN.208).

The gear ladder is the research notes' (Environmental.Extremes.and.Protection.Research.md,
Thermal.Respiratory.Protection.Research.md): modular gear (clothing, masks), intermediate
systems (counterpressure and chemical-barrier suits) and fully enclosed suits, one ladder
per stressor, with the limits of `protection_limits` (one table, the more permissive value
wherever the notes and the built table differ). A suit built for one stressor is not one built
for another: a heated sealed suit does not cool, a cooled one does not shield from radiation,
and solar-storm doses need a shelter, not a suit.

Everything here is a pure function of a `habitability.World`.
"""

from planetgen.physics import habitability
from planetgen.physics.habitability import (
    ARMSTRONG_KPA, GAS_LIMITS_KPA, MASK_MIN_KPA, TIERS, gas_inputs, inspired_o2_kpa, tier, water_inputs, wet_bulb_c)
from planetgen.physics.protection_limits import limit

SUIT_RANK = {"none": 0, "mask": 1, "scrubber": 2, "garment": 3, "suit": 4, "hardened": 5}
"""dict: How heavy a class of gear is, lightest first: nothing, a breathing mask, a mask with a
scrubber, clothing or garments for the weather, a sealed suit, a hardened habitat (no suit helps)."""

GAS_NAMES = {"CO2": "carbon dioxide", "CO": "carbon monoxide", "H2S": "hydrogen sulfide", "SO2": "sulfur dioxide"}


def _num(value):
    if abs(value) >= 1000.0:
        return f"{value:,.0f}"
    return f"{value:.3g}" if abs(value) >= 0.001 or value == 0 else f"{value:.1e}"


def _join(items, conjunction="and"):
    """`a`, `a and b`, `a, b and c`."""
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + f" {conjunction} " + items[-1]


def _kpa(value):
    return f"{_num(value)} kPa"


def _celsius(value):
    return f"{_num(value)} C"


def _dose(sv_yr):
    return f"{_num(sv_yr)} Sv a year" if sv_yr >= 1.0 else f"{_num(sv_yr * 1000.0)} mSv a year"


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
        return "suit", "a sealed suit with active heating", why
    if t < limit("cold_clothing_c"):
        return "garment", "an insulated suit", why
    if t < limit("comfort_low_c"):
        return "garment", "warm clothing", why
    if t > limit("heat_shielded_c"):
        return "suit", "a heat-shielded suit with a liquid cooling garment", why
    if t > limit("heat_garment_c"):
        return "suit", "an actively cooled sealed suit", why
    bulb = wet_bulb_c(t, world.relative_humidity) if -20.0 <= t <= 50.0 else None
    if t > limit("comfort_high_c") or (bulb is not None and bulb >= limit("wet_bulb_c")):
        if bulb is not None and bulb >= limit("wet_bulb_c"):
            why = f"{why} and the wet bulb reaches {_celsius(bulb)}"
        return "garment", "cooling garments", why
    return "none", None, why


def air_need(world):
    """`(rank, gear, why)` for the air: the worst of its oxygen and the gases a filter can take out."""
    if not habitability.has_atmosphere(world):
        return "suit", "a sealed suit with its own air", f"there is no atmosphere ({world.pressure_kpa:.1f} kPa)"
    need = ("none", None, "the oxygen and the other gases are within their long-term limits")
    o2 = inspired_o2_kpa(world)
    if o2 < limit("inspired_o2_min_kpa"):
        reaches = "no oxygen reaches the lungs" if o2 < 0.05 else f"only {_kpa(o2)} of oxygen reaches the lungs"
        need = ("mask", "a breathing mask with supplied oxygen", reaches)
    elif o2 > limit("inspired_o2_chronic_kpa"):
        need = ("suit", "a sealed suit with its own air", f"{_kpa(o2)} of oxygen reaches the lungs, more than months of breathing allow")
    too_much = [(gas, kpa) for gas, kpa, band in gas_inputs(world) if habitability.band_score(kpa, band, log=True) < 1.0]
    if too_much:
        names = " and ".join(GAS_NAMES[gas] for gas, _kpa_value in too_much)
        worst_score = min(habitability.band_score(kpa, band, log=True) for _g, kpa, band in gas_inputs(world))
        if worst_score < 0.4:
            return "suit", "a sealed suit with its own air", f"{names} is far past what a filter can take out"
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


SHORT_GEAR = {
    "a full gas-pressurised suit": "a pressure suit", "a counterpressure or pressurised suit": "a pressure suit",
    "a hyperbaric pressure suit": "a hyperbaric suit", "a rigid armoured pressure suit": "an armoured suit",
    "a sealed suit with active heating": "an environment suit", "an actively cooled sealed suit": "an environment suit",
    "a heat-shielded suit with a liquid cooling garment": "an environment suit",
    "a breathing mask with supplied oxygen": "a breathing mask", "a sealed suit with its own air": "a sealed suit",
    "a shielded shelter for long stays": "a shielded shelter", "a hardened, shielded habitat": "a hardened habitat",
}
"""dict: The gear in the short words a description uses ("so a pressure suit is needed")."""

OVERALL_GEAR = {"none": "no equipment", "mask": "a breathing mask", "scrubber": "a mask with a scrubber",
                "suit": "a sealed suit", "hardened": "a hardened habitat"}
"""dict: What a visitor needs overall by the heaviest class of gear any stressor asks for (the lighter
classes name their own gear)."""


def _short(gear):
    return SHORT_GEAR.get(gear, gear)


def kit_sentence(world, needs=None):
    """The gear a visitor needs overall, as a phrase: the heaviest class of gear any stressor asks for."""
    needs = needs or requirements(world)
    body = [need for name, need in needs.items() if name in ("Pressure", "Temperature", "Air", "Radiation")]
    top = max(SUIT_RANK[need[0]] for need in body)
    rank = next(name for name, value in SUIT_RANK.items() if value == top)
    if rank in OVERALL_GEAR:
        return OVERALL_GEAR[rank]
    return _join(dict.fromkeys(_short(need[1]) for need in body if SUIT_RANK[need[0]] == top and need[1]))


def paragraph(world, bold, compact_host=None):
    """
    The habitability description of a scored rocky body as one paragraph (GEN.208, GEN.216): the PHI-6
    value and what a visitor needs overall, then each of the six categories with its colour and value,
    one short sentence of the number and its consequence. `bold` wraps a phrase in the page format's
    emphasis. Nothing is said beyond the categories and the gear.
    """
    scores = habitability.phi4(world)
    needs = requirements(world)
    parts = [f"{bold('PHI-6 habitability')} is {scores['PHI-4']:.2f}, and a visitor needs {kit_sentence(world, needs)}."]
    for name in habitability.DOMAINS:
        value, colour = scores[name]
        sentence = _suitability_clause(world) if name == "Suitability" else _clause(name, world, needs[name], compact_host)
        parts.append(f"{bold(f'{name}: {colour}')} ({value:.2f}). {sentence}")
    return " ".join(parts)


def _clause(name, world, need, compact_host):
    """One category's sentence: the number, then the gear it asks for."""
    rank, gear, why = need
    text = why[0].upper() + why[1:]
    if name == "Air" and not habitability.has_atmosphere(world):
        gear = None
    elif name == "Water" and world.solvent != "water":
        gear = None
    elif name == "Radiation" and rank != "hardened":
        gear = None
    if gear:
        text += f", so {_short(gear)} is needed" if name != "Water" else f", so {gear}"
    elif name in ("Pressure", "Temperature") and rank == "none":
        text += ", so no gear is needed"
    text += "."
    if name == "Radiation" and compact_host:
        text += f" A planet of a {compact_host} is rated as a lethal dose."
    return text


def _suitability_clause(world):
    """Suitability's sentence: which of food, water and shelter are readily available."""
    rungs = {"water": habitability.water_access_score(world), "food": habitability.food_score(world),
             "shelter": habitability.shelter_score(world)}
    none = [name for name in ("food", "water", "shelter") if tier(rungs[name]) == "Red"]
    limited = [name for name in ("food", "water", "shelter") if tier(rungs[name]) == "Yellow"]
    good = [name for name in ("food", "water", "shelter") if tier(rungs[name]) in ("Blue", "Green")]
    if not none and not limited:
        return "Water, food and shelter are all readily available."
    if none and limited:
        text = f"There is no {_join(none, 'or')} readily available, and only limited {_join(limited)}"
    elif none:
        text = f"There is no {_join(none, 'or')} readily available"
    else:
        text = f"There is only limited {_join(limited)} readily available"
    text += f", but there is {_join(good)} available." if good else "."
    ceiling = habitability.suitability_ceiling(habitability.core_scores(world))
    held = tier(habitability.suitability_score(world)) == ceiling != "Blue" and tier(min(rungs.values())) != ceiling
    if held:
        text += f" The other categories hold it to {ceiling}."
    return text
