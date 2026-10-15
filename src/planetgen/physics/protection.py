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


def _join(items):
    """`a`, `a and b`, `a, b and c`."""
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _kpa(value):
    return f"{_num(value)} kPa"


def _celsius(value):
    return f"{_num(value)} C"


def _dose(sv_yr):
    return f"{_num(sv_yr)} Sv a year" if sv_yr >= 1.0 else f"{_num(sv_yr * 1000.0)} mSv a year"


def pressure_need(world):
    """`(rank, gear, why)` for the pressure: the gear class the ambient pressure needs, or `None` gear."""
    p = world.pressure_kpa
    if p < ARMSTRONG_KPA:
        return "suit", "a full gas-pressurised suit", f"at {_kpa(p)} water boils at body temperature"
    if p < MASK_MIN_KPA:
        return "suit", "a counterpressure or pressurised suit", (
            f"at {_kpa(p)} even a pure-oxygen mask cannot give enough oxygen")
    if p < limit("pressure_no_gear_kpa"):
        return "mask", "a breathing mask", f"at {_kpa(p)} the oxygen is thin and a mask makes up the difference"
    if p <= limit("pressure_narcosis_kpa"):
        return "none", None, f"at {_kpa(p)} no gear is needed"
    if p <= limit("pressure_hyperbaric_max_kpa"):
        return "suit", "a hyperbaric pressure suit", f"at {_kpa(p)} the gas is dense enough to narcotise and overload breathing"
    return "suit", "a rigid armoured pressure suit", f"at {_kpa(p)} the crushing pressure and gas density defeat any soft suit"


def temperature_need(world):
    """`(rank, gear, why)` for the temperature, with the wet bulb when it matters."""
    t = world.temperature_c
    if t < limit("cold_suit_c"):
        return "suit", "a sealed suit with active heating", f"at {_celsius(t)} insulation alone cannot hold body heat"
    if t < limit("cold_clothing_c"):
        return "garment", "an insulated suit", f"at {_celsius(t)} cold-weather clothing is not enough"
    if t < limit("comfort_low_c"):
        return "garment", "warm clothing", f"at {_celsius(t)} ordinary warm clothing is needed"
    if t > limit("heat_shielded_c"):
        return "suit", "a heat-shielded suit with a liquid cooling garment", (
            f"at {_celsius(t)} the skin burns and a mask's seal fails without one")
    if t > limit("heat_garment_c"):
        return "suit", "an actively cooled sealed suit", f"at {_celsius(t)} the body cannot shed heat"
    bulb = wet_bulb_c(t, world.relative_humidity) if -20.0 <= t <= 50.0 else None
    if t > limit("comfort_high_c") or (bulb is not None and bulb >= limit("wet_bulb_c")):
        reason = (f"at {_celsius(t)} the wet bulb reaches {_celsius(bulb)}, where sweat no longer cools"
                  if bulb is not None and bulb >= limit("wet_bulb_c") else f"at {_celsius(t)} it is hot")
        return "garment", "cooling garments", reason
    return "none", None, f"at {_celsius(t)} no gear is needed"


def air_need(world):
    """`(rank, gear, why)` for the air: the worst of its oxygen and the gases a filter can take out."""
    need = ("none", None, "the oxygen and the gases are within their long-term limits")
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
            return "scrubber", "a mask with a scrubber", f"{names} is above its long-term limit and must be filtered out"
    return need


def water_need(world):
    """`(rank, gear, why)` for the water: whether it can be drunk as it is."""
    if world.solvent != "water":
        return "none", "it has to be brought or made from ice", "there is no liquid water"
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
    dose = world.dose_sv_yr
    colour = tier(habitability.radiation_score(world))
    if colour == "Blue":
        return "none", None, f"{_dose(dose)} is under the {_dose(limit('dose_blue_sv_yr'))} limit"
    if colour == "Green":
        return "none", None, f"{_dose(dose)} is over the {_dose(limit('dose_blue_sv_yr'))} limit, which a visit allows"
    if colour == "Yellow":
        return "garment", "a shielded shelter for long stays", (
            f"{_dose(dose)} is past the {_dose(limit('dose_green_sv_yr'))} a visit allows")
    return "hardened", "a hardened, shielded habitat", (
        f"{_dose(dose)} is over the {_dose(limit('dose_yellow_sv_yr'))} that causes radiation sickness and no suit shields from it")


def requirements(world):
    """`{category: (rank, gear, why)}` for pressure, temperature, air, water and radiation."""
    return {"Pressure": pressure_need(world), "Temperature": temperature_need(world), "Air": air_need(world),
            "Water": water_need(world), "Radiation": radiation_need(world)}


def kit_sentence(world, needs=None):
    """The one gear sentence: the heaviest class any stressor needs, and a clause for the others that
    make the lighter ones moot (a sealed suit carries its own air, so air in Yellow does not add a mask)."""
    needs = needs or requirements(world)
    body = {name: need for name, need in needs.items() if name in ("Pressure", "Temperature", "Air", "Radiation")}
    top = max(body.values(), key=lambda need: SUIT_RANK[need[0]])
    if SUIT_RANK[top[0]] == 0:
        return "A visitor needs no equipment."
    heavy = [name for name, need in body.items() if SUIT_RANK[need[0]] == SUIT_RANK[top[0]] and need[1]]
    gear = [body[name][1] for name in heavy]
    unique = list(dict.fromkeys(gear))
    what = _join(unique)
    because = _join(name.lower() for name in heavy)
    return f"A visitor needs {what}, because of the {because}."


def paragraph(world, bold, compact_host=None):
    """
    The habitability description of a scored rocky body as one paragraph (GEN.208): its overall score
    and colour, then each of the six categories with its colour, its value and why it falls there and the
    gear it asks for, and last what the visitor needs overall. A category that a heavier suit makes moot
    says so (when the temperature is Red the air is moot, because the suit has its own), and a factor
    that is fine while another makes the surface impossible is named. `bold` wraps a phrase in the page
    format's emphasis. Nothing is said beyond the categories and the gear.
    """
    scores = habitability.phi4(world)
    needs = requirements(world)
    parts = [f"{bold('Habitability')} is {scores['PHI-4']:.2f} ({tier(scores['PHI-4'])})."]
    covering = [name for name in ("Pressure", "Temperature", "Radiation") if SUIT_RANK[needs[name][0]] >= SUIT_RANK["suit"]]
    for name in habitability.DOMAINS:
        colour = scores[name][1]
        clause = _suitability_clause(world) if name == "Suitability" else _clause(name, needs[name], compact_host)
        if name == "Air" and needs["Air"][1] and SUIT_RANK[needs["Air"][0]] < SUIT_RANK["suit"] and covering:
            clause += (f" The {covering[0].lower()} already calls for a sealed suit with its own air, "
                       "so this changes nothing.")
        parts.append(f"{bold(name + ': ' + colour)}. {clause}")
    parts.append(kit_sentence(world, needs))
    impossible = [name for name in ("Pressure", "Temperature", "Air", "Radiation") if scores[name][1] == "Red"]
    fine = [name for name in ("Pressure", "Temperature", "Air", "Water", "Radiation") if scores[name][1] in ("Blue", "Green")]
    if impossible and fine:
        parts.append(f"The {_join(name.lower() for name in impossible)} "
                     f"{'makes' if len(impossible) == 1 else 'make'} the surface impossible without protection, "
                     f"even though the {_join(name.lower() for name in fine)} "
                     f"{'is' if len(fine) == 1 else 'are'} fine.")
    return " ".join(parts)


def _clause(name, need, compact_host):
    """One category's sentence: why it falls where it does and the gear it asks for."""
    _rank, gear, why = need
    text = why[0].upper() + why[1:]
    if gear:
        text += f", so {gear}" if name == "Water" else f", so a visitor needs {gear}"
    text += "."
    if name == "Radiation" and compact_host:
        text += f" A planet of a {compact_host} is rated as a lethal dose."
    return text


def _suitability_clause(world):
    """Suitability's sentence: the weakest of water within reach, food and shelter, and the ceiling the other categories set."""
    rungs = {"water within reach": habitability.water_access_score(world), "food": habitability.food_score(world),
             "shelter": habitability.shelter_score(world)}
    weakest = min(rungs, key=rungs.get)
    ceiling = habitability.suitability_ceiling(habitability.core_scores(world))
    held = tier(habitability.suitability_score(world)) == ceiling != "Blue" and tier(min(rungs.values())) != ceiling
    if rungs[weakest] >= 1.0:
        text = "Water, food and shelter are all within easy reach."
    else:
        stage = world.life_stage if world.life_on_surface else "none"
        detail = {
            "water within reach": {"liquid": "there is liquid water at the surface", "ice": "the water is ice or lies under ice",
                                   "hydrated": "water is only in minerals", "vapour": "water is only vapour",
                                   None: "there is no water"}[world.water_source],
            "food": {"none": "there is no native food on the surface", "microbial": "the native life is microbial only",
                     "macroscopic": "the native life is the first plants and animals",
                     "biosphere": "the native life is a full biosphere"}[stage],
            "shelter": "the radiation dose sets how much shelter is needed",
        }[weakest]
        text = f"The weakest of water, food and shelter is {weakest}: {detail}."
    if held:
        text += f" The other categories hold it to {ceiling}."
    return text
