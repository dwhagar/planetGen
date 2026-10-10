# planetgen/physics/habitability_explain.py

"""
The words behind a body's PHI-4 colours (UX.91): for each of the four
factors (pressure, temperature, chemistry, radiation) its colour and score,
the inputs that make it up with this world's value for each, and which input
sets the colour. It reads the same `habitability.World` the stored numbers
come from (`habitability_world.world_from`) and the same band edges
(`habitability.PRESSURE_BAND` and the others), so it can never disagree with
the chip on the system list; the wording matches the Habitability levels page
(`web/lib/habitability_ref.py`, UX.90).
"""

from planetgen.physics import habitability
from planetgen.physics.habitability import (
    ARMSTRONG_KPA, DOSE_BAND, GAS_LIMITS_KPA, INSPIRED_O2_MIN_KPA, O2_HIGH_BAND, PRESSURE_BAND, TEMPERATURE_BAND,
    WET_BULB_CRIT_C, band_score, gas_inputs, inspired_o2_kpa, tier, water_inputs, wet_bulb_c)

DOMAIN_ORDER = ("Pressure", "Temperature", "Chemistry", "Radiation")
"""tuple: The four factors, in the order the description lists them."""

GAS_NAMES = {"CO2": "carbon dioxide", "CO": "carbon monoxide", "H2S": "hydrogen sulfide", "SO2": "sulfur dioxide"}
"""dict: What each filtered gas is called in a sentence."""

KIT_TEXT = {
    "ideal": "no equipment (conditions are ideal)",
    "breathing mask": "a breathing mask",
    "mask with scrubber": "a mask with a scrubber",
    "sealed suit": "a sealed suit",
    "full life support with radiation hardening": "full life support with radiation hardening",
}
"""dict: The equipment levels (`habitability.equipment`) as a visitor's needs in a sentence."""


def _num(value):
    """Three significant figures, whole numbers with commas from 1,000 up (never scientific notation)."""
    if abs(value) >= 1000.0:
        return f"{value:,.0f}"
    return f"{value:.3g}" if abs(value) >= 0.001 or value == 0 else f"{value:.2g}".replace("e-0", "e-").replace("e-", " x 10^-")


def _kpa(value):
    return f"{_num(value)} kPa"


def _celsius(value):
    return f"{_num(value)} C"


def _sv(value):
    """A dose in Sv a year, in mSv below 1 Sv."""
    return f"{_num(value)} Sv a year" if abs(value) >= 1.0 else f"{_num(value * 1000.0)} mSv a year"


def _plain(value):
    return _num(value)


def _span(low, high, fmt):
    """`low` to `high` with an open end (`None`) as "up to" or "or more"."""
    if low is None and high is None:
        return "any value"
    if low is None:
        return f"up to {fmt(high)}"
    if high is None:
        return f"{fmt(low)} or more"
    first, last = fmt(low), fmt(high)
    unit = last.partition(" ")[2]
    if unit and first.endswith(" " + unit):
        first = first[:-len(unit) - 1]
    return f"{first} to {last}"


def band_phrase(colour, band, fmt):
    """Where a value of colour `colour` sits against `band` (see the band
    constants in `habitability`), as the page words it: the Blue range, then
    the wider Green and Yellow ranges a value of that colour or worse falls
    outside of. A band with no Green range of its own (the same edges as
    Blue) skips it."""
    blue = (band[3], band[4])
    green = (band[2], band[5])
    yellow = (band[1], band[6])
    if colour == "Blue":
        return f"inside the Blue band ({_span(*blue, fmt)})"
    parts = [f"outside the Blue band ({_span(*blue, fmt)})"]
    if colour == "Green":
        parts.append(f"inside the Green band ({_span(*green, fmt)})")
        return ", ".join(parts)
    if green != blue:
        parts.append(f"outside the Green band ({_span(*green, fmt)})")
    if colour == "Yellow":
        parts.append(f"inside the Yellow band ({_span(*yellow, fmt)})")
    else:
        parts.append(f"beyond the Yellow band ({_span(*yellow, fmt)})")
    return ", ".join(parts)


def _item(label, value_text, score, phrase):
    return {"label": label, "value": value_text, "score": score, "tier": tier(score), "phrase": phrase}


def _banded(label, value, band, fmt, log=False):
    score = band_score(value, band, log=log)
    return _item(label, fmt(value), score, band_phrase(tier(score), band, fmt))


def _oxygen_item(world):
    o2 = inspired_o2_kpa(world)
    high = band_score(o2, O2_HIGH_BAND, log=True)
    if o2 < INSPIRED_O2_MIN_KPA:
        return _item("oxygen", f"{_kpa(o2)} reaching the lungs", min(0.75, high),
                     f"under the {_kpa(INSPIRED_O2_MIN_KPA)} a person needs, so a breathing mask is required "
                     "(Green: a mask makes up the difference where one works)")
    if high >= 1.0:
        return _item("oxygen", f"{_kpa(o2)} reaching the lungs", 1.0,
                     f"inside the Blue band ({_kpa(INSPIRED_O2_MIN_KPA)} to {_kpa(O2_HIGH_BAND[4])})")
    return _item("oxygen", f"{_kpa(o2)} reaching the lungs", high,
                 f"above the {_kpa(O2_HIGH_BAND[4])} a person can breathe for months, "
                 f"{'inside' if tier(high) == 'Yellow' else 'beyond'} the Yellow band "
                 f"(up to {_kpa(O2_HIGH_BAND[6])})")


def _pressure_items(world):
    item = _banded("air pressure", world.pressure_kpa, PRESSURE_BAND, _kpa, log=True)
    if world.pressure_kpa < PRESSURE_BAND[2]:
        item["phrase"] += (f"; under {_kpa(PRESSURE_BAND[2])} even a pure-oxygen mask cannot give a person enough "
                           "oxygen, so a pressure suit is needed")
    elif world.pressure_kpa < ARMSTRONG_KPA:
        item["phrase"] += "; water boils at body temperature here"
    return [item]


def _temperature_items(world):
    items = [_banded("air temperature", world.temperature_c, TEMPERATURE_BAND, _celsius)]
    if -20.0 <= world.temperature_c <= 50.0:
        bulb = wet_bulb_c(world.temperature_c, world.relative_humidity)
        if bulb >= WET_BULB_CRIT_C:
            items.append(_item("wet-bulb temperature (heat and humidity together)", _celsius(bulb), 0.75 - 1e-9,
                               f"at or above {_celsius(WET_BULB_CRIT_C)}, where a body can no longer shed heat, "
                               "so the colour is Yellow at best"))
        else:
            items.append(_item("wet-bulb temperature (heat and humidity together)", _celsius(bulb), 1.0,
                               f"under {_celsius(WET_BULB_CRIT_C)}, where a body can still shed heat"))
    return items


def _chemistry_items(world):
    items = [_oxygen_item(world)]
    for gas, kpa, band in gas_inputs(world):
        item = _banded(GAS_NAMES[gas], kpa, band, _kpa, log=True)
        chronic = GAS_LIMITS_KPA[gas][0]
        if item["tier"] == "Blue":
            item["phrase"] = f"under its long-term limit of {_kpa(chronic)}"
        items.append(item)
    for label, value, band in water_inputs(world):
        name = {"pH": "ocean pH", "water activity": "water activity", "chaotropicity": "chaotropicity"}[label]
        items.append(_banded(name, value, band, _plain))
    return items


def _radiation_items(world, compact_host=None):
    item = _banded("radiation dose", world.dose_sv_yr, DOSE_BAND, _sv, log=True)
    if compact_host:
        item["phrase"] += f"; a planet of a {compact_host} is rated as a lethal dose"
    return [item]


def factors(world, compact_host=None):
    """
    The four factors of `world` as dicts, in `DOMAIN_ORDER`: `domain`,
    `score`, `tier`, `items` (each `label`, `value`, `score`, `tier`,
    `phrase`) and `limiting`, the item that sets the colour (the lowest
    score, the first of equals).
    """
    scores = habitability.phi4(world)
    builders = {"Pressure": _pressure_items, "Temperature": _temperature_items, "Chemistry": _chemistry_items,
                "Radiation": lambda w: _radiation_items(w, compact_host)}
    out = []
    for domain in DOMAIN_ORDER:
        score, colour = scores[domain]
        items = builders[domain](world)
        limiting = min(items, key=lambda entry: entry["score"])
        out.append({"domain": domain, "score": score, "tier": colour, "items": items, "limiting": limiting})
    return out


def paragraphs(world, bold, compact_host=None):
    """
    The explanation as paragraphs of text for a planet or moon's description
    (UX.91): one introducing the score and the equipment, then one per factor
    naming its colour, the input that sets it with this world's value, and
    the other inputs with their colours. `bold` wraps a phrase in the page
    format's emphasis (Markdown or wikitext).
    """
    world_scores = habitability.scores(world)
    kit = world_scores["equipment"]
    lead = (f"{bold('PHI-4 habitability')} is {world_scores['PHI-4']:.2f}, and a human visitor needs "
            f"{KIT_TEXT[kit]}. "
            "The score is the geometric mean of four colour factors, pressure, temperature, chemistry and "
            "radiation; each factor takes the colour of its worst input, from Blue (best) through Green and "
            "Yellow to Red.")
    out = [lead]
    for factor in factors(world, compact_host):
        limiting = factor["limiting"]
        text = (f"{bold(factor['domain'] + ': ' + factor['tier'])} ({factor['score']:.2f}). "
                f"The {limiting['label']} is {limiting['value']}, {limiting['phrase']}.")
        others = [entry for entry in factor["items"] if entry is not limiting]
        if others:
            text += " Also considered: " + "; ".join(
                f"the {entry['label']}, {entry['value']} ({entry['tier']})" for entry in others) + "."
        out.append(text)
    return out
