# planetgen/web/lib/habitability_ref.py

"""
The text of the habitability levels page (`/classes/habitability`, UX.90):
the five equipment levels, the six PHI-6 domains with their Blue, Green,
Yellow and Red thresholds, the compact-host rule and the chips a planet row
shows. The equipment labels come from `planetgen.physics.habitability_world`;
the thresholds are read from `planetgen.physics.habitability` and `protection_limits`, where they are named with
their sources (`domains`), so the tables cannot disagree with the score.
"""

from planetgen.physics import atmosphere, habitability
from planetgen.physics.habitability_world import EQUIPMENT_LABELS, LETHAL_DOSE_SV_YR
from planetgen.physics.protection_limits import limit
from planetgen.web.lib import habitability_symbols

EQUIPMENT_MEANINGS = (
    "Ordinary weather clothing is enough: the air has the right pressure and oxygen, no harmful gas, "
    "a comfortable temperature and a low radiation dose.",
    "A breathing mask: oxygen is low or the pressure is a little off, but a mask makes the air safe.",
    "A mask with a scrubber: a gas such as carbon dioxide, carbon monoxide, hydrogen sulfide or sulfur dioxide "
    "is above its long-term limit and has to be filtered out.",
    "An environment suit: the pressure is too low for a mask to work, or the temperature or the air's chemistry "
    "is outside what a person can stand.",
    "Full life support with radiation hardening: the radiation dose is beyond what any suit can shield, "
    "so a visitor needs a hardened habitat.",
)
"""tuple: What each equipment level means, in the order of `EQUIPMENT_LABELS`."""

DOMAIN_TIERS = (
    ("Blue", "The best band: a visitor notices nothing wrong."),
    ("Green", "Workable with ordinary care."),
    ("Yellow", "Needs equipment or limits the time spent outside."),
    ("Red", "Beyond what a person can survive without heavy protection."),
)
"""tuple: `(colour, meaning)` of a PHI-6 domain's score, best first."""

def _n(value):
    """A limit as the page prints it: whole numbers with thousands separators, others as they are."""
    return f"{int(value):,}" if float(value).is_integer() else f"{value:g}"


def _between(low, high, unit):
    return f"{_n(low)} to {_n(high)} {unit}".strip()


def _up_to(value, unit):
    return f"up to {_n(value)} {unit}"


def _under(value, unit):
    return f"under {_n(value)} {unit}"


def _gas_rows(gases):
    """One row per gas a filter handles: Blue under the long-term limit, Yellow up to ten times the short-term one."""
    rows = []
    for gas, _kpa, band in gases:
        name = atmosphere.SPECIES_NAMES[gas.lower()].capitalize()
        rows.append((f"{name} (a filter takes it out)", _under(band[4], "kPa"), "", _up_to(band[6], "kPa")))
    return rows


def domains(bands=habitability.HUMAN):
    """
    `(domain, [(variable, blue, green, yellow)])`, beyond the Yellow column a variable is Red. Every number is
    read from `bands` (the limits in `protection_limits`) and the constants of `habitability`, so the page cannot
    drift from the score.
    """
    pressure, temperature, dose = bands.pressure, bands.temperature, bands.dose
    gases = habitability.gas_inputs(habitability.World(pressure_kpa=0.0, temperature_c=0.0, gases={}))
    co2, filtered = gases[0], gases[1:]
    water = {label: band for label, _value, band in habitability.water_inputs(
        habitability.World(pressure_kpa=0.0, temperature_c=0.0, gases={}, ph=7.0))}
    ph, activity, chaotropicity = water["pH"], water["water activity"], water["chaotropicity"]
    return (
        ("Pressure", [
            ("Air pressure", _between(pressure[3], pressure[4], "kPa"), _between(pressure[2], pressure[5], "kPa"),
             _between(pressure[1], pressure[6], "kPa")),
        ]),
        ("Temperature", [
            ("Air temperature", _between(temperature[3], temperature[4], "C"),
             _between(temperature[2], temperature[5], "C"), _between(temperature[1], temperature[6], "C")),
            ("Wet-bulb temperature (heat and humidity together)", _under(bands.wet_bulb, "C"),
             _under(bands.wet_bulb, "C"), f"{_n(bands.wet_bulb)} C or more (a body can no longer shed heat)"),
        ]),
        ("Air", [
            (f"Oxygen the lungs get (the pressure minus {_n(bands.armstrong)} kPa of water vapour), too little",
             f"{_n(bands.o2_min)} kPa or more", f"{_under(bands.o2_min, 'kPa')} (a mask works)", ""),
            ("Oxygen the lungs get, too much", _up_to(bands.o2_high[4], "kPa"), "",
             _up_to(bands.o2_high[6], "kPa")),
            ("Carbon dioxide", _under(co2[2][4], "kPa"), _up_to(co2[2][5], "kPa"), _up_to(co2[2][6], "kPa")),
            *_gas_rows(filtered),
        ]),
        ("Water", [
            ("Liquid water to drink", "present", "present", "present (none at all is Red)"),
            ("Water acidity (pH)", _between(ph[3], ph[4], ""), _between(ph[2], ph[5], ""),
             _between(ph[1], ph[6], "")),
            ("Water activity (how freely water is available)", f"{_n(activity[3])} or more",
             f"{_n(activity[2])} or more", f"{_n(activity[1])} or more"),
            ("Chaotropicity (how strongly dissolved salts pull the water's structure apart)",
             _under(chaotropicity[4], "kJ/kg"), "", ""),
        ]),
        ("Radiation", [
            ("Dose per year", _under(dose[4], "Sv"), _up_to(dose[5], "Sv"), _up_to(dose[6], "Sv")),
        ]),
        ("Suitability", [
            ("Water within reach (the lowest of the three rungs scores)", "fresh liquid water",
             "rivers, lakes and rain", "ice or an ocean under ice"),
            ("Food", "forests and animals", "algae, mosses and the first land plants",
             "microbial mats and plankton only"),
            ("Shelter", "the surface dose is Blue", "the surface dose is Green",
             "the surface dose is Yellow (caves or burial)"),
        ]),
    )


SCORE_RULE = (
    "Each variable scores 1 in Blue and falls to 0.75 at the edge of Green, 0.4 at the edge of Yellow and 0 "
    "at the outer edge of Red. A domain takes the score of its worst variable: Blue at 1, Green from 0.75, "
    "Yellow from 0.4, Red below. The PHI-6 score on a "
    "planet's tooltip is the geometric mean of the six domains, held under the colour gates: Radiation and "
    "Temperature cap it at their own colour, Pressure in Yellow or Red caps it there (unless a full suit is "
    "needed anyway), Air, Water or a Red Suitability cap it at Yellow, and it is Blue only if all six are. "
    "Suitability, the long-term domain, is the lowest of water, food and shelter, and is never better than "
    "the colour the other five allow."
)

EQUIPMENT_RULE = (
    "The equipment level comes from the worst of pressure, temperature, air and radiation: radiation in the Red band needs full life support; "
    "pressure below Green, or any domain in the Red band, needs an environment suit; a gas above its long-term limit "
    "needs a mask with a scrubber; low oxygen or a pressure below Blue needs a mask; anything else is ideal. "
    f"Air hotter than {limit('heat_shielded_c'):g} C burns a mask-only visitor's skin at the seal, so it needs a "
    "heat-shielded suit with a liquid cooling garment."
)

COMPACT_HOST_RULE = (
    f"A planet of a pulsar, a neutron star or a black hole is rated as if its surface took "
    f"{LETHAL_DOSE_SV_YR:,.0f} Sv a year. That puts radiation in Red, so the planet needs full life support "
    "and its life scores are 0. Its note names the host."
)

CHIPS = (
    ("Habitable class", "The planet's class can carry life (an Earth-like or similar world). It is not a rating "
                        "and not a promise that anything lives there: the equipment level beside it is the rating."),
    ("Habitable-class moon", "A moon of this planet is of a habitable class in that sense."),
    ("Inhabited", "A species lives here."),
)
"""tuple: The chips a planet or moon row can show beside its equipment level."""


SYMBOL_INTRO = (
    "Lists show a rating as a small shape drawn in its colour. The shape says which domain it is, and its tooltip "
    "and its screen-reader text give the domain and the colour in words, so the colour is never the only sign."
)

SYMBOL_ROWS = tuple((name, habitability_symbols.LABELS[name], habitability_symbols.SHAPES[name])
                    for name in habitability_symbols.SYMBOL_NAMES)
"""tuple: `(name, words, what the shape looks like)` of each symbol, the six domains and then the overall one."""


def page():
    """Everything the page shows."""
    return {
        "symbol_intro": SYMBOL_INTRO,
        "symbol_rows": SYMBOL_ROWS,
        "tier_numbers": tuple(range(len(habitability_symbols.TIER_CLASSES))),
        "equipment": list(zip(EQUIPMENT_LABELS, EQUIPMENT_MEANINGS)),
        "domain_tiers": DOMAIN_TIERS,
        "domains": domains(),
        "score_rule": SCORE_RULE,
        "equipment_rule": EQUIPMENT_RULE,
        "compact_host_rule": COMPACT_HOST_RULE,
        "chips": CHIPS,
    }
