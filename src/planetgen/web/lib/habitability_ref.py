# planetgen/web/lib/habitability_ref.py

"""
The text of the habitability page (`/help/habitability`, UX.90, UX.180):
the survival equipment items, the six PHI-6 domains with their Blue, Green,
Yellow and Red thresholds, the compact-host rule and the chips a planet row
shows. The equipment items come from `planetgen.physics.equipment`;
the thresholds are read from `planetgen.physics.habitability` and `protection_limits`, where they are named with
their sources (`domains`), so the tables cannot disagree with the score.
"""

from planetgen.physics import atmosphere, constants, habitability, orbits
from planetgen.physics.habitability_world import LETHAL_DOSE_SV_YR
from planetgen.physics.protection_limits import limit
from planetgen.web.lib import habitability_symbols

EQUIPMENT_ITEMS = (
    ("Clothing", "Cold-weather, cooling or sun-cover clothing: the temperature is a little outside comfort or the "
                 "ultraviolet is strong."),
    ("Thermal suit", "An insulated, liquid-cooled or heated diving suit: the temperature is well outside what "
                     "clothing handles, or the pressure needs a heated diving suit."),
    ("Sealed suit", "A fully enclosed suit with its own climate: heated for the deep cold, heat-shielded for the "
                    "hot, chemical-barrier where a corrosive gas or liquid is about."),
    ("Pressure suit", "The ambient pressure is too low for a mask to give enough oxygen."),
    ("Armoured suit", "A rigid suit that holds one atmosphere inside: crushing pressure or heat past 400 C."),
    ("Vehicle only", "No suit exists for it: heat past 500 C, pressure past 100 MPa or gravity past 3 g."),
    ("Powered frame", "The gravity is past 1.5 g (recommended from 1.25 g)."),
    ("Dust filter or seals", "Fine regolith dust wrecks open seals (a rocky body with almost no air)."),
    ("Gas filter, full-face or powered", "A gas such as carbon monoxide, hydrogen sulfide or sulfur dioxide is above "
                                         "its long-term limit and has to be filtered out."),
    ("CO2 scrubber", "Carbon dioxide is above its long-term limit but short of what is acutely dangerous."),
    ("Gas-mix mask", "Oxygen is low or high, or the pressure is dense: an open mask fixes the mix."),
    ("Self-contained air", "A sealed suit or mask on its own air: the ambient gas cannot be breathed (too little "
                           "pressure, too much oxygen or carbon dioxide, a toxic gas past its immediate danger "
                           "limit, or air too cold or hot)."),
    ("Mixed-gas rig and isolated air", "Heliox or trimix for pressure from 4 to 100 bar; air carried inside a rigid "
                                       "suit beyond that."),
    ("Shielded garment, storm vest, hardened suit", "Three steps of radiation shielding worn on the body, from a "
                                                    "thin shielded layer to a hardened suit with storm warning."),
    ("Shielded shelter and hardened habitat", "Thick shielding to stay in: water or regolith walls, or a buried "
                                              "habitat when cosmic rays are the problem."),
)
"""tuple: `(item, what makes it needed)` for the survival equipment line (GEN.655; UX.153 builds the page from the
limits table)."""

EQUIPMENT_ROW_IDS = {
    "Clothing": "equipment-clothing",
    "Thermal suit": "equipment-thermal-suit",
    "Sealed suit": "equipment-sealed-suit",
    "Pressure suit": "equipment-pressure-suit",
    "Armoured suit": "equipment-armoured-suit",
    "Vehicle only": "equipment-vehicle",
    "Powered frame": "equipment-powered-frame",
    "Dust filter or seals": "equipment-dust",
    "Gas filter, full-face or powered": "equipment-gas-filter",
    "CO2 scrubber": "equipment-co2-scrubber",
    "Gas-mix mask": "equipment-gas-mix-mask",
    "Self-contained air": "equipment-self-contained-air",
    "Mixed-gas rig and isolated air": "equipment-mixed-gas",
    "Shielded garment, storm vest, hardened suit": "equipment-radiation-garment",
    "Shielded shelter and hardened habitat": "equipment-shelter",
}
"""dict: The id of each `EQUIPMENT_ITEMS` row on the page, the target of the equipment chips (UX.181). The rows an item
name reaches are `habitability_symbols.EQUIPMENT_ANCHORS`."""

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
    "A planet's equipment line lists only what a visitor needs, never what is not needed. Each item is required "
    "(past the physical limit) or recommended (in the band before it, where a person is impaired or at risk, or "
    "where only a region or the night asks for it); a world that needs nothing shows no line. The line follows the "
    "physics, so it can be stricter than a domain's colour: a gas past its immediate danger limit needs sealed "
    "air even where the score only says Yellow. Radiation goes up in steps: a shielded garment, a storm vest, a "
    "hardened suit, then a shielded shelter and a hardened habitat, each worked out from how much of the dose is "
    f"cosmic rays and how much is storm particles. Air hotter than {limit('heat_shielded_c'):g} C burns a mask-only "
    "visitor's skin at the seal, so it needs a heat-shielded suit."
)

COMPACT_HOST_RULE = (
    f"A planet of a pulsar, a neutron star or a black hole is rated as if its surface took "
    f"{LETHAL_DOSE_SV_YR:,.0f} Sv a year. That puts radiation in Red, so the planet needs a hardened habitat "
    "and its life scores are 0. Its note names the host."
)

CHIPS = (
    ("Habitable class", "The planet's class can carry life (an Earth-like or similar world). It is not a rating "
                        "and not a promise that anything lives there: the equipment level beside it is the rating."),
    ("Habitable-class moon", "A moon of this planet is of a habitable class in that sense."),
    ("Inhabited", "A species lives here."),
)
"""tuple: The chips a planet or moon row can show beside its equipment level."""

CHIP_IDS = {chip: habitability_symbols.CHIP_ANCHORS[chip.lower()] for chip, _meaning in CHIPS}
"""dict: The id of each `CHIPS` entry on the page, the target of the flag chips (UX.181)."""


SYMBOL_INTRO = (
    "Lists show a rating as a small shape drawn in its colour. The shape says which domain it is, and its tooltip "
    "and its screen-reader text give the domain and the colour in words, so the colour is never the only sign."
)

LIFE_INTRO = (
    "Four life scores sit beside the six symbols. Each is a whole number from 0 (none) to 9 (abundant), set by the "
    "surface conditions and the stage life has reached on that world, and is drawn as a shape with the number next to it."
)

LIFE_ROWS = tuple((habitability_symbols.LIFE_NAMES[index], habitability_symbols.LIFE_LABELS[name],
                   habitability_symbols.LIFE_SHAPES[name])
                  for index, name in enumerate(habitability_symbols.LIFE_NAMES))
"""tuple: `(name, words, what the shape looks like)` of each life score."""

SYMBOL_ROWS = tuple((name, habitability_symbols.LABELS[name], habitability_symbols.SHAPES[name])
                    for name in habitability_symbols.SYMBOL_NAMES)
"""tuple: `(name, words, what the shape looks like)` of each symbol, the six domains and then the overall one."""


WHERE_INTRO = (
    "The rating is shown the same way wherever a planet or moon appears: the overall symbol first, then the six "
    "domain symbols. A world with no score yet (a gas giant has no surface, and a body not yet rated) shows none."
)

WHERE_ROWS = (
    ("Planets and Moons lists", "/lists/planet and /lists/moon: a Habitability column before the Sector column "
                               "(it is not sortable)."),
    ("A planet or moon page", "The seven symbols sit under the title of its page."),
    ("A system page", "Each planet and moon row carries the seven symbols in place of the old notes."),
    ("Search", "The planet and moon panels of Search show the same symbols beside each result."),
)
"""tuple: `(place, what it shows)` for the places the symbols appear (UX.100, UX.155, UX.156)."""


def zone_text():
    """
    What the habitable zone is and where its edges lie, with the numbers read from `physics/orbits` (GEN.146) so the
    page cannot drift from the generator.
    """
    low, high = orbits.KOPPARAPU_TEFF_RANGE_K
    inner, outer = orbits.calculate_habitable_zone(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K)
    return (
        "A star's habitable zone is the ring of orbits where a planet like Earth could keep liquid water on its "
        "surface. Its inner edge is the runaway greenhouse limit, where a planet loses its water to space, and its "
        "outer edge is the maximum greenhouse limit, beyond which a carbon dioxide atmosphere can no longer keep the "
        "surface warm (Kopparapu et al., 2014). Both edges are a flux, in units of the sunlight Earth gets, that "
        "depends on the star's temperature: a cool star's light is redder, so a planet takes more of it at the same "
        "flux and the zone lies nearer in than the star's brightness alone says. "
        f"For the Sun the edges are {inner:.2f} and {outer:.2f} AU. The fits cover stars from {low:,.0f} to "
        f"{high:,.0f} K; a star outside that range (a brown dwarf, a white dwarf, a remnant) takes the nearest end's "
        "edges. The zone is a guide to where to look and is not the rating: a planet is rated by the six domains "
        "above, from its own air, water, temperature and radiation, wherever it orbits."
    )


def page():
    """Everything the page shows."""
    return {
        "symbol_intro": SYMBOL_INTRO,
        "symbol_rows": SYMBOL_ROWS,
        "life_intro": LIFE_INTRO,
        "life_rows": LIFE_ROWS,
        "life_scores": (0, 3, 9),
        "where_intro": WHERE_INTRO,
        "where_rows": WHERE_ROWS,
        "zone_text": zone_text(),
        "tier_numbers": tuple(range(len(habitability_symbols.TIER_CLASSES))),
        "equipment": EQUIPMENT_ITEMS,
        "equipment_ids": EQUIPMENT_ROW_IDS,
        "domain_tiers": DOMAIN_TIERS,
        "domains": domains(),
        "score_rule": SCORE_RULE,
        "equipment_rule": EQUIPMENT_RULE,
        "compact_host_rule": COMPACT_HOST_RULE,
        "chips": CHIPS,
        "chip_ids": CHIP_IDS,
        "domain_ids": habitability_symbols.DOMAIN_ANCHORS,
        "life_ids": habitability_symbols.LIFE_ANCHORS,
        "overall_id": habitability_symbols.SYMBOL_ANCHORS[habitability_symbols.OVERALL],
    }
