# planetgen/web/lib/habitability_ref.py

"""
The text of the habitability page (`/help/habitability`, UX.90, UX.180):
the survival equipment items, the six PHI-6 domains with their Blue, Green,
Yellow and Red thresholds, the compact-host rule and the chips a planet row
shows. The equipment items come from `planetgen.physics.equipment`;
the thresholds are read from `planetgen.physics.habitability` and `protection_limits`, where they are named with
their sources (`domains`), so the tables cannot disagree with the score.
"""

from planetgen import tuning
from planetgen.physics import atmosphere, constants, habitability, latitude_bands, life_classes, orbits
from planetgen.physics.habitability_world import LETHAL_DOSE_SV_YR
from planetgen.physics.protection_limits import limit
from planetgen.web.lib import habitability_symbols, help_registry

EQUIPMENT_ITEMS = (
    ("Clothing", "Cold-weather, cooling or sun-cover clothing: the temperature is a little outside comfort or the "
                 "ultraviolet is strong."),
    ("Thermal suit", "An insulated, liquid-cooled or heated diving suit: the temperature is well outside what "
                     "clothing handles, or the pressure needs a heated diving suit."),
    ("Sealed suit", "A fully enclosed suit with its own climate: heated for the deep cold, heat-shielded for the "
                    "hot, chemical-barrier where a corrosive gas or liquid is about."),
    ("Pressure suit", "The ambient pressure is too low for a mask to give enough oxygen."),
    ("Armoured suit", f"A rigid suit that holds one atmosphere inside: crushing pressure or heat past {limit('armoured_heat_c'):g} C."),
    ("Vehicle only", f"No suit exists for it: heat past {limit('vehicle_heat_c'):g} C, pressure past {limit('vehicle_pressure_kpa') / 1000:g} MPa "
                      f"or gravity past {limit('vehicle_gravity_g'):g} g."),
    ("Powered frame", f"The gravity is past {limit('frame_gravity_g'):g} g (recommended from {limit('frame_gravity_rec_g'):g} g)."),
    ("Dust filter or seals", "Fine regolith dust wrecks open seals (a rocky body with almost no air), or the dust in the "
                            f"air passes {limit('dust_recommended_mg_m3'):g} mg/m3 (recommended) or "
                            f"{limit('dust_required_mg_m3'):g} mg/m3 (required)."),
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
        name = atmosphere.species_name(gas.lower()).capitalize()
        rows.append((f"{name} (a filter takes it out)", _under(band[4], "kPa"), "", _up_to(band[6], "kPa")))
    return rows


def domains(bands=habitability.HUMAN):
    """
    `(domain, [(variable, blue, green, yellow)])`, beyond the Yellow column a variable is Red. Every number is
    read from `bands` (the limits in `protection_limits`) and the constants of `habitability`, so the page cannot
    drift from the score.
    """
    pressure, temperature, dose, uv = bands.pressure, bands.temperature, bands.dose, bands.uv
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
            ("Ultraviolet light at the ground (DNA-damaging dose)", _under(uv[4], "times Earth's"), _up_to(uv[5], "times Earth's"),
             _up_to(uv[6], "times Earth's")),
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


_BLUE, _GREEN, _YELLOW = habitability._EDGE_SCORES[3], habitability._EDGE_SCORES[2], habitability._EDGE_SCORES[1]
"""float: The domain score at the edge of Blue, Green and Yellow (`habitability._EDGE_SCORES`)."""

SCORE_RULE = (
    f"Each variable scores {_BLUE:g} in Blue and falls to {_GREEN:g} at the edge of Green, {_YELLOW:g} at the edge of Yellow and 0 "
    f"at the outer edge of Red. A domain takes the score of its worst variable: Blue at {_BLUE:g}, Green from {_GREEN:g}, "
    f"Yellow from {_YELLOW:g}, Red below. The PHI-6 score on a "
    "planet's tooltip is the geometric mean of the six domains, held under the colour gates: Radiation and "
    "Temperature cap it at their own colour, Pressure in Yellow or Red caps it there (unless a full suit is "
    "needed anyway), Air, Water or a Red Suitability cap it at Yellow, and it is Blue only if all six are. "
    "Suitability, the long-term domain, is the lowest of water, food and shelter, and is never better than "
    "the colour the other five allow."
)

EQUIPMENT_RULE = (
    "A planet's equipment line lists only what a visitor needs, never what is not needed. Each item is required "
    "(past the physical limit) or recommended (in the band before it, where a person is impaired or at risk, or "
    "where only a region or the night asks for it); a body gives a Required line and a Recommended line, and leaves off one with nothing in it. The line follows the "
    "physics, so it can be stricter than a domain's colour: a gas past its immediate danger limit needs sealed "
    "air even where the score only says Yellow. Radiation goes up in steps: a shielded garment, a storm vest, a "
    "hardened suit, then a shielded shelter and a hardened habitat, each worked out from how much of the dose is "
    "cosmic rays and how much is storm particles. Where no lighter gear brings the year under the Green edge (cosmic "
    "rays on Mars or the Moon) the hardened suit is required and the hardened habitat recommended. "
    f"Air hotter than {limit('heat_shielded_c'):g} C burns a mask-only "
    "visitor's skin at the seal, so it needs a heat-shielded suit."
)

COMPACT_HOST_RULE = (
    f"A planet of a pulsar, a neutron star or a black hole is rated as if its surface took "
    f"{LETHAL_DOSE_SV_YR:,.0f} Sv a year. That puts radiation in Red, so the planet needs a hardened habitat "
    "and its life scores are 0. Its note names the host."
)

CHIPS = (
    ("Habitable", "The world's overall PHI-6 colour is Blue or Green, so it needs at most a mask. It is not a promise "
                  "that anything lives there. A habitable world is shown as this chip alone, not with Habitable class as well."),
    ("Habitable moon", "A moon of this planet is habitable in that sense (shown instead of Habitable-class moon when both hold)."),
    ("Habitable class", "The planet's class can carry life (an Earth-like or similar world). It is not a rating "
                        "and not a promise that anything lives there: the equipment level beside it is the rating."),
    ("Habitable-class moon", "A moon of this planet is of a habitable class in that sense."),
    ("Inhabited", "A species lives here. Such a world is habitable and of a habitable class, so it shows this chip alone, not those two as well."),
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

LIFE_NEEDS = (
    ("microbial", "Any step of life at all, from abiogenesis on. How much comes from how kind the surface is to microbes."),
    ("flora", (f"From photosynthesis. Until life is multicellular it is mats and algae and stays at "
               f"{tuning.LIFE_SCORE_SIMPLE_FLORA_MAX} or less; after that it follows how kind the surface is to complex life.")),
    ("fauna", ("From multicellular life. How much comes from how kind the surface is to complex life, and never more "
               "than the Flora the animals feed on.")),
    ("sentient", ("Only where a technological civilization has arisen. It takes Fauna's number and scales it by how long "
                  f"the civilization has existed, in full after {tuning.LIFE_SCORE_SENTIENT_FULL_GYR * 1000:g} million "
                  "years; never below 1 once there is one.")),
)
"""tuple: `(life score name, what it needs and what sets how much)`, the numbers read from `planetgen.tuning`."""

LIFE_FOOD = (
    "A 0 means none of that kind. A world with a life chemistry always has a timeline, so these four read where it "
    "got to. Food for a visitor comes from the same numbers: any Flora is good food, Flora together with Fauna is "
    "better, Fauna on its own counts, and microbes alone are poor food. Ice over the life puts it out of reach. "
    f"The scale tops out at {tuning.LIFE_SCORE_MAX}; a surface too hostile to a kind of life scores it 0, whatever stage "
    "the world reached."
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
    ("The equipment chips", "A system page row and the Equipment column of Search show one chip per item of the "
                            "equipment line: solid for a required item, outlined for a recommended one, in the colour "
                            "of the domain it comes from."),
)
"""tuple: `(place, what it shows)` for the places the symbols appear (UX.100, UX.155, UX.156)."""


def _class_ids():
    """The first life class ids, as the map names them (`1a, 1b, 1c, 2, 3a`)."""
    return ", ".join(row.id for row in life_classes.CLASSES[:5])


MAP_INTRO = (
    "The Galaxy Map shows where habitable worlds are, in two ways that you choose. A habitable world here means a "
    "planet or moon whose overall PHI-6 colour is Blue or Green; Yellow and Red worlds are not counted."
)

MAP_ROWS = (
    ("Color by: Habitable worlds", "In the map's Menu, Color by paints each sector block by how many habitable worlds "
                                   "its systems hold, from none to the most, with a legend; the choice stays in the "
                                   "address as color=habitable."),
    ("Marker: Habitable worlds", "In the Menu's marker list (off until you turn it on, and kept in the address as "
                                 "show=marker-habitable) every star system with at least one habitable world gets a "
                                 "marker. Hovering names the system and how many habitable worlds it has; a click "
                                 "opens its summary."),
    ("Marker: Life", "The Menu's Points of interest group also has three markers for worlds with life: Life (any of the four "
                     "life scores is 1 or more), Sentient life and Civilizations. They are off until you turn them on "
                     "(show=marker-life, marker-sentient, marker-civilization in the address); a system with several is "
                     "marked once, by the strongest."),
    ("Marker: World classes", "The Menu's World classes group has a marker for each planet class, all off until you turn "
                              "one on (kept in the address as show=marker-class-K for class K). A system with a world of "
                              "the class gets a disc with the class letter, once the map is close enough to have its "
                              "sector; hovering names the system and how many worlds of the class it holds."),
    ("Marker: Life classes", f"The Menu's Life classes group has a marker for each of the {len(life_classes.CLASSES)} life classes ({_class_ids()} and "
                             "so on), all off until you turn one on (kept in the address as show=marker-lifeclass-1a for "
                             "class 1a). A system with a world on which that class of life lives gets a green disc with the "
                             "class id, once the map is close enough to have its sector; hovering names the system and how "
                             "many worlds have the class."),
)
"""tuple: `(control, what it does)` for the habitable-world views of the Galaxy Map (MAP.143, MAP.220)."""

_NUMBER_WORDS = {3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight"}


def temperature_text():
    """The temperature paragraph, with the bands read from `physics/latitude_bands.py` (their number, their width and each one's share of the surface)."""
    count = latitude_bands.BAND_COUNT
    shares = [latitude_bands._area_share(high, low) for high, low in latitude_bands._bounds()]
    middle = count // 2
    return (
        "The surface temperature of a world is its mean, and the ratings use the mean. The Temperature row of a planet's "
        f"page also says how the mean is spread: a world that turns on its own has an annual mean for each of {_NUMBER_WORDS.get(count, count)} bands "
        f"of latitude, each {latitude_bands.BAND_WIDTH_DEG:.1f} degrees wide. The equator band reaches {latitude_bands.BAND_WIDTH_DEG / 2:.1f} degrees "
        f"each side of the equator and holds {shares[middle] * 100:.0f} percent of the surface; the bands on each side of it hold "
        f"{', '.join(f'{share * 100:.0f}' for share in shares[:middle][::-1][:-1])} and, for the one that takes the pole, "
        f"{shares[0] * 100:.0f} percent. Each is found from the sunlight it gets over a year at the "
        "world's axial tilt, with a thick atmosphere carrying heat toward the poles. It also gives the day and night swing "
        "of each band, which a slow turn on bare ground makes as large as the mean itself and a thick atmosphere shrinks "
        "to a few degrees. A world that is tidally locked has no day and night; its page gives rings from the day side, through the "
        "terminator ring, to the night side. Where a band or the night needs heavier gear than the mean does, the row names "
        "the place and the gear as recommended."
    )


def radiation_text():
    """
    Where the Radiation rating comes from and what sets each term, with the limits read from `protection_limits` so the page
    cannot drift from the score (GEN.445, DOC.75).
    """
    return (
        "The Radiation rating is the worse of two scores: the yearly dose and the ground ultraviolet. The dose adds four terms, "
        "each set by the world's make-up. Cosmic rays are cut by the air above the ground (its column, the pressure over the "
        "gravity, in grams per square centimetre) and by the dipole magnetic field; a body with no air takes them whole. "
        "Stellar particles come from the star's flares, fall with the square of the distance, and are cut by the same air and "
        "field. The ground gives the gamma rays and the radon of the potassium, thorium and uranium in its rock: the same "
        "stored content that heats the interior, so a world rich in them is hotter inside and carries a larger ground dose. "
        "The ultraviolet is the DNA-damaging dose at the ground as a multiple of Earth's: the star's own light through the "
        "world's air, its ozone (made from its oxygen), haze, dust and cloud, and the unfiltered light for a world with no air. "
        f"The dose is Blue under {limit('dose_blue_sv_yr') * 1000:g} mSv a year, Green to {limit('dose_green_sv_yr') * 1000:g}, "
        f"Yellow to {limit('dose_yellow_sv_yr') * 1000:g} mSv, and the ultraviolet is Blue under {limit('uv_cover_index'):g} times "
        f"Earth's, Green to {limit('uv_suit_index'):g} times and Yellow to {limit('uv_lethal_index'):g} times. A planet's page lists "
        "each term with the air, the field, the ozone and the rock that set it, and says which of the two sets the rating."
    )


def zone_text():
    """
    What the habitable zone is and where its edges lie, with the numbers read from `physics/orbits` (GEN.146) so the
    page cannot drift from the generator.
    """
    low, high = orbits.KOPPARAPU_TEFF_RANGE_K
    inner, outer = orbits.calculate_habitable_zone(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K)
    light_mass, heavy_mass = orbits.KOPPARAPU_MASS_RANGE_EARTH
    light_inner = orbits.hz_limits(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K, light_mass)[0]
    heavy_inner = orbits.hz_limits(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K, heavy_mass)[0]
    return (
        "A star's habitable zone is the ring of orbits where a planet like Earth could keep liquid water on its "
        "surface. Its inner edge is the runaway greenhouse limit, where a planet loses its water to space, and its "
        "outer edge is the maximum greenhouse limit, beyond which a carbon dioxide atmosphere can no longer keep the "
        "surface warm [R2447]. Both edges are a flux, in units of the sunlight Earth gets, that "
        "depends on the star's temperature: a cool star's light is redder, so a planet takes more of it at the same "
        "flux and the zone lies nearer in than the star's brightness alone says. "
        f"For the Sun the edges are {inner:.2f} and {outer:.2f} AU. The inner edge moves with the planet's mass, because a "
        f"heavier planet holds its water to a higher flux: for the Sun it is {light_inner:.2f} AU at {light_mass:g} Earth "
        f"masses and {heavy_inner:.2f} AU at {heavy_mass:g}; the star's own zone is for one Earth mass. A cool star is "
        "also brightest when young, so a planet in its zone today may have crossed the runaway limit early and lost "
        "water, and its oceans are smaller for it. "
        f"The fits cover stars from {low:,.0f} to "
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
        "life_needs": dict(LIFE_NEEDS),
        "life_food": LIFE_FOOD,
        "where_intro": WHERE_INTRO,
        "where_rows": WHERE_ROWS,
        "map_intro": MAP_INTRO,
        "map_rows": MAP_ROWS,
        "temperature_text": temperature_text(),
        "radiation_text": radiation_text(),
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


help_registry.register_category("worlds", "Worlds and life", "How a visitor's gear is judged on a world and how life on it is described.", order=20)
help_registry.register_topic("habitability_levels", "Habitability", "worlds",
                             "What the equipment levels, the Blue to Red domain colors and the Habitable chips mean.", order=20)
