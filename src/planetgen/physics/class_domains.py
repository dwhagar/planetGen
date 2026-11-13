# planetgen/physics/class_domains.py

"""
The class-by-domain matrix (GEN.1221; design: reports/class-type-and-domain-scoring-design-2026-10-31.md section 3 and
reports/class-families-and-positions-design-2026-10-31.md)
==========================================================================================================================

Every class of `physics/data/planet_class_fit.json` is judged on the same physical properties, grouped into domains:

| domain | quantities |
|---|---|
| `structure` | the mass (and with GEN.1219 the body type, which is a hard gate) |
| `gravity` | the surface gravity |
| `thermal` | the mean surface temperature |
| `air_amount` | the surface pressure (the airless gate at 10 Pa) |
| `air_chemistry` | the shares of oxygen, carbon dioxide, hydrogen, sulfur dioxide and the reducing gases |
| `water` | the ocean cover, the water mass fraction, the relative humidity, the groundwater and the ammonia of the water |
| `regional` | the hottest and the coldest latitude band (annual means of the seven bands) |
| `day` | the hottest day (the hottest band's mean plus half its swing) |
| `night` | the coldest night (the coldest band's mean less half its swing) |
| `wet_bulb` | the wet-bulb temperature at the hottest day and the stored humidity (`habitability.wet_bulb_c`) |
| `phase` | whether the family's liquid is a liquid at the body's pressure and temperature (`phase_gate`) |

The cell of a class in a domain is the conditions of the fit table whose quantity is in that domain (`matrix`): a core edge and, past it, a shoulder
(the width of the quantity in `fudge_table.json`, GEN.1314). A class with no condition in a domain has the full range there (score 1). `domain_scores` grades a body in each
domain as the mean membership of the cell's conditions; `overall` is the weighted geometric mean of the six main domains (`weights` in the data file, equal to start).

Quantities a body does not store (a row made before GEN.1108 has no bands, before GEN.1221 no humidity) are `None` and score 0, as everywhere in the fit table.

Use:

    class_domains.derived(values)                   # the regional, day, night, wet-bulb, humidity, groundwater and phase quantities of a body's stored values
    class_domains.matrix(cls)                       # {"thermal": [conditions], "water": [...], ...}
    class_domains.domain_scores(cls, found)         # {"thermal": 1.0, "water": 0.4, ...}
"""

import math

from planetgen.physics import atmosphere, constants, habitability, phase_gate

DOMAINS = ("structure", "gravity", "thermal", "air_amount", "air_chemistry", "water")
"""tuple: The six main domains, over which the overall score is taken."""

EXTRA_DOMAINS = ("regional", "day", "night", "wet_bulb", "phase")
"""tuple: The domains of the regional, day, night and wet-bulb temperatures and of the phase gate of a family."""

ALL_DOMAINS = DOMAINS + EXTRA_DOMAINS

PHASE_PREFIX = "phase_"

WATER_TABLE_COVER_MAX = 0.2
"""float: Over none and under this ocean cover a world whose surface holds liquid (`WATER_TABLE_SURFACES`) has land with a water table under it (GEN.1638). It is the cover class H
(a desert world with some water) is cut at in `planet_class_fit.json`; the aquifer model (`aquifers`) only covers worlds with no liquid at the surface, so a world with
a little liquid at the surface never had a stored aquifer and H could not be reached. A test holds it equal to that cut. The reach is not a new number: a table under land with liquid at the surface is as reachable as that liquid."""

WATER_TABLE_SURFACES = ("surface ocean", "brine seep")
"""tuple: The hydrosphere states that hold liquid water at the surface where a land share can have a water table (a hycean or ice-covered ocean has no such land)."""

QUANTITY_DOMAIN = {
    "mass_earth": "structure", "gravity_g": "gravity", "temperature_k": "thermal", "pressure_pa": "air_amount",
    "ocean_fraction": "water", "water_mass_fraction": "water", "relative_humidity": "water", "groundwater_present": "water",
    "hottest_region_k": "regional", "coldest_region_k": "regional", "day_temperature_k": "day", "night_temperature_k": "night",
    "wet_bulb_c": "wet_bulb", "reducing_share": "air_chemistry", "o2_partial_pa": "air_chemistry", "co2_partial_pa": "air_chemistry", "ammonia_water_mass_fraction": "water",
}
"""dict: The domain of each quantity of the fit table that is not a gas share (`share_<gas>` is air chemistry, `phase_<family>` is the phase gate)."""

BAND_COUNT = constants.BAND_COUNT


def domain_of(quantity):
    """The domain a quantity of the fit table belongs to."""
    if quantity.startswith("share_"):
        return "air_chemistry"
    if quantity.startswith(PHASE_PREFIX):
        return "phase"
    return QUANTITY_DOMAIN[quantity]


def _bands(values):
    """`(mean, swing)` of each stored latitude band, or an empty list when the body has none."""
    bands = []
    for number in range(1, BAND_COUNT + 1):
        mean, swing = values.get(f"band{number}_mean_k"), values.get(f"band{number}_swing_k")
        if mean is None:
            return []
        bands.append((mean, swing or 0.0))
    return bands


def relative_humidity(values):
    """The stored relative humidity of a body, else the one its air and temperature give (`habitability_world.relative_humidity`); `None` for a body with no temperature."""
    stored = values.get("relative_humidity")
    if stored is not None:
        return stored
    temperature = values.get("surface_temperature_k")
    if temperature is None:
        return None
    from planetgen.physics import habitability_world
    return habitability_world.relative_humidity(values)


def derived(values):
    """The quantities of the matrix that are not in `planet_classes.quantities`, from a body's stored `values`:
    `hottest_region_k`, `coldest_region_k`, `day_temperature_k`, `night_temperature_k`, `relative_humidity`, `wet_bulb_c`, `groundwater_present` (1.0 or 0.0; 1.0 too for a liquid-surface world under `WATER_TABLE_COVER_MAX` of ocean, GEN.1638, which leaves the stored groundwater and so the habitability scores as they were), `ammonia_water_mass_fraction` (the ammonia over the water's mass, `None` when no ice share is stored) and
    `phase_<family>` (1.0 where the family's liquid is a liquid at the body's temperature and pressure, else 0.0)."""
    found = {}
    bands = _bands(values)
    found["hottest_region_k"] = max(mean for mean, _ in bands) if bands else None
    found["coldest_region_k"] = min(mean for mean, _ in bands) if bands else None
    found["day_temperature_k"] = max(mean + swing / 2.0 for mean, swing in bands) if bands else None
    found["night_temperature_k"] = min(mean - swing / 2.0 for mean, swing in bands) if bands else None
    humidity = relative_humidity(values)
    found["relative_humidity"] = humidity
    day = found["day_temperature_k"]
    found["wet_bulb_c"] = None if day is None or humidity is None else habitability.wet_bulb_c(day - 273.15, humidity)
    from planetgen.physics import surface_liquid_world
    found["ammonia_water_mass_fraction"] = surface_liquid_world.ammonia_water_mass_fraction(values)      # GEN.1258: the stored percent of the liquid, else the body's ice and water
    groundwater = values.get("groundwater")
    if groundwater is None and values.get("surface_temperature_k") is not None:
        from planetgen.physics import habitability_world
        groundwater = habitability_world.groundwater(values)
    found["groundwater_present"] = None if groundwater is None else (0.0 if groundwater == "none" else 1.0)
    cover = values.get("ocean_fraction")
    if (found["groundwater_present"] == 0.0 and values.get("hydrosphere") in WATER_TABLE_SURFACES      # GEN.1638: the land of a mostly dry liquid-surface world
            and cover is not None and 0.0 < cover < WATER_TABLE_COVER_MAX):
        found["groundwater_present"] = 1.0
    temperature, pressure = values.get("surface_temperature_k"), values.get("atmospheric_pressure_pa") or 0.0
    for family in phase_gate.FAMILIES:
        found[PHASE_PREFIX + family] = None if temperature is None else (1.0 if phase_gate.holds(family, temperature, pressure) else 0.0)
    return found


def matrix(cls):
    """`{domain: [conditions]}` for every domain, the conditions of the class `cls` (a `class_fit.FitClass`) that lie in it; an empty list is the full range."""
    cells = {domain: [] for domain in ALL_DOMAINS}
    for condition in cls.when:
        cells[domain_of(condition["quantity"])].append(condition)
    return cells


def domain_scores(cls, found, membership):
    """`{domain: score}` of a body for `cls`: the mean membership of the conditions of the domain (1 for a domain with none). `membership(condition, found)` is
    `class_fit.membership`."""
    scores = {}
    for domain, conditions in matrix(cls).items():
        scores[domain] = 1.0 if not conditions else sum(membership(condition, found) for condition in conditions) / len(conditions)
    return scores


def overall(scores, weights):
    """The weighted geometric mean of the six main domain scores (0 when any is 0)."""
    total = sum(weights[domain] for domain in DOMAINS)
    if any(scores[domain] <= 0.0 for domain in DOMAINS):
        return 0.0
    return math.exp(sum(weights[domain] * math.log(scores[domain]) for domain in DOMAINS) / total)
