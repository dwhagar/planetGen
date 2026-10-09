"""
Mantle redox and atmosphere species for planets and moons (GEN.85).

Each body with air stores the partial pressures, kPa, of ten species
(`SPECIES`) and its mantle's redox state; the free-text `atmosphere` is
written from them. The design is in docs/design/habitability-index.md
section 4 (rows 1 and 2) and docs/design/atmospheres-retention-and-classes.md
section 3.3 (the cold side).

1. **Mantle redox.** The upper mantle's oxygen fugacity, in log units above
   the iron-wustite buffer (`delta_iw`), rises with mass: a big rocky planet
   disproportionates ferrous iron in its lower mantle and sinks the metal,
   leaving an oxidized upper mantle (Earth sits near QFM, about IW+3.5),
   while a small one stays near or below IW (the Moon about IW-1). The mean
   is `REDOX_EARTH_DELTA_IW + REDOX_SLOPE * log10(M / M_earth)` with a
   normal scatter `REDOX_SIGMA` [D: fitted by eye to the Moon, Mars and
   Earth]. Reduced at or below IW, oxidized from `REDOX_OXIDIZED_DELTA_IW`.
   Gas giants have no mantle redox (`None`).
2. **Mixing ratios.** Each class has a mole-fraction template
   (`CLASS_MIXES`, read from today's class atmosphere text); every fraction
   gets a log-normal scatter. A gas outside the ten (helium, ammonia,
   sodium) is the class's `balance`, named in the text but not stored.
3. **Redox shift.** A reduced mantle outgasses H2, CO, CH4 and H2S in place
   of H2O, CO2 and SO2; in an anoxic atmosphere (O2 under 0.1 percent) part
   of each oxidized gas moves to its reduced partner (`REDOX_SHIFT`). An
   oxygenated atmosphere, made by life, keeps its mix.
4. **Cold side.** No species can exceed its vapour pressure at the surface
   temperature (Clausius-Clapeyron through the design doc's 100 kPa and
   0.001 kPa points, `CONDENSATION_K`): the excess is ice, and the total
   surface pressure drops by it. Gas giants have no surface and skip this.
"""

import math

from planetgen.physics import constants
from planetgen.util import draw

SPECIES = ("o2", "co2", "co", "n2", "ar", "h2", "h2o", "ch4", "h2s", "so2")
"""tuple: The stored species, in column order (`p_<species>_kpa`)."""

SPECIES_NAMES = {
    "o2": "oxygen", "co2": "carbon dioxide", "co": "carbon monoxide", "n2": "nitrogen",
    "ar": "argon", "h2": "hydrogen", "h2o": "water vapor", "ch4": "methane",
    "h2s": "hydrogen sulfide", "so2": "sulfur dioxide",
}

PARTIAL_PRESSURE_FIELDS = tuple(f"p_{gas}_kpa" for gas in SPECIES)
"""tuple: The attributes (and columns) holding each species' partial pressure."""

ATMOSPHERE_FIELDS = ("mantle_redox", "mantle_delta_iw") + PARTIAL_PRESSURE_FIELDS
"""tuple: Every attribute (and column) this module sets."""

REDOX_EARTH_DELTA_IW = 3.5
REDOX_SLOPE = 2.0
REDOX_SIGMA = 1.0
REDOX_RANGE = (-6.0, 6.0)
REDOX_OXIDIZED_DELTA_IW = 2.5
"""float: Above this (about QFM - 1) a mantle counts as oxidized [D]."""

MIX_SCATTER_SIGMA = 0.25
"""float: Log-normal sigma on every template fraction."""

ANOXIC_O2_FRACTION = 1e-3

REDOX_SHIFT = {
    # Fraction of each oxidized gas that moves to its reduced partners.
    "reduced": 0.6,
    "intermediate": 0.15,
    "oxidized": 0.0,
}
REDUCED_PARTNERS = {
    "co2": (("co", 0.8), ("ch4", 0.2)),
    "so2": (("h2s", 1.0),),
    "h2o": (("h2", 1.0),),
}

CONDENSATION_K = {
    # (T at 100 kPa, T at 0.001 kPa), atmospheres-retention-and-classes.md 3.3
    "h2": (20.0, 6.0), "n2": (77.0, 33.0), "co": (81.0, 36.0), "ar": (87.0, 38.0),
    "o2": (90.0, 40.0), "ch4": (112.0, 48.0), "co2": (195.0, 112.0),
    "so2": (263.0, 131.0), "h2o": (373.0, 212.0),
}
"""dict: H2S has no row; it condenses near SO2 and is capped the same way."""
CONDENSATION_K["h2s"] = CONDENSATION_K["so2"]

CLASS_MIXES = {
    # class: (mole fractions of the ten species, balance gas name or None)
    "A": ({"co2": 0.60, "so2": 0.35, "n2": 0.05}, None),
    "B": ({"o2": 0.10}, ("helium", "sodium")),
    "E": ({"h2o": 0.55, "ch4": 0.25, "n2": 0.05, "h2": 0.05}, "ammonia"),
    "F": ({"co2": 0.60, "ch4": 0.20, "n2": 0.10}, "ammonia"),
    "G": ({"co2": 0.50, "n2": 0.35, "o2": 0.10, "h2o": 0.04, "ar": 0.01}, None),
    "H": ({"n2": 0.70, "o2": 0.25, "ar": 0.04, "h2o": 0.01}, None),
    "I": ({"h2": 0.86, "ch4": 3e-3, "h2o": 1e-3, "h2s": 1e-4}, "helium"),
    "J": ({"h2": 0.86, "ch4": 3e-3, "h2o": 1e-3, "h2s": 1e-4}, "helium"),
    "K": ({"co2": 0.95, "n2": 0.028, "ar": 0.019, "o2": 1.5e-3, "co": 7e-4, "h2o": 3e-4}, None),
    "L": ({"ar": 0.60, "o2": 0.30, "n2": 0.08, "co2": 0.02}, None),
    "M": ({"n2": 0.7703, "o2": 0.21, "ar": 9.3e-3, "h2o": 0.01, "co2": 4e-4, "ch4": 1.9e-6}, None),
    "N": ({"co2": 0.965, "n2": 0.035, "so2": 1.5e-4, "h2s": 1e-4, "ar": 7e-5, "co": 1.7e-5, "h2o": 2e-5}, None),
    "O": ({"n2": 0.74, "o2": 0.22, "h2o": 0.03, "ar": 9e-3, "co2": 1e-3}, None),
    "P": ({"n2": 0.75, "o2": 0.22, "ar": 0.02, "co2": 4e-3, "h2o": 1e-3}, None),
    "Q": ({"n2": 0.78, "o2": 0.15, "ar": 0.05, "co2": 0.01, "h2o": 0.01}, None),
    "T": ({"h2": 0.84, "ch4": 0.02, "h2o": 1e-3}, "helium"),
    "V": ({"co2": 0.85, "o2": 0.08, "n2": 0.04, "h2": 0.01}, "helium"),
}
"""dict: Each class's mix, from its `tuning.PLANET_CLASSES` atmosphere text
and the solar-system analogue it is modelled on (K Mars, M Earth, N Venus,
I/J/T Jupiter-like) [D]. Fractions plus the balance sum to 1."""

DEFAULT_MIX = ({"n2": 0.78, "o2": 0.2, "ar": 0.01, "co2": 0.01}, None)

THIN_KPA = 1.0
DENSE_KPA = 500.0
MAJOR_FRACTION = 0.01
TRACE_FRACTION = 1e-4
HUMID_H2O_FRACTION = 0.02


def redox_class(delta_iw):
    """'reduced' (at or below IW), 'oxidized' (from
    `REDOX_OXIDIZED_DELTA_IW`) or 'intermediate'."""
    if delta_iw <= 0.0:
        return "reduced"
    if delta_iw >= REDOX_OXIDIZED_DELTA_IW:
        return "oxidized"
    return "intermediate"


def draw_mantle_delta_iw(mass_kg):
    """A rocky body's upper-mantle oxygen fugacity, log units above IW."""
    mass_earth = max(mass_kg / constants.EARTH_MASS_TO_KG, 1e-6)
    mean = REDOX_EARTH_DELTA_IW + REDOX_SLOPE * math.log10(mass_earth)
    low, high = REDOX_RANGE
    return min(high, max(low, draw.gauss(mean, REDOX_SIGMA)))


def vapour_pressure_kpa(gas, temperature_k):
    """The most of `gas` that can be vapour at `temperature_k`, kPa
    (Clausius-Clapeyron through `CONDENSATION_K`'s two points)."""
    t_100, t_low = CONDENSATION_K[gas]
    latent_over_r = math.log(1e5) / (1.0 / t_low - 1.0 / t_100)
    return 100.0 * math.exp(-latent_over_r * (1.0 / max(temperature_k, 1.0) - 1.0 / t_100))


def draw_mixing_ratios(planet_class, redox):
    """
    The ten species' mole fractions and the balance's, for a body of
    `planet_class` over a mantle of `redox` (`None` for a gas giant).

    Returns:
        tuple: (dict gas -> fraction, balance fraction, balance name or None)
    """
    template, balance_name = CLASS_MIXES.get(planet_class, DEFAULT_MIX)
    fractions = {gas: 0.0 for gas in SPECIES}
    for gas, value in template.items():
        fractions[gas] = value * math.exp(draw.gauss(0.0, MIX_SCATTER_SIGMA))
    balance = 0.0
    if balance_name is not None:
        balance = (1.0 - sum(template.values())) * math.exp(draw.gauss(0.0, MIX_SCATTER_SIGMA))
    total = sum(fractions.values()) + balance
    fractions = {gas: value / total for gas, value in fractions.items()}
    balance /= total

    shift = REDOX_SHIFT.get(redox, 0.0)
    if shift and fractions["o2"] < ANOXIC_O2_FRACTION:
        for gas, partners in REDUCED_PARTNERS.items():
            moved = fractions[gas] * shift
            fractions[gas] -= moved
            for partner, share in partners:
                fractions[partner] += moved * share
    return fractions, balance, balance_name


def describe(partials_kpa, balance_kpa, balance_name):
    """The free-text atmosphere, e.g. "a mix of nitrogen, oxygen, and
    argon, with traces of water vapor and carbon dioxide"."""
    total = sum(partials_kpa.values()) + balance_kpa
    if total <= 0.0:
        return "None"
    parts = [(value / total, SPECIES_NAMES[gas]) for gas, value in partials_kpa.items() if value > 0.0]
    if balance_name is not None and balance_kpa > 0.0:
        names = (balance_name,) if isinstance(balance_name, str) else balance_name
        parts += [(balance_kpa / total / len(names), name) for name in names]
    parts.sort(key=lambda part: -part[0])
    majors = [name for fraction, name in parts if fraction >= MAJOR_FRACTION][:4]
    traces = [name for fraction, name in parts if TRACE_FRACTION <= fraction < MAJOR_FRACTION][:2]

    qualifiers = []
    if total < THIN_KPA:
        qualifiers.append("thin")
    elif total > DENSE_KPA:
        qualifiers.append("dense")
    if partials_kpa.get("h2o", 0.0) / total >= HUMID_H2O_FRACTION:
        qualifiers.append("humid")
    reduced = sum(partials_kpa.get(gas, 0.0) for gas in ("h2", "co", "ch4", "h2s"))
    oxidized = sum(partials_kpa.get(gas, 0.0) for gas in ("o2", "co2", "so2"))
    if partials_kpa.get("h2", 0.0) / total < 0.5 and reduced > oxidized:
        qualifiers.append("reducing")
    lead = "a " + ", ".join(qualifiers) + " " if qualifiers else "a "

    if len(majors) == 1:
        text = f"{lead}mix of mostly {majors[0]}"
    else:
        text = f"{lead}mix of {_join(majors)}"
    if traces:
        text += f", with traces of {_join(traces)}"
    return text


def _join(names):
    if len(names) <= 2:
        return " and ".join(names)
    return ", ".join(names[:-1]) + ", and " + names[-1]


def generate_atmosphere(planet):
    """
    Sets `planet`'s `ATMOSPHERE_FIELDS` and its `atmosphere` text from its
    class, mass, surface temperature and `atmospheric_pressure` (Pa), and
    lowers `atmospheric_pressure` by any gas the cold side freezes out. An
    airless body (`atmosphere == "None"`) keeps its redox and gets zero
    partial pressures.
    """
    giant = planet.body_type == "g"
    if giant:
        planet.mantle_delta_iw = None
        planet.mantle_redox = None
    else:
        planet.mantle_delta_iw = draw_mantle_delta_iw(planet.mass)
        planet.mantle_redox = redox_class(planet.mantle_delta_iw)

    if planet.atmosphere == "None" or not planet.atmospheric_pressure:
        for field in PARTIAL_PRESSURE_FIELDS:
            setattr(planet, field, 0.0)
        return

    total_kpa = planet.atmospheric_pressure / 1000.0
    fractions, balance, balance_name = draw_mixing_ratios(planet.planet_class, planet.mantle_redox)
    partials = {gas: fraction * total_kpa for gas, fraction in fractions.items()}
    balance_kpa = balance * total_kpa
    if not giant:
        for gas, value in partials.items():
            if value > 0.0:
                partials[gas] = min(value, vapour_pressure_kpa(gas, planet.surface_temperature))
        planet.atmospheric_pressure = (sum(partials.values()) + balance_kpa) * 1000.0
    for gas, value in partials.items():
        setattr(planet, f"p_{gas}_kpa", value)
    planet.atmosphere = describe(partials, balance_kpa, balance_name)


def partial_pressures_kpa(body):
    """`body`'s stored partial pressures as a dict gas -> kPa (zeros for a
    body saved before GEN.85)."""
    return {gas: getattr(body, f"p_{gas}_kpa", None) or 0.0 for gas in SPECIES}


def hydrogen_kpa(body):
    """The hydrogen in `body`'s air as an H2-equivalent partial pressure,
    kPa: H2 + H2O + 2 CH4 + H2S. The water-loss clock's inventory."""
    p = partial_pressures_kpa(body)
    return p["h2"] + p["h2o"] + 2.0 * p["ch4"] + p["h2s"]
