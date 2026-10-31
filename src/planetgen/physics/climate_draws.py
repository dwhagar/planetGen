# planetgen/physics/climate_draws.py

"""
The albedo and the greenhouse quantity a body's temperature is made of (GEN.412).

`planets.calculate_atmospheric_conditions` draws a Bond albedo A from the range of the body's regime (`albedo_range`, GEN.684) and a greenhouse quantity G (the
CO2-likeness of the air), and the surface temperature is exactly what they give:

    T_e = ((1 - A) S / (4 sigma))^(1/4)         the temperature with no air
    T_s = T_e (1 + G)^(1/4)                      G = (T_s / T_e)^4 - 1

`S` is the starlight at the body's orbit (W/m2; a moon's is its parent planet's) and nothing is colder than the cosmic
microwave background, so `T = max(2.725 K, formula)`. The generator now keeps both numbers on the body (`bond_albedo`,
`greenhouse_factor`) instead of the temperature alone, and `albedo_basis` says where they came from. The design is
docs/design/composition/atmosphere-colour.md sections 5 and 8.

A galaxy made before the columns existed keeps its stored scale height, whose barometric law holds the no-air
temperature, so `from_stored` recovers both draws exactly (`BASIS_INVERSE`) and says when it cannot, and
`backfill_draw` then takes a seeded draw from the range of the row's regime with G worked out from the stored temperature, so the
stored temperature never moves (`BASIS_BACKFILL_DRAW`).
"""

import math

from planetgen import tuning
from planetgen.physics import albedo_color, constants
from planetgen.util import draw

CLIMATE_FIELDS = ("bond_albedo", "greenhouse_factor", "albedo_basis")
"""tuple: The attributes (and `planets`/`moons` columns) this module sets."""

BASIS_GENERATED = 0
BASIS_INVERSE = 1
BASIS_BACKFILL_DRAW = 2
BASIS_NAMES = {BASIS_GENERATED: "generated", BASIS_INVERSE: "exact inverse", BASIS_BACKFILL_DRAW: "backfill draw"}
"""dict: `albedo_basis`: the generator recorded the draws (0); worked out exactly from a stored row (1); a seeded draw
from the range of the body's regime because the row could not be inverted (2)."""

DEFAULT_ALBEDO_RANGE = (0.12, 0.35)
"""tuple: The Bond albedo range of a class that declares none (the generator's own default)."""

GIANT_BOND_ALBEDO_RANGE = (0.29, 0.503)
"""tuple: The Bond albedo a generated giant draws from (GEN.1030): the lowest and the highest of the four measured in the Solar System, Neptune 0.29 (Pearl and Conrath 1991), Uranus 0.30 (Pearl et al. 1990), Saturn 0.342 and Jupiter 0.503 (Li et al. 2018). Registry row CC9."""

GIANT_ALBEDO_RANGE = (0.02, 0.91)
"""tuple: The open bracket for a giant with no class range (Sudarsky 2000: about 0.02 to 0.8; Cahoy 2010: Bond up to
0.914), used only to backfill a row (docs/design/composition/atmosphere-colour.md section 8)."""

GREENHOUSE_RANGE = (0.0, 500.0)
"""tuple: The greenhouse quantity an inversion may give; outside it the row is not one this code built."""

PRESSURE_TOLERANCE = 1e-6
"""float: How far above the barometric law's pressure a stored one may sit before the row is called edited."""

FLOOR_TOLERANCE = 1e-9
"""float: A temperature this close to the 2.725 K floor counts as on it (the equation is then an inequality)."""


def starlight_w_m2(luminosity_w, distance_km):
    """The starlight at `distance_km` from a star of `luminosity_w`, W/m2 (the generator's own form: W per km2, 1e-6)."""
    return (luminosity_w / (4.0 * math.pi * distance_km ** 2)) / 1e6


def equilibrium_temperature_k(albedo, flux_w_m2):
    """T_e, K: the temperature the starlight gives a body of Bond albedo `albedo` with no air, floored at the cosmic
    microwave background."""
    return max(constants.COSMIC_BACKGROUND_TEMPERATURE_K,
               ((1.0 - albedo) * flux_w_m2 / (4.0 * constants.STEFAN_BOLTZMANN_CONSTANT)) ** 0.25)


def surface_temperature_k(albedo, greenhouse_factor, flux_w_m2):
    """T_s, K: `T_e (1 + G)^(1/4)` with the same floor. `greenhouse_factor` 0 is no air."""
    return max(constants.COSMIC_BACKGROUND_TEMPERATURE_K,
               ((1.0 - albedo) * flux_w_m2 * (1.0 + greenhouse_factor)
                / (4.0 * constants.STEFAN_BOLTZMANN_CONSTANT)) ** 0.25)


def greenhouse_from_temperature(surface_k, albedo, flux_w_m2):
    """G = (T_s / T_e)^4 - 1 for a stored `surface_k`, `albedo` and starlight (the inverse of `surface_temperature_k`
    for a temperature above the floor). `None` when the starlight is zero."""
    if not flux_w_m2 or flux_w_m2 <= 0.0 or albedo >= 1.0:
        return None
    return 4.0 * constants.STEFAN_BOLTZMANN_CONSTANT * surface_k ** 4 / ((1.0 - albedo) * flux_w_m2) - 1.0


def record(planet, albedo, greenhouse_factor):
    """Sets the three fields on a body just drawn."""
    planet.bond_albedo = albedo
    planet.greenhouse_factor = greenhouse_factor
    planet.albedo_basis = BASIS_GENERATED


def from_stored(row, flux_w_m2):
    """
    The exact albedo and greenhouse quantity of a stored body, or `None` with the reason when it has no solution.

    `row` holds `surface_temperature_k`, `atmosphere` (`"None"` for no air), `scale_height_km`, `atm_molar_density` (kg/mol),
    `gravity_g`, `atm_density` and `atmospheric_pressure_pa`; `flux_w_m2` is the starlight at the body's orbit. The scale
    height keeps T_e (`H = R T_e / (M g)`), so `1 - A = 4 sigma T_e^4 / S` and `G = (T_s / T_e)^4 - 1` follow, with no draw
    (airless: `T_e = T_s`, `G = 0`).

    Returns:
        tuple: `(albedo, greenhouse_factor, None)`, or `(None, None, reason)`.
    """
    if not flux_w_m2 or flux_w_m2 <= 0.0:
        return None, None, "no starlight to work from"
    surface_k = row.get("surface_temperature_k")
    if not surface_k or surface_k <= 0.0:
        return None, None, "no surface temperature"
    floor = constants.COSMIC_BACKGROUND_TEMPERATURE_K * (1.0 + FLOOR_TOLERANCE)
    if surface_k <= floor:
        return None, None, "the surface sits on the background floor"
    sigma = constants.STEFAN_BOLTZMANN_CONSTANT
    airless = row.get("atmosphere") in (None, "None")
    if airless:
        equilibrium = surface_k
        greenhouse = 0.0
    else:
        height_m = (row.get("scale_height_km") or 0.0) * 1000.0
        molar = row.get("atm_molar_density") or 0.0
        gravity = (row.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
        if height_m <= 0.0 or molar <= 0.0 or gravity <= 0.0:
            return None, None, "the scale height, molar mass or gravity is missing"
        equilibrium = height_m * molar * gravity / constants.R
        if equilibrium <= floor:
            return None, None, "the no-air temperature sits on the background floor"
        pressure = row.get("atmospheric_pressure_pa")
        density = row.get("atm_density")
        if pressure is not None and density:
            ceiling = density * gravity * height_m * (row.get("gravity_g") or 0.0)
            if pressure > ceiling * (1.0 + PRESSURE_TOLERANCE):
                return None, None, "the pressure is above the barometric law (an edited row)"
        greenhouse = (surface_k / equilibrium) ** 4 - 1.0
    albedo = 1.0 - 4.0 * sigma * equilibrium ** 4 / flux_w_m2
    if not 0.0 <= albedo < 1.0:
        return None, None, f"the albedo {albedo:.3g} is outside 0 to 1"
    if not GREENHOUSE_RANGE[0] - 1e-9 <= greenhouse <= GREENHOUSE_RANGE[1]:
        return None, None, f"the greenhouse quantity {greenhouse:.3g} is outside {GREENHOUSE_RANGE}"
    return albedo, max(0.0, greenhouse), None


def albedo_range(body_type, water_mass_fraction, flux_w_m2, star_temperature_k=None):
    """
    The Bond albedo range a body draws from before its air is solved (GEN.684): `GIANT_BOND_ALBEDO_RANGE` for a giant (GEN.1030),
    `DEFAULT_ALBEDO_RANGE` for a dry body and a body with no starlight; a body with surface water (`tuning.ALBEDO_WET_MIN_WATER_FRACTION` or more) is ice-covered and
    bright (`tuning.ALBEDO_FROZEN_RANGE`) where its no-air equilibrium temperature at `tuning.ALBEDO_REGIME_REFERENCE` is under
    `tuning.ALBEDO_FROZEN_BELOW_K`, else ocean and cloud (`tuning.ALBEDO_WET_RANGE`). The two bright ranges are scaled to the colour of the
    star of `star_temperature_k` (`albedo_color.shifted_range`, GEN.912); with no star temperature they are the Sun's.
    """
    if body_type == "g":
        return GIANT_BOND_ALBEDO_RANGE
    if (not water_mass_fraction or water_mass_fraction < tuning.ALBEDO_WET_MIN_WATER_FRACTION
            or not flux_w_m2 or flux_w_m2 <= 0.0):
        return DEFAULT_ALBEDO_RANGE
    if equilibrium_temperature_k(tuning.ALBEDO_REGIME_REFERENCE, flux_w_m2) < tuning.ALBEDO_FROZEN_BELOW_K:
        return albedo_color.shifted_range(tuning.ALBEDO_FROZEN_RANGE, "frozen", star_temperature_k)
    return albedo_color.shifted_range(tuning.ALBEDO_WET_RANGE, "wet", star_temperature_k)


def backfill_albedo_range(row, flux_w_m2):
    """The range a stored row's albedo is backfilled from: the open giant bracket for a giant, else `albedo_range`."""
    if row.get("body_type") == "g":
        return GIANT_ALBEDO_RANGE
    return albedo_range(row.get("body_type"), row.get("water_mass_fraction"), flux_w_m2)


def backfill_draw(row, flux_w_m2, seed):
    """
    A seeded albedo from the range of the row's regime (`backfill_albedo_range`) for a row `from_stored` could not invert, and the greenhouse quantity the
    stored temperature then needs, so the stored temperature does not move. The albedo is held at or above
    `1 - 4 sigma T_s^4 / S` (G is never negative). `seed` keys the draw (the body's table and id).

    Returns:
        tuple: `(albedo, greenhouse_factor)`, the factor `None` when there is no starlight.
    """
    stream = draw.Stream(f"albedo-backfill:{seed}")
    low, high = backfill_albedo_range(row, flux_w_m2)
    surface_k = row.get("surface_temperature_k") or 0.0
    if flux_w_m2 and flux_w_m2 > 0.0 and surface_k > 0.0:
        low = max(low, 1.0 - 4.0 * constants.STEFAN_BOLTZMANN_CONSTANT * surface_k ** 4 / flux_w_m2)
    high = max(high, low)
    albedo = min(stream.uniform(low, high), 0.999)
    return albedo, greenhouse_from_temperature(surface_k, albedo, flux_w_m2) if surface_k else None
