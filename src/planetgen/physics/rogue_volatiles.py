# planetgen/physics/rogue_volatiles.py

"""
What air and water a thrown-out planet carries (GEN.1112, RP3)
==============================================================

A rogue rock drew its hydrogen envelope from a coin (`ROGUE_HYDROGEN_ENVELOPE_CHANCE` 0.1 and 0.5), its water from another (`ROGUE_WATER_RICH_CHANCE`), its frozen air from a mass cut
(`ROGUE_FROZEN_ATMOSPHERE_MIN_MASS_EARTH`) and its gas-giant class from a mass threshold (`ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER`), none of them tied to the system the planet came from
(`rogue-planet-review.md` rows 5, 6, 8 and 9). A planet that was ejected carries what its formation gave it less what its star took before it left, so this module asks the
code that already decides those for in-system planets, with the ejected planet's own record:

* **The envelope** the planet kept (`envelope_at_ejection`) is `envelope_escape.envelope_left_earth` at the planet's orbit, for the time it spent lit by its star (the ejection delay of `rogue_age`).
  That loss law starts at 1 Gyr (`envelope_escape.ONSET_GY`; what is lost earlier is in the starting envelope), so a planet thrown out in its first Gyr leaves with its starting envelope whole. After the
  ejection there is no star, so nothing more is lost to XUV.
* **The envelope's base pressure** (`envelope_base_pressure_bar`) is the weight of the envelope on the planet's surface, `G M M_env / (4 pi R^4)`: 1,100 bar for an Earth-mass planet with a tenth of a percent
  of its mass in hydrogen, inside Stevenson's 10^2 to 10^3 bar for rogue Earths. It replaces the log-uniform draw.
* **A gas giant** (`is_gas_giant`) is a planet whose envelope makes its radius at least `tuning.GIANT_MIN_RADIUS_KM`, the radius the generator already calls a giant (GEN.681,
  `planet_kind.kind_of_given`), from the radius law of Lopez and Fortney 2014 (`bulk_structure.envelope_radius_earth`) with no stellar heating. The mass threshold becomes a result: at 3 Gyr the
  15,000 km is reached by 5 Earth masses with 5 percent of their mass in hydrogen, by 10 with 3 percent and by 17 with 1.5 percent (`smallest_giant_envelope_fraction`); above about 30 a rock core alone is that big.
* **The water** (`water_mass_fraction`) is `delivered_water.sample_fraction` at the distance the planet formed at in units of its star's snow line: the law the in-system rocky planets take their water from.
* **An air frozen on the surface** (`had_air`) is a planet that, at the orbit it left, held an air by the shoreline test of `atmosphere_retention` (a verdict of `HELD_AIR_VERDICTS`).

Everything is a pure function of the ejection record's numbers; the random numbers come from the caller. Nothing calls it yet (GEN.1111 creates the rogue rows). Stand-ins are in
`docs/design/research/gap-filling-plan.md` section 6.111.
"""

import math

from planetgen import tuning
from planetgen.physics import atmosphere_retention, bulk_structure, constants, delivered_water, envelope_escape, hydrosphere, unit_registry as units
from planetgen.util import draw
from planetgen.util.checks import finite_domain

_mearth_to_kg = units.converter("Mearth", "kg")

HELD_AIR_VERDICTS = ("dense", "moderate", "tenuous")
"""tuple: The shoreline verdicts (`atmosphere_retention.verdict_of`) that count as a body having held an air of its own; `trace` and `none` are bare rock. Source tier T6, tweak level L0."""
NO_STELLAR_FLUX_EARTH = 1.0e-6
"""float: The stellar flux on a rogue's envelope, Earth units, small enough that the radius law of Lopez and Fortney sees no heating (numerical; the law is fitted to irradiated planets)."""
ENVELOPE_LAW_MASS_RANGE_EARTH = (1.0, 30.0)
"""tuple: The planet masses the radius law of Lopez and Fortney was fitted to, Earth masses (`envelope_escape` uses it for 1 to 30); `envelope_radius_km` holds a mass outside it at the nearest edge. Source tier T1, tweak level L0."""
_BAR_PA = 1.0e5


class RogueVolatilesError(ValueError):
    """An input the model cannot use."""


@finite_domain()
def envelope_at_ejection(core_mass_earth, starting_envelope_fraction, insolation_earth, host_mass_solar, luminosity_solar, delay_gy):
    """
    The envelope a planet leaves its system with, Earth masses, after the XUV of its star has worked on it for the ejection delay.

    Args:
        core_mass_earth (float): The planet's mass (the rows of `envelope_escape` pass the planet's mass here).
        starting_envelope_fraction (float): The envelope as a fraction of the planet's mass the formation gave it (0 for a bare rock).
        insolation_earth (float): The bolometric flux at its orbit, Earth units.
        host_mass_solar, luminosity_solar (float): The star's.
        delay_gy (float): How long after formation it was thrown out (`rogue_age`).
    """
    if not core_mass_earth > 0.0 or not insolation_earth > 0.0 or not delay_gy >= 0.0 or not 0.0 <= starting_envelope_fraction < 1.0:
        raise RogueVolatilesError("the mass and the insolation must be positive, the delay not negative, the envelope fraction from 0 to under 1")
    return envelope_escape.envelope_left_earth(core_mass_earth, starting_envelope_fraction, insolation_earth, host_mass_solar, luminosity_solar, delay_gy)


@finite_domain()
def envelope_base_pressure_bar(mass_earth, radius_earth, envelope_earth):
    """The pressure at the bottom of a hydrogen envelope, bar: its weight on the surface, `G M M_env / (4 pi R^4)`."""
    if not mass_earth > 0.0 or not radius_earth > 0.0 or not envelope_earth >= 0.0:
        raise RogueVolatilesError("the mass and radius must be positive and the envelope not negative")
    mass = _mearth_to_kg(mass_earth)
    radius = radius_earth * constants.EARTH_RADIUS_KM * 1.0e3
    return constants.G * mass * envelope_earth * constants.EARTH_MASS_TO_KG / (4.0 * math.pi * radius ** 4) / _BAR_PA


@finite_domain()
def envelope_radius_km(mass_earth, envelope_fraction, age_gy):
    """The radius of a planet of `mass_earth` Earth masses holding `envelope_fraction` of its mass as hydrogen and helium, with no star to heat it, km (Lopez and Fortney 2014, for 1 to 30 Earth masses; a mass outside is held at the nearest edge)."""
    if not mass_earth > 0.0 or not 0.0 <= envelope_fraction < 1.0 or not age_gy > 0.0:
        raise RogueVolatilesError("the mass and age must be positive and the envelope fraction from 0 to under 1")
    mass = min(max(mass_earth, ENVELOPE_LAW_MASS_RANGE_EARTH[0]), ENVELOPE_LAW_MASS_RANGE_EARTH[1])
    return bulk_structure.envelope_radius_earth(mass, envelope_fraction, NO_STELLAR_FLUX_EARTH, age_gy) * constants.EARTH_RADIUS_KM


def is_gas_giant(mass_earth, envelope_fraction, age_gy):
    """Whether the planet's envelope makes it a giant: its radius is at least `tuning.GIANT_MIN_RADIUS_KM`."""
    return envelope_radius_km(mass_earth, envelope_fraction, age_gy) >= tuning.GIANT_MIN_RADIUS_KM


def smallest_giant_envelope_fraction(mass_earth, age_gy, tolerance=1.0e-4):
    """The envelope fraction at which a planet of `mass_earth` Earth masses becomes a giant, or `None` when even half its mass in hydrogen leaves it smaller (not reached inside the law's mass range)."""
    low, high = 0.0, 0.5
    if not is_gas_giant(mass_earth, high, age_gy):
        return None
    while high - low > tolerance:
        mid = 0.5 * (low + high)
        low, high = (low, mid) if is_gas_giant(mass_earth, mid, age_gy) else (mid, high)
    return high


def water_mass_fraction(mass_earth, formation_au, snow_line_au, rng=draw):
    """
    The water share of a rocky planet's mass, from where it formed: `delivered_water.sample_fraction` at `formation_au / snow_line_au` snow lines. Two draws from `rng`.

    Args:
        mass_earth (float): The planet's mass.
        formation_au (float): The distance it formed at.
        snow_line_au (float): Its star's snow line then.
    """
    if not mass_earth > 0.0 or not formation_au > 0.0 or not snow_line_au > 0.0:
        raise RogueVolatilesError("the mass and both distances must be positive")
    return delivered_water.sample_fraction(mass_earth, formation_au / snow_line_au, rng)


def is_water_rich(water_fraction):
    """Whether a water share makes an icy body (`hydrosphere.ICY_MIN_WATER_FRACTION`)."""
    return water_fraction >= hydrosphere.ICY_MIN_WATER_FRACTION


def had_air(mass_earth, radius_earth, insolation_earth, **kwargs):
    """
    Whether the planet held an air of its own at the orbit it left: its shoreline verdict (`atmosphere_retention.retention`) is one of `HELD_AIR_VERDICTS`.

    Args:
        mass_earth, radius_earth (float): The planet's.
        insolation_earth (float): The flux at the orbit it left, Earth units.
        kwargs: The rest of `atmosphere_retention.retention` (the star's XUV index and flux).
    """
    if not mass_earth > 0.0 or not radius_earth > 0.0 or not insolation_earth > 0.0:
        raise RogueVolatilesError("the mass, radius and insolation must be positive")
    return atmosphere_retention.retention(mass_earth, radius_earth, insolation_earth, **kwargs).verdict in HELD_AIR_VERDICTS
