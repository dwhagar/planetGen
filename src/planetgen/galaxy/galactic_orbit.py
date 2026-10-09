# planetgen/galaxy/galactic_orbit.py

"""
Galactic Orbit
==============

A star system's circular orbit around the galactic center: its speed,
period and phase, and the text every "Galactic Orbit" row shows.
"""

import math

from planetgen.physics import constants as physical_constants
from planetgen.physics.orbits import minimum_update_interval_years
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw
from planetgen.util.checks import finite_domain
from planetgen.util.format import format_period_years, format_speed_kms


# -inf is just "not positive" (the documented (0.0, 0.0)); +inf fails the result check.
@finite_domain(allow_inf=("distance_ly",))
def calculate_galactic_orbit(distance_ly):
    """
    Estimates a star system's circular orbital speed and orbital period
    around the galactic center, given its distance from it.

    Uses `physical_constants.GALACTIC_ROTATION_FLAT_VELOCITY_KMS`/
    `GALACTIC_ROTATION_CORE_RADIUS_PC`'s simple rotation-curve model (see
    that module's comment for the physical justification and calibration
    against Sol's own distance) rather than a Keplerian point-mass orbit
    around `physical_constants.MILKY_WAY_MASS` -- the latter would put a
    Sol-distance orbit at nearly 800 km/s, ~4x the real value, since most
    of the galaxy's mass isn't actually enclosed within that radius the way
    a naive point-mass calculation assumes. This function is deliberately
    independent of the orbiting body's own mass (unlike
    `calculate_hill_sphere`), matching real orbital mechanics at galactic
    scale: essentially every star's mass is negligible next to the
    galaxy's, so orbital speed at a given radius is the same for any star
    there, not a function of that star's own mass.

    Args:
        distance_ly (float): The star system's distance from the galactic
                             center, in light-years.

    Returns:
        tuple: `(orbital_speed_kms, orbital_period_gy)` -- circular orbital
              speed in km/s, and orbital period in billions of years (Gy),
              the same unit `Star.age`/`lifespan` already use. Both `0.0`
              for a system placed exactly at the galactic center (r = 0,
              where a circular orbit is degenerate).
    """
    if distance_ly <= 0:
        return 0.0, 0.0

    distance_pc = ly_to_pc(distance_ly)
    core_radius_pc = physical_constants.GALACTIC_ROTATION_CORE_RADIUS_PC
    orbital_speed_kms = (
        physical_constants.GALACTIC_ROTATION_FLAT_VELOCITY_KMS
        * distance_pc / math.sqrt(distance_pc ** 2 + core_radius_pc ** 2)
    )

    circumference_m = 2 * math.pi * distance_ly * physical_constants.LY_TO_M
    orbital_period_s = circumference_m / (orbital_speed_kms * physical_constants.KM_TO_M_FACTOR)
    orbital_period_gy = (orbital_period_s / physical_constants.SECONDS_PER_YEAR) / 1e9

    return orbital_speed_kms, orbital_period_gy


def generate_galactic_orbit_fields(galactic_center_dist_ly=None, galactic_orbital_phase_deg=None):
    """
    Generates the four `galactic_orbital_*` fields every gravitationally-
    bound-to-the-galaxy body this generator tracks needs: a `Star`
    (`Star.__init__`), a compact remnant (`compactRemnant.CompactRemnant.
    _finish_init`), and every standalone exotic phenomenon (`nebulaData.
    Nebula`, `supernovaRemnantData.SupernovaRemnant`, `roguePlanetData.
    RoguePlanet`/`InterstellarComet`, `asteroidFieldData.AsteroidField`) --
    a rogue planet or a comet passing through is unbound from any specific
    star, not from the galaxy itself, so it still orbits the galactic
    center on the same timescale a lone star does, via the exact same
    mass-independent `calculate_galactic_orbit` formula. Extracted here so
    every caller computes this identically instead of re-deriving it.

    Args:
        galactic_center_dist_ly (float, optional): This body's actual
            distance from the galactic center, in light-years -- see
            `Star.calculate_system_perimeter`'s docstring for the same
            parameter/fallback convention. `None` (the default) uses the
            fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY` constant.
        galactic_orbital_phase_deg (float, optional): This body's current
            angular position around its galactic orbit, in degrees -- see
            `Star.__init__`'s docstring for the same parameter. `None`
            (the default) rolls a fresh random value in `[0, 360)`.

    Returns:
        tuple: `(galactic_orbital_speed_kms, galactic_orbital_period_gy,
              galactic_orbital_phase_deg, galactic_min_update_interval_years)`.
    """
    if galactic_center_dist_ly is None:
        galactic_center_dist_ly = physical_constants.GALACTIC_CENTER_DISTANCE_LY
    speed_kms, period_gy = calculate_galactic_orbit(galactic_center_dist_ly)
    phase_deg = galactic_orbital_phase_deg if galactic_orbital_phase_deg is not None else draw.uniform(0, 360)
    min_update_interval_years = minimum_update_interval_years(period_gy * 1e9)
    return speed_kms, period_gy, phase_deg, min_update_interval_years


def format_galactic_orbit(speed_kms, period_gy):
    """
    Formats a `(galactic_orbital_speed_kms, galactic_orbital_period_gy)`
    pair as the display string every "Galactic Orbit" table row uses (e.g.
    `Star.get_table_properties`, `doubleStar.BinaryStarProxy.
    get_table_properties`, `compactRemnant.BlackHole`/`NeutronStar`, and
    every standalone exotic phenomenon) -- extracted so all of them render
    it identically instead of re-deriving the same f-string.

    Args:
        speed_kms (float): Circular orbital speed, km/s.
        period_gy (float): Orbital period, billions of years.

    Returns:
        str: e.g. `"206 km/s (236 My per orbit)"`, both on the shared
             ladders (`format_speed_kms`, `format_period_years`).
    """
    return f"{format_speed_kms(speed_kms)} ({format_period_years(period_gy * 1e9)} per orbit)"
