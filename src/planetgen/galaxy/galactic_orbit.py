# planetgen/galaxy/galactic_orbit.py

"""
Galactic Orbit
==============

A star system's circular orbit around the galactic center: its speed,
period and phase.
"""

import math

from planetgen.galaxy import potential, scales
from planetgen.physics.orbits import minimum_update_interval_years
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw
from planetgen.util.checks import finite_domain
from planetgen.physics import unit_registry as units

_ly_to_m = units.converter("ly", "m")
_km_to_m = units.converter("km", "m")
_s_to_yr = units.converter("s", "yr")
_yr_to_gyr = units.converter("yr", "Gyr")


# -inf is just "not positive" (the documented (0.0, 0.0)); +inf fails the result check.
@finite_domain(allow_inf=("distance_ly",))
def calculate_galactic_orbit(distance_ly):
    """
    A star system's circular orbital speed and orbital period around the galactic center, given its distance from it, in the potential of the galaxy the
    running unit belongs to (`galaxy.potential.current`: the Hernquist spheres, Miyamoto-Nagai disks, NFW halo and central black hole its density model and
    mass chain give, GEN.338; the fallback seed's galaxy outside any unit). The speed is the potential's own, `sqrt(R |a_R|)` in the plane, so it differs from
    galaxy to galaxy and with radius and holds no Milky Way value (GEN.339: the 220 km/s curve with a 3 kpc core this used to give every galaxy is gone). This
    function is independent of the orbiting body's own mass, as essentially every star's mass is negligible next to the galaxy's.

    Args:
        distance_ly (float): The star system's distance from the galactic center, in light-years.

    Returns:
        tuple: `(orbital_speed_kms, orbital_period_gy)` -- circular orbital speed in km/s, and orbital period in billions of years (Gy), the same unit
              `Star.age`/`lifespan` already use. Both `0.0` for a system placed exactly at the galactic center (r = 0, where a circular orbit is degenerate).
    """
    if distance_ly <= 0:
        return 0.0, 0.0

    distance_pc = ly_to_pc(distance_ly)
    orbital_speed_kms = potential.current_circular_speed_kms(distance_pc)
    if orbital_speed_kms <= 0.0:
        return 0.0, 0.0

    circumference_m = _ly_to_m(2 * math.pi * distance_ly)
    orbital_period_s = circumference_m / _km_to_m(orbital_speed_kms)
    orbital_period_gy = _yr_to_gyr(_s_to_yr(orbital_period_s))

    return orbital_speed_kms, orbital_period_gy


def generate_galactic_orbit_fields(galactic_center_dist_ly=None, galactic_orbital_phase_deg=None):
    """
    Generates the four `galactic_orbital_*` fields every gravitationally-
    bound-to-the-galaxy body this generator tracks needs: a `Star`
    (`Star.__init__`), a compact remnant (`phenomena.compact_remnant.CompactRemnant.
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
            parameter/fallback convention. `None` (the default) uses the radius that holds half of the galaxy's disc stars (`galaxy.scales`).
        galactic_orbital_phase_deg (float, optional): This body's current
            angular position around its galactic orbit, in degrees -- see
            `Star.__init__`'s docstring for the same parameter. `None`
            (the default) rolls a fresh random value in `[0, 360)`.

    Returns:
        tuple: `(galactic_orbital_speed_kms, galactic_orbital_period_gy,
              galactic_orbital_phase_deg, galactic_min_update_interval_years)`.
    """
    if galactic_center_dist_ly is None:
        galactic_center_dist_ly = scales.default_galactocentric_ly()
    speed_kms, period_gy = calculate_galactic_orbit(galactic_center_dist_ly)
    phase_deg = galactic_orbital_phase_deg if galactic_orbital_phase_deg is not None else draw.uniform(0, 360)
    min_update_interval_years = minimum_update_interval_years(period_gy * 1e9)
    return speed_kms, period_gy, phase_deg, min_update_interval_years

