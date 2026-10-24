# planetgen/physics/body_state.py

"""
A planet's or moon's state vector and how it moves (GEN.289).

A planet or moon is stored as the state of a point: its position and velocity relative to what it goes round (the
`position_*_km` and `velocity_*_kms` columns), its mass, the epoch the state is true at, and its parent. There is no
stored phase. Motion is `kepler.universal_step` on that state, so a body's orbit is whatever the state says it is:
circular when generated, and free to change when a pull from another body changes the vector (GEN.109).

`distance_km` and `period_years` stay beside the vector as the orbit's semi-major axis and period, and the pull of the
body it goes round is `mu = (2 pi / period)^2 distance^3`, the same Kepler's-third-law figure the generator worked its
period out with (`Planet.orbit_mu_au3_per_year2`), so a state advanced here and a circular phase agree to rounding.

The browser twin is `static/orbitpositions.js`; `tests/test_js_unit.py` holds the two together.
"""

import math

from planetgen.physics import unit_registry as units
from planetgen.physics.kepler import universal_step


_yr_to_s = units.converter("yr", "s")


def mu_km3_s2(distance_km, period_years):
    """`G M` of what a body goes round, km^3/s^2, from its orbit's semi-major axis and period (Kepler's third law)."""
    omega = 2.0 * math.pi / _yr_to_s(period_years)
    return omega * omega * distance_km ** 3


def argument_of_latitude_deg(position, inclination_deg, ascending_node_deg):
    """
    How far along its orbit plane (measured from the ascending node) a position is, in [0, 360). The plane is the
    one `orbital_position_au` puts a body in; a position off the plane is read by its projection onto it.
    """
    i = math.radians(inclination_deg)
    node = math.radians(ascending_node_deg)
    along = (math.cos(node), math.sin(node), 0.0)
    across = (-math.sin(node) * math.cos(i), math.cos(node) * math.cos(i), math.sin(i))
    x = sum(a * b for a, b in zip(position, along))
    y = sum(a * b for a, b in zip(position, across))
    return math.degrees(math.atan2(y, x)) % 360.0


def advance(position_km, velocity_kms, distance_km, period_years, seconds):
    """
    The state `seconds` later (negative for earlier), by `kepler.universal_step`.

    Returns:
        tuple: `(position_km, velocity_kms)`.

    Raises:
        ValueError: As `universal_step`.
    """
    return universal_step(position_km, velocity_kms, mu_km3_s2(distance_km, period_years), seconds)
