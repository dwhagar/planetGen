# planetgen/web/maps/moonorbits.py

"""
The shape of a moon's orbit about its planet, for the System Map's moon system map (MAP.308).

A moon is stored as a state vector about its planet (`position_*_km`, `velocity_*_kms`, GEN.289). `moon_orbit` works
the osculating ellipse out of it (`physics.state_vectors.elements_from_state`, the planet's and the moon's masses for
the gravitational parameter) and samples it in the plane the map is drawn in (x, y, the top-down view), in kilometres
about the planet. A moon with no usable state (no velocity, no mass, or an orbit that is not a closed ellipse) falls
back to the circle its `distance_km` names, which is what the 3D view draws (`systemview3d.js`).

Pure functions of the rows: no database, no drawing. Tests: tests/test_moon_orbits.py.
"""

import math

from planetgen.physics import constants
from planetgen.physics.state_vectors import elements_from_state, state_from_elements

_G_KM3_PER_KG_S2 = constants.G * 1.0e-9
"""float: The gravitational constant in km^3 / (kg s^2), for positions in km and speeds in km/s."""

SAMPLES = 96
"""int: Points sampled round an orbit; the path closes on the first."""


def _vector(row, prefix, suffix):
    values = [row.get(f"{prefix}{axis}{suffix}") for axis in "xyz"]
    if any(v is None or not math.isfinite(v) for v in values):
        return None
    return tuple(float(v) for v in values)


def _circle(radius_km):
    points = [(radius_km * math.cos(2 * math.pi * n / SAMPLES), radius_km * math.sin(2 * math.pi * n / SAMPLES))
              for n in range(SAMPLES)]
    return {"points": points, "periapsis_km": radius_km, "apoapsis_km": radius_km, "eccentricity": 0.0, "circular": True}


def moon_orbit(planet, moon):
    """
    The orbit of `moon` about `planet`.

    Args:
        planet (dict): The planet's row (`mass_kg`).
        moon (dict): The moon's row (`position_x/y/z_km`, `velocity_x/y/z_kms`, `mass_kg`, `distance_km`).

    Returns:
        dict: `points` (x, y) in km about the planet round the orbit projected on the map's plane, `periapsis_km` and
        `apoapsis_km` (true distances from the planet, not the projection), `eccentricity` and `circular` (whether it is
        the fallback circle). The points are empty only for a moon with neither a state nor a distance.
    """
    position = _vector(moon, "position_", "_km")
    velocity = _vector(moon, "velocity_", "_kms")
    planet_mass = planet.get("mass_kg") or 0.0
    mu = _G_KM3_PER_KG_S2 * (planet_mass + (moon.get("mass_kg") or 0.0))
    if position is not None and velocity is not None and mu > 0.0:
        try:
            elements = elements_from_state(position, velocity, mu)
        except ValueError:
            elements = None
        if elements is not None and elements["kind"] == "elliptical":
            periapsis = elements["periapsis_distance"]
            eccentricity = elements["eccentricity"]
            points = []
            for n in range(SAMPLES):
                nu = 2 * math.pi * n / SAMPLES
                where, _speed = state_from_elements(mu, periapsis, eccentricity, elements["inclination"],
                                                    elements["ascending_node"], elements["arg_periapsis"], nu)
                points.append((where[0], where[1]))
            return {"points": points, "periapsis_km": periapsis, "apoapsis_km": periapsis * (1 + eccentricity) / (1 - eccentricity),
                    "eccentricity": eccentricity, "circular": False}
    radius = math.hypot(*position) if position is not None else 0.0
    radius = radius or moon.get("distance_km") or 0.0
    if radius > 0.0:
        return _circle(radius)
    return {"points": [], "periapsis_km": 0.0, "apoapsis_km": 0.0, "eccentricity": 0.0, "circular": True}
