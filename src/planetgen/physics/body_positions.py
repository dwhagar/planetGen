# planetgen/physics/body_positions.py

"""
Where every body of a system is at any time (MAP.70), from the scene
`GET /api/systems/<id>/scene` returns (`web/maps/systemscene.py`).

`positions_at(scene, years)` moves each body `years` past the scene's
`epoch` along its stored orbit and returns `{ref: [x, y, z]}` in
kilometres in the scene's frame (the system barycenter, or the first star
for a wide binary, at the origin):

- A planet, moon or the second star of a pair goes round on a circular
  orbit: its stored phase plus 360 degrees for each period
  (`orbits.orbital_position_au`).
- A comet goes round on its Kepler orbit: an elliptical one's mean anomaly
  advances by 360 degrees a period (`kepler.comet_orbital_state`), a
  parabolic one's parabolic mean anomaly by `sqrt(mu / 2q^3)` a year.
- A close pair's stars sit either side of the barycenter in proportion to
  the other's mass (`secondary_mass_fraction`), a wide pair's second star
  orbits the first.

`static/orbitpositions.js` is the browser twin; `tests/test_js_unit.py`
checks the two against each other. NAV.6 plans courses on this one.
"""

import math

from planetgen.physics import constants
from planetgen.physics.kepler import comet_orbital_state, gravitational_parameter_au3_yr2
from planetgen.physics.orbits import orbital_position_au

AU_KM = constants.AU_TO_KM


def _relative_circular_km(orbit, years):
    """A circular orbit's position relative to what it goes round, in km."""
    period = orbit["period_years"]
    phase = orbit["phase_deg"] + (360.0 * years / period if period else 0.0)
    x, y, z = orbital_position_au(orbit["distance_km"] / AU_KM, orbit["inclination_deg"],
                                  orbit["ascending_node_deg"], phase % 360.0)
    return [x * AU_KM, y * AU_KM, z * AU_KM]


def _comet_relative_km(orbit, years):
    """A comet's position relative to its star, in km."""
    k = orbit["kepler"]
    q_au = k["perihelion_distance_km"] / AU_KM
    if orbit["type"] == "elliptical":
        anomaly = math.radians(k["mean_anomaly_deg"] + 360.0 * years / k["period_years"])
        state = comet_orbital_state(
            "elliptical", q_au, k["eccentricity"], k["inclination_deg"], k["arg_periapsis_deg"],
            k["ascending_node_deg"], k["primary_mass_solar"], mean_anomaly_rad=anomaly,
            orbital_period_years=k["period_years"])
    else:
        mu = gravitational_parameter_au3_yr2(k["primary_mass_solar"])
        anomaly = k["parabolic_mean_anomaly"] + math.sqrt(mu / (2 * q_au ** 3)) * years
        state = comet_orbital_state(
            "parabolic", q_au, 1.0, k["inclination_deg"], k["arg_periapsis_deg"],
            k["ascending_node_deg"], k["primary_mass_solar"], parabolic_mean_anomaly_value=anomaly)
    return [state["position_x_au"] * AU_KM, state["position_y_au"] * AU_KM, state["position_z_au"] * AU_KM]


def _add(a, b):
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def positions_at(scene, years):
    """
    Every star, planet, moon and comet's position `years` after the
    scene's epoch (negative for before), in km, as `{ref: [x, y, z]}`.
    Belts have no single point and are left out.
    """
    out = {}
    stars = scene["stars"]
    for star in stars:
        orbit = star["orbit"]
        if orbit is None:
            out[star["ref"]] = [0.0, 0.0, 0.0]
    for star in stars:
        orbit = star["orbit"]
        if orbit is None:
            continue
        relative = _relative_circular_km(orbit, years)
        if orbit["around"] == "barycenter":
            fraction = orbit["secondary_mass_fraction"]
            out[star["ref"]] = [c * (1.0 - fraction) for c in relative]
            for other in stars:
                if other["orbit"] is None:
                    out[other["ref"]] = [-c * fraction for c in relative]
        else:
            out[star["ref"]] = _add(out[orbit["around"]], relative)

    def origin(around):
        return [0.0, 0.0, 0.0] if around == "barycenter" else out[around]

    for planet in scene["planets"]:
        out[planet["ref"]] = _add(origin(planet["orbit"]["around"]), _relative_circular_km(planet["orbit"], years))
        for moon in planet["moons"]:
            out[moon["ref"]] = _add(out[planet["ref"]], _relative_circular_km(moon["orbit"], years))
    for comet in scene["comets"]:
        out[comet["ref"]] = _add(origin(comet["orbit"]["around"]), _comet_relative_km(comet["orbit"], years))
    return out
