# planetgen/generation/pair_orbit.py

"""
The orbit of a system's own pair of stars, read from its orbit tree (GEN.640 slice b).

A system with two or more stars has one pair, the two stars `primary_star` and `secondary_star` (the node both hang in, or the lowest
node above them). The pair's orbit is that node's: its axis, eccentricity, plane, phase and period are the node's, held in
`StarSystem.orbit_nodes`. `PairOrbit` is the read-only view of it that a reader wants (the numbers the old `WideBinaryPair` carried: the
period and speed, the closest and farthest approach, the secondary's place about the primary and each star's place about the pair's
center of mass), worked out on each read, so there is no second copy of the orbit to keep in step with the node.
"""

import dataclasses

from planetgen.physics import unit_registry as units
from planetgen.physics.orbits import circular_orbital_speed_kms, minimum_update_interval_years
from planetgen.physics.system_tree import node_relative_au

_km_to_au = units.converter("km", "au")


@dataclasses.dataclass(frozen=True)
class PairOrbit:
    """
    The pair's orbit as the pair node holds it. Lengths are AU, angles degrees, times years.

    Attributes:
        separation_au (float): The orbit's semi-major axis.
        eccentricity (float): Its eccentricity (0 for a close pair, whose node is a circle).
        period_years (float): The mutual orbital period.
        speed_kms (float): The circular-orbit speed at the axis (`orbits.circular_orbital_speed_kms`).
        periapsis_au (float): The closest approach, `separation_au * (1 - eccentricity)`.
        apoapsis_au (float): The farthest separation, `separation_au * (1 + eccentricity)`.
        inclination_deg (float): The plane's tilt.
        ascending_node_deg (float): The plane's longitude of the ascending node.
        phase_deg (float): The node's mean anomaly at its own epoch.
        min_update_interval_years (float): `orbits.minimum_update_interval_years` of the period.
        position_x_au/position_y_au/position_z_au (float): The secondary's place relative to the primary at the node's epoch (on its ellipse).
        primary_position_*_au/secondary_position_*_au (float): Each star's place about the pair's center of mass.
        secondary_mass_fraction (float): The lighter star's share of the pair's mass (the pair's primary is the heavier star).
    """
    separation_au: float
    eccentricity: float
    period_years: float
    speed_kms: float
    periapsis_au: float
    apoapsis_au: float
    inclination_deg: float
    ascending_node_deg: float
    phase_deg: float
    min_update_interval_years: float
    position_x_au: float
    position_y_au: float
    position_z_au: float
    primary_position_x_au: float
    primary_position_y_au: float
    primary_position_z_au: float
    secondary_position_x_au: float
    secondary_position_y_au: float
    secondary_position_z_au: float
    secondary_mass_fraction: float


def pair_orbit(node, primary_mass, secondary_mass, turned=False):
    """
    The `PairOrbit` of the pair node `node` (an `OrbitNode`) for stars of masses `primary_mass` and `secondary_mass` (any one unit).

    The node is the orbit of its slot-1 child about its slot-0 child and the pair is the second star about the first, so `turned` (the first
    star hangs under slot 1 of the node) turns the node's orbit half way round. A close pair's node holds a circle with the angle from the
    ascending node as its mean anomaly, which is the phase the pair is tracked by.
    """
    if turned:
        node = dataclasses.replace(node, arg_periapsis_deg=(node.arg_periapsis_deg + 180.0) % 360.0)
    separation_au = _km_to_au(node.separation_km)
    lighter, heavier = sorted((primary_mass, secondary_mass))
    fraction = lighter / (lighter + heavier)
    x, y, z = node_relative_au(node, node.epoch_unix)
    return PairOrbit(
        separation_au=separation_au, eccentricity=node.eccentricity, period_years=node.period_years,
        speed_kms=circular_orbital_speed_kms(separation_au, node.period_years),
        periapsis_au=separation_au * (1 - node.eccentricity), apoapsis_au=separation_au * (1 + node.eccentricity),
        inclination_deg=node.inclination_deg, ascending_node_deg=node.ascending_node_deg, phase_deg=node.mean_anomaly_deg,
        min_update_interval_years=minimum_update_interval_years(node.period_years),
        position_x_au=x, position_y_au=y, position_z_au=z,
        primary_position_x_au=-fraction * x, primary_position_y_au=-fraction * y, primary_position_z_au=-fraction * z,
        secondary_position_x_au=(1.0 - fraction) * x, secondary_position_y_au=(1.0 - fraction) * y, secondary_position_z_au=(1.0 - fraction) * z,
        secondary_mass_fraction=fraction)
