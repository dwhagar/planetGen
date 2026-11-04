# generation/rogue_landing.py

"""
Where a rogue planet is now (GEN.1111, part 4, docs/design/systems/region-object-store.md section 8): the donor's place in the galaxy, the kick and the
time since the ejection give the rogue's position and velocity and so the sector it is stored in.

The velocity of the body relative to its donor is the planet's orbital velocity at the ejection plus the scattering kick (`ejection_events`). The donor
goes round the galaxy on a circular orbit (`galactic_orbit`), so the rogue's motion from it is the epicyclic drift of `rogue_drift`: the offset along the
orbit is an angle round the galactic centre (`ring_drift`), which is how a rogue that has gone far lands on the other side of the galaxy's rotation
rather than on a straight line. The axes are the galaxy frame's own (`galaxy/frame.py`): outward from the axis, along the rotation (counterclockwise from
galactic north) and up. The epicyclic constants are those of `rogue_drift.SOLAR_NEIGHBOURHOOD` (a stand-in for every donor's radius, registered with
its rows in `docs/design/research/gap-filling-plan.md`).
"""

import collections
import math

from planetgen.galaxy import frame, galactic_orbit, spherical
from planetgen.physics import rogue_drift, rogue_kick
from planetgen.physics.units import pc_to_ly

Landing = collections.namedtuple("Landing", ("position_pc", "velocity_km_s", "address", "separation_pc"))
"""Where a rogue is: its galaxy-frame position (pc) and velocity (km/s), the `(shell, band, wedge)` of the sector that holds it, and how far it is from its donor (pc)."""


def _axes(position):
    return frame.sector_orientation(position)


def _along(vector, axes):
    return tuple(sum(vector[k] * axes[k][i] for k in range(3)) for i in range(3))


def landing(donor_position_pc, event, edge_pc, circular_speed_km_s, galaxy=rogue_drift.SOLAR_NEIGHBOURHOOD):
    """
    The landing of the rogue an ejection event made.

    Args:
        donor_position_pc (sequence): The donor's galaxy-frame position, pc.
        event (ejection_events.EjectionEvent): What was thrown out (`kick_km_s`, `orbital_velocity_km_s`, `since_ejection_gy`).
        edge_pc (float): The galaxy's sector edge, pc (`sectors` shape).
        circular_speed_km_s (float): The donor's circular speed round the galaxy (`galactic_orbit.calculate_galactic_orbit`), which the rogue's own
            velocity is measured against; zero leaves the velocity relative to the local standard of rest.
        galaxy: The epicyclic constants.

    Returns:
        Landing: The position, velocity, sector address and distance from the donor.
    """
    x, y, z = donor_position_pc
    radius = math.hypot(x, y)
    kick = tuple(k + o for k, o in zip(event.kick_km_s, event.orbital_velocity_km_s))
    time_myr = event.since_ejection_gy * 1000.0
    if radius < 1e-9:      # on the galaxy's axis the epicyclic frame has no direction: the offset is taken in the galaxy axes
        here = rogue_drift.drift(kick, time_myr, galaxy)
        position = (x + here.x_pc, y + here.y_pc, z + here.z_pc)
        axes = _axes(position)
    else:
        ring = rogue_drift.ring_drift(kick, time_myr, radius, galaxy)
        new_radius = radius + ring.radial_pc
        angle = math.atan2(y, x) + ring.azimuth_rad
        position = (new_radius * math.cos(angle), new_radius * math.sin(angle), z + ring.vertical_pc)
        axes = _axes(position)
        here = rogue_drift.drift(kick, time_myr, galaxy)
    relative = rogue_kick.rogue_velocity((0.0, 0.0, 0.0), kick, time_myr, galaxy)   # outward, along the rotation, up, against the circular speed
    local = (relative[0], relative[1] + circular_speed_km_s, relative[2])
    velocity = _along(local, axes)
    separation = math.dist(position, donor_position_pc)
    return Landing(position, velocity, spherical.sector_address_at(position, edge_pc), separation)


def placement_of(landing_):
    """The `placement` dict the store takes (`center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc`) for a `Landing`."""
    x, y, z = landing_.position_pc
    return {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "galactic_radius_pc": math.sqrt(x * x + y * y + z * z)}


def apply_motion(planet, landing_):
    """
    Gives a rogue planet the motion the store keeps: its galaxy-frame velocity (`ejecta_velocity_kms`, which `store._free_vector` stores) and the
    galactic orbit fields of a free object (`galactic_orbit.generate_galactic_orbit_fields` at its distance from the centre, the phase its own
    azimuth, so no draw is made).
    """
    x, y, _z = landing_.position_pc
    distance_ly = pc_to_ly(placement_of(landing_)["galactic_radius_pc"])
    speed, period, phase, interval = galactic_orbit.generate_galactic_orbit_fields(distance_ly, math.degrees(math.atan2(y, x)) % 360.0)
    planet.ejecta_velocity_kms = tuple(landing_.velocity_km_s)
    planet.galactic_orbital_speed_kms = speed
    planet.galactic_orbital_period_gy = period
    planet.galactic_orbital_phase_deg = phase
    planet.galactic_min_update_interval_years = interval
    return planet
