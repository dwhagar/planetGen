# generation/ejection_events.py

"""
The ejection events of a donor system (GEN.1111, part 2, docs/design/systems/region-object-store.md section 8.1).

A rogue planet exists only because a generated system threw it out. The planet architecture (GEN.728, GEN.1119) records what it threw out
(`StarSystem.ejected_losses`, each a `planet_architecture.Loss` with the share that ends ejected); this module turns each such loss into the events
a rogue is built from: whether the body ended ejected, when, how fast and in which direction, how old it is now, and how wide a moon system it
kept. Every random number comes from the stream the caller passes (the donor's own, `draw` by default) in a fixed order, so the set of rogues of a
donor is a pure function of its seed and its record, the same whatever else is charted.

The velocity of the body relative to its donor is the planet's orbital velocity at the ejection plus the kick of the scattering (Boss, 5:29 PM EDT
2026-10-30: isotropic, built from the origin planet's actual motion, not assumed perpendicular to the galactic motion). A planet's orbital plane and
phase are random against the galaxy, so the direction of its orbital velocity is uniform on the sphere; its size is the circular speed at its
distance. The star's own velocity through the galaxy is added by the caller (`rogue_kick.rogue_velocity`).
"""

import collections

from planetgen.generation.phenomena import starless_host
from planetgen.physics import rogue_age, rogue_kick
from planetgen.util import draw

EjectionEvent = collections.namedtuple("EjectionEvent", (
    "donor_uid", "index", "mass_kg", "core_kg", "formation_au", "distance_au", "perturber_au", "perturber_mass_kg", "ice_rich",
    "delay_gy", "age_gy", "since_ejection_gy", "kick_km_s", "orbital_velocity_km_s", "moon_limit_km"))
"""One body a donor threw out: the donor's uid and the event's index in the donor (the key a rogue row carries, so saving a donor again replaces exactly
the rogues it made), the body's mass and core mass, where it formed and where it was thrown from, the scatterer's orbit and mass, whether it formed beyond
the snow line, the ejection delay after formation, the age now (the donor's) and the time since the ejection (all Gy), the kick of the scattering and the
planet's orbital velocity then (both `(outward, along the rotation, up)` in km/s, relative to the donor), and the outer edge of the moon system it kept (km)."""


def events(donor_uid, donor, losses, rng=draw):
    """
    The events of the bodies a donor threw out.

    Args:
        donor_uid: The donor system's uid (any value that is the same each time the donor is built).
        donor: The star the bodies were thrown from (`StarSystem.star`, or the star `StarSystem.ejected_losses` names).
        losses (iterable): `planet_architecture.Loss` records of kind `thrown` (`StarSystem.ejected_losses()` gives them with their star).
        rng: Anything with `random`; `draw` by default.

    Returns:
        list: An `EjectionEvent` for each body that ended ejected, in the order of `losses`; the index counts every body of every loss, ejected or
            not, so it stays the same if the shares change a little.
    """
    found = []
    index = 0
    for loss in losses:
        for _body in range(loss.count):
            mass_kg = loss.mass_kg / loss.count
            core_kg = loss.core_kg / loss.count
            u_fate, u_delay, u_speed, u_cos, u_phi, u_orbit_cos, u_orbit_phi, u_moons = (rng.random() for _ in range(8))
            event_index, index = index, index + 1
            if u_fate >= loss.ejected:
                continue
            delay_gy = rogue_age.ejection_delay_gy(u_delay, donor.age)
            ages = rogue_age.from_record(donor.age, delay_gy)
            speed = rogue_kick.ejection_speed_km_s(u_speed, loss.perturber_mass_kg, loss.perturber_au, mass_kg)
            kick = rogue_kick.kick(speed, rogue_kick.isotropic_direction(u_cos, u_phi))
            circular = rogue_kick.circular_speed_km_s(donor.mass, loss.distance_au)
            orbital = rogue_kick.kick(circular, rogue_kick.isotropic_direction(u_orbit_cos, u_orbit_phi))
            limit = starless_host.ejected_moon_limit_km(loss, mass_kg, donor.mass, u_moons)
            found.append(EjectionEvent(donor_uid, event_index, mass_kg, core_kg, loss.formation_au, loss.distance_au, loss.perturber_au,
                                       loss.perturber_mass_kg, loss.ice_rich, delay_gy, ages.age_gy, ages.since_ejection_gy, tuple(kick), tuple(orbital), limit))
    return found
