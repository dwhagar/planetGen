# generation/rogue_topup_build.py

"""
The rogue planets of a sector's top-up, built (GEN.1120; docs/design/systems/region-object-store.md section 8.5, stage 3).

`physics/rogue_topup.plan` says how many rogues a sector needs, their masses, places, velocities and apparent origins. This builds each one the way a rogue
a system threw out is built (`rogue_yield`): through the one planet chain against a starless host (`starless_host`), with the mass it was given and a
donor star of the sector (the star whose system gave the mass, never a claimed ejector: the origin is marked "apparent" and no donor is linked). The body's
age, its time since the ejection and its speed are the plan's; the rest of its history is a synthetic ejection event of those numbers, with a scatterer
drawn at the planner's range (`rogue_topup.SCATTERER_AU_RANGE`) for the moons it kept.

    made = build(topup.rogues, hosts, system_config, "region-key:5", address, edge_pc)   # [ToppedRogue(rogue, planet, landing, note)]
"""

import collections
import math

from planetgen.galaxy import frame, galactic_orbit, spherical
from planetgen.generation import ejection_events, rogue_landing
from planetgen.generation.phenomena import starless_host
from planetgen.physics import constants, planet_architecture, rogue_moons, rogue_topup
from planetgen.physics.units import pc_to_ly
from planetgen.util import draw

ToppedRogue = collections.namedtuple("ToppedRogue", ("rogue", "planet", "landing", "note"))
"""One top-up rogue: the plan's `TopUpRogue`, the planet through the chain, where it is now (`rogue_landing.Landing`) and the sentence that says where it
seems to have come from."""

ORIGIN_LABEL = rogue_topup.ORIGIN


def note_of(rogue):
    """The rogue note of a top-up rogue: it seems to have come from another system, never a real one."""
    origin = rogue.apparent_origin
    return (f"Thrown out of a system about {origin.distance_pc:.0f} pc from here {origin.since_ejection_gy:.2f} Gy ago, at {origin.speed_km_s:.1f} km/s "
            f"(an {ORIGIN_LABEL} origin: no system is claimed).")


def _event(index, rogue, donor, formation_au, scatter_au, uid_text):
    """The synthetic `EjectionEvent` of a top-up rogue: only the mass, the formation distance and the moon limit are read by the chain."""
    mass_kg = rogue.mass_earth * constants.EARTH_MASS_TO_KG
    origin = rogue.apparent_origin
    jupiter_kg = constants.JUPITER_MASS_TO_KG
    hill_au = formation_au * (mass_kg / (3.0 * donor.mass)) ** (1.0 / 3.0)
    loss = planet_architecture.Loss("thrown", 1, mass_kg, mass_kg, formation_au, formation_au, 0.0, scatter_au, jupiter_kg, None, hill_au,
                                    bool(formation_au > getattr(donor, "snow_line_au", math.inf)), 1.0, 0.0, 0.0)
    stream = draw.Stream(f"topup-moons:{uid_text}:{index}")
    limit = starless_host.ejected_moon_limit_km(loss, mass_kg, donor.mass, stream.random())
    delay_gy = max(origin.age_gy - origin.since_ejection_gy, 0.0)
    return ejection_events.EjectionEvent(
        None, index, mass_kg, mass_kg, formation_au, formation_au, scatter_au, jupiter_kg, loss.ice_rich, delay_gy, origin.age_gy, origin.since_ejection_gy,
        (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), limit)


def _along(vector, axes):
    """A `(outward, along the rotation, up)` vector in the galaxy frame, given the three unit axes of a sector."""
    return tuple(sum(vector[k] * axes[k][i] for k in range(3)) for i in range(3))


def _landing(rogue, address, edge_pc, circular_speed_km_s):
    """Where the rogue is: its place in the sector (the plan's offset from the corner, along the sector's own axes) and its velocity through the Galaxy."""
    center = spherical.sector_position_pc(*address, edge_pc)
    axes = frame.sector_orientation(center)
    half = 0.5 * edge_pc
    offset = _along(tuple(c - half for c in rogue.position_pc), axes)
    position = tuple(c + o for c, o in zip(center, offset))
    u, v, w = rogue.velocity_km_s
    velocity = _along((u, v + circular_speed_km_s, w), frame.sector_orientation(position))
    return rogue_landing.Landing(position, velocity, tuple(address), rogue.apparent_origin.distance_pc)


def build(rogues, hosts, system_config, stream_key, address, edge_pc):
    """
    Every rogue of a top-up, built and placed.

    Args:
        rogues (iterable): The plan's `TopUpRogue`s.
        hosts (sequence): The sector's own stars that can donate a mass (`star.mass` kg, `age`, `snow_line_au`), at least one when `rogues` is not empty. A
            rogue's donor is the one the stream picks for it.
        system_config (SystemConfig): The shared settings.
        stream_key (str): A key unique to this sector's top-up (its public ID), which the draws are seeded from: the same sector and state give the same rogues.
        address (tuple): The `(shell, band, wedge)` of the sector.
        edge_pc (float): The grid's sector edge, pc.

    Returns:
        list: A `ToppedRogue` for each plan rogue, in order.
    """
    rogues = list(rogues)
    if rogues and not hosts:
        raise rogue_topup.RogueTopUpError("a top-up needs at least one host star to donate the mass")
    center = spherical.sector_position_pc(*address, edge_pc)
    circular_speed, _period = galactic_orbit.calculate_galactic_orbit(pc_to_ly(math.sqrt(sum(c * c for c in center))))
    low, high = rogue_topup.SCATTERER_AU_RANGE
    made = []
    for index, rogue in enumerate(rogues):
        stream = draw.Stream(f"topup-body:{stream_key}:{index}")
        donor = hosts[min(int(stream.random() * len(hosts)), len(hosts) - 1)]
        scatter_au = low * (high / low) ** stream.random()
        formation_au = low * (high / low) ** stream.random()
        event = _event(index, rogue, donor, formation_au, scatter_au, stream_key)
        host = starless_host.StarlessHost(donor, rogue.apparent_origin.age_gy)
        with draw.bound(f"topup-planet:{stream_key}:{index}"):
            planet = starless_host.build_rogue_body(system_config, host, event)
        where = _landing(rogue, address, edge_pc, circular_speed)
        rogue_landing.apply_motion(planet, where)
        made.append(ToppedRogue(rogue, planet, where, note_of(rogue)))
    return made
