# generation/rogue_yield.py

"""
The rogue planets a generated system made (GEN.1111, part 7, docs/design/systems/region-object-store.md section 8): the ejection events of its planet
architecture, each built through the one planet chain against a starless host and put where it has drifted to.

A rogue exists only because a system threw it out, so the set of a donor's rogues is a function of the donor alone: its uid (the stream the events and
the bodies are drawn from), its record of losses (`StarSystem.ejected_losses`), its place in the galaxy and the galaxy's potential. Charting a region in
a different order, or charting a neighbour, never changes a donor's rogues. The save of the donor (`db/store.save_donor_rogues`) writes them, replacing
the ones the donor made before.

    made = rogues_of(donor_uid, system, donor_position_pc, edge_pc)     # [Rogue(event, landing, planet, note)]
"""

import collections

from planetgen.galaxy import galactic_orbit
from planetgen.generation import ejection_events, rogue_landing
from planetgen.generation.phenomena import starless_host
from planetgen.physics.units import pc_to_ly
from planetgen.util import draw

Rogue = collections.namedtuple("Rogue", ("event", "landing", "planet", "note"))
"""One rogue planet: the event that made it, where it is now, the planet through the chain and the sentence that says where it came from."""


def note_of(system_name, event):
    """The rogue note: the system it was thrown out of and how long ago, in the units of the site (Gy)."""
    return (f"Thrown out of the {system_name} system {event.since_ejection_gy:.2f} Gy ago, {event.delay_gy:.2f} Gy after that system formed, "
            f"and it has drifted {{separation_pc:.0f}} pc from there since.")


def rogues_of(donor_uid, system, donor_position_pc, edge_pc, system_name=None):
    """
    Every rogue planet `system` made, built and placed.

    Args:
        donor_uid (bytes): The donor's uid (12 bytes); the events and the bodies are drawn from streams seeded from it.
        system (StarSystem): The freshly generated donor (a system loaded from storage has no loss record and gives none).
        donor_position_pc (sequence): The donor's galaxy-frame position, pc.
        edge_pc (float): The sector edge of the galaxy's grid, pc.
        system_name (str, optional): For the note; `system.name` by default.

    Returns:
        list: A `Rogue` for every body that ended ejected, in event order.
    """
    pairs = system.ejected_losses()
    if not pairs:
        return []
    uid_text = donor_uid.hex() if isinstance(donor_uid, (bytes, bytearray)) else str(donor_uid)
    distance_ly = pc_to_ly((donor_position_pc[0] ** 2 + donor_position_pc[1] ** 2 + donor_position_pc[2] ** 2) ** 0.5)
    circular_speed, _period = galactic_orbit.calculate_galactic_orbit(distance_ly)
    name = system_name if system_name is not None else system.name
    stars = []
    for star, loss in pairs:
        if not stars or stars[-1][0] is not star:
            stars.append((star, []))
        stars[-1][1].append(loss)
    stream = draw.Stream(f"ejection-events:{uid_text}")
    made, first_index = [], 0
    for star, losses in stars:
        for event in ejection_events.events(uid_text, star, losses, rng=stream, first_index=first_index):
            where = rogue_landing.landing(donor_position_pc, event, edge_pc, circular_speed)
            host = starless_host.host_of(star, event)
            with draw.bound(f"rogue-body:{uid_text}:{event.index}"):
                planet = starless_host.build_rogue_body(system.system_config, host, event)
            rogue_landing.apply_motion(planet, where)
            made.append(Rogue(event, where, planet, note_of(name, event).format(separation_pc=where.separation_pc)))
        first_index += sum(loss.count for loss in losses)
    return made
