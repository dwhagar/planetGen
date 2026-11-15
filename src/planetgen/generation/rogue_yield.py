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
from planetgen.physics import rogue_pairs
from planetgen.physics.constants import JUPITER_MASS_TO_KG
from planetgen.physics.units import pc_to_ly
from planetgen.util import detmath, draw

Pair = collections.namedtuple("Pair", ("partner_index", "partner_mass_kg", "separation_au", "period_yr"))
"""The companion a rogue left its donor with (GEN.332, slice D): the event index of the other body (the link both rows carry), its mass, and the circular orbit of the two (AU, years)."""

Rogue = collections.namedtuple("Rogue", ("event", "landing", "planet", "note", "pair"), defaults=(None,))
"""One rogue planet: the event that made it, where it is now, the planet through the chain, the sentence that says where it came from and its `Pair` (`None` for a rogue that left alone)."""


def note_of(system_name, event):
    """The rogue note: the system it was thrown out of (or pulled off by its birth cluster) and how long ago, in the units of the site (Gy)."""
    if event.cause == "cluster":
        return (f"Unbound by its birth cluster from the {system_name} system {event.since_ejection_gy:.2f} Gy ago, {event.delay_gy:.2f} Gy after that system formed, "
                f"and it has drifted {{separation_pc:.0f}} pc from there since.")
    return (f"Thrown out of the {system_name} system {event.since_ejection_gy:.2f} Gy ago, {event.delay_gy:.2f} Gy after that system formed, "
            f"and it has drifted {{separation_pc:.0f}} pc from there since.")


def _pair_note(pair):
    """The sentence a bound rogue's note ends with: its companion, how far and how fast round."""
    return (f" It left bound to a companion of {pair.partner_mass_kg / JUPITER_MASS_TO_KG:.1f} Jupiter masses, {pair.separation_au:.0f} AU away, "
            f"which the two go round in {pair.period_yr:,.0f} years.")


def bind_pairs(uid_text, found):
    """
    Decides which ejected bodies of a donor left as a pair (GEN.332, slice D; `physics/rogue_pairs.py`).

    Each body not yet taken is a system of free planets that is a bound pair with the chance `rogue_pairs.bound_share` of its age since the ejection; its companion is the next body
    the donor ejected from the same star by the same cause. The companion left at the same moment and with the same velocity, so it takes its partner's delay, age and kick, and
    both land together. The draws come from the donor's own stream (two for every body not yet taken, in event order), so the pairs of a donor are the same whatever else is charted.

    Args:
        uid_text (str): The donor's uid as text.
        found (list): `(star, EjectionEvent)` in event order.

    Returns:
        tuple: The events (the companions changed) and `{event index: Pair}` for both members of every pair.
    """
    events = [event for _star, event in found]
    pairs, taken = {}, set()
    stream = draw.Stream(f"rogue-pairs:{uid_text}")
    for at, (star, first) in enumerate(found):
        if first.index in taken:
            continue
        u_bound, u_separation = stream.random(), stream.random()
        if at + 1 >= len(found) or u_bound >= rogue_pairs.bound_share(first.since_ejection_gy):
            continue
        next_star, second = found[at + 1]
        if next_star is not star or second.cause != first.cause:
            continue
        separation = rogue_pairs.separation_au(u_separation, first.mass_kg + second.mass_kg)
        if separation is None:
            continue
        period = rogue_pairs.period_yr(separation, first.mass_kg, second.mass_kg)
        events[at + 1] = second._replace(delay_gy=first.delay_gy, age_gy=first.age_gy, since_ejection_gy=first.since_ejection_gy, kick_km_s=first.kick_km_s,
                                         orbital_velocity_km_s=first.orbital_velocity_km_s)
        pairs[first.index] = Pair(second.index, second.mass_kg, separation, period)
        pairs[second.index] = Pair(first.index, first.mass_kg, separation, period)
        taken.update((first.index, second.index))
    return events, pairs


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
    distance_ly = pc_to_ly(detmath.hypot(*donor_position_pc))
    circular_speed, _period = galactic_orbit.calculate_galactic_orbit(distance_ly)
    name = system_name if system_name is not None else system.name
    stars = []
    for star, loss in pairs:
        if not stars or stars[-1][0] is not star:
            stars.append((star, []))
        stars[-1][1].append(loss)
    stream = draw.Stream(f"ejection-events:{uid_text}")
    found, first_index = [], 0
    for star, losses in stars:
        found.extend((star, event) for event in ejection_events.events(uid_text, star, losses, rng=stream, first_index=first_index))
        first_index += sum(loss.count for loss in losses)
    events, pairs = bind_pairs(uid_text, found)
    made = []
    for (star, _first), event in zip(found, events):
        where = rogue_landing.landing(donor_position_pc, event, edge_pc, circular_speed)
        host = starless_host.host_of(star, event)
        with draw.bound(f"rogue-body:{uid_text}:{event.index}"):
            planet = starless_host.build_rogue_body(system.system_config, host, event)
        rogue_landing.apply_motion(planet, where)
        pair = pairs.get(event.index)
        note = note_of(name, event) + (_pair_note(pair) if pair is not None else "")
        made.append(Rogue(event, where, planet, note.format(separation_pc=where.separation_pc), pair))
    return made
