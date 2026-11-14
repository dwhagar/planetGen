# planetgen/generation/belt_planets.py

"""
The heaviest round members of a belt become planets (GEN.1510; Boss, 2026-11-09, 6:08 PM EST and his answers of 2026-11-10; the design is `reports/belt-bodies-as-planets-2026-11-09.md`)
=====================================================================================================================================================================================

A belt stores the bodies of its ledger as members (`AsteroidBelt.members`, GEN.1070), each with its mass, orbit, composition and a status (`body_status`: `minor body`, `dwarf planet`, `planet`). A body that is round has the
conditions of a planet of Class D at least, and Boss's rule is that it is made the way a planet is made, not a second way. The generator cannot make every round member of a belt (a belt stores up to 300, and a
planet takes a tenth of a second), so each belt promotes its `tuning.BELT_PLANETS_PER_BELT` heaviest round members (ten); the others stay members.

What a promoted member is:

* **A planet of its star** (`Planet`), made by the planet generator with the member's mass, orbit (eccentricity, plane, node, periapsis, place on the orbit) and shares of rock, metal, water and ice given and everything
  else drawn as for any planet: the air it holds (none, mostly), its interior, surface, spin, ring, field, moons, class and life. The planet draws from a stream of its own, keyed on the belt and the member's rank, so
  a promoted body never moves the draws of the system round it. It keeps the pair `belt_order` and `belt_rank`, the belt's place among its star's belts and the member's rank in it, which is how a body caught from a belt is
  named too (`irregular_moons`), so the member's designation and page still resolve.
* **In distance order among its star's planets.** The planets and belts of a star are one list in order of distance; the promoted bodies are put in by their semi-major axes, so the planets after them are numbered after them
  (their names, and the order they are stored in, move; `names.bodies.name_bodies` numbers them in the order they stand). They are exempt from the spacing of the slots (`validation.slot_planets`): they orbit inside the
  belt, on the orbits the belt gave them.
* **Out of the belt.** The member leaves the belt's members (its rank is not used again) and its mass leaves the belt's (`mass_kg`), so the belt and the planets made of it add up to the belt that was drawn. The belt's
  largest body (`largest_body_km`) and its other numbers are as they were drawn.

It runs once, when the planets, belts and ages of the system are final and before the irregular moons are given (`StarSystem._add_irregular_moons`), so a promoted planet can catch moons by the same rule as any planet, and
no body that a planet caught from the belt (status `moon`) is promoted. A belt that was read back from a store holds the members it was saved with, minus the promoted ones; promoting is only done by generation.
"""

from planetgen import tuning
from planetgen.generation import irregular_moons, validation
from planetgen.generation.planet import Planet
from planetgen.physics import body_status
from planetgen.util import draw, log

TRIES_PER_PLANET = 2
"""int: How many members a belt tries, at most, for each planet it may promote (`promote`). A stand-in search bound (it limits the work, it is not a physical quantity)."""

ROUND_STATUSES = (body_status.DWARF_PLANET, body_status.PLANET)
"""tuple: The statuses of a member that is round (`body_status`): the ones that can be promoted."""


def chosen(belt, limit=None):
    """
    The members of `belt` to promote: the `limit` (default `tuning.BELT_PLANETS_PER_BELT`) heaviest that are round, heaviest first. A member a planet caught (status `moon`) is a moon and is not one.
    """
    limit = tuning.BELT_PLANETS_PER_BELT if limit is None else limit
    round_ones = [member for member in belt.members if member["status"] in ROUND_STATUSES]
    return sorted(round_ones, key=lambda member: (-member["mass_kg"], member["rank"]))[:max(limit, 0)]


def planet_of(system, star, belt, order, member):
    """
    The planet the member `member` of `belt` (the `order`-th of the belts of its star) becomes: made by the planet generator from a stream of its own, with the member's mass, orbit and makeup.
    """
    stream = draw.Stream(f"belt-planet:{belt.distance!r}:{belt.lower_limit!r}:{star.mass!r}:{member['rank']}")
    with draw.bound(stream):
        return Planet(system.system_config, star, star.habitable_zone, member["distance_au"], mass=member["mass_kg"], belt_member=member, belt_order=order)


def promote(system, star, bodies):
    """
    Promotes the round members of every belt in `bodies` (one star's planets and belts, in orbital order): makes their planets, takes them out of their belts' members and masses and puts them in `bodies` by distance.
    A member whose planet the planet generator does not call round (it has no class: its own roundness ratio under 1, as for a member only just over the round mass, or its orbit reaches beyond the farthest stable orbit of the star, `validation.orbit_ceiling_au`, as a member on a long ellipse at the belt's edge can) stays a member, and the next heaviest is tried,
    to at most twice the limit of tries a belt, so every promoted body is a planet with a class.

    Returns:
        list: The new planets, in the order they were made (belt by belt, heaviest first).
    """
    made = []
    limit = tuning.BELT_PLANETS_PER_BELT
    ceiling = validation.orbit_ceiling_au(star)
    for order, belt in irregular_moons.belts_of(bodies).items():
        taken = 0
        for member in chosen(belt, limit * TRIES_PER_PLANET):
            if taken >= limit:
                break
            planet = planet_of(system, star, belt, order, member)
            if planet.planet_class is None or validation.reach_au(planet) > ceiling * (1 + validation.RELATIVE_TOLERANCE):
                continue
            taken += 1
            made.append(planet)
            belt.members.remove(member)
            belt.mass_kg = max(belt.mass_kg - member["mass_kg"], 0.0)
            log.choice("Belt member made a planet", f"rank {member['rank']} of the belt at {belt.distance:.4g} AU",
                       f"{member['mass_kg']:.3g} kg at {member['distance_au']:.4g} AU, {member['status']}, class {planet.planet_class}")
    if made:
        bodies[:] = sorted([*bodies, *made], key=lambda body: body.distance)
    return made
