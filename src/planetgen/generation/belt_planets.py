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
  largest body (`largest_body_km`) is the largest member it still holds (GEN.1561: it was the drawn body's, a planet now); its other numbers are as they were drawn.
* **Every member over the oligarch mass is a planet** (GEN.1561). The ten is a cost limit, not a physical one: a member heavier than the isolation mass of the stretch the belt lies in (`oligarch_cap_kg`) would have
  cleared it, so it is promoted whatever the count, and the belt keeps its mass in the members below the cap. A real belt holds no body near that size (the main belt and the Kuiper belt
  with its largest objects hold 4.0e-4 and 2.0e-2 Earth masses, R3539), while a debris disc's belt is a mass of small planetesimals (R5248), so the cut leaves the mass and takes only the worlds out.

It runs once, when the planets, belts and ages of the system are final and before the irregular moons are given (`StarSystem._add_irregular_moons`), so a promoted planet can catch moons by the same rule as any planet, and
no body that a planet caught from the belt (status `moon`) is promoted. The first settle pass (`StarSystem._settle_to_final_age`) is over by then, so a promoted planet is given its own (`StarSystem._settle_promoted`, GEN.1668): the tides of its orbit, its envelope and end stages, then the tilt range and chaotic-tilt check (`settle_obliquity`) and the locked moons' Cassini states (`settle_moon_tilts`) against the modes of all its star's planets. The secular exchange and the Kozai pumping are not run again. A belt that was read back from a store holds the members it was saved with, minus the promoted ones; promoting is only done by generation.
"""

from planetgen import tuning
from planetgen.generation import irregular_moons, validation
from planetgen.generation.planet import Planet
from planetgen.physics import belt_origin, body_status, bulk_draw, orbits, unit_registry as units
from planetgen.util import draw, log

TRIES_PER_PLANET = 2
"""int: How many members a belt tries, at most, for each planet it may promote (`promote`). A stand-in search bound (it limits the work, it is not a physical quantity)."""

ROUND_STATUSES = (body_status.DWARF_PLANET, body_status.PLANET)
"""tuple: The statuses of a member that is round (`body_status`): the ones that can be promoted."""


def oligarch_cap_kg(star, belt):
    """
    The mass above which a member of `belt` is a planet and not a belt member (GEN.1561; Boss, 10:06 PM EST on 2026-11-12: "belt members above the oligarch mass go to the planet pass and the belt keeps
    its mass"): the isolation mass of the stretch the belt lies in, the same figure the belt's origin gives the most solids a planetless stretch can hold (`belt_origin.isolation_cap_earth`, GEN.1261). A body
    of that mass would have cleared the stretch and left no belt.
    """
    return units.convert(belt_origin.isolation_cap_earth(belt.lower_limit, belt.upper_limit, star.mass, star.snow_line_au, bulk_draw.metal_ratio_to_solar(star)), "Mearth", "kg")


def in_ring(belt, member):
    """
    True when the whole orbit of `member` lies between the limits of `belt`: its closest and farthest points, `a (1 - e)` and `a (1 + e)` AU, are inside `lower_limit` and `upper_limit`. A member on a long ellipse that leaves the ring crosses the planets
    and belts beside it, and the ring is what the system generator has cleared of them, so it stays a member of the belt (a population that crosses what is round it) and is not made a planet.
    """
    distance, eccentricity = member["distance_au"], member.get("eccentricity") or 0.0
    return distance * (1.0 - eccentricity) >= belt.lower_limit and distance * (1.0 + eccentricity) <= belt.upper_limit


def clear_of_planets(star, member, order, planets):
    """
    True when `member` stands `tuning.MUTUAL_HILL_RADII_SEPARATION` mutual Hill radii (Gladman 1993, `orbits.mutual_hill_radius_au`) from every planet of `planets` that did not come from the belt `order`, as the planets of the star stand from each other
    (`validation.min_distance_after_au`): a body of the ring keeps that room from the planets beside it, and one that is closer stays a member. The planets of its own belt are the ring and are not asked.
    """
    for planet in planets:
        if planet.body_type == 'a' or planet.belt_order == order:
            continue
        hill = orbits.mutual_hill_radius_au(member["mass_kg"], planet.mass, member["distance_au"], planet.distance, star.mass)
        if abs(member["distance_au"] - planet.distance) < tuning.MUTUAL_HILL_RADII_SEPARATION * hill:
            return False
    return True


def chosen(belt, limit=None):
    """
    The members of `belt` to promote: the `limit` (default `tuning.BELT_PLANETS_PER_BELT`) heaviest that are round and whose whole orbit lies inside the belt's ring, heaviest first. A member a planet caught (status `moon`) is a moon and is not one.
    """
    limit = tuning.BELT_PLANETS_PER_BELT if limit is None else limit
    round_ones = [member for member in belt.members if member["status"] in ROUND_STATUSES and in_ring(belt, member)]
    return sorted(round_ones, key=lambda member: (-member["mass_kg"], member["rank"]))[:max(limit, 0)]


def largest_left_km(belt):
    """The diameter in km of the largest member `belt` still holds (a member a planet caught, status `moon`, is not one), 0.0 for a belt with none left."""
    return max((2.0 * member["radius_km"] for member in belt.members if member["status"] != body_status.MOON), default=0.0)


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
    ceiling = validation.orbit_ceiling_au(star)
    for order, belt in irregular_moons.belts_of(bodies).items():
        taken = 0
        cap = oligarch_cap_kg(star, belt)
        # the heaviest `BELT_PLANETS_PER_BELT` round members, and every round member over the oligarch mass, whatever their number (GEN.1561)
        limit = 0 if tuning.BELT_PLANETS_PER_BELT <= 0 else max(tuning.BELT_PLANETS_PER_BELT, sum(1 for member in belt.members if member["status"] in ROUND_STATUSES and member["mass_kg"] > cap))
        for member in chosen(belt, limit * TRIES_PER_PLANET):
            if taken >= limit:
                break
            if not clear_of_planets(star, member, order, [*bodies, *made]):
                continue
            planet = planet_of(system, star, belt, order, member)
            if planet.planet_class is None or validation.reach_au(planet) > ceiling * (1 + validation.RELATIVE_TOLERANCE):
                continue
            taken += 1
            made.append(planet)
            belt.members.remove(member)
            belt.mass_kg = max(belt.mass_kg - member["mass_kg"], 0.0)
            log.choice("Belt member made a planet", f"rank {member['rank']} of the belt at {belt.distance:.4g} AU",
                       f"{member['mass_kg']:.3g} kg at {member['distance_au']:.4g} AU, {member['status']}, class {planet.planet_class}")
        if taken:
            belt.largest_body_km = largest_left_km(belt)       # GEN.1561: the belt's largest body is the largest it still holds
    if made:
        bodies[:] = sorted([*bodies, *made], key=lambda body: body.distance)
    return made
