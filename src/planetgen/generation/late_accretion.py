# planetgen/generation/late_accretion.py

"""
The air reads the system's late accretion (GEN.1493 slice 1, docs/design/composition/late-accretion-from-the-system.md)
=========================================================================================================================

A rocky planet's air is made when the planet is, with one late-accretion number drawn from a law (`impact_history.late_accreted_fraction`, GEN.739) because the belts of its system do not exist yet.
Once the system is whole, `apply` takes that number from the system instead: the mass of the early belts that lands on the planet, past the other planets that share it
(`physics/late_supply.py`), and makes the planet's air, water and ground again on the draws they began with (`planets.remake_air`), so only what the late accretion changes is different.

* **The supply.** `late_supply.veneer_fraction` over the reservoirs of the planet's own star: its belts, the stretches no belt is left of, and the leftover planetesimals of its rocky zone (`leftover_of`, GEN.1493 slice 3). The planetesimals are followed past
  the system's bodies (`late_supply.fates`), so a giant between a belt and the planets takes what would reach them, and the leftovers of the planets' own zone, which lies inside the giants, carry the veneer.
* **The floor.** A planet no belt reaches (31 percent of 581 sampled planets) still took the planetesimals that never formed a belt: it gets at least `FLOOR_SHARE` of the mean the law gives its mass
  (`impact_history.late_accreted_percent_mean`). [Stand-in 6.255 LA1, tier T6, tweak L3.]
* **The ceiling.** The old cap `impact_history.LATE_ACCRETED_PERCENT_MAX` stands.
* **The spread.** There is no scatter of its own: the planets of different systems differ because their belts and neighbours do (the law's log-normal scatter is replaced by that). [Stand-in 6.255 LA2, tier T6, tweak L3.]

The air is made again for every rocky planet of a star (a moon, a giant and a planet with no star are left), from a stream of its own, so no draw of the run's stream moves; the cost is one more
solve of the air for each rocky planet.
"""

from planetgen.physics import constants, impact_history, late_supply, planets as planet_physics, planet_type
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

FLOOR_SHARE = 0.05
"""float: The least a planet takes, as a share of the mean the present law gives its mass: the planetesimals that never formed a belt. The lowest published veneers (a few hundredths of a percent for the Moon against
about half a percent for the Earth) are about 0.05 of the mean of a body of their mass, recalled, not registered. [Stand-in 6.255 LA1, tier T6, tweak L3.]"""

LEFTOVER_SHARE = 0.007
"""float: The leftover planetesimals of a star's rocky zone (the ones that rain on its planets for the first hundred million years, and whose fall is the veneer of the Earth), as a share of the mass of the rocky planets that grew from the zone: the planets
took the planetesimals at an efficiency that does not change much from one disk to the next, so what was left goes with what was taken. Set so the Sun's own terrestrial planets (1.8 Earth masses), with the giants beyond them, put 0.2 to 0.7 percent of the Earth's
mass on the Earth (test). [Stand-in 6.255 LA4, tier T6, tweak L4.]"""

LEFTOVER_ZONE = (0.8, 1.2)
"""tuple: The ring of the leftovers, as a multiple of the distance of the innermost and of the outermost rocky planet of the star (the planets' own feeding zones reach a fifth either side). [Stand-in 6.255 LA5, tier T6, tweak L1.]"""

EARTH_KG = 5.9722e24


def _primordial_kg(belt):
    return (getattr(belt, "primordial_mass_earth", None) or 0.0) * EARTH_KG


def early_of_belt(belt):
    """The belt as it was early (`late_supply.Early`): its primordial mass, its ring, and the share of it sent inward by what the giants and companions did to it (`late_supply.inward_share_of`)."""
    return late_supply.Early("belt", _primordial_kg(belt), belt.lower_limit, belt.upper_limit,
                             inward_share=late_supply.inward_share_of(getattr(belt, "scatter_fates", None), getattr(belt, "stirred", None)))


def early_of_origin(origin, lower_au, upper_au, fates):
    """
    A stretch a system made as a belt and then lost (the giants cleared it of everything, or a belt of what is left under the least a belt keeps) as it was early (`belt_birth.build` gives these to the
    system). It sent inward what the giants threw onto crossing orbits before it went, so the planets inside it keep their share of that (GEN.1493 slice 2); what the planets or the star took from it is not a leftover and counts for nothing.

    Args:
        origin: The stretch (`belt_origin.Origin`): its primordial mass, whether it is stirred.
        lower_au, upper_au (float): Its ring, cut to the orbits the star holds.
        fates (dict or None): Its shares after the outermost giant (`belt_birth._fates`).
    """
    return late_supply.Early("belt", origin.primordial_mass_earth * EARTH_KG, lower_au, upper_au,
                             inward_share=late_supply.inward_share_of(fates, origin.stirred, gone=True))


def leftover_of(star, siblings):
    """
    The leftover planetesimals of the zone a star's rocky planets stand in (`late_supply.Early`, all of it sent inward), or `None` for a star with no rocky planet: `LEFTOVER_SHARE` of the mass of the rocky planets between `LEFTOVER_ZONE` times the innermost and
    the outermost of them, the zone cut at the star's snow line (past which a rocky planet's zone is not the one the planets grew in).

    Args:
        star: The star the planets orbit.
        siblings: The bodies of the star (a rocky planet is a planet of type `t` that is not a moon).
    """
    rocky = [body for body in siblings if body.body_type == "t" and not body.is_moon and body.mass and body.distance > 0.0]
    if not rocky:
        return None
    nearest, farthest = min(body.distance for body in rocky), max(body.distance for body in rocky)
    lower = LEFTOVER_ZONE[0] * nearest
    upper = min(LEFTOVER_ZONE[1] * farthest, max(star.snow_line_au, lower * LEFTOVER_ZONE[1] / LEFTOVER_ZONE[0]))
    grown = sum(body.mass for body in rocky if body.distance <= upper)
    return late_supply.Early("leftover", LEFTOVER_SHARE * grown, lower, upper, inward_share=1.0)


def percent_for(planet, siblings, star_mass_solar, lost=()):
    """The late accretion of `planet`, percent of its mass, from the planets and belts `siblings` of its star (the planet among them) and the `lost` reservoirs of the star (`late_supply.Early`: stretches no belt is left of)."""
    mass_earth = planet.mass / EARTH_KG
    bodies = [(body.radius, body.distance, body.mass) for body in siblings
              if body.body_type in ("t", "g") and not body.is_moon and body.mass and body.radius]
    reservoirs = [early_of_belt(belt) for belt in siblings if belt.body_type == "a" and _primordial_kg(belt) > 0.0] + list(lost)
    leftover = leftover_of(planet.star, siblings)
    if leftover is not None:
        reservoirs.append(leftover)
    supply = 100.0 * late_supply.veneer_fraction((planet.radius, planet.distance, planet.mass), bodies, star_mass_solar, reservoirs)
    floor = FLOOR_SHARE * impact_history.late_accreted_percent_mean(mass_earth)
    return min(max(supply, floor), impact_history.LATE_ACCRETED_PERCENT_MAX)


def is_rocky_of_a_star(planet):
    """Whether the late accretion reaches the planet's air at all: a rocky planet of a star, with ground (a moon, a giant, a belt and a planet with no star are left)."""
    star = getattr(planet, "star", None)
    return (planet.body_type == "t" and not planet.is_moon and star is not None and not getattr(star, "starless", False)
            and bool(planet.mass) and bool(planet.radius) and planet_type.has_ground(planet))


def applies_to(planet):
    """Whether the planet's air is made again: a rocky planet of a star, with ground, whose air was made here."""
    return is_rocky_of_a_star(planet) and planet.build.air_replay is not None


def percent_in(planet, planets, lost=None):
    """The late accretion of `planet` from the planets and belts `planets` of the system (`StarSystem.all_planets`): those of its own star, and the reservoirs lost for it (`lost`: the `{star id: [Early]}` of the system)."""
    siblings = [body for body in planets if body.body_type == "a" and _same_star(body, planet) or body.body_type != "a" and body.star is planet.star]
    return percent_for(planet, siblings, _kg_to_msun(planet.star.mass), (lost or {}).get(id(planet.star), ()))


def give_first_air(planets, lost=None):
    """
    Gives every rocky planet whose air is still to be made (`Planet.finish_climate`, GEN.1548) the late accretion the system gives it now (PERF.271), so the air is made once with a number near the final one
    instead of with the law's draw and again at `apply`: the belts and the planets' orbits are known by then, and what is left to `apply` is what the belts promoted to planets add. The number is the same
    function (`percent_in`) the later pass uses; `planets.remake_air` still gates the second pass by how far the number moved the carbon and the nitrogen.
    """
    for planet in planets:
        if planet.body_type != "a" and planet.build.deferred_climate is not None and is_rocky_of_a_star(planet):
            planet.build.late_accreted_percent = percent_in(planet, planets, lost)


def apply(planets, lost=None, gate=planet_physics.LATE_ACCRETION_GATE):
    """Makes the air again, with the late accretion the system gives, of every rocky planet of `planets` (every planet and belt of the system, `StarSystem.all_planets`). A planet whose new late accretion moves its carbon and nitrogen by under `gate` of what its volcanoes gave keeps its air (`planets.remake_air`). `lost` is the `{star id: [Early]}` of the system's lost reservoirs. Returns how many were made again."""
    made = 0
    for planet in planets:
        if planet.body_type == "a" or not applies_to(planet):
            continue
        made += planet_physics.remake_air(planet, percent_in(planet, planets, lost), gate)
    return made


def _same_star(belt, planet):
    from planetgen.generation import belt as belt_module

    return belt_module._STARS.get(belt) is planet.star
