# planetgen/generation/late_accretion.py

"""
The air reads the system's late accretion (GEN.1493 slice 1, docs/design/composition/late-accretion-from-the-system.md)
=========================================================================================================================

A rocky planet's air is made when the planet is, with one late-accretion number drawn from a law (`impact_history.late_accreted_fraction`, GEN.739) because the belts of its system do not exist yet.
Once the system is whole, `apply` takes that number from the system instead: the mass of the early belts that lands on the planet, past the other planets that share it
(`physics/late_supply.py`), and makes the planet's air, water and ground again on the draws they began with (`planets.remake_air`), so only what the late accretion changes is different.

* **The supply.** `late_supply.veneer_fraction` over the planets and belts of the planet's own star.
* **The floor.** A planet no belt reaches (31 percent of 581 sampled planets) still took the planetesimals that never formed a belt: it gets at least `FLOOR_SHARE` of the mean the law gives its mass
  (`impact_history.late_accreted_percent_mean`). [Stand-in 6.255 LA1, tier T6, tweak L3.]
* **The ceiling.** The old cap `impact_history.LATE_ACCRETED_PERCENT_MAX` stands.
* **The spread.** There is no scatter of its own: the planets of different systems differ because their belts and neighbours do (the law's log-normal scatter is replaced by that). [Stand-in 6.255 LA2, tier T6, tweak L3.]

The air is made again for every rocky planet of a star (a moon, a giant and a planet with no star are left), from a stream of its own, so no draw of the run's stream moves; the cost is one more
solve of the air for each rocky planet.
"""

from planetgen.physics import constants, impact_history, influx, late_supply, planets as planet_physics, planet_type
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

FLOOR_SHARE = 0.05
"""float: The least a planet takes, as a share of the mean the present law gives its mass: the planetesimals that never formed a belt. The lowest published veneers (a few hundredths of a percent for the Moon against
about half a percent for the Earth) are about 0.05 of the mean of a body of their mass, recalled, not registered. [Stand-in 6.255 LA1, tier T6, tweak L3.]"""

EARTH_KG = 5.9722e24


def _primordial_kg(belt):
    return (getattr(belt, "primordial_mass_earth", None) or 0.0) * EARTH_KG


def percent_for(planet, siblings, star_mass_solar):
    """The late accretion of `planet`, percent of its mass, from the planets and belts `siblings` of its star (the planet among them)."""
    mass_earth = planet.mass / EARTH_KG
    bodies = [(body.radius, body.distance, body.mass) for body in siblings
              if body.body_type in ("t", "g") and not body.is_moon and body.mass and body.radius]
    reservoirs = [influx.Reservoir("belt", _primordial_kg(belt), belt.lower_limit, belt.upper_limit) for belt in siblings
                  if belt.body_type == "a" and _primordial_kg(belt) > 0.0]
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
    return is_rocky_of_a_star(planet) and getattr(planet, "_air_replay", None) is not None


def percent_in(planet, planets):
    """The late accretion of `planet` from the planets and belts `planets` of the system (`StarSystem.all_planets`): those of its own star."""
    siblings = [body for body in planets if body.body_type == "a" and _same_star(body, planet) or body.body_type != "a" and body.star is planet.star]
    return percent_for(planet, siblings, _kg_to_msun(planet.star.mass))


def give_first_air(planets):
    """
    Gives every rocky planet whose air is still to be made (`Planet.finish_climate`, GEN.1548) the late accretion the system gives it now (PERF.271), so the air is made once with a number near the final one
    instead of with the law's draw and again at `apply`: the belts and the planets' orbits are known by then, and what is left to `apply` is what the belts promoted to planets add. The number is the same
    function (`percent_in`) the later pass uses; `planets.remake_air` still gates the second pass by how far the number moved the carbon and the nitrogen.
    """
    for planet in planets:
        if planet.body_type != "a" and getattr(planet, "_deferred_climate", None) is not None and is_rocky_of_a_star(planet):
            planet._late_accreted_percent = percent_in(planet, planets)


def apply(planets, gate=planet_physics.LATE_ACCRETION_GATE):
    """Makes the air again, with the late accretion the system gives, of every rocky planet of `planets` (every planet and belt of the system, `StarSystem.all_planets`). A planet whose new late accretion moves its carbon and nitrogen by under `gate` of what its volcanoes gave keeps its air (`planets.remake_air`). Returns how many were made again."""
    made = 0
    for planet in planets:
        if planet.body_type == "a" or not applies_to(planet):
            continue
        made += planet_physics.remake_air(planet, percent_in(planet, planets), gate)
    return made


def _same_star(belt, planet):
    from planetgen.generation import belt as belt_module

    return belt_module._STARS.get(belt) is planet.star
