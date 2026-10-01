# stellarObjects/adminEdits.py

"""
An admin's edits to a stored star system, made on the `StarSystem` object
`_db.load_star_system` returns and then written back with
`editStore.save_system_edits` (TODO ADM.1):

- ADM.8: regenerate or remove one planet, moon or asteroid belt.

After every edit the system goes through `validation.stabilize_star_system`,
which moves bodies (never removes them) until it validates, and the
returned `EditResult` says what moved and what is still wrong.
"""

import random
from collections import namedtuple

from . import planetPhysics, program_constants, validation
from .asteroidData import AsteroidBelt
from .bodyNames import moon_letters
from .planetData import Planet

EditResult = namedtuple("EditResult", ["summary", "moved", "reclassified", "removed", "warnings"])
"""What an edit did: a one-line `summary`, the names of the bodies the
stabilize pass `moved`, `reclassified` and `removed`, and `warnings`
(each a plain sentence) for whatever still doesn't validate."""

BODY_KINDS = ("planet", "moon", "belt")
"""tuple: The kinds of body these edits act on."""


class BodyNotFound(LookupError):
    """No body of that kind and id in the system."""


def find_body(system, kind, body_id):
    """
    Finds a loaded body by kind and row id.

    Returns:
        tuple: `(body, owner)`: for a planet or belt, `owner` is the list
        that holds it (one star's `planets`); for a moon, its planet.

    Raises:
        BodyNotFound: If the system has no such body.
    """
    for _star, planets in validation.star_lists(system):
        for body in planets:
            if body.body_type == 'a':
                if kind == "belt" and getattr(body, "db_id", None) == body_id:
                    return body, planets
                continue
            if kind == "planet" and getattr(body, "db_id", None) == body_id:
                return body, planets
            if kind == "moon":
                for moon in body.moons:
                    if getattr(moon, "db_id", None) == body_id:
                        return moon, body
    raise BodyNotFound(f"no such {kind}: {body_id}")


def _name_moons(planet):
    """Names a planet's moons `<planet name><letter>`, in orbital order."""
    planet.moons.sort(key=lambda m: m.distance)
    for index, moon in enumerate(planet.moons):
        moon.name = f"{planet.name}{moon_letters(index)}"


def _finish(system, summary, pinned=()):
    """Stabilizes the edited system and wraps up the `EditResult`."""
    report = validation.stabilize_star_system(system, pinned=pinned)
    warnings = [f"{p.body}: {p.message}." for p in report.problems]
    return EditResult(summary, report.moved, report.reclassified, report.removed, warnings)


def regenerate_planet(system, planet, owner):
    """
    Replaces a planet with a freshly generated one at the same orbit: a
    new class (any that fits its zone), size, atmosphere, life and moons.
    It keeps its name and its database row, so facilities on it stay.
    """
    fresh = Planet(system.system_config, planet.star, planet.star.habitable_zone, planet.distance)
    if fresh.distance != planet.distance:
        # An ecosphere class places itself within the habitable zone;
        # the regenerated planet stays on the old orbit.
        fresh.distance = planet.distance
        validation.reconcile_moved_planet(fresh, keep_class=True)
    fresh.name = planet.name
    fresh.db_id = planet.db_id
    _name_moons(fresh)
    for body in [fresh] + fresh.moons:
        validation.reapply_life(body)
    owner[owner.index(planet)] = fresh
    return _finish(system, f"Regenerated {planet.name}: now class {fresh.planet_class} "
                           f"with {len(fresh.moons)} moon{'' if len(fresh.moons) == 1 else 's'}.")


def moon_classes(planet):
    """The classes a moon of `planet` may roll, as `generate_moons` picks
    them: terrestrial, not barred from being a moon, valid in the planet's
    zone, and light and small enough for the planet to hold."""
    max_mass = planet.mass / 10
    max_radius = validation.max_moon_radius_km(planet)
    return [c for c, (_low, high) in planetPhysics.planet_mass_ranges.items()
            if program_constants.PLANET_CLASSES[c][planet.zone]
            and program_constants.PLANET_CLASSES[c]["type"] == 't'
            and c not in program_constants.MOON_BLACKLIST
            and high <= max_mass
            and program_constants.PLANET_CLASSES[c]["radius_range"][1] <= max_radius]


def make_moon(planet, distance_au, moon_class):
    """A new moon of `planet` of `moon_class` at `distance_au` from it."""
    max_radius = validation.max_moon_radius_km(planet)
    low, high = program_constants.PLANET_CLASSES[moon_class]["radius_range"]
    radius = planetPhysics._sample_class_radius(moon_class, low, min(high, max_radius))
    return Planet(planet.system_config, planet.star, planet.habitable_zone, distance_au,
                  radius=radius, planet_class=moon_class, zone_override=planet.zone,
                  distance_override=planet.distance, is_moon=True, primary_mass_kg=planet.mass)


def regenerate_moon(system, moon, planet):
    """
    Replaces a moon with a freshly generated one at the same orbit around
    its planet, keeping its name and row.

    Raises:
        ValueError: If no moon class fits this planet.
    """
    classes = moon_classes(planet)
    if not classes:
        raise ValueError(f"no moon class fits {planet.name}")
    fresh = make_moon(planet, moon.distance, random.choice(classes))
    fresh.name = moon.name
    fresh.db_id = moon.db_id
    validation.reapply_life(fresh)
    planet.moons[planet.moons.index(moon)] = fresh
    validation.refresh_planet_reflex(planet)
    return _finish(system, f"Regenerated {moon.name}: now class {fresh.planet_class}.")


def regenerate_belt(system, belt, owner):
    """Replaces an asteroid belt with a fresh roll (density and
    composition) over the same span, keeping its row."""
    fresh = AsteroidBelt(system.system_config, belt.distance, belt.lower_limit, belt.upper_limit)
    fresh.db_id = belt.db_id
    owner[owner.index(belt)] = fresh
    return _finish(system, f"Regenerated the asteroid belt: now {fresh.density}, "
                           f"{fresh.get_composition_summary()}.")


def remove_body(system, body, owner):
    """
    Removes a planet (with its moons), a moon or an asteroid belt from the
    system. The other bodies keep their names and orbits.
    """
    if isinstance(owner, list):
        owner.remove(body)
    else:
        owner.moons.remove(body)
        validation.refresh_planet_reflex(owner, force=True)
    label = body.name if getattr(body, "name", None) else "the asteroid belt"
    return _finish(system, f"Removed {label}.")


def regenerate_body(system, kind, body, owner):
    """`regenerate_planet`, `regenerate_moon` or `regenerate_belt` by kind."""
    if kind == "planet":
        return regenerate_planet(system, body, owner)
    if kind == "moon":
        return regenerate_moon(system, body, owner)
    return regenerate_belt(system, body, owner)
