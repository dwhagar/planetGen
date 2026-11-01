# planetgen/admin/edits.py

"""
An admin's edits to a stored star system, made on the `StarSystem` object
`store.load_star_system` returns and then written back with
`db.edits.save_system_edits` (TODO ADM.1):

- ADM.8: regenerate or remove one planet, moon or asteroid belt.
- (ADM.6, changing a body's class, ended with the class cut-over: a class is derived from the finished body's conditions, so there
  is nothing to set; the edit of the inputs it follows from belongs to the 8.6 edit chain.)
- ADM.7: change a single-star system's star.

After every edit the system goes through `validation.stabilize_star_system`,
which moves bodies (never removes them) until it validates, and the
returned `EditResult` says what moved and what is still wrong.
"""

from collections import namedtuple

from planetgen.generation import validation
from planetgen.physics import constants, daylight, magnetism, planet_kind, planets as planetPhysics, radiation
from planetgen.generation.belt import AsteroidBelt
from planetgen.names.bodies import moon_letters
from planetgen.generation.planet import Planet
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

_au_to_km = units.converter("au", "km")

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


def _finish(system, summary):
    """Stabilizes the edited system and wraps up the `EditResult`."""
    report = validation.stabilize_star_system(system)
    warnings = [f"{p.body}: {p.message}." for p in report.problems]
    return EditResult(summary, report.moved, report.reclassified, report.removed, warnings)


def regenerate_planet(system, planet, owner):
    """
    Replaces a planet with a freshly generated one at the same orbit: a
    new mass, make, atmosphere, life and moons (and so, derived from them, class).
    It keeps its name and its database row, so facilities on it stay.
    """
    fresh = Planet(system.system_config, planet.star, planet.star.habitable_zone, planet.distance)
    fresh.name = planet.name
    fresh.db_id = planet.db_id
    _name_moons(fresh)
    for moon in fresh.moons:
        radiation.apply_parent_belts(moon, fresh)  # GEN.983
    for body in [fresh] + fresh.moons:
        validation.reapply_life(body)
    magnetism.generate_stellar_induction(fresh)  # GEN.981
    for moon in fresh.moons:
        magnetism.generate_induced_field(moon, fresh)
    owner[owner.index(planet)] = fresh
    return _finish(system, f"Regenerated {planet.name}: now class {fresh.planet_class} "
                           f"with {len(fresh.moons)} moon{'' if len(fresh.moons) == 1 else 's'}.")


def make_moon(planet, distance_au):
    """A new moon of `planet` at `distance_au` from it: a mass from the moon size law (`planet_kind.draw_moon_mass_kg`), or `None`
    when the planet can hold none."""
    mass_kg = planet_kind.draw_moon_mass_kg(planet.mass, planet.body_type == "g", planetPhysics.moon_size_limits(planet)[1])
    if mass_kg is None:
        return None
    moon = Planet(planet.system_config, planet.star, planet.habitable_zone, distance_au,
                  mass=mass_kg, zone_override=planet.zone,
                  distance_override=planet.distance, is_moon=True, primary_mass_kg=planet.mass)
    magnetism.generate_induced_field(moon, planet)
    radiation.apply_parent_belts(moon, planet)  # GEN.983
    daylight.generate(moon, planet)
    return moon


def regenerate_moon(system, moon, planet):
    """
    Replaces a moon with a freshly generated one at the same orbit around
    its planet, keeping its name and row.

    Raises:
        ValueError: If the planet can hold no moon.
    """
    fresh = make_moon(planet, moon.distance)
    if fresh is None:
        raise ValueError(f"{planet.name} can hold no moon")
    fresh.name = moon.name
    fresh.db_id = moon.db_id
    validation.reapply_life(fresh)
    planet.moons[planet.moons.index(moon)] = fresh
    validation.refresh_planet_reflex(planet)
    return _finish(system, f"Regenerated {moon.name}: now class {fresh.planet_class}.")


def regenerate_belt(system, belt, owner):
    """Replaces an asteroid belt with a fresh roll (density and
    composition) over the same span, keeping its row."""
    star = system.star if owner is system.planets else system.secondary_star
    fresh = AsteroidBelt(system.system_config, star, belt.distance, belt.lower_limit, belt.upper_limit)
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


# ---------------------------------------------------------------------
# ADM.7: a system's star
# ---------------------------------------------------------------------

def _scale_belt(belt, factor):
    belt.distance *= factor
    belt.lower_limit *= factor
    belt.upper_limit *= factor


def _rehost(body, star):
    """Points a planet (and its moons) at a new star and its habitable zone."""
    body.star = star
    body.habitable_zone = star.habitable_zone
    for moon in body.moons:
        moon.star = star
        moon.habitable_zone = star.habitable_zone


def _rescale_comet(comet, factor, mass_ratio, new_mass_solar):
    """Moves a comet's orbit out by `factor` around a star `mass_ratio`
    times as heavy (Kepler: period ~ a^1.5 / sqrt(M), speed ~
    sqrt(M / a))."""
    comet.perihelion_distance_au *= factor
    comet.distance_au *= factor
    comet.set_position_au(comet.position_x_au * factor, comet.position_y_au * factor, comet.position_z_au * factor)
    if comet.orbital_speed_kms is not None:
        comet.orbital_speed_kms *= (mass_ratio / factor) ** 0.5
    if comet.orbital_period_years is not None:
        comet.orbital_period_years *= factor ** 1.5 / mass_ratio ** 0.5
    comet.primary_mass_solar = new_mass_solar


def change_star(system, star_type):
    """
    Replaces a single-star system's star with a newly generated one of
    `star_type` (e.g. "K2V") (ADM.7). Every planet, moon and belt stays,
    with its mass and make (its air, temperature and class are made again for the new light):
    their orbits are scaled with the square root of the
    luminosity, so each keeps its place relative to the habitable zone,
    then the innermost is moved clear of the new star's engulfment radius
    and everything is re-spaced for the new star's mass. Bodies left
    beyond the new star's farthest stable orbit are removed, outermost
    first, and named in the result. Comets follow the same scaling.

    Returns:
        tuple: `(EditResult, new_star)`; `new_star` carries the old one's
        `db_id`, name and galactic orbit.

    Raises:
        ValueError: For a binary system, a black hole or neutron star, or
            a malformed star type.
    """
    from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
    from planetgen.generation.star import STAR_TYPE_PATTERN, Star
    from planetgen.generation.system import StarSystem

    if system.binary_type is not None:
        raise ValueError("only a single star can be changed; regenerate a binary system instead")
    old = system.star
    if isinstance(old, (BlackHole, NeutronStar)):
        raise ValueError("a black hole or neutron star can't be changed into a star")
    star_type = (star_type or "").strip().upper()
    if not STAR_TYPE_PATTERN.fullmatch(star_type):
        raise ValueError(f"{star_type!r} is not a spectral type: expected a class letter (OBAFGKM), a subclass "
                         "digit (0-9) and a Yerkes class (0, IA+, IA, IAB, IB, II, III, IV, V, VI, VII or D), "
                         "e.g. G2V")

    config = system.system_config
    forced_type = config.STAR_TYPE
    config.STAR_TYPE = star_type
    try:
        new = Star(config, name=old.name, galactic_orbital_phase_deg=old.galactic_orbital_phase_deg)
    finally:
        config.STAR_TYPE = forced_type
    new.db_id = getattr(old, "db_id", None)
    # Same place in the galaxy: the galactic orbit stays, and the star's
    # Hill sphere against the galaxy scales with the cube root of its mass.
    for field in ("galactic_orbital_speed_kms", "galactic_orbital_period_gy",
                  "galactic_orbital_phase_deg", "galactic_min_update_interval_years"):
        setattr(new, field, getattr(old, field))
    new.system_perimeter = old.system_perimeter * (new.mass / old.mass) ** (1 / 3)

    factor = (new.luminosity / old.luminosity) ** 0.5
    mass_ratio = new.mass / old.mass
    planets = system.planets
    for body in planets:
        if body.body_type == 'a':
            _scale_belt(body, factor)
        else:
            body.distance *= factor
            _rehost(body, new)
    system.star = system.primary_star = new
    system.stars = [new]

    floor_au = StarSystem._engulfment_radius_au(new)
    if planets and validation.inner_edge_au(planets[0]) < floor_au:
        shift = floor_au * 1.01 - validation.inner_edge_au(planets[0])
        for body in planets:
            if body.body_type == 'a':
                body.distance += shift
                body.lower_limit += shift
                body.upper_limit += shift
            else:
                body.distance += shift

    dropped_moons = []
    for body in planets:
        if body.body_type == 'a':
            continue
        validation.reconcile_moved_planet(body)
        # A planet moved in has a smaller Hill sphere: the moons it can no
        # longer hold are lost, like the bodies past the new ceiling.
        before = list(body.moons)
        planetPhysics.drop_unstable_moons(body)
        dropped_moons += [m.name for m in before if m not in body.moons]
        validation.refresh_planet_reflex(body, force=True)
        for item in [body] + body.moons:
            validation.reapply_life(item)

    for comet in system.comets:
        _rescale_comet(comet, factor, mass_ratio, _kg_to_msun(new.mass))

    report = validation.stabilize_star_system(system, allow_removal=True)
    warnings = [f"{p.body}: {p.message}." for p in report.problems]
    removed = report.removed + dropped_moons
    summary = f"The star is now {new.type}; orbits were scaled by {factor:.3g}."
    if removed:
        summary += f" No room for {len(removed)} bod{'y' if len(removed) == 1 else 'ies'}: " \
                   f"{', '.join(removed)} removed."
    return EditResult(summary, report.moved, report.reclassified, removed, warnings), new
