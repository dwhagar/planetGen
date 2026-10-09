# planetgen/admin/edits.py

"""
An admin's edits to a stored star system, made on the `StarSystem` object
`_db.load_star_system` returns and then written back with
`editStore.save_system_edits` (TODO ADM.1):

- ADM.8: regenerate or remove one planet, moon or asteroid belt.
- ADM.6: change a planet's or moon's class, to a recommended one (fits
  where it is) or a forced one (anything).
- ADM.7: change a single-star system's star.

After every edit the system goes through `validation.stabilize_star_system`,
which moves bodies (never removes them) until it validates, and the
returned `EditResult` says what moved and what is still wrong.
"""

import math
from collections import namedtuple

from planetgen.generation import validation
from planetgen.physics import constants, planets as planetPhysics
from planetgen import tuning
from planetgen.generation.belt import AsteroidBelt
from planetgen.names.bodies import moon_letters
from planetgen.generation.planet import Planet
from planetgen.util import draw

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
    them (`planetPhysics.moon_class_options`): terrestrial, not barred
    from being a moon, valid in the planet's zone, and able to fit under
    the size and mass the planet can hold."""
    return list(planetPhysics.moon_class_options(planet, planet.zone))


def make_moon(planet, distance_au, moon_class):
    """A new moon of `planet` of `moon_class` at `distance_au` from it."""
    ceiling = planetPhysics.moon_class_options(planet, planet.zone)[moon_class]
    low = tuning.PLANET_CLASSES[moon_class]["radius_range"][0]
    radius = planetPhysics._sample_class_radius(moon_class, low, ceiling)
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
    fresh = make_moon(planet, moon.distance, draw.choice(classes))
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


# ---------------------------------------------------------------------
# ADM.6: a planet's or moon's class
# ---------------------------------------------------------------------

def class_fits_mass(planet_class, mass_kg):
    """Whether a body of `mass_kg` can be `planet_class` (ADM.27): its
    mass lies in the class's mass range (`planetPhysics.planet_mass_ranges`),
    since a class change keeps the body's mass."""
    low, high = planetPhysics.planet_mass_ranges[planet_class]
    return low * (1 - 1e-9) <= mass_kg <= high * (1 + 1e-9)


def _mass_refusal(body, planet_class):
    """The message refusing `planet_class` for `body` because of its mass."""
    low, high = planetPhysics.planet_mass_ranges[planet_class]
    earth = constants.EARTH_MASS_TO_KG
    return (f"{body.name} is {body.mass / earth:.3g} Earth masses; class {planet_class} needs "
            f"{low / earth:.3g} to {high / earth:.3g}, so its class wasn't changed")


def _hill_radius_au(distance_au, mass_kg, star_mass_kg):
    return distance_au * (mass_kg / (3 * star_mass_kg)) ** (1 / 3)


def _fits_between(body, owner, mass_kg):
    """Whether a planet of `mass_kg` at `body`'s orbit keeps the spacing
    `validation.space_orbits` asks of its inner and outer neighbors."""
    index = owner.index(body)
    star_mass = body.star.mass
    stand_in = _StandIn(body.distance, mass_kg, body.star,
                        _hill_radius_au(body.distance, mass_kg, star_mass))
    if index > 0:
        inner = owner[index - 1]
        if body.distance < validation.min_distance_after_au(stand_in, inner) * (1 - validation.RELATIVE_TOLERANCE):
            return False
    if index + 1 < len(owner):
        outer = owner[index + 1]
        if outer.body_type == 'a':
            if outer.distance - body.distance < stand_in.min_orbit_distance:
                return False
        elif outer.distance < validation.mutual_min_distance_au(outer, stand_in) * (1 - validation.RELATIVE_TOLERANCE):
            return False
    return True


def _can_hold_moons(body, planet_class):
    """Whether some planet of `planet_class` could hold `body`'s moons:
    big enough for the largest. Its mass, and so its Hill sphere, stays
    what it is (ADM.27)."""
    if not body.moons:
        return True
    largest = max(m.radius for m in body.moons)
    low, high = tuning.PLANET_CLASSES[planet_class]["radius_range"]
    return (low * high) ** 0.5 / (10 ** (1 / 3)) >= largest


def _roll_fits(body, owner):
    """Whether `body`'s current roll needs nothing moved: a planet clears
    its neighbors and holds its moons where they are; a moon fits its
    planet."""
    if body.is_moon:
        return body.radius <= validation.max_moon_radius_km(owner) and body.mass <= owner.mass / 10
    if not _fits_between(body, owner, body.mass) or not body.moons:
        return _fits_between(body, owner, body.mass)
    # Moons may still be nudged outward past a grown planet, but they must
    # fit and stay inside its Hill sphere.
    _low_km, high_km = planetPhysics.moon_orbit_bounds_km(body)
    return (max(m.radius for m in body.moons) <= validation.max_moon_radius_km(body)
            and max(m.distance for m in body.moons) * constants.AU_TO_KM <= high_km)


RECOMMENDED_ROLLS = 40
"""int: How many rolls of a recommended class are tried for one that
needs nothing moved, before keeping the last one."""


class _StandIn:
    """A planet-shaped stand-in for spacing arithmetic: `distance` (AU),
    `mass` (kg), `star`, `hill_radius` (km) and `min_orbit_distance` (AU)."""
    body_type = 't'

    def __init__(self, distance, mass, star, hill_radius_au):
        self.distance = distance
        self.mass = mass
        self.star = star
        self.hill_radius = hill_radius_au * constants.AU_TO_KM
        self.min_orbit_distance = 5 * hill_radius_au


def recommended_classes(system, body, owner):
    """
    The classes `body` could take without moving any other planet ("recommended",
    ADM.6), best first by how common they are: valid in its zone, not a
    habitable class where the system rules those out, a class its mass
    fits (`class_fits_mass`; the change keeps it, ADM.27), and -- for a
    planet -- big enough to hold its moons; for a moon, a class its planet
    can hold (`moon_classes`). Its current class is left out.
    """
    if body.is_moon:
        candidates = moon_classes(owner)
    else:
        zone = validation.zone_for(body.habitable_zone, body.distance)
        candidates = [
            c for c, data in tuning.PLANET_CLASSES.items()
            if data[zone] and _can_hold_moons(body, c)
        ]
        if planetPhysics._habitable_classes_barred(body, zone):
            candidates = [c for c in candidates if c not in tuning.HABITABLE_PLANET_CLASSES]
    candidates = [c for c in candidates if c != body.planet_class and class_fits_mass(c, body.mass)]
    return sorted(candidates, key=lambda c: -tuning.PLANET_CLASS_PROBABILITIES.get(c, 0))


def class_options(system):
    """`{"planet:<id>" or "moon:<id>": recommended classes}` for every
    loaded planet and moon (`recommended_classes`)."""
    options = {}
    for _star, planets in validation.star_lists(system):
        for body in planets:
            if body.body_type == 'a':
                continue
            options[f"planet:{body.db_id}"] = recommended_classes(system, body, planets)
            for moon in body.moons:
                options[f"moon:{moon.db_id}"] = recommended_classes(system, moon, body)
    return options


def _keep_mass(body, mass_kg):
    """Puts `body` back at `mass_kg` after a re-roll drew its own: a giant
    on the mass-radius relation takes the relation's radius for it; any
    other class keeps its drawn density where the radius that gives stays
    in the class's range, else the nearest end of that range and the
    density the mass then makes (in range too, since the mass fits the
    class)."""
    body.mass = mass_kg
    if planetPhysics.uses_giant_mass_radius(body.planet_class):
        body.radius = None
        planetPhysics._apply_giant_mass_and_radius(body, radius_given=False, mass_given=True)
    else:
        low, high = tuning.PLANET_CLASSES[body.planet_class]["radius_range"]
        radius_m = (3 * mass_kg / (4 * math.pi * body.density * 1000)) ** (1 / 3)
        body.radius = min(max(radius_m / constants.KM_TO_M_FACTOR, low), high)
        radius_m = body.radius * constants.KM_TO_M_FACTOR
        body.density = mass_kg / ((4 / 3) * math.pi * radius_m ** 3) / 1000
    body.volume = (4 / 3) * math.pi * body.radius ** 3


def _apply_class(body, owner, planet_class):
    """Re-generates `body` as `planet_class` in place (ADM.27): radius,
    density, composition, atmosphere, temperature, pressure, life and
    everything derived, keeping its orbit, mass and name. A class not
    valid in the body's zone is generated as if it were in a zone where it
    is, then put back in its real zone (a forced change)."""
    mass = body.mass
    distance = body.distance
    parent_distance = owner.distance if body.is_moon else None
    real_zone = validation.zone_for(body.habitable_zone, parent_distance if body.is_moon else distance)
    data = tuning.PLANET_CLASSES[planet_class]
    gen_zone = real_zone if data[real_zone] else next(z for z in "ech" if data[z])
    config = body.system_config
    habitable_rule = config.HABITABLE_WORLD
    config.HABITABLE_WORLD = None  # an admin's choice overrides the system's no-habitable-world rule
    try:
        body.planet_class = planet_class
        body.radius = None
        body.mass = None
        planetPhysics.generate_planet_properties(body, zone_override=gen_zone)
    finally:
        config.HABITABLE_WORLD = habitable_rule
    _keep_mass(body, mass)
    body.distance = distance
    body.zone = real_zone
    primary_mass = owner.mass if body.is_moon else body.star.mass
    body.period = planetPhysics.calculate_orbital_period_years(body.distance, primary_mass)
    planetPhysics.calculate_surface_gravity(body)
    planetPhysics.calculate_atmospheric_conditions(body, parent_distance)
    planetPhysics.generate_orbital_motion_properties(body, primary_mass)
    planetPhysics.update_hill_sphere(body)
    validation.reapply_life(body)
    if not body.is_moon:
        for moon in body.moons:
            validation.refresh_moon_orbit(moon, body)
        validation.refresh_planet_reflex(body, force=True)
    else:
        validation.refresh_planet_reflex(owner)


def change_class(system, body, owner, planet_class, force=False):
    """
    Changes a planet's or moon's class (ADM.6), re-generating its surface
    conditions as that class and keeping its orbit, mass and name
    (ADM.27). Without `force` only a recommended class
    (`recommended_classes`) is accepted; with it, any known class its mass
    fits, even one that can't exist where the body orbits. The body
    keeps its new class no matter what (it is pinned); the rest of the
    system is then re-spaced from the moons outward until it validates,
    and whatever still doesn't (a forced class out of its zone, a moon too
    large for its planet) comes back as warnings. Nothing is removed.

    Raises:
        ValueError: For an unknown class, a class the body's mass doesn't
            fit, or a class that isn't recommended without `force`.
            Nothing is changed.
    """
    if planet_class not in tuning.PLANET_CLASSES:
        raise ValueError(f"unknown class: {planet_class}")
    if not class_fits_mass(planet_class, body.mass):
        raise ValueError(_mass_refusal(body, planet_class))
    if not force and planet_class not in recommended_classes(system, body, owner):
        raise ValueError(f"class {planet_class} doesn't fit {body.name} where it is; force it to set it anyway")
    old_class = body.planet_class
    for _ in range(1 if force else RECOMMENDED_ROLLS):
        _apply_class(body, owner, planet_class)
        if force or _roll_fits(body, owner):
            break
    result = _finish(system, f"{body.name} changed from class {old_class} to class {planet_class}.", pinned=[body])
    return result


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
    with its class: their orbits are scaled with the square root of the
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
        validation.reconcile_moved_planet(body, keep_class=True)
        # A planet moved in has a smaller Hill sphere: the moons it can no
        # longer hold are lost, like the bodies past the new ceiling.
        before = list(body.moons)
        planetPhysics.drop_unstable_moons(body)
        dropped_moons += [m.name for m in before if m not in body.moons]
        validation.refresh_planet_reflex(body, force=True)
        for item in [body] + body.moons:
            validation.reapply_life(item)

    for comet in system.comets:
        _rescale_comet(comet, factor, mass_ratio, new.mass / constants.SOLAR_MASS_TO_KG)

    # Every planet keeps its class, as promised above: re-spacing must not
    # reclassify one it nudged into another zone (TEST.80).
    kept = [body for body in system.planets if body.body_type != 'a']
    report = validation.stabilize_star_system(system, pinned=kept, allow_removal=True)
    warnings = [f"{p.body}: {p.message}." for p in report.problems]
    removed = report.removed + dropped_moons
    summary = f"The star is now {new.type}; orbits were scaled by {factor:.3g}."
    if removed:
        summary += f" No room for {len(removed)} bod{'y' if len(removed) == 1 else 'ies'}: " \
                   f"{', '.join(removed)} removed."
    return EditResult(summary, report.moved, report.reclassified, removed, warnings), new
