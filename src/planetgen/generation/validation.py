# planetgen/generation/validation.py

"""
One place to validate a planet, a lunar system (a planet and its moons)
and a star system (TODO ADM.5).

Three kinds of function live here:

- **Checks** (`check_planet`, `check_lunar_system`, `check_star_system`)
  read a body or system and return a list of `Problem`s; they never change
  anything. A freshly generated system almost never has any.
- **Spacing** (`space_orbits`, `reconcile_moved_planet`,
  `trim_to_orbit_ceiling`, `cross_star_clearance`, and the distance
  helpers they use) is the orbit-spacing pass generation has always run.
  `StarSystem.validate_system`, `_trim_to_orbit_ceiling` and
  `_validate_cross_star_clearance` (`generation/system.py`) now call these, with
  the same results.
- **Stabilizing** (`stabilize_lunar_system`, `stabilize_star_system`) is
  for a system someone edited after generation (an admin's star
  override, ADM.7): it re-spaces moons and orbits outward until
  the checks pass, and reports what moved, what changed class, what
  was removed and what is still wrong.
"""

import math
from collections import namedtuple

from planetgen.generation import life as planetLife
from planetgen.physics import constants, daylight, ocean_tides, orbit_elements, planet_classes, planets as planetPhysics, rings
from planetgen import tuning
from planetgen.util import draw
from planetgen.physics.orbits import calculate_reflex_offset, mutual_hill_radius_au
from planetgen.physics import unit_registry as units

_km_to_au = units.converter("km", "au")
_au_to_km = units.converter("au", "km")

RELATIVE_TOLERANCE = 1e-9
"""float: Slack the checks allow on a distance comparison, so a body the
spacing pass put exactly on its limit (then saved and reloaded through km)
isn't reported for a rounding error."""

Problem = namedtuple("Problem", ["body", "message"])
"""A failed check: `body` is the name of the body (or system) it is
about, `message` says what is wrong, in plain words."""

StabilizeReport = namedtuple("StabilizeReport", ["moved", "reclassified", "removed", "problems"])
"""What `stabilize_star_system` did: the names of the bodies it `moved`,
`reclassified` (their class, derived from their conditions, changed because they moved into another
zone or a different light) and `removed`, and the `Problem`s still left (empty when the system
is stable)."""


def _label(body):
    """A body's name, or a stand-in for one that has none yet."""
    name = getattr(body, "name", None)
    if name:
        return name
    if getattr(body, "body_type", None) == 'a':
        return f"asteroid belt at {body.distance:.4g} AU"
    if hasattr(body, "perihelion_distance_au"):
        return f"comet with perihelion at {body.perihelion_distance_au:.4g} AU"
    return f"body at {body.distance:.4g} AU"


def _below(value, limit):
    """`value < limit` beyond `RELATIVE_TOLERANCE`."""
    return value < limit - abs(limit) * RELATIVE_TOLERANCE


def edge_au(obj):
    """The outer edge of a planet (its orbit) or a belt (`upper_limit`), in AU."""
    return obj.upper_limit if obj.body_type == 'a' else obj.distance


def inner_edge_au(obj):
    """The inner edge of a planet (its orbit) or a belt (`lower_limit`), in AU."""
    return obj.lower_limit if obj.body_type == 'a' else obj.distance


def periapsis_au(planet):
    """The closest a planet or moon comes to what it orbits, in AU (`distance` is its semi-major axis, GEN.727)."""
    return orbit_elements.periapsis_au(planet.distance, planet.eccentricity)


def apoapsis_au(planet):
    """The farthest a planet or moon goes from what it orbits, in AU."""
    return orbit_elements.apoapsis_au(planet.distance, planet.eccentricity)


def reach_au(obj):
    """How far out a body's influence reaches, in AU: a belt's
    `upper_limit`, or a planet's apoapsis plus its own `min_orbit_distance`
    (5 Hill radii), so its whole Hill sphere, not just its center, must
    stay inside a star's orbit ceiling."""
    return obj.upper_limit if obj.body_type == 'a' else apoapsis_au(obj) + obj.min_orbit_distance


def exceeds_orbit_ceiling(obj, ceiling_au):
    """True if `obj` reaches past `ceiling_au` (`reach_au`)."""
    return reach_au(obj) > ceiling_au


def orbit_ceiling_au(star):
    """
    The AU distance beyond which no body can stably orbit `star`: its Hill
    sphere against the galaxy (`system_perimeter`), or, for one star of a
    wide binary, the tighter of that and its Holman & Wiegert limit
    against the companion (`a_crit_au`).
    """
    if star.a_crit_au is None:
        return star.system_perimeter
    return min(star.system_perimeter, star.a_crit_au)


# ---------------------------------------------------------------------
# Spacing distances
# ---------------------------------------------------------------------

def min_distance_past_belt_au(planet, belt):
    """
    The closest `planet` can sit outside `belt`: its own Hill sphere
    must clear the belt's outer edge by the same 5 Hill radii
    (`min_orbit_distance`) a belt placed after a planet must keep, plus
    `MIN_ASTEROID_BELT_SEPARATION`.

    The Hill radius scales linearly with distance (`r_H = x * c`, with
    `c = (m / 3M)^(1/3)`), and the planet's periapsis `x (1 - e)` is what
    must clear the belt (GEN.727), so this solves `x (1 - e) - 5*c*x >=
    edge + gap` in closed form: `x >= (edge + gap) / (1 - e - 5c)`. `5c`
    is clamped to 0.9, and the divisor to 0.1, so a body heavy or
    eccentric enough to make that unsolvable still gets a large but
    finite distance.
    """
    c = _km_to_au(planet.hill_radius) / planet.distance
    k = min(5 * c, 0.9)
    return (belt.upper_limit + tuning.MIN_ASTEROID_BELT_SEPARATION) / max(1 - k - (planet.eccentricity or 0.0), 0.1)


def mutual_min_distance_au(planet, last_planet):
    """
    The closest `planet` (the later/farther of the pair) can stably
    sit to `last_planet` (fixed -- it's already been placed and
    won't move again this pass), via their *mutual* Hill radius
    (`orbits.mutual_hill_radius_m`) rather than either one's own
    individual Hill radius alone -- see
    `tuning.MUTUAL_HILL_RADII_SEPARATION`'s docstring for
    the stability-literature basis.

    Solved in closed form for `planet`'s own distance rather than
    evaluated once at the current (pre-correction) positions and
    added on top: the mutual Hill radius depends on the *average* of
    the two distances, so a naive "compute now, then push `planet`
    out by that much" approximation understates the requirement --
    moving `planet` out raises the average, which raises the
    requirement further, and `MUTUAL_HILL_RADII_SEPARATION`'s margin
    (10x) is large enough that this understatement is measurable, not
    just a rounding-level slip (confirmed by a real
    `assert_no_orbital_overlap` failure before this closed form
    replaced the naive version).

    Derivation: let d = `last_planet.distance` (fixed this pass), x =
    `planet`'s target distance, and kappa =
    `MUTUAL_HILL_RADII_SEPARATION * ((m_planet + m_last) / (3 *
    M_star)) ** (1/3)`. The stability requirement `x - d >= kappa *
    (x + d) / 2` rearranges to `x >= d * (1 + kappa/2) / (1 -
    kappa/2)`. With eccentric orbits (GEN.727) the gap that must hold the
    separation is the one between the inner planet's apoapsis and the
    outer one's periapsis, `x (1 - e2) - d (1 + e1) >= kappa (x + d) / 2`,
    which gives `x >= d * (1 + e1 + kappa/2) / (1 - e2 - kappa/2)` (the
    same at e = 0; the divisor is held at 0.1 or more).

    `M_star` is `planet.star.mass` rather than the system's `star` --
    for a single star or a P-type (close) binary's merged proxy both
    planets' own `.star` already *is* the system's star, but for an
    S-type (wide) binary's secondary list the system's `star` is still
    the *primary* -- using it here would derive the secondary's own
    planet spacing against the wrong star's mass entirely.

    Args:
        planet (Planet): The planet whose distance may need raising.
        last_planet (Planet): The fixed reference planet, closer to
                              the star.

    Returns:
        float: The minimum distance (AU) `planet` can sit at,
              measured from the star -- not a gap.
    """
    kappa = tuning.MUTUAL_HILL_RADII_SEPARATION * (
        (planet.mass + last_planet.mass) / (3 * planet.star.mass)
    ) ** (1 / 3)
    # A pair whose combined mass is a large enough fraction of the
    # star's own that kappa approaches/exceeds 2 has no finite stable
    # separation under this linear model at all -- clamp well short
    # of that so the formula below always returns a large-but-finite,
    # rather than negative or infinite, distance.
    kappa = min(kappa, 1.8)
    inner_e = last_planet.eccentricity or 0.0
    outer_e = planet.eccentricity or 0.0
    return last_planet.distance * (1 + inner_e + kappa / 2) / max(1 - outer_e - kappa / 2, 0.1)


def min_distance_after_au(body, last):
    """
    The closest `body` can sit outside `last` (its inner neighbor), in AU,
    by the same rules `space_orbits` enforces: for a belt, the inner edge
    (`lower_limit`) against `last`'s outer edge; for a planet, its orbit.
    """
    if body.body_type == 'a':
        if last.body_type == 'a':
            return last.upper_limit + tuning.MIN_ASTEROID_BELT_SEPARATION
        return apoapsis_au(last) + last.min_orbit_distance
    if last.body_type == 'a':
        return min_distance_past_belt_au(body, last)
    return mutual_min_distance_au(body, last)


# ---------------------------------------------------------------------
# Spacing pass (generation's own)
# ---------------------------------------------------------------------

def reconcile_moved_planet(planet):
    """
    Call right after `space_orbits` moves a top-level planet's `distance` to resolve an orbital overlap.

    The body keeps its mass, radius and composition (where it sits does not change what it is made of), but its climate is the
    one of its actual final position: its zone is re-derived from the new `distance` (`planetPhysics.reconcile_zone`), its air,
    water and temperature are made again there (`calculate_atmospheric_conditions`) and its class follows (`apply_class`, derived from
    those conditions), the same thing real orbital migration does to a body's conditions. Its Hill sphere, period, position and
    spin follow the distance too.

    Every moon is reconciled the same way: a moon's own zone is always its parent's (see `planetPhysics.generate_moons`'
    `zone_override`), and its climate depends on the parent's distance (its `distance_override`), so a moon whose planet moved
    needs the identical treatment even though the moon's own orbit around its planet never changed.

    Args:
        planet (Planet): The planet that moved.
    """
    planetPhysics.reconcile_zone(planet, planet.star.mass)
    # The Hill sphere grows with distance, so a moved planet's own clearance (and the next body's spacing) must follow it.
    planetPhysics.update_hill_sphere(planet)
    planetPhysics.calculate_atmospheric_conditions(planet)
    planetPhysics.apply_class(planet)
    planet.period = planetPhysics.calculate_orbital_period_years(planet.distance, planet.star.mass)
    planetPhysics.update_orbital_position(planet)
    planetPhysics.generate_spin(planet, planet.star.mass)

    for moon in planet.moons:
        planetPhysics.reconcile_zone(moon, planet.mass, distance_override=planet.distance)
        planetPhysics.calculate_atmospheric_conditions(moon, planet.distance, planet.eccentricity)
        planetPhysics.apply_class(moon)
        # GEN.704: the moon's star distance changed, so the dose, UV and score `generate_spin` ends with are re-run as a moved planet's are.
        planetPhysics.generate_spin(moon, planet.mass)

    # v20: planet.reflex_offset_x/y/z (this planet's own wobble from its moons -- see Planet.reflex_offset_x's docstring) depends on
    # every moon's own mass/position, any of which the reconciliation above may just have moved. Recomputed unconditionally here,
    # the same "cheap, always recompute from current state" treatment store.advance_orbital_phases gives this same value.
    refresh_planet_reflex(planet)

    # GEN.904: the spin, the orbital period or the class above can be new, and the day and night shares were worked out from the
    # old ones (or went with a regenerated body): they follow the geometry the planet and its moons have now (no draws).
    daylight.generate(planet)
    for moon in planet.moons:
        daylight.generate(moon, planet)
    # GEN.1034: and the tide its sea raises, from the orbits and spins as they are now
    ocean_tides.generate(planet)
    for moon in planet.moons:
        ocean_tides.generate(moon, planet)


def refresh_planet_reflex(planet, force=False):
    """Recomputes a planet's wobble from its moons (`reflex_offset_*`);
    with no moons it is zeroed when `force` is set (they were just
    removed) and otherwise left alone."""
    if planet.moons:
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons]
        )
    elif force:
        planet.reflex_offset_x = planet.reflex_offset_y = planet.reflex_offset_z = 0.0


def _shift_belt(belt, amount):
    belt.distance += amount
    belt.upper_limit += amount
    belt.lower_limit += amount


def space_orbits(planets):
    """
    Validates and adjusts the distances of `planets` (one star's ordered
    planets and asteroid belts) so no two orbits are too close, pushing
    each body outward past its inner neighbor where needed. This is
    generation's spacing pass (`StarSystem.validate_system` calls it).

    Two adjacent planets' minimum separation is enforced via their
    *mutual* Hill radius (`mutual_min_distance_au`), not either one's own
    individual Hill radius alone -- see
    `tuning.MUTUAL_HILL_RADII_SEPARATION`'s docstring for why.
    A belt has no mass/Hill-radius concept of its own, so any correction
    involving one uses the single real planet's own `min_orbit_distance`
    (5 Hill radii) on whichever side of the belt the planet is
    (`min_distance_past_belt_au` for a planet after a belt), or the fixed
    `MIN_ASTEROID_BELT_SEPARATION` between two belts.

    A moved planet is reconciled to its new distance
    (`reconcile_moved_planet`): its climate and class are made again at the
    new position, along with its moons'.

    Args:
        planets (list): One star's `Planet`/`AsteroidBelt` list, in
            orbital order. Changed in place.

    Returns:
        list: The bodies that moved.
    """
    moved = []
    if len(planets) < 2:
        return moved

    for i in range(1, len(planets)):
        planet = planets[i]
        last_planet = planets[i - 1]
        before = planet.distance

        if last_planet.body_type == 'a':
            distance_to_last = planet.distance - last_planet.upper_limit
        else:
            distance_to_last = planet.distance - apoapsis_au(last_planet)

        if distance_to_last < 0:
            # abs(distance_to_last) alone is exactly enough to cancel
            # the negative gap and land `planet` right on
            # `last_planet`'s own edge (upper_limit for a belt,
            # distance for a planet) -- the real minimum separation
            # (MIN_ASTEROID_BELT_SEPARATION or min_orbit_distance) is
            # then added on top of THIS corrected position by the
            # branches below. An earlier version of this line also
            # added `last_planet.distance`, double-counting that
            # offset on top of the already-correct cancellation and
            # roughly doubling the corrected distance instead of
            # nudging it just past the obstruction.
            additional_correction = abs(distance_to_last)
        else:
            additional_correction = 0

        if planet.body_type == 'a':
            if last_planet.body_type == 'a':
                if _below(distance_to_last, tuning.MIN_ASTEROID_BELT_SEPARATION):
                    _shift_belt(planet, tuning.MIN_ASTEROID_BELT_SEPARATION + additional_correction)
                    draw.redrawn("belt.spacing_push")
            elif _below(distance_to_last, last_planet.min_orbit_distance):
                _shift_belt(planet, last_planet.min_orbit_distance + additional_correction)
                draw.redrawn("belt.spacing_push")
        else:
            if last_planet.body_type == 'a':
                min_distance = min_distance_past_belt_au(planet, last_planet)
            else:
                min_distance = mutual_min_distance_au(planet, last_planet)
            if _below(planet.distance, min_distance):
                planet.distance = min_distance
                reconcile_moved_planet(planet)
                draw.redrawn("planet.hill_spacing_push" if last_planet.body_type != 'a' else "planet.belt_spacing_push")
        if planet.distance != before:
            moved.append(planet)
    return moved


def trim_to_orbit_ceiling(planets, ceiling_au):
    """
    Removes trailing planets/belts from `planets` that `space_orbits`'
    own overlap correction pushed past `ceiling_au`.

    Generation already keeps every body it places within the ceiling,
    but `space_orbits` (called afterward, to resolve spacing collisions
    between neighbors) only ever pushes a body's `distance` further
    OUTWARD, never closer -- for the galactic-Hill-sphere ceiling that
    correction was never practically reachable (negligible next to a
    light-year-scale boundary), but a much tighter S-type `a_crit_au`
    ceiling makes it a real possibility. Trims from the end (outermost
    first) since an outward-only correction cannot newly push an
    *interior* body past the ceiling without its outer neighbor
    already having been pushed past it first.

    Args:
        planets (list): The list to trim in place.
        ceiling_au (float): The boundary to check against, in AU (see
            `orbit_ceiling_au`).

    Returns:
        list: The removed bodies, outermost first.
    """
    removed = []
    while planets and exceeds_orbit_ceiling(planets[-1], ceiling_au):
        removed.append(planets.pop())
        draw.redrawn("planet.orbit_ceiling_trim")
    return removed


def comet_apoapsis_au(comet):
    """How far from its star a comet's orbit reaches, in AU: its
    farthest point, `q * (1 + e) / (1 - e)`; infinite for a parabolic
    comet, which leaves the star for good."""
    if comet.orbit_type != "elliptical":
        return math.inf
    return comet.perihelion_distance_au * (1 + comet.eccentricity) / (1 - comet.eccentricity)


def comet_fits_ceiling(comet, ceiling_au):
    """True if `comet`'s whole orbit stays inside `ceiling_au`."""
    return comet_apoapsis_au(comet) <= ceiling_au * (1 + RELATIVE_TOLERANCE)


def trim_comets_to_orbit_ceiling(comets, ceiling_au):
    """
    Removes from `comets`, in place, every comet whose orbit leaves
    `ceiling_au`: past a wide binary's star's stable limit
    (`orbit_ceiling_au`) a comet is on the way to the companion star,
    where its path would cross the companion's planets (GEN.228).

    Returns:
        list: The removed comets.
    """
    removed = [comet for comet in comets if not comet_fits_ceiling(comet, ceiling_au)]
    comets[:] = [comet for comet in comets if comet_fits_ceiling(comet, ceiling_au)]
    return removed


def comet_lists(system):
    """`(star, comets)` for each star whose comets must stay inside its
    stable limit: both stars of a wide binary, none otherwise."""
    if getattr(system, "binary_type", None) != "wide":
        return []
    return [(system.primary_star, system.comets), (system.secondary_star, system.secondary_comets)]


def _cross_star_threshold_au(system, outer_p, outer_s):
    """The worst-case gap two wide-binary stars' outermost bodies must
    keep: Gladman's mutual Hill criterion for two planets, else the
    fixed belt separation."""
    if outer_p.body_type == 'a' or outer_s.body_type == 'a':
        return tuning.MIN_ASTEROID_BELT_SEPARATION
    central_mass_kg = system.primary_star.mass + system.secondary_star.mass
    r_h_mutual_au = mutual_hill_radius_au(
        outer_p.mass, outer_s.mass, outer_p.distance, outer_s.distance, central_mass_kg
    )
    return constants.GLADMAN_MUTUAL_HILL_STABILITY_FACTOR * r_h_mutual_au


def cross_star_clearance(system):
    """
    For an S-type (wide) binary, prunes the outermost body from whichever
    star's list gravitationally encroaches on the other star's own
    outermost body, until the worst-case gap between them (the two pointed
    straight at each other with the stars at their closest approach,
    `periapsis - edge_p - edge_s`) clears Gladman's mutual Hill criterion (or, with a belt
    involved, `MIN_ASTEROID_BELT_SEPARATION`). The star whose outermost
    body has the least room to its own `a_crit_au` loses it. See
    `StarSystem._validate_cross_star_clearance` for the full reasoning.

    Returns:
        list: The removed bodies. Empty unless `system.binary_type` is
             "wide" and both lists are non-empty.
    """
    removed = []
    if system.binary_type != "wide" or not system.planets or not system.secondary_planets:
        return removed

    closest_au = system.wide_binary.periapsis_au
    while system.planets and system.secondary_planets:
        outer_p = max(system.planets, key=edge_au)
        outer_s = max(system.secondary_planets, key=edge_au)
        worst_case_gap_au = closest_au - edge_au(outer_p) - edge_au(outer_s)
        if worst_case_gap_au >= _cross_star_threshold_au(system, outer_p, outer_s):
            break
        draw.redrawn("planet.cross_star_trim")
        primary_room = system.primary_star.a_crit_au - edge_au(outer_p)
        secondary_room = system.secondary_star.a_crit_au - edge_au(outer_s)
        if primary_room <= secondary_room:
            system.planets.remove(outer_p)
            removed.append(outer_p)
        else:
            system.secondary_planets.remove(outer_s)
            removed.append(outer_s)
    return removed


# ---------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------

def check_planet(body, parent=None):
    """
    Checks one planet or moon on its own: its class is the one its conditions give (GEN.689: every body has exactly one,
    derived by `planet_classes.classify`, so a stored class that is not that one is a stale or edited row), its zone matches where
    it orbits (a moon's is its parent's), and it has a mass and a radius.

    Args:
        body (Planet): The planet or moon.
        parent (Planet, optional): A moon's planet.

    Returns:
        list: `Problem`s, empty when it validates.
    """
    problems = []
    name = _label(body)
    if body.planet_class not in planet_classes.LETTERS:
        return [Problem(name, f"class {body.planet_class!r} is not a known planet class")]
    derived = planet_classes.classify(body)
    if derived != body.planet_class:
        problems.append(Problem(name, f"its class is {body.planet_class} but its conditions give class {derived}"))
    distance = parent.distance if parent is not None else body.distance
    zone = planetPhysics.zone_for(body.habitable_zone, distance)
    if body.zone != zone:
        problems.append(Problem(name, f"its zone is {body.zone!r} but it orbits in zone {zone!r}"))
    if body.radius is None or not body.radius > 0:
        problems.append(Problem(name, f"its radius ({body.radius} km) is not positive"))
    if body.mass is None or not body.mass > 0:
        problems.append(Problem(name, "it has no mass"))
    return problems


def max_moon_radius_km(planet):
    """The largest moon `planet` can hold: `radius / 10**(1/3)`, the cap
    `planetPhysics.generate_moons` uses."""
    return planet.radius / (10 ** (1 / 3))


def check_lunar_system(planet):
    """
    Checks a planet and its moons: the planet itself (`check_planet`),
    then every moon -- its own class and zone, that it is no larger than
    the planet can hold, that it orbits inside the planet's stable range
    (`planetPhysics.moon_orbit_bounds_km`), and that each moon clears the
    one inside it by that moon's own `min_orbit_distance`.

    Returns:
        list: `Problem`s, empty when it validates.
    """
    problems = check_planet(planet)
    if not planet.moons:
        return problems
    low_km, high_km = planetPhysics.moon_orbit_bounds_km(planet)
    largest_km = max_moon_radius_km(planet)
    last = None
    for moon in sorted(planet.moons, key=lambda m: m.distance):
        name = _label(moon)
        problems.extend(check_planet(moon, parent=planet))
        if moon.radius is not None and moon.radius > largest_km * (1 + RELATIVE_TOLERANCE):
            problems.append(Problem(name, f"it is too large for {_label(planet)} to hold "
                                          f"({moon.radius:.4g} km, at most {largest_km:.4g} km)"))
        distance_km = _au_to_km(moon.distance)
        if _below(distance_km, low_km) or _au_to_km(apoapsis_au(moon)) > high_km * (1 + RELATIVE_TOLERANCE):
            problems.append(Problem(name, f"its orbit ({distance_km:.4g} km) is outside {_label(planet)}'s stable "
                                          f"range for moons ({low_km:.4g}-{high_km:.4g} km)"))
        if last is not None and _below(moon.distance, last.distance + last.min_orbit_distance):
            problems.append(Problem(name, f"its orbit is too close to {_label(last)}'s"))
        last = moon
    return problems


def check_orbit_spacing(planets):
    """
    Checks one star's planets and belts: in orbital order, each clears
    its inner neighbor by the distance `space_orbits` would enforce
    (`min_distance_after_au`).

    Returns:
        list: `Problem`s, empty when it validates.
    """
    problems = []
    for last, body in zip(planets, planets[1:]):
        if _below(inner_edge_au(body), inner_edge_au(last)):
            problems.append(Problem(_label(body), f"it is listed after {_label(last)} but orbits inside it"))
            continue
        needed = min_distance_after_au(body, last)
        if _below(inner_edge_au(body), needed):
            problems.append(Problem(_label(body), f"its orbit is too close to {_label(last)}'s "
                                                  f"(at {inner_edge_au(body):.4g} AU, needs at least {needed:.4g} AU)"))
    return problems


def _member_hosts(system):
    """`[(star, planets)]` of the further stars of the tree, none for a system that has no such stars."""
    return system.member_hosts() if hasattr(system, "member_hosts") else []


def star_lists(system):
    """`(star, planets)` for each star that holds its own planet list:
    the system's star (a single star, or a close pair's merged proxy),
    then a wide pair's secondary, then each further star of the tree (GEN.332)."""
    lists = [(system.star, system.planets)]
    if getattr(system, "binary_type", None) == "wide":
        lists.append((system.secondary_star, system.secondary_planets))
    lists.extend(_member_hosts(system))
    return lists


def check_star_system(system):
    """
    Checks a whole star system: for each star, the spacing of its planets
    and belts (`check_orbit_spacing`), that none lies beyond the star's
    orbit ceiling (`orbit_ceiling_au`), and every planet with its moons
    (`check_lunar_system`); for a wide binary, that the two stars'
    outermost bodies keep clear of each other (`cross_star_clearance`'s
    rule).

    Returns:
        list: `Problem`s, empty when it validates.
    """
    problems = []
    for star, planets in star_lists(system):
        problems.extend(check_orbit_spacing(planets))
        ceiling = orbit_ceiling_au(star)
        for body in planets:
            if reach_au(body) > ceiling * (1 + RELATIVE_TOLERANCE):
                problems.append(Problem(_label(body), f"it orbits beyond the farthest stable orbit "
                                                      f"({ceiling:.4g} AU)"))
            if body.body_type != 'a':
                problems.extend(check_lunar_system(body))
    for star, comets in comet_lists(system):
        ceiling = orbit_ceiling_au(star)
        for comet in comets:
            if not comet_fits_ceiling(comet, ceiling):
                problems.append(Problem(_label(comet), f"its orbit reaches beyond the farthest stable orbit "
                                                       f"({ceiling:.4g} AU)"))
    if getattr(system, "binary_type", None) == "wide" and system.planets and system.secondary_planets:
        outer_p = max(system.planets, key=edge_au)
        outer_s = max(system.secondary_planets, key=edge_au)
        gap = system.wide_binary.periapsis_au - edge_au(outer_p) - edge_au(outer_s)
        needed = _cross_star_threshold_au(system, outer_p, outer_s)
        if _below(gap, needed):
            problems.append(Problem(_label(outer_p), f"it and {_label(outer_s)} come too close to each other "
                                                     f"across the two stars"))
    return problems


# ---------------------------------------------------------------------
# Stabilizing an edited system
# ---------------------------------------------------------------------

def refresh_moon_orbit(moon, planet):
    """A moon's orbit-derived values (period, position, speed, Hill
    sphere, climate) after its own distance or its planet's mass or
    distance changed."""
    moon.period = planetPhysics.calculate_orbital_period_years(moon.distance, planet.mass)
    # The same Hill sphere generation gives a moon (`Planet.__init__`), about its planet.
    moon.primary_mass_kg = planet.mass
    planetPhysics.update_hill_sphere(moon)
    planetPhysics.calculate_atmospheric_conditions(moon, planet.distance, planet.eccentricity)
    planetPhysics.generate_orbital_motion_properties(moon, planet.mass)


def stabilize_lunar_system(planet):
    """
    Re-spaces a planet's moons after the planet or one of its moons
    changed: in orbital order, each moon moves out (never in) until it is
    past the planet's innermost stable moon orbit and clears the moon
    inside it by that moon's `min_orbit_distance`. Every moon's orbit is
    then refreshed for the planet's current mass and distance. Nothing is
    removed: a moon pushed past the stable range, or too large for the
    planet, is reported by `check_lunar_system`.

    Returns:
        list: The moons that moved.
    """
    if not planet.moons:
        return []
    planet.moons.sort(key=lambda m: m.distance)
    low_km, _high_km = planetPhysics.moon_orbit_bounds_km(planet)
    floor_au = _km_to_au(low_km)
    moved = []
    last = None
    for moon in planet.moons:
        refresh_moon_orbit(moon, planet)
        needed = floor_au if last is None else max(floor_au, last.distance + last.min_orbit_distance)
        if moon.distance < needed:
            moon.distance = needed
            refresh_moon_orbit(moon, planet)
            moved.append(moon)
        last = moon
    refresh_planet_reflex(planet)
    return moved


def refresh_star_reflex(system):
    """Each star's wobble from the planets orbiting it (a close pair's is
    shared by the pair), as `StarSystem.__init__` sets it."""
    def children(planets):
        return [(p.mass, p.position_x, p.position_y, p.position_z) for p in planets if p.body_type != 'a']

    if system.binary_type == "close":
        system.binary_planetary_wobble_x, system.binary_planetary_wobble_y, system.binary_planetary_wobble_z = \
            calculate_reflex_offset(system.star.mass, children(system.planets))
        return
    star = system.primary_star
    star.reflex_offset_x, star.reflex_offset_y, star.reflex_offset_z = \
        calculate_reflex_offset(star.mass, children(system.planets))
    if system.binary_type == "wide":
        star = system.secondary_star
        star.reflex_offset_x, star.reflex_offset_y, star.reflex_offset_z = \
            calculate_reflex_offset(star.mass, children(system.secondary_planets))
    for star, planets in _member_hosts(system):
        star.reflex_offset_x, star.reflex_offset_y, star.reflex_offset_z = calculate_reflex_offset(star.mass, children(planets))


def stabilize_star_system(system, allow_removal=False):
    """
    Brings an edited star system back to a stable layout, from the inside
    out: every lunar system first (`stabilize_lunar_system`), then each
    star's planets and belts, in their stored order, re-spaced outward
    (`space_orbits`; a moved planet whose class no longer fits its new
    zone is reclassified unless it is `pinned`), then the orbit ceiling
    and, for a wide binary, the clearance between the two stars.

    Bodies are only moved, never removed, unless `allow_removal` is set:
    then the bodies left beyond a star's orbit ceiling or crowding the
    companion star are dropped, outermost first, as generation does
    (`trim_to_orbit_ceiling`, `cross_star_clearance`). Without it they
    stay and are reported in `problems`.

    Args:
        system (StarSystem): Changed in place.
        allow_removal (bool): Drop what can't fit instead of reporting it.

    Returns:
        StabilizeReport: What moved, was reclassified or removed, and the
        problems still left (`check_star_system`).
    """
    before = {}
    for _star, planets in star_lists(system):
        for body in planets:
            before[id(body)] = (body.distance, getattr(body, "planet_class", None))
            for moon in getattr(body, "moons", ()):
                before[id(moon)] = (moon.distance, moon.planet_class)

    removed = []
    for star, planets in star_lists(system):
        for body in planets:
            if body.body_type != 'a':
                stabilize_lunar_system(body)
        space_orbits(planets)
        if allow_removal:
            removed.extend(trim_to_orbit_ceiling(planets, orbit_ceiling_au(star)))
    if allow_removal:
        removed.extend(cross_star_clearance(system))
        for star, comets in comet_lists(system):
            removed.extend(trim_comets_to_orbit_ceiling(comets, orbit_ceiling_au(star)))
    refresh_star_reflex(system)

    moved, reclassified = [], []
    for _star, planets in star_lists(system):
        for body in planets:
            bodies = [body] + list(getattr(body, "moons", ()))
            for item in bodies:
                distance, planet_class = before.get(id(item), (None, None))
                if distance is not None and item.distance != distance:
                    moved.append(_label(item))
                if planet_class is not None and getattr(item, "planet_class", None) != planet_class:
                    reclassified.append(_label(item))
    return StabilizeReport(moved, reclassified, [_label(body) for body in removed], check_star_system(system))


def reapply_life(body):
    """A body's life class, evolutionary pace and timeline, worked out again after its class or climate
    changed (`planetLife.apply_life_data`). The class and the pace it had stay when the new climate still allows them
    (DB.100); only one that no longer fits is drawn again."""
    planetLife.apply_life_data(body, keep_stored=True)
