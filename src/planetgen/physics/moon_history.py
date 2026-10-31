# planetgen/physics/moon_history.py

"""
The history of a planet's moons under its tide, applied to the bodies (GEN.959)
===============================================================================

`planet_kind` and `generate_moons` give a planet's moons the distances they were *born* at (a regular moon's disc, a captured one's catch). This pass
runs the history that `moon_orbits.evolve` integrates from those distances and the planet's birth spin to the age of the star, with the star's tide and
the thermal tide of the air on the planet's spin in the same sum (`body_spin`, GEN.957), and puts the result on the bodies:

* **The planet's spin** is where the exchange leaves it: a planet that turns faster than its moons go round gives them angular momentum and slows (the
  Earth and the Moon), one the star has locked or that turns slower than the moons pulls them in, and a moon that falls takes its angular momentum to
  the planet's spin. A planet within `SYNCHRONOUS_TOLERANCE` of the star's orbital rate is locked to it; the resonance windows of `body_spin` are for a planet
  that has no moons (TM5).
* **Each moon** that stays is moved to the distance it reached (`planets.move_orbit`: period, Hill sphere, position and speed follow, the plane, the
  shape and where it is on its orbit do not), and its own spin is run again at that distance when it moved by `RESPIN_FRACTION` (a ten-thousandth) or more.
* **The lunar system stays one a generated one could be:** `_settle` pushes a moon out to clear the one inside it, and loses one below the innermost orbit (ring material) or whose apoapsis passes the stable range.
* **A moon that falls** inside its Roche limit is removed and its mass and the time it fell are kept (`fallen_moons`, `ring_origin.FallenMoon`s, which `rings.apply` reads
  through `ring_origin.ring_origin`: the share of it that stays in orbit is a ring, spread for the time since, GEN.1057); a moon dragged below the innermost orbit is kept the same way
  at the age of the planet (the history does not say when); one pushed out of the stable part of the Hill sphere is lost.
* **The planet's heat** reads the moons that stay: the power their tides dissipate in a rocky planet is part of `tidal_heating.Forcing` (`moons`), found
  from the stored distances, masses and spin by `interior_heat`, at the temperature of the mantle.

The pass draws nothing except the spin axis' precession angle that `spin.set_spin` places again, and it can be repeated: it starts from the bodies
`_born_moons` (kept on the planet, with each moon's `_birth_distance_au`) and the planet's `_birth_spin`. A planet restored from the database has neither
and is left as it is. A history no integrator can finish (`moon_orbits.IntegrationFailed`) leaves the moons as born (TM8). A planet cannot be spun past `spin.BREAKUP_FRACTION` of its breakup rate by the moons that fall (TM9). Stand-ins are TM5 to TM9 of `docs/design/research/gap-filling-plan.md` section 6.25.
"""

import math

from planetgen.physics import body_spin, constants, moon_orbits, ring_origin, spin, spin_torques
from planetgen.physics.orbits import calculate_reflex_offset

RESPIN_FRACTION = 1.0e-4
"""float: A moon whose distance changed by this share of itself or more has its spin history run again at the new distance (its period moves by one and a half times as
much, and a locked moon's turns per orbit must stay within `spin.SPIN_LOCK_TOLERANCE` of one). Source tier T6, tweak level L0."""
SYNCHRONOUS_TOLERANCE = 1.0e-3
"""float: A planet with moons that ends within this share of its star's orbital rate is locked to the star. Source tier T6, tweak level L0."""
RETROGRADE_INCLINATION_DEG = 90.0
"""float: A moon whose orbit is inclined more than this goes round against its planet's spin."""


def retrograde(moon):
    """Whether a moon goes round the other way from its planet's spin."""
    return (getattr(moon, "orbital_inclination_deg", 0.0) or 0.0) > RETROGRADE_INCLINATION_DEG


def density_kg_m3(body):
    """A body's mean density, kg/m3."""
    return body.mass / (4.0 / 3.0 * math.pi * (body.radius * 1000.0) ** 3)


def born_moons(planet):
    """The moons the planet was generated with, in order, each with its birth distance (`_birth_distance_au`); the current moons the first time."""
    born = getattr(planet, "_born_moons", None)
    if born is None:
        born = list(planet.moons)
        for moon in born:
            moon._birth_distance_au = moon.distance
        planet._born_moons = born
    return born


def _settle(planet, survivors):
    """
    The moons that stay, in order of distance, held to the rules a generated lunar system keeps (`validation.check_lunar_system`): a moon below the
    innermost orbit (inside the planet's air and the room the largest moon needs) is lost to drag and its mass is ring material, a moon is pushed out
    (never in) to clear the one inside it by that one's `min_orbit_distance`, and one whose apoapsis then lies beyond the stable range is lost (TM4: moons
    do not meet, so the tide's drift is the only thing that moves them, and this keeps them apart). Returns `(moons, masses)`, the masses (kg) of the moons lost to drag.
    """
    from planetgen.physics import planets
    low_km, high_km = planets.moon_orbit_bounds_km(planet)
    kept, dragged, last = [], [], None
    for moon in survivors:
        if moon.distance * constants.AU_TO_M / 1000.0 < low_km:
            dragged.append(moon.mass)
            continue
        if last is not None and moon.distance < last.distance + last.min_orbit_distance:
            planets.move_orbit(moon, last.distance + last.min_orbit_distance, planet.mass)
        if moon.distance * (1.0 + (moon.eccentricity or 0.0)) * constants.AU_TO_M / 1000.0 > high_km * (1.0 + 1.0e-9):
            continue
        kept.append(moon)
        last = moon
    return kept, dragged


def apply(planet, primary_mass_kg):
    """
    Runs the history of `planet`'s moons and its spin together and sets them on the bodies (the module docstring). Called by `planets.generate_spin` for a planet with
    moons, which then sets the spin state returned here. A planet with no birth spin on record (restored from the database) is left as it is.
    """
    from planetgen.physics import planets
    birth = getattr(planet, "_birth_spin", None)
    if birth is None:
        return None
    born = born_moons(planet)
    age_s = planet.star.age * 1.0e9 * constants.SECONDS_PER_YEAR
    semi_major_axis_m = planet.distance * constants.AU_TO_M
    body = body_spin.planet_of(planet)
    star = spin_torques.Raiser(primary_mass_kg, semi_major_axis_m, spin_torques.orbital_rate(primary_mass_kg, semi_major_axis_m))
    moons = [moon_orbits.Moon(m.mass, m._birth_distance_au * constants.AU_TO_M, density_kg_m3(m), retrograde(m)) for m in born]
    try:
        result = moon_orbits.evolve(body, star, body_spin.dissipation_of(planet, body), moons, birth[0], age_s, body_spin.thermal_tide_of(planet, star),
                                    density_kg_m3(planet), (planet.hill_radius or 0.0) * 1000.0 or None)
    except moon_orbits.IntegrationFailed:
        # TM8: no integrator finished (a very stiff planet): the moons stay as they were born and the planet's spin is the one it would have without them
        planet.fallen_moons = ()
        for moon in born:
            if abs(moon.distance / moon._birth_distance_au - 1.0) > 1.0e-12:
                planets.move_orbit(moon, moon._birth_distance_au, planet.mass)
        planet.moons[:] = born
        found = body_spin.evolve(planet, primary_mass_kg, semi_major_axis_m, birth[0], age_s, getattr(planet, "eccentricity", 0.0) or 0.0, capture_units=birth[2])
        planets.generate_spin(planet, primary_mass_kg, birth=birth, found=found)
        return None
    entry = {id(m): m.distance for m in born}
    survivors = []
    for moon, fate in zip(born, result.moons):
        if fate.fate == "stays":
            target_au = fate.distance_m / constants.AU_TO_M
            if abs(target_au / moon.distance - 1.0) > 1.0e-12:
                planets.move_orbit(moon, target_au, planet.mass)
            survivors.append(moon)
    survivors.sort(key=lambda m: m.distance)      # a moon that drifted past another keeps its place by distance
    survivors, dragged = _settle(planet, survivors)
    for moon in survivors:
        if abs(moon.distance / entry[id(moon)] - 1.0) >= RESPIN_FRACTION:
            planets.generate_spin(moon, planet.mass, birth=getattr(moon, "_birth_spin", None))
    planet.moons[:] = survivors
    age_yr = planet.star.age * 1.0e9
    planet.fallen_moons = tuple(
        [ring_origin.FallenMoon(fate.moon.mass_kg, fate.time_s / constants.SECONDS_PER_YEAR) for fate in result.moons if fate.fate == "fell"]
        + [ring_origin.FallenMoon(mass, age_yr) for mass in dragged])
    if survivors:
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in survivors])
    else:
        planet.reflex_offset_x = planet.reflex_offset_y = planet.reflex_offset_z = 0.0
    breakup = spin.BREAKUP_FRACTION * spin.critical_spin_rad_s(planet.mass, planet.radius)
    rate = math.copysign(min(abs(result.spin_rate), breakup), result.spin_rate)      # the angular momentum of falling moons cannot spin a planet past its breakup (TM9)
    if abs(rate - star.orbital_rate) <= SYNCHRONOUS_TOLERANCE * star.orbital_rate:
        found = spin_torques.SpinResult(star.orbital_rate, 1.0, "synchronous")
    else:
        found = spin_torques.SpinResult(rate, None, "free")
    planets.generate_spin(planet, primary_mass_kg, birth=birth, found=found)
    return None
