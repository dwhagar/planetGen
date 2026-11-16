# planetgen/physics/moon_history.py

"""
The history of a planet's moons under its tide, applied to the bodies (GEN.959)
===============================================================================

`planet_kind` and `generate_moons` give a planet's moons the distances they were *born* at (a regular moon's disk, a captured one's catch). This pass
runs the history that `moon_orbits.evolve` integrates from those distances and the planet's birth spin to the age of the star, with the star's tide and
the thermal tide of the air on the planet's spin in the same sum (`body_spin`, GEN.957), and puts the result on the bodies:

* **The planet's spin** is where the exchange leaves it: a planet that turns faster than its moons go round gives them angular momentum and slows (the
  Earth and the Moon), one the star has locked or that turns slower than the moons pulls them in, and a moon that falls takes its angular momentum to
  the planet's spin. A planet within `SYNCHRONOUS_TOLERANCE` of the star's orbital rate is locked to it; the resonance windows of `body_spin` are for a planet
  that has no moons (TM5).
* **Each moon** that stays is moved to the distance it reached (`planets.move_orbit`: period, Hill sphere, position and speed follow, the plane, the
  shape and where it is on its orbit do not), and its own spin is run again at that distance when it moved by `RESPIN_FRACTION` (a ten-thousandth) or more.
* **Two moons that swept a resonance** on the way out are held in it, the inner pushing the outer along (`moon_resonances.capture`, GEN.1449 slice 2); the eccentricity the resonance forces is given once the orbits are final.
* **A moon the evection resonance caught** (the planet's bulge turns its orbit once a year of the star's, `evection`, GEN.1750) is held at it as the planet slows and given the eccentricity it raises
  (`evection.held_eccentricity`, kept as a floor by `settle_orbit` like a moon pair's forced one); one the resonance only crossed is untouched.
* **The lunar system stays one a generated one could be:** `_settle` pushes a moon out to clear the one inside it, and loses one below the innermost orbit (ring material) or whose apoapsis passes the stable range.
* **A moon that falls** inside its Roche limit is removed and its mass and the time it fell are kept (`fallen_moons`, `ring_origin.FallenMoon`s, which `rings.apply` reads
  through `ring_origin.ring_origin`: the share of it that stays in orbit is a ring, spread for the time since, GEN.1057); a moon dragged below the innermost orbit is kept the same way
  at the age of the planet (the history does not say when); one pushed out of the stable part of the Hill sphere is lost.
* **The planet's heat** reads the moons that stay: the power their tides dissipate in a rocky planet is part of `tidal_heating.Forcing` (`moons`), found
  from the stored distances, masses and spin by `interior_heat`, at the temperature of the mantle.

The pass draws nothing except the spin axis' precession angle that `spin.set_spin` places again, and it can be repeated: it starts from the bodies
`_born_moons` (kept on the planet, with each moon's `_birth_distance_au`) and the planet's `_birth_spin`. A planet restored from the database has them from its rows (GEN.1753: `body_history`, the `lost_moons` table
and `restore_born`), the moons the history lost as `LostMoon`s; one saved before that starts from the moons it has, each at the orbit it has now. A history no integrator can finish (`moon_orbits.IntegrationFailed`) leaves the moons as born (TM8). A planet cannot be spun past `spin.BREAKUP_FRACTION` of its breakup rate by the moons that fall (TM9). Stand-ins are TM5 to TM9 of `docs/design/research/gap-filling-plan.md` section 6.25.
"""

import math

from planetgen.physics import body_spin, constants, evection, hill_perturbers, moon_orbits, moon_resonances, ring_origin, spin, spin_torques
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


class LostMoon:
    """
    A moon the history of a stored planet took away, kept for the history to run again (GEN.1753, `body_history.lost_rows`): the four numbers `moon_orbits.Moon` reads (its mass, birth orbit, density and
    direction) and nothing else. It is never a body of the system: `apply` lets it take part and drops it again, and one that the history now lets stay is counted with the moons that fell.

    Attributes:
        mass (float): kg.
        radius (float): km, from the density.
        distance (float): AU, the birth orbit; a moon that was lost is nowhere else.
        orbital_inclination_deg (float): 180 for a moon that went round against its planet's spin, else 0.
    """

    def __init__(self, mass_kg, birth_distance_au, density_kg_m3, retrograde_flag):
        self.mass = mass_kg
        self.radius = (3.0 * mass_kg / (4.0 * math.pi * density_kg_m3)) ** (1.0 / 3.0) / 1000.0
        self.distance = birth_distance_au
        self._birth_distance_au = birth_distance_au
        self.orbital_inclination_deg = 180.0 if retrograde_flag else 0.0


def restore_born(planet, kept, lost):
    """
    Gives a planet restored from the database the list of the moons it was born with (`_born_moons`) that a run keeps (GEN.1753): its moons that have a stored birth orbit and place, and the lost ones as `LostMoon`s,
    in the order they were born in. Nothing is set unless every moon the planet has comes with a place (a planet saved before the lost moons were kept starts its history from the moons it has, as before), and
    the planet needs its birth spin (`_birth_spin`) for the history to be run.

    Args:
        planet: A restored planet.
        kept (list[tuple]): `(born_index, moon)` for each moon of the planet that has a stored place.
        lost (list[dict]): The planet's `lost_moons` rows (`body_history.LOST_FIELDS`).

    Returns:
        bool: Whether the born list was set.
    """
    if getattr(planet, "_birth_spin", None) is None or len(kept) != len(getattr(planet, "moons", ())) or not kept and not lost:
        return False
    placed = {place: moon for place, moon in kept}
    for row in lost:
        placed[row["born_index"]] = LostMoon(row["mass_kg"], row["birth_distance_au"], row["density_kg_m3"], row["retrograde"])
    if len(placed) != len(kept) + len(lost) or sorted(placed) != list(range(len(placed))):
        return False
    planet._born_moons = [placed[place] for place in range(len(placed))]
    return True


def born_moons(planet):
    """The moons the planet was generated with, in order, each with its birth distance (`_birth_distance_au`); the current moons the first time."""
    born = getattr(planet, "_born_moons", None)
    if born is None:
        born = list(planet.moons)
        for moon in born:
            if getattr(moon, "_birth_distance_au", None) is None:      # a restored moon brings the orbit it was born on (GEN.1753)
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
    kept, dragged = [], []
    for moon in survivors:
        if moon.distance * constants.AU_TO_M / 1000.0 < low_km:
            dragged.append(moon.mass)
            continue
        clearance = max((inner.distance + inner.min_orbit_distance for inner in kept), default=0.0)       # every moon kept inside it, not the last alone (GEN.1460)
        if kept and moon.distance < clearance:
            planets.move_orbit(moon, clearance, planet.mass)
        if moon.distance * (1.0 + (moon.eccentricity or 0.0)) * constants.AU_TO_M / 1000.0 > high_km * (1.0 + 1.0e-9):
            continue
        kept.append(moon)
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
                                    density_kg_m3(planet), (hill_perturbers.tidal_hill_radius_km(planet) or 0.0) * 1000.0 or None,      # GEN.1449: against every body that pulls on the planet
                                    resonance=evection.Resonance(getattr(planet, "body_type", None) == "g", star.orbital_rate))
    except moon_orbits.IntegrationFailed:
        # TM8: no integrator finished (a very stiff planet): the moons stay as they were born and the planet's spin is the one it would have without them
        planet.fallen_moons = ()
        born = [moon for moon in born if not isinstance(moon, LostMoon)]      # a stand-in is no body of the system
        for moon in born:
            if abs(moon.distance / moon._birth_distance_au - 1.0) > 1.0e-12:
                planets.move_orbit(moon, moon._birth_distance_au, planet.mass)
        planet.moons[:] = born
        found = body_spin.evolve(planet, primary_mass_kg, semi_major_axis_m, birth[0], age_s, getattr(planet, "eccentricity", 0.0) or 0.0, capture_units=birth[2])
        planets.generate_spin(planet, primary_mass_kg, birth=birth, found=found)
        return None
    entry = {id(m): m.distance for m in born}
    survivors = []
    stayed = []      # lost moons of a stored planet that the history, run again, lets stay: the system has no body for them, so they are counted as lost
    for moon, fate in zip(born, result.moons):
        if fate.fate == "stays" and isinstance(moon, LostMoon):
            stayed.append(moon.mass)
        elif fate.fate == "stays":
            target_au = fate.distance_m / constants.AU_TO_M
            if abs(target_au / moon.distance - 1.0) > 1.0e-12:
                planets.move_orbit(moon, target_au, planet.mass)
            survivors.append(moon)
            if fate.held:      # GEN.1750: caught in the evection resonance, held there with the eccentricity it raises (`evection.held_eccentricity`)
                held_e = evection.held_eccentricity(fate.distance_m, planet.radius * 1000.0)
                if held_e - (moon.eccentricity or 0.0) > planets.ECCENTRICITY_SETTLE_TOLERANCE:
                    moon.eccentricity = held_e
                    planets.update_orbital_position(moon)
                    moon._resonance_e = held_e
                    if getattr(moon, "_age_drawn", None) is not None:
                        moon._age_drawn = 0.0      # its air, class and heat are made again once, on the orbit it ends with (`planets.stale_air`)
    survivors.sort(key=lambda m: m.distance)      # a moon that drifted past another keeps its place by distance
    if len(survivors) > 1:        # GEN.1449 slice 2: a pair that swept a first-order resonance on the way out is held in it (`moon_resonances.capture`)
        held, _ = moon_resonances.capture([m.mass for m in survivors], [m._birth_distance_au for m in survivors], [m.distance for m in survivors])
        for moon, target_au in zip(survivors, held):
            if abs(target_au / moon.distance - 1.0) > 1.0e-12:
                planets.move_orbit(moon, target_au, planet.mass)
    survivors, dragged = _settle(planet, survivors)
    for moon in survivors:
        if abs(moon.distance / entry[id(moon)] - 1.0) >= RESPIN_FRACTION and not planets.is_bare(moon):     # a moon built bare has no spin yet, it takes its first at its final orbit (GEN.1548)
            planets.generate_spin(moon, planet.mass, birth=getattr(moon, "_birth_spin", None))
    planet.moons[:] = survivors
    age_yr = planet.star.age * 1.0e9
    planet.fallen_moons = tuple(
        [ring_origin.FallenMoon(fate.moon.mass_kg, fate.time_s / constants.SECONDS_PER_YEAR) for fate in result.moons if fate.fate == "fell"]
        + [ring_origin.FallenMoon(mass, age_yr) for mass in dragged + stayed])
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
