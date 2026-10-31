# planetgen/physics/planets.py

"""
Planet Physical and Orbital Generation
=======================================

This module determines what a planet or moon physically *is*: its orbital
zone, kind, mass, radius, density, composition, atmosphere, surface
gravity, and atmospheric conditions, plus the generation of any moons it
has. No class letter is read: the mass comes from the formation law
(`planet_kind.draw_formation`) or a moon's size law, the composition from the
star and where the body formed, and the class is derived from the finished
body (`apply_class`, `planet_classes.classify`). It operates on `Planet` instances (from `planet`) passed in as
arguments, rather than as methods on the class, so that generation logic is
kept separate from the `Planet` object's own definition and presentation.

Life chemistry and evolutionary data are intentionally NOT determined here —
see `life.apply_life_data`, which `StarSystem` applies to every planet
and moon in a system after all of them have been generated.
"""

import math
import re

from planetgen.physics import moon_history, body_state, bulk_draw, constants, giant_draw, giant_interior, giant_tilt, kepler, metallicity, moon_composition, orbit_elements, planet_classes, planet_kind, planet_mass, rings
from planetgen.physics import aerosols, atmosphere, body_colour, body_spin, daylight, spin_torques, climate_draws, habitability_world, hosting_gate, hydrosphere, interior_heat, magnetism, ocean_tides, radiation, radiogenic, relations, spin, surface_air, surface_composition, tidal_spin
from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.physics.orbits import calculate_hill_sphere, calculate_reflex_offset, minimum_update_interval_years
from planetgen.physics.state_vectors import state_from_elements
from planetgen.util.checks import finite_domain
from planetgen.util.random import sample_bounded_bell, sample_power_law
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_msun = units.converter("kg", "Msun")
_mearth_to_kg = units.converter("Mearth", "kg")

_km_to_m = units.converter("km", "m")
_m_to_km = units.converter("m", "km")
_au_to_m = units.converter("au", "m")
_au_to_km = units.converter("au", "km")
_km_to_au = units.converter("km", "au")
_gcm3_to_kgm3 = units.converter("g_cm3", "kg_m3")
_h_to_s = units.converter("h", "s")
_yr_to_s = units.converter("yr", "s")
_gyr_to_yr = units.converter("Gyr", "yr")
_yr_to_h = units.converter("yr", "h")
_au_yr_to_kms = units.converter("au_yr", "km_s")
_TWO_PI_RAD = 2.0 * math.pi


def uses_giant_mass_radius(planet):
    """Whether `planet` gets its mass and radius from the giant-planet mass-radius relation (GEN.34): a giant."""
    return planet.body_type == "g"


def _giant_transition_earth():
    """`(mass, radius)` in Earth units where the Neptunian and Jovian
    branches of the giant mass-radius relation meet."""
    m_n, r_n, s_n = constants.GIANT_NEPTUNIAN_MASS_RADIUS
    m_j, r_j, s_j = constants.GIANT_JOVIAN_MASS_RADIUS
    # r_n * (M/m_n)^s_n == r_j * (M/m_j)^s_j, solved in log space.
    log_m = (math.log(r_j / r_n) + s_n * math.log(m_n) - s_j * math.log(m_j)) / (s_n - s_j)
    mass = math.exp(log_m)
    return mass, r_n * (mass / m_n) ** s_n


GIANT_TRANSITION_MASS_EARTH, GIANT_TRANSITION_RADIUS_EARTH = _giant_transition_earth()


def giant_regime(mass_kg):
    """`"neptunian"` below the relation's transition mass, else `"jovian"`."""
    mass_earth = _kg_to_mearth(mass_kg)
    return "neptunian" if mass_earth < GIANT_TRANSITION_MASS_EARTH else "jovian"


def giant_radius_km(mass_kg):
    """The giant mass-radius relation's median radius (km) for `mass_kg`
    (`constants.GIANT_NEPTUNIAN_MASS_RADIUS`), with no scatter."""
    mass_earth = _kg_to_mearth(mass_kg)
    if giant_regime(mass_kg) == "neptunian":
        m0, r0, slope = constants.GIANT_NEPTUNIAN_MASS_RADIUS
    else:
        m0, r0, slope = constants.GIANT_JOVIAN_MASS_RADIUS
    return r0 * (mass_earth / m0) ** slope * constants.EARTH_RADIUS_KM


def _sample_giant_mass_kg(low_kg, high_kg):
    """A giant's mass in [low_kg, high_kg], from dN/dlogM ~
    M^-GIANT_MASS_FUNCTION_SLOPE."""
    return sample_power_law(low_kg, high_kg, tuning.GIANT_MASS_FUNCTION_SLOPE)


def sample_giant_radius_km(mass_kg):
    """`giant_radius_km(mass_kg)` with the relation's own scatter
    (`GIANT_RADIUS_SCATTER`, cut at 3 sigma). Bound giants
    (`_sample_giant_radius_km`) and rogue gas giants (GEN.60) both draw
    their radius here."""
    sigma = constants.GIANT_RADIUS_SCATTER[giant_regime(mass_kg)]
    factor = 1 + max(-3.0, min(3.0, draw.gauss(0.0, 1.0))) * sigma
    return giant_radius_km(mass_kg) * factor


def _giant_mass_for_radius_kg(radius_km, star_mass_kg):
    """
    A mass (kg) for a giant whose radius was fixed first (a caller's radius): the Neptunian branch's inverse below the transition
    radius, never under the core-critical mass, else a draw from the runaway mass function (the Jovian branch is nearly flat, so
    radius barely constrains mass there), never under the transition mass.
    """
    transition_kg = _mearth_to_kg(GIANT_TRANSITION_MASS_EARTH)
    radius_earth = radius_km / constants.EARTH_RADIUS_KM
    if radius_earth < GIANT_TRANSITION_RADIUS_EARTH:
        m0, r0, slope = constants.GIANT_NEPTUNIAN_MASS_RADIUS
        mass_kg = m0 * (radius_earth / r0) ** (1 / slope) * constants.EARTH_MASS_TO_KG
        return max(mass_kg, _mearth_to_kg(tuning.GIANT_CORE_CRITICAL_MASS_EARTH))
    return max(transition_kg, planet_mass.runaway_mass_kg(star_mass_kg))


def apply_class(planet):
    """
    Sets the class of a finished planet or moon from its stored conditions (GEN.689, `planet_classes.classify`), with the description
    text of that class (the composition text is written from the stored fractions, GEN.430). Called when the body is built (so a
    request can test it) and again once the system is final (the seasonal swing of its starlight, `flux_variation`, is known only then).
    """
    planet.planet_class = planet_classes.classify(planet)
    data = tuning.PLANET_CLASSES[planet.planet_class]
    planet.description = data["description"]
    if planet.is_moon:
        # Class descriptions are authored generically (e.g. "shields the inner planets"); a moon knows it's a moon.
        planet.description = re.sub(r'\bplanets\b', 'moons', planet.description)
        planet.description = re.sub(r'\bplanet\b', 'moon', planet.description)


@finite_domain()
def calculate_orbital_period_years(distance_au, primary_mass_kg):
    """
    Kepler's third law: T(years) = sqrt(a(AU)^3 / M_primary(Msun)).

    "Primary" is whatever body this orbit is actually around -- the host
    star for an ordinary planet, but the parent planet for a moon (see
    `Planet.__init__`'s `primary_mass_kg` parameter). A single shared
    helper so this formula is computed the same way wherever an orbital
    distance is set or changed -- also `StarSystem.validate_system`, which
    adjusts `distance` after generation to resolve overlapping orbits and
    must keep `period` in sync with it.

    Args:
        distance_au (float): Orbital distance from the primary, in AU.
        primary_mass_kg (float): The primary's mass, in kg.

    Returns:
        float: Orbital period in years.

    Raises:
        ValueError: If `distance_au` or `primary_mass_kg` isn't positive --
            every real caller always has a positive orbital distance and
            primary mass, so this guards against a `ZeroDivisionError`/
            `math domain error` escaping from a coding mistake upstream
            (e.g. an orbital-overlap correction pushing `distance` to zero
            or negative) with the same clean-`ValueError` treatment
            `calculate_surface_gravity` already gives a non-positive
            gravity, instead of a confusing low-level exception.
    """
    if distance_au <= 0 or primary_mass_kg <= 0:
        raise ValueError(
            f"calculate_orbital_period_years: distance_au ({distance_au}) and primary_mass_kg "
            f"({primary_mass_kg}) must both be positive."
        )
    primary_mass_sol = _kg_to_msun(primary_mass_kg)
    return math.sqrt(distance_au ** 3 / primary_mass_sol)


def formation_distance_au(planet):
    """Where the body formed, in AU: its own distance, or its parent's for a moon (`Planet.formation_distance_au`)."""
    distance = getattr(planet, "formation_distance_au", None)
    return planet.distance if distance is None else distance


def formed_beyond_snow_line(planet):
    """Whether `planet` formed at or beyond its star's snow line (False for a star with none)."""
    line = getattr(getattr(planet, "star", None), "snow_line_au", None)
    return line is not None and formation_distance_au(planet) >= line


def generate_planet_properties(planet, zone_override=None):
    """
    Generates a planet's physical and orbital properties: zone, kind, mass, radius, density, composition, and the interior
    (GEN.646, no class letter is read).

    A planet draws its kind and mass from the formation law at its distance (`planet_kind.draw_formation`); a moon's mass is
    given by its parent (`generate_moons`, the moon size law); a caller may give a mass, a radius or both (kept). A terrestrial
    body then draws its bulk composition from its star and its origin (`bulk_draw.draw_terrestrial`) and a giant its heavy
    elements and radius (`giant_draw.draw_giant`). The air, the water and the class follow (`calculate_atmospheric_conditions`,
    `apply_class`).

    Args:
        planet (Planet): The planet to generate properties for. Any of `planet.radius`/`planet.mass` already set are kept; unset
                         ones are filled in.
        zone_override (str, optional): A character ('h', 'c', 'e') to manually set the planet's zone, overriding the
                                       calculation based on distance.
    """
    # Determine the planet's zone (hot, ecosphere, or cold)
    zone = zone_for(planet.habitable_zone, planet.distance)

    if zone_override and zone_override.lower() in "hce":
        zone = zone_override.lower()
    planet.zone = zone

    radius_given = planet.radius is not None
    mass_given = planet.mass is not None
    beyond = formed_beyond_snow_line(planet)
    _check_given_bounds(planet.mass if mass_given else None, planet.radius if radius_given else None)

    if planet.is_moon:
        kind = "terrestrial"
        if not mass_given and not radius_given:
            raise ValueError("A moon is built from a mass or a radius its parent gives (generate_moons).")
    elif mass_given or radius_given:
        kind = planet_kind.kind_of_given(planet.mass if mass_given else None, planet.radius if radius_given else None, beyond)
    else:
        formation = planet_kind.draw_formation(planet.distance, planet.star.mass, planet.star.snow_line_au,
                                               _star_fe_h(planet.star), bulk_draw.metal_ratio_to_solar(planet.star))
        kind = formation.kind
        planet.mass = max(formation.mass_kg, tuning.PLANET_MIN_MASS_KG)
        mass_given = True
        log.choice("Planet kind", kind, f"formation at {planet.distance:.4g} AU gave {planet.mass:.4g} kg "
                   f"(core {formation.core_kg:.4g} kg, runaway {formation.runaway}, beyond the snow line {formation.beyond_snow_line})")

    if mass_given and radius_given:
        _check_given_density(planet.mass, planet.radius, kind)

    planet.body_type = "g" if kind == "giant" else "t"
    rings.clear(planet)  # a new draw draws its ring again (GEN.732)

    if uses_giant_mass_radius(planet):
        apply_giant(planet, radius_given, mass_given)
    else:
        # The composition comes first and the radius, density and gravity follow from it (GEN.428, GEN.429; `bulk_draw`).
        bulk_draw.draw_terrestrial(planet, radius_given, mass_given)
        giant_interior.clear(planet)

    # The air is made later, with the temperature it gives (`calculate_atmospheric_conditions`, GEN.411 and GEN.422): whether the
    # body holds one is the solver's call. Clear what a previous draw left behind -- `reconcile_zone` regenerates a body
    # in place -- and the class, which is derived once the body is finished (`apply_class`).
    planet.atmosphere = None
    planet.atm_density = None
    planet.atm_molar_density = None
    planet.scale_height = None
    planet.planet_class = None
    planet.description = None

    planet.volume = (4 / 3) * math.pi * planet.radius ** 3
    # GEN.459: the potassium, thorium, uranium and aluminium-26 the interior heat and the ground's dose both read (a stored draw of the
    # body, so the later recalculation of its air and water, which moves nothing in its rock, does not redraw it).
    # GEN.428: Th and U go with the rock's Al2O3 (part E 2.8b), the body's over the bulk silicate Earth's.
    radiogenic.generate_content(planet, refractory_enrichment=radiogenic.rock_of(planet).aluminium / radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    # GEN.699: the rheology of the mantle that the heat it was born with leaves through (the heat itself: `internal_heat_flux_w_m2`).
    interior_heat.generate_interior(planet)

    update_hill_sphere(planet)


ECCENTRICITY_SETTLE_TOLERANCE = 1e-4
"""float: A planet or captured moon whose eccentricity at the star's final age is within this of the one drawn keeps its orbit (`settle_orbit`);
below it the orbit cannot be told apart. A stand-in, registry row OE14."""


def settle_orbit(planet):
    """
    Takes the tides of a planet's (or captured moon's) orbit again at the star's final age (GEN.985): `generate_orbit_shape` wore the
    birth eccentricity down for the age the star had when the planet was drawn, and the age rises when the planets are in. The birth
    eccentricity is kept with the body for the generation (`_orbit_birth`, not stored), and the circularisation time is read again from the
    body's radius as it is now. Ages only rise, so the eccentricity only falls (it is never raised above what the body has, which the spacing of the orbits may have cut), and a lower one clears the star and the neighbours better.
    The position and velocity are made again on the new ellipse (`update_orbital_position`, the body keeps its place in its plane) and
    a planet's moons are told (`host_eccentricity`). The air, the tide's heat and the spin of a rocky body are read again by the
    re-finish of GEN.992, not here; the star's XUV flux, exposure and flares at the new mean light are set here, so the stored flux follows the stored orbit.

    Returns:
        bool: Whether the eccentricity moved by more than `ECCENTRICITY_SETTLE_TOLERANCE`.
    """
    birth = getattr(planet, "_orbit_birth", None)
    if birth is None:
        return False
    time_yr = orbit_elements.circularisation_time_yr(planet, planet.primary_mass_kg, planet.distance)
    # never above what the planet has: the spacing of the orbits may have cut it (`validation`) and ages only rise
    eccentricity = min(planet.eccentricity, orbit_elements.circularised(birth, _gyr_to_yr(planet.star.age), time_yr))
    if planet.eccentricity - eccentricity <= ECCENTRICITY_SETTLE_TOLERANCE:
        return False
    planet.eccentricity = eccentricity
    for moon in getattr(planet, "moons", ()):
        moon.host_eccentricity = eccentricity
    update_orbital_position(planet)
    # the star's XUV and flares at the mean light of the new orbit (no draws), so the stored flux follows the stored ellipse; the air waits for the re-finish
    light_au = orbit_elements.flux_distance_au(planet.distance, eccentricity)
    for body in (planet, *getattr(planet, "moons", ())):
        if getattr(body, "_star_distance_au", None):
            magnetism.update_exposure(body, light_au)
    return True


AGE_SETTLE_TOLERANCE = 0.02
"""float: A rocky planet or moon made when the star's age was within this share of its final age keeps its air, class and dose
(`refinish_rocky`); the age the decay clock and the heat history read is then within 2 percent of the right one. A stand-in, registry
row OE15."""


def age_moved(body):
    """Whether the star's age has risen by more than `AGE_SETTLE_TOLERANCE` of the age `body` was made at."""
    drawn = getattr(body, "_age_drawn", None)
    return drawn is not None and body.star.age - drawn > AGE_SETTLE_TOLERANCE * max(drawn, 0.1)


def refinish_rocky(planet):
    """
    Makes a rocky planet's end stages again at the star's final age, and its moons' (GEN.992, found by GEN.986): the air, the water and the
    surface (`calculate_atmospheric_conditions`), the class, the spin at that age (its own birth spin worn by the tide again, `generate_spin` with `birth`, GEN.961) and what ends it (the field, the dose, the colour and the
    habitability score; `generate_spin`), the same sequence a moved planet gets (`validation.reconcile_moved_planet`), but from a stream of
    the body's own (mass, distance, star and the age it settles to), so the run's stream does not move. A body keeps its mass, radius,
    composition and orbit. The forced requests are judged on the classes this makes (`StarSystem`), and the moons' induced fields, rings and
    life follow in later passes. A planet or moon whose star's age rose by under `AGE_SETTLE_TOLERANCE` of the age it was made at is left.

    Args:
        planet (Planet): A rocky planet (`body_type` "t", not a moon) with its moons; a giant's moons are made again too.

    Returns:
        int: How many bodies were made again (the planet and its moons).
    """
    made = 0
    with_planet = planet.body_type == "t" and not planet.is_moon and age_moved(planet)

    def bodies():
        if with_planet:
            yield planet
        yield from list(planet.moons)      # read after the planet is made again: its moons' history (GEN.959) can replace a moon with a new body

    for body in bodies():
        if not age_moved(body):
            continue
        stream = draw.Stream(f"refinish-rocky:{float(body.mass)!r}:{float(body.distance)!r}:{float(planet.distance)!r}:"
                             f"{float(body.star.mass)!r}:{float(body.star.age)!r}:{float(body.radius)!r}:{int(body.is_moon)}")
        with draw.bound(stream):
            if body.is_moon:
                calculate_atmospheric_conditions(body, planet.distance, planet.eccentricity)
                apply_class(body)
                generate_spin(body, planet.mass, birth=getattr(body, "_birth_spin", None))
            else:
                calculate_atmospheric_conditions(body)
                apply_class(body)
                generate_spin(body, body.primary_mass_kg, birth=getattr(body, "_birth_spin", None))
        body._age_drawn = body.star.age
        made += 1
    if made:  # GEN.767: the spin and period changed, the moons' day and night shares follow (no draws)
        for moon in planet.moons:
            daylight.generate(moon, planet)
    if made:
        # GEN.1034: the sea and the spin are new, and the tide worked out from the old ones follows them (no draws)
        ocean_tides.generate(planet)
        for moon in planet.moons:
            ocean_tides.generate(moon, planet)
    return made


def settle_giant(planet):
    """
    Settles a giant sized from the radius tables once the star's age and the orbit are final (GEN.961, `giant_draw.settle`): when its
    radius moves, the volume and the gravity follow, and the air is made again from a stream of its own (the gravity sets its scale
    height; the run's stream is not touched), then the class and the spin (which ends with the field, the dose, the colour and the habitability score).

    Returns:
        bool: Whether the giant's radius moved.
    """
    if not giant_draw.settle(planet):
        return False
    planet.volume = (4 / 3) * math.pi * planet.radius ** 3
    calculate_surface_gravity(planet)
    stream = draw.Stream(f"settle-giant-air:{float(planet.mass)!r}:{float(planet.distance)!r}:{float(planet.star.mass)!r}:"
                         f"{float(planet.star.age)!r}:{float(planet.radius)!r}")
    with draw.bound(stream):
        calculate_atmospheric_conditions(planet)
        apply_class(planet)
        generate_spin(planet, planet.primary_mass_kg, birth=getattr(planet, "_birth_spin", None))  # the birth spin again worn to the final age, then the field, the dose, the colour and the score
    if drop_unstable_moons(planet):  # a smaller giant holds a smaller moon and a smaller range of orbits
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons])
    return True


def _check_given_bounds(mass_kg, radius_km):
    """The mass and radius a caller gave must be inside `tuning.GIVEN_MASS_KG_BOUNDS` and `tuning.GIVEN_RADIUS_KM_BOUNDS`.

    Raises:
        ValueError: One is outside its bounds (or not a number).
    """
    for name, value, unit, (low, high) in (("mass", mass_kg, "kg", tuning.GIVEN_MASS_KG_BOUNDS),
                                           ("radius", radius_km, "km", tuning.GIVEN_RADIUS_KM_BOUNDS)):
        if value is not None and not low <= value <= high:
            raise ValueError(f"The given {name} ({value:.4g} {unit}) is outside {low:.3g} to {high:.3g} {unit}.")


def _check_given_density(mass_kg, radius_km, kind):
    """A mass and a radius a caller gave together must make a body of a density that kind of body can have.

    Raises:
        ValueError: The mean density is outside the bounds.
    """
    density = mass_kg / ((4 / 3) * math.pi * (radius_km * constants.KM_TO_M_FACTOR) ** 3) / 1000.0
    low, high = tuning.GIANT_DENSITY_BOUNDS_G_CM3 if kind == "giant" else bulk_draw.DENSITY_BOUNDS_G_CM3
    if not low <= density <= high:
        raise ValueError(f"A {kind} body of {mass_kg:.4g} kg and {radius_km:.4g} km has a density of {density:.4g} g/cm3; "
                         f"a {kind} body has {low:g} to {high:g}.")


def _star_fe_h(star):
    """The star's [Fe/H] for the formation law (0 for a star without one, a remnant)."""
    fe_h = getattr(star, "metallicity_fe_h", None)
    return 0.0 if fe_h is None else fe_h


def apply_giant(planet, radius_given, mass_given):
    """
    Sets a giant's mass, radius, density and composition: the mass first (`_apply_giant_mass`), then its heavy elements, radius
    and interior (`giant_draw.draw_giant`, GEN.460 and GEN.642), with the anchored mass-radius relation (GEN.34) for a body the
    radius tables do not hold.
    """
    _apply_giant_mass(planet, radius_given, mass_given)
    giant_draw.draw_giant(planet, radius_given, sample_giant_radius_km)


def _apply_giant_mass(planet, radius_given, mass_given):
    """
    Sets a giant's mass: the formation law's or the caller's, kept; with only a radius, the anchored relation's mass for it
    (`_giant_mass_for_radius_kg`).
    """
    if not mass_given:
        planet.mass = _giant_mass_for_radius_kg(planet.radius, planet.star.mass)


def update_hill_sphere(planet):
    """
    Sets `planet.hill_radius` (km) and `planet.min_orbit_distance` (5 Hill
    radii, in AU: the clearance the next body out must keep) from the
    planet's current distance and mass about `planet.primary_mass_kg` (the
    star for a planet, the parent planet for a moon). Call again whenever
    either changes, e.g. after `StarSystem.validate_system` moves the planet.
    """
    distance_m = _au_to_m(planet.distance)
    # A moon's Hill sphere is within its planet's pull, not the star's (GEN.138).
    primary_mass_kg = getattr(planet, "primary_mass_kg", None) or planet.star.mass
    planet.hill_radius = _m_to_km(calculate_hill_sphere(distance_m, planet.mass, primary_mass_kg))
    planet.min_orbit_distance = _km_to_au(5 * planet.hill_radius)


def calculate_surface_gravity(planet):
    """
    Calculates the surface gravity of the planet in g's (multiples of Earth's gravity).

    This computes the surface gravity based on the planet's mass and
    radius. The result is normalized to Earth's gravity (g's). It includes
    special adjustments for certain planet classes to ensure realistic values.

    Args:
        planet (Planet): The planet to calculate surface gravity for. Its
                         `gravity` attribute is set in place.

    Raises:
        ValueError: If the computed gravity is zero or negative.
    """
    radius_meters = planet.radius * 1000
    surface_gravity = relations.surface_gravity(planet.mass, radius_meters)
    surface_gravity_g = surface_gravity / constants.EARTH_GRAVITY
    if surface_gravity_g <= 0:
        raise ValueError('Invalid value for gravity.')
    planet.gravity = surface_gravity_g


def calculate_atmospheric_conditions(planet, distance_override=None, host_eccentricity=None):
    """
    Makes the air of the planet or moon and the surface temperature it gives (GEN.411, GEN.422; `surface_air`).

    The Bond albedo is drawn from the range of the body's regime (`climate_draws.albedo_range`); the heat the interior has carried out (`interior_heat.update_heat`) says how much
    of the body's volatiles has reached the surface (`surface_pressure`); the water it holds is the reservoir its air's water vapour
    stands in contact with; and the mixture, the pressure with the freeze-out, the density, the scale height, the greenhouse and the
    temperature are solved together (`surface_air.solve`). A body that holds no air, and a body whose air is lost, have the temperature the
    starlight alone gives.

    Args:
        planet (Planet): The planet to calculate atmospheric conditions for. Its `surface_temperature`, `atmospheric_pressure`,
                         `atmosphere`, `atm_density`, `atm_molar_density`, `scale_height`, gases, hydrosphere and aerosols are set
                         in place.
        distance_override (float, optional): An override for the planet's distance, used for special cases like moons.
        host_eccentricity (float, optional): For a moon, its parent planet's eccentricity, which the light it gets goes by (GEN.727);
                                             kept on the moon (`Planet.host_eccentricity`) and left as it was when not given.
    """
    if host_eccentricity is not None:
        planet.host_eccentricity = host_eccentricity
    distance = float(distance_override) if distance_override is not None else float(planet.distance)
    # GEN.727: the light of an eccentric orbit, averaged over it, is that of a circle of radius a (1 - e^2)^(1/4)
    distance = orbit_elements.flux_distance_au(distance, planet.host_eccentricity if planet.is_moon else planet.eccentricity)
    orbital_radius_km = _au_to_km(distance)
    output_area = 4 * math.pi * orbital_radius_km ** 2
    solar_output_at_orbit = (planet.star.luminosity / output_area) / 1e6
    # The range of the body's regime (a frozen or a wet world reflects more than dry rock): it reads only what the body is made of and
    # the starlight, so it can be drawn before the air it will be solved with.
    albedo = draw.uniform(*climate_draws.albedo_range(planet.body_type, planet.water_mass_fraction, solar_output_at_orbit,
                                                      getattr(planet.star, "temperature", None)))
    # Floored at the cosmic microwave background: starlight alone can put a
    # body around a very dim star (or far out) below ~2.7 K, but nothing in
    # space is colder than the background it sits in.
    surface_temperature_no_atmosphere = climate_draws.equilibrium_temperature_k(albedo, solar_output_at_orbit)

    # GEN.699: the heat the interior has carried out (no draws), which the degassing of the air follows; GEN.88: the surface water, drawn
    # before the air because the ocean is the reservoir the air's water vapour and its carbonate cycle stand in contact with.
    magnetism.update_exposure(planet, distance)     # GEN.86: the star's XUV and flares at this distance (no draws); the water loss reads it
    heat = interior_heat.update_heat(planet)
    giant = planet.body_type == "g"
    surface_w = None if giant else hydrosphere.draw_surface_water(planet)
    planet._surface_w = surface_w                   # GEN.839: kept for `refresh_at_final_age`, not stored
    water_lost = (not giant) and hydrosphere.water_loss_kg(planet) > 0.0

    planet.surface_temperature = surface_temperature_no_atmosphere
    planet.atmospheric_pressure = 0.0
    climate_draws.record(planet, albedo, 0.0)       # the mixture subject reads the albedo
    setting = surface_air.setting_of(planet, albedo=albedo, starlight_w_m2=solar_output_at_orbit, star_distance_au=distance,
                                     heat=heat, surface_w=surface_w, water_lost=water_lost)
    air = surface_air.solve(setting, surface_temperature_no_atmosphere)
    if air is None or not air.partials_pa:
        planet.atm_density = planet.atm_molar_density = planet.scale_height = None
        atmosphere.store_air(planet, {}, None if giant else _mantle_redox(planet, air, distance, water_lost))
    else:
        planet.surface_temperature = air.temperature_k
        planet.atmospheric_pressure = air.pressure_pa
        planet.atm_molar_density = air.state.mean_molar_mass_kg_mol
        planet.atm_density = air.state.density_kg_m3
        planet.scale_height = _m_to_km(air.state.scale_height_m)
        factor = air.greenhouse_factor
        if giant:
            # GEN.412 reads the temperature as starlight alone times (1 + G)^(1/4); a giant also radiates its own heat, which the air's factor
            # leaves out of T_e, so the stored G is taken against the starlight alone (GEN.1030) and the pair reproduces the temperature
            factor = climate_draws.greenhouse_from_temperature(air.temperature_k, albedo, solar_output_at_orbit)
            factor = air.greenhouse_factor if factor is None else max(factor, 0.0)
        climate_draws.record(planet, albedo, factor)     # GEN.412: the two numbers the temperature is made of
        atmosphere.store_air(planet, air.partials_pa, None if giant else _mantle_redox(planet, air, distance, water_lost))

    # GEN.88: water, land, ice and the ocean's chemistry.
    hydrosphere.generate_hydrosphere(planet, 0.0 if heat is None else heat.flux_w_m2, surface_w)
    # GEN.721: stagnant or mobile lid, which the magnetic field, the heat and the ground below all read (the hydrosphere above
    # read the heat of a stagnant lid, the regime not being drawn yet; it needs the hydrosphere). The heat is made again with it.
    interior_heat.generate_lid(planet)
    interior_heat.update_heat(planet)
    # GEN.441: the minerals, ices and ocean on the ground, from the bulk fractions, the temperature, the air and the water.
    surface_composition.generate_surface(planet)
    # GEN.414: cloud, haze and dust depths, on the body's own stream (no draw of the run's moves).
    aerosols.generate_aerosols(planet)


def refresh_at_final_age(body):
    """
    GEN.839: makes again, at the heat of the star's final age, what read the heat of the provisional one: the ice shell and the rest of
    the hydrosphere (`hydrosphere.generate_hydrosphere`), the ground with its volcanism gate (`surface_composition.generate_surface`),
    the ground dose when the ground changed (`radiation.update_ground_dose`), and the habitability score that reads them all. Called in the life pass, after
    `interior_heat.update_heat` has run again at the final age. The draws these make (the brine, the ocean chemistry, the ground and its
    gates) come from the body's own stream (`util.draw.Stream`, seeded from what is final about it), so no draw of the run's stream moves;
    a body with no water drawn (a giant, a belt, or one loaded from the database) is left as it is.
    """
    surface_w = getattr(body, "_surface_w", None)
    if body.body_type in ("g", "a") or surface_w is None:
        return
    seed = f"final-age:{body.is_moon}:{body.planet_class}:{body.mass!r}:{body.radius!r}:{body.distance!r}:{getattr(body.star, 'age', None)!r}"
    flux = internal_heat_flux_w_m2(body)            # the heat at the star's age now (no draws; the life pass has run it already)
    ground = surface_composition.mix_of(body)
    with draw.bound(draw.Stream(seed)):
        hydrosphere.generate_hydrosphere(body, flux, surface_w)
        surface_composition.generate_surface(body)
        if surface_composition.mix_of(body) != ground:
            radiation.update_ground_dose(body)
    habitability_world.generate(body, flux)


def _mantle_redox(planet, air, distance, water_lost):
    """The mantle redox (log10 above iron-wuestite) of a rocky body: the air mixture's when it has one, else (no air, or an envelope's
    mixture, which has none) the recipe's redox rows at the body's own mass and age (`surface_air.mantle_redox_diw`)."""
    if air is not None and air.mixture.redox_diw is not None:
        return air.mixture.redox_diw
    return surface_air.mantle_redox_diw(planet, distance, water_lost)


def internal_heat_flux_w_m2(planet):
    """A rocky body's own heat flow, W/m^2, at its star's age now: radiogenic plus primordial heat (`interior_heat.update_heat`, GEN.699),
    stored with the body (its Urey ratio is an output), and with it the flow of its tide from its orbit (GEN.737, `tidal_flux_w_m2`, not in
    this figure: the habitability energy term adds the two). 0 for a body with no stored content (a giant, or a row from before GEN.699)."""
    found = interior_heat.update_heat(planet)
    return found.own_flux_w_m2 if found is not None else 0.0


def generate_orbit_shape(planet, primary_mass_kg, planned_objects=None, parent_distance_au=None):
    """
    Draws `planet.eccentricity` (GEN.727, `physics.orbit_elements`): a planet's from the law of the system it is in (`planned_objects`
    is the number of planets and belts planned round its star, one for a single planet; `None` counts as several), cut where its
    periapsis would fall inside the star's Roche limit; a moon's from how it came to be there (`planet.moon_origin`), cut where its
    apoapsis would leave its parent's stable Hill sphere (`parent_distance_au` is the parent's distance from the star). The tides of the
    orbit then take what they can in the age of the star, from a planet's and a captured moon's eccentricity, which nothing renews; a
    regular moon's is forced (by its neighbours and the resonances among them, as Io's is by Europa and Ganymede) and stays, which is what
    makes its tide a lasting heat source (GEN.737). Draws the eccentricity before the climate, which goes by the mean light of
    the orbit; the plane and the angles come with the rest of the motion (`generate_orbital_motion_properties`).

    Args:
        planet (Planet): The planet or moon; needs `distance`, `mass`, `radius` and `star`.
        primary_mass_kg (float): The mass (kg) of what it orbits.
        planned_objects (int, optional): As above.
        parent_distance_au (float, optional): For a moon, the parent's distance from the star, AU.
    """
    if planet.is_moon:
        hill_au = None
        if parent_distance_au is not None and primary_mass_kg and planet.star.mass:
            hill_au = _km_to_au(_m_to_km(calculate_hill_sphere(_au_to_m(parent_distance_au), primary_mass_kg, planet.star.mass)))
        e_max = orbit_elements.max_moon_eccentricity(planet.distance, hill_au)
        birth = orbit_elements.draw_moon_eccentricity(planet.moon_origin, e_max)
    else:
        e_max = orbit_elements.max_planet_eccentricity(planet.distance, planet.star, planet)
        birth = orbit_elements.draw_planet_eccentricity(None if planned_objects is None else planned_objects > 1, e_max)
    if planet.is_moon and planet.moon_origin != "captured":
        planet.eccentricity = birth
        return
    time_yr = orbit_elements.circularisation_time_yr(planet, primary_mass_kg, planet.distance)
    planet.eccentricity = orbit_elements.circularised(birth, _gyr_to_yr(planet.star.age), time_yr)
    planet._orbit_birth = birth  # not stored: `settle_orbit` takes the tides again once the star's age is final (GEN.985)


def generate_orbital_motion_properties(planet, primary_mass_kg):
    """
    Sets a planet's (or moon's) 3D orbital orientation, position, orbital
    speed, and rotation period.

    Together with `planet.distance` (its semi-major axis) and
    `planet.eccentricity` (`generate_orbit_shape`), `orbital_inclination_deg`
    (the tilt of the orbital plane) and `orbital_ascending_node_deg` (where
    that plane crosses its primary's reference plane) fully orient the
    orbit in 3D, with a periapsis drawn here (uniform round the plane); the
    mean anomaly drawn here is where the body starts around it (uniform
    in time, so a body is found nearer apoapsis most often), and neither
    is kept: the body is its state vector (`position_x/y/z`, `velocity_*_kms`,
    see `update_orbital_position`), which `planetgen.cli.orbits` advances
    over time on `planet.period` with `kepler.universal_step` (GEN.289).

    `position_x/y/z` (AU, via `update_orbital_position`) are this body's
    Cartesian position *relative to its orbital anchor* -- the star (or,
    for a binary system, the `BinaryStarProxy` standing in for the system's
    combined center) for a planet, the parent planet for a moon -- the same
    "each body positioned relative to its immediate primary, not some
    absolute frame" convention `docs/design/galaxy/galaxy-coordinate-system.md`
    already uses one level up for sectors/systems relative to the galactic
    center.

    `rotation_period_hours` is a separate, purely descriptive "day length"
    stat (this generator doesn't track rotational phase -- nothing consumes
    "which side currently faces the primary"). `generate_spin` draws the
    body's birth spin and wears it down by the tide of what it orbits for the
    age of its star (`body_spin`, GEN.731 and GEN.957): the moon ends up locked
    (`rotation_period_hours` == its own orbital period, converted to hours)
    or in a spin-orbit resonance when the tide has had time to bring it there,
    which is the norm for large or close-in moons, not a rare special case,
    but genuinely not universal: a large moon far from a low-mass primary can
    easily have a locking time longer than the system itself, which is
    exactly why not every moon comes out tidally locked here either.

    Args:
        planet (Planet): The planet or moon to set these properties on.
                         Must already have `distance`, `period`, `mass`,
                         `radius`, and `body_type` set (see
                         `Planet.__init__`).
        primary_mass_kg (float): The mass (kg) of the body this one
                                 actually orbits -- the star for an
                                 ordinary planet, the parent planet for a
                                 moon (see `Planet.__init__`'s parameter of
                                 the same name, which this is threaded
                                 straight through from).
    """
    if planet.is_moon and getattr(planet, "moon_origin", None) is not None:
        # GEN.714: the orbit follows from how the moon came to be there (a captured one can be retrograde)
        planet.orbital_inclination_deg = moon_composition.draw_inclination_deg(planet.moon_origin)
    elif planet.is_moon:
        planet.orbital_inclination_deg = draw.uniform(0, constants.MOON_ORBITAL_INCLINATION_MAX_DEG)
    else:
        planet.orbital_inclination_deg = orbit_elements.draw_planet_inclination_deg()  # GEN.727: Rayleigh, as the Kepler multiples are
    planet.orbital_ascending_node_deg = draw.uniform(0, 360)
    arg_periapsis_deg = orbit_elements.draw_argument_of_periapsis_deg()
    mean_anomaly_rad = draw.uniform(0, _TWO_PI_RAD)

    generate_spin(planet, primary_mass_kg)
    true_anomaly_rad, _ = kepler.true_anomaly_and_distance_elliptical(mean_anomaly_rad, planet.eccentricity, planet.distance)
    update_orbital_position(planet, arg_periapsis_deg + math.degrees(true_anomaly_rad), arg_periapsis_deg)


def draw_birth_spin(planet):
    """The spin a planet or moon is born with, `(rate rad/s, tilt degrees, capture units)`: a rocky planet's from the giant impacts that built it
    (`spin.accretion_spin`), a moon's or giant's at a rate from `spin.planet_day_hours` with a Rayleigh tilt (GEN.104, GEN.979). The capture units are one uniform
    for each resonance the spin can pass without being caught (`spin_torques.CAPTURE_RATIOS`, GEN.1172), drawn here so that running the history again decides the same."""
    if planet.body_type == "t" and not planet.is_moon:
        rate, tilt = spin.accretion_spin(planet.mass, planet.radius)
    else:
        rate, tilt = 2.0 * math.pi / _h_to_s(spin.planet_day_hours(planet.body_type)), spin.rayleigh_tilt_deg()
        if planet.body_type == "g" and not planet.is_moon:
            tilt = giant_tilt.tilt_after_merger_deg(planet.mass, planet.radius, tilt)    # GEN.1165: an ice giant's merger can lay it on its side
    return rate, tilt, tuple(draw.random() for _ in spin_torques.CAPTURE_RATIOS)


def generate_spin(planet, primary_mass_kg, birth=None, found=None):
    """
    Draws a planet's or moon's `rotation_period_hours` and spin axis (GEN.104, GEN.731; `physics.spin.SPIN_FIELDS`). It is born spinning
    (`draw_birth_spin`, kept on the body as `_birth_spin` so that a later pass can run its history again without drawing it) and the tide
    of what it orbits then wears the spin toward the orbit's rate for the age of its star (`body_spin.evolve`: the gravitational tide with
    the body's own frequency-dependent rheology and, for a rocky planet with an air, the thermal tide of the air, GEN.957; the orbit's
    eccentricity from `generate_orbit_shape`). A planet with moons runs the history of the moons with its spin
    (`moon_history.apply`, GEN.959), which calls this again with the `found` state it reached. A body
    the tide has brought all the way down turns once an orbit (locked) or in a spin-orbit resonance such as Mercury's 3:2
    (`spin.spin_state`), with its axis in the first Cassini state (a small share of the orbit's inclination off the orbit normal, `tidal_spin.cassini_tilt_deg`); one the thermal and gravitational tides hold in balance away from synchronism turns at that
    rate, upright or upside down; any other keeps its tilt, a rocky planet's from the impacts and a moon's or giant's a Rayleigh one. Called again when a move changes the body's period, since that can change the state.

    Args:
        planet (Planet): The planet or moon.
        primary_mass_kg (float): What it orbits.
        birth (tuple, optional): `(rate, tilt, capture units)` to start from instead of drawing it.
        found (spin_torques.SpinResult, optional): The state to set instead of running the history (`moon_history.apply`'s).
    """
    if birth is None:
        birth = draw_birth_spin(planet)
    planet._birth_spin = birth
    spin_rate, birth_tilt, capture_units = birth
    if found is None:
        if not planet.is_moon and (getattr(planet, "moons", None) or getattr(planet, "_born_moons", None)):
            return moon_history.apply(planet, primary_mass_kg)
        found = body_spin.evolve(planet, primary_mass_kg, _au_to_m(planet.distance), spin_rate, _yr_to_s(_gyr_to_yr(planet.star.age)),
                                 getattr(planet, "eccentricity", 0.0) or 0.0, capture_units=capture_units)
    period_hours = planet.period * _yr_to_h(1.0)
    if found.ratio is not None:
        planet.rotation_period_hours = period_hours / found.ratio
        tilt = tidal_spin.cassini_tilt_deg(planet.orbital_inclination_deg, draw.random())    # GEN.1165: the first Cassini state, not exactly 0
    else:
        planet.rotation_period_hours = 2.0 * math.pi / abs(found.rate) / 3600.0
        if found.kind == "asynchronous":
            tilt = 0.0 if found.rate > 0.0 else 180.0      # the tides have brought it to its balance: upright, or upside down (Venus)
        else:
            tilt = birth_tilt
            if found.rate < 0.0:
                tilt = 180.0 - tilt
    spin.set_spin(planet, spin.orbit_normal(planet.orbital_inclination_deg,
                                            planet.orbital_ascending_node_deg), tilt)
    # GEN.767: the day and night shares follow the rotation and tilt just set (a planet's; a moon's need its parent's axis, `refinish_rocky`)
    if not planet.is_moon:
        daylight.generate(planet)
    # GEN.86: the field depends on the rotation just set.
    magnetism.generate_field(planet, None if planet.body_type == "g" else interior_heat.core_of(planet))
    # GEN.87: the surface dose, from the air, water and field above.
    radiation.generate_dose(planet)
    # DB.96: the atmosphere colour, from the air, the ozone and the aerosols above.
    body_colour.generate(planet)
    # GEN.89: the habitability score, from all of the above.
    habitability_world.generate(planet, internal_heat_flux_w_m2(planet))



def move_orbit(body, distance_au, primary_mass_kg):
    """
    Sets the semi-major axis of a planet's or moon's orbit and what follows from it (period, Hill sphere, position and speed), keeping the
    plane, the shape and where the body is along the orbit (no draws): the move of a moon by the tide of its planet (`moon_history`, GEN.959).
    """
    orbit = body.orbit_from_vector()
    arg_periapsis_deg = math.degrees(orbit["arg_periapsis"])
    phase_deg = body_state.argument_of_latitude_deg((body.position_x, body.position_y, body.position_z),
                                                    body.orbital_inclination_deg, body.orbital_ascending_node_deg)
    body.distance = distance_au
    body.period = calculate_orbital_period_years(distance_au, primary_mass_kg)
    update_hill_sphere(body)
    update_orbital_position(body, phase_deg, arg_periapsis_deg)


def update_orbital_position(planet, phase_deg=None, arg_periapsis_deg=0.0):
    """
    Recomputes `planet.position_x/y/z` (AU, relative to its orbital
    anchor), its velocity, `planet.orbital_speed_kms`, and
    `planet.min_update_interval_years` from its current `distance` (the
    semi-major axis), `eccentricity`, `period`, and orbital plane
    (`orbital_inclination_deg`/`orbital_ascending_node_deg`), as the
    ellipse `state_vectors.state_from_elements` makes of them (a circle
    when the eccentricity is 0).

    Called both at initial generation (`generate_orbital_motion_properties`,
    which passes the angles it rolled) and again by
    `StarSystem.validate_system` whenever it corrects a top-level planet's
    `distance` post-hoc to resolve an orbital overlap -- the same
    "recompute anything that depends on distance" treatment `period`
    already gets there (see that method's docstring). There the body keeps
    where it is along its plane and the orbit keeps its periapsis: no
    angle is passed and the ones the current state has are used.
    Position/speed/`min_update_interval_years` would
    otherwise silently go stale relative to the corrected `distance`/
    `period` the same way `period` itself used to before that fix.

    `planetgen.cli.orbits`'s periodic time-based advancement is a separate
    path (`planetgen.db.store.advance_orbital_phases`) that moves the
    stored state vector in the database rather than through this function
    -- this generator has no live "simulation loop" over in-memory
    objects, only one-shot generation (this function) followed by periodic
    database updates (see that module's own docstring).
    `min_update_interval_years` depends only on `period`, fixed once
    generated, so the update never touches it -- it only *reads* the
    column, to decide whether a given call's `elapsed_years` is even worth
    writing for this row.

    Args:
        planet (Planet): The planet or moon to update in place. Must
                         already have `distance`, `eccentricity`, `period`,
                         `orbital_inclination_deg` and
                         `orbital_ascending_node_deg` set, and, without
                         `phase_deg`, its position and velocity.
        phase_deg (float, optional): Where along its orbit plane the body
                         is (the argument of latitude, from the ascending
                         node, which is the argument of periapsis plus the
                         true anomaly); the one its current position has
                         when not given.
        arg_periapsis_deg (float, optional): The angle from the ascending
                         node to periapsis (meaningless for a circle);
                         the one its current state has when `phase_deg` is
                         not given.
    """
    if phase_deg is None:
        orbit = planet.orbit_from_vector()
        arg_periapsis_deg = math.degrees(orbit["arg_periapsis"])
        phase_deg = body_state.argument_of_latitude_deg(
            (planet.position_x, planet.position_y, planet.position_z),
            planet.orbital_inclination_deg, planet.orbital_ascending_node_deg)
    position, velocity = state_from_elements(
        (2.0 * math.pi / planet.period) ** 2 * planet.distance ** 3, orbit_elements.periapsis_au(planet.distance, planet.eccentricity), planet.eccentricity,
        math.radians(planet.orbital_inclination_deg), math.radians(planet.orbital_ascending_node_deg), math.radians(arg_periapsis_deg),
        math.radians(phase_deg - arg_periapsis_deg))
    planet.set_position_au(*position)
    to_kms = _au_yr_to_kms(1.0)
    planet.set_velocity_kms(*(v * to_kms for v in velocity))
    planet.orbital_speed_kms = math.sqrt(sum(v * v for v in velocity)) * to_kms
    planet.min_update_interval_years = minimum_update_interval_years(planet.period)


def reconcile_zone(planet, primary_mass_kg, distance_override=None):
    """
    Re-derives `planet`'s zone from its *current* orbital distance and records the mass of the body it orbits.

    Meant to be called after something external to normal generation -- `StarSystem.validate_system`'s orbital-overlap
    correction -- has moved `planet.distance` (for an ordinary planet) or its parent's distance (for a moon, via
    `distance_override`). The mass, radius and composition of the body do not change with where it sits; its climate does, so the
    caller makes the air again (`calculate_atmospheric_conditions`) and the class follows (`apply_class`).

    Args:
        planet (Planet): The planet or moon to reconcile in place. Must already have `distance`, `zone` and `habitable_zone` set.
        primary_mass_kg (float): The mass (kg) of the body this one actually orbits -- the star for an ordinary planet, the parent
            planet for a moon.
        distance_override (float, optional): The distance (AU) to determine the zone from, when it differs from `planet.distance`
            -- a moon's own `distance` is its orbit around its *parent planet*, so a moon passes its parent's `distance` here.

    Returns:
        bool: True if the zone changed.
    """
    distance = distance_override if distance_override is not None else planet.distance
    new_zone = zone_for(planet.habitable_zone, distance)
    zone_changed = new_zone != planet.zone
    planet.zone = new_zone
    planet.primary_mass_kg = primary_mass_kg  # what its Hill sphere is measured against (GEN.138)
    return zone_changed


def zone_for(habitable_zone, distance_au):
    """The zone letter ('h' hot, 'e' ecosphere, 'c' cold) of `distance_au` against `habitable_zone` (inner, outer) in AU."""
    inner_bound, outer_bound = habitable_zone
    if distance_au < inner_bound:
        return 'h'
    if distance_au > outer_bound:
        return 'c'
    return 'e'


def hosting_population(planet):
    """
    What the hosting gate (`physics.hosting_gate`) calls `planet`: a small body (under `tuning.HOSTING_SMALL_BODY_MAX_RADIUS_KM` in
    radius) is a `trans-neptunian-cold` or `trans-neptunian-hot` body past `tuning.HOSTING_TRANS_NEPTUNIAN_MIN_SNOW_LINES`
    snow lines of its star (cold under `tuning.HOSTING_COLD_MAX_INCLINATION_DEG` of inclination, Stephens and Noll 2006), a
    `belt-small-body` in the cold zone short of that, and an `inner-small-body` in the hot or habitable zone; every other
    planet is a `planet`, which the gate's share table has no rows for (its host-class table decides).
    """
    if planet.radius >= tuning.HOSTING_SMALL_BODY_MAX_RADIUS_KM:
        return "planet"
    if planet.zone != "c":
        return "inner-small-body"
    if planet.distance < tuning.HOSTING_TRANS_NEPTUNIAN_MIN_SNOW_LINES * planet.star.snow_line_au:
        return "belt-small-body"
    inclination = planet.orbital_inclination_deg or 0.0
    return "trans-neptunian-cold" if inclination < tuning.HOSTING_COLD_MAX_INCLINATION_DEG else "trans-neptunian-hot"


def hosting_host(planet):
    """`planet` as the hosting gate reads it (`hosting_gate.Host`): its population, size, mass and orbit about what it
    orbits (the planets of this generator are on circular orbits); a gas or ice giant is fluid."""
    return hosting_gate.Host(
        hosting_population(planet), 2.0 * _km_to_m(planet.radius), planet.mass, parent_mass_kg=planet.primary_mass_kg,
        orbit_radius_m=_au_to_m(planet.distance), fluid=planet.body_type != "t")


def satellite_share(planet):
    """The chance that `planet` has satellites from the gate's `satellite-share` table, or `None` for a planet the table has
    no rows for (the host-class chance of the draw decides)."""
    return hosting_gate.satellite_share(hosting_host(planet))


def moon_orbit_bounds_km(planet):
    """
    The range of orbits, in km from `planet`'s center, where a moon of it
    can sit.

    Innermost: clear of the planet's body and the largest moon it could
    have (`planet.radius / 10**(1/3)`, so no moon touches it; this is also
    past the rigid-body Roche limit of ~1.26 planet radii for similar
    densities), plus 15 atmospheric scale heights (or 100 km with no
    atmosphere) of drag-free margin. Outermost: the prograde stability
    limit, `MOON_PROGRADE_STABLE_HILL_FRACTION` of the Hill radius; past it
    the star strips the moon.

    Returns:
        tuple: `(low_km, high_km)`; `low_km >= high_km` means no room.
    """
    max_moon_radius = planet.radius / (10 ** (1 / 3))
    atmosphere_margin_km = planet.scale_height * 15 if planet.scale_height else 100
    low_km = planet.radius + max_moon_radius + atmosphere_margin_km
    high_km = planet.hill_radius * tuning.MOON_PROGRADE_STABLE_HILL_FRACTION
    return low_km, high_km


def drop_unstable_moons(planet):
    """
    Removes the moons of `planet` that no longer fit after its class (and
    so its radius and mass) was regenerated: one orbiting outside
    `moon_orbit_bounds_km`, or one larger than the planet could hold
    (`moon_size_limits`: `planet.radius / 10**(1/3)` and a tenth of its
    mass). A real planet that lost mass this way
    would lose those moons to the star or to a collision.

    Returns:
        int: How many moons were dropped.
    """
    low_km, high_km = moon_orbit_bounds_km(planet)
    max_moon_radius, max_moon_mass = moon_size_limits(planet)
    kept = [moon for moon in planet.moons
            if low_km <= _au_to_km(moon.distance) and _au_to_km(orbit_elements.apoapsis_au(moon.distance, moon.eccentricity)) <= high_km
            and moon.radius <= max_moon_radius and moon.mass <= max_moon_mass]
    dropped = len(planet.moons) - len(kept)
    planet.moons[:] = kept
    return dropped


def moon_size_limits(parent):
    """The largest moon `parent` can hold: `(radius_km, mass_kg)`, at most
    `radius / 10**(1/3)` and a tenth of its mass."""
    return parent.radius / (10 ** (1 / 3)), parent.mass / 10


MOON_REJECTION_LIMIT = 5
"""int: How many moons in a row `generate_moons` may discard for being larger than their parent can hold before it stops."""


def generate_moons(planet, moon_count=None):
    """
    Generates a system of moons for the given planet.

    This procedurally generates moons for the planet. It determines the
    number, size, and orbital distance of the moons based on the planet's
    properties, such as its mass and Hill radius. The generated moons are
    themselves `Planet` instances, with the `is_moon` flag set, appended to
    `planet.moons`.

    Args:
        planet (Planet): The planet to generate moons for.
        moon_count (int, optional): An exact number of moons to attempt to
            generate. If None (the default), moons are generated until no
            more orbital room is available, as before. If given, generation
            stops once `moon_count` moons exist, even if more room remains;
            if the planet doesn't have room for `moon_count` moons, fewer
            than requested may be generated.
    """
    if moon_count == 0:
        return
    # GEN.761: a body asks the hosting gate first (a satellite must fit, a stable orbit must exist, and the body's
    # population and size must have a chance of satellites at all); the share of the gated table is the draw's, in
    # `Planet._draw_moons`
    if not hosting_gate.can_host(hosting_host(planet)).allowed:
        return
    # The lightest moon must fit under the largest the parent can hold (`moon_size_limits`).
    limit_radius_km, limit_mass_kg = moon_size_limits(planet)
    if limit_mass_kg <= tuning.MOON_MIN_MASS_KG:
        return
    parent_is_giant = planet.body_type == "g"

    low_orbit, high_orbit = moon_orbit_bounds_km(planet)
    total_orbit_distance = low_orbit
    parent_kind = moon_composition.parent_kind(planet.envelope_mass_fraction or 0.0, planet.mass)
    # GEN.466: a giant's regular moons share one circumplanetary disk (the lifetime, the ice line, the ice ramp)
    disk_system = moon_composition.draw_disk_system(planet.mass) if parent_kind == moon_composition.GIANT else None
    if disk_system is not None:
        moon_composition.store_disk(planet, disk_system)

    # Deferred import: planet imports this module at load time, so Planet
    # can't be imported here at module level without a circular import.
    from planetgen.generation.planet import Planet

    rejected = 0
    while total_orbit_distance < high_orbit and total_orbit_distance < _au_to_km(planet.distance):
        if moon_count is not None and len(planet.moons) >= moon_count:
            break

        moon_mass = planet_kind.draw_moon_mass_kg(planet.mass, parent_is_giant, limit_mass_kg)
        if moon_mass is None:
            break

        # Log-uniform, not linear-uniform: [total_orbit_distance, high_orbit]
        # can span many orders of magnitude (high_orbit reaches out to about
        # half the planet's own Hill radius, which for a large planet is tens
        # of millions of km -- far beyond where any real large moon
        # actually orbits, e.g. our Moon at ~384,400 km), and a plain
        # random.uniform over that range spends almost all its density in the
        # single largest order of magnitude, so nearly every moon landed
        # implausibly far out. Real moon systems are much closer to
        # log-spaced -- e.g. the Galilean moons run 421,700 / 671,100 /
        # 1,070,400 / 1,882,700 km, each roughly 1.5-1.6x the last, not a
        # near-flat distribution across the whole possible range. Sampling
        # log-uniformly gives every order of magnitude equal weight instead,
        # which both matches that real spacing pattern better and -- since
        # tidal-locking timescale scales with distance^6
        # (`body_spin.evolve`) -- stops real tidal-locking
        # physics from calling almost every moon unlocked purely because the
        # old distribution pushed it implausibly far from its primary.
        moon_distance_km = math.exp(draw.uniform(math.log(total_orbit_distance), math.log(high_orbit)))
        moon_distance = _km_to_au(moon_distance_km)
        # GEN.714: the origin comes first and the orbit follows from it
        origin = moon_composition.draw_origin(parent_kind).row.name
        new_moon = Planet(planet.system_config, planet.star, planet.habitable_zone, moon_distance,
                          mass=moon_mass, zone_override=planet.zone,
                          distance_override=planet.distance, is_moon=True, primary_mass_kg=planet.mass,
                          moon_origin=origin, moon_disk_system=disk_system, host_eccentricity=planet.eccentricity)
        if new_moon.radius > limit_radius_km:    # denser or lighter than the parent could hold; its slot stays open
            draw.redrawn("moon.too_large")
            rejected += 1
            if rejected >= MOON_REJECTION_LIMIT:
                break
            continue
        rejected = 0
        planet.moons.append(new_moon)
        total_orbit_distance = _au_to_km(new_moon.distance) + _au_to_km(new_moon.min_orbit_distance)

    if planet.moons:
        # This planet's own reflex-offset "wobble" (schema v20) from the
        # combined pull of its own moons -- a proper two-body treatment
        # alongside each moon's own unchanged position_x/y/z (relative to
        # this planet); see Planet.reflex_offset_x's own docstring and
        # orbits.calculate_reflex_offset. Left at its 0.0 default when this
        # planet ends up with no moons (every early-return path above).
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons]
        )
