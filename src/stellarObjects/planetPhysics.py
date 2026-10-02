# stellarObjects/planetPhysics.py

"""
Planet Physical and Orbital Generation
=======================================

This module determines what a planet or moon physically *is*: its orbital
zone, class, radius, mass, density, composition, atmosphere, surface
gravity, and atmospheric conditions, plus the generation of any moons it
has. It operates on `Planet` instances (from `planetData`) passed in as
arguments, rather than as methods on the class, so that generation logic is
kept separate from the `Planet` object's own definition and presentation.

Life chemistry and evolutionary data are intentionally NOT determined here —
see `planetLife.apply_life_data`, which `StarSystem` applies to every planet
and moon in a system after all of them have been generated.
"""

import math
import random
import re

from . import log, physical_constants, program_constants
from .utils import (calculate_object_mass, calculate_hill_sphere, calculate_reflex_offset,
                    circular_orbital_speed_kms, minimum_update_interval_years,
                    finite_domain, orbital_position_au, sample_bounded_bell,
                    sample_power_law)


def _sample_class_radius(cls, min_radius, max_radius):
    """
    Draws a radius in [min_radius, max_radius] using `cls`'s declared
    `size_mode` (see `program_constants.PLANET_CLASSES`), via
    `utils.sample_bounded_bell` -- a bell-curve draw peaking at the
    class's statistically-most-common size, rather than a flat uniform
    draw across its whole declared range. `min_radius`/`max_radius` are
    passed separately from `cls`'s own declared `radius_range` (rather
    than read directly here) since a moon's actually-available range can
    be narrower (capped by its parent's Hill sphere/mass) -- `size_mode`'s
    "how far through the available range" reading is evaluated against
    whatever range is actually being drawn from for this call.
    """
    size_mode = program_constants.PLANET_CLASSES[cls].get("size_mode", 0.5)
    return sample_bounded_bell(min_radius, max_radius, size_mode)


def uses_giant_mass_radius(cls):
    """Whether class `cls` gets its mass and radius from the giant-planet
    mass-radius relation (GEN.34): every gas-giant class without its own
    `"density_range"`."""
    data = program_constants.PLANET_CLASSES[cls]
    return data["type"] == "g" and "density_range" not in data


def giant_mass_range_kg(cls):
    """Class `cls`'s giant mass range (`GIANT_CLASS_MASS_RANGE_EARTH`), in kg."""
    low, high = program_constants.GIANT_CLASS_MASS_RANGE_EARTH.get(
        cls, program_constants.GIANT_DEFAULT_MASS_RANGE_EARTH)
    return low * physical_constants.EARTH_MASS_TO_KG, high * physical_constants.EARTH_MASS_TO_KG


def _giant_transition_earth():
    """`(mass, radius)` in Earth units where the Neptunian and Jovian
    branches of the giant mass-radius relation meet."""
    m_n, r_n, s_n = physical_constants.GIANT_NEPTUNIAN_MASS_RADIUS
    m_j, r_j, s_j = physical_constants.GIANT_JOVIAN_MASS_RADIUS
    # r_n * (M/m_n)^s_n == r_j * (M/m_j)^s_j, solved in log space.
    log_m = (math.log(r_j / r_n) + s_n * math.log(m_n) - s_j * math.log(m_j)) / (s_n - s_j)
    mass = math.exp(log_m)
    return mass, r_n * (mass / m_n) ** s_n


GIANT_TRANSITION_MASS_EARTH, GIANT_TRANSITION_RADIUS_EARTH = _giant_transition_earth()


def giant_regime(mass_kg):
    """`"neptunian"` below the relation's transition mass, else `"jovian"`."""
    mass_earth = mass_kg / physical_constants.EARTH_MASS_TO_KG
    return "neptunian" if mass_earth < GIANT_TRANSITION_MASS_EARTH else "jovian"


def giant_radius_km(mass_kg):
    """The giant mass-radius relation's median radius (km) for `mass_kg`
    (`physical_constants.GIANT_NEPTUNIAN_MASS_RADIUS`), with no scatter."""
    mass_earth = mass_kg / physical_constants.EARTH_MASS_TO_KG
    if giant_regime(mass_kg) == "neptunian":
        m0, r0, slope = physical_constants.GIANT_NEPTUNIAN_MASS_RADIUS
    else:
        m0, r0, slope = physical_constants.GIANT_JOVIAN_MASS_RADIUS
    return r0 * (mass_earth / m0) ** slope * physical_constants.EARTH_RADIUS_KM


def _sample_giant_mass_kg(low_kg, high_kg):
    """A giant's mass in [low_kg, high_kg], from dN/dlogM ~
    M^-GIANT_MASS_FUNCTION_SLOPE."""
    return sample_power_law(low_kg, high_kg, program_constants.GIANT_MASS_FUNCTION_SLOPE)


def _sample_giant_radius_km(cls, mass_kg):
    """`giant_radius_km(mass_kg)` with the relation's own scatter
    (`GIANT_RADIUS_SCATTER`, cut at 3 sigma), kept inside class `cls`'s
    radius range."""
    sigma = physical_constants.GIANT_RADIUS_SCATTER[giant_regime(mass_kg)]
    factor = 1 + max(-3.0, min(3.0, random.gauss(0.0, 1.0))) * sigma
    low, high = program_constants.PLANET_CLASSES[cls]["radius_range"]
    return min(high, max(low, giant_radius_km(mass_kg) * factor))


def _giant_mass_for_radius_kg(cls, radius_km):
    """
    A mass (kg) for a giant of class `cls` whose radius was fixed first (a
    caller's radius): the Neptunian branch's inverse below the transition
    radius, else a draw from the class's masses past the transition (the
    Jovian branch is nearly flat, so radius barely constrains mass there).
    Always inside the class's mass range.
    """
    low_kg, high_kg = giant_mass_range_kg(cls)
    transition_kg = GIANT_TRANSITION_MASS_EARTH * physical_constants.EARTH_MASS_TO_KG
    radius_earth = radius_km / physical_constants.EARTH_RADIUS_KM
    if radius_earth < GIANT_TRANSITION_RADIUS_EARTH or high_kg <= transition_kg:
        m0, r0, slope = physical_constants.GIANT_NEPTUNIAN_MASS_RADIUS
        mass_kg = m0 * (radius_earth / r0) ** (1 / slope) * physical_constants.EARTH_MASS_TO_KG
    else:
        mass_kg = _sample_giant_mass_kg(max(low_kg, transition_kg), high_kg)
    return min(high_kg, max(low_kg, mass_kg))


def get_planet_mass_ranges():
    """
    Calculates and returns a dictionary of valid mass ranges (in kg) for each planet class.

    This function pre-calculates the minimum and maximum possible mass for each
    planet class based on its radius and density ranges. This is used for
    validation when generating planets with a given mass. The calculations
    account for the different compositions of terrestrial planets and gas giants.

    Returns:
        dict: A dictionary where keys are planet classes (e.g., 'M', 'N') and
              values are tuples of (min_mass, max_mass) in kilograms.
    """
    mass_ranges = {}
    for planet_class, data in program_constants.PLANET_CLASSES.items():
        # radius_range is in kilometers (as used everywhere else -- e.g. class
        # M's 5000-10000 matches Earth's ~6371 km radius); convert to meters
        # to match the kg/m^3 densities used below.
        min_radius, max_radius = data["radius_range"]
        min_radius *= physical_constants.KM_TO_M_FACTOR
        max_radius *= physical_constants.KM_TO_M_FACTOR
        planet_type = data["type"]

        # Class-specific density range if declared (e.g. a brown-dwarf-like
        # sub-stellar class, far denser than an ordinary gas giant -- see
        # program_constants.PLANET_CLASSES), else the default range shared
        # by every other class of this body type. Same per-class-override
        # pattern as atm_molar_density_range.
        min_density, max_density = data.get("density_range", physical_constants.PLANET_DENSITY[planet_type])  # g/cm^3

        # Convert density from g/cm³ to kg/m³ for mass calculation
        min_density *= 1000
        max_density *= 1000

        if uses_giant_mass_radius(planet_class):
            min_mass, max_mass = giant_mass_range_kg(planet_class)
        else:  # Terrestrial, or a gas giant with its own density_range
            min_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density
            max_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density

        mass_ranges[planet_class] = (min_mass, max_mass)
    return mass_ranges


planet_mass_ranges = get_planet_mass_ranges()


def _choose_weighted_planet_class(valid_classes):
    """
    Draws a single planet class from `program_constants.PLANET_CLASS_PROBABILITIES`,
    restricted to `valid_classes`.

    This re-weights the distribution to only the eligible classes and draws
    once, rather than repeatedly sampling the full distribution and
    rejecting draws that fall outside `valid_classes`.

    Args:
        valid_classes (iterable): The planet class codes eligible to be chosen.

    Returns:
        str: The chosen planet class code.
    """
    valid_classes = set(valid_classes)
    eligible = [c for c in program_constants.PLANET_CLASS_PROBABILITIES if c in valid_classes]
    weights = [program_constants.PLANET_CLASS_PROBABILITIES[c] for c in eligible]
    chosen = random.choices(eligible, weights=weights, k=1)[0]
    log.choice("Planet class", chosen,
               f"weighted draw among {len(eligible)} eligible classes {eligible} "
               f"(weights {weights}) out of {len(program_constants.PLANET_CLASS_PROBABILITIES)} total")
    return chosen


def _habitable_classes_barred(planet, zone):
    """
    Whether a randomly chosen class for `planet` must skip the habitable
    classes: the system disallows habitable worlds and this is the
    ecosphere, or its star is younger than `LIFE_MIN_STAR_AGE_GY` (too
    young for a crust and oceans, let alone life). An explicitly requested
    class is only held to the first rule (`_validate_no_habitable_world`).
    """
    if planet.system_config.HABITABLE_WORLD is False and zone == 'e':
        return True
    star_age = getattr(getattr(planet, "star", None), "age", None)
    return star_age is not None and star_age < program_constants.LIFE_MIN_STAR_AGE_GY


def _validate_no_habitable_world(planet, zone):
    """
    Raises if the system disallows habitable worlds and this planet's
    class is a habitable one being placed in the ecosphere.

    Args:
        planet (Planet): The planet being generated.
        zone (str): The zone ('h', 'c', 'e') the planet would occupy.

    Raises:
        ValueError: If `planet.system_config.HABITABLE_WORLD` is False,
                   `zone` is the ecosphere ('e'), and `planet.planet_class`
                   is a habitable class.
    """
    if planet.system_config.HABITABLE_WORLD is False and zone == 'e' and planet.planet_class in program_constants.HABITABLE_PLANET_CLASSES:
        raise ValueError(f"Cannot generate habitable planet class {planet.planet_class} in ecosphere when HABITABLE_WORLD is False.")


def _validate_planet_class(planet, zone):
    """
    Validates if the planet's class is valid for its zone.

    Args:
        planet (Planet): The planet being generated.
        zone (str): The zone ('h', 'c', 'e') of the planet.

    Raises:
        ValueError: If the planet class is not valid for the given zone.
    """
    if planet.planet_class not in program_constants.PLANET_CLASSES or not program_constants.PLANET_CLASSES[planet.planet_class][zone]:
        raise ValueError("Invalid planet class for this zone")


def _validate_radius(planet):
    """
    Validates if the planet's radius is within the allowed range for its class.

    Args:
        planet (Planet): The planet being generated.

    Raises:
        ValueError: If the radius is outside the valid range for the planet's class.
    """
    min_radius, max_radius = program_constants.PLANET_CLASSES[planet.planet_class]["radius_range"]
    if not (min_radius <= planet.radius <= max_radius):
        raise ValueError("Invalid radius for planet class")


def _validate_mass(planet):
    """
    Validates if the planet's mass is within the allowed range for its class.

    Args:
        planet (Planet): The planet being generated.

    Raises:
        ValueError: If the mass is outside the valid range for the planet's class.
    """
    min_mass, max_mass = planet_mass_ranges[planet.planet_class]
    if not (min_mass <= planet.mass <= max_mass):
        raise ValueError("Invalid mass for planet class")


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
    primary_mass_sol = primary_mass_kg / physical_constants.SOLAR_MASS_TO_KG
    return math.sqrt(distance_au ** 3 / primary_mass_sol)


def generate_planet_properties(planet, zone_override=None):
    """
    Generates a planet's physical and orbital properties: zone, class,
    radius, mass, density, composition, and atmosphere.

    This is the core of the planet generation logic. It determines the
    planet's properties such as class, composition, and atmosphere. The
    generation can be fully random or guided by inputs already set on
    `planet` (a specific radius, mass, or planet class). It ensures the
    generated properties are consistent with each other and the planet's
    orbital zone.

    Args:
        planet (Planet): The planet to generate properties for. Any of
                         `planet.planet_class`/`planet.radius`/`planet.mass`
                         already set constrain the generation; unset ones
                         are filled in.
        zone_override (str, optional): A character ('h', 'c', 'e') to manually
                                       set the planet's zone, overriding the
                                       calculation based on distance.
    """
    # Determine the planet's zone (hot, ecosphere, or cold)
    inner_bound, outer_bound = planet.habitable_zone
    if planet.distance < inner_bound:
        zone = 'h'
    elif planet.distance > outer_bound:
        zone = 'c'
    else:
        zone = 'e'

    if zone_override and zone_override.lower() in "hce":
        zone = zone_override.lower()
    planet.zone = zone

    # --- Input Validation and Random Generation ---
    # This section handles the logic for generating planet properties based on
    # the inputs provided. It can generate a fully random planet, or generate
    # properties based on a given class, radius, or mass.

    radius_given = planet.radius is not None
    mass_given = planet.mass is not None

    if planet.planet_class is None and planet.radius is None and planet.mass is None:
        # Fully random generation
        valid_classes = [c for c, data in program_constants.PLANET_CLASSES.items() if data[zone]]
        if _habitable_classes_barred(planet, zone):
            valid_classes = [c for c in valid_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]

        planet.planet_class = _choose_weighted_planet_class(valid_classes)
        min_radius, max_radius = program_constants.PLANET_CLASSES[planet.planet_class]["radius_range"]
        planet.radius = _sample_class_radius(planet.planet_class, min_radius, max_radius)

    elif planet.planet_class is not None and planet.radius is None and planet.mass is None:
        # Class given, generate radius
        _validate_planet_class(planet, zone)
        _validate_no_habitable_world(planet, zone)
        min_radius, max_radius = program_constants.PLANET_CLASSES[planet.planet_class]["radius_range"]
        planet.radius = _sample_class_radius(planet.planet_class, min_radius, max_radius)

    elif planet.planet_class is None and planet.radius is not None and planet.mass is None:
        # Radius given, determine possible classes
        possible_classes = [c for c, data in program_constants.PLANET_CLASSES.items()
                            if data[zone] and data["radius_range"][0] <= planet.radius <= data["radius_range"][1]]
        if _habitable_classes_barred(planet, zone):
            possible_classes = [c for c in possible_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]
        if not possible_classes:
            raise ValueError("No valid planet class for the given radius in this zone")
        planet.planet_class = random.choice(possible_classes)
        _validate_radius(planet)

    elif planet.planet_class is None and planet.radius is None and planet.mass is not None:
        # Mass given, determine possible classes
        possible_classes = [c for c, data in program_constants.PLANET_CLASSES.items()
                            if planet_mass_ranges[c][0] <= planet.mass <= planet_mass_ranges[c][1] and data[zone]]
        if _habitable_classes_barred(planet, zone):
            possible_classes = [c for c in possible_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]
        if not possible_classes:
            raise ValueError("No valid planet class for the given mass in this zone")
        planet.planet_class = random.choice(possible_classes)
        _validate_mass(planet)
        # Everything downstream needs a radius, so draw one for the class,
        # the same as the class+mass branch below.
        min_radius, max_radius = program_constants.PLANET_CLASSES[planet.planet_class]["radius_range"]
        planet.radius = _sample_class_radius(planet.planet_class, min_radius, max_radius)

    elif planet.planet_class is not None and planet.radius is not None and planet.mass is None:
        # Class and radius given, validate
        _validate_planet_class(planet, zone)
        _validate_no_habitable_world(planet, zone)
        _validate_radius(planet)

    elif planet.planet_class is not None and planet.radius is None and planet.mass is not None:
        # Class and mass given, validate and generate radius
        _validate_planet_class(planet, zone)
        _validate_no_habitable_world(planet, zone)
        _validate_mass(planet)
        min_radius, max_radius = program_constants.PLANET_CLASSES[planet.planet_class]["radius_range"]
        planet.radius = _sample_class_radius(planet.planet_class, min_radius, max_radius)

    elif planet.planet_class is None and planet.radius is not None and planet.mass is not None:
        # Radius and mass given, determine possible classes
        possible_classes = []
        for c, data in program_constants.PLANET_CLASSES.items():
            min_mass, max_mass = planet_mass_ranges[c]
            min_radius, max_radius = data["radius_range"]
            if min_mass <= planet.mass <= max_mass and min_radius <= planet.radius <= max_radius and data[zone]:
                possible_classes.append(c)
        if _habitable_classes_barred(planet, zone):
            possible_classes = [c for c in possible_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]
        if not possible_classes:
            raise ValueError("No valid planet class for the given radius/mass in this zone")
        planet.planet_class = random.choice(possible_classes)
        _validate_radius(planet)
        _validate_mass(planet)

    else:
        # All inputs provided, fully validate
        _validate_planet_class(planet, zone)
        _validate_no_habitable_world(planet, zone)
        _validate_radius(planet)
        _validate_mass(planet)

    class_data = program_constants.PLANET_CLASSES[planet.planet_class]

    # Ecosphere-zone classes with a declared "zone_position_mode" (Venus/
    # Earth/Mars-analog-style classes -- see program_constants.PLANET_CLASSES)
    # get placed at a class-appropriate distance within the zone instead of
    # wherever the caller's initial estimate happened to land -- the zone
    # itself doesn't change (still 'e'), only the position within it, so
    # this runs after `zone`/`planet.zone` are already settled above. Moons
    # are excluded: `planet.distance` for a moon is its orbit around the
    # *parent planet*, not an AU-scale position within the star's habitable
    # zone `planet.habitable_zone` describes, so redrawing it here would
    # corrupt it, not correct it (see generate_moons/zone_override).
    if zone == 'e' and not planet.is_moon and "zone_position_mode" in class_data:
        inner_bound, outer_bound = planet.habitable_zone
        planet.distance = sample_bounded_bell(inner_bound, outer_bound, class_data["zone_position_mode"])

    planet.composition = class_data["composition"]
    planet.description = class_data["description"]
    if planet.is_moon:
        # Class descriptions are authored generically (e.g. "shields the inner
        # planets"); a moon knows it's a moon, so correct the wording once here
        # rather than scrubbing every rendered paragraph for it later.
        planet.description = re.sub(r'\bplanets\b', 'moons', planet.description)
        planet.description = re.sub(r'\bplanet\b', 'moon', planet.description)
    planet.body_type = class_data["type"]

    if uses_giant_mass_radius(planet.planet_class):
        _apply_giant_mass_and_radius(planet, radius_given, mass_given)
    else:
        # Class-specific density range if declared (e.g. a brown-dwarf-like
        # sub-stellar class -- see get_planet_mass_ranges above and
        # program_constants.PLANET_CLASSES), else the default range shared by
        # every other class of this body type.
        default_density = physical_constants.PLANET_DENSITY[planet.body_type]
        min_density, max_density = class_data.get("density_range", default_density)
        planet.density = random.uniform(min_density, max_density)

    if class_data["atmosphere"] is None:
        planet.atmosphere = "None"
        # Clear any atmosphere a previous class left behind --
        # `reconcile_zone_and_class` regenerates a body in place, so a body
        # moved from an atmosphered class into an airless one would
        # otherwise keep (and save) its old class's densities.
        planet.atm_density = None
        planet.atm_molar_density = None
        planet.scale_height = None
    else:
        planet.atmosphere = class_data["atmosphere"]
        # Class-specific atmosphere-density/molar-density ranges if declared
        # (e.g. Class N's Venus-like dense, CO2-heavy atmosphere), else the
        # default range shared by every other class of this body type. This
        # replaces the old Class N hardcoded special case with the same
        # general per-class-override mechanism Class P's albedo_range uses,
        # so every class's atmosphere composition can be tuned independently
        # (see program_constants.PLANET_CLASSES and
        # docs/analysis/habitability-atmosphere-sanity-review.md).
        default_a_density = physical_constants.ATMOSPHERE_DENSITY[planet.body_type]
        min_a_density, max_a_density = class_data.get("atm_density_range", default_a_density)
        planet.atm_density = random.uniform(min_a_density, max_a_density)
        default_am_density = physical_constants.ATMOSPHERIC_MOLAR_DENSITY[planet.body_type]
        min_am_density, max_am_density = class_data.get("atm_molar_density_range", default_am_density)
        planet.atm_molar_density = random.uniform(min_am_density, max_am_density)

    planet.volume, planet.mass = calculate_object_mass(planet.planet_class, planet.radius, program_constants.PLANET_CLASSES, physical_constants.PLANET_DENSITY,
                                              planet.density)

    update_hill_sphere(planet)


def _apply_giant_mass_and_radius(planet, radius_given, mass_given):
    """
    Sets a gas giant's mass, radius and density from the giant mass-radius
    relation (GEN.34) instead of a radius and a density drawn apart: with
    neither given, a mass from the class's mass function (at most
    `GIANT_MAX_STAR_MASS_RATIO` of its star's mass) and a radius from the
    relation; with only a mass, the relation's radius for it; with only
    a radius, the relation's mass for it (`_giant_mass_for_radius_kg`);
    with both, those two. The density is whatever the mass and radius make,
    so `calculate_object_mass` gives the same mass back.
    """
    if not mass_given:
        if radius_given:
            planet.mass = _giant_mass_for_radius_kg(planet.planet_class, planet.radius)
        else:
            low_kg, high_kg = giant_mass_range_kg(planet.planet_class)
            star_mass = getattr(getattr(planet, "star", None), "mass", None)
            if star_mass:
                high_kg = max(low_kg, min(high_kg, program_constants.GIANT_MAX_STAR_MASS_RATIO * star_mass))
            planet.mass = _sample_giant_mass_kg(low_kg, high_kg)
    if not radius_given:
        planet.radius = _sample_giant_radius_km(planet.planet_class, planet.mass)
    volume_m3 = (4 / 3) * math.pi * (planet.radius * physical_constants.KM_TO_M_FACTOR) ** 3
    planet.density = planet.mass / volume_m3 / 1000  # g/cm^3


def update_hill_sphere(planet):
    """
    Sets `planet.hill_radius` (km) and `planet.min_orbit_distance` (5 Hill
    radii, in AU: the clearance the next body out must keep) from the
    planet's current distance and mass. Call again whenever either changes,
    e.g. after `StarSystem.validate_system` moves the planet.
    """
    distance_m = planet.distance * physical_constants.AU_TO_M
    planet.hill_radius = calculate_hill_sphere(distance_m, planet.mass, planet.star.mass) / 1000  # Convert to km
    planet.min_orbit_distance = (5 * planet.hill_radius) / physical_constants.AU_TO_KM


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
    surface_gravity = (physical_constants.G * planet.mass) / (radius_meters ** 2)
    surface_gravity_g = surface_gravity / physical_constants.EARTH_GRAVITY
    if surface_gravity_g <= 0:
        raise ValueError('Invalid value for gravity.')
    planet.gravity = surface_gravity_g


@finite_domain()
def _atmosphere_retention_factor(gravity_g):
    """
    A gravity-based atmosphere-retention scaling factor, normalized to 1.0
    at Earth gravity (`gravity_g == 1.0`).

    Substituting `scale_height_m` into the old `atmospheric_pressure =
    atm_density * surface_gravity_ms2 * scale_height_m` line shows gravity
    cancels out exactly (`atmospheric_pressure = atm_density * R * T /
    atm_molar_density`), making pressure unphysically independent of a
    planet's gravity -- higher gravity should let a planet retain more
    atmosphere and support higher surface pressure. This factor is applied
    to an *effective* atmospheric density used only in the pressure
    calculation (see `calculate_atmospheric_conditions`), not to
    `planet.atm_density` itself, since that value also feeds the gas-giant
    density blend in `generate_planet_properties` -- modifying it here would
    create a circular dependency, since gravity is itself derived partly
    from density for gas giants.

    k=1 (linear scaling) is a starting point, not a value derived from a
    real atmospheric-retention model; tune here if generated pressure
    distributions warrant a different curve.

    Args:
        gravity_g (float): Surface gravity in Earth g's.

    Returns:
        float: The retention factor (1.0 at gravity_g == 1.0).
    """
    return gravity_g ** 1  # k=1 linear scaling as the default/starting point


def calculate_atmospheric_conditions(planet, distance_override=None):
    """
    Calculates the atmospheric conditions of the planet, including surface
    temperature and atmospheric pressure.

    This models the planet's atmospheric conditions. It calculates the
    surface temperature considering the star's output, the planet's distance,
    and the greenhouse effect of its atmosphere. It also estimates the
    atmospheric pressure based on the atmospheric mass and planet's gravity.

    Args:
        planet (Planet): The planet to calculate atmospheric conditions for.
                         Its `surface_temperature`, `atmospheric_pressure`,
                         and (if it has an atmosphere) `scale_height`
                         attributes are set in place.
        distance_override (float, optional): An override for the planet's
                                             distance, used for special cases
                                             like moons.
    """
    distance = float(distance_override) if distance_override is not None else float(planet.distance)
    orbital_radius_km = distance * physical_constants.AU_TO_KM
    output_area = 4 * math.pi * orbital_radius_km ** 2
    solar_output_at_orbit = (planet.star.luminosity / output_area) / 1e6
    class_data = program_constants.PLANET_CLASSES.get(planet.planet_class, {})
    # Class-specific albedo range if declared (e.g. Class P's icy/glaciated
    # surface reflects more than the default rocky/Earth-like range), else
    # the default range used for every other class.
    albedo_range = class_data.get("albedo_range", (0.12, 0.35))
    albedo = random.uniform(*albedo_range)
    # Floored at the cosmic microwave background: starlight alone can put a
    # body around a very dim star (or far out) below ~2.7 K, but nothing in
    # space is colder than the background it sits in.
    surface_temperature_no_atmosphere = max(
        physical_constants.COSMIC_BACKGROUND_TEMPERATURE_K,
        ((1 - albedo) * solar_output_at_orbit / (4 * physical_constants.STEFAN_BOLTZMANN_CONSTANT)) ** (1 / 4),
    )

    if planet.atmosphere == "None":
        planet.surface_temperature = surface_temperature_no_atmosphere
        planet.atmospheric_pressure = 0.0
    else:
        scale_height_m = (physical_constants.R * surface_temperature_no_atmosphere) / (
                    planet.atm_molar_density * planet.gravity * physical_constants.EARTH_GRAVITY)
        # planet.radius (and everything derived from it below) is in km, so
        # convert the scale height -- dimensionally meters, per the R*T/(M*g)
        # formula -- to km to match before it's combined with radius.
        scale_height = scale_height_m / physical_constants.KM_TO_M_FACTOR
        planet.scale_height = scale_height
        # Closed-form barometric formula for an isothermal, hydrostatic
        # atmosphere: integrating dP/dz = -rho*g from the surface to
        # infinity gives P_surface = rho_surface * g * H, where H is the
        # scale height in meters. This replaces a previous shell-integration
        # loop that summed km^3 shell volumes against a kg/m^3 density
        # without converting units, undercounting atmospheric_mass by
        # roughly 9 orders of magnitude and requiring an arbitrary "* 7500"
        # fudge factor to partially compensate.
        # Applied to an effective atmospheric density used only for this
        # pressure calculation (not to planet.atm_density itself -- see
        # _atmosphere_retention_factor's docstring for why) so higher-gravity
        # planets retain more atmosphere and support higher surface pressure,
        # reintroducing the gravity dependence that otherwise cancels out of
        # this formula algebraically.
        effective_atm_density = planet.atm_density * _atmosphere_retention_factor(planet.gravity)
        surface_gravity_ms2 = planet.gravity * physical_constants.EARTH_GRAVITY
        atmospheric_pressure = effective_atm_density * surface_gravity_ms2 * scale_height_m

        # atm_molar_density is the only atmosphere-composition signal in the
        # data model today (no explicit CO2-fraction field exists), so
        # base_ratio measures how CO2-like this draw's composition is
        # (>=1.0 at/above CO2's own molar mass). On its own that ratio can't
        # calibrate both a thin-but-CO2-heavy atmosphere (e.g. real Mars,
        # ~43.3 g/mol, negligible greenhouse effect) and a dense CO2-heavy
        # one (real Venus, ~43.45 g/mol -- almost the same molar mass, ~100x
        # the greenhouse forcing): composition alone doesn't capture
        # quantity/potency. greenhouse_multiplier_range is the per-class
        # knob for that second axis (see program_constants.PLANET_CLASSES),
        # independent of how heavy/light the class's own atmosphere is.
        # CO2_MAX_GREENHOUSE_FACTOR is now a generous safety ceiling rather
        # than the value most classes hit -- real Venus's own ratio is
        # ~101, so it only guards against a badly-configured future class.
        base_ratio = planet.atm_molar_density / physical_constants.CO2_BASE_MOLAR_DENSITY
        greenhouse_multiplier_range = class_data.get("greenhouse_multiplier_range", (1.0, 1.0))
        greenhouse_multiplier = random.uniform(*greenhouse_multiplier_range)
        greenhouse_factor = min(program_constants.CO2_MAX_GREENHOUSE_FACTOR, base_ratio * greenhouse_multiplier)
        surface_temperature_atmosphere = ((1 - albedo) * solar_output_at_orbit * (1 + greenhouse_factor) / (4 * physical_constants.STEFAN_BOLTZMANN_CONSTANT)) ** (1 / 4)
        planet.surface_temperature = max(physical_constants.COSMIC_BACKGROUND_TEMPERATURE_K, surface_temperature_atmosphere)
        planet.atmospheric_pressure = atmospheric_pressure

        # Class M/P's forced pressure/temperature clamps are disabled.
        # Several underlying model bugs they may have been compensating for
        # (the inverted greenhouse factor above, and Class M/P sharing
        # identical atmosphere sampling with no differentiation between
        # them) have since been fixed, but that isn't a guarantee every
        # generated value now lands in a narrow realistic band -- it's a
        # statistical question the physical-plausibility tooling
        # (stellarObjects/plausibility.py) is better suited to monitor than
        # a hard clamp. Commented out rather than deleted in case it needs
        # restoring.
        # if planet.planet_class == "M":
        #     if planet.atmospheric_pressure < 90000 or planet.atmospheric_pressure > 112000:
        #         planet.atmospheric_pressure = random.uniform(90000, 112000)
        #     if planet.surface_temperature < 283 or planet.surface_temperature > 290:
        #         planet.surface_temperature = random.uniform(283, 290)
        # elif planet.planet_class == "P" and planet.surface_temperature >= 283:
        #     # If surface_temperature_no_atmosphere is already above 283, we need a different approach
        #     # to ensure the P class planet remains cold.
        #     if surface_temperature_no_atmosphere < 283:
        #         planet.surface_temperature = random.uniform(surface_temperature_no_atmosphere, 283)
        #     else:
        #         planet.surface_temperature = random.uniform(200, 283) # A reasonable cold range for P class


def _tidal_locking_timescale_seconds(moon, primary_mass_kg, initial_rotation_period_hours):
    """
    Estimates how long tidal forces would take to lock `moon`'s rotation
    to its orbital period, in seconds, via the standard simplified
    tidal-despinning formula (Murray & Dermott, "Solar System Dynamics"):

        t_lock = omega0 * a^6 * I * Q / (3 * G * M_primary^2 * k2 * R_moon^5)

    with the uniform-sphere moment of inertia I = (2/5) * m_moon * R^2
    substituted in, which simplifies to:

        t_lock = (2*Q / (15*k2)) * (omega0 * a^6 * m_moon) / (G * M_primary^2 * R_moon^3)

    `Q`/`k2` are fixed representative values for a rocky/icy body
    (`physical_constants.MOON_TIDAL_DISSIPATION_Q`/`_LOVE_NUMBER_K2` --
    see that constant's own comment for real-example verification), not
    modeled per body.

    Args:
        moon (Planet): The moon (`is_moon=True`) -- must already have
                       `distance` (AU), `mass` (kg), and `radius` (km) set.
        primary_mass_kg (float): The parent planet's mass, in kg -- the
                                 actual body raising the tide, not the
                                 grandparent star `moon.star` points at.
        initial_rotation_period_hours (float): The moon's assumed
                                               pre-locking rotation period,
                                               in hours -- despinning takes
                                               longer starting from a
                                               faster spin (more angular
                                               momentum to shed).

    Returns:
        float: Estimated locking timescale, in seconds.
    """
    omega0 = 2 * math.pi / (initial_rotation_period_hours * 3600)
    distance_m = moon.distance * physical_constants.AU_TO_M
    radius_m = moon.radius * 1000
    q_over_k2 = physical_constants.MOON_TIDAL_DISSIPATION_Q / physical_constants.MOON_TIDAL_LOVE_NUMBER_K2
    numerator = omega0 * (distance_m ** 6) * moon.mass
    denominator = physical_constants.G * (primary_mass_kg ** 2) * (radius_m ** 3)
    return (2 * q_over_k2 / 15) * (numerator / denominator)


def generate_orbital_motion_properties(planet, primary_mass_kg):
    """
    Sets a planet's (or moon's) 3D orbital orientation, position, orbital
    speed, and rotation period.

    Together with `planet.distance` (this generator only ever models
    circular orbits -- no eccentricity), `orbital_inclination_deg` (the
    tilt of the orbital plane) and `orbital_ascending_node_deg` (where that
    plane crosses its primary's reference plane) fully orient the orbit in
    3D; `orbital_phase_deg` is where the body currently sits around it.
    Only `orbital_phase_deg` (and the `position_x/y/z`/`orbital_speed_kms`
    derived from it, see `update_orbital_position`) ever change after
    generation -- `updateOrbits.py` advances phase (and position in
    lockstep) over time based on `planet.period` -- the plane itself is
    fixed for the body's lifetime, the same way its `distance` is.

    `position_x/y/z` (AU, via `update_orbital_position`) are this body's
    Cartesian position *relative to its orbital anchor* -- the star (or,
    for a binary system, the `BinaryStarProxy` standing in for the system's
    combined center) for a planet, the parent planet for a moon -- the same
    "each body positioned relative to its immediate primary, not some
    absolute frame" convention `docs/design/galaxy-coordinate-system.md`
    already uses one level up for sectors/systems relative to the galactic
    center.

    `rotation_period_hours` is a separate, purely descriptive "day length"
    stat (this generator doesn't track rotational phase -- nothing consumes
    "which side currently faces the primary"). A candidate (pre-locking)
    rotation period is always drawn first from a `body_type`-appropriate
    range; for a moon, that candidate is kept only if real tidal physics
    (`_tidal_locking_timescale_seconds`) says there *hasn't* been enough
    time (`planet.star.age`, the best available proxy for the system's
    age -- planets/moons don't carry an independent age of their own) to
    lock it yet. Otherwise the moon ends up actually locked
    (`rotation_period_hours` == its own orbital period, converted to
    hours) -- which real physics makes the norm for large or close-in
    moons, not a rare special case, but genuinely not universal: a large
    moon far from a low-mass primary can easily have a locking timescale
    longer than the system itself, which is exactly why not every moon
    comes out tidally locked here either.

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
    inclination_max = (
        physical_constants.MOON_ORBITAL_INCLINATION_MAX_DEG if planet.is_moon
        else physical_constants.PLANET_ORBITAL_INCLINATION_MAX_DEG
    )
    planet.orbital_inclination_deg = random.uniform(0, inclination_max)
    planet.orbital_ascending_node_deg = random.uniform(0, 360)
    planet.orbital_phase_deg = random.uniform(0, 360)

    min_hours, max_hours = physical_constants.ROTATION_PERIOD_RANGE_HOURS[planet.body_type]
    candidate_rotation_period_hours = random.uniform(min_hours, max_hours)

    is_locked = False
    if planet.is_moon:
        lock_timescale_s = _tidal_locking_timescale_seconds(
            planet, primary_mass_kg, candidate_rotation_period_hours
        )
        system_age_s = planet.star.age * 1e9 * physical_constants.SECONDS_PER_YEAR
        is_locked = lock_timescale_s < system_age_s

    if is_locked:
        planet.rotation_period_hours = planet.period * (physical_constants.SECONDS_PER_YEAR / 3600)
    else:
        planet.rotation_period_hours = candidate_rotation_period_hours

    update_orbital_position(planet)


def update_orbital_position(planet):
    """
    Recomputes `planet.position_x/y/z` (AU, relative to its orbital
    anchor), `planet.orbital_speed_kms`, and
    `planet.min_update_interval_years` from its current `distance`,
    `period`, and orbital elements (`orbital_inclination_deg`/
    `orbital_ascending_node_deg`/`orbital_phase_deg`), via
    `utils.orbital_position_au`/`circular_orbital_speed_kms`/
    `minimum_update_interval_years`.

    Called both at initial generation (`generate_orbital_motion_properties`,
    right after `orbital_phase_deg` is rolled) and again by
    `StarSystem.validate_system` whenever it corrects a top-level planet's
    `distance` post-hoc to resolve an orbital overlap -- the same
    "recompute anything that depends on distance" treatment `period`
    already gets there (see that method's docstring). Position/speed/
    `min_update_interval_years` would otherwise silently go stale relative
    to the corrected `distance`/`period` the same way `period` itself
    used to before that fix.

    `updateOrbits.py`'s periodic time-based advancement is a separate,
    SQL-only path (`stellarObjects._db.advance_orbital_phases`) that
    recomputes position directly in the database rather than through this
    function -- this generator has no live "simulation loop" over
    in-memory objects, only one-shot generation (this function) followed
    by periodic database updates (see that module's own docstring).
    `min_update_interval_years` doesn't depend on phase at all (only
    `period`, fixed once generated), so `advance_orbital_phases` never
    needs to touch it -- it only *reads* the column, to decide whether a
    given call's `elapsed_years` is even worth writing for this row.

    Args:
        planet (Planet): The planet or moon to update in place. Must
                         already have `distance`, `period`,
                         `orbital_inclination_deg`,
                         `orbital_ascending_node_deg`, and
                         `orbital_phase_deg` set.
    """
    planet.position_x, planet.position_y, planet.position_z = orbital_position_au(
        planet.distance, planet.orbital_inclination_deg,
        planet.orbital_ascending_node_deg, planet.orbital_phase_deg,
    )
    planet.orbital_speed_kms = circular_orbital_speed_kms(planet.distance, planet.period)
    planet.min_update_interval_years = minimum_update_interval_years(planet.period)


def reconcile_zone_and_class(planet, primary_mass_kg, distance_override=None, parent=None):
    """
    Re-derives `planet`'s zone from its *current* orbital distance and, if
    the class it already has isn't valid there, regenerates its class and
    every class-derived physical property from scratch for the corrected
    zone -- rather than silently reporting a stale class (e.g. "a
    terrestrial Earth-like world") at a distance that class doesn't
    physically support.

    Meant to be called only after something external to normal generation
    -- `StarSystem.validate_system`'s orbital-overlap correction -- has
    moved `planet.distance` (for an ordinary planet) or its parent's
    distance (for a moon, via `distance_override`) out from under an
    already-chosen class. A class is a statement about the conditions at a
    body's *actual* final position, not an independently-revisable label:
    once something else changes that position, the class has to be
    re-derived from it, the same way real orbital migration changes a
    body's real climate, not just its assumed one.

    Args:
        planet (Planet): The planet or moon to reconcile in place. Must
                         already have `distance`, `planet_class`, `zone`,
                         and `habitable_zone` set (i.e. already fully
                         generated once).
        primary_mass_kg (float): The mass (kg) of the body this one
                                 actually orbits -- the star for an
                                 ordinary planet, the parent planet for a
                                 moon (same meaning as everywhere else this
                                 parameter appears).
        distance_override (float, optional): The distance (AU) to
                                              determine the zone from and
                                              feed to
                                              `calculate_atmospheric_conditions`,
                                              when it differs from
                                              `planet.distance` -- a moon's
                                              own `distance` is its orbit
                                              around its *parent planet*,
                                              not an AU-scale position in
                                              the star's own zone, so a
                                              moon always passes its
                                              parent's (corrected)
                                              `distance` here instead.
        parent (Planet, optional): A moon's planet. A regenerated moon
            only gets a class and size a moon of it may have
            (`choose_moon_class_and_radius`; GEN.25, GEN.36), never a gas
            giant, a blacklisted class or one too large for its planet.

    Returns:
        bool: True if the class (and everything derived from it) was
             actually regenerated; False if the existing class remained
             valid for the (possibly still new) zone, in which case only
             `planet.zone` itself was refreshed.
    """
    inner_bound, outer_bound = planet.habitable_zone
    distance = distance_override if distance_override is not None else planet.distance
    if distance < inner_bound:
        new_zone = 'h'
    elif distance > outer_bound:
        new_zone = 'c'
    else:
        new_zone = 'e'

    zone_changed = new_zone != planet.zone
    planet.zone = new_zone

    if not zone_changed or program_constants.PLANET_CLASSES[planet.planet_class][new_zone]:
        return False

    moon_class = moon_radius = None
    if parent is not None:
        moon_class, moon_radius = choose_moon_class_and_radius(parent, new_zone)
        if moon_class is None:
            # Unreachable in practice (class D fits any planet, anywhere);
            # leave the moon for check_lunar_system to report.
            return False

    planet.planet_class = moon_class
    planet.radius = moon_radius
    planet.mass = None
    # The caller already decided where this body sits (validate_system
    # pushed it clear of its inner neighbor), so keep that distance:
    # generate_planet_properties redraws an ecosphere class with a
    # "zone_position_mode" anywhere in the habitable zone, which could
    # drop it back inside the neighbor it was just moved past (an
    # asteroid belt's span, in a giant star's wide habitable zone).
    distance_before = planet.distance
    generate_planet_properties(planet, zone_override=new_zone)
    planet.distance = distance_before
    # The Hill sphere generate_planet_properties just set is for the
    # redrawn distance, not this one.
    update_hill_sphere(planet)
    planet.period = calculate_orbital_period_years(planet.distance, primary_mass_kg)
    calculate_surface_gravity(planet)
    calculate_atmospheric_conditions(planet, distance_override)
    generate_orbital_motion_properties(planet, primary_mass_kg)
    return True


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
    high_km = planet.hill_radius * program_constants.MOON_PROGRADE_STABLE_HILL_FRACTION
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
            if low_km <= moon.distance * physical_constants.AU_TO_KM <= high_km
            and moon.radius <= max_moon_radius and moon.mass <= max_moon_mass]
    dropped = len(planet.moons) - len(kept)
    planet.moons[:] = kept
    return dropped


def moon_size_limits(parent):
    """The largest moon `parent` can hold: `(radius_km, mass_kg)`, at most
    `radius / 10**(1/3)` and a tenth of its mass."""
    return parent.radius / (10 ** (1 / 3)), parent.mass / 10


def moon_class_options(parent, zone):
    """
    Every class a moon of `parent` may have in `zone`, mapped to the
    largest radius (km) a moon of that class may be (GEN.35, GEN.36):
    terrestrial, not in `MOON_BLACKLIST`, valid in `zone`, not a habitable
    class where `_habitable_classes_barred` says so, and able to fit under
    `moon_size_limits` at all. The ceiling is the class's own top radius,
    the largest moon radius, or the radius at which the class's densest
    rock would pass the mass limit, whichever is least. A class only has
    to fit at its smallest size, not across its whole range: requiring the
    whole range left a rocky planet only Class D moons.
    """
    max_radius_km, max_mass_kg = moon_size_limits(parent)
    barred = _habitable_classes_barred(parent, zone)
    options = {}
    for cls, data in program_constants.PLANET_CLASSES.items():
        if (not data[zone] or data["type"] != 't' or cls in program_constants.MOON_BLACKLIST
                or (barred and cls in program_constants.HABITABLE_PLANET_CLASSES)):
            continue
        max_density_kg_m3 = data.get("density_range", physical_constants.PLANET_DENSITY['t'])[1] * 1000
        mass_radius_km = (max_mass_kg / ((4 / 3) * math.pi * max_density_kg_m3)) ** (1 / 3) / physical_constants.KM_TO_M_FACTOR
        ceiling_km = min(data["radius_range"][1], max_radius_km, mass_radius_km)
        if data["radius_range"][0] <= ceiling_km:
            options[cls] = ceiling_km
    return options


def choose_moon_class_and_radius(parent, zone):
    """
    Draws a moon's class (weighted like any planet's, among
    `moon_class_options`) and a radius that fits under its ceiling.

    Returns:
        tuple: `(class, radius_km)`, or `(None, None)` when no class fits.
    """
    options = moon_class_options(parent, zone)
    if not options:
        return None, None
    moon_class = _choose_weighted_planet_class(options)
    low = program_constants.PLANET_CLASSES[moon_class]["radius_range"][0]
    return moon_class, _sample_class_radius(moon_class, low, options[moon_class])


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
    # A gas giant placed in 'e' with HABITABLE_WORLD=False must not roll a
    # habitable-class moon either: moon_class_options applies the same
    # _habitable_classes_barred rule as generate_planet_properties.
    if not moon_class_options(planet, planet.zone):
        return

    low_orbit, high_orbit = moon_orbit_bounds_km(planet)
    total_orbit_distance = low_orbit

    # Deferred import: planetData imports this module at load time, so Planet
    # can't be imported here at module level without a circular import.
    from .planetData import Planet

    while total_orbit_distance < high_orbit and total_orbit_distance < (planet.distance * physical_constants.AU_TO_KM):
        if moon_count is not None and len(planet.moons) >= moon_count:
            break

        moon_class, moon_radius = choose_moon_class_and_radius(planet, planet.zone)

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
        # (_tidal_locking_timescale_seconds) -- stops real tidal-locking
        # physics from calling almost every moon unlocked purely because the
        # old distribution pushed it implausibly far from its primary.
        moon_distance_km = math.exp(random.uniform(math.log(total_orbit_distance), math.log(high_orbit)))
        moon_distance = moon_distance_km / physical_constants.AU_TO_KM
        new_moon = Planet(planet.system_config, planet.star, planet.habitable_zone, moon_distance,
                          radius=moon_radius, planet_class=moon_class, zone_override=planet.zone,
                          distance_override=planet.distance, is_moon=True, primary_mass_kg=planet.mass)
        planet.moons.append(new_moon)
        total_orbit_distance = (new_moon.distance * physical_constants.AU_TO_KM) + (new_moon.min_orbit_distance * physical_constants.AU_TO_KM)

    if planet.moons:
        # This planet's own reflex-offset "wobble" (schema v20) from the
        # combined pull of its own moons -- a proper two-body treatment
        # alongside each moon's own unchanged position_x/y/z (relative to
        # this planet); see Planet.reflex_offset_x's own docstring and
        # utils.calculate_reflex_offset. Left at its 0.0 default when this
        # planet ends up with no moons (every early-return path above).
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons]
        )
