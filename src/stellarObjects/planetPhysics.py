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
import secrets

from . import physical_constants, program_constants
from .utils import (calculate_object_mass, calculate_hill_sphere, circular_orbital_speed_kms,
                    orbital_position_au, position_change_interval_hours, reseed_rng,
                    sample_bounded_bell)


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

        if planet_type == "t":  # Terrestrial planet
            min_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density
            max_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density
        else:  # Gas giant
            min_atm_density, max_atm_density = physical_constants.ATMOSPHERE_DENSITY[planet_type]
            min_core_ratio, max_core_ratio = program_constants.GAS_GIANT_CORE_ATMOSPHERE_RATIO

            min_core_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density * min_core_ratio
            max_core_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density * max_core_ratio

            min_atm_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_atm_density * (1 - min_core_ratio)
            max_atm_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_atm_density * (1 - max_core_ratio)

            min_mass = min_core_mass + min_atm_mass
            max_mass = max_core_mass + max_atm_mass

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
    return random.choices(eligible, weights=weights, k=1)[0]


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
    """
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
    reseed_rng()
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

    if planet.planet_class is None and planet.radius is None and planet.mass is None:
        # Fully random generation
        valid_classes = [c for c, data in program_constants.PLANET_CLASSES.items() if data[zone]]
        if planet.system_config.HABITABLE_WORLD is False and zone == 'e':
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
        if planet.system_config.HABITABLE_WORLD is False and zone == 'e':
            possible_classes = [c for c in possible_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]
        if not possible_classes:
            raise ValueError("No valid planet class for the given radius in this zone")
        planet.planet_class = secrets.choice(possible_classes)
        _validate_radius(planet)

    elif planet.planet_class is None and planet.radius is None and planet.mass is not None:
        # Mass given, determine possible classes
        possible_classes = [c for c, data in program_constants.PLANET_CLASSES.items()
                            if planet_mass_ranges[c][0] <= planet.mass <= planet_mass_ranges[c][1] and data[zone]]
        if planet.system_config.HABITABLE_WORLD is False and zone == 'e':
            possible_classes = [c for c in possible_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]
        if not possible_classes:
            raise ValueError("No valid planet class for the given mass in this zone")
        planet.planet_class = secrets.choice(possible_classes)
        _validate_mass(planet)

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
        if planet.system_config.HABITABLE_WORLD is False and zone == 'e':
            possible_classes = [c for c in possible_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]
        if not possible_classes:
            raise ValueError("No valid planet class for the given radius/mass in this zone")
        planet.planet_class = secrets.choice(possible_classes)
        _validate_radius(planet)
        _validate_mass(planet)

    else:
        # All inputs provided, fully validate
        _validate_planet_class(planet, zone)
        _validate_no_habitable_world(planet, zone)
        _validate_radius(planet)
        _validate_mass(planet)

    class_data = program_constants.PLANET_CLASSES[planet.planet_class]
    planet.composition = class_data["composition"]
    planet.description = class_data["description"]
    if planet.is_moon:
        # Class descriptions are authored generically (e.g. "shields the inner
        # planets"); a moon knows it's a moon, so correct the wording once here
        # rather than scrubbing every rendered paragraph for it later.
        planet.description = re.sub(r'\bplanets\b', 'moons', planet.description)
        planet.description = re.sub(r'\bplanet\b', 'moon', planet.description)
    planet.body_type = class_data["type"]

    # Class-specific density range if declared (e.g. a brown-dwarf-like
    # sub-stellar class -- see get_planet_mass_ranges above and
    # program_constants.PLANET_CLASSES), else the default range shared by
    # every other class of this body type.
    default_density = physical_constants.PLANET_DENSITY[planet.body_type]
    min_density, max_density = class_data.get("density_range", default_density)
    planet.density = random.uniform(min_density, max_density)

    if class_data["atmosphere"] is None:
        planet.atmosphere = "None"
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

    if planet.body_type == 'g' and "density_range" not in class_data:
        # Core/envelope density blend -- for ORDINARY gas giants only (no
        # class-specific density_range declared). A class that declares its
        # own density_range (e.g. a brown-dwarf-like sub-stellar class, see
        # program_constants.PLANET_CLASSES) skips this blend entirely and
        # keeps the density already drawn above as final: real brown dwarfs
        # don't have a meaningfully separate light envelope over a denser
        # core the way an ordinary gas giant does -- electron degeneracy
        # pressure keeps the whole body close to uniformly dense throughout,
        # so blending toward a light "puffy" envelope value here would just
        # dilute the real, elevated density_range right back down (the
        # blend is a harmonic mean, dominated by whichever term is smaller
        # almost regardless of mass fraction).
        core_to_atmosphere_ratio = random.uniform(*program_constants.GAS_GIANT_CORE_ATMOSPHERE_RATIO)
        # core_to_atmosphere_ratio is a MASS fraction (core mass / total
        # mass), not a volume/density-averaging weight. The physically
        # correct way to combine the core and atmosphere densities via a
        # mass fraction is the mass-weighted harmonic mean (1/density_total
        # = mass_fraction/density_a + (1-mass_fraction)/density_b) -- an
        # arithmetic mean of the two densities (the previous formula) has no
        # physical basis and let atmosphere-heavy blends drag the whole
        # planet's density down far below either component's own range.
        # The envelope side of the blend uses physical_constants.
        # GAS_ENVELOPE_BULK_DENSITY (a real, g/cm^3-scale "puffy gas giant"
        # bulk density), NOT planet.atm_density -- that value is ~1000x
        # lighter (it's the thin surface/pressure-layer density used by
        # calculate_atmospheric_conditions, a different physical layer), and
        # a harmonic mean is dominated by whichever term is smaller almost
        # regardless of mass fraction, so using it here collapsed every gas
        # giant's density to a near-zero, physically meaningless value no
        # matter how dense its core (see GAS_ENVELOPE_BULK_DENSITY's
        # docstring for the real-world numbers this replaced it with).
        envelope_density_gcm3 = random.uniform(*physical_constants.GAS_ENVELOPE_BULK_DENSITY)
        planet.density = 1 / (core_to_atmosphere_ratio / planet.density + (1 - core_to_atmosphere_ratio) / envelope_density_gcm3)

    planet.volume, planet.mass = calculate_object_mass(planet.planet_class, planet.radius, program_constants.PLANET_CLASSES, physical_constants.PLANET_DENSITY,
                                              planet.density)

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
    reseed_rng()
    radius_meters = planet.radius * 1000
    surface_gravity = (physical_constants.G * planet.mass) / (radius_meters ** 2)
    surface_gravity_g = surface_gravity / physical_constants.EARTH_GRAVITY
    if surface_gravity_g <= 0:
        raise ValueError('Invalid value for gravity.')
    planet.gravity = surface_gravity_g


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
    reseed_rng()
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
    surface_temperature_no_atmosphere = (
                                                (1 - albedo) * solar_output_at_orbit / (4 * physical_constants.STEFAN_BOLTZMANN_CONSTANT)) ** (
                                                    1 / 4)

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
        planet.surface_temperature = surface_temperature_atmosphere
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
    reseed_rng()
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
    `planet.position_change_interval_hours` from its current `distance`,
    `period`, `radius`, and orbital elements (`orbital_inclination_deg`/
    `orbital_ascending_node_deg`/`orbital_phase_deg`), via
    `utils.orbital_position_au`/`circular_orbital_speed_kms`/
    `position_change_interval_hours`.

    Called both at initial generation (`generate_orbital_motion_properties`,
    right after `orbital_phase_deg` is rolled) and again by
    `StarSystem.validate_system` whenever it corrects a top-level planet's
    `distance` post-hoc to resolve an orbital overlap -- the same
    "recompute anything that depends on distance" treatment `period`
    already gets there (see that method's docstring). Position/speed/
    notice-interval would otherwise silently go stale relative to the
    corrected `distance` the same way `period` used to before that fix.

    `updateOrbits.py`'s periodic time-based advancement is a separate,
    SQL-only path (`stellarObjects._db.advance_orbital_phases`) that
    recomputes position directly in the database rather than through this
    function -- this generator has no live "simulation loop" over
    in-memory objects, only one-shot generation (this function) followed
    by periodic database updates (see that module's own docstring).
    `position_change_interval_hours` doesn't depend on phase at all
    (only `distance`/`radius`/`period`, all fixed once generated), so
    `advance_orbital_phases` never needs to touch it.

    Args:
        planet (Planet): The planet or moon to update in place. Must
                         already have `distance`, `period`, `radius`,
                         `orbital_inclination_deg`,
                         `orbital_ascending_node_deg`, and
                         `orbital_phase_deg` set.
    """
    planet.position_x, planet.position_y, planet.position_z = orbital_position_au(
        planet.distance, planet.orbital_inclination_deg,
        planet.orbital_ascending_node_deg, planet.orbital_phase_deg,
    )
    planet.orbital_speed_kms = circular_orbital_speed_kms(planet.distance, planet.period)
    planet.position_change_interval_hours = position_change_interval_hours(
        planet.radius, planet.orbital_speed_kms
    )


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
    reseed_rng()
    if moon_count == 0:
        return
    max_moon_mass = planet.mass / 10
    max_moon_radius = planet.radius / (10 ** (1 / 3))
    possible_classes = [c for c, data in planet_mass_ranges.items()
                        if program_constants.PLANET_CLASSES[c][planet.zone] and program_constants.PLANET_CLASSES[c]["type"] == 't' and c not in program_constants.MOON_BLACKLIST
                        and data[1] <= max_moon_mass and program_constants.PLANET_CLASSES[c]['radius_range'][1] <= max_moon_radius]
    # Mirrors generate_planet_properties' own "Fully random generation"
    # habitable-world filtering: a moon is subject to the same
    # HABITABLE_WORLD=False rule as any other body, but wasn't checked here
    # -- unreachable in practice until gas giants could be placed in zone
    # 'e' (see PLANET_CLASSES["J"]'s zone rework), since no non-gas-giant
    # planet generates moons of its own zone's habitable classes any
    # differently. A gas giant now placed in 'e' with HABITABLE_WORLD=False
    # must not roll a habitable-class moon.
    if planet.system_config.HABITABLE_WORLD is False and planet.zone == 'e':
        possible_classes = [c for c in possible_classes if c not in program_constants.HABITABLE_PLANET_CLASSES]
    if not possible_classes:
        return

    low_orbit = planet.scale_height * 15 if planet.scale_height else 100
    high_orbit = planet.min_orbit_distance * physical_constants.AU_TO_KM
    total_orbit_distance = low_orbit

    # Deferred import: planetData imports this module at load time, so Planet
    # can't be imported here at module level without a circular import.
    from .planetData import Planet

    while total_orbit_distance < high_orbit and total_orbit_distance < (planet.distance * physical_constants.AU_TO_KM):
        if moon_count is not None and len(planet.moons) >= moon_count:
            break

        moon_class = _choose_weighted_planet_class(possible_classes)

        radius_limit = program_constants.PLANET_CLASSES[moon_class]['radius_range'][1] if max_moon_radius > \
                                                                        program_constants.PLANET_CLASSES[moon_class]['radius_range'][
                                                                            1] else max_moon_radius
        # Log-uniform, not linear-uniform: [total_orbit_distance, high_orbit]
        # can span many orders of magnitude (high_orbit reaches out to 1/5 of
        # the planet's own Hill radius, which for a large planet is tens to
        # hundreds of millions of km -- far beyond where any real large moon
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
        # radius_limit may be narrower than moon_class's own declared
        # radius_range ceiling (capped by the parent's Hill sphere/mass
        # above) -- _sample_class_radius's size_mode reading is evaluated
        # against this actually-available [min, radius_limit] window, not
        # necessarily the class's full range.
        moon_radius = _sample_class_radius(moon_class, program_constants.PLANET_CLASSES[moon_class]['radius_range'][0], radius_limit)

        new_moon = Planet(planet.system_config, planet.star, planet.habitable_zone, moon_distance,
                          radius=moon_radius, planet_class=moon_class, zone_override=planet.zone,
                          distance_override=planet.distance, is_moon=True, primary_mass_kg=planet.mass)
        planet.moons.append(new_moon)
        total_orbit_distance = (new_moon.distance * physical_constants.AU_TO_KM) + (new_moon.min_orbit_distance * physical_constants.AU_TO_KM)
