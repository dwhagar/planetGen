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
from .utils import calculate_object_mass, calculate_hill_sphere, reseed_rng


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

        min_density, max_density = physical_constants.PLANET_DENSITY[planet_type]  # g/cm^3

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
        planet.radius = random.uniform(min_radius, max_radius)

    elif planet.planet_class is not None and planet.radius is None and planet.mass is None:
        # Class given, generate radius
        _validate_planet_class(planet, zone)
        _validate_no_habitable_world(planet, zone)
        min_radius, max_radius = program_constants.PLANET_CLASSES[planet.planet_class]["radius_range"]
        planet.radius = random.uniform(min_radius, max_radius)

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
        planet.radius = random.uniform(min_radius, max_radius)

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

    min_density, max_density = physical_constants.PLANET_DENSITY[planet.body_type]
    planet.density = random.uniform(min_density, max_density)

    if class_data["atmosphere"] is None:
        planet.atmosphere = "None"
    else:
        planet.atmosphere = class_data["atmosphere"]
        if planet.planet_class == 'N':
            planet.atm_density = 65
            min_am_density, max_am_density = physical_constants.ATMOSPHERIC_MOLAR_DENSITY[planet.body_type]
            planet.atm_molar_density = max_am_density
        else:
            min_a_density, max_a_density = physical_constants.ATMOSPHERE_DENSITY[planet.body_type]
            planet.atm_density = random.uniform(min_a_density, max_a_density)
            min_am_density, max_am_density = physical_constants.ATMOSPHERIC_MOLAR_DENSITY[planet.body_type]
            planet.atm_molar_density = random.uniform(min_am_density, max_am_density)

    if planet.body_type == 'g':
        core_to_atmosphere_ratio = random.uniform(*program_constants.GAS_GIANT_CORE_ATMOSPHERE_RATIO)
        planet.density = planet.density * core_to_atmosphere_ratio + (
                    1 - core_to_atmosphere_ratio) * (planet.atm_density / 1000)

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
    if planet.planet_class == "M" and (surface_gravity_g < 0.75 or surface_gravity_g > 1.25):
        surface_gravity_g = random.uniform(0.75, 1.25)
    planet.gravity = surface_gravity_g


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
    albedo = random.uniform(0.12, 0.35)
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
        atmosphere_thickness = scale_height * random.uniform(5, 7)
        planet_volume = (4 * math.pi * planet.radius ** 3) / 3
        atmosphere_volume = (4 * math.pi * (planet.radius + atmosphere_thickness) ** 3) / 3 - planet_volume
        atmospheric_mass = 0
        num_zones = round(random.uniform(5, 7)) # Consistent with atmosphere_thickness calculation
        for zone in range(num_zones):
            zone_volume = atmosphere_volume + planet_volume - (
                        4 * math.pi * (planet.radius + (zone * scale_height)) ** 3) / 3
            zone_density = planet.atm_density / (zone * 2.7) if zone >= 1 else planet.atm_density
            atmospheric_mass += zone_volume * zone_density
        atmospheric_force = atmospheric_mass * (planet.gravity * physical_constants.EARTH_GRAVITY)
        planet_surface_area = 4 * math.pi * (planet.radius * 1000) ** 2
        atmospheric_pressure = (atmospheric_force / planet_surface_area) * 7500

        greenhouse_factor = abs((planet.atm_molar_density - physical_constants.CO2_BASE_MOLAR_DENSITY) / physical_constants.CO2_BASE_MOLAR_DENSITY * program_constants.CO2_MAX_GREENHOUSE_FACTOR)
        surface_temperature_atmosphere = ((1 - albedo) * solar_output_at_orbit * (1 + greenhouse_factor) / (4 * physical_constants.STEFAN_BOLTZMANN_CONSTANT)) ** (1 / 4)
        planet.surface_temperature = surface_temperature_atmosphere
        planet.atmospheric_pressure = atmospheric_pressure

        if planet.planet_class == "M":
            if planet.atmospheric_pressure < 90000 or planet.atmospheric_pressure > 112000:
                planet.atmospheric_pressure = random.uniform(90000, 112000)
            if planet.surface_temperature < 283 or planet.surface_temperature > 290:
                planet.surface_temperature = random.uniform(283, 290)
        elif planet.planet_class == "P" and planet.surface_temperature >= 283:
            # If surface_temperature_no_atmosphere is already above 283, we need a different approach
            # to ensure the P class planet remains cold.
            if surface_temperature_no_atmosphere < 283:
                planet.surface_temperature = random.uniform(surface_temperature_no_atmosphere, 283)
            else:
                planet.surface_temperature = random.uniform(200, 283) # A reasonable cold range for P class


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
        moon_distance = random.uniform(total_orbit_distance, high_orbit) / physical_constants.AU_TO_KM
        moon_radius = random.uniform(program_constants.PLANET_CLASSES[moon_class]['radius_range'][0], radius_limit)

        new_moon = Planet(planet.system_config, planet.star, planet.habitable_zone, moon_distance,
                          radius=moon_radius, planet_class=moon_class, zone_override=planet.zone,
                          distance_override=planet.distance, is_moon=True)
        planet.moons.append(new_moon)
        total_orbit_distance = (new_moon.distance * physical_constants.AU_TO_KM) + (new_moon.min_orbit_distance * physical_constants.AU_TO_KM)
