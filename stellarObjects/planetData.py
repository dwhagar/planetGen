# stellarObjects/planetData.py

"""
Planet and Asteroid Belt Generation
===================================

This module contains the `Planet` and `AsteroidBelt` classes, which are used
to generate and represent celestial bodies within a star system. The `Planet`
class is a comprehensive model that includes physical properties, atmospheric
conditions, and orbital characteristics. It can also generate its own system
of moons. The `AsteroidBelt` class provides a simpler representation for
asteroid belts within the system.
"""

import math
import random, secrets

from .config import SystemConfig # Updated import
from stellarObjects import constants
from .evolution import get_evolutionary_timeline
from .names import (MOON_NAMES, MOON_PREFIXES, MOON_SUFFIXES, PLANET_NAMES,
                    PLANET_PREFIXES, PLANET_SUFFIXES)
from .utils import (calc_object_mass, calculate_hill_sphere, _format_age_string,
                    generate_phoneme_salad_name, properties_to_string,
                    to_paragraph, to_scientific_notation, years_to_time_string)


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
    for planet_class, data in constants.PLANET_CLASSES.items():
        min_radius, max_radius = data["radius_range"]
        planet_type = data["type"]

        min_density, max_density = constants.PLANET_DENSITY[planet_type]  # g/cm^3

        # Convert density from g/cm³ to kg/m³ for mass calculation
        min_density *= 1000
        max_density *= 1000

        if planet_type == "t":  # Terrestrial planet
            min_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density
            max_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density
        else:  # Gas giant
            min_atm_density, max_atm_density = constants.ATMOSPHERE_DENSITY[planet_type]
            min_core_ratio, max_core_ratio = constants.GAS_GIANT_CORE_ATMOSPHERE_RATIO

            min_core_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density * min_core_ratio
            max_core_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density * max_core_ratio

            min_atm_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_atm_density * (1 - min_core_ratio)
            max_atm_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_atm_density * (1 - max_core_ratio)

            min_mass = min_core_mass + min_atm_mass
            max_mass = max_core_mass + max_atm_mass

        mass_ranges[planet_class] = (min_mass, max_mass)
    return mass_ranges


planet_mass_ranges = get_planet_mass_ranges()


class Planet:
    """
    A class representing a single planet or moon and all of its properties.

    This class encapsulates the generation and properties of a celestial body,
    including its physical characteristics, orbital data, and atmospheric
    conditions. It can represent both planets and moons, and includes logic
    for generating a system of moons for a planet, influenced by settings in
    the global `config` object.

    The generation process is flexible, allowing for fully random creation or
    generation based on specified parameters like radius, mass, or class.

    Attributes:
        is_moon (bool): True if the object is a moon, False otherwise.
        moons (list): A list of `Planet` objects representing the moons of this planet.
        zone (str): The orbital zone of the planet ('h' for hot, 'e' for ecosphere, 'c' for cold).
        description (str): A text description of the planet's class.
        star_radius (float): The radius of the star in meters.
        star_output (float): The total energy output of the star in Watts.
        star_temperature (float): The surface temperature of the star in Kelvin.
        atm_molar_density (float): The molar density of the atmosphere.
        gravity (float): The surface gravity of the planet in g's.
        atm_density (float): The density of the atmosphere.
        surface_temperature (float): The surface temperature of the planet in Kelvin.
        density (float): The overall density of the planet in g/cm³.
        atmospheric_pressure (float): The atmospheric pressure at the surface in Pascals.
        mass (float): The mass of the planet in kilograms.
        atmosphere (str): The atmospheric composition.
        composition (str): The chemical composition of the planet.
        radius (float): The radius of the planet in kilometers.
        planet_class (str): The classification of the planet (e.g., 'M', 'N').
        distance (float): The distance from the star in AU.
        type (str): The type of planet ('t' for terrestrial, 'g' for gas giant').
        scale_height (float): The atmospheric scale height in kilometers.
        star_mass (float): The mass of the star in kilograms.
        hill_radius (float): The Hill radius of the planet in kilometers.
        min_orbit_distance (float): The minimum stable orbit distance for a satellite in AU.
        name (str): The generated name of the planet or moon.
        life_chemical (str): The primary chemical basis for any potential life.
        hab (tuple): The inner and outer bounds of the habitable zone in AU.
        star_type (str): The spectral type of the star (e.g., 'G', 'M').
        star (Star): The Star object this planet orbits.
        evolutionary_data (list): A list of strings describing the evolutionary timeline.
        flavor_text (str): A randomly selected flavor text for the planet.
        flavor_text_count (int): The number of flavor texts added to this planet.
    """

    def __init__(self, system_config: SystemConfig, star, hab_zone, distance, star_type, star_output, star_radius, star_temperature, star_mass,
                 radius=None, planet_class=None, mass=None, zone_override=None, distance_override=None,
                 is_moon=False):
        """
        Initializes a Planet object with its properties and orbital context.

        Args:
            system_config (SystemConfig): The shared SystemConfig object for the system.
            star (Star): The Star object this planet orbits.
            hab_zone (tuple): A tuple containing the inner and outer bounds of the
                              habitable zone in AU.
            distance (float): The planet's distance from its star in AU.
            star_type (str): The spectral type of the star (e.g., 'G', 'M').
            star_output (float): The total energy output of the star in Watts.
            star_radius (float): The radius of the star in kilometers.
            star_temperature (float): The surface temperature of the star in Kelvin.
            star_mass (float): The mass of the star in kilograms.
            radius (float, optional): The radius of the planet in kilometers.
            planet_class (str, optional): The class of the planet (e.g., 'M', 'N').
            mass (float, optional): The mass of the planet in kilograms.
            zone_override (str, optional): A single character ('h', 'c', 'e') to
                                           override the planet's calculated zone.
            distance_override (float, optional): An override for the planet's
                                                 distance, used in specific calculations.
            is_moon (bool, optional): Flag indicating if the object is a moon.
        """
        self.system_config = system_config # Store SystemConfig
        self.is_moon = is_moon
        self.moons = []
        self.zone = None
        self.description = None
        self.star_radius = star_radius * 1000  # in meters
        self.star_output = star_output  # in Watts
        self.star_temperature = star_temperature  # in Kelvin
        self.atm_molar_density = None
        self.gravity = None
        self.atm_density = None
        self.surface_temperature = None
        self.density = None
        self.atmospheric_pressure = None
        self.mass = mass
        self.atmosphere = None
        self.composition = None
        self.radius = radius
        self.planet_class = planet_class
        self.distance = distance
        self.type = None
        self.scale_height = None
        self.star_mass = star_mass
        self.hill_radius = None
        self.min_orbit_distance = None
        self.name = None
        self.life_chemical = None
        self.evolution = None
        self.reflection_spectrum_visible = None
        self.reflection_spectrum_non_visible = None
        self.evolutionary_data = [] # Initialize evolutionary data
        self.flavor_text = None # Initialize flavor text
        self.flavor_text_count = 0 # Initialize flavor text count for this planet

        # From the star, should not be changed.
        self.hab = hab_zone
        self.star_type = star_type
        self.star = star # Store the Star object

        if self.is_moon:
            self.name = generate_phoneme_salad_name(MOON_NAMES, MOON_PREFIXES, MOON_SUFFIXES)
        else:
            self.name = generate_phoneme_salad_name(PLANET_NAMES, PLANET_PREFIXES, PLANET_SUFFIXES)

        # Calculate additional properties.
        self.generate_planet(zone_override)
        self.volume = (4 / 3) * math.pi * self.radius ** 3  # Calculate volume in km^3
        self.period = math.sqrt(self.distance ** 3)
        self.calculate_surface_gravity()
        self.calculate_atmospheric_conditions(distance_override)

        # Generate evolutionary timeline data if the planet is habitable and not a moon
        if self.zone == 'e' and not self.is_moon: # Only for habitable planets
            self.evolutionary_data = get_evolutionary_timeline(self.star)

        if (self.system_config.FORCE_MOONS or secrets.randbelow(2) == 1) and not self.is_moon: # Use self.system_config
            self.generate_moons()

    def generate_planet(self, zone_override=None):
        """
        Generates random planet properties based on its location and any provided constraints.

        This method is the core of the planet generation logic. It determines the
        planet's properties such as class, composition, and atmosphere. The generation
        can be fully random or guided by inputs like a specific radius, mass, or
        planet class. It ensures that the generated properties are consistent with
        each other and the planet's orbital zone.

        Args:
            zone_override (str, optional): A character ('h', 'c', 'e') to manually
                                           set the planet's zone, overriding the
                                           calculation based on distance.
        """
        random.seed(secrets.randbits(128)) # Re-seed at the start of the function
        # Determine the planet's zone (hot, ecosphere, or cold)
        inner_bound, outer_bound = self.hab
        if self.distance < inner_bound:
            zone = 'h'
        elif self.distance > outer_bound:
            zone = 'c'
        else:
            zone = 'e'

        if zone_override and zone_override.lower() in "hce":
            zone = zone_override.lower()
        self.zone = zone

        # --- Input Validation and Random Generation ---
        # This section handles the logic for generating planet properties based on
        # the inputs provided. It can generate a fully random planet, or generate
        # properties based on a given class, radius, or mass.

        if self.planet_class is None and self.radius is None and self.mass is None:
            # Fully random generation
            valid_classes = [c for c, data in constants.PLANET_CLASSES.items() if data[zone]]
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e': # Use self.system_config
                valid_classes = [c for c in valid_classes if c not in constants.HABITABLE_PLANET_CLASSES]
            
            classes = list(constants.PLANET_CLASS_PROBABILITIES.keys())
            probabilities = list(constants.PLANET_CLASS_PROBABILITIES.values())
            class_valid = False
            while not class_valid:
                class_choice = random.choices(classes, weights=probabilities, k=1)[0]
                if class_choice in valid_classes:
                    self.planet_class = class_choice
                    class_valid = True
            min_radius, max_radius = constants.PLANET_CLASSES[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        elif self.planet_class is not None and self.radius is None and self.mass is None:
            # Class given, generate radius
            self._validate_planet_class(zone)
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e' and self.planet_class in constants.HABITABLE_PLANET_CLASSES: # Use self.system_config
                raise ValueError(f"Cannot generate habitable planet class {self.planet_class} in ecosphere when NO_HABITABLE_WORLD is True.")
            min_radius, max_radius = constants.PLANET_CLASSES[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        elif self.planet_class is None and self.radius is not None and self.mass is None:
            # Radius given, determine possible classes
            possible_classes = [c for c, data in constants.PLANET_CLASSES.items()
                                if data[zone] and data["radius_range"][0] <= self.radius <= data["radius_range"][1]]
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e': # Use self.system_config
                possible_classes = [c for c in possible_classes if c not in constants.HABITABLE_PLANET_CLASSES]
            if not possible_classes:
                raise ValueError("No valid planet class for the given radius in this zone")
            self.planet_class = secrets.choice(possible_classes)
            self._validate_radius()

        elif self.planet_class is None and self.radius is None and self.mass is not None:
            # Mass given, determine possible classes
            possible_classes = [c for c, data in constants.PLANET_CLASSES.items()
                                if planet_mass_ranges[c][0] <= self.mass <= planet_mass_ranges[c][1] and data[zone]]
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e': # Use self.system_config
                possible_classes = [c for c in possible_classes if c not in constants.HABITABLE_PLANET_CLASSES]
            if not possible_classes:
                raise ValueError("No valid planet class for the given mass in this zone")
            self.planet_class = secrets.choice(possible_classes)
            self._validate_mass()

        elif self.planet_class is not None and self.radius is not None and self.mass is None:
            # Class and radius given, validate
            self._validate_planet_class(zone)
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e' and self.planet_class in constants.HABITABLE_PLANET_CLASSES: # Use self.system_config
                raise ValueError(f"Cannot generate habitable planet class {self.planet_class} in ecosphere when NO_HABITABLE_WORLD is True.")
            self._validate_radius()

        elif self.planet_class is not None and self.radius is None and self.mass is not None:
            # Class and mass given, validate and generate radius
            self._validate_planet_class(zone)
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e' and self.planet_class in constants.HABITABLE_PLANET_CLASSES: # Use self.system_config
                raise ValueError(f"Cannot generate habitable planet class {self.planet_class} in ecosphere when NO_HABITABLE_WORLD is True.")
            self._validate_mass()
            min_radius, max_radius = constants.PLANET_CLASSES[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        elif self.planet_class is None and self.radius is not None and self.mass is not None:
            # Radius and mass given, determine possible classes
            possible_classes = []
            for c, data in constants.PLANET_CLASSES.items():
                min_mass, max_mass = planet_mass_ranges[c]
                min_radius, max_radius = data["radius_range"]
                if min_mass <= self.mass <= max_mass and min_radius <= self.radius <= max_radius and data[zone]:
                    possible_classes.append(c)
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e': # Use self.system_config
                possible_classes = [c for c in possible_classes if c not in constants.HABITABLE_PLANET_CLASSES]
            if not possible_classes:
                raise ValueError("No valid planet class for the given radius/mass in this zone")
            self.planet_class = secrets.choice(possible_classes)
            self._validate_radius()
            self._validate_mass()

        else:
            # All inputs provided, fully validate
            self._validate_planet_class(zone)
            if self.system_config.NO_HABITABLE_WORLD and zone == 'e' and self.planet_class in constants.HABITABLE_PLANET_CLASSES: # Use self.system_config
                raise ValueError(f"Cannot generate habitable planet class {self.planet_class} in ecosphere when NO_HABITABLE_WORLD is True.")
            self._validate_radius()
            self._validate_mass()

        class_data = constants.PLANET_CLASSES[self.planet_class]
        self.composition = class_data["composition"]
        self.description = class_data["description"]
        self.type = class_data["type"]

        min_density, max_density = constants.PLANET_DENSITY[self.type]
        self.density = random.uniform(min_density, max_density)

        if class_data["atmosphere"] is None:
            self.atmosphere = "None"
        else:
            self.atmosphere = class_data["atmosphere"]
            if self.planet_class == 'N':
                self.atm_density = 65
                min_am_density, max_am_density = constants.ATMOSPHERIC_MOLAR_DENSITY[self.type]
                self.atm_molar_density = max_am_density
            else:
                min_a_density, max_a_density = constants.ATMOSPHERE_DENSITY[self.type]
                self.atm_density = random.uniform(min_a_density, max_a_density)
                min_am_density, max_am_density = constants.ATMOSPHERIC_MOLAR_DENSITY[self.type]
                self.atm_molar_density = random.uniform(min_am_density, max_am_density)

        self.volume, self.mass = calc_object_mass(self.planet_class, self.radius, constants.PLANET_CLASSES, constants.PLANET_DENSITY,
                                                  self.density)

        if self.type == 'g':
            core_to_atmosphere_ratio = random.uniform(*constants.GAS_GIANT_CORE_ATMOSPHERE_RATIO)
            self.density = self.density * core_to_atmosphere_ratio + (
                        1 - core_to_atmosphere_ratio) * (self.atm_density / 1000)

        self.volume, self.mass = calc_object_mass(self.planet_class, self.radius, constants.PLANET_CLASSES, constants.PLANET_DENSITY,
                                                  self.density)
        
        # Get the dictionary of viable chemicals and their weights
        viable_chems = self.get_viable_life_chemicals()
        
        if viable_chems:
            # random.choices returns a list, so we grab the first [0] element
            self.life_chemical = random.choices(
                population=list(viable_chems.keys()),
                weights=list(viable_chems.values()),
                k=1
            )[0]
            life_chem_data = constants.LIFE_CHEMICALS.get(self.life_chemical, {})
            self.reflection_spectrum_visible = life_chem_data.get("reflection_spectrum_visible")
            self.reflection_spectrum_non_visible = life_chem_data.get("reflection_spectrum_non_visible")
        else:
            self.life_chemical = None

        # Determine the evolutionary speed based on the star and chosen chemical
        self.evolution = self.get_evolutionary_speed()

        distance_m = self.distance * constants.AU_TO_M
        self.hill_radius = calculate_hill_sphere(distance_m, self.mass, self.star_mass) / 1000  # Convert to km
        self.min_orbit_distance = (5 * self.hill_radius) / constants.AU_TO_KM

    def calculate_surface_gravity(self):
        """
        Calculates the surface gravity of the planet in g's (multiples of Earth's gravity).

        This method computes the surface gravity based on the planet's mass and
        radius. The result is normalized to Earth's gravity (g's). It includes
        special adjustments for certain planet classes to ensure realistic values.
        """
        random.seed(secrets.randbits(128)) # Re-seed at the start of the function
        radius_meters = self.radius * 1000
        surface_gravity = (constants.G * self.mass) / (radius_meters ** 2)
        surface_gravity_g = surface_gravity / constants.EARTH_GRAVITY
        if surface_gravity_g <= 0:
            raise ValueError('Invalid value for gravity.')
        if self.planet_class == "M" and (surface_gravity_g < 0.75 or surface_gravity_g > 1.25):
            surface_gravity_g = random.uniform(0.75, 1.25)
        self.gravity = surface_gravity_g

    def get_viable_life_chemicals(self):
        """
        Determines the viable life chemicals for the planet by finding the
        intersection of chemicals supported by both the planet's class and
        the star's spectral type. Normalizes the probabilities to sum to 1.0.

        Returns:
            dict: A dictionary mapping viable life chemical strings to their
                  normalized float probabilities (e.g., {"Melanin": 0.6}).
        """
        # Retrieve the base list of possible chemicals for this specific planet class
        planet_data = constants.PLANET_CLASSES.get(self.planet_class)
        if not planet_data or not planet_data.get("life_chemical"):
            return {}

        planet_chems = planet_data["life_chemical"]

        # Determine the star's spectral class (e.g., 'G' from 'G2V')
        if not self.star_type:
            return {}

        spectral_class = self.star_type[0].upper()

        # Retrieve the potentially viable chemicals for the star's evolutionary scale
        star_data = constants.STAR_EVOLUTION.get(spectral_class)
        if not star_data or not star_data.get("supported_evolutionary_scales"):
            return {}

        star_chems = star_data["potentially_viable_chemicals"]

        raw_probabilities = {}
        total_raw_prob = 0

        # Intersect the lists using substring matching and collect raw probabilities
        for p_chem in planet_chems:
            if any(p_chem in s_chem for s_chem in star_chems):
                chem_prob = constants.LIFE_CHEMICALS.get(p_chem, {}).get(
                    "star_spectra_probabilities", {}).get(spectral_class, 0)

                if chem_prob > 0:
                    raw_probabilities[p_chem] = chem_prob
                    total_raw_prob += chem_prob

        # Normalize the probabilities proportionally so they sum to 1.0
        normalized_chemicals = {}
        if total_raw_prob > 0:
            for chem, prob in raw_probabilities.items():
                normalized_chemicals[chem] = prob / total_raw_prob

        return normalized_chemicals

    def get_evolutionary_speed(self):
        """
        Determines the evolutionary speed for the planet's biosphere.
        It finds the intersection of speeds supported by the star's lifespan
        and the speed required by the planet's specific life chemical, then
        chooses one randomly.

        Returns:
            str: The chosen evolutionary speed (e.g., 'fast', 'normal', 'slow'),
                 or None if data is missing.
        """
        if not self.star_type:
            return None

        spectral_class = self.star_type[0].upper()

        # 1. Get the speeds supported by the star
        star_data = constants.STAR_EVOLUTION.get(spectral_class)
        if not star_data or not star_data.get("supported_evolutionary_scales"):
            return None

        star_speeds = star_data["supported_evolutionary_scales"]

        # 2. Get the speeds supported by the chemical (if one is assigned)
        if self.life_chemical:
            chem_data = constants.LIFE_CHEMICALS.get(self.life_chemical)
            if chem_data and chem_data.get("evolutionary_time_scale"):
                chem_speeds = chem_data["evolutionary_time_scale"]

                # Ensure it's a list for intersection logic, even though
                # constants.py currently stores it as a single string
                if isinstance(chem_speeds, str):
                    chem_speeds = [chem_speeds]

                # Find the intersection
                valid_speeds = [speed for speed in star_speeds if speed in chem_speeds]

                if valid_speeds:
                    return secrets.choice(valid_speeds)

        # 3. Fallback: If no chemical is set (or if there's somehow no overlap),
        # just pick a random speed supported by the star.
        if star_speeds:
            return secrets.choice(star_speeds)

        return None

    def calculate_atmospheric_conditions(self, distance_override=None):
        """
        Calculates the atmospheric conditions of the planet, including surface
        temperature and atmospheric pressure.

        This method models the planet's atmospheric conditions. It calculates the
        surface temperature considering the star's output, the planet's distance,
        and the greenhouse effect of its atmosphere. It also estimates the
        atmospheric pressure based on the atmospheric mass and planet's gravity.

        Args:
            distance_override (float, optional): An override for the planet's
                                                 distance, used for special cases
                                                 like moons.
        """
        random.seed(secrets.randbits(128)) # Re-seed at the start of the function
        distance = float(distance_override) if distance_override is not None else float(self.distance)
        orbital_radius_km = distance * constants.AU_TO_KM
        output_area = 4 * math.pi * orbital_radius_km ** 2
        solar_output_at_orbit = (self.star_output / output_area) / 1e6
        albedo = random.uniform(0.12, 0.35)
        surface_temperature_no_atmosphere = (
                                                    (1 - albedo) * solar_output_at_orbit / (4 * constants.STEFAN_BOLTZMANN_CONSTANT)) ** (
                                                        1 / 4)

        if self.atmosphere == "None":
            self.surface_temperature = surface_temperature_no_atmosphere
            self.atmospheric_pressure = 0.0
        else:
            scale_height = (constants.R * surface_temperature_no_atmosphere) / (
                        self.atm_molar_density * self.gravity * constants.EARTH_GRAVITY)
            self.scale_height = scale_height
            atmosphere_thickness = scale_height * random.uniform(5, 7)
            planet_volume = (4 * math.pi * self.radius ** 3) / 3
            atmosphere_volume = (4 * math.pi * (self.radius + atmosphere_thickness) ** 3) / 3 - planet_volume
            atmospheric_mass = 0
            num_zones = round(random.uniform(5, 7)) # Consistent with atmosphere_thickness calculation
            for zone in range(num_zones):
                zone_volume = atmosphere_volume + planet_volume - (
                            4 * math.pi * (self.radius + (zone * scale_height)) ** 3) / 3
                zone_density = self.atm_density / (zone * 2.7) if zone >= 1 else self.atm_density
                atmospheric_mass += zone_volume * zone_density
            atmospheric_force = atmospheric_mass * (self.gravity * constants.EARTH_GRAVITY)
            planet_surface_area = 4 * math.pi * (self.radius * 1000) ** 2
            atmospheric_pressure = (atmospheric_force / planet_surface_area) * 7500
            
            greenhouse_factor = abs((self.atm_molar_density - constants.CO2_BASE_MOLAR_DENSITY) / constants.CO2_BASE_MOLAR_DENSITY * constants.CO2_MAX_GREENHOUSE_FACTOR)
            surface_temperature_atmosphere = ((1 - albedo) * solar_output_at_orbit * (1 + greenhouse_factor) / (4 * constants.STEFAN_BOLTZMANN_CONSTANT)) ** (1 / 4)
            self.surface_temperature = surface_temperature_atmosphere
            self.atmospheric_pressure = atmospheric_pressure

            if self.planet_class == "M":
                if self.atmospheric_pressure < 90000 or self.atmospheric_pressure > 112000:
                    self.atmospheric_pressure = random.uniform(90000, 112000)
                if self.surface_temperature < 283 or self.surface_temperature > 290:
                    self.surface_temperature = random.uniform(283, 290)
            elif self.planet_class == "P" and self.surface_temperature >= 283:
                # If surface_temperature_no_atmosphere is already above 283, we need a different approach
                # to ensure the P class planet remains cold.
                if surface_temperature_no_atmosphere < 283:
                    self.surface_temperature = random.uniform(surface_temperature_no_atmosphere, 283)
                else:
                    self.surface_temperature = random.uniform(200, 283) # A reasonable cold range for P class

    def generate_moons(self):
        """
        Generates a system of moons for the given planet.

        This method procedurally generates moons for the planet. It determines the
        number, size, and orbital distance of the moons based on the planet's
        properties, such as its mass and Hill radius. The generated moons are
        themselves instances of the `Planet` class, with the `is_moon` flag set.
        """
        random.seed(secrets.randbits(128)) # Re-seed at the start of the function
        max_moon_mass = self.mass / 10
        max_moon_radius = self.radius / (10 ** (1 / 3))
        possible_classes = [c for c, data in planet_mass_ranges.items()
                            if constants.PLANET_CLASSES[c][self.zone] and constants.PLANET_CLASSES[c]["type"] == 't' and c not in constants.MOON_BLACKLIST
                            and data[1] <= max_moon_mass and constants.PLANET_CLASSES[c]['radius_range'][1] <= max_moon_radius]
        if not possible_classes:
            return

        low_orbit = self.scale_height * 15 if self.scale_height else 100
        high_orbit = self.min_orbit_distance * constants.AU_TO_KM
        total_orbit_distance = low_orbit

        while total_orbit_distance < high_orbit and total_orbit_distance < (self.distance * constants.AU_TO_KM):
            classes = list(constants.PLANET_CLASS_PROBABILITIES.keys())
            probabilities = list(constants.PLANET_CLASS_PROBABILITIES.values())
            moon_class = random.choices(classes, weights=probabilities, k=1)[0]
            while moon_class not in possible_classes:
                moon_class = random.choices(classes, weights=probabilities, k=1)[0]

            radius_limit = constants.PLANET_CLASSES[moon_class]['radius_range'][1] if max_moon_radius > \
                                                                            constants.PLANET_CLASSES[moon_class]['radius_range'][
                                                                                1] else max_moon_radius
            moon_distance = random.uniform(total_orbit_distance, high_orbit) / constants.AU_TO_KM
            moon_radius = random.uniform(constants.PLANET_CLASSES[moon_class]['radius_range'][0], radius_limit)

            new_moon = Planet(self.system_config, self.star, self.hab, moon_distance, self.star_type, self.star_output, self.star_radius, # Pass system_config
                              self.star_temperature, self.star_mass, radius=moon_radius,
                              planet_class=moon_class, zone_override=self.zone,
                              distance_override=self.distance, is_moon=True)
            self.moons.append(new_moon)
            total_orbit_distance = (new_moon.distance * constants.AU_TO_KM) + (new_moon.min_orbit_distance * constants.AU_TO_KM)

    def _validate_planet_class(self, zone):
        """
        Validates if the planet's class is valid for its zone.

        Args:
            zone (str): The zone ('h', 'c', 'e') of the planet.

        Raises:
            ValueError: If the planet class is not valid for the given zone.
        """
        if self.planet_class not in constants.PLANET_CLASSES or not constants.PLANET_CLASSES[self.planet_class][zone]:
            raise ValueError("Invalid planet class for this zone")

    def _validate_radius(self):
        """
        Validates if the planet's radius is within the allowed range for its class.

        Raises:
            ValueError: If the radius is outside the valid range for the planet's class.
        """
        min_radius, max_radius = constants.PLANET_CLASSES[self.planet_class]["radius_range"]
        if not (min_radius <= self.radius <= max_radius):
            raise ValueError("Invalid radius for planet class")

    def _validate_mass(self):
        """
        Validates if the planet's mass is within the allowed range for its class.

        Raises:
            ValueError: If the mass is outside the valid range for the planet's class.
        """
        min_mass, max_mass = planet_mass_ranges[self.planet_class]
        if not (min_mass <= self.mass <= max_mass):
            raise ValueError("Invalid mass for planet class")

    def _generate_life_and_flavor_paragraphs(self, object_type_desc, sentences):
        """
        Generates paragraphs related to life, evolution, and flavor text for the planet.

        Args:
            object_type_desc (str): Description of the object type ("planet" or "moon").
            sentences (list): A list of sentences to append life-related descriptions to.

        Returns:
            list: A list of generated paragraphs.
        """
        random.seed(secrets.randbits(128)) # Re-seed at the start of the function
        life_paragraphs = []

        if self.life_chemical:
            if self.evolution:
                sentences.append(
                    f"The {object_type_desc}'s conditions are suitable for the development of life based on {self.life_chemical.lower()}, which is expected to evolve at a {self.evolution} pace.")
            else:
                sentences.append(
                    f"The {object_type_desc}'s conditions are suitable for the development of life based on {self.life_chemical.lower()}.")
            if self.reflection_spectrum_visible:
                visible_spectrum = ", ".join(self.reflection_spectrum_visible)
                sentences.append(
                    f"The visible reflection spectrum of the life on this {object_type_desc} is characterized by {visible_spectrum.lower()}.")
            if self.reflection_spectrum_non_visible:
                non_visible_spectrum = ", ".join(self.reflection_spectrum_non_visible)
                sentences.append(f"In the non-visible spectrum, it exhibits {non_visible_spectrum.lower()}.")

        life_paragraphs.append(to_paragraph(sentences))

        # Add the evolutionary timeline data if available
        if self.evolutionary_data and self.planet_class in constants.HABITABLE_PLANET_CLASSES:
            life_paragraphs.extend(self.evolutionary_data)

        # Add flavor text if the random chance passes and limits are not exceeded
        if random.random() < constants.FLAVOR_CHANCE_PLANET and self.system_config.system_flavor_count < constants.MAX_FLAVOR_TOTAL:
            selected_flavor = None
            # Check for habitable and multicellular/technological life
            is_habitable = self.planet_class in constants.HABITABLE_PLANET_CLASSES
            has_multicellular_life = False
            if is_habitable and self.evolutionary_data:
                for stage_paragraph in self.evolutionary_data:
                    if "multicellularity" in stage_paragraph.lower() or "technological_civilization" in stage_paragraph.lower():
                        has_multicellular_life = True
                        break

            if is_habitable and has_multicellular_life:
                selected_flavor = secrets.choice(constants.HABITABLE_FLAVOR)
            elif self.type == "t" and self.planet_class != "A":
                selected_flavor = secrets.choice(constants.PLANET_FLAVOR)
            elif self.type == "g" or self.planet_class == "A":
                selected_flavor = secrets.choice(constants.ORBITAL_FLAVOR)

            if selected_flavor:
                self.flavor_text = selected_flavor
                life_paragraphs.append(f"Sensors show {self.flavor_text}")
                self.system_config.system_flavor_count += 1
                self.flavor_text_count += 1
        
        return life_paragraphs

    def to_paragraph_list(self):
        """
        Returns a list of strings, where each string is a paragraph describing
        the planet.
        """
        object_type_desc = "moon" if self.is_moon else "planet"

        if self.is_moon:
            # Moons orbit their parent planet, so their distance is from the planet, not the star.
            # This distance is typically much smaller and best represented in kilometers.
            distance_text = f"{to_scientific_notation(self.system_config, self.distance * constants.AU_TO_KM, 4)} km" # Pass system_config
            header_level = '###' if self.system_config.MARKDOWN else '===' # Use self.system_config
        else:
            # For planets, the distance is from the star. We check if this distance
            # is large enough to warrant using light-years.
            distance_ly = self.distance * constants.AU_TO_LY
            if distance_ly < constants.LY_THRESHOLD:
                # For "normal" sized systems, display in AU.
                # If less than 1 AU, also show in km for context.
                distance_text = f"{to_scientific_notation(self.system_config, self.distance * constants.AU_TO_KM, 1)} km ({self.distance:.3f} AU)" if self.distance < 1 else f"{self.distance:.3f} AU" # Pass system_config
            else:
                # For very large systems, display in light-years.
                distance_text = f"{distance_ly:.4f} light-years"
            header_level = '##' if self.system_config.MARKDOWN else '==' # Use self.system_config

        output_paragraphs = []
        header = f"{header_level} {self.name} {header_level if not self.system_config.MARKDOWN else ''}" # Use self.system_config
        output_paragraphs.append(header)

        radius_string = f"{round(self.radius, 2):,} km" if self.radius <= 100000 else f"{to_scientific_notation(self.system_config, self.radius, 2)} km" # Pass system_config

        planet_properties = {
            "class": self.planet_class,
            "distance": distance_text,
            "period": years_to_time_string(self.period),
            "radius": radius_string,
            "gravity": f"{round(self.gravity, 3)} g",
        }

        output_paragraphs.append(properties_to_string(self.system_config, planet_properties, "Planet Data")) # Pass system_config

        sentences = []
        if not self.is_moon:
            if len(self.moons) > 1:
                sentences.append(f"There are {len(self.moons)} moons orbiting this planet.")
            elif len(self.moons) == 1:
                sentences.append("There is 1 moon orbiting this planet.")
            else:
                sentences.append("There are no moons orbiting this planet.")

        if self.type == "t":
            if self.atmosphere != "None":
                sentences.append(
                    f"This {object_type_desc} has a surface pressure of {self.atmospheric_pressure / 1000:.1f} kPa or {self.atmospheric_pressure / 101300:.2f} atmospheres and a temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
                sentences.append(
                    f"It is {self.description.lower()} with an atmosphere of {self.atmosphere.lower()} and a composition of {self.composition.lower()}.")
            else:
                sentences.append(
                    f"This {object_type_desc} has no atmosphere and a surface temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
                sentences.append(f"It is {self.description.lower()} with a composition of {self.composition.lower()}.")
        else:
            sentences.append(
                f"This gas giant has an internal pressure of {self.atmospheric_pressure / 1000:.1f} kPa or {self.atmospheric_pressure / 101300:.1f} atmospheres and a temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
            sentences.append(f"It is {self.description.lower()} with a composition of {self.composition.lower()}.")

        # Call the new method to get life and flavor text paragraphs
        output_paragraphs.extend(self._generate_life_and_flavor_paragraphs(object_type_desc, sentences))

        if self.moons:
            for moon in self.moons:
                output_paragraphs.extend(moon.to_paragraph_list())  # Recursively get moon paragraphs

        # Post-processing: Replace "planet" with "moon" if the object is a moon
        if self.is_moon:
            processed_paragraphs = []
            for paragraph in output_paragraphs:
                paragraph = paragraph.replace("planet", "moon")
                paragraph = paragraph.replace("Planet", "Moon")
                paragraph = paragraph.replace("{{Moon Data", "{{Class Data")
                processed_paragraphs.append(paragraph)
            output_paragraphs = processed_paragraphs

        return output_paragraphs

    def __str__(self):
        """
        Returns a string representation of the planet, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())