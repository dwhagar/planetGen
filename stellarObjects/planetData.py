# stellarObjects/planetData.py

"""
Planet and Asteroid Belt Generation
===================================

This module contains the `Planet` and `Asteroid_Belt` classes, which are used
to generate and represent celestial bodies within a star system.
"""

import math
import random
from .utils import to_scientific_notation, years_to_time_string, calc_object_mass, generate_phoneme_salad_name
from .constants import (EARTH_RADIUS_KM, EARTH_GRAVITY, AU_TO_KM, G, R,
                        STEFAN_BOLTZMANN_CONSTANT, PLANET_DENSITY, ATMOSPHERE_DENSITY,
                        ATMOSPHERIC_MOLAR_DENSITY, GAS_GIANT_CORE_ATMOSPHERE_RATIO,
                        PLANET_CLASSES, PLANET_CLASS_PROBABILITIES)
from .names import (PLANET_NAMES, MOON_NAMES, PLANET_PREFIXES, PLANET_SUFFIXES, 
                    MOON_PREFIXES, MOON_SUFFIXES)

def get_planet_mass_ranges():
    """
    Calculates and returns a dictionary of valid mass ranges (in kg) for each planet class.

    This function pre-calculates the minimum and maximum possible mass for each
    planet class based on its radius and density ranges. This is used for
    validation when generating planets with a given mass.

    Returns:
        dict: A dictionary where keys are planet classes and values are tuples
              of (min_mass, max_mass).
    """
    mass_ranges = {}
    for planet_class, data in PLANET_CLASSES.items():
        min_radius, max_radius = data["radius_range"]
        planet_type = data["type"]

        min_density, max_density = PLANET_DENSITY[planet_type]  # g/cm^3

        # Convert density from g/cm³ to kg/m³ for mass calculation
        min_density *= 1000
        max_density *= 1000

        if planet_type == "t":  # Terrestrial planet
            min_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density
            max_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density
        else:  # Gas giant
            min_atm_density, max_atm_density = ATMOSPHERE_DENSITY[planet_type]
            min_core_ratio, max_core_ratio = GAS_GIANT_CORE_ATMOSPHERE_RATIO

            min_core_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density * min_core_ratio
            max_core_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density * max_core_ratio

            min_atm_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_atm_density * (1 - min_core_ratio)
            max_atm_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_atm_density * (1 - max_core_ratio)

            min_mass = min_core_mass + min_atm_mass
            max_mass = max_core_mass + max_atm_mass

        mass_ranges[planet_class] = (min_mass, max_mass)
    return mass_ranges

planet_mass_ranges = get_planet_mass_ranges()

class Asteroid_Belt:
    """
    A basic class to store information for an asteroid belt.
    """
    def __init__(self, distance, lower_limit, upper_limit):
        self.distance = distance
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.type = 'a'

    def __str__(self):
        return f"== Asteroid Belt ==\nAn asteroid belt orbits roughly between {self.lower_limit:.3f} AU and {self.upper_limit:.3f} AU."

class Planet:
    """
    A class representing a single planet or moon and all of its properties.
    """
    def __init__(self, hab_zone, distance, star_output, star_radius, star_temperature, star_mass,
                 radius=None, planet_class=None, mass=None, zone_override=None, distance_override=None,
                 is_moon=False, force_moons=False):
        """
        Initializes a Planet object with its radius and the spectral class of its host star.
        """
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

        # From the star, should not be changed.
        self.hab = hab_zone

        if self.is_moon:
            self.name = generate_phoneme_salad_name(MOON_NAMES, MOON_PREFIXES, MOON_SUFFIXES)
            force_moons = False
        else:
            self.name = generate_phoneme_salad_name(PLANET_NAMES, PLANET_PREFIXES, PLANET_SUFFIXES)

        # Calculate additional properties.
        self.generate_planet(zone_override)
        self.volume = (4/3) * math.pi * self.radius**3  # Calculate volume in km^3
        self.period = math.sqrt(self.distance**3)
        self.calculate_surface_gravity()
        self.calculate_atmospheric_conditions(distance_override)

        if (force_moons or random.randint(0, 1) == 1) and not self.is_moon:
            self.generate_moons()

    def generate_planet(self, zone_override=None):
        """
        Generates random planet properties (class, composition, atmosphere, mass)
        based on distance and optional radius/class/mass inputs.
        """
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
            valid_classes = [c for c, data in PLANET_CLASSES.items() if data[zone]]
            classes = list(PLANET_CLASS_PROBABILITIES.keys())
            probabilities = list(PLANET_CLASS_PROBABILITIES.values())
            class_valid = False
            while not class_valid:
                class_choice = random.choices(classes, weights=probabilities, k=1)[0]
                if class_choice in valid_classes:
                    self.planet_class = class_choice
                    class_valid = True
            min_radius, max_radius = PLANET_CLASSES[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        elif self.planet_class is not None and self.radius is None and self.mass is None:
            # Class given, generate radius
            self._validate_planet_class(zone)
            min_radius, max_radius = PLANET_CLASSES[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        elif self.planet_class is None and self.radius is not None and self.mass is None:
            # Radius given, determine possible classes
            possible_classes = [c for c, data in PLANET_CLASSES.items()
                                if data[zone] and data["radius_range"][0] <= self.radius <= data["radius_range"][1]]
            if not possible_classes:
                raise ValueError("No valid planet class for the given radius in this zone")
            self.planet_class = random.choice(possible_classes)
            self._validate_radius()

        elif self.planet_class is None and self.radius is None and self.mass is not None:
            # Mass given, determine possible classes
            possible_classes = [c for c, data in PLANET_CLASSES.items()
                                if planet_mass_ranges[c][0] <= self.mass <= planet_mass_ranges[c][1] and data[zone]]
            if not possible_classes:
                raise ValueError("No valid planet class for the given mass in this zone")
            self.planet_class = random.choice(possible_classes)
            self._validate_mass()

        elif self.planet_class is not None and self.radius is not None and self.mass is None:
            # Class and radius given, validate
            self._validate_planet_class(zone)
            self._validate_radius()

        elif self.planet_class is not None and self.radius is None and self.mass is not None:
            # Class and mass given, validate and generate radius
            self._validate_planet_class(zone)
            self._validate_mass()
            min_radius, max_radius = PLANET_CLASSES[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)
            self._validate_radius()

        elif self.planet_class is None and self.radius is not None and self.mass is not None:
            # Radius and mass given, determine possible classes
            possible_classes = []
            for c, data in PLANET_CLASSES.items():
                min_mass, max_mass = planet_mass_ranges[c]
                min_radius, max_radius = data["radius_range"]
                if min_mass <= self.mass <= max_mass and min_radius <= self.radius <= max_radius and data[zone]:
                    possible_classes.append(c)
            if not possible_classes:
                raise ValueError("No valid planet class for the given radius/mass in this zone")
            self.planet_class = random.choice(possible_classes)
            self._validate_radius()
            self._validate_mass()

        else:
            # All inputs provided, fully validate
            self._validate_planet_class(zone)
            self._validate_radius()
            self._validate_mass()

        class_data = PLANET_CLASSES[self.planet_class]
        self.composition = class_data["composition"]
        self.description = class_data["description"]
        self.type = class_data["type"]

        min_density, max_density = PLANET_DENSITY[self.type]
        self.density = random.uniform(min_density, max_density)

        if class_data["atmosphere"] is None:
            self.atmosphere = "None"
        else:
            self.atmosphere = class_data["atmosphere"]
            if self.planet_class == 'N':
                self.atm_density = 65
                min_am_density, max_am_density = ATMOSPHERIC_MOLAR_DENSITY[self.type]
                self.atm_molar_density = max_am_density
            else:
                min_a_density, max_a_density = ATMOSPHERE_DENSITY[self.type]
                self.atm_density = random.uniform(min_a_density, max_a_density)
                min_am_density, max_am_density = ATMOSPHERIC_MOLAR_DENSITY[self.type]
                self.atm_molar_density = random.uniform(min_am_density, max_am_density)

        self.volume, self.mass = calc_object_mass(self.planet_class, self.radius, PLANET_CLASSES, PLANET_DENSITY, self.density)

        if self.type == 'g':
            core_to_atmosphere_ratio = random.uniform(*GAS_GIANT_CORE_ATMOSPHERE_RATIO)
            self.density = self.density * core_to_atmosphere_ratio + (1 - core_to_atmosphere_ratio) * (self.atm_density / 1000)

        self.volume, self.mass = calc_object_mass(self.planet_class, self.radius, PLANET_CLASSES, PLANET_DENSITY, self.density)

        self.hill_radius = (self.distance * AU_TO_KM) * (self.mass / (3 * self.star_mass)) ** (1 / 3)
        self.min_orbit_distance = (5 * self.hill_radius) / AU_TO_KM

    def calculate_surface_gravity(self):
        """
        Calculates the surface gravity of the planet in g's (Earth's gravity).
        """
        radius_meters = self.radius * 1000
        surface_gravity = (G * self.mass) / (radius_meters ** 2)
        surface_gravity_g = surface_gravity / EARTH_GRAVITY
        if surface_gravity_g <= 0:
            raise ValueError('Invalid value for gravity.')
        if self.planet_class == "M" and (surface_gravity_g < 0.75 or surface_gravity_g > 1.25):
            surface_gravity_g = random.uniform(0.75, 1.25)
        self.gravity = surface_gravity_g

    def calculate_atmospheric_conditions(self, distance_override=None):
        """
        Calculates the atmospheric conditions of the planet, including surface
        temperature and atmospheric pressure.
        """
        distance = float(distance_override) if distance_override is not None else float(self.distance)
        orbital_radius_km = distance * AU_TO_KM
        output_area = 4 * math.pi * orbital_radius_km ** 2
        solar_output_at_orbit = (self.star_output / output_area) / 1e6
        albedo = random.uniform(0.12, 0.35)
        surface_temperature_no_atmosphere = ((1 - albedo) * solar_output_at_orbit / (4 * STEFAN_BOLTZMANN_CONSTANT)) ** (1 / 4)

        if self.atmosphere == "None":
            self.surface_temperature = surface_temperature_no_atmosphere
            self.atmospheric_pressure = 0.0
        else:
            scale_height = (R * surface_temperature_no_atmosphere) / (self.atm_molar_density * self.gravity * EARTH_GRAVITY)
            self.scale_height = scale_height
            atmosphere_thickness = scale_height * random.uniform(5, 7)
            planet_volume = (4 * math.pi * self.radius ** 3) / 3
            atmosphere_volume = (4 * math.pi * (self.radius + atmosphere_thickness) ** 3) / 3 - planet_volume
            atmospheric_mass = 0
            for zone in range(round(random.uniform(5, 7))):
                zone_volume = atmosphere_volume + planet_volume - (4 * math.pi * (self.radius + (zone * scale_height)) ** 3) / 3
                zone_density = self.atm_density / (zone * 2.7) if zone >= 1 else self.atm_density
                atmospheric_mass += zone_volume * zone_density
            atmospheric_force = atmospheric_mass * (self.gravity * EARTH_GRAVITY)
            planet_surface_area = 4 * math.pi * (self.radius * 1000) ** 2
            atmospheric_pressure = (atmospheric_force / planet_surface_area) * 7500
            CO2_BASE_MOLAR_DENSITY = 0.04345
            CO2_MAX_GREENHOUSE_FACTOR = 5
            greenhouse_factor = abs((self.atm_molar_density - CO2_BASE_MOLAR_DENSITY) / CO2_BASE_MOLAR_DENSITY * CO2_MAX_GREENHOUSE_FACTOR)
            surface_temperature_atmosphere = ((1 - albedo) * solar_output_at_orbit * (1 + greenhouse_factor) / (4 * STEFAN_BOLTZMANN_CONSTANT)) ** (1 / 4)
            self.surface_temperature = surface_temperature_atmosphere
            self.atmospheric_pressure = atmospheric_pressure

            if self.planet_class == "M":
                if self.atmospheric_pressure < 90000 or self.atmospheric_pressure > 112000:
                    self.atmospheric_pressure = random.uniform(90000, 112000)
                if self.surface_temperature < 283 or self.surface_temperature > 290:
                    self.surface_temperature = random.uniform(283, 290)
            elif self.planet_class == "P" and self.surface_temperature >= 283:
                self.surface_temperature = random.uniform(surface_temperature_no_atmosphere, 283) if surface_temperature_no_atmosphere < 283 else random.uniform(-surface_temperature_no_atmosphere, 283)

    def generate_moons(self):
        """
        Generates a system of moons for the given planet.
        """
        moon_blacklist = ['Q', 'R', 'V', 'W', 'X', 'Y']
        max_moon_mass = self.mass / 10
        max_moon_radius = self.radius / (10 ** (1 / 3))
        possible_classes = [c for c, data in planet_mass_ranges.items()
                            if PLANET_CLASSES[c][self.zone] and PLANET_CLASSES[c]["type"] == 't' and c not in moon_blacklist
                            and data[1] <= max_moon_mass and PLANET_CLASSES[c]['radius_range'][1] <= max_moon_radius]
        if not possible_classes:
            return

        low_orbit = self.scale_height * 15 if self.scale_height else 100
        high_orbit = self.min_orbit_distance * AU_TO_KM
        total_orbit_distance = low_orbit

        while total_orbit_distance < high_orbit and total_orbit_distance < (self.distance * AU_TO_KM):
            classes = list(PLANET_CLASS_PROBABILITIES.keys())
            probabilities = list(PLANET_CLASS_PROBABILITIES.values())
            moon_class = random.choices(classes, weights=probabilities, k=1)[0]
            while moon_class not in possible_classes:
                moon_class = random.choices(classes, weights=probabilities, k=1)[0]
            
            radius_limit = PLANET_CLASSES[moon_class]['radius_range'][1] if max_moon_radius > PLANET_CLASSES[moon_class]['radius_range'][1] else max_moon_radius
            moon_distance = random.uniform(total_orbit_distance, high_orbit) / AU_TO_KM
            moon_radius = random.uniform(PLANET_CLASSES[moon_class]['radius_range'][0], radius_limit)

            new_moon = Planet(self.hab, moon_distance, self.star_output, self.star_radius,
                              self.star_temperature, self.star_mass, radius=moon_radius,
                              planet_class=moon_class, zone_override=self.zone,
                              distance_override=self.distance, is_moon=True)
            self.moons.append(new_moon)
            total_orbit_distance = (new_moon.distance + new_moon.min_orbit_distance) * AU_TO_KM

    def _validate_planet_class(self, zone):
        if self.planet_class not in PLANET_CLASSES or not PLANET_CLASSES[self.planet_class][zone]:
            raise ValueError("Invalid planet class for this zone")

    def _validate_radius(self):
        min_radius, max_radius = PLANET_CLASSES[self.planet_class]["radius_range"]
        if not (min_radius <= self.radius <= max_radius):
            raise ValueError("Invalid radius for planet class")

    def _validate_mass(self):
        min_mass, max_mass = planet_mass_ranges[self.planet_class]
        if not (min_mass <= self.mass <= max_mass):
            raise ValueError("Invalid mass for planet class")

    def __str__(self):
        """
        Returns the wiki template text for this object.
        """
        if self.is_moon:
            distance_text = f"|distance={to_scientific_notation(self.distance * AU_TO_KM, 5)} km"
            header_level = '==='
        else:
            distance_text = f"|distance={to_scientific_notation(self.distance * AU_TO_KM, 1)} km ({round(self.distance, 3)} AU)" if self.distance < 1 else f"|distance={round(self.distance, 3)} AU"
            header_level = '=='

        header = f"{header_level} {self.name} {header_level}"
        radius_string = f"|radius={round(self.radius, 2):,} km" if self.radius <= 100000 else f"|radius={to_scientific_notation(self.radius, 2)} km"

        output = [
            header, "{{Planet Data", f"|name={self.name}", f"|class={self.planet_class}",
            distance_text, f"|period={years_to_time_string(self.period)}",
            radius_string, f"|gravity={round(self.gravity, 3)} g", "}}",
        ]

        if not self.is_moon:
            if len(self.moons) > 1:
                output.append(f"There are {len(self.moons)} moons orbiting this planet.")
            elif len(self.moons) == 1:
                output.append("There is 1 moon orbiting this planet.")
            else:
                output.append("There are no moons orbiting this planet.")

        if self.type == "t":
            if self.atmosphere != "None":
                output.append(f"Surface conditions are an average of {self.atmospheric_pressure / 1000:.1f} kPa or {self.atmospheric_pressure / 101300:.2f} atmospheres and an average surface temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
                output.append(self.atmosphere)
            else:
                output.append(f"There is no atmosphere and the surface has an average temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
        else:
            output.append(f"The average internal conditions of this gas giant are an average of {self.atmospheric_pressure / 1000:.1f} kPa or {self.atmospheric_pressure / 101300:.1f} atmospheres and {self.surface_temperature - 273.15:.1f} degrees C.")
            output.append(f"The atmospheric pressure drops by a third to a half for every {self.scale_height / 1000:.1f} km from the core.")

        output.append(self.description)
        if output[-1] != self.composition:
            output.append(self.composition)

        if self.moons:
            for moon in self.moons:
                output.append('\n' + str(moon))

        return '\n'.join(output)