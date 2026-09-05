# stellarObjects/doubleStar.py

"""
Binary Star Proxy
=================

This module defines the `BinaryStarProxy` class, which acts as a unified
representation for a binary star system. It allows the `StarSystem` to treat
a binary pair as a single, effective central star for the purpose of planet
generation, particularly for P-type (circumbinary) orbits.

The proxy aggregates properties like mass, luminosity, and habitable zone from
its constituent stars, providing a consistent interface while encapsulating
the complexity of a multi-star system.
"""

import random

from .config import SystemConfig
from . import physical_constants, program_constants
from .starData import Star
from .utils import (format_age_string, calculate_habitable_zone, calculate_hill_sphere,
                    format_relative_to_sol, properties_to_string, to_scientific_notation)

class BinaryStarProxy(Star):
    """
    A proxy class representing a binary star system as a single effective star.

    This class encapsulates two `Star` objects (primary and secondary) and
    exposes their combined properties (mass, luminosity, habitable zone, etc.)
    to the `StarSystem` class. This allows the existing planet generation logic
    to function with minimal changes, treating the binary pair as a single
    gravitational and energetic source for circumbinary (P-type) planets.

    Attributes:
        _primary (Star): The primary star of the binary system.
        _secondary (Star): The secondary star of the binary system.
        stars (list): A list containing both primary and secondary Star objects.
        _effective_mass (float): The combined mass of both stars in kilograms.
        _effective_luminosity (float): The combined luminosity of both stars in Watts.
        _binary_separation_au (float): The orbital separation between the two stars in AU.
    """

    def __init__(self, system_config: SystemConfig, primary_star: Star, secondary_star: Star):
        """
        Initializes a BinaryStarProxy object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object for the system.
            primary_star (Star): The primary star instance.
            secondary_star (Star): The secondary star instance.
        """
        # Initialize the base Star class with _skip_property_init=True
        # The name will be overridden, and other properties will be handled by getters.
        super().__init__(system_config, name=f"{primary_star.name} Binary System", _skip_property_init=True)

        if secondary_star.mass > primary_star.mass:
            temp_star = secondary_star
            secondary_star = primary_star
            primary_star = temp_star

        self._primary = primary_star
        self._secondary = secondary_star

        self.age = max(self._primary.age, self._secondary.age)
        self.lifespan = max(self._primary.lifespan, self._secondary.lifespan)

        self.stars = [self._primary, self._secondary]

        # Calculate combined properties
        self._effective_mass = sum(s.mass for s in self.stars)
        self._effective_luminosity = sum(s.luminosity for s in self.stars)

        # Generate binary separation (0.05 to 0.25 AU for P-type systems)
        # Add the radii of both stars to the separation to ensure they don't overlap
        # and to account for their physical size in the orbital distance.
        base_separation = random.uniform(0.05, 0.25)
        # Convert star radii from kilometers to AU before adding to separation
        self._binary_separation_au = base_separation + \
                                     (self._primary.radius / physical_constants.AU_TO_KM) + \
                                     (self._secondary.radius / physical_constants.AU_TO_KM)

        # Override base Star properties with effective values
        # Changed this line to only use the primary star's name for the system name
        self.name = f"{primary_star.name} Binary System"
        self.type = f"Binary ({primary_star.type.split(' ')[0]}/{secondary_star.type.split(' ')[0]})" # Simplified type
        self.temperature = (primary_star.temperature + secondary_star.temperature) / 2 # Simple average
        self.radius = max(primary_star.radius, secondary_star.radius) # Use the larger radius for approximation
        # Note: self.lifespan was already set above to max(primary, secondary) lifespan;
        # do not overwrite it here with the primary-only value.

        # Recalculate derived properties based on effective values
        self.habitable_zone = calculate_habitable_zone(self._effective_luminosity)
        self.system_perimeter = self._calculate_system_perimeter_static(self._effective_mass)
        # For heliosphere, we need an effective radius and type for the static method.
        # This is a simplification, as binary heliospheres are complex.
        self.heliosphere_radius = Star._calculate_heliosphere_radius_static(
            self._effective_mass,
            self._effective_luminosity,
            self.radius, # Using the larger radius as an approximation
            self._primary.type, # Using primary's type for wind model selection
            self._primary.yerkes_class # Using primary's yerkes class for wind model selection
        )

    @property
    def mass(self):
        """Returns the combined mass of the binary system."""
        return self._effective_mass

    @property
    def luminosity(self):
        """Returns the combined luminosity of the binary system."""
        return self._effective_luminosity

    @property
    def binary_separation_au(self):
        """Returns the orbital separation between the two stars in AU."""
        return self._binary_separation_au

    def adjust_age_for_planets(self, planets):
        """
        Adjusts the age and lifespan of both constituent stars to accommodate
        the evolutionary requirements of the system's planets, then recomputes
        the proxy's combined `age` and `lifespan` as the maximum of the two
        (now-adjusted) stars' values.

        Args:
            planets (list): A list of `Planet` objects in the system.
        """
        for star_obj in self.stars:
            star_obj.adjust_age_for_planets(planets)
        # After adjusting individual stars, update the proxy's age to reflect the primary's
        self.age = max(self._primary.age, self._secondary.age)
        self.lifespan = max(self._primary.lifespan, self._secondary.lifespan)

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the binary star system's
        combined properties. It does NOT include individual star details.

        Returns:
            list: A two-element list: [combined data block, combined age sentence].
        """
        paragraphs = []

        # 1. Combined Binary System Data Block
        mass_string = format_relative_to_sol(self.system_config, self.mass, physical_constants.SOLAR_MASS_TO_KG, "kg")
        lum_string = format_relative_to_sol(self.system_config, self.luminosity, physical_constants.SOLAR_LUMINOSITY, "W", low_percent_precision=4)

        hab_lower = str(round(self.habitable_zone[0], program_constants.ROUND_HABITABLE_ZONE_AU))
        hab_upper = str(round(self.habitable_zone[1], program_constants.ROUND_HABITABLE_ZONE_AU))

        # Calculate separation in kilometers and format with scientific notation
        separation_km = self.binary_separation_au * physical_constants.AU_TO_KM
        separation_km_scientific = to_scientific_notation(self.system_config, separation_km)
        separation_string = f"{separation_km_scientific} km ({self.binary_separation_au:.2f} AU)"

        binary_properties = {
            "type": self.type,
            "mass": mass_string,
            "lum": lum_string,
            "hab": f"Between {hab_lower} and {hab_upper} AU",
            "separation": separation_string,
            "loc": f"{self._primary.name} & {self._secondary.name} Binary System" # Use full name for location
        }
        markdown_key_map = {
            "type": "Type", "mass": "Mass", "lum": "Luminosity",
            "hab": "Habitable Zone", "separation": "Stellar Separation", "loc": "Location"
        }
        paragraphs.append(properties_to_string(self.system_config, binary_properties, "Binary System Data", markdown_key_map=markdown_key_map))

        # 2. Combined Age and Evolutionary Notes (simplified for binary)
        age_sentence = f"The binary system is approximately {format_age_string(self.age)} old."
        paragraphs.append(age_sentence)

        return paragraphs

    @staticmethod
    def _calculate_system_perimeter_static(mass):
        """
        Calculates the system perimeter (Hill sphere relative to the galaxy)
        for an arbitrary mass, mirroring `Star.calculate_system_perimeter` but
        usable without a bound `Star` instance.

        Args:
            mass (float): The mass to use for the calculation, in kilograms
                          (typically the binary system's combined mass).

        Returns:
            float: The radius of the Hill sphere in Astronomical Units (AU).
        """
        galactic_center_dist_m = physical_constants.GALACTIC_CENTER_DISTANCE_LY * physical_constants.LY_TO_M
        hill_radius_m = calculate_hill_sphere(galactic_center_dist_m, mass, physical_constants.MILKY_WAY_MASS)
        return hill_radius_m / physical_constants.AU_TO_M