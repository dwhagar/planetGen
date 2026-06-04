# stellarObjects/starData.py

"""
Star Generation
===============

This module contains the `Star` class, which is used to generate and represent
a single star within a star system.
"""

import math
import random
from .utils import to_scientific_notation, calculate_habitable_zone, calculate_stellar_radius, generate_phoneme_salad_name
from .constants import (STEFAN_BOLTZMANN_CONSTANT, SOLAR_MASS_TO_KG, SOLAR_LUMINOSITY,
                        LUMINOSITY_RANGES, TEMP_RANGES)
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES

class Star:
    """
    A class representing a star and its properties.
    """

    def __init__(self, spectral_class=None, temperature=None, force_large=False, absurd=False):
        """
        Initializes a Star object with default values or user-provided values.

        Args:
            spectral_class (str, optional): The spectral class of the star. Defaults to None.
            temperature (int, optional): The temperature of the star in Kelvin. Defaults to None.
            force_large (bool, optional): Whether to force the generation of a large star. Defaults to False.
            absurd (bool, optional): Whether to generate an absurdly large star. Defaults to False.
        """
        self.name = generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.luminosity = None
        self.temperature = None
        self.mass = None
        self.radius = None
        self.type = None
        self.habitable_zone = None
        self.generate_star(force_large=force_large, spectral_class=spectral_class, temperature=temperature, absurd=absurd)
        self.habitable_zone = calculate_habitable_zone(self.luminosity)

    def __str__(self):
        """
        Returns a string representation of the star in the format of a wiki template.
        """
        if round(self.habitable_zone[0], 2) == round(self.habitable_zone[1], 2):
            hab_lower = str(round(self.habitable_zone[0], 5))
            hab_upper = str(round(self.habitable_zone[1], 5))
        else:
            if self.habitable_zone[0] < 0.01:
                hab_lower = str(round(self.habitable_zone[0], 5))
            else:
                hab_lower = str(round(self.habitable_zone[0], 2))

            if self.habitable_zone[1] < 0.01:
                hab_upper = str(round(self.habitable_zone[1], 5))
            else:
                hab_upper = str(round(self.habitable_zone[1], 2))

        sol_mass = round(self.mass / SOLAR_MASS_TO_KG, 1)
        sol_lum = round(self.luminosity / SOLAR_LUMINOSITY, 1)

        if sol_mass <= 0:
            sol_mass = round(self.mass / SOLAR_MASS_TO_KG * 100, 4)
        if sol_lum <= 0:
            sol_lum = round(self.luminosity / SOLAR_LUMINOSITY * 100, 4)

        if sol_mass <= 2:
            mass_string = f"|mass={to_scientific_notation(self.mass)} kg ({sol_mass * 100}% of Sol)"
        elif sol_mass > 100:
            exponent = int(math.floor(math.log10(abs(sol_mass))))
            mass_string = f"|mass={to_scientific_notation(self.mass)} kg (10<sup>{exponent}</sup>x Sol)"
        else:
            mass_string = f"|mass={to_scientific_notation(self.mass)} kg ({sol_mass}x Sol)"

        if sol_lum <= 2:
            lum_string = f"|lum={to_scientific_notation(self.luminosity)} W ({sol_lum * 100}% of Sol)"
        elif sol_lum > 100:
            exponent = int(math.floor(math.log10(abs(sol_lum))))
            lum_string = f"|lum={to_scientific_notation(self.luminosity)} W (10<sup>{exponent}</sup>x Sol)"
        else:
            lum_string = f"|lum={to_scientific_notation(self.luminosity)} W ({sol_lum}x Sol)"

        if self.radius <= 100000:
            radius_string = f"|radius={round(self.radius, 2):,} km"
        else:
            radius_string = f"|radius={to_scientific_notation(self.radius, 2)} km"

        output = ["{{Star Data", 
                  f"|name={self.name}",
                  f"|type={self.type}",
                  radius_string,
                  mass_string,
                  f"|temp={self.temperature} K",
                  lum_string,
                  f"|hab=Between {hab_lower} and {hab_upper} AU",
                  "}}"]
        return '\n'.join(output)

    def set_radius_bounds(self, luminosity, temperature, spectral_class):
        """
        Checks if a luminosity is reasonable for a given temperature and spectral class.

        Args:
            luminosity (float): The luminosity of the star in solar luminosities.
            temperature (float): The temperature of the star in Kelvin.
            spectral_class (str): The spectral class of the star ('O', 'B', 'A', 'F', 'G', 'K', or 'M').

        Returns:
            tuple: A tuple containing the minimum and maximum allowed radius in meters.
        """
        luminosity_watts = luminosity * SOLAR_LUMINOSITY
        radius_meters = calculate_stellar_radius(luminosity_watts, temperature)
        min_radius, max_radius = LUMINOSITY_RANGES[spectral_class]
        min_radius_meters = min_radius * 6.957e8
        max_radius_meters = max_radius * 6.957e8
        return min_radius_meters, max_radius_meters

    def generate_star(self, spectral_class=None, temperature=None, force_large=False, absurd=False):
        """
        Generates a random star's properties, optionally taking spectral class
        and temperature as input.
        """
        if absurd:
            spectral_probabilities = {'O': 100, 'B': 0, 'A': 0, 'F': 0, 'G': 0, 'K': 0, 'M': 0}
        elif force_large:
            spectral_probabilities = {'O': 10, 'B': 20, 'A': 30, 'F': 30, 'G': 10, 'K': 0, 'M': 0}
        else:
            spectral_probabilities = {'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45}

        if spectral_class is None:
            spectral_class = random.choices(list(spectral_probabilities.keys()), weights=spectral_probabilities.values(), k=1)[0]
        elif spectral_class not in spectral_probabilities:
            raise ValueError("Invalid spectral class")

        valid_temp_range = TEMP_RANGES.get(spectral_class)
        if valid_temp_range is None:
            raise ValueError("Invalid spectral class")

        if temperature is None:
            temperature = int(round(random.uniform(*valid_temp_range), -2)) if not absurd else valid_temp_range[1]
        elif not (valid_temp_range[0] <= temperature <= valid_temp_range[1]):
            raise ValueError("Temperature out of range for the given spectral class")

        min_luminosity, max_luminosity = LUMINOSITY_RANGES[spectral_class]
        luminosity = random.uniform(min_luminosity, max_luminosity) if not absurd else max_luminosity

        radius_min, radius_max = self.set_radius_bounds(luminosity, temperature, spectral_class)

        min_temp, max_temp = TEMP_RANGES[spectral_class]
        temp_range_size = max_temp - min_temp
        subclass = 9 - round((temperature - min_temp) / temp_range_size * 9)

        luminosity_watts = luminosity * SOLAR_LUMINOSITY
        radius = math.sqrt(luminosity_watts / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature ** 4)) / 1000
        if radius > radius_max:
            radius = radius_max
        elif radius < radius_min:
            radius = radius_min

        mass = (luminosity**(1/3.5) * SOLAR_MASS_TO_KG)

        if luminosity > 100000:  
            yerkes_class, yerkes_type = "0", "Hypergiant"
        elif luminosity > 30000: 
            yerkes_class, yerkes_type = "Ia+", "Luminous Supergiant"
        elif luminosity > 10000:       
            yerkes_class, yerkes_type = "Ia", "Supergiant"
        elif luminosity > 1000:
            yerkes_class, yerkes_type = "Ib", "Less Luminous Supergiant"
        elif luminosity > 25:        
            yerkes_class, yerkes_type = "II", "Bright Giant"
        elif luminosity > 5:           
            yerkes_class, yerkes_type = "III", "Giant"
        elif luminosity > 1.5:          
            yerkes_class, yerkes_type = "IV", "Subgiant"
        elif luminosity > 0.08:        
            yerkes_class, yerkes_type = "V", "Main Sequence"
        else:
            yerkes_class, yerkes_type = "D", "Dwarf"

        color_descriptions = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}
        star_type = f"{spectral_class}{subclass}{yerkes_class} {color_descriptions[spectral_class]} {yerkes_type} Star"

        self.type = star_type
        self.radius = radius
        self.mass = mass
        self.temperature = temperature
        self.luminosity = luminosity_watts