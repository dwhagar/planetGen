# stellarObjects/starData.py

"""
Star Generation and Properties
==============================

This module defines the `Star` class, the central component for creating and
managing the properties of a star within a procedurally generated star system.
The class handles the generation of fundamental stellar characteristics such as
mass, radius, temperature, and luminosity based on astrophysical models and
probabilistic distributions.

The generation process can be random or guided by specific parameters like
spectral class. Once the core properties are established, the class calculates
derived attributes, including the habitable zone, the gravitational extent of
the system (Hill sphere), and the size of the star's stellar wind bubble
(heliosphere).
"""

import math
import random
from .utils import to_scientific_notation, calculate_habitable_zone, calculate_stellar_radius, generate_phoneme_salad_name, calculate_hill_sphere, to_paragraph
from .constants import (STEFAN_BOLTZMANN_CONSTANT, SOLAR_MASS_TO_KG, SOLAR_LUMINOSITY, LUMINOSITY_RANGES, TEMP_RANGES, G,
                        MILKY_WAY_MASS, GALACTIC_CENTER_DISTANCE_LY, LY_TO_M, AU_TO_M, ISM_PRESSURE, SOLAR_RADIUS_M,
                        SOLAR_ESCAPE_VELOCITY, SOLAR_WIND_VELOCITY, SOLAR_MASS_LOSS_RATE)
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import config

class Star:
    """
    Represents a single star, encapsulating its physical and orbital properties.

    This class orchestrates the procedural generation of a star, from its
    spectral type and luminosity to its mass and radius. It serves as the
    gravitational and energetic center of a `StarSystem`, defining the context
    for all orbiting planets and other celestial bodies.

    A class representing a star and its properties.
    """

    def calculate_system_perimeter(self):
        """
        Calculates the Hill sphere for the star system relative to the galaxy.

        The Hill sphere is the region around a celestial body (in this case, the
        star) where its own gravity is the dominant force for attracting satellites.
        Objects orbiting within this sphere are gravitationally bound to the star,
        not the larger body it orbits (the galactic center).

        This method models the star's position relative to the Milky Way's center
        to determine the ultimate gravitational boundary of the star system. It
        serves as a practical "end" for the system, beyond which interstellar
        space truly begins. The result is converted to Astronomical Units (AU) for
        consistency with other system-scale measurements.

        Calculates the Hill sphere for the star system relative to the galaxy.
        """
        galactic_center_dist_m = GALACTIC_CENTER_DISTANCE_LY * LY_TO_M
        hill_radius_m = calculate_hill_sphere(galactic_center_dist_m, self.mass, MILKY_WAY_MASS)
        return hill_radius_m / AU_TO_M

    def calculate_heliosphere(self):
        """
        Estimates the radius of the star's heliosphere (astrosphere) based on its physical properties.

        The heliosphere is a vast "bubble" in space created by the stellar wind, a
        stream of charged particles flowing outward from the star. The edge of this
        bubble, called the heliopause, is the point where the pressure of the
        stellar wind is balanced by the pressure of the surrounding interstellar
        medium (ISM).

        The calculation involves several steps:
        1.  **Determine Stellar Wind Properties**: The model first calculates the
            star's mass-loss rate (how much mass it sheds via wind per year) and
            the wind's terminal velocity. These properties vary significantly
            between hot, massive stars (like O- and B-types) and cooler, less
            massive stars.
        2.  **Calculate Momentum Flux**: The momentum carried by the wind per second
            is calculated (mass_loss_rate * wind_velocity).
        3.  **Balance Pressures**: The heliopause radius is found by setting the
            dynamic pressure of the wind equal to the static pressure of the ISM.

        This provides a physically-grounded estimate for one of the key boundaries
        of the star system.
        Returns:
            float: The estimated radius of the heliosphere in AU.
        """
        # --- 1. Get Fundamental Stellar Properties ---
        radius_m = self.radius * 1000  # Convert radius from km to meters
        lum_sol = self.luminosity / SOLAR_LUMINOSITY
        mass_sol = self.mass / SOLAR_MASS_TO_KG
        radius_sol = radius_m / SOLAR_RADIUS_M

        # --- 2. Calculate Mass-Loss Rate (M-dot) and Wind Velocity (v_inf) ---
        # A single formula for mass loss is insufficient. We use a tiered system based on
        # the star's Yerkes luminosity class (evolutionary stage) and spectral type.

        escape_velocity = math.sqrt(2 * G * self.mass / radius_m)
        yerkes_class = self.yerkes_class

        # TIER 1: Hypergiants (Class 0)
        # These stars have extreme radiation-driven winds. We use a power-law that is
        # more stable than other models at this high-luminosity extreme.
        if yerkes_class == "0":
            # M-dot ~ L^1.5. This avoids the runaway effect of steeper power laws.
            mass_loss_rate_smyr = (10**-5.0) * (lum_sol**1.5)
            # Wind velocity is a high multiple of escape velocity.
            wind_velocity = 2.6 * escape_velocity

        # TIER 2: Supergiants (Class I) and Bright Giants (Class II)
        # For these highly evolved stars, the classic Reimers' Law provides a stable
        # and physically appropriate model for their powerful stellar winds.
        elif yerkes_class in ["Ia+", "Ia", "Ib", "II"]:
            eta = 2.0  # Higher efficiency factor for these very luminous stars.
            # Reimers' Law: M-dot = 4e-13 * η * (L*R/M)
            mass_loss_rate_smyr = (4e-13) * eta * (lum_sol * radius_sol / mass_sol)
            # Wind velocity is a smaller fraction of escape velocity for these cooler giants.
            wind_velocity = 0.3 * escape_velocity

        # TIER 3: Giants (Class III) and Subgiants (Class IV)
        # These are less luminous evolved stars. Reimers' Law is still the best model,
        # but with a standard efficiency factor.
        elif yerkes_class in ["III", "IV"]:
            eta = 1.0  # Standard eta for giants.
            mass_loss_rate_smyr = (4e-13) * eta * (lum_sol * radius_sol / mass_sol)
            wind_velocity = 0.3 * escape_velocity

        # TIER 4: Main Sequence (Class V) and Dwarfs (Class D)
        # For sun-like stars, mass loss is very low. We scale directly from the Sun's
        # known properties for a stable and physically grounded result.
        else:
            # For cool main-sequence stars, mass loss is driven by magnetic activity, not
            # radiation pressure. The Reimers' law scaling (L*R/M) is inaccurate here.
            # We use a more modern, empirically-derived formula from Wood et al. (2005),
            # which relates mass-loss rate to surface X-ray flux.
            # M-dot ~ R^2 * F_x^1.34, and F_x ~ L^-1.5
            # This results in a scaling relationship of M-dot ~ R^2 * L^-0.5
            # We normalize this to the Sun's known mass-loss rate.
            # The mass term is added to account for gravitational binding.
            scaling_factor = (radius_sol**2) * (lum_sol**-0.5) * (mass_sol**-1.0)
            mass_loss_rate_smyr = (2e-14) * scaling_factor # 2e-14 is the Sun's M-dot in M_sol/yr

            # Scale wind velocity based on the star's escape velocity relative to the Sun's.
            wind_velocity = SOLAR_WIND_VELOCITY * (escape_velocity / SOLAR_ESCAPE_VELOCITY)

        # --- 3. Convert and Calculate Final Radius ---

        # Convert M-dot from (solar masses/year) to (kg/s)
        mass_loss_rate_kgs = mass_loss_rate_smyr * SOLAR_MASS_TO_KG / (365.25 * 24 * 3600)

        # The heliopause radius is where the stellar wind's momentum flux balances the ISM pressure.
        # R = sqrt( (M-dot * v_inf) / (4 * pi * P_ism) )
        momentum_flux = mass_loss_rate_kgs * wind_velocity

        # Add a small floor value to prevent math domain errors for stars with near-zero wind.
        # This represents a minimal, residual pressure.
        momentum_flux = max(momentum_flux, 1e15)

        try:
            heliopause_radius_m = math.sqrt(momentum_flux / (4 * math.pi * ISM_PRESSURE))
        except ValueError:
            # This should not happen with the floor value, but as a final safety net.
            heliopause_radius_m = 0

        # Convert the final radius from meters to Astronomical Units (AU) for output.
        return heliopause_radius_m / AU_TO_M

    def __init__(self, spectral_class=None, temperature=None, force_large=False, absurd=False):
        """
        Initializes a Star object, generating its properties based on specified constraints.

        This constructor orchestrates the creation of a star. It can generate a
        star randomly or be guided by specific parameters like spectral class or
        temperature. Flags can be used to force the generation of particularly
        large or even physically absurd stars for creative purposes.

        Once the core properties (mass, radius, luminosity, temperature) are
        established, it calculates derived properties such as the habitable zone,
        system perimeter (Hill sphere), and the heliosphere radius.

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
        self.yerkes_class = None
        self.habitable_zone = None
        self.system_perimeter = None
        self.heliosphere_radius = None
        self.generate_star(force_large=force_large, spectral_class=spectral_class, temperature=temperature, absurd=absurd)
        self.habitable_zone = calculate_habitable_zone(self.luminosity)
        self.system_perimeter = self.calculate_system_perimeter()
        self.heliosphere_radius = self.calculate_heliosphere()

    def __str__(self):
        """
        Returns a string representation of the star in the format of a wiki template
        or a markdown table.

        The output format is determined by the `config.MARKDOWN` flag. This method
        formats the star's key properties—such as its name, type, radius, mass,
        temperature, luminosity, and habitable zone—into a clean, readable
        summary suitable for display in either a markdown table or a wikitext
        template.

        It includes logic to format very large or small numbers into scientific
        notation or percentages of solar values for improved readability.

        Returns:
            str: A formatted string representing the star's data.

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

        # Refactored mass and luminosity string generation for clarity and correctness
        sol_mass_val = self.mass / SOLAR_MASS_TO_KG
        if sol_mass_val < 0.01: # Less than 1%
            mass_string = f"{to_scientific_notation(self.mass)} kg ({sol_mass_val * 100:.2f}% of Sol)"
        elif sol_mass_val < 2: # Between 1% and 200%
            mass_string = f"{to_scientific_notation(self.mass)} kg ({sol_mass_val * 100:.1f}% of Sol)"
        else: # Greater than 200%
            mass_string = f"{to_scientific_notation(self.mass)} kg ({sol_mass_val:.1f}× Sol)"

        sol_lum_val = self.luminosity / SOLAR_LUMINOSITY
        if sol_lum_val < 0.01: # Less than 1%
            lum_string = f"{to_scientific_notation(self.luminosity)} W ({sol_lum_val * 100:.4f}% of Sol)"
        elif sol_lum_val < 2: # Between 1% and 200%
            lum_string = f"{to_scientific_notation(self.luminosity)} W ({sol_lum_val * 100:.1f}% of Sol)"
        else: # Greater than 200%
            lum_string = f"{to_scientific_notation(self.luminosity)} W ({sol_lum_val:.1f}× Sol)"

        if self.radius <= 100000:
            radius_string = f"{round(self.radius, 2):,} km"
        else:
            radius_string = f"{to_scientific_notation(self.radius, 2)} km"

        data = {
            "name": self.name,
            "type": self.type,
            "radius": radius_string,
            "mass": mass_string,
            "temp": f"{self.temperature} K",
            "lum": lum_string,
            "hab": f"Between {hab_lower} and {hab_upper} AU"
        }

        if config.MARKDOWN:
            output = ["# Star: " + data['name'],
                      "| Property | Value |",
                      "|---|---|",
                      f"| Type | {data['type']} |",
                      f"| Radius | {data['radius']} |",
                      f"| Mass | {data['mass']} |",
                      f"| Temperature | {data['temp']} |",
                      f"| Luminosity | {data['lum']} |",
                      f"| Habitable Zone | {data['hab']} |"]
        else:
            output = ["{{Star Data",
                      f"|name={data['name']}",
                      f"|type={data['type']}",
                      f"|radius={data['radius']}",
                      f"|mass={data['mass']}",
                      f"|temp={data['temp']}",
                      f"|lum={data['lum']}",
                      f"|hab={data['hab']}",
                      "}}"]
        
        return '\n'.join(output)

    def set_radius_bounds(self, luminosity, temperature, spectral_class):
        """
        Calculates the theoretical minimum and maximum radius for a star of a given class.

        This function is not currently used to constrain the radius but serves as a
        utility to understand the expected size range for a star based on its
        spectral class's luminosity range. It uses the Stefan-Boltzmann law to
        relate luminosity, temperature, and radius.

        Note: The current implementation calculates radius directly from luminosity
        and temperature, making this function more of a validation or reference tool.


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
        Generates the physical properties of the star based on procedural rules.

        This is the core generation method for the `Star` class. It follows these steps:
        1.  **Determine Spectral Class**: Selects a spectral class (O, B, A, etc.)
            based on probabilities that can be modified by flags like `force_large`.
        2.  **Determine Temperature**: Assigns a surface temperature within the
            range of the chosen spectral class.
        3.  **Determine Luminosity**: Assigns a luminosity within the range of the
            spectral class.
        4.  **Calculate Radius**: Computes the star's radius using the Stefan-Boltzmann
            law from its luminosity and temperature.
        5.  **Determine Luminosity Class**: Assigns a Yerkes classification (e.g.,
            V for Main Sequence, III for Giant) based on the luminosity.
        6.  **Calculate Mass**: Estimates the star's mass using an appropriate
            mass-luminosity relation for its class.
        7.  **Set Final Properties**: Assembles the final descriptive type string and
            assigns all calculated properties to the `Star` object.
        Generates a random star's properties, optionally taking spectral class
        and temperature as constraints.

        Args:
            spectral_class (str, optional): A specific spectral class ('O', 'B', etc.) to generate.
            temperature (int, optional): A specific temperature in Kelvin to generate.
            force_large (bool, optional): If True, biases generation towards larger, more massive stars.
            absurd (bool, optional): If True, generates an extremely large 'O' type star.
        """
        # Set spectral class probabilities based on generation flags
        if absurd:
            spectral_probabilities = {'O': 100, 'B': 0, 'A': 0, 'F': 0, 'G': 0, 'K': 0, 'M': 0}
        elif force_large:
            spectral_probabilities = {'O': 10, 'B': 20, 'A': 30, 'F': 30, 'G': 10, 'K': 0, 'M': 0}
        else:
            spectral_probabilities = {'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45}

        # Validate or generate spectral class
        if spectral_class is None:
            spectral_class = random.choices(list(spectral_probabilities.keys()), weights=spectral_probabilities.values(), k=1)[0]
        elif spectral_class not in spectral_probabilities:
            raise ValueError("Invalid spectral class")

        # Validate or generate temperature
        valid_temp_range = TEMP_RANGES.get(spectral_class)
        if valid_temp_range is None:
            raise ValueError("Invalid spectral class")

        if temperature is None:
            temperature = int(round(random.uniform(*valid_temp_range), -2)) if not absurd else valid_temp_range[1]
        elif not (valid_temp_range[0] <= temperature <= valid_temp_range[1]):
            raise ValueError("Temperature out of range for the given spectral class")

        # Generate luminosity
        min_luminosity, max_luminosity = LUMINOSITY_RANGES[spectral_class]
        luminosity = random.uniform(min_luminosity, max_luminosity) if not absurd else max_luminosity

        # Validate radius
        radius_min, radius_max = self.set_radius_bounds(luminosity, temperature, spectral_class)

        # Calculate spectral subclass
        min_temp, max_temp = TEMP_RANGES[spectral_class]
        temp_range_size = max_temp - min_temp
        subclass = 9 - round((temperature - min_temp) / temp_range_size * 9)

        # Calculate final radius and mass
        luminosity_watts = luminosity * SOLAR_LUMINOSITY
        radius = math.sqrt(luminosity_watts / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature ** 4)) / 1000
        if radius > radius_max:
            radius = radius_max
        elif radius < radius_min:
            radius = radius_min

        # Determine Yerkes spectral classification
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

        # Calculate mass based on a more robust mass-luminosity relation
        lum_sol = luminosity  # luminosity is already in solar units here
        if yerkes_class == "V":  # Main Sequence
            if lum_sol < 0.23:
                mass = (0.23 * (lum_sol ** 2) + 0.75 * lum_sol + 0.02) * SOLAR_MASS_TO_KG
            else:
                mass = (lum_sol ** (1 / 3.5)) * SOLAR_MASS_TO_KG
        elif yerkes_class in ["III", "IV"]:  # Giants and Subgiants
            # Giants have higher luminosity for a given mass than main sequence stars
            mass = (lum_sol ** (1 / 3.0)) * SOLAR_MASS_TO_KG
        else:  # Supergiants, Hypergiants
            # Mass-luminosity relation is less steep for very massive stars
            mass = (lum_sol ** (1 / 1.5)) * SOLAR_MASS_TO_KG

        # Set star type string
        color_descriptions = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}
        star_type = f"{spectral_class}{subclass}{yerkes_class} {color_descriptions[spectral_class]} {yerkes_type} Star"

        # Set final star properties
        self.type = star_type
        self.yerkes_class = yerkes_class
        self.radius = radius
        self.mass = mass
        self.temperature = temperature
        self.luminosity = luminosity * SOLAR_LUMINOSITY