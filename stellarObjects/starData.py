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
import re
from .utils import to_scientific_notation, calculate_habitable_zone, calculate_stellar_radius, generate_phoneme_salad_name, calculate_hill_sphere, to_paragraph
from .constants import (STEFAN_BOLTZMANN_CONSTANT, SOLAR_MASS_TO_KG, SOLAR_LUMINOSITY,
                        SPECTRAL_LUMINOSITY_RANGES, YERKES_LUMINOSITY_RANGES, YERKES_MASS_CONSTRAINTS,
                        WHITE_DWARF_BASE_RADIUS_KM, CHANDRASEKHAR_LIMIT_SOL,
                        TEMP_RANGES, G, MILKY_WAY_MASS, GALACTIC_CENTER_DISTANCE_LY, LY_TO_M,
                        AU_TO_M, ISM_PRESSURE, SOLAR_RADIUS_M, SOLAR_ESCAPE_VELOCITY,
                        SOLAR_WIND_VELOCITY, SOLAR_MASS_LOSS_RATE, STAR_EVOLUTION)
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import config

class Star:
    """
    Represents a single star, encapsulating its physical and orbital properties.

    This class orchestrates the procedural generation of a star, from its
    spectral type and luminosity to its mass and radius. It serves as the
    gravitational and energetic center of a `StarSystem`, defining the context
    for all orbiting planets and other celestial bodies.

    The generation process can be random or guided by specific parameters, such
    as forcing a large star or specifying a precise star type (e.g., 'G2V').
    Once the core properties are established, the class calculates derived
    attributes, including the habitable zone, the system's gravitational
    perimeter, and the heliosphere's radius.
    """

    def get_star_age(self):
        """
        Calculates the star's age and lifespan based on its spectral and Yerkes class.

        This method determines the age and expected lifespan of the star by using
        pre-defined evolutionary data from the `STAR_EVOLUTION` constant. The
        lifespan is primarily based on the star's spectral class, but is then
        adjusted based on its Yerkes luminosity class to account for different
        evolutionary stages (e.g., giants, dwarfs).

        For most stars, the age is randomly determined within a range appropriate
        for its class. For example, giant stars are assumed to be near the end of
        their lives, while main-sequence stars can be of any age. White dwarfs
        are a special case, with a practically infinite lifespan representing their
        cooling period.

        Returns:
            tuple: A tuple containing the star's age and lifespan in billions of years.
                   The lifespan can be `float('inf')` for white dwarfs.
        """
        spectral_class = self.type[0]
        star_info = STAR_EVOLUTION.get(spectral_class, {})
        
        base_lifespan = star_info.get("lifespan_gy", 10) # Default to 10 GY for G-class stars
        
        # Adjust lifespan based on Yerkes class
        if self.yerkes_class in ["0", "IA+", "IA", "IAB", "IB", "II", "III"]: # Giants/Supergiants
            lifespan = base_lifespan * 0.1 # Significantly shorter lifespan
            age = random.uniform(lifespan * 0.8, lifespan) # Likely to be near the end of their life
        elif self.yerkes_class == "IV": # Subgiants
            lifespan = base_lifespan * 0.5
            age = random.uniform(lifespan * 0.8, lifespan)
        elif self.yerkes_class == "VI": # Subdwarfs
            lifespan = base_lifespan * 1.2 # Longer lifespan
            age = random.uniform(0.1, lifespan)
        elif self.yerkes_class in ["VII", "D"]: # White Dwarfs
            lifespan = float('inf') # Effectively infinite lifespan for cooling
            age = random.uniform(0.1, 10) # Age represents cooling time
        else: # Main Sequence (V)
            lifespan = base_lifespan
            age = random.uniform(0.1, lifespan)
            
        return age, lifespan

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

        Returns:
            float: The radius of the Hill sphere in Astronomical Units (AU).
        """
        galactic_center_dist_m = GALACTIC_CENTER_DISTANCE_LY * LY_TO_M
        hill_radius_m = calculate_hill_sphere(galactic_center_dist_m, self.mass, MILKY_WAY_MASS)
        return hill_radius_m / AU_TO_M

    def calculate_heliosphere(self):
        """
        Estimates the radius of the star's heliosphere (astrosphere).

        The heliosphere is a "bubble" in space created by the stellar wind, a
        stream of charged particles flowing from the star. The edge of this
        bubble, the heliopause, is where the stellar wind's pressure is balanced
        by the interstellar medium (ISM).

        The calculation is based on the star's mass-loss rate and wind velocity,
        which are determined by its spectral type and luminosity class. The
        method uses a tiered approach to model these properties, from the intense
        winds of hypergiants to the more gentle outflows of sun-like stars.

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
        elif yerkes_class in ["IA+", "IA", "IAB", "IB", "II"]:
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

        # TIER 4: Main Sequence (Class V), Subdwarfs (VI), and White Dwarfs (VII)
        # For sun-like stars and stellar remnants, mass loss is very low. We scale
        # directly from the Sun's known properties for a stable, physically grounded result.
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

    def __init__(self):
        """
        Initializes a Star object, generating its properties based on specified constraints.

        This constructor orchestrates the creation of a star. It can generate a
        star randomly or be guided by specific parameters. The `star_type`
        parameter allows for the creation of a star with a specific spectral
        class, subclass, and Yerkes classification (e.g., 'G2V').

        Once the core properties are established, it calculates derived
        properties such as the habitable zone, system perimeter (Hill sphere),
        and the heliosphere radius. These calculations are essential for the
        subsequent generation of planets and other celestial bodies in the
        `StarSystem` class.
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
        self.generate_star()
        self.age, self.lifespan = self.get_star_age()
        self.habitable_zone = calculate_habitable_zone(self.luminosity)
        self.system_perimeter = self.calculate_system_perimeter()
        self.heliosphere_radius = self.calculate_heliosphere()

    def __str__(self):
        """
        Returns a string representation of the star's data.

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
            mass_string = f"{to_scientific_notation(self.mass)} kg ({sol_mass_val:.1f}% of Sol)"
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
        
        sentences = []
        if self.lifespan == float('inf'):
            sentences.append(f"The star is approximately {self.age:.2f} billion years old and is now a white dwarf, which will cool for trillions of years.")
        else:
            sentences.append(f"The star is approximately {self.age:.2f} billion years old, with an expected lifespan of {self.lifespan:.2f} billion years.")
        
        output.append(to_paragraph(sentences))
        return '\n'.join(output)

    def generate_star(self):
        """
        Generates the physical properties of the star based on a hierarchical,
        physically-constrained procedural model.

        The generation process is fundamentally different for user-specified stars
        (via `config.STAR_TYPE`) versus randomly generated ones.

        For specified stars:
        1.  It parses the type (e.g., 'G2V') to get spectral class, subclass, and Yerkes class.
        2.  Temperature is calculated from the spectral class and subclass.
        3.  A valid luminosity range is found by taking the *intersection* of the
            allowed luminosities for the spectral and Yerkes classes.
        4.  A final luminosity is chosen from this valid range.
        5.  Mass is then determined based on physical constraints for that Yerkes class.
        6.  Radius is calculated last, using specific physics for the star type (e.g.,
            white dwarf mass-radius relation vs. Stefan-Boltzmann law for others).

        For random stars:
        1.  A spectral class is chosen based on galactic population probabilities.
            If `config.FORCE_HABITABLE_WORLD` is True, the choice is restricted
            to star types that can support complex life (A, F, G, K, M).
        2.  A luminosity is randomly chosen from that spectral class's typical range.
        3.  The Yerkes class is *determined* by this luminosity.
        4.  Temperature, mass, and radius then follow from these properties.

        This ensures that all generated stars, whether specified or random, adhere
        to the fundamental physical relationships between their core properties.
        """
        yerkes_lookup = {
            "0": "Hypergiant", "IA+": "Luminous Supergiant", "IA": "Supergiant",
            "IAB": "Intermediate-size Luminous Supergiant", "IB": "Less Luminous Supergiant",
            "II": "Bright Giant", "III": "Giant", "IV": "Subgiant", "V": "Main Sequence",
            "VI": "Subdwarf", "VII": "White Dwarf", "D": "White Dwarf"
        }

        if config.STAR_TYPE:
            # --- GENERATE STAR FROM SPECIFIED TYPE ---

            # 1. Parse and validate the specified star type string.
            match = re.match(r"([OBAFGKM])([0-9])(IA\+|IAB|VII|III|IA|IB|II|IV|VI|0|V|D)", config.STAR_TYPE.upper())
            if not match:
                raise ValueError("Invalid star type format. Expected format is e.g., G2V.")
            spectral_class, subclass_str, yerkes_class_str = match.groups()
            subclass = int(subclass_str)
            self.yerkes_class = yerkes_class_str
            yerkes_type = yerkes_lookup[yerkes_class_str]

            # 2. Calculate Temperature from spectral class and subclass.
            min_temp, max_temp = TEMP_RANGES[spectral_class]
            temp_range_size = max_temp - min_temp
            temperature = min_temp + (9 - subclass) * (temp_range_size / 9)
            temperature = int(round(temperature, -2))

            # 3. Determine a physically valid Luminosity.
            # Find the overlapping luminosity range between the spectral and Yerkes classes.
            spec_min_lum, spec_max_lum = SPECTRAL_LUMINOSITY_RANGES[spectral_class]
            yerkes_min_lum, yerkes_max_lum = YERKES_LUMINOSITY_RANGES[yerkes_class_str]
            
            # Special case for hot, young white dwarfs, which can be temporarily very luminous.
            if yerkes_class_str in ["VII", "D"] and spectral_class in ["O", "B"]:
                yerkes_min_lum, yerkes_max_lum = 0.1, 100

            valid_min_lum = max(spec_min_lum, yerkes_min_lum)
            valid_max_lum = min(spec_max_lum, yerkes_max_lum)
            
            if valid_min_lum > valid_max_lum:
                # If there's no overlap, the combination is physically unlikely.
                # We will prioritize the Yerkes class range, as it's the primary
                # driver of luminosity for evolved stars.
                valid_min_lum, valid_max_lum = yerkes_min_lum, yerkes_max_lum

            luminosity = random.uniform(valid_min_lum, valid_max_lum)

        else:
            # --- GENERATE STAR RANDOMLY ---

            # 1. Generate Spectral Class based on galactic population.
            if config.ABSURD:
                spectral_probabilities = {'O': 100, 'B': 0, 'A': 0, 'F': 0, 'G': 0, 'K': 0, 'M': 0}
            elif config.FORCE_LARGE_STAR:
                spectral_probabilities = {'O': 10, 'B': 20, 'A': 30, 'F': 30, 'G': 10, 'K': 0, 'M': 0}
            elif config.FORCE_HABITABLE_WORLD:
                # If forcing a habitable world, restrict to stars that can support "normal" or "slow" evolution.
                spectral_probabilities = {'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45}
            else:
                spectral_probabilities = {'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45}
            spectral_class = random.choices(list(spectral_probabilities.keys()), weights=spectral_probabilities.values(), k=1)[0]

            # 2. Generate Luminosity from the spectral class's typical range.
            min_luminosity, max_luminosity = SPECTRAL_LUMINOSITY_RANGES[spectral_class]
            luminosity = max_luminosity if config.ABSURD else random.uniform(min_luminosity, max_luminosity)

            # 3. Determine Yerkes Class from the resulting luminosity.
            if luminosity > YERKES_LUMINOSITY_RANGES["0"][0]:
                self.yerkes_class, yerkes_type = "0", "Hypergiant"
            elif luminosity > YERKES_LUMINOSITY_RANGES["IA"][0]:
                self.yerkes_class, yerkes_type = "IA", "Supergiant"
            elif luminosity > YERKES_LUMINOSITY_RANGES["IAB"][0]:
                self.yerkes_class, yerkes_type = "IAB", "Intermediate-size Luminous Supergiant"
            elif luminosity > YERKES_LUMINOSITY_RANGES["IB"][0]:
                self.yerkes_class, yerkes_type = "IB", "Less Luminous Supergiant"
            elif luminosity > YERKES_LUMINOSITY_RANGES["II"][0]:
                self.yerkes_class, yerkes_type = "II", "Bright Giant"
            elif luminosity > YERKES_LUMINOSITY_RANGES["III"][0]:
                self.yerkes_class, yerkes_type = "III", "Giant"
            elif luminosity > YERKES_LUMINOSITY_RANGES["IV"][0]:
                self.yerkes_class, yerkes_type = "IV", "Subgiant"
            elif luminosity > YERKES_LUMINOSITY_RANGES["VI"][0]: # Check Subdwarf before Main Sequence
                 self.yerkes_class, yerkes_type = "VI", "Subdwarf"
            elif luminosity > SPECTRAL_LUMINOSITY_RANGES["M"][0]: # Check against dimmest main sequence
                self.yerkes_class, yerkes_type = "V", "Main Sequence"
            else:
                self.yerkes_class, yerkes_type = "VII", "White Dwarf"

            # 4. Calculate Temperature and Subclass.
            min_temp, max_temp = TEMP_RANGES[spectral_class]
            temperature = int(round(random.uniform(min_temp, max_temp), -2))
            temp_range_size = max_temp - min_temp
            subclass = 9 - round((temperature - min_temp) / temp_range_size * 9)

        # --- CALCULATE FINAL PROPERTIES (COMMON TO BOTH PATHS) ---

        # 5. Calculate Mass based on Yerkes class constraints.
        min_mass, max_mass = YERKES_MASS_CONSTRAINTS[self.yerkes_class]
        
        if self.yerkes_class == "V":
            # For main-sequence stars, the mass-luminosity relation is strong.
            mass_sol = luminosity ** (1 / 3.5)
        elif self.yerkes_class in ["VII", "D"]:
            # For white dwarfs, mass is tightly constrained. Hotter (younger) ones
            # are typically more massive, closer to the Chandrasekhar limit.
            if spectral_class in ["O", "B"]:
                mass_sol = random.uniform(1.1, CHANDRASEKHAR_LIMIT_SOL)
            else:
                mass_sol = random.uniform(min_mass, 1.0)
        else:
            # For giants and supergiants, mass is less predictable from luminosity alone.
            # We choose a random mass within the physically allowed range for the class.
            mass_sol = random.uniform(min_mass, max_mass)
            
        # Ensure the calculated mass is within the absolute physical bounds for its class.
        mass_sol = max(min(mass_sol, max_mass), min_mass)
        mass = mass_sol * SOLAR_MASS_TO_KG

        # 6. Calculate Radius based on the star's type.
        if self.yerkes_class in ["VII", "D"]:
            # White dwarf radius follows an inverse mass-radius relationship.
            # R ∝ M^(-1/3). A 1 solar mass WD is the base.
            radius = WHITE_DWARF_BASE_RADIUS_KM * (mass_sol ** (-1/3))
        else:
            # For all other stars, radius is calculated from luminosity and temperature
            # using the Stefan-Boltzmann law.
            luminosity_watts = luminosity * SOLAR_LUMINOSITY
            radius = math.sqrt(luminosity_watts / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature ** 4)) / 1000

        # 7. Set final star properties.
        color_descriptions = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}
        star_type_str = f"{spectral_class}{subclass}{self.yerkes_class} {color_descriptions[spectral_class]} {yerkes_type} Star"

        self.type = star_type_str
        self.radius = radius
        self.mass = mass
        self.temperature = temperature
        self.luminosity = luminosity * SOLAR_LUMINOSITY