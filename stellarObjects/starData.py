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
from .utils import to_scientific_notation, calculate_habitable_zone, calculate_stellar_radius, generate_phoneme_salad_name, calculate_hill_sphere, to_paragraph, properties_to_string
from .constants import (STEFAN_BOLTZMANN_CONSTANT, SOLAR_MASS_TO_KG, SOLAR_LUMINOSITY,
                        SPECTRAL_LUMINOSITY_RANGES, YERKES_LUMINOSITY_RANGES, YERKES_MASS_CONSTRAINTS,
                        WHITE_DWARF_BASE_RADIUS_KM, CHANDRASEKHAR_LIMIT_SOL,
                        TEMP_RANGES, G, MILKY_WAY_MASS, GALACTIC_CENTER_DISTANCE_LY, LY_TO_M,
                        AU_TO_M, ISM_PRESSURE, SOLAR_RADIUS_M, SOLAR_ESCAPE_VELOCITY,
                        SOLAR_WIND_VELOCITY, SOLAR_MASS_LOSS_RATE, STAR_EVOLUTION,
                        PLANET_CLASSES, EVOLUTIONARY_TIMELINES)
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
    Once the core properties are established, it calculates derived
    attributes, including the habitable zone, the system's gravitational
    perimeter, and the heliosphere's radius.
    """

    def _calculate_initial_star_age_and_lifespan(self):
        """
        Calculates the star's initial age and lifespan based on its spectral class
        and Yerkes class, using data from `STAR_EVOLUTION`.

        The lifespan is primarily determined by the star's spectral class, with
        adjustments made for its Yerkes luminosity class (e.g., giants have
        shorter lifespans than main-sequence stars of the same spectral type).
        For white dwarfs, the lifespan is considered effectively infinite as they
        slowly cool.

        The initial age is randomly chosen to be a fraction of its total lifespan,
        ensuring the star is in a plausible evolutionary stage. This age can be
        further adjusted later by `adjust_age_for_planets` to accommodate the
        evolutionary requirements of its planets.

        Returns:
            tuple: A tuple containing:
                - age (float): The initial age of the star in billions of years (GY).
                - lifespan (float): The total expected lifespan of the star in
                                    billions of years (GY), or `float('inf')` for white dwarfs.
        """
        spectral_class_char = self.type[0]
        star_info = STAR_EVOLUTION.get(spectral_class_char, {})

        min_lifespan, max_lifespan = star_info["lifespan_gy"]
        lifespan = random.uniform(min_lifespan, max_lifespan)

        min_age = 0.1
        max_age = lifespan * 0.9

        if config.AGE == "old":
            min_age = lifespan * 0.5
        elif config.AGE == "young":
            max_age = lifespan * (1/3)

        age = random.uniform(min_age, max_age) # Ensure age is less than lifespan

        return age, lifespan

    def adjust_age_for_planets(self, planets):
        """
        Adjusts the star's age to ensure it is at least as old as the minimum
        age required for the most "evolved" planet in its system. This method
        should be called by the `StarSystem` class after planets have been generated.

        The star's age is a crucial factor in planetary evolution. This function
        iterates through all generated planets, identifies their minimum required
        ages based on their `planet_class` and the star's `supported_evolutionary_scales`,
        and then ensures the star's age meets or exceeds the highest of these
        minimum requirements. This prevents scenarios where a young star might
        host an "old" planet (e.g., a terraformed world requiring billions of years
        of development).

        If the star's initial lifespan is too short to accommodate the required
        planetary age, the star's age is set to be near the end of its lifespan,
        reflecting a system that has evolved rapidly.

        Args:
            planets (list): A list of `Planet` objects in the system.
        """
        spectral_class_char = self.type[0]
        star_evolution_data = STAR_EVOLUTION.get(spectral_class_char, {})
        supported_scales = star_evolution_data.get("supported_evolutionary_scales", [])

        min_required_age_for_system = 0.0

        for planet in planets:
            # Assuming planet has a 'planet_class' attribute
            if hasattr(planet, 'planet_class') and planet.planet_class in PLANET_CLASSES:
                planet_class_data = PLANET_CLASSES[planet.planet_class]
                age_ranges = planet_class_data.get("age_ranges", {})

                for scale in supported_scales:
                    if scale in age_ranges:
                        min_planet_age_for_scale, _ = age_ranges[scale]
                        min_required_age_for_system = max(min_required_age_for_system, min_planet_age_for_scale)

        # Ensure the star's age is at least the minimum required by its planets
        if self.age < min_required_age_for_system:
            # If the star's lifespan is too short to support the required age,
            # we cap the age at a reasonable fraction of its lifespan.
            if self.lifespan != float('inf') and self.lifespan < min_required_age_for_system:
                # If lifespan is too short, set age to be near the end of its short life
                self.age = random.uniform(min_required_age_for_system * 0.8, self.lifespan * 0.95)
                if self.age > self.lifespan: # Final check
                    self.age = self.lifespan * 0.95
            else:
                # Otherwise, set age to be between the required minimum and near end of lifespan
                self.age = random.uniform(min_required_age_for_system, self.lifespan * 0.95 if self.lifespan != float('inf') else min_required_age_for_system + 5) # Add 5 GY for WD if no upper bound

        # Ensure age doesn't exceed lifespan (unless lifespan is infinite)
        if self.lifespan != float('inf') and self.age >= self.lifespan:
            self.age = self.lifespan * 0.95 # Star is near end of life

    def calculate_system_perimeter(self):
        """
        Calculates the Hill sphere for the star system relative to the galaxy.

        The Hill sphere defines the gravitational sphere of influence of the star
        within the larger galactic environment. Any object orbiting the star
        beyond this radius would be more strongly perturbed by the galaxy's
        gravity than by the star's, making it unlikely to be gravitationally
        bound to the star system.

        This method uses the star's mass, the estimated mass of the Milky Way,
        and the star's distance from the galactic center to compute this boundary.
        The result is converted from meters to Astronomical Units (AU) for
        convenience and consistency with other system-scale measurements.

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

        This function employs a tiered approach to model the stellar wind's
        mass-loss rate and velocity, which are highly dependent on the star's
        evolutionary stage (Yerkes luminosity class) and spectral type.
        Different physical models (e.g., power-laws for hypergiants, Reimers'
        Law for giants, and empirical scalings for main-sequence stars) are
        applied to ensure physically plausible results across the stellar spectrum.

        The heliopause radius is then calculated based on the balance between
        the stellar wind's momentum flux and the pressure of the interstellar medium.

        Returns:
            float: The estimated radius of the heliosphere in Astronomical Units (AU).
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
        momentum_flux = max(mass_loss_rate_kgs * wind_velocity, 1e-15) # Ensure momentum_flux is not zero or negative

        try:
            heliopause_radius_m = math.sqrt(momentum_flux / (4 * math.pi * ISM_PRESSURE))
        except ValueError:
            # This should not happen with the floor value, but as a final safety net.
            heliopause_radius_m = 0

        # Convert the final radius from meters to Astronomical Units (AU) for output.
        return heliopause_radius_m / AU_TO_M

    def __init__(self, name=None):
        """
        Initializes a Star object, generating its properties based on specified constraints.

        This constructor orchestrates the creation of a star. It can generate a
        star randomly or be guided by specific parameters provided through the
        `config` module. The `config.STAR_TYPE` parameter allows for the creation
        of a star with a specific spectral class, subclass, and Yerkes
        classification (e.g., 'G2V'). If `config.FORCE_LARGE_STAR` or `config.ABSURD`
        are set, these influence the star's initial generation to favor larger,
        more massive stars.

        If a `name` is provided, it will be used for the star; otherwise, a
        random name will be generated using `generate_phoneme_salad_name`.

        Once the core properties (`type`, `radius`, `mass`, `temperature`, `luminosity`)
        are established by `generate_star()`, this constructor calculates derived
        properties such as the `age`, `lifespan`, `habitable_zone`,
        `system_perimeter` (Hill sphere), and `heliosphere_radius`.

        Args:
            name (str, optional): The name to assign to the star. If None, a
                                  random name will be generated. Defaults to None.
        """
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
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
        self.age, self.lifespan = self._calculate_initial_star_age_and_lifespan()
        self.habitable_zone = calculate_habitable_zone(self.luminosity)
        self.system_perimeter = self.calculate_system_perimeter()
        self.heliosphere_radius = self.calculate_heliosphere()

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs about the star's properties.

        This method formats the star's physical characteristics (type, radius,
        mass, temperature, luminosity, habitable zone) into a structured block
        using `properties_to_string`. It also constructs a sentence describing
        the star's age and expected lifespan, incorporating evolutionary notes
        from `STAR_EVOLUTION` if available.

        The output is designed to be human-readable and can be formatted either
        as Wikitext or Markdown based on the `config.MARKDOWN` flag.

        Returns:
            list: A list of strings, where each string represents a paragraph
                  or a formatted data block describing the star.
        """
        paragraphs = [] # This will contain the final list of "paragraphs" for the star.

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

        star_properties = {
            "type": self.type,
            "radius": radius_string,
            "mass": mass_string,
            "temp": f"{self.temperature} K",
            "lum": lum_string,
            "hab": f"Between {hab_lower} and {hab_upper} AU",
            "loc": self.name # Adding the star's name as location
        }

        markdown_key_map = {
            "type": "Type",
            "radius": "Radius",
            "mass": "Mass",
            "temp": "Temperature",
            "lum": "Luminosity",
            "hab": "Habitable Zone",
            "loc": "Location"
        }
        
        star_block = properties_to_string(star_properties, "Star Data", markdown_key_map=markdown_key_map)
        paragraphs.append(star_block)

        # Construct the age and evolutionary notes sentence
        age_sentence_base = ""
        if self.lifespan == float('inf'):
            age_sentence_base = f"The star is approximately {self.age:.2f} billion years old and is now a white dwarf, which will cool for trillions of years"
        else:
            age_sentence_base = f"The star is approximately {self.age:.2f} billion years old, with an expected lifespan of {self.lifespan:.2f} billion years"

        spectral_class_char = self.type[0]
        star_info = STAR_EVOLUTION.get(spectral_class_char, {})
        
        full_age_and_notes_sentence = age_sentence_base
        if "evolutionary_constraint_notes" in star_info:
            # Append notes directly as they are a continuation, ensuring a space and period at the end
            full_age_and_notes_sentence += ". " + star_info["evolutionary_constraint_notes"]
        full_age_and_notes_sentence += "." # Ensure the entire combined sentence ends with a period

        paragraphs.append(full_age_and_notes_sentence)

        return paragraphs

    def __str__(self):
        """
        Returns a complete string representation of the star, suitable for display.

        This method compiles all the descriptive paragraphs generated by
        `to_paragraph_list()` and joins them with double newlines to create
        a well-formatted, multi-paragraph string.

        Returns:
            str: A formatted string describing the star and its properties.
        """
        return '\n\n'.join(self.to_paragraph_list())

    def generate_star(self):
        """
        Generates the physical properties of the star based on a hierarchical,
        physically-constrained procedural model.

        The generation process is fundamentally different for user-specified stars
        (via `config.STAR_TYPE`) versus randomly generated ones.

        For specified stars (when `config.STAR_TYPE` is set):
        1.  **Parsing Input**: The method first parses the `config.STAR_TYPE` string
            (e.g., 'G2V') to extract the spectral class (O, B, A, F, G, K, M),
            subclass (0-9), and Yerkes luminosity class (e.g., V for Main Sequence).
            It raises a `ValueError` if the format is invalid.
        2.  **Temperature Calculation**: The star's surface temperature is calculated
            based on its spectral class and subclass, interpolating within the
            `TEMP_RANGES` data.
        3.  **Luminosity Determination**: A valid luminosity range is established
            by finding the *intersection* of the allowed luminosity ranges for
            both the spectral class (`SPECTRAL_LUMINOSITY_RANGES`) and the
            Yerkes class (`YERKES_LUMINOSITY_RANGES`). A random luminosity is
            then chosen from this overlapping range. Special handling is included
            for white dwarfs, which can have temporarily high luminosities.
        4.  **Mass Calculation**: The star's mass is determined based on the
            physical constraints defined for its Yerkes class in `YERKES_MASS_CONSTRAINTS`.
            For main-sequence stars, a mass-luminosity relation is used. For white
            dwarfs, an inverse mass-radius relation is considered, with mass
            constrained by the Chandrasekhar limit. For giants and supergiants,
            mass is randomly selected within their class-specific bounds.
        5.  **Radius Calculation**: The star's radius is calculated last. For white
            dwarfs, an inverse mass-radius relationship is applied. For all other
            star types, the Stefan-Boltzmann law is used, deriving radius from
            luminosity and temperature.

        For random stars (when `config.STAR_TYPE` is None):
        1.  **Spectral Class Selection**: A spectral class is chosen probabilistically
            based on their prevalence in the galactic population. `config.ABSURD`
            or `config.FORCE_LARGE_STAR` can bias this selection towards hotter,
            more massive stars.
        2.  **Luminosity Generation**: A luminosity is randomly chosen from the
            typical range for the selected spectral class. `config.ABSURD` forces
            the maximum luminosity for the chosen class.
        3.  **Yerkes Class Determination**: The Yerkes luminosity class (e.g.,
            Supergiant, Main Sequence, White Dwarf) is *derived* from the
            generated luminosity, using thresholds defined in `YERKES_LUMINOSITY_RANGES`.
        4.  **Temperature and Subclass Calculation**: The star's surface temperature
            is randomly selected within the `TEMP_RANGES` for its spectral class,
            and the subclass is then derived from this temperature.
        5.  **Mass and Radius Calculation**: Similar to the specified star path,
            mass is determined based on Yerkes class constraints and luminosity,
            and radius is calculated using the appropriate physical laws.

        This comprehensive approach ensures that all generated stars, whether
        user-specified or randomly created, adhere to the fundamental physical
        relationships between their core properties, resulting in astrophysically
        plausible stellar objects.
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