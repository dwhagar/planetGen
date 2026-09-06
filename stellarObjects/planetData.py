# stellarObjects/planetData.py

"""
Planet Generation
==================

This module defines the `Planet` class, the representation of a single
planet or moon within a star system. `Planet` itself holds the object's
state and its presentation (`to_paragraph_list`); the logic that determines
what a planet physically *is* lives in `planetPhysics`, and the logic that
determines whether/how life arises on it lives in `planetLife` — both
operate on `Planet` instances passed in as arguments, rather than as
methods on the class, to keep generation concerns separate from the object
itself.

A `Planet` is constructed *without* life data (see `planetPhysics`'s module
docstring); `StarSystem` applies `planetLife.apply_life_data` to every
planet and moon once the whole system has been generated.

The `AsteroidBelt` class, a simpler representation of an asteroid belt, is
defined separately in `asteroidData`.
"""

import math
import random
import secrets

from .config import SystemConfig
from .names import (MOON_NAMES, MOON_PREFIXES, MOON_SUFFIXES, PLANET_NAMES,
                    PLANET_PREFIXES, PLANET_SUFFIXES)
from . import physical_constants, planetPhysics, program_constants
from .utils import (format_length_km, generate_phoneme_salad_name, properties_to_string,
                    reseed_rng, to_paragraph, to_scientific_notation, years_to_time_string)


class Planet:
    """
    A class representing a single planet or moon and all of its properties.

    This class encapsulates the physical characteristics, orbital data, and
    atmospheric conditions of a celestial body. It can represent both
    planets and moons, and includes a list of any moons it has.

    The generation process is flexible, allowing for fully random creation or
    generation based on specified parameters like radius, mass, or class;
    this generation work is delegated to `planetPhysics.generate_planet_properties`
    and related functions in that module, called from `__init__`. Life
    chemistry and evolutionary data are NOT set during construction — see
    `planetLife.apply_life_data`, applied later by `StarSystem`.

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
        life_chemical (str): The primary chemical basis for any potential life
                             (set by `planetLife.apply_life_data`, None until then).
        hab (tuple): The inner and outer bounds of the habitable zone in AU.
        star_type (str): The spectral type of the star (e.g., 'G', 'M').
        star (Star): The Star object this planet orbits.
        evolutionary_data (list): A list of strings describing the evolutionary timeline
                                  (set by `planetLife.apply_life_data`, empty until then).
        flavor_text (str): A randomly selected flavor text for the planet.
        flavor_text_count (int): The number of flavor texts added to this planet.
    """

    def __init__(self, system_config: SystemConfig, star, hab_zone, distance, star_type, star_output, star_radius, star_temperature, star_mass,
                 radius=None, planet_class=None, mass=None, zone_override=None, distance_override=None,
                 is_moon=False, moon_count=None):
        """
        Initializes a Planet object with its physical and orbital properties.

        Life chemistry and evolutionary data are intentionally not set here;
        call `planetLife.apply_life_data(planet)` separately once the whole
        system's planets and moons have been generated (this is done
        automatically by `StarSystem.__init__`).

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
            moon_count (int, optional): An exact number of moons to generate
                                        for this planet (0 for none). Ignored
                                        if `is_moon` is True. If None, moon
                                        generation falls back to
                                        `system_config.MOONS`/random chance.
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
        self.evolutionary_data = [] # Populated later by planetLife.apply_life_data
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

        # Calculate physical/orbital properties (no life data yet).
        planetPhysics.generate_planet_properties(self, zone_override)
        self.volume = (4 / 3) * math.pi * self.radius ** 3  # Calculate volume in km^3
        self.period = math.sqrt(self.distance ** 3)
        planetPhysics.calculate_surface_gravity(self)
        planetPhysics.calculate_atmospheric_conditions(self, distance_override)

        if not self.is_moon:
            if moon_count is not None:
                planetPhysics.generate_moons(self, moon_count=moon_count)
            elif self.system_config.MOONS is not False and (self.system_config.MOONS is True or secrets.randbelow(2) == 1):
                planetPhysics.generate_moons(self)

    def _generate_life_and_flavor_paragraphs(self, object_type_desc, sentences):
        """
        Generates paragraphs related to life, evolution, and flavor text for the planet.
        It also ensures that recently used flavor texts are not repeated.

        Reads (but does not compute) `self.life_chemical`, `self.evolution`,
        and `self.evolutionary_data`, which are set by
        `planetLife.apply_life_data` before this is ever called (`to_paragraph_list`
        is only invoked once a `StarSystem` is fully generated).

        Args:
            object_type_desc (str): Description of the object type ("planet" or "moon").
            sentences (list): A list of sentences to append life-related descriptions to.

        Returns:
            list: A list of generated paragraphs.
        """
        reseed_rng()
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
        if self.evolutionary_data and self.planet_class in program_constants.HABITABLE_PLANET_CLASSES:
            life_paragraphs.extend(self.evolutionary_data)

        # Add flavor text if the random chance passes and limits are not exceeded
        if random.random() < program_constants.FLAVOR_CHANCE_PLANET and self.system_config.system_flavor_count < program_constants.MAX_FLAVOR_TOTAL:
            selected_flavor = None

            # Filter out recently used flavor texts
            available_habitable_flavor = [f for f in program_constants.HABITABLE_FLAVOR if f not in self.system_config.recent_flavor_texts]
            available_planet_flavor = [f for f in program_constants.PLANET_FLAVOR if f not in self.system_config.recent_flavor_texts]
            available_orbital_flavor = [f for f in program_constants.ORBITAL_FLAVOR if f not in self.system_config.recent_flavor_texts]

            # If all flavors of a category have been recently used, reset to the full list
            if not available_habitable_flavor:
                available_habitable_flavor = program_constants.HABITABLE_FLAVOR
            if not available_planet_flavor:
                available_planet_flavor = program_constants.PLANET_FLAVOR
            if not available_orbital_flavor:
                available_orbital_flavor = program_constants.ORBITAL_FLAVOR

            # Check for habitable and multicellular/technological life
            is_habitable = self.planet_class in program_constants.HABITABLE_PLANET_CLASSES
            has_multicellular_life = False
            if is_habitable and self.evolutionary_data:
                for stage_paragraph in self.evolutionary_data:
                    if "multicellularity" in stage_paragraph.lower() or "technological civilization" in stage_paragraph.lower():
                        has_multicellular_life = True
                        break

            if is_habitable and has_multicellular_life:
                selected_flavor = secrets.choice(available_habitable_flavor)
            elif self.type == "t" and self.planet_class != "A":
                selected_flavor = secrets.choice(available_planet_flavor)
            elif self.type == "g" or self.planet_class == "A":
                selected_flavor = secrets.choice(available_orbital_flavor)

            if selected_flavor:
                self.flavor_text = selected_flavor
                life_paragraphs.append(f"Sensors show {self.flavor_text}")
                self.system_config.system_flavor_count += 1
                self.flavor_text_count += 1

                # Add the selected flavor text to the recent list
                self.system_config.recent_flavor_texts.append(selected_flavor)
                # Keep the recent_flavor_texts list to a maximum size
                if len(self.system_config.recent_flavor_texts) > program_constants.MAX_RECENT_FLAVOR_TEXTS:
                    self.system_config.recent_flavor_texts.pop(0) # Remove the oldest entry

        return life_paragraphs

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs about the planet or moon,
        including its data block, atmospheric/surface conditions, life and
        flavor text, and (recursively) the paragraphs for any moons it has.

        Returns:
            list: A list of strings, where each string is a paragraph describing
                  the planet, or the moon if `self.is_moon` is set (moon-specific
                  wording is resolved via `object_type_desc` and the `is_moon`-aware
                  `self.description` set by `planetPhysics.generate_planet_properties`,
                  not by post-processing the rendered text here).
        """
        object_type_desc = "moon" if self.is_moon else "planet"

        if self.is_moon:
            # Moons orbit their parent planet, so their distance is from the planet, not the star.
            # This distance is typically much smaller and best represented in kilometers.
            distance_text = f"{to_scientific_notation(self.system_config, self.distance * physical_constants.AU_TO_KM, 4)} km" # Pass system_config
            header_level = '###' if self.system_config.MARKDOWN else '===' # Use self.system_config
        else:
            # For planets, the distance is from the star. We check if this distance
            # is large enough to warrant using light-years.
            distance_ly = self.distance * physical_constants.AU_TO_LY
            if distance_ly < program_constants.LY_THRESHOLD:
                # For "normal" sized systems, display in AU.
                # If less than 1 AU, also show in km for context.
                distance_text = f"{to_scientific_notation(self.system_config, self.distance * physical_constants.AU_TO_KM, 1)} km ({self.distance:.3f} AU)" if self.distance < 1 else f"{self.distance:.3f} AU" # Pass system_config
            else:
                # For very large systems, display in light-years.
                distance_text = f"{distance_ly:.4f} light-years"
            header_level = '##' if self.system_config.MARKDOWN else '==' # Use self.system_config

        output_paragraphs = []
        header = f"{header_level} {self.name} {header_level if not self.system_config.MARKDOWN else ''}".rstrip() # Use self.system_config
        output_paragraphs.append(header)

        radius_string = format_length_km(self.system_config, self.radius, 100000, 2) # Pass system_config

        planet_properties = {
            "class": self.planet_class,
            "distance": distance_text,
            "period": years_to_time_string(self.period),
            "radius": radius_string,
            "gravity": f"{round(self.gravity, 3)} g",
        }

        # Moons use a distinct wiki template name ("Class Data") instead of "Planet Data",
        # set directly here rather than via text substitution further down.
        template_name = "Class Data" if self.is_moon else "Planet Data"
        output_paragraphs.append(properties_to_string(self.system_config, planet_properties, template_name)) # Pass system_config

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

        return output_paragraphs

    def __str__(self):
        """
        Returns a string representation of the planet, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
