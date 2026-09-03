# stellarObjects/systemData.py

"""
Star System Generation
======================

This module contains the `StarSystem` class, which is used to generate and
represent a full star system, including a central star and a list of planets.
The `StarSystem` class orchestrates the creation of the star and its planets,
applying a set of procedural generation rules to create a diverse and
scientifically grounded star system.

The generation process can be customized through a variety of parameters,
allowing for the creation of specific types of systems, such as those with
habitable worlds, asteroid belts, or a particular star type. The class also
includes methods for validating the system's orbital mechanics, counting the
number of celestial objects, and generating a detailed string representation of
the system.
"""

import random
import math
import copy # Added import for deepcopy
from .planetData import Planet
from . import planetPhysics
from . import planetLife
from .doubleStar import BinaryStarProxy # New import for binary star systems
from .asteroidData import AsteroidBelt
from .starData import Star
from stellarObjects import constants
from .utils import to_paragraph
from .config import SystemConfig # Updated import

class StarSystem:
    """
    A class representing a star system, containing a central star and a list of planets.

    The `StarSystem` class is the main entry point for generating a complete star
    system. It initializes a central star and then procedurally generates a list
    of planets and other celestial objects orbiting it. The generation process
    can be influenced by a variety of flags, allowing for fine-tuned control over
    the final system's characteristics.

    The class handles the complex logic of placing planets in stable orbits,
    ensuring that the system is physically plausible. It also includes methods
    for counting the number of different types of objects in the system and
    generating a detailed, human-readable summary of the system's properties.

    Attributes:
        star (Star): The central star of the system.
        planets (list): A list of `Planet` and `AsteroidBelt` objects orbiting the star.
        planet_count (int): The total number of planets in the system.
        belt_count (int): The total number of asteroid belts in the system.
        moon_count (int): The total number of moons in the system.
        hab_count (int): The total number of potentially habitable worlds.
        m_count (int): The total number of Class M worlds.
        system_flavor_count (int): Tracks the total number of flavor texts added across the system.
    """

    def __init__(self, system_config: SystemConfig): # Updated signature to accept system_config
        """
        Initializes a StarSystem object, generating a star and its planets.

        This constructor orchestrates the entire star system generation process.
        It begins by creating a `Star` instance, which can be customized with
        flags to control its size and type. It then estimates the number of
        celestial objects the system can support and proceeds to generate them,
        placing them in orbit around the star.

        The generation logic can be influenced by several parameters, allowing for
        fine-tuned control over the final system's characteristics. For example,
        the `force_hab` flag ensures that at least one habitable planet is generated,
        while the `no_planets` flag can be used to create a star with no orbiting bodies.
        The placement of planets is done sequentially, with each new planet's
        orbit being determined based on the position of the previous one to ensure
        a degree of realism in orbital spacing.

        If `config.NAME` is provided, it will be used as the name for the star
        system, overriding the default random name generation.

        Every `Planet` (and moon) is generated without life data — see
        `planetPhysics`'s module docstring. Once all planets/moons exist,
        orbits are validated, and the star's age has been finalized via
        `Star.adjust_age_for_planets`, this constructor makes one pass over
        every planet and moon and applies `planetLife.apply_life_data` to
        each, so evolutionary timelines reflect the star's final age.
        """
        self.system_config = system_config # Assign the passed SystemConfig instance
        self.star = Star(self.system_config, name=self.system_config.NAME) # Pass system_config and use its NAME
        self.primary_star = self.star # For single star systems, the primary is the star
        self.planets = []
        self.stars = [self.primary_star] # Keep track of individual stars

        if self.system_config.IS_BINARY_SYSTEM:
            # Create a copy of the system_config for the secondary star
            secondary_star_config = copy.deepcopy(self.system_config)
            # Ensure FORCE_LARGE_STAR is False for the secondary star
            secondary_star_config.FORCE_LARGE_STAR = False

            # Logic to generate a secondary star (e.g., random mass relative to primary)
            # This is new generation logic, but contained.
            # For simplicity, secondary mass is a fraction of primary mass.
            secondary_mass_factor = random.uniform(0.1, 0.8)
            secondary_mass = self.primary_star.mass * secondary_mass_factor
            # Create secondary star, potentially with a different name or type if desired
            self.secondary_star = Star(secondary_star_config, name=f"{self.primary_star.name} B", mass_override=secondary_mass)
            self.stars.append(self.secondary_star)
            self.star = BinaryStarProxy(self.system_config, self.primary_star, self.secondary_star) # self.star now points to the proxy

        # Removed: self.system_flavor_count = 0 # Initialize system flavor count
        system_objects = self.estimate_num_objects()
        star_factor = self.star.mass / constants.SOLAR_MASS_TO_KG

        required_objects = 0
        if self.system_config.FORCE_HABITABLE_WORLD: # Use self.system_config
            required_objects += 1
        if self.system_config.FORCE_ASTEROID_BELT: # Use self.system_config
            required_objects += 1

        if system_objects < required_objects:
            system_objects = required_objects

        if system_objects > 0:
            belt_index = random.randint(0, system_objects - 1) if self.system_config.FORCE_ASTEROID_BELT else -1 # Use self.system_config
            found_hab = False
            i = -1

            while i < system_objects - 1:
                i += 1
                last_asteroid = False
                if i > 0:
                    last_planet = self.planets[i - 1]
                    random_buffer = random.uniform(0, star_factor)
                    if last_planet.type == 'a':
                        estimated_distance = last_planet.upper_limit + random_buffer * 2
                        last_asteroid = True
                    else:
                        estimated_distance = (last_planet.distance + last_planet.min_orbit_distance) + random_buffer
                else:
                    estimated_distance = constants.INITIAL_PLANET_DISTANCE_FACTOR * star_factor

                hz = self.star.habitable_zone[0] < estimated_distance < self.star.habitable_zone[1]

                if self.system_config.FORCE_HABITABLE_WORLD and not found_hab: # Use self.system_config
                    if not hz and i == 0:
                        if (estimated_distance > self.star.habitable_zone[1] or
                                0 < self.star.habitable_zone[0] - estimated_distance < 0.2 or system_objects == 1):
                            estimated_distance = random.uniform(self.star.habitable_zone[0],
                                                                self.star.habitable_zone[1])
                            hz = True
                    elif not hz and i > 0:
                        last_planet = self.planets[i - 1]
                        beyond_hz = last_planet.type == 'a' and last_planet.upper_limit > self.star.habitable_zone[1] or \
                                    last_planet.type != 'a' and (last_planet.distance + last_planet.min_orbit_distance >
                                                                 self.star.habitable_zone[1])

                        if beyond_hz:
                            estimated_distance = random.uniform(self.star.habitable_zone[0],
                                                                self.star.habitable_zone[1])
                            planet = Planet(self.system_config, self.star, self.star.habitable_zone, estimated_distance, self.star.type[0], # Pass system_config
                                            self.star.luminosity,
                                            self.star.radius, self.star.temperature, self.star.mass,
                                            planet_class="M")
                            self.planets[i - 1] = planet
                            i -= 1
                            found_hab = True
                            continue
                        elif i == system_objects - 1:
                            estimated_distance = random.uniform(self.star.habitable_zone[0],
                                                                self.star.habitable_zone[1])
                            hz = True

                    if hz:
                        planet = Planet(self.system_config, self.star, self.star.habitable_zone, estimated_distance, self.star.type[0], # Pass system_config
                                        self.star.luminosity,
                                        self.star.radius, self.star.temperature, self.star.mass,
                                        planet_class="M")
                        found_hab = True
                        self.planets.append(planet)
                        continue

                if (random.random() < constants.ASTEROID_BELT_PROBABILITY or i == belt_index) and not last_asteroid and not hz:
                    min_distance = estimated_distance
                    max_distance = estimated_distance * random.uniform(constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN, constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)
                    self.planets.append(AsteroidBelt(self.system_config, estimated_distance, min_distance, max_distance)) # Pass system_config
                else:
                    planet = Planet(self.system_config, self.star, self.star.habitable_zone, estimated_distance, self.star.type[0], # Pass system_config
                                    self.star.luminosity,
                                    self.star.radius, self.star.temperature, self.star.mass)
                    if planet.planet_class == "M":
                        found_hab = True
                    self.planets.append(planet)

        self.validate_system()
        self.star.adjust_age_for_planets(self.planets)

        # All planets and moons are generated without life data (see planetPhysics's
        # module docstring). Apply it now, in one pass over the finished system, so
        # evolutionary timelines are computed against the star's final, planet-adjusted
        # age rather than its provisional pre-adjustment one.
        for obj in self.planets:
            if obj.type == 'a': # Skip asteroid belts; they carry no life data.
                continue
            planetLife.apply_life_data(obj)
            for moon in obj.moons:
                planetLife.apply_life_data(moon)

        self.planet_count, self.belt_count, self.moon_count = self.count_objects()
        self.hab_count, self.m_count = self.count_habitable()

    def count_objects(self):
        """
        Counts the number of planets, asteroid belts, and moons in the system.

        This method iterates through the list of celestial objects in the system
        and tallies the number of each type. It distinguishes between planets and
        asteroid belts, and also counts the total number of moons orbiting the
        planets. This information is used to provide a summary of the system's
        composition in the `__str__` method.

        Returns:
            tuple: A tuple containing the total number of planets, asteroid belts,
                   and moons in the system, in that order.
        """
        planet_counter, moon_counter, belt_counter = 0, 0, 0
        for planet in self.planets:
            if planet.type == 'a':
                belt_counter += 1
            else:
                planet_counter += 1
                moon_counter += len(planet.moons)
        return planet_counter, belt_counter, moon_counter

    def count_habitable(self):
        """
        Counts the number of potentially habitable worlds in the system.

        This method iterates through all planets and their moons, checking their
        classification to determine if they are potentially habitable. It counts
        worlds classified as Class H, K, L, M, O, or P, and also keeps a separate
        count of Class M worlds, which are considered the most Earth-like. This
        is used for the system summary output.

        Returns:
            tuple: A tuple containing the total number of potentially habitable
                   worlds and the total number of Class M worlds, in that order.
        """
        habitable_classes = "HKLMOP"
        hab_count, m_count = 0, 0
        for planet in self.planets:
            if planet.type != 'a':
                if planet.planet_class in habitable_classes:
                    hab_count += 1
                if planet.planet_class == "M":
                    m_count += 1
                for moon in planet.moons:
                    if moon.planet_class in habitable_classes:
                        hab_count += 1
                    if moon.planet_class == "M":
                        m_count += 1
        return hab_count, m_count

    def estimate_num_objects(self):
        """
        Estimates the number of objects in a star system based on the star's mass.

        This method calculates the potential number of celestial objects a star
        can support based on its mass. The number of objects scales with the
        star's mass, with more massive stars being able to support more objects.
        The calculation is tempered by a logarithmic function to prevent an
        excessive number of objects for very massive stars. The final number can
        be a random value up to the calculated maximum, or the maximum itself if
        `FORCE_MAX_PLANETS` is enabled.

        Returns:
            int: The estimated number of objects to be generated in the system.
                 Returns 0 if `NO_PLANETS` is enabled in the configuration.
        """
        if self.system_config.NO_PLANETS: # Use self.system_config
            return 0
        solar_masses = self.star.mass / constants.SOLAR_MASS_TO_KG

        # This provides a continuous scaling factor based on mass.
        # For a 1 solar mass star, this is 1.
        # For smaller stars, it's < 1; for larger stars, it's > 1.
        # The logarithm helps temper the explosive growth for very massive stars.
        scaling_factor = 1 + math.log10(solar_masses) if solar_masses >= 1 else solar_masses

        # Base number of objects for a 1 solar mass star is 15.
        max_objects = constants.BASE_MAX_SYSTEM_OBJECTS * scaling_factor
        if max_objects > constants.ABSOLUTE_MAX_SYSTEM_OBJECTS:
            max_objects = constants.ABSOLUTE_MAX_SYSTEM_OBJECTS
        return math.ceil(max_objects) if self.system_config.FORCE_MAX_PLANETS else random.randint(0, math.ceil(max_objects)) # Use self.system_config

    def validate_system(self):
        """
        Validates and adjusts the distances of stellar objects to prevent orbital overlap.

        This method iterates through the generated planets and other celestial
        objects, checking for any orbital overlaps. If two objects are too close
        to each other, it adjusts their orbits to ensure a safe distance is
        maintained. This is crucial for creating a physically plausible and stable
        star system.

        The method accounts for the different types of objects, such as planets
        and asteroid belts, and applies appropriate corrections to ensure that
        their orbits are not just non-overlapping, but also realistically spaced.
        If an adjustment is made, the planet's atmospheric conditions are
        recalculated to reflect its new orbital distance.
        """
        if len(self.planets) < 2:
            return

        for i in range(1, len(self.planets)):
            planet = self.planets[i]
            last_planet = self.planets[i - 1]

            if last_planet.type == 'a':
                distance_to_last = planet.distance - last_planet.upper_limit
            else:
                distance_to_last = planet.distance - last_planet.distance

            if distance_to_last < 0:
                additional_correction = abs(distance_to_last) + last_planet.distance
            else:
                additional_correction = 0

            if planet.type == 'a':
                if last_planet.type == 'a':
                    if distance_to_last < constants.MIN_ASTEROID_BELT_SEPARATION:
                        planet.distance += constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        planet.upper_limit += constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        planet.lower_limit += constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                elif distance_to_last < last_planet.min_orbit_distance:
                    planet.distance += last_planet.min_orbit_distance + additional_correction
                    planet.upper_limit += last_planet.min_orbit_distance + additional_correction
                    planet.lower_limit += last_planet.min_orbit_distance + additional_correction
            else:
                if last_planet.type == 'a':
                    if distance_to_last < constants.MIN_ASTEROID_BELT_SEPARATION:
                        planet.distance += constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        planetPhysics.calculate_atmospheric_conditions(planet)
                else:
                    min_orbit = max(planet.min_orbit_distance, last_planet.min_orbit_distance)
                    if distance_to_last < min_orbit:
                        planet.distance += min_orbit + additional_correction
                        planetPhysics.calculate_atmospheric_conditions(planet)

    def __str__(self):
        """
        Generates a string output for the system data, including a summary and
        details for each celestial body.

        This method compiles a comprehensive summary of the entire star system,
        including details about the central star and each of its orbiting
        objects. The output is formatted as a human-readable string, suitable for
        display in a console or for writing to a file.

        The method provides a high-level overview of the system, including the
        total number of planets, asteroid belts, and moons, as well as the number
        of potentially habitable worlds. It then lists each celestial body in
        order, providing a detailed description of its properties. The final
        output may also include a category tag for wiki-based systems.

        Returns:
            str: A formatted string representing the entire star system.
        """
        all_output_parts = []

        # Add level 1 header for the system name
        if self.system_config.MARKDOWN:
            all_output_parts.append(f"# {self.star.name}\n\n")
        else:
            all_output_parts.append(f"= {self.star.name} =\n\n")

        # Generate system summary sentences (needed before star details for binary systems)
        system_summary_sentences = []
        segments = []
        if 0 < self.planet_count != self.m_count and self.hab_count != self.planet_count:
            segments.append(f"{self.planet_count} planet{'s' if self.planet_count > 1 else ''}")
        if self.belt_count > 0:
            segments.append(f"{self.belt_count} asteroid belt{'s' if self.belt_count > 1 else ''}")
        if self.moon_count > 0:
            segments.append(f"{self.moon_count} moon{'s' if self.moon_count > 1 else ''}")

        if not segments:
            system_summary_sentences.append("There are no stellar objects in this system.")
        else:
            system_string = "This system contains " + ", ".join(segments)
            if len(segments) == 2:
                system_string = system_string.replace(f", {segments[-1]}", f" and {segments[-1]}")
            elif len(segments) > 2:
                system_string = system_string.replace(f", {segments[-1]}", f", and {segments[-1]}")
            system_string += "."
            system_summary_sentences.append(system_string)

        if self.m_count == 1:
            m_string = "1 of which is class M"
        elif self.m_count > 1:
            m_string = f"{self.m_count} of which are class M"
        else:
            m_string = "none of which are class M"

        if self.hab_count == 1 and self.m_count < 1:
            system_summary_sentences.append("There is 1 potentially habitable world in the system.")
        elif self.hab_count == 1 and self.m_count == 1:
            system_summary_sentences.append("There is 1 class M world in the system.")
        elif self.hab_count > 1 and self.hab_count == self.m_count:
            system_summary_sentences.append(f"There are {self.hab_count} class M worlds in the system.")
        elif self.hab_count > 1:
            system_summary_sentences.append(f"There are {self.hab_count} potentially habitable worlds ({m_string})")
        else:
            system_summary_sentences.append("There are no potentially habitable worlds in this system.")

        perimeter_ly = self.star.system_perimeter * constants.AU_TO_LY
        heliosphere_ly = self.star.heliosphere_radius * constants.AU_TO_LY
        if heliosphere_ly < 0.1:
            heliosphere_text = f"{self.star.heliosphere_radius:.4f} AU"
        else:
            heliosphere_text = f"{heliosphere_ly:.4f} light-years"

        system_summary_sentences.append(
            f"The star's stellar wind creates a bubble, known as the heliosphere, which extends out to approximately {heliosphere_text}.")
        system_summary_sentences.append(
            f"Beyond this, the star's gravitational influence extends out to a distance of {perimeter_ly:.2f} light-years, marking the ultimate edge of the system.")

        combined_system_summary_paragraph = to_paragraph(system_summary_sentences)


        if isinstance(self.star, BinaryStarProxy):
            # Get combined binary system data and age from the proxy
            # BinaryStarProxy.to_paragraph_list() now returns [data_block, "", age_sentence, ""]
            binary_proxy_paragraphs = self.star.to_paragraph_list()

            # 1. Append the combined binary system data block
            all_output_parts.append(binary_proxy_paragraphs[0]) # Binary System Data
            all_output_parts.append('\n\n')

            # 2. Append the system summary
            all_output_parts.append(combined_system_summary_paragraph)
            all_output_parts.append('\n\n')

            # 3. Append the combined age sentence for the binary system
            all_output_parts.append(binary_proxy_paragraphs[1]) # Binary System Age sentence

            # Add flavor text if the random chance passes and total flavor text limit is not exceeded
            if random.random() < constants.FLAVOR_CHANCE_SYSTEM and self.system_config.system_flavor_count < constants.MAX_FLAVOR_TOTAL:
                flavor_text = random.choice(constants.SYSTEM_FLAVOR)
                all_output_parts.append(f"\n\nSensors show {flavor_text}")
                self.system_config.system_flavor_count += 1

            # all_output_parts.append('\n\n')

            # 4. Append individual star details for primary and secondary stars
            for star_obj in self.stars: # self.stars contains primary and secondary Star objects
                all_output_parts.append('\n\n') # Blank line after age sentence
                header_level = '===' if not self.system_config.MARKDOWN else '###'
                all_output_parts.append(f"{header_level} {star_obj.name} {header_level if not self.system_config.MARKDOWN else ''}")
                all_output_parts.append('\n') # Add a newline after the header

                # Each individual star's to_paragraph_list() returns [data_block, "", age_sentence, ""]
                individual_star_details = star_obj.to_paragraph_list()
                all_output_parts.append(individual_star_details[0]) # Individual Star Data block
                all_output_parts.append('\n\n') # Blank line after data block
                all_output_parts.append(individual_star_details[1]) # Individual Star Age sentence

        else:
            # For single star:
            # Get single star data and age from the star object
            # Star.to_paragraph_list() returns [data_block, "", age_sentence, ""]
            single_star_paragraphs = self.star.to_paragraph_list()

            # 1. Append the single star data block
            all_output_parts.append(single_star_paragraphs[0])
            all_output_parts.append('\n\n')
            # 2. Append the single star age sentence
            all_output_parts.append(single_star_paragraphs[1])
            all_output_parts.append('\n\n')
            # 3. Append the system summary
            all_output_parts.append(combined_system_summary_paragraph)
            # all_output_parts.append('\n\n')

            # Add flavor text if the random chance passes and total flavor text limit is not exceeded
            if random.random() < constants.FLAVOR_CHANCE_SYSTEM and self.system_config.system_flavor_count < constants.MAX_FLAVOR_TOTAL:
                flavor_text = random.choice(constants.SYSTEM_FLAVOR)
                all_output_parts.append(f"\n\nSensors show {flavor_text}")
                self.system_config.system_flavor_count += 1

        # Add planet/belt paragraphs, each separated by a double newline from the previous.
        if self.planets:
            for planet in self.planets:
                all_output_parts.append('\n\n' + '\n\n'.join(planet.to_paragraph_list()))

        # Add category tag if not markdown
        if not self.system_config.MARKDOWN:
            all_output_parts.append('\n\n' + '[[Category:Star Systems]]')

        return ''.join(all_output_parts)