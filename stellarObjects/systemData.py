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

import copy
import math
import random

from .asteroidData import AsteroidBelt
from .config import SystemConfig
from .doubleStar import BinaryStarProxy
from . import physical_constants, planetLife, planetPhysics, program_constants
from .planetData import Planet
from .starData import Star
from .utils import to_paragraph

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
        `config.HABITABLE_WORLD = True` ensures that at least one habitable planet
        is generated, while `config.PLANETS = False` can be used to create a star
        with no orbiting bodies. `config.NUM_ORBITS` and `config.SLOTS` allow the
        number of orbits and the exact contents of specific orbital slots to be
        specified explicitly, overriding the random placement logic for those
        slots. The placement of planets is done sequentially, with each new
        planet's orbit being determined based on the position of the previous one
        to ensure a degree of realism in orbital spacing.

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

        if self.system_config.BINARY_SYSTEM:
            # Create a copy of the system_config for the secondary star
            secondary_star_config = copy.deepcopy(self.system_config)
            # Ensure LARGE_STAR is not forced for the secondary star
            secondary_star_config.LARGE_STAR = False

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
        star_factor = self.star.mass / physical_constants.SOLAR_MASS_TO_KG

        required_objects = 0
        if self.system_config.HABITABLE_WORLD is True:
            required_objects += 1
        if self.system_config.ASTEROID_BELT is True:
            required_objects += 1

        if system_objects < required_objects:
            system_objects = required_objects

        slots = self.system_config.SLOTS or []

        if system_objects > 0:
            belt_index = random.randint(0, system_objects - 1) if self.system_config.ASTEROID_BELT is True else -1
            found_hab = False
            i = -1

            while i < system_objects - 1:
                i += 1
                last_asteroid = False
                slot_spec = slots[i] if i < len(slots) else None
                prev_slot_explicit = i > 0 and (i - 1) < len(slots) and slots[i - 1] is not None

                if i > 0:
                    last_planet = self.planets[i - 1]
                    random_buffer = random.uniform(0, star_factor)
                    if last_planet.type == 'a':
                        estimated_distance = last_planet.upper_limit + random_buffer * 2
                        last_asteroid = True
                    else:
                        estimated_distance = (last_planet.distance + last_planet.min_orbit_distance) + random_buffer
                else:
                    estimated_distance = program_constants.INITIAL_PLANET_DISTANCE_FACTOR * star_factor

                hz = self.star.habitable_zone[0] < estimated_distance < self.star.habitable_zone[1]

                # An explicit per-slot specification takes priority over all of the
                # normal random/forced generation logic below.
                if slot_spec is not None:
                    obj = self.generate_slot_object(slot_spec, estimated_distance)
                    if getattr(obj, 'planet_class', None) in program_constants.HABITABLE_PLANET_CLASSES:
                        found_hab = True
                    self.planets.append(obj)
                    continue

                if self.system_config.HABITABLE_WORLD is True and not found_hab:
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

                        if beyond_hz and not prev_slot_explicit:
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

                if self.system_config.ASTEROID_BELT is not False and (random.random() < program_constants.ASTEROID_BELT_PROBABILITY or i == belt_index) and not last_asteroid and not hz:
                    min_distance = estimated_distance
                    max_distance = estimated_distance * random.uniform(program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN, program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)
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

    def generate_slot_object(self, slot_spec, estimated_distance):
        """
        Builds the celestial object explicitly requested for one orbital slot
        by `system_config.SLOTS`.

        Args:
            slot_spec (dict): The slot specification, with a required "type"
                              key ("planet" or "asteroid_belt") and optional
                              "planet_class"/"moons" keys (see `SystemConfig.SLOTS`).
            estimated_distance (float): The orbital distance (in AU) computed
                                        for this slot by the generation loop.

        Returns:
            Planet or AsteroidBelt: The generated object for this slot.

        Raises:
            ValueError: If `slot_spec["type"]` is not "planet" or "asteroid_belt".
        """
        slot_type = slot_spec.get("type", "planet")

        if slot_type == "asteroid_belt":
            min_distance = estimated_distance
            max_distance = estimated_distance * random.uniform(program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN, program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)
            return AsteroidBelt(self.system_config, estimated_distance, min_distance, max_distance)

        if slot_type == "planet":
            planet_class = slot_spec.get("planet_class")
            distance = self.calculate_distance_for_class(planet_class, estimated_distance)
            return Planet(self.system_config, self.star, self.star.habitable_zone, distance, self.star.type[0],
                         self.star.luminosity, self.star.radius, self.star.temperature, self.star.mass,
                         planet_class=planet_class, moon_count=slot_spec.get("moons"))

        raise ValueError(f"Invalid slot type '{slot_type}'; expected 'planet' or 'asteroid_belt'.")

    def calculate_distance_for_class(self, planet_class, estimated_distance):
        """
        Adjusts an orbital distance so that it falls in a zone (hot, cold, or
        ecosphere) that actually supports a requested `planet_class`.

        The generation loop picks `estimated_distance` before knowing what
        the slot will contain, so a user-requested class (e.g. "M", which
        only exists in the ecosphere) may not be valid at that distance. This
        nudges the distance into a supporting zone, the same way the existing
        forced-habitable-world logic snaps a planet's distance into the
        habitable zone, so an explicit `planet_class` request doesn't fail
        with an "Invalid planet class for this zone" error just because of
        where its slot happened to land in the orbit sequence.

        Args:
            planet_class (str or None): The requested planet class, or None
                                        if the slot doesn't specify one.
            estimated_distance (float): The orbital distance (in AU) computed
                                        for this slot by the generation loop.

        Returns:
            float: `estimated_distance`, or an adjusted distance (in AU) that
                  falls within a zone supporting `planet_class`.
        """
        class_data = program_constants.PLANET_CLASSES.get(planet_class)
        if class_data is None:
            return estimated_distance

        inner, outer = self.star.habitable_zone
        if estimated_distance < inner:
            zone = 'h'
        elif estimated_distance > outer:
            zone = 'c'
        else:
            zone = 'e'

        if class_data.get(zone):
            return estimated_distance

        if class_data.get('e'):
            return random.uniform(inner, outer)
        if class_data.get('h'):
            return random.uniform(inner * 0.05, inner * 0.95)
        if class_data.get('c'):
            return outer * random.uniform(1.05, 3.0)

        # No zone supports this class; leave the distance as-is and let
        # planetPhysics raise its usual, clearer validation error.
        return estimated_distance

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
        classification against `program_constants.HABITABLE_PLANET_CLASSES` to
        determine if they are potentially habitable, and also keeps a separate
        count of Class M worlds, which are considered the most Earth-like. This
        is used for the system summary output.

        Returns:
            tuple: A tuple containing the total number of potentially habitable
                   worlds and the total number of Class M worlds, in that order.
        """
        habitable_classes = program_constants.HABITABLE_PLANET_CLASSES
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
        be a random value between the minimum and calculated maximum, or the
        maximum/minimum itself if `MAX_PLANETS` is True/False.

        `NUM_ORBITS`, when set, bypasses this estimate entirely and is
        returned as-is. `PLANETS` set to False forces zero objects; set to
        True, it raises the minimum considered to 1.

        Returns:
            int: The estimated number of objects to be generated in the system.
        """
        if self.system_config.PLANETS is False:
            return 0

        if self.system_config.NUM_ORBITS is not None:
            return self.system_config.NUM_ORBITS

        solar_masses = self.star.mass / physical_constants.SOLAR_MASS_TO_KG

        # This provides a continuous scaling factor based on mass.
        # For a 1 solar mass star, this is 1.
        # For smaller stars, it's < 1; for larger stars, it's > 1.
        # The logarithm helps temper the explosive growth for very massive stars.
        scaling_factor = 1 + math.log10(solar_masses) if solar_masses >= 1 else solar_masses

        # Base number of objects for a 1 solar mass star is 15.
        max_objects = program_constants.BASE_MAX_SYSTEM_OBJECTS * scaling_factor
        if max_objects > program_constants.ABSOLUTE_MAX_SYSTEM_OBJECTS:
            max_objects = program_constants.ABSOLUTE_MAX_SYSTEM_OBJECTS
        max_objects = math.ceil(max_objects)

        min_objects = 1 if self.system_config.PLANETS is True else 0
        max_objects = max(max_objects, min_objects)

        if self.system_config.MAX_PLANETS is True:
            return max_objects
        if self.system_config.MAX_PLANETS is False:
            return min_objects
        return random.randint(min_objects, max_objects)

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
                    if distance_to_last < program_constants.MIN_ASTEROID_BELT_SEPARATION:
                        planet.distance += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        planet.upper_limit += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        planet.lower_limit += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                elif distance_to_last < last_planet.min_orbit_distance:
                    planet.distance += last_planet.min_orbit_distance + additional_correction
                    planet.upper_limit += last_planet.min_orbit_distance + additional_correction
                    planet.lower_limit += last_planet.min_orbit_distance + additional_correction
            else:
                if last_planet.type == 'a':
                    if distance_to_last < program_constants.MIN_ASTEROID_BELT_SEPARATION:
                        planet.distance += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
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

        if segments:
            system_string = "This system contains " + ", ".join(segments)
            if len(segments) == 2:
                system_string = system_string.replace(f", {segments[-1]}", f" and {segments[-1]}")
            elif len(segments) > 2:
                system_string = system_string.replace(f", {segments[-1]}", f", and {segments[-1]}")
            system_string += "."
            system_summary_sentences.append(system_string)
        elif self.planet_count == 0 and self.belt_count == 0 and self.moon_count == 0:
            # No segment was built and the system is genuinely empty (as opposed to having
            # its planet count omitted above to avoid redundancy with the habitability
            # sentence that follows).
            system_summary_sentences.append("There are no stellar objects in this system.")

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

        perimeter_ly = self.star.system_perimeter * physical_constants.AU_TO_LY
        heliosphere_ly = self.star.heliosphere_radius * physical_constants.AU_TO_LY
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
            # BinaryStarProxy.to_paragraph_list() returns [data_block, age_sentence]
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
            if random.random() < program_constants.FLAVOR_CHANCE_SYSTEM and self.system_config.system_flavor_count < program_constants.MAX_FLAVOR_TOTAL:
                flavor_text = random.choice(program_constants.SYSTEM_FLAVOR)
                all_output_parts.append(f"\n\nSensors show {flavor_text}")
                self.system_config.system_flavor_count += 1

            # 4. Append individual star details for primary and secondary stars
            for star_obj in self.stars: # self.stars contains primary and secondary Star objects
                all_output_parts.append('\n\n') # Blank line after age sentence
                header_level = '===' if not self.system_config.MARKDOWN else '###'
                all_output_parts.append(f"{header_level} {star_obj.name} {header_level if not self.system_config.MARKDOWN else ''}")
                all_output_parts.append('\n') # Add a newline after the header

                # Each individual star's to_paragraph_list() returns [data_block, age_sentence]
                individual_star_details = star_obj.to_paragraph_list()
                all_output_parts.append(individual_star_details[0]) # Individual Star Data block
                all_output_parts.append('\n\n') # Blank line after data block
                all_output_parts.append(individual_star_details[1]) # Individual Star Age sentence

        else:
            # For single star:
            # Get single star data and age from the star object
            # Star.to_paragraph_list() returns [data_block, age_sentence]
            single_star_paragraphs = self.star.to_paragraph_list()

            # 1. Append the single star data block
            all_output_parts.append(single_star_paragraphs[0])
            all_output_parts.append('\n\n')
            # 2. Append the single star age sentence
            all_output_parts.append(single_star_paragraphs[1])
            all_output_parts.append('\n\n')
            # 3. Append the system summary
            all_output_parts.append(combined_system_summary_paragraph)

            # Add flavor text if the random chance passes and total flavor text limit is not exceeded
            if random.random() < program_constants.FLAVOR_CHANCE_SYSTEM and self.system_config.system_flavor_count < program_constants.MAX_FLAVOR_TOTAL:
                flavor_text = random.choice(program_constants.SYSTEM_FLAVOR)
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