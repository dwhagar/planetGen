# stellarObjects/systemData.py

"""
Star System Generation
======================

This module contains the `StarSystem` class, which is used to generate and
represent a full star system, including a central star and a list of planets.
"""

import random
import math
from .planetData import Planet, Asteroid_Belt
from .starData import Star
from .constants import SOLAR_MASS_TO_KG
from .utils import to_paragraph
from . import config

class StarSystem:
    """
    A class representing a star system, containing a central star and a list of planets.
    """

    def __init__(self, force_hab=False, force_belt=False, force_large=False, force_moons=False,
                 force_planets=False, absurd=False):
        """
        Initializes a StarSystem object.

        Args:
            force_hab (bool, optional): Whether to force the generation of a habitable world. Defaults to False.
            force_belt (bool, optional): Whether to force the generation of an asteroid belt. Defaults to False.
            force_large (bool, optional): Whether to force the generation of a large star. Defaults to False.
            force_moons (bool, optional): Whether to force the generation of lots of moons. Defaults to False.
            force_planets (bool, optional): Whether to force the system to have the maximum number of planets. Defaults to False.
            absurd (bool, optional): Whether to generate an absurdly large system. Defaults to False.
        """
        self.star = Star(force_large=force_large, absurd=absurd)
        self.planets = []
        system_objects = self.estimate_num_objects(force_max=force_planets)
        star_factor = self.star.mass / SOLAR_MASS_TO_KG

        required_objects = 0
        if force_hab:
            required_objects += 1
        if force_belt:
            required_objects += 1

        if system_objects < required_objects:
            system_objects = required_objects

        if system_objects > 0:
            belt_index = random.randint(0, system_objects - 1) if force_belt else -1
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
                    estimated_distance = 0.55 * star_factor

                hz = self.star.habitable_zone[0] < estimated_distance < self.star.habitable_zone[1]

                if force_hab and not found_hab:
                    if not hz and i == 0:
                        if (estimated_distance > self.star.habitable_zone[1] or
                                0 < self.star.habitable_zone[0] - estimated_distance < 0.2 or system_objects == 1):
                            estimated_distance = random.uniform(self.star.habitable_zone[0], self.star.habitable_zone[1])
                            hz = True
                    elif not hz and i > 0:
                        last_planet = self.planets[i - 1]
                        beyond_hz = last_planet.upper_limit > self.star.habitable_zone[1] if last_planet.type == 'a' else \
                            (last_planet.distance + last_planet.min_orbit_distance > self.star.habitable_zone[1])

                        if beyond_hz:
                            estimated_distance = random.uniform(self.star.habitable_zone[0], self.star.habitable_zone[1])
                            planet = Planet(self.star.habitable_zone, estimated_distance, self.star.luminosity,
                                            self.star.radius, self.star.temperature, self.star.mass,
                                            planet_class="M", force_moons=force_moons)
                            self.planets[i - 1] = planet
                            i -= 1
                            found_hab = True
                            continue
                        elif i == system_objects - 1:
                            estimated_distance = random.uniform(self.star.habitable_zone[0], self.star.habitable_zone[1])
                            hz = True

                    if hz:
                        planet = Planet(self.star.habitable_zone, estimated_distance, self.star.luminosity,
                                        self.star.radius, self.star.temperature, self.star.mass,
                                        planet_class="M", force_moons=force_moons)
                        found_hab = True
                        self.planets.append(planet)
                        continue

                if (random.random() < 0.1 or i == belt_index) and not last_asteroid and not hz:
                    min_distance = estimated_distance
                    max_distance = estimated_distance * random.uniform(1.1, 2)
                    self.planets.append(Asteroid_Belt(estimated_distance, min_distance, max_distance))
                else:
                    planet = Planet(self.star.habitable_zone, estimated_distance, self.star.luminosity,
                                    self.star.radius, self.star.temperature, self.star.mass,
                                    force_moons=force_moons)
                    if planet.planet_class == "M":
                        found_hab = True
                    self.planets.append(planet)

        self.validate_system()
        self.planet_count, self.belt_count, self.moon_count = self.count_objects()
        self.hab_count, self.m_count = self.count_habitable()

    def count_objects(self):
        """
        Counts the number of planets, asteroid belts, and moons in the system.
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
        Counts the number of potentially habitable (Class H, K, L, M, O, P) and
        Class M worlds in the system.
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

    def estimate_num_objects(self, force_max=False):
        """
        Estimates the number of objects in a star system based on the star's mass.
        """
        solar_masses = self.star.mass / SOLAR_MASS_TO_KG

        # This provides a continuous scaling factor based on mass.
        # For a 1 solar mass star, this is 1.
        # For smaller stars, it's < 1; for larger stars, it's > 1.
        # The logarithm helps temper the explosive growth for very massive stars.
        scaling_factor = 1 + math.log10(solar_masses) if solar_masses >= 1 else solar_masses

        # Base number of objects for a 1 solar mass star is 15.
        max_objects = 15 * scaling_factor
        if max_objects > 500:
            max_objects = 500
        return math.ceil(max_objects) if force_max else random.randint(0, math.ceil(max_objects))

    def validate_system(self):
        """
        Validates and adjusts the distances of stellar objects in a system to
        prevent orbital overlap.
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
                    if distance_to_last < 0.05:
                        planet.distance += 0.05 + additional_correction
                        planet.upper_limit += 0.05 + additional_correction
                        planet.lower_limit += 0.05 + additional_correction
                elif distance_to_last < last_planet.min_orbit_distance:
                    planet.distance += last_planet.min_orbit_distance + additional_correction
                    planet.upper_limit += last_planet.min_orbit_distance + additional_correction
                    planet.lower_limit += last_planet.min_orbit_distance + additional_correction
            elif last_planet.type == 'a':
                if distance_to_last < 0.05:
                    planet.distance += 0.05 + additional_correction
                    planet.calculate_atmospheric_conditions()
            else:
                min_orbit = max(planet.min_orbit_distance, last_planet.min_orbit_distance)
                if distance_to_last < min_orbit:
                    planet.distance += min_orbit + additional_correction
                    planet.calculate_atmospheric_conditions()

    def __str__(self):
        """
        Generates a string output for the system data, including a summary and
        details for each celestial body.
        """
        output = [str(self.star)]
        
        data = {
            "planets": self.planet_count,
            "belts": self.belt_count,
            "moons": self.moon_count,
            "habitable": self.hab_count,
            "m_class": self.m_count
        }

        segments = []
        if data['planets'] > 0:
            segments.append(f"{data['planets']} planet{'s' if data['planets'] > 1 else ''}")
        if data['belts'] > 0:
            segments.append(f"{data['belts']} asteroid belt{'s' if data['belts'] > 1 else ''}")
        if data['moons'] > 0:
            segments.append(f"{data['moons']} moon{'s' if data['moons'] > 1 else ''}")

        if not segments:
            system_string = "There are no stellar objects in this system."
        else:
            system_string = "This system contains " + ", ".join(segments)
            if len(segments) > 1:
                system_string = system_string.replace(f", {segments[-1]}", f" and {segments[-1]}")
            system_string += "."

        if data['m_class'] == 1:
            m_string = "1 of which is class M"
        elif data['m_class'] > 1:
            m_string = f"{data['m_class']} of which are class M"
        else:
            m_string = "none of which are class M"

        if data['habitable'] == 1 and data['m_class'] < 1:
            habitable_string = "There is 1 potentially habitable world in the system."
        elif data['habitable'] == 1 and data['m_class'] == 1:
            habitable_string = "There is 1 class M world in the system."
        elif data['habitable'] > 1 and data['habitable'] == data['m_class']:
            habitable_string = f"There are {data['habitable']} class M worlds in the system."
        elif data['habitable'] > 1:
            habitable_string = f"There are {data['habitable']} potentially habitable worlds ({m_string})"
        else:
            habitable_string = "There are no potentially habitable worlds in this system."

        sentences = [system_string, habitable_string]

        # Convert AU to Light-Years for the descriptive text (1 LY = 63241.1 AU)
        perimeter_ly = self.star.system_perimeter / 63241.1

        # Convert heliosphere radius to light-years to check the condition
        heliosphere_ly = self.star.heliosphere_radius / 63241.1
        if heliosphere_ly < 0.1:
            heliosphere_text = f"{self.star.heliosphere_radius:.4f} AU"
        else:
            heliosphere_text = f"{heliosphere_ly:.4f} light-years"

        sentences.append(f"The star's stellar wind creates a bubble, known as the heliosphere, which extends out to approximately {heliosphere_text}.")
        sentences.append(f"Beyond this, the star's gravitational influence extends out to a distance of {perimeter_ly:.2f} light-years, marking the ultimate edge of the system.")

        output.append(to_paragraph(sentences) + "\n")

        if self.planets:
            for planet in self.planets:
                output.append(str(planet) + '\n')
        if not config.MARKDOWN:
            output.append('\n[[Category:Star Systems]]')
        return '\n'.join(output)