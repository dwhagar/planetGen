import random
import math
from .planetData import Planet, Asteroid_Belt
from .starData import Star

SOLAR_MASS = 1.989e30  # kg
SOLAR_RADIUS = 6.9634e8  # m
SOLAR_GRAVITY = 28 # g's
GRAVITATIONAL_CONSTANT = 6.674e-11  # N(m/kg)^2
EARTH_G = 9.80665 # m/s^2

class StarSystem:
    """
    A class representing a star system, containing a central star and a list of planets.
    """

    def __init__(self, force_hab = False, force_belt = False, force_large = False):
        """
        Initializes a StarSystem object.
        """

        self.outer_limit = None
        self.inner_limit = None
        self.star = Star(force_large = force_large)
        self.planets = []
        system_objects = self.estimate_num_objects()
        star_factor = self.star.mass / SOLAR_MASS

        # Distances here are arbitrary and the idea is to create a system that does not need correction later
        # this is not scientific, it is highly probable I'm doing this wrong.

        required_objects = 0
        if force_hab:
            required_objects += 1
        if force_belt:
            required_objects += 1

        if system_objects < required_objects:
            system_objects += required_objects

        if system_objects > 0:
            if force_belt:
                belt_index = random.randint(0, system_objects - 1)
            else:
                belt_index = -1

            found_hab = False

            for i in range(system_objects):
                last_asteroid = False # Was the last system an asteroid belt?

                if i > 0:
                    last_planet = self.planets[i - 1]
                    random_buffer = random.uniform(0, star_factor)
                    if last_planet.type == 'a':
                        estimated_distance = last_planet.upper_limit + random_buffer * 2
                        last_asteroid = True
                    else:
                        estimated_distance = (last_planet.distance + last_planet.min_orbit_distance) + random_buffer
                else:
                    estimated_distance = 0.35 * star_factor

                hz = self.star.habitable_zone[0] < estimated_distance < self.star.habitable_zone[1]

                if force_hab and not found_hab:
                    if not hz and i == 0:
                        if (estimated_distance > self.star.habitable_zone[1] or
                                0 < self.star.habitable_zone[0] - estimated_distance < 0.2 or system_objects == 1):
                            estimated_distance = random.uniform(self.star.habitable_zone[0],
                                                                self.star.habitable_zone[1])
                            hz = self.star.habitable_zone[0] < estimated_distance < self.star.habitable_zone[1]
                    elif not hz and i > 0:
                        beyond_hz =\
                            (self.planets[i - 1].distance + self.planets[i - 1].min_orbit_distance >
                             self.star.habitable_zone[1])
                        if beyond_hz:
                            estimated_distance = random.uniform(self.star.habitable_zone[0],
                                                                self.star.habitable_zone[1])
                            planet = Planet(self.star.habitable_zone, estimated_distance,
                                            self.star.luminosity, self.star.radius, self.star.temperature,
                                            self.star.mass,
                                            planet_class="M")
                            self.planets[i - 1] = planet
                            i -= 1 # Walk in the index back since we've replaced the last planet
                            found_hab = True
                            continue

                        elif i == system_objects - 1:
                            estimated_distance = random.uniform(self.star.habitable_zone[0],
                                                                self.star.habitable_zone[1])

                    if hz:
                        planet = Planet(self.star.habitable_zone, estimated_distance,
                                        self.star.luminosity, self.star.radius, self.star.temperature, self.star.mass,
                                        planet_class="M")

                        found_hab = True
                        self.planets.append(planet)
                        continue

                if (random.random() < 0.1 or i == belt_index) and not last_asteroid:
                    if star_factor < 1:
                        min_distance = estimated_distance
                        max_distance = estimated_distance / star_factor
                    else:
                        min_distance = estimated_distance
                        max_distance = estimated_distance * star_factor

                    self.planets.append(Asteroid_Belt(estimated_distance, min_distance, max_distance))
                else:
                    planet = Planet(self.star.habitable_zone, estimated_distance,
                                    self.star.luminosity, self.star.radius, self.star.temperature, self.star.mass)

                    if planet.planet_class == "M":
                        found_hab = True

                    self.planets.append(planet)

        self.validate_system()

    def estimate_num_objects(self):
        """
        Estimates the number of objects in a star system based on the star's mass.
        """
        solar_masses = self.star.mass / SOLAR_MASS
        max_objects = 25 * solar_masses

        # Add some randomness for variation
        num_objects = random.randint(0, math.ceil(max_objects))

        return num_objects

    def validate_system(self):
        """
        This function will validate (and adjust) the distances of stellar objects in a system.
        """
        # We are going to use the Hill Radius to determine how far out from a planet before
        # the planet's gravity doesn't dominate over the stars own gravity.

        # If there aren't at least 2 things in the system, no need to do this.
        if len(self.planets) < 2:
            return

        for idx in range(len(self.planets)):
            # If it's an asteroid belt it'll be a string.
            planet = self.planets[idx]

            if idx > 0:
                last_planet = self.planets[idx - 1]
                if last_planet.type == 'a' and planet.type == 'a':
                    distance_to_last = planet.lower_limit - last_planet.upper_limit
                elif last_planet.type == 'a':
                    distance_to_last = planet.distance - last_planet.upper_limit
                else:
                    distance_to_last = planet.distance - last_planet.distance
            else:
                continue

            # For large star systems it's possible this algorithm will move planets out of order
            # this should correct that by adding a correction to all the values so that they will
            # ensure the distance is corrected if negative.
            if distance_to_last < 0:
                distance_to_last = abs(distance_to_last)
                additional_correction = last_planet.distance
            else:
                additional_correction = 0

            ## TODO: Asteroid belt orbital distances just don't work, need to find a way to validate

            # If both are asteroid belts (highly unlikely)
            if planet.type == 'a' and last_planet.type == 'a':
                if distance_to_last < 0.05:
                    planet.distance += 0.05 + additional_correction
                    planet.upper_limit += 0.05 + additional_correction
                    planet.lower_limit += 0.05 + additional_correction
            # If the current planet is an asteroid belt.
            elif planet.type == 'a':
                if distance_to_last < last_planet.min_orbit_distance:
                    planet.distance += last_planet.min_orbit_distance + additional_correction
                    planet.upper_limit += last_planet.min_orbit_distance + additional_correction
                    planet.lower_limit += last_planet.min_orbit_distance + additional_correction
            # If the last planet is an asteroid belt.
            elif last_planet.type == 'a':
                if distance_to_last < 0.05:
                    planet.distance += 0.05 + additional_correction
                    planet.calculate_atmospheric_conditions()
            # If both are planets (most likely)
            else:
                if planet.min_orbit_distance < last_planet.min_orbit_distance:
                    min_orbit = planet.min_orbit_distance
                else:
                    min_orbit = last_planet.min_orbit_distance

                if distance_to_last < min_orbit:
                    planet.distance += min_orbit + additional_correction
                    planet.calculate_atmospheric_conditions()

    def __str__(self):
        """
        Generates a string output for the system data.
        """
        output = [str(self.star)]

        if len(self.planets) == 1:
            output.append(f"There is 1 stellar object in the system, not including any moons.\n")
            for planet in self.planets:
                output.append(str(planet) + '\n')
        elif len(self.planets) > 1:
            output.append(f"There are {len(self.planets)} stellar objects in the system, not including any moons.\n")
            for planet in self.planets:
                output.append(str(planet) + '\n')
        else:
            output.append("There are no planets or asteroid belts in this system.")

        output.append('\n[[Category:Star Systems]]')
        
        return '\n'.join(output)