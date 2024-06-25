import random
import math
from .planetData import Planet
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

    def __init__(self):
        """
        Initializes a StarSystem object.
        """

        self.outer_limit = None
        self.inner_limit = None
        self.star = Star()
        self.planets = []
        system_objects = self.estimate_num_objects()
        star_factor = self.star.mass / SOLAR_MASS

        if system_objects > 0:
            for i in range(system_objects):
                # 1 in 10 chance of an asteroid belt
                normalized_distance = i / system_objects
                estimated_distance =  0.72 * math.exp(0.18 * normalized_distance) * star_factor

                if random.random() < 0.1:
                    min_distance = round(estimated_distance * 0.79, 3)
                    max_distance = round(estimated_distance * 1.21, 3)
                    self.planets.append(f"An asteroid belt orbits roughly between {min_distance} AU and {max_distance} AU.")
                else:
                    planet = Planet(self.star.habitable_zone, estimated_distance, self.star.luminosity, self.star.radius, self.star.temperature)
                    self.planets.append(planet)

    def estimate_num_objects(self):
        """
        Estimates the number of objects in a star system based on the star's mass.
        """
        solar_masses = self.star.mass / SOLAR_MASS
        max_objects = 25 * solar_masses

        # Add some randomness for variation
        num_objects = random.randint(0, math.ceil(max_objects))

        return num_objects
    
    def __str__(self):
        """
        Generates a string output for the system data.
        """
        output = [str(self.star)]

        if len(self.planets) > 0:
            for planet in self.planets:
                output.append(str(planet))
        else:
            output.append("There are no planets or asteroid belts in this system.")
        
        return '\n'.join(output)