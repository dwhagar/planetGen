import random
from .planetData import Planet
from .starData import Star

class StarSystem:
    """
    A class representing a star system, containing a central star and a list of planets.
    """

    def __init__(self):
        """
        Initializes a StarSystem object.
        """

        self.star = Star()
        self.planets = []

        self.system_bounds()
        system_objects = self.define_system_size()

        for i in range(system_objects):
            normalized_distance = i / (system_objects - 1)  # Calculate normalized distance
            estimated_distance = self.estimate_distance(normalized_distance)
            estimated_distance *= random.uniform(0.5, 2.0)  # Random multiplier for variation

            planet = Planet(self.star.habitable_zone, estimated_distance)
            self.planets.append(planet)

    def system_bounds(self):
        """
        Estimates the inner and outer bounds of a planetary system based on the star's mass.
        """

        # Constants based on our solar system
        INNER_LIMIT_FACTOR = 0.05  # Mercury's distance in AU relative to the Sun's mass
        OUTER_LIMIT_FACTOR = 30.0  # Neptune's distance in AU relative to the Sun's mass

        # Calculate inner and outer limits
        inner_limit = INNER_LIMIT_FACTOR * self.star.mass
        outer_limit = OUTER_LIMIT_FACTOR * self.star.mass

        self.inner_limit = inner_limit
        self.outer_limit = outer_limit

    def define_system_size(self):
        """
        Based on size of the system in AU's, it generates the size of the system by number of
        stellar objects.
        """
        # Density of stellar objects within our solar system as a function of AU per object
        # a random factor is included to make it less predictable.
        object_density = (9 / 30) * random.uniform(0.5, 2.0)

        total_space = self.outer_limit - self.inner_limit
        total_objects = object_density * total_space
        print(total_objects)
        exit()
        return int(total_objects)

    def estimate_distance(self, normalized_distance):
        """
        Estimates a planet's distance from the Sun based on a normalized distance.

        Args:
            normalized_distance (float): A value between 0 inner-most planet and 1 outer-most planet.

        Returns:
            float: The estimated distance from the main star in AU.
        """
        return 0.35 * 1.11**normalized_distance
    
    def __str__(self):
        """
        Generates a string output for the system data.
        """
        output = []
        output.append(str(self.star))

        for planet in self.planets:
            output.append(str(planet))
        
        return '\n'.join(output)