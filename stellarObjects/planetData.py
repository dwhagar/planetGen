import math
import random

EARTH_RADIUS_KM = 6371

planet_density = { # g per cubic centimeter
    "Terrestrial": (3.93, 5.51), # Range from Mars to Earth
    "Gas Giant": (0.69, 1.64)    # Range from Saturn to Neptune
}

atmosphere_density = { # kg per cubic meter
    "Terrestrial": (0.020, 67),  # Range from Mars to Venus (kg/m³)
    "Gas Giant": (0.69, 1.33),   # Approximate range using Jupiter and Saturn's overall densities
}

planet_classes = {
    "A": {
        "description": "Small, barren, volcanic",
        "composition": "Igneous silica, basalt",
        "radius_range": (500, 5000), # In meters
        "h": True,  # Hot Zone
        "e": True,  # Ecosphere
        "c": True,  # Cold Zone
        "atmosphere": "Sulfur dioxide, carbon dioxide",
        "type": "t"
    },
    "B": {
        "description": "Small, molten, thin atmosphere",
        "composition": "Iron, potassium, silicon",
        "radius_range": (500, 5000),
        "h": True, 
        "e": False,
        "c": False,
        "atmosphere": "Helium, sodium, oxygen",
        "type": "t"
    },
    "C": {
        "description": "Dead worlds, no atmosphere",
        "composition": "Anthracite, basalt, hydrocarbons",
        "radius_range": (500, 10000),
        "h": True,  # Hot Zone
        "e": True,  # Ecosphere
        "c": True,  # Cold Zone
        "atmosphere": None,
        "type": "t"
    },
    "D": {
        "description": "Small icy bodies (not true planets)",
        "composition": "Frozen hydrocarbons, ice",
        "radius_range": (50, 3000),
        "h": True,
        "e": True,
        "c": True,
        "atmosphere": None,  # Or very tenuous
        "type": "t"
    },
    "E": {
        "description": "Molten core & crust, thin atmosphere",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": True,
        "e": True,
        "c": True,
        "atmosphere": "Hydrogen compounds",
        "type": "t"
    },
    "F": {
        "description": "Volcanic, shallow seas, bacterial life",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": True,
        "e": True,
        "c": True,
        "atmosphere": "Carbon dioxide, ammonia, methane",
        "type": "t"
    },
    "G": {
        "description": "Rocky, barren, simple life",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": True,
        "e": True,
        "c": True,
        "atmosphere": "Carbon dioxide, oxygen, nitrogen",
        "type": "t"
    },
    "H": {
        "description": "Desert worlds, minimal water",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Oxygen, nitrogen, argon, metals",
        "type": "t"
    },
    "I": {
        "description": "Ice giant, tilted magnetic field",
        "composition": "Rock, ice, methane, ammonia",
        "radius_range": (15000, 50000),
        "h": False,
        "e": False,
        "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "J": {
        "description": "Gas giant, turbulent atmosphere, rings",
        "composition": "Hydrogen, helium",
        "radius_range": (25000, 250000),
        "h": False,
        "e": False,
        "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "K": {
        "description": "Adaptable, thin atmosphere",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (2500, 7500),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Carbon dioxide",
        "type": "t"
    },
    "L": {
        "description": "Marginal, varied atmosphere, vegetation",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 7500),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Argon, oxygen, trace elements",
        "type": "t"
    },
    "M": {
        "description": "Terrestrial, Earth-like",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Oxygen, nitrogen, argon",
        "type": "t"
    },
    "N": {
        "description": "Reducing, hot, dense atmosphere",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Carbon dioxide, sulfides",
        "type": "t"
    },
    "O": {
        "description": "Pelagic, mostly water",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Oxygen, nitrogen, argon",
        "type": "t"
    },
    "P": {
        "description": "Glaciated, cold",
        "composition": "Silicon, iron, magnesium, ice",
        "radius_range": (5000, 10000),
        "h": False,
        "e": True,
        "c": True,
        "atmosphere": "Oxygen, nitrogen, argon (thinning with age)",
        "type": "t"
    },
    "Q": {
        "description": "Eccentric orbit, extreme variations",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (2000, 7500),
        "h": True,
        "e": True,
        "c": True,
        "atmosphere": "Varies widely (nitrogen, oxygen, argon, thin to dense)",
        "type": "t"
    },
    "R": {
        "description": "Ejected, geologically active",
        "composition": "Silicate compounds, iron",
        "radius_range": (7500, 10000),
        "h": False,  # Technically can be in any zone due to ejection
        "e": False, 
        "c": False, 
        "atmosphere": "Volcanic outgassing",
        "type": "t"
    },
    "S": {
        "description": "Supergiant, shield for inner planets",
        "composition": "Hydrogen, helium",
        "radius_range": (250000, 50000000),
        "h": False,
        "e": False,
        "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "T": {
        "description": "Gas dwarf, thick atmosphere",
        "composition": "Hydrogen, helium, hydrocarbons",
        "radius_range": (250000, 25000000),
        "h": False,
        "e": False,
        "c": True,
        "atmosphere": "Hydrogen, helium, hydrocarbons",
        "type": "g"
    },
    "U": {
        "description": "Ultragiant, can become a star",
        "composition": "Hydrogen, helium",
        "radius_range": (25000000, 60000000),
        "h": False,
        "e": False,
        "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "V": {
        "description": "Super-Earth, high gravity",
        "composition": "Iron, iridium, tungsten, nickel",
        "radius_range": (10000, 15000),
        "h": True,
        "e": True,
        "c": True,
        "atmosphere": "Carbon dioxide, oxygen, hydrogen, helium",
        "type": "t"
    },
    "W": {
        "description": "Tidally locked, extreme variations",
        "composition": "Iron, potassium, silicon",
        "radius_range": (500, 10000),
        "h": True,
        "e": True,
        "c": False,
        "atmosphere": "Oxygen, sodium, hydrogen",
        "type": "t"
    },
    "X": {
        "description": "Stripped core, no atmosphere",
        "composition": "Molten iron",
        "radius_range": (500, 5000),
        "h": True,
        "e": False,
        "c": False,
        "atmosphere": None,
        "type": "t"
    },
    "Y": {
        "description": "Small, rocky planet", 
        "composition": "Silicates, metals",
        "radius_range": (3000, 5000), 
        "h": True,
        "e": True,
        "c": True,
        "atmosphere": "Usually thin or absent",
        "type": "t"
    },
    "Z": {
        "description": "Demon class, toxic atmosphere",
        "composition": "Molten iron, sulfur, deuterium",
        "radius_range": (5000, 7500),
        "h": True,
        "e": False,
        "c": False,
        "atmosphere": "Turbulent, toxic, radiation",
        "type": "t"
    }
}

class Planet:
    """
    A Class representing a single planet and all of it's properties.
    """

    def __init__(self, hab_zone, distance, radius=None, planet_class=None):
        """
        Initializes a Planet object with its radius and the spectral class of its host star.

        Args:
            hab_zone: A tuple with where a star's habitable zone is in AU.
            distance: A float with the distance from the star in AU.
            radius: The radius of the planet in Earth km.
            planet_class: The classification of the planet.
        """

        self.radius = radius
        self.planet_class = planet_class
        self.distance = distance

        # From the star, should not be changed.
        self.hab = hab_zone

        # Calculate additional properties.
        self.volume = (4/3) * math.pi * self.radius_km**3  # Calculate volume in km^3
        self.mass = self.volume * self.density * (10**-12) # Convert density from g/cm^3 to kg/km^3 and calculate mass in Earth masses

    def generate_planet(self, stellar_mass):
        """
        Generates random planet properties (class, composition, atmosphere) based on distance and radius.
        """

        # Habitable Zone
        inner_bound = self.hab[0]  
        outer_bound = self.hab[1]

        # Determine zone based on distance
        if self.orbital_distance < inner_bound:
            zone = 'h'
        elif self.orbital_distance > outer_bound:
            zone = 'c'
        else:
            zone = 'e'

        # If planet_class and radius are None, generate them randomly
        if self.planet_class is None and self.radius is None:
            # Filter possible planet classes based on zone
            possible_classes = [
                c for c, data in planet_classes.items() if data[zone]
            ]

            # Choose a random planet class
            self.planet_class = random.choice(possible_classes)

            # Generate random radius within the allowed range for the chosen class
            min_radius, max_radius = planet_classes[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius) / EARTH_RADIUS_KM # convert km radius to Earth radii

        # If planet_class is present but radius is None, generate radius based on class
        elif self.planet_class is not None and self.radius is None:
            # Validate planet_class
            if self.planet_class not in planet_classes:
                raise ValueError("Invalid planet class")

            # Generate random radius within the allowed range for the chosen class
            min_radius, max_radius = planet_classes[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius) / EARTH_RADIUS_KM

        # If radius is present but planet_class is None, determine possible classes and choose randomly
        elif self.planet_class is None and self.radius is not None:
            # Filter possible planet classes based on zone and radius
            possible_classes = [
                c for c, data in planet_classes.items()
                if data[zone] and data["radius_range"][0] <= self.radius * EARTH_RADIUS_KM <= data["radius_range"][1]
            ]

            if not possible_classes:
                raise ValueError("No valid planet class for the given radius in this zone")

            # Choose a random planet class from the filtered list
            self.planet_class = random.choice(possible_classes)
            
        # If both planet_class and radius are present, validate them
        else:
            # Validate planet_class
            if self.planet_class not in planet_classes:
                raise ValueError("Invalid planet class")

            # Validate radius
            min_radius, max_radius = planet_classes[self.planet_class]["radius_range"]
            if not (min_radius <= self.radius * EARTH_RADIUS_KM <= max_radius):
                raise ValueError("Radius out of range for the given planet class")


        # Get the properties for the chosen planet class
        class_data = planet_classes[self.planet_class]

        # Update planet properties based on class data
        self.composition = class_data["composition"]

        # Determine planet type (Terrestrial or Gas Giant) based on class data
        planet_type = class_data["type"]  # Assuming the planet type is "t" or "g" in class_data

        # Generate random density based on planet type
        min_density, max_density = planet_density[planet_type]
        p_density = random.uniform(min_density, max_density)

        # Handle atmosphere (set to "None" if it's None in class_data)
        if class_data["atmosphere"] is None:
            self.atmosphere = "None"
            a_density = 0.0  # No atmosphere, density is 0
        else:
            self.atmosphere = class_data["atmosphere"]

            # Approximate atmospheric density based on planet type (This can be refined later)
            min_a_density, max_a_density = atmosphere_density[planet_type]
            a_density = random.uniform(min_a_density, max_a_density)

        # Recalculate mass based on the new density
        self.volume = (4/3) * math.pi * self.radius_km**3  # Calculate volume in km^3
        self.mass = self.volume * p_density * (10**-12)  # Update mass calculation
        self.density = p_density
        self.atm_density = a_density