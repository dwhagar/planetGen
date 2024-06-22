import math
import random

STEFAN_BOLTZMANN_CONSTANT = 5.67e-8  # W/m²/K⁴
SOLAR_MASS_TO_KG = 1.989e30

class Star:
    """
    A class representing a star and it's properties.
    """

    def __init__(self, spectral_class=None, temperature=None):
        """
        Initializes a Star object with default values or user-provided values.
        """
        self.generate_star(spectral_class=spectral_class, temperature=temperature)
        self.calculate_habitable_zone()

    def __str__(self):
        """
        Returns a string representation of the star.
        """
        output = []
        output.append("{{Star Data")
        output.append(f"|type={self.type}")
        output.append(f"|radius={self.radius} km")
        output.append(f"|mass={self.mass} kg")
        output.append(f"|temp={self.temperature} K")
        output.append(f"|lum={self.luminosity} W")
        output.append(f"|hab=Between {self.habitable_zone[0]} and {self.habitable_zone[1]} AU")
        output.append("}}")

        return '\n'.join(output)
    
    def calculate_habitable_zone(self):
        """
        Calculates the inner and outer boundaries of the habitable zone for a star based on its
        luminosity measured as multiples of the sun's luminosity.  Using the Kopparapu et al. (2013)
        model for conservative habitable zone boundaries measured in AU.
        """
        inner_radius = math.sqrt(self.luminosity / 1.1)  
        outer_radius = math.sqrt(self.luminosity / 0.53)
        self.habitable_zone = (inner_radius, outer_radius)

    def generate_star(self, spectral_class=None, temperature=None):
        """
        Generates a random star's properties, optionally taking spectral class
        and temperature as input.
        """
        # Spectral class probabilities (adjust as needed)
        spectral_probabilities = {
            'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45
        }

        # Temperature range based on spectral class (Kelvin)
        temp_ranges = {
            'O': (30000, 60000), 'B': (10000, 30000), 'A': (7500, 10000),
            'F': (6000, 7500), 'G': (5200, 6000), 'K': (3700, 5200), 'M': (2400, 3700)
        }

        # Validate or generate spectral_class
        if spectral_class is None:
            spectral_class = random.choices(list(spectral_probabilities.keys()), weights=spectral_probabilities.values(), k=1)[0]
        elif spectral_class not in spectral_probabilities:
            raise ValueError("Invalid spectral class")

        # Validate or generate temperature
        valid_temp_range = temp_ranges.get(spectral_class)
        if valid_temp_range is None:
            raise ValueError("Invalid spectral class")

        if temperature is None:
            temperature = random.uniform(*valid_temp_range)
        elif not (valid_temp_range[0] <= temperature <= valid_temp_range[1]):
            raise ValueError("Temperature out of range for the given spectral class")

        # Luminosity-Radius-Temperature Relation & Mass-Luminosity Relation approximations
        luminosity = temperature**4  # Approximate Stefan-Boltzmann law
        radius = math.sqrt(luminosity / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature**4)) / 1000
        mass = (luminosity**(1/3.5) * SOLAR_MASS_TO_KG) / 1000 # Approximate Mass-Luminosity Relation

        # Yerkes spectral classification based on luminosity and radius
        if luminosity > 10000:
            yerkes_class = "Ia+"
            yerkes_type = "Luminous Supergiant"
        elif luminosity > 1000:
            yerkes_class = "Ia"
            yerkes_type = "Supergiant"
        elif luminosity > 100:
            yerkes_class = "II"
            yerkes_type = "Bright Giant"
        elif luminosity > 10:
            yerkes_class = "III"
            yerkes_type = "Giant"
        elif luminosity > 1:
            yerkes_class = "IV"
            yerkes_type = "Subgiant"
        else:
            if radius > 10:
                yerkes_class = "V"
                yerkes_type = "Main Sequence"
            else:
                yerkes_class = "VI"
                yerkes_type = "Subdwarf"

        # Create type string based on spectral class, luminosity class, and color
        color_descriptions = {
            'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 
            'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'
        }
        star_type = f"{spectral_class}{yerkes_class} {color_descriptions[spectral_class]} "
        star_type += yerkes_type + " Star"

        # Set the star's properties
        self.type = star_type
        self.radius = radius
        self.mass = mass
        self.temperature = temperature
        self.luminosity = luminosity