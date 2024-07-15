import math
import random

STEFAN_BOLTZMANN_CONSTANT = 5.67e-8  # W/m²/K⁴
SOLAR_MASS_TO_KG = 1.989e30
SOLAR_LUMINOSITY = 3.82e26

luminosity_ranges = {
    'O': (30000, 1000000),
    'B': (25, 30000),
    'A': (5, 25),
    'F': (1.5, 5),
    'G': (0.6, 1.5),
    'K': (0.08, 0.6),
    'M': (0.00008, 0.08)
}

temp_ranges = {
    'O': (30000, 60000),
    'B': (10000, 30000),
    'A': (7500, 10000),
    'F': (6000, 7500),
    'G': (5200, 6000),
    'K': (3700, 5200),
    'M': (2400, 3700)
}

def to_scientific_notation(number, precision=2):
    """
    Converts a number to scientific notation with the specified precision.

    Args:
        number (float): The number to convert.
        precision (int, optional): The number of decimal places to show. Defaults to 2.

    Returns:
        str: The number in scientific notation (e.g., "1.23e+04").
    """
    if number == 0:
        return "0"
    exponent = int(math.floor(math.log10(abs(number))))
    coefficient = number / 10**exponent
    output = f"Exp|{coefficient:.{precision}f}|{exponent:d}"
    return "{{" + output + "}}"

class Star:
    """
    A class representing a star and it's properties.
    """

    def __init__(self, spectral_class=None, temperature=None, force_large=False, absurd=False):
        """
        Initializes a Star object with default values or user-provided values.
        """
        self.luminosity = None
        self.temperature = None
        self.mass = None
        self.radius = None
        self.type = None
        self.habitable_zone = None
        self.generate_star(force_large=force_large, spectral_class=spectral_class, temperature=temperature, absurd=absurd)
        self.calculate_habitable_zone()

    def __str__(self):
        """
        Returns a string representation of the star.
        """
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

        sol_mass = round(self.mass / SOLAR_MASS_TO_KG, 1)
        sol_lum = round(self.luminosity / SOLAR_LUMINOSITY, 1)

        if sol_mass <= 0:
            sol_mass = round(self.mass / SOLAR_MASS_TO_KG * 100, 4)
        if sol_lum <= 0:
            sol_lum = round(self.luminosity / SOLAR_LUMINOSITY * 100, 4)

        if sol_mass <= 2:
            mass_string = f"|mass={to_scientific_notation(self.mass)} kg ({sol_mass * 100}% of Sol)"
        elif sol_mass > 100:
            exponent = int(math.floor(math.log10(abs(sol_mass))))
            mass_string = f"|mass={to_scientific_notation(self.mass)} kg (10<sup>{exponent}</sup>x Sol)"
        else:
            mass_string = f"|mass={to_scientific_notation(self.mass)} kg ({sol_mass}x Sol)"

        if sol_lum <= 2:
            lum_string = f"|lum={to_scientific_notation(self.luminosity)} W ({sol_lum * 100}% of Sol)"
        elif sol_lum > 100:
            exponent = int(math.floor(math.log10(abs(sol_lum))))
            lum_string = f"|lum={to_scientific_notation(self.luminosity)} W (10<sup>{exponent}</sup>x Sol)"
        else:
            lum_string = f"|lum={to_scientific_notation(self.luminosity)} W ({sol_lum}x Sol)"

        if self.radius <= 100000:
            radius_string = f"|radius={round(self.radius, 2):,} km"
        else:
            radius_string = f"|radius={to_scientific_notation(self.radius, 2)} km"

        output = ["{{Star Data", f"|type={self.type}",
                  radius_string,
                  mass_string,
                  f"|temp={self.temperature} K",
                  lum_string,
                  f"|hab=Between {hab_lower} and {hab_upper} AU",
                  "}}"]
        return '\n'.join(output)
    
    def calculate_habitable_zone(self):
        """
        Calculates the inner and outer boundaries of the habitable zone for a star based on its
        luminosity measured as multiples of the sun's luminosity.  Using the Kopparapu et al. (2013)
        model for conservative habitable zone boundaries measured in AU.
        """
        solar_lum = self.luminosity / SOLAR_LUMINOSITY

        inner_radius = math.sqrt(solar_lum / 1.1)  
        outer_radius = math.sqrt(solar_lum / 0.53)
        self.habitable_zone = (inner_radius, outer_radius)

    def calculate_stellar_radius(self, luminosity, temperature):
        """
        Calculates the radius of a star in meters given its luminosity in watts and temperature in Kelvin.
        @param luminosity: Luminosity in Watts
        @param temperature: Temperature in Kelvin
        @return: Radius in Meters
        """
        # Calculate radius in meters
        radius_meters = math.sqrt(luminosity / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature ** 4))

        return radius_meters

    def set_radius_bounds(self, luminosity, temperature, spectral_class):
        """
        Checks if a luminosity is reasonable for a given temperature and spectral class.

        Args:
            luminosity (float): The luminosity of the star in solar luminosities.
            temperature (float): The temperature of the star in Kelvin.
            spectral_class (str): The spectral class of the star ('O', 'B', 'A', 'F', 'G', 'K', or 'M').

        Returns:
            tuple: minimum radius in meters, maximum radius in meters
        """

        # Convert luminosity to watts (using solar luminosity)
        luminosity_watts = luminosity * SOLAR_LUMINOSITY

        # Calculate radius in meters using Stefan-Boltzmann Law
        radius_meters = self.calculate_stellar_radius(luminosity_watts, temperature)

        # Get allowed radius range for the spectral class
        min_radius, max_radius = luminosity_ranges[spectral_class]

        # Convert allowed radius range from solar radii to meters
        min_radius_meters = min_radius * 6.957e8  # 1 solar radius = 6.957e8 meters
        max_radius_meters = max_radius * 6.957e8

        return min_radius_meters, max_radius_meters

    def generate_star(self, spectral_class=None, temperature=None, force_large=False, absurd=False):
        """
        Generates a random star's properties, optionally taking spectral class
        and temperature as input.
        """
        # Spectral class probabilities (adjust as needed)
        if absurd:
            # If set, the flag to force a large system should emphasize stars larger than Earth's Sun.
            spectral_probabilities = {
                'O': 100, 'B': 0, 'A': 0, 'F': 0, 'G': 0, 'K': 0, 'M': 0
            }
        elif force_large:
            # If set, the flag to force a large system should emphasize stars larger than Earth's Sun.
            spectral_probabilities = {
                'O': 10, 'B': 20, 'A': 30, 'F': 30, 'G': 10, 'K': 0, 'M': 0
            }
        else:
            spectral_probabilities = {
                'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45
            }

        # Validate or generate spectral_class
        if spectral_class is None:
            # noinspection PyTypeChecker
            spectral_class = random.choices(list(spectral_probabilities.keys()), weights=spectral_probabilities.values(), k=1)[0]
        elif spectral_class not in spectral_probabilities:
            raise ValueError("Invalid spectral class")

        # Validate or generate temperature
        valid_temp_range = temp_ranges.get(spectral_class)
        if valid_temp_range is None:
            raise ValueError("Invalid spectral class")

        if temperature is None:
            if absurd:
                temperature = valid_temp_range[1]
            else:
                temperature = int(round(random.uniform(*valid_temp_range), -2))
        elif not (valid_temp_range[0] <= temperature <= valid_temp_range[1]):
            raise ValueError("Temperature out of range for the given spectral class")

        # Generate the Luminosity
        min_luminosity, max_luminosity = luminosity_ranges[spectral_class]

        if absurd:
            luminosity = max_luminosity
        else:
            luminosity = random.uniform(min_luminosity, max_luminosity)

        # Validate Luminosity and Temperature
        radius_min, radius_max = self.set_radius_bounds(luminosity, temperature, spectral_class)

        # Calculate spectral subclass (0-9) based on temperature
        min_temp, max_temp = temp_ranges[spectral_class]
        temp_range_size = max_temp - min_temp
        subclass = 9 - round((temperature - min_temp) / temp_range_size * 9)  # Inverted scale (0 is hottest)

        luminosity_watts = luminosity * SOLAR_LUMINOSITY
        radius = math.sqrt(luminosity_watts / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature ** 4)) / 1000
        if radius > radius_max:
            radius = radius_max
        elif radius < radius_min:
            radius = radius_min

        mass = (luminosity**(1/3.5) * SOLAR_MASS_TO_KG) # Approximate Mass-Luminosity Relation

        # Yerkes spectral classification based on luminosity and radius
        if luminosity > 100000:  
            yerkes_class = "0"         # Hypergiant
            yerkes_type = "Hypergiant"
        elif luminosity > 30000: 
            yerkes_class = "Ia+"       # Luminous Supergiant
            yerkes_type = "Luminous Supergiant"
        elif luminosity > 10000:       
            yerkes_class = "Ia"        # Supergiant
            yerkes_type = "Supergiant"
        elif luminosity > 1000:
            yerkes_class = "Ib"         # Less Luminous Supergiant
            yerkes_type = "Less Luminous Supergiant"
        elif luminosity > 25:        
            yerkes_class = "II"        # Bright Giant
            yerkes_type = "Bright Giant"
        elif luminosity > 5:           
            yerkes_class = "III"       # Giant
            yerkes_type = "Giant"
        elif luminosity > 1.5:          
            yerkes_class = "IV"        # Subgiant
            yerkes_type = "Subgiant"
        elif luminosity > 0.08:        
            yerkes_class = "V"         # Dwarf (Main Sequence)
            yerkes_type = "Main Sequence"
        else:
            yerkes_class = "D"         # White Dwarf
            yerkes_type = "Dwarf"

        # Create type string based on spectral class, luminosity class, and color
        color_descriptions = {
            'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 
            'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'
        }
        star_type = f"{spectral_class}{subclass}{yerkes_class} {color_descriptions[spectral_class]} "
        star_type += yerkes_type + " Star"

        # Set the star's properties
        self.type = star_type
        self.radius = radius
        self.mass = mass
        self.temperature = temperature
        self.luminosity = luminosity_watts