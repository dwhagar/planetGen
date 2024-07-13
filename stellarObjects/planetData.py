import math
import random
from .starData import to_scientific_notation

EARTH_RADIUS_KM = 6371
EARTH_GRAVITY = 9.807 # in m/s^2
AU_TO_KM = 1.496e8
# Gravitational constant in m^3/kg/s^2
G = 6.6743e-11
R = 8.314
BOLTZMANN = 1.381e-23 # in J/K
STEFAN_BOLTZMANN = 5.67e-8 # in W/m^2 / K^4

planet_density = { # g per cubic centimeter
    "t": (3.93, 5.51), # Range from Mars to Earth
    "g": (0.69, 1.64)  # Range from Saturn to Neptune
}

atmosphere_density = { # kg per cubic meter
    "t": (0.02, 1.2),    # Range from Mars to Venus (kg/m³)
    "g": (0.69, 1.33),   # Approximate range using Jupiter and Saturn's overall densities
}

atmospheric_molar_density = { # kg per mol
    "t": (0.02897, 0.04347),  # Range from Earth to Venus (kg/mol)
    "g": (0.00226, 0.00416),   # Range from Jupiter to Neptune (kg/mol)
}

gas_giant_core_atmosphere_ratio = (0.03, 0.6) # average ratio of core to atmosphere

planet_classes = {
    "A": {
        "description": "Small, barren, volcanic",
        "composition": "Igneous silica, basalt",
        "radius_range": (500, 5000), # In meters
        "h": True,  # Hot Zone
        "e": False,  # Ecosphere
        "c": False,  # Cold Zone
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
        "description": "Small icy bodies (not true planets), often called dwarf or pseudo planets",
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
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Hydrogen compounds",
        "type": "t"
    },
    "F": {
        "description": "Volcanic, shallow seas, bacterial life",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Carbon dioxide, ammonia, methane",
        "type": "t"
    },
    "G": {
        "description": "Rocky, barren, simple life",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False,
        "e": True,
        "c": False,
        "atmosphere": "Carbon dioxide, oxygen, nitrogen",
        "type": "t"
    },
    "H": {
        "description": "Desert worlds, minimal water, must be less than 10% liquid water",
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
        "description": "Pelagic, mostly water, ocean worlds with > 90% of the surface is liquid water",
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
        "c": False,
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
        "h": False,  # Technically can't be in any zone due to ejection
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
        "h": False,
        "e": True,
        "c": False,
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

planet_class_probabilities = {
    'A': 0.1399,
    'B': 0.0722,
    'C': 0.2365,
    'D': 0.0142,
    'E': 0.0239,
    'F': 0.0335,
    'G': 0.0432,
    'H': 0.0915,
    'I': 0.0722,
    'J': 0.0529,
    'K': 0.0142,
    'L': 0.0335,
    'M': 0.0915,
    'N': 0.0239,
    'O': 0.0045,
    'P': 0.0046,
    'Q': 0.0001,
    'R': 0.0000,
    'S': 0.0001,
    'T': 0.0001,
    'U': 0.0001,
    'V': 0.0045,
    'W': 0.0001,
    'X': 0.0002,
    'Y': 0.0432
}

def get_planet_mass_ranges():
    """
    Calculates and returns a dictionary of valid mass ranges (in kg) for each planet class.

    Returns:
        dict: A dictionary where keys are planet classes and values are tuples of (min_mass, max_mass)
    """

    mass_ranges = {}
    for planet_class, data in planet_classes.items():
        min_radius, max_radius = data["radius_range"]
        planet_type = data["type"]

        min_density, max_density = planet_density[planet_type]  # g/cm^3

        # Convert density from g/cm³ to kg/m³ for mass calculation
        min_density *= 1000
        max_density *= 1000

        if planet_type == "t":  # Terrestrial planet
            min_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density
            max_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density

        else:  # Gas giant
            # Get atmosphere density
            min_atm_density, max_atm_density = atmosphere_density[planet_type]  # kg/m³

            # Core-to-atmosphere calculations (similar to generate_planet)
            min_core_ratio, max_core_ratio = gas_giant_core_atmosphere_ratio

            # Calculate core and atmosphere masses separately for min/max
            min_core_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_density * min_core_ratio
            max_core_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_density * max_core_ratio

            min_atm_mass = (4 / 3) * math.pi * (min_radius ** 3) * min_atm_density * (1 - min_core_ratio)
            max_atm_mass = (4 / 3) * math.pi * (max_radius ** 3) * max_atm_density * (1 - max_core_ratio)

            min_mass = min_core_mass + min_atm_mass
            max_mass = max_core_mass + max_atm_mass

        mass_ranges[planet_class] = (min_mass, max_mass)  # Store as tuple

    return mass_ranges

planet_mass_ranges  = get_planet_mass_ranges()

def years_to_time_string(years):
    """
    Converts a decimal number of years into a human-readable string
    like "x years y days z hours m minutes" (omitting any components with zero values).
    """
    total_minutes = round(years * 365.25 * 24 * 60)  # Approximate total minutes in a year

    # Calculate individual time components
    years = total_minutes // (365 * 24 * 60)
    remaining_minutes = total_minutes % (365 * 24 * 60)
    days = remaining_minutes // (24 * 60)
    remaining_minutes %= 24 * 60
    hours = remaining_minutes // 60
    minutes = remaining_minutes % 60

    time_parts = []
    if years > 0:
        time_parts.append(f"{years} year{'s' if years > 1 else ''}")
    if days > 0:
        time_parts.append(f"{days} day{'s' if days > 1 else ''}")
    if hours > 0:
        time_parts.append(f"{hours} hour{'s' if hours > 1 else ''}")
    if minutes > 0:
        time_parts.append(f"{minutes} minute{'s' if minutes > 1 else ''}")

    # Join the non-zero time parts with "and"
    if len(time_parts) > 1:
        time_parts[-1] = f"and {time_parts[-1]}"

    return " ".join(time_parts)

def calc_object_mass(object_class, object_radius, object_density = None):
    """
    Calculates the mass of an object in kg.

    @param object_class: Class of the object (String)
    @param object_radius: Radius of the object in km (float)
    @param object_density: Density of the object in g/cm^3 (float)
    @return: Tuple with the volume and the mass (volume, mass) in m^3 and kg
    """
    if object_density is None:
        min_density, max_density = planet_density[planet_classes[object_class]['type']]
        p_density = random.uniform(min_density, max_density)
    else:
        p_density = object_density

    volume = (4 / 3) * math.pi * (object_radius * 1000) ** 3  # Calculate volume in m^3
    mass = volume * p_density * 1000

    return volume, mass

class Asteroid_Belt:
    """
    A basic class to store information for an asteroid belt.
    """

    def __init__(self, distance, lower_limit, upper_limit):
        self.distance = distance
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.type = 'a'

    def __str__(self):
        return f"== Asteroid Belt ==\nAn asteroid belt orbits roughly between {self.lower_limit:.3f} AU and {self.upper_limit:.3f} AU."

class Planet:
    """
    A Class representing a single planet and all of its properties.
    """

    def __init__(self, hab_zone, distance, star_output, star_radius, star_temperature, star_mass,
                 radius=None, planet_class=None, mass=None, zone_override = None, distance_override = None,
                 is_moon = False, force_moons = False):
        """
        Initializes a Planet object with its radius and the spectral class of its host star.
        """
        self.is_moon = is_moon
        self.moons = []
        self.zone = None
        self.description = None
        self.star_radius = star_radius * 1000 # in meters
        self.star_output = star_output # in Watts
        self.star_temperature = star_temperature # in Kelvin
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

        # From the star, should not be changed.
        self.hab = hab_zone

        if self.is_moon:
            id_number = random.randint(100, 999)
            self.id_string = f"{id_number}"
            force_moons = False
        else:
            id_number = random.randint(1000, 9999)
            id_letters = "".join(random.choices("ABCDEFGHIJKLMNOPQRSTUVWXYZ", k=2))
            self.id_string = f"{id_number}-{id_letters}"

        # Calculate additional properties.
        self.generate_planet(zone_override)
        self.volume = (4/3) * math.pi * self.radius**3  # Calculate volume in km^3
        self.period = math.sqrt(self.distance**3)
        self.calculate_surface_gravity()
        self.calculate_atmospheric_conditions(distance_override)

        if (force_moons or random.randint(0,1) == 1) and not self.is_moon:
            self.generate_moons()

    def generate_planet(self, zone_override = None):
        """
        Generates random planet properties (class, composition, atmosphere, mass) based on distance and optional radius/class/mass inputs.
        """

        inner_bound = self.hab[0]
        outer_bound = self.hab[1]

        if self.distance < inner_bound:
            zone = 'h'
        elif self.distance > outer_bound:
            zone = 'c'
        else:
            zone = 'e'

        if zone_override:
            if zone_override.lower() in "hce":
                zone = zone_override.lower()

        self.zone = zone

        # --- Input Validation and Random Generation ---

        if self.planet_class is None and self.radius is None and self.mass is None:
            # Fully random generation (no inputs provided)
            valid_classes = [c for c, data in planet_classes.items() if data[zone]]
            classes = list(planet_class_probabilities.keys())
            probabilities = list(planet_class_probabilities.values())
            class_valid = False
            while not class_valid:
                class_choice = random.choices(classes, weights=probabilities, k=1)[0]
                if class_choice in valid_classes:
                    self.planet_class = class_choice
                    class_valid = True

            min_radius, max_radius = planet_classes[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        elif self.planet_class is not None and self.radius is None and self.mass is None:
            # Planet class given, generate radius
            self._validate_planet_class(zone)
            min_radius, max_radius = planet_classes[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        elif self.planet_class is None and self.radius is not None and self.mass is None:
            # Radius given, determine possible classes
            possible_classes = [
                c for c, data in planet_classes.items()
                if data[zone] and data["radius_range"][0] <= self.radius <= data["radius_range"][1]
            ]

            if not possible_classes:
                raise ValueError("No valid planet class for the given radius in this zone")

            self.planet_class = random.choice(possible_classes)
            self._validate_radius()

        elif self.planet_class is None and self.radius is None and self.mass is not None:
            # Mass given, determine possible classes and radii
            possible_classes = []

            for planet_class, class_data in planet_classes.items():
                min_mass, max_mass = planet_mass_ranges[planet_class]

                if min_mass <= self.mass <= max_mass and class_data[zone]:
                    possible_classes.append(planet_class)

            self._validate_mass()

        elif self.planet_class is not None and self.radius is not None and self.mass is None:
            # Class and radius given, validate
            self._validate_planet_class(zone)
            self._validate_radius()

        elif self.planet_class is not None and self.radius is None and self.mass is not None:
            # Class and mass given, validate and generate radius
            self._validate_planet_class(zone)
            self._validate_mass()
            min_radius, max_radius = planet_classes[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)
            self._validate_radius()

        elif self.planet_class is None and self.radius is not None and self.mass is not None:
            # Radius and mass given, validate and determine possible classes
            possible_classes = []

            for planet_class, class_data in planet_classes.items():
                min_mass, max_mass = planet_mass_ranges[planet_class]
                min_radius, max_radius = class_data["radius_range"]

                mass_valid = min_mass <= self.mass <= max_mass and class_data[zone]
                radius_valid = min_radius <= self.radius <= max_radius and class_data[zone]

                if mass_valid and radius_valid:
                    possible_classes.append(planet_class)

            if not possible_classes:
                raise ValueError("No valid planet class for the given radius / mass in this zone")

            self.planet_class = random.choice(possible_classes)
            self._validate_radius()
            self._validate_mass()

        else:
            # All inputs provided, fully validate
            self._validate_planet_class(zone)
            self._validate_radius()
            self._validate_mass()

            # Get the properties for the chosen planet class
        class_data = planet_classes[self.planet_class]

        # Update planet properties based on class data
        self.composition = class_data["composition"]
        self.description = class_data["description"]

        # Determine planet type (Terrestrial or Gas Giant) based on class data
        planet_type = class_data["type"]  # Assuming the planet type is "t" or "g" in class_data
        self.type = planet_type

        # Generate random density based on planet type
        min_density, max_density = planet_density[planet_type]
        self.density = random.uniform(min_density, max_density)

        # Handle atmosphere (set to "None" if it's None in class_data)
        if class_data["atmosphere"] is None:
            self.atmosphere = "None"
        else:
            self.atmosphere = class_data["atmosphere"]

            if self.planet_class == 'N': # Class N worlds are basically life genus.
                self.atm_density = 65
                min_am_density, max_am_density = atmospheric_molar_density[planet_type]
                self.atm_molar_density =  max_am_density
            else:
                # Approximate atmospheric density based on planet type (This can be refined later)
                min_a_density, max_a_density = atmosphere_density[planet_type]
                self.atm_density = random.uniform(min_a_density, max_a_density)
                min_am_density, max_am_density = atmospheric_molar_density[planet_type]
                self.atm_molar_density = random.uniform(min_am_density, max_am_density)

        self.volume, self.mass = calc_object_mass(self.planet_class, self.radius, self.density)

        if self.type == 'g':
             core_to_atmosphere_ratio = random.uniform(*gas_giant_core_atmosphere_ratio)
             self.density = self.density * core_to_atmosphere_ratio + (core_to_atmosphere_ratio - 1) * (self.atm_density / 1000)

        self.volume, self.mass = calc_object_mass(self.planet_class, self.radius, self.density)

        self.hill_radius = (self.distance * AU_TO_KM) * (self.mass / (3 * self.star_mass)) ** (1 / 3) # in km!
        self.min_orbit_distance = (5 * self.hill_radius) / AU_TO_KM # Now back to AU

    def calculate_surface_gravity(self):
        """
        Calculates the surface gravity of the planet in g's (Earth's gravity).
        """
        # Conversion from kilometers to meters
        radius_meters = self.radius * 1000

        # Calculate surface gravity using Newton's law of universal gravitation
        surface_gravity = (G * self.mass) / (radius_meters ** 2)
        surface_gravity_g = surface_gravity / EARTH_GRAVITY  # Convert to g's

        if surface_gravity_g <= 0:
            raise ValueError('Invalid value for gravity.')

        self.gravity = surface_gravity_g

    def calculate_atmospheric_conditions(self, distance_override = None):
        """
        Calculates the atmospheric conditions of the planet in bars.
        """
        # This will give us the option to override the distance to generate atmospheric conditions for a moon
        # based on the orbital distance of the planet.  This value will be the OFFSET for this object's distance
        # and the distance you want it to be.  This function does not know the difference between a moon and a planet
        # and will treat it like a planet of the appropriate class.
        if not distance_override is None:
            distance = float(distance_override)
        else:
            distance = float(self.distance)

        planet_volume = (4 * math.pi * self.radius ** 3) / 3
        orbital_radius_km = distance * AU_TO_KM
        output_area = 4 * math.pi * orbital_radius_km ** 2     # Output of the star at orbit in square kilometers
        # Start using watts per square meter here, as that is what the equations expect.
        solar_output_at_orbit = (self.star_output / output_area) / 1e6   # in Watts per square meter <----
        albedo = random.uniform(0.12, 0.35)
        surface_temperature_no_atmosphere = ((1 - albedo) * solar_output_at_orbit / (4 * STEFAN_BOLTZMANN)) ** (1 / 4)

        # Calculate effective solar input at the surface (adjusted for atmospheric pressure)\
        if self.atmosphere == "None":
            self.surface_temperature = surface_temperature_no_atmosphere  # This figure does not require an atmosphere.
            self.atmospheric_pressure = 0.0
        else:
            # Calculate the Pressure
            scale_height = \
                (R * surface_temperature_no_atmosphere) / (self.atm_molar_density * self.gravity * EARTH_GRAVITY)
            scale_factor = random.uniform(5, 7) # How many scale heights to calculate for
            self.scale_height = scale_height

            atmosphere_thickness = scale_height * random.uniform(5, 7)
            atmosphere_volume = (4 * math.pi * (self.radius + atmosphere_thickness) ** 3) / 3 - planet_volume

            atmospheric_mass = 0
            for zone in range(round(scale_factor)):
                zone_volume = atmosphere_volume + planet_volume - (4 * math.pi * (self.radius + (zone * scale_height)) ** 3) / 3
                if zone < 1:
                    zone_density = self.atm_density
                else:
                    zone_density = self.atm_density / (zone * 2.7)
                atmospheric_mass += zone_volume * zone_density
            atmospheric_force = atmospheric_mass * (self.gravity * EARTH_GRAVITY)  # Force in Newtons (N)
            planet_surface_area = 4 * math.pi * (self.radius * 1000) ** 2  # Surface area in square meters (m^2)
            atmospheric_pressure = (atmospheric_force / planet_surface_area) * 7500  # Pressure in Pascals (Pa)

            # Assuming a linear relationship between CO₂ abundance and greenhouse it causes
            CO2_BASE_MOLAR_DENSITY = 0.04345  # kg/mol (approximate for Mars)
            CO2_MAX_GREENHOUSE_FACTOR = 5  # Venus's greenhouse factor is much higher than Earth's

            greenhouse_factor =\
                abs((self.atm_molar_density - CO2_BASE_MOLAR_DENSITY) / CO2_BASE_MOLAR_DENSITY * CO2_MAX_GREENHOUSE_FACTOR)
            surface_temperature_atmosphere =\
                ((1 - albedo) * solar_output_at_orbit * (1 + greenhouse_factor) / (4 * STEFAN_BOLTZMANN)) ** (1 / 4)

            self.surface_temperature = surface_temperature_atmosphere
            self.atmospheric_pressure = atmospheric_pressure

    def generate_moons(self):
        """
        Generates a system of moons for the given planet.  This function does not deal with probability.
        """

        # There is a cubic relationship between radius (volume ha radius cubed) and mass.
        moon_blacklist = ['Q', 'R', 'V', 'W', 'X', 'W', 'X', 'Y']
        max_moon_mass = self.mass / 10
        max_moon_radius = self.radius / 10 ** (1 / 3) # Here the reduction factor is 1/10 of the mass

        possible_classes = []

        for c, data in planet_mass_ranges.items():
            if planet_classes[c][self.zone] and planet_classes[c]["type"] == 't' and not c in moon_blacklist:
                if data[1] <= max_moon_mass and planet_classes[c]['radius_range'][1] <= max_moon_radius:
                    possible_classes.append(c)

        if not possible_classes:
            return

        # Leaving these as kilometers to avoid precision errors for something so close to 0 as lunar distances in AU
        if self.scale_height:
            low_orbit = self.scale_height * 15
        else:
            low_orbit = 100
        high_orbit = self.min_orbit_distance * AU_TO_KM
        total_orbit_distance = low_orbit

        while True:
            if ((total_orbit_distance >= high_orbit
                 or total_orbit_distance >= (self.distance * AU_TO_KM))
                    or total_orbit_distance >= (self.min_orbit_distance * AU_TO_KM)):
                break

            # This is about choosing the distance, which must in AU as that's what everything else works in.
            # last_orbital_min is going to be the minimum orbital distance of the last object and should ensure that
            # no two orbits overlap.

            classes = list(planet_class_probabilities.keys())
            probabilities = list(planet_class_probabilities.values())
            class_valid = False
            moon_class = None
            while not class_valid:
                moon_class = random.choices(classes, weights=probabilities, k=1)[0]
                if moon_class in possible_classes:
                    class_valid = True

            if max_moon_radius > planet_classes[moon_class]['radius_range'][1]:
                radius_limit = planet_classes[moon_class]['radius_range'][1]
            else:
                radius_limit = max_moon_radius

            # This will be the distance of the moon from the planet.
            moon_distance = random.uniform(total_orbit_distance, high_orbit) / AU_TO_KM # Needs to be in AU
            moon_radius = random.uniform(planet_classes[moon_class]['radius_range'][0], radius_limit)

            new_moon = Planet(self.hab, moon_distance,
                              self.star_output, self.star_radius,
                              self.star_temperature, self.star_mass,
                              radius=moon_radius, planet_class=moon_class,
                              zone_override=self.zone, distance_override=self.distance, is_moon=True)
            self.moons.append(new_moon)

            total_orbit_distance = (new_moon.distance + new_moon.min_orbit_distance) * AU_TO_KM

        return

    def _validate_planet_class(self, zone):
        if self.planet_class not in planet_classes or not planet_classes[self.planet_class][zone]:
            raise ValueError("Invalid planet class for this zone")

    def _validate_radius(self):
        radius_range = planet_classes[self.planet_class]["radius_range"]
        valid = (radius_range[0] <= self.radius <= radius_range[1])

        if not valid:
            raise ValueError("Invalid radius for planet class")

    def _validate_mass(self):
        mass_range = planet_mass_ranges[self.planet_class]
        valid = (mass_range[0] <= self.mass <= mass_range[1])

        if not valid:
            raise ValueError("Invalid mass for planet class")

    def __str__(self):
        """
        Returns the wiki template text for this object.
        """
        if self.is_moon:
            distance_text = f"|distance={to_scientific_notation(self.distance * AU_TO_KM, 5)} km"
            header_level = '==='
        else:
            if self.distance < 1:
                distance_text = f"|distance={to_scientific_notation(self.distance * AU_TO_KM, 1)} km ({round(self.distance, 3)} AU)"
            else:
                distance_text = f"|distance={round(self.distance, 3)} AU"
            header_level = '=='

        random_id = f"{header_level} {self.id_string} {header_level}"

        output = [
            random_id,
            "{{Planet Data",
            f"|class={self.planet_class}",
            distance_text,
            f"|period={years_to_time_string(self.period)}",
            f"|radius={round(self.radius, 1):,} km",
            f"|gravity={round(self.gravity, 2)} g",
            "}}",
        ]

        # Check if atmosphere exists before adding information
        if len(self.moons) > 1 and not self.is_moon:
            output.append(f"There are {len(self.moons)} moons orbiting this planet.")
        elif len(self.moons) == 1 and not self.is_moon:
            output.append("There is 1 moon orbiting this planet.")
        elif not self.is_moon:
            output.append("There are no moons orbiting this planet.")

        if self.type == "t":
            if self.atmosphere != "None":
                output.append(
                    f"Surface conditions are an average of {self.atmospheric_pressure / 1000:.1f} kPa or {self.atmospheric_pressure / 101300:.2f} atmospheres and an average surface temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
                output.append(self.atmosphere)
            else:
                output.append(f"There is no atmosphere and the surface has an average temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
        else:
            output.append(f"The average internal conditions of this gas giant are an average of {self.atmospheric_pressure / 1000:.1f} kPa or {self.atmospheric_pressure / 101300:.1f} atmospheres and {self.surface_temperature - 273.15:.1f} degrees C.")
            output.append(f"The atmospheric pressure drops by a third to a half for every {self.scale_height / 1000:.1f} km from the core.")

        output.append(self.description)

        if not output[-1] == self.composition:
            output.append(self.composition)

        if not self.moons is None and not self.is_moon:
            for moon in self.moons:
                output.append('\n' + str(moon))

        return '\n'.join(output)
