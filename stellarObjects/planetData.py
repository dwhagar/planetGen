import math
import random
from .starData import to_scientific_notation

EARTH_RADIUS_KM = 6371
SOLAR_CONSTANT = 1361  # W/m²
AU_TO_KM = 1.496e8
# Gravitational constant in m^3/kg/s^2
G = 6.6743e-11
R = 8.31

planet_density = { # g per cubic centimeter
    "t": (3.93, 5.51), # Range from Mars to Earth
    "g": (0.69, 1.64)    # Range from Saturn to Neptune
}

atmosphere_density = { # kg per cubic meter
    "t": (0.020, 67),  # Range from Mars to Venus (kg/m³)
    "g": (0.69, 1.33),   # Approximate range using Jupiter and Saturn's overall densities
}

atmospheric_molar_density = { # kg per mol
    "t": (0.02897, 0.04347),  # Range from Earth to Venus (kg/mol)
    "g": (0.00226, 0.00416),   # Range from Jupiter to Neptune (kg/mol)
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
        self.atm_molar_density = None
        self.gravity = None
        self.atm_density = None
        self.surface_temperature = None
        self.density = None
        self.atmospheric_pressure = None
        self.mass = None
        self.atmosphere = None
        self.composition = None
        self.radius = radius
        self.planet_class = planet_class
        self.distance = distance

        # From the star, should not be changed.
        self.hab = hab_zone

        # Calculate additional properties.
        self.generate_planet()
        self.volume = (4/3) * math.pi * self.radius**3  # Calculate volume in km^3
        self.period = math.sqrt(self.distance**3)
        self.calculate_surface_gravity()
        self.calculate_atmospheric_pressure()
        self.calculate_atmospheric_temperature()

    def generate_planet(self):
        """
        Generates random planet properties (class, composition, atmosphere) based on distance and radius.
        """

        # Habitable Zone
        inner_bound = self.hab[0]  
        outer_bound = self.hab[1]

        # Determine zone based on distance
        if self.distance < inner_bound:
            zone = 'h'
        elif self.distance > outer_bound:
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
            self.radius = random.uniform(min_radius, max_radius)

        # If planet_class is present but radius is None, generate radius based on class
        elif self.planet_class is not None and self.radius is None:
            # Validate planet_class
            if self.planet_class not in planet_classes:
                raise ValueError("Invalid planet class")

            # Generate random radius within the allowed range for the chosen class
            min_radius, max_radius = planet_classes[self.planet_class]["radius_range"]
            self.radius = random.uniform(min_radius, max_radius)

        # If radius is present but planet_class is None, determine possible classes and choose randomly
        elif self.planet_class is None and self.radius is not None:
            # Filter possible planet classes based on zone and radius
            possible_classes = [
                c for c, data in planet_classes.items()
                if data[zone] and data["radius_range"][0] <= self.radius <= data["radius_range"][1]
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
            if not (min_radius <= self.radius <= max_radius):
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
        else:
            self.atmosphere = class_data["atmosphere"]
            # Approximate atmospheric density based on planet type (This can be refined later)
            min_a_density, max_a_density = atmosphere_density[planet_type]
            self.atm_density = random.uniform(min_a_density, max_a_density)
            min_am_density, max_am_density = atmospheric_molar_density[planet_type]
            self.atm_molar_density = random.uniform(min_am_density, max_am_density)

        # Recalculate mass based on the new density
        self.volume = (4/3) * math.pi * self.radius**3  # Calculate volume in km^3
        self.mass = self.volume * p_density * (10**-12)  # Update mass calculation
        self.density = p_density

    def calculate_surface_gravity(self):
        """
        Calculates the surface gravity of the planet in g's (Earth's gravity).
        """
        # Conversion from kilometers to meters
        radius_meters = self.radius * 1000

        # Calculate surface gravity using Newton's law of universal gravitation
        surface_gravity = (G * self.mass * 5.972e24) / (radius_meters ** 2)
        surface_gravity_g = surface_gravity / 9.80665  # Convert to g's

        self.gravity = surface_gravity_g

    def calculate_atmospheric_pressure(self):
        """
        Calculates the atmospheric conditions of the planet in bars.
        """

        if self.atmosphere == "None":
            self.atmospheric_pressure = 0.0
            return  # No atmosphere, no pressure

        # Calculate the orbital area (in square kilometers)
        orbital_radius_km = self.distance * AU_TO_KM
        orbital_area = 4 * math.pi * orbital_radius_km**2

        # Calculate atmospheric mass (in kilograms) using scale height (approximate for simplicity)
        solar_input_no_atm = SOLAR_CONSTANT * (1 / self.distance ** 2)
        scale_height = (R * solar_input_no_atm * 1000) / (self.atm_molar_density * self.gravity * 9.81)
        atmosphere_volume = (4/3) * math.pi * scale_height**3
        planet_volume = (4/3) * math.pi * self.radius**3
        atm_mass_kg = (atmosphere_volume - planet_volume) * self.atm_density * 10**9  # Convert density to kg/km^3

        # Calculate atmospheric pressure (in bars)
        # This is a simplification, as atmospheric pressure is complex and depends on many factors
        atmospheric_pressure = (atm_mass_kg * self.gravity) / orbital_area

        self.atmospheric_pressure = atmospheric_pressure

    def calculate_atmospheric_temperature(self):
        """
        Calculates the solar input at the planet's surface in Watts per square meter (W/m²),
        taking into account atmospheric pressure.
        """
        # Calculate solar input without atmosphere (inverse square law)
        solar_input_no_atm = SOLAR_CONSTANT * (1 / self.distance**2)

        # Calculate effective solar input at the surface (adjusted for atmospheric pressure)
        if self.atmosphere == "None":
            solar_input = solar_input_no_atm
        else:
            # Calculate atmospheric pressure in Pascals (Pa)
            atm_pressure_pa = self.atmospheric_pressure * 1e5  # 1 bar = 1e5 Pa

            # Apply atmospheric attenuation factor (simplified approximation)
            # Note: This is a very rough estimate, actual attenuation is more complex
            attenuation_factor = 1 - (0.2 * math.log10(atm_pressure_pa))
            solar_input = solar_input_no_atm * attenuation_factor

        # Calculate surface radiation in W/m^2 adjusted for atmospheric pressure
        surface_temperature = (solar_input / (5.67e-8))**(1/4)

        self.surface_temperature = surface_temperature

    def __str__(self):
        """
        Returns the wiki template text for this object.
        """
        output = [
            "{{Planet Data",
            f"|class={self.planet_class}",
            f"|distance={round(self.distance, 3)} AU",
            f"|period={years_to_time_string(self.period)}",
            f"|radius={round(self.radius, 1)} km",
            f"|gravity={round(self.gravity, 2)} G",
            "}}",
        ]

        # Check if atmosphere exists before adding information
        if self.atmosphere != "None":
            output.append(
                f"Atmospheric Conditions are an average of {self.atmospheric_pressure * 100:.1f} kPa and an average surface temperature of {self.surface_temperature - 273.15:.1f} degrees C.")
            output.append(self.atmosphere)
        else:
            output.append(f"There is no atmosphere and the surface has an average temperature of {self.surface_temperature - 273.15:.1f} degrees C.")

        if not output[-1] == self.composition:
            output.append(self.composition)

        return '\n'.join(output)
