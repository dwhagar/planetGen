# stellarObjects/constants.py

"""
Physical and Astronomical Constants
===================================

This module contains a collection of physical and astronomical constants used
for calculations within the planetGen package. All constants are in SI units
unless otherwise specified.
"""

# --- Physical Constants ---
EARTH_RADIUS_KM = 6371  # Earth's mean radius in kilometers
EARTH_GRAVITY = 9.807  # Standard Earth gravity in m/s^2
AU_TO_KM = 1.496e8  # Astronomical Unit to kilometers conversion factor
G = 6.6743e-11  # Gravitational constant in m^3/kg/s^2
R = 8.314  # Ideal gas constant in J/(mol·K)
BOLTZMANN = 1.381e-23  # Boltzmann constant in J/K
STEFAN_BOLTZMANN_CONSTANT = 5.67e-8  # Stefan-Boltzmann constant in W/m²/K⁴

# --- Astronomical Constants ---
SOLAR_MASS_TO_KG = 1.989e30  # Solar mass to kilograms conversion factor
SOLAR_LUMINOSITY = 3.82e26  # Solar luminosity in Watts

# --- Planet Generation Parameters ---

# Density ranges for terrestrial and gas giant planets in g/cm³
PLANET_DENSITY = {
    "t": (3.93, 5.51),  # Terrestrial: Range from Mars to Earth
    "g": (0.69, 1.64)   # Gas Giant: Range from Saturn to Neptune
}

# Atmospheric density ranges for terrestrial and gas giant planets in kg/m³
ATMOSPHERE_DENSITY = {
    "t": (0.02, 1.2),    # Terrestrial: Range from Mars to Venus
    "g": (0.69, 1.33),   # Gas Giant: Approximate range using Jupiter and Saturn's overall densities
}

# Atmospheric molar density ranges for terrestrial and gas giant planets in kg/mol
ATMOSPHERIC_MOLAR_DENSITY = {
    "t": (0.02897, 0.04347),  # Terrestrial: Range from Earth to Venus
    "g": (0.00226, 0.00416),   # Gas Giant: Range from Jupiter to Neptune
}

# The average ratio of a gas giant's core mass to its total mass
GAS_GIANT_CORE_ATMOSPHERE_RATIO = (0.03, 0.6)

# --- Planet Classification Data ---

# A dictionary defining the properties of different planet classes.
# Each class has a description, composition, radius range (in meters),
# habitable zone compatibility, atmosphere type, and planet type ('t' for terrestrial, 'g' for gas giant).
PLANET_CLASSES = {
    "A": {
        "description": "Small, barren, volcanic",
        "composition": "Igneous silica, basalt",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": "Sulfur dioxide, carbon dioxide",
        "type": "t"
    },
    "B": {
        "description": "Small, molten, thin atmosphere",
        "composition": "Iron, potassium, silicon",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": "Helium, sodium, oxygen",
        "type": "t"
    },
    "C": {
        "description": "Dead worlds, no atmosphere",
        "composition": "Anthracite, basalt, hydrocarbons",
        "radius_range": (500, 10000),
        "h": True, "e": True, "c": True,
        "atmosphere": None,
        "type": "t"
    },
    "D": {
        "description": "Small icy bodies (not true planets), often called dwarf or pseudo planets",
        "composition": "Frozen hydrocarbons, ice",
        "radius_range": (50, 3000),
        "h": True, "e": True, "c": True,
        "atmosphere": None,
        "type": "t"
    },
    "E": {
        "description": "Molten core & crust, thin atmosphere",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Hydrogen compounds",
        "type": "t"
    },
    "F": {
        "description": "Volcanic, shallow seas, bacterial life",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Carbon dioxide, ammonia, methane",
        "type": "t"
    },
    "G": {
        "description": "Rocky, barren, simple life",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Carbon dioxide, oxygen, nitrogen",
        "type": "t"
    },
    "H": {
        "description": "Desert worlds, minimal water, must be less than 10% liquid water",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Oxygen, nitrogen, argon, metals",
        "type": "t"
    },
    "I": {
        "description": "Ice giant, tilted magnetic field",
        "composition": "Rock, ice, methane, ammonia",
        "radius_range": (15000, 50000),
        "h": False, "e": False, "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "J": {
        "description": "Gas giant, turbulent atmosphere, rings",
        "composition": "Hydrogen, helium",
        "radius_range": (25000, 250000),
        "h": False, "e": False, "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "K": {
        "description": "Adaptable, thin atmosphere",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (2500, 7500),
        "h": False, "e": True, "c": False,
        "atmosphere": "Carbon dioxide",
        "type": "t"
    },
    "L": {
        "description": "Marginal, varied atmosphere, vegetation",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 7500),
        "h": False, "e": True, "c": False,
        "atmosphere": "Argon, oxygen, trace elements",
        "type": "t"
    },
    "M": {
        "description": "Terrestrial, Earth-like",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Oxygen, nitrogen, argon",
        "type": "t"
    },
    "N": {
        "description": "Reducing, hot, dense atmosphere",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Carbon dioxide, sulfides",
        "type": "t"
    },
    "O": {
        "description": "Pelagic, mostly water, ocean worlds with > 90% of the surface is liquid water",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Oxygen, nitrogen, argon",
        "type": "t"
    },
    "P": {
        "description": "Glaciated, cold",
        "composition": "Silicon, iron, magnesium, ice",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Oxygen, nitrogen, argon (thinning with age)",
        "type": "t"
    },
    "Q": {
        "description": "Eccentric orbit, extreme variations",
        "composition": "Silicon, iron, magnesium",
        "radius_range": (2000, 7500),
        "h": True, "e": True, "c": True,
        "atmosphere": "Varies widely (nitrogen, oxygen, argon, thin to dense)",
        "type": "t"
    },
    "R": {
        "description": "Ejected, geologically active",
        "composition": "Silicate compounds, iron",
        "radius_range": (7500, 10000),
        "h": False, "e": False, "c": False,
        "atmosphere": "Volcanic outgassing",
        "type": "t"
    },
    "S": {
        "description": "Supergiant, shield for inner planets",
        "composition": "Hydrogen, helium",
        "radius_range": (250000, 50000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "T": {
        "description": "Gas dwarf, thick atmosphere",
        "composition": "Hydrogen, helium, hydrocarbons",
        "radius_range": (250000, 25000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "Hydrogen, helium, hydrocarbons",
        "type": "g"
    },
    "U": {
        "description": "Ultragiant, can become a star",
        "composition": "Hydrogen, helium",
        "radius_range": (25000000, 60000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "Hydrogen, helium",
        "type": "g"
    },
    "V": {
        "description": "Super-Earth, high gravity",
        "composition": "Iron, iridium, tungsten, nickel",
        "radius_range": (10000, 15000),
        "h": False, "e": True, "c": False,
        "atmosphere": "Carbon dioxide, oxygen, hydrogen, helium",
        "type": "t"
    },
    "W": {
        "description": "Tidally locked, extreme variations",
        "composition": "Iron, potassium, silicon",
        "radius_range": (500, 10000),
        "h": True, "e": True, "c": False,
        "atmosphere": "Oxygen, sodium, hydrogen",
        "type": "t"
    },
    "X": {
        "description": "Stripped core, no atmosphere",
        "composition": "Molten iron",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": None,
        "type": "t"
    },
    "Y": {
        "description": "Demon class, toxic atmosphere",
        "composition": "Molten iron, sulfur, deuterium",
        "radius_range": (5000, 7500),
        "h": True, "e": False, "c": False,
        "atmosphere": "Turbulent, toxic, radiation",
        "type": "t"
    }
}

# Probabilities for each planet class to be generated
PLANET_CLASS_PROBABILITIES = {
    'A': 0.1399, 'B': 0.0722, 'C': 0.2365, 'D': 0.0142, 'E': 0.0239, 'F': 0.0335,
    'G': 0.0432, 'H': 0.0915, 'I': 0.0722, 'J': 0.0529, 'K': 0.0142, 'L': 0.0335,
    'M': 0.1345, 'N': 0.0239, 'O': 0.0045, 'P': 0.0046, 'Q': 0.0001, 'R': 0.0000,
    'S': 0.0001, 'T': 0.0001, 'U': 0.0001, 'V': 0.0045, 'W': 0.0001, 'X': 0.0002,
    'Y': 0.0002
}

# --- Star Generation Parameters ---

# Luminosity ranges for each spectral class in solar luminosities
LUMINOSITY_RANGES = {
    'O': (30000, 1000000), 'B': (25, 30000), 'A': (5, 25), 'F': (1.5, 5),
    'G': (0.6, 1.5), 'K': (0.08, 0.6), 'M': (0.00008, 0.08)
}

# Temperature ranges for each spectral class in Kelvin
TEMP_RANGES = {
    'O': (30000, 60000), 'B': (10000, 30000), 'A': (7500, 10000), 'F': (6000, 7500),
    'G': (5200, 6000), 'K': (3700, 5200), 'M': (2400, 3700)
}
