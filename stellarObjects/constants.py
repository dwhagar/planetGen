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
MILKY_WAY_MASS = 1.15e12 * SOLAR_MASS_TO_KG  # Mass of the Milky Way in kg
GALACTIC_CENTER_DISTANCE_LY = 25800  # Distance from Sol to the Galactic Center in light-years
LY_TO_M = 9.461e+15  # Light-year to meters conversion factor
AU_TO_M = 1.496e+11  # Astronomical Unit to meters conversion factor
ISM_PRESSURE = 2.5e-13  # Pressure of the local interstellar medium in Pascals (N/m^2)
SOLAR_RADIUS_M = 6.957e8  # Radius of the Sun in meters
SOLAR_ESCAPE_VELOCITY = 617.7 * 1000  # Sun's escape velocity in m/s
SOLAR_WIND_VELOCITY = 400 * 1000  # Average solar wind velocity in m/s
SOLAR_MASS_LOSS_RATE = 2e-14 * SOLAR_MASS_TO_KG * (365.25 * 24 * 3600)  # Sun's mass loss rate in kg/s

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
        "description": "a small, barren, and volcanic world",
        "composition": "igneous silica and basalt",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": "a mix of sulfur dioxide and carbon dioxide",
        "type": "t"
    },
    "B": {
        "description": "a small, molten world with a thin atmosphere",
        "composition": "iron, potassium, and silicon",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": "a mix of helium, sodium, and oxygen",
        "type": "t"
    },
    "C": {
        "description": "a dead world",
        "composition": "anthracite, basalt, and hydrocarbons",
        "radius_range": (500, 10000),
        "h": True, "e": True, "c": True,
        "atmosphere": None,
        "type": "t"
    },
    "D": {
        "description": "a small icy body",
        "composition": "frozen hydrocarbons and ice",
        "radius_range": (50, 3000),
        "h": True, "e": True, "c": True,
        "atmosphere": None,
        "type": "t"
    },
    "E": {
        "description": "a world with a molten core and crust, and a thin atmosphere",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "hydrogen compounds",
        "type": "t"
    },
    "F": {
        "description": "a volcanic world with shallow seas and bacterial life",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide, ammonia, and methane",
        "type": "t"
    },
    "G": {
        "description": "a rocky, barren world with simple life",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide, oxygen, and nitrogen",
        "type": "t"
    },
    "H": {
        "description": "a desert world with minimal water (less than 10% liquid water)",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, argon, and metals",
        "type": "t"
    },
    "I": {
        "description": "an ice giant with a tilted magnetic field",
        "composition": "rock, ice, methane, and ammonia",
        "radius_range": (15000, 50000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g"
    },
    "J": {
        "description": "a gas giant with a turbulent atmosphere and rings",
        "composition": "hydrogen and helium",
        "radius_range": (25000, 250000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g"
    },
    "K": {
        "description": "an adaptable world with a thin atmosphere",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (2500, 7500),
        "h": False, "e": True, "c": False,
        "atmosphere": "carbon dioxide",
        "type": "t"
    },
    "L": {
        "description": "a marginally habitable world with vegetation",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 7500),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of argon, oxygen, and trace elements",
        "type": "t"
    },
    "M": {
        "description": "a terrestrial Earth-like world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, and argon",
        "type": "t"
    },
    "N": {
        "description": "a hot world with a dense, reducing atmosphere",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide and sulfides",
        "type": "t"
    },
    "O": {
        "description": "a pelagic (ocean) world with greater than 90% of its surface covered in liquid water",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, and argon",
        "type": "t"
    },
    "P": {
        "description": "a cold, glaciated world",
        "composition": "silicon, iron, magnesium, and ice",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, and argon (thinning with age)",
        "type": "t"
    },
    "Q": {
        "description": "a world with an eccentric orbit and extreme temperature variations",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (2000, 7500),
        "h": True, "e": True, "c": True,
        "atmosphere": "a variable atmosphere (thin to dense) of nitrogen, oxygen, and argon",
        "type": "t"
    },
    "R": {
        "description": "an ejected, geologically active world",
        "composition": "silicate compounds and iron",
        "radius_range": (7500, 10000),
        "h": False, "e": False, "c": False,
        "atmosphere": "volcanic outgassing",
        "type": "t"
    },
    "S": {
        "description": "a supergiant that shields the inner planets",
        "composition": "hydrogen and helium",
        "radius_range": (250000, 50000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g"
    },
    "T": {
        "description": "a gas dwarf with a thick atmosphere",
        "composition": "hydrogen, helium, and hydrocarbons",
        "radius_range": (250000, 25000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen, helium, and hydrocarbons",
        "type": "g"
    },
    "U": {
        "description": "an ultragiant that could become a star",
        "composition": "hydrogen and helium",
        "radius_range": (25000000, 60000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g"
    },
    "V": {
        "description": "a Super-Earth with high gravity",
        "composition": "iron, iridium, tungsten, and nickel",
        "radius_range": (10000, 15000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide, oxygen, hydrogen, and helium",
        "type": "t"
    },
    "W": {
        "description": "a tidally locked world with extreme temperature variations",
        "composition": "iron, potassium, and silicon",
        "radius_range": (500, 10000),
        "h": True, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, sodium, and hydrogen",
        "type": "t"
    },
    "X": {
        "description": "a stripped core from a gas giant with no atmosphere",
        "composition": "molten iron",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": None,
        "type": "t"
    },
    "Y": {
        "description": "a 'demon' class world with a toxic atmosphere",
        "composition": "molten iron, sulfur, and deuterium",
        "radius_range": (5000, 7500),
        "h": True, "e": False, "c": False,
        "atmosphere": "a turbulent, toxic, and irradiated atmosphere",
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