# stellarObjects/constants.py

"""
Physical and Astronomical Constants
===================================

This module contains a collection of physical and astronomical constants used
for calculations within the planetGen package. All constants are in SI units
unless otherwise specified.
"""

import math

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

# --- Star System Generation Parameters ---
INITIAL_PLANET_DISTANCE_FACTOR = 0.55
ASTEROID_BELT_PROBABILITY = 0.1
ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN = 1.1
ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX = 2
BASE_MAX_SYSTEM_OBJECTS = 15
ABSOLUTE_MAX_SYSTEM_OBJECTS = 500
MIN_ASTEROID_BELT_SEPARATION = 0.05

# Renamed AU_TO_LIGHT_YEAR to LY_TO_AU for clarity and consistency.
LY_TO_AU = 63241.1
"""
float: Conversion factor from Light-Years (LY) to Astronomical Units (AU).
1 Light-Year is approximately 63241.1 AU.
"""
AU_TO_LY = 1 / LY_TO_AU
"""
float: Conversion factor from Astronomical Units (AU) to Light-Years (LY).
1 AU is approximately 1/63241.1 Light-Years.
"""

# Asteroid Belt Configuration
ASTEROID_COMPONENTS = [
    "carbon", "silicon", "magnesium", "aluminum", "calcium",
    "sulfur", "phosphorus", "iron", "nickel",
    "iridium", "palladium", "platinum", "gold", "osmium", "ruthenium",
    "rhodium", "olivine", "pyroxene", "plagioclase feldspars", "kamacite",
    "taenite", "troilite", "schreibersite", "cohenite", "serpentine",
    "magnetite", "hematite", "chromite", "silicon carbide",
    # Rare earth element compounds
    "bastnasite", "monazite", "xenotime", "cerite", "gadolinite", "samarskite",
    "fergusonite", "euxenite",
    # Platinum-group metal compounds
    "osmiridium", "sperrylite", "cooperite", "braggite",
    "laurite", "vysotskite",
    # Naturally occurring radioactive material compounds
    "uraninite", "thorianite", "carnotite", "autunite", "brannerite",
    "torbernite", "coffinite",
    # Titanium and titanium compounds
    "titanium", "rutile", "ilmenite", "titanite", "perovskite",
    "brookite", "anatase",
    # Rock, crystal, and gem compounds
    "silicon dioxide", # Quartz, amethyst, citrine
    "aluminum oxide", # Corundum, ruby, sapphire
    "beryllium aluminum cyclosilicate", # Beryl, emerald, aquamarine
    "aluminum silicate fluoride hydroxide", # Topaz
    "magnesium aluminum oxide", # Spinel
    "zirconium silicate", # Zircon
    "borosilicate", # Tourmaline
    "potassium aluminum silicate", # Orthoclase
    "calcium carbonate", # Calcite, aragonite
    "sodium chloride", # Halite
    "calcium fluoride", # Fluorite
    "calcium fluorophosphate", # Apatite
    "hydrated copper aluminum phosphate", # Turquoise
    "copper carbonate hydroxide", # Malachite
    "almandine", # Iron aluminum silicate (garnet)
    "pyrope", # Magnesium aluminum silicate (garnet)
    "spessartine", # Manganese aluminum silicate (garnet)
    "grossular", # Calcium aluminum silicate (garnet)
    "muscovite", # Hydrated potassium aluminum silicate (mica)
    "biotite", # Iron magnesium potassium aluminum silicate (mica)
    "chrysoberyl", # Beryllium aluminum oxide
    # Solid Hydrocarbons (Polycyclic Aromatic Hydrocarbons)
    "naphthalene", "anthracene", "phenanthrene", "pyrene",
    "coronene", "fluoranthene"
]
"""
list: A comprehensive list of various components that can be found in asteroids.
These components are used to generate the composition of asteroid belts.
"""

LY_THRESHOLD = 1.0
"""
float: The distance in Light-Years (LY) at which the display format for
asteroid belt distances switches from AU to LY for better readability.
"""

# --- Star Calculation Parameters ---
HABITABLE_ZONE_BUFFER_AU = 0.2
HELIOSPHERE_DISPLAY_THRESHOLD_LY = 0.1
KM_TO_M_FACTOR = 1000
ESCAPE_VELOCITY_CONSTANT = 2
HYPERGIANT_MASS_LOSS_RATE_FACTOR = 10**-5.0
HYPERGIANT_MASS_LOSS_RATE_EXPONENT = 1.5
RADIUS_SOL_EXPONENT = 2
LUMINOSITY_SOL_EXPONENT = -0.5
MASS_SOL_EXPONENT = -1.0
SECONDS_PER_YEAR = 365.25 * 24 * 3600
HELIOPAUSE_RADIUS_DEFAULT_M = 0
ROUND_HABITABLE_ZONE_AU = 2
ROUND_HABITABLE_ZONE_AU_SMALL = 5
ROUND_RADIUS_KM = 2
SCIENTIFIC_NOTATION_DECIMAL_PLACES = 2
ROUND_TEMPERATURE_NEAREST_HUNDRED = -2

# --- Star Generation Parameters ---
MIN_INITIAL_STAR_AGE_GY = 0.1
MAX_INITIAL_STAR_AGE_LIFESPAN_RATIO = 0.9
OLD_STAR_AGE_LIFESPAN_RATIO = 0.5
YOUNG_STAR_AGE_LIFESPAN_RATIO = 1/3
MIN_PLANET_AGE_ADJUSTMENT_FACTOR = 0.8
MAX_PLANET_AGE_ADJUSTMENT_FACTOR = 0.95
WHITE_DWARF_AGE_ADDITION_GY = 5
HYPERGIANT_WIND_VELOCITY_FACTOR = 2.6
REIMERS_LAW_CONSTANT = 4e-13
SUPERGIANT_REIMERS_ETA = 2.0
GIANT_WIND_VELOCITY_FACTOR = 0.3
STANDARD_REIMERS_ETA = 1.0
SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR = 2e-14
MIN_MOMENTUM_FLUX = 1e-15
FOUR_PI = 4 * math.pi
PERCENT_SOL_THRESHOLD_LOW = 0.01
PERCENT_SOL_THRESHOLD_HIGH = 2
RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD = 100000
PERCENT_MULTIPLIER = 100
SUBCLASS_MAX_VALUE = 9
MAIN_SEQUENCE_MASS_LUMINOSITY_EXPONENT = 3.5
HOT_WHITE_DWARF_MIN_MASS_SOL = 1.1
COOL_WHITE_DWARF_MAX_MASS_SOL = 1.0
WHITE_DWARF_MASS_RADIUS_EXPONENT = -1/3
SPECTRAL_CLASS_COLORS = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}
SPECTRAL_PROBABILITIES_ABSURD = {'O': 100, 'B': 0, 'A': 0, 'F': 0, 'G': 0, 'K': 0, 'M': 0}
SPECTRAL_PROBABILITIES_LARGE_STAR = {'O': 10, 'B': 20, 'A': 30, 'F': 30, 'G': 10, 'K': 0, 'M': 0}
SPECTRAL_PROBABILITIES_NORMAL = {'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45}

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
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.005),
            "normal": (0.0, 0.5),
            "slow": (0.0, 3.0)
        }
    },
    "B": {
        "description": "a small, molten world with a thin atmosphere",
        "composition": "iron, potassium, and silicon",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": "a mix of helium, sodium, and oxygen",
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.005),
            "normal": (0.0, 0.5),
            "slow": (0.0, 3.0)
        }
    },
    "C": {
        "description": "a dead world",
        "composition": "anthracite, basalt, and hydrocarbons",
        "radius_range": (500, 10000),
        "h": True, "e": True, "c": True,
        "atmosphere": None,
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "D": {
        "description": "a small icy body",
        "composition": "frozen hydrocarbons and ice",
        "radius_range": (50, 500),
        "h": True, "e": True, "c": True,
        "atmosphere": None,
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "E": {
        "description": "a world with a molten core and crust, and a thin atmosphere",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "hydrogen compounds",
        "type": "t",
        "life_chemical": ["Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.005, 0.015),
            "normal": (0.5, 1.5),
            "slow": (3.0, 8.0)
        }
    },
    "F": {
        "description": "a volcanic world with shallow seas and bacterial life",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide, ammonia, and methane",
        "type": "t",
        "life_chemical": ["Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.005, 0.015),
            "normal": (0.5, 1.5),
            "slow": (3.0, 8.0)
        }
    },
    "G": {
        "description": "a rocky, barren world with simple life",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide, oxygen, and nitrogen",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.015, 0.05),
            "normal": (1.5, 4.0),
            "slow": (8.0, 25.0)
        }
    },
    "H": {
        "description": "a desert world with minimal water (less than 10% liquid water)",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, argon, and metals",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "I": {
        "description": "an ice giant with a tilted magnetic field",
        "composition": "rock, ice, methane, and ammonia",
        "radius_range": (15000, 50000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "J": {
        "description": "a gas giant with a turbulent atmosphere and rings",
        "composition": "hydrogen and helium",
        "radius_range": (25000, 250000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "K": {
        "description": "an adaptable world with a thin atmosphere",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (2500, 7500),
        "h": False, "e": True, "c": False,
        "atmosphere": "carbon dioxide",
        "type": "t",
        "life_chemical": ["Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.015, 0.05),
            "normal": (1.5, 4.0),
            "slow": (8.0, 25.0)
        }
    },
    "L": {
        "description": "a marginally habitable world with vegetation",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 7500),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of argon, oxygen, and trace elements",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "M": {
        "description": "a terrestrial Earth-like world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, and argon",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "N": {
        "description": "a hot world with a dense, reducing atmosphere",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide and sulfides",
        "type": "t",
        "life_chemical": ["Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.005, 0.015),
            "normal": (0.5, 1.5),
            "slow": (3.0, 8.0)
        }
    },
    "O": {
        "description": "a pelagic (ocean) world with greater than 90% of its surface covered in liquid water",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, and argon",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "P": {
        "description": "a cold, glaciated world",
        "composition": "silicon, iron, magnesium, and ice",
        "radius_range": (5000, 10000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, nitrogen, and argon (thinning with age)",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.08, 0.1),
            "normal": (8.0, 10.0),
            "slow": (45.0, 100.0)
        }
    },
    "Q": {
        "description": "a world with an eccentric orbit and extreme temperature variations",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (2000, 7500),
        "h": True, "e": True, "c": True,
        "atmosphere": "a variable atmosphere (thin to dense) of nitrogen, oxygen, and argon",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "R": {
        "description": "an ejected, geologically active world",
        "composition": "silicate compounds and iron",
        "radius_range": (7500, 10000),
        "h": False, "e": False, "c": False,
        "atmosphere": "volcanic outgassing",
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "S": {
        "description": "a supergiant that shields the inner planets",
        "composition": "hydrogen and helium",
        "radius_range": (250000, 50000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "T": {
        "description": "a gas dwarf with a thick atmosphere",
        "composition": "hydrogen, helium, and hydrocarbons",
        "radius_range": (250000, 25000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen, helium, and hydrocarbons",
        "type": "g",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "U": {
        "description": "an ultragiant that could become a star",
        "composition": "hydrogen and helium",
        "radius_range": (25000000, 60000000),
        "h": False, "e": False, "c": True,
        "atmosphere": "a mix of hydrogen and helium",
        "type": "g",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "V": {
        "description": "a Super-Earth with high gravity",
        "composition": "iron, iridium, tungsten, and nickel",
        "radius_range": (10000, 15000),
        "h": False, "e": True, "c": False,
        "atmosphere": "a mix of carbon dioxide, oxygen, hydrogen, and helium",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.015, 0.1),
            "normal": (1.5, 9.0),
            "slow": (8.0, 50.0)
        }
    },
    "W": {
        "description": "a tidally locked world with extreme temperature variations",
        "composition": "iron, potassium, and silicon",
        "radius_range": (500, 10000),
        "h": True, "e": True, "c": False,
        "atmosphere": "a mix of oxygen, sodium, and hydrogen",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.015, 0.1),
            "normal": (1.5, 9.0),
            "slow": (8.0, 50.0)
        }
    },
    "X": {
        "description": "a stripped core from a gas giant with no atmosphere",
        "composition": "molten iron",
        "radius_range": (500, 5000),
        "h": True, "e": False, "c": False,
        "atmosphere": None,
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.01, 0.1),
            "normal": (1.0, 10.0),
            "slow": (10.0, 100.0)
        }
    },
    "Y": {
        "description": "a 'demon' class world with a toxic atmosphere",
        "composition": "molten iron, sulfur, and deuterium",
        "radius_range": (5000, 7500),
        "h": True, "e": False, "c": False,
        "atmosphere": "a turbulent, toxic, and irradiated atmosphere",
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.005),
            "normal": (0.0, 0.5),
            "slow": (0.0, 3.0)
        }
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

# The Chandrasekhar limit for white dwarf mass in solar masses.
CHANDRASEKHAR_LIMIT_SOL = 1.44

# A base radius for a 1 solar mass white dwarf, in kilometers.
WHITE_DWARF_BASE_RADIUS_KM = 5800

# Temperature ranges for each spectral class in Kelvin
TEMP_RANGES = {
    'O': (30000, 60000), 'B': (10000, 30000), 'A': (7500, 10000), 'F': (6000, 7500),
    'G': (5200, 6000), 'K': (3700, 5200), 'M': (2400, 3700)
}

# Luminosity ranges based on a star's spectral class (color/temperature).
# This is primarily useful for Main Sequence (V) stars, where the relationship is strong.
SPECTRAL_LUMINOSITY_RANGES = {
    'O': (30000, 1000000), 'B': (25, 30000), 'A': (5, 25), 'F': (1.5, 5),
    'G': (0.6, 1.5), 'K': (0.08, 0.6), 'M': (0.0001, 0.08)
}

# Approximate luminosity ranges for each Yerkes luminosity class, in solar luminosities.
# This is the primary determinant of a star's energy output.
YERKES_LUMINOSITY_RANGES = {
    "0": (500000, 2000000),      # Hypergiants
    "IA+": (500000, 2000000),    # Luminous Supergiants (upper end)
    "IA": (50000, 500000),       # Luminous Supergiants
    "IAB": (10000, 100000),      # Intermediate-size Luminous Supergiants
    "IB": (1000, 50000),         # Less Luminous Supergiants
    "II": (100, 10000),          # Bright Giants
    "III": (50, 1000),           # Normal Giants
    "IV": (2, 50),               # Subgiants
    "V": (0.0001, 1000000),      # Main-sequence (very broad)
    "VI": (0.00001, 0.1),        # Subdwarfs (dimmer than main sequence)
    "VII": (0.0001, 0.1),        # White Dwarfs (standard)
    "D": (0.0001, 0.1)           # Alias for White Dwarfs
}

# Physically allowed mass ranges for each Yerkes luminosity class, in solar masses.
# This enforces realistic constraints on generated stars.
YERKES_MASS_CONSTRAINTS = {
    "0": (20, 150), "IA+": (20, 150), "IA": (10, 40), "IAB": (8, 20), "IB": (8, 20),
    "II": (2, 15),
    "III": (0.8, 8),
    "IV": (1, 5),
    "V": (0.08, 150),
    "VI": (0.1, 0.8),
    "VII": (0.5, CHANDRASEKHAR_LIMIT_SOL),
    "D": (0.5, CHANDRASEKHAR_LIMIT_SOL)
}

# --- Life and Photosynthesis Data ---

LIFE_CHEMICALS = {
    "Retinal": {
        "description": "a primitive photosensitive chemical that represents a young/early biosphere precursor, operating via proton-motive gradients rather than electron transport chains in pre-oxygenated, reducing atmospheres",
        "evolutionary_time_scale": "fast",
        "absorption_spectrum": ["Green", "Yellow", "500-650 nm"],
        "reflection_spectrum_visible": ["Blue", "Red", "Purple", "Magenta"],
        "reflection_spectrum_non_visible": ["Green edge (inverse of the standard Vegetation Red Edge)"],
        "star_spectra_probabilities": {
            "O": 10, "B": 10, "A": 15, "F": 20, "G": 35, "K": 35, "M": 20
        }
    },
    "Melanin": {
        "description": "a radiotrophic chemical that thrives in environments with lethal ionizing radiation or stellar flaring, and is viable on rogue planets subject to cosmic ray bombardment without a thick protective atmosphere",
        "evolutionary_time_scale": "normal",
        "absorption_spectrum": ["Extreme Ultraviolet (EUV)", "X-rays", "Gamma rays", "Ionizing radiation"],
        "reflection_spectrum_visible": ["None (Charcoal-black)"],
        "reflection_spectrum_non_visible": ["Extremely low albedo across all bands"],
        "star_spectra_probabilities": {
            "O": 90, "B": 85, "A": 60, "F": 30, "G": 5, "K": 5, "M": 75
        }
    },
    "Chlorophyll a": {
        "description": "a standard porphyrin that drives oxygenic photosynthesis, requiring an oxygenated atmosphere as a climax ecology and stable, moderate radiation environments",
        "evolutionary_time_scale": "slow",
        "absorption_spectrum": ["Red (644-746 nm)", "Blue (468-476 nm)"],
        "reflection_spectrum_visible": ["Green"],
        "reflection_spectrum_non_visible": ["Near-infrared edge (Canonical Vegetation Red Edge / VRE)"],
        "star_spectra_probabilities": {
            "O": 0, "B": 0, "A": 5, "F": 25, "G": 90, "K": 60, "M": 5
        }
    },
    "Blue-Optimized Porphyrins": {
        "description": "a porphyrin that requires simultaneous evolution of biofluorescence for UV shielding, down-converting lethal high-energy light into safe visible wavelengths and creating transient, global optical fluorescence aligned with stellar flares",
        "evolutionary_time_scale": "slow",
        "absorption_spectrum": ["Intense Blue (468-476 nm)", "Ultraviolet (UV)"],
        "reflection_spectrum_visible": ["Red", "Orange", "Yellow", "Biofluorescent Green", "Biofluorescent Red"],
        "reflection_spectrum_non_visible": ["Blue-shifted vegetation edge", "Strong UV absorption"],
        "star_spectra_probabilities": {
            "O": 5, "B": 10, "A": 40, "F": 85, "G": 15, "K": 0, "M": 0
        }
    },
    "Bacteriochlorophylls": {
        "description": "a chemical that drives anoxygenic photosynthesis, does not split water, and requires chemical electron donors like Hydrogen Sulfide (H2S), Ferrous Iron (Fe2+), or Molecular Hydrogen (H2), resulting in reducing atmospheres lacking Oxygen (O2) or Ozone (O3)",
        "evolutionary_time_scale": "slow",
        "absorption_spectrum": ["Broad Visible Spectrum", "Near-Infrared (805-890 nm)", "Deep Infrared (987-1050 nm)"],
        "reflection_spectrum_visible": ["None (Black)"],
        "reflection_spectrum_non_visible": ["Deep infrared edge (> 1.0 µm)"],
        "star_spectra_probabilities": {
            "O": 0, "B": 0, "A": 0, "F": 5, "G": 10, "K": 40, "M": 95
        }
    },
    "Zinc-Bacteriochlorophyll": {
        "description": "a secondary adaptation derived from standard magnesium-chelatase pathways that requires highly acidic global oceans (pH 1.5 - 3.0) and planetary heavy metal abundance and mobilization, often tied to severe volcanic outgassing",
        "evolutionary_time_scale": "slow",
        "absorption_spectrum": ["Broad Visible Spectrum", "Near-Infrared (793 nm, 853 nm)"],
        "reflection_spectrum_visible": ["None (Black)"],
        "reflection_spectrum_non_visible": ["Blue-shifted near-infrared edge"],
        "star_spectra_probabilities": {
            "O": 5, "B": 5, "A": 10, "F": 10, "G": 10, "K": 20, "M": 40
        }
    }
}

STAR_EVOLUTION = {
    "O": {
        "lifespan_gy": (0.001, 0.01), # 1 Million - 10 Million years
        "supported_evolutionary_scales": ["fast"],
        "potentially_viable_chemicals": ["Retinal"],
        "evolutionary_constraint_notes": "has an extremely short lifespan, as the star exhausts its fuel and goes supernova before planets can sufficiently cool to form stable liquid oceans, meaning only the absolute fastest, most primitive precursor photochemistry could theoretically occur"
    },
    "B": {
        "lifespan_gy": (0.01, 0.1), # 10 Million - 100 Million years
        "supported_evolutionary_scales": ["fast"],
        "potentially_viable_chemicals": ["Retinal"],
        "evolutionary_constraint_notes": "is highly volatile and short-lived, providing barely enough time for planetary cooling and the emergence of single-celled life utilizing simple proton-motive gradients"
    },
    "A": {
        "lifespan_gy": (0.1, 1.0), # 100 Million - 1 Billion years
        "supported_evolutionary_scales": ["fast", "normal"],
        "potentially_viable_chemicals": ["Retinal", "Melanin"],
        "evolutionary_constraint_notes": "allows enough time for the development of early biospheres and radiotrophic organisms adapted to high-energy radiation, but likely dies before complex, slow-evolving porphyrin-based photosynthesis can be perfected"
    },
    "F": {
        "lifespan_gy": (2.0, 4.0), # 2 Billion - 4 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Blue-Optimized Porphyrins (+ GFPs)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "possesses a main-sequence habitable zone lifespan spanning between 2 and 4 billion years, which provides sufficient time for complex biospheres to emerge, including advanced UV-shielding and biofluorescent adaptations"
    },
    "G": {
        "lifespan_gy": (8.0, 12.0), # 8 Billion - 12 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Chlorophyll a (Standard Porphyrins)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "is the solar standard, providing stable, long-term conditions ideal for the slow evolution of highly complex, oxygenic photosynthesis based on the tetrapyrrole chlorin ring"
    },
    "K": {
        "lifespan_gy": (15.0, 40.0), # 15 Billion - 40 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Chlorophyll a (Standard Porphyrins)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "is exceptionally stable and long-lived, offering tens of billions of years for slow-evolving biospheres to reach climax ecologies and adapt to slight red-shifts in the stellar spectrum"
    },
    "M": {
        "lifespan_gy": (100.0, 1000.0), # 100 Billion - 1 Trillion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Bacteriochlorophylls (BChls)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "is the longest-lived star in the universe, and while early flaring requires rapid adaptation (e.g., radiotrophic melanin) during the initial hundreds of millions of years, their trillions of years of stability allow for deep-infrared anoxygenic photosynthesis to dominate permanently"
    }
}

EVOLUTIONARY_TIMELINES = {
    "normal": { # Changed from "norm" to "normal" for consistency
        "star_lifespan": 10.0, # Billion Years
        "evolutionary_pace": "Standard",
        "abiogenesis": 0.5, # Billion Years
        "photosynthesis": 1.5, # Billion Years
        "complex_cells": 2.5, # Billion Years
        "multicellularity": 4.0, # Billion Years
        "technological_civilization": 4.5 # Billion Years
    },
    "fast": {
        "star_lifespan": 0.1, # Billion Years (100 Million Years)
        "evolutionary_pace": "Hyper-Accelerated",
        "abiogenesis": 0.005, # Billion Years (5 Million Years)
        "photosynthesis": 0.015, # Billion Years (15 Million Years)
        "complex_cells": 0.030, # Billion Years (30 Million Years)
        "multicellularity": 0.050, # Billion Years (50 Million Years)
        "technological_civilization": 0.080 # Billion Years (80 Million Years)
    },
    "slow": {
        "star_lifespan": 1000.0, # Billion Years (1 Trillion Years)
        "evolutionary_pace": "Decelerated",
        "abiogenesis": 3.0, # Billion Years
        "photosynthesis": 8.0, # Billion Years
        "complex_cells": 15.0, # Billion Years
        "multicellularity": 25.0, # Billion Years
        "technological_civilization": 35.0 # Billion Years
    }
}

EVOLUTIONARY_TEXT = {
    "abiogenesis": (
        "Life is in its absolute infancy, consisting entirely of simple, single-celled organisms "
        "without a nucleus. There is no true flora or fauna. The biosphere is largely limited to "
        "chemotrophs and extremophiles thriving in nutrient-rich primordial pools or around "
        "hydrothermal vents, feeding on chemical reactions."
    ),
    "photosynthesis": (
        "The biosphere is dominated by simple, light-harvesting microbes. While macroscopic flora "
        "and fauna still do not exist, vast microbial mats and stony biological structures (similar "
        "to stromatolites) line the shallow waters and coastlines. These organisms are actively "
        "terraforming the planet by releasing oxygen into the atmosphere."
    ),
    "complex_cells": (
        "Microscopic life has evolved distinct internal structures and nuclei (eukaryotes). This stage "
        "introduces the first single-celled protozoa-analogues (early microscopic 'fauna') and "
        "planktonic autotrophs (early microscopic 'flora'). The land remains barren, but the oceans "
        "teem with a complex microscopic food web of tiny predators and prey."
    ),
    "multicellularity": (
        "Cells have cooperated to form complex, macroscopic organisms, marking the true arrival of "
        "flora and fauna. The oceans are populated by macro-algae, sponges, and early invertebrates "
        "like arthropod or jellyfish analogues. On land, pioneering flora—such as moss, ferns, and "
        "primitive vascular plants—have taken root, supporting emerging land-dwelling animals."
    ),
    "technological_civilization": (
        "The planet boasts a fully mature, highly complex biosphere featuring diverse ecosystems, "
        "topped by at least one sapient, tool-using species. Complex flora forms sprawling forests "
        "and varied biomes, while highly evolved fauna occupy vast ecological niches. The sapient "
        "population actively reshapes the world through agriculture, architecture, and industry."
    )
}

# --- Planet Generation Specific Constants ---
HABITABLE_PLANET_CLASSES = ['E', 'F', 'G', 'H', 'K', 'L', 'M', 'O', 'P', 'V', 'W']
"""
list: A list of planet class codes that are considered habitable.
"""

MOON_BLACKLIST = ['Q', 'R', 'V', 'W', 'X', 'Y']
"""
list: A list of planet class codes that cannot be generated as moons.
"""

CO2_BASE_MOLAR_DENSITY = 0.04345
"""
float: The base molar density of CO2 used in atmospheric calculations.
"""

CO2_MAX_GREENHOUSE_FACTOR = 5
"""
int: The maximum greenhouse effect factor for CO2 in atmospheric calculations.
"""