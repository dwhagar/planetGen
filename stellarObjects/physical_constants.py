# stellarObjects/physical_constants.py

"""
Physical, Mathematical, and Astronomical Constants
=====================================================

This module holds constants that are true independent of this program's
design: universal physical/mathematical constants, unit conversions, real
astronomical reference values (solar/galactic figures, real density and
temperature ranges, stellar classification ranges), and the empirical
coefficients used by the stellar-wind/heliosphere physics model. Nothing
here is a game-design or generation-tuning choice — see `program_constants`
for those. Nothing in this module has side effects; it is imported wherever
these values are needed.
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

# --- Real Physical Density Ranges (used for planet/atmosphere generation) ---

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

# --- Physics-Model Constants (heliosphere / stellar wind) ---
KM_TO_M_FACTOR = 1000
ESCAPE_VELOCITY_CONSTANT = 2
HYPERGIANT_MASS_LOSS_RATE_FACTOR = 10**-5.0
HYPERGIANT_MASS_LOSS_RATE_EXPONENT = 1.5
RADIUS_SOL_EXPONENT = 2
LUMINOSITY_SOL_EXPONENT = -0.5
MASS_SOL_EXPONENT = -1.0
SECONDS_PER_YEAR = 365.25 * 24 * 3600
HELIOPAUSE_RADIUS_DEFAULT_M = 0
HYPERGIANT_WIND_VELOCITY_FACTOR = 2.6
REIMERS_LAW_CONSTANT = 4e-13
SUPERGIANT_REIMERS_ETA = 2.0
GIANT_WIND_VELOCITY_FACTOR = 0.3
STANDARD_REIMERS_ETA = 1.0
SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR = 2e-14
MIN_MOMENTUM_FLUX = 1e-15
FOUR_PI = 4 * math.pi

# --- Real Stellar Classification Data ---
SUBCLASS_MAX_VALUE = 9  # Spectral subclasses run 0-9 by astronomical convention
MAIN_SEQUENCE_MASS_LUMINOSITY_EXPONENT = 3.5
SPECTRAL_CLASS_COLORS = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}

# --- White Dwarf Physics ---

# The Chandrasekhar limit for white dwarf mass in solar masses.
CHANDRASEKHAR_LIMIT_SOL = 1.44

# A base radius for a 1 solar mass white dwarf, in kilometers.
WHITE_DWARF_BASE_RADIUS_KM = 5800

HOT_WHITE_DWARF_MIN_MASS_SOL = 1.1
COOL_WHITE_DWARF_MAX_MASS_SOL = 1.0
WHITE_DWARF_MASS_RADIUS_EXPONENT = -1/3

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

CO2_BASE_MOLAR_DENSITY = 0.04345
"""
float: The real molar density of CO2 (kg/mol), used as a reference point in
atmospheric greenhouse-effect calculations.
"""
