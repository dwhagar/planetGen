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

# Average local stellar density in the solar neighborhood, in systems per
# cubic light-year. Cross-checked from two independent real surveys (see
# `spaceSector`'s module docstring for the worked conversions):
#   - Mamajek's stellar-density review (~0.0984 stars/pc^3;
#     https://www.pas.rochester.edu/~emamajek/memo_star_dens.html), converted
#     via 1 pc^3 = 34.706 ly^3, gives ~0.00284 stars/ly^3.
#   - RECONS' 10-parsec census (~270-414 objects within 10 pc, i.e. within
#     145,120 ly^3; https://chview.nova.org/solcom/stars/pc10.htm) gives a
#     consistent ~0.0019-0.0029 objects/ly^3.
LOCAL_STELLAR_DENSITY_LY3 = 0.00284

# --- Real Physical Density Ranges (used for planet/atmosphere generation) ---

# Density ranges for terrestrial and gas giant planets in g/cm³
PLANET_DENSITY = {
    "t": (3.93, 5.51),  # Terrestrial: Range from Mars to Earth
    "g": (0.69, 1.64)   # Gas Giant: Range from Saturn to Neptune
}

# Atmospheric density ranges for terrestrial and gas giant planets in kg/m³
ATMOSPHERE_DENSITY = {
    # Terrestrial: real surface air density spans Mars (~0.02 kg/m^3) to
    # Earth (~1.225 kg/m^3) here -- NOT Venus, whose real surface air
    # density (~65 kg/m^3) is roughly 50x this range's own upper bound; the
    # comment previously named Venus, which was incorrect and made this
    # range look far more permissive than it actually is (see
    # docs/analysis/habitability-atmosphere-sanity-review.md for how this
    # undershoot shows up in generated Class M output).
    "t": (0.02, 1.2),
    # Gas Giant: the previous (0.69, 1.33) range was Saturn's and Jupiter's
    # *bulk/overall* density (the same g/cm^3-scale figures used correctly
    # in PLANET_DENSITY above) -- the average density of the whole planet
    # including its dense core, not the density of the gas at the 1-bar
    # reference "surface" level this table is meant to describe. Derived
    # here instead from the ideal gas law (rho = P*M/(R*T), P = 1 bar) using
    # each planet's real 1-bar temperature and mean atmospheric molar mass:
    # Jupiter ~0.17, Saturn ~0.19, Uranus ~0.42, Neptune ~0.44 kg/m^3.
    "g": (0.16, 0.45),   # Gas Giant: Range from Jupiter to Neptune at the 1-bar reference level
}

# Atmospheric molar density ranges for terrestrial and gas giant planets in kg/mol
ATMOSPHERIC_MOLAR_DENSITY = {
    "t": (0.02897, 0.04347),  # Terrestrial: Range from Earth to Venus
    "g": (0.00207, 0.00266),   # Gas Giant: Range from Saturn to Neptune's real mean molar mass (~2.07-2.66 g/mol)
}

# Bulk density of a gas giant's light, puffy H/He envelope, in g/cm^3 -- used
# only by planetPhysics.generate_planet_properties' core/envelope density
# blend (and its mirror, plausibility.theoretical_gravity_bounds_g), NOT the
# same quantity as ATMOSPHERE_DENSITY["g"] above (that one feeds the
# atmospheric-pressure/scale-height calculation, a different physical layer).
# Deliberately a separate constant: ATMOSPHERE_DENSITY["g"]'s values are
# kg/m^3 gas densities at the 1-bar reference level (~0.16-0.45), two to
# three orders of magnitude smaller than a g/cm^3 bulk density -- reusing
# them directly as a kg/m^3 input to the density blend (which divides by
# 1000 to convert to g/cm^3, correct for the "t" case but not this one)
# would understate the envelope density by ~1000x, which the harmonic-mean
# blend below is
# highly sensitive to (the smaller of the two blended densities dominates
# the result almost regardless of mass fraction) -- collapsing every gas
# giant's overall density to a near-zero, physically meaningless value
# regardless of PLANET_CLASSES' own (or a class's density_range override's)
# core density. Range grounded in real measured "puffy" gas giants: WASP-193b
# (~0.06 g/cm^3, the lowest confirmed bulk density known) up through a
# representative light-envelope ceiling comfortably below Saturn's own real
# 0.69 g/cm^3 (the lightest actual solar-system planet), so the envelope
# term stays legitimately the *lighter* of the two blended components
# without collapsing the result the way the old ~0.0007-0.0013 g/cm^3 draw
# did.
GAS_ENVELOPE_BULK_DENSITY = (0.06, 0.3)

# n=1 polytrope (Lane-Emden) central-pressure formula for gas giants, used
# by planetData's flavor text as a genuinely-derived interior estimate
# rather than a curve-fit.
#
# A polytrope assumes P = K*rho^Gamma; for Gamma=2 (polytropic index n=1),
# the Lane-Emden equation (1/xi^2) d/dxi(xi^2 dtheta/dxi) = -theta^n has the
# closed-form solution theta(xi) = sin(xi)/xi, with its first zero (the
# body's surface) at xi_1 = pi. Combined with the definitions of a
# polytrope's total mass and radius, that solution reduces to a single
# closed-form relation between mass, radius, and central pressure -- no
# separately-fit constant needed:
#
#     P_core = pi * G * M^2 / (8 * R^4)
#
# (planetData.py computes this directly rather than exposing it as a named
# constant here, since it needs per-planet M and R). This is exactly
# pi^2/3 (~3.29x) higher than the cruder uniform-density-sphere estimate
# (P_c = 3*G*M^2/(8*pi*R^4)), because a real self-gravitating body's density
# rises toward the center rather than staying uniform. n=1 is the classic,
# simplest polytrope with an exact analytic solution, and is a standard
# textbook approximation for cold, dense hydrogen/helium (Zapolsky &
# Salpeter 1969 used the closely related n=3/2 for fully degenerate
# matter). Checked against Jupiter's real numbers (M=1.898e27 kg,
# R=71,492 km): this formula alone gives ~3,620 GPa, within ~10% of the
# ~4,000 GPa commonly cited for Jupiter's actual core pressure, with no
# tuning at all.
#
# Still not a real equation of state: no molecular-to-metallic hydrogen
# phase transition, no core/envelope differentiation, no temperature
# dependence. A rigorous calculation needs tabulated real EOS data (e.g.
# Saumon-Chabrier-van Horn) and a two-point boundary-value integration of
# hydrostatic equilibrium, well beyond flavor-text scope. (For reference,
# real estimates of Jupiter's actual core pressure run close to 4,000 GPa.)
# See planetData.py's gas giant text generation for the fuller caveats.

# Adiabatic index (ratio of specific heats, Cp/Cv) for a hydrogen/helium
# atmosphere in the molecular (non-dissociated, non-metallic) regime, used
# only to extend the already-computed surface_temperature (the same
# equilibrium-temperature calculation used for every planet, terrestrial or
# gas giant; see planetPhysics.calculate_atmospheric_conditions) into a
# convective-adiabat temperature-with-depth estimate for planetData's gas
# giant flavor text: T(P) = T_ref * (P/P_ref)^((gamma-1)/gamma). Real gas
# giant atmospheres measure Gamma1 ~= 1.4-1.44 in this molecular regime;
# 1.4 (the textbook diatomic-ideal-gas value) is used as a representative
# constant rather than a per-planet fit.
ADIABATIC_INDEX_H2_HE = 1.4

# Approximate pressure (Pa) at which hydrogen stops behaving as a molecular
# gas/fluid and becomes a liquid metallic state, used as a "where does the
# gas end" milestone for the adiabatic estimate above. Real estimates for
# this transition vary by model and temperature, roughly 1-3 Mbar; 1 Mbar
# (this value) is a commonly cited round figure. Applying the same
# isothermal barometric law used for the shallow atmosphere all the way out
# to this pressure is itself an extrapolation well past where that law is
# strictly valid, so this is a narrative milestone, not a claim of
# precision.
HYDROGEN_METALLIZATION_PRESSURE_PA = 1e11

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

AU_PER_PARSEC = 206264.80625
"""
float: The IAU-defined parsec, in Astronomical Units (AU). Used only as the
conversion path for the database persistence layer's sector-position
storage (star_systems.position_x/y/z_mpc, sectors.edge_mpc -- see
stellarObjects/schema.sql) -- not used anywhere in generation/physics code,
which keeps its own native light-year units for sector geometry.
"""
AU_PER_MILLIPARSEC = AU_PER_PARSEC / 1000
"""
float: One milliparsec (1/1000 parsec), in Astronomical Units (AU).
"""

# --- Physics-Model Constants (heliosphere / stellar wind) ---
KM_TO_M_FACTOR = 1000
ESCAPE_VELOCITY_CONSTANT = 2
RADIUS_SOL_EXPONENT = 2
LUMINOSITY_SOL_EXPONENT = -0.5
MASS_SOL_EXPONENT = -1.0
SECONDS_PER_YEAR = 365.25 * 24 * 3600
HELIOPAUSE_RADIUS_DEFAULT_M = 0
HYPERGIANT_WIND_VELOCITY_FACTOR = 2.6
GIANT_WIND_VELOCITY_FACTOR = 0.3
SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR = 2e-14
# Evolved-star mass loss (giants, subgiants, bright giants, supergiants,
# hypergiants): Nieuwenhuijzen & de Jager (1990), Mdot = C * L^a * M^b * R^c
# (L, M, R in solar units, Mdot in Msun/yr) -- a single empirical fit to 247
# stars spanning the whole HR diagram. This one formula replaces what used
# to be three independently hand-fit tiers stitched together at Yerkes-class
# boundaries (a separate hypergiant power-law, plus Reimers' Law with two
# different eta values for "giants" vs "supergiants") -- exactly the
# structure that let the old hypergiant constant drift ~9 orders of
# magnitude out of calibration unnoticed. Verified against measured rates:
# gives ~factor-of-4 agreement for ordinary giants/supergiants, but
# systematically overestimates by ~10x above log(L/Lsun) > 5 -- the
# supergiant/hypergiant regime -- per Mauron & Josselin (2011)'s comparison
# against 40 Galactic red supergiants, hence the correction below.
NDJ_MASS_LOSS_COEFFICIENT = 9.63e-15
NDJ_LUMINOSITY_EXPONENT = 1.42
NDJ_MASS_EXPONENT = 0.16
NDJ_RADIUS_EXPONENT = 0.81
NDJ_HIGH_LUMINOSITY_THRESHOLD_LSUN = 1e5
NDJ_HIGH_LUMINOSITY_CORRECTION_FACTOR = 0.1
# Kept as its own tier rather than folded into the Nieuwenhuijzen & de Jager
# formula above: stellar-evolution comparisons (e.g. against Vink et al.
# 2000/2001 radiation-driven-wind models, the standard for hot dwarfs) show
# N&dJ overestimates main-sequence O-star mass loss by up to ~20x, since it
# was fit mostly to evolved/luminous stars. Calibrated instead against
# measured O/B dwarf rates (e.g. Krticka 2014; Vink et al. 2000), which put
# late-O dwarfs (~1e5 Lsun) around 3-7e-8 Msun/yr. The previous value
# (3.2e-20) was ~3-4 orders of magnitude too low, making hot dwarf winds
# weaker than a cool M dwarf's and understating their heliosphere size.
OB_DWARF_MASS_LOSS_RATE_FACTOR = 1.7e-16
OB_DWARF_MASS_LOSS_RATE_EXPONENT = 1.75
OB_DWARF_WIND_VELOCITY_FACTOR = 2.0
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

# Typical main-sequence mass ranges by spectral class (letter), in solar
# masses, per standard stellar-classification references. Unlike
# YERKES_MASS_CONSTRAINTS above (which spans every luminosity class, so its
# "V" entry alone covers this whole range), this is specifically dwarfs
# broken out by spectral letter -- used to scale a star's galactic Hill
# radius (see `spaceSector.hill_radius_ly`) across spectral types.
SPECTRAL_MASS_RANGES = {
    'O': (15, 90), 'B': (2, 16), 'A': (1.4, 2.1), 'F': (1.0, 1.4),
    'G': (0.89, 1.07), 'K': (0.6, 0.9), 'M': (0.08, 0.45)
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

# --- Orbital motion (planetPhysics.generate_orbital_motion_properties) ---
# Real solar-system planets sit close to the ecliptic (Mercury, the most
# tilted, is ~7 degrees) -- 10 degrees gives a little headroom without
# implying a genuinely different orbital-plane population. Moons range much
# wider in reality (close-in regular moons are typically near-equatorial,
# but captured/irregular moons can be tilted or even retrograde) -- 25
# degrees is a simple, larger-but-still-modest scatter that doesn't require
# separately modeling regular vs. irregular moon populations.
PLANET_ORBITAL_INCLINATION_MAX_DEG = 10.0
MOON_ORBITAL_INCLINATION_MAX_DEG = 25.0

# Axial rotation ("day length") ranges, by body_type, in hours. Terrestrial:
# real solar-system terrestrial bodies span Earth/Mars-fast (~24-25h) to
# Mercury/Venus-slow (~1400-5800h) -- capped at 1400h (Mercury's own,
# real 58.6-day rotation) rather than reaching all the way to Venus' still
# slower, retrograde ~5800h, so the range stays "slow" without a separate
# retrograde-rotation concept. Gas giant: real gas giants all spin fast
# (Jupiter ~9.9h, Saturn ~10.7h, Uranus ~17.2h, Neptune ~16.1h).
ROTATION_PERIOD_RANGE_HOURS = {
    "t": (10.0, 1400.0),
    "g": (8.0, 20.0),
}

# Tidal-locking (despinning) timescale, per the standard simplified
# formula (Murray & Dermott, "Solar System Dynamics"; the uniform-sphere
# moment of inertia I = (2/5) m R^2 folds Q/k2 into a single 2Q/(15*k2)
# factor -- see planetPhysics._tidal_locking_timescale_seconds):
#
#   t_lock = (2*Q / (15*k2)) * (omega0 * a^6 * m_moon) / (G * M_primary^2 * R_moon^3)
#
# Q (tidal dissipation factor) and k2 (Love number) aren't modeled per
# body -- real values are only known precisely for a handful of solar-
# system bodies -- so these are single representative values for a
# rocky/icy moon (order-of-magnitude figures spanning the Moon's real
# Q~30-100/k2~0.024 and similar estimates for other rocky/icy satellites).
# Verified directly against real examples: with a 10-hour candidate
# initial rotation period, this formula gives ~47 million years for the
# real Earth-Moon system (real estimates: tens of millions of years,
# consistent with the Moon's long-since-locked observed state), ~90 years
# for Mars/Deimos (consistent with Deimos being locked given its tiny
# size and short distance), and ~1 billion years for Saturn/Iapetus
# (consistent with real estimates of a billion-year-plus despinning time
# for that unusually slow case) -- three real systems spanning many
# orders of magnitude in outcome, all landing in the right ballpark.
MOON_TIDAL_DISSIPATION_Q = 100.0
MOON_TIDAL_LOVE_NUMBER_K2 = 0.03
