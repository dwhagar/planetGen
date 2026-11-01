# planetgen/physics/constants.py

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

# Values that are measured or defined are written out in SI floats, from the
# [R4361] and IAU 2015 tables (GEN.66, GEN.211): astropy supplied them
# until GEN.211 dropped it at run time (it cost 280 ms and 29 MB per process).
# The astropy constants test holds each against astropy's own value.

# --- Physical Constants ---
EARTH_RADIUS_KM = 6371  # Earth's mean radius in kilometers
EARTH_GRAVITY = 9.80665  # Standard Earth gravity in m/s^2 (9.80665)
# --- Distance ladder (Boss's exact values, in meters) ---
# Every displayed distance goes through `format.format_distance_m`, which picks
# the largest of these units the value is at least 1 of: km < AU < mpc <
# cpc < ly < pc < kpc < Mpc < Gpc. Every other length conversion below is
# derived from these, so none of them is rounded on its own.
KM_M = 1e3
AU_M = 149_597_870_700  # IAU 2012, exact
LIGHTYEAR_M = 9_460_730_472_580_800  # c x 365.25 days, exact
PARSEC_M = 3.085677581491367e+16  # 648000/pi AU
MILLIPARSEC_M = PARSEC_M * 1e-3
CENTIPARSEC_M = PARSEC_M * 1e-2
KILOPARSEC_M = PARSEC_M * 1e3
MEGAPARSEC_M = PARSEC_M * 1e6
GIGAPARSEC_M = PARSEC_M * 1e9
AU_TO_KM = AU_M / KM_M  # Astronomical Unit to kilometers conversion factor
G = 6.6743e-11  # Gravitational constant in m^3/kg/s^2
SPEED_OF_LIGHT_M_S = 299_792_458  # Speed of light in a vacuum, m/s, exact (the SI definition; LIGHTYEAR_M uses the same c) -- used for Schwarzschild radius (phenomena.compact_remnant.BlackHole)
CELSIUS_ZERO_K = 273.15  # 0 °C in K (format.format_temperature_k)
STANDARD_ATMOSPHERE_PA = 101_325.0  # 1 atm, exact, Pa (format.format_pressure_pa)
INCH_M = 0.0254  # The international inch, m (exact; the Customary unit preference)
FOOT_M = 0.3048  # The international foot, m (exact; the Customary unit preference)
MILE_M = 1_609.344  # The international mile, m (exact)
POUND_KG = 0.45359237  # The avoirdupois pound, kg (exact)
PSI_PA = 6_894.757293168361  # 1 lbf/in², Pa (format.format_pressure_pa)
SPEED_OF_LIGHT_KMS = SPEED_OF_LIGHT_M_S / KM_M  # Exact (the SI definition), km/s -- the speed ladder's c (format.format_speed_kms, static/unitladder.js)
R = 8.31446261815324  # Ideal gas constant in J/(mol·K)
BOLTZMANN = 1.380649e-23  # Boltzmann constant in J/K
REDUCED_PLANCK = 1.0545718176461565e-34  # hbar, J*s -- Hawking radiation (phenomena.compact_remnant)
STEFAN_BOLTZMANN_CONSTANT = 5.6703744191844314e-08  # Stefan-Boltzmann constant in W/m²/K⁴
COSMIC_BACKGROUND_TEMPERATURE_K = 2.725  # Cosmic microwave background, K -- the floor for any body's surface temperature

# --- Astronomical Constants ---
SOLAR_MASS_TO_KG = 1.988409870698051e+30  # Solar mass to kilograms conversion factor
SOLAR_LUMINOSITY = 3.828e+26  # Solar luminosity in Watts
MIN_TIDAL_DISTANCE_LY = 4.89  # 1.5 pc, a point's mean distance in the 2 pc core ball: the nearest to the center a system's Hill sphere is worked out at (at 0 it would be 0, leaving no room for a planet)
EARTH_MASS_TO_KG = 5.972167867791379e+24  # Earth mass to kilograms conversion factor
MOON_MASS_KG = 7.342e+22  # The Moon's mass, kg: project-nominal (the IAU defines no lunar mass), for the mass ladder's lunar rung
JUPITER_MASS_TO_KG = 1.8981245973360505e+27  # Jupiter mass to kilograms conversion factor
JUPITER_RADIUS_KM = 71_492.0  # Jupiter's mean equatorial radius in kilometers

# Minimum Mass Solar Nebula [R1854] *solid* surface density at
# 1 AU, and its power-law falloff exponent: Sigma_s(a) = SOL_GCM2 *
# (a/1AU)^EXPONENT, valid inside the snow line. This is the *dust/ice*
# surface density (what actually accretes into planetesimals/embryos), not
# the much larger *gas* surface density Hayashi's model is often quoted
# with (~1700 g/cm^2 at 1 AU) -- using the gas value here would overstate
# every isolation mass below by roughly (1700/7)^1.5 ~= 3800x. See
# `formation.mmsn_surface_density_gcm2`/`formation.isolation_mass_kg`, used by
# `planet_mass.oligarch_core_kg`.
MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2 = 7.0
MMSN_SURFACE_DENSITY_EXPONENT = -1.5

# Rotation-curve model for a star's circular orbit around the galactic
# center (see `galactic_orbit.calculate_galactic_orbit`): v(r) = FLAT_VELOCITY * r /
# sqrt(r^2 + CORE_RADIUS^2) -- a pseudo-isothermal-halo-style curve that
# rises ~linearly (near-solid-body) inside the core radius and flattens to
# the asymptotic velocity well outside it, the same qualitative shape real
# spiral-galaxy rotation curves show (steep bulge rise, flat disk plateau)
# and the standard explanation for why "flat rotation curve" is cited as
# dark-matter evidence in the first place. At the Milky Way's 8.128 kpc (`mathcheck.SUN_REFERENCE_RADIUS_PC`, a test
# reference only) it gives 206 km/s and a 242 million-year orbit: slower and longer than the measured ~230 km/s and
# ~203 to 220 million years (the Sun's angular speed from the Sgr A* proper motion, and the circular period of the
# potential model, docs/design/galaxy/galactic-potential.md section 2). A stand-in that the potential's circular speed replaces
# (GEN.115), not a calibration to the real Sun.
GALACTIC_ROTATION_FLAT_VELOCITY_KMS = 220.0  # Asymptotic circular velocity, km/s -- the IAU (1985) standard local-standard-of-rest value, still within the range of modern estimates (~220-236 km/s)
GALACTIC_ROTATION_CORE_RADIUS_PC = 3000  # Turnover radius of the rotation curve (rise -> plateau), in parsecs -- approximates the real Milky Way's bulge/inner-disk transition

LY_TO_M = LIGHTYEAR_M  # Light-year to meters conversion factor
AU_TO_M = AU_M  # Astronomical Unit to meters conversion factor
ISM_PRESSURE = 2.5e-13  # Pressure of the local interstellar medium in Pascals (N/m^2)
# Inside a nebula or supernova remnant the cloud's own pressure replaces
# ISM_PRESSURE at the heliopause: its ram pressure (rho * v^2, with the
# star moving through the cloud) plus its thermal pressure (n k T). See
# star.compressed_heliosphere_radius.
HYDROGEN_ATOM_MASS_KG = 1.67353286206015e-27  # Mass of a hydrogen atom in kg
ISM_MASS_PER_HYDROGEN = 1.4  # Gas mass per hydrogen atom, in hydrogen masses (helium adds ~40%)
STAR_CLOUD_RELATIVE_SPEED_MS = 26.0 * 1000  # A star's typical speed through its surrounding gas, m/s (the Sun's through the local cloud)
CM3_TO_M3 = 1e6  # Particles per cm^3 to particles per m^3
SOLAR_RADIUS_M = 695_700_000.0  # Radius of the Sun in meters
SOLAR_ESCAPE_VELOCITY = math.sqrt(2 * G * SOLAR_MASS_TO_KG / SOLAR_RADIUS_M)  # Sun's escape velocity in m/s (about 617.5 km/s)
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

# Giant-planet mass-radius relation (GEN.34), the shape of Chen & Kipping
# (2017, ApJ 834:17): below about 0.4 Jupiter masses a volatile-rich giant
# grows as R ~ M^0.52 (their "Neptunian" regime, anchored here on Saturn);
# above it electron degeneracy holds the radius nearly flat, R ~ M^-0.04
# (their "Jovian" regime, anchored on Jupiter), up to the deuterium-burning
# limit. The two meet at about 140 Earth masses. This gives Neptune ~1.6,
# Saturn ~0.7, Jupiter ~1.3 and a 10 Jupiter-mass super-Jupiter ~16 g/cm^3,
# where drawing a radius and a density separately (the old way) gave a
# median of 0.25 g/cm^3 and no super-Jupiters at all. Each anchor is
# (mass in Earth masses, radius in Earth radii, exponent); the scatter is
# the relation's intrinsic fractional spread in radius at fixed mass
# (theirs is about 15% and 7%; tightened a little here), cut at 3 sigma.
GIANT_NEPTUNIAN_MASS_RADIUS = (95.16, 9.449, 0.519)
GIANT_JOVIAN_MASS_RADIUS = (317.8, 11.21, -0.04)
GIANT_RADIUS_SCATTER = {"neptunian": 0.10, "jovian": 0.07}

# n=1 polytrope (Lane-Emden) central-pressure formula for gas giants, used
# by planet's flavor text as a genuinely-derived interior estimate
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
# (generation/planet.py computes this directly rather than exposing it as a named
# constant here, since it needs per-planet M and R). This is exactly
# pi^2/3 (~3.29x) higher than the cruder uniform-density-sphere estimate
# (P_c = 3*G*M^2/(8*pi*R^4)), because a real self-gravitating body's density
# rises toward the center rather than staying uniform. n=1 is the classic,
# simplest polytrope with an exact analytic solution, and is a standard
# textbook approximation for cold, dense hydrogen/helium ([R4867]
# used the closely related n=3/2 for fully degenerate
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
# See generation/planet.py's gas giant text generation for the fuller caveats.

# Adiabatic index (ratio of specific heats, Cp/Cv) for a hydrogen/helium
# atmosphere in the molecular (non-dissociated, non-metallic) regime, used
# only to extend the already-computed surface_temperature (the same
# equilibrium-temperature calculation used for every planet, terrestrial or
# gas giant; see physics.planets.calculate_atmospheric_conditions) into a
# convective-adiabat temperature-with-depth estimate for planet's gas
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
LY_TO_AU = LIGHTYEAR_M / AU_M
"""
float: Conversion factor from Light-Years (LY) to Astronomical Units (AU).
1 Light-Year is about 63241.08 AU.
"""
AU_TO_LY = 1 / LY_TO_AU
"""
float: Conversion factor from Astronomical Units (AU) to Light-Years (LY).
1 AU is about 1/63241.08 Light-Years.
"""

AU_PER_PARSEC = PARSEC_M / AU_M
"""
float: The IAU-defined parsec, in Astronomical Units (AU). Used only as the
conversion path for the database persistence layer's sector-position
storage (star_systems.position_x/y/z_mpc, sectors.edge_mpc -- see
planetgen/db/schema.sql) -- not used anywhere in generation/physics code,
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
SECONDS_PER_YEAR = 31_557_600.0  # the Julian year, 365.25 days
ISOTROPIC_TOTAL_OVER_ONE_AXIS = math.sqrt(3.0)  # an isotropic velocity dispersion's three-axis total over its one-axis value
HELIOPAUSE_RADIUS_DEFAULT_M = 0
HYPERGIANT_WIND_VELOCITY_FACTOR = 2.6
GIANT_WIND_VELOCITY_FACTOR = 0.3
SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR = 2e-14
# Evolved-star mass loss (giants, subgiants, bright giants, supergiants,
# hypergiants): the empirical fit, Mdot = C * L^a * M^b * R^c
# (L, M, R in solar units, Mdot in Msun/yr) -- a single empirical fit to 247
# stars spanning the whole HR diagram. This one formula replaces what used
# to be three independently hand-fit tiers stitched together at Yerkes-class
# boundaries (a separate hypergiant power-law, plus Reimers' Law with two
# different eta values for "giants" vs "supergiants") -- exactly the
# structure that let the old hypergiant constant drift ~9 orders of
# magnitude out of calibration unnoticed. Verified against measured rates:
# gives ~factor-of-4 agreement for ordinary giants/supergiants, but
# systematically overestimates by ~10x above log(L/Lsun) > 5 -- the
# supergiant/hypergiant regime -- per [R2958]'s comparison
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
# measured O/B dwarf rates (e.g. [R4557]), which put
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
SPECTRAL_CLASS_COLORS = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}

# --- White Dwarf Physics ---

# The Chandrasekhar limit for white dwarf mass in solar masses.
CHANDRASEKHAR_LIMIT_SOL = 1.44

# A base radius for a 1 solar mass white dwarf, in kilometers.
WHITE_DWARF_BASE_RADIUS_KM = 5800

# A white dwarf cooler than this is a cool one (the hot ones are the cores of planetary nebulae and the youngest remnants).
COOL_WHITE_DWARF_BELOW_K = 10000
WHITE_DWARF_MASS_RADIUS_EXPONENT = -1/3

# The temperatures between which the nearest spectral type is of each letter, in Kelvin: the midpoints between the last type of one
# letter and the first of the next in the Pecaut and Mamajek dwarf sequence (GEN.943; `mamajek_sequence.letter_range_k` reads them
# from the file and a test keeps these equal). The map shades a star within its letter's band with them.
TEMP_RANGES = {
    'O': (32350, 62500), 'B': (10200, 32350), 'A': (7310, 10200), 'F': (5990, 7310),
    'G': (5325, 5990), 'K': (3890, 5325), 'M': (2285, 3890)
}

# Typical main-sequence mass ranges by spectral class (letter), in solar
# masses, per standard stellar-classification references: dwarfs broken
# out by spectral letter -- used to scale a star's galactic Hill
# radius (see `spaceSector.hill_radius_ly`) across spectral types.
SPECTRAL_MASS_RANGES = {
    'O': (15, 90), 'B': (2, 16), 'A': (1.4, 2.1), 'F': (1.0, 1.4),
    'G': (0.89, 1.07), 'K': (0.58, 0.9), 'M': (0.08, 0.58)
}

# --- Orbital motion (physics.planets.generate_orbital_motion_properties) ---
# Real solar-system planets sit close to the ecliptic (Mercury, the most
# tilted, is ~7 degrees) -- 10 degrees gives a little headroom without
# implying a genuinely different orbital-plane population. Moons range much
# wider in reality (close-in regular moons are typically near-equatorial,
# but captured/irregular moons can be tilted or even retrograde) -- 25
# degrees is a simple, larger-but-still-modest scatter that doesn't require
# separately modeling regular vs. irregular moon populations.
PLANET_ORBITAL_INCLINATION_MAX_DEG = 10.0
MOON_ORBITAL_INCLINATION_MAX_DEG = 25.0

# Axial rotation ("day length") ranges, by body_type, in hours, for a body
# no tide has slowed (`spin.planet_day_hours`). Terrestrial: giant impacts
# leave a newly formed rocky planet spinning in a few hours to about two
# days ([R2431]:L21), drawn log-uniform; the slow days
# of Mercury and Venus come from tides, which the tidal-locking rule
# covers (Boss, 2026-10-10: realistic days). Gas giant: real gas giants
# all spin fast (Jupiter ~9.9h, Saturn ~10.7h, Uranus ~17.2h, Neptune
# ~16.1h), drawn uniform.
ROTATION_PERIOD_RANGE_HOURS = {
    "t": (8.0, 48.0),
    "g": (8.0, 20.0),
}

# Tidal-locking (despinning) timescale, per the standard simplified
# formula (Murray & Dermott, "Solar System Dynamics"; the uniform-sphere
# moment of inertia I = (2/5) m R^2 folds Q/k2 into a single 2Q/(15*k2)
# factor -- see physics.planets._tidal_locking_timescale_seconds):
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

# --- S-type (wide) binary stability (orbits.holman_wiegert_critical_semimajor_axis,
# orbits.mutual_hill_radius_au, system.StarSystem) ---

# [R1990], AJ 117:621, "Long-Term Stability of Planets in
# Binary Systems" -- their numerically-fit S-type critical semi-major axis
# formula was tested over companion mass fraction mu in [0.1, 0.9] and
# binary eccentricity e in [0.0, 0.7-0.8]; inputs are clamped to this range
# (orbits.holman_wiegert_critical_semimajor_axis) rather than extrapolated.
HOLMAN_WIEGERT_MU_RANGE = (0.1, 0.9)
HOLMAN_WIEGERT_ECCENTRICITY_RANGE = (0.0, 0.8)

# [R1583], "Dynamical stability of the outer solar system and the
# delivery of comets," Icarus 106, 247 -- for two coplanar, circular-orbit
# planets sharing one central mass, a mutual separation of at least
# 2*sqrt(3) (~3.464) mutual Hill radii guarantees Hill stability (no close
# encounters ever occur). Used by
# system.StarSystem._validate_cross_star_clearance as a physically-
# motivated extension of this criterion across two planets that orbit
# *different* stars of a wide binary rather than one shared star -- see
# that method's docstring and orbits.mutual_hill_radius_au for the caveat.
GLADMAN_MUTUAL_HILL_STABILITY_FACTOR = 2 * math.sqrt(3)

ELECTRON_VOLT_J = 1.602176634e-19  # One electron volt in joules (the elementary charge times 1 V)

# Earth's XUV flux today, 4.6 erg/cm^2/s at 1 AU [R3753], W/m^2:
# the unit of a planet's `xuv_flux_earth` (GEN.86).
EARTH_XUV_FLUX_W_M2 = 4.6e-3
