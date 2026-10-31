# planetgen/tuning.py

"""
Program and Procedural-Generation Constants
==============================================

This module holds constants that are specific to this program's design
rather than physical fact: generation-tuning knobs and probabilities,
display/rounding/formatting thresholds, and the large data tables that
drive procedural content — planet class definitions, life chemistry,
evolutionary timelines and narrative text, asteroid belt components, and
flavor text pools. See `physical_constants` for real physical/astronomical
constants. Nothing in this module has side effects; it is imported wherever
these values or tables are needed.

Rules for a new constant (`tests/test_fuzz_utils_and_config.py` checks them, so
a breach shows as a red test): a number is finite; a pair is written
`(low, high)` and `low <= high`; a `_MIN`/`_MAX`, `_LOW`/`_HIGH` or
`_INNER`/`_OUTER` pair of names is ordered the same way; and a scalar at the
top level of the module is positive. A value that can be negative is kept in a
dict or a list, not as a bare module constant; is named for what it is
(`EXPONENT`, `ROUND_`, `_OFFSET`, `_DELTA`, `_SLOPE` and `_INTERCEPT` are the
names the test lets through); or is stored as a positive "fall" (how far
down, not minus how far) and subtracted where it is used.
"""

import math as _math

from planetgen.generation.phenomena import registry as _phenomena_registry
from planetgen.physics import constants as _physical_constants
from planetgen.physics import mass_boundaries as _boundaries
from planetgen.physics import unit_registry as units

_myr_to_gyr = units.converter("Myr", "Gyr")
_mearth_to_kg = units.converter("Mearth", "kg")
_kg_to_mjup = units.converter("kg", "Mjup")

# --- Generation directives (GEN.96, GEN.180) ---
DIRECTIVE_ATTEMPTS = 200
"""int: K, the draws made at most for a sector that must meet a directive's minimums (`--directive`); a habitable-worlds
minimum still missed after K draws takes the forced fill (`generation.directives`, status `met_forced`)."""

PEAK_RELATIVE_DENSITY = 160.0
"""float: A bound above the densest relative stellar density the galaxy model gives (154.2 at the centre, measured with
`density.relative_density`), for judging up front what a `galaxy` or `sectors` run's densest sector could hold."""

DIRECTIVE_IMPOSSIBLE_BELOW = 1e-9
"""float: A directive whose minimum a sector reaches with a probability below this, and that the forced fill cannot reach
either, is refused before any drawing (`generation.directives.impossible`)."""

# --- Planet Generation Parameters ---
# Domingos, Winter & Yokoyama (2006), MNRAS 373:1227, "Stable satellites
# around extrasolar giant planets" -- a prograde moon on a circular orbit
# stays bound out to about 0.4895 of its planet's Hill radius; beyond that
# the star strips it (physics.planets.generate_moons' outer limit).
MOON_PROGRADE_STABLE_HILL_FRACTION = 0.4895


# Giant masses are drawn from dN/dlogM ~ M^-GIANT_MASS_FUNCTION_SLOPE inside
# their class's range: Cumming et al. (2008, PASP 120:531) measured about
# 0.31 for radial-velocity giants, so lighter giants are commoner but a
# third of class J giants are still past a Jupiter mass.
GIANT_MASS_FUNCTION_SLOPE = 0.31

# A giant drawn around a star is at most this fraction of the star's mass
# (but never below its class's floor): a disk holds about 1% of its star's
# mass, so a red dwarf can't build a 10 Jupiter-mass planet (the Sun's
# limit is about 10 Jupiter masses).
GIANT_MAX_STAR_MASS_RATIO = 0.01

# The planet mass law once no class supplies the mass (GEN.650, `physics/planet_mass.py`) is the planet architecture's law for one
# oligarch (GEN.728, GEN.730, `physics/planet_architecture.py`): the isolation mass of the feeding zone, `OLIGARCH_SEPARATION_MEDIAN`
# mutual Hill radii wide, in a disc scattered by `DISC_SOLIDS_SCATTER_DEX`. Each value is a planetGen stand-in (registry section 6.36),
# tuned against the observed occurrence (`physics/data/planet_occurrence.json`).
# A core that reaches this many Earth masses before the gas disc is gone can run away into a gas giant (core accretion:
# Pollack et al. 1996; Rafikov 2006 give 5 to 20 Earth masses, recalled).
GIANT_CORE_CRITICAL_MASS_EARTH = 10.0
# The share of runaway-capable cores that still find gas left in the disc around a star of the disc's metallicity (scaled for [Fe/H] by
# `metallicity.giant_occurrence_factor` and `GIANT_GAS_SURVIVAL_METALLICITY_EXPONENT`); the rest stay ice giants. Tuned (tier T6, tweak L3)
# to the giants at 2 to 8 AU round a Sun-like star (Fulton et al. 2021, arXiv:2105.11584, read: 0.14), their fall with the star's mass
# (Johnson et al. 2010, arXiv:1005.3084, read) and the fall of the hot Jupiters round M dwarfs (Gan et al. 2023, read).
GIANT_GAS_SURVIVAL_SHARE = 0.766
# Beyond this distance the share falls as the square of the ratio: the cores of the outer disc grow slowly, so few reach
# the critical mass before the gas is gone (direct imaging finds giants beyond 10 AU around about 1 percent of stars).
GIANT_GAS_SURVIVAL_CUT_AU = 9.022
# A core that finds no gas stops growing near an ice giant's mass (Uranus and Neptune are 14.5 and 17.1 Earth masses): the
# late solids beyond the snow line are not all added to it. Earth masses.
GIANT_NO_GAS_CORE_MAX_EARTH = 20.0
# The most a solid core inside the snow line reaches without gas, in Earth masses: heavier rock-and-iron cores bind a hydrogen
# envelope and either run away or are rare (the hot-Neptune desert; super-Earths top out at 10 to 20 Earth masses: Fulton et al.
# 2017, Owen and Wu 2013, recalled). A planetGen stand-in (tier L3), tuned with the pull of the outer-disc taper below.
PLANET_INNER_CORE_MAX_EARTH = 12.0
# How sharply a core is pressed under its cap (`planet_mass.saturate_kg`): the larger, the closer to its own mass until the cap.
PLANET_CORE_SATURATION_SHARPNESS = 4.0
# A runaway giant is at least this heavy: below it the envelope has not yet outweighed the core (Earth masses).
GIANT_RUNAWAY_MIN_MASS_EARTH = 30.0
# How many times the minimum-mass nebula's solids the mass law assumes: below one, because the oligarchs are spaced wide and the
# isolation mass goes as the spacing to the 3/2 (`OLIGARCH_SEPARATION_MEDIAN`); the Kepler systems need several times the nebula
# in all (Chiang and Laughlin 2013, recalled). A planetGen stand-in (tier T6, tweak L3), tuned to the observed occurrence.
PLANET_DISC_SOLIDS_ENHANCEMENT = 0.4647
# How the solids the planets are built from scale with the star's mass, as a power of the mass in suns. Shallower than the
# dust mass of the disc (`DISK_MASS_STELLAR_MASS_EXPONENT`, which sets how many objects a disc can hold): the planets of M dwarfs
# are about as many as those of Sun-like stars (Dressing and Charbonneau 2015; Hsu et al. 2020, recalled), as pebble drift
# concentrates solids close to a cool star. A planetGen stand-in (tier L3), tuned to that occurrence.
PLANET_SOLIDS_STELLAR_MASS_EXPONENT = 1.2
# Inside 1 AU the solids rise faster than the nebula's -1.5 power: the Kepler systems' minimum-mass extrasolar nebula falls as
# about the -1.6 to -2 power (Chiang & Laughlin 2013, recalled), and solids drifting in feed the inner disc. The surface density
# is multiplied by (distance / 1 AU) to the minus this power inside 1 AU. Tuned (tier L3) with `PLANET_DISC_SOLIDS_ENHANCEMENT`
# so that a Sun-like star has about 0.8 planets of 1 to 20 Earth masses inside 0.45 AU and an M dwarf about 1 (Petigura et al.
# 2018; Dressing & Charbonneau 2015: about 2.5 planets of 0.5 to 4 Earth radii within 50 days; recalled).
PLANET_INNER_DISC_STEEPENING = 0.75
# Beyond this distance the cores of the outer disc have not had time to assemble before the gas and the planetesimals are gone:
# orbits are slow and the disc is thin (Levison et al. 2010 grow Uranus and Neptune in a few Myr only by pebble accretion,
# recalled), so the assembled solids fall as (this distance / the distance) to `PLANET_OUTER_GROWTH_EXPONENT`. A planetGen
# stand-in (tier T6, tweak L3) tuned with the giants' occurrence (`GIANT_GAS_SURVIVAL_SHARE`) and the count of the Neptune-mass planets beyond the
# snow line (Suzuki et al. 2016: about 0.5 per star, recalled).
PLANET_OUTER_GROWTH_AU = 11.36
PLANET_OUTER_GROWTH_EXPONENT = 2.889
# A planet's solids are never more than this fraction of its star's mass (a disc holds about 1% of its star).
PLANET_MAX_STAR_MASS_RATIO = GIANT_MAX_STAR_MASS_RATIO

# How much more a heavier-metal star's cores find gas, as 10 to this power times [Fe/H], beside `metallicity.giant_occurrence_factor`: with the
# solids the giants' occurrence goes as 10 to 1.2 times [Fe/H] (Johnson et al. 2010, arXiv:1005.3084, read). A stand-in (tier T6, tweak L3).
GIANT_GAS_SURVIVAL_METALLICITY_EXPONENT = 0.6498

# Planet architecture (GEN.728, GEN.729, GEN.730, `physics/planet_architecture.py`, docs/design/composition/planet-architecture.md):
# how a disc builds a chain of oligarchs, how dynamical instability thins it and what is left, from the star and its disc alone.
# The solids of a disc scatter about the median of its star's mass with a root-mean-square 0.7 dex (Andrews et al. 2013, arXiv:1305.5262,
# read: the dispersion of the dust mass about the host-mass relation). A disc light by one dex builds no planet of an Earth radius or more
# inside 1 AU; a heavy one builds a chain of sub-Neptunes and, with the metallicity, the giants. Whether a star has a Kepler-like inner system at all is the sourced draw below.
DISC_SOLIDS_SCATTER_DEX = 0.7
# The spacing at which oligarchs form, in mutual Hill radii: Pu and Wu 2015 (arXiv:1502.05449, read) find the observed spacings of high
# multiples clustered near 12 and propose that systems formed with an effective spacing of 8.1 or more and were thinned by instability. The
# median and the width (a log-normal) are planetGen stand-ins (tier T6, tweak L3), tuned with the sculpting below so that the spacing of the survivors
# peaks near 20 mutual Hill radii and 93 percent of pairs are over 10 apart (Weiss et al. 2018, arXiv:1706.06204, read).
OLIGARCH_SEPARATION_MEDIAN = 15.07
OLIGARCH_SEPARATION_SIGMA_LN = 0.2
OLIGARCH_SEPARATION_MIN = 4.0
# The inner edge of the chain, as an orbital period in days: the occurrence of planets falls inward of a 10-day period at every spectral
# type (Mulders et al. 2015, arXiv:1406.7356, read), the edge a disc truncates at being the co-rotation radius. The first oligarch sits
# at a period drawn log-normally about this median with this width in dex. A planetGen stand-in (tier T6, tweak L3) tuned to the rise of the
# Kepler occurrence with period (Zhu and Dong 2021, arXiv:2103.02127, read).
OLIGARCH_INNER_EDGE_PERIOD_DAYS = 5.0
# The dust sublimation temperature, kelvin: no solid body builds where a blackbody at its distance would be hotter than the silicates
# evaporate, so the chain of a luminous star starts no closer than the distance of this temperature (a = (278.3 K / T)^2 sqrt(L / Lsun) AU).
# About 1500 K (Dullemond and Monnier 2010, arXiv:1002.3056, secondary). A planetGen stand-in (tier T3, tweak L1).
DUST_SUBLIMATION_K = 1500.0
# A belt takes the place of a planet in its slot only where the disc built a body lighter than this, in Earth masses: a belt is what a
# feeding zone leaves when it did not build a planet, so a giant or a body of an Earth mass or more is never swapped for one. A
# planetGen stand-in (tier T5, tweak L0).
BELT_REPLACES_BELOW_EARTH = 1.0
OLIGARCH_INNER_EDGE_SIGMA_DEX = 0.2
# The age to which the chain is thinned by instability, years. Pu and Wu's survival law (arXiv:1502.05449, read) loses 0.7 in spacing
# per decade of time, so a factor of ten in age matters little; the chain is thinned to this age whatever the star's, and the later
# change of orbits is the orbital update (GEN.105). A planetGen stand-in (tier T6, tweak L0).
SCULPT_AGE_YR = 1.0e9
# The inner system (a planet of an Earth radius or more with a period under 400 days). Whether a star has one: the share of stars with
# Kepler-like planets, 30 percent of Sun-like stars (Zhu et al. 2018, arXiv:1802.09526, read), about 75 percent of stars under 5000 K and 35
# percent over 6500 K (Yang et al. 2020, arXiv:2002.02840, read): the masses (in suns) and the shares at them, log-linear in the mass between.
INNER_SYSTEM_FRACTION_MASSES = (0.15, 0.6, 1.0, 1.4)
INNER_SYSTEM_FRACTION_SHARES = (0.85, 0.75, 0.30, 0.30)
# How many planets it starts with before the thinning: 1 + a Poisson draw of this mean. A planetGen stand-in (tier T6, tweak L3), tuned so the
# survivors average 3.0 per Sun-like system and 2.8 per M dwarf (Zhu et al. 2018; Yang et al. 2020).
INNER_INITIAL_COUNT_EXTRA = 5.619
# The spread of that mean between systems (the width in ln of a log-normal factor), so the counts are broader than a Poisson: the survey's
# many single transiting planets need systems that start with few. A planetGen stand-in (tier T6, tweak L3), tuned to the transit multiplicity.
INNER_COUNT_SIGMA_LN = 0.8
# The mass of a planet of the inner system, Earth masses, for a Sun-like star and a median disc: the planets Kepler counts are 1 to 4 Earth
# radii, 1 to 10 Earth masses (Zhu and Dong 2021, arXiv:2103.02127, read). The scale goes as the star's mass to the exponent (a planet-to-star
# mass ratio about constant: Mulders et al. 2015), as the disc's scatter to the exponent (a heavy disc builds heavy planets and giants:
# Zhu and Wu 2018) and as 10 to the exponent times [Fe/H] (the occurrence of sub-Neptunes doubles over 0.8 dex, Petigura et al. 2018,
# arXiv:1712.04042, read). Planetgen stand-ins (tier T6, tweak L3).
INNER_PLANET_MASS_EARTH = 3.631
INNER_PLANET_MASS_STELLAR_EXPONENT = 0.6357
INNER_PLANET_MASS_DISC_EXPONENT = 0.01653
INNER_PLANET_MASS_METALLICITY_EXPONENT = 0.7565
# Neighbours of a system are alike in mass: the spread of the log10 of a planet's mass about its system's scale, and how much heavier each next
# planet is (the outer planet is the larger in 65 percent of pairs: Weiss et al. 2018, arXiv:1706.06204, read). Stand-ins (tier T6, tweak L3).
INNER_NEIGHBOUR_SIGMA_DEX = 0.15
INNER_OUTWARD_MASS_DEX_PER_PLANET = 0.04
# The spacing at which the inner planets start, in mutual Hill radii, before the thinning: Pu and Wu's effective 8.1 or more
# (arXiv:1502.05449, read), and the spread of one pair about its system's. A stand-in (tier T6, tweak L3).
INNER_SEPARATION_MEDIAN = 5.983
INNER_SPACING_SIGMA_LN = 0.25
# The inner system ends at this period (days); the outer chain starts this many times beyond it (as a period) or past the last inner planet (as a distance).
INNER_SYSTEM_PERIOD_LIMIT_DAYS = 400.0
OUTER_CHAIN_START_RATIO = 1.3

# The lightest body that is a planet, Earth masses; lighter oligarchs are the debris of a belt (Mercury is 0.055). A stand-in (tier T6, tweak L0).
PLANET_MIN_MASS_EARTH = 0.1
# Beyond the snow line a body this light is a dwarf planet of the belt beyond, not a planet: the ice makes the planets there heavier (Earth
# masses). A stand-in (tier T6, tweak L3).
PLANET_MIN_MASS_OUTER_EARTH = 0.6016
# The share of the mass of two oligarchs that collide which stays in the planet; the rest is debris (Pu and Wu 2015 expect a destructive
# collision at the relative speeds of the compact systems and cannot predict the outcome). A planetGen stand-in (tier T4, tweak L0).
COLLISION_MASS_RETENTION = 0.9
# The share of giants that migrate in to a hot Jupiter orbit, and how its rate rises with [Fe/H] beyond the giants' (10 to the 3.4
# for hot Jupiters, Petigura et al. 2018, arXiv:1712.04042, read, against 1.2 for giants, Johnson et al. 2010, arXiv:1005.3084, read):
# a giant formed beyond the snow line migrates in with this probability at solar metallicity. Tuned (tier T6, tweak L3) to 0.6 percent of Sun-like stars
# (Zhu and Dong 2021, arXiv:2103.02127, read) and 0.27 percent of early M dwarfs (Gan et al. 2023, arXiv:2210.08313, read).
GIANT_MIGRATION_SHARE = 0.0092
GIANT_MIGRATION_METALLICITY_EXTRA_EXPONENT = 2.2
# Where a migrated giant ends, as an orbital period in days: log-uniform between these (the hot Jupiters pile up near 3 days).
GIANT_MIGRATION_PERIOD_RANGE_DAYS = (1.0, 10.0)
# A giant that does not become a hot Jupiter still moves in from where it formed: the occurrence of giants rises at 1 AU by a factor of 4 and
# peaks at 1 to 3 AU, inside the snow line (Fulton et al. 2021, arXiv:2105.11584, read; Cumming et al. 2008, arXiv:0803.3357, read). It
# moves in by the factor e to minus the absolute value of a normal of this width in ln. A planetGen stand-in (tier T6, tweak L3).
GIANT_INWARD_MIGRATION_SIGMA_LN = 0.7212
# A giant that moved in scatters the planets this close, in mutual Hill radii, and a hot Jupiter every planet inside its orbit up to this
# distance in AU (planets inside a hot orbit are gone). Stand-ins (tier T6, tweak L0).
GIANT_CLEARING_SPACING = 8.0
GIANT_HOT_ORBIT_AU = 0.1

# The kind of body a caller's mass or radius alone names (GEN.681, `physics/planet_kind.kind_of_given`), when the formation law
# was not asked: a radius from this many km up is a giant (the smallest radius a Neptune-like envelope body reaches; the old
# lowest giant class radius). Tier "from needed output", level L1.
GIANT_MIN_RADIUS_KM = 15000.0
# The mean density a given mass and radius may imply for a giant, g/cm3: a hot Saturn-mass giant swells to under 0.1, thirteen Jupiter masses
# at 0.8 of Jupiter's radius is about 40. A terrestrial body's bounds are `bulk_draw.DENSITY_BOUNDS_G_CM3`.
GIANT_DENSITY_BOUNDS_G_CM3 = (0.05, 60.0)
# The mass (kg) and radius (km) a caller may give a planet: the lightest body bulk_draw builds (1e12 kg) to the deuterium-burning limit,
# 13 Jupiter masses (a heavier body is a brown dwarf), and from half a kilometre to the largest inflated giant.
GIVEN_MASS_KG_BOUNDS = (1.0e12, 2.47e28)
GIVEN_RADIUS_KM_BOUNDS = (0.5, 200000.0)

# How many planets a request that names a class (a slot's `planet_class`, `FORCE_CLASS`, a forced habitable world) draws at its distance
# before it keeps the last (GEN.688; `StarSystem.build_planet`).
REQUEST_MAX_PLANET_ATTEMPTS = 40
# A request for a habitable world draws this many planets freely, then draws planets of an Earth-like mass (log-uniform between the two bounds,
# Earth masses) for the rest: a planetGen stand-in (registry row CC10) found by the sweep, where a free draw in the ecosphere of an M5 dwarf
# gave a life class in none of 30 and an Earth-like mass in 20 of 30.
HABITABLE_REQUEST_FREE_ATTEMPTS = 8
HABITABLE_REQUEST_MASS_EARTH = (0.5, 2.0)
CLASS_REQUEST_FREE_ATTEMPTS = 8
CLASS_REQUEST_MASS_EARTH = (0.05, 1000.0)
"""int, tuple: A request for a planet class (a slot's `planet_class`, `FORCE_CLASS`) is first left to chance for this many draws; the draws
after it take a mass log-uniform over this range in Earth masses, so a class the disc rarely makes at that slot (an ice giant at a cold
one) is still reached inside `REQUEST_MAX_PLANET_ATTEMPTS`. A planetGen stand-in (tier T9, tweak L3; registry row PA14)."""
# What a request for a habitable world asks of a finished planet or moon (GEN.916, `generation/habitable_request.py`): conditions, not a
# class letter. The surface temperature range is the liquid-water range, the pressure floor is the air a world must hold to count as having
# any (10 kPa), and the water a surface ocean; the target is the middle of the range, where a draw that missed by its temperature aims the
# next (GEN.911).
REQUEST_TEMPERATURE_K = (273.15, 323.15)
REQUEST_MIN_PRESSURE_PA = 10000.0
REQUEST_HYDROSPHERES = ("surface ocean",)
REQUEST_TARGET_K = 298.15
# A habitable request leaves its first draws at the distance it was given, then aims each next draw by the temperature the last one solved.
HABITABLE_REQUEST_AIM_AFTER_ATTEMPTS = 2

# How many times a planet or moon is drawn again when its finished class is one that can carry life and the system bars habitable worlds
# from its ecosphere (`HABITABLE_WORLD` False); the last draw is kept.
BARRED_CLASS_MAX_ATTEMPTS = 30

# The lightest planet the formation law may give (GEN.646): the heaviest a heap of ice and rock can be before its own gravity rounds it
# (`bulk_structure.SMALL_BODY_MASS_KG`, about 540 km radius); the old lowest class radius was 500 km. Lighter bodies are the small bodies
# (belts, comets), not planets. Tier "from needed output", level L1.
PLANET_MIN_MASS_KG = 1.0e21

# The moon size law once no class supplies the radius (GEN.646, `physics/planet_kind.draw_moon_mass_kg`). All four values are
# planetGen stand-ins (tier "educated guess from physics", level L3) until a moon-formation law is wired:
# moon masses follow dN/dlogM ~ M^-MOON_MASS_FUNCTION_SLOPE between the lightest moon the generator makes and a share of the
# parent's mass. Regular satellites of the giants are at most a few 1e-4 of the parent (Titan 2.4e-4 of Saturn; Ganymede 7.8e-5
# of Jupiter), so the giants' cap is 1e-3; a rocky parent can hold a moon of a tenth of its mass (the Moon 0.012, Charon 0.12 of
# Pluto), which is also the largest `physics/planets.moon_size_limits` allows.
MOON_MASS_FUNCTION_SLOPE = 0.5
# The lightest moon generated (about an Enceladus, 250 km across): the generator leaves the swarm of smaller captured rocks to the belts and
# rings. A planetGen stand-in (tier L3): the old class draw built moons of 500 km and more, 1.1 per planet.
MOON_MIN_MASS_KG = 1.0e20
MOON_GIANT_PARENT_MAX_MASS_RATIO = 1.0e-3
MOON_TERRESTRIAL_PARENT_MAX_MASS_RATIO = 0.1

# The water a rocky terrestrial body is born with (GEN.849; `physics/delivered_water`): its class (dry, low, high, rich) comes from the Bern
# population synthesis law the atmosphere recipes already read (`data/ngpps_formation_water.json`), by its mass and the snow lines it formed
# at; inside a class the share of its whole mass is log-uniform. The class edges are the model's. The dry class reaches down to
# `DELIVERED_WATER_DRY_MIN` (Venus and Mercury hold far less than Earth's 2.3e-4 in the oceans), the rich class up to
# `DELIVERED_WATER_RICH_MAX`, the largest water share a rocky body keeps and still has the rock to be a planet (Thiabaud et al. 2014
# find at most about 0.5 for the bodies the model formed). Both ends are stand-ins (registry rows DW2 and DW3, tier "educated guess
# from physics", level L1). A body that formed beyond the snow line of a class that is not rocky is ice and rock: its share is
# drawn from `ICY_WATER_MASS_FRACTION_FLOOR` up to the water share of the solids that the disk temperature at its formation distance condenses (`condensate_shares`, Lodders 2003 Table 11: 0.54 to 0.38, GEN.950; Europa 8 percent water, Ganymede about 45, Marboeuf et al. 2014: 5 to 50 percent).
DELIVERED_WATER_DRY_MIN = 1.0e-6
DELIVERED_WATER_RICH_MAX = 0.5
ICY_WATER_MASS_FRACTION_FLOOR = 0.05

# A terrestrial body's Bond albedo, drawn before its air is solved so it can only read what the body is made of and the starlight
# it gets (GEN.684; `physics/climate_draws.albedo_range`): dry rock keeps the generator's own range (`climate_draws.DEFAULT_ALBEDO_RANGE`); a body with surface water
# (a water mass fraction of at least `ALBEDO_WET_MIN_WATER_FRACTION`) and an equilibrium temperature under
# `ALBEDO_FROZEN_BELOW_K` (at the middle albedo) is ice-covered and bright; one above it is ocean and cloud. A Venus-like
# cloud deck (albedo 0.75 to 0.9) is not reached until the clouds are drawn before the temperature (GEN.414 is drawn after).
# Stand-ins, tier "educated guess from physics", level L3 (hand-tuned against the old class albedo ranges; the Earth, at 255 K with
# no air, and the outer edge of the habitable zone, at 223 K, must stay on the ocean side of 230 K or an Earth twin freezes).
ALBEDO_WET_RANGE = (0.25, 0.40)
ALBEDO_FROZEN_RANGE = (0.45, 0.70)
ALBEDO_WET_MIN_WATER_FRACTION = 1.0e-4
ALBEDO_FROZEN_BELOW_K = 230.0
ALBEDO_REGIME_REFERENCE = 0.30

# --- Star System Generation Parameters ---
ASTEROID_BELT_PROBABILITY = 0.1
ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN = 1.1
ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX = 2
ABSOLUTE_MAX_SYSTEM_OBJECTS = 500
MIN_ASTEROID_BELT_SEPARATION = 0.05

# Minimum stable orbital separation between two adjacent planets, in units
# of their *mutual* Hill radius (orbits.mutual_hill_radius_m -- the pair's
# combined mass and average distance, not either body's own individual
# Hill radius alone). The analytically rigorous minimum for guaranteed
# two-planet Hill stability in the circular, coplanar case is 2*sqrt(3)
# ~= 3.46 mutual Hill radii (Gladman 1993); long-term (10^8-10^9 orbit)
# N-body integrations of systems with more than two planets recommend a
# larger safety margin, commonly cited around 8-10 mutual Hill radii
# (Chambers, Wetherill & Boslough 1996; Smith & Lissauer 1999/2009). This
# generator uses the upper/safer end of that range -- see
# `StarSystem.validate_system`'s planet-planet spacing check, which is
# the only place this applies (a moon's own orbital limit around its
# parent is a single-body Hill-sphere question, not a mutual one -- see
# `Planet.min_orbit_distance`/`physics.planets.generate_moons`).
MUTUAL_HILL_RADII_SEPARATION = 10

# Real protoplanetary disks are observed (sub-mm/ALMA continuum surveys)
# to carry more solid mass around more massive stars, and steeply so --
# not a flat "1 + log10(solar_masses)" scaling. Disk dust-mass-vs-
# stellar-mass surveys (Andrews et al. 2013, Taurus; Pascucci et al. 2016,
# multi-region) find M_dust ~ M_star^1.8 in young (~1-3 Myr) star-forming
# regions, steepening further (~M_star^2.7) in older ones. This scales
# `physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2` for any given
# star relative to the Sun -- see `formation.disk_surface_density_scale`,
# used by `planet_mass.solids_surface_density_gcm2`. The
# younger-region exponent is used since this generator has no notion of a
# system's disk-formation age (only its current, post-formation age).
DISK_MASS_STELLAR_MASS_EXPONENT = 1.8

# Real disks are truncated far short of a star's own galactic-tidal Hill
# sphere (`Star.system_perimeter`, tens to hundreds of thousands of AU) --
# viscous spreading and photoevaporation cut them off at tens to a few
# hundred AU (ALMA disk-size surveys). Rather than a flat AU figure, this
# scales with the same star-dependent quantity that sets where solids can
# even condense in the first place: the snow line
# (`formation.snow_line_au`, drawn per star and stored: `Star.snow_line_au`) --
# our own Solar System's own giant-planet/Kuiper-belt region extends to
# roughly this same multiple (~18x) of its own ~2.7 AU snow line. See
# `StarSystem._plan_architecture`.
DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER = 18

# The scatter of a star's own snow line around the blackbody position for its
# luminosity (`formation.draw_snow_line_au`, GEN.458): 0.2 dex in log10, the
# 0.4 dex scatter of accretion rates carried through the steady-disc snow-line
# law (Mulders et al. 2015, arXiv:1505.03516). The draw is cut at
# `SNOW_LINE_SCATTER_MAX_SIGMA` sigma, a factor of 4 either way, inside the
# factor 2 to 3 the observed lines move.
SNOW_LINE_SCATTER_DEX = 0.2
SNOW_LINE_SCATTER_MAX_SIGMA = 3.0

# The farthest a planet zone in a hierarchy of stars reaches (GEN.314,
# `physics.planet_zones`): a protoplanetary disk around a star or a pair is
# not larger than this, whatever the stability limit allows.
PLANET_ZONE_DISK_CAP_AU = 100.0

# How many times StarSystem.__init__ retries its whole placement loop (fresh
# object count, positions, and validate_system pass, same star) before
# giving up when HABITABLE_WORLD=True/ASTEROID_BELT=True still isn't
# satisfied afterward. Needed because validate_system's own orbital-overlap
# correction can, in rare cases, push a deliberately-placed guaranteed body
# (e.g. the forced Class M world) into a zone its class no longer supports --
# physics.planets.reconcile_zone correctly redraws it away
# rather than reporting a physically inconsistent class, which can undo the
# very guarantee that placement was satisfying. Retrying with a fresh
# placement is simpler and more robust than trying to reshuffle every
# neighbor's spacing to protect one body's zone in place.
MAX_SYSTEM_GENERATION_ATTEMPTS = 8

# How many trees of stars a system draws (GEN.631) while looking for one that is what its directives and its place need: the
# kind of pair `WIDE_BINARY` or a `wide_binary` prevalence asks for, a pair at all under `BINARY_SYSTEM`, a close pair whose
# circumbinary limit leaves room for the habitable world `HABITABLE_WORLD` forces. After this many a forced pair keeps the
# last tree and any other system is its star alone. A forced close pair is the rare draw (about 1 in 14 pairs has an axis
# under 0.25 AU), so the count is large; the first draw is kept whenever nothing is asked.
COMPANION_DRAW_ATTEMPTS = 200

# How many whole systems (star and all) `planetgen system` builds while
# looking for one that meets every forced option (`StarSystem.
# unmet_requirements`) before it gives up with an error and saves nothing.
# Some stars rarely allow a forced body: around an O5V star a habitable
# world fails about two systems in five.
SINGLE_SYSTEM_GENERATION_ATTEMPTS = 5

# --- Binary System Generation Parameters ---

# Holman & Wiegert (1999), AJ 117:621 -- their P-type (circumbinary) fit,
# for a planet orbiting both stars of a close pair
# (orbits.holman_wiegert_circumbinary_a_crit_au), was tested over mu (the
# lighter star's mass fraction) in [0.1, 0.5] and e in [0.0, 0.7]; inputs
# are clamped to this range rather than extrapolated.
HOLMAN_WIEGERT_P_TYPE_MU_RANGE = (0.1, 0.5)
HOLMAN_WIEGERT_P_TYPE_ECCENTRICITY_RANGE = (0.0, 0.7)

# The companions of a star by orbital period, for a primary of each mass (Moe & Di Stefano 2017, ApJS 230:15, Table 13,
# arXiv 1606.05347, "after correcting for observational selection effects"): the companions per decade of log10 of the
# period in days (mass ratios above 0.1) at log P 1, 3, 5 and 7 (each across a decade centred there, a = 0.1 to 10^3
# AU for a pair of two solar masses), for the five mass bins the table has: solar-type (0.8 to 1.2 Msun), A and
# late B (2 to 5), mid-B (5 to 9), early-B (9 to 16) and O (above 16). `generation.companion_classes` draws a
# companion's class and axis from the straight lines between the four points (GEN.448). The bins' centres are the
# geometric means of their edges, and 25 Msun for the open-ended O bin (a design default: the IMF-weighted middle
# of 16 to 150 Msun). The 1 sigma errors (0.009 to 0.11) are reference values for TEST.196, not inputs.
COMPANION_LOG_PERIODS_DAYS = (1.0, 3.0, 5.0, 7.0)
COMPANIONS_PER_DECADE_BY_PRIMARY_MASS = (
    (1.0, (0.027, 0.057, 0.095, 0.075)),
    (3.16, (0.07, 0.12, 0.13, 0.09)),
    (6.71, (0.14, 0.22, 0.20, 0.11)),
    (12.0, (0.19, 0.26, 0.23, 0.13)),
    (25.0, (0.29, 0.32, 0.30, 0.18)),
)

# S-type binary separation is drawn from the companions-per-decade lines above (GEN.448), from this lower bound,
# the top of the close/P-type BinaryStarProxy separation range (0.05-0.25 AU), so the two configurations'
# separations never overlap.
# The upper bound is not a constant: a pair survives only inside the
# widest axis its place allows (GEN.103, `generation.place_terms.
# wide_axis_limit_au`: passing stars and the galactic tide, from the place's
# density, velocity dispersion and the pair's age), and a place that allows
# none under this lower bound keeps no wide pair.
WIDE_BINARY_SEPARATION_MIN_AU = 0.25
WIDE_BINARY_SEPARATION_CEILING_AU = 1.0e6
"""float: Only where a place's density is zero is the survival limit infinite; no pair is a pair past a few parsecs
(the Jacobi radius at the halo's density is under 10 pc), and this keeps the log-uniform draw finite."""

HOSTING_SMALL_BODY_MAX_RADIUS_KM = 500.0
"""float: A body under this radius is a small body to the hosting gate (`planets.hosting_population`, GEN.761); the old smallest-class
(class D) ceiling. A planet starts at `PLANET_MIN_MASS_KG`, so few bodies are under it."""

HOSTING_TRANS_NEPTUNIAN_MIN_SNOW_LINES = 10.0
"""float: A small body in the cold zone farther than this many snow lines from its star is a trans-neptunian body to the
hosting gate (`planets.hosting_population`, GEN.761); nearer, a belt body. STAND-IN [P]: the Solar System's main belt lies at
about 1 snow line (2.2 to 3.3 AU against 2.7) and the Kuiper belt starts at Neptune's 30 AU, about 11; bracket 5 to 15."""

HOSTING_COLD_MAX_INCLINATION_DEG = 5.0
"""float: A trans-neptunian body under this inclination is `trans-neptunian-cold`, the others `trans-neptunian-hot` (the cut of
Stephens and Noll 2006 behind the gate's binary fractions, 22 and 5.5 percent; GEN.761)."""

# Cap on the sampled binary orbital eccentricity for an S-type pair (see
# binary.WideBinaryPair.generate's thermal-distribution sampling).
# Wide binaries never tidally circularize the way the close/P-type pair
# does, and real wide pairs are broadly consistent with a "thermal"
# eccentricity distribution (f(e) = 2e) -- but this is capped short of the
# full [0, 1) thermal range at physical_constants.HOLMAN_WIEGERT_ECCENTRICITY_RANGE's
# own upper bound, since the Holman & Wiegert (1999) stability fit this
# eccentricity feeds isn't validated past there.
WIDE_BINARY_ECCENTRICITY_MAX = 0.8

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
This is a vocabulary of names the materials table resolves (a belt's or field's composition is drawn from asteroid taxa, GEN.442 and GEN.477).
"""

# --- Space Sector Generation Parameters ---

# The standard sector edge length, in whole parsecs: the galaxy grid's
# ring width, layer height and (near enough) slot arc, and the cube edge of
# a standalone sector. A design choice, not a physical constant. 4 pc
# (~13.05 ly, ~2,220 ly^3 of volume, ~6.3 systems at local density -- see
# `physical_constants.LOCAL_STELLAR_DENSITY_LY3`) puts the default Milky
# Way shape's one-star-per-sector edge at ~50,000 ly, the real disk's
# radius; see `docs/design/galaxy/galaxy-coordinate-system.md`, "Sector size".
DEFAULT_SECTOR_EDGE_PC = 4

# `DEFAULT_SECTOR_EDGE_PC` in light-years (~13.046), for the code that
# sizes sectors in light-years.
DEFAULT_SECTOR_EDGE_LY = DEFAULT_SECTOR_EDGE_PC * _physical_constants.AU_PER_PARSEC / _physical_constants.LY_TO_AU

# Maximum offset (on each axis) from dead-center when placing a sector's
# "home" system -- see `spaceSector.SpaceSector.add_home_system`.
DEFAULT_HOME_JITTER_LY = 1.0

# Retry cap for rejection-sampling a new system's random position against
# every existing system's minimum-separation requirement before giving up.
SECTOR_MAX_PLACEMENT_ATTEMPTS = 100

# --- Space Sector Growth Parameters (Poisson-disk placement) ---

# Candidate attempts tried around an active point (Bridson's Fast Poisson
# Disk Sampling active-list scheme -- "Fast Poisson Disk Sampling in
# Arbitrary Dimensions", SIGGRAPH 2007 sketch) before giving up on it and
# removing it from the active list -- see
# `spaceSector.SpaceSector.grow_from_seed`. This program's own tuning
# choice, not a physical constant; Bridson's own sketch uses k=30 for its
# 2D/3D examples, kept as-is here since nothing about this sector's scale
# (a handful to a few dozen systems) argues for a different value.
SECTOR_GROWTH_POISSON_DISK_K = 30

# Maximum number of "push away from whichever star turned out to be
# nearest" nudges `SpaceSector.grow_from_seed`'s fine-tuning step will try
# for a single candidate position before giving up on it (moving away from
# one violating neighbor can bring the candidate closer to a different one,
# so this can take more than one nudge to converge). This program's own
# tuning choice; the sector's realistic scale (a handful of systems) means
# this should converge in only one or two nudges the vast majority of the
# time, so this cap is generous headroom rather than an expected typical
# count.
SECTOR_GROWTH_FINE_TUNE_MAX_ITERATIONS = 20

# Floating-point slack, in light-years, below which a fine-tuning deficit
# (required separation minus actual distance) is treated as "resolved"
# rather than a real violation -- see
# `spaceSector.SpaceSector._fine_tune_position`. Without this, nudging a
# position to exactly its required distance can leave a residual
# floating-point error (sqrt/division noise) of a few parts in 1e-15, which
# would otherwise register as a still-positive deficit forever and nudge
# the same point back to the same place every remaining iteration instead
# of ever converging. A numerical-precision safety margin, not a physical
# or generation-design choice -- 1e-9 ly is still astronomically negligible
# next to any Hill-sphere-scale distance this module works with.
SECTOR_GROWTH_FLOATING_POINT_TOLERANCE_LY = 1e-9

# Roman-numeral Octant (see `spaceSector.classify_octant`) labels for the
# 8 sign-combinations of an (x, y, z) position relative to a sector's
# center -- displayed as "Octant" (not "Quadrant") specifically to stay
# distinct from planetgen/web/maps/galaxymap.py's own, unrelated galaxy-scale
# Quadrant concept (4 azimuthal regions spanning many sectors); see
# spaceSector.py's module docstring, "Named locations (octants)". NOT a
# universal mathematical standard -- unlike the 2D I-IV quadrant
# convention, there is no single authoritative numbering for 3D octants
# (Wikipedia's "Octant (solid geometry)" article recommends explicit
# sign-tuple notation instead, precisely because no such standard exists).
# This table adopts a
# commonly *taught* (not ISO-standardized) extension of the 2D pattern:
# octants with z >= 0 are numbered I-IV in the same counterclockwise sign
# pattern as the 2D quadrants, then z < 0 continues the same x/y pattern as
# V-VIII. Each tuple is (roman_numeral, x_is_non_negative,
# y_is_non_negative, z_is_non_negative); a coordinate of exactly 0.0 is
# treated as non-negative (an arbitrary but consistent tie-break).
SECTOR_OCTANT_LABELS = [
    ("I",    True,  True,  True),
    ("II",   False, True,  True),
    ("III",  False, False, True),
    ("IV",   True,  False, True),
    ("V",    True,  True,  False),
    ("VI",   False, True,  False),
    ("VII",  False, False, False),
    ("VIII", True,  False, False),
]

# Decimal places used when formatting a named location's magnitudes -- see
# `spaceSector.format_named_location`.
SECTOR_LOCATION_DECIMAL_PLACES = 2

# --- Display / Formatting Parameters ---
SCIENTIFIC_NOTATION_DECIMAL_PLACES = 2
ROUND_TEMPERATURE_NEAREST_HUNDRED = -2

# --- Star Generation Parameters ---
# The current age of the universe (~13.8 Gy, standard cosmological
# consensus) -- an absolute ceiling on any star's age: no star is older, however
# long its own track could keep it shining (an M dwarf's runs to trillions of years).
UNIVERSE_AGE_GY = 13.8
UNREACHABLE_PLANET_AGE_MIN_LIFESPAN_RATIO = 0.85
OLD_STAR_AGE_LIFESPAN_RATIO = 1/3
YOUNG_STAR_AGE_LIFESPAN_RATIO = 1/3
MAX_PLANET_AGE_ADJUSTMENT_FACTOR = 0.95
WHITE_DWARF_AGE_ADDITION_GY = 5
# The heavy-element mass fraction the lifetime's [Fe/H] is measured from:
# the proto-solar Z of the composition scale, Asplund 2021 (stellar-
# composition.md D1). A solar mass lasts 10.2 Gy at it.
MS_LIFETIME_REFERENCE_Z = 0.0154
MS_LIFETIME_Z_RANGE = (1.0e-4, 0.03)
"""tuple: The metallicities (mass fraction of metals, Z) Hurley, Pols and Tout 2000 fitted their formulae over (the
paper's abstract: 0.0001 to 0.03); `ms_lifetime.zeta` holds a star's to it (GEN.651)."""
# The spectral types hotter than O3, which the Pecaut and Mamajek sequence does not hold (GEN.943): the
# effective temperatures of O2, O1 and O0 (stand-ins, an extrapolation of the sequence's step of about 2,000 K a
# type to 5,000 K, the hot end of the observed O2 stars near 50,000 K, Walborn et al. 2002; Martins et al. 2005).
HOT_O_SUBCLASS_TEMPERATURES_K = {0: 60000.0, 1: 55000.0, 2: 50000.0}
# A main-sequence star this far below the luminosity of a solar-metallicity dwarf of its temperature (log10 of the
# luminosity) is a subdwarf, Yerkes VI (GEN.943). Subdwarfs lie one to two magnitudes below the dwarf sequence in
# visual magnitude (Gizis 1997; Yerkes/MK class sd); in bolometric terms the model's metal-poor stars are 0.4 dex dimmer
# (F and G), 0.6 (K) and 1.0 (M) at [Fe/H] -1.5, while a solar-metallicity M dwarf of the model sits 0.17 dex under the
# sequence (95 percent within 0.23), so 0.3 labels the first and not the second.
SUBDWARF_LUMINOSITY_DEFICIT_DEX = 0.3
# A forced star type (GEN.941, generation/forced_star.py) is a star the model
# makes. A subdwarf (Yerkes VI) is a main-sequence star under-luminous for its
# temperature (SUBDWARF_LUMINOSITY_DEFICIT_DEX), which a metal-poor star of the
# halo is: Gizis 1997 define the sd classes as metal-poor dwarfs, [Fe/H] below
# about -1; the range runs to the halo's lower edge. For the hot letters it is
# also a hot horizontal-branch star (sdB and sdO stars are helium-burning cores,
# Heber 2016) from SUBDWARF_HOT_MIN_TEMPERATURE_K, the extreme horizontal
# branch's floor (HORIZONTAL_BRANCH_COLOUR_TEMPERATURE_RANGE_K).
SUBDWARF_FE_H_RANGE = (-2.0, -1.0)
SUBDWARF_HOT_MIN_TEMPERATURE_K = 20000.0
# The scan the forced-type table is built from: masses log-spaced across the
# mass function, and per mass the ages of the main sequence, of the stretch
# after it and the white dwarf cooling (log-spaced), each count a resolution of
# the scan, not a physical value. A cell is jittered by up to half its size
# before the star is evolved, tried FORCED_STAR_JITTER_TRIES times before the
# cell's own star is taken. A type the model never makes takes the nearest
# model star of its Yerkes class by temperature, among those within
# FORCED_STAR_NEAREST_SLACK (in ln T) of the nearest.
FORCED_STAR_MASS_STEPS = 140
FORCED_STAR_AGE_STEPS = (60, 300, 80)
FORCED_STAR_JITTER_TRIES = 40
FORCED_STAR_NEAREST_SLACK = 0.03

# Real stellar-multiplicity surveys consistently find binary/multiple
# companionship rising with primary mass, not a single flat rate: low-mass
# M dwarfs are the *least* likely to have a companion, while the most
# massive O stars are very nearly certain to. Keyed by spectral letter
# (`Star.type[0]`) the same way SPECTRAL_PROBABILITIES_NORMAL above is,
# and used the same way when `SystemConfig.BINARY_SYSTEM` is left at its
# default `None` -- see `StarSystem._companion`. Anchor
# points, in solar-mass order:
#   M: 0.26  -- Duchene & Kraus 2013 (Annual Review of Astronomy and
#               Astrophysics 51:269-310), low-mass-star multiplicity
#               fraction 26 +/- 3%.
#   G: 0.44  -- Raghavan et al. 2010 (ApJS 190:1), solar-type (F6-K3)
#               multiplicity fraction 46% (Duchene & Kraus's own review
#               cites this population as 44 +/- 2%; G anchors the middle
#               of that same F/G/K grouping here).
#   F, K: 0.47/0.40 -- interpolated either side of the G anchor along the
#               same "solar-type" grouping (F/G/K aren't broken out
#               separately in the literature above), consistent with the
#               broader monotonic-with-mass trend every other anchor here
#               shows.
#   A, B: 0.55/0.65 -- Duchene & Kraus 2013 describe intermediate-mass
#               (A/B) multiplicity as ">=50%"; split across the two
#               letters along the same mass trend, both comfortably at or
#               above that floor. Moe & Di Stefano 2017 Table 13 puts the
#               A / late-B bin (2 to 5 solar masses) at 0.59 +/- 0.11, so the
#               letter's 0.65 stands for late B (9 to 16 and 5 to 9 solar
#               masses are the two rows below).
#   B-mid, B-early: 0.76/0.84 -- Table 13's mid-B (5 to 9 solar masses) and
#               early-B (9 to 16 solar masses) single-star fractions, 0.24
#               and 0.16 (`MULTIPLICITY_BY_PRIMARY_MASS_BIN`), taken from
#               the primary's birth mass (GEN.449).
#   O: 0.90  -- Moe & Di Stefano 2017 (ApJS 230:15) revise O-star
#               multiplicity up to 94 +/- 14%; kept just under certainty
#               given that real uncertainty rather than treating O stars
#               as *always* binary.
MULTIPLICITY_BY_PRIMARY_MASS_BIN = {
    "solar": {"masses": (0.8, 1.2), "multiplicity": 0.50, "close": 0.15, "single": 0.60, "binary": 0.30, "multiple": 0.10},
    "A / late B": {"masses": (2.0, 5.0), "multiplicity": 0.84, "close": 0.37, "single": 0.41, "binary": 0.37, "multiple": 0.22},
    "mid B": {"masses": (5.0, 9.0), "multiplicity": 1.3, "close": 0.63, "single": 0.24, "binary": 0.36, "multiple": 0.40},
    "early B": {"masses": (9.0, 16.0), "multiplicity": 1.6, "close": 0.8, "single": 0.16, "binary": 0.32, "multiple": 0.52},
    "O": {"masses": (16.0, 40.0), "multiplicity": 2.1, "close": 1.0, "single": 0.06, "binary": 0.21, "multiple": 0.73},
}
"""dict: Moe and Di Stefano 2017 (ApJS 230, 15), Table 13, multiplicity statistics after correction for selection effects,
by primary mass bin (solar masses). `multiplicity`: companions per primary with q above 0.1; `close`: companions
per primary inside log P (days) 3.7; `single`, `binary`, `multiple`: the fractions of primaries with no companion, one,
and two or more (triples and quadruples). The chance of a companion is 1 - `single`. The rows of
`generation.phenomena.hosts.STAR_COUNT_PROBABILITIES` for mid and early B stars are rederived from these (GEN.449)."""

# --- Star population model (stellarEvolution.py) ---
# A random star is drawn as physics, not from a letter table: a mass from
# the initial mass function, an age from the star-formation history, and
# the star's present state (main sequence, subgiant, giant, supergiant or
# white dwarf) from how far that age is into its own lifetime. Its letter,
# subclass and Yerkes class follow from the resulting temperature and
# luminosity. SPECTRAL_PROBABILITIES_NORMAL above is kept as the census the
# model is tested against, not as a draw.

# The single-star initial mass function (GEN.927): Chabrier (2003), PASP 115:763, Table 2, the Galactic disk, individual stars:
# a log-normal in log mass below 1 Msun and a power law above it, joined continuously. All four are recalled, not read in the
# paper (tier Confident Unverified, level L0; docs/design/galaxy/stellar-mass-function.md).
IMF_BELL_CENTRE_SOL = 0.079
"""float: The solar mass at which the log-normal part of the single-star mass function peaks per logarithm of mass (Chabrier
2003, Table 2, disk, individual stars: `m_c` = 0.079; the systems' own function, 0.22, is derived here, not set)."""
IMF_BELL_WIDTH_DEX = 0.69
"""float: The width of that log-normal in dex of mass (Chabrier 2003: `sigma` = 0.69)."""
IMF_BELL_TOP_SOL = 1.0
"""float: The mass at which the log-normal gives way to the power law above it (Chabrier 2003: 1 solar mass)."""
IMF_HIGH_SLOPE = 2.3
"""float: The power law above `IMF_BELL_TOP_SOL`, `dN/dM ~ M^-2.3` (Chabrier 2003: `x` = 1.3 per logarithm of mass; Kroupa 2001:
2.3 above 0.5). A birth function: the steeper 2.7 of the field (Kroupa, Tout and Gilmore 1993) counts stars still alive and
is what the disk history reproduces from this and the ages (docs/design/galaxy/stellar-mass-function.md)."""
IMF_NODE_SPACING_DEX = 0.05
"""float: The spacing of the nodes the log-normal is cut into, each pair joined by a power law through the function at both
(the function is kept as a continuous piecewise power law, `physics/imf.py`); at 0.05 dex the chord is within 0.07% of
the log-normal."""


def _bell_rows():
    """`(breaks, slopes)` of the log-normal below `IMF_BELL_TOP_SOL` cut into power-law segments at `IMF_NODE_SPACING_DEX`
    (a segment's slope is the chord of `ln dN/dM` against `ln M` between its nodes), then the high power law up to 150 Msun."""
    low = _boundaries.sol("hydrogen_burning_limit")
    cells = int(_math.ceil(_math.log10(IMF_BELL_TOP_SOL / low) / IMF_NODE_SPACING_DEX))
    breaks = [low * (IMF_BELL_TOP_SOL / low) ** (i / cells) for i in range(cells + 1)]

    def log_density(mass):
        return -(_math.log10(mass) - _math.log10(IMF_BELL_CENTRE_SOL)) ** 2 / (2.0 * IMF_BELL_WIDTH_DEX ** 2) - _math.log(mass)

    slopes = [-(log_density(b) - log_density(a)) / _math.log(b / a) for a, b in zip(breaks, breaks[1:])]
    return tuple(breaks) + (150.0,), tuple(slopes) + (IMF_HIGH_SLOPE,)


IMF_BREAKS_SOL, IMF_SLOPES = _bell_rows()
"""tuple: The segment edges (solar masses) and each segment's `alpha` (`dN/dM ~ M^-alpha`) of the single-star function."""

# The default star-formation history: ages uniform over the thin disk's
# 0-10 Gy. Population-dependent ages by position come on top of this (the
# sector fill passes an age instead of letting the star draw one).
STAR_FORMATION_AGE_RANGE_GY = (0.0, 10.0)

# Post-main-sequence phases (physics/post_main_sequence.py, GEN.935). The main-sequence lifetime is
# physics/ms_lifetime.py's. A star lighter than HELIUM_FLASH_MASS_SOL burns helium in its core after the red giant branch
# (the horizontal branch) and then climbs the asymptotic giant branch before it becomes a white dwarf (GEN.286). For it
# the subgiant phase ends at, and the red giant branch ends at, the times of physics/giant_branch.py (GEN.502, Hurley,
# Pols and Tout 2000: 1.05 and 1.12 of the main-sequence lifetime at a solar mass, from the core-mass-luminosity
# relation of the branch), and HORIZONTAL_BRANCH_DURATION_GY and AGB_MS_FRACTION are added for it. A star from
# HELIUM_FLASH_MASS_SOL up to SUPERGIANT_MIN_MASS_SOL has the same two times (0.6 and 1.2% of the main-sequence
# lifetime at 3 solar masses), then burns helium in its core and climbs the asymptotic giant branch for as long as the
# MIST track (Choi et al. 2016) does, as a share of the main-sequence lifetime (24% and 2% at 3 solar masses, 14% and 1%
# at 5); a star of SUPERGIANT_MIN_MASS_SOL or more has one supergiant phase as long as the MIST track's, from the end
# of the main sequence to the end of core helium burning (about 11% at 15 solar masses). The luminosity, temperature and
# radius in every phase are the MIST track's. [V]
HELIUM_FLASH_MASS_SOL = 2.2
# Core helium burning lasts the same 100 Myr in every star that ignites it in a
# flash, since the degenerate core is about 0.5 solar masses whatever the star's
# mass (Salaris et al. 2004, the 1.44 +- 0.06 horizontal-branch stars to each red
# giant brighter than the branch in globular clusters, against the 69 Myr the
# giant branch of physics/giant_branch.py spends there) [R].
HORIZONTAL_BRANCH_DURATION_GY = 0.100
AGB_MS_FRACTION = 0.002

# Stars at or above this initial mass end as supergiants and then neutron
# stars or black holes (core collapse); below it, as giants and white dwarfs.
SUPERGIANT_MIN_MASS_SOL = _boundaries.sol("core_collapse_limit")
# Giants at or above this initial mass are bright giants (II), not III.
BRIGHT_GIANT_MIN_MASS_SOL = 4.0
# A giant or supergiant on its track is a classical Cepheid inside the instability strip (a luminous star between a red
# edge near 5,000 K and a blue edge near 7,000 K), at least this luminous (the shortest-period classical Cepheids,
# about one day, are M_V near -2, a few hundred to a thousand Lsun) and born with at most this mass (classical Cepheids
# are 3 to 13 solar masses; the yellow stars of 20 to 45 are supergiants and hypergiants, not pulsators of that kind)
# [R, recalled, unconfirmed: the instability strip of Bono et al. 2000 and Turner 2010; the star class audit of
# 2026-10-29 section 1, H1]. The horizontal branch's own RR Lyrae stars have their colour
# (HORIZONTAL_BRANCH_COLOUR_TEMPERATURE_RANGE_K).
CEPHEID_TEMPERATURE_RANGE_K = (5000.0, 7000.0)
CEPHEID_MIN_LUMINOSITY_SOL = 500.0
CEPHEID_MAX_INITIAL_MASS_SOL = 13.0
# The mass a star has lost by a point on its track (GEN.940). MIST (Choi et al. 2016, section 3.7.1) takes Reimers (1975)
# mass loss on the red giant branch with eta_R = 0.1, tuned to the white dwarf initial-final mass relation, and a star
# that loses that little is 0.954 solar masses at the tip from 1; globular clusters want a horizontal branch of 0.55 to
# 0.7. McDonald and Zijlstra 2015 (MNRAS 446, 2226) fit eta_R = 0.477 +- 0.070 to the horizontal-branch morphology of 22
# clusters, so the loss of the first ascent is the track's times their eta over MIST's. Above SUPERGIANT_MIN_MASS_SOL the
# track's winds (Vink et al., Nugis and Lamers) are kept as they are.
RED_GIANT_BRANCH_REIMERS_ETA = 0.477
MIST_REIMERS_ETA = 0.1
# The horizontal branch (core helium burning): one luminosity range, whatever
# the colour; the colour (temperature) is drawn from [Fe/H] and age
# (HORIZONTAL_BRANCH_COLOURS below). Red clump and red horizontal-branch stars
# are the metal-rich, RR Lyrae the middle, blue horizontal-branch the
# metal-poor and extreme horizontal-branch (hot subdwarfs) a few percent of
# old stars (docs/design/galaxy/realistic-galactic-core.md 5.3) [R, recalled, unconfirmed:
# Lee, Demarque and Zinn 1994; Catelan 2009]. The luminosity is the RR Lyrae
# level, 36 to 47 Lsun (M_V 0.56 to 0.80 at [Fe/H] -1.6, Catelan 2009; Gould and
# Popowski 1998) and the red clump's 30 to 40 [P, stellar tables check]; a hot
# subdwarf is dimmer, 10 to 40 Lsun, and skips the AGB, which is not modelled
# (a few percent of old stars). The strip is 6,100 to 7,500 K (blue edge 7,400
# to 7,600 K, red edge 6,100 to 6,250 K); red horizontal-branch stars fill the
# 4,600 to 6,100 K below it; 20,000 K is the usual extreme boundary [P].
HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL = (35.0, 50.0)
HORIZONTAL_BRANCH_COLOUR_TEMPERATURE_RANGE_K = {
    "red": (4600.0, 6100.0), "rr-lyrae": (6100.0, 7500.0),
    "blue": (7500.0, 20000.0), "extreme": (20000.0, 35000.0),
}
# Colour shares (blue, RR Lyrae, red, extreme) at [Fe/H] anchors, linear
# between them and held beyond the ends: nearly all blue at -2.2, a mix at
# -1.5 (the intermediate halo), nearly all red above -0.5. The anchors are on
# the Zinn and West (1984) scale, the one the 57 globular clusters of Salaris et
# al. 2004, table 2, are listed on: their mean HB type, (blue - red) / all, is
# +0.75 at -2.2, +0.50 at -1.5 and -0.92 at -0.7 (the Carretta and Gratton 1997
# scale, 0.2 dex higher at -1.5, gives +0.72 there). The RR Lyrae and extreme
# shares and the -0.1 row are unverified [R].
HORIZONTAL_BRANCH_COLOURS = (
    (-2.2, (0.82, 0.07, 0.02, 0.09)),
    (-1.5, (0.55, 0.25, 0.13, 0.07)),
    (-0.7, (0.06, 0.06, 0.86, 0.02)),
    (-0.1, (0.0, 0.0, 0.98, 0.02)),
)
# The second parameter: a younger population has a redder branch, as if its
# [Fe/H] were higher by this many dex for each Gy it is younger than the
# reference age. No published figure was found; the papers imply 0.12 to 0.23
# dex per Gy indirectly (Lee et al. 1999, 2001; Catelan 2009) [R].
HORIZONTAL_BRANCH_AGE_SHIFT_DEX_PER_GY = 0.1
HORIZONTAL_BRANCH_REFERENCE_AGE_GY = 12.0
# Supergiant Yerkes class by luminosity (Lsun, lowest first): below the
# first threshold IB, then IAB, IA, and 0 (hypergiant) at the top.
SUPERGIANT_YERKES_THRESHOLDS_SOL = {"IAB": 10000.0, "IA": 50000.0, "0": 500000.0}

# Wolf-Rayet stars (GEN.243, docs/design/galaxy/intergalactic-and-halo-contents.md 5.11): a very massive star
# stripped to its hot helium core in the last stretch of its life, in place of the end of its supergiant phase.
WR_MIN_INITIAL_MASS_SOL = 25.0
"""float: The least initial mass of a star that ends as a Wolf-Rayet star at the Sun's metallicity, solar masses
(the Galactic Wolf-Rayet census against the Galactic birth rate gives about 25 [R]; the count check in
`tests/test_wolf_rayet.py` is the test of it)."""
WR_THRESHOLD_METALLICITY_EXPONENT = 0.5
"""float: The threshold rises as Z^-0.5: the winds that strip the star weaken with metals [R], so 25 solar masses at
the Sun's metallicity is 35 at half and 56 at a fifth of it [I]."""
WR_PHASE_MYR = 0.3
"""float: How long a Wolf-Rayet star stays one, Myr (the last stretch of its life) [R]."""
WR_MAX_SUPERGIANT_SHARE = 0.5
"""float: The phase never takes more than this share of the supergiant window it replaces, so the heaviest stars,
whose whole post-main-sequence life is short, still pass through a supergiant phase first."""
WR_SUBTYPE_PHASE_SHARES = (("WN", 0.6), ("WC", 0.3), ("WO", 0.1))
"""tuple: The subtype by the share of the phase gone: nitrogen-rich first, then carbon-rich, then the oxygen-rich last
tenth (only for stars of `WO_MIN_INITIAL_MASS_SOL` or more; lighter ones stay WC to the end)."""
WO_MIN_INITIAL_MASS_SOL = 40.0
"""float: The least initial mass that reaches the oxygen-rich WO stage [I]."""
WR_TEMPERATURE_RANGE_K = {"WN": (30000.0, 100000.0), "WC": (40000.0, 120000.0), "WO": (150000.0, 200000.0)}
"""dict: The effective temperature of each subtype, K, drawn uniformly: the whole 30,000 to 200,000 K of the design
note, split by subtype [R, I]."""
WR_SUBCLASS_RANGE = {"WN": (2, 9), "WC": (4, 9), "WO": (1, 4)}
"""dict: The subclass numbers of each subtype, hottest first (WN2 to WN9, WC4 to WC9, WO1 to WO4)."""
WR_LUMINOSITY_RANGE_SOL = (1.0e5, 3.0e6)
"""tuple: The range of the luminosity of a Wolf-Rayet star, solar luminosities: the helium-star relation of its present mass
(`wolf_rayet.state`, GEN.939) cut at both ends. It was 1e5 to 1e6 for a star of any mass; a present mass of 40 gives
2e6 [R]."""
WR_COMPANION_CHANCE = {"WN": 0.52, "WC": 0.90, "WO": 0.90}
"""dict: The share of Wolf-Rayet stars with a companion, by subtype (Dsilva et al., a spectroscopic multiplicity survey of
Galactic Wolf-Rayet stars: 52% of WN, 96% of WC, capped at 0.9 as in the phenomena catalogue; WO takes WC's)."""

# The mass range the MIST main-sequence grid holds (GEN.652): a mass outside it is carried out with the power laws at
# its edge (`physics/main_sequence.py`). The generator's own piecewise luminosity and radius laws that stood here are
# gone: a main-sequence star's luminosity, temperature and radius are the MIST row of its mass, age and [Fe/H].
MAIN_SEQUENCE_GRID_MASS_RANGE_SOL = (0.1, 240.0)
# The "map luminosity" relation (GEN.200): the luminosity a main-sequence star
# of a given mass would have, used to rank a faint but massive object (a
# neutron star, a black hole) for the Galaxy Map by the star of its mass
# (Astronomy Notes, "Mass-Luminosity Relation Explained"; Varsity Tutors,
# "Stellar Mass & Lifetime"). The star generators no longer
# read a piecewise law (GEN.652); this is the map's alone, with Boss's table's constants and
# breaks (1.5 M^3.5 from 2 to 20 Msun, not 1.4 M^3.5 to 55), and above 20 Msun
# L is proportional to M, its coefficient fixed by continuity with the band
# below (1.5 * 20^2.5, about 2,683 Lsun per Msun). So 20 Msun is about
# 53,700 Lsun and a 4 million Msun black hole about 1e10 Lsun.
MAP_MASS_LUMINOSITY_PIECES = (
    {"max_mass_sol": 0.43, "coeff": 0.23, "exponent": 2.3},
    {"max_mass_sol": 2.0, "coeff": 1.0, "exponent": 4.0},
    {"max_mass_sol": 20.0, "coeff": 1.5, "exponent": 3.5},
    {"max_mass_sol": None, "coeff": 1.5 * 20.0 ** 2.5, "exponent": 1.0},
)
# MAP.208: a black hole, neutron star or quasar is drawn on the maps as bright
# as a star of its mass (map_luminosity_sol), a black hole this share of it
# (a little dimmer, so its purple does not out-shine a star); an object with no
# mass on record is drawn as a star of this mass.
BLACK_HOLE_MAP_BRIGHTNESS = 0.8
MASSLESS_MARKER_MASS_SOL = 1.0
# The Sun's effective temperature (IAU 2015 nominal), anchoring
# Stefan-Boltzmann in solar units: T = T_sun * (L / R^2)^(1/4).
SUN_EFFECTIVE_TEMPERATURE_K = 5772.0

# White dwarfs: final mass from the initial-final mass relation of
# Cummings et al. (2018), ApJ 866:21 (arXiv 1809.01673), the MIST-based fit, equations 4 to 6, read; luminosity from Mestel
# cooling, L = WD_COOLING_L0 * (M/0.6) * t_cool^-1.4 (Lsun, t_cool in Gy), clamped.
WD_IFMR_PIECES = (
    (0.83, 2.85, 0.080, 0.489),
    (2.85, 3.60, 0.187, 0.184),
    (3.60, 7.20, 0.107, 0.471),
)
"""tuple: The white dwarf initial-final mass relation, `(low, high, slope, intercept)` per piece: the final mass is
`intercept + slope * M_i` solar masses for an initial mass between `low` and `high` (Cummings et al. 2018, equations 4, 5
and 6, the fit to the MIST tracks; the pieces join at 2.85 and 3.60). The first piece is carried down to the lightest
progenitor (`WHITE_DWARF_PROGENITOR_SOL`) and the last up to the core collapse limit. Verified (read in the paper); one
relation, read by the star model, the member draw, the disk history and `mass_loss`. The relation of Kalirai et al. (2008)
the generator used before runs 0.04 to 0.06 under it, and Weidemann (2000), the disk history's, 0.1 to 0.16 under at 5 to 7
solar masses (the audit of the star classes, M2)."""
WD_MASS_RANGE_SOL = (0.5, 1.35)
WD_COOLING_L0_SOL = 1.0e-3
WD_COOLING_EXPONENT = -1.4
WD_LUMINOSITY_RANGE_SOL = (1.0e-5, 50.0)
WD_MIN_COOLING_AGE_GY = 1.0e-4

# `+large_star` draws from the same model with the IMF truncated below
# this mass, and an age inside the star's own lifetime so it's still alive.
LARGE_STAR_MIN_MASS_SOL = 1.4

# A random primary that has already collapsed (a neutron star or black
# hole) is redrawn: isolated remnants are generated as phenomena, so their
# rates stay in one place. Caps the redraws.
STAR_MODEL_MAX_REDRAWS = 1000

# Planets by the star's age and history. A star younger than
# PLANET_MIN_STAR_AGE_GY still has only a debris disk (planet formation
# takes ~10 Myr), so it keeps belts at most; no habitable-class world is
# drawn around a star younger than LIFE_MIN_STAR_AGE_GY (a crust and
# oceans take ~0.1 Gy). A giant engulfs planets inside
# GIANT_ENGULFMENT_RADIUS_FACTOR times its present radius (tidal capture
# reaches a few stellar radii), and a white dwarf's progenitor, on its
# asymptotic giant branch, engulfed or drove off everything inside
# the distance a lookup through WD_ENGULFMENT_POINTS_AU gives for a progenitor
# of M solar masses (GEN.1012, GM1/GM2: the three Jovian-planet points of
# Mustill and Villaver 2012, ApJ 761, 121, Figure 7, with their Zahn tides:
# 2.7 AU for a solar-mass star, 3.8 AU at 2.5 and 10 AU at 5; straight lines
# between them, the last slope carried on to 8). A white dwarf without a
# model progenitor mass (a forced type, a planetary nebula's central star)
# takes the mass that leaves it
# (WHITE_DWARF_REMNANT_MASS_SOL inverted), kept inside
# WD_PROGENITOR_MASS_RANGE_SOL: nothing under a solar mass has had time to
# become one, and above 8 the relation hands over to neutron stars.
PLANET_MIN_STAR_AGE_GY = 0.01
LIFE_MIN_STAR_AGE_GY = 0.1
GIANT_ENGULFMENT_RADIUS_FACTOR = 2.0
WD_ENGULFMENT_POINTS_AU = ((1.0, 2.7), (2.5, 3.8), (5.0, 10.0))
"""tuple: `(progenitor solar masses, AU)` points of the engulfment radius, Mustill and Villaver 2012 (Jovian planets). Source tier T3, tweak level L1 (GM1, GM2)."""
WD_PROGENITOR_MASS_RANGE_SOL = (1.0, _boundaries.sol("core_collapse_limit"))

# A planetary nebula's central star is a close pair in PLANETARY_NEBULA_COMPANION_CHANCE of nebulae (a post-common-envelope
# pair); such a pair still inside its ejecta holds a circumbinary disk and no planets (GEN.318, design note 1.3). A white
# dwarf at least this hot (the coolest central star type, B0VII) belongs to such a pair; the disk's mass is a stand-in
# bracket, the two figures the note gives (a disk approaching 0.1 solar masses can trigger the streaming instability, one of 0.01
# cannot), and its outer edge the young stars' range, since these disks look like young protoplanetary ones.
PLANETARY_NEBULA_CENTRAL_MIN_TEMPERATURE_K = 30000.0
POST_AGB_DISK_MASS_RANGE_SOLAR = (0.01, 0.1)

# The pre-main-sequence state (GEN.242, GEN.529, physics/pre_main_sequence.py).
# A star below PRE_MAIN_SEQUENCE_MAX_MASS_SOL younger than its contraction time
# is a protostar (Class 0 and I, under PROTOSTAR_MAX_AGE_YR) or a
# pre-main-sequence star (Class II and III), however old that time is: the
# track takes 0.1 solar masses about a billion years, and a red dwarf still
# shrinking at 100 Myr keeps its pre-main-sequence class. The track is
# derived: a convective star (polytrope n = 3/2) holds 3/7 G M^2 / R of
# energy and a radiative one (n = 3) 3/4 G M^2 / R, it shrinks at the
# luminosity it radiates, and on its Hayashi line (convective) the effective
# temperature stays put, so R^-3 grows in step with time. Three numbers are
# read from the published tracks, not recalled (GEN.529): Baraffe, Homeier,
# Allard and Chabrier 2015 (A&A 577, A42, arXiv:1503.04107, the BHAC15 tracks
# table at perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir, solar metallicity,
# 0.1 to 1.4 solar masses, from 0.5 Myr) [V]. The Hayashi temperature is the
# effective temperature at 1 Myr, 4357 K at one solar mass and rising as
# m^0.1815 (18 masses, rms 0.004 dex); the birthline radius is the radius at
# the table's first age, 3.19 solar radii at one solar mass and rising as
# m^0.4335 (rms 0.03 dex); both are extended past 1.4 solar masses by their
# power law. The energy factor is the one number the table fixes only through
# its end: contraction alone, at the luminosity the star radiates, is too
# quick, and deuterium burning and the heat of the interior pay for part of
# the light, so the time is stretched by 2.0, the value at which the model
# reaches the main sequence at the age the table's radius settles to within
# 2% of its minimum, rms 0.12 dex over the 7 masses of the test (every one within
# 0.19 dex; 3.0 before GEN.652, when the model ended on the generator's own,
# brighter and larger main sequence). The model now ends on the MIST main
# sequence's first row at the star's mass and [Fe/H] (`physics/main_sequence.py`,
# GEN.652), the same table the star's own record reads, so the track and the star
# it becomes agree. At 1 to 10 Myr the luminosity is within 0.25 dex of the table up
# to 0.3 solar masses and up to 0.35 dex too bright from 0.5 to 1.4, where the
# derived Hayashi phase falls more slowly than the table's.
PRE_MAIN_SEQUENCE_MAX_MASS_SOL = _boundaries.sol("core_collapse_limit")
# No star contracts longer than this: the track takes the lightest star (0.074 solar masses, the hydrogen-burning limit at solar
# metallicity, GEN.936) 4.1 Gyr, against 0.94 Gyr for the 0.08 of before (GEN.652; the main sequence under 0.1 solar masses now is
# Baraffe's, which arrives later). A shortcut for the old stars, which are nearly all of them; the test checks it against the track.
PRE_MAIN_SEQUENCE_LONGEST_YR = 4.5e9
PROTOSTAR_MAX_AGE_YR = 0.5e6
HAYASHI_TEMPERATURE_1_SOL_K = 4357.0
HAYASHI_TEMPERATURE_MASS_EXPONENT = 0.1815
BIRTHLINE_RADIUS_1_SOL_RSUN = 3.187
BIRTHLINE_RADIUS_MASS_EXPONENT = 0.4335
CONVECTIVE_BINDING_COEFFICIENT = 3.0 / 7.0
RADIATIVE_BINDING_COEFFICIENT = 3.0 / 4.0
PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR = 2.0
# The disk (VANDAM and the TODO text, GEN.242): Class 0 under 0.15 Myr, I under
# 0.5, II under 3, and none after; mass 0.01 to 0.3 of the star and outer
# radius 10 to 300 AU, each log-uniform; the inner edge is where a blackbody
# grain sublimates (DISK_DUST_SUBLIMATION_K).
DISK_CLASS_AGES_YR = (("0", 0.15e6), ("I", 0.5e6), ("II", 3.0e6))
RUNAWAY_HISTORY_START_YR = 30.0e6
"""float: The age the water-loss history starts at (GEN.639): 30 Myr, when the last giant impacts of the rocky planets are
done (Earth's Moon-forming impact, 30 to 100 Myr after the disk; Kleine et al. 2009, Jacobson et al. 2014) and a surface
with its water is bare to the star. A planet that was still growing earlier had no oceans to lose; the choice keeps the
Sun's own habitable zone dry-proof (the young Sun was past the runaway limit at 1 AU only until about 15 Myr). DESIGN."""

RUNAWAY_HISTORY_STEPS = 64
"""int: Log-spaced steps of the water-loss history integral from `RUNAWAY_HISTORY_START_YR` to the star's age."""

ENERGY_LIMITED_EFFICIENCY = 0.01
"""float: The share of the absorbed XUV that drives hydrogen away in the energy-limited escape,
Mdot = eps pi R^3 F_xuv / (G M) (Watson et al. 1981; Luger and Barnes 2015 take 0.1, other estimates run from about 0.01 to
0.3). 0.01, a tenth of the usual value, so that only the worst red dwarf zones lose their oceans (Boss, 7:09 PM EDT
2026-10-23; DESIGN stand-in, GEN.639)."""

MANTLE_OUTGAS_SHARE = 0.3
"""float: The share of a planet's mantle water that volcanism brings to the surface after an early runaway has stripped the
surface reservoir, because escape reaches only the air and the oceans, not the mantle (GEN.639; DESIGN stand-in)."""

AQUIFER_CONDUCTIVITY_W_M_K = 2.0
"""float: The thermal conductivity of a world's crust, W/m/K, that turns its heat flow into a geotherm (Clifford 1993, JGR 98,
10973, column average about 2; secondary source, DESIGN stand-in, GEN.654)."""

AQUIFER_POROSITY_SURFACE = 0.2
"""float: The pore share of the crust at the surface (the lunar and Martian megaregolith value; DESIGN stand-in, GEN.654)."""

AQUIFER_POROSITY_MIN = 0.01
"""float: The pore share under which the crust holds no usable water, so the water-bearing layer ends (DESIGN stand-in, GEN.654)."""

AQUIFER_POROSITY_SCALE_KM_MARS = 2.82
"""float: The e-folding depth of the crust's porosity, phi(z) = phi0 exp(-z/K), on Mars (Clifford 1993; secondary source). On
another world it scales with the inverse of the gravity (DESIGN stand-in, GEN.654)."""

AQUIFER_MARS_GRAVITY_M_S2 = 3.71
"""float: The gravity `AQUIFER_POROSITY_SCALE_KM_MARS` is for."""

AQUIFER_LIFE_MAX_K = 395.15
"""float: The hottest water life can use, 122 C (Takai et al. 2008, PNAS 105, 10949, Methanopyrus kandleri); the water-bearing
layer stops where the geotherm reaches it."""

AQUIFER_MIN_THICKNESS_KM = 0.01
"""float: A water-bearing layer thinner than this (10 m) is none."""

AQUIFER_REACH_KM = 5.0
"""float: How deep a base can drill for water (the deepest mines are about 4 km; DESIGN stand-in, GEN.654): an aquifer with its
top above this depth is a water source for the Water and Suitability rungs, a deeper one only shelters microbes."""

WATER_PER_HYDROGEN_MASS = 8.94
"""float: Kilograms of water gone for each kilogram of hydrogen that escapes (H2O 18.015 g/mol over its two hydrogens
2.016 g/mol; the oxygen is dragged along or left behind to oxidise the crust, and neither returns the water)."""

DISK_MASS_FRACTION_RANGE = (0.01, 0.3)
DISK_OUTER_RADIUS_RANGE_AU = (10.0, 300.0)
DISK_DUST_SUBLIMATION_K = 1500.0
# The cold debris disk of a main-sequence star (GEN.787, physics/debris_disk.py). Sibthorpe et al. 2018 (MNRAS 475, 3046, the
# Herschel DEBRIS survey; arXiv:1803.00072, sections 2.5 and 3 and Tables 1 and 2 read 2026-10-30): the share of stars with a
# detected cold disk is F 23.9, G 14.3 and K 13.0 percent (47 discs round 275), M 2.2 percent (89 stars) and, raw, 24 percent for
# A stars. `DEBRIS_DISK_INCIDENCE` is the share for each mass band, the bands ending at the masses of
# `DEBRIS_DISK_INCIDENCE_MASS_LIMITS_SOL` (the last band has no upper mass). The masses are the spectral types'
# (M under 0.45, K to 0.8, G to 1.04, F to 1.4, A and hotter above); Matthews et al. 2014 find the dust shows no clear dependence on
# type from late B to M, so a hotter star takes the A figure.
DEBRIS_DISK_INCIDENCE_MASS_LIMITS_SOL = (0.45, 0.8, 1.04, 1.4)
DEBRIS_DISK_INCIDENCE = (0.022, 0.130, 0.143, 0.239, 0.24)
BLACKBODY_AT_ONE_AU_K = 278.3
# The dust temperature of a detected disk: the 47 discs of Sibthorpe 2018 Table 1 have blackbody radii from 0.63 to 280 AU with a
# median of 23 AU (7 to 40 AU for most), which at one solar luminosity is a median of 58 K; the log-normal width 0.45 puts the
# middle two thirds of Sun-like stars at about 7 to 40 AU (a stand-in fitted to that statement, not a fit to the table).
DEBRIS_DISK_TEMPERATURE_MEDIAN_K = 58.0
DEBRIS_DISK_TEMPERATURE_LOG_SIGMA = 0.45
# The belt lies further out than a blackbody grain would: the paper's 7 to 40 AU blackbody radii "scale to about 40 AU for realistic
# dust properties", about 1.7 times the median 23 AU. The edges at 0.7 and 1.5 of the belt radius are a design choice (a ring about
# as wide as it is far out), not read from a source.
DEBRIS_RADIUS_OVER_BLACKBODY = 1.74
DEBRIS_DISK_EDGES_OVER_RADIUS = (0.7, 1.5)
# Companions of a protostar: the multiplicity fraction of Perseus' Class 0
# (0.57) and Class I (0.23) protostars, Tobin et al. 2016 (ApJ 818:73,
# arXiv:1601.00692, abstract read 2026-10-17).
PROTOSTAR_MULTIPLICITY_BY_CLASS = {"0": 0.57, "I": 0.23}

# Stellar populations by age (galaxy.density.population_densities,
# star_population). Each has its own age range (uniform within it); the
# three disk populations share the disk's stars in proportion to the star
# formation that fell in their age range (GEN.469, galaxy.disk_history: the
# disk's own gas history; it was the range's length over
# STAR_FORMATION_AGE_RANGE_GY's when the history was flat). Young stars sit close to the plane and
# crowd the spiral arms; old stars are puffed up by billions of years of
# scattering and spread evenly in azimuth. Scale heights are relative to
# the galaxy's own thin-disk scale height, all exponential scale heights
# (`galaxy.density._vertical`): the old thin disk's ~300 pc in the Milky
# Way (BHG16), young ~50 pc (OB stars, Reed 2000) and intermediate
# ~150 pc (A to F dwarfs, Bovy 2017). The bulge is old
# (8-12 Gy), and so is the thick disk, whose stars count as "old"
# (`galaxy.density.population_densities`). The halo and the globular clusters are older still, 10 to 13 Gy
# (GEN.234). A satellite galaxy's stars (GEN.237) are three populations, split by age like the disk's and by the same windows
# (GEN.526): `satellite_young` and `satellite_intermediate` are the shares its star formation history puts in the last
# 0.1 and 3 Gy (`galaxy.satellite_history`; 0 for a satellite that no longer forms stars and holds none in its
# history's last 3 Gy), and `satellite` is the rest, 3 Gy to the age of the Universe. The nuclear
# populations (GEN.280) are windows with the mean ages of the observed mixes (the cluster 80% 5 to 12 Gyr and 20%
# 0.1 to 5, a mean of 7.3 Gyr; the disc 70% 8 to 12 and 30% 0 to 8, a mean of 8.2) [recalled, to verify].
THIN_DISC_OLD_WINDOW_GY = (3.0, 10.0)
"""tuple: The thin disc's share of the `old` population's ages (GEN.1025): the thin disc's own history ends at 10 Gy, so its law, its
share of the disc, its heating and its arm contrast keep this window while `old` itself reaches the thick disc's 13."""
STELLAR_POPULATION_AGE_RANGES_GY = {
    "young": (0.0, 0.1), "intermediate": (0.1, 3.0), "old": (3.0, 13.0), "bulge": (8.0, 12.0),
    "halo": (10.0, 13.0), "cluster": (10.0, 13.0), "satellite": (3.0, 13.0),
    "nuclear_cluster": (3.0, 12.0), "nuclear_disc": (4.0, 12.0),
    "satellite_young": (0.0, 0.1), "satellite_intermediate": (0.1, 3.0),
}
# The young open clusters (GEN.244, `generation/open_clusters.py`) are coeval, so their stars cannot share one age
# window: a cluster's stars belong to the population of the quarter-decade of age it falls in, 1 Myr to 1 Gyr
# (`open_cluster_00` is 1 to 1.78 Myr, `open_cluster_11` 562 Myr to 1 Gyr), and a star's age is drawn from that window.
OPEN_CLUSTER_AGE_EDGES_MYR = tuple(10.0 ** (k / 4.0) for k in range(13))
OPEN_CLUSTER_POPULATIONS = tuple(f"open_cluster_{k:02d}" for k in range(len(OPEN_CLUSTER_AGE_EDGES_MYR) - 1))
STELLAR_POPULATION_AGE_RANGES_GY.update({
    name: (_myr_to_gyr(OPEN_CLUSTER_AGE_EDGES_MYR[k]), _myr_to_gyr(OPEN_CLUSTER_AGE_EDGES_MYR[k + 1]))
    for k, name in enumerate(OPEN_CLUSTER_POPULATIONS)})
# The three thin-disc populations' scale heights over the old one's come from the galaxy's own dispersions
# (`galaxy.structure.population_height_ratios`, GEN.532), not from a table here.
DISK_AGE_LAW_POPULATIONS = ("young", "intermediate", "old")
"""tuple: The populations whose ages inside their window follow the disk's star formation rate (GEN.546,
`galaxy.disk_history.population_age_laws`). The others (bulge, halo, clusters, satellites) draw uniformly: their history
is not the disk's."""
AGE_LAW_POINTS = 33
"""int: Points of the rate curve an age law keeps across a window (the curve is smooth, and 33 points hold a window's
cumulative to better than 1e-3 of the window)."""

# How much of the thin disk's recent star formation a place holds (GEN.134, `galaxy.star_formation`): derived from the
# galaxy's shape, with the Milky Way's profile only as a test target.
STAR_FORMATION_LAW_EXPONENT = 1.4
"""float: The Kennicutt-Schmidt exponent: the rate per area goes as the gas per area to this power (Kennicutt 1998:
1.4 +/- 0.15 over 61 spirals and 36 starbursts; other fits run from 1.3 to 1.65)."""
ARM_GAS_TO_INTERARM_RATIO = 4.0
"""float: The molecular gas on a spiral arm's crest over the gas midway between arms (CO arm to interarm contrasts run about 2
to 5; Koda et al. 2009, ApJ 700:L132). `NEBULA_FIELD_ARM_AMPLITUDE` is this ratio's amplitude, `(r - 1) / (r + 1)` (GEN.931)."""
ARM_CONTRAST_OLD_DISK = 0.10
"""float: The arm contrast, `1 + A cos(arm phase)`, of the old disk's stars: the surface mass density contrast between arm and
inter-arm is about +-10% at the solar radius (Eilers et al. 2020, ApJ 900:186, red giants, abstract), 14% in Siebert et al.
2012, and the near-infrared arm to inter-arm ratio is 1.2 in J and 1.32 in K (Drimmel and Spergel 2001, ApJ 556:181, section
4.3), amplitudes 0.09 to 0.14 (GEN.931)."""
ARM_CONTRAST_AGE_YOUNG_GY = 0.05
"""float: The age the young population's arm contrast belongs to: the middle of its window, 0 to 0.1 Gy (GEN.931)."""


def _arm_amplitude(ratio):
    return (ratio - 1.0) / (ratio + 1.0)


def _stellar_population_arm_amplitudes():
    """The arm contrast of each thin-disk population (GEN.931). Young stars (under 0.1 Gy) are where they were born, so they
    follow the star formation, which goes as the gas to the Kennicutt-Schmidt index: the gas ratio to that power. The old
    disk is the observed one. Between, the contrast falls with the log of the age from the young value to the old one,
    at each population's geometric mean age (a stand-in: stars leave their arm in a few hundred Myr, and the observed
    contrasts of A and F stars, between the two, are not published as a number)."""
    young = _arm_amplitude(ARM_GAS_TO_INTERARM_RATIO ** STAR_FORMATION_LAW_EXPONENT)
    old = ARM_CONTRAST_OLD_DISK
    windows = {"intermediate": STELLAR_POPULATION_AGE_RANGES_GY["intermediate"], "old": THIN_DISC_OLD_WINDOW_GY}
    ages = {name: _math.sqrt(window[0] * window[1]) for name, window in windows.items()}
    span = _math.log(ages["old"] / ARM_CONTRAST_AGE_YOUNG_GY)
    fraction = _math.log(ages["intermediate"] / ARM_CONTRAST_AGE_YOUNG_GY) / span
    return {"young": young, "intermediate": young + (old - young) * fraction, "old": old}


STELLAR_POPULATION_ARM_AMPLITUDE = _stellar_population_arm_amplitudes()
"""dict: Each thin-disk population's arm contrast `A` in `1 + A cos(arm phase)`: about 0.75, 0.42 and 0.10 for the young,
intermediate and old (GEN.931). Derived from `ARM_GAS_TO_INTERARM_RATIO`, `STAR_FORMATION_LAW_EXPONENT` and
`ARM_CONTRAST_OLD_DISK`; the shares of a place's mix shift, never its total density."""
STAR_FORMATION_LOW_DENSITY_EXPONENT = 2.1
"""float: The index the star-formation law steepens to where the gas thins (GEN.450): resolved observations give 1.5 to 2.7 for
the total gas (Kennicutt and Evans), and 2.1 is the middle of that range. Open: no source gives the index of the outer
discs of a single galaxy (Bigiel et al. 2008 find the total-gas relation weak where atomic gas dominates)."""
STAR_FORMATION_SATURATION_GAS_MSUN_PC2 = 9.0
"""float: The gas surface density, solar masses per pc2, below which the index rises (GEN.450): atomic gas saturates near
9 solar masses per pc2 (Bigiel et al. 2008), so the gas that forms stars at a constant efficiency, the molecular gas,
is what lies above it; read as a 'half-way' scale of the smooth rise, not a step."""
GAS_TO_STELLAR_SCALE_LENGTH_RATIO = 1.0
"""float: The gas disk's scale length over the stellar disk's. An assumption, not a measurement: the observed gas
disk is more extended, which would make the outer disk younger; 1 is the plain reading, and the profile test shows
what it gives against the observed fall."""
BAR_GAP_RADIUS_RATIO = 2.2
"""float: The radius inside which the bar has depleted the star-forming gas, over the bulge's scale radius (the bar's
length scale). The molecular gas has a gap from 0.5 to 3.5 kpc in a 2026 reconstruction whose outer edge is taken to match the
bar's half-length (arXiv 2603.17341, read; another estimate of the half-length is 5.0 kpc), against the bulge's 1.58 kpc scale radius:
2.2 (a model constant of the shape, like `BULGE_AXIS_RATIO_Y`, not derived)."""
BAR_GAP_STAR_FORMATION_FLOOR = 0.03
"""float: What is left of the star formation at the centre of the bar gap, as a share of the unbarred rate. The
simulations remove the gas "on 1 Gyr timescales" after the bar forms (arXiv 1908.11119, read), so 0.03 is what 3.5 Gyr
leaves (`exp(-3.5)`); the observed deficit depends on the reconstruction (design note, 'the bar gap'), so this is a default [open]: below about 0.05
the rate per area peaks at the gap's edge, the ring a bar leaves."""
BAR_GAP_EDGE_POWER = 4.0
"""float: How sharp the gap's edge is: the rate rises as `u^p / (1 + u^p)` with `u` the radius over the gap radius.
The bar's gas is swept in over its length and what lies beyond is untouched, so the edge is steep; 2 gives a ramp
that never leaves a ring, 4 a ring at the bar's end [a model default, open]."""

# The thin disk's star formation history over time (GEN.469, `galaxy.disk_history`,
# docs/design/galaxy/disk-star-formation-history.md): gas continuity per ring, the inflow from the halo's growth, the
# Kennicutt-Schmidt law above as the depletion. The Milky Way's own history is only a test target.
# The cosmology is the mass chain's (`galaxy.mass_chain.HUBBLE_H`, `OMEGA_MATTER`: Planck-like, 13.8 Gy old).
HALO_GROWTH_MASS_EXPONENT = 1.1
"""float: How the mean growth rate of a halo's mass goes with its mass, `M^1.1` (Fakhouri, Ma and Boylan-Kolchin 2010,
MNRAS 406, 2267, arXiv 1001.2304, eq. 2, read: fitted to the Millennium simulations)."""
HALO_GROWTH_REDSHIFT_SLOPE = 1.11
"""float: The `1 + 1.11 z` in the same fit's mean growth rate, `46.1 Msun/yr (M / 1e12)^1.1 (1 + 1.11 z) E(z)`."""
STELLAR_TO_HALO_EVOLUTION = {"log_m1": 1.195, "norm": -0.0247, "beta": -0.826, "gamma": 0.329}
"""dict: How the stellar-to-halo mass relation moves with redshift: each of its four parameters is its present value
(`galaxy.mass_chain.LOG_M1`, `N_NORM`, `BETA`, `GAMMA`) plus this times `z / (z + 1)` (Moster, Naab and White 2013,
MNRAS 428, 3121, arXiv 1205.5807, eqs. 11 to 14 and Table 1, read: the multi-epoch abundance-matching fit with errors
in the observed masses). It gives the stars a halo of a given mass and age holds, so the disk's gas inflow
follows the stellar mass the halo's growth builds."""
HALO_GROWTH_RATE_AT_1E12_MSUN_YR = 46.1
"""float: That fit's amplitude, solar masses per year for a halo of 1e12 solar masses at z = 0. It only sets how the
halo's mass went back in time; the disk's inflow is scaled to the stars it must build, so its size drops out."""
NEUTRON_STAR_RETURN_REMNANT_SOL = 1.4
"""float: The remnant of a star of 8 to `BLACK_HOLE_RETURN_MASS_SOL` solar masses (the same source's choice)."""
BLACK_HOLE_RETURN_MASS_SOL = 40.0
"""float: Above this birth mass a star keeps all its mass, collapsing to a black hole (the same source's choice; few
stars form above it)."""
RECYCLING_MASS_SOL = 1.0
"""float: Stars above this mass return their gas within the history's step, stars below it never (the instantaneous
recycling approximation, Madau and Dickinson 2014, who take 1 solar mass; it fails only between 0.9 and 2.5)."""
DISK_HISTORY_RINGS = 40
"""int: The rings, over 0 to `DISK_HISTORY_MAX_RADIUS_SCALE_LENGTHS` scale lengths, each with its own gas history."""
DISK_HISTORY_MAX_RADIUS_SCALE_LENGTHS = 10.0
"""float: The outermost ring, in disk scale lengths (the disk holds 0.05% of its stars beyond)."""
DISK_HISTORY_TIME_STEPS = 800
"""int: Steps of the gas history over the disk's age (12.5 Myr each at 10 Gyr)."""

# Metallicity ([Fe/H], dex; GEN.158, docs/design/galaxy/globular-clusters.md 5). The
# disc populations follow a mild radial gradient of -0.06 dex per kpc [recalled,
# to verify; written as a fall, so positive] that is zero at the Sun's distance from the centre; each population
# has its own offset and scatter (the young disc is the most uniform, the bulge
# the widest). The halo is -1.6 inside 25 kpc and -2.2 outside, scatter 0.5
# (docs/design/galaxy/intergalactic-and-halo-contents.md 5.3). A cluster carries its own
# value. `METALLICITY_RANGE_FE_H` clamps every draw.
METALLICITY_DISC_FALL_DEX_PER_KPC = 0.06
METALLICITY_BY_POPULATION = {
    "young": {"offset": 0.0, "sigma": 0.10, "gradient": True},
    "intermediate": {"offset": 0.0, "sigma": 0.15, "gradient": True},
    # GEN.1025: the old population's ages run 3 to 13 Gyr now (the thick disc's 8 to 13 Gy added), and the relation's mix over that window averages 0.125 dex
    # lower than over 3 to 10 Gy (-0.255 against -0.130, `age_metallicity`); the offset was -0.2 and is -0.325 so the population's [Fe/H] is still the
    # relation's mix over its window and the age tilt does not crowd it toward the young. [Stand-in, tier T6, tweak L2: tuned to the relation.]
    "old": {"offset": -0.325, "sigma": 0.25, "gradient": True},
    "bulge": {"offset": -0.05, "sigma": 0.40, "gradient": False},
    # The nuclear cluster and disc are metal rich, +0.1 to +0.3 [recalled, GEN.280], and the cluster's spread is broad.
    "nuclear_cluster": {"offset": 0.20, "sigma": 0.40, "gradient": False},
    "nuclear_disc": {"offset": 0.20, "sigma": 0.30, "gradient": False},
}
METALLICITY_HALO = {"inner": -1.6, "outer": -2.2, "break_kpc": 25.0, "sigma": 0.5}
METALLICITY_RANGE_FE_H = (-3.5, 0.6)
# Solid surface density of a protoplanetary disk scales as 10^[Fe/H] (solids are
# a fixed share of the metals), and giant-planet occurrence as 10^(2 [Fe/H]),
# capped at 1 (Fischer and Valenti 2005 form, recalled).
GIANT_PLANET_METALLICITY_EXPONENT = 2.0
# GEN.236: the planet cull. The generator's planet counts are those of a disk-population star ([Fe/H] about -0.2 to 0),
# so occurrence is left alone at and above the floor (`PLANET_OCCURRENCE_FLOOR_DEX_BELOW_SOLAR` below the Sun) and falls
# below it: planets of any kind as
# `10^(PLANET_OCCURRENCE_METALLICITY_EXPONENT ([Fe/H] - floor))` (small planets keep a weak dependence, so a gentler
# slope than the giants'), giants as the Fischer and Valenti form above, taken from the floor too. So a metal-poor
# halo or cluster star hosts a system of one or two small planets at most and almost never a giant (the floor of -0.5
# is recalled as where Kepler's small-planet occurrence starts to fall, Buchhave 2012; not looked up this session).
PLANET_OCCURRENCE_FLOOR_DEX_BELOW_SOLAR = 0.5
"""float: Dex below the Sun's [Fe/H] where planet occurrence starts to fall (the floor is [Fe/H] -0.5), kept as a
positive depth."""
PLANET_OCCURRENCE_METALLICITY_EXPONENT = 1.0

# Stellar abundances (GEN.439, `physics/stellar_abundances.py`, docs/design/composition/stellar-composition.md and part A of
# docs/design/composition/body-composition-research.md). The reference scale is Asplund, Amarsi and Grevesse 2021 (arXiv:2105.01661, Table 2
# log eps on H = 12; Boss 2026-10-17 1:24 AM EDT), moved to proto-solar by A21's diffusion correction of +0.064 dex for the
# metals and +0.070 dex for helium; it is a scale and never a preset. Verification states are in `STELLAR_RELATION_STATES`.
STELLAR_SOLAR_LOG_EPS = {
    "C": 8.46, "N": 7.83, "O": 8.69, "Ne": 8.06, "Na": 6.22, "Mg": 7.55, "Al": 6.43, "Si": 7.51, "S": 7.12, "Ar": 6.38,
    "K": 5.07, "Ca": 6.30, "Fe": 7.46, "Ni": 6.20, "Eu": 0.52, "Th": 0.03, "U": -0.54}
"""dict: Element symbol to A21 photospheric log eps (H = 12); the 14 elements that make 99.7 percent of Z, and Eu, Th and U."""
STELLAR_PROTO_METAL_OFFSET_DEX = 0.064
STELLAR_PROTO_X = 0.7121
STELLAR_PROTO_Z = 0.0154
STELLAR_PRIMORDIAL_Y = 0.24672
"""float: Primordial helium mass fraction, Planck 2018 VI eq. 72a (arXiv:1807.06209). The proto-solar X, Z and Y (0.2725 by
difference) are A21 re-summed with the offsets above."""
STELLAR_THICK_START_GY = 8.0
"""float: Stars younger than this are thin-disc; older ones draw the branch from the logistic. Haywood et al. 2013 (the
thick disc formed from about 13 to about 8 Gy ago, the age scale uncertain by 1 to 1.5 Gy; Boss's papers); range 6.5 to
9.5. The first draft's 3.0 was the project's intermediate/old boundary and the paper contradicts it (Boss approved
8.0, 2026-10-17)."""

THICK_DISC_AGE_RANGE_GY = (STELLAR_THICK_START_GY, 13.0)
"""tuple: The ages, Gy, over which the thick disc formed its stars: from where its branch starts (`STELLAR_THICK_START_GY`) to the oldest the old
population holds (GEN.1025, `galaxy.thick_ages`). [Stand-in TD1, tier T3, tweak L0: Haywood et al. 2013 (the thick disc formed from about 13 to
about 8 Gy ago).]"""

STELLAR_HIGH_ALPHA_LOGISTIC = (-0.350, 0.132)
"""tuple: `(f0, w)`: P(high-alpha) = 1 / (1 + exp(([Fe/H] - f0) / w)) for a star old enough to be thick-disc (Bensby, Feltzing
and Oey 2014, fit to the chemically high-alpha fraction; their paper is not in the repository's papers, so partly verified)."""
STELLAR_FE_H_FIT_RANGE = {"thin": (-0.9, 0.5), "high-alpha": (-1.6, 0.1)}
"""dict: The [Fe/H] each branch's [X/Fe] fits were made over; outside it the fit is held at its edge (flagged)."""
# [X/Fe] means by branch, as ("linear", a, b): a + b F; ("min", c, b, k): c + b min(F - k, 0), the thin disc's plateau above k
# and rise below; ("max", p, s, k): p + s max(F - k, 0), the high-alpha plateau below k. Fitted to Bensby 2014 (Mg, Si, Ca,
# Al), Nissen 2014 NLTE (C, O) and Battistini and Bensby 2016 (Eu), part A 2.7. The thin-disc Mg and Si are refitted to Brewer and
# Fischer 2016 (GEN.989, and the thin-disc O of GEN.991; `scripts/fit_stellar_xfe_refit.py`): their stars sit 0.06 dex under the Sun in
# Mg, 0.03 in Si and 0.03 in O at [Fe/H] above solar, where Bensby's plateau is at the Sun and Nissen's O line goes on falling.
STELLAR_XFE_FITS = {
    "thin": {
        "Mg": ("min", -0.067, -0.284, 0.06), "Si": ("min", -0.047, -0.298, -0.07), "Ca": ("min", -0.041, -0.163, 0.42),
        "O": ("min", -0.028, -0.414, 0.07), "C": ("linear", -0.030, -0.163), "Al": ("min", 0.032, -0.029, 0.06),
        "Eu": ("linear", 0.060, -0.404)},
    "high-alpha": {
        "Mg": ("max", 0.357, -0.373, -0.73), "Si": ("max", 0.253, -0.199, -0.99), "Ca": ("max", 0.282, -0.251, -0.99),
        "O": ("max", 0.606, -0.495, -0.91), "C": ("max", 0.191, -0.521, -0.10), "Al": ("max", 0.225, -0.298, -0.48),
        "Eu": ("max", 0.446, -0.290, -0.93)},
}
STELLAR_XFE_SIGMA = {
    "thin": {"C": 0.095, "O": 0.102, "Mg": 0.046, "Si": 0.072, "Ca": 0.032, "Eu": 0.077},
    "high-alpha": {"C": 0.108, "O": 0.128, "Mg": 0.051, "Si": 0.126, "Ca": 0.035, "Eu": 0.084}}
"""dict: Residual scatter (dex) of the fits above, by branch. C, O, Mg, Si and Ca are refitted to Brewer and Fischer 2016 (1,615 dwarfs,
`scripts/fit_stellar_xfe_refit.py`, GEN.988 and GEN.990): the scale of each that puts the generated variance of every populous [Fe/H] bin
on the table's robust one, the means and the branch mixture included; the ratio between the branches is the old one, and calcium's is
magnesium's. Eu is Battistini and Bensby 2016's."""
STELLAR_XFE_CORRELATION = (
    (1.00, 0.47, 0.57, 0.64, 0.35),
    (0.47, 1.00, 0.26, 0.42, 0.41),
    (0.57, 0.26, 1.00, 0.68, 0.02),
    (0.64, 0.42, 0.68, 1.00, 0.30),
    (0.35, 0.41, 0.02, 0.30, 1.00))
"""tuple: Residual correlation of C, O, Mg, Si and Ca (Brewer and Fischer 2016, 1,558 dwarfs of the full table with all five, the residual
being [X/Fe] less its median in the 0.1 dex bin of [Fe/H]), so Mg/Si, a ratio of two correlated terms, has its narrow real scatter. The
844-dwarf table of the first fit gave 0.56, 0.69, 0.71, 0.36, 0.37, 0.77 for the C, O, Mg, Si pairs; calcium is new (GEN.988)."""
STELLAR_CEMP = {"fe_h_below": -2.0, "probability": 0.25, "carbon_fe_min": 0.7, "carbon_fe_span": 1.3, "rich_probability_below": -4.0}
"""dict: Carbon-enhanced metal-poor stars (Frebel and Norris 2015): below [Fe/H] -2 a fraction `probability` have [C/Fe] of
`carbon_fe_min` plus a uniform share of `carbon_fe_span`. The literature's rise to 1 below -4 is unreachable under the stored
[Fe/H] floor of -3.5 (`METALLICITY_RANGE_FE_H`)."""
STELLAR_NEON_FOLLOWS = "O"
"""str: Which element's [X/Fe] neon takes: "O" (the GEN.464 filing; oxygen from the NLTE fit) or "Mg" (part A 3.2); the two
differ in Z/Z0 by 8 percent at most (stellar-composition.md D11)."""
STELLAR_TH_EU = {"log_ratio_mean": 0.014, "log_ratio_sigma": 0.045}
"""dict: [Th/Eu] at zero age (Botelho et al. 2019, arXiv:1810.10413 Table 4); Th is measured near solar [Fe/H] only, so Th
elsewhere is inferred through Eu."""
STELLAR_SOLAR_SYSTEM_AGE_GY = 4.568
"""float: How far back the solar scale's Th and U are decayed to give their birth abundances (the age of the oldest solids;
the photospheric Th and U of A21 are meteoritic)."""
# Surface changes of the first dredge-up on the red-giant branch (Boothroyd and Sackmann 1999, arXiv:astro-ph/9512121 Tables 2 and
# 4): per initial mass, the rise in helium mass fraction and the C/Z and N/Z at the end, for Z = 0.02 (initial Y 0.28, C/Z
# 0.1733, N/Z 0.0531) and Z = 0.001 (initial Y 0.24, C/Z 0.0765; N not extracted at low Z).
STELLAR_DREDGE_UP_MASSES_SOL = (1.0, 1.5, 2.0, 2.5)
STELLAR_DREDGE_UP_Z_HIGH = {
    "z": 0.02, "initial_y": 0.28, "initial_c_over_z": 0.1733, "initial_n_over_z": 0.0531,
    "final_y": (0.3025, 0.2957, 0.2888, 0.2915), "final_c_over_z": (0.1530, 0.1273, 0.1167, 0.1143),
    "final_n_over_z": (0.0771, 0.1073, 0.1197, 0.1324)}
STELLAR_DREDGE_UP_Z_LOW = {
    "z": 0.001, "initial_y": 0.24, "initial_c_over_z": 0.0765,
    "final_y": (0.2598, 0.2606, 0.2545, None), "final_c_over_z": (0.0654, 0.0518, 0.0449, None)}
# Wolf-Rayet surface composition (Crowther 2007, arXiv:astro-ph/0610356 section 3.5): mass fractions.
STELLAR_WOLF_RAYET = {
    "WN-early": {"H": 0.0, "N": 0.01, "C": 0.0005},
    "WN-late": {"H": (0.05, 0.25), "N": 0.01},
    "WN/WC": {"C": 0.05, "N": 0.01},
    "WC": {"C": (0.20, 0.55)},
    "WO": {"C": (0.30, 0.55), "O": (0.10, 0.40)}}
"""dict: Surface mass fractions by Wolf-Rayet subtype; a pair is a range. WC carbon is the refined 20 to 55 percent (C/He 0.1 to
0.5 by number), inside the wider 10 to 60 percent of the literature; the rest of each star is helium (and some oxygen in WN and WC,
the split open). WO is a stand-in (registry row BC24): Crowther 2007 names the oxygen-rich stage without mass fractions, so the carbon and
oxygen ranges are argued from helium burning (carbon first made, then turned to oxygen by alpha capture; carbon still the larger,
oxygen well above WC's, the rest helium) and the draw takes both at one position, so the two never sum past 0.95 [I]."""
WR_WN_EARLY_MAX_SUBCLASS = 5
"""int: A WN star of this subclass number or lower is an early one (WNE, WN2 to WN5, hydrogen-free); a higher number is a late one
(WNL, WN6 to WN9, with the hydrogen the late row gives). The usual convention [R]."""
# White dwarfs: core type by mass (Althaus et al. 2010, arXiv:1007.2659), atmosphere type shares (Xu, Rogers and Blouin 2024,
# arXiv:2404.15425).
STELLAR_WD_HELIUM_CORE_MAX_SOL = 0.45
"""float: A white dwarf below this mass has a helium core. Recalled, open (no source read for the number)."""
STELLAR_WD_ONE_CORE_MIN_SOL = 1.05
"""float: A white dwarf above this mass has an oxygen-neon core (Althaus et al. 2010)."""
STELLAR_WD_DA_SHARE = 0.80
STELLAR_WD_POLLUTED_SHARE = (0.25, 0.50)
# Neutron stars (Ozel and Freire 2016, arXiv:1603.02698; Romani 2022; Rezzolla et al. 2018; Miller et al. 2021).
STELLAR_NS_MASS_GAUSSIANS = (("slow", 1.49, 0.19), ("recycled", 1.54, 0.23))
"""tuple: `(population, mean, sigma)` of the observed neutron-star masses in solar masses."""
STELLAR_NS_MASS_FLOOR_SOL = 1.17
STELLAR_NS_MASS_HEAVIEST_SOL = 2.35
STELLAR_NS_MASS_HEAVIEST_SIGMA_SOL = 0.17
STELLAR_NS_MAXIMUM_NONROTATING_MASS_SOL = (2.01, 2.16)
STELLAR_NS_RADIUS_BAND_KM = (11.6, 13.1)
"""tuple: The hull of Miller et al. 2021's 1-sigma radii at 1.4 and 2.08 solar masses; nearly mass independent, not a cap
(J0740+6620 alone is 13.7 +2.6 -1.5 km). GEN.463 owns the radius function."""
# Brown dwarfs: limits in Jupiter masses (Spiegel, Burrows and Milsom 2011, arXiv:1008.5150 Table 1) and in solar masses
# (Burrows et al. 2001, arXiv:astro-ph/0103383).
STELLAR_BROWN_DWARF = {
    "deuterium_mj": 13.55, "deuterium_solar_z": 0.0189, "deuterium_z_slope_mj_per_dex": -0.9,
    "deuterium_y_slope_mj_per_0p01": -0.2, "deuterium_reference_y": 0.25, "deuterium_range_mj": (11.0, 16.3),
    "lithium_sol": 0.06, "hydrogen_solar_z_sol": (0.070, 0.074), "hydrogen_zero_z_sol": 0.092}

# Globular clusters (GEN.159, docs/design/galaxy/globular-clusters.md). A cluster is stored as a density component, not as
# stars: its profile is a King (1966) model of one of ten concentrations `KING_W0_RANGE` (inclusive), solved once
# (`physics/king.py`) on `KING_TABLE_POINTS` log-spaced radii and scaled by length. A cluster's systems number its
# mass over the mean mass of its systems, worked out from its age (`generation/system_mass.py`, GEN.598). The tidal
# radius a cluster may not exceed is the Jacobi radius in an isothermal halo of circular speed
# `CLUSTER_HALO_CIRCULAR_KMS`.
KING_W0_RANGE = (3, 12)
KING_TABLE_POINTS = 2000
# Satellite galaxies (GEN.237, docs/design/galaxy/intergalactic-and-halo-contents.md 5.5) are stored the same way: the Clouds
# and Sagittarius as an exponential profile (density falling as exp(-r / a)), the dwarf spheroidals as a Plummer sphere
# (density as (1 + (r / a)^2)^-2.5), each cut off at a tidal radius that is `SATELLITE_TIDAL_OVER_HALF_MASS` times
# its half-mass radius (an exponential at 16 scale lengths holds all but 1e-5 of its mass, a Plummer sphere at
# 7.8 scale lengths 98%; the density is normalised so the cut profile holds the whole mass). The rows are real
# objects, but their sizes are recalled from McConnachie (2012) and not checked (the sources were out of reach).
SATELLITE_TIDAL_OVER_HALF_MASS = 6.0
CLUSTER_HALO_CIRCULAR_KMS = 220.0

# A generated galaxy's globular cluster system (`generation/clusters.py`, docs/design/galaxy/globular-clusters.md section 8):
# the count from the halo mass (`CLUSTER_SYSTEM_FRACTION_RANGE`); masses lognormal about `CLUSTER_MASS_PEAK_SOL` with `CLUSTER_MASS_SIGMA_DEX`, truncated to
# `CLUSTER_MASS_RANGE_SOL`; half-mass radii lognormal about `CLUSTER_HALF_MASS_RADIUS_PC` (width in dex, growing with
# distance from the centre as `(R / reference) ** CLUSTER_RADIUS_GROWTH` (the reference radius in the galaxy's own
# scale lengths, `CLUSTER_RADIUS_REFERENCE_SCALE_LENGTHS`), limited to `CLUSTER_HALF_MASS_RADIUS_RANGE_PC`);
# a fraction `CLUSTER_CORE_COLLAPSED_FRACTION` is core-collapsed (W0 at the top of the range) and the rest draw their
# concentration c from a normal about `CLUSTER_CONCENTRATION_MEAN` (all recalled from Harris 1996 and Baumgardt and
# Hilker 2018, good to the shape of the distribution). [Fe/H] has two groups, `(fraction, mean, sigma)`: the
# metal-poor halo clusters and the metal-rich ones of the bulge and disk (the two peaks near -1.5 and -0.5), the
# metal-rich group sitting on flattened orbits about the centre (median `CLUSTER_BULGE_MEDIAN_PC`, rotating at
# `CLUSTER_HALO_CIRCULAR_KMS`) and the halo group on a density falling as r**-3.5 between `CLUSTER_HALO_RADIUS_RANGE_PC`.
CLUSTER_SYSTEM_FRACTION_RANGE = (2.9e-5, 6.0e-5)
"""float: The globular cluster system's mass over its halo's mass, `eta` (GEN.403): Harris, Blakeslee and Harris 2017
(2.9e-5) to Harris, Harris and Alessi 2013 (6e-5), Hudson 2014 (4e-5) between them. The three disagree by a factor of two,
so each galaxy draws its `eta` log-uniform over the whole range (the standing rule: two in question, the widest range)."""
CLUSTER_SYSTEM_SCATTER_DEX = 0.2
"""float: The scatter of a system's mass about `eta` times the halo, dex (Hudson et al. 2014: at most 0.2 dex intrinsic)."""
CLUSTER_MASS_PEAK_SOL = 2.0e5
CLUSTER_MASS_SIGMA_DEX = 0.5
CLUSTER_MASS_RANGE_SOL = (1.0e4, 3.0e6)
CLUSTER_HALF_MASS_RADIUS_PC = 5.0
CLUSTER_HALF_MASS_RADIUS_SIGMA_DEX = 0.3
CLUSTER_HALF_MASS_RADIUS_RANGE_PC = (1.5, 40.0)
CLUSTER_RADIUS_GROWTH = 0.3
CLUSTER_RADIUS_REFERENCE_SCALE_LENGTHS = 8000.0 / 2600.0
"""float: The radius clusters have their typical size at, in the disk's scale lengths (the Milky Way's 8 kpc over its
2.6 kpc, scaled to the galaxy; GEN.404 took the fixed 8 kpc out)."""
CLUSTER_CORE_COLLAPSED_FRACTION = 0.2
CLUSTER_CONCENTRATION_MEAN = 1.45
CLUSTER_CONCENTRATION_SIGMA = 0.4
CLUSTER_AGE_RANGE_GY = (10.0, 13.0)
CLUSTER_FE_H_GROUPS = ((0.65, -1.55, 0.3), (0.35, -0.55, 0.25))
CLUSTER_BULGE_MEDIAN_PC = 3000.0
CLUSTER_BULGE_FLATTENING = 0.5
CLUSTER_HALO_RADIUS_RANGE_PC = (1000.0, 120000.0)
CLUSTER_HALO_DISPERSION_KMS = 110.0
CLUSTER_BULGE_DISPERSION_KMS = 60.0

# A cluster's density in the sector gate (`galaxy/component_field.py`, GEN.160): the density a cell gets is the cluster's
# systems inside the cell over the systems a cell holds at density 1, so the cube is integrated: on
# `CLUSTER_CELL_SUBSAMPLES` points a side within `CLUSTER_NEAR_EDGES` cell edges of a cluster's centre (where its core
# varies faster than a cell) and at the centre alone beyond. The skeleton's bound grows a cluster by
# `CLUSTER_BOUND_SLACK_EDGES` cell edges so a band that holds part of one still qualifies.
CLUSTER_CELL_SUBSAMPLES = 6
CLUSTER_NEAR_EDGES = 3.0
CLUSTER_BOUND_SLACK_EDGES = 2.0
CLUSTER_SMALL_REACH_PC = 256.0
"""float: A component whose tidal sphere is no wider than this (the young open clusters, GEN.244) is found by a window in a
sorted radius in `ComponentField.bound` and `high_band`, not tested one by one: tens of thousands of them are drawn."""
CLUSTER_GRID_PC = 1024.0

# A generated galaxy's young open clusters (GEN.244, `generation/open_clusters.py`, `physics/open_cluster_evolution.py`,
# docs/design/galaxy/intergalactic-and-halo-contents.md 5.14). Every number is from a paper read for this item, none from the
# Milky Way's catalogue of clusters; the sources are in the design note's table.
# Formed: a share of the disk's star formation (the cluster formation efficiency, Kruijssen 2012, MNRAS 426, 3008,
# arXiv:1208.2963, his equation 45: `Gamma = 1 / (1.15 + 0.6 S^-0.4 + 0.05 S^-1)` for the star formation rate density
# `S` in solar masses per year per kpc2) goes into bound clusters whose initial masses fall as M^-2 between
# `OPEN_CLUSTER_MASS_RANGE_SOL` (Lamers et al. 2005, A&A 441, 117, arXiv:astro-ph/0505558: slope -2, a physical upper
# limit as high as 1e6; Lada and Lada 2003, ARA&A 41, 57, arXiv:astro-ph/0301540: 90% of embedded stars are in groups
# of 100 or more members, masses over 50).
OPEN_CLUSTER_CFE_COEFFICIENTS = (1.15, 0.6, -0.4, 0.05, -1.0)
OPEN_CLUSTER_MASS_RANGE_SOL = (50.0, 1.0e6)
OPEN_CLUSTER_MASS_FUNCTION_SLOPE = 2.0
OPEN_CLUSTER_MAX_COUNT = 100_000
"""int: The most young open clusters a galaxy draws, a program limit (GEN.244). A galaxy shaped to hold its star formation
in a few parsecs would form millions (the cluster formation efficiency rises to 0.7); the draw is cut by raising the low end
of the initial masses until this many are expected, so the clusters left are the heavier ones (the mass function falls as
M^-2, so the lighter ones hold little of the stars). The galaxies the generator draws form 5,000 to 92,000."""
OPEN_CLUSTER_AGE_RANGE_MYR = (OPEN_CLUSTER_AGE_EDGES_MYR[0], OPEN_CLUSTER_AGE_EDGES_MYR[-1])
# Lost: the cluster loses mass to stellar evolution (`OPEN_CLUSTER_EVOLUTION_LOSS`, Lamers et al. 2005 Table 1, rows of
# `(Z, a, b, c)`: the share of the initial mass lost is `10 ** ((log10 t - a) ** b + c)` for `t` in years, above
# `OPEN_CLUSTER_EVOLUTION_LOSS_START_MYR`) and to tidal disruption on a timescale that falls as the cluster loses mass,
# `t_dis = t0 (M / Msun) ** gamma` with `gamma = 0.62`, `t0 = Cenv 1e-4 gamma (rho_amb / (Msun / pc3)) ** -0.5` (Lamers et
# al. 2005 and Lamers, Gieles and Portegies Zwart 2005, A&A 429, 173, arXiv:astro-ph/0408235, where `Cenv` is 300 to
# 800 Myr and `rho_amb` the mean density where the cluster lies). `Cenv` is drawn once per galaxy, from the seed. A
# cluster whose mass falls below `OPEN_CLUSTER_MASS_RANGE_SOL`'s low end is gone. The solar metallicity of the table
# is Z = 0.02 (its own convention); a cluster's Z is that times 10 ** [Fe/H].
OPEN_CLUSTER_DISRUPTION_GAMMA = 0.62
OPEN_CLUSTER_DISRUPTION_CENV_MYR = (300.0, 800.0)
OPEN_CLUSTER_EVOLUTION_LOSS_START_MYR = 12.5
OPEN_CLUSTER_EVOLUTION_LOSS = (
    (0.0004, 7.06, 0.265, -1.790), (0.0040, 7.06, 0.260, -1.800), (0.0080, 7.03, 0.260, -1.800),
    (0.0200, 7.00, 0.255, -1.805), (0.0500, 7.00, 0.250, -1.820),
)
OPEN_CLUSTER_SOLAR_Z = 0.02
# Sized: a young cluster's half-mass density is 10^(3 +- 1) solar masses per pc3 at ages under 10 Myr and its radius grows
# by a factor of 5 to 10 in its first 100 Myr, then stays (Portegies Zwart, McMillan and Gieles 2010, ARA&A 48, 431,
# arXiv:1002.1961, section 2.3 and figures 8 and 9). Its half-mass radius is `(3 M / (8 pi rho)) ** (1/3)` at birth;
# a tidally limited one keeps a sixth of its Jacobi radius at most [Henon's tidally filling cluster,
# recalled, to verify]: the profile is cut at six half-mass radii (`SATELLITE_TIDAL_OVER_HALF_MASS`), so the cap is
# a sixth of the Jacobi radius.
OPEN_CLUSTER_HALF_MASS_DENSITY_LOG10 = {"mean": 3.0, "sigma": 1.0}
OPEN_CLUSTER_EXPANSION_FACTOR_RANGE = (5.0, 10.0)
OPEN_CLUSTER_EXPANSION_AGE_GY = 0.1
# Placed: where stars of the cluster's age lie (`galaxy.density.population_densities`, the `young` population under
# 0.1 Gy and the `intermediate` one over it), on the circular orbit of the disk with that population's dispersion.
# A cluster's stars share one [Fe/H] to within `COMPONENT_FE_H_SIGMA["open_cluster"]`, which stands in for the
# homogeneity of open clusters [stand-in, to verify]; its value is the disk's mean at its radius, drawn with the
# population's scatter.
OPEN_CLUSTER_DISK_POPULATION_EDGE_GY = 0.1

# The galaxy density model's fixed Milky Way structure
# (`galaxy.density`, GEN.118, GEN.119): everything that isn't a
# `GalaxyShape` field, as ratios to those fields so it scales with the
# galaxy. Sources: Bland-Hawthorn & Gerhard 2016 (ARA&A 54:529, "BHG16"),
# Dwek et al. 1995 (ApJ 445:716), Wegg & Gerhard 2013 (MNRAS 435:1874).
#
# The reference radius, over the thin disk's scale length: the radius the populations' motions are characterised at,
# the nebula and star-formation profiles are quoted at, and the bar and the satellites' frame are angled from. A
# stand-in ratio from the Milky Way's R0 = 8.2 kpc and L = 2.6 kpc (BHG16), scaled to whatever galaxy this is; it
# fixes no position in the galaxy (GEN.404 removed the home position and the density's calibration point).
REFERENCE_RADIUS_SCALE_LENGTHS = 8200.0 / 2600.0
# The thick disk's mass, length and height, the bulge's axes and angle, and the stellar halo's share are drawn
# from the galaxy seed (`galaxy.structure_draw`, GEN.532, GEN.533); the Milky Way's published ones (BHG16: a thick
# disk of 900 pc and 2.0 kpc against the thin disk's 300 pc and 2.6 kpc; Dwek et al. 1995's G2 bulge fit with
# y0/x0 = 0.392 and z0/x0 = 0.272) are `galaxy.structure.REFERENCE_STRUCTURE`, a check target. The bulge is the bar's
# boxy/peanut core, as COBE/DIRBE sees it edge-on: rho = exp(-r_s^2 / 2) with r_s^4 =
# ((x/x0)^2 + (y/y0)^2)^2 + (z/z0)^4; `bulge_scale_radius_pc` is x0 (along the bar).

# The nuclear star cluster and the nuclear stellar disc (GEN.280, docs/design/galaxy/realistic-galactic-core.md 3), the
# two structures at the real centre, are density terms of the model. Both are ratios to the galaxy's own bulge, so
# a toy galaxy gets a proportionate centre, and `GalaxyShape.nuclear_amplitude` (1 is the observed centre, 0 is
# none) scales the pair together. Calibrated on the default Milky Way-scale model, whose bulge holds 3.9e10 solar
# masses at 0.4 solar masses a system (the real bulge's 1.5e10 is in the model's units 2.6 times that, as every
# component's is; the totals are tuned to star counts, not to mass). The cluster's mass is drawn from the galaxy's
# stellar mass (`galaxy.structure_draw`, GEN.402: Neumayer 2020 and Georgiev 2016, or none); the disc's from the stellar mass
# too (GEN.628: Gadotti and de Sa-Freitas 2025), 65% of an exponential disc of scale length 90 pc lying inside 200 pc.
# The size ratios below are STAND-INS (gap plan 6) from the Milky Way's measured centre, kept until a size law is read.
NUCLEAR_CLUSTER_SCALE_RATIO = 2.07 / 1580.0
"""float: The cluster's Hernquist scale radius over the bulge scale radius (2.07 pc against 1,580 pc: a half-mass
radius of 5 pc, a projected half-light radius of 3.8 to 4.2 pc, against the observed 4.2)."""
NUCLEAR_DISC_LENGTH_RATIO = 90.0 / 1580.0
NUCLEAR_DISC_HEIGHT_RATIO = 35.0 / 1580.0
"""float: The nuclear disc's exponential scale length and `sech^2` scale height over the bulge scale radius (90 and
35 pc; the observed height is 30 to 45 pc, recalled)."""
NUCLEAR_RESOLUTION_PC = 2.0
"""float: The radius of the grid's core ball (shell 0, half the 4 pc shell depth). The cluster's cusp is held below
the radius at which its density equals the ball's mean, so the ball's cell holds exactly the cluster's mass inside
2 pc and no cell is asked for a density the grid cannot hold (`galaxy.density.model_terms`)."""

# The galactic centre's gas and young stars (GEN.282, `galaxy.core_gas`, docs/design/galaxy/realistic-galactic-core.md 3.4).
# A bar carries the gas of its gap inward (`BAR_GAP_RADIUS_RATIO`) and piles it into a ring at the centre, the central
# molecular zone; the ring forms stars at a small share of what comes in. All three are from the review of Henshaw et al.
# 2023 (PPVII, arXiv:2203.11223, read): the bar's gap gas leaves on a Gy (the simulations behind `BAR_GAP_STAR_FORMATION_FLOOR`),
# the inflow to the ring is 0.8 solar masses a year for the Milky Way (Hatchfield et al. 2021, 30% efficiency of the
# dust lanes' flow) against a star formation rate of 0.07 (Table 1: source counts 0.07, integrated light 0.07 to 0.09, more
# or less constant over the last 5 Myr to 100 Myr), by continuity `inflow = star formation + outflow + growth` (their eq. 1).
CORE_GAS_INFLOW_TIME_GY = 1.0
"""float: The time the bar takes to carry the gas of its gap to the centre, Gy."""
CORE_STAR_FORMATION_PER_INFLOW = 0.07 / 0.8
"""float: The ring's star formation rate over the gas inflow to it (the Milky Way's measured 0.07 over 0.8 solar masses a
year; one galaxy measured, so a stand-in for the ratio in others) [R, stand-in for other galaxies]."""
CORE_RING_RADIUS_BULGE_RATIO = 100.0 / 1580.0
CORE_RING_WIDTH_BULGE_RATIO = 50.0 / 1580.0
CORE_RING_HEIGHT_BULGE_RATIO = 30.0 / 1580.0
"""float: The ring's radius, full width at half maximum and vertical scale height over the bulge's scale radius, as the
nuclear disc's are (`NUCLEAR_DISC_LENGTH_RATIO`): 100 pc, 50 pc and 30 pc for the Milky Way's 1.58 kpc [recalled: the
review gives the ring at 100 to 200 pc, so these are placeholders tied to the nuclear disc, to verify]."""
CORE_STAR_FORMATION_SPAN_GY = 1.0
"""float: How long the ring is taken to have formed stars at its present rate (open clusters of the ring are drawn over
this span; the review has the rate constant for at least 5 Myr and 0.09 over 100 Myr) [placeholder]."""
CORE_RING_GAS_RESIDENCE_MYR = 50.0
"""float: How long the gas the bar brings stays in the central ring, Myr: the ring's gas is the inflow times this, over the ring's
area its gas surface density (`galaxy.core_gas`, GEN.1195). The Milky Way's ring holds about 4e7 solar masses of gas (recalled,
2e7 to 6e7) against the inflow of 0.8 a year [stand-in, one galaxy; the mass of gas of the ring is the sourced quantity]."""

# The most massive cluster a ring can form and the clouds that disrupt it (GEN.1195, `physics/cluster_maximum_mass.py`,
# `physics/open_cluster_evolution.cloud_shock_mass_lost_sol`). The maximum is Reina-Campos and Kruijssen 2017 (MNRAS 469,
# 1282, arXiv:1704.00732, read Eqs 1 to 11, section 5 and footnote 9): `M_cl,max = eps Gamma M_T f_coll`, with the Toomre
# mass `M_T = 4 pi^5 G^2 Sigma_g^3 / kappa^4` and the collapsed fraction `min(1, t_fb / t_ff,2D) ** 4`.
CLUSTER_MAX_MASS_EFFICIENCY = 0.1
"""float: The star formation efficiency of a collapsing cloud (`eps`, Eq 10; Lada and Lada 2003 give 0.1 to 0.3 for embedded
clusters, 0.03 to 0.06 for nearby clouds; the paper takes 0.1 and a factor of three either way)."""
CLUSTER_MAX_MASS_CEILING_SOL = 1.0e5
"""float: A hard ceiling on the maximum cluster mass of a ring: the most massive cluster seen in the Milky Way's central ring
implies a maximum of a few 1e4 and "certainly under 1e5" (Reina-Campos and Kruijssen 2017, footnote 9, by the size-of-sample
argument of Gieles 2009) [guard: the relation goes as Sigma_g cubed, so a factor of two in the gas is a factor of eight]."""
CORE_RING_GAS_DISPERSION_KMS = 5.0
DISC_GAS_DISPERSION_KMS = 10.0
"""float: The one-dimensional velocity dispersion of the gas of the central ring (5 km/s, Henshaw et al. 2016) and of a disc
(10 km/s, Heiles and Troland 2003), the values Reina-Campos and Kruijssen 2017 take for the Milky Way (Section 4.1); it sets
the Toomre parameter `Q = kappa sigma / (pi G Sigma_g)` (Eq 5) [stand-in: one galaxy's values for every galaxy]."""
CLUSTER_FEEDBACK_FIRST_SUPERNOVA_MYR = 3.0
CLUSTER_FEEDBACK_ENERGY_RATE_CM2_S3 = 0.16
CLUSTER_STAR_FORMATION_PER_FREE_FALL = 0.012
CLUSTER_MIDPLANE_GRAVITY_FACTOR = 3.0
"""The feedback time of Eq 4 to 7: the time of the first supernova (3 Myr), the energy injection rate per unit stellar mass
(0.16 cm2 s-3), the star formation per free-fall time of the gas (0.012; Elmegreen 2002, Krumholz et al. 2012) and the factor
for the gravity of the stars in the mid-plane density (`phi_P` = 3)."""
CLUSTER_SHOCK_GIELES_COEFFICIENT_GYR = 2.0
CLUSTER_SHOCK_GIELES_NUMERATOR = 5.1
CLUSTER_SHOCK_REFERENCE_MASS_SOL = 1.0e4
CLUSTER_SHOCK_REFERENCE_RADIUS_PC = 3.75
"""Disruption by passing molecular clouds (Gieles et al. 2006, MNRAS 371, 793, arXiv:astro-ph/0606451, Eq 24):
`t_dis = 2.0 S (M / 1e4) (3.75 pc / r_h) ** 3` Gyr with `S = 5.1 / (Sigma_n rho_n)`, `Sigma_n` the surface density of the
clouds and `rho_n` their density (5.1 is the solar neighbourhood's 170 times 0.03 as the paper writes it)."""
CORE_RING_CLOUD_SURFACE_DENSITY_MSUN_PC2 = 5.0e3
CORE_RING_CLOUD_DENSITY_MSUN_PC3 = 150.0
"""float: The clouds of the central ring that disrupt its clusters: the Brick, the template cloud, at 5e3 solar masses per pc2,
and the mean molecular density of the ring at 150 solar masses per pc3 (Kruijssen et al. 2014, MNRAS 440, 3370, arXiv:1303.6286,
sections 3.1.2.6 and 3.1.3: 5.7 Myr for a cluster of 1e4 solar masses and half-mass radius 0.5 pc; they call it an upper
limit of the disruptive potential, and the ring's own density from the generated ring is 17, a factor of nine lower)
[stand-in: one template cloud for every galaxy]."""
# The bar's long axis leads the Sun-center line by an angle drawn from the seed (`galaxy.structure_draw.draw_bar`);
# the Milky Way's is 27 degrees (its near end at positive galactic longitude; Wegg & Gerhard 2013, BHG16: 28-33).

# The stellar halo (GEN.234, docs/design/galaxy/intergalactic-and-halo-contents.md 2.1 and 4) replaces the flat
# floor GEN.78 put under the density (0.1% of the local density everywhere, half the real halo at the Sun
# and 200 to 340,000 times too much at 50 to 262 kpc). It is a spheroid with a broken power law,
# `A * (max(m, core) / break) ** -n` with `m = sqrt(R**2 + (z / q)**2)`, `n` the inner slope inside the break and
# the outer slope beyond it. Radii are ratios to the thin disk's scale length (2.6 kpc in the Milky Way), so a
# galaxy of any size gets a halo that scales with it. Sources: Deason et al. 2019 (1.5e9 solar masses), SDSS
# fits (inner slope 2.3 to 3.3, outer 4.6 to 4.8, break 17 to 27 kpc, flattening 0.6 to 0.8), all in the design
# note.
HALO_SCALE_LENGTH_PC = 2600.0
HALO_BREAK_RATIO = 25000.0 / HALO_SCALE_LENGTH_PC
"""float: The break radius over the disk scale length (25 kpc, 9.6 scale lengths)."""
HALO_CORE_RATIO = 500.0 / HALO_SCALE_LENGTH_PC
"""float: The radius inside which the halo is flat, over the disk scale length (500 pc; the bulge dominates there)."""
HALO_INNER_SLOPE = 2.5
HALO_OUTER_SLOPE = 4.5
HALO_AXIS_RATIO = 0.7
# The halo's amount is the galaxy's own (`galaxy.structure.stellar_halo_system_fraction`, GEN.532): its mass over
# the galaxy's stellar mass drawn from the placeholder law of `galaxy.structure_draw`; the amplitude is solved from it
# like `k_norm` (`galaxy.density.halo_terms`).

# The sector profile (GEN.235, docs/design/galaxy/intergalactic-and-halo-contents.md 5.3): what a place sets for the
# rules that never change by place. A first model, not a fit, in the design note's order. The thin disc's, the
# bulge's and the halo's velocity dispersions are the galaxy's own (`galaxy.structure`, GEN.532); the rest are the
# middle of the note's ranges (a cluster 1 to 15 km/s inside and its orbit outside); `rotation` says what the
# population's net motion follows. [Recalled, to verify.]
# These are one-axis dispersions, the kind the observed ranges are quoted in (a line-of-sight value): the profile states
# them as three-axis totals, `sqrt(3)` times these, like the disc's, bulge's and halo's (GEN.284), and the pair and
# planet limits read the one-axis value back (`generation.place_terms.from_profile`).
STELLAR_POPULATION_DISPERSION_KMS = {
    "cluster": 8.0,
    "satellite": 25.0,       # the dwarfs' and Clouds' own 20 to 30 km/s (design note 5.3)
    "satellite_young": 25.0, "satellite_intermediate": 25.0,       # the same satellite, so the same motion (GEN.526)
    "nuclear_cluster": 150.0, "nuclear_disc": 80.0,
}
# A young open cluster's stars (GEN.244) take the thin disc's age heating of their own bin, from the galaxy's structure
# (`galaxy.structure.disc_dispersion_kms`), not a value here.
STELLAR_POPULATION_ROTATION = {
    "young": "disk", "intermediate": "disk", "old": "disk", "bulge": "bar", "halo": "none", "cluster": "orbit",
    "satellite": "orbit", "satellite_young": "orbit", "satellite_intermediate": "orbit",
    "nuclear_cluster": "none", "nuclear_disc": "disk",
}
STELLAR_POPULATION_ROTATION.update({name: "disk" for name in OPEN_CLUSTER_POPULATIONS})
# A cluster carries its own [Fe/H] (`physics.metallicity`); until the component table (GEN.237) gives each one
# its own, the profile of a cluster cell shows the two peaks of the Milky Way's clusters (-1.5 and -0.5,
# docs/design/galaxy/globular-clusters.md), as `(share, mean, sigma)`; the shares are recalled (about two thirds are the
# metal-poor ones).
METALLICITY_CLUSTER_PEAKS = ((0.65, -1.5, 0.3), (0.35, -0.5, 0.3))
# Where a cell holds components (GEN.237) the peaks of its `cluster` and `satellite` populations are the components'
# own [Fe/H], each with this scatter (dex), weighted by the systems each puts in the cell; the peaks below stand in
# for a place no component reaches (a sample drawn without a profile). A satellite's stand-in is the Clouds'
# and dwarfs' middle (LMC about -0.4, SMC about -0.9, the dwarf spheroidals -1.0 to -2.1, recalled).
COMPONENT_FE_H_SIGMA = {"cluster": 0.15, "satellite": 0.3, "open_cluster": 0.05}
METALLICITY_SATELLITE_PEAKS = ((1.0, -1.2, 0.5),)

# A satellite's star formation history (GEN.526, `galaxy.satellite_history`, docs/design/galaxy/intergalactic-and-halo-contents.md
# 5.5). Worked out from what its row holds (stellar mass, half-light radius, mean stellar age), with no value from the
# Clouds' own histories (those are the test targets of `tests/test_satellite_history.py`, never inputs).
SATELLITE_HISTORY_SPAN_GY = STELLAR_POPULATION_AGE_RANGES_GY["satellite"][1]
"""float: How long ago a satellite's first stars formed, Gy: the top of the age window, which is the age of the
Universe in this model. Its history runs from there to the present."""
STAR_FORMATION_LAW_COEFFICIENT = 2.5e-4
"""float: The Schmidt-Kennicutt law's coefficient: the star formation rate per area is this times the gas per area (solar
masses per pc2) to `STAR_FORMATION_LAW_EXPONENT`, in solar masses per year per kpc2 (Kennicutt 1998, ApJ 498, 541, eq. 4:
(2.5 +/- 0.7) e-4, read in arXiv astro-ph/9712213)."""
SATELLITE_GAS_TO_STELLAR_SCALE_LENGTH_RATIO = 1.8
"""float: A satellite's gas disc scale length over its stars' (GEN.1192). The HI of a late-type dwarf reaches 1.8 +/- 0.8 times
the optical diameter, the extent defined as 6.4 disc scale lengths (Swaters et al. 2002, A&A 390, 829, arXiv astro-ph/0204525,
abstract read); the extent ratio is read as the scale length ratio [I], and the Clouds' gas (Large Cloud 4.8e8 solar masses,
Staveley-Smith et al. 2003, MNRAS 339, 87, arXiv astro-ph/0210501, abstract read) is the check, within a factor of 2.5."""
SATELLITE_STAR_FORMATION_THRESHOLD_MSUN_PC2 = 3.0
"""float: The gas surface density below which a satellite forms no stars, solar masses per pc2. Schaye (2004, ApJ 609,
667, read in arXiv astro-ph/0205125) puts the empirical threshold at 3 to 10 and the point where the gas turns cold
at about 10 (some H II regions form at 2 to 4); the lower end is taken [I], so a dwarf with a trace of gas still forms
stars."""

# How deep one dense cell is filled (GEN.281, docs/design/galaxy/realistic-galactic-core.md 4.2). A cell that expects more
# systems than the cell star ceiling fills only the stars at or above the fill luminosity cap; the rest is
# unresolved, its count and mass computed from the address and never stored. Both are galaxy settings (`planetgen
# plan --cell-star-ceiling` and `--fill-luminosity-cap`, stored on `galaxy_shape`). Ordinary space never reaches the
# ceiling (a cell of the Milky Way's bar holds about a thousand systems).
CELL_STAR_CEILING = 5_000.0
"""float: The default cell star ceiling, systems a cell expects before it is filled only down to the cap (14
minutes of generation, 2.4 GB, at the measured 16.4 ms and 48 KiB a system)."""
CELL_SECONDS_PER_SYSTEM = 0.02
"""float: What generating one system of a crowded cell takes, seconds (measured with the galactic core's cells, GEN.364:
about 20 ms), for the estimate beside the cell star ceiling on the Generate page."""
CELL_BYTES_PER_SYSTEM = 48 * 1024
"""int: What one generated system takes in the database, bytes (GEN.281's measurement, 48 KiB)."""
FILL_LUMINOSITY_CAP_SOL = 100.0
"""float: The default fill luminosity cap, solar luminosities: the dimmest star a cell above the ceiling resolves.
The lowest a cap can be is the brightest white dwarf (`WD_LUMINOSITY_RANGE_SOL`), the dimmest threshold the star
model takes."""
TIER_DOMINANCE_SHARE = 0.5
"""float: The share of a band's density a term must reach for the band to take that term's tier
(`galaxy.scatter_limits`): the nuclear terms for the core tier, the halo for the halo tier (GEN.333)."""

CORE_TIER_FACTOR = 1.10
"""float: How much higher than the baseline the core's scatter limits default (docs/design/galaxy/realistic-galactic-core.md
5.1, Boss): the core luminosity floor and mass limit follow the galaxy's baseline times this until they are set."""
HALO_MIN_LUMINOSITY_SOL = 500.0
"""float: The default halo luminosity floor, solar luminosities (GEN.283, Boss 1:09 AM EDT 2026-10-16)."""
HALO_MIN_LUMINOSITY_RANGE_SOL = (50.0, 4.0e6)
"""tuple: The lowest and highest halo luminosity floor `plan` accepts (Boss, 10:06 PM EDT 2026-10-15: lowest 50)."""

# Gas and dark matter in a cell (docs/design/galaxy/intergalactic-and-halo-contents.md 2.4 and 2.6). The ambient medium
# is whichever of three is densest: the disc's interstellar medium (about 1 H per cm3 at the plane near the Sun
# [recalled], scaled by the cold-gas fraction), the hot halo (1e-4 per cm3 at 50 kpc falling as r^-1.5 [recalled]),
# and the intergalactic medium (the cosmic mean 1.9e-7 H per cm3, the Local Group 25 times that). Radii scale
# with the galaxy's disk scale length like the stellar halo's. 0.0326 solar masses per pc3 per H per cm3 is
# hydrogen plus helium (the note's 1.9e-7 H per cm3 is 6.2e-9 solar masses per pc3).
DISK_ISM_DENSITY_CM3 = 1.0
HOT_HALO_DENSITY_CM3 = 1.0e-4
HOT_HALO_REFERENCE_PC = 50000.0
HOT_HALO_SLOPE = 1.5
HOT_HALO_CORE_PC = 1000.0
COSMIC_MEAN_DENSITY_CM3 = 1.9e-7
LOCAL_GROUP_GAS_OVERDENSITY = 25.0
GAS_SOLAR_MASS_PER_PC3_PER_CM3 = 0.0326
# The NFW dark-matter halo of docs/design/galaxy/galactic-potential.md (M_h = 5.4e11 solar masses, r_h = 16 kpc): 0.0091
# solar masses per pc3 at the Sun. Both scale with the disk scale length (a flat rotation curve: mass over radius
# is constant).
DARK_MATTER_SCALE_MASS_SOLAR = 5.4e11
DARK_MATTER_SCALE_RADIUS_PC = 16000.0

# Pre-placed bright stars (star_population.sample_bright_stars): the
# mass grid, in log-spaced cells over the IMF's range, on which the chance
# of a star being at least the threshold's brightness is tabulated.
BRIGHT_STAR_MASS_GRID_CELLS = 4000

BRIGHT_TABLE_FE_H_STEP_DEX = 0.25
"""float: The width, dex, of the [Fe/H] bins the bright-star sampling tables are built at (GEN.651): a star's lifetimes
and giant branch follow its birth [Fe/H], so a table's cell bounds are valid for the stars of its own bin only."""

BRIGHT_TABLE_BIN_GRID_CELLS = 500
"""int: The mass cells of a sampling table of an [Fe/H] bin other than the reference composition's
(`BRIGHT_STAR_MASS_GRID_CELLS` has the reference one's, which works out the share the counts use): such a table's
bounds already stretch over the bin's lifetimes, by more than a coarse cell's width, and it takes an eighth of the time
to build, which the worker warm pays for every bin of every population (GEN.651)."""

BRIGHT_TABLE_BOUND_MARGIN = 1.1
"""float: What a table's cell bound is multiplied by so a star anywhere in its [Fe/H] bin is under it (calibrated by
`tests/test_star_population.py`, which draws masses and compositions across every bin and checks the ratio)."""
BRIGHT_TABLE_EVOLVED_MARGIN = 2.0
"""float: What the greatest measure of the windows after the main sequence read at a cell's masses is multiplied by to bound
every mass of the cell (GEN.935): the pieces of a track as spiky as the thermal pulses fall differently for each mass, a
few tens of percent either way. `tests/test_star_population.py` draws masses and compositions across every bin and checks it."""
BRIGHT_TABLE_MAIN_SEQUENCE_SLACK = 1.1
"""float: How much brighter than its threshold a cell's bound reads the main sequence (GEN.652): the MIST rows step by a
few percent over a few tenths of a percent of mass, so the end of a cell's main sequence is read a tenth brighter."""

OBJECT_FIRST_CAPACITIES = (1, 2, 4, 8, 16)
"""tuple: The object counts a sector can take in the object-first scatter (PERF.58): its capacity is the smallest
of these that its expected count rarely exceeds (`OBJECT_FIRST_TAIL`); an object past it is dropped."""

OBJECT_FIRST_TAIL = 1.0e-4
"""float: The chance a sector's Poisson count may exceed its capacity (PERF.58)."""

SCATTER_COVERS_REACH = True
"""bool: The bright-star and mass passes draw the grid's whole reach (262 kpc by default), the bands beyond the outline
included (GEN.322). Always true in a run; the tests turn it off, because a toy galaxy's fixed-size halo would fill the
reach with stars."""

OBJECT_FIRST_FAR_PIECES = 16
"""int: Steps of colatitude each hemisphere of a far run of shells is cut into when the shells too faint to bound one
by one are drawn together (`object_first.ReachGroups.far`, GEN.322)."""

OBJECT_FIRST_FAR_STEPS = 128
"""int: A far run of shells is as thick as 1/this of its distance from the centre (`ReachGroups.far_blocks`)."""

OBJECT_FIRST_FAR_MEAN = 4.0
"""float: A shell beyond the outline is drawn on its own, in groups of bands, while its coarse bound expects at least
this many stars; the shells past the last such one are drawn together (GEN.322). A waste of this many candidates a shell
costs about a hundred microseconds, so it is not set to keep the loss under 10%."""

OBJECT_FIRST_GROUP_BANDS = 64
"""int: Beyond the galaxy's outline the object-first scatter bounds the density over groups of this many bands (GEN.322),
the group's bound being taken at its nearest point; a pass to the 262 kpc reach would otherwise bound 6.75e9 bands."""

POISSON_REJECTION_MEAN = 10.0
"""float: From this mean up `util/poisson.py` draws a Poisson count by transformed rejection (constant time);
below it, by Knuth's method (PERF.60)."""

# Binary mass ratio q = M2/M1, uniform (Moe & Di Stefano 2017, ApJS 230:15,
# find it close to flat); the secondary shares the primary's age.
BINARY_MASS_RATIO_RANGE = (0.1, 1.0)

BROWN_DWARF_COMPANION_CHANCE = 0.03
"""float: The chance a star heavier than the hydrogen-burning limit has a brown dwarf companion besides its stars (GEN.899).
The `star-count` table counts stars, so a companion drawn under the limit used up one of its places (an M2V primary lost
nearly half of them); the brown dwarfs are an extra member in this share of the systems instead. Stand-in tier
"Unverified" (second-hand: Metchev and Hillenbrand 2009 find 3.2% (+3.1, -2.7) of solar-type stars with a companion of
0.013 to 0.075 solar masses at 28 to 1590 AU, the "brown dwarf desert"), tweak scale L0."""

# --- Close-binary types (GEN.270, `generation/multistar.py`, docs/design/multi-star-systems.md 13.2) ---
SPECTROSCOPIC_MAX_AXIS_AU = 0.26
"""float: A pair closer than this is a spectroscopic binary (the catalogue's "close detached binary, a under 0.26 AU",
design note 13.1); `MemberTree.spectroscopic_nodes`."""

CONTACT_TERTIARY_KEEP_SHARE = 0.26
"""float: The share of the contact binaries that hold the placed star and have a tertiary which keep it (GEN.270). The
member-tree draw gives 94% of contact binaries a tertiary (84% before GEN.652 read the stars' radii from MIST), because
the stability screen leaves the inner pair of a hierarchy only short periods; stripping the others from 74% of those where
the pair holds the placed star brings the share to the catalogue's 60% (measured over 120,000 seeded disc systems: 48% with
none kept, 94% with all kept; a stand-in until GEN.631's fill measures it again)."""

YOUNG_GROUP_AGE_MAX_YR = 3.0e5
"""float: A massive system (O or early B, a birth mass of `mass_boundaries.B_SUBCLASS_LOWER_MASS_SOL`'s early-B limit and up) older than this is
never drawn as a young group (the design note's 0.3 Myr, GEN.327). The share that is a group is the `young-group` row of
`SYSTEM_TYPE_RATES`."""

YOUNG_GROUP_DECAY_AGE_RANGE_YR = (1.0e4, 4.0e4)
"""tuple: The range a group's `decay_age_yr` is drawn from, log-uniform: trapezia decay in 1e4 to 4e4 years (arXiv 1809.03537,
1701.03440, 1908.11437; multi-star-systems.md 3.3). A group younger than its decay age is stored as a group, an older one is
drawn as its remnant."""

YOUNG_GROUP_SIZE_WEIGHTS = {3: 0.3, 4: 0.4, 5: 0.3}
"""dict: The weight of each member count of a young group. The Orion Trapezium has 4 to 5 members and the catalogued
trapezium-type multiples run from 3; the weights are a stand-in [P, bracket 3 to 5 members] until a count by size is read."""

YOUNG_GROUP_RADIUS_AU_RANGE = (3.0e3, 3.0e4)
"""tuple: The radius within which a young group's members are drawn, in AU, log-uniform (a stand-in [P]: the Trapezium's
core is some 0.03 pc, 6e3 AU; the bracket is a factor of 3 either side of it)."""

YOUNG_GROUP_VELOCITY_DISPERSION_KMS = 2.5
"""float: The one-dimensional velocity dispersion of a young group's members, km/s (a stand-in [P]: the Orion Trapezium
cluster is some 2 to 3 km/s, recalled, not read in this repository). Each component of a member's velocity is a draw of it."""

YOUNG_GROUP_REMNANT_MAX_AXIS_AU = 300.0
"""float: The widest axis of the hierarchy a decayed group leaves, AU (a stand-in [P]): the lightest member is ejected and the
rest harden, so the remnant is tighter than the field's draw (multi-star-systems.md, `decayed-trapezium-remnant`)."""

MERGER_MASSIVE_MIN_MASS_SOL = 9.0
"""float: The birth mass from which a star counts as massive for the merger-product and mass-gainer shares (the
early-B and O stars, `mass_boundaries.B_SUBCLASS_LOWER_MASS_SOL`; de Mink et al. 2014 give the 8% and 20% for O and early-B stars)."""

YOUNG_POPULATION_MAX_AGE_GY = STELLAR_POPULATION_AGE_RANGES_GY["young"][1]
"""float: The age under which a star belongs to the young population, Gy: the gate of the Be stars, the massive merger
products and the young-only types of the catalogue (the generator knows a star's age, not its population)."""

OLD_POPULATION_MIN_AGE_GY = STELLAR_POPULATION_AGE_RANGES_GY["old"][0]
"""float: The age from which a star belongs to the old population, Gy: the gate of the R Coronae Borealis and FK Comae
stars and the Thorne-Zytkow objects, and from which an Algol no longer forms (the catalogue's `algol` row)."""

BE_STRIPPED_COMPANION_MASS_RANGE_SOL = (0.3, 0.5)
"""tuple: The mass of the stripped helium star that is a Be star's companion (design note 13.2, as the hot subdwarf row:
0.3 to 0.5 solar masses)."""

BE_STRIPPED_PROGENITOR_RATIO_RANGE = (1.2, 2.0)
"""tuple: A stripped companion's birth mass over the Be star's birth mass, uniform. A stand-in: the donor of the transfer
that spun the Be star up was the heavier of the pair, so the ratio is over 1; the range is a design default."""

STRIPPED_HELIUM_RADIUS_SOL_AT_0_47 = 0.2
"""float: The radius of a 0.47 solar mass stripped helium star or hot subdwarf, solar radii (a typical sdB, Heber 2016,
PASP 128:082001); the radius of another mass scales as the mass to the power -1/3 (the degenerate-envelope trend)."""

STRIPPED_PAIR_PERIOD_DAYS = (30.0, 1000.0)
"""tuple: The periods a Be star and its stripped companion are drawn over, log-uniform, days. A stand-in: such pairs
are detached after the transfer and have periods of tens to hundreds of days; the range is a design default, not a
measured one."""

OLD_GIANT_SYSTEM_SHARE = 0.005
"""float: The share of all systems whose placed star is a giant of the old population. A stand-in, used to turn a type's
Galaxy count (R Coronae Borealis, FK Comae, Thorne-Zytkow) into the chance among old giants; remeasure from a fill."""

# --- White-dwarf binary types (GEN.271, `generation/multistar.py`, docs/design/systems/multi-star-systems.md 13.3) ---
DOUBLE_WHITE_DWARF_AXIS_AU = (0.001, 4.0)
"""tuple: The axis of a close double white dwarf, log-uniform, AU (Maoz et al.: within 4 AU); the Peters rule
(`stability.survives_gravitational_waves`) drops a pair that would have merged by the system's age."""

POST_COMMON_ENVELOPE_AXIS_AU = (0.005, 2.0)
"""tuple: The axis of a post-common-envelope white dwarf and companion, log-uniform, AU (design note 13.3)."""

CATACLYSMIC_PERIOD_HOURS = (1.3, 10.0)
"""tuple: The orbital period of a cataclysmic variable, hours (design note 13.3)."""

CATACLYSMIC_MAGNETIC_SHARE = 0.36
"""float: The share of cataclysmic variables whose white dwarf is magnetic (a polar or an intermediate polar)."""

CATACLYSMIC_DONOR_MASS_RANGE_SOL = (0.08, 0.7)
"""tuple: The mass of a cataclysmic variable's main-sequence donor, uniform, solar masses. A stand-in: the donor fills its
lobe, so its mass sets the period, and this range gives the periods of `CATACLYSMIC_PERIOD_HOURS`; a draw whose period
falls outside is made again."""

AM_CVN_PERIOD_MINUTES = (5.0, 65.0)
"""tuple: The orbital period of an AM Canum Venaticorum binary, minutes (design note 13.3)."""

AM_CVN_DONOR_MASS_RANGE_SOL = (0.009, 0.13)
"""tuple: The mass of an AM CVn helium donor, log-uniform, solar masses. A stand-in that gives the periods of
`AM_CVN_PERIOD_MINUTES` with the radius law below; a draw whose period falls outside is made again."""

AM_CVN_DONOR_RADIUS_SOL_AT_0_05 = 0.03
"""float: The radius of a 0.05 solar mass degenerate helium donor, solar radii; other masses scale as the mass to the power
-1/3 (a stand-in for the semi-degenerate donors of AM CVn systems, Nelemans et al. 2001)."""

ACCRETING_FILL_FACTOR = 1.05
"""float: How far a cataclysmic variable's or AM CVn's donor fills its Roche lobe (1 is exactly full); a pair at this
filling has `interaction` accreting."""

HOT_SUBDWARF_POST_CE_SHARE = 0.05
"""float: The share of post-common-envelope pairs whose companion is a stripped helium star or hot subdwarf (kind
`stripped-helium`, 0.3 to 0.5 solar masses). A stand-in: the observed sdB and sdO binaries are a few percent of the
white dwarf plus star pairs, and the note gives no share."""

BARIUM_GIANT_MAX_INITIAL_MASS_SOL = 3.0
"""float: A giant born lighter than this is a G or K giant for the barium and symbiotic shares (the giants of the
intermediate and old populations are nearly all of such mass; a stand-in for the spectral class, which the member does
not carry)."""

BARIUM_PERIOD_DAYS = (80.0, 1.0e4)
"""tuple: The period of a barium or CH star's white dwarf companion, log-uniform, days (Jorissen et al. 2019: tens of
days to tens of years)."""

SYMBIOTIC_PERIOD_DAYS = (200.0, 2000.0)
"""tuple: The period of a symbiotic star's white dwarf, log-uniform, days (Mikolajewska 2012: 200 to 1,500 days)."""

# Accreting compact binaries (GEN.272, docs/design/systems/multi-star-systems.md 13.4). Every figure below is a stand-in
# unless it names a source; the shares of the registry rows (`high-mass-xray-binary`, `low-mass-xray-binary`,
# `microquasar`, `millisecond-pulsar-companion`) are in `generation/phenomena/system_types.py`.
SUPERGIANT_XRAY_PERIOD_DAYS = (1.4, 100.0)
"""tuple: The period of a wind-fed (supergiant) high-mass X-ray binary, log-uniform, days (Cyg X-1 5.6 d, Vela X-1 9 d,
GX 301-2 41 d: 1 to 100 d). A star whose lobe the span would put inside takes the span out to `XRAY_SPAN_ROOM` times the
closest the pair may be."""

BE_XRAY_PERIOD_DAYS = (10.0, 300.0)
"""tuple: The period of a Be/X-ray binary, log-uniform, days (Corbet 1986 and Reig 2011: 10 to 300 d; the orbit is eccentric)."""

LOW_MASS_XRAY_PERIOD_DAYS = (0.17, 1.0)
"""tuple: The period of a low-mass X-ray binary with a main-sequence donor, days (4 to 24 hours, most of the 349 known)."""

LOW_MASS_XRAY_DONOR_MASS_RANGE_SOL = (0.1, 1.5)
"""tuple: The initial mass of the low-mass star that a neutron star or a black hole accretes from, solar masses (a stand-in
for the range of donors of the low-mass X-ray binaries)."""

XRAY_SPAN_ROOM = 4.0
"""float: How far past the closest allowed axis an X-ray binary's span reaches when its period range would put the donor
inside its lobe (see `SUPERGIANT_XRAY_PERIOD_DAYS`)."""

SPIDER_PERIOD_DAYS = (0.1, 1.0)
"""tuple: The orbital period of a black widow or a redback, log-uniform, days (hours to a day, arXiv 2210.16295)."""

BLACK_WIDOW_MASS_RANGE_SOL = (0.01, 0.05)
"""tuple: A black widow's companion mass, solar masses (the registry row: 0.01 to 0.05)."""

REDBACK_MASS_RANGE_SOL = (0.2, 0.4)
"""tuple: A redback's companion mass, solar masses (the registry row: 0.2 to 0.4; 0.07 to 0.1 is empty)."""

MICROQUASAR_SHARE_OF_XRAY_BINARIES = 0.10
"""float: The share of accreting neutron stars and black holes that launch jets (a microquasar). A stand-in: several dozen
microquasars (the registry row's 50) against the 165 high-mass and 349 low-mass X-ray binaries known."""

ACCRETION_EFFICIENCY_BLACK_HOLE = 0.1
"""float: The share of the rest-mass energy of accreted gas that a black hole radiates (a stand-in for a moderate spin; the
neutron star's is its surface potential, `physics/accretion.py`)."""

WIND_MASS_LOSS_MSUN_YR_AT_1E5_LSUN = 1.0e-6
WIND_MASS_LOSS_LUMINOSITY_EXPONENT = 1.5
"""A hot massive star's wind loses `WIND_MASS_LOSS_MSUN_YR_AT_1E5_LSUN` solar masses a year at 1e5 solar luminosities,
scaling as the luminosity to this power (stand-ins within the range of Vink et al. 2001 for O and B supergiants)."""

WIND_SPEED_OVER_ESCAPE = 1.5
"""float: The terminal wind speed of a hot massive star over its surface escape speed (1.5 to 3 for OB stars; the low end)."""

BE_DISC_MASS_LOSS_MSUN_YR = 1.0e-8
BE_DISC_OUTFLOW_KMS = 20.0
"""A Be star's decretion disc supplies the neutron star at periastron: the rate and the slow outflow speed of the disc gas
(stand-ins; Be/X-ray outbursts average 1e36 to 1e37 erg/s)."""

OVERFLOW_TIMESCALE_YR = {"main-sequence": 5.0e8, "giant": 1.0e7}
"""dict: A lobe-filling donor gives up its mass over this many years, by the kind of donor (a stand-in for the angular
momentum loss of a main-sequence donor and the expansion of a giant)."""

# The chance of a companion inside 10 AU against the star's own [Fe/H] (GEN.103; Boss, 2026-10-16, 6:12 to 6:15 PM
# EDT: one pooled trend, every star mass). Moe, Kratter & Badenes 2019, ApJ 875:61 (arXiv 1808.02116), for solar-type
# primaries (0.6 to 1.5 Msun, periods under 10^4 d): 53 +-12% at [Fe/H] -3.0, 40 +-6% at -1.0, 24 +-4% at -0.2 and
# 10 +-3% at +0.5. `generation.close_binary_trend.close_companion_factor` reads the straight line through them, held flat
# past the ends, as `F(fe_h) / F(0)` on the odds of the generator's own close share (Boss, 2026-10-16: the odds, not the chance). The errors are reference values for the check
# (TEST.196), not inputs.
CLOSE_BINARY_FRACTION_BY_FE_H = {-3.0: 0.53, -1.0: 0.40, -0.2: 0.24, 0.5: 0.10}

# The semi-major axis at which half of the pairs of age t survive passing stars is
# `ENCOUNTER_LIMIT_COEFFICIENT * (M1 + M2) * sigma / (G * rho_2 * t)` (Bahcall, Hut & Tremaine 1985, ApJ 290:15, as
# Jiang & Tremaine 2010, MNRAS 401:977 restate it in their equation (1), arXiv 0907.2952), with sigma the dispersion of
# the relative velocity (root 2 times a population's one-dimensional dispersion) and rho_2 the second moment of the
# perturbers' mass function (GEN.447, `generation.place_terms`). The coefficient reproduces their worked 3.1e4 AU.
ENCOUNTER_LIMIT_COEFFICIENT = 0.002

# A system drawn with no galaxy placement (a single system from the command line) has no cell to read a density and
# a dispersion from; it takes the Sun's neighbourhood (a stand-in until the place's own density is read, like the
# distance from the centre `galaxy.scales.default_galactocentric_pc` gives it). Recalled round values (a mass density of 0.1 Msun/pc^3 and a field dispersion of
# 30 km/s; Jiang & Tremaine 2010 use 0.045 Msun/pc^3 for the observed stars and 29 km/s, their equation 47 and section 5),
# to verify.
SOLAR_NEIGHBOURHOOD_MASS_DENSITY_SOLAR_PC3 = 0.1
SOLAR_NEIGHBOURHOOD_DISPERSION_KMS = 30.0

# The orbit beyond which a planet has been stripped by passing systems (GEN.325, `physics.planet_survival`): the largest
# axis whose expected number of ionising passes in `PLANET_SURVIVAL_AGE_GY` is one, for a host of
# `PLANET_SURVIVAL_HOST_MASS_SOLAR` (docs/design/galaxy/globular-clusters.md section 7 and realistic-galactic-core.md
# section 7: 12 Gy and 0.5 solar masses, a typical host of the old populations). A star younger than the age has been
# exposed to its own age only (`generation.place_terms.planet_limit_au`).
PLANET_SURVIVAL_AGE_GY = 12.0
PLANET_SURVIVAL_HOST_MASS_SOLAR = 0.5

# --- Planet Classification Data ---

# The words that describe each class, and its kind ('t' for terrestrial, 'g' for gas giant). A class is derived from a finished body's
# conditions (`physics/planet_classes.py`, `physics/data/planet_classes.json`); nothing here says what a body of the class
# looks like for generation: no radius, zone, density, albedo or probability.
PLANET_CLASSES = {
    'A': {
        "description": 'a small, barren, and volcanic world, occasionally scarred by intense, toxic eruptions',
        "composition": 'igneous silica and basalt',
        "type": 't',
    },
    'B': {
        "description": 'a small, molten world, occasionally the stripped core of a former gas giant',
        "composition": 'iron, potassium, silicon, and sulfur',
        "type": 't',
    },
    'C': {
        "description": 'a dead world',
        "composition": 'anthracite, basalt, and hydrocarbons',
        "type": 't',
    },
    # GEN.1091: the class is a rule over the finished body's conditions (`physics/data/planet_classes.json`): surface pressure under 10 Pa and at
    # least 1 percent of the mass water or ice. Size is not in the rule; the intended scale is a dwarf planet or a moon (radius 50 to 500 km), and
    # the class covers belt members and dwarf planets as well as moons (docs/design/systems/belt-origins.md).
    'D': {
        "description": 'a small icy body',
        "composition": 'frozen hydrocarbons and ice',
        "type": 't',
    },
    'E': {
        "description": 'a deep-frozen world',
        "composition": 'silicon, iron, magnesium, and ice',
        "type": 't',
    },
    'F': {
        "description": 'a volcanic world with shallow seas',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'G': {
        "description": 'a rocky, barren world',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'H': {
        "description": 'a desert world with little liquid water',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'I': {
        "description": 'an ice giant with a tilted magnetic field',
        "composition": 'rock, ice, methane, and ammonia',
        "type": 'g',
    },
    'J': {
        "description": 'a gas giant with a turbulent atmosphere and rings',
        "composition": 'hydrogen and helium',
        "type": 'g',
    },
    'K': {
        "description": 'an adaptable world',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'L': {
        "description": 'a marginally habitable world',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'M': {
        "description": 'a temperate terrestrial world',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'N': {
        "description": 'a hot world',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'O': {
        "description": 'an ocean world with most of its surface covered in liquid water',
        "composition": 'silicon, iron, and magnesium',
        "type": 't',
    },
    'P': {
        "description": 'a cold, glaciated world',
        "composition": 'silicon, iron, magnesium, and ice',
        "type": 't',
    },
    'S': {
        "description": 'a barren, rocky Super-Earth, its surface molten when it orbits close to its star',
        "composition": 'silicate rock over a large iron core',
        "type": 't',
    },
    'T': {
        "description": 'a gas dwarf with a thick atmosphere',
        "composition": 'hydrogen, helium, and hydrocarbons',
        "type": 'g',
        "also_type": 't',       # GEN.1142: a terrestrial-kind body under about 10 Earth masses with a hydrogen-dominated air is a gas dwarf too
    },
    'V': {
        "description": 'a Super-Earth with high gravity',
        "composition": 'iron, iridium, tungsten, and nickel',
        "type": 't',
    },
}


# --- Life Data ---
# What lives on a world is a life class (`physics/data/life_classes.json`, GEN.889), drawn from the conditions the world stores.

# Not a world's schedule since GEN.750: a world draws its own five ages (`physics.life_timeline`, `LIFE_STEP_MEANS_GYR`). This table
# names the three paces (the words `life_facts.pace_label` writes), gives the star-level helpers the paces a star's lifespan allows
# (`stellar_object.supported_evolutionary_scales`, `age_needs`; GEN.942) and is what migration 0156 read to fill the ages of an older database.
EVOLUTIONARY_TIMELINES = {
    "normal": { # Changed from "norm" to "normal" for consistency
        "star_lifespan": 10.0, # Billion Years
        "abiogenesis": 0.5, # Billion Years
        "photosynthesis": 1.5, # Billion Years
        "complex_cells": 2.5, # Billion Years
        "multicellularity": 4.0, # Billion Years
        "technological_civilization": 4.5 # Billion Years
    },
    "fast": {
        "star_lifespan": 0.1, # Billion Years (100 Million Years)
        "abiogenesis": 0.005, # Billion Years (5 Million Years)
        "photosynthesis": 0.015, # Billion Years (15 Million Years)
        "complex_cells": 0.030, # Billion Years (30 Million Years)
        "multicellularity": 0.050, # Billion Years (50 Million Years)
        "technological_civilization": 0.080 # Billion Years (80 Million Years)
    },
    "slow": {
        "star_lifespan": 1000.0, # Billion Years (1 Trillion Years)
        "abiogenesis": 3.0, # Billion Years
        "photosynthesis": 8.0, # Billion Years
        "complex_cells": 15.0, # Billion Years
        "multicellularity": 25.0, # Billion Years
        "technological_civilization": 35.0 # Billion Years
    }
}

EVOLUTIONARY_PACE_LABELS = {"fast": "hyper-accelerated", "normal": "standard", "slow": "decelerated"}
"""dict: How each `EVOLUTIONARY_TIMELINES` pace is worded in the life sentence (DB.100: the words are not in the milestone
table, so a page is written from a world's stored pace and needs no table of ages)."""

# --- Planet Generation Specific Constants ---
HABITABLE_PLANET_CLASSES = ['E', 'F', 'G', 'H', 'K', 'L', 'M', 'O', 'P', 'V']
"""
list: A list of planet class codes that are considered habitable.
"""

"""
dict[str, str]: Maps a habitable planet class code to the highest
`EVOLUTIONARY_TIMELINES`-milestone key its own generated description text
is consistent with -- see `evolution.get_evolutionary_state`'s
`planet_class` parameter, the only reader of this mapping. A class not
present here is uncapped (full range, up to and including a technological
civilization).
"""

HABITABILITY_LIFE_SCORE_MIN = 0.005
"""
float: A planet's `phi_bio` (microbial) or `phi_cpx` (complex life) likelihood
below this counts as none (GEN.206, and a life score of 0, GEN.649): the habitability
index allows no life of that kind, so the description must not claim it.
"""

LIFE_MILESTONES = ("abiogenesis", "photosynthesis", "complex_cells", "multicellularity", "technological_civilization")
"""
tuple[str]: The evolutionary milestones (`EVOLUTIONARY_TIMELINES` keys) in ascending order of life
complexity. `evolution.MILESTONE_KEYS` and `physics.life_scores` read it.
"""

LIFE_SCHEDULE_FIELDS = tuple(f"life_age_{key}_gyr" for key in LIFE_MILESTONES)
"""tuple[str]: The body columns holding a world's milestone schedule (DB.100), one age in Gyr for each `LIFE_MILESTONES`
milestone in order: `life_age_abiogenesis_gyr` ... `life_age_technological_civilization_gyr`."""

LIFE_STEP_MEANS_GYR = (0.5, 1.5, 2.5, 3.0, 60.0)
"""
tuple[float]: The mean waiting time of each hard step of `LIFE_MILESTONES`, Gyr (GEN.750; docs/design/habitability/life-timeline-and-oxygen.md,
reports/best-guess-standins-2026-10-24.md section 14): abiogenesis, oxygenic-type photosynthesis, complex cells, multicellular animals,
civilisation. Each step takes an exponential waiting time with this mean once the one before it is done (the hard-step model of Carter 1983 and
Watson 2008; the early means are shorter than the later because the first steps came early on the one world we know, which is a test point and
not a preset). Boss chose these on 10-28 (9:03 PM EDT) from a Monte Carlo of the timeline: complex cells 2.5, animals 3 and a civilization 60
Gyr make animals common and a civilization rare. Stand-in T7, L3.
"""

LIFE_STEP_SCALE = 1.0
"""
float: A factor on every `LIFE_STEP_MEANS_GYR` mean, the one knob that scales how common every stage is (GEN.750). It is 1.0: the means are the
ones Boss chose. Measured by drawing 40,000 timelines (`life_timeline.schedule`) against a star of any age from 3 to 13 Gyr, oxygenic life at the
reference temperature and oxygen 2 Gyr after the first photosynthesis: complex cells are reached in 65 in 100, animals in 41 in 100 and a
civilization in 2 in 100 (oxygen at 0.7 Gyr: 76, 51 and 2.9). For stars of 1 to 10 Gyr (G) and 1 to 13.8 Gyr (K, M) animals are 21, 33 and 23 in
100 and a civilization 0.6, 1.5 and 0.7. Stand-in T9, L3: chosen so that animals are common and a civilization rare; Boss can retune it.
"""

LIFE_START_AGE_GYR = 0.1
"""float: The age of a world, Gyr, at which the waiting for abiogenesis starts: its surface has cooled and its first ocean has settled (stand-in
T6, L0; Earth's oceans were there by about 4.4 Ga, 0.15 Gyr after it formed)."""

LIFE_ACTIVATION_ENERGY_EV = 0.30
"""float: The activation energy, eV, of the rate of every hard step: a step goes faster on a warm world by a Boltzmann factor. Metabolic theory
gives 0.65 eV for the rate of a cell's chemistry (Gillooly et al. 2001; Allen and Gillooly 2006 for the rate of evolution itself; from memory);
a hard step is a rare event of mutation and selection, whose temperature dependence the geological record cannot pin, and 0.65 eV makes a world at
26 C three times faster than the reference, so that nearly half of the warm worlds of a Sun-like star would have a civilization by 8 Gyr. Half of it
gives 1.6. Stand-in T7, L3: lowered to keep the stages rare on a warm world, inside the range of 0.2 to 0.8 eV of the same papers."""

LIFE_REFERENCE_TEMPERATURE_K = 288.0
"""float: The surface temperature, K, at which a step takes its mean time (the temperature of the world that gave the first test point)."""

LIFE_RATE_FACTOR_RANGE = (0.25, 2.5)
"""tuple[float]: The least and most the temperature changes the rate of a step (stand-in T6, L1: the factor of a few over the liquid range that
the Boltzmann law gives at `LIFE_ACTIVATION_ENERGY_EV`, kept from running away at the edges of it)."""

LIFE_PACE_FAST_ABOVE = 1.5
LIFE_PACE_SLOW_BELOW = 0.67
"""float: A world whose rate factor is above `LIFE_PACE_FAST_ABOVE` has the pace "fast", below `LIFE_PACE_SLOW_BELOW` "slow", else "normal"
(the words the life sentence uses; the pace is a label of how warm the world is, not a table)."""

LIFE_FORCED_TRIES = 4000
"""int: How many timelines a forced technological civilization (`INTELLIGENT_LIFE`) tries before it falls back to squeezing the last one into the
world's age (`life_timeline.forced_waits`)."""

LIFE_SCORE_MAX = 9
"""int: The top of the four life scores (Microbial, Flora, Fauna, Sentient): 0 is none, 9 is abundant (GEN.649)."""

LIFE_SCORE_SIMPLE_FLORA_MAX = 3
"""int: The most flora a world can show before multicellular life: mats and algae (stage photosynthesis or complex
cells) never rise above this, however kind the surface (DESIGN, GEN.649)."""

LIFE_SCORE_SENTIENT_FULL_GYR = 0.01
"""float: Gyr a technological civilization needs after its milestone to fill its surface (Sentient 9 at the best
conditions); a younger one scores in proportion, never less than 1 (DESIGN, GEN.649: 10 Myr)."""

LIFE_NPP = {
    "chemosynthesis_share": 1.0e-3,
    "light_ratio_max": 3.0,
    "tie_factor": 2.0,
}
"""
dict: The net primary production of a life class, the number that picks the class that lives on a world (GEN.889,
docs/design/habitability/life-classes.md 8.2). A light class's is the usable light of the surface over Earth's, times its pigment band's share of
the star's photons over the Sun's share of 400 to 700 nm (at most `light_ratio_max`); a chemical class's is `chemosynthesis_share` times its heat
flow in Earth's. Classes whose production is within `tie_factor` of the largest are tied and the draw between them is weighted by production.
`chemosynthesis_share`: chemosynthesis is about a thousandth of the primary production of Earth's ocean (vent and seep estimates of 0.02 to 0.2
Gt C a year against 50; recalled), stand-in T4, L3. `light_ratio_max` keeps a hot star's photon surplus from running away, T9, L3. `tie_factor` is
the design's "within a factor of 2" (life-classes.md 8.2), T9, L3.
"""

LIFE_COMPLEXITY = {
    "thresholds": {"complex_cells": 1.0e-5, "multicellularity": 3.0e-5, "A1": 1.0e-4, "A2": 1.0e-3, "A3": 1.0e-1},
    "oxygen_floor_kpa": {"A1": 0.1, "A2": 2.0, "A3": 6.0},
    "oxygen_full_yield_kpa": 2.0,
    "anaerobic_yield": 0.11,
    "full_complexity": 0.3,
    "fauna_cap": (0, 3, 6, 9),
    "anoxic_animal_chance": 0.03,
    "anoxic_animal_donor_min": 0.3,
}
"""
dict: The complexity index and the animal grades (GEN.890, GEN.877, GEN.882 to GEN.884; life-classes.md 7.3). A class's index C is its net primary
production over Earth's (`LIFE_NPP`, `life_classes.ClassScore.npp`) times the yield of what it breathes (1 for aerobic respiration, `anaerobic_yield`,
a ninth, for none; an oxygen breather's follows the air, rising from the anaerobic yield to 1 at `oxygen_full_yield_kpa`) times the temperature rate
factor of its clock (`life_timeline.rate_factor`). `thresholds` are the C a class needs for complex cells, multicellular bodies and the three animal
grades A1 (sessile), A2 (small mobile) and A3 (large active; a civilization needs it). `oxygen_floor_kpa` is the oxygen partial pressure, kPa, an
oxygen-breathing class needs for each grade: about 0.1 to 0.8 kPa for sponge-grade animals (Mills et al. 2014; recalled), rising with body size by the
diffusion argument (Catling et al. 2005; recalled). `full_complexity` is the C at which the complex-life score is full. `fauna_cap` is the most
Fauna each grade shows (0 to 9). `anoxic_animal_chance` is the share of anoxic worlds with a donor supply of `anoxic_animal_donor_min` or more that
get animals that have lost aerobic respiration (route R4; contested: Danovaro 2010, Yahalomi 2020). Stand-ins T9, L3 (the floors T2 for A1, T7 above;
the R4 chance T4); set from the first census, not fixed.
"""

FLAVOR_CHANCE_SYSTEM = 0.05 # The chance flavor text will be added to a system.
FLAVOR_CHANCE_PLANET = 0.05 # The chance flavor text will be added to a planet.
MAX_FLAVOR_TOTAL = 3 # The maximum flavor items total in a system.
MAX_RECENT_FLAVOR_TEXTS = 5 # Maximum number of recent flavor texts to keep track of.

SYSTEM_FLAVOR = [ # System-only flavor text.
    "a derelict chain of deep-space platforms firing concentrated energy beams into the void for an unknown purpose.",
    "an automated buoy equipped with a mass-driver that catches physical data pods and flings them toward inner-system hubs.",
    "a massive, ultra-thin sheet of reflective foil drifting through the system, originally deployed to shade a now-dead world.",
    "thousands of small, passive reflectors scattered across the system that reflect electromagnetic waves back and forth.",
    "an ancient, nuclear-powered navigational beacon at the system's edge, driven and maintained by a local artificial intelligence.",
    "a drifting, chaotic cloud of uncoordinated mirror tiles arrayed like a fragmented automated radiation collection swarm.",
    "a system-wide network of navigational satellites that constantly shift and optimize their orbits using a decentralized swarm intelligence.",
    "an invisible mass sitting in an empty orbit that clearly registers on navigation sensors, warping the local space-time.",
    "a spherical dead-zone in deep space where all electromagnetic scans, communication signals, and radiation instantly drop to absolute zero.",
    "a jagged, glowing fracture in space-time that continuously vents electromagnetic energy into the surrounding void.",
    "a decentralized network of automated defensive satellites arrayed across the system.",
    "an array of laser transceivers locked in continuous, high-bandwidth optical communication with one another.",
    "a swarm of automated solar collectors skimming the star's corona, routing steady power into massive, fully charged battery arrays.",
    "a partially built ring of solar collectors encircling the star, its gaps suggesting the project was abandoned mid-construction.",
    "a cluster of dormant self-replicating probes holding formation, their fabrication cycles halted mid-sequence.",
    "a derelict generation ship on a slow, decaying orbit, its hull dark and its drive cold after generations adrift.",
    "a stellar mirror array angled to nudge the star's own light pressure into a gentle directional thrust, now drifting loose.",
    "a pair of gravitational-lensing relays aimed at a target far beyond the system, still patiently listening.",
    "a broken ring segment holding a stable orbit on its own, as if waiting for the rest of the structure to return.",
    "a scattering of paired relay buoys still exchanging signals with a partner that fell silent long ago.",
    "an automated probe holding station at the system's edge, its greeting message looping to no one.",
    "a magnetic siphon rig drawing material from the star's outer layers, still running long after its operators left.",
    "a shell of solar-collecting satellites in a wide, half-finished orbit, their construction drones nowhere to be found."
]

PLANET_FLAVOR = [ # Flavor text for any planet.
    "a massive hollow metallic cylinder embedded in the crust.",
    "the crashed hull of a derelict pre-warp ship.",
    "the crashed hull of a derelict primitive warp ship.",
    "a solar-powered beacon looping a brief audio distress call.",
    "an array of automated radio dishes tracking empty coordinates in deep space.",
    "a vast, precise grid of ancient cargo containers showing faint energy readings.",
    "an abandoned solar farm attached to overfilled battery capacitors that periodically discharge themselves in violent plasma arcs.",
    "a vast, perfectly smooth glass crater with refined materials fused into its center.",
    "an automated surface rover actively sorting stones by size and composition into a grid.",
    "the collapsed base of a space elevator, its severed cable coiled across the surrounding terrain for kilometers.",
    "a sealed cryogenic vault built into the bedrock, its status lights still cycling after some uncounted span of time.",
    "a field of dormant terraforming seeders, half-buried and inert, that should have blanketed the world in engineered soil long ago.",
    "a buried genetic seed vault with a collapsed access shaft, its contents unknown.",
    "a standing obelisk of unidentified alloy on an otherwise featureless plain, untouched by weathering.",
    "an abandoned colony's foundation grid, its structural bones outlining streets and buildings that were never raised.",
    "the fallen anchor of an orbital elevator, its counterweight long gone, embedded deep in the crust.",
    "an automated excavation still cutting the same exhausted mineral seam decades after it ran dry.",
    "a scorched landing zone cleared by a controlled burn that was, by every visible sign, never used.",
    "a shallow crater ringed with fused glass, the unmistakable signature of an old kinetic impact test."
]

ORBITAL_FLAVOR = [ # Flavor text for orbital space only.
    "several dead delivery drones and cargo pods in stable orbits.",
    "a dense debris ring consisting of shattered hull wreckage, dormant escape pods, and thousands of desiccated corpses in identical hazard suits.",
    "a non-sentient drone swarm operating via networked intelligence that actively intercepts passing vessels it flags as requiring extensive repairs.",
    "a primitive interstellar probe caught in orbit.",
    "a derelict light cruiser venting plasma slowly into the surrounding void.",
    "a derelict transport ship whose severed bow trails forty kilometers behind its own stern.",
    "a dense ring of cobalt-ferrite dust, in a thick fog of ionized sodium-potassium gas broadcasting radio signals across a variety of frequencies.",
    "a monolithic stone sculpture of a humanoid hand orbiting with its index finger pointing at magnetic north.",
    "a constellation of large solar mirrors suspended in a rough orbital grid.",
    "a ring of cryostasis pods tumbling freely, their occupants' status unreadable at this range.",
    "the wreckage of a one-sided battle, hulls scattered with no intact vessel among them.",
    "a stationary listening post bristling with dish arrays, all aimed at a patch of empty sky.",
    "a dormant ring-shaped gate of unfamiliar design, its inner surface dark and its function unconfirmed.",
    "a cluster of sealed seed pods tumbling in a slow, decaying orbit, clearly not of local origin.",
    "an automated checkpoint platform still broadcasting docking instructions in a language no database recognizes.",
    "a field of armed kinetic interceptors, live and waiting, with no record of who deployed them.",
    "a solitary dark-hulled relay emitting a tight beam toward the star and nothing else.",
    "the fused, burned-out frame of an old drive test rig, its containment ring melted solid.",
    "a loose formation of repair drones still patrolling a ship that no longer exists."
]

HABITABLE_FLAVOR = [ # Flavor text for any habitable world.
    "microscopic airborne organisms that cause nighttime clouds to glow with blue bioluminescence.",
    "armored, tortoise-like animals that eat raw rock and excrete exceptionally fine glass sand.",
    "high-altitude predators with hyper-thin membranous wings that ride the solar wind.",
    "jagged, forest-like continental structures built from the hard mineral shells of colonial insects.",
    "massive, lighter-than-air floating beasts that use natural gas bladders to filter-feed on atmospheric plankton.",
    "migratory animals that navigate along the electromagnetic hum of buried, abandoned colonial power grids.",
    "tall, iron-rich trees that act as natural lightning rods to pop open their heat-activated seed pods.",
    "lowland basins heavily blanketed in a thick yellow fog caused by massive fungal spore releases.",
    "pack hunters communicating via low-frequency infrasound that induces unexplained dread and vertigo in explorers.",
    "dense subterranean vines that register rapid, targeted growth vectors toward concentrated heat sources.",
    "large, eyeless amphibious predators observed actively tracking the dissolved chemical residue of engine exhaust.",
    "radiotrophic fungi that thrive in irradiated caverns, growing darker and denser near hot mineral veins.",
    "colonial filter-feeders anchored to reef structures, pulsing in loose synchrony with the tides.",
    "burrowing, desiccation-proof invertebrates that survive years dormant between rare rainfalls.",
    "grazing animals with slowly crystallizing shells that eventually fuse into the rock formations they feed on.",
    "blind cave-dwellers that navigate total darkness by sensing faint bioelectric fields in the walls around them.",
    "lichen-like organisms that ride migrating herds, trading camouflage for transport.",
    "herd animals that bury themselves and estivate for years to survive multi-year droughts.",
    "predators whose skin patterns mimic local flora so precisely they vanish at rest.",
    "colonial organisms that briefly link into a single distributed nervous system during seasonal migration.",
    "subterranean dwellers that cultivate bacterial mats in lightless tunnels as a stable food source across generations."
]

# --- Exotic Stellar Phenomena (phenomenonGen.py, generation/phenomena/compact_remnant.py,
# nebulaData.py, supernovaRemnantData.py, roguePlanetData.py) ---
#
# These five phenomena are generated by a separate, rarer, on-demand
# mode (phenomenonGen.py) rather than StarSystem._generate_planets's
# normal per-slot rolls -- see that module's docstring.

# Stellar-mass black hole mass function: real X-ray binary mass
# measurements cluster in a ~5-20 Msun range with a paucity of
# compact objects between the heaviest neutron stars and the lightest
# black holes (the "lower mass gap", Ozel et al. 2010, ApJ 725:1918;
# Farr et al. 2011, ApJ 741:103). Modeled here as a simple uniform draw
# over the well-populated range; INTERMEDIATE_MASS_CHANCE occasionally
# instead draws an intermediate-mass black hole (IMBH), log-uniform over
# ~1e2-1e5 Msun -- the range between stellar and supermassive black
# holes (Greene, Strader & Ho 2020, ARA&A 58:257). Boss (2026-09-30)
# wants them "a smattering (rare)": retuned 2026-10-09 from 2% to 0.1%, which
# gave about 1e5 in a Milky Way-sized galaxy against the handful of observed
# candidates, and on 2026-10-16 (Boss, decision card, GEN.259) to 3e-6, so that
# a Milky Way-sized galaxy, with about 1e8 black holes, holds about 300 of them
# ("hundreds, as the observations suggest") [R]. GEN.563 (Boss, card BH4, 12:34 AM EDT on 2026-10-18) holds that count: the
# stellar-mass black holes now follow the stars' metallicity, but these are not remnants of failed supernovae, so
# they keep the rate they had, 3e-6 of 7e-5 per pc^3 over `REFERENCE_STELLAR_DENSITY_PC3`.
BLACK_HOLE_MASS_RANGE_SOLAR = (5.0, 20.0)
BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR = (1e2, 1e5)
BLACK_HOLE_INTERMEDIATE_MASS_CHANCE = 3.0e-6
BLACK_HOLE_INTERMEDIATE_RATE_PER_STAR = BLACK_HOLE_INTERMEDIATE_MASS_CHANCE * 7.0e-5 / 0.14

# A supermassive black hole sits at the center of every galaxy (generate.add_galactic_nucleus), quiescent or an
# active quasar. Its mass follows the galaxy's bulge or stellar mass (`galaxy.structure.black_hole_mass_solar`,
# GEN.402). The class edge is the usual definition, a million solar masses.
BLACK_HOLE_SUPERMASSIVE_MIN_MASS_SOLAR = 1e6

# A quiescent nucleus accretes far below its Eddington limit: Sgr A*
# shines at ~1e-9 of it (Genzel, Eisenhauer & Gillessen 2010, RvMP
# 82:3121), low-luminosity nuclei up to ~1e-5. Drawn log-uniform; it
# always keeps a (faint) accretion flow.
BLACK_HOLE_SUPERMASSIVE_EDDINGTON_RATIO_RANGE = (1e-9, 1e-6)

# A supermassive black hole's sphere of influence, r = G*M / sigma^2, stands in for the Hill sphere a stellar-mass
# one uses; sigma is the bulge dispersion of the M-sigma relation for its mass (`galaxy.structure`, GEN.402).

BLACK_HOLE_MASS_CLASSES = ("stellar", "intermediate", "supermassive")
"""tuple: `black_holes.mass_class` values (schema v36)."""

# Dimensionless spin parameter a* = c*J/(G*M^2), physically bounded to
# [0, 1) (a*=1 is the extremal Kerr limit). Population-synthesis and
# X-ray-binary spin measurements span this whole range with no strong
# preferred value (Fragos & McClintock 2015, ApJ 800:17), so a uniform
# draw is used rather than a peaked distribution.
BLACK_HOLE_SPIN_RANGE = (0.0, 0.998)

# Chance a generated black hole retains a visible accretion disk (most
# stellar-mass black holes are quiescent/isolated; only a minority are
# actively accreting from a companion or interstellar medium enough to
# be luminous, e.g. Cygnus X-1-like systems).
BLACK_HOLE_ACCRETION_DISK_CHANCE = 0.001

# Where remnants are against where stars are (GEN.132; see
# docs/design/systems/compact-remnant-regions.md). Estimates from the trends there.
# A remnant's layer is this many times as thick as the thin disk: compact
# remnants about 3x the stars' height (Sweeney et al., 1260 pc against
# 306-420 pc), black holes 2.5x (786 pc against 306 pc).
REMNANT_SCALE_HEIGHT_RATIO = {"neutron-star": 3.0, "black-hole": 2.5}
REMNANT_VERTICAL_FACTOR_MAX = 4.0
# Black holes' extra share toward the core (smaller kicks keep them near
# where they formed): this much more at the center, over this scale.
BLACK_HOLE_CORE_EXCESS = 1.0
BLACK_HOLE_CORE_SCALE_PC = 2000.0
# Active radio pulsars' radial profile, Lorimer et al. 2006 (a = 1.9, b = 5.0),
# against the disk's own stars, with the radius it is quoted at (`REFERENCE_RADIUS_SCALE_LENGTHS`) and the disk's scale
# length both the galaxy's own (`galaxy.scales`).
PULSAR_RADIAL_A = 1.9
PULSAR_RADIAL_B = 5.0

# Neutron star mass range in solar masses -- real measured masses cluster
# tightly around ~1.4 Msun with a hard floor near the Chandrasekhar-like
# collapse threshold and a soft ceiling near the maximum mass general
# relativity allows before further collapse to a black hole (Ozel &
# Freire 2016, ARA&A 54:401, compiling radio pulsar timing masses).
NEUTRON_STAR_MASS_RANGE_SOLAR = (1.1, 2.2)

# The initial mass of the stars that end as a white dwarf. A star of 8 solar masses or more ends as a neutron star or a
# black hole, by the windows of `EXPLOSION_FAILURE_ISLANDS_SOL` (GEN.554, `physics.remnant_progenitors`).
WHITE_DWARF_PROGENITOR_SOL = (_boundaries.sol("white_dwarf_progenitor_floor"), _boundaries.sol("core_collapse_limit"))

EXPLOSION_FAILURE_ISLANDS_SOL = ((23.0, 25.0), (35.0, 45.0))
"""tuple: The windows of initial mass, solar masses, in which a star of solar metallicity fails to explode and leaves a black
hole, when little of the massive stars fail: O'Connor and Ott 2011 (ApJ 730, 70) for the Woosley and Heger 2007 progenitors,
"~23-25 Msun and ~35-45 Msun" (their 0 to 7% of the massive stars), the islands of Sukhbold et al. 2016 (ApJ 821, 38),
"islands of 'explodability' in a sea of black hole formation" (GEN.560, docs/design/systems/black-hole-formation-sources.md).
Only these two windows are sourced. A galaxy that fails a smaller share shrinks them about their centres; one that fails a
larger share keeps them and adds the top of the mass range down to the mass that makes up the share (the widest verified
share, 45%, reaches down to about 14.5 solar masses, where Ertl et al. 2016, quoting Ugliano et al. 2012, put the lowest
island) [S, stand-in: the nesting of the windows is not sourced]."""

# Pulsar spin period, in milliseconds, split into the two real observed
# populations (Lorimer 2008, "Binary and Millisecond Pulsars", Living
# Reviews in Relativity 11:8, reviewing the ATNF pulsar catalog): young,
# high-field pulsars spin slower (tens of ms to a few seconds) while
# recycled millisecond pulsars (spun up by past accretion from a binary
# companion) spin far faster (~1-10 ms).
PULSAR_SPIN_PERIOD_MS_RANGE_YOUNG = (16.0, 2000.0)
PULSAR_SPIN_PERIOD_MS_RANGE_MILLISECOND = (1.4, 10.0)

# Surface magnetic field strength, in Gauss, for the same two
# populations (Lorimer 2008) -- young pulsars carry the strong
# fossil fields inherited from core collapse, while millisecond
# pulsars' fields have decayed/been buried over a much longer active
# lifetime.
PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_YOUNG = (1e11, 1e13)
PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_MILLISECOND = (1e8, 1e9)

# Whether a neutron star is actively pulsing follows its age (GEN.153): a
# young radio pulsar stays active about 1e7 years, so the chance it still
# pulses is exp(-age / PULSAR_ACTIVE_LIFETIME_YR); older stars that pulse
# are recycled millisecond pulsars, spun up by a companion, at a small
# fixed share (about 1e5 known beamed radio pulsars among about 5e8
# isolated neutron stars; Lorimer 2008, Living Rev. Relativ. 11:8). The
# result is scaled by PULSAR_RATE_SCALE (1.0 is the physical rate).
PULSAR_ACTIVE_LIFETIME_YR = 1.0e7
PULSAR_RECYCLED_CHANCE = 1.0e-4
PULSAR_RATE_SCALE = 1.0

# Planets of a millisecond pulsar (GEN.317; docs/design/systems/exotic-environments-planets-and-compact-binaries.md section
# 2.1). Whether a pulsar holds any is the remnant satellite table's (0.7% of millisecond pulsars, a young pulsar none);
# these are the origins of the ones that do, from the four known systems (low confidence): the fallback disc of a
# disrupted white dwarf (PSR B1257+12), a captured circumbinary giant (PSR B1620-26 b, needs a globular cluster) and
# an ablated companion (PSR J1719-1438 b, J2322-2650 b). Only the disc origin is built; the other two draw and give
# the pulsar no planets until their kinds exist (the giant GEN.163, the sub-stellar companion with the pair builder).
PULSAR_PLANET_ORIGIN_WEIGHTS = {"disc": 25.0, "captured-giant": 25.0, "ablated-companion": 50.0}
PULSAR_PLANET_BUILT_ORIGINS = ("disc",)

# The disc origin: 1 to 3 planets in 0.1 to 0.6 AU, masses log-uniform from 0.01 to 10 Earth masses (the least
# massive known pulsar planet is 0.02, the most 4.3), each period near a low-order ratio to the one inside it (PSR
# B1257+12's outer pair sits near 3:2). A planet's period is its ratio times the inner one's, scattered by a normal
# deviation of this fraction; a mass over what the carbon-rich airless class holds is not drawn.
PULSAR_DISC_PLANET_COUNT = (1, 3)
PULSAR_DISC_DISTANCE_AU_RANGE = (0.1, 0.6)
PULSAR_DISC_MASS_EARTH_RANGE = (0.01, 10.0)
PULSAR_DISC_PERIOD_RATIOS = (4.0 / 3.0, 3.0 / 2.0, 5.0 / 3.0, 2.0)
PULSAR_DISC_PERIOD_RATIO_SCATTER = 0.01
PULSAR_DISC_MAX_RADIUS_KM = 10000.0  # the largest disc planet, at the mean density of a terrestrial body (the old airless class ceiling)

# The pulsar wind's momentum flux on a body is XI * Edot / (4 pi r^2 c) for a spin-down power Edot and a distance r
# (the same order of unity the design note takes).
PULSAR_WIND_XI = 1.0

# Magnetars (GEN.153), a subtype of neutron star: B of 1e14 to 1e15 G and a
# spin of 2 to 12 s, active only about 1e4 years (an older one has decayed
# to an ordinary high-field star), born in 1 to 10% of core collapses,
# default 10% (Olausen & Kaspi 2014, ApJS 212:6; Gill & Heyl 2007). The
# birth lottery applies to the cores of young supernova remnants; the
# free-standing ones are scatter rows (below).
MAGNETAR_ACTIVE_LIFETIME_YR = 1.0e4
MAGNETAR_BIRTH_FRACTION = 0.10
MAGNETAR_MAGNETIC_FIELD_GAUSS_RANGE = (1e14, 1e15)
MAGNETAR_SPIN_PERIOD_MS_RANGE = (2000.0, 12000.0)

# Magnetars as a scatter class of neutron stars (GEN.278, design note
# later-phenomena-scatter-rules.md section 3): about 100 active in the
# Galaxy against about 1e9 neutron stars (Olausen & Kaspi 2014; 30 to 600),
# so 1e-7 of them. They follow the young population, not the stars, so the
# rate per star is divided by the young population's share of all stars in
# a Milky Way-like galaxy (the integral of the young density over that of
# every component, 0.16%: since GEN.469 the disk's populations have their shares from its gas history (a galaxy of
# 5e10 solar masses; a lighter one forms its stars later and holds more young ones) and since GEN.544 the exponential
# height (it was 0.19% over the isothermal one); test_gen_magnetar_scatter checks it); a galaxy
# with a different young share holds proportionally more or fewer.
MAGNETAR_PER_NEUTRON_STAR = 1.0e-7
MAGNETAR_YOUNG_STAR_SHARE = 0.0016
# A scattered magnetar's age: from the first century to the end of its active life.
MAGNETAR_AGE_RANGE_GY = (1.0e-7, MAGNETAR_ACTIVE_LIFETIME_YR / 1e9)

# Neutron star surface temperature range in Kelvin, from thermal
# X-ray emission measurements of young-to-middle-aged isolated
# neutron stars (Potekhin et al. 2020, A&A Rev 28:3, cooling-curve
# review) -- hotter shortly after formation, cooling over ~10^5-10^6
# years toward the low end of this range.
NEUTRON_STAR_SURFACE_TEMPERATURE_RANGE_K = (5e4, 3e6)

# Age since the core-collapse supernova that formed the compact
# remnant, in billions of years -- capped well below
# UNIVERSE_AGE_GY since these are always described as galactic
# (Population I/II), not primordial, remnants. The floor is a magnetar's
# active life (MAGNETAR_ACTIVE_LIFETIME_YR): a free neutron star is never
# young enough to be one (GEN.278); the magnetars are scatter rows and the
# cores of young supernova remnants, which carry their own ages.
COMPACT_REMNANT_AGE_RANGE_GY = (1.0e-5, 10.0)

# The time the node phases of a scatter drawn before any orbit update hold at (GEN.262): 2000-01-01 00:00 UTC. The
# orbit clock does not exist until the first update, and a fixed origin keeps one galaxy seed's members identical in
# every run; the nightly update never writes a node, and a member's place at time t advances from its phase here.
SCATTER_PHASE_ORIGIN_UNIX = 946684800.0

# --- Nebulae and supernova remnants (nebulaData.Nebula,
# supernovaRemnantData.SupernovaRemnant) ---
#
# Nebulae fall into five families (Osterbrock & Ferland 2006, "Astrophysics
# of Gaseous Nebulae and Active Galactic Nuclei", 2nd ed.): diffuse
# interstellar gas, emission (H II) regions ionized by hot young stars,
# reflection nebulae scattering a nearby star's light off dust, planetary
# nebulae expelled by a dying low/intermediate-mass star, and dark
# (molecular) clouds seen in silhouette. `nebulae.nebula_type` stores the
# family; the class letter below says what is actually in the cloud.
NEBULA_FAMILIES = {
    "diffuse": {
        "composition": "thin, warm interstellar gas, mostly hydrogen and helium",
        "formation_cause": "the galaxy's general interstellar medium, stirred and heated by starlight and old supernova blasts",
    },
    "emission": {
        "composition": "ionized hydrogen (H II) glowing under ultraviolet radiation from nearby hot young stars",
        "formation_cause": "ultraviolet radiation from newly-formed O and B class stars ionizing the surrounding hydrogen cloud",
    },
    "reflection": {
        "composition": "fine interstellar dust scattering the blue light of an adjacent bright star",
        "formation_cause": "a bright star passing through or forming within a dense dust cloud, illuminating it by reflection rather than ionization",
    },
    "planetary": {
        "composition": "ionized gas shells of hydrogen, helium, oxygen, and nitrogen expelled by a dying star",
        "formation_cause": "a low-to-intermediate-mass star shedding its outer envelope at the end of its asymptotic giant branch phase, exposing a hot white dwarf core",
    },
    "dark": {
        "composition": "dense molecular hydrogen and cold dust grains, opaque to visible light",
        "formation_cause": "a cold, dense molecular cloud that has not yet collapsed to form stars, visible only in silhouette against background starlight",
    },
}

# Nebula and supernova remnant classes, one letter each like
# PLANET_CLASSES (GEN.11; reasoning, sources and the family tables in
# docs/design/systems/nebula-and-asteroid-field-classes.md). I and O are unused so
# they aren't read as 1 and 0; Y and Z are reserved. A-Q and X (the Wolf-Rayet ring, GEN.243) are nebulae, R-W are
# supernova remnants ("stellar remnants" means supernova remnants only,
# Boss 2026-09-30; compact objects keep their own tables).
#
# Each class: `name`, `family` (a NEBULA_FAMILIES key, or
# "supernova-remnant"), `species` (dominant contents), `density_range_cm3`
# (particle density nH, log-uniform), `temperature_range_k` (log-uniform),
# `extinction_range_av` (optical extinction in magnitudes, log-uniform
# unless the low end is 0), `center` (the central-object rule of GEN.10
# builds on), and `frequency` (relative weight within its family group).
# A class with `host_only` is never drawn free (by family or frequency): only a star of its kind lights one (X).
# Nebulae also carry `radius_range_ly`; a remnant's radius comes from its
# age (Sedov-Taylor, below), so remnants carry `morphology` (the Vink 2012
# shell/plerion/composite shape), `age_range_years` and `compact` (which
# compact remnants the class allows: "neutron_star", "black_hole", None).
NEBULA_CLASSES = {
    "A": {"name": "Diffuse neutral cloud", "family": "diffuse",
          "species": "neutral hydrogen and helium with trace ions",
          "radius_range_ly": (10.0, 150.0), "density_range_cm3": (0.1, 10.0),
          "temperature_range_k": (6000.0, 10000.0), "extinction_range_av": (0.0, 0.1),
          "center": "none", "frequency": 3.0},
    "B": {"name": "Diffuse ionized gas", "family": "diffuse",
          "species": "ionized hydrogen and free electrons",
          "radius_range_ly": (20.0, 200.0), "density_range_cm3": (0.1, 1.0),
          "temperature_range_k": (8000.0, 10000.0), "extinction_range_av": (0.0, 0.1),
          "center": "none nearby", "frequency": 2.0},
    "C": {"name": "Compact H II region", "family": "emission",
          "species": "ionized gas still inside its dusty birth cloud",
          "radius_range_ly": (0.3, 3.0), "density_range_cm3": (1e3, 1e4),
          "temperature_range_k": (8000.0, 12000.0), "extinction_range_av": (1.0, 10.0),
          "center": "one young O or early-B star", "frequency": 1.0},
    "D": {"name": "Classical H II region", "family": "emission",
          "species": "H+, e-, O2+, N+ and S+",
          "radius_range_ly": (10.0, 100.0), "density_range_cm3": (10.0, 1e3),
          "temperature_range_k": (8000.0, 12000.0), "extinction_range_av": (0.0, 1.0),
          "center": "a small O/B cluster", "frequency": 2.0},
    "E": {"name": "Giant H II complex", "family": "emission",
          "species": "H+, e-, O2+, N+ and S+ in many nested shells",
          "radius_range_ly": (50.0, 200.0), "density_range_cm3": (10.0, 1e3),
          "temperature_range_k": (8000.0, 12000.0), "extinction_range_av": (0.0, 1.0),
          "center": "a young massive cluster", "frequency": 0.5},
    "F": {"name": "Reflection nebula", "family": "reflection",
          "species": "silicate grains, PAHs, carbon soot and ices in neutral gas",
          "radius_range_ly": (1.0, 20.0), "density_range_cm3": (1e2, 1e3),
          "temperature_range_k": (10.0, 100.0), "extinction_range_av": (0.5, 3.0),
          "center": "a B or A star", "frequency": 2.0},
    "G": {"name": "Emission-reflection nebula", "family": "emission",
          "species": "ionized core inside a dusty rim",
          "radius_range_ly": (2.0, 20.0), "density_range_cm3": (1e2, 1e3),
          "temperature_range_k": (50.0, 10000.0), "extinction_range_av": (0.5, 3.0),
          "center": "an early-B star", "frequency": 1.0},
    "X": {"name": "Wolf-Rayet ring nebula", "family": "emission",
          "species": "ionized hydrogen and helium in a thin swept-up shell around a hollow interior",
          "radius_range_ly": (6.5, 114.0), "density_range_cm3": (10.0, 1e3),
          "temperature_range_k": (8000.0, 12000.0), "extinction_range_av": (0.0, 1.0),
          "center": "a Wolf-Rayet star (required)", "frequency": 1.0, "host_only": True},
    "H": {"name": "Young planetary nebula", "family": "planetary",
          "species": "ionized hydrogen, carbon, nitrogen, oxygen and neon",
          "radius_range_ly": (0.1, 0.5), "density_range_cm3": (1e4, 1e5),
          "temperature_range_k": (10000.0, 20000.0), "extinction_range_av": (0.0, 0.3),
          "center": "a hot central star", "frequency": 1.0},
    "J": {"name": "Evolved planetary nebula", "family": "planetary",
          "species": "thinning ionized hydrogen, carbon, nitrogen, oxygen and neon",
          "radius_range_ly": (0.5, 3.0), "density_range_cm3": (1e2, 1e3),
          "temperature_range_k": (10000.0, 20000.0), "extinction_range_av": (0.0, 0.3),
          "center": "a white dwarf", "frequency": 2.0},
    "K": {"name": "Carbon-rich planetary nebula", "family": "planetary",
          "species": "ionized gas with carbon-rich dust and PAHs",
          "radius_range_ly": (0.1, 2.0), "density_range_cm3": (1e2, 1e4),
          "temperature_range_k": (10000.0, 20000.0), "extinction_range_av": (0.0, 0.3),
          "center": "a central star from a 1.5-3 Msun progenitor", "frequency": 1.0},
    "L": {"name": "Nitrogen-rich bipolar planetary nebula", "family": "planetary",
          "species": "nitrogen- and helium-enriched gas in bipolar lobes",
          "radius_range_ly": (0.2, 3.0), "density_range_cm3": (1e2, 1e4),
          "temperature_range_k": (10000.0, 20000.0), "extinction_range_av": (0.0, 0.3),
          "center": "a central star from a 3-8 Msun progenitor, often binary", "frequency": 1.0},
    "M": {"name": "Giant molecular cloud", "family": "dark",
          "species": "H2, He, CO, PAHs, silicates",
          "radius_range_ly": (20.0, 150.0), "density_range_cm3": (1e2, 1e6),
          "temperature_range_k": (10.0, 30.0), "extinction_range_av": (10.0, 100.0),
          "center": "none; embedded young clusters possible", "frequency": 1.0},
    "N": {"name": "Dark cloud", "family": "dark",
          "species": "H2, CO, cold dust",
          "radius_range_ly": (1.0, 50.0), "density_range_cm3": (1e3, 1e5),
          "temperature_range_k": (10.0, 30.0), "extinction_range_av": (5.0, 50.0),
          "center": "none", "frequency": 3.0},
    "P": {"name": "Bok globule", "family": "dark",
          "species": "H2, CO, organics and dust",
          "radius_range_ly": (0.3, 3.0), "density_range_cm3": (1e4, 1e5),
          "temperature_range_k": (10.0, 30.0), "extinction_range_av": (5.0, 50.0),
          "center": "none or one protostar", "frequency": 2.0},
    "Q": {"name": "Star-forming core", "family": "dark",
          "species": "H2 with outflows, Herbig-Haro jets and masers",
          "radius_range_ly": (0.1, 1.0), "density_range_cm3": (1e5, 1e6),
          "temperature_range_k": (10.0, 30.0), "extinction_range_av": (10.0, 100.0),
          "center": "embedded protostars", "frequency": 1.0},
    "R": {"name": "Young ejecta-dominated remnant", "family": "supernova-remnant",
          "species": "Fe, Si, S and O ejecta",
          "morphology": "shell", "age_range_years": (50.0, 3000.0),
          "compact": ("neutron_star", "black_hole", None),
          "density_range_cm3": (0.1, 100.0),
          "temperature_range_k": (1e6, 1e7), "extinction_range_av": (0.0, 0.1),
          "center": "a neutron star or black hole", "frequency": 1.0},
    "S": {"name": "Shell remnant", "family": "supernova-remnant",
          "species": "shocked interstellar gas and ejecta",
          "morphology": "shell", "age_range_years": (50.0, 100000.0),
          "compact": ("neutron_star", "black_hole", None),
          "density_range_cm3": (0.1, 10.0),
          "temperature_range_k": (1e6, 1e7), "extinction_range_av": (0.0, 0.1),
          "center": "a neutron star, black hole or nothing", "frequency": 3.0},
    "T": {"name": "Pulsar wind nebula (plerion)", "family": "supernova-remnant",
          "species": "relativistic electrons in a magnetic field",
          "morphology": "plerion", "age_range_years": (50.0, 20000.0),
          "compact": ("neutron_star",),
          "density_range_cm3": (0.01, 1.0),
          "temperature_range_k": (1e4, 1e6), "extinction_range_av": (0.0, 0.1),
          "center": "a pulsar (required)", "frequency": 1.0},
    "U": {"name": "Composite remnant", "family": "supernova-remnant",
          "species": "a shocked shell around a pulsar wind nebula",
          "morphology": "composite", "age_range_years": (1000.0, 30000.0),
          "compact": ("neutron_star",),
          "density_range_cm3": (0.1, 10.0),
          "temperature_range_k": (1e6, 1e7), "extinction_range_av": (0.0, 0.1),
          "center": "a pulsar", "frequency": 1.0},
    "V": {"name": "Old radiative remnant", "family": "supernova-remnant",
          "species": "a cooling shell merging with the interstellar medium",
          "morphology": "shell", "age_range_years": (20000.0, 100000.0),
          "compact": ("neutron_star", "black_hole", None),
          "density_range_cm3": (1.0, 100.0),
          "temperature_range_k": (1e4, 1e6), "extinction_range_av": (0.0, 0.1),
          "center": "a neutron star far off-center, or nothing", "frequency": 2.0},
    "W": {"name": "Thermonuclear remnant", "family": "supernova-remnant",
          "species": "iron-rich ejecta with no hydrogen",
          "morphology": "shell", "age_range_years": (50.0, 100000.0),
          "compact": (None,),
          "density_range_cm3": (0.1, 10.0),
          "temperature_range_k": (1e6, 1e7), "extinction_range_av": (0.0, 0.1),
          "center": "nothing (a Type Ia supernova leaves no core)", "frequency": 1.0},
}

SUPERNOVA_REMNANT_MORPHOLOGIES = ("shell", "plerion", "composite")
"""tuple: The remnant shapes `supernova_remnants.morphology` allows (Vink
2012, A&A Rev 20:49): a limb-brightened shell, a plerion (filled by a
pulsar wind nebula, like the Crab), or a composite of both. Each remnant
class in `NEBULA_CLASSES` fixes its own."""

# Sedov-Taylor phase expansion, R = C * (E51 / n)^(1/5) * t^(2/5): the
# self-similar blast-wave solution for a remnant that has swept up enough
# interstellar mass to decelerate from its initial free-expansion phase
# (Taylor 1950, Proc. Roy. Soc. A 201:159; Sedov 1959, "Similarity and
# Dimensional Methods in Mechanics"). C = 1.03 ly for E51 = 1 (an explosion
# energy of 1e51 erg), n in cm^-3 and t in years: 5.0 pc at 1,000 yr in
# n = 1 (GEN.151). The radius stops growing at the end of the Sedov phase,
# where radiative cooling begins and the snowplow starts, R_pds = 45.7 ly
# * E51^(2/7) * n^(-3/7) (14.0 pc, Cioffi, McKee & Bertschinger 1988, ApJ
# 334:252).
SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY = 1.03
SEDOV_TAYLOR_TIME_EXPONENT = 2 / 5
SEDOV_TAYLOR_DENSITY_EXPONENT = -1 / 5
SEDOV_TAYLOR_ENERGY_E51 = 1.0
SNOWPLOW_ONSET_RADIUS_COEFFICIENT_LY = 45.66
SNOWPLOW_ONSET_ENERGY_EXPONENT = 2 / 7
SNOWPLOW_ONSET_DENSITY_EXPONENT = -3 / 7

# Supernova remnant age range, in years -- from newly-formed (a visible
# shell needs at least some expansion) out to the point (order 10^5
# years) real remnants fade into the general interstellar medium and
# are no longer identifiable as such (Vink 2012).
SUPERNOVA_REMNANT_AGE_RANGE_YEARS = (50.0, 100000.0)

# A Type Ia (thermonuclear disruption of a white dwarf, no compact
# remnant left behind) vs. core-collapse (massive-star death, leaves a
# neutron star or black hole) progenitor, and the chance a core-collapse
# remnant's compact core is still detectable within it -- real core-
# collapse supernovae outnumber Type Ia roughly 3:1 in star-forming
# galaxies (Li et al. 2011, MNRAS 412:1441).
SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE = 0.25
SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE = 0.6

# Of a core-collapse supernova's visible compact remnant, the chance it is
# specifically a black hole rather than a neutron star -- most core-
# collapse progenitors (roughly 8-20 Msun) leave a neutron star, while only
# the most massive (>~20-25 Msun) leave a black hole (Heger et al. 2003,
# ApJ 591:288, "How Massive Single Stars End Their Life").
SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE = 0.15

SUPERNOVA_KICK_SPEED_RANGE_KMS = {"neutron_star": (100.0, 700.0), "black_hole": (20.0, 200.0)}
"""dict: Birth-kick speed of a core-collapse remnant's compact core, km/s,
drawn log-uniformly -- pulsars move at ~400 km/s on average (Hobbs et al.
2005, MNRAS 360:974); black holes get smaller kicks (Repetto et al. 2012,
MNRAS 425:2799). Times the remnant's age, it puts the core off the
remnant's center: a young remnant's core sits near the middle, an old
remnant's can have left it entirely (class V's "far off-center")."""

NEBULA_HOST_RULES = (
    ("O", 0, 9, ("C", "D", "E"), 0.25),
    ("B", 0, 2, ("C", "G"), 0.25),
    ("B", 3, 9, ("F",), 0.05),
    ("A", 0, 9, ("F",), 0.01),
    ("W", 0, 9, ("X",), 0.4),
)
"""tuple: Nebulae a sector grows around its own stars (GEN.10), as
`(spectral letter, lowest subclass, highest subclass, classes, chance)`
for a main-sequence primary. Only O and early-B stars put out enough
ultraviolet below 91.2 nm to ionize hydrogen (Osterbrock & Ferland 2006),
and a catalogued H II region is about one per four O stars: about 8,000
are known against some 270,000 O stars (GEN.204; the nebula density study,
docs/design/systems/nebula-density-vs-reality.md). The `chance` is that the star
still sits in its natal cloud (the one number the catalogues fix, 0.25 for
O and B0-B2 stars), and whether it hosts a class is then decided by the
star's own age against `NEBULA_CLASS_MAX_AGE_MYR`: a B1V of 20 Myr is too
old for a class G region, an O star always young enough. The rows used to
carry a chance per spectral type (B0-B2 0.1); the age test now thins B0-B2
to about 0.06 of them with no number of its own (GEN.150). A single B0-B2
star cannot ionize 10 to 100 ly (one B2V ionizes about 0.8 ly), so class D
is no longer one of their classes, and G needs an early-B star, so it left
the B3-B9 row. Later B and A stars light dust without ionizing it, a
reflection nebula (F) when a dusty cloud happens to be near (van den Bergh
1966, AJ 71:990 catalogs ~150 within a few kpc); their chance is that of
the cloud. The class is drawn among the classes the star is young enough
for, by NEBULA_CLASSES frequency."""

NEBULA_CLASS_MAX_AGE_MYR = {"C": 1.0, "D": 10.0, "E": 20.0, "F": 50.0, "G": 5.0, "X": 10.0}
"""dict: The oldest a host star can be and still sit in a nebula of the class, Myr
(design note, "Lifetimes by class": C 0.1 to 1, D 3 to 10, E 5 to 20, F 1 to 50, G 1
to 5; X: the longest a Wolf-Rayet star's life runs, 8.5 Myr at the threshold mass, GEN.243).
Only the upper limits are used: the lower ones come out of the physics, a
young region has not grown to the class's radius (`physics.hii`)."""

NEBULA_CLUSTER_MASS_RANGE_SOL = {"D": (300.0, 1.0e4), "E": (1.0e4, 1.0e5)}
"""dict: The mass of the cluster behind a class D or E region, solar masses, log-uniform:
a small O/B cluster of 300 to 10,000 Msun (log Q 49.2 to 50.7) and a young massive
cluster of 10,000 to 100,000 (50.7 to 51.7). The other ionized classes are one star's
(design note, "Q needed by class" and the cluster mass table)."""

NEBULA_SEED_MASS_RANGE_SOL = {
    "C": (11.0, 60.0), "D": (18.0, 150.0), "E": (18.0, 150.0), "F": (2.0, 8.5), "G": (11.0, 25.0),
}
"""dict: The star a nebula built with no host star of its own is given (a nebula the
phenomenon mode makes), solar masses, from the IMF: one O or early-B star for
class C, an O star for D and E, an early B star for G, a B or A star for F."""

PLANETARY_NEBULA_COMPANION_CHANCE = 0.20
"""float: The share of planetary nebulae whose central star has a close companion (the class table says one is
"often" there; GEN.320)."""

NEBULA_CLOUD_PROTOSTARS = {"P": {"chance": 0.30}, "Q": {"mean": 2.0}}
"""dict: The protostars a field cloud holds (GEN.320): a Bok globule (P) one in `chance` of globules, a star-forming
core (Q) a zero-truncated Poisson count of mean `mean` (one to four are the usual; the draw is capped at 4)."""
NEBULA_CLOUD_PROTOSTAR_MAX = 4
NEBULA_CLOUD_CLUMP_FRACTION = 1.0 / 3.0
"""float: The protostars of a core are scattered about its centre with a Gaussian of this fraction of its radius
[I: "clumped near the centre"], redrawn if one falls outside the cloud."""

NEBULA_CLUSTER_HALF_MASS_PC = {"D": 1.0, "E": 2.0}
"""dict: The half-mass radius, parsecs, of the Plummer sphere the ionizing members of a class D or E region are placed in
(GEN.321; a small O and B cluster, a young massive one; the design note's 1 and 2 pc), clipped to the nebula's radius."""

HABING_FLUX_ERG_CM2_S = 1.6e-3
"""float: One Habing unit (G0), the far-ultraviolet flux that is the yardstick of external photoevaporation
(nebula-and-asteroid-field-classes.md, "Photoevaporation and G0"; arXiv 1909.04093)."""

CLOUD_FUV_FRACTION_BY_TEFF_K = {12000.0: 0.16, 20000.0: 0.46, 30000.0: 0.54, 40000.0: 0.45}
"""dict: Teff, kelvin, to f_FUV, the share of a hot star's light in the 6 to 13.6 eV band, linear between the rows and
held at the end rows beyond them (the design note's four values, [C]); the FUV flux of a host is `f_FUV L / (4 pi d^2)`."""

CLOUD_G0_TIER_EDGES = (1.0e2, 1.0e3, 1.0e4)
"""tuple: The G0 at which the radiation tier rises by one (design note, "Rule table": under 1e2 adds nothing; 1e3 is
where mass loss is significant and 1e4 is where outer discs are destroyed in about 1 Myr [R, Johnstone 1998; Haworth 2018])."""

CLOUD_HOST_MIN_DISTANCE_PC = 0.01
"""float: A host star is never taken closer than this to a system, parsecs, so the sum over a cluster placed at the
cloud's centre stays finite for a system near it. A stand-in: the hosts of a seeded cloud have no stored positions yet."""

CLOUD_CLASS_RADIATION_TIER = {
    "C": 1, "G": 1,
    "H": 3, "J": 2, "K": 2, "L": 2,
    "R": 3, "S": 2, "T": 3, "U": 2, "V": 1, "W": 1,
}
"""dict: The radiation tier a cloud of the class adds before any host's G0 (design note, "Rule table"): the planetary
nebulae's ultraviolet and the young remnants' cosmic rays and X-rays. A class not listed adds 0; D, E and X rise by their
hosts' G0 alone and M, N and F by the heliopause tier, which is the habitability index's. The tier a system gets is
the larger of this and the G0 tier. Design values [R, design]."""

CLOUD_ENRICHMENT_CLASSES = ("D", "E", "R", "S", "T", "U", "V", "W")
"""tuple: The classes whose systems draw an 26Al enrichment `k` (design note, "Rule table": a system inside D or E, or in
a remnant)."""

CLOUD_ENRICHMENT_FIELD_MEDIAN = 0.1
"""float: The median of `k`, the 26Al/27Al ratio in units of `CLOUD_ENRICHMENT_REFERENCE_RATIO`, for a system without a
recent nearby source: the galactic mean 6e-6 against the Solar System's 5.25e-5 [C, design note]."""

CLOUD_ENRICHMENT_RICH_MEDIAN = 1.0
"""float: The median of `k` for a system a massive star or supernova enriched."""

CLOUD_ENRICHMENT_SIGMA_DEX = 0.5
"""float: The spread of `k`, dex, log-normal (design values)."""

CLOUD_ENRICHMENT_RICH_CHANCE = {"D": 0.15, "E": 0.3, "R": 0.5, "S": 0.5, "U": 0.5}
"""dict: The chance a system in a cloud of the class has the rich median rather than the field one (design values)."""

PLANETARY_NEBULA_CENTRAL_STAR_TYPES = ("O3VII", "O5VII", "O7VII", "O9VII", "B0VII")
"""tuple: The central star a planetary nebula is generated around (GEN.10):
the exposed hot core of a dying 0.8-8 Msun star, 30,000 K and
up, already a white dwarf in the Yerkes scheme (class VII)."""

# --- Rogue Planets & Interstellar Comets (roguePlanetData) ---
#
# Free-floating ("rogue"/nomad) planet mass range, in Jupiter masses --
# microlensing surveys (Sumi et al. 2011, Nature 473:349) found
# Jupiter-mass free-floating objects roughly twice as common as main-
# sequence stars in the galaxy, with subsequent analysis (Mroz et al.
# 2017, Nature 548:183) revising the abundance down but still finding a
# genuine population spanning sub-Earth to super-Jupiter masses
# (Strigari et al. 2012, MNRAS 423:1856, order-of-magnitude population
# estimate). This range covers the terrestrial-to-giant span those
# surveys probe.
# Mass bins (Boss's research, 2026-09-30; docs/design/systems/
# interstellar-object-rates.md): a rogue planet draws a bin by its
# per-star rate, then a mass inside it. Since GEN.45 the rates and the draw
# inside a bin both follow one mass function, dN/dlogM ~ M^-0.65 from 0.1
# Earth masses to 13 Jupiter masses (Boss's research, 2026-10-01, which
# writes it dN/dM ~ M^-0.65; read per unit mass that would make about 87%
# of rogues gas giants, against its own aim of keeping them rare, so it is
# read per logarithm of mass). The total is ROGUE_PLANET_RATE_PER_STAR, 5.8 per
# star (the bins sum to it). That gives:
#   - terrestrial 0.1-2 Earth masses, about 5.0 per star (2-10; Johnson
#     et al. 2020, AJ 160:123; Mroz et al. 2020, ApJL 903:L11);
#   - sub-Neptune / ice giant 2-20, about 0.64 per star (Sumi et al. 2023
#     find Neptune-mass candidates but give no rate);
#   - Saturn-class 20 Earth masses to 1 Jupiter mass, about 0.15;
#   - Jupiter-mass 1-13 Jupiter masses, about 0.025 per star (under Mroz
#     et al. 2017's upper limit of 0.25, Nature 548:183).
# So about 96% of rogues are not gas giants and 4% are (past
# ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER). The old hand-set rates
# (5, 1, 0.25, 0.25) with a log-uniform mass in each bin gave 91% and 9%,
# twice the reference's share of gas giants.
# Above 13 Jupiter masses is a brown dwarf (ROGUE_BROWN_DWARF_*).
ROGUE_PLANET_RATE_PER_STAR = 5.8
ROGUE_PLANET_MASS_FUNCTION_SLOPE = 0.65
_ROGUE_PLANET_BIN_EDGES_EARTH = {
    "terrestrial": (0.1, 2.0),
    "sub-neptune": (2.0, 20.0),
    "saturn": (20.0, 317.8),
    "jupiter": (317.8, _boundaries.jupiter("deuterium_burning_limit") * 317.8),
}


def _rogue_planet_mass_bins():
    slope = ROGUE_PLANET_MASS_FUNCTION_SLOPE
    low_all = min(low for low, _high in _ROGUE_PLANET_BIN_EDGES_EARTH.values())
    high_all = max(high for _low, high in _ROGUE_PLANET_BIN_EDGES_EARTH.values())
    total = low_all ** -slope - high_all ** -slope
    return {name: (low, high, ROGUE_PLANET_RATE_PER_STAR * (low ** -slope - high ** -slope) / total)
            for name, (low, high) in _ROGUE_PLANET_BIN_EDGES_EARTH.items()}


ROGUE_PLANET_MASS_BINS = _rogue_planet_mass_bins()
"""dict: bin name -> `(min_mass_earth, max_mass_earth, per_star_rate)`.
The rates' sum is the rogue-planet rate per star (`ROGUE_PLANET_RATE_PER_STAR`)."""

ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER = (_boundaries.jupiter("deuterium_burning_limit"), _boundaries.jupiter("brown_dwarf_upper_limit"))
"""tuple: A free-floating brown dwarf, between the deuterium-burning
(~13 Mjup) and hydrogen-burning (~80 Mjup) limits. Stored as a
`rogue_planets` row with `mass_bin = 'brown-dwarf'` (schema v37)."""

ROGUE_PLANET_MASS_BIN_CHOICES = tuple(ROGUE_PLANET_MASS_BINS) + ("brown-dwarf",)
"""tuple: Every `rogue_planets.mass_bin` value."""

ROGUE_LEVEL_CEILINGS = {"None": 0.0, "Few": 2.5, "Normal": 4.3, "Many": 10.0}
"""dict: The Rogues level -> the most rogue planets per star a sector may hold after its top-up (GEN.1125; Boss accepted these numbers 5:39 PM EDT 2026-10-30; the page shows only the
words). Source tier T9, tweak level L3 (GEN.725 row 6.86:RT1)."""

ROGUE_TOPUP_PLACEMENT = 1.0
"""float: The share of the sector's edge, centred on the sector, inside which the top-up places a rogue uniformly (1 is the whole sector). Boss wants to play with it
(`region-object-store.md` 8.5). Source tier T9, tweak level L3 (GEN.725 row 6.86:RT2)."""

# Above this mass (in Jupiter masses, 10 Earth masses), a generated rogue planet is a giant rather than terrestrial: the mass
# a core must reach to hold a hydrogen and helium envelope (`GIANT_CORE_CRITICAL_MASS_EARTH`), so a rogue is a giant exactly where
# a bound planet is (GEN.1133; a bound rocky planet never passes it, and a rogue's rock does not either). Between this and
# `GIANT_RUNAWAY_MIN_MASS_EARTH` it is an ice giant (class I by its conditions); the rogue's old line of 0.05 (about 16 Earth masses,
# "where the Solar System's ice giants sit") let 10 to 16 Earth-mass rock come out at up to 17,600 km, past class S.
ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER = (
    _kg_to_mjup(_mearth_to_kg(GIANT_CORE_CRITICAL_MASS_EARTH)))

# --- Rogue planet surface conditions (rogueSurface; Boss's research,
# 2026-10-01: Stevenson 1999, Burrows et al. 2001, Guillot 2005, Abbot &
# Switzer 2011, Turcotte & Schubert 2014, Pierrehumbert 2010) ---
#
# With no star, a rogue planet's only heat is its own: radioactive decay
# and leftover formation heat in a rocky body, slow Kelvin-Helmholtz
# contraction in a giant. Values marked "a default" are this generator's
# picks where the research gives no number.

ROGUE_PLANET_AGE_RANGE_GY = (0.1, 10.0)
"""tuple: A rogue planet's age, drawn uniformly: the thin disk has formed
stars (and so ejected planets) at a roughly steady rate for ~10 Gyr (a
default)."""

ROGUE_EARTH_REFERENCE_AGE_GY = 4.5
"""float: The Earth's age, the reference a giant's cooling law is scaled from (`ROGUE_GIANT_COOLING_TIME_EXPONENT`)."""

RADIOGENIC_SILICATE_MASS_FRACTION = 0.675
"""float: The heat-producing (bulk silicate) share of a rocky body's mass: the Earth's 1 minus its core mass fraction 0.325 (part E
2.8e). The stand-in for a body with no stored composition; `radiogenic.rock_of` reads the drawn one (GEN.428)."""

GCR_SOLAR_MODULATION_RANGE = (0.5, 1.64)
"""tuple: The cosmic-ray dose as a fraction of its mean over the solar cycle, drawn log-uniformly per body (GEN.453): the
Mars dose equivalent runs from 0.35 to 1.15 mSv/day around a mean of 0.7 (Guo 2015, Hassler 2014; a total spread of
3.3). The same field at another star differs in the same way, so the draw stands for the host's modulation level."""


RADIOGENIC_OFFSET_SIGMA_DEX = 0.22
"""float: The scatter of a host-less body's thorium and uranium against the solar scale, dex: chosen so the 1 to 99 percent range of the
Hypatia catalogue's Eu/Mg, plus or minus 0.5 dex, is 2.33 sigma (Nimmo et al. 2020; radiation-model.md R15). The stand-in for a body
whose host star has no stored composition (a rogue planet)."""

RADIOGENIC_OFFSET_LIMIT_DEX = 0.5
"""float: The truncation of that offset, dex, host-less or from the host star (radiation-model.md R15)."""

RADIOGENIC_K_DEPLETION_RANGE = (0.03, 1.0)
"""tuple: The potassium depletion against uranium as the CI chondrites have it (K/U 67407), drawn log-uniformly. The Earth's 0.177 lies
inside it; the range itself is a placeholder tied to the large K differences of Solar System bodies (Nimmo et al. 2020), not a measured
spread (radiation-model.md R15, part E 3.3)."""

INTERIOR_VISCOSITY_PREFACTOR_RANGE_PA_S = (9.3e7, 1.9e11)
"""tuple: The mantle viscosity prefactor `eta0` of `eta = eta0 exp(Ea / R T)`, Pa s, drawn uniformly (GEN.1185). Guimond, Rudge and Shorttle 2022 take 2.6e10 to 5.3e13 with an
activation energy of 200 to 300 kJ/mol (arXiv:2201.05636 section 2.2.3; READ); the activation range moved up to the diffusion-creep values of the laboratory
(`INTERIOR_ACTIVATION_RANGE_KJ_MOL`), and these bounds are theirs multiplied by `exp(-(325 - 250) kJ/mol / (R 1600 K))`, 3.56e-3, so the viscosity of a mantle at 1,600 K in the middle of the
ranges holds (1.7e20 Pa s) and the corners of the viscosity at 1,600 K go from 9e16 and 3e23 Pa s to 3e16 and 1e24 Pa s. Source tier T4, tweak level L3 (registry GEN.725, row IHA1)."""

INTERIOR_ACTIVATION_RANGE_KJ_MOL = (260.0, 390.0)
"""tuple: The viscosity activation energy `Ea`, kJ/mol, drawn uniformly (GEN.1185): the diffusion creep a tide and a convecting mantle follow, 261 to 375 kJ/mol dry and 335 to 387 wet
(Korenaga and Karato 2008 Table 3; READ). Guimond et al. 2022 take 200 (wet enstatite) to 300 (dry olivine), below every laboratory value for olivine. Source tier T1, tweak level L1 (row IHA2)."""

INTERIOR_LIQUIDUS_COEFFICIENT_K_GPA = 310.317
"""float: `a_l` of the mantle's liquidus `T_l(p) = a_l (p + p_l)^c_l`, K per GPa^c_l: chondritic material, Manske, Ruedas, Plesa et al. 2025
(arXiv:2506.18718, Table 1 and Eq. 7, a Simon and Glatzel form fitted to experiments; READ). 1,874 K at the surface and 4,717 K at the
Earth's core-mantle boundary (135 GPa)."""

INTERIOR_LIQUIDUS_OFFSET_GPA = 29.102
"""float: `p_l` of the liquidus, GPa (the same table)."""

INTERIOR_LIQUIDUS_EXPONENT = 0.5337
"""float: `c_l` of the liquidus (the same table)."""

ROGUE_TIDAL_FLUX_RANGE_W_M2 = (1e-4, 1e-2)
"""tuple: Tidal heating, W/m^2, drawn log-uniformly for a rocky rogue
that kept a moon (a default: Earth's own tidal dissipation is ~0.007
W/m^2; a captured moon on an eccentric orbit can add more)."""

ROGUE_GIANT_COOLING_TIME_EXPONENT = -1.3
"""float: Giant planets and brown dwarfs fade as L ~ t^-1.3 (Burrows &
Liebert 1993; Burrows et al. 2001, Rev. Mod. Phys. 73:719)."""

ROGUE_JUPITER_INTERNAL_FLUX_W_M2 = 5.4
"""float: Jupiter's internal heat flux at 4.5 Gyr (T_int ~ 99 K), the
anchor for a rogue giant's own heat."""

ROGUE_NEPTUNE_INTERNAL_TEMPERATURE_K = 53.0
"""float: Neptune's internal temperature (Pearl & Conrath 1991). With
Jupiter it fixes how a sub-Jupiter giant's heat scales with mass
(L ~ M^1.65)."""

# Brown dwarfs as stars (GEN.252, GEN.549, GEN.948, docs/design/systems/brown-dwarf-state.md): the cooling track (the Baraffe table,
# `datasets/baraffe_isochrones.py`), the mass function and the map from temperature to spectral type.
BROWN_DWARF_TRACK_END_GY = 10.0
"""float: The age the Baraffe table is cut at (10 Gy): a mass the table still holds at its last age has settled on its plateau."""

BROWN_DWARF_TEMPERATURE_AGE_EXPONENT = -0.32
"""float: Burrows et al. 2001 (Rev. Mod. Phys. 73:719, eq. 2, read, arXiv astro-ph/0103383): the effective temperature of
an old brown dwarf falls as `t^-0.32`. It continues a mass past the last age the Baraffe table holds it at (the models stop where an
object has cooled below their range, about 1,500 K)."""

BROWN_DWARF_MASS_FUNCTION_SLOPE = 0.6
"""float: The substellar mass function `dN/dM ~ M^-0.6` (Kirkpatrick et al. 2021, ApJS 253:7, arXiv 2011.11616, read:
alpha = 0.6 +- 0.1 from the 20 pc census of 525 L, T and Y dwarfs). It is drawn between
`ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER`; the generator's own stellar mass function is not continued below
`IMF_BREAKS_SOL[0]`, since continuing it would give 0.7 brown dwarfs per star against the census's 0.2 (the
free-floating rate `PHENOMENON_DENSITY_PC3["brown-dwarf"]`)."""

ROGUE_GIANT_PHOTOSPHERE_PRESSURE_BAR = 0.3
"""float: Where a Jupiter-gravity giant's photosphere (optical depth 2/3)
sits. It fixes one Rosseland mean opacity (kappa = (2/3) g / P) for every
giant, so a heavier one's photosphere sits deeper, and puts a 1 Mjup,
4.5 Gyr rogue at ~140 K at 1 bar (Boss's research: 130-170 K)."""

ROGUE_ENVELOPE_PHOTOSPHERE_PRESSURE_BAR = 0.1
"""float: A rocky rogue's hydrogen envelope radiates from about 0.1 bar;
below that, collision-induced absorption makes it opaque (Stevenson 1999;
Pierrehumbert 2010)."""

ROGUE_HYDROGEN_ENVELOPE_CHANCE = {"terrestrial": 0.1, "sub-neptune": 0.5}
"""dict: Chance a rocky rogue kept a primordial hydrogen envelope, by
mass bin -- one ejected before its star could strip it (Stevenson
1999). Defaults: rare for Earth-mass rogues, common for sub-Neptunes."""

ROGUE_HYDROGEN_ENVELOPE_PRESSURE_BAR = {"terrestrial": (10.0, 1000.0), "sub-neptune": (100.0, 10000.0)}
"""dict: The envelope's base pressure, bar, drawn log-uniformly (Stevenson:
10^2-10^3 bar; deeper for sub-Neptunes, a default)."""

ROGUE_MIN_ESCAPE_PARAMETER = 15.0
"""float: The Jeans escape parameter (G M m / k T R) above which a gas
stays bound for billions of years; below it the envelope is lost."""

ROGUE_WATER_RICH_CHANCE = {"terrestrial": 0.35, "sub-neptune": 0.8}
"""dict: Chance a rocky rogue is water-rich, by mass bin (defaults: ejected
planets often formed beyond a snow line)."""

ROGUE_WATER_MASS_FRACTION_RANGE = (1e-4, 0.1)
"""tuple: A water-rich rogue's water share of its mass, drawn log-uniformly
(Earth's oceans are 2.3e-4)."""

ROGUE_FROZEN_ATMOSPHERE_MIN_MASS_EARTH = 0.3
"""float: A rocky rogue at least this heavy once held an air (like Earth or
Venus, unlike the Moon) that has frozen onto its surface as nitrogen and
carbon dioxide frost (a default); a lighter, dry one is bare rock."""

ROGUE_ACTIVE_HEAT_FLUX_W_M2 = 0.04
"""float: A rocky rogue whose heat flow reaches this is still geologically
active (`has_internal_heat`): Earth flows 0.087 W/m^2, Mars about 0.02
(a default between them). A giant always is."""

# Interstellar object hyperbolic excess speed range, in km/s -- both
# confirmed interstellar visitors, 1I/'Oumuamua and 2I/Borisov, were
# measured on unbound hyperbolic trajectories with excess speeds in
# this range relative to the Sun (Jewitt & Seligman 2023, Annual Review
# of Astronomy and Astrophysics 61:197, review of both objects).
INTERSTELLAR_OBJECT_SPEED_KMS_RANGE = (10.0, 90.0)

# Interstellar comet nucleus diameter range, in km -- 'Oumuamua was
# estimated at ~0.1-0.4 km and Borisov at ~0.4-1 km (Jewitt & Seligman
# 2023), widened modestly at both ends for generation variety.
INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM = (0.05, 5.0)

# Chance a generated interstellar comet is currently active (showing a
# coma/tail from sublimating ices, as Borisov did) rather than an inert,
# 'Oumuamua-like body with no detected outgassing.
INTERSTELLAR_COMET_ACTIVE_CHANCE = 0.5

# A comet nucleus's makeup is drawn by `physics/small_body_composition.draw_comet` (GEN.442): the refractory-to-ice ratio, the
# organic share, the porosity and the volatiles, with the density they make.

# --- Star-Bound Comets (cometData.Comet) ---
#
# Unlike InterstellarComet (unbound, hyperbolic, encountered once, always
# standalone), a Comet here is bound to a star's own system, propagated
# via real two-body Kepler/Barker orbital mechanics (see keplerMotion.py)
# rather than a fixed hyperbolic excess speed. See
# docs/design/systems/comet-orbital-realism.md for the full research/design
# writeup this implements.

# Bound short-period/long-period comet nuclei run larger, on average,
# than the two confirmed interstellar visitors this generator also models
# (INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM, ~0.1-1 km) -- e.g. Halley
# is ~15x8 km (Keller et al. 1986, ESA Giotto results) and Jupiter-family
# comet nuclei are typically observed in the 1-10 km range (Lamy et al.
# 2004, "Sizes, Shapes, Albedos, and Colors of Cometary Nuclei"). A
# separate, wider range reflects that real difference rather than reusing
# the interstellar range as-is.
BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM = (0.5, 20.0)

# Perihelion distance range, in AU, for a star-bound comet -- wide enough
# to cover both a sungrazer-like close pass and an outer-system comet that
# never gets very active. COMET_ACTIVITY_PERIHELION_THRESHOLD_AU (below)
# is what actually drives whether a given perihelion produces visible
# activity, not this range's own bounds.
COMET_PERIHELION_DISTANCE_RANGE_AU = (0.05, 5.0)

# Real ices (dominated by water ice) begin sublimating noticeably inside
# roughly 2.5-3 AU of a Sun-like star -- the classic "ices turn on" comet
# activity threshold (e.g. Meech & Svoren 2004, "Physical and chemical
# evolution of cometary nuclei", in "Comets II"). Used by
# `Comet._roll_activity` to scale activity chance by perihelion distance
# rather than a flat roll (contrast INTERSTELLAR_COMET_ACTIVE_CHANCE,
# where an interstellar comet's arbitrary/often-irrelevant perihelion
# makes a flat roll the more honest choice).
COMET_ACTIVITY_PERIHELION_THRESHOLD_AU = 3.0

# Activity-chance floor/ceiling `Comet._roll_activity` interpolates
# between across COMET_ACTIVITY_PERIHELION_THRESHOLD_AU: near-certain
# activity for a sungrazer-close perihelion, down to a small residual
# chance (a comet can still show faint activity, or have shed enough
# volatiles over past passages, well beyond the nominal sublimation
# threshold) at or beyond it.
COMET_ACTIVITY_MAX_CHANCE = 0.9
COMET_ACTIVITY_MIN_CHANCE = 0.05

COMET_PERIOD_CLASSES = {
    # Real dynamical comet families (see docs/design/systems/comet-orbital-realism.md's
    # taxonomy table) -- each entry's period_range_years/eccentricity_range/
    # inclination_max_deg are drawn from, in that order: Jupiter-family
    # comets (repeatedly perturbed by Jupiter into short, low-inclination,
    # moderately eccentric orbits -- e.g. 2P/Encke, 67P/Churyumov-
    # Gerasimenko), Halley-type comets (intermediate period, can retain a
    # high or even retrograde inclination from their less-processed
    # Oort-Cloud-adjacent origin -- 1P/Halley itself is inclined 162 deg,
    # i.e. retrograde), and long-period comets (e very close to 1,
    # isotropic inclination -- true Oort Cloud origin, only weakly
    # dynamically processed). `weight` is this family's relative share of
    # generated elliptical comets, not a physical population estimate.
    "jupiter_family": {
        "period_range_years": (3.3, 20.0),
        "eccentricity_range": (0.2, 0.75),
        "inclination_max_deg": 30.0,
        "weight": 0.5,
    },
    "halley_type": {
        "period_range_years": (20.0, 200.0),
        "eccentricity_range": (0.5, 0.97),
        "inclination_max_deg": 180.0,
        "weight": 0.3,
    },
    "long_period": {
        "period_range_years": (200.0, 4_000_000.0),
        "eccentricity_range": (0.9, 0.999),
        "inclination_max_deg": 180.0,
        "weight": 0.2,
    },
}
"""
dict: Elliptical star-bound comet subtype table, keyed by `period_class`
(see `cometData.Comet`'s own docstring) -- flavor/plausibility metadata
only (`keplerMotion.py`'s propagation is identical for every subtype; only
the *sampled ranges* differ per family), not a separate physics model.
"""

# A "parabolic" Comet's stored eccentricity is drawn from just under 1.0
# rather than fixed at exactly 1.0 -- realistic near-parabolic long-period/
# single-apparition comets are never mathematically exact parabolas -- but
# `keplerMotion.py` always propagates orbit_type == "parabolic" via the
# exact Barker's-equation solution regardless of the specific value drawn
# here (the numerical difference between e=1.0 and e=0.998 over a single
# apparition is negligible, and the elliptical Kepler-equation solver
# converges poorly as e -> 1 anyway -- see keplerMotion.solve_eccentric_anomaly).
PARABOLIC_COMET_ECCENTRICITY_RANGE = (0.995, 1.0)

# Chance a star-bound comet generates as parabolic (single-apparition,
# escapes after this one perihelion passage) rather than elliptical
# (periodic, returns every orbit) -- see docs/design/systems/comet-orbital-realism.md's
# taxonomy table. Most bound comets observed are periodic; long-period/
# parabolic "great comets" are the rarer, more dramatic case.
COMET_PARABOLIC_CHANCE = 0.3

# Inclination is drawn independently of COMET_PERIOD_CLASSES'
# inclination_max_deg for a parabolic comet (which has no period_class of
# its own) -- isotropic, like a true long-period/Oort-Cloud-origin comet,
# since a single-apparition object hasn't been dynamically flattened into
# a low-inclination orbit the way a repeatedly-perturbed Jupiter-family
# comet has.
PARABOLIC_COMET_INCLINATION_MAX_DEG = 180.0

# Chance a given star (per `StarSystem._generate_comets`, one roll per
# star -- the primary, and again for a wide binary's independently-rolled
# secondary) has any native comets at all, when
# `SystemConfig.COMETS` is left at its default `None` (random chance) --
# most real stars aren't known to host an observed comet population, so
# this stays well under 1.0 rather than guaranteeing every system gets
# one.
SYSTEM_COMET_CHANCE = 0.3

# --- Prevalence (generation/prevalence.py, GEN.52) ---
# The share of systems that come out with each feature that no single draw
# decides, measured over 20,000 systems generated with default options
# (four seeds of 5,000; 2026-10-07, remeasured 2026-10-17 after GEN.447 to
# GEN.449 changed the companions and 2026-10-18 after GEN.404 moved a system with
# no place to the galaxy's half-mass radius; habitable_world alone remeasured
# 2026-10-23 after GEN.146 widened the zone to the Kopparapu edges, the same
# 20,000 systems: 0.4216 where it was 0.253; all four remeasured 2026-10-23 over
# 20,000 systems, four seeds of 5,000, after GEN.631 drew a system's stars as a
# tree: its pairs are closer, so fewer have room for planets; habitable_world alone
# remeasured 2026-10-29 over 6,000 systems, four seeds of 1,500, after GEN.855 let
# the air follow its life: 0.354 where it was 0.402; and 2026-10-29 over 3,000 systems, two seeds of 1,500, after GEN.646 derived the class
# from the finished body: 0.722, most of it classes Q, P, V and K, which the feature counts as able to carry life while few are habitable
# by the score; all four remeasured 2026-10-30 over 2,000 systems, seed 11, after trains 3 and 4: habitable_world 0.598, asteroid_belt 0.039 where
# it was 0.549 (the belt mass cut keeps only the belts that clear it), large_star 0.050, planets 0.809). A sector or galaxy run's
# `--prevalence feature=+50` forces the feature on in enough systems to
# make that share 1.5 times as large. Remeasure after a change to
# generation moves one of them: tests/test_prevalence.py checks they still
# hold within a few points.
PREVALENCE_BASE_SHARES = {
    "habitable_world": 0.598,   # a planet or moon of a class that can carry life (see above)
    "asteroid_belt": 0.039,     # an asteroid belt
    "large_star": 0.050,        # a primary born at LARGE_STAR_MIN_MASS_SOL or more
    "planets": 0.809,           # at least one planet or belt
}

# The usual share of each feature one draw decides (ADM.16: the Generate
# page starts each prevalence field at its usual share), measured over
# 10,000 systems with default options (four seeds of 2,500; 2026-10-08;
# comets, binary_system and wide_binary remeasured 2026-10-17 over 20,000
# systems, four seeds of 5,000: GEN.448 puts 93% of the binaries in the wide
# class, the close class being the 0.05 to 0.25 AU range only; binary_system and
# wide_binary remeasured 2026-10-23 over 20,000 systems, four seeds of 5,000,
# after GEN.631 drew a system's stars as a tree: 0.239 have a pair, 0.810 of the
# pairs have an axis over 0.25 AU).
# Each is a share of the things in its comment, not always of systems.
# Remeasure with PREVALENCE_BASE_SHARES.
PREVALENCE_DRAW_SHARES = {
    "comets": 0.321,            # systems with comets (each star of a wide pair draws)
    "binary_system": 0.315,     # systems with a pair (the two stars of a system of more than two count)
    "wide_binary": 0.810,       # pairs that are wide (S-type: an axis of 0.25 AU or more)
    "moons": 0.429,             # planets with at least one moon
    "max_planets": 0.037,       # stars that get the most orbits they can hold
    "multi_star": 0.0044,       # systems with five or more stars (see below)
}
# `multi_star` is not measured from generated systems: `StarSystem` builds only a
# single star, a close pair or a wide pair until the member-tree host lands
# (GEN.630), so no run yet has five stars. It is the `star-count` table's
# five-or-more columns (`generation.multistar.star_count_weights`) averaged over the
# placed star's primary mass function (`imf.PRIMARY`, GEN.928; 400 log bins
# from 0.08 to 150 solar masses), before any age or remnant, 0.442% of systems.
# Remeasure it from generated systems once GEN.630 merges.

# The most stars a system may hold (`SystemConfig.MAX_STARS`, GEN.267; the
# member-tree draw's cap, multi-star-systems.md 11 item 11). 12 is Boss's default
# of 2026-10-15 6:08 PM EDT; the ceiling a user may raise the cap to is a stand-in:
# a chain of axis ratio 3 fits a dozen levels (section 3 item 1), so a few dozen stars
# is as many as branching could hold in the 1e6 span of axes, and a stored tree of
# more would be unreadable on a system map. Gap plan registry: "stars ceiling".
MAX_STARS_DEFAULT = 12
MAX_STARS_CEILING = 64

# How many comets a star that DOES get any (per SYSTEM_COMET_CHANCE, or
# `SystemConfig.COMETS = True` forcing at least this many) actually
# generates -- a small handful, not a full population (this generator
# only models a system's few most notable comets, the same way
# ASTEROID_BELT models zero-or-one belts, not an exhaustive minor-body
# catalog).
SYSTEM_COMET_COUNT_RANGE = (1, 3)

# A star of a wide binary keeps every comet inside its stable limit
# (a_crit_au), or the comet's orbit would cross the companion's (GEN.228).
# Each comet is re-rolled this many times to find an orbit that fits, then
# dropped: close wide pairs leave little room, so most comets there don't.
COMET_CEILING_ATTEMPTS = 25

# --- Asteroid Fields (asteroidFieldData.AsteroidField) ---

# Radius range for a standalone asteroid field drifting in open
# interstellar space (as opposed to AsteroidBelt, which orbits a star),
# in light-years. Lower bound is well above our own Kuiper Belt's real
# scale (~30-50 AU, ~0.0005-0.0008 ly) -- a field with no central star to
# hold it together can plausibly be spread far wider by galactic tidal
# shear over billions of years -- while the upper bound stays below a
# planetary nebula's own radius range (NEBULA_CLASSES "H"-"L", 0.1-3
# ly) so the two remain visually/narratively distinct phenomena.
ASTEROID_FIELD_RADIUS_RANGE_LY = (0.001, 1.0)

# Asteroid field classes (GEN.14, Boss 2026-09-30: "a digit and part
# of the class"). The letter comes from composition and density, following
# asteroid taxonomy (C carbonaceous, S stony, M metallic, D/P icy primitive,
# V basaltic; Bus & Binzel 2002, Icarus 158:146; DeMeo et al. 2009, Icarus
# 202:160); the digit is the field's size, floor(log10(radius in AU)), so a
# 0.001-1 ly field is 1-4. A full class reads like "C3". V-Z are reserved.
#
# Each composition family: `letters` keyed by density (`sparse`, `typical`,
# `dense`), `description`, and `frequency` (relative weight; what its bodies are
# made of is `small_body_composition.FIELD_FAMILIES`).
ASTEROID_FIELD_COMPOSITIONS = {
    "carbonaceous": {
        "letters": {"sparse": "A", "typical": "B", "dense": "C"},
        "description": "dark carbon-rich bodies with hydrated clays",
        "frequency": 30.0,
    },
    "stony": {
        "letters": {"sparse": "D", "typical": "E", "dense": "F"},
        "description": "silicate rock with flecks of metal",
        "frequency": 20.0,
    },
    "metallic": {
        "letters": {"sparse": "G", "typical": "H", "dense": "J"},
        "description": "iron-nickel cores of shattered bodies",
        "frequency": 8.0,
    },
    "icy": {
        "letters": {"sparse": "K", "typical": "L", "dense": "M"},
        "description": "volatile-rich primitive bodies, ices bound in dark crusts",
        "frequency": 20.0,
    },
    "basaltic": {
        "letters": {"sparse": "N", "typical": "N", "dense": "P"},
        "description": "basalt from a differentiated body's crust",
        "frequency": 4.0,
    },
    "mixed": {
        "letters": {"sparse": "Q", "typical": "R", "dense": "S"},
        "description": "a mix of rock, metal and carbonaceous bodies",
        "frequency": 10.0,
    },
    "dust": {
        "letters": {"sparse": "T", "typical": "T", "dense": "T"},
        "description": "mostly dust and gravel, with few large bodies",
        "frequency": 5.0,
    },
    "collisional": {
        "letters": {"sparse": "U", "typical": "U", "dense": "U"},
        "description": "fragments of one parent body broken apart in a collision",
        "frequency": 3.0,
    },
}
"""dict: See the comment above. What a family is made of is its mixture of asteroid taxa
(`physics/small_body_composition.FIELD_FAMILIES`, GEN.477)."""

# --- Collisions (GEN.345, docs/design/orbits/collisions-and-mergers.md section 5.1) ---

ALWAYS_DESTROY_TERRESTRIAL = False
"""bool: Boss's literal rule for two rocky or icy bodies (2026-10-03): any contact destroys both and leaves one field of
the total mass. Off, the outcome table of section 5.1 decides (merge, hit-and-run, erosion or disruption by speed,
angle and mass ratio)."""

COLLISION_MERGE_SPEED_RATIO = 1.2
"""float: Two terrestrial bodies that meet slower than this many mutual escape speeds merge (a design value, section 5.1)."""

COLLISION_HIT_AND_RUN_SPEED_RATIO = 1.5
"""float: A grazing pass faster than this many escape speeds is a hit-and-run, slower a graze-and-merge (Genda, Kokubo
and Ida 2012; a design value)."""

COLLISION_EROSION_X = 0.1
"""float: Below this share of the catastrophic energy a head-on pass of two terrestrial bodies is a clean merger; above it
(and under 1) a partial erosion."""

COLLISION_FIELD_MIN_SHARE = 0.01
"""float: A hit-and-run's eroded mass becomes an asteroid field only when it is at least this share of the pair's mass,
else it is only logged."""

GIANT_MERGER_LOSS_RANGE = (0.01, 0.5)
"""tuple: The share of a merged giant's mass lost as gas, `Q_R / (2 Q*_RD)` held between these (SPH finds at least 97 percent
kept for quick Jupiter-Jupiter mergers, section 4.2)."""

COLLISION_BORDERLINE_SHARE = 0.03
"""float: A merger within this share of the hydrogen-burning limit is flagged borderline in the report (section 4.3)."""

ASTEROID_FIELD_SIZE_DIGIT_RANGE = (1, 4)
"""tuple: The size digit's clamp -- floor(log10(radius in AU)) over
`ASTEROID_FIELD_RADIUS_RANGE_LY` (63 AU to 63,241 AU)."""

# --- Quasars (quasarData.Quasar) ---
#
# A quasar is not a free-floating object: it is a galaxy's own central
# supermassive black hole caught in a phase of near-Eddington accretion,
# its accretion disk outshining every star in the host galaxy combined.
# There is exactly one such nucleus per galaxy, at its dynamical center,
# so the generator only ever places a quasar there (see
# `generate.add_galactic_nucleus`), never scattered
# through ordinary sectors the way `PHENOMENON_DENSITY_PC3`'s
# types are.

QUASAR_ACTIVE_NUCLEUS_CHANCE = 0.1
"""
float: The chance a generated galaxy's nucleus is active -- i.e. that the
one sector that hosts the galactic center gets a quasar at all. Real
quasar activity is rarer still in today's universe (quasars peaked around
redshift ~2, some 10 billion years ago -- Richards et al. 2006, AJ
131:2766 -- and the Milky Way's own Sagittarius A* is quiescent), so this
is a deliberately generous "rare but possible" setting rather than a
derived rate.
"""

QUASAR_EDDINGTON_RATIO_RANGE = (0.1, 1.0)
"""
tuple: Bolometric luminosity as a fraction of the Eddington limit, sampled
log-uniformly -- luminous quasars cluster between ~0.1 and 1 with a
median near 0.25 (Kollmeier et al. 2006, ApJ 648:128).
"""

EDDINGTON_LUMINOSITY_W_PER_SOLAR_MASS = 1.26e31
"""
float: The Eddington luminosity per solar mass of accretor, in watts
(`4 pi G M m_p c / sigma_T` for ionized hydrogen, ~1.26e38 erg/s per
solar mass).
"""

QUASAR_RADIATIVE_EFFICIENCY = 0.1
"""
float: Fraction of accreted rest-mass energy radiated away, `L = eta *
Mdot * c^2` -- the standard ~10% thin-disk value (Soltan 1982, MNRAS
200:115; Yu & Tremaine 2002, MNRAS 335:965).
"""

QUASAR_RADIO_LOUD_CHANCE = 0.1
"""
float: The chance a quasar is radio-loud (launches relativistic jets) --
roughly 10% of optically selected quasars are (Ivezic et al. 2002, AJ
124:2364).
"""

QUASAR_JET_LENGTH_RANGE_LY = (3e4, 3e6)
"""
tuple: Jet/radio-lobe extent range for a radio-loud quasar, in
light-years, sampled log-uniformly -- from ~10 kpc to the ~1 Mpc of giant
radio sources. Far larger than the host galaxy, which is why it is
described rather than drawn at map scale.
"""

QUASAR_ACTIVE_AGE_RANGE_YEARS = (1e6, 1e8)
"""
tuple: How long this episode of quasar activity has been running, in
years, sampled log-uniformly -- quasar lifetimes are constrained to
~1e6-1e8 years (Martini 2004, "QSO Lifetimes", Carnegie Obs. Astrophys.
Ser. 1).
"""

QUASAR_BOLOMETRIC_CORRECTION_5100 = 9.26
"""
float: Bolometric luminosity divided by the 5100-angstrom monochromatic
luminosity (`L_bol / lambda L_5100`) for a typical quasar spectrum
(Richards et al. 2006, ApJS 166:470).
"""

QUASAR_BLR_RADIUS_LIGHT_DAYS_AT_1E44 = 33.65
QUASAR_BLR_RADIUS_LUMINOSITY_SLOPE = 0.533
"""
Broad-line-region radius-luminosity relation from reverberation mapping,
`R_BLR = 33.65 light-days * (lambda L_5100 / 1e44 erg/s)^0.533` (Bentz et
al. 2013, ApJ 767:149).
"""

MILKY_WAY_STELLAR_LUMINOSITY_W = 2.5e10 * 3.828e26
"""
float: The combined starlight of a Milky-Way-sized galaxy, in watts
(~2.5e10 solar luminosities -- Licquia, Newman & Brinchmann 2015, ApJ
809:96), which a quasar's description compares itself against.
"""

# --- Phenomenon Generation Mode (phenomenonGen.py) ---

PHENOMENON_TYPE_CHOICES = _phenomena_registry.phenomenon_types()
"""
tuple: The valid `--type` values `phenomenonGen.py` accepts, read from the
phenomena registry (GEN.250). Omitting
`--type` picks uniformly at random among `RANDOM_PHENOMENON_TYPE_CHOICES`.
"""

RANDOM_PHENOMENON_TYPE_CHOICES = tuple(
    choice for choice in PHENOMENON_TYPE_CHOICES if choice != "quasar"
)
"""
tuple: `PHENOMENON_TYPE_CHOICES` minus `"quasar"`, which only exists at a
galaxy's center (see `QUASAR_ACTIVE_NUCLEUS_CHANCE`) and so is only ever
generated when asked for by name.
"""

# --- Sector-Level Exotic Phenomena (generation/run_sector.py, generation/run_galaxy.py) ---
#
# Every generated sector seeds a population of exotic phenomena, sampled
# independently per type via a Poisson draw (spaceSector.
# _sample_poisson_count, generate.generate_sector_phenomena).
#
# Rates come from Boss's research of 2026-09-30 (docs/design/systems/
# interstellar-object-rates.md): each type has a local number density
# n_i0 (objects per pc^3) at the solar neighborhood's stellar density
# n_*0 = REFERENCE_STELLAR_DENSITY_PC3, and scales with the local stellar
# density, n_i = n_i0 * n_* / n_*0. A sector's star count already follows
# n_* * V, so that is a fixed rate per STAR, n_i0 / n_*0, times the
# sector's own star count (`phenomenon_rate_per_star`).

REFERENCE_STELLAR_DENSITY_PC3 = 0.14
"""float: Stars per pc^3 in the solar neighborhood (0.10-0.14 in the
literature), the density every PHENOMENON_DENSITY_PC3 figure is quoted at."""

PHENOMENON_DENSITY_PC3 = {
    # Sum of ROGUE_PLANET_MASS_BINS' per-star rates (5.8) at n_*0.
    # Terrestrial alone: about 0.78 pc^-3 (0.5-1.4); Jupiter-mass about
    # 0.004 (under Mroz et al. 2017's limit of 0.035).
    "rogue-planet": sum(rate for _lo, _hi, rate in ROGUE_PLANET_MASS_BINS.values()) * REFERENCE_STELLAR_DENSITY_PC3,
    # 0.025-0.035 (Kirkpatrick et al. 2021, ApJS 253:7), ~1 per 4.7 stars.
    "brown-dwarf": 0.03,
    # Stars moving > 30 km/s, 1-2% of stars (Tauris 2015, MNRAS 448:L6).
    # A flag on an ordinary generated system, not a phenomenon row.
    "runaway-star": 2.1e-3,
    # Stars moving > 500 km/s, ejected by the central black hole, at
    # HVS_REFERENCE_SCALE_LENGTHS; scales as r_GC^-2 (Brown 2015, ARA&A
    # 53:15). The research's 5e-9 implies ~4e5 out to 100 kpc against a
    # stated 1e3-1e4 total; 1e-10 matches the total. Also a system flag.
    "hypervelocity-star": 1e-10,
    # ~1e9 isolated neutron stars and ~1e8 black holes galaxy-wide
    # (Sartore et al. 2010, A&A 510:A23; Olejak et al. 2020, A&A 638:A94;
    # Sahu et al. 2022, ApJ 933:83). Tuned to the star model's own remnant
    # shares (docs: galaxy-studies star-fix spec, 2026-09-30): a Kroupa IMF
    # with a 10 Gy thin disk leaves ~0.5% of stars as neutron stars and
    # ~0.1% as black holes, i.e. 0.005 and 0.001 per star at n*0 = 0.14.
    # GEN.100 retune (Boss, 2026-10-09): the central observed values, 1e9
    # neutron stars (0.4% of ~2.5e11 stars) and 1e8 black holes (0.05%).
    "neutron-star": 5.6e-4,
    # The stellar-mass black holes are not a rate of their own: they follow from the stars' metallicity
    # (`galaxy.remnant_density`, GEN.563). The old 7e-5 (1e8 in the Galaxy, 0.05% of ~2.5e11 stars) is a test target.
    # Giant molecular clouds: 1e-6 to 1e-5 (Kennicutt & Evans 2012, ARA&A
    # 50:531): cloud centers, generated as dark-family nebulae (classes
    # M-Q). GMC_ARM_FILLING_FACTOR is the matching fraction of arm volume
    # inside a cloud, for reference.
    "molecular-cloud": 5e-6,
    # ~20,000 planetary nebulae galaxy-wide (Frew & Parker 2010, PASA
    # 27:129).
    # 1.4e-8 is 1e-7 per star (20,000 over ~2e11 stars).
    "planetary-nebula": 1.4e-8,
    # About 3,000 supernova remnants galaxy-wide (GEN.203): the supernova
    # rate gives 400 to 3,000 that are still visible (Green 2019, 310
    # catalogued; Ranasinghe & Leahy 2022), 1.5e-9 per pc^3 is the top of
    # that range. It was 1e-8, 20,000 remnants (n_* * rho_gas, Draine 2011),
    # which counted every faint one.
    "supernova-remnant": 1.5e-9,
    # Notable interstellar comets kept as rows, a design rate (today's
    # 0.05 per star). The real population of comets and planetesimals,
    # INTERSTELLAR_DEBRIS_DENSITY_PC3, is a computed sector figure.
    "comet": 0.007,
    # Free asteroid fields disperse in 1e6-1e7 years (Raymond et al.
    # 2020, ApJL 894:L22), so none are generated; the type stays for
    # hand-made fields and facilities.
    "asteroid-field": 0.0,
}
"""dict: Local number density per pc^3 at REFERENCE_STELLAR_DENSITY_PC3,
the research value for each kind (see the design doc's table)."""

WEIRDNESS_DIAL = 1.0
"""float: The world-level weirdness dial (GEN.1036, docs/design/systems/weird-worlds.md 4.4): it multiplies the rate of every oddity row
(`physics/oddities.py`); 0 turns every oddity off, 1 gives the literature rates. Default 1.0 (design; Boss asked for a dial,
4:15 PM EDT 2026-10-30)."""

WEIRDNESS_RATE_SCALE = {}
"""dict: A multiplier per oddity id (`{"W10": 0.5}`, default 1.0 for any id not named), applied with `WEIRDNESS_DIAL` to the rate of that row."""

PHENOMENON_RATE_SCALE = {kind: 1.0 for kind in _phenomena_registry.rated_kinds()}
"""dict: A multiplier per kind (default 1.0) to dial a rate up or down
without touching its research value -- e.g. 0.1 for "rogue-planet" gives
about 6 rogues per local sector instead of about 58 (Boss chose the full
rate on 2026-09-30)."""

GMC_ARM_FILLING_FACTOR = 0.015
"""float: Share of spiral-arm volume inside a giant molecular cloud
(0.01-0.02), for GEN.10's cloud placement. Nothing places clouds from it;
the nebula audit test (TEST.136) holds the field's arm-crest filling to its band."""

NEBULA_FIELD_KEEP = {"M": 0.095, "N": 0.5, "P": 1.0, "Q": 1.0}
"""dict: The share of each field class's drawn clouds that the galaxy keeps (GEN.152). The field draws 700,000
clouds at `PHENOMENON_DENSITY_PC3["molecular-cloud"]`, class M ten times the 8,107 to 9,710 catalogued giant
molecular clouds and the volume inside the clouds' bounding spheres 2 to 24% against 0.5 to 2% observed. Each
cloud keeps or loses by a hash of its cell and draw number (`nebula_field.cloud_kept`), so the draw, and every
survivor's position and object ID, are as they were. Class N (dark clouds) has the lowest confidence; P and Q are
in range. Study: docs/design/systems/nebula-density-vs-reality.md."""

GMC_GAS_DENSITY_EXPONENT = 1.4
"""float: Molecular clouds scale with gas density to this power
(Schmidt-Kennicutt), for GEN.10."""

NEBULA_FIELD_CELL_PC = 50.0
"""float: Edge of the galaxy's molecular cloud field cells, parsecs
(GEN.47, `nebulaField`): each cell draws its own clouds from its own
seed. About a cloud's own size, so a sector reads a few dozen cells."""

NEBULA_FIELD_ARM_AMPLITUDE = 0.6
"""float: The molecular gas's arm contrast, `1 + A cos(arm phase)` like
`STELLAR_POPULATION_ARM_AMPLITUDE`: 0.6 puts four times as much gas on an
arm's crest as midway between arms (CO arm to interarm contrasts run
about 2 to 5; Koda et al. 2009, ApJ 700:L132)."""

NEBULA_FIELD_MAX_GAS_FACTOR = 3.0
"""float: The most molecular gas the cloud field assumes anywhere,
relative to the solar circle's average (`nebulaField.gas_factor`). The
young-star tracer climbs twentyfold toward the inner disk; uncapped, the
inner galaxy would be solid cloud. At 3, an arm crest or the inner disk
has about a quarter of its volume inside a dark cloud."""

SUPERBUBBLE_CELL_PC = 250.0
"""float: Edge of the superbubble field's columns, parsecs (GEN.254, `superbubble_field`). The design note said 50,
the molecular cloud field's cell, but a bubble reaches up to `SUPERBUBBLE_RADIUS_RANGE_PC[1]` and a point would
read 20 cells in each direction (8,000 columns); at 250 it reads 5."""

SUPERBUBBLE_RADIUS_RANGE_PC = (50.0, 500.0)
"""tuple: The radius of a superbubble, parsecs, drawn log-uniform: 100 to 1,000 pc across, the Local Bubble's 100 to
200 pc (Zucker et al. 2022, Nature 601:334) at the small end, the large superbubbles of other galaxies (Chu 2008,
IAU 255) at the other."""

SUPERBUBBLE_FILL_FRACTION = 0.12
"""float: The share of the plane at the solar circle inside a superbubble. Superbubbles fill about 10 to 15% of the
disk (Heiles 1979; McKee & Ostriker 1977 put the hot phase at 10 to 30%); the arms and inner disk hold more, as the
molecular gas does (`nebulaField.gas_factor`)."""

SUPERBUBBLE_CENTRE_HEIGHT_SIGMA_PC = 60.0
"""float: The spread of a superbubble's centre above and below the plane, parsecs, about the young stars' layer: they
are blown by the supernovae of OB associations, which are born in it."""

SUPERBUBBLE_AMBIENT_DENSITY_RANGE_CM3 = (0.001, 0.01)
"""tuple: The gas inside a cavity, H atoms per cm3, log-uniform; the Local Bubble holds about 0.005 (Snowden et al.
1990, Lallement et al. 2003) against the disk's 1. Remnants inside swell to match (GEN.151)."""

SUPERBUBBLE_TEMPERATURE_K = 1.0e6
"""float: The gas inside a cavity: ionised, about a million kelvin (the soft X-ray background)."""

HVS_REFERENCE_SCALE_LENGTHS = 8000.0 / 2600.0
"""float: Galactic radius the "hypervelocity-star" density is quoted at, in the galaxy's thin disk scale lengths (the Milky
Way's 8 kpc over its 2.6 kpc, scaled to whatever galaxy this is; `galaxy.scales`, GEN.404)."""
DISK_HALF_MASS_SCALE_LENGTHS = 1.6783469900170382
"""float: The radius that holds half of an exponential disc's stars, in scale lengths (the root of `1 - (1 + x) exp(-x)
= 1/2`): where a system with no place in the galaxy is taken (`galaxy.scales.default_galactocentric_pc`)."""

INTERSTELLAR_DEBRIS_DENSITY_PC3 = 1e15
"""float: Interstellar objects of about 100 m and up per pc^3 at REFERENCE_STELLAR_DENSITY_PC3: 0.1 to 0.2 per au^3,
which is 7e14 to 1.8e15 per pc^3 (Do, Tucker and Tonry 2018 from 'Oumuamua; Jewitt and Seligman 2023, ARA&A), 1e15 the
middle of the range; smaller bodies are more numerous still (`docs/design/small-body-size-laws.md`). Far too many for
rows, so a sector shows it as a figure (`planetgen.db.query.interstellar_debris_count`)."""

RUNAWAY_STAR_SPEED_RANGE_KMS = (30.0, 200.0)
"""tuple: A runaway star's speed relative to its neighbors, km/s
(log-uniform; Tauris 2015)."""

HYPERVELOCITY_STAR_SPEED_RANGE_KMS = (500.0, 1000.0)
"""tuple: A hypervelocity star's speed, km/s (Brown 2015)."""

HYPERVELOCITY_STARS_PER_GALAXY = (1e3, 1e4)
"""tuple: How many hypervelocity stars a galaxy's scatter places out to `HYPERVELOCITY_COUNT_RADIUS_KPC`, drawn
log-uniformly between the two ends (Boss 2026-10-08: "about 1e3 to 1e4 galaxy-wide"; Brown 2015). Each starts at the
central black hole and moves outward at `HYPERVELOCITY_STAR_SPEED_RANGE_KMS`; their density falls as r^-2, so the
count out to a reach is this number times the reach over that radius (GEN.238)."""

NEUTRON_STAR_KICK_SIGMA_KMS = 265.0
"""float: The one-dimensional dispersion of a neutron star's natal kick, km/s, a Maxwellian in three dimensions (Hobbs et
al. 2005's pulsar birth velocities, 265 km/s, as `docs/design/systems/multistar-and-compact-systems.md` records it) [R, recalled in
the repository's notes and not re-read: Boss, 2026-10-21, 1:17 PM EDT, "assume already included information is valid"]. With
the galaxy's escape speed (`galaxy.natal_kick`) it sets the share of the neutron stars a galaxy ejects (GEN.554)."""

NEUTRON_STAR_KICK_BIRTH_RADII = 32
"""int: How many equal shares of the thin disc's massive stars the neutron stars' birth radius is averaged over
(`galaxy.natal_kick`); the integral over a smooth tail, so more changes the share in the third digit."""

NEUTRON_STAR_UNBOUND_FRACTION_TARGET_RANGE = (0.2, 0.4)
"""tuple: The share of the Milky Way's neutron stars the kicks send out of it, a test target and no input (GEN.554): 20% to
40% (the old per-galaxy draw: Hobbs et al. 2005's kick distribution against the escape speed [S, recalled]; "about 40%" in
one population synthesis, `docs/design/systems/compact-remnant-regions.md`). A galaxy of the Milky Way's masses lands inside it."""

BLACK_HOLE_UNBOUND_FRACTION = 0.02
"""float: The share of a galaxy's black holes that are ejected (GEN.238; about 2%) [S, recalled; to verify]."""

CORE_COLLAPSE_MASS_RANGE_SOL = (SUPERGIANT_MIN_MASS_SOL, IMF_BREAKS_SOL[-1])
"""tuple: The initial masses that end as a neutron star or a black hole, solar masses: 8 (the supergiants' edge) to the
top of the IMF, 150 (GEN.560). 0.713% of the stars born (GEN.927)."""

BLACK_HOLE_FAILED_SHARE_SOLAR_RANGE = (0.0, 0.45)
"""tuple: The share of the massive stars of solar metallicity that fail to explode and leave a black hole: the widest
verified range, 0 (O'Connor and Ott 2011's lowest) to 45% (Sukhbold et al. 2016's highest; Ebinger et al. 2020 give 5 to
21%). A galaxy draws its own, uniform over the range (only the range is known) and cut by the neutron star budget
(`physics.black_hole_formation.share_cap`). GEN.563, docs/design/systems/black-hole-formation-sources.md."""

BLACK_HOLE_FAILED_SHARE_METAL_FREE_RANGE = (0.15, 0.32)
"""tuple: The same share for a metal-free population, 15% (O'Connor and Ott 2011) to 32% (Ebinger et al. 2020's
highest; Sukhbold et al. give solar only). Drawn uniform over the range, cut at the neutron star budget and not below
the galaxy's solar share. GEN.563."""

REFERENCE_BLACK_HOLE_DRAWS = (0.55, 0.5)
"""tuple: The two uniform draws (solar share first, then metal-free) of a galaxy built with no seed, the reference structure that
stands for the Milky Way (`physics.black_hole_formation.reference_failed_shares`). The metal-free draw is the middle of its range.
The solar draw was the middle (0.5) until the halo of the 5e10 solar masses of stars got lighter (GEN.1154: 1.6e12 against 2.0e12
solar masses), which lifted the neutron stars' unbound share from 0.159 to 0.208, cut the black hole budget (`share_cap`) from 0.35
to 0.32 and so the middle share from 0.175 to 0.161: 9.3e7 black holes in a galaxy of 5e10 solar masses of stars, under the 1e8 to 1e9
of the older estimates Olejak et al. 2020 list (their own count is 1.3e8). 0.55 gives a share of 0.177 (0.175 before) and 1.02e8, back
inside that range at its low edge, and keeps the second moment the design note records (the draw that would give Olejak's 1.3e8,
0.70, moves that moment by 13 percent). A stand-in tuned to the edge of one published range (tier T1, level L2, row BR1 of the registry,
GEN.725)."""

BLACK_HOLE_FORMED_MASS_STEPS = 64
BLACK_HOLE_FORMED_AGE_STEPS = 64
"""int: Log-spaced progenitor masses and ages across a window the share of formed remnants is averaged over (numerical,
`physics.black_hole_formation.remnant_formed_share`)."""

BLACK_HOLE_FAILED_MEAN_TABLE_STEPS = 1025
"""int: Points of the table the failed share's mean over a metallicity peak is read from (it is smooth in the peak's mean)."""
BLACK_HOLE_GALAXY_MEAN_RADIAL_STEPS = 72
BLACK_HOLE_GALAXY_MEAN_POLAR_STEPS = 24
BLACK_HOLE_GALAXY_MEAN_AZIMUTH_STEPS = 6
"""int: The log-spaced radii, equal-area polar angles and azimuths of the grid a galaxy's mean black hole rate per star
is averaged over (`galaxy.remnant_density.mean_rate_per_star`, GEN.563)."""

EJECTA_SPEED_RANGE_KMS = (300.0, 500.0)
"""tuple: The speed an ejected stellar black hole has far out, above the galaxy's escape speed, km/s (GEN.238) [I]: a
stand-in, since no kick model gives the black holes' (`BLACK_HOLE_UNBOUND_FRACTION` is a placeholder too). The population's
density takes the middle of it; each object draws its own, uniform. An ejected neutron star's speed follows the kick law
(`galaxy.natal_kick.ejecta_speeds_kms`)."""

EJECTA_SPEED_QUANTILES = 16
"""int: How many equal-share quantiles of an ejected neutron star's far speed a galaxy's structure holds
(`Structure.neutron_star_ejecta_speeds_kms`); the density takes their mean and each object draws between two of them."""

EJECTA_AGE_GY = 10.0
"""float: How long the ejecta of a galaxy have been leaving, Gyr: with `EJECTA_SPEED_RANGE_KMS` it sets where the
population ends, at `v T` (3 to 5 Mpc, far beyond any grid's reach) (GEN.238) [I]."""

EJECTA_MAX_PER_CELL = 50.0
"""float: The most ejected objects of one kind a cell draws (GEN.238), a program limit: the draw is O(mean), and a galaxy
shaped to hold astronomically many systems in a few cells would otherwise draw more than a cell has room for. The real
galaxy's densest cells draw far fewer."""

HYPERVELOCITY_COUNT_RADIUS_KPC = 100.0
"""float: The radius the count `HYPERVELOCITY_STARS_PER_GALAXY` is read to: "galaxy-wide" is taken as out to 100 kpc,
the halo's extent (docs/design/galaxy/intergalactic-and-halo-contents.md 2.3) [I]."""


PHENOMENON_MIN_MASS_PRESETS = (8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0)
"""tuple: The mass limits the user may pick (GEN.183, Boss 2026-10-10): the
Generate page's slider and `--phenomenon-min-mass` accept only these."""

COMPACT_MIN_MASS_PRESETS = (1.0, 2.0, 4.0, 6.0)
"""tuple: The compact-object limits the user may pick (GEN.195, Boss 2026-10-10): the
lowest mass of a neutron star or black hole the scatter places. `--compact-min-mass`
and the Generate page's compact-object control accept only these, or "star" to use the
stellar mass limit."""

COMPACT_MIN_MASS_STAR = "star"
"""str: `--compact-min-mass`'s word for "use the stellar mass limit" (the default)."""

PHENOMENON_MIN_MASS_SOLAR = 14.0
"""float: The phenomenon scatter's lowest mass, the default of
`PHENOMENON_MIN_MASS_PRESETS` (GEN.167, `planetgen plan
--phenomenon-min-mass`; docs/design/systems/phenomenon-scatter-mass-cut.md, Boss
2026-10-09). The scatter places only the neutron stars and black holes at
or above it; a sector draws the rest when it is filled (GEN.168). An
object's tidal sphere (1.44 pc per cube root of a solar mass) reaches one
4 pc sector edge at about 21 solar masses, so below the cut an object's
pull stays within its own and the neighbouring sectors. 20 first (only the
intermediate-mass black holes scattered: about 2.7e5 rows at a 0.1%
share, about 800 at the 3e-6 of GEN.259, instead of 1.17e9); 8 on 2026-10-10 (Boss, the bottom preset)."""


def phenomenon_rate_per_star(kind):
    """Expected count of `kind` (a PHENOMENON_DENSITY_PC3 key) per star:
    its density over REFERENCE_STELLAR_DENSITY_PC3, times its
    PHENOMENON_RATE_SCALE."""
    return PHENOMENON_DENSITY_PC3[kind] * PHENOMENON_RATE_SCALE.get(kind, 1.0) / REFERENCE_STELLAR_DENSITY_PC3


# --- Galaxy Random-Start Generation (generation/run_galaxy.py) ---

GALAXY_RADIUS_PC = 15000.0
"""
float: The Milky Way's real approximate radius, in parsecs (commonly cited
~15 kpc). Only a fallback for the Galaxy Map's camera range before any
`planetgen plan` has run; once a galaxy is planned, its stored outline
(`galaxy_layer`) is the only bound anything uses -- generation never picks
or accepts an address outside it.
"""

DEFAULT_GENERATE_RADIUS_PC = 12.0
"""
float: The default generate-around sphere, in parsecs: `planetgen
galaxy`'s no-argument random start, and `POST
/api/sectors/<id>/generate-neighborhood` (the Sector page's and the
Generate page's neighborhood) when no radius is given, fill every
not-yet-generated sector this close to their center. Boss (GEN.23,
2026-10-01): "about 10 pc, rounded up" -- to the next whole 4 pc sector,
12 pc. It was 100 ly (about 30.7 pc) before; that sphere now gets only
its bright stars (`BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`).
"""

PROGRESS_BAR_SECONDS = 15.0
"""
float: A step of a generation or maintenance run draws its own progress bar
when it is predicted to take longer than this (Boss, 2026-10-09: "a progress
bar on all substeps that automatically activates on the start of that
sub-item if it's predicted to take longer than 15 seconds"), and the moment
it runs past it if it was not (`planetgen.generation.steps`).
"""

PROGRESS_LOG_SECONDS = 30.0
"""
float: A long step that works in many units writes a log line about every this many seconds (PERF.189), so a log
read without the bar (a job's stage log, `tail -f`) never goes quiet longer than this.
"""

BRIGHT_STAR_BACKFILL_RING_MASSES_SOL = (1.0, 2.0, 5.0, 8.0)
"""
tuple: The bright-star backfill around the generated sectors (GEN.187;
Boss, GitHub issue #952, 2026-10-10), as the initial mass in solar masses
each ring is filled down to, nearest first. The first ring is every
unfilled sector a face away (no diagonals) from a sector the run
generated and gets every living star born with at least 1 solar mass; the
second, a face past the first, gets 2; the third 5; the fourth 8. Beyond
the fourth only the scatter's own mass pass (`galaxy_shape`'s mass limit)
holds. A sector keeps the lightest mass it has been filled to
(`sector_stats.bright_mass_sol`), so a sector a nearer one reaches later
is topped up with only the band it lacks.
"""

RANDOM_START_MAX_PLACEMENT_ATTEMPTS = 1000
"""
int: Retry cap for picking a random, not-yet-occupied sector address before
`generation/run_galaxy.py`'s random-start mode gives up -- generous, since even a
fairly well-populated galaxy leaves overwhelmingly more addresses empty
than occupied (see docs/design/galaxy/galaxy-coordinate-system.md section 9's
storage-analysis addendum: ~320 billion addressable sector slots total).
"""

# --- Navigation Parameters ---

WARP_FACTORS_FOR_NAV = (1, 2, 4, 8, 9, 9.5, 9.9, 9.995)
"""
The warp factors NAV output reports travel time at (see
`planetgen.galaxy.navigation.warp_travel_times`) -- the rows of Boss's warp
table (docs/design/orbits/navigation-frames.md, "Travel speeds"), from warp 1
(exactly c) up to 9.995, where the curve climbs steeply toward its
asymptote at warp 10.
"""

WARP_VELOCITY_EXPONENT = 10 / 3
"""
The exponent of the warp curve's base term, `warp_factor **
WARP_VELOCITY_EXPONENT` (in multiples of light-speed) -- "warp factor to
the 3.33...". Below about warp 9 the whole curve is effectively this term;
see `planetgen.galaxy.navigation.warp_speed_c` for the full formula.
"""

WARP_TRANSITION_STEEPNESS = 9.3575
"""
The logistic steepness `k` in the warp curve's blend term, `1 / (1 +
e^(-k (w - WARP_TRANSITION_MIDPOINT)))`, which hands the curve over from
the plain `w^(10/3)` term to the asymptotic high-warp term around
`WARP_TRANSITION_MIDPOINT`.
"""

WARP_TRANSITION_MIDPOINT = 9.5
"""
The warp factor where the warp curve's blend term is exactly one half
(see `WARP_TRANSITION_STEEPNESS`).
"""

WARP_ASYMPTOTE_COEFFICIENT = 198.9
"""
The numerator of the warp curve's asymptotic term, `198.9 / (10 - w) **
WARP_ASYMPTOTE_EXPONENT`, which grows without bound as `w` nears
`WARP_FACTOR_LIMIT`.
"""

WARP_ASYMPTOTE_EXPONENT = 0.75
"""
The exponent on `(10 - w)` in the warp curve's asymptotic term (see
`WARP_ASYMPTOTE_COEFFICIENT`).
"""

WARP_HIGH_WARP_OFFSET = 1721.7
"""
The constant added to the asymptotic term inside the warp curve's blend,
in multiples of light-speed: `(198.9 / (10 - w)^0.75 + 1721.7 - w^(10/3))`.
"""

WARP_FACTOR_LIMIT = 10
"""
The unreachable top of the warp scale: the asymptotic term divides by
`(WARP_FACTOR_LIMIT - w)`, so a warp factor must be below this.
"""

FOLD_FACTORS_FOR_NAV = (4, 5, 6, 6.5, 7, 7.5, 8, 8.5)
"""
The dimensional fold factors NAV output reports travel time at (see
`planetgen.galaxy.navigation.fold_travel_times`) -- the rows of Boss's fold
table.
"""

FOLD_SPEED_COEFFICIENT = 6
"""
The coefficient in the dimensional fold curve, speed in c = `6 F^4 / (10 -
F)` (see `planetgen.galaxy.navigation.fold_speed_c`).
"""

FOLD_SPEED_EXPONENT = 4
"""
The exponent on the fold factor `F` in the dimensional fold curve.
"""

FOLD_FACTOR_LIMIT = 10
"""
The unreachable top of the fold scale: the fold curve divides by
`(FOLD_FACTOR_LIMIT - F)`, so a fold factor must be below this.
"""

NAV_ADJACENCY_K = 6
"""
How many nearest neighbors each system is connected to when building the
NAV adjacency graph (`planetgen.galaxy.nav_graph.build_knn_adjacency`) that
optimal-route pathfinding runs over. Symmetrized after building (see that
function's docstring), so a system can end up connected to more than `k`
neighbors if others chose it as one of theirs.
"""

NAV_ISLAND_LINKS = 6
"""
How many nearest islands each island of the NAV route graph is linked to
(`planetgen.galaxy.nav_graph.join_islands`, NAV.34). The 6-nearest graph
splits into islands wherever separately generated areas sit apart; in
the hop-length study (2026-10-02), linking each island to its 6 nearest
joined all 714 islands of 2,000 scattered sectors and gave a
cross-galaxy route 1.29 times the direct distance (3 left 2 islands, 1
left 208). Rounds repeat until one island is left, so any value of 1 or
more always joins the graph; a bigger value gives shorter detours.
"""

NAV_CORRIDOR_FRACTION = 0.25
"""float: A route is first searched among the systems within this fraction
of the direct distance of the straight line between its ends (NAV.10),
clamped to `NAV_CORRIDOR_MIN_LY` to `NAV_CORRIDOR_START_MAX_LY`."""

NAV_CORRIDOR_MIN_LY = 25.0
NAV_CORRIDOR_START_MAX_LY = 100.0
NAV_CORRIDOR_MAX_LY = 6400.0
"""float: The corridor's half-width starts between the first two and
doubles up to the last while the route found still has a hop longer than
half the half-width (a hop that long may be hugging the corridor's edge,
with stepping stones just outside it)."""

NAV_CORRIDOR_MAX_SYSTEMS = 300000
"""int: A corridor holding more systems than this is not widened further."""

NAV_BYPASS_MAX_EDGES = 60000
"""int: Most distinct hops NAV.48's bypass test (a route search that may not cross an
uncharted sector) looks at before it gives up and says the answer is unknown."""

NAV_CHART_CONFIRM_SECTORS = 5000
"""int: Charting the sectors of a course (NAV.48) asks for an explicit confirmation past this many sectors
(Boss, 2026-10-10); the PERF.3 disk refusal applies on top."""

NAV_STONE_KINDS = ("black-hole", "neutron-star", "quasar", "hypervelocity-star")
"""tuple: The kinds of scattered phenomenon a jump through unknown space may stop at (NAV.47, Boss,
2026-10-03), beside the scattered bright stars. Planetary nebulae and supernova remnants are not stops
(default decided 2026-10-10)."""

NAV_STONE_CORRIDOR_MIN_PC = 8.0
NAV_STONE_CORRIDOR_FRACTION = 0.1
"""float: The stops of a jump through unknown space are read from the cells within this half-width of
the straight jump: the larger of `NAV_STONE_CORRIDOR_MIN_PC` and this fraction of its length."""

NAV_STONE_MAX_CELLS = 3000
"""int: Most uncharted cells read for the stops of one jump (the corridor is narrowed to fit)."""

NAV_STONE_MAX_STONES = 5000
"""int: Most candidate stops one jump considers (the nearest to the straight jump are kept)."""

NAV_UNKNOWN_JUMP_PENALTY = 0.5
"""float: A chain of stops replaces a jump through unknown space only when it is no longer than the
jump times one plus this (the direct jump is kept in reserve at that cost)."""


# --- Facilities (schema v42) ---

FACILITY_KINDS = {
    "colony": "A settlement people live in; makes its world inhabited.",
    "outpost": "A small crewed post: a research, relay or survey station.",
    "mining-colony": "A settlement built to work an asteroid belt.",
    "station": "A crewed space station.",
    "starbase": "A large station: a port, yard and base in one.",
}
"""dict: Every facility kind (`facilities.kind`), with a one-line meaning."""

FACILITY_RULES = {
    # (placement, host type): the kinds allowed there. Boss, 2026-09-30:
    # gas giants take orbital facilities only; terrestrial worlds and moons
    # take colonies and orbital facilities; asteroid belts take outposts and
    # mining colonies, asteroid fields outposts; a star can have an outpost
    # (or a station) in orbit around it; stand-alone ones park in space.
    ("terrestrial", "planet"): ("colony", "outpost"),
    ("terrestrial", "moon"): ("colony", "outpost"),
    ("orbital", "star"): ("outpost", "station", "starbase"),
    ("orbital", "planet"): ("outpost", "station", "starbase"),
    ("orbital", "moon"): ("outpost", "station", "starbase"),
    ("asteroid", "asteroid_belt"): ("outpost", "mining-colony"),
    ("asteroid", "asteroid_field"): ("outpost",),
    ("standalone", "space"): ("outpost", "station", "starbase"),
}
"""dict: Which facility kinds (`FACILITY_KINDS`) may go where, keyed by
`(placement, host_type)`. A terrestrial facility also needs a terrestrial
(`body_type = 't'`) planet or moon -- `facilities.check_facility`."""

FACILITY_RATING_GREEN = 0.66
"""float: A suggested rating score at or above this is Green (POP.10, `facility_suggest`)."""

FACILITY_RATING_YELLOW = 0.33
"""float: At or above this (and below `FACILITY_RATING_GREEN`) is Yellow; below it, Red."""

FACILITY_RATING_NOISE = 0.1
"""float: A suggested score gets a stable draw of up to this much either way, seeded from the facility."""

FACILITY_RATING_UNKNOWN = 0.5
"""float: What a missing input (an index nobody drew, a host with no habitability) counts as."""

FACILITY_RATING_ORBITAL_HABITABILITY = 0.6
"""float: The habitability an orbital or free-flying facility counts as: sealed, off the surface."""

FACILITY_RATING_UNAFFILIATED_PENALTY = 0.2
"""float: Taken off a facility's crime score when no polity polices it."""

FACILITY_KIND_SIZE = {"colony": 1.0, "starbase": 0.8, "mining-colony": 0.7, "station": 0.6, "outpost": 0.4}
"""dict: How much room a kind of facility has to house people, 0 to 1 (the housing score)."""

FACILITY_HOST_RICHNESS = {"asteroid_belt": 0.9, "asteroid_field": 0.9, "planet": 0.6, "moon": 0.6, "star": 0.5,
                          "space": 0.2}
"""dict: What a host offers a facility to work with, 0 to 1 (the resources score)."""

FACILITY_DEFAULT_ORBIT_RADII = 3.0
"""float: An orbital facility around a planet or moon with no distance
given orbits at this many of its host's radii."""

FACILITY_ORBIT_FLOOR = 1.01
"""float: The lowest orbit the facility form offers, in host radii: just
above the surface (`facilities.orbit_limits`)."""

FACILITY_ORBIT_FALLBACK_RADII = 1000.0
"""float: The highest orbit the facility form offers, in host radii, for
a host with no stored sphere of influence (a planet or moon without a
`hill_radius_km`, a star without a heliosphere)."""

FACILITY_ORBIT_STEPS = 1000
"""int: The facility form's orbit slider runs from 0 to this; each step
is the same factor of distance (`facilities.distance_from_step`)."""

# --- Galaxy pre-placement (schema v43) ---

BRIGHT_STAR_MIN_LUMINOSITY_SOL = 9000.0
"""float: Every star at least this bright (solar luminosities) is generated
and placed galaxy-wide right after `planetgen plan`, before any sector is
filled (`bright_stars`, schema v43). Its sector is still generated later,
around it. 500 on 2026-10-01, then 1000 (GEN.30), with the dimmer stars
filled in near generated sectors by the backfill (luminosity tiers, then the
mass rings of GEN.187); 3000 on 2026-10-10 (GEN.184, Boss), then 5000 the same day. The user
picks it from the ladder in `generation/luminosity_floor.py`, from
`BRIGHT_STAR_FLOOR_MIN_SOL` to `BRIGHT_STAR_FLOOR_MAX_SOL`. A galaxy
already scattered keeps the level it was scattered at. The value a scatter
used is stored in `galaxy_shape.bright_star_min_luminosity_sol`, and
filling reads that, not this, so retuning it can't make a fill
double-count or skip stars. See
/mnt/project-files/galaxy-studies/bright-star-preplacement-plan.md."""

BRIGHT_STAR_FLOOR_MIN_SOL = 2500.0
"""float: The lowest luminosity floor a user may choose (GEN.184, Boss
2026-10-10 03:01Z)."""

BRIGHT_STAR_FLOOR_MAX_SOL = 4_000_000.0
"""float: The highest luminosity floor a user may choose (GEN.184)."""

BRIGHT_STAR_FLOOR_MIN_STEP_SOL = 100.0
"""float: The gap between presets at the lowest floor, for the best control
there (GEN.184)."""

BRIGHT_STAR_FLOOR_MAX_STEP_SOL = 500_000.0
"""float: The gap between presets near the highest floor (GEN.184)."""

# --- Population and politics (schema v44, POP.1 to POP.4) ---
# See docs/design/habitability/population-and-politics.md.

LIFE_WORLD_STAGES = ("multicellularity", "technological_civilization")
"""tuple: The evolutionary milestones (`evolution.MILESTONE_KEYS`) that
make a planet a life world with a named dominant species. Simpler
biospheres keep their life chemistry but get no `species` row."""

CIVILIZATION_CHANCE = 0.001
"""float: The chance that a planet whose timeline reached the
`technological_civilization` milestone (which only says the window for one
has opened -- about one generated system in seven reaches it) actually has
a civilization now. 1 in 1,000 gives about one civilization per 7,000
systems, so spacefaring neighbors sit roughly 150-200 ly apart and their
territories (`TERRITORY_REACH_CAP_LY`) just about meet. A system generated
with intelligent life forced on (`system_configs.intelligent_life = 1`)
always has one. The rest still get a named dominant species, with no
civilization."""

CIVILIZATION_MIN_AGE_YEARS = 100.0
"""float: The youngest a technological civilization can be. Its age is a
log-uniform draw between this and the time since its homeworld's
`technological_civilization` milestone (`population.civilization_age`)."""

CIVILIZATION_ERAS = (
    ("Industrial", 0.0, False),
    ("Interplanetary", 300.0, False),
    ("Interstellar", 2_000.0, True),
    ("Established", 50_000.0, True),
    ("Ancient", 1_000_000.0, True),
    ("Elder", 100_000_000.0, True),
)
"""tuple: `(era, starts at age in years, spacefaring)` in ascending order.
A civilization is in the last era whose start its age has reached. Tuning
these needs no migration: every population pass recomputes `species.era`
and `species.spacefaring` from the stored age."""

TECH_DOMAINS = ("energy", "materials", "information", "medical", "propulsion", "defense")
"""tuple: The six tech domains of the Technological Assessment, in the order of `TECH_WEIGHTS`
(Energy, Materials, Information, Medical, Propulsion, Defense)."""

TECH_WEIGHTS = (0.25, 0.20, 0.20, 0.15, 0.10, 0.10)
"""tuple: The weights of the six domain indices in a tech level (Boss's Technological Assessment, "Domain
Weighting System"); they sum to 1. TL = 0.25 EI + 0.20 MI + 0.20 II + 0.15 MeI + 0.10 PI + 0.10 DI."""

TECH_DOMAIN_OFFSETS = (0.066, 0.066, -0.171, -0.053, 0.066, -0.053)
"""tuple: Each domain's shift in decades of log10 effective age (negative is earlier), solved so the median
index vector at 266 years is Earth's 2026 (E 3.0, M 3.0, I 3.4, Me 3.2, P 3.0, D 3.2)."""

TECH_CURVE_STEEPNESS = 0.97
"""float: `k`, the index gained per decade of log10 age at the middle of a domain's logistic."""

TECH_CURVE_CENTRE = 2.655
"""float: `m`, the centre of the logistic in log10 years (about 450 years, where an index is 3.5)."""

TECH_PACE_SIGMA = 0.30
"""float: The spread of a species' shared pace, in decades (a factor of 2 in time per sigma)."""

TECH_DOMAIN_SIGMA = 0.09
"""float: The scale of each domain's own t(4) noise, in decades (standard deviation about 0.13)."""

TECH_CEILING_MEAN = 6.2
TECH_CEILING_SD = 0.6
TECH_CEILING_DOMAIN_SD = 0.4
TECH_CEILING_SOFTNESS = 5.0
"""float: A species' ceiling is a normal draw of `TECH_CEILING_MEAN` and `TECH_CEILING_SD` clipped to 4 to 7,
each domain's `TECH_CEILING_DOMAIN_SD` away from it; the index is a soft minimum of the curve and the ceiling
with this sharpness."""

TECH_STAGNATION_CHANCE = 0.12
TECH_STAGNATION_LOG10_YEARS = (2.7, 6.0)
"""A species freezes at a plateau age (log10 years, drawn evenly between these) with this chance."""

TECH_COLLAPSE_RATE_PER_DECADE = 0.10
TECH_COLLAPSE_KEEP = (0.01, 0.30)
"""Collapses are a Poisson process of this rate per decade of log10 age past `CIVILIZATION_MIN_AGE_YEARS`;
at each one the effective age is multiplied by a draw between these."""

TECH_SPIKE_CHANCE = 0.03
TECH_SPIKE_SIZE = (1.5, 3.0)
"""A lopsided species has one domain raised by a draw between these after the prerequisite clamps."""

TECH_PREREQUISITES = (
    ("information", ("materials",), 1.5),
    ("medical", ("information", "materials"), 1.5),
    ("propulsion", ("energy", "materials", "information"), 1.0),
    ("defense", ("energy", "materials"), 1.0),
)
"""tuple: `(domain, the domains it needs, margin)`: the domain is at most the lowest of those plus the margin,
applied in this order."""

TECH_BANDS = (
    (0.0, "Primal"), (0.5, "Agrarian"), (1.5, "Industrial"), (2.5, "Atomic"), (3.5, "Fusion"),
    (4.5, "Gravitic"), (5.5, "Antimatter"), (6.5, "Transcendent"),
)
"""tuple: `(tech level where the band starts, name)` in ascending order."""

TERRITORY_BASE_REACH_LY = 5.0
"""float: A polity's territory reach (light years) the moment its species
goes interstellar. Reach grows with the square root of age from there:
`TERRITORY_BASE_REACH_LY * sqrt(age / interstellar start)`."""

TERRITORY_REACH_CAP_LY = 100.0
"""float: The farthest any polity's territory reaches (light years),
reached at about 800,000 years."""

GOVERNMENT_FORMS = (
    "Hegemony", "Concord", "Republic", "Union", "Directorate", "Dominion", "Assembly",
    "Commonwealth", "Collective", "Compact", "Federation", "Sovereignty", "Ascendancy", "League",
)
"""tuple: The form of government a polity's name ends with ("the Voranthi
Concord")."""


WEB_LOCK_WAIT_SECONDS = 5
"""int: How long a web request waits for a row lock in the control schema (the account, session and audit tables)
before it gives up, where the server's own limit is 50 s (ADM.53). The app has five threads and a 60 s request
limit, so a handful of requests stuck behind one long transaction would otherwise stop every page, the login
included."""

# --- Species environment limits (POP.12) -------------------------------------------------------

SPECIES_LIMIT_MARGIN = {"pressure_factor": 1.5, "temperature_c": 10.0, "wet_bulb_c": 2.0, "dose_factor": 3.0,
                        "o2_factor": 1.25}
"""dict: How far inside Blue a species' homeworld sits when its limits are generated from it
(`species_limits.generated`): the homeworld's pressure is 1.5 times inside the Blue pressure edges, its
temperature cycle 10 degrees inside the Blue temperature edges, a wet bulb limit 2 degrees above its wet bulb,
its dose a third of the Blue dose edge, and the chronic inspired O2 limit a quarter above its own. planetGen
defaults (Boss's default of 2026-10-15: a species passes its own homeworld)."""

SPECIES_EQUIPMENT_TECH = {
    "cold-weather clothing": 0.0, "cooling clothing": 0.0, "sun-cover clothing": 0.0, "weather clothing": 0.0,
    "insulated suit": 2.0, "liquid-cooled suit": 2.0, "heated diving suit": 2.0, "full-cover suit": 2.0,
    "thermal suit": 2.0,
    "heated sealed suit": 3.0, "heat-shielded suit": 3.0, "chemical-barrier suit": 3.0, "sealed suit": 3.0,
    "pressure suit": 3.0, "armoured suit": 4.0, "vehicle only": 99.0, "powered frame": 3.0,
    "dust filter": 2.0, "dust-tight seals": 2.0, "gas filter": 2.5, "full-face gas filter": 2.5,
    "powered gas filter": 3.0, "CO2 scrubber": 2.5,
    "gas-mix mask": 2.0, "self-contained air": 3.0, "mixed-gas rig": 3.5, "isolated air": 4.0,
    "shielded garment": 1.0, "storm vest": 2.5, "hardened suit": 3.5, "shielded shelter": 2.0, "hardened habitat": 4.0,
}
"""dict: The least tech level (0 to 7) a species needs to make each item of the survival equipment
(`physics/equipment.py`, by the item's name; a "hardened" suit needs the larger of the hardened suit's and its own), for
the "could live there with equipment" list (POP.12, POP.15). It replaces the five-tier tuple (0, 2, 2.5, 3, 4). Vehicle
only is never within reach (99). A planetGen default."""

COLD_END_LABORATORY_CORRECTION = True
"""bool: Whether the cold-end vapour pressures of nitrogen, carbon monoxide, methane, oxygen and ammonia are the 2026
Grundy preprint's set (True), which multiplies the 2024 laboratory fits by 1.19 +/- 0.03, or the 2024 papers' (False)
(`physics.saturation`, GEN.497). A switch until the preprint is read in print. A planetGen default."""

ORBIT_UPDATE_ROW_CEILING = 10_000_000
"""int: The most rows one orbital-update run moves (PERF.103), counted across every table it moves. When more rows
are due, the run moves this many (the oldest due of each table first, a share of it per table by how many are due),
warns with the backlog and leaves the rest due for the next night. 0 moves everything that is due, however many.
`planetgen.cli.orbits --row-ceiling` sets it for one run. A planetGen default."""

# --- Slow page requests (PERF.105) -------------------------------------------------------------

SLOW_PAGE_ROW_BUDGET = 2_000_000
"""int: A page request that reads more rows than this is slow (`web/lib/slowpages.py`): the rows are the
database's own estimate (`information_schema`) of the tables the route reads in full, and are used only until a
build of that exact request has been timed. A planetGen default."""

SLOW_PAGE_SECONDS = 3.0
"""float: A page request whose last build took longer than this is slow, whatever its row estimate said. A
measured time always beats the estimate. A planetGen default."""

SLOW_PAGE_ROWS_PER_SECOND = 250_000
"""int: How many rows a build reads per second, for the first estimate of how long a request will take, before one
has been timed. A planetGen default."""

SLOW_PAGE_PATH_SECONDS = 600
"""int: How long the slowest recent build of a route (any request on it) counts for a request of that route that
has not been timed itself (`slowpages.cost`). A planetGen default."""

SLOW_PAGE_JOB_SECONDS = 3600
"""int: How long a slow request's job (its id, progress and finished state) is kept in Redis. A planetGen default."""

SLOW_PAGE_FAILED_SECONDS = 30
"""int: How long a failed or cancelled job is shown to the visitors who ask for the same request, before the next
one queues a new job. A planetGen default."""

SLOW_PAGE_DELIVER_SECONDS = 300
"""int: How long a finished page job keeps its answer for the visitors waiting on it when the page cache would not take
it (data changed while it was built, or Redis had no room), so the page that asked is answered instead of queuing the
same job again (UX.234). It is served only while the page cache has no copy, and only to a page that waited for the job.
A planetGen default."""

SLOW_PAGE_ESTIMATE_SECONDS = 60
"""int: How long a web process keeps the table row estimates before it asks the API again. A planetGen default."""

SLOW_PAGE_TIMING_SECONDS = 7 * 86400
"""int: How long the build time of a request is remembered in Redis. A planetGen default."""

# --- Nightly New galaxy smoke run (PERF.183) ---------------------------------------------------

SMOKE_RADIUS_PC = 12.0
"""float: The radius, in parsecs, of the one neighborhood the nightly smoke run (`planetgen.cli.smoke`) fills, with
the scatters on: small enough to finish in minutes, large enough to run every step of a New galaxy job. A planetGen
default."""

SMOKE_BUDGET_MB = {"mysql": 1800, "workers": 1400, "step": 600, "host": 5000}
"""dict: The most memory, in megabytes, the smoke run may peak at in any one step of the job, by group (`memwatch`):
`mysql` (every mysqld), `workers` (every planetgen worker process), `step` (the step's own process tree) and `host`
(the whole host's memory in use). A step over any of them fails the run, so a change that makes a step hold far more
memory is caught the night it merges and not on the first big galaxy. Set about twice the first measured run (906 MB
MySQL, 661 MB workers, 271 MB a step, 2,535 MB host in the worst step, 10-24). A planetGen default."""

SMOKE_STEP_BUDGET_MB = {"Scatter the phenomena": {"mysql": 3600, "workers": 1500, "step": 4200, "host": 9700}}
"""dict: Budgets of single steps (by the step's label on the Generate page) that replace `SMOKE_BUDGET_MB` for it. The
phenomena scatter held 6,391 MB in its own process on the first measured run (a 12 pc neighborhood, host 8,655 MB),
far above every other step; PERF.187 made its block light totals write each block as the shells pass it and its nebula
block totals write one level at a time, and the run then measured 2,100 MB in the step (the nebula block totals pass,
which still holds the cover of the level-3 blocks as dicts), 1,819 MB MySQL, 644 MB workers and 4,853 MB host
(10-24). Set about twice that. A planetGen default."""
