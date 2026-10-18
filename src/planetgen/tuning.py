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

from planetgen.generation.phenomena import registry as _phenomena_registry
from planetgen.physics import constants as _physical_constants

# --- Planet Generation Parameters ---
# Domingos, Winter & Yokoyama (2006), MNRAS 373:1227, "Stable satellites
# around extrasolar giant planets" -- a prograde moon on a circular orbit
# stays bound out to about 0.4895 of its planet's Hill radius; beyond that
# the star strips it (physics.planets.generate_moons' outer limit).
MOON_PROGRADE_STABLE_HILL_FRACTION = 0.4895

# Each gas-giant class's mass range in Earth masses (GEN.34): a giant's
# mass is drawn from it first and its radius follows from the mass-radius
# relation (physical_constants.GIANT_NEPTUNIAN_MASS_RADIUS). Ice giants
# (Uranus 14.5, Neptune 17.1); gas giants from 0.1 Jupiter masses to the
# 13 Jupiter-mass deuterium-burning limit (Saturn 95, Jupiter 318); gas
# dwarfs from the smallest H/He-envelope planets that fill class T's
# 15,000 km floor. A class with its own "density_range" (none today) keeps
# drawing a radius and a density instead.
GIANT_CLASS_MASS_RANGE_EARTH = {"I": (10.0, 60.0), "J": (30.0, 13 * 317.8), "T": (6.0, 40.0)}
GIANT_DEFAULT_MASS_RANGE_EARTH = (10.0, 13 * 317.8)

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

# --- Star System Generation Parameters ---
INITIAL_PLANET_DISTANCE_FACTOR = 0.55
# Where a system's first slot sits, as a fraction of its star's habitable
# zone inner edge, and how far out each next slot sits, as a ratio of the
# last one's distance (GEN.37). Adjacent real planets sit at period ratios
# of about 1.5-2.5 (Kepler multis; the solar system's distance ratios run
# 1.4-3.4), i.e. distance ratios of about 1.3-1.8, so placement is
# geometric and scales with the star's light; slots stop at the disk's
# outer edge (StarSystem._disk_outer_edge_au). Before, the first slot sat
# at 0.55 AU per solar mass and each next one a fixed number of AU
# further, which put about 96% of planets in the cold zone; now about 32%
# are hot, 6% in the ecosphere and 62% cold (the solar system: 25%, 12-25%
# and 50-62%).
FIRST_PLANET_HZ_FRACTION_RANGE = (0.05, 0.4)
PLANET_SPACING_RATIO_RANGE = (1.3, 1.8)
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
# not the flat "1 + log10(solar_masses)" scaling the old
# StarSystem.estimate_num_objects formula used. Disk dust-mass-vs-
# stellar-mass surveys (Andrews et al. 2013, Taurus; Pascucci et al. 2016,
# multi-region) find M_dust ~ M_star^1.8 in young (~1-3 Myr) star-forming
# regions, steepening further (~M_star^2.7) in older ones. This scales
# `physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2` for any given
# star relative to the Sun -- see `formation.disk_surface_density_scale`,
# used by `StarSystem._estimate_max_objects_from_disk_physics`. The
# younger-region exponent is used since this generator has no notion of a
# system's disk-formation age (only its current, post-formation age).
DISK_MASS_STELLAR_MASS_EXPONENT = 1.8

# Real disks are truncated far short of a star's own galactic-tidal Hill
# sphere (`Star.system_perimeter`, tens to hundreds of thousands of AU) --
# viscous spreading and photoevaporation cut them off at tens to a few
# hundred AU (ALMA disk-size surveys). Rather than a flat AU figure, this
# scales with the same star-dependent quantity that sets where solids can
# even condense in the first place: the snow line
# (`physical_constants.SNOW_LINE_AU_AT_1_LSUN`/`formation.snow_line_au`) --
# our own Solar System's own giant-planet/Kuiper-belt region extends to
# roughly this same multiple (~18x) of its own 2.7 AU snow line. See
# `StarSystem._estimate_max_objects_from_disk_physics`.
DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER = 18

# The farthest a planet zone in a hierarchy of stars reaches (GEN.314,
# `physics.planet_zones`): a protoplanetary disk around a star or a pair is
# not larger than this, whatever the stability limit allows.
PLANET_ZONE_DISK_CAP_AU = 100.0

# The disk-physics walk (`StarSystem._estimate_max_objects_from_disk_physics`)
# first counts how many isolation-mass "oligarchs" (Kokubo & Ida
# 2000/2002) the disk's solid budget can support, spaced by their own
# mutual Hill radius -- but not every oligarch survives as a final planet.
# N-body integrations of the subsequent giant-impact phase (Chambers 2001)
# show a large fraction of oligarchs merge or get ejected before the
# system settles down; only a minority survive as the final planet
# count. This is the multiplicative attrition factor applied to the raw
# oligarch count. Tuned toward the middle of that literature's range
# (rather than one precise figure) so a solar-mass star's typical final
# count stays in the same well-tested, playable range this generator
# already verified via repeated full-suite runs (see CHANGELOG.md).
GIANT_IMPACT_SURVIVAL_FRACTION = 0.4

# How many times StarSystem.__init__ retries its whole placement loop (fresh
# object count, positions, and validate_system pass, same star) before
# giving up when HABITABLE_WORLD=True/ASTEROID_BELT=True still isn't
# satisfied afterward. Needed because validate_system's own orbital-overlap
# correction can, in rare cases, push a deliberately-placed guaranteed body
# (e.g. the forced Class M world) into a zone its class no longer supports --
# physics.planets.reconcile_zone_and_class correctly reclassifies it away
# rather than reporting a physically inconsistent class, which can undo the
# very guarantee that placement was satisfying. Retrying with a fresh
# placement is simpler and more robust than trying to reshuffle every
# neighbor's spacing to protect one body's zone in place.
MAX_SYSTEM_GENERATION_ATTEMPTS = 8

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
These components are used to generate the composition of asteroid belts.
"""

# --- Space Sector Generation Parameters ---

# The standard sector edge length, in whole parsecs: the galaxy grid's
# ring width, layer height and (near enough) slot arc, and the cube edge of
# a standalone sector. A design choice, not a physical constant. 4 pc
# (~13.05 ly, ~2,220 ly^3 of volume, ~6.3 systems at local density -- see
# `physical_constants.LOCAL_STELLAR_DENSITY_LY3`) puts the default Milky
# Way shape's one-star-per-sector edge at ~50,000 ly, the real disk's
# radius; see `docs/design/galaxy-coordinate-system.md`, "Sector size".
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
HABITABLE_ZONE_BUFFER_AU = 0.2
ROUND_RADIUS_KM = 2
SCIENTIFIC_NOTATION_DECIMAL_PLACES = 2
ROUND_TEMPERATURE_NEAREST_HUNDRED = -2

# --- Star Generation Parameters ---
# The current age of the universe (~13.8 Gy, standard cosmological
# consensus) -- an absolute ceiling on any star's age, regardless of how
# long its own spectral class could theoretically keep burning. Without
# this, extremely long-lived classes (M dwarfs: STAR_EVOLUTION["M"]'s
# lifespan_gy runs up to 5500 Gy, a real theoretical figure -- no red dwarf
# has ever actually died of old age) could roll an "age" of hundreds of
# billions of years. That age is still less than the star's own enormous
# lifespan (the existing age <= lifespan invariant holds), but far older
# than any star that could actually exist yet, since the universe itself
# isn't old enough. White dwarfs already have their own bounded cooling-age
# range (WHITE_DWARF_MIN/MAX_AGE_GY) and don't need this cap.
UNIVERSE_AGE_GY = 13.8
MIN_INITIAL_STAR_AGE_GY = 0.1
MIN_INITIAL_STAR_AGE_LIFESPAN_RATIO = 0.01
MAX_INITIAL_STAR_AGE_LIFESPAN_RATIO = 0.9
UNREACHABLE_PLANET_AGE_MIN_LIFESPAN_RATIO = 0.85
OLD_STAR_AGE_LIFESPAN_RATIO = 1/3
YOUNG_STAR_AGE_LIFESPAN_RATIO = 1/3
MAX_PLANET_AGE_ADJUSTMENT_FACTOR = 0.95
WHITE_DWARF_AGE_ADDITION_GY = 5
# White dwarf cooling-age bounds: a WD's spectral letter reflects its current
# temperature, not its progenitor's mass, so its age can't be drawn as a
# fraction of a STAR_EVOLUTION main-sequence lifespan the way other classes
# are. The oldest known white dwarfs (e.g. in globular cluster M4) have
# cooling ages of ~12-13 Gy, close to the age of the universe.
WHITE_DWARF_MIN_AGE_GY = 0.1
WHITE_DWARF_MAX_AGE_GY = 12.0
# For any evolved-star Yerkes class other than main-sequence (V), white
# dwarf (VII/D), or subdwarf (VI, see below) -- giants, subgiants, bright
# giants, supergiants, hypergiants -- the current spectral letter reflects
# only present-day temperature, not a lifespan (see
# get_star_evolutionary_profile in utils.py for the same issue on the
# planet-life side). Instead, the star's own already-generated mass is run
# through the main-sequence lifetime of physics/ms_lifetime.py (the fit of
# Hurley, Pols and Tout 2000 to detailed models, GEN.467; the STAR_EVOLUTION
# table above was derived from the older t = 10 Gy * (M/Msun)^-2.5 law),
# then extended by the standard rule-of-thumb that a star spends about 90%
# of its total lifetime on the main sequence -- so the star must already be
# at least as old as that main-sequence lifespan (it has to have finished
# that phase to be observed as an evolved class), and at most its total
# (MS + post-MS) lifespan.
# The heavy-element mass fraction the lifetime's [Fe/H] is measured from:
# the proto-solar Z of the composition scale, Asplund 2021 (stellar-
# composition.md D1). A solar mass lasts 10.2 Gy at it.
MS_LIFETIME_REFERENCE_Z = 0.0154
MS_LIFESPAN_FRACTION_OF_TOTAL = 0.9
# Subdwarfs (Yerkes VI, e.g. real sdB/sdO stars) are excluded from the
# mass-derived model above -- not just from its mass-sampling rejection
# check (_sample_evolved_star_mass_sol), but from the model itself. Real
# subdwarfs are thought to form via binary mass-stripping near the tip of a
# lower/intermediate-mass progenitor's red-giant branch, not through
# ordinary single-star post-main-sequence evolution, so a star's *current*
# (post-strip) mass has no direct relationship to its progenitor's own
# main-sequence lifespan the way it does for an ordinary giant/supergiant --
# using it that way is what let a generated subdwarf's age come out older
# than the universe (formerly tracked in TODO.md's Future Ideas). Real
# subdwarf progenitors are typically old, low-mass population stars (a
# short-lived, higher-mass star wouldn't have had time to reach the RGB
# tip and get stripped), so age is instead drawn directly from an
# old-population-biased range, independent of this star's own mass. The
# core-helium-burning subdwarf phase itself is short relative to that age
# (real sdB stars: roughly 0.05-0.3 Gy) -- see
# Star._calculate_initial_star_age_and_lifespan's Yerkes-VI branch.
SUBDWARF_MIN_AGE_GY = 1.0
SUBDWARF_MAX_AGE_GY = 13.5
SUBDWARF_REMAINING_PHASE_MIN_GY = 0.05
SUBDWARF_REMAINING_PHASE_MAX_GY = 0.3
# A sub-solar-mass progenitor (roughly under ~0.93 Msun, per the
# ms_lifetime.py lifetime above) already has
# a main-sequence lifespan longer than the universe itself -- it couldn't
# possibly have evolved off the main sequence yet, so it must never be
# accepted as an evolved-class star's progenitor mass in the first place
# (see Star._sample_evolved_star_mass_sol, used by Star.generate_star).
# This caps the reject-and-resample loop that enforces that, the same way
# SECTOR_MAX_PLACEMENT_ATTEMPTS caps SpaceSector's own reject-and-resample
# loop: most classes' allowed mass range sits entirely clear of the cutoff
# (so the very first draw is accepted) and only the class III (Giant) low
# end (0.8-~0.93 Msun, out of an 0.8-8 Msun range) ever needs a resample,
# so 100 attempts leaves effectively zero chance of exhausting this for any
# class it's actually applied to.
EVOLVED_STAR_MASS_MAX_RESAMPLE_ATTEMPTS = 100
RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD = 100000
SPECTRAL_PROBABILITIES_LARGE_STAR = {'O': 10, 'B': 20, 'A': 30, 'F': 30, 'G': 10, 'K': 0, 'M': 0}
SPECTRAL_PROBABILITIES_NORMAL = {'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45}

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

B_SUBCLASS_LOWER_MASS_SOL = (("B-early", 9.0), ("B-mid", 5.0))
"""tuple: The birth mass at which the B-star sub-classes begin, heaviest first (Table 13's bins); a B star below the
last is late B and keeps the letter's row."""

BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS = {
    'O': 0.90, 'B': 0.65, 'A': 0.55, 'F': 0.47, 'G': 0.44, 'K': 0.40, 'M': 0.26,
    'B-mid': round(1.0 - MULTIPLICITY_BY_PRIMARY_MASS_BIN["mid B"]["single"], 2),
    'B-early': round(1.0 - MULTIPLICITY_BY_PRIMARY_MASS_BIN["early B"]["single"], 2),
}

# --- Star population model (stellarEvolution.py) ---
# A random star is drawn as physics, not from a letter table: a mass from
# the initial mass function, an age from the star-formation history, and
# the star's present state (main sequence, subgiant, giant, supergiant or
# white dwarf) from how far that age is into its own lifetime. Its letter,
# subclass and Yerkes class follow from the resulting temperature and
# luminosity. SPECTRAL_PROBABILITIES_NORMAL above is kept as the census the
# model is tested against, not as a draw.

# Kroupa (2001), MNRAS 322:231 -- the initial mass function as a broken
# power law dN/dM ~ M^-alpha: alpha = 1.3 from 0.08 to 0.5 Msun and 2.3
# above, up to the ~150 Msun upper limit of real stars. IMF_BREAKS_SOL are
# the segment edges and IMF_SLOPES each segment's alpha.
IMF_BREAKS_SOL = (0.08, 0.5, 150.0)
IMF_SLOPES = (1.3, 2.3)

# The default star-formation history: ages uniform over the thin disk's
# 0-10 Gy. Population-dependent ages by position come on top of this (the
# sector fill passes an age instead of letting the star draw one).
STAR_FORMATION_AGE_RANGE_GY = (0.0, 10.0)

# Post-main-sequence phases, as multiples of the main-sequence lifetime
# (physics/ms_lifetime.py, above). A star lighter than HELIUM_FLASH_MASS_SOL
# burns helium in its core after the red giant branch (the horizontal branch)
# and then climbs the asymptotic giant branch before it becomes a white dwarf
# (GEN.286). For it the subgiant phase ends at, and the red giant branch ends
# at, the times of physics/giant_branch.py (GEN.502, Hurley, Pols and Tout
# 2000: 1.05 and 1.12 of the main-sequence lifetime at a solar mass, from the
# core-mass-luminosity relation of the branch), and HORIZONTAL_BRANCH_DURATION_GY
# and AGB_MS_FRACTION are added for it. A heavier star keeps the two figures
# below, a subgiant phase up to 1.1 and a giant or supergiant phase up to 1.2,
# then a remnant (real post-MS phases take roughly 10-20% of the
# main-sequence time) [R]: its core helium burning is not modelled.
SUBGIANT_PHASE_END_MS_FRACTION = 1.1
GIANT_PHASE_END_MS_FRACTION = 1.2
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
SUPERGIANT_MIN_MASS_SOL = 8.0
# Giants at or above this initial mass are bright giants (II), not III.
BRIGHT_GIANT_MIN_MASS_SOL = 4.0
# Giant-branch luminosities (Lsun) and temperatures (K). A giant lighter than
# HELIUM_FLASH_MASS_SOL is on the red giant branch of physics/giant_branch.py
# (GEN.502): the helium core grows by the hydrogen shell burning at the star's
# luminosity, dMc/dt = A_H L, and the luminosity climbs with the core mass as
# L = min(B Mc^q, D Mc^p) (Hurley, Pols and Tout 2000, eq. 37, p = 6 and q = 3
# at low mass), so the time at each luminosity, the branch's duration and its
# tip (where the core helium ignites, about 2,700 Lsun at the reference
# composition, 1,900 at Z 0.0001 against Serenelli et al. 2017's 1,650 to
# 1,780) follow from the relation, the slope flattening from -0.83 to -0.67
# along the branch [V]. It reproduces the paper's Table 1 durations
# and the observed ratio of horizontal-branch stars to red giants brighter
# than the horizontal branch, 1.44 +- 0.06 for the 41 halo clusters of Salaris
# et al. 2004, table 2 (1.51 +- 0.05 for all 57) [P]. The shares of the branch
# spent above 50, 100, 500 and 1,000 Lsun (8.7, 4.7, 1.0 and 0.42%) are 10 to
# 30% under the 9.6, 5.5, 1.3 and 0.59% of the fit of arXiv 2503.17772
# (eq. 1, one data point [P]), whose slope flattens further, to -0.56 at the tip.
# The two figures below are what is left of the older single-exponent law
# (GEN.287): the intermediate-mass giants, lighter than
# BRIGHT_GIANT_MIN_MASS_SOL but not lighter than HELIUM_FLASH_MASS_SOL, rise from
# SUBGIANT_MAX_BRIGHTENING times the main-sequence luminosity to the tip, the
# time at each luminosity falling as L^-GIANT_BRANCH_TIME_SLOPE per unit of ln L,
# and the asymptotic giant branch (AGB_LUMINOSITY_RANGE_SOL) takes the same slope
# [R, recalled, unconfirmed: Refsdal and Weigert 1970; Salaris and Cassisi 2005].
GIANT_BRANCH_TIP_LUMINOSITY_SOL = 2500.0
GIANT_BRANCH_TIME_SLOPE = 0.78
BRIGHT_GIANT_LUMINOSITY_RANGE_SOL = (1000.0, 10000.0)  # log-uniform
GIANT_TEMPERATURE_RANGE_K = (3500.0, 5000.0)
# The horizontal branch (core helium burning): one luminosity range, whatever
# the colour; the colour (temperature) is drawn from [Fe/H] and age
# (HORIZONTAL_BRANCH_COLOURS below). Red clump and red horizontal-branch stars
# are the metal-rich, RR Lyrae the middle, blue horizontal-branch the
# metal-poor and extreme horizontal-branch (hot subdwarfs) a few percent of
# old stars (docs/design/realistic-galactic-core.md 5.3) [R, recalled, unconfirmed:
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
# The asymptotic giant branch: the same time-at-luminosity law, between these.
AGB_LUMINOSITY_RANGE_SOL = (100.0, 2500.0)
# A subgiant brightens to about three and a half times its main-sequence luminosity (the
# zero-age relation under-counts the turn-off star's own brightening, and the base
# of the red giant branch of an old 0.85 solar mass star is near 2 to 3 Lsun)
# and cools toward the base of the giant branch while crossing.
SUBGIANT_MAX_BRIGHTENING = 3.5
SUBGIANT_END_TEMPERATURE_K = 5000.0
# Massive stars evolve at nearly constant luminosity; about a third of
# observed supergiants are red (K/M), the rest blue or white (O/B/A).
SUPERGIANT_LUMINOSITY_GROWTH_RANGE = (1.0, 1.5)
RED_SUPERGIANT_FRACTION = 1 / 3
RED_SUPERGIANT_TEMPERATURE_RANGE_K = (3500.0, 4500.0)
BLUE_SUPERGIANT_TEMPERATURE_RANGE_K = (9000.0, 35000.0)
# Supergiant Yerkes class by luminosity (Lsun, lowest first): below the
# first threshold IB, then IAB, IA, and 0 (hypergiant) at the top.
SUPERGIANT_YERKES_THRESHOLDS_SOL = {"IAB": 10000.0, "IA": 50000.0, "0": 500000.0}

# Wolf-Rayet stars (GEN.243, docs/design/intergalactic-and-halo-contents.md 5.11): a very massive star
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
WR_PRESENT_MASS_FRACTION_RANGE = (0.4, 0.6)
"""tuple: The present mass of a Wolf-Rayet star over its initial mass, drawn uniformly [R]."""
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
WR_LUMINOSITY_RANGE_SOL = (1.0e5, 1.0e6)
"""tuple: The luminosity of a Wolf-Rayet star, solar luminosities, log-uniform [R]."""
WR_COMPANION_CHANCE = {"WN": 0.52, "WC": 0.90, "WO": 0.90}
"""dict: The share of Wolf-Rayet stars with a companion, by subtype (Dsilva et al., a spectroscopic multiplicity survey of
Galactic Wolf-Rayet stars: 52% of WN, 96% of WC, capped at 0.9 as in the phenomena catalogue; WO takes WC's)."""

# Piecewise main-sequence mass-luminosity relation, L = coeff * M^exponent
# (Lsun, Msun), each piece for M below its max_mass_sol and the last one
# open-ended (Duric 2004; Salaris & Cassisi 2005).
MS_MASS_LUMINOSITY_PIECES = (
    {"max_mass_sol": 0.43, "coeff": 0.23, "exponent": 2.3},
    {"max_mass_sol": 2.0, "coeff": 1.0, "exponent": 4.0},
    {"max_mass_sol": 55.0, "coeff": 1.4, "exponent": 3.5},
    {"max_mass_sol": None, "coeff": 32000.0, "exponent": 1.0},
)
# The "map luminosity" relation (GEN.200): the luminosity a main-sequence star
# of a given mass would have, used to rank a faint but massive object (a
# neutron star, a black hole) for the Galaxy Map by the star of its mass
# (Astronomy Notes, "Mass-Luminosity Relation Explained"; Varsity Tutors,
# "Stellar Mass & Lifetime"). Same form as MS_MASS_LUMINOSITY_PIECES, whose
# star generators it does not replace: it has Boss's table's constants and
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
# Main-sequence mass-radius relation, R = M^exponent (Rsun, Msun), with the
# exponent changing at 1 Msun (Demircan & Kahraman 1991).
MS_RADIUS_EXPONENT_BELOW_1_SOL = 0.8
MS_RADIUS_EXPONENT_ABOVE_1_SOL = 0.57
# The Sun's effective temperature (IAU 2015 nominal), anchoring
# Stefan-Boltzmann in solar units: T = T_sun * (L / R^2)^(1/4).
SUN_EFFECTIVE_TEMPERATURE_K = 5772.0

# White dwarfs: final mass from the initial-final mass relation of
# Kalirai et al. (2008), ApJ 676:594 (M_f = 0.109 M_i + 0.394), clamped to
# the observed range; luminosity from Mestel cooling,
# L = WD_COOLING_L0 * (M/0.6) * t_cool^-1.4 (Lsun, t_cool in Gy), clamped.
WD_IFMR_SLOPE = 0.109
WD_IFMR_INTERCEPT_SOL = 0.394
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
# reaches a few stellar radii), and a population-model white dwarf's
# progenitor, on its asymptotic giant branch, engulfed or drove off
# everything inside WD_PROGENITOR_ENGULFMENT_AU.
PLANET_MIN_STAR_AGE_GY = 0.01
LIFE_MIN_STAR_AGE_GY = 0.1
GIANT_ENGULFMENT_RADIUS_FACTOR = 2.0
WD_PROGENITOR_ENGULFMENT_AU = 1.5

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
# the light, so the time is stretched by 3.0, the value at which the model
# reaches the main sequence at the age the table's radius settles to within
# 2% of its minimum, rms 0.106 dex over the 18 masses (0.7 to 1.4 solar masses
# within 0.06 dex, 0.1 to 0.2 up to 0.21 dex too slow). The model ends on the
# generator's own main sequence (`stellar_evolution`), which is 0.05 to 0.28
# dex brighter than the table's 1 Gyr star below 1.2 solar masses (the table
# is solar-composition at 1 Gyr, the generator's relation is anchored on the
# Sun), so a state against the table is that much too bright at the end of the
# contraction. At 1 to 10 Myr the luminosity is within 0.25 dex of the table up
# to 0.3 solar masses and up to 0.35 dex too bright from 0.5 to 1.4, where the
# derived Hayashi phase falls more slowly than the table's.
PRE_MAIN_SEQUENCE_MAX_MASS_SOL = 8.0
# No star contracts longer than this: the track takes the lightest star (0.08 solar masses) 0.94 Gyr. A shortcut for
# the old stars, which are nearly all of them; the test checks it against the track.
PRE_MAIN_SEQUENCE_LONGEST_YR = 1.0e9
PROTOSTAR_MAX_AGE_YR = 0.5e6
HAYASHI_TEMPERATURE_1_SOL_K = 4357.0
HAYASHI_TEMPERATURE_MASS_EXPONENT = 0.1815
BIRTHLINE_RADIUS_1_SOL_RSUN = 3.187
BIRTHLINE_RADIUS_MASS_EXPONENT = 0.4335
CONVECTIVE_BINDING_COEFFICIENT = 3.0 / 7.0
RADIATIVE_BINDING_COEFFICIENT = 3.0 / 4.0
PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR = 3.0
# The disk (VANDAM and the TODO text, GEN.242): Class 0 under 0.15 Myr, I under
# 0.5, II under 3, and none after; mass 0.01 to 0.3 of the star and outer
# radius 10 to 300 AU, each log-uniform; the inner edge is where a blackbody
# grain sublimates (DISK_DUST_SUBLIMATION_K).
DISK_CLASS_AGES_YR = (("0", 0.15e6), ("I", 0.5e6), ("II", 3.0e6))
DISK_MASS_FRACTION_RANGE = (0.01, 0.3)
DISK_OUTER_RADIUS_RANGE_AU = (10.0, 300.0)
DISK_DUST_SUBLIMATION_K = 1500.0
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
STELLAR_POPULATION_AGE_RANGES_GY = {
    "young": (0.0, 0.1), "intermediate": (0.1, 3.0), "old": (3.0, 10.0), "bulge": (8.0, 12.0),
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
    name: (OPEN_CLUSTER_AGE_EDGES_MYR[k] / 1000.0, OPEN_CLUSTER_AGE_EDGES_MYR[k + 1] / 1000.0)
    for k, name in enumerate(OPEN_CLUSTER_POPULATIONS)})
# The three thin-disc populations' scale heights over the old one's come from the galaxy's own dispersions
# (`galaxy.structure.population_height_ratios`, GEN.532), not from a table here.
STELLAR_POPULATION_ARM_AMPLITUDE = {"young": 0.9, "intermediate": 0.4, "old": 0.0}

# How much of the thin disk's recent star formation a place holds (GEN.134, `galaxy.star_formation`): derived from the
# galaxy's shape, with the Milky Way's profile only as a test target.
STAR_FORMATION_LAW_EXPONENT = 1.4
"""float: The Kennicutt-Schmidt exponent: the rate per area goes as the gas per area to this power (Kennicutt 1998:
1.4 +/- 0.15 over 61 spirals and 36 starbursts; other fits run from 1.3 to 1.65)."""
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
# docs/design/disk-star-formation-history.md): gas continuity per ring, the inflow from the halo's growth, the
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
WHITE_DWARF_REMNANT_MASS_SOL = {"intercept": 0.444, "slope": 0.084}
"""dict: A star of 1 to 7 solar masses leaves a white dwarf of `intercept + slope * m` (Weidemann 2000's initial-final masses
as Madau and Dickinson 2014, ARA&A 52, 415, arXiv 1403.0007, quote them, read)."""
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

# Metallicity ([Fe/H], dex; GEN.158, docs/design/globular-clusters.md 5). The
# disc populations follow a mild radial gradient of -0.06 dex per kpc [recalled,
# to verify; written as a fall, so positive] that is zero at the Sun's distance from the centre; each population
# has its own offset and scatter (the young disc is the most uniform, the bulge
# the widest). The halo is -1.6 inside 25 kpc and -2.2 outside, scatter 0.5
# (docs/design/intergalactic-and-halo-contents.md 5.3). A cluster carries its own
# value. `METALLICITY_RANGE_FE_H` clamps every draw.
METALLICITY_DISC_FALL_DEX_PER_KPC = 0.06
METALLICITY_BY_POPULATION = {
    "young": {"offset": 0.0, "sigma": 0.10, "gradient": True},
    "intermediate": {"offset": 0.0, "sigma": 0.15, "gradient": True},
    "old": {"offset": -0.2, "sigma": 0.25, "gradient": True},
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
PLANET_OCCURRENCE_REFERENCE_OBJECTS = 16
"""
int: The most objects a metal-rich system's slots hold in practice. The disk-physics ceiling of
`StarSystem._estimate_max_objects_from_disk_physics` is far above what the orbit spacing lets fit (60 to 200 for a G
dwarf; the most that were placed in 60 draws was 15), so the cull cuts the count by `planet_occurrence_factor` from this
number, not from the ceiling: `ceil(16 x factor)` objects at most, which is 2 at [Fe/H] -1.65 and 16 at the floor.
"""

# Stellar abundances (GEN.439, `physics/stellar_abundances.py`, docs/design/stellar-composition.md and part A of
# docs/design/body-composition-research.md). The reference scale is Asplund, Amarsi and Grevesse 2021 (arXiv:2105.01661, Table 2
# log eps on H = 12; Boss 2026-10-17 1:24 AM EDT), moved to proto-solar by A21's diffusion correction of +0.064 dex for the
# metals and +0.070 dex for helium; it is a scale and never a preset. Verification states are in `STELLAR_RELATION_STATES`.
STELLAR_SOLAR_LOG_EPS = {
    "C": 8.46, "N": 7.83, "O": 8.69, "Ne": 8.06, "Na": 6.22, "Mg": 7.55, "Al": 6.43, "Si": 7.51, "S": 7.12, "Ar": 6.38,
    "K": 5.07, "Ca": 6.30, "Fe": 7.46, "Ni": 6.20, "Eu": 0.52, "Th": 0.03, "U": -0.54}
"""dict: Element symbol to A21 photospheric log eps (H = 12); the 14 elements that make 99.7 percent of Z, and Eu, Th and U."""
STELLAR_PROTO_METAL_OFFSET_DEX = 0.064
STELLAR_PROTO_HELIUM_OFFSET_DEX = 0.070
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
STELLAR_HIGH_ALPHA_LOGISTIC = (-0.350, 0.132)
"""tuple: `(f0, w)`: P(high-alpha) = 1 / (1 + exp(([Fe/H] - f0) / w)) for a star old enough to be thick-disc (Bensby, Feltzing
and Oey 2014, fit to the chemically high-alpha fraction; their paper is not in the repository's papers, so partly verified)."""
STELLAR_FE_H_FIT_RANGE = {"thin": (-0.9, 0.5), "high-alpha": (-1.6, 0.1)}
"""dict: The [Fe/H] each branch's [X/Fe] fits were made over; outside it the fit is held at its edge (flagged)."""
# [X/Fe] means by branch, as ("linear", a, b): a + b F; ("min", c, b, k): c + b min(F - k, 0), the thin disc's plateau above k
# and rise below; ("max", p, s, k): p + s max(F - k, 0), the high-alpha plateau below k. Fitted to Bensby 2014 (Mg, Si, Ca,
# Al), Nissen 2014 NLTE (C, O) and Battistini and Bensby 2016 (Eu), part A 2.7.
STELLAR_XFE_FITS = {
    "thin": {
        "Mg": ("min", 0.006, -0.343, 0.03), "Si": ("min", 0.017, -0.166, -0.07), "Ca": ("min", -0.041, -0.163, 0.42),
        "O": ("linear", -0.011, -0.373), "C": ("linear", -0.030, -0.163), "Al": ("min", 0.032, -0.029, 0.06),
        "Eu": ("linear", 0.060, -0.404)},
    "high-alpha": {
        "Mg": ("max", 0.357, -0.373, -0.73), "Si": ("max", 0.253, -0.199, -0.99), "Ca": ("max", 0.282, -0.251, -0.99),
        "O": ("max", 0.606, -0.495, -0.91), "C": ("max", 0.191, -0.521, -0.10), "Al": ("max", 0.225, -0.298, -0.48),
        "Eu": ("max", 0.446, -0.290, -0.93)},
}
STELLAR_XFE_SIGMA = {
    "thin": {"C": 0.044, "O": 0.052, "Mg": 0.066, "Si": 0.042, "Eu": 0.077},
    "high-alpha": {"C": 0.050, "O": 0.065, "Mg": 0.073, "Si": 0.073, "Eu": 0.084}}
"""dict: Residual scatter (dex) of the fits above, by branch."""
STELLAR_XFE_CORRELATION = (
    (1.00, 0.56, 0.69, 0.71),
    (0.56, 1.00, 0.36, 0.37),
    (0.69, 0.36, 1.00, 0.77),
    (0.71, 0.37, 0.77, 1.00))
"""tuple: Residual correlation of C, O, Mg and Si (Brewer and Fischer 2016, 844 dwarfs), so Mg/Si, a ratio of two correlated
terms, has its narrow real scatter."""
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
STELLAR_U_TH_UNCERTAINTY = 0.06
"""float: The relative spread on the birth U/Th (0.572 for the solar scale): the 2003 Palme and O'Neill cross-check, decayed
back, sits 6 percent from it (stellar-composition.md step 10). Carried as the uncertainty, not applied."""
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
    "WC": {"C": (0.20, 0.55)}}
"""dict: Surface mass fractions by Wolf-Rayet subtype; a pair is a range. WC carbon is the refined 20 to 55 percent (C/He 0.1 to
0.5 by number), inside the wider 10 to 60 percent of the literature; the rest of each star is helium (and oxygen in WO), the
split open."""
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

# Globular clusters (GEN.159, docs/design/globular-clusters.md). A cluster is stored as a density component, not as
# stars: its profile is a King (1966) model of one of ten concentrations `KING_W0_RANGE` (inclusive), solved once
# (`physics/king.py`) on `KING_TABLE_POINTS` log-spaced radii and scaled by length. A cluster's systems number its
# mass over `CLUSTER_MEAN_SYSTEM_MASS_SOL` (a Kroupa mean of 0.30 plus remnants and binaries, recalled, good to about
# a factor 1.5). The tidal radius a cluster may not exceed is the Jacobi radius in an isothermal halo of circular
# speed `CLUSTER_HALO_CIRCULAR_KMS`.
KING_W0_RANGE = (3, 12)
KING_TABLE_POINTS = 2000
# Satellite galaxies (GEN.237, docs/design/intergalactic-and-halo-contents.md 5.5) are stored the same way: the Clouds
# and Sagittarius as an exponential profile (density falling as exp(-r / a)), the dwarf spheroidals as a Plummer sphere
# (density as (1 + (r / a)^2)^-2.5), each cut off at a tidal radius that is `SATELLITE_TIDAL_OVER_HALF_MASS` times
# its half-mass radius (an exponential at 16 scale lengths holds all but 1e-5 of its mass, a Plummer sphere at
# 7.8 scale lengths 98%; the density is normalised so the cut profile holds the whole mass). The rows are real
# objects, but their sizes are recalled from McConnachie (2012) and not checked (the sources were out of reach).
SATELLITE_TIDAL_OVER_HALF_MASS = 6.0
CLUSTER_MEAN_SYSTEM_MASS_SOL = 0.4
CLUSTER_HALO_CIRCULAR_KMS = 220.0

# A generated galaxy's globular cluster system (`generation/clusters.py`, docs/design/globular-clusters.md section 8):
# the Milky Way's count; masses lognormal about `CLUSTER_MASS_PEAK_SOL` with `CLUSTER_MASS_SIGMA_DEX`, truncated to
# `CLUSTER_MASS_RANGE_SOL`; half-mass radii lognormal about `CLUSTER_HALF_MASS_RADIUS_PC` (width in dex, growing with
# distance from the centre as `(R / 8 kpc) ** CLUSTER_RADIUS_GROWTH`, limited to `CLUSTER_HALF_MASS_RADIUS_RANGE_PC`);
# a fraction `CLUSTER_CORE_COLLAPSED_FRACTION` is core-collapsed (W0 at the top of the range) and the rest draw their
# concentration c from a normal about `CLUSTER_CONCENTRATION_MEAN` (all recalled from Harris 1996 and Baumgardt and
# Hilker 2018, good to the shape of the distribution). [Fe/H] has two groups, `(fraction, mean, sigma)`: the
# metal-poor halo clusters and the metal-rich ones of the bulge and disk (the two peaks near -1.5 and -0.5), the
# metal-rich group sitting on flattened orbits about the centre (median `CLUSTER_BULGE_MEDIAN_PC`, rotating at
# `CLUSTER_HALO_CIRCULAR_KMS`) and the halo group on a density falling as r**-3.5 between `CLUSTER_HALO_RADIUS_RANGE_PC`.
CLUSTER_COUNT_MILKY_WAY = 157
CLUSTER_MASS_PEAK_SOL = 2.0e5
CLUSTER_MASS_SIGMA_DEX = 0.5
CLUSTER_MASS_RANGE_SOL = (1.0e4, 3.0e6)
CLUSTER_HALF_MASS_RADIUS_PC = 5.0
CLUSTER_HALF_MASS_RADIUS_SIGMA_DEX = 0.3
CLUSTER_HALF_MASS_RADIUS_RANGE_PC = (1.5, 40.0)
CLUSTER_RADIUS_GROWTH = 0.3
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
# docs/design/intergalactic-and-halo-contents.md 5.14). Every number is from a paper read for this item, none from the
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
# The Sun's distance from the center over the thin disk's scale length:
# R0 = 8.2 kpc, L = 2.6 kpc (BHG16). The calibration point (relative
# density 1.0) sits here by default, and the bar is angled from it.
GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH = 8200.0 / 2600.0
# The thick disk's mass, length and height, the bulge's axes and angle, and the stellar halo's share are drawn
# from the galaxy seed (`galaxy.structure_draw`, GEN.532, GEN.533); the Milky Way's published ones (BHG16: a thick
# disk of 900 pc and 2.0 kpc against the thin disk's 300 pc and 2.6 kpc; Dwek et al. 1995's G2 bulge fit with
# y0/x0 = 0.392 and z0/x0 = 0.272) are `galaxy.structure.REFERENCE_STRUCTURE`, a check target. The bulge is the bar's
# boxy/peanut core, as COBE/DIRBE sees it edge-on: rho = exp(-r_s^2 / 2) with r_s^4 =
# ((x/x0)^2 + (y/y0)^2)^2 + (z/z0)^4; `bulge_scale_radius_pc` is x0 (along the bar).

# The nuclear star cluster and the nuclear stellar disc (GEN.280, docs/design/realistic-galactic-core.md 3), the
# two structures at the real centre, are density terms of the model. Both are ratios to the galaxy's own bulge, so
# a toy galaxy gets a proportionate centre, and `GalaxyShape.nuclear_amplitude` (1 is the observed centre, 0 is
# none) scales the pair together. Calibrated on the default Milky Way-scale model, whose bulge holds 3.9e10 solar
# masses at 0.4 solar masses a system (the real bulge's 1.5e10 is in the model's units 2.6 times that, as every
# component's is; the totals are tuned to star counts, not to mass). The cluster's mass is drawn from the galaxy's
# stellar mass (`galaxy.structure_draw`, GEN.402: Neumayer 2020 and Georgiev 2016, or none); the disc is a multiple of it.
# The size ratios below are STAND-INS (gap plan 6) from the Milky Way's measured centre, kept until a size law is read.
NUCLEAR_CLUSTER_SCALE_RATIO = 2.07 / 1580.0
"""float: The cluster's Hernquist scale radius over the bulge scale radius (2.07 pc against 1,580 pc: a half-mass
radius of 5 pc, a projected half-light radius of 3.8 to 4.2 pc, against the observed 4.2)."""
NUCLEAR_DISC_CLUSTER_MASS_RATIO = 61.0
"""float: The nuclear disc's total mass over the cluster's (1.5e9 against 2.5e7). 65% of an exponential disc of
scale length 90 pc lies inside 200 pc, which makes the 1.0e9 inside 200 pc the observations give (recalled)."""
NUCLEAR_DISC_LENGTH_RATIO = 90.0 / 1580.0
NUCLEAR_DISC_HEIGHT_RATIO = 35.0 / 1580.0
"""float: The nuclear disc's exponential scale length and `sech^2` scale height over the bulge scale radius (90 and
35 pc; the observed height is 30 to 45 pc, recalled)."""
NUCLEAR_RESOLUTION_PC = 2.0
"""float: The radius of the grid's core ball (shell 0, half the 4 pc shell depth). The cluster's cusp is held below
the radius at which its density equals the ball's mean, so the ball's cell holds exactly the cluster's mass inside
2 pc and no cell is asked for a density the grid cannot hold (`galaxy.density.model_terms`)."""

# The galactic centre's gas and young stars (GEN.282, `galaxy.core_gas`, docs/design/realistic-galactic-core.md 3.4).
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
# The bar's long axis leads the Sun-center line by an angle drawn from the seed (`galaxy.structure_draw.draw_bar`);
# the Milky Way's is 27 degrees (its near end at positive galactic longitude; Wegg & Gerhard 2013, BHG16: 28-33).

# The stellar halo (GEN.234, docs/design/intergalactic-and-halo-contents.md 2.1 and 4) replaces the flat
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
HALO_LOCAL_DENSITY_FRACTION = 5.0e-3
"""float: The alternative normalisation the note keeps (Juric et al. 2008): the halo at the Sun is 0.5% of the
local density, about 2.7 times what `HALO_SYSTEM_FRACTION` gives (0.185%). Not used by the generator;
`halo_local_density_scale` is the factor that would switch to it."""

# The sector profile (GEN.235, docs/design/intergalactic-and-halo-contents.md 5.3): what a place sets for the
# rules that never change by place. A first model, not a fit, in the design note's order. The thin disc's, the
# bulge's and the halo's velocity dispersions are the galaxy's own (`galaxy.structure`, GEN.532); the rest are the
# middle of the note's ranges (a cluster 1 to 15 km/s inside and its orbit outside); `rotation` says what the
# population's net motion follows. [Recalled, to verify.]
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
# docs/design/globular-clusters.md), as `(share, mean, sigma)`; the shares are recalled (about two thirds are the
# metal-poor ones).
METALLICITY_CLUSTER_PEAKS = ((0.65, -1.5, 0.3), (0.35, -0.5, 0.3))
# Where a cell holds components (GEN.237) the peaks of its `cluster` and `satellite` populations are the components'
# own [Fe/H], each with this scatter (dex), weighted by the systems each puts in the cell; the peaks below stand in
# for a place no component reaches (a sample drawn without a profile). A satellite's stand-in is the Clouds'
# and dwarfs' middle (LMC about -0.4, SMC about -0.9, the dwarf spheroidals -1.0 to -2.1, recalled).
COMPONENT_FE_H_SIGMA = {"cluster": 0.15, "satellite": 0.3, "open_cluster": 0.05}
METALLICITY_SATELLITE_PEAKS = ((1.0, -1.2, 0.5),)

# A satellite's star formation history (GEN.526, `galaxy.satellite_history`, docs/design/intergalactic-and-halo-contents.md
# 5.5). Worked out from what its row holds (stellar mass, half-light radius, mean stellar age), with no value from the
# Clouds' own histories (those are the test targets of `tests/test_satellite_history.py`, never inputs).
SATELLITE_HISTORY_SPAN_GY = STELLAR_POPULATION_AGE_RANGES_GY["satellite"][1]
"""float: How long ago a satellite's first stars formed, Gy: the top of the age window, which is the age of the
Universe in this model. Its history runs from there to the present."""
STAR_FORMATION_LAW_COEFFICIENT = 2.5e-4
"""float: The Schmidt-Kennicutt law's coefficient: the star formation rate per area is this times the gas per area (solar
masses per pc2) to `STAR_FORMATION_LAW_EXPONENT`, in solar masses per year per kpc2 (Kennicutt 1998, ApJ 498, 541, eq. 4:
(2.5 +/- 0.7) e-4, read in arXiv astro-ph/9712213)."""
SATELLITE_STAR_FORMATION_THRESHOLD_MSUN_PC2 = 3.0
"""float: The gas surface density below which a satellite forms no stars, solar masses per pc2. Schaye (2004, ApJ 609,
667, read in arXiv astro-ph/0205125) puts the empirical threshold at 3 to 10 and the point where the gas turns cold
at about 10 (some H II regions form at 2 to 4); the lower end is taken [I], so a dwarf with a trace of gas still forms
stars."""

# How deep one dense cell is filled (GEN.281, docs/design/realistic-galactic-core.md 4.2). A cell that expects more
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
"""float: How much higher than the baseline the core's scatter limits default (docs/design/realistic-galactic-core.md
5.1, Boss): the core luminosity floor and mass limit follow the galaxy's baseline times this until they are set."""
HALO_MIN_LUMINOSITY_SOL = 500.0
"""float: The default halo luminosity floor, solar luminosities (GEN.283, Boss 1:09 AM EDT 2026-10-16)."""
HALO_MIN_LUMINOSITY_RANGE_SOL = (50.0, 4.0e6)
"""tuple: The lowest and highest halo luminosity floor `plan` accepts (Boss, 10:06 PM EDT 2026-10-15: lowest 50)."""
MEAN_SYSTEM_MASS_SOLAR = 0.4
"""float: The mass of an average system, remnants and companions included (docs/design/globular-clusters.md; recalled,
within a factor of 1.5): what the mass of the unresolved part of a cell is worked out from."""

# Gas and dark matter in a cell (docs/design/intergalactic-and-halo-contents.md 2.4 and 2.6). The ambient medium
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
# The NFW dark-matter halo of docs/design/galactic-potential.md (M_h = 5.4e11 solar masses, r_h = 16 kpc): 0.0091
# solar masses per pc3 at the Sun. Both scale with the disk scale length (a flat rotation curve: mass over radius
# is constant).
DARK_MATTER_SCALE_MASS_SOLAR = 5.4e11
DARK_MATTER_SCALE_RADIUS_PC = 16000.0

# Pre-placed bright stars (star_population.sample_bright_stars): the
# mass grid, in log-spaced cells over the IMF's range, on which the chance
# of a star being at least the threshold's brightness is tabulated.
BRIGHT_STAR_MASS_GRID_CELLS = 4000

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

# The chance of a companion inside 10 AU against the star's own [Fe/H] (GEN.103; Boss, 2026-10-16, 6:12 to 6:15 PM
# EDT: one pooled trend, every star mass). Moe, Kratter & Badenes 2019, ApJ 875:61 (arXiv 1808.02116), for solar-type
# primaries (0.6 to 1.5 Msun, periods under 10^4 d): 53 +-12% at [Fe/H] -3.0, 40 +-6% at -1.0, 24 +-4% at -0.2 and
# 10 +-3% at +0.5. `generation.place_terms.close_companion_factor` reads the straight line through them, held flat
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
# a dispersion from; it takes the Sun's neighbourhood, the way `constants.GALACTIC_CENTER_DISTANCE_LY` stands for its
# distance from the centre. Recalled round values (a mass density of 0.1 Msun/pc^3 and a field dispersion of
# 30 km/s; Jiang & Tremaine 2010 use 0.045 Msun/pc^3 for the observed stars and 29 km/s, their equation 47 and section 5),
# to verify.
SOLAR_NEIGHBOURHOOD_MASS_DENSITY_SOLAR_PC3 = 0.1
SOLAR_NEIGHBOURHOOD_DISPERSION_KMS = 30.0

# --- Planet Classification Data ---

# A dictionary defining the properties of different planet classes.
# Each class has a description, composition, radius range (in meters),
# habitable zone compatibility, atmosphere type, and planet type ('t' for terrestrial, 'g' for gas giant).
# Zone flags: "h" hot zone, "e" ecosphere, "c" cold zone, and "r" (GEN.8)
# whether a free-floating rogue planet may have the class -- "r" is the
# rogue planet's zone. With no star to warm it, a rogue is a dead (C),
# icy (D), rocky super-Earth (S) or giant (I, J, T) world: every class with life needs starlight
# (each one's life_chemical is a light-harvesting pigment), and the hot,
# molten and volcanic classes need a star's heat.
PLANET_CLASSES = {
    "A": {
        # Radius ceiling raised to 7500 to absorb a formerly-separate
        # "demon world" class's size range after its removal, with its
        # toxic/irradiated flavor folded in here as a variant rather than
        # kept as a separate class. See CHANGELOG.md for the full removal
        # rationale.
        "description": "a small, barren, and volcanic world, occasionally scarred by intense, toxic eruptions",
        "composition": "igneous silica and basalt",
        "radius_range": (500, 7500),
        # Mercury (2,439.7km) sits ~28% through this range.
        "size_mode": 0.28,
        "h": True, "e": False, "c": False, "r": False,
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
        # "with a thin atmosphere" moved off the description and onto the
        # atmosphere field itself ("a thin mix of...") -- the description
        # used to collide with the render template's own "with an
        # atmosphere of {atmosphere}" clause (generation/planet.py's
        # to_paragraph_list), producing "...world with a thin atmosphere
        # with an atmosphere of...".
        # Radius ceiling raised to 7500 to absorb the same formerly-separate
        # "demon world" class's range as Class A above, and a second,
        # formerly-separate "stripped core from a gas giant" class (no
        # atmosphere) folded in as an alternate origin story for the same
        # small-molten-world physical envelope, rather than kept as a
        # separate atmosphere-less class -- real Mercury-analog worlds this
        # close to their star already have only a negligible exosphere, so
        # B's existing "thin" atmosphere already covers that "no atmosphere"
        # identity closely enough. "and sulfur" added to the composition for
        # the shared volcanic/irradiated theme. See CHANGELOG.md for the
        # full removal rationale.
        "description": "a small, molten world, occasionally the stripped core of a former gas giant",
        "composition": "iron, potassium, silicon, and sulfur",
        "radius_range": (500, 7500),
        # Slightly below Class A's own Mercury anchor -- a freshly molten
        # or newly-stripped world skews a bit smaller.
        "size_mode": 0.25,
        "h": True, "e": False, "c": False, "r": False,
        "atmosphere": "a thin mix of helium, sodium, and oxygen",
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
        # Real small-body populations (asteroids, KBOs) follow a size-
        # frequency distribution strongly weighted toward smaller objects,
        # and this class's own broadest-of-any range skews the same way.
        "size_mode": 0.20,
        "h": True, "e": True, "c": True, "r": True,
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
        # Same real small-body size-frequency reasoning as Class C -- Ceres
        # (469.7km) sits near this range's own ceiling, but the population
        # as a whole skews toward its smaller end.
        "size_mode": 0.25,
        "h": True, "e": True, "c": True, "r": True,
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
        # "and a thin atmosphere" moved off the description and onto the
        # atmosphere field ("a thin mix of...") -- collided with the render
        # template's own "with an atmosphere of {atmosphere}" clause
        # otherwise (see Class B's note above).
        "description": "a world with a molten core and crust",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale range shared with F/G/M/N/O/P; no single real analog
        # for a young, molten-crust world, so anchored near Earth's own
        # (6,371km, ~27% through this range) like the rest of that group.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False, "r": False,
        # Hottest of the habitable classes (see below), so pushed toward
        # the inner (hotter) edge of the zone -- see Class K's own
        # `zone_position_mode` note.
        "zone_position_mode": 0.15,
        # "hydrogen compounds" made concrete as a real Hadean/Archean-analog
        # reducing mix (water vapor, ammonia, methane).
        "atmosphere": "a thin mix of water vapor, ammonia, and methane",
        "type": "t",
        # Hottest of the habitable (life-bearing) classes -- youngest,
        # volcanic, "barely supports life." Dark volcanic rock/minimal ice
        # keeps albedo low; the greenhouse_multiplier is the highest of the
        # E/F/G progression, reflecting real methane/ammonia's outsized
        # per-molecule greenhouse potency versus CO2. Verified via
        # climate_tuning_cli.py --class E: mean surface_temperature ~399K
        # over a 400-sample run (up from ~374K before `zone_position_mode`
        # above started placing E near the zone's hotter inner edge -- still
        # the top of the E->F->G cooling progression toward M/O/K/L/N
        # below).
        "albedo_range": (0.10, 0.18),
        "atm_molar_density_range": (0.0290, 0.0310),
        "atm_density_range": (0.3, 1.2),
        "greenhouse_multiplier_range": (4.5, 7.5),
        "life_chemical": ["Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.005, 0.015),
            "normal": (0.5, 1.5),
            "slow": (3.0, 8.0)
        }
    },
    "F": {
        "description": "a volcanic world with shallow seas",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False, "r": False,
        # Middle of the E->F->G progression positionally too -- see Class
        # K's own `zone_position_mode` note.
        "zone_position_mode": 0.30,
        "atmosphere": "a mix of carbon dioxide, ammonia, and methane",
        "type": "t",
        # Middle step of the E->F->G cooling progression: cooler than E
        # (higher albedo, lighter greenhouse_multiplier) but still hotter
        # than G/M. Verified via climate_tuning_cli.py --class F: mean
        # surface_temperature ~342K over a 300-sample run.
        "albedo_range": (0.15, 0.22),
        "atm_molar_density_range": (0.0295, 0.0315),
        "atm_density_range": (0.3, 1.0),
        "greenhouse_multiplier_range": (2.5, 4.0),
        "life_chemical": ["Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.005, 0.015),
            "normal": (0.5, 1.5),
            "slow": (3.0, 8.0)
        }
    },
    "G": {
        "description": "a rocky, barren world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False, "r": False,
        # Converges near M's own position too -- see Class K's own
        # `zone_position_mode` note.
        "zone_position_mode": 0.44,
        "atmosphere": "a mix of carbon dioxide, oxygen, and nitrogen",
        "type": "t",
        # Final step of the E->F->G cooling progression, converging near
        # M/O's own Earth-like range (the point of the progression --
        # "moving toward an M, O, K, L, or N"). Verified via
        # climate_tuning_cli.py --class G: mean surface_temperature ~297K
        # over a 300-sample run (vs Class M's ~286K, Class O's ~293K).
        "albedo_range": (0.20, 0.28),
        "atm_molar_density_range": (0.0300, 0.0320),
        "atm_density_range": (0.2, 0.8),
        "greenhouse_multiplier_range": (1.3, 2.0),
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.015, 0.05),
            "normal": (1.5, 4.0),
            "slow": (8.0, 25.0)
        }
    },
    "H": {
        "description": "a desert world with little liquid water",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False, "r": False,
        # Hot/dry -> pushed toward the zone's hotter inner half -- see
        # Class K's own `zone_position_mode` note.
        "zone_position_mode": 0.33,
        # "metals" replaced with "mineral dust" -- real deserts loft
        # particulate, not metal vapor (that's a magma-ocean/ultra-hot-rocky
        # -exoplanet phenomenon, not a match for a class that still has some
        # liquid water).
        "atmosphere": "a mix of oxygen, nitrogen, argon, and mineral dust",
        "type": "t",
        # Tuned hot/dry: lower albedo (dark exposed rock, minimal ice/cloud
        # cover) and a heavier, more CO2-loaded molar density than M/O drive
        # the heat; lower atm_density than M/O keeps it drier/thinner
        # (minimal water -> less retained humidity). Verified via
        # climate_tuning_cli.py --class H: mean surface_temperature ~335K
        # (vs Class M's ~286K), mean atmospheric_pressure ~38kPa (vs Class
        # M's ~99kPa) over a 400-sample run.
        "albedo_range": (0.18, 0.26),
        "atm_molar_density_range": (0.0325, 0.0345),
        "atm_density_range": (0.5, 0.9),
        "greenhouse_multiplier_range": (2.6, 3.2),
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
        # Uranus (25,362km) / Neptune (24,622km) sit ~28-30% through this
        # range.
        "size_mode": 0.28,
        # Added "e" (warm-Neptune analog -- real Neptune-mass planets in or
        # near a star's temperate zone are a common, well-documented
        # exoplanet category). Deliberately NOT given "h": real close-in,
        # Neptune-mass planets are rare -- the observed "hot Neptune desert"
        # -- because a star's X-ray/EUV irradiation photoevaporates a
        # Neptune-mass H/He envelope down to a bare rocky/metal core well
        # before it could stay class I; that outcome is already represented
        # by Class B's "occasionally the stripped core of a former gas
        # giant" (see Class B's note).
        "h": False, "e": True, "c": True, "r": True,
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
        # Jupiter (69,911km) / Saturn (58,232km) sit ~15-20% through this
        # range.
        "size_mode": 0.18,
        # Added "h" and "e": real Jupiter/Saturn-mass gas giants are
        # routinely found close to their star ("hot Jupiters", orbital
        # period < 10 days) and at intermediate distances ("warm Jupiters",
        # 10-365 days, sometimes within or near the habitable zone) as well
        # as at Jupiter/Saturn-like wide separations ("cold Jupiters") --
        # this is a standard three-way real observational classification,
        # not a stretch. A warm/cold Jupiter placed in zone 'e' can also
        # generate ordinary terrestrial moons via the existing moon-
        # generation path, including habitable-class ones -- the
        # "habitable exomoon around a giant planet" trope.
        "h": True, "e": True, "c": True, "r": True,
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
        # "with a thin atmosphere" dropped -- already collided with the
        # render template's own "with an atmosphere of {atmosphere}" clause,
        # and the atmosphere field below already says "a thin mix of...".
        "description": "an adaptable world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (2500, 7500),
        # Mars (3,389.5km) sits ~18% through this range.
        "size_mode": 0.18,
        "h": False, "e": True, "c": False, "r": False,
        # How far through the ecosphere zone's own [inner, outer] AU range
        # this class is generated, via a `random.sample_bounded_bell` draw
        # centered here instead of the zone's full width being equally
        # likely (see `physics.planets.generate_planet_properties`'s
        # "zone_position_mode" handling) -- Mars sits much farther from the
        # Sun than Earth, so K is pushed toward the outer (colder) edge of
        # the zone rather than sharing Class M's own ~0.5 (center) position.
        "zone_position_mode": 0.90,
        "atmosphere": "a thin mix of carbon dioxide and nitrogen",
        "type": "t",
        # Mars analog. Real Mars and Venus have almost identical mean
        # atmospheric molar mass (~43.3 vs 43.45 g/mol) but differ by ~100x
        # in greenhouse forcing -- composition (atm_molar_density) alone
        # can't tell them apart, so K keeps a realistically heavy/CO2-like
        # molar density (like N/Venus below) but gets a tiny
        # greenhouse_multiplier instead of N's huge one: same composition,
        # utterly different quantity/potency. Low atm_density keeps it
        # genuinely thin (Mars' real ~0.020 kg/m^3). Retuned after
        # `zone_position_mode` above started actually placing K near the
        # outer (colder) edge of the zone instead of sharing Class M's
        # midpoint position -- see docs/TODO.md's now-resolved "Open items"
        # entry for the history of the old, distance-blind version of this
        # class (mean surface_temperature ~231K/+9.9%, mean
        # atmospheric_pressure ~540Pa/-11.6% vs real Mars, "as close as
        # achievable without a zone change"). Verified via
        # climate_tuning_cli.py --class K: mean surface_temperature ~214K
        # (real Mars ~210K, +1.9%), mean atmospheric_pressure ~611Pa (real
        # Mars ~610Pa, +0.2%) over a 1000-sample run -- the zone change.
        "albedo_range": (0.36, 0.44),
        "atm_molar_density_range": (0.0420, 0.0433),
        "atm_density_range": (0.022, 0.042),
        "greenhouse_multiplier_range": (0.02, 0.05),
        "life_chemical": ["Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.015, 0.05),
            "normal": (1.5, 4.0),
            "slow": (8.0, 25.0)
        }
    },
    "L": {
        "description": "a marginally habitable world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 7500),
        # No single real analog; skewed toward the lower-mid range --
        # smaller worlds retain a thinner, more "marginal" atmosphere more
        # easily, consistent with this class's own tuning (see below).
        "size_mode": 0.35,
        "h": False, "e": True, "c": False, "r": False,
        # Between M and K positionally too, closer to K -- see Class K's
        # own `zone_position_mode` note.
        "zone_position_mode": 0.70,
        "atmosphere": "a mix of argon, oxygen, and trace elements",
        "type": "t",
        # K + "usually has vegetation" -> a modestly thicker, warmer, more
        # retained atmosphere than K (still CO2/N2-leaning composition, but
        # a meaningfully larger greenhouse_multiplier and atm_density than
        # K's near-zero values) -- enough to support vegetation without
        # approaching M's Earth-like identity. Verified via
        # climate_tuning_cli.py --class L: mean surface_temperature ~248K
        # (vs Class K's ~214K), mean atmospheric_pressure ~2.0kPa (vs
        # Class K's ~611Pa) over a 300-sample run.
        "albedo_range": (0.24, 0.30),
        "atm_molar_density_range": (0.0400, 0.0430),
        "atm_density_range": (0.03, 0.08),
        "greenhouse_multiplier_range": (0.25, 0.45),
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "M": {
        "description": "a temperate terrestrial world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth (6,371km) sits ~27% through this range -- exact.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False, "r": False,
        # The reference/baseline position every other ecosphere class's own
        # `zone_position_mode` is described relative to (see Class K's
        # note) -- dead center of the zone, same as the exact midpoint this
        # class was already tuned against before `zone_position_mode`
        # existed, so this is a no-op for M specifically.
        "zone_position_mode": 0.50,
        "atmosphere": "a mix of oxygen, nitrogen, and argon",
        "type": "t",
        # Tuned to real Earth: albedo ~0.29-0.31 (Earth's own Bond albedo,
        # ~0.3), atm_molar_density ~0.0288-0.0292 kg/mol (real Earth air,
        # ~0.02897), atm_density 1.3-1.55 kg/m^3 (near/above Earth's real
        # ~1.225 -- the corner of the pressure formula that actually reaches
        # ~1 atm, see docs/analysis/habitability-atmosphere-sanity-review.md),
        # greenhouse_multiplier 1.65-1.85 (calibrates base_ratio -- the
        # composition-only greenhouse proxy -- up to Earth's real ~33K
        # greenhouse effect; see physics.planets.calculate_atmospheric_conditions).
        # Verified via climate_tuning_cli.py --class M: mean surface_temperature
        # 286K (-0.6% vs 288K), mean atmospheric_pressure ~98.6kPa (-2.7% vs
        # 101,325 Pa) over a 400-sample run across the full host-star grid.
        "albedo_range": (0.29, 0.31),
        "atm_molar_density_range": (0.0288, 0.0292),
        "atm_density_range": (1.3, 1.55),
        "greenhouse_multiplier_range": (1.65, 1.85),
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "N": {
        # "dense, reducing atmosphere" moved onto the atmosphere field
        # itself -- collided with the render template's own "with an
        # atmosphere of {atmosphere}" clause otherwise.
        "description": "a hot world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Venus (6,051.8km) sits ~21% through this range -- exact.
        "size_mode": 0.21,
        "h": False, "e": True, "c": False, "r": False,
        # See Class K's own `zone_position_mode` note. Venus sits much
        # closer to the Sun than Earth, so N is pushed toward the inner
        # (hotter) edge of the ecosphere zone.
        "zone_position_mode": 0.05,
        "atmosphere": "a dense, reducing mix of carbon dioxide and sulfides",
        "type": "t",
        # Venus analog -- tuned to Venus's real surface temperature (737K)
        # and pressure (~9.2MPa). atm_molar_density near Venus's real
        # ~0.04345 kg/mol (near-pure CO2); albedo 0.75-0.90 matches Venus's
        # real highly-reflective cloud deck. Retuned down after
        # `zone_position_mode` above started actually placing N near the
        # inner (hotter) edge of the zone instead of sharing Class M's
        # midpoint position: the old greenhouse_multiplier (370-420) was
        # far above Venus's own real ratio (~101) specifically to compensate
        # for the wrong, too-cold midpoint distance (see docs/TODO.md's
        # now-resolved "Open items" entry for that history); the real
        # distance now does most of the work, so a much smaller multiplier
        # reaches the same target. atm_density is still well above Venus's
        # real ~65 kg/m^3 surface air density: this model's scale-height
        # formula (calculate_atmospheric_conditions) uses the pre-greenhouse
        # airless-equilibrium temperature rather than the final,
        # greenhouse-boosted surface temperature, which understates scale
        # height (and therefore pressure, P = density * g * H) regardless of
        # orbital distance -- compensated for here via atm_density rather
        # than by changing the shared scale-height formula, which affects
        # every class. Verified via climate_tuning_cli.py --class N: mean
        # surface_temperature ~737K (+0.0% vs real Venus), mean
        # atmospheric_pressure ~9.17MPa (-0.3% vs real Venus) over a
        # 1000-sample run.
        "albedo_range": (0.75, 0.90),
        "atm_molar_density_range": (0.0433, 0.0435),
        "atm_density_range": (300, 350),
        "greenhouse_multiplier_range": (260, 295),
        # No life_chemical (GEN.147): at Venus's 737 K and 9 MPa nothing lives on the
        # surface, and the class is not on HABITABLE_PLANET_CLASSES.
        "age_ranges": {
            "fast": (0.005, 0.015),
            "normal": (0.5, 1.5),
            "slow": (3.0, 8.0)
        }
    },
    "O": {
        "description": "an ocean world with most of its surface covered in liquid water",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False, "r": False,
        # Slightly hotter than M positionally too (see below) -- see Class
        # K's own `zone_position_mode` note.
        "zone_position_mode": 0.42,
        # Text now distinct from Class M's identical-before-this "oxygen,
        # nitrogen, and argon" -- water vapor is a real, tracked constituent
        # here (and is *lighter* than N2/O2, hence O's atm_molar_density_range
        # sitting below M's, not above -- see the tuning note below).
        "atmosphere": "a humid mix of oxygen, nitrogen, and water vapor",
        "type": "t",
        # Tuned warm/wet: a higher albedo than M (more cloud cover over an
        # ocean-dominated surface -- a real waterworld-climate-literature
        # finding) is more than offset by a stronger greenhouse_multiplier
        # (water vapor's real greenhouse contribution), netting *warmer*
        # than M despite the higher albedo. atm_molar_density is lighter
        # than M's (water vapor's molar mass, 18g/mol, is below N2/O2's) --
        # physically correct even though counterintuitive. Verified via
        # climate_tuning_cli.py --class O: mean surface_temperature ~298K
        # (vs Class M's ~286K), mean atmospheric_pressure ~83kPa over a
        # 400-sample run.
        "albedo_range": (0.28, 0.35),
        "atm_molar_density_range": (0.0270, 0.0285),
        "atm_density_range": (1.2, 1.5),
        "greenhouse_multiplier_range": (2.1, 2.4),
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
        # No single real analog; skewed slightly above the Earth anchor --
        # icy worlds carry proportionally more low-density ice content,
        # plausibly running a bit larger for a given mass than a pure
        # rock/iron world would.
        "size_mode": 0.30,
        "h": False, "e": True, "c": False, "r": False,
        # Cold/glaciated -> pushed toward the zone's colder outer edge, same
        # as K -- see Class K's own `zone_position_mode` note.
        "zone_position_mode": 0.90,
        "atmosphere": "a mix of oxygen, nitrogen, and argon (thinning with age)",
        "type": "t",
        # Icy/glaciated surfaces reflect far more sunlight than the default
        # rocky/Earth-like range (0.12, 0.35) -- real ice/snow Bond albedo is
        # roughly 0.5-0.9 (e.g. Europa ~0.68, Enceladus ~0.81). Gives P a
        # genuine cold bias from the unclamped physics instead of relying on
        # a post-hoc temperature clamp (see physics.planets.calculate_atmospheric_conditions).
        "albedo_range": (0.5, 0.7),
        # Given the same atm_molar_density_range/atm_density_range/
        # greenhouse_multiplier_range treatment the M/O/H/K/L/N/E/F/G/V pass
        # (CHANGELOG.md [5.3.7]) gave the other habitable classes -- P had
        # already stopped there with just its own albedo_range, per that
        # pass's note. No single real analog, so this isn't chasing a target
        # delta the way M/K/N are; instead the composition/quantity/potency
        # levers are set to make the class's own "cold, glaciated"/"thinning
        # with age" flavor text physically real rather than incidental.
        # atm_molar_density stays near Earth's own real ~0.02897 kg/mol
        # (the class's atmosphere text names oxygen/nitrogen/argon, not a
        # heavier CO2-like mix the way K/N's do) -- composition isn't what
        # makes this class cold. atm_density is set thin (between Mars'
        # ~0.02 and Earth's ~1.2 kg/m^3) for "thinning with age", and
        # greenhouse_multiplier is set weak (well below M's own 1.65-1.85)
        # so the cold comes from genuine physics, on top of the high albedo
        # above, rather than albedo alone. Verified via
        # climate_tuning_cli.py --class P: mean surface_temperature ~204K
        # (well below freezing, clearly colder than Class M's ~286K -- see
        # test_class_p_is_colder_on_average_than_class_m), mean
        # atmospheric_pressure ~9.5kPa (~0.1 atm, genuinely thin) over a
        # 400-sample run across the full host-star grid.
        "atm_molar_density_range": (0.0285, 0.0300),
        "atm_density_range": (0.05, 0.35),
        "greenhouse_multiplier_range": (0.3, 0.6),
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
        # No single real analog; moderate lower-mid skew.
        "size_mode": 0.35,
        # Was h/e/c all True -- but Q carries a life_chemical (it's a
        # life-bearing class), and habitable/life-bearing classes are
        # restricted to the ecosphere zone only (every other
        # life_chemical-bearing class -- E/F/G/H/K/L/M/N/O/P/V -- is
        # already e-only; see test_life_bearing_classes_are_ecosphere_only in
        # test_planets.py). Its "eccentric orbit" flavor still
        # holds fully confined to zone e -- a highly eccentric orbit *within*
        # the habitable zone still swings meaningfully between its own
        # perihelion and aphelion.
        "h": False, "e": True, "c": False, "r": False,
        "atmosphere": "a variable mix (thin to dense) of nitrogen, oxygen, and argon",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
        }
    },
    "S": {
        # Rocky super-Earth (GEN.28's first class, built with GEN.38): 1.2-1.8
        # Earth radii and 2-10 Earth masses, barren, in every zone. V stays
        # the life-bearing super-Earth. Compressed rock is denser than
        # Earth's own 5.51 g/cm3, so the density range is set so the class's
        # masses run 2-10 Earth masses across its radii (1.2 Earth radii at
        # 6.38 g/cm3 is 2 Earth masses, 1.8 at 9.45 is 10). Rogue-eligible:
        # without it, rocky rogues past Class C's 10,000 km ceiling had no
        # class that fit their size (GEN.38).
        "description": "a barren, rocky Super-Earth, its surface molten when it orbits close to its star",
        "composition": "silicate rock over a large iron core",
        "radius_range": (7645, 11468),
        # Smaller super-Earths are commoner (Fulton et al. 2017's radius
        # distribution falls from 1.2 toward 1.8 Earth radii).
        "size_mode": 0.30,
        "density_range": (6.38, 9.45),
        "h": True, "e": True, "c": True, "r": True,
        "atmosphere": None,
        "type": "t",
        "life_chemical": None,
        "age_ranges": {
            "fast": (0.0, 0.1),
            "normal": (0.0, 10.0),
            "slow": (0.0, 100.0)
        }
    },
    "T": {
        # Radius corrected 250,000-25,000,000km -> 15,000-55,000km -- real
        # brown-dwarf/giant-planet physics (electron degeneracy pressure
        # keeps sub-stellar objects within ~15% of Jupiter's own radius,
        # ~69,911km, regardless of mass) ruled out the old range regardless.
        # Now spans ice-giant-to-Saturn scale (Neptune 24,622km, Saturn
        # 58,232km), staying meaningfully below Class J's Jupiter-and-up
        # range -- a genuine "dwarf" relative to it.
        "description": "a gas dwarf with a thick atmosphere",
        "composition": "hydrogen, helium, and hydrocarbons",
        "radius_range": (15000, 55000),
        # Neptune (24,622km) sits ~24% through this range.
        "size_mode": 0.24,
        "h": False, "e": False, "c": True, "r": True,
        "atmosphere": "a mix of hydrogen, helium, and hydrocarbons",
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
        # The real rocky-to-gaseous transition radius (where planets start
        # retaining a significant H/He envelope) is commonly cited around
        # 1.5-1.6 Earth radii (~9,500-10,200km) -- right at this range's own
        # floor, so peaked low to stay on the rocky side of that boundary.
        "size_mode": 0.25,
        "h": False, "e": True, "c": False, "r": False,
        # Hot mainly via its high gravity/thick atmosphere rather than a
        # real-analog orbital position (unlike N/K), but still one of the
        # hotter habitable classes -- see Class K's own `zone_position_mode`
        # note.
        "zone_position_mode": 0.25,
        # Resolves the composition fork research flagged (H/He-retained
        # sub-Neptune-like reading vs. CO2-retained reading) toward the
        # latter -- "thick atmosphere with high surface temperature and
        # pressure" needs the heavier, more CO2-leaning composition; trace
        # hydrogen/helium kept as a minor primordial remnant, not the bulk.
        "atmosphere": "a thick mix of carbon dioxide, oxygen, and trace hydrogen and helium",
        "type": "t",
        # High gravity (from the larger radius_range above, mean ~1.7g)
        # retains a genuinely thick, hot atmosphere: heavier molar density
        # and a stronger greenhouse_multiplier than M/O/H, plus an
        # atm_density range well above every other terrestrial class except
        # N. Verified via climate_tuning_cli.py --class V: mean
        # surface_temperature ~383K, mean atmospheric_pressure ~279kPa
        # (~2.8 atm) over a 300-sample run -- clearly hot and thick, short
        # of N/Venus's full extreme.
        "albedo_range": (0.20, 0.30),
        "atm_molar_density_range": (0.0380, 0.0420),
        "atm_density_range": (2.0, 5.0),
        "greenhouse_multiplier_range": (3.5, 6.0),
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Bacteriochlorophylls", "Zinc-Bacteriochlorophyll", "Retinal", "Melanin"],
        "age_ranges": {
            "fast": (0.015, 0.1),
            "normal": (1.5, 9.0),
            "slow": (8.0, 50.0)
        }
    },
}

# Probabilities for each planet class to be generated. Six classes were
# removed over time: two small hot-zone rocky variants merged into A/B, two
# brown-dwarf-like sub-stellar classes cut entirely once their radius ranges
# were corrected to real physics and turned out to be redundant
# near-duplicates of each other, Class R (never reachable -- h/e/c were all
# False) cut entirely rather than left as permanent dead weight, and Class W
# ("tidally locked world with extreme temperature variations") cut entirely
# -- its day/night-split identity can't be produced from a single global
# surface_temperature scalar under this generator's climate model (see
# docs/TODO.md's now-resolved "Investigate Further" entry) -- (see
# CHANGELOG.md). Each removed class's weight was folded into the class(es)
# that absorbed its concept rather than just dropped; R already carried a
# weight of 0.0000, so nothing needed redistributing, and these weights are
# just relative (this dict isn't required to sum to 1.0 -- see
# physics.planets._choose_weighted_planet_class), so dropping W's tiny 0.0001
# share needed no redistribution either.
PLANET_CLASS_PROBABILITIES = {
    'A': 0.1400, 'B': 0.0725, 'C': 0.2365, 'D': 0.0142, 'E': 0.0239, 'F': 0.0335,
    'G': 0.0432, 'H': 0.0915, 'I': 0.0722, 'J': 0.0531, 'K': 0.0142, 'L': 0.0335,
    'M': 0.1345, 'N': 0.0239, 'O': 0.0045, 'P': 0.0046, 'Q': 0.0001,
    'S': 0.0300, 'T': 0.0001, 'V': 0.0045
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

# Main-sequence lifespan ranges, derived from t = 10 Gy * (M/Msun)^-2.5 (the
# standard mass-luminosity-based lifetime scaling, anchored to the Sun at 1
# Msun/10 Gy) applied to each spectral class's real mass range (O: >=16 Msun,
# B: 2.1-16, A: 1.4-2.1, F: 1.04-1.4, G: 0.8-1.04, K: 0.45-0.8, M: 0.08-0.45
# Msun -- the same boundaries implied by SPECTRAL_LUMINOSITY_RANGES). Each
# class's range now meets its neighbors' at the shared mass boundary instead
# of jumping (the previous table, e.g., capped A at 1.0 Gy while F started at
# 2.0 Gy, silently skipping the ~1.5-2.3 Gy lifespans of real late-A stars).
# O keeps a literature-cited floor (~1-10 Myr) rather than extrapolating the
# power law out to 150 Msun, which understates real mass-loss effects at
# extreme masses and would predict lifespans of only tens of thousands of
# years for the most massive O stars -- shorter than any observed star.
STAR_EVOLUTION = {
    "O": {
        "lifespan_gy": (0.001, 0.01), # 1 Million - 10 Million years
        "supported_evolutionary_scales": ["fast"],
        "potentially_viable_chemicals": ["Retinal"],
        "evolutionary_constraint_notes": "It burns out and explodes before planets can cool enough to hold oceans, so only the simplest light-driven chemistry could arise."
    },
    "B": {
        "lifespan_gy": (0.01, 1.5), # 10 Million - 1.5 Billion years
        "supported_evolutionary_scales": ["fast"],
        "potentially_viable_chemicals": ["Retinal"],
        "evolutionary_constraint_notes": "It is volatile and short-lived, which leaves barely enough time for a planet to cool and for single-celled life to appear."
    },
    "A": {
        "lifespan_gy": (1.5, 4.3), # 1.5 Billion - 4.3 Billion years
        "supported_evolutionary_scales": ["fast", "normal"],
        "potentially_viable_chemicals": ["Retinal", "Melanin"],
        "evolutionary_constraint_notes": "It lasts long enough for early biospheres and radiation-tolerant life, but probably not for slow, complex photosynthesis to mature."
    },
    "F": {
        "lifespan_gy": (4.3, 9.1), # 4.3 Billion - 9.1 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Blue-Optimized Porphyrins (+ GFPs)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "It lasts long enough for complex biospheres to emerge, including UV shielding and fluorescent adaptations."
    },
    "G": {
        "lifespan_gy": (9.1, 17.5), # 9.1 Billion - 17.5 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Chlorophyll a (Standard Porphyrins)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "It is steady for a long time, which suits the slow evolution of complex, oxygen-making photosynthesis."
    },
    "K": {
        "lifespan_gy": (17.5, 73.5), # 17.5 Billion - 73.5 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Chlorophyll a (Standard Porphyrins)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "It is very stable and lives for tens of billions of years, so slow-evolving biospheres can reach a climax and adapt to its slight red shift."
    },
    "M": {
        "lifespan_gy": (73.5, 5500.0), # 73.5 Billion - 5.5 Trillion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Bacteriochlorophylls (BChls)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "Early flares demand quick adaptation, such as radiotrophic melanin, but trillions of years of stability let deep-infrared photosynthesis last for good."
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
HABITABLE_PLANET_CLASSES = ['E', 'F', 'G', 'H', 'K', 'L', 'M', 'O', 'P', 'Q', 'V']
"""
list: A list of planet class codes that are considered habitable.
"""

PLANET_CLASS_MAX_LIFE_STAGE = {
    # Each habitable class's OWN "description" text above already commits
    # to a life-complexity ceiling -- `evolution.get_evolutionary_timeline`
    # used to have no idea what `planet_class` even was, so a Class G
    # planet (description: "simple life") could still roll all the way to
    # "Technological Civilization" if the star was old/fast-evolving
    # enough, contradicting its own generated text. Keys are
    # `EVOLUTIONARY_TIMELINES`/`EVOLUTIONARY_TEXT`'s own milestone keys;
    # a class not listed here (H, K, M, O, P, V) is uncapped -- its
    # description doesn't commit to any life-complexity ceiling at all
    # ("an adaptable world", "a temperate terrestrial world", etc.), so
    # the full milestone range -- including a technological civilization,
    # same as real Earth -- remains fair game.
    "E": "abiogenesis",         # "barely supports life" (this class's own
                                 # PLANET_CLASSES comment) -- the hottest,
                                 # youngest of the E->F->G progression;
                                 # capped at the most minimal stage,
                                 # "chemotrophs and extremophiles" only.
    "F": "photosynthesis",      # description: "...bacterial life" --
                                 # prokaryotic/microbial, matching
                                 # EVOLUTIONARY_TEXT['photosynthesis']'s
                                 # "simple, light-harvesting microbes...
                                 # microbial mats" (no nucleated cells yet).
    "G": "photosynthesis",      # description: "...simple life" -- same
                                 # prokaryotic tier as Class F.
    "Q": "photosynthesis",      # "extreme temperature variations" (GEN.147): an
                                 # eccentric orbit swings it between freezing and
                                 # scorching, so only microbial life survives.
    "L": "multicellularity",    # description: "...with vegetation" --
                                 # matches EVOLUTIONARY_TEXT
                                 # ['multicellularity']'s own "pioneering
                                 # flora...moss, ferns...vascular plants",
                                 # but not yet a sapient civilization.
}
"""
dict[str, str]: Maps a habitable planet class code to the highest
`EVOLUTIONARY_TIMELINES`-milestone key its own generated description text
is consistent with -- see `evolution.get_evolutionary_timeline`'s
`planet_class` parameter, the only reader of this mapping. A class not
present here is uncapped (full range, up to and including a technological
civilization).
"""

HABITABILITY_LIFE_SCORE_MIN = 0.005
"""
float: A planet's `phi_bio` (microbial) or `phi_cpx` (complex life) score
below this shows as 0.00 and counts as none (GEN.206): the habitability
index allows no life of that kind, so the description must not claim it.
"""

MOON_BLACKLIST = ['Q', 'S', 'V']
"""
list: A list of planet class codes that cannot be generated as moons.
"""

CO2_MAX_GREENHOUSE_FACTOR = 500
"""
int: A generous safety ceiling on greenhouse_factor (physics/planets.py's
calculate_atmospheric_conditions), not the per-class calibration knob --
that's PLANET_CLASSES[cls]["greenhouse_multiplier_range"]. Real Venus's own
airless-equilibrium-to-surface ratio is ~101, still well below N's (the
Venus analog) own tuned range even after that range was cut roughly in
half once N's `zone_position_mode` started placing it near the zone's
real, close-in inner edge instead of sharing Class M's midpoint position
(see PLANET_CLASSES["N"]'s tuning note) -- the atmospheric-pressure
scale-height understatement that note describes accounts for the rest of
the gap. 500 leaves generous headroom above N's own tuned range while
still guarding against a badly-configured future class producing a
runaway/non-finite temperature.
"""

FLAVOR_CHANCE_SYSTEM = 0.05 # The chance flavor text will be added to a system.
FLAVOR_CHANCE_PLANET = 0.05 # The chance flavor text will be added to a planet.
MAX_FLAVOR_TOTAL = 3 # The maximum flavor items total in a system.
MAX_FLAVOR_PLANET = 1 # The maximum flavor text for a single planet.
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
# docs/design/compact-remnant-regions.md). Estimates from the trends there.
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
# against the Milky Way's own stars: the Sun's radius and the disk's scale length.
PULSAR_RADIAL_A = 1.9
PULSAR_RADIAL_B = 5.0
GALAXY_SUN_RADIUS_PC = 8200.0
GALAXY_DISK_SCALE_LENGTH_PC = 2600.0

# Neutron star mass range in solar masses -- real measured masses cluster
# tightly around ~1.4 Msun with a hard floor near the Chandrasekhar-like
# collapse threshold and a soft ceiling near the maximum mass general
# relativity allows before further collapse to a black hole (Ozel &
# Freire 2016, ARA&A 54:401, compiling radio pulsar timing masses).
NEUTRON_STAR_MASS_RANGE_SOLAR = (1.1, 2.2)

# Neutron star radius (GEN.463, physics/neutron_star.py): the centre and the
# spread, km, at two masses, from Miller et al. 2021 (ApJL 918:L28,
# arXiv:2105.06979): R(1.4 Msun) = 12.45 +- 0.65 and R(2.08 Msun) = 12.35 +-
# 0.75, the full range of the 68% intervals over three equation-of-state
# frameworks; a star draws a normal deviation about the line between them,
# clipped at NEUTRON_STAR_RADIUS_CLIP_SIGMA deviations (the 95% band).
NEUTRON_STAR_RADIUS_ANCHORS = ((1.4, 12.45, 0.65), (2.08, 12.35, 0.75))
NEUTRON_STAR_RADIUS_CLIP_SIGMA = 2.0

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
# every component, 0.19% since GEN.469 gave the disk's populations their
# shares from its gas history (a galaxy of 5e10 solar masses; a lighter one forms its stars later and holds more young ones): test_gen_magnetar_scatter checks it); a galaxy
# with a different young share holds proportionally more or fewer.
MAGNETAR_PER_NEUTRON_STAR = 1.0e-7
MAGNETAR_YOUNG_STAR_SHARE = 0.0019
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
# docs/design/nebula-and-asteroid-field-classes.md). I and O are unused so
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
docs/design/nebula-density-vs-reality.md). The `chance` is that the star
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
# Mass bins (Boss's research, 2026-09-30; docs/design/
# interstellar-object-rates.md): a rogue planet draws a bin by its
# per-star rate, then a mass inside it. Since GEN.45 the rates and the draw
# inside a bin both follow one mass function, dN/dlogM ~ M^-0.65 from 0.1
# Earth masses to 13 Jupiter masses (Boss's research, 2026-10-01, which
# writes it dN/dM ~ M^-0.65; read per unit mass that would make about 87%
# of rogues gas giants, against its own aim of keeping them rare, so it is
# read per logarithm of mass). The total stays 6.5 per star. That gives:
#   - terrestrial 0.1-2 Earth masses, about 5.6 per star (2-10; Johnson
#     et al. 2020, AJ 160:123; Mroz et al. 2020, ApJL 903:L11);
#   - sub-Neptune / ice giant 2-20, about 0.72 per star (Sumi et al. 2023
#     find Neptune-mass candidates but give no rate);
#   - Saturn-class 20 Earth masses to 1 Jupiter mass, about 0.18;
#   - Jupiter-mass 1-13 Jupiter masses, about 0.028 per star (under Mroz
#     et al. 2017's upper limit of 0.25, Nature 548:183).
# So about 96% of rogues are terrestrial and 4% gas giants (past
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
    "jupiter": (317.8, 13 * 317.8),
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

ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER = (13.0, 80.0)
"""tuple: A free-floating brown dwarf, between the deuterium-burning
(~13 Mjup) and hydrogen-burning (~80 Mjup) limits. Stored as a
`rogue_planets` row with `mass_bin = 'brown-dwarf'` (schema v37)."""

ROGUE_PLANET_MASS_BIN_CHOICES = tuple(ROGUE_PLANET_MASS_BINS) + ("brown-dwarf",)
"""tuple: Every `rogue_planets.mass_bin` value."""

# Above this mass (in Jupiter masses, ~16 Earth masses), a generated rogue
# planet is treated as a gas giant rather than terrestrial -- roughly
# where the solar system's own ice giants (Uranus ~14.5, Neptune ~17
# Earth masses) sit, a common rough dividing line between rocky/icy
# terrestrial-scale bodies and volatile-dominated giants.
ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER = 0.05

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

ROGUE_RADIOGENIC_ISOTOPES = (
    # (isotope, present mantle concentration kg/kg, heat W/kg, half-life Gyr)
    ("U-238", 31e-9 * 0.9928, 9.46e-5, 4.468),
    ("U-235", 31e-9 * 0.0071, 5.69e-4, 0.704),
    ("Th-232", 124e-9, 2.64e-5, 14.05),
    ("K-40", 310e-6 * 1.19e-4, 2.92e-5, 1.248),
)
"""tuple: Earth's mantle heat sources today (Turcotte & Schubert 2014,
Geodynamics 3rd ed., tables 4.2-4.3: U 31 ppb, Th 124 ppb, K 310 ppm,
about 7.4e-12 W/kg). A rogue of age t carries Earth's concentrations
as they stood at age t (each scaled by exp(lambda * (4.5 Gyr - t)))."""

ROGUE_EARTH_REFERENCE_AGE_GY = 4.5
"""float: The age the concentrations above describe (Earth's)."""

ROGUE_MANTLE_MASS_FRACTION = 0.67
"""float: The heat-producing (silicate) share of a rocky rogue's mass,
Earth's (4.0e24 of 5.97e24 kg)."""

ROGUE_RADIOGENIC_ABUNDANCE_RANGE = (0.5, 2.0)
"""tuple: A rogue's radioactive abundance relative to Earth's, drawn
log-uniformly -- stars, and so their planets, differ in U, Th and K (a
default)."""

ROGUE_UREY_RATIO = 0.67
"""float: Radiogenic heat over total heat flow. The rest is leftover
formation heat (secular cooling), kept in Earth's proportion: Earth's
mantle radiogenic ~0.058 W/m^2 over 0.67 gives ~0.087 W/m^2, Earth's
measured flow and the 0.08-0.1 W/m^2 Boss's research quotes."""

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

ROGUE_MAX_EFFECTIVE_TEMPERATURE_K = 2800.0
"""float: The hottest a rogue gets. The cooling law above assumes an old,
shrunken body; a brown dwarf younger than a few hundred million years is
still puffed up, so the law (with today's radius) overshoots. Young brown
dwarfs near the hydrogen-burning limit are late-M dwarfs, ~2,800 K."""

# Brown dwarfs as stars (GEN.252, GEN.549, docs/design/brown-dwarf-state.md): the cooling track, the mass function and the
# map from temperature to spectral type.
BROWN_DWARF_TRACKS = (
    (0.01, ((0.01, 1570, 0.156),)),
    (0.015, ((0.01, 2261, 0.189), (0.03, 2139, 0.171), (0.05, 1760, 0.146))),
    (0.02, ((0.01, 2519, 0.244), (0.03, 2184, 0.159), (0.05, 1838, 0.137))),
    (0.03, ((0.01, 2682, 0.283), (0.03, 2448, 0.167), (0.1, 1912, 0.124), (0.12, 1779, 0.12))),
    (0.04, ((0.01, 2768, 0.292), (0.03, 2632, 0.185), (0.1, 2311, 0.13), (0.3, 1667, 0.108))),
    (0.05, ((0.01, 2843, 0.315), (0.03, 2755, 0.2), (0.1, 2525, 0.139), (0.3, 2063, 0.109), (0.625, 1555, 0.099))),
    (0.06, ((0.01, 2900, 0.335), (0.03, 2860, 0.217), (0.1, 2670, 0.148), (0.3, 2331, 0.112), (1.0, 1584, 0.094))),
    (0.07, ((0.01, 2946, 0.358), (0.03, 2925, 0.233), (0.1, 2794, 0.158), (0.3, 2519, 0.117), (1.0, 2047, 0.094),
            (3.0, 1650, 0.09), (10.0, 1626, 0.09))),
    (0.072, ((0.01, 2954, 0.363), (0.03, 2936, 0.236), (0.1, 2813, 0.16), (0.3, 2552, 0.118), (1.0, 2109, 0.095),
             (3.0, 1793, 0.091), (10.0, 1790, 0.091))),
    (0.075, ((0.01, 2967, 0.369), (0.03, 2952, 0.241), (0.1, 2839, 0.162), (0.3, 2594, 0.119), (1.0, 2243, 0.097),
             (3.0, 2052, 0.092), (10.0, 2037, 0.091))),
    (0.08, ((0.01, 2987, 0.383), (0.03, 2977, 0.249), (0.1, 2879, 0.167), (0.3, 2660, 0.122), (1.0, 2400, 0.102),
            (3.0, 2343, 0.099), (10.0, 2345, 0.099))),
)
"""tuple: The cooling track of a brown dwarf, one row per mass: `(mass in solar masses, ((age in Gy, effective
temperature in K, radius in solar radii), ...))`, oldest last. These are values of the solar-metallicity models of
Baraffe, Homeier, Allard and Chabrier 2015 (A&A 577:A42, arXiv 1503.04107, read; the 2MASS isochrone file of
perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir, whose BT-Settl atmospheres reach down to about 1,500 K), thinned to nine
ages and rounded as printed. Their tables stop where an object has cooled below the models' range (the last age of a
row), which is why Phillips et al. 2020 (A&A 637:A38, arXiv 2003.13717, read) limit their own ATMO 2020 grid to below
about 2,000 K and send the hotter objects to these models. A row that ends before 10 Gy continues at
`BROWN_DWARF_TEMPERATURE_AGE_EXPONENT`; a row that reaches 10 Gy holds its last values (the near-limit objects settle on a
plateau)."""

BROWN_DWARF_TRACK_END_GY = 10.0
"""float: The age the tabulated tracks are cut at: a `BROWN_DWARF_TRACKS` row whose last age is this one has settled."""

BROWN_DWARF_TEMPERATURE_AGE_EXPONENT = -0.32
"""float: Burrows et al. 2001 (Rev. Mod. Phys. 73:719, eq. 2, read, arXiv astro-ph/0103383): the effective temperature of
an old brown dwarf falls as `t^-0.32`. It continues a `BROWN_DWARF_TRACKS` row past its last age."""

BROWN_DWARF_MASS_FUNCTION_SLOPE = 0.6
"""float: The substellar mass function `dN/dM ~ M^-0.6` (Kirkpatrick et al. 2021, ApJS 253:7, arXiv 2011.11616, read:
alpha = 0.6 +- 0.1 from the 20 pc census of 525 L, T and Y dwarfs). It is drawn between
`ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER`; the generator's own stellar mass function is not continued below
`IMF_BREAKS_SOL[0]`, since continuing it would give 0.7 brown dwarfs per star against the census's 0.2 (the
free-floating rate `PHENOMENON_DENSITY_PC3["brown-dwarf"]`)."""

BROWN_DWARF_SPECTRAL_TYPE_FITS = (
    (0.0, 8.75, (2237.5, -144.96, 4.0301), 134.0),
    (8.75, 14.75, (1437.9, -18.309, 0.0), 79.0),
    (14.75, 22.0, (5141.3, -368.65, 6.7301), 79.0),
)
"""tuple: Kirkpatrick et al. 2021 Table 13 (read): `Teff` (K) as a polynomial in the spectral subclass number
(0 = L0, 10 = T0, 20 = Y0), per range of that number: `(low, high, coefficients c0 + c1 x + c2 x^2, rms K)` (T0 to T6 is a straight line, so its c2 is 0). The
scatter is the fit's, the measurement's, and is not drawn. Y0 falls at 461 K, inside the 371 to 510 K the paper
measured."""

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

COMET_COMPOSITION = [
    "water ice", "carbon dioxide ice", "carbon monoxide ice",
    "methane ice", "ammonia ice", "silicate dust", "amorphous carbon",
    "complex organic compounds", "hydrogen cyanide", "formaldehyde",
]
"""
list: Common cometary nucleus components (real, spectroscopically-
identified cometary ices/dust, e.g. as surveyed in A'Hearn et al. 1995,
Icarus 118:223), used the same way ASTEROID_COMPONENTS is -- a random
subset sampled per object for descriptive composition text. Shared by
`InterstellarComet` and `cometData.Comet` alike -- nucleus chemistry
doesn't depend on whether the comet is bound to a star.
"""

# --- Star-Bound Comets (cometData.Comet) ---
#
# Unlike InterstellarComet (unbound, hyperbolic, encountered once, always
# standalone), a Comet here is bound to a star's own system, propagated
# via real two-body Kepler/Barker orbital mechanics (see keplerMotion.py)
# rather than a fixed hyperbolic excess speed. See
# docs/design/comet-orbital-realism.md for the full research/design
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
COMET_NUCLEUS_DENSITY_KG_M3 = 533.0
"""float: Bulk density of a comet nucleus, kg/m^3 -- 67P/Churyumov-Gerasimenko, measured by Rosetta (Paetzold et al.
2016, 533 +/- 6), the one nucleus whose mass and volume are both known. A nucleus's mass is this times its sphere's
volume (`physics/bulk_mass.comet_nucleus_mass_kg`)."""

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
    # Real dynamical comet families (see docs/design/comet-orbital-realism.md's
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
# (periodic, returns every orbit) -- see docs/design/comet-orbital-realism.md's
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
# GEN.449 changed the companions). A sector or galaxy run's
# `--prevalence feature=+50` forces the feature on in enough systems to
# make that share 1.5 times as large. Remeasure after a change to
# generation moves one of them: tests/test_prevalence.py checks they still
# hold within a few points.
PREVALENCE_BASE_SHARES = {
    "habitable_world": 0.257,   # a habitable planet or moon
    "asteroid_belt": 0.607,     # an asteroid belt
    "large_star": 0.050,        # a primary born at LARGE_STAR_MIN_MASS_SOL or more
    "planets": 0.961,           # at least one planet or belt
}

# The usual share of each feature one draw decides (ADM.16: the Generate
# page starts each prevalence field at its usual share), measured over
# 10,000 systems with default options (four seeds of 2,500; 2026-10-08;
# comets, binary_system and wide_binary remeasured 2026-10-17 over 20,000
# systems, four seeds of 5,000: GEN.448 puts 93% of the binaries in the wide
# class, the close class being the 0.05 to 0.25 AU range only).
# Each is a share of the things in its comment, not always of systems.
# Remeasure with PREVALENCE_BASE_SHARES.
PREVALENCE_DRAW_SHARES = {
    "comets": 0.321,            # systems with comets (each star of a wide pair draws)
    "binary_system": 0.260,     # systems with two stars
    "wide_binary": 0.928,       # binaries that are wide pairs
    "moons": 0.429,             # planets with at least one moon
    "max_planets": 0.037,       # stars that get the most orbits they can hold
}

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
# `dense`), `components` (what its bodies are made of, sampled for the
# composition text), `description`, and `frequency` (relative weight).
ASTEROID_FIELD_COMPOSITIONS = {
    "carbonaceous": {
        "letters": {"sparse": "A", "typical": "B", "dense": "C"},
        "components": ["carbon", "serpentine", "magnetite", "silicon carbide", "sulfur",
                       "phosphorus", "olivine", "troilite"],
        "description": "dark carbon-rich bodies with hydrated clays",
        "frequency": 30.0,
    },
    "stony": {
        "letters": {"sparse": "D", "typical": "E", "dense": "F"},
        "components": ["olivine", "pyroxene", "plagioclase feldspars", "silicon", "magnesium",
                       "troilite", "iron", "silicon dioxide"],
        "description": "silicate rock with flecks of metal",
        "frequency": 20.0,
    },
    "metallic": {
        "letters": {"sparse": "G", "typical": "H", "dense": "J"},
        "components": ["iron", "nickel", "kamacite", "taenite", "schreibersite", "cohenite",
                       "iridium", "platinum", "osmiridium"],
        "description": "iron-nickel cores of shattered bodies",
        "frequency": 8.0,
    },
    "icy": {
        "letters": {"sparse": "K", "typical": "L", "dense": "M"},
        "components": ["water ice", "carbon dioxide ice", "ammonia ice", "amorphous carbon",
                       "complex organic compounds", "silicate dust", "serpentine"],
        "description": "volatile-rich primitive bodies, ices bound in dark crusts",
        "frequency": 20.0,
    },
    "basaltic": {
        "letters": {"sparse": "N", "typical": "N", "dense": "P"},
        "components": ["pyroxene", "plagioclase feldspars", "olivine", "ilmenite", "chromite",
                       "magnesium", "calcium"],
        "description": "basalt from a differentiated body's crust",
        "frequency": 4.0,
    },
    "mixed": {
        "letters": {"sparse": "Q", "typical": "R", "dense": "S"},
        "components": None,
        "description": "a mix of rock, metal and carbonaceous bodies",
        "frequency": 10.0,
    },
    "dust": {
        "letters": {"sparse": "T", "typical": "T", "dense": "T"},
        "components": ["silicon dioxide", "carbon", "olivine", "pyroxene", "magnetite"],
        "description": "mostly dust and gravel, with few large bodies",
        "frequency": 5.0,
    },
    "collisional": {
        "letters": {"sparse": "U", "typical": "U", "dense": "U"},
        "components": None,
        "description": "fragments of one parent body broken apart in a collision",
        "frequency": 3.0,
    },
}
"""dict: See the comment above. The icy family borrows cometary ices
(`COMET_COMPOSITION`); `components` None means sample from all of
`ASTEROID_COMPONENTS` (a mixed field, or a collisional family whose
parent body could have been anything)."""

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
# Rates come from Boss's research of 2026-09-30 (docs/design/
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
    # Sum of ROGUE_PLANET_MASS_BINS' per-star rates (6.5) at n_*0.
    # Terrestrial alone: about 0.78 pc^-3 (0.5-1.4); Jupiter-mass about
    # 0.004 (under Mroz et al. 2017's limit of 0.035).
    "rogue-planet": sum(rate for _lo, _hi, rate in ROGUE_PLANET_MASS_BINS.values()) * REFERENCE_STELLAR_DENSITY_PC3,
    # 0.025-0.035 (Kirkpatrick et al. 2021, ApJS 253:7), ~1 per 4.7 stars.
    "brown-dwarf": 0.03,
    # Stars moving > 30 km/s, 1-2% of stars (Tauris 2015, MNRAS 448:L6).
    # A flag on an ordinary generated system, not a phenomenon row.
    "runaway-star": 2.1e-3,
    # Stars moving > 500 km/s, ejected by the central black hole, at
    # HVS_REFERENCE_RADIUS_PC; scales as r_GC^-2 (Brown 2015, ARA&A
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
    # (`galaxy.black_hole_density`, GEN.563). The old 7e-5 (1e8 in the Galaxy, 0.05% of ~2.5e11 stars) is a test target.
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
in range. Study: docs/design/nebula-density-vs-reality.md."""

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

HVS_REFERENCE_RADIUS_PC = 8000.0
"""float: Galactic radius the "hypervelocity-star" density is quoted at."""

INTERSTELLAR_DEBRIS_DENSITY_PC3 = 1e12
"""float: Interstellar comets and planetesimals (m-km) per pc^3 at
REFERENCE_STELLAR_DENSITY_PC3, 1e11-1e14 (Engelhardt et al. 2017, AJ
153:133; Seligman & Laughlin 2020, ApJL 896:L8) -- far too many for rows,
so a sector shows it as a figure (`planetgen.db.query.interstellar_debris_count`)."""

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

NEUTRON_STAR_UNBOUND_FRACTION_RANGE = (0.2, 0.4)
"""tuple: The share of a galaxy's neutron stars a birth kick sends out of it, drawn once per galaxy from its seed
(GEN.238; 20% to 40% of the isolated ones: Hobbs et al. 2005's kick distribution against the escape speed)
[S, recalled; to verify]."""

BLACK_HOLE_UNBOUND_FRACTION = 0.02
"""float: The share of a galaxy's black holes that are ejected (GEN.238; about 2%) [S, recalled; to verify]."""

CORE_COLLAPSE_MASS_RANGE_SOL = (SUPERGIANT_MIN_MASS_SOL, IMF_BREAKS_SOL[-1])
"""tuple: The initial masses that end as a neutron star or a black hole, solar masses: 8 (the supergiants' edge) to the
top of the IMF, 150 (GEN.560). 0.637% of the stars born."""

BLACK_HOLE_FAILED_SHARE_SOLAR_RANGE = (0.0, 0.45)
"""tuple: The share of the massive stars of solar metallicity that fail to explode and leave a black hole: the widest
verified range, 0 (O'Connor and Ott 2011's lowest) to 45% (Sukhbold et al. 2016's highest; Ebinger et al. 2020 give 5 to
21%). A galaxy draws its own, uniform over the range (only the range is known) and cut by the neutron star budget
(`physics.black_hole_formation.share_cap`). GEN.563, docs/design/black-hole-formation-sources.md."""

BLACK_HOLE_FAILED_SHARE_METAL_FREE_RANGE = (0.15, 0.32)
"""tuple: The same share for a metal-free population, 15% (O'Connor and Ott 2011) to 32% (Ebinger et al. 2020's
highest; Sukhbold et al. give solar only). Drawn uniform over the range, cut at the neutron star budget and not below
the galaxy's solar share. GEN.563."""

BLACK_HOLE_FORMED_MASS_STEPS = 64
BLACK_HOLE_FORMED_AGE_STEPS = 64
"""int: Log-spaced progenitor masses and ages across a window the share of formed remnants is averaged over (numerical,
`physics.black_hole_formation.remnant_formed_share`)."""

BLACK_HOLE_GALAXY_MEAN_RADIAL_STEPS = 72
BLACK_HOLE_GALAXY_MEAN_POLAR_STEPS = 24
BLACK_HOLE_GALAXY_MEAN_AZIMUTH_STEPS = 6
"""int: The log-spaced radii, equal-area polar angles and azimuths of the grid a galaxy's mean black hole rate per star
is averaged over (`galaxy.black_hole_density.mean_rate_per_star`, GEN.563)."""

EJECTA_SPEED_RANGE_KMS = (300.0, 500.0)
"""tuple: The speed a kicked neutron star or black hole has far out, above the galaxy's escape speed, km/s (GEN.238)
[I]. The population's density takes the middle of it; each object draws its own, uniform."""

EJECTA_AGE_GY = 10.0
"""float: How long the ejecta of a galaxy have been leaving, Gyr: with `EJECTA_SPEED_RANGE_KMS` it sets where the
population ends, at `v T` (3 to 5 Mpc, far beyond any grid's reach) (GEN.238) [I]."""

EJECTA_MAX_PER_CELL = 50.0
"""float: The most ejected objects of one kind a cell draws (GEN.238), a program limit: the draw is O(mean), and a galaxy
shaped to hold astronomically many systems in a few cells would otherwise draw more than a cell has room for. The real
galaxy's densest cells draw far fewer."""

HYPERVELOCITY_COUNT_RADIUS_KPC = 100.0
"""float: The radius the count `HYPERVELOCITY_STARS_PER_GALAXY` is read to: "galaxy-wide" is taken as out to 100 kpc,
the halo's extent (docs/design/intergalactic-and-halo-contents.md 2.3) [I]."""


PHENOMENON_MIN_MASS_PRESETS = (8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0)
"""tuple: The mass limits the user may pick (GEN.183, Boss 2026-10-10): the
Generate page's slider and `--phenomenon-min-mass` accept only these."""

COMPACT_MIN_MASS_PRESETS = (1.0, 2.0, 4.0, 6.0)
"""tuple: The compact-object limits the user may pick (GEN.195, Boss 2026-10-10): the
lowest mass of a neutron star or black hole the scatter places. `--compact-min-mass`
and the Generate page's compact-object control accept only these, or "star" to use the
stellar mass limit."""

REDO_SCATTERS = ("mass", "luminosity", "phenomena")
REDO_NEBULA_FIELD_MODES = ("dry-run", "apply")
"""tuple: The scatters `planetgen plan --redo-scatters` can redo (GEN.196), in the order they are listed."""

COMPACT_MIN_MASS_STAR = "star"
"""str: `--compact-min-mass`'s word for "use the stellar mass limit" (the default)."""

PHENOMENON_MIN_MASS_SOLAR = 14.0
"""float: The phenomenon scatter's lowest mass, the default of
`PHENOMENON_MIN_MASS_PRESETS` (GEN.167, `planetgen plan
--phenomenon-min-mass`; docs/design/phenomenon-scatter-mass-cut.md, Boss
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
than occupied (see docs/design/galaxy-coordinate-system.md section 9's
storage-analysis addendum: ~320 billion addressable sector slots total).
"""

# --- Navigation Parameters ---

WARP_FACTORS_FOR_NAV = (1, 2, 4, 8, 9, 9.5, 9.9, 9.995)
"""
The warp factors NAV output reports travel time at (see
`planetgen.galaxy.navigation.warp_travel_times`) -- the rows of Boss's warp
table (docs/design/navigation-frames.md, "Travel speeds"), from warp 1
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

NAV_COURSE_DECIMAL_PLACES = 2
"""
Decimal places used when formatting NAV distance/azimuth/altitude output.
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
# See docs/design/population-and-politics.md.

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

SPECIES_EQUIPMENT_MIN_TECH = (0.0, 2.0, 2.5, 3.0, 4.0)
"""tuple: The least tech level (0 to 7) a species needs to make the kit of each equipment tier
(`habitability_world.EQUIPMENT_NAMES`: ideal, breathing mask, mask with scrubber, environment suit, full life
support with radiation hardening), for the "could live there with equipment" list. A planetGen default."""

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

SLOW_PAGE_ESTIMATE_SECONDS = 60
"""int: How long a web process keeps the table row estimates before it asks the API again. A planetGen default."""

SLOW_PAGE_TIMING_SECONDS = 7 * 86400
"""int: How long the build time of a request is remembered in Redis. A planetGen default."""
