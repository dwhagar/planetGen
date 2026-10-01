# stellarObjects/program_constants.py

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
"""

from . import physical_constants as _physical_constants

# --- Planet Generation Parameters ---
# Domingos, Winter & Yokoyama (2006), MNRAS 373:1227, "Stable satellites
# around extrasolar giant planets" -- a prograde moon on a circular orbit
# stays bound out to about 0.4895 of its planet's Hill radius; beyond that
# the star strips it (planetPhysics.generate_moons' outer limit).
MOON_PROGRADE_STABLE_HILL_FRACTION = 0.4895

# The average ratio of a gas giant's core mass to its total mass
GAS_GIANT_CORE_ATMOSPHERE_RATIO = (0.03, 0.6)

# --- Star System Generation Parameters ---
INITIAL_PLANET_DISTANCE_FACTOR = 0.55
ASTEROID_BELT_PROBABILITY = 0.1
ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN = 1.1
ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX = 2
ABSOLUTE_MAX_SYSTEM_OBJECTS = 500
MIN_ASTEROID_BELT_SEPARATION = 0.05

# Minimum stable orbital separation between two adjacent planets, in units
# of their *mutual* Hill radius (utils.mutual_hill_radius_m -- the pair's
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
# `Planet.min_orbit_distance`/`planetPhysics.generate_moons`).
MUTUAL_HILL_RADII_SEPARATION = 10

# Real protoplanetary disks are observed (sub-mm/ALMA continuum surveys)
# to carry more solid mass around more massive stars, and steeply so --
# not the flat "1 + log10(solar_masses)" scaling the old
# StarSystem.estimate_num_objects formula used. Disk dust-mass-vs-
# stellar-mass surveys (Andrews et al. 2013, Taurus; Pascucci et al. 2016,
# multi-region) find M_dust ~ M_star^1.8 in young (~1-3 Myr) star-forming
# regions, steepening further (~M_star^2.7) in older ones. This scales
# `physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2` for any given
# star relative to the Sun -- see `utils.disk_surface_density_scale`,
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
# (`physical_constants.SNOW_LINE_AU_AT_1_LSUN`/`utils.snow_line_au`) --
# our own Solar System's own giant-planet/Kuiper-belt region extends to
# roughly this same multiple (~18x) of its own 2.7 AU snow line. See
# `StarSystem._estimate_max_objects_from_disk_physics`.
DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER = 18

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
# planetPhysics.reconcile_zone_and_class correctly reclassifies it away
# rather than reporting a physically inconsistent class, which can undo the
# very guarantee that placement was satisfying. Retrying with a fresh
# placement is simpler and more robust than trying to reshuffle every
# neighbor's spacing to protect one body's zone in place.
MAX_SYSTEM_GENERATION_ATTEMPTS = 8

# --- Binary System Generation Parameters ---

# Holman & Wiegert (1999), AJ 117:621 -- their P-type (circumbinary) fit,
# for a planet orbiting both stars of a close pair
# (utils.holman_wiegert_circumbinary_a_crit_au), was tested over mu (the
# lighter star's mass fraction) in [0.1, 0.5] and e in [0.0, 0.7]; inputs
# are clamped to this range rather than extrapolated.
HOLMAN_WIEGERT_P_TYPE_MU_RANGE = (0.1, 0.5)
HOLMAN_WIEGERT_P_TYPE_ECCENTRICITY_RANGE = (0.0, 0.7)

# Probability an S-type (wide) binary is chosen over a P-type (close) one
# when BINARY_SYSTEM is True and WIDE_BINARY is left at None. A plain
# coin-flip, not an attempt to model true field-star multiplicity
# demographics (real wide pairs vastly outnumber sub-0.25-AU pairs) --
# BINARY_SYSTEM's own binary-vs-none split is already a game-design toggle
# rather than a demographic model, so this keeps the two forceable binary
# configurations symmetric with each other.
WIDE_BINARY_DEFAULT_CHANCE = 0.5

# S-type binary separation is sampled log-uniformly between these bounds
# (see doubleStar.WideBinaryPair.generate) rather than with plain
# random.uniform -- real wide-binary separations are observed roughly
# log-uniform/log-normal over several decades (Duquennoy & Mayor 1991;
# Raghavan et al. 2010), and a flat-uniform draw would spend almost all its
# density in the single largest order of magnitude (the same reasoning
# `planetPhysics.generate_moons` already applies to moon spacing). The lower
# bound sits comfortably above the close/P-type BinaryStarProxy separation
# range (0.05-0.25 AU) so the two configurations' separations never overlap;
# the upper bound keeps generated pairs within the separation regime that
# typically survives Galactic tidal shear and passing-star perturbations
# over a stellar lifetime (real wide pairs increasingly disrupt beyond
# roughly 0.1-0.2 pc, ~20,000-41,000 AU -- Jiang & Tremaine 2010), without
# this generator needing to model that disruption process directly.
WIDE_BINARY_SEPARATION_MIN_AU = 50.0
WIDE_BINARY_SEPARATION_MAX_AU = 10000.0

# Cap on the sampled binary orbital eccentricity for an S-type pair (see
# doubleStar.WideBinaryPair.generate's thermal-distribution sampling).
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
# distinct from html/lib/galaxymap.py's own, unrelated galaxy-scale
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
# through the standard mass-luminosity-based main-sequence lifetime scaling
# (the same t = 10 Gy * (M/Msun)^-2.5 relation used to derive the
# STAR_EVOLUTION table above, anchored to the Sun), then extended by the
# standard rule-of-thumb that a star spends about 90% of its total lifetime
# on the main sequence -- so the star must already be at least as old as
# that main-sequence lifespan (it has to have finished that phase to be
# observed as an evolved class), and at most its total (MS + post-MS)
# lifespan.
SOLAR_MS_LIFESPAN_GY = 10.0
MS_LIFESPAN_MASS_EXPONENT = -2.5
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
# A sub-solar-mass progenitor (roughly under ~0.88 Msun, per the
# SOLAR_MS_LIFESPAN_GY/MS_LIFESPAN_MASS_EXPONENT formula above) already has
# a main-sequence lifespan longer than the universe itself -- it couldn't
# possibly have evolved off the main sequence yet, so it must never be
# accepted as an evolved-class star's progenitor mass in the first place
# (see Star._sample_evolved_star_mass_sol, used by Star.generate_star).
# This caps the reject-and-resample loop that enforces that, the same way
# SECTOR_MAX_PLACEMENT_ATTEMPTS caps SpaceSector's own reject-and-resample
# loop: most classes' allowed mass range sits entirely clear of the cutoff
# (so the very first draw is accepted) and only the class III (Giant) low
# end (0.8-~0.88 Msun, out of an 0.8-8 Msun range) ever needs a resample,
# so 100 attempts leaves effectively zero chance of exhausting this for any
# class it's actually applied to.
EVOLVED_STAR_MASS_MAX_RESAMPLE_ATTEMPTS = 100
PERCENT_SOL_THRESHOLD_LOW = 0.01
PERCENT_SOL_THRESHOLD_HIGH = 2
RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD = 100000
PERCENT_MULTIPLIER = 100
SPECTRAL_PROBABILITIES_LARGE_STAR = {'O': 10, 'B': 20, 'A': 30, 'F': 30, 'G': 10, 'K': 0, 'M': 0}
SPECTRAL_PROBABILITIES_NORMAL = {'O': 0.0001, 'B': 0.12, 'A': 0.6, 'F': 3.0, 'G': 7.6, 'K': 12.1, 'M': 76.45}

# Real stellar-multiplicity surveys consistently find binary/multiple
# companionship rising with primary mass, not a single flat rate: low-mass
# M dwarfs are the *least* likely to have a companion, while the most
# massive O stars are very nearly certain to. Keyed by spectral letter
# (`Star.type[0]`) the same way SPECTRAL_PROBABILITIES_NORMAL above is,
# and used the same way when `SystemConfig.BINARY_SYSTEM` is left at its
# default `None` -- see `StarSystem._should_generate_binary`. Anchor
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
#               above that floor.
#   O: 0.90  -- Moe & Di Stefano 2017 (ApJS 230:15) revise O-star
#               multiplicity up to 94 +/- 14%; kept just under certainty
#               given that real uncertainty rather than treating O stars
#               as *always* binary.
BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS = {
    'O': 0.90, 'B': 0.65, 'A': 0.55, 'F': 0.47, 'G': 0.44, 'K': 0.40, 'M': 0.26,
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
# (SOLAR_MS_LIFESPAN_GY * M^MS_LIFESPAN_MASS_EXPONENT above): subgiant up
# to 1.1, giant or supergiant up to 1.2, then a remnant. Real post-MS
# phases take roughly 10-20% of the main-sequence time.
SUBGIANT_PHASE_END_MS_FRACTION = 1.1
GIANT_PHASE_END_MS_FRACTION = 1.2

# Stars at or above this initial mass end as supergiants and then neutron
# stars or black holes (core collapse); below it, as giants and white dwarfs.
SUPERGIANT_MIN_MASS_SOL = 8.0
# Giants at or above this initial mass are bright giants (II), not III.
BRIGHT_GIANT_MIN_MASS_SOL = 4.0
# Giant-branch luminosity (log-uniform, Lsun) and temperature (K) ranges.
GIANT_LUMINOSITY_RANGE_SOL = (50.0, 1000.0)
BRIGHT_GIANT_LUMINOSITY_RANGE_SOL = (1000.0, 10000.0)
GIANT_TEMPERATURE_RANGE_K = (3500.0, 5000.0)
# A subgiant brightens to about twice its main-sequence luminosity and
# cools toward the base of the giant branch while crossing.
SUBGIANT_MAX_BRIGHTENING = 2.0
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

# Piecewise main-sequence mass-luminosity relation, L = coeff * M^exponent
# (Lsun, Msun), each piece for M below its max_mass_sol and the last one
# open-ended (Duric 2004; Salaris & Cassisi 2005).
MS_MASS_LUMINOSITY_PIECES = (
    {"max_mass_sol": 0.43, "coeff": 0.23, "exponent": 2.3},
    {"max_mass_sol": 2.0, "coeff": 1.0, "exponent": 4.0},
    {"max_mass_sol": 55.0, "coeff": 1.4, "exponent": 3.5},
    {"max_mass_sol": None, "coeff": 32000.0, "exponent": 1.0},
)
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
WD_LUMINOSITY_RANGE_SOL = (1.0e-5, 100.0)
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

# Stellar populations by age (galaxyDensity.population_densities,
# stellarPopulation). Each has its own age range (uniform within it); the
# three disk populations share the disk's star formation evenly in time,
# so each one's share of disk stars is its age range's length over
# STAR_FORMATION_AGE_RANGE_GY's. Young stars sit close to the plane and
# crowd the spiral arms; old stars are puffed up by billions of years of
# scattering and spread evenly in azimuth. Scale heights are relative to
# the galaxy's own disk scale height (the old disk's, ~350 pc in the Milky
# Way: young ~80 pc, intermediate ~200 pc). The bulge is old (8-12 Gy).
STELLAR_POPULATION_AGE_RANGES_GY = {
    "young": (0.0, 0.1), "intermediate": (0.1, 3.0), "old": (3.0, 10.0), "bulge": (8.0, 12.0),
}
STELLAR_POPULATION_SCALE_HEIGHT_RATIO = {"young": 80 / 350, "intermediate": 200 / 350, "old": 1.0}
STELLAR_POPULATION_ARM_AMPLITUDE = {"young": 0.9, "intermediate": 0.4, "old": 0.0}

# Pre-placed bright stars (stellarPopulation.sample_bright_stars): the
# mass grid, in log-spaced cells over the IMF's range, on which the chance
# of a star being at least the threshold's brightness is tabulated.
BRIGHT_STAR_MASS_GRID_CELLS = 4000

# Binary mass ratio q = M2/M1, uniform (Moe & Di Stefano 2017, ApJS 230:15,
# find it close to flat); the secondary shares the primary's age.
BINARY_MASS_RATIO_RANGE = (0.1, 1.0)

# --- Planet Classification Data ---

# A dictionary defining the properties of different planet classes.
# Each class has a description, composition, radius range (in meters),
# habitable zone compatibility, atmosphere type, and planet type ('t' for terrestrial, 'g' for gas giant).
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
        # "with a thin atmosphere" moved off the description and onto the
        # atmosphere field itself ("a thin mix of...") -- the description
        # used to collide with the render template's own "with an
        # atmosphere of {atmosphere}" clause (planetData.py's
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
        "h": True, "e": False, "c": False,
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
        # Same real small-body size-frequency reasoning as Class C -- Ceres
        # (469.7km) sits near this range's own ceiling, but the population
        # as a whole skews toward its smaller end.
        "size_mode": 0.25,
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
        "h": False, "e": True, "c": False,
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
        "description": "a volcanic world with shallow seas and bacterial life",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False,
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
        "description": "a rocky, barren world with simple life",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False,
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
        "description": "a desert world with minimal water (less than 10% liquid water)",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False,
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
        "h": False, "e": True, "c": True,
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
        "h": True, "e": True, "c": True,
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
        "h": False, "e": True, "c": False,
        # How far through the ecosphere zone's own [inner, outer] AU range
        # this class is generated, via a `utils.sample_bounded_bell` draw
        # centered here instead of the zone's full width being equally
        # likely (see `planetPhysics.generate_planet_properties`'s
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
        "description": "a marginally habitable world with vegetation",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 7500),
        # No single real analog; skewed toward the lower-mid range --
        # smaller worlds retain a thinner, more "marginal" atmosphere more
        # easily, consistent with this class's own tuning (see below).
        "size_mode": 0.35,
        "h": False, "e": True, "c": False,
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
        "description": "a terrestrial Earth-like world",
        "composition": "silicon, iron, and magnesium",
        "radius_range": (5000, 10000),
        # Earth (6,371km) sits ~27% through this range -- exact.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False,
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
        # greenhouse effect; see planetPhysics.calculate_atmospheric_conditions).
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
        "h": False, "e": True, "c": False,
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
        # Earth-scale anchor, see Class E's note.
        "size_mode": 0.27,
        "h": False, "e": True, "c": False,
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
        "h": False, "e": True, "c": False,
        # Cold/glaciated -> pushed toward the zone's colder outer edge, same
        # as K -- see Class K's own `zone_position_mode` note.
        "zone_position_mode": 0.90,
        "atmosphere": "a mix of oxygen, nitrogen, and argon (thinning with age)",
        "type": "t",
        # Icy/glaciated surfaces reflect far more sunlight than the default
        # rocky/Earth-like range (0.12, 0.35) -- real ice/snow Bond albedo is
        # roughly 0.5-0.9 (e.g. Europa ~0.68, Enceladus ~0.81). Gives P a
        # genuine cold bias from the unclamped physics instead of relying on
        # a post-hoc temperature clamp (see planetPhysics.calculate_atmospheric_conditions).
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
        "h": False, "e": True, "c": False,
        "atmosphere": "a variable mix (thin to dense) of nitrogen, oxygen, and argon",
        "type": "t",
        "life_chemical": ["Chlorophyll a", "Blue-Optimized Porphyrins", "Melanin"],
        "age_ranges": {
            "fast": (0.05, 0.1),
            "normal": (4.0, 9.0),
            "slow": (25.0, 50.0)
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
    "V": {
        "description": "a Super-Earth with high gravity",
        "composition": "iron, iridium, tungsten, and nickel",
        "radius_range": (10000, 15000),
        # The real rocky-to-gaseous transition radius (where planets start
        # retaining a significant H/He envelope) is commonly cited around
        # 1.5-1.6 Earth radii (~9,500-10,200km) -- right at this range's own
        # floor, so peaked low to stay on the rocky side of that boundary.
        "size_mode": 0.25,
        "h": False, "e": True, "c": False,
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
# planetPhysics._choose_weighted_planet_class), so dropping W's tiny 0.0001
# share needed no redistribution either.
PLANET_CLASS_PROBABILITIES = {
    'A': 0.1400, 'B': 0.0725, 'C': 0.2365, 'D': 0.0142, 'E': 0.0239, 'F': 0.0335,
    'G': 0.0432, 'H': 0.0915, 'I': 0.0722, 'J': 0.0531, 'K': 0.0142, 'L': 0.0335,
    'M': 0.1345, 'N': 0.0239, 'O': 0.0045, 'P': 0.0046, 'Q': 0.0001,
    'T': 0.0001, 'V': 0.0045
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
        "evolutionary_constraint_notes": "burns through its fuel and goes supernova before planets can sufficiently cool to form stable liquid oceans, meaning only the absolute fastest, most primitive precursor photochemistry could theoretically occur"
    },
    "B": {
        "lifespan_gy": (0.01, 1.5), # 10 Million - 1.5 Billion years
        "supported_evolutionary_scales": ["fast"],
        "potentially_viable_chemicals": ["Retinal"],
        "evolutionary_constraint_notes": "is highly volatile and short-lived, providing barely enough time for planetary cooling and the emergence of single-celled life utilizing simple proton-motive gradients"
    },
    "A": {
        "lifespan_gy": (1.5, 4.3), # 1.5 Billion - 4.3 Billion years
        "supported_evolutionary_scales": ["fast", "normal"],
        "potentially_viable_chemicals": ["Retinal", "Melanin"],
        "evolutionary_constraint_notes": "allows enough time for the development of early biospheres and radiotrophic organisms adapted to high-energy radiation, but likely dies before complex, slow-evolving porphyrin-based photosynthesis can be perfected"
    },
    "F": {
        "lifespan_gy": (4.3, 9.1), # 4.3 Billion - 9.1 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Blue-Optimized Porphyrins (+ GFPs)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "provides a long enough main-sequence window for complex biospheres to emerge, including advanced UV-shielding and biofluorescent adaptations"
    },
    "G": {
        "lifespan_gy": (9.1, 17.5), # 9.1 Billion - 17.5 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Chlorophyll a (Standard Porphyrins)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "is the solar standard, providing stable, long-term conditions ideal for the slow evolution of highly complex, oxygenic photosynthesis based on the tetrapyrrole chlorin ring"
    },
    "K": {
        "lifespan_gy": (17.5, 73.5), # 17.5 Billion - 73.5 Billion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Chlorophyll a (Standard Porphyrins)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "is exceptionally stable and long-lived, offering tens of billions of years for slow-evolving biospheres to reach climax ecologies and adapt to slight red-shifts in the stellar spectrum"
    },
    "M": {
        "lifespan_gy": (73.5, 5500.0), # 73.5 Billion - 5.5 Trillion years
        "supported_evolutionary_scales": ["fast", "normal", "slow"],
        "potentially_viable_chemicals": ["Retinal", "Melanin", "Bacteriochlorophylls (BChls)", "Zinc-Bacteriochlorophyll"],
        "evolutionary_constraint_notes": "is the longest-lived star in the universe; while early flaring requires rapid adaptation (e.g., radiotrophic melanin) during the initial hundreds of millions of years, its trillions of years of stability allow deep-infrared anoxygenic photosynthesis to dominate permanently"
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
HABITABLE_PLANET_CLASSES = ['E', 'F', 'G', 'H', 'K', 'L', 'M', 'O', 'P', 'V']
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
    # ("an adaptable world", "a terrestrial Earth-like world", etc.), so
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

MOON_BLACKLIST = ['Q', 'V']
"""
list: A list of planet class codes that cannot be generated as moons.
"""

CO2_MAX_GREENHOUSE_FACTOR = 500
"""
int: A generous safety ceiling on greenhouse_factor (planetPhysics.py's
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

# --- Exotic Stellar Phenomena (phenomenonGen.py, compactRemnant.py,
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
# wants them "a smattering (rare)": 2% of black holes, tweak here.
BLACK_HOLE_MASS_RANGE_SOLAR = (5.0, 20.0)
BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR = (1e2, 1e5)
BLACK_HOLE_INTERMEDIATE_MASS_CHANCE = 0.02

# A quiescent supermassive black hole sits at the center of every galaxy
# whose nucleus isn't an active quasar (generate.add_galactic_nucleus).
# Sagittarius A* is ~4.3e6 Msun (GRAVITY Collaboration 2019, A&A
# 625:L10); quiescent nuclei in Milky-Way-like disks span roughly
# 1e6-1e8 Msun, below the 1e8-1e10 of QUASAR_BLACK_HOLE_MASS_RANGE_SOLAR.
# Drawn log-uniform.
BLACK_HOLE_SUPERMASSIVE_MASS_RANGE_SOLAR = (1e6, 1e8)

# A quiescent nucleus accretes far below its Eddington limit: Sgr A*
# shines at ~1e-9 of it (Genzel, Eisenhauer & Gillessen 2010, RvMP
# 82:3121), low-luminosity nuclei up to ~1e-5. Drawn log-uniform; it
# always keeps a (faint) accretion flow.
BLACK_HOLE_SUPERMASSIVE_EDDINGTON_RATIO_RANGE = (1e-9, 1e-6)

# Bulge velocity dispersion for a supermassive black hole's sphere of
# influence, r = G*M / sigma^2 (~2 pc for Sgr A* at ~100 km/s), which
# stands in for the Hill sphere a stellar-mass one uses.
BLACK_HOLE_SUPERMASSIVE_VELOCITY_DISPERSION_KMS = 100.0

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
BLACK_HOLE_ACCRETION_DISK_CHANCE = 0.15

# Neutron star mass range in solar masses -- real measured masses cluster
# tightly around ~1.4 Msun with a hard floor near the Chandrasekhar-like
# collapse threshold and a soft ceiling near the maximum mass general
# relativity allows before further collapse to a black hole (Ozel &
# Freire 2016, ARA&A 54:401, compiling radio pulsar timing masses).
NEUTRON_STAR_MASS_RANGE_SOLAR = (1.1, 2.2)

# Neutron star radius range in km, from NICER x-ray pulse-profile
# modeling of PSR J0030+0451 and PSR J0740+6620 (Miller et al. 2019,
# ApJL 887:L24; Miller et al. 2021, ApJL 918:L28), which converged on
# radii clustering near 11-13 km across the measured mass range.
NEUTRON_STAR_RADIUS_RANGE_KM = (10.0, 13.0)

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

# Chance a generated neutron star is actively pulsing (vs. a quiescent/
# non-pulsing remnant, e.g. an old isolated neutron star whose beam no
# longer sweeps past this vantage point), and the chance a pulsing one
# is specifically a recycled millisecond pulsar rather than a young one.
NEUTRON_STAR_PULSAR_CHANCE = 0.7
PULSAR_MILLISECOND_CHANCE = 0.3

# Neutron star surface temperature range in Kelvin, from thermal
# X-ray emission measurements of young-to-middle-aged isolated
# neutron stars (Potekhin et al. 2020, A&A Rev 28:3, cooling-curve
# review) -- hotter shortly after formation, cooling over ~10^5-10^6
# years toward the low end of this range.
NEUTRON_STAR_SURFACE_TEMPERATURE_RANGE_K = (5e4, 3e6)

# Age since the core-collapse supernova that formed the compact
# remnant, in billions of years -- capped well below
# UNIVERSE_AGE_GY since these are always described as galactic
# (Population I/II), not primordial, remnants.
COMPACT_REMNANT_AGE_RANGE_GY = (0.001, 10.0)

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
# PLANET_CLASSES (TODO item 28; reasoning, sources and the family tables in
# docs/design/nebula-and-asteroid-field-classes.md). I and O are unused so
# they aren't read as 1 and 0; X-Z are reserved. A-Q are nebulae, R-W are
# supernova remnants ("stellar remnants" means supernova remnants only,
# Boss 2026-09-30; compact objects keep their own tables).
#
# Each class: `name`, `family` (a NEBULA_FAMILIES key, or
# "supernova-remnant"), `species` (dominant contents), `density_range_cm3`
# (particle density nH, log-uniform), `temperature_range_k` (log-uniform),
# `extinction_range_av` (optical extinction in magnitudes, log-uniform
# unless the low end is 0), `center` (the central-object rule item 27
# builds on), and `frequency` (relative weight within its family group).
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

# Sedov-Taylor phase expansion coefficient and exponent: R(t) = C *
# t^(2/5), the self-similar blast-wave solution for a remnant that has
# swept up enough interstellar mass to decelerate from its initial
# free-expansion phase (Taylor 1950, Proc. Roy. Soc. A 201:159; Sedov
# 1959, "Similarity and Dimensional Methods in Mechanics"). C is
# calibrated so a ~1,000-year-old remnant (roughly Tycho's/Cassiopeia
# A's real age) comes out to a few light-years across, matching their
# real observed sizes.
SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY = 0.35
SEDOV_TAYLOR_TIME_EXPONENT = 2 / 5

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
    ("O", 0, 9, ("C", "D", "E"), 1.0),
    ("B", 0, 2, ("C", "D", "G"), 0.5),
    ("B", 3, 9, ("F", "G"), 0.05),
    ("A", 0, 9, ("F",), 0.01),
)
"""tuple: Nebulae a sector grows around its own stars (TODO item 27), as
`(spectral letter, lowest subclass, highest subclass, classes, chance)`
for a main-sequence primary. Only O and early-B stars put out enough
ultraviolet below 91.2 nm to ionize hydrogen, so every O star sits in an
H II region (C-E) and half the B0-B2 stars do (Osterbrock & Ferland
2006); later B and A stars light dust without ionizing it, a reflection
nebula (F, or G for the hotter ones) when a dusty cloud happens to be near
(van den Bergh 1966, AJ 71:990 catalogs ~150 within a few kpc). The class
is drawn among `classes` by NEBULA_CLASSES frequency."""

PLANETARY_NEBULA_CENTRAL_STAR_TYPES = ("O3VII", "O5VII", "O7VII", "O9VII", "B0VII")
"""tuple: The central star a planetary nebula is generated around (TODO
item 27): the exposed hot core of a dying 0.8-8 Msun star, 30,000 K and
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
# per-star rate, then a mass log-uniformly inside it. Low-mass rogues
# dominate -- disk scattering ejects small bodies while giants stay bound
# -- and no single power law fits both ends, hence bins.
#   - terrestrial 0.1-2 Earth masses, 5 per star (2-10; Johnson et al.
#     2020, AJ 160:123; Mroz et al. 2020, ApJL 903:L11);
#   - sub-Neptune / ice giant 2-20, 1 per star (a default: Sumi et al.
#     2023 find Neptune-mass candidates but give no rate);
#   - Saturn-class 20 Earth masses to 1 Jupiter mass, 0.25 per star (a
#     default filling the gap);
#   - Jupiter-mass 1-13 Jupiter masses, 0.25 per star (Mroz et al. 2017,
#     Nature 548:183, superseding Sumi et al. 2011's 1.8).
# Above 13 Jupiter masses is a brown dwarf (ROGUE_BROWN_DWARF_*).
ROGUE_PLANET_MASS_BINS = {
    "terrestrial": (0.1, 2.0, 5.0),
    "sub-neptune": (2.0, 20.0, 1.0),
    "saturn": (20.0, 317.8, 0.25),
    "jupiter": (317.8, 13 * 317.8, 0.25),
}
"""dict: bin name -> `(min_mass_earth, max_mass_earth, per_star_rate)`.
The rates' sum is the rogue-planet rate per star (6.5)."""

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

# Chance a generated rogue planet retains detectable internal heat
# (radiogenic/primordial, same mechanism warming Earth's mantle or
# Jupiter's interior) worth describing, vs. one long since frozen solid
# with no internal activity.
ROGUE_PLANET_INTERNAL_HEAT_CHANCE = 0.4
ROGUE_PLANET_MOON_CHANCE = 0.2

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

# How many comets a star that DOES get any (per SYSTEM_COMET_CHANCE, or
# `SystemConfig.COMETS = True` forcing at least this many) actually
# generates -- a small handful, not a full population (this generator
# only models a system's few most notable comets, the same way
# ASTEROID_BELT models zero-or-one belts, not an exhaustive minor-body
# catalog).
SYSTEM_COMET_COUNT_RANGE = (1, 3)

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

# Asteroid field classes (TODO item 31, Boss 2026-09-30: "a digit and part
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

QUASAR_BLACK_HOLE_MASS_RANGE_SOLAR = (1e8, 1e10)
"""
tuple: Central black hole mass range, in solar masses, sampled
log-uniformly -- the span of virial mass estimates across the SDSS quasar
catalog (Shen et al. 2011, ApJS 194:45). Below ~1e8 a black hole can't
reach quasar luminosity even at its Eddington limit.
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

PHENOMENON_TYPE_CHOICES = (
    "black-hole", "neutron-star", "nebula", "supernova-remnant",
    "rogue-planet", "comet", "asteroid-field", "quasar",
)
"""
tuple: The valid `--type` values `phenomenonGen.py` accepts. Omitting
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

# --- Sector-Level Exotic Phenomena (sectorGen.py/galaxyGen.py) ---
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
    # Terrestrial alone: 0.7 pc^-3 (0.5-1.4), Jupiter-mass 0.035.
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
    "neutron-star": 7e-4,
    "black-hole": 1.4e-4,
    # Giant molecular clouds: 1e-6 to 1e-5 (Kennicutt & Evans 2012, ARA&A
    # 50:531): cloud centers, generated as dark-family nebulae (classes
    # M-Q). GMC_ARM_FILLING_FACTOR is the matching fraction of arm volume
    # inside a cloud, for reference.
    "molecular-cloud": 5e-6,
    # ~20,000 planetary nebulae galaxy-wide (Frew & Parker 2010, PASA
    # 27:129).
    "planetary-nebula": 3e-8,
    # 1e-8 to 1e-7, really n_* * rho_gas (Draine 2011); 1e-8 counts
    # faint remnants too (~1,000 are detectable, ~1.5e-9).
    "supernova-remnant": 1e-8,
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

PHENOMENON_RATE_SCALE = {kind: 1.0 for kind in PHENOMENON_DENSITY_PC3}
"""dict: A multiplier per kind (default 1.0) to dial a rate up or down
without touching its research value -- e.g. 0.1 for "rogue-planet" gives
about 6 rogues per local sector instead of about 58 (Boss chose the full
rate on 2026-09-30)."""

GMC_ARM_FILLING_FACTOR = 0.015
"""float: Share of spiral-arm volume inside a giant molecular cloud
(0.01-0.02), for TODO item 27's cloud placement."""

GMC_GAS_DENSITY_EXPONENT = 1.4
"""float: Molecular clouds scale with gas density to this power
(Schmidt-Kennicutt), for TODO item 27."""

HVS_REFERENCE_RADIUS_PC = 8000.0
"""float: Galactic radius the "hypervelocity-star" density is quoted at."""

INTERSTELLAR_DEBRIS_DENSITY_PC3 = 1e12
"""float: Interstellar comets and planetesimals (m-km) per pc^3 at
REFERENCE_STELLAR_DENSITY_PC3, 1e11-1e14 (Engelhardt et al. 2017, AJ
153:133; Seligman & Laughlin 2020, ApJL 896:L8) -- far too many for rows,
so a sector shows it as a figure (`queryDb.interstellar_debris_count`)."""

RUNAWAY_STAR_SPEED_RANGE_KMS = (30.0, 200.0)
"""tuple: A runaway star's speed relative to its neighbors, km/s
(log-uniform; Tauris 2015)."""

HYPERVELOCITY_STAR_SPEED_RANGE_KMS = (500.0, 1000.0)
"""tuple: A hypervelocity star's speed, km/s (Brown 2015)."""


def phenomenon_rate_per_star(kind):
    """Expected count of `kind` (a PHENOMENON_DENSITY_PC3 key) per star:
    its density over REFERENCE_STELLAR_DENSITY_PC3, times its
    PHENOMENON_RATE_SCALE."""
    return PHENOMENON_DENSITY_PC3[kind] * PHENOMENON_RATE_SCALE.get(kind, 1.0) / REFERENCE_STELLAR_DENSITY_PC3


# --- Galaxy Random-Start Generation (galaxyGen.py) ---

GALAXY_RADIUS_PC = 15000.0
"""
float: The Milky Way's real approximate radius, in parsecs (commonly cited
~15 kpc). Only a fallback for the Galaxy Map's camera range before any
`generate.py plan` has run; once a galaxy is planned, its stored outline
(`galaxy_layer`) is the only bound anything uses -- generation never picks
or accepts an address outside it.
"""

RANDOM_START_NEIGHBORHOOD_RADIUS_LY = 100.0
"""
float: The neighborhood radius, in light-years, `galaxyGen.py`'s
no-argument "random start" mode generates around its randomly chosen
starting sector -- every not-yet-generated sector within this radius, in
every direction, per this feature's own request.
"""

RANDOM_START_MAX_PLACEMENT_ATTEMPTS = 1000
"""
int: Retry cap for picking a random, not-yet-occupied sector address before
`galaxyGen.py`'s random-start mode gives up -- generous, since even a
fairly well-populated galaxy leaves overwhelmingly more addresses empty
than occupied (see docs/design/galaxy-coordinate-system.md section 9's
storage-analysis addendum: ~320 billion addressable sector slots total).
"""

# --- Navigation Parameters ---

WARP_FACTORS_FOR_NAV = (1, 2, 4, 8, 9, 9.5, 9.9, 9.995)
"""
The warp factors NAV output reports travel time at (see
`stellarObjects.navigation.warp_travel_times`) -- the rows of Boss's warp
table (docs/design/navigation-frames.md, "Travel speeds"), from warp 1
(exactly c) up to 9.995, where the curve climbs steeply toward its
asymptote at warp 10.
"""

WARP_VELOCITY_EXPONENT = 10 / 3
"""
The exponent of the warp curve's base term, `warp_factor **
WARP_VELOCITY_EXPONENT` (in multiples of light-speed) -- "warp factor to
the 3.33...". Below about warp 9 the whole curve is effectively this term;
see `stellarObjects.navigation.warp_speed_c` for the full formula.
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
`stellarObjects.navigation.fold_travel_times`) -- the rows of Boss's fold
table.
"""

FOLD_SPEED_COEFFICIENT = 6
"""
The coefficient in the dimensional fold curve, speed in c = `6 F^4 / (10 -
F)` (see `stellarObjects.navigation.fold_speed_c`).
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
NAV adjacency graph (`stellarObjects.navGraph.build_knn_adjacency`) that
optimal-route pathfinding runs over. Symmetrized after building (see that
function's docstring), so a system can end up connected to more than `k`
neighbors if others chose it as one of theirs.
"""

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

FACILITY_DEFAULT_ORBIT_RADII = 3.0
"""float: An orbital facility around a planet or moon with no distance
given orbits at this many of its host's radii."""

# --- Galaxy pre-placement (schema v43) ---

BRIGHT_STAR_MIN_LUMINOSITY_SOL = 100.0
"""float: Every star at least this bright (solar luminosities) is generated
and placed galaxy-wide right after `generate.py plan`, before any sector is
filled (`bright_stars`, schema v43). Its sector is still generated later,
around it. Boss, 2026-09-30 (500), lowered to 100 on 2026-10-01;
`--bright-star-min-luminosity` raises it for a quick test galaxy. It can't
go below the brightest white dwarf (`WD_LUMINOSITY_RANGE_SOL`). The value a
scatter used is stored in `galaxy_shape.bright_star_min_luminosity_sol`, and
filling reads that, not this, so retuning it can't make a fill
double-count or skip stars. See
/mnt/project-files/galaxy-studies/bright-star-preplacement-plan.md."""
