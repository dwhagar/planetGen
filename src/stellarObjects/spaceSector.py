# stellarObjects/spaceSector.py

"""
Space Sector Storage
=====================

This module contains the `SpaceSector` and `SectorSystemEntry` classes, which
store a named collection of star systems laid out in three-dimensional space.

Where `StarSystem` represents the contents of a single system, `SpaceSector`
represents a whole cubic region of space containing many of them, each placed
at its own `(x, y, z)` position in light-years relative to the sector's
center (its origin), so that distances and neighbors between systems can be
reasoned about.

Sector shape and density
-------------------------
A sector is a cube `program_constants.DEFAULT_SECTOR_EDGE_LY` light-years on
a side (~11.5 ly, ~1,521 cubic ly), spanning `[-edge/2, +edge/2]` on each
axis around the origin. `physical_constants.LOCAL_STELLAR_DENSITY_LY3` --
how many systems, on average, occupy one cubic light-year -- is not a
made-up number; it is derived from real observational surveys of the Sun's
own neighborhood:

    - Mamajek's stellar-density review
      (https://www.pas.rochester.edu/~emamajek/memo_star_dens.html) gives
      ~0.0984 +/- 0.0068 stars/pc^3 in the solar neighborhood (star-forming
      objects only; excludes brown dwarfs, neutron stars, black holes).
      Converting with 1 pc = 3.26156 ly (so 1 pc^3 = 34.706 ly^3) gives
      ~0.00284 stars/ly^3.
    - RECONS' 10-parsec census (https://chview.nova.org/solcom/stars/pc10.htm)
      independently counts on the order of 270-414 stellar objects within 10
      pc of the Sun -- a sphere of ~4,189 pc^3 (145,120 ly^3) -- giving a
      consistent ~0.0019-0.0029 objects/ly^3.

Both land in the same neighborhood, so this module uses ~0.00284 systems per
cubic light-year: on average, one system per ~352 cubic light-years. This is
deliberately sparse -- it is why the Sun's own nearest neighbor, Alpha
Centauri, is 4.37 ly away -- so a single sector this size will typically hold
only a handful of systems (`expected_system_count`), which is realistic
rather than a bug.

Minimum separation (Hill spheres)
-----------------------------------
Random placement (`add_system` with no explicit `position`) keeps two systems
from landing so close together that they'd gravitationally interfere with
each other, using physics this package already computes for every system:
each `Star`/`BinaryStarProxy` has a `system_perimeter` (its Hill sphere
*relative to the Milky Way* -- see `Star.calculate_system_perimeter` --
i.e. the radius beyond which the galaxy's own gravity dominates over the
star's, roughly the same concept that sets the outer edge of the Oort
cloud). Two systems are required to be at least the *sum* of their two Hill
radii apart (`required_separation_ly`), so those radii can touch but not
overlap.

This isn't a new number invented for sector placement -- it's the same
`system_perimeter` already shown in every generated system's output, and it
checks out against real astronomy: with this package's own constants, the
Sun's `system_perimeter` comes out to ~1.71 ly (~108,000 AU), squarely inside
the ~100,000-150,000 AU commonly cited for the Oort cloud's outer boundary.
Scaled by (mass/solar mass)^(1/3) across spectral types (using the typical
main-sequence mass ranges in `physical_constants.SPECTRAL_MASS_RANGES`), it
runs from ~0.7 ly for a small M dwarf up to ~7-8 ly for a massive O star.
Two Sun-like stars' Hill spheres summed come to ~3.4 ly -- which lines up
neatly with the ~3.9 ly mean nearest-neighbor distance the ~0.00284/ly^3
density above implies for a random (Poisson) distribution of systems. In
other words, real stars are spaced, on average, almost exactly as close as
their mutual gravity allows without one's Oort cloud bleeding into the
other's -- and this module's default placement reproduces that.

Position sampling uses `secrets.SystemRandom`, a `random.Random` subclass
backed by the OS's own CSPRNG (`os.urandom`) rather than the deterministic,
seedable Mersenne Twister behind the plain `random` module. It cannot be
seeded to reproduce a past run -- there is no "true random" seed to save --
and it is a separate generator entirely, so sector placement is never
accidentally made predictable by code elsewhere reseeding the global
`random` module (as the rest of this package's star/planet generation does;
see `systemGen.main`).

Distances
---------
`distance_between` (and `SectorSystemEntry.distance_to`) compute straight-line
3D Euclidean distance. At the scale of a single sector (light-years), this is
exact -- no relativistic or cosmological-expansion correction applies; those
only become relevant at scales many, many orders of magnitude larger.

Growth (sample, then fine-tune)
---------------------------------
`SpaceSector.grow_from_seed` fills a sector outward from an already-placed
system, borrowing the "active list" structure from Bridson's Fast Poisson
Disk Sampling algorithm ("Fast Poisson Disk Sampling in Arbitrary
Dimensions", SIGGRAPH 2007 sketch --
https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf) --
maintain a list of points that might still have room for a neighbor,
repeatedly pick a random one, retire it after
`program_constants.SECTOR_GROWTH_POISSON_DISK_K` failed attempts -- but
placing each candidate in two distinct phases rather than pure
reject-and-resample:

1. **Initial sample**: a candidate `StarSystem` is generated (via a
   caller-supplied factory -- a position can't be chosen until the
   candidate's own stellar mass, and hence its own Hill radius, is known),
   then placed at a random spherical-coordinate offset from the active
   parent: distance uniform-by-volume between `hill_radius_ly(parent)` (the
   parent's own Hill radius alone -- a quick, cheap floor) and
   `mean_nearest_neighbor_ly()` (the typical spacing implied by this
   sector's target density -- see below), in a uniformly random direction.
2. **Fine-tune**: the *actual* gravitational check is against whichever
   system turns out to be closest to the candidate's current position --
   not necessarily the parent it was sampled from. If that nearest
   neighbor's and the candidate's combined Hill spheres
   (`required_separation_ly`) overlap, the candidate is pushed directly
   away from that neighbor, out to exactly the required distance, and the
   check repeats (moving away from one neighbor can bring the candidate
   closer to a different one) up to
   `program_constants.SECTOR_GROWTH_FINE_TUNE_MAX_ITERATIONS` times.

A candidate is accepted once fine-tuning converges on a position that's
both clear of every existing system's Hill sphere and still inside the
sector's cube; otherwise this attempt is discarded and, if the active point
has attempts left, a fresh independent candidate is generated and sampled
again from scratch.

`mean_nearest_neighbor_ly()` is the mean nearest-neighbor distance for a 3D
Poisson process at `physical_constants.LOCAL_STELLAR_DENSITY_LY3` --
`Γ(4/3) * (3 / (4*pi*density))^(1/3)`, ~3.9 ly with this module's current
constants -- i.e. a single density-derived value, independent of which
specific stars are involved (unlike the Hill-sphere-based floor, which
depends on the parent's own mass).

By default, growth stops once the sector's system count reaches a
Poisson-distributed target drawn around `expected_system_count()` -- matching
how a real spatial Poisson process's count in a fixed volume is itself
Poisson-distributed, rather than always packing sectors to the same density.
It also stops early (without error) if the active list empties first, i.e.
there simply wasn't room for more at the required spacing.

Named locations (quadrants)
------------------------------
`classify_octant`/`format_named_location` (and
`SectorSystemEntry.named_location`) express a raw `(x, y, z)` position as a
Roman-numeral "quadrant" -- properly an *octant* in 3D, 8 regions rather than
4, kept as "quadrant" here to match this project's own terminology -- plus
the position's three positive magnitudes. Unlike 2D quadrants (I-IV is a
genuine, universal mathematical standard), there is no single authoritative
Roman-numeral numbering for 3D octants -- Wikipedia's "Octant (solid
geometry)" article says so explicitly, recommending plain sign-tuple
notation instead for exactly that reason. `program_constants.SECTOR_OCTANT_LABELS`
adopts a commonly *taught* (not ISO-standardized) extension of the 2D
pattern instead; see that constant's own comment for the exact convention.
This is purely a derived display label, generated on demand from the stored
raw position -- there is no reverse conversion, since raw `(x, y, z)` is the
only form actually persisted.

Persistence
-----------
A sector can be saved to (and reloaded from) JSON. Since Phase 1.5 (see
TODO.md), each entry's `to_dict()` writes TWO things: `config`, the
`SystemConfig` its `StarSystem` was built from (the sector's *recipe* --
star types, forced features, names, orbit counts, and so on), and
`generated`, the full object graph itself (`StarSystem.to_dict()`,
Phase 1's serialization). `from_dict` prefers `generated` when present,
giving an EXACT, byte-for-byte reload of the original system -- not a
fresh roll -- since a fair amount of this package's generation (planet
class rolls, moon coin-flips, name generation, flavor text) is drawn from
the `secrets` module specifically because it should be unpredictable, so
there's no seed that could play it back exactly any other way.

A file written before this addition (or a hand-authored one supplying only
`config`) has no `generated` key, so `from_dict` falls back to the
original behavior instead: a fresh, independently random system is built
from the same recipe -- not the same system come back. The one detail
pinned across THAT fallback path's reload is the star's name: if a
system's config didn't already force one, the name it was actually given
gets baked into the saved config so it doesn't change out from under a
saved sector.
"""

import json
import math
import secrets

from . import physical_constants, program_constants
from .config import SystemConfig
from .systemData import StarSystem

# A CSPRNG-backed generator (os.urandom under the hood), used for every
# position placed in a sector instead of the deterministic, seedable `random`
# module -- see "Position sampling" in the module docstring.
_rng = secrets.SystemRandom()


def distance_between(a, b):
    """
    Computes the exact Euclidean distance, in light-years, between two
    positions (or entries) -- see "Distances" in the module docstring for why
    straight-line 3D distance is exact at this scale.

    Args:
        a (SectorSystemEntry or tuple): The first entry, or a raw `(x, y, z)`
                                        position.
        b (SectorSystemEntry or tuple): The second entry, or a raw
                                        `(x, y, z)` position.

    Returns:
        float: The distance in light-years.
    """
    position_a = a.position if isinstance(a, SectorSystemEntry) else a
    position_b = b.position if isinstance(b, SectorSystemEntry) else b
    return math.dist(position_a, position_b)


def hill_radius_ly(star_system):
    """
    Returns a system's Hill-sphere radius relative to the galaxy, in
    light-years -- see "Minimum separation (Hill spheres)" in the module
    docstring.

    Every generated `Star` (and `BinaryStarProxy`, for binary systems)
    already computes this at generation time as `system_perimeter`, in AU
    (see `Star.calculate_system_perimeter`); this just reads it back
    converted to light-years.

    Args:
        star_system (StarSystem): The system to measure.

    Returns:
        float: The Hill-sphere radius in light-years.
    """
    return star_system.star.system_perimeter * physical_constants.AU_TO_LY


def required_separation_ly(system_a, system_b):
    """
    The minimum distance, in light-years, at which two systems' galactic
    Hill spheres (`hill_radius_ly`) stop overlapping: the sum of each
    system's own Hill radius. Placing two systems closer than this puts each
    one, at least partly, inside the other's zone of gravitational
    dominance -- their Oort-cloud-scale outskirts would be stripped,
    perturbed, or (at the extreme) drawn into the other system rather than
    staying bound to their own star.

    Args:
        system_a (StarSystem): The first system.
        system_b (StarSystem): The second system.

    Returns:
        float: The minimum non-overlapping separation, in light-years.
    """
    return hill_radius_ly(system_a) + hill_radius_ly(system_b)


def mean_nearest_neighbor_ly():
    """
    The mean nearest-neighbor distance, in light-years, for a 3D Poisson
    process at `physical_constants.LOCAL_STELLAR_DENSITY_LY3` -- the
    statistically expected typical spacing between neighboring systems at
    this sector's target density, independent of which specific stars are
    involved (unlike `hill_radius_ly`, which depends on a star's own mass).
    See "Growth (sample, then fine-tune)" in the module docstring for how
    this is used.

    For a 3D Poisson process with intensity (density) `lambda`, the
    nearest-neighbor distance follows a Weibull distribution (shape 3) whose
    mean works out to `Γ(4/3) * (3 / (4*pi*lambda))^(1/3)`
    (`math.gamma(4/3)` for the Gamma function). With this module's current
    density constant, this comes out to ~3.9 ly -- close to the Sun's own
    real nearest neighbor, Alpha Centauri, at 4.37 ly.

    Returns:
        float: The mean nearest-neighbor distance in light-years.
    """
    density = physical_constants.LOCAL_STELLAR_DENSITY_LY3
    return math.gamma(4 / 3) * (3 / (4 * math.pi * density)) ** (1 / 3)


def _sample_poisson_count(mean, rng=_rng):
    """
    Draws a Poisson-distributed non-negative integer with the given mean,
    using only uniform draws from `rng` -- Knuth's simple multiplicative
    algorithm: repeatedly multiply `rng.random()` draws together until the
    running product drops below `exp(-mean)`, counting how many draws that
    took. This is exact (not an approximation), and it stays on this
    module's own `_rng` (true OS entropy, not the seedable global `random`
    module) rather than pulling in `numpy.random.poisson`, which would add a
    new dependency this project doesn't otherwise have.

    This runs in O(mean) draws, so it's only efficient for small means. The
    one place this module uses it (`SpaceSector.grow_from_seed`'s default
    `target_count`, via `expected_system_count`) is realistically single
    digits for the default sector size -- comfortably within where this
    approach is appropriate. Do not reuse this for a mean much above a few
    dozen without switching to a rejection-based algorithm instead.

    Args:
        mean (float): The Poisson distribution's mean (lambda). Values <= 0
                      always return 0.
        rng (random.Random): The generator to draw from.

    Returns:
        int: A Poisson-distributed sample.
    """
    if mean <= 0:
        return 0

    threshold = math.exp(-mean)
    count = 0
    product = 1.0
    while True:
        product *= rng.random()
        if product <= threshold:
            return count
        count += 1


def _random_unit_direction(rng=_rng):
    """
    Samples a uniformly random direction on the unit sphere.

    Uses `theta = uniform(0, 2*pi)` and `phi = acos(uniform(-1, 1))` --
    NOT `phi = uniform(0, pi)`, which would bias samples toward the poles,
    since a fixed-size angular step near the poles covers far less surface
    area than the same step near the equator.

    Args:
        rng (random.Random): The generator to draw from.

    Returns:
        tuple: A unit-length `(dx, dy, dz)` direction vector.
    """
    theta = rng.uniform(0, 2 * math.pi)
    phi = math.acos(rng.uniform(-1.0, 1.0))
    sin_phi = math.sin(phi)
    return (sin_phi * math.cos(theta), sin_phi * math.sin(theta), math.cos(phi))


def _random_point_in_annulus(center, inner_radius, outer_radius, rng=_rng):
    """
    Samples a point uniformly *by volume* (not by radius) within the
    spherical annulus between `inner_radius` and `outer_radius` around
    `center` -- part of `SpaceSector.grow_from_seed`'s initial-sample phase;
    see "Growth (sample, then fine-tune)" in the module docstring.

    Sampling the distance directly as `uniform(inner_radius, outer_radius)`
    would bias points toward the inner radius, since a thin shell near the
    inner radius has less volume than an equally thin shell farther out; the
    cube-root form below corrects for that.

    Args:
        center (tuple): The `(x, y, z)` point the annulus is centered on.
        inner_radius (float): The annulus's inner radius.
        outer_radius (float): The annulus's outer radius.
        rng (random.Random): The generator to draw from.

    Returns:
        tuple: The sampled `(x, y, z)` position.
    """
    u = rng.random()
    distance = (inner_radius ** 3 + u * (outer_radius ** 3 - inner_radius ** 3)) ** (1 / 3)
    dx, dy, dz = _random_unit_direction(rng)
    return (center[0] + distance * dx, center[1] + distance * dy, center[2] + distance * dz)


def _nudge_away(anchor_position, position, target_distance, rng=_rng):
    """
    Moves `position` directly away from `anchor_position`, out to exactly
    `target_distance` -- part of `SpaceSector.grow_from_seed`'s fine-tuning
    phase; see "Growth (sample, then fine-tune)" in the module docstring.

    If `position` coincides exactly with `anchor_position` (the direction
    between them is undefined), a fresh random direction is used instead
    (see `_random_unit_direction`) -- astronomically this never happens by
    chance with continuous random coordinates, but it keeps this function
    well-defined rather than dividing by zero.

    Args:
        anchor_position (tuple): The point being moved away from.
        position (tuple): The current `(x, y, z)` position to nudge.
        target_distance (float): The desired distance from
                                 `anchor_position` after nudging.
        rng (random.Random): The generator to draw a fallback direction
                             from.

    Returns:
        tuple: The nudged `(x, y, z)` position.
    """
    dx = position[0] - anchor_position[0]
    dy = position[1] - anchor_position[1]
    dz = position[2] - anchor_position[2]
    current_distance = math.sqrt(dx * dx + dy * dy + dz * dz)

    if current_distance == 0:
        dx, dy, dz = _random_unit_direction(rng)
    else:
        dx, dy, dz = dx / current_distance, dy / current_distance, dz / current_distance

    return (
        anchor_position[0] + dx * target_distance,
        anchor_position[1] + dy * target_distance,
        anchor_position[2] + dz * target_distance,
    )


def classify_octant(position):
    """
    Classifies a raw `(x, y, z)` position (relative to a sector's center)
    into one of the 8 sign-combination "quadrants" and its three positive
    magnitudes -- see "Named locations (quadrants)" in the module docstring
    and `program_constants.SECTOR_OCTANT_LABELS` for the convention and its
    caveats.

    A coordinate of exactly 0.0 is treated as non-negative, matching the
    convention's own tie-break.

    Args:
        position (tuple): Raw `(x, y, z)` in light-years, relative to
                          center.

    Returns:
        tuple: `(roman_numeral_label, (abs_x, abs_y, abs_z))`.
    """
    x, y, z = position
    signs = (x >= 0, y >= 0, z >= 0)
    for label, sign_x, sign_y, sign_z in program_constants.SECTOR_OCTANT_LABELS:
        if (sign_x, sign_y, sign_z) == signs:
            return label, (abs(x), abs(y), abs(z))

    raise AssertionError("unreachable: SECTOR_OCTANT_LABELS covers all 8 sign combinations")


def format_named_location(position):
    """
    Formats a raw `(x, y, z)` position as a human-readable "named location"
    string, e.g. `"Quadrant III (2.10, 4.40, 1.05 ly from center)"`. See
    "Named locations (quadrants)" in the module docstring.

    Args:
        position (tuple): Raw `(x, y, z)` in light-years, relative to
                          center.

    Returns:
        str: The formatted label.
    """
    label, (abs_x, abs_y, abs_z) = classify_octant(position)
    places = program_constants.SECTOR_LOCATION_DECIMAL_PLACES
    return (
        f"Quadrant {label} ({abs_x:.{places}f}, {abs_y:.{places}f}, "
        f"{abs_z:.{places}f} ly from center)"
    )


class SectorSystemEntry:
    """
    One star system's placement within a `SpaceSector`.

    Attributes:
        star_system (StarSystem): The generated system occupying this slot.
        position (tuple): `(x, y, z)` coordinates in light-years, relative to
                          the sector's center.
        system_config (SystemConfig): The config `star_system` was generated
                                      from, kept so the sector can be saved
                                      and later rebuilt from the same recipe.
    """

    def __init__(self, star_system, position, system_config=None):
        self.star_system = star_system
        self.position = tuple(position)
        self.system_config = system_config if system_config is not None else star_system.system_config

    def distance_to(self, other):
        """
        Computes the Euclidean distance, in light-years, to another entry (or
        a raw `(x, y, z)` position). See `distance_between`.

        Args:
            other (SectorSystemEntry or tuple): The other entry, or a raw
                                                `(x, y, z)` position.

        Returns:
            float: The distance in light-years.
        """
        return distance_between(self, other)

    def named_location(self):
        """
        Returns this entry's position formatted as a "named location"
        string. See `format_named_location`.

        Returns:
            str: The formatted label.
        """
        return format_named_location(self.position)

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this entry.

        The saved config's `name` is pinned to the star's actual generated
        name when the config itself didn't already force one, so a reloaded
        sector keeps the same names even though other details are freshly
        randomized when rebuilt from `config` alone (see the module
        docstring's "Persistence" section) -- the `generated` key below is
        what makes an *exact* reload possible instead, when present.

        Returns:
            dict: With `name`, `position`, `config`, and `generated` (this
                 entry's `star_system.to_dict()`, Phase 1's full
                 object-graph serialization) keys.
        """
        config_dict = self.system_config.to_dict()
        if config_dict.get("name") is None:
            config_dict["name"] = self.star_system.star.name

        return {
            "name": self.star_system.star.name,
            "position": list(self.position),
            "config": config_dict,
            "generated": self.star_system.to_dict(),
        }


class SpaceSector:
    """
    A named collection of star systems laid out within a cubic region of
    three-dimensional space.

    Attributes:
        name (str): The sector's name.
        edge_ly (float): The length, in light-years, of the cubic volume's
                         edge. New systems are randomly placed within
                         `[-edge_ly / 2, edge_ly / 2]` on each axis when no
                         explicit position is given.
        entries (list): The `SectorSystemEntry` instances the sector contains.
    """

    def __init__(self, name, edge_ly=program_constants.DEFAULT_SECTOR_EDGE_LY):
        self.name = name
        self.edge_ly = edge_ly
        self.entries = []

    def __len__(self):
        return len(self.entries)

    @property
    def volume_ly3(self):
        """float: The sector's volume in cubic light-years (`edge_ly ** 3`)."""
        return self.edge_ly ** 3

    def expected_system_count(self):
        """
        Estimates how many systems a sector this size would realistically
        contain, based on real local stellar density (see the module
        docstring for sourcing). This is an average, not a guarantee --
        actual space is clumpy, not uniform -- so treat it as a reasonable
        target to generate around, not an exact count to enforce.

        Returns:
            float: The expected number of systems (`volume_ly3 *
                  physical_constants.LOCAL_STELLAR_DENSITY_LY3`).
        """
        return self.volume_ly3 * physical_constants.LOCAL_STELLAR_DENSITY_LY3

    def _random_position(self, star_system, min_separation_ly):
        """
        Chooses a random `(x, y, z)` point uniformly within the sector's
        cubic volume, far enough from every system already placed.

        Args:
            star_system (StarSystem): The system being placed, needed (when
                `min_separation_ly` is None) to work out its Hill radius.
            min_separation_ly (float or None): The minimum allowed distance
                to every existing entry. If None, the requirement is worked
                out per-neighbor from Hill spheres (`required_separation_ly`)
                instead of a single flat distance -- see "Minimum separation
                (Hill spheres)" in the module docstring.

        Returns:
            tuple: The chosen `(x, y, z)` position.

        Raises:
            ValueError: If no such point could be found within
                       `program_constants.SECTOR_MAX_PLACEMENT_ATTEMPTS` tries.
        """
        half_edge = self.edge_ly / 2
        for _ in range(program_constants.SECTOR_MAX_PLACEMENT_ATTEMPTS):
            candidate = tuple(_rng.uniform(-half_edge, half_edge) for _ in range(3))
            if all(
                distance_between(candidate, entry.position) >= (
                    min_separation_ly if min_separation_ly is not None
                    else required_separation_ly(star_system, entry.star_system)
                )
                for entry in self.entries
            ):
                return candidate

        raise ValueError(
            f"Could not place a new system far enough from every existing system within a "
            f"{self.edge_ly} ly sector after {program_constants.SECTOR_MAX_PLACEMENT_ATTEMPTS} "
            f"attempts; the sector may be too full or too small."
        )

    def add_system(self, star_system, position=None, system_config=None, min_separation_ly=None):
        """
        Adds an already-generated `StarSystem` to the sector.

        Args:
            star_system (StarSystem): The system to place.
            position (tuple, optional): An explicit `(x, y, z)` position in
                light-years. If omitted, a random position is chosen within
                the sector's cube, far enough from every system already
                placed (see `min_separation_ly`).
            system_config (SystemConfig, optional): The config `star_system`
                was generated from, so the sector can later be saved and
                rebuilt from the same recipe. Defaults to
                `star_system.system_config`.
            min_separation_ly (float, optional): A flat minimum allowed
                distance to every existing system when auto-placing,
                overriding the default physics-based check. If omitted (the
                default), each existing system is instead required to be at
                least `required_separation_ly` away -- the sum of the two
                systems' own Hill-sphere radii -- so their zones of
                gravitational dominance never overlap; see "Minimum
                separation (Hill spheres)" in the module docstring. Ignored
                entirely if `position` is given explicitly.

        Returns:
            SectorSystemEntry: The newly created entry.
        """
        if position is None:
            position = self._random_position(star_system, min_separation_ly)

        entry = SectorSystemEntry(star_system, position, system_config=system_config)
        self.entries.append(entry)
        return entry

    def add_home_system(self, star_system, system_config=None,
                        jitter_ly=program_constants.DEFAULT_HOME_JITTER_LY):
        """
        Adds a system near the sector's center -- meant for the sector's
        first, "home" system, since a sector otherwise has no reason to favor
        any particular location. The position is randomized within
        `+/- jitter_ly` of dead-center on each axis (clamped to stay inside
        the sector) rather than placed exactly at the origin, since nothing
        about the origin itself is physically special.

        Args:
            star_system (StarSystem): The system to place.
            system_config (SystemConfig, optional): As in `add_system`.
            jitter_ly (float): The maximum offset from center, on each axis.

        Returns:
            SectorSystemEntry: The newly created entry.
        """
        half_edge = self.edge_ly / 2
        jitter_ly = min(jitter_ly, half_edge)
        position = tuple(_rng.uniform(-jitter_ly, jitter_ly) for _ in range(3))
        return self.add_system(star_system, position=position, system_config=system_config)

    def _fine_tune_position(self, position, candidate_system,
                            max_iterations=program_constants.SECTOR_GROWTH_FINE_TUNE_MAX_ITERATIONS):
        """
        Nudges `position` until `candidate_system`, placed there, clears
        every existing system's Hill sphere -- not necessarily just the
        point it was originally sampled around -- by repeatedly finding the
        worst violation (if any) and pushing directly away from it
        (`_nudge_away`), out to exactly `required_separation_ly`. Moving
        away from one neighbor can bring the position closer to (or newly
        inside) a different one, so this repeats until a position clear of
        every neighbor is found or `max_iterations` is exhausted. See
        "Growth (sample, then fine-tune)" in the module docstring.

        Every entry is re-checked on every iteration, not just whichever one
        is geometrically nearest by raw distance: `required_separation_ly`
        depends on each system's own Hill radius, so a farther-away but more
        massive system can be violated while a closer, smaller one isn't --
        checking only the nearest-by-distance system would miss that.

        A deficit within `program_constants.SECTOR_GROWTH_FLOATING_POINT_TOLERANCE_LY`
        of zero is treated as resolved rather than a real violation: nudging
        a position to exactly its required distance can leave a residual
        floating-point error of a few parts in 1e-15 from the sqrt/division
        in `_nudge_away`, which without this tolerance would register as a
        still-positive deficit forever, nudging the same point back to the
        same place every remaining iteration instead of ever converging.

        Args:
            position (tuple): The starting `(x, y, z)` position to fine-tune.
            candidate_system (StarSystem): The system being placed.
            max_iterations (int): The maximum number of nudges to attempt.

        Returns:
            tuple or None: The fine-tuned position, or None if no position
                           clear of every neighbor was found within
                           `max_iterations`.
        """
        tolerance = program_constants.SECTOR_GROWTH_FLOATING_POINT_TOLERANCE_LY

        for _ in range(max_iterations):
            worst_entry = None
            worst_required = None
            worst_deficit = tolerance

            for entry in self.entries:
                required = required_separation_ly(candidate_system, entry.star_system)
                deficit = required - distance_between(position, entry.position)
                if deficit > worst_deficit:
                    worst_deficit = deficit
                    worst_entry = entry
                    worst_required = required

            if worst_entry is None:
                return position

            position = _nudge_away(worst_entry.position, position, worst_required)

        return None

    def grow_from_seed(self, seed, system_factory, target_count=None,
                       k=program_constants.SECTOR_GROWTH_POISSON_DISK_K):
        """
        Grows the sector outward from an already-placed system -- see
        "Growth (sample, then fine-tune)" in the module docstring for the
        full two-phase mechanism (initial spherical-coordinate sample from
        the active parent, then fine-tuning against whichever system turns
        out to actually be nearest).

        Args:
            seed (SectorSystemEntry): An entry already present in
                `self.entries` to grow outward from -- typically the
                sector's home system, from `add_home_system`.
            system_factory (callable): A zero-argument callable returning a
                freshly generated `StarSystem` for each placement attempt.
                This module deliberately does not generate systems itself
                (it doesn't import `systemGen.py`'s CLI-oriented generation
                logic) -- the caller supplies it.
            target_count (int, optional): The sector's desired total system
                count once growth stops (including systems already
                present). If None, a Poisson-distributed count with mean
                `self.expected_system_count()` is drawn instead
                (`_sample_poisson_count`), so per-sector variance matches a
                real spatial Poisson process. If the drawn or given target
                is at or below the sector's current size, growth simply
                does nothing.
            k (int): Independent placement attempts per active point before
                giving up on it and retiring it from the active list.

        Returns:
            list: The newly created `SectorSystemEntry` instances, in
                 placement order. Does not include `seed`.

        Raises:
            ValueError: If `seed` is not already one of `self.entries`.
        """
        if seed not in self.entries:
            raise ValueError("seed must already be an entry in this sector (e.g. via add_home_system)")

        if target_count is None:
            target_count = _sample_poisson_count(self.expected_system_count())

        half_edge = self.edge_ly / 2
        active_list = [seed]
        new_entries = []

        while active_list and len(self.entries) < target_count:
            parent = _rng.choice(active_list)
            placed = False

            for _ in range(k):
                candidate_system = system_factory()

                inner_radius = hill_radius_ly(parent.star_system)
                outer_radius = max(mean_nearest_neighbor_ly(), inner_radius)
                position = _random_point_in_annulus(parent.position, inner_radius, outer_radius)

                position = self._fine_tune_position(position, candidate_system)
                if position is None:
                    continue
                if not all(-half_edge <= coordinate <= half_edge for coordinate in position):
                    continue

                entry = self.add_system(candidate_system, position=position)
                active_list.append(entry)
                new_entries.append(entry)
                placed = True
                break

            if not placed:
                active_list.remove(parent)

        return new_entries

    @staticmethod
    def distance_between(a, b):
        """
        Computes the distance, in light-years, between any two entries (or
        raw `(x, y, z)` positions) in the sector. See the module-level
        `distance_between`.
        """
        return distance_between(a, b)

    @staticmethod
    def required_separation_ly(system_a, system_b):
        """
        The minimum non-overlapping distance, in light-years, between two
        systems' Hill spheres. See the module-level `required_separation_ly`.
        """
        return required_separation_ly(system_a, system_b)

    @staticmethod
    def mean_nearest_neighbor_ly():
        """
        The mean nearest-neighbor distance, in light-years, at this
        module's target stellar density. See the module-level
        `mean_nearest_neighbor_ly`.
        """
        return mean_nearest_neighbor_ly()

    @staticmethod
    def classify_octant(position):
        """
        Classifies a raw `(x, y, z)` position into its "quadrant" label and
        positive magnitudes. See the module-level `classify_octant`.
        """
        return classify_octant(position)

    @staticmethod
    def format_named_location(position):
        """
        Formats a raw `(x, y, z)` position as a "named location" string. See
        the module-level `format_named_location`.
        """
        return format_named_location(position)

    def nearest_neighbors(self, entry, count=1):
        """
        Finds the systems closest to a given entry.

        Args:
            entry (SectorSystemEntry): The entry to measure distances from.
            count (int): The maximum number of neighbors to return.

        Returns:
            list: Up to `count` other `SectorSystemEntry` instances, nearest
                 first.
        """
        others = [other for other in self.entries if other is not entry]
        others.sort(key=entry.distance_to)
        return others[:count]

    def to_dict(self):
        """
        Returns a JSON-serializable dict of the whole sector.

        Returns:
            dict: With `name`, `edge_ly`, and `systems` (a list of each
                 entry's `to_dict()`) keys.
        """
        return {
            "name": self.name,
            "edge_ly": self.edge_ly,
            "systems": [entry.to_dict() for entry in self.entries],
        }

    def save(self, path):
        """
        Writes this sector to a JSON file.

        Args:
            path (str): The file path to write to.
        """
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, data):
        """
        Rebuilds a `SpaceSector` from `to_dict()`'s output.

        Two-tier behavior per entry, per the module docstring's
        "Persistence" section: when a `generated` key is present (Phase
        1.5's addition, written by every `to_dict()` since), the entry is
        rebuilt exactly via `StarSystem.from_dict` -- a faithful replay of
        the original, not a fresh roll. When it's absent (a file written
        before this addition, or a hand-authored one supplying only
        `config`), this falls back to the original behavior: rebuilding a
        fresh `StarSystem` from the stored `SystemConfig` recipe, freshly
        and independently randomized rather than byte-identical to the
        original.

        Args:
            data (dict): A dict as produced by `to_dict()`.

        Returns:
            SpaceSector: The reconstructed sector.
        """
        sector = cls(data["name"], edge_ly=data.get("edge_ly", program_constants.DEFAULT_SECTOR_EDGE_LY))

        for system_data in data["systems"]:
            if "generated" in system_data:
                star_system = StarSystem.from_dict(system_data["generated"])
                config = star_system.system_config
            else:
                config = SystemConfig.from_dict(system_data["config"])
                star_system = StarSystem(system_config=config)
            sector.entries.append(SectorSystemEntry(
                star_system, system_data["position"], system_config=config,
            ))

        return sector

    @classmethod
    def load(cls, path):
        """
        Reads a sector back from a JSON file written by `save()`.

        Args:
            path (str): The file path to read from.

        Returns:
            SpaceSector: The reconstructed sector.
        """
        with open(path, 'r') as f:
            return cls.from_dict(json.load(f))
