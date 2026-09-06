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

Persistence
-----------
A sector can be saved to (and reloaded from) JSON: each entry keeps the
`SystemConfig` its `StarSystem` was built from, so the sector's *recipe* --
star types, forced features, names, orbit counts, and so on -- survives a
save/reload cycle. Reloading does NOT reproduce byte-identical systems: a
fair amount of the generation in this package (planet class rolls, moon
coin-flips, name generation, flavor text) is drawn from the `secrets` module
specifically because it should be unpredictable, so there's no seed that
could play it back exactly. A reloaded system is a fresh, independently
random system built from the same recipe -- not the same system come back.
The one detail pinned across reload is the star's name: if a system's config
didn't already force one, the name it was actually given gets baked into the
saved config so it doesn't change out from under a saved sector.
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

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this entry.

        The saved config's `name` is pinned to the star's actual generated
        name when the config itself didn't already force one, so a reloaded
        sector keeps the same names even though other details are freshly
        randomized (see the module docstring).

        Returns:
            dict: With `name`, `position`, and `config` keys.
        """
        config_dict = self.system_config.to_dict()
        if config_dict.get("name") is None:
            config_dict["name"] = self.star_system.star.name

        return {
            "name": self.star_system.star.name,
            "position": list(self.position),
            "config": config_dict,
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
        Rebuilds a `SpaceSector` from `to_dict()`'s output, rebuilding each
        system from its stored `SystemConfig`. Rebuilt systems are freshly,
        independently randomized from that same recipe -- not byte-identical
        to the originals -- see the module docstring for why.

        Args:
            data (dict): A dict as produced by `to_dict()`.

        Returns:
            SpaceSector: The reconstructed sector.
        """
        sector = cls(data["name"], edge_ly=data.get("edge_ly", program_constants.DEFAULT_SECTOR_EDGE_LY))

        for system_data in data["systems"]:
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
