# planetgen/generation/binary.py

"""
The host of a pair's circumbinary planets
=========================================

This module defines the `PairHost` class (GEN.640 slice c, which replaced `BinaryStarProxy`), the two stars of an orbit-tree node standing as one
effective central star for the planets that circle both (P-type orbits): the system's own close pair, and each inner `close` node of a bigger tree
that holds circumbinary planets (GEN.332).

The host aggregates the mass, light, snow line, habitable zone and galactic place of its two stars. The pair's orbit is not kept on it: the host
holds the orbit node of its pair and `pair_orbit` reads the `PairOrbit` from it on each use, as `StarSystem.pair_orbit` does for the system's own
pair, so there is no second copy of the mutual orbit to keep in step with the node, and the stored `star_systems` row has no column for it.
"""

import math

from planetgen.generation import pair_orbit
from planetgen.generation.config import SystemConfig
from planetgen.physics import constants
from planetgen.util import detmath, draw
from planetgen.util.serialization import fields_from_dict
from planetgen.generation.star import Star
from planetgen.generation.mass_object import GALACTIC_MOTION_UNITS
from planetgen.generation.stellar_object import StellarObject
from planetgen.galaxy import enclosed_mass, scales
from planetgen.physics.orbits import calculate_habitable_zone, calculate_hill_sphere, minimum_update_interval_years
from planetgen.util.checks import finite_domain
from planetgen.physics import unit_registry as units

_km_to_au = units.converter("km", "au")
_ly_to_m = units.converter("ly", "m")
_m_to_pc = units.converter("m", "pc")
_m_to_au = units.converter("m", "au")

class PairHost(StellarObject):
    """
    The host of a pair's circumbinary planets: two stars of an orbit-tree node standing as one effective star.

    This class encapsulates two `Star` objects (primary and secondary) and
    exposes their combined properties (mass, luminosity, habitable zone, etc.)
    to the `StarSystem` class, so the planet generation logic treats the pair
    as a single gravitational and energetic source for circumbinary (P-type)
    planets. It draws nothing of the pair's orbit: the orbit is the node's.

    Attributes:
        _primary (Star): The primary (heavier) star of the pair.
        _secondary (Star): The secondary star of the pair.
        stars (list): A list containing both primary and secondary Star objects.
        _effective_mass (float): The combined mass of both stars in kilograms.
        _effective_luminosity (float): The combined luminosity of both stars in Watts.
        orbit_node (OrbitNode): The orbit node of the pair, the one place its orbit is kept (`pair_orbit` reads it).
        turned (bool): Whether the primary hangs under slot 1 of the node, which turns the orbit half way round (`pair_orbit.pair_orbit`).
    """

    FIELD_UNITS = {
        "name": "text", "type": "text", "temperature": "K", "radius": "km", "age": "Gyr", "lifespan": "Gyr",
        "habitable_zone": "au", "system_perimeter": "au", "heliosphere_radius": "au", **GALACTIC_MOTION_UNITS,
        "_effective_mass": "kg", "_effective_luminosity": "W",
    }
    host_node_no = None
    """int or None: The node of the system's orbit tree this host stands for when it hosts the circumbinary planets of an inner
    pair (GEN.332, `StarSystem._make_node_hosts`); `None` for the host of the system's own close pair, whose node the system
    knows (`StarSystem.pair_node_no`)."""

    SERIALIZABLE_FIELDS = list(FIELD_UNITS)
    """
    Mirrors `Star.SERIALIZABLE_FIELDS`, minus `yerkes_class` (not set on a
    host) and `mass`/`luminosity` (read-only properties -- their
    underscore-prefixed backing attributes, `_effective_mass`/
    `_effective_luminosity`, are listed instead, since `from_dict`'s plain
    `setattr` can't target a property name).
    `primary`/`secondary` are handled separately in `to_dict`/`from_dict`
    (nested full `Star` dicts), not via this list. The pair's orbit is not
    among them: it is the orbit node's, which the system stores once.
    """

    def __init__(self, system_config: SystemConfig, primary_star: Star, secondary_star: Star, orbit_node,
                 galactic_center_dist_ly=None, galactic_orbital_phase_deg=None, turned=False):
        """
        Initializes a PairHost object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object for the system.
            primary_star (Star): The primary star instance, the heavier of the two (the tree puts it in slot 0 of its node).
            secondary_star (Star): The secondary star instance.
            orbit_node (OrbitNode): The pair's node of the system's orbit tree: its axis, plane, phase and period are the pair's, and nothing is
                drawn for them (`pair_orbit` reads them from it).
            galactic_center_dist_ly (float, optional): This system's actual
                distance from the galactic center, in light-years, threaded
                into `_calculate_system_perimeter_static` -- see
                `Star.calculate_system_perimeter`'s docstring for the same
                parameter. `None` (the default) uses the radius that holds half of the galaxy's disk stars (`galaxy.scales`).
            galactic_orbital_phase_deg (float, optional): This system's
                current angular position around its galactic orbit, in
                degrees -- see `Star.__init__`'s docstring for the same
                parameter. `StarSystem.__init__` rolls this once and passes
                the identical value here and to both `primary_star`/
                `secondary_star`, so all three agree (a binary pair moves
                around the galaxy together, not independently). `None` (the
                default) rolls a fresh random value in `[0, 360)`.
            turned (bool): Whether the primary hangs under slot 1 of the node, which turns the node's orbit half way round (`pair_orbit.pair_orbit`).
        """
        # The host stands in for the whole pair: its own `.name` is what `store.insert_star_system` stores as `star_systems.name`, so it takes the
        # primary's own bare name -- what a wide (S-type) pair's system name already is.
        super().__init__(system_config, name=primary_star.name)

        self._primary = primary_star
        self._secondary = secondary_star
        self.orbit_node = orbit_node
        self.turned = turned

        self.age = max(self._primary.age, self._secondary.age)
        self.lifespan = max(self._primary.lifespan, self._secondary.lifespan)
        # Both formed in one cloud: the pair's disk is the primary's (the planet cull and the solids read it).
        self.metallicity_fe_h = self._primary.metallicity_fe_h
        # GEN.458: the two stars' lights add (luminosity, so the squares of their snow lines do: the line follows
        # the square root of the luminosity), each with the scatter it drew.
        self.snow_line_au = detmath.hypot(self._primary.snow_line_au, self._secondary.snow_line_au)

        self.stars = [self._primary, self._secondary]

        # Calculate combined properties
        self._effective_mass = sum(s.mass for s in self.stars)
        self._effective_luminosity = sum(s.luminosity for s in self.stars)

        lead = primary_star.leaf_stars[0] if isinstance(primary_star, PairHost) else primary_star    # a pair of pairs (GEN.332) goes by its heaviest star
        self.type = f"Binary ({primary_star.type.split(' ')[0]}/{secondary_star.type.split(' ')[0]})" # Simplified type
        self.temperature = (primary_star.temperature + secondary_star.temperature) / 2 # Simple average
        self.radius = max(primary_star.radius, secondary_star.radius) # Use the larger radius for approximation
        # Note: self.lifespan was already set above to max(primary, secondary) lifespan;
        # do not overwrite it here with the primary-only value.

        # Recalculate derived properties based on effective values
        # The zone follows the colour of the light: the pair's temperature weighted by each star's luminosity (GEN.146),
        # not the plain average `self.temperature` keeps.
        light_primary, light_secondary = primary_star.luminosity, secondary_star.luminosity
        light_temperature = ((light_primary * primary_star.temperature + light_secondary * secondary_star.temperature)
                             / (light_primary + light_secondary))
        self.habitable_zone = calculate_habitable_zone(self._effective_luminosity, light_temperature)
        self.system_perimeter = self._calculate_system_perimeter_static(
            self._effective_mass, galactic_center_dist_ly
        )
        # Galactic orbit doesn't depend on mass (see
        # `Star.calculate_galactic_orbit`'s docstring), so the inherited
        # instance method is reused directly rather than needing a static
        # variant the way `_calculate_system_perimeter_static` does.
        self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy = \
            self.calculate_galactic_orbit(galactic_center_dist_ly)
        self.galactic_orbital_phase_deg = (
            galactic_orbital_phase_deg if galactic_orbital_phase_deg is not None else draw.uniform(0, 360)
        )
        self.galactic_min_update_interval_years = \
            minimum_update_interval_years(self.galactic_orbital_period_gy * 1e9)
        # For heliosphere, we need an effective radius and type for the static method.
        # This is a simplification, as binary heliospheres are complex.
        self.heliosphere_radius = Star._calculate_heliosphere_radius_static(
            self._effective_mass,
            self._effective_luminosity,
            self.radius, # Using the larger radius as an approximation
            lead.type, # Using primary's type for wind model selection
            lead.yerkes_class, # Using primary's yerkes class for wind model selection
            self.log_lx_lbol # The pair's combined X-ray share, which a cool dwarf's coronal wind is made from (GEN.1625)
        )

    @property
    def leaf_stars(self):
        """list of Star: The stars the host is made of, in order, a pair of pairs (an inner node holding a sub-node, GEN.332) flattened."""
        return [leaf for part in self.stars for leaf in (part.leaf_stars if isinstance(part, PairHost) else (part,))]

    @property
    def mass(self):
        """Returns the combined mass of the binary system."""
        return self._effective_mass

    @property
    def luminosity(self):
        """Returns the combined luminosity of the binary system."""
        return self._effective_luminosity

    def _summed(self, field):
        values = [getattr(star, field, None) for star in self.stars]
        return None if all(value is None for value in values) else sum(value or 0.0 for value in values)

    @property
    def l_xuv_w(self):
        """Both stars' XUV, W (GEN.86): a circumbinary planet gets both."""
        return self._summed("l_xuv_w")

    @property
    def xuv_fluence_j(self):
        """Both stars' XUV over their lives, J (GEN.86)."""
        return self._summed("xuv_fluence_j")

    @property
    def lethal_event_rate_per_gyr(self):
        """The pair's galactic hazard (GEN.87): both stars share one place."""
        return getattr(self.stars[0], "lethal_event_rate_per_gyr", None)

    @property
    def flare_n33_per_yr(self):
        """Both stars' flares above 1e33 erg a year (GEN.86)."""
        return self._summed("flare_n33_per_yr")

    @property
    def log_lx_lbol(self):
        """The pair's combined `L_X/L_bol`, log10 (GEN.86): what drives
        the wind a circumbinary planet sits in."""
        l_x = sum(10 ** star.log_lx_lbol * star.luminosity for star in self.stars
                  if getattr(star, "log_lx_lbol", None) is not None)
        return math.log10(l_x / self.luminosity) if l_x > 0.0 and self.luminosity else None

    @property
    def pair_orbit(self):
        """The orbit of the pair (`pair_orbit.PairOrbit`), read from its node on each use: the axis, plane, period, phase and each star's place about the centre of mass."""
        return pair_orbit.pair_orbit(self.orbit_node, self._primary.mass, self._secondary.mass, self.turned)

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this host's snapshot state, per
        `SERIALIZABLE_FIELDS`, plus nested `primary`/`secondary` star dicts.

        The combined/derived fields (`_effective_mass`,
        `_effective_luminosity`, `habitable_zone`, `system_perimeter`,
        `heliosphere_radius`, etc.) are stored exactly as generated rather
        than recomputed from `primary`/`secondary` on load -- they were
        computed once, at generation time, from whichever `random`
        rolls happened to occur, and recomputing them later from the
        (faithfully reloaded) constituent stars would only risk drift if
        the derivation formulas ever change between save and load. The pair's
        orbit is the system's `orbit_nodes` entry, written once with the system.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`, plus
                 `primary`/`secondary` (each a full `Star.to_dict()`).
        """
        data = super().to_dict()
        data["primary"] = self._primary.to_dict()
        data["secondary"] = self._secondary.to_dict()
        return data

    @classmethod
    def from_dict(cls, data, system_config, orbit_node, turned=False):
        """
        Reconstructs a `PairHost` from a dict in the shape
        `to_dict()` produces, without re-running generation (`__init__` is
        bypassed via `object.__new__`) and without recomputing any derived
        field -- every stored value is assigned back as-is.

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The system's shared config,
                                          threaded into both `primary`'s and
                                          `secondary`'s `Star.from_dict`
                                          calls as well as the host itself
                                          (see `StarSystem.from_dict`).
            orbit_node (OrbitNode): The pair's node of the system's stored orbit tree.
            turned (bool): Whether the primary hangs under slot 1 of the node (`__init__`).

        Returns:
            PairHost: The reconstructed host.
        """
        host = object.__new__(cls)
        host.system_config = system_config
        primary = Star.from_dict(data["primary"], system_config)
        secondary = Star.from_dict(data["secondary"], system_config)
        host._primary = primary
        host._secondary = secondary
        host.stars = [primary, secondary]
        host.orbit_node = orbit_node
        host.turned = turned
        fields_from_dict(host, data, cls.SERIALIZABLE_FIELDS)
        host.snow_line_au = detmath.hypot(primary.snow_line_au, secondary.snow_line_au)  # derived from the stored pair, as in `__init__`
        host.habitable_zone = tuple(data["habitable_zone"])
        host.lifespan = float('inf') if data["lifespan"] is None else data["lifespan"]
        return host

    @staticmethod
    @finite_domain()
    def _calculate_system_perimeter_static(mass, galactic_center_dist_ly=None):
        """
        Calculates the system perimeter (Hill sphere relative to the galaxy)
        for an arbitrary mass, mirroring `Star.calculate_system_perimeter` but
        usable without a bound `Star` instance.

        Args:
            mass (float): The mass to use for the calculation, in kilograms
                          (typically the binary system's combined mass).
            galactic_center_dist_ly (float, optional): This system's actual
                distance from the galactic center, in light-years -- see
                `Star.calculate_system_perimeter`'s docstring. `None` (the
                default) uses the radius that holds half of the galaxy's disk stars (`galaxy.scales`).

        Returns:
            float: The radius of the Hill sphere in Astronomical Units (AU).
        """
        if galactic_center_dist_ly is None:
            galactic_center_dist_ly = scales.default_galactocentric_ly()
        galactic_center_dist_m = _ly_to_m(max(galactic_center_dist_ly, constants.MIN_TIDAL_DISTANCE_LY))
        hill_radius_m = calculate_hill_sphere(galactic_center_dist_m, mass, enclosed_mass.current_enclosed_mass_kg(_m_to_pc(galactic_center_dist_m)))
        return _m_to_au(hill_radius_m)

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)
object_timing.track(PairHost, "binary-star")
