# planetgen/generation/mass_object.py

"""
Mass Object
===========

The core every body with a mass shares (GEN.565, the class the unified
star-like objects design calls `MassObject`): its motion and the one `to_dict`/`from_dict` pair that each family used
to write out by hand.

Motion (Boss, 5:20 PM EDT on 2026-10-18: every object with mass follows the same rules) is a state relative to the
body's parent: the position and velocity (`spatial`, a `SpatialPosition3D` whose "system" frame is the offset from the
primary: the star, or the planet for a moon) and the pull of the primary (`orbit_mu_au3_per_year2`, which each family
sets), moved by `advance` on `kepler.universal_step`; or, for a body that goes round the galaxy, the circular
galactic-orbit fields (`set_galactic_motion`), moved by `advance` along the orbit's phase. One function moves a planet,
a moon, a comet, a star and a wanderer.

A family class (`Star`, `Planet`, `Comet`, `RoguePlanet`, ...) lists its
own `SERIALIZABLE_FIELDS` and, where a value is not plain JSON, overrides
`_encode` (object to dict) and `_decode` (dict back onto the object). Nothing
here draws a random number except `set_galactic_motion`, which draws what
the older copies drew, in the same place in each constructor.
"""

import math

from planetgen.galaxy.galactic_orbit import generate_galactic_orbit_fields
from planetgen.physics import constants as physical_constants
from planetgen.physics import state_vectors
from planetgen.physics.kepler import universal_step
from planetgen.physics.position import SpatialPosition3D
from planetgen.util.serialization import fields_from_dict, fields_to_dict

GALACTIC_MOTION_FIELDS = [
    "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
    "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
]
"""The four fields `set_galactic_motion` sets: the circular speed and period about the galactic centre, the phase
along that orbit, and the shortest interval between orbital updates."""


class MassObject:
    """The core of a body with a mass: galactic motion and serialization from a field list. A family whose name is
    a salad drawn on first read lists `LazySaladName` before it."""

    SERIALIZABLE_FIELDS = []
    """The attributes a family persists; each family sets its own."""

    spatial = None
    """SpatialPosition3D or None: where the body is, with its "system" frame the offset from its primary (a family
    that holds one declares its axes with `position.axis_property` and `velocity_axis_property`)."""
    _staged_au = None
    _staged_velocity_kms = None

    # --- Motion relative to a primary --------------------------------------------------------------------------

    def set_position_au(self, x, y, z):
        """Puts the body at `(x, y, z)` AU from its primary; creates its
        position (anchors at the origin until it is placed in a sector)
        on first use, with the mass and mu it has by then."""
        if self.spatial is None:
            self.spatial = SpatialPosition3D((x, y, z), (0.0, 0.0, 0.0), star_center_galactic=(0.0, 0.0, 0.0),
                                             length_unit_m=physical_constants.AU_M)
        else:
            self.spatial.set_system_cartesian(x, y, z)
        self._staged_au = None
        mass = getattr(self, "mass", None)
        if isinstance(mass, (int, float)) and not isinstance(mass, bool) and mass >= 0:
            self.spatial.set_mass(mass)
        if self._staged_velocity_kms is not None and all(v is not None for v in self._staged_velocity_kms):
            self.set_velocity_kms(*self._staged_velocity_kms)

    def set_velocity_kms(self, vx, vy, vz):
        """Sets the velocity relative to the primary, km/s. A body with no
        position yet keeps it until it has one."""
        if self.spatial is None:
            self._staged_velocity_kms = (vx, vy, vz)
            return
        self.spatial.set_velocity_cartesian(vx * 1000.0, vy * 1000.0, vz * 1000.0, frame="system")
        self._staged_velocity_kms = None

    def orbit_mu_au3_per_year2(self):
        """The gravitational parameter, AU^3/yr^2, of what the body orbits;
        each family that holds an orbit position defines it."""
        raise NotImplementedError

    def orbit_from_vector(self):
        """
        The orbit the body is on right now, worked out from its position
        and velocity relative to its primary (`state_vectors.elements_from_state`,
        AU and years): periapsis distance, eccentricity, inclination,
        ascending node, argument of periapsis, true anomaly (radians),
        semi-major axis, period, `kind` and energy. Never stored, so it
        cannot fall out of step with the vector: anything that moves the
        body or changes its velocity changes this orbit with it.

        Raises:
            ValueError: With no position or velocity yet, or for a body
                moving straight at or away from its primary.
        """
        if self.spatial is None or self.spatial.get_velocity_vector("system") is None:
            raise ValueError("the body has no position and velocity yet")
        to_au_per_year = physical_constants.SECONDS_PER_YEAR / physical_constants.AU_M
        position = self.spatial.get_coordinates("system", "cartesian")
        velocity = tuple(v * to_au_per_year for v in self.spatial.get_velocity_vector("system"))
        return state_vectors.elements_from_state(position, velocity, self.orbit_mu_au3_per_year2())

    def projected_orbit_au(self, count=128):
        """
        `count` points, AU from the primary, round the closed ellipse the
        body is on now (`orbit_from_vector`), starting at periapsis: the
        course it follows when nothing disturbs it. Wobble and other
        perturbations in the vector are not part of it.

        Raises:
            ValueError: For an orbit that does not close (parabolic or
                hyperbolic), or as `orbit_from_vector`.
        """
        orbit = self.orbit_from_vector()
        if orbit["kind"] != "elliptical":
            raise ValueError(f"a {orbit['kind']} orbit does not close")
        return state_vectors.closed_orbit_points(
            orbit["periapsis_distance"], orbit["eccentricity"], orbit["inclination"], orbit["ascending_node"],
            orbit["arg_periapsis"], count)

    def advance(self, seconds):
        """
        Moves the body `seconds` along its trajectory (negative: back). A body with a position and velocity relative
        to a primary is stepped on that state by `kepler.universal_step` with the pull of the primary
        (`orbit_mu_au3_per_year2`), so an orbit another body's pull has changed is carried on as it now is; a body
        that goes round the galaxy on a circular orbit has its phase advanced by the orbit's period, the arithmetic
        `store.advance_galactic_positions` does in bulk.

        Raises:
            ValueError: For a body with neither a state vector nor a galactic orbit.
        """
        years = seconds / physical_constants.SECONDS_PER_YEAR
        if self.spatial is not None and self.spatial.get_velocity_vector("system") is not None:
            to_au_per_year = physical_constants.SECONDS_PER_YEAR / physical_constants.AU_M
            position = self.spatial.get_coordinates("system", "cartesian")
            velocity = tuple(v * to_au_per_year for v in self.spatial.get_velocity_vector("system"))
            position, velocity = universal_step(position, velocity, self.orbit_mu_au3_per_year2(), years)
            self.set_position_au(*position)
            self.set_velocity_kms(*(v / to_au_per_year / 1000.0 for v in velocity))
            return
        period_gy = getattr(self, "galactic_orbital_period_gy", None)
        if period_gy:
            period_years = period_gy * 1.0e9
            turn = math.fmod(years, period_years) / period_years
            self.galactic_orbital_phase_deg = (self.galactic_orbital_phase_deg + 360.0 * turn) % 360.0
            return
        raise ValueError("the body has no state vector and no galactic orbit to move along")

    def set_galactic_motion(self, galactic_center_dist_ly=None, galactic_orbital_phase_deg=None):
        """Sets the four `GALACTIC_MOTION_FIELDS` from the body's distance from the galactic centre (none: the
        generator's default distance) and a given phase (none: drawn)."""
        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields(galactic_center_dist_ly, galactic_orbital_phase_deg)

    def to_dict(self):
        """A JSON-serializable dict of the body's `SERIALIZABLE_FIELDS`, with `_encode` applied."""
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        self._encode(data)
        return data

    def _encode(self, data):
        """Turns the values of `data` that are not plain JSON into ones that are; the family overrides."""

    @classmethod
    def _rebuild(cls, data, system_config, **context):
        """The body for `data`, without running `__init__` (that would draw it again): its shared config, any
        `context` attributes (what a family threads in rather than stores, such as a planet's star), the
        `SERIALIZABLE_FIELDS` set, and `_decode` applied."""
        body = object.__new__(cls)
        body.system_config = system_config
        body.__dict__.update(context)
        fields_from_dict(body, data, cls.SERIALIZABLE_FIELDS)
        body._decode(data)
        return body

    def _decode(self, data):
        """Turns the values `_encode` changed back; the family overrides."""

    @classmethod
    def from_dict(cls, data, system_config):
        """Reconstructs the body from a dict in the shape `to_dict` produces."""
        return cls._rebuild(data, system_config)
