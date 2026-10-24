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

import dataclasses
import math

from planetgen.galaxy.galactic_orbit import generate_galactic_orbit_fields
from planetgen.physics import constants as physical_constants
from planetgen.physics import spin, state_vectors, unit_registry
from planetgen.physics.kepler import universal_step
from planetgen.physics.position import SpatialPosition3D, update_interval_s
from planetgen.physics.quantity import Quantity
from planetgen.util.serialization import fields_from_dict, fields_to_dict

_KM_S = unit_registry.unit("km_s").factor
"""float: Metres a second in a kilometre a second."""

PLAIN = ("text", "ratio", "count", "index", "flag", "dex", "unit-vector", "id")
"""tuple: What a field declares in place of a unit when it holds no quantity: a name or a label (`text`), a plain
number (`ratio`), a whole number (`count`, `index`, `id`), a yes or no (`flag`), a logarithmic scale (`dex`) or one
component of a unit vector. The same words the column registry (`db.column_units`) uses for a column with no unit."""

GALACTIC_MOTION_UNITS = {
    "galactic_orbital_speed_kms": "km_s", "galactic_orbital_period_gy": "Gyr",
    "galactic_orbital_phase_deg": "deg", "galactic_min_update_interval_years": "yr",
}
"""dict: The four fields `set_galactic_motion` sets, with their native units: the circular speed and period about the
galactic centre, the phase along that orbit, and the shortest interval between orbital updates."""

GALACTIC_MOTION_FIELDS = list(GALACTIC_MOTION_UNITS)
"""list: The names of `GALACTIC_MOTION_UNITS`."""


@dataclasses.dataclass(frozen=True)
class KinematicRecord:
    """The kinematic record of a mass object (GEN.565, card S16 of docs/design/architecture/unified-system-barycentre.md): what the
    stored point of every body holds, all relative to its parent in the shared non-rotating (galactic) axes. A field is
    `None` where the body has not got it yet (no position, an object saved before the spin was drawn)."""

    position_au: tuple = None
    """tuple or None: Where the body is, AU from its primary."""
    velocity_kms: tuple = None
    """tuple or None: Its velocity relative to its primary, km/s."""
    spin_axis: tuple = None
    """tuple or None: The unit vector of its angular velocity."""
    spin_rate_rad_s: float = None
    """float or None: Its rotation rate, rad/s; zero when it is not rotating."""
    mass_kg: float = None
    """float or None: Its mass, kg."""
    mu_m3_s2: float = None
    """float or None: G times its mass, m^3/s^2."""
    epoch_unix: float = None
    """float or None: When the position and velocity hold, Unix seconds."""
    next_due_unix: float = None
    """float or None: When the orbit update next has to move it, Unix seconds."""


class MassObject:
    """The core of a body with a mass: galactic motion and serialization from a field list. A family whose name is
    a salad drawn on first read lists `LazySaladName` before it."""

    FIELD_UNITS = {}
    """dict: Each attribute a family persists, with the unit the instance holds it in (GEN.570): a key of
    `unit_registry.UNITS` for a quantity (a plain float in that native unit, so a kernel pays for nothing), or a word of
    `PLAIN` for a field that holds none. A family sets its own, in the order it persists them; `quantity` builds the
    `Quantity` of one at the edges, and `tests/test_field_units.py` fails when a persisted field is not declared here or
    is declared in another kind than its column."""

    SERIALIZABLE_FIELDS = []
    """The attributes a family persists; each family sets its own (the names of its `FIELD_UNITS`)."""

    spatial = None
    """SpatialPosition3D or None: where the body is, with its "system" frame the offset from its primary (a family
    that holds one declares its axes with `position.axis_property` and `velocity_axis_property`)."""
    _staged_au = None
    _staged_velocity_kms = None

    # --- Quantities ----------------------------------------------------------------------------------------------

    @classmethod
    def field_unit(cls, name):
        """The unit key or `PLAIN` word `FIELD_UNITS` declares for the attribute `name`; a `KeyError` for an
        undeclared one."""
        return cls.FIELD_UNITS[name]

    def quantity(self, name):
        """
        The attribute `name` as a `Quantity` in its declared native unit (a pair, such as the habitable zone, as a
        tuple of them), every other unit of its kind one read away: `star.quantity("mass").Msun`.

        Raises:
            KeyError: For an attribute `FIELD_UNITS` does not declare.
            ValueError: For one that holds no quantity (`PLAIN`), or that has no value yet.
        """
        unit = self.FIELD_UNITS[name]
        if unit in PLAIN:
            raise ValueError(f"{name} holds no quantity (it is {unit})")
        value = getattr(self, name)
        if value is None:
            raise ValueError(f"{name} has no value yet")
        kind = unit_registry.unit(unit).kind
        if isinstance(value, (tuple, list)):
            return tuple(Quantity(kind, part, unit) for part in value)
        return Quantity(kind, value, unit)

    # --- Motion relative to a primary --------------------------------------------------------------------------

    def set_position_au(self, x, y, z):
        """Puts the body at `(x, y, z)` AU from its primary; creates its
        position (anchors at the origin until it is placed in a sector)
        on first use, with the mass and mu it has by then."""
        if self.spatial is None:
            self.spatial = SpatialPosition3D((x, y, z), (0.0, 0.0, 0.0), star_center_galactic=(0.0, 0.0, 0.0),
                                             length_unit="au")
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
        self.spatial.set_velocity_cartesian(vx * _KM_S, vy * _KM_S, vz * _KM_S, frame="system")
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
            self.set_velocity_kms(*(v / to_au_per_year / _KM_S for v in velocity))
            return
        period_gy = getattr(self, "galactic_orbital_period_gy", None)
        if period_gy:
            period_years = period_gy * 1.0e9
            turn = math.fmod(years, period_years) / period_years
            self.galactic_orbital_phase_deg = (self.galactic_orbital_phase_deg + 360.0 * turn) % 360.0
            return
        raise ValueError("the body has no state vector and no galactic orbit to move along")

    @property
    def spin_axis(self):
        """`(x, y, z)` unit vector of the body's angular velocity, or `None` for a body with no spin drawn."""
        axis = spin.spin_values(self)[:3]
        return None if any(value is None for value in axis) else axis

    @property
    def spin_rate_rad_s(self):
        """The body's rotation rate, rad/s (`2 pi` over its rotation period), or `None` for a body with no period."""
        period_hours = getattr(self, "rotation_period_hours", None)
        if not period_hours:
            return None
        return 2.0 * math.pi / (period_hours * 3600.0)

    def speed_kms(self):
        """The body's speed, km/s: the length of its velocity relative to its primary, or, for a body that goes round the
        galaxy with no such vector, its circular galactic speed; `None` with neither."""
        if self.spatial is not None and self.spatial.get_velocity_vector("system") is not None:
            return self.spatial.get_speed("system") / _KM_S
        return getattr(self, "galactic_orbital_speed_kms", None)

    def axial_tilt_from_state(self):
        """
        The angle, degrees, between the spin axis and the normal of the orbit about the primary (`r x v`), or, for a
        body with no state vector, the galaxy's pole (`physics.spin`): what `axial_tilt_deg` stores, worked out from the
        axis and the motion. `None` for a body with no spin axis.
        """
        axis = self.spin_axis
        if axis is None:
            return None
        normal = spin.GALACTIC_POLE
        if self.spatial is not None and self.spatial.get_velocity_vector("system") is not None:
            normal = spin.orbit_normal_from_state(self.spatial.get_coordinates("system", "cartesian"),
                                                  self.spatial.get_velocity_vector("system"))
        return spin.obliquity_between(axis, normal)

    def kinematic_record(self, scale="system"):
        """
        The body's `KinematicRecord`, read from its point and its spin. `scale` ("galactic", "system" or "planetary",
        `position.THRESHOLDS_M`) picks the distance the next-due time is the time to move.
        """
        position = velocity = mass = mu = epoch = next_due = None
        point = self.spatial
        if point is not None:
            position = point.get_coordinates("system", "cartesian")
            if point.get_velocity_vector("system") is not None:
                velocity = tuple(v / _KM_S for v in point.get_velocity_vector("system"))
            mass, mu, epoch = point.mass_kg, point.mu, point.epoch_unix
        speed = self.speed_kms()
        if epoch is not None and speed is not None:
            next_due = epoch + update_interval_s(speed * _KM_S, scale)
        return KinematicRecord(position, velocity, self.spin_axis, self.spin_rate_rad_s, mass, mu, epoch, next_due)

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
