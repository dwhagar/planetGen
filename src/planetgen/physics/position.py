# planetgen/physics/position.py

"""
One point in space (GEN.74)
===========================

`SpatialPosition3D` keeps where one body is in three frames at once, each
in Cartesian, cylindrical and spherical form, and keeps all of them in step:
setting any coordinate, in any frame and any form, recomputes the rest.

Frames
    galactic   the galaxy's own axes, origin at the galactic center.
    sector     offset from the center of the body's sector.
    system     offset from the center of the body's nearest star (absent
               for a star itself and for a body with no star reference).

Forms
    cartesian    (x, y, z)
    cylindrical  (r, theta, z): theta is the azimuth, in radians, from +x
                 toward +y, in (-pi, pi].
    spherical    (r, theta, phi): theta the same azimuth, phi the polar
                 angle from +z, in [0, pi].

Lengths are in `length_unit_m` metres (default 1: metres), so an object
kept in light-years or AU holds its own numbers exactly, with no round trip
through metres; speeds are always metres a second and masses kilograms. The
storage layer converts (km, mpc) at its edge (`planetgen.physics.units`).

Precision: the position that was last set is kept as given, in the frame it
was given in, and the other frames are derived from it. A moon placed by its
system offset therefore keeps its metres exactly even though the same point
in galactic coordinates (about 1e20 m from the center) cannot.

Sector address: with the sector edge known (`sector_edge_pc`), the position
also knows which cell of the galaxy's sector grid it is in (`sector_address`:
ring, layer, slot; `planetgen.galaxy.geometry`). It is worked out from the
galactic position whenever the position changes, by any setter or anchor
move, so it is never stale. `set_sector_address` moves the body to another
cell the way `carry_sector_center` does, keeping its offset from the
sector's center.

Velocity: one vector in two frames. The galactic velocity (the sector frame
has the same axes and does not move) is the velocity in the galaxy's axes;
the system velocity is relative to the nearest star, so a planet keeps its
30 km/s round its star separate from the star's own 220 km/s round the
galaxy. The velocity that was last set is kept as given, in its frame, and the
other is derived: galactic = the star's velocity (`star_velocity`) + system.
`epoch_unix` is when the position and velocity hold (Unix seconds; `None`
when unknown), so a stored position and velocity can be advanced to another
time.

Observable movement: how long a body at its speed takes to move the distance
that has to be re-stored (`THRESHOLDS_M`, from docs/design/orbital-updates.md
section 3), capped at `MAX_UPDATE_INTERVAL_S`. A body at rest never comes
due.

Limits (where the maths stops holding): a coordinate that is not finite, a
negative radius, a polar angle outside [0, pi], and a speed at or above the
speed of light are refused with `ValueError`; at the origin and on the polar
axis the angles are not defined and read as 0.
"""

import math

from planetgen.galaxy import geometry
from planetgen.physics import constants as physical_constants

FRAMES = ("galactic", "sector", "system")
"""tuple: The frames a position is kept in."""

FORMS = ("cartesian", "cylindrical", "spherical")
"""tuple: The coordinate forms each frame is kept in."""

VELOCITY_FRAMES = ("galactic", "system")
"""tuple: The frames a velocity is kept in; "sector" reads as "galactic"."""

THRESHOLDS_M = {
    "galactic": 0.01 * physical_constants.AU_PER_MILLIPARSEC * physical_constants.AU_M,
    "system": 0.01 * physical_constants.AU_M,
    "planetary": 1.0e8,
}
"""dict: The least move worth storing for each scale of object, in metres:
0.01 mpc (about 2 AU) for stars and other bodies outside a system, 0.01 AU
for planets and companion stars, 100,000 km for moons."""

MAX_UPDATE_INTERVAL_S = 1.0e9 * 365.25 * 86400.0
"""float: The longest a next-due time runs, a billion years in seconds, so a
very slow body still gets looked at."""

SPEED_OF_LIGHT_MS = 299_792_458.0
"""float: The speed no body reaches, m/s."""


def _finite(values, what):
    for value in values:
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{what} must be finite numbers, got {tuple(values)!r}")


def cartesian_to_cylindrical(x, y, z):
    """(x, y, z) -> (r, theta, z); theta is 0 on the axis."""
    return (math.hypot(x, y), math.atan2(y, x), z)


def cartesian_to_spherical(x, y, z):
    """(x, y, z) -> (r, theta, phi); both angles are 0 at the origin and phi
    is 0 (or pi) on the polar axis."""
    r = math.sqrt(x * x + y * y + z * z)
    if r == 0.0:
        return (0.0, 0.0, 0.0)
    return (r, math.atan2(y, x), math.acos(max(-1.0, min(1.0, z / r))))


def cylindrical_to_cartesian(r, theta, z):
    """(r, theta, z) -> (x, y, z). Raises `ValueError` for a negative `r`."""
    _finite((r, theta, z), "cylindrical coordinates")
    if r < 0.0:
        raise ValueError(f"cylindrical radius must not be negative, got {r!r}")
    return (r * math.cos(theta), r * math.sin(theta), z)


def spherical_to_cartesian(r, theta, phi):
    """(r, theta, phi) -> (x, y, z). Raises `ValueError` for a negative `r`
    or a `phi` outside [0, pi]."""
    _finite((r, theta, phi), "spherical coordinates")
    if r < 0.0:
        raise ValueError(f"spherical radius must not be negative, got {r!r}")
    if not 0.0 <= phi <= math.pi:
        raise ValueError(f"spherical polar angle must be in [0, pi], got {phi!r}")
    return (r * math.sin(phi) * math.cos(theta), r * math.sin(phi) * math.sin(theta), r * math.cos(phi))


def mu_of(mass_kg):
    """The gravitational parameter G * M (m^3/s^2) of `mass_kg`."""
    return physical_constants.G * mass_kg


def mass_of(mu):
    """The mass (kg) whose gravitational parameter is `mu`."""
    return mu / physical_constants.G


class SpatialPosition3D:
    """
    One body's position, velocity and point-mass data, with every frame and
    form kept in step.

    Args:
        galactic_cartesian (tuple): Galactic (x, y, z), metres.
        sector_center_galactic (tuple): The sector center, galactic metres.
        velocity_vector_cartesian (tuple): Velocity, galactic axes, m/s.
            (Set a velocity relative to the star with
            `set_velocity_cartesian(..., frame="system")`.)
        star_center_galactic (tuple or None): The nearest star's center,
            galactic metres; `None` for a body with no star reference.
        is_star (bool): A star has no system frame.
        mass_kg (float or None): The body's mass; its `mu` follows.
        length_unit_m (float): Metres in one unit of every coordinate and
            anchor (1 for metres, `LY_TO_M` for light-years, `AU_M` for AU).
        star_velocity_galactic (tuple): The nearest star's velocity, galactic
            axes, m/s: what the system frame's velocity is measured from.
        epoch_unix (float or None): When the position and velocity hold.
        sector_edge_pc (float or None): The sector grid's edge length, parsecs.
            Given, the position keeps its `sector_address`; `None` leaves
            the address unknown.
    """

    def __init__(self, galactic_cartesian, sector_center_galactic, velocity_vector_cartesian=(0.0, 0.0, 0.0),
                 star_center_galactic=None, is_star=False, mass_kg=None, length_unit_m=1.0,
                 star_velocity_galactic=(0.0, 0.0, 0.0), epoch_unix=None, sector_edge_pc=None):
        _finite((length_unit_m,), "length unit")
        if length_unit_m <= 0.0:
            raise ValueError(f"the length unit must be positive, got {length_unit_m!r}")
        self.length_unit_m = float(length_unit_m)
        _finite(galactic_cartesian, "galactic position")
        _finite(sector_center_galactic, "sector center")
        if star_center_galactic is not None:
            _finite(star_center_galactic, "star center")
        self._sector_center = tuple(float(v) for v in sector_center_galactic)
        self._star_center = None if star_center_galactic is None else tuple(float(v) for v in star_center_galactic)
        self._is_star = bool(is_star)
        _finite(star_velocity_galactic, "star velocity")
        self._star_velocity = tuple(float(v) for v in star_velocity_galactic)
        self._epoch_unix = None
        self.set_epoch_unix(epoch_unix)
        if sector_edge_pc is not None:
            _finite((sector_edge_pc,), "sector edge")
            if sector_edge_pc <= 0.0:
                raise ValueError(f"the sector edge must be positive, got {sector_edge_pc!r}")
            sector_edge_pc = float(sector_edge_pc)
        self._sector_edge_pc = sector_edge_pc
        self._sector_address = None
        self._truth = ("galactic", tuple(float(v) for v in galactic_cartesian))
        self._mass_kg = None
        self._mu = None
        if mass_kg is not None:
            self.set_mass(mass_kg)
        self._velocity_truth = ("galactic", (0.0, 0.0, 0.0))
        self.set_velocity_cartesian(*velocity_vector_cartesian)
        self._sync()

    # --- Frames ---------------------------------------------------------------------

    def _anchor(self, frame):
        if frame == "galactic":
            return (0.0, 0.0, 0.0)
        if frame == "sector":
            return self._sector_center
        return self._star_center

    @property
    def sector_center(self):
        """tuple: The sector's center, galactic, in this position's unit."""
        return self._sector_center

    @property
    def has_system_frame(self):
        """Whether the system frame exists (not a star, and a star is known)."""
        return not self._is_star and self._star_center is not None

    @property
    def sector_edge_pc(self):
        """float or None: The sector grid's edge length, parsecs."""
        return self._sector_edge_pc

    @property
    def sector_address(self):
        """tuple or None: The `(ring, layer, slot)` of the sector cell the
        position is in; `None` when the sector edge is not known."""
        return self._sector_address

    def set_sector_edge_pc(self, sector_edge_pc):
        """Sets (or, with `None`, forgets) the sector grid's edge length,
        parsecs; the address is worked out again."""
        if sector_edge_pc is not None:
            _finite((sector_edge_pc,), "sector edge")
            if sector_edge_pc <= 0.0:
                raise ValueError(f"the sector edge must be positive, got {sector_edge_pc!r}")
            sector_edge_pc = float(sector_edge_pc)
        self._sector_edge_pc = sector_edge_pc
        self._sector_address = None
        self._sync()

    def set_sector_address(self, ring_index, layer_index, slot_index):
        """Moves the body to the sector cell `(ring, layer, slot)`, keeping
        its offset from the sector's center (`carry_sector_center` to that
        cell's center). Raises `ValueError` without a sector edge or for a
        slot the ring does not have."""
        if self._sector_edge_pc is None:
            raise ValueError("the sector edge is not known, so there is no sector address to set")
        center_pc = geometry.sector_position_pc(ring_index, layer_index, slot_index, self._sector_edge_pc)
        from_pc = physical_constants.PARSEC_M / self.length_unit_m
        self.carry_sector_center(tuple(c * from_pc for c in center_pc))

    def _galactic_cartesian(self):
        frame, xyz = self._truth
        anchor = self._anchor(frame)
        return tuple(a + b for a, b in zip(anchor, xyz))

    def _sync(self):
        gal = self._galactic_cartesian()
        self._coords = {}
        for frame in FRAMES:
            if frame == "system" and not self.has_system_frame:
                continue
            if frame == self._truth[0]:
                xyz = self._truth[1]
            else:
                anchor = self._anchor(frame)
                xyz = tuple(g - a for g, a in zip(gal, anchor))
            self._coords[frame] = {
                "cartesian": xyz,
                "cylindrical": cartesian_to_cylindrical(*xyz),
                "spherical": cartesian_to_spherical(*xyz),
            }
        if self._sector_edge_pc is not None:
            to_pc = self.length_unit_m / physical_constants.PARSEC_M
            self._sector_address = geometry.sector_address_at(tuple(c * to_pc for c in gal), self._sector_edge_pc)

    def _set(self, frame, xyz):
        _finite(xyz, "coordinates")
        if frame == "system" and not self.has_system_frame:
            raise ValueError("a star, or a body with no star reference, has no system coordinates")
        self._truth = (frame, tuple(float(v) for v in xyz))
        self._sync()

    # --- Setters: any frame, any form -------------------------------------------------

    def set_galactic_cartesian(self, x, y, z):
        self._set("galactic", (x, y, z))

    def set_galactic_cylindrical(self, r, theta, z):
        self._set("galactic", cylindrical_to_cartesian(r, theta, z))

    def set_galactic_spherical(self, r, theta, phi):
        self._set("galactic", spherical_to_cartesian(r, theta, phi))

    def set_sector_cartesian(self, x, y, z):
        self._set("sector", (x, y, z))

    def set_sector_cylindrical(self, r, theta, z):
        self._set("sector", cylindrical_to_cartesian(r, theta, z))

    def set_sector_spherical(self, r, theta, phi):
        self._set("sector", spherical_to_cartesian(r, theta, phi))

    def set_system_cartesian(self, x, y, z):
        self._set("system", (x, y, z))

    def set_system_cylindrical(self, r, theta, z):
        self._set("system", cylindrical_to_cartesian(r, theta, z))

    def set_system_spherical(self, r, theta, phi):
        self._set("system", spherical_to_cartesian(r, theta, phi))

    def set_coordinates(self, frame, form, values):
        """Sets the position from `values` in `frame` and `form` (the
        `set_<frame>_<form>` methods by name)."""
        if frame not in FRAMES or form not in FORMS:
            raise KeyError(f"unknown frame {frame!r} or form {form!r}")
        getattr(self, f"set_{frame}_{form}")(*values)

    # --- Anchors --------------------------------------------------------------------

    def set_sector_center(self, sector_center_galactic):
        """Moves the sector's center; the body stays where it is."""
        _finite(sector_center_galactic, "sector center")
        gal = self._galactic_cartesian()
        self._sector_center = tuple(float(v) for v in sector_center_galactic)
        self._rebase(gal)

    def set_nearest_star_center(self, star_center_galactic):
        """Sets, moves or (`None`) drops the nearest star's center; the body
        stays where it is."""
        if star_center_galactic is not None:
            _finite(star_center_galactic, "star center")
        gal = self._galactic_cartesian()
        galactic_velocity = self._velocity_galactic()
        self._star_center = None if star_center_galactic is None else tuple(float(v) for v in star_center_galactic)
        if not self.has_system_frame:
            self._velocity_truth = ("galactic", galactic_velocity)
        self._rebase(gal)

    def _rebase(self, gal):
        if self._truth[0] == "system" and not self.has_system_frame:
            self._truth = ("galactic", gal)
        elif self._truth[0] != "galactic":
            self._truth = (self._truth[0], tuple(g - a for g, a in zip(gal, self._anchor(self._truth[0]))))
        else:
            self._truth = ("galactic", gal)
        self._sync()

    # --- Velocity, epoch and the next-due time ----------------------------------------

    def _velocity_frame(self, frame):
        key = str(frame).lower()
        if key == "sector":
            key = "galactic"
        if key not in VELOCITY_FRAMES:
            raise KeyError(f"Invalid velocity frame {frame!r}. Valid frames: {', '.join(VELOCITY_FRAMES)}, sector.")
        if key == "system" and not self.has_system_frame:
            raise ValueError("a star, or a body with no star reference, has no system velocity")
        return key

    def _velocity_galactic(self):
        frame, vector = self._velocity_truth
        if frame == "galactic":
            return vector
        return tuple(a + b for a, b in zip(self._star_velocity, vector))

    def set_velocity_cartesian(self, vx, vy, vz, frame="galactic"):
        """Sets the velocity, m/s, in `frame` ("galactic" or "system"; the
        system velocity is relative to the nearest star). Raises
        `ValueError` for a velocity (or, in the system frame, the galactic
        velocity it makes) at or above the speed of light."""
        key = self._velocity_frame(frame)
        _finite((vx, vy, vz), "velocity")
        vector = (float(vx), float(vy), float(vz))
        galactic = vector if key == "galactic" else tuple(a + b for a, b in zip(self._star_velocity, vector))
        for candidate in (vector, galactic):
            if math.sqrt(sum(v * v for v in candidate)) >= SPEED_OF_LIGHT_MS:
                raise ValueError("a body cannot move at or above the speed of light")
        self._velocity_truth = (key, vector)

    def get_velocity_vector(self, frame="galactic"):
        """The velocity, m/s, in `frame` ("galactic", "sector" or "system");
        the system frame is `None` for a star or a body with no star
        reference."""
        key = str(frame).lower()
        if key == "system" and not self.has_system_frame:
            return None
        key = self._velocity_frame(frame)
        if key == "galactic":
            return self._velocity_galactic()
        return self._velocity_truth[1] if self._velocity_truth[0] == "system" else tuple(
            g - s for g, s in zip(self._velocity_galactic(), self._star_velocity))

    def get_speed(self, frame="galactic"):
        velocity = self.get_velocity_vector(frame)
        return None if velocity is None else math.sqrt(sum(v * v for v in velocity))

    def get_velocity_direction(self, frame="galactic"):
        """The unit vector the body moves along in `frame`; (0, 0, 0) at rest."""
        velocity = self.get_velocity_vector(frame)
        if velocity is None:
            return None
        speed = math.sqrt(sum(v * v for v in velocity))
        if speed == 0.0:
            return (0.0, 0.0, 0.0)
        return tuple(v / speed for v in velocity)

    @property
    def star_velocity(self):
        """tuple: The nearest star's velocity, galactic axes, m/s."""
        return self._star_velocity

    def set_star_velocity(self, star_velocity_galactic):
        """Sets the nearest star's velocity; the body keeps its galactic
        velocity (its system velocity changes)."""
        _finite(star_velocity_galactic, "star velocity")
        galactic = self._velocity_galactic()
        self._star_velocity = tuple(float(v) for v in star_velocity_galactic)
        self._velocity_truth = ("galactic", galactic)

    def carry_star_velocity(self, star_velocity_galactic):
        """Sets the nearest star's velocity and the body's with it: its
        system velocity stays, its galactic velocity follows. Needs a system
        frame."""
        _finite(star_velocity_galactic, "star velocity")
        if not self.has_system_frame:
            raise ValueError("a body with no system frame has no relative velocity to carry")
        keep = self.get_velocity_vector("system")
        self._star_velocity = tuple(float(v) for v in star_velocity_galactic)
        self._velocity_truth = ("system", keep)

    @property
    def epoch_unix(self):
        """float or None: When the position and velocity hold, Unix seconds."""
        return self._epoch_unix

    def set_epoch_unix(self, epoch_unix):
        """Sets (or, with `None`, clears) the time the position and velocity hold at."""
        if epoch_unix is not None:
            _finite((epoch_unix,), "epoch")
            epoch_unix = float(epoch_unix)
        self._epoch_unix = epoch_unix

    def get_time_to_observable_movement(self, scale="system"):
        """Seconds until the body has moved the `THRESHOLDS_M[scale]`
        distance ("galactic", "system" or "planetary"), capped at
        `MAX_UPDATE_INTERVAL_S`; the cap for a body at rest. A "galactic"
        move is measured by the galactic velocity, the others by the
        velocity relative to the star (the galactic one for a body that has
        no system frame)."""
        if scale not in THRESHOLDS_M:
            raise KeyError(f"unknown scale {scale!r}; use one of {sorted(THRESHOLDS_M)}")
        frame = "galactic" if scale == "galactic" or not self.has_system_frame else "system"
        speed = self.get_speed(frame)
        if speed == 0.0:
            return MAX_UPDATE_INTERVAL_S
        return min(THRESHOLDS_M[scale] / speed, MAX_UPDATE_INTERVAL_S)

    # --- Carrying the anchors ---------------------------------------------------------

    def carry_sector_center(self, sector_center_galactic):
        """Moves the sector's center and the body with it: its sector
        coordinates stay, its galactic ones follow. For placing a sector
        whose bodies were made before its place in the galaxy was known."""
        _finite(sector_center_galactic, "sector center")
        keep = self._coords["sector"]["cartesian"]
        self._sector_center = tuple(float(v) for v in sector_center_galactic)
        self._truth = ("sector", keep)
        self._sync()

    def carry_anchors(self, sector_center_galactic, star_center_galactic):
        """Moves the sector's center and the nearest star's center together,
        and the body with them: its system offset stays, its sector and
        galactic coordinates follow. Needs a system frame."""
        _finite(sector_center_galactic, "sector center")
        _finite(star_center_galactic, "star center")
        if not self.has_system_frame:
            raise ValueError("a body with no system frame has no offset to carry")
        keep = self._coords["system"]["cartesian"]
        self._sector_center = tuple(float(v) for v in sector_center_galactic)
        self._star_center = tuple(float(v) for v in star_center_galactic)
        self._truth = ("system", keep)
        self._sync()

    def carry_star_center(self, star_center_galactic):
        """Moves the nearest star's center and the body with it: its system
        coordinates stay, its galactic and sector ones follow. Needs a
        system frame."""
        _finite(star_center_galactic, "star center")
        if self._is_star:
            raise ValueError("a star has no system frame to carry")
        keep = self._coords["system"]["cartesian"] if self.has_system_frame else None
        self._star_center = tuple(float(v) for v in star_center_galactic)
        self._truth = ("system", keep if keep is not None else (0.0, 0.0, 0.0))
        self._sync()

    # --- Mass and mu ------------------------------------------------------------------

    def set_mass(self, mass_kg):
        """Sets the mass (kg); `mu` follows. Raises `ValueError` if it is
        negative or not finite."""
        _finite((mass_kg,), "mass")
        if mass_kg < 0.0:
            raise ValueError(f"mass must not be negative, got {mass_kg!r}")
        self._mass_kg = float(mass_kg)
        self._mu = mu_of(self._mass_kg)

    def set_mu(self, mu):
        """Sets the gravitational parameter (m^3/s^2); the mass follows."""
        _finite((mu,), "mu")
        if mu < 0.0:
            raise ValueError(f"mu must not be negative, got {mu!r}")
        self._mu = float(mu)
        self._mass_kg = mass_of(self._mu)

    @property
    def mass_kg(self):
        """The mass, kg (`None` if never set)."""
        return self._mass_kg

    @property
    def mu(self):
        """G * mass, m^3/s^2 (`None` if the mass was never set)."""
        return self._mu

    # --- Reading ----------------------------------------------------------------------

    def get_coordinates(self, frame, form):
        """The position in `frame` ("galactic", "sector", "system") and
        `form` ("cartesian", "cylindrical", "spherical"); `None` for the
        system frame of a star or a body with no star reference."""
        frame_key, form_key = str(frame).lower(), str(form).lower()
        if frame_key not in FRAMES:
            raise KeyError(f"Invalid frame {frame!r}. Valid frames: {', '.join(FRAMES)}.")
        if form_key not in FORMS:
            raise KeyError(f"Invalid coordinate type {form!r}. Valid types: {', '.join(FORMS)}.")
        entry = self._coords.get(frame_key)
        return None if entry is None else entry[form_key]

    def is_star(self):
        return self._is_star


def axis_property(index):
    """
    A property for one Cartesian axis (0 x, 1 y, 2 z), in AU, of a body that
    holds a `SpatialPosition3D` as `spatial` and mixes in `HoldsOrbitPosition`.
    Reading gives the body's offset from its primary (the "system" frame);
    writing before the body has a position stages the value until all three
    axes are given, so `fields_from_dict` and `None` defaults keep working.
    """
    def read(self):
        if self.spatial is not None:
            return self.spatial.get_coordinates("system", "cartesian")[index]
        staged = self._staged_au
        return None if staged is None else staged[index]

    def write(self, value):
        staged = list(self._staged_au or (None, None, None))
        staged[index] = value
        self._staged_au = tuple(staged)
        if all(v is not None for v in staged):
            self.set_position_au(*staged)

    return property(read, write)


def velocity_axis_property(index):
    """
    `axis_property`'s counterpart for the velocity: a property for one
    Cartesian axis (0 x, 1 y, 2 z), in km/s, of a body that holds a
    `SpatialPosition3D` as `spatial` and mixes in `HoldsOrbitPosition`.
    Reading gives the body's velocity relative to its primary (the "system"
    frame); writing before the body has a position stages the value until
    all three axes are given.
    """
    def read(self):
        if self.spatial is not None:
            return self.spatial.get_velocity_vector("system")[index] / 1000.0
        staged = self._staged_velocity_kms
        return None if staged is None else staged[index]

    def write(self, value):
        staged = list(self._staged_velocity_kms or (None, None, None))
        staged[index] = value
        self._staged_velocity_kms = tuple(staged)
        if all(v is not None for v in staged):
            self.set_velocity_kms(*staged)

    return property(read, write)


class HoldsOrbitPosition:
    """
    Mixin for a planet, moon or comet: its position is a `SpatialPosition3D`
    in AU (`spatial`), whose "system" frame is the offset from the body's
    primary (the star, or the parent planet for a moon). The subclass
    declares its three attributes with `axis_property`; `set_position_au`
    is the one way to move it. Its velocity relative to the primary is the
    "system" frame velocity of `spatial`, in km/s (`velocity_axis_property`);
    `set_velocity_kms` sets it.
    """

    spatial = None
    _staged_au = None
    _staged_velocity_kms = None

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
