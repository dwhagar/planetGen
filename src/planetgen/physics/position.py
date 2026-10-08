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

All lengths are metres, speeds metres a second and masses kilograms; the
storage layer converts (km, mpc, ly) at its edge (`planetgen.physics.units`).

Precision: the position that was last set is kept as given, in the frame it
was given in, and the other frames are derived from it. A moon placed by its
system offset therefore keeps its metres exactly even though the same point
in galactic coordinates (about 1e20 m from the center) cannot.

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

from planetgen.physics import constants as physical_constants

FRAMES = ("galactic", "sector", "system")
"""tuple: The frames a position is kept in."""

FORMS = ("cartesian", "cylindrical", "spherical")
"""tuple: The coordinate forms each frame is kept in."""

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
        star_center_galactic (tuple or None): The nearest star's center,
            galactic metres; `None` for a body with no star reference.
        is_star (bool): A star has no system frame.
        mass_kg (float or None): The body's mass; its `mu` follows.
    """

    def __init__(self, galactic_cartesian, sector_center_galactic, velocity_vector_cartesian=(0.0, 0.0, 0.0),
                 star_center_galactic=None, is_star=False, mass_kg=None):
        _finite(galactic_cartesian, "galactic position")
        _finite(sector_center_galactic, "sector center")
        if star_center_galactic is not None:
            _finite(star_center_galactic, "star center")
        self._sector_center = tuple(float(v) for v in sector_center_galactic)
        self._star_center = None if star_center_galactic is None else tuple(float(v) for v in star_center_galactic)
        self._is_star = bool(is_star)
        self._truth = ("galactic", tuple(float(v) for v in galactic_cartesian))
        self._mass_kg = None
        self._mu = None
        if mass_kg is not None:
            self.set_mass(mass_kg)
        self._velocity = (0.0, 0.0, 0.0)
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
    def has_system_frame(self):
        """Whether the system frame exists (not a star, and a star is known)."""
        return not self._is_star and self._star_center is not None

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
        self._star_center = None if star_center_galactic is None else tuple(float(v) for v in star_center_galactic)
        self._rebase(gal)

    def _rebase(self, gal):
        if self._truth[0] == "system" and not self.has_system_frame:
            self._truth = ("galactic", gal)
        elif self._truth[0] != "galactic":
            self._truth = (self._truth[0], tuple(g - a for g, a in zip(gal, self._anchor(self._truth[0]))))
        else:
            self._truth = ("galactic", gal)
        self._sync()

    # --- Velocity and the next-due time -----------------------------------------------

    def set_velocity_cartesian(self, vx, vy, vz):
        """Sets the velocity (galactic axes, m/s). Raises `ValueError` for
        a speed at or above the speed of light."""
        _finite((vx, vy, vz), "velocity")
        if math.sqrt(vx * vx + vy * vy + vz * vz) >= SPEED_OF_LIGHT_MS:
            raise ValueError("a body cannot move at or above the speed of light")
        self._velocity = (float(vx), float(vy), float(vz))

    def get_velocity_vector(self):
        return self._velocity

    def get_speed(self):
        return math.sqrt(sum(v * v for v in self._velocity))

    def get_velocity_direction(self):
        """The unit vector the body moves along; (0, 0, 0) at rest."""
        speed = self.get_speed()
        if speed == 0.0:
            return (0.0, 0.0, 0.0)
        return tuple(v / speed for v in self._velocity)

    def get_time_to_observable_movement(self, scale="system"):
        """Seconds until the body has moved the `THRESHOLDS_M[scale]`
        distance ("galactic", "system" or "planetary"), capped at
        `MAX_UPDATE_INTERVAL_S`; the cap for a body at rest."""
        if scale not in THRESHOLDS_M:
            raise KeyError(f"unknown scale {scale!r}; use one of {sorted(THRESHOLDS_M)}")
        speed = self.get_speed()
        if speed == 0.0:
            return MAX_UPDATE_INTERVAL_S
        return min(THRESHOLDS_M[scale] / speed, MAX_UPDATE_INTERVAL_S)

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
