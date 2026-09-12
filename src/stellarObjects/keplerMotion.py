# stellarObjects/keplerMotion.py

"""
Kepler Orbital Motion
=====================

Two-body Kepler/Barker-equation propagation for eccentric (elliptical,
`0 <= e < 1`) and parabolic (`e == 1`) orbits -- the physically correct
alternative to `utils.orbital_position_au`'s uniform-angular-speed
circular-orbit model, needed for any body whose orbit isn't nearly
circular. `cometData.Comet` is the motivating case here (see
docs/design/comet-orbital-realism.md): a Halley-like bound comet reaches
eccentricities around 0.97, and a single-apparition long-period comet is,
for practical purposes, parabolic. `planetPhysics.py`'s planets/moons
stay on the existing circular-orbit model (`orbital_phase_deg` advancing
linearly) -- real planetary eccentricities are small enough that the
difference is negligible, and introducing time-varying angular speed
there would be a much larger, unrelated change.

Kepler's second law (equal areas in equal time) means a body on an
eccentric orbit moves far faster near perihelion than near aphelion --
angular position is NOT a linear function of time the way it is for a
circular orbit, so a "mean anomaly" (a fictitious angle that DOES advance
linearly with time) has to be converted to the real "true anomaly" (the
body's actual angular position) via Kepler's equation (elliptical) or
Barker's equation (parabolic) first.

All functions here use this generator's existing AU/years/solar-mass unit
convention (see `planetPhysics.calculate_orbital_period_years`), under
which the two-body gravitational parameter mu = G*M simplifies to
`4*pi^2 * M_solar_masses` -- exactly `4*pi^2` for a one-solar-mass primary,
since Kepler's third law (`T^2 = a^3/M` in these units) reduces to
`T = 2*pi*sqrt(a^3/mu)` with `mu = 4*pi^2*M`.

Once a true anomaly and radius are known, 3D placement reuses
`utils.orbital_position_au` as-is: that function's rotation math takes an
arbitrary radius and "argument of latitude" angle -- it happens to be
called elsewhere with a *fixed* radius and phase-as-argument-of-latitude
(a circular orbit has no periapsis, so phase alone plays that role there)
-- here it's called with a *varying* radius (from Kepler's/Barker's
equation) and `argument_of_periapsis + true_anomaly` standing in for that
angle instead. No new rotation math is needed.
"""

import math

from . import physical_constants
from .utils import orbital_position_au

TWO_PI = 2 * math.pi

AU3_PER_YR2_PER_SOLAR_MASS = 4 * math.pi ** 2
"""
float: The two-body gravitational parameter `mu = G*M`, in AU^3/yr^2, for
a one-solar-mass primary -- see this module's docstring for the Kepler's-
third-law derivation. Multiply by a primary's mass in solar masses to get
its own `mu` in these units.
"""


def gravitational_parameter_au3_yr2(primary_mass_solar):
    """
    `mu = G*M`, in AU^3/yr^2, for a primary of the given mass.

    Args:
        primary_mass_solar (float): The primary's mass, in solar masses.

    Returns:
        float: `mu`, in AU^3/yr^2.
    """
    return AU3_PER_YR2_PER_SOLAR_MASS * primary_mass_solar


def mean_motion_per_year(semi_major_axis_au, primary_mass_solar):
    """
    An elliptical orbit's mean motion `n = 2*pi/P` (radians/year), via
    Kepler's third law (`planetPhysics.calculate_orbital_period_years`'s
    same `T = sqrt(a^3/M)` formula) -- the constant rate a *fictitious*
    mean anomaly advances at, standing in for the body's real (non-
    uniform) angular speed.

    Args:
        semi_major_axis_au (float): Orbital semi-major axis, in AU.
        primary_mass_solar (float): The primary's mass, in solar masses.

    Returns:
        float: Mean motion, in radians/year.
    """
    period_years = math.sqrt(semi_major_axis_au ** 3 / primary_mass_solar)
    return TWO_PI / period_years


def solve_eccentric_anomaly(mean_anomaly_rad, eccentricity, tolerance=1e-10, max_iterations=100):
    """
    Newton-Raphson solution of Kepler's equation `M = E - e*sin(E)` for
    the eccentric anomaly `E`, given the mean anomaly `M` (radians) and
    eccentricity `e` in `[0, 1)`.

    Kepler's equation has no closed-form solution -- this is the standard
    numerical approach (e.g. Meeus, "Astronomical Algorithms", ch. 30),
    converging quadratically from an initial guess of `M` itself for
    low-to-moderate eccentricity, or `pi` for a high-eccentricity orbit
    (where `M` alone is a poor starting guess -- see e.g. Danby, "Fundamentals
    of Celestial Mechanics"). `tolerance`/`max_iterations` bound how far
    this iterates; for any physically realistic comet eccentricity (this
    generator caps elliptical comets below `PARABOLIC_COMET_ECCENTRICITY_RANGE`'s
    floor -- see `program_constants.COMET_PERIOD_CLASSES`), convergence to
    machine precision typically takes under 10 iterations.

    Args:
        mean_anomaly_rad (float): Mean anomaly, in radians (any real
            value -- wrapped into `[0, 2*pi)` internally).
        eccentricity (float): Orbital eccentricity, `0 <= e < 1`.
        tolerance (float, optional): Convergence tolerance, in radians.
        max_iterations (int, optional): Maximum Newton-Raphson iterations.

    Returns:
        float: The eccentric anomaly `E`, in radians.
    """
    m = mean_anomaly_rad % TWO_PI
    eccentric_anomaly = m if eccentricity < 0.8 else math.pi
    for _ in range(max_iterations):
        delta = (eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly) - m) / (
            1 - eccentricity * math.cos(eccentric_anomaly)
        )
        eccentric_anomaly -= delta
        if abs(delta) < tolerance:
            break
    return eccentric_anomaly


def true_anomaly_and_distance_elliptical(mean_anomaly_rad, eccentricity, semi_major_axis_au):
    """
    The true anomaly and current orbital radius for a bound (`e < 1`)
    orbit at a given mean anomaly.

    `true_anomaly` uses the `atan2`-based half-angle formula (rather than
    the more commonly quoted `tan(nu/2) = sqrt((1+e)/(1-e))*tan(E/2)`),
    which is well-behaved at every eccentric anomaly (no division blowing
    up as `E` approaches an odd multiple of pi, where a plain tangent
    would).

    Args:
        mean_anomaly_rad (float): Mean anomaly, in radians.
        eccentricity (float): Orbital eccentricity, `0 <= e < 1`.
        semi_major_axis_au (float): Orbital semi-major axis, in AU.

    Returns:
        tuple: `(true_anomaly_rad, distance_au)`.
    """
    eccentric_anomaly = solve_eccentric_anomaly(mean_anomaly_rad, eccentricity)
    true_anomaly_rad = 2 * math.atan2(
        math.sqrt(1 + eccentricity) * math.sin(eccentric_anomaly / 2),
        math.sqrt(1 - eccentricity) * math.cos(eccentric_anomaly / 2),
    )
    distance_au = semi_major_axis_au * (1 - eccentricity * math.cos(eccentric_anomaly))
    return true_anomaly_rad, distance_au


def _real_cube_root(value):
    """
    The real cube root of `value`, including negative `value` -- Python's
    `**(1/3)` returns a complex result for a negative base with a
    fractional exponent, which `solve_barker_equation`'s closed-form
    solution needs to avoid (its second Cardano term is always negative).

    Args:
        value (float): The value to take the cube root of.

    Returns:
        float: The real cube root of `value`.
    """
    return math.copysign(abs(value) ** (1 / 3), value)


def solve_barker_equation(parabolic_mean_anomaly):
    """
    Closed-form (Cardano) solution of Barker's equation, `D^3 + 3*D =
    3*Mp`, for `D = tan(true_anomaly / 2)` -- the parabolic (`e == 1`)
    analog of Kepler's equation. Unlike the elliptical case, this depressed
    cubic has exactly one real root (its discriminant is always positive,
    since `D^3 + 3D` is strictly increasing), so an exact closed form
    exists and no iteration is needed (see e.g. Danby, "Fundamentals of
    Celestial Mechanics", ch. 6, or Meeus ch. 34).

    Derivation: for the depressed cubic `t^3 + p*t + q = 0` with `p = 3`,
    `q = -3*Mp`, Cardano's formula gives
    `t = cbrt(-q/2 + sqrt((q/2)^2 + (p/3)^3)) + cbrt(-q/2 - sqrt((q/2)^2 + (p/3)^3))`,
    which simplifies (with `(p/3)^3 = 1`) to the `s`/`w` form below.

    Args:
        parabolic_mean_anomaly (float): `Mp = sqrt(mu / (2*q^3)) * (t -
            t_perihelion)` (see `parabolic_mean_anomaly`) -- negative
            before perihelion, zero at perihelion, positive after.

    Returns:
        float: `D = tan(true_anomaly / 2)`.
    """
    w = 1.5 * parabolic_mean_anomaly
    s = math.sqrt(w * w + 1)
    return _real_cube_root(w + s) + _real_cube_root(w - s)


def parabolic_mean_anomaly(time_since_perihelion_years, perihelion_distance_au, primary_mass_solar):
    """
    The parabolic "mean anomaly" `Mp = sqrt(mu / (2*q^3)) * dt` --
    Barker's equation's own fictitious linear-in-time angle, playing the
    same role for a parabolic orbit that the ordinary mean anomaly `M =
    n*dt` plays for an elliptical one (see `mean_motion_per_year`), just
    without a period to normalize by (a parabolic orbit's period is
    infinite).

    Args:
        time_since_perihelion_years (float): Time since (or, if negative,
            until) perihelion passage, in years.
        perihelion_distance_au (float): Perihelion distance `q`, in AU.
        primary_mass_solar (float): The primary's mass, in solar masses.

    Returns:
        float: The parabolic mean anomaly `Mp` (dimensionless).
    """
    mu = gravitational_parameter_au3_yr2(primary_mass_solar)
    return math.sqrt(mu / (2 * perihelion_distance_au ** 3)) * time_since_perihelion_years


def true_anomaly_and_distance_parabolic(parabolic_mean_anomaly_value, perihelion_distance_au):
    """
    The true anomaly and current orbital radius for a parabolic (`e ==
    1`) orbit at a given parabolic mean anomaly.

    `distance_au = q*(1 + D^2)` follows from the conic equation `r =
    p/(1 + cos(true_anomaly))` with semi-latus-rectum `p = 2*q` (true for
    any parabola) and the half-angle identity `1 + cos(true_anomaly) =
    2*cos^2(true_anomaly/2)`.

    Args:
        parabolic_mean_anomaly_value (float): `Mp`, see
            `parabolic_mean_anomaly`.
        perihelion_distance_au (float): Perihelion distance `q`, in AU.

    Returns:
        tuple: `(true_anomaly_rad, distance_au)`.
    """
    d = solve_barker_equation(parabolic_mean_anomaly_value)
    true_anomaly_rad = 2 * math.atan(d)
    distance_au = perihelion_distance_au * (1 + d * d)
    return true_anomaly_rad, distance_au


def vis_viva_speed_kms(distance_au, semi_major_axis_au, primary_mass_solar):
    """
    The vis-viva equation, `v = sqrt(mu * (2/r - 1/a))` -- a body's own
    orbital speed at a given distance `r` from its primary, general for
    any conic section. `semi_major_axis_au = math.inf` collapses the `1/a`
    term to `0`, giving the parabolic-orbit special case `v =
    sqrt(2*mu/r)` for free, with no separate formula needed.

    Args:
        distance_au (float): Current distance from the primary, in AU.
        semi_major_axis_au (float): Orbital semi-major axis, in AU, or
            `math.inf` for a parabolic orbit.
        primary_mass_solar (float): The primary's mass, in solar masses.

    Returns:
        float: Orbital speed, in km/s.
    """
    mu = gravitational_parameter_au3_yr2(primary_mass_solar)
    inv_a = 0.0 if math.isinf(semi_major_axis_au) else 1 / semi_major_axis_au
    speed_au_per_year = math.sqrt(mu * (2 / distance_au - inv_a))
    au_per_year_to_kms = physical_constants.AU_TO_KM / physical_constants.SECONDS_PER_YEAR
    return speed_au_per_year * au_per_year_to_kms


def comet_orbital_state(orbit_type, perihelion_distance_au, eccentricity, inclination_deg,
                        arg_periapsis_deg, ascending_node_deg, primary_mass_solar,
                        mean_anomaly_rad=None, parabolic_mean_anomaly_value=None,
                        orbital_period_years=None):
    """
    A star-bound comet's full orbital state -- 3D position, distance from
    its primary, and orbital speed -- at whatever point along its orbit
    `mean_anomaly_rad`/`parabolic_mean_anomaly_value` (whichever applies
    to `orbit_type`) describes.

    Dispatches on `orbit_type` ("elliptical" uses Kepler's equation and a
    finite semi-major axis derived from `perihelion_distance_au`/
    `eccentricity`; "parabolic" uses Barker's equation and an infinite
    semi-major axis) and reuses `utils.orbital_position_au` for the actual
    3D placement -- see this module's own docstring for why that function,
    built for a fixed-radius circular orbit, applies unchanged here.

    Args:
        orbit_type (str): `"elliptical"` or `"parabolic"`.
        perihelion_distance_au (float): Perihelion distance `q`, in AU.
        eccentricity (float): Orbital eccentricity (`0 <= e < 1` for
            elliptical; ignored for parabolic, which always uses the
            exact `e == 1` Barker solution regardless of the specific
            near-1 value a `Comet` may still store for descriptive text --
            see `program_constants.PARABOLIC_COMET_ECCENTRICITY_RANGE`).
        inclination_deg (float): Orbital plane tilt, in degrees.
        arg_periapsis_deg (float): Argument of periapsis, in degrees.
        ascending_node_deg (float): Longitude of the ascending node, in
            degrees.
        primary_mass_solar (float): The host star's mass, in solar
            masses.
        mean_anomaly_rad (float, optional): Required for `orbit_type ==
            "elliptical"` -- the current mean anomaly, in radians.
        parabolic_mean_anomaly_value (float, optional): Required for
            `orbit_type == "parabolic"` -- see `parabolic_mean_anomaly`.
        orbital_period_years (float, optional): The orbit's own period, in
            years -- only meaningful (and only used) for `orbit_type ==
            "elliptical"`.

    Returns:
        dict: `position_x_au`, `position_y_au`, `position_z_au`,
             `distance_au`, `orbital_speed_kms`.

    Raises:
        ValueError: If `orbit_type` isn't `"elliptical"` or `"parabolic"`,
            or the anomaly required for that type wasn't given.
    """
    if orbit_type == "elliptical":
        if mean_anomaly_rad is None:
            raise ValueError("comet_orbital_state: mean_anomaly_rad is required for an elliptical orbit.")
        semi_major_axis_au = perihelion_distance_au / (1 - eccentricity)
        true_anomaly_rad, distance_au = true_anomaly_and_distance_elliptical(
            mean_anomaly_rad, eccentricity, semi_major_axis_au
        )
    elif orbit_type == "parabolic":
        if parabolic_mean_anomaly_value is None:
            raise ValueError("comet_orbital_state: parabolic_mean_anomaly_value is required for a parabolic orbit.")
        semi_major_axis_au = math.inf
        true_anomaly_rad, distance_au = true_anomaly_and_distance_parabolic(
            parabolic_mean_anomaly_value, perihelion_distance_au
        )
    else:
        raise ValueError(f"comet_orbital_state: unknown orbit_type '{orbit_type}'.")

    argument_of_latitude_deg = math.degrees(math.radians(arg_periapsis_deg) + true_anomaly_rad) % 360
    position_x_au, position_y_au, position_z_au = orbital_position_au(
        distance_au, inclination_deg, ascending_node_deg, argument_of_latitude_deg
    )
    orbital_speed_kms = vis_viva_speed_kms(distance_au, semi_major_axis_au, primary_mass_solar)

    return {
        "position_x_au": position_x_au,
        "position_y_au": position_y_au,
        "position_z_au": position_z_au,
        "distance_au": distance_au,
        "orbital_speed_kms": orbital_speed_kms,
    }
