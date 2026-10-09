# planetgen/physics/kepler.py

"""
Kepler Orbital Motion
=====================

Two-body Kepler/Barker-equation propagation for eccentric (elliptical,
`0 <= e < 1`) and parabolic (`e == 1`) orbits -- the physically correct
alternative to `orbits.orbital_position_au`'s uniform-angular-speed
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
`orbits.orbital_position_au` as-is: that function's rotation math takes an
arbitrary radius and "argument of latitude" angle -- it happens to be
called elsewhere with a *fixed* radius and phase-as-argument-of-latitude
(a circular orbit has no periapsis, so phase alone plays that role there)
-- here it's called with a *varying* radius (from Kepler's/Barker's
equation) and `argument_of_periapsis + true_anomaly` standing in for that
angle instead. No new rotation math is needed.
"""

import math
import warnings

import numpy as np
from scipy import optimize

from planetgen.physics import constants
from planetgen.physics.orbits import orbital_position_au
from planetgen.physics.state_vectors import state_from_elements
from planetgen.util.checks import finite_domain

TWO_PI = 2 * math.pi

KEPLER_XTOL = 1e-14
"""float: How close, in radians, Kepler's equation is solved to (about 50 times machine epsilon at 2*pi)."""

AU3_PER_YR2_PER_SOLAR_MASS = 4 * math.pi ** 2
"""
float: The two-body gravitational parameter `mu = G*M`, in AU^3/yr^2, for
a one-solar-mass primary -- see this module's docstring for the Kepler's-
third-law derivation. Multiply by a primary's mass in solar masses to get
its own `mu` in these units.
"""


@finite_domain()
def gravitational_parameter_au3_yr2(primary_mass_solar):
    """
    `mu = G*M`, in AU^3/yr^2, for a primary of the given mass.

    Args:
        primary_mass_solar (float): The primary's mass, in solar masses.

    Returns:
        float: `mu`, in AU^3/yr^2.
    """
    return AU3_PER_YR2_PER_SOLAR_MASS * primary_mass_solar


@finite_domain()
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


@finite_domain()
def solve_eccentric_anomaly(mean_anomaly_rad, eccentricity):
    """
    Solves Kepler's equation `M = E - e*sin(E)` for the eccentric anomaly
    `E`, given the mean anomaly `M` (radians) and eccentricity `e` in
    `[0, 1)`, with `scipy.optimize.brentq`.

    Kepler's equation has no closed-form solution. For `e < 1` the left
    side minus the right, `E - e*sin(E) - M`, rises monotonically from
    `-M <= 0` at `E = 0` to `2*pi - M > 0` at `E = 2*pi`, so that interval
    always brackets the one root and Brent's method always converges,
    including as `e` approaches 1 where a Newton iteration from `M` alone
    wanders.

    Args:
        mean_anomaly_rad (float): Mean anomaly, in radians (any real
            value -- wrapped into `[0, 2*pi)` internally).
        eccentricity (float): Orbital eccentricity, `0 <= e < 1`.

    Returns:
        float: The eccentric anomaly `E`, in radians, in `[0, 2*pi]`.

    Raises:
        ValueError: If `eccentricity` is outside `[0, 1)`.
    """
    if not (0 <= eccentricity < 1):
        raise ValueError(
            f"solve_eccentric_anomaly: eccentricity ({eccentricity}) must be in [0, 1) for Kepler's equation."
        )
    m = mean_anomaly_rad % TWO_PI
    return float(optimize.brentq(_kepler_residual, 0.0, TWO_PI, args=(eccentricity, m), xtol=KEPLER_XTOL))


def _kepler_residual(eccentric_anomaly, eccentricity, mean_anomaly):
    return eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly


def _kepler_slope(eccentric_anomaly, eccentricity, mean_anomaly):
    return 1 - eccentricity * np.cos(eccentric_anomaly)


def _kepler_curvature(eccentric_anomaly, eccentricity, mean_anomaly):
    return eccentricity * np.sin(eccentric_anomaly)


def solve_eccentric_anomalies(mean_anomalies_rad, eccentricities):
    """
    `solve_eccentric_anomaly` for many orbits at once: one vectorised
    Halley iteration (`scipy.optimize.newton` with arrays) for the whole
    batch, and `brentq` for any orbit it did not settle (very eccentric
    ones, mostly).

    Args:
        mean_anomalies_rad (array-like): Mean anomalies, radians (any real values, finite).
        eccentricities (array-like): Eccentricities, each in `[0, 1)`; same length.

    Returns:
        numpy.ndarray: The eccentric anomalies, radians, each in `[0, 2*pi]`.

    Raises:
        ValueError: If any eccentricity is outside `[0, 1)`, any mean
            anomaly is not finite, or the lengths differ.
    """
    m = np.asarray(mean_anomalies_rad, dtype=float)
    e = np.asarray(eccentricities, dtype=float)
    if m.shape != e.shape or m.ndim != 1:
        raise ValueError("solve_eccentric_anomalies: mean anomalies and eccentricities must be 1-D and the same length.")
    if not np.all(np.isfinite(m)):
        raise ValueError("solve_eccentric_anomalies: mean anomalies must be finite numbers.")
    if not np.all((e >= 0) & (e < 1)):
        raise ValueError("solve_eccentric_anomalies: every eccentricity must be in [0, 1) for Kepler's equation.")
    if m.size < 2:  # scipy's newton() takes its scalar path for a single start value
        return np.array([solve_eccentric_anomaly(mean, ecc) for mean, ecc in zip(m, e)], dtype=float)
    m = m % TWO_PI
    start = np.where(e < 0.8, m, np.pi)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # the orbits it leaves unsettled go to brentq below
        solved, converged, _zero_derivative = optimize.newton(
            _kepler_residual, start, fprime=_kepler_slope, args=(e, m), fprime2=_kepler_curvature,
            tol=KEPLER_XTOL, maxiter=50, full_output=True, disp=False)
    solved = np.array(solved, dtype=float)
    settled = np.asarray(converged) & (solved >= 0) & (solved <= TWO_PI)
    for index in np.flatnonzero(~settled):
        solved[index] = solve_eccentric_anomaly(m[index], e[index])
    return solved


@finite_domain()
def true_anomaly_and_distance_elliptical(mean_anomaly_rad, eccentricity, semi_major_axis_au, eccentric_anomaly_rad=None):
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
        eccentric_anomaly_rad (float, optional): The eccentric anomaly for
            this mean anomaly when the caller already has it (from
            `solve_eccentric_anomalies` for a whole batch); solved here
            when left out.

    Returns:
        tuple: `(true_anomaly_rad, distance_au)`.

    Raises:
        ValueError: If `eccentricity` is outside `[0, 1)` -- outside that
            range the `sqrt(1 - eccentricity)` term below goes complex (a
            silent `nan`, since Python's `**`/`math.sqrt` on a negative
            float raises for `math.sqrt` but this uses the two-argument
            `sqrt`-under-`atan2` form) or the orbit isn't actually
            elliptical (`e >= 1` is `comet_orbital_state`'s "parabolic"
            case instead, via `true_anomaly_and_distance_parabolic`) --
            without this check, an out-of-range `e` silently returns a
            plausible-looking but physically meaningless number instead of
            failing loudly.
    """
    if not (0 <= eccentricity < 1):
        raise ValueError(
            f"true_anomaly_and_distance_elliptical: eccentricity ({eccentricity}) must be in [0, 1) "
            f"for a bound elliptical orbit -- e >= 1 is a parabolic/hyperbolic orbit, not this function."
        )
    eccentric_anomaly = (solve_eccentric_anomaly(mean_anomaly_rad, eccentricity)
                         if eccentric_anomaly_rad is None else eccentric_anomaly_rad)
    true_anomaly_rad = 2 * math.atan2(
        math.sqrt(1 + eccentricity) * math.sin(eccentric_anomaly / 2),
        math.sqrt(1 - eccentricity) * math.cos(eccentric_anomaly / 2),
    )
    # a (1 - e cos E) written as q + 2 a e sin^2(E/2) (GEN.108): near
    # perihelion of a near-parabolic orbit, 1 - e cos E subtracts two
    # nearly equal numbers and the radius loses digits; this form doesn't.
    half_sin = math.sin(eccentric_anomaly / 2)
    distance_au = semi_major_axis_au * ((1 - eccentricity) + 2 * eccentricity * half_sin * half_sin)
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


@finite_domain()
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
    # sqrt(w^2 + 1) == |w| to double precision long before w*w overflows
    # (which would make `w - s` an inf - inf = NaN).
    s = math.sqrt(w * w + 1) if abs(w) < 1e150 else abs(w)
    return _real_cube_root(w + s) + _real_cube_root(w - s)


@finite_domain()
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


@finite_domain()
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


@finite_domain(allow_inf=("semi_major_axis_au",))
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
    au_per_year_to_kms = constants.AU_TO_KM / constants.SECONDS_PER_YEAR
    return speed_au_per_year * au_per_year_to_kms


UNIVERSAL_TOLERANCE = 1e-13
"""float: The relative error in the universal anomaly where
`universal_step`'s iteration stops."""

UNIVERSAL_MAX_ITERATIONS = 60
"""int: `universal_step`'s iteration cap; reaching it raises (logged by the
caller), never returns a guess."""


def stumpff_c2_c3(z):
    """
    The Stumpff functions `c2(z) = (1 - cos sqrt z) / z` and
    `c3(z) = (sqrt z - sin sqrt z) / z^(3/2)` (hyperbolic forms for `z < 0`),
    with their series near `z = 0`, where the closed forms cancel to
    nothing: the near-parabolic case.

    Returns:
        tuple: `(c2, c3)`.
    """
    if abs(z) < 1e-3:
        # Six terms of each series: the next is below 1e-20 at |z| < 1e-3.
        c2 = 1 / 2 - z / 24 + z * z / 720 - z ** 3 / 40320 + z ** 4 / 3628800 - z ** 5 / 479001600
        c3 = 1 / 6 - z / 120 + z * z / 5040 - z ** 3 / 362880 + z ** 4 / 39916800 - z ** 5 / 6227020800
        return c2, c3
    if z > 0:
        root = math.sqrt(z)
        return (1 - math.cos(root)) / z, (root - math.sin(root)) / (root ** 3)
    root = math.sqrt(-z)
    return (math.cosh(root) - 1) / -z, (math.sinh(root) - root) / (root ** 3)


def universal_step(position, velocity, mu, dt):
    """
    Moves a two-body state `dt` forward (or back, for negative `dt`) on
    its conic, whatever the conic: the universal-variable form of Kepler's
    equation (docs/design/orbital-updates.md section 10.6, with the
    document's double-counted term corrected), the GEN.108 guard for
    near-parabolic orbits, where the elliptical and hyperbolic forms both
    lose precision, and for steps of any length (a whole number of orbits
    drops out of the bracket below).

    The equation in the universal anomaly `chi`,
    `F = s0 chi^2 c2 + (1 - alpha r0) chi^3 c3 + r0 chi - sqrt(mu) dt`,
    with `alpha = 2 / r0 - v0^2 / mu` and `z = alpha chi^2`, is solved by
    Halley's method kept inside a bracket that always holds the root
    (`F` rises with `chi`, since `F' = r > 0`), falling back to bisection
    whenever a Halley step leaves it: it cannot diverge.

    Args:
        position, velocity (sequence): The state relative to the primary,
            in any consistent units (`L`, `L/T`).
        mu (float): The primary's `G M`, in `L^3/T^2`. Positive.
        dt (float): The time to move, in `T`.

    Returns:
        tuple: `(position, velocity)` after `dt`, each an (x, y, z) tuple.

    Raises:
        ValueError: For a non-positive `mu`, a non-finite input, a body at
            the primary, or no convergence within
            `UNIVERSAL_MAX_ITERATIONS`.
    """
    if not mu > 0 or not math.isfinite(mu):
        raise ValueError(f"universal_step: mu must be positive, got {mu!r}")
    r0v = [float(c) for c in position]
    v0v = [float(c) for c in velocity]
    if not all(math.isfinite(c) for c in r0v + v0v + [dt]):
        raise ValueError("universal_step: the state and step must be finite")
    if dt == 0:
        return tuple(r0v), tuple(v0v)
    r0 = math.sqrt(sum(c * c for c in r0v))
    if r0 == 0:
        raise ValueError("universal_step: the body sits on its primary")
    v0_sq = sum(c * c for c in v0v)
    sqrt_mu = math.sqrt(mu)
    s0 = sum(a * b for a, b in zip(r0v, v0v)) / sqrt_mu
    alpha = 2 / r0 - v0_sq / mu

    def residual(chi):
        z = alpha * chi * chi
        c2, c3 = stumpff_c2_c3(z)
        f = s0 * chi * chi * c2 + (1 - alpha * r0) * chi ** 3 * c3 + r0 * chi - sqrt_mu * dt
        df = s0 * chi * (1 - z * c3) + (1 - alpha * r0) * chi * chi * c2 + r0
        d2f = s0 * (1 - z * c2) + (1 - alpha * r0) * chi * (1 - z * c3)
        return f, df, d2f, c2, c3

    # First guess (Vallado): the mean motion for an ellipse; the signed
    # log form for a hyperbola (the document drops sign(dt)); r0-scaled
    # for a parabola.
    if alpha > 1e-12 / r0:
        chi = sqrt_mu * dt * alpha
        period = 2 * math.pi / (sqrt_mu * alpha ** 1.5)
        if abs(dt) > period:
            # A whole number of orbits changes nothing: step the rest.
            dt = math.fmod(dt, period)
            chi = sqrt_mu * dt * alpha
    elif alpha < -1e-12 / r0:
        a = 1 / alpha
        sign = 1.0 if dt > 0 else -1.0
        inner = (-2 * mu * alpha * dt) / (s0 * sqrt_mu + sign * math.sqrt(-mu * a) * (1 - r0 * alpha))
        chi = sign * math.sqrt(-a) * math.log(inner) if inner > 0 else sign * math.sqrt(-a)
    else:
        chi = sqrt_mu * dt / r0

    # A bracket: F(0) = -sqrt(mu) dt, and F grows at least as fast as
    # r_min * chi, so widen until the sign changes.
    low, high = (0.0, abs(chi) or 1.0) if dt > 0 else (-(abs(chi) or 1.0), 0.0)
    for _ in range(200):
        if dt > 0 and residual(high)[0] < 0:
            low, high = high, high * 2
        elif dt < 0 and residual(low)[0] > 0:
            low, high = low * 2, low
        else:
            break
    if not low <= chi <= high:
        chi = 0.5 * (low + high)

    for _ in range(UNIVERSAL_MAX_ITERATIONS):
        f, df, d2f, c2, c3 = residual(chi)
        if f > 0:
            high = chi
        else:
            low = chi
        step = 2 * f * df / (2 * df * df - f * d2f) if (2 * df * df - f * d2f) != 0 else f / df
        new = chi - step
        if not low < new < high:
            new = 0.5 * (low + high)  # the Halley step left the bracket: bisect
        if abs(new - chi) <= UNIVERSAL_TOLERANCE * max(1.0, abs(new)) or high - low <= 4e-16 * max(1.0, abs(new)):
            chi = new
            break
        chi = new
    else:
        raise ValueError(f"universal_step: no convergence in {UNIVERSAL_MAX_ITERATIONS} iterations")

    z = alpha * chi * chi
    c2, c3 = stumpff_c2_c3(z)
    f = 1 - chi * chi * c2 / r0
    g = dt - chi ** 3 * c3 / sqrt_mu
    r_vec = [f * a + g * b for a, b in zip(r0v, v0v)]
    r = math.sqrt(sum(c * c for c in r_vec))
    g_dot = 1 - chi * chi * c2 / r
    f_dot = sqrt_mu * chi * (z * c3 - 1) / (r * r0)
    v_vec = [f_dot * a + g_dot * b for a, b in zip(r0v, v0v)]
    return tuple(r_vec), tuple(v_vec)


KEPLER_MAX_ECCENTRICITY = 0.999
"""float: The most eccentric ellipse Kepler's equation is trusted for
(GEN.108). Past it, `M = E - e sin E` subtracts two nearly equal numbers
near perihelion and the anomaly loses digits (at e = 1 - 1e-9 it is off
by percent), so `comet_orbital_state` moves the body from perihelion
with `universal_step` instead. Generated elliptical comets stop at 0.999
(`tuning.COMET_PERIOD_CLASSES`); an edit can go further."""


def _near_parabolic_anomaly_and_distance(mean_anomaly_rad, eccentricity, perihelion_distance_au,
                                         primary_mass_solar):
    """The true anomaly and radius of a near-parabolic ellipse at
    `mean_anomaly_rad`, by `universal_step` from perihelion over the time
    the mean anomaly stands for (M / n, M taken in [-pi, pi])."""
    mu = gravitational_parameter_au3_yr2(primary_mass_solar)
    semi_major_axis_au = perihelion_distance_au / (1 - eccentricity)
    mean = math.remainder(mean_anomaly_rad, TWO_PI)  # exact, unlike (M + pi) % 2 pi - pi for a tiny M
    since_perihelion = mean / math.sqrt(mu / semi_major_axis_au ** 3)
    position, velocity = state_from_elements(mu, perihelion_distance_au, eccentricity, 0.0, 0.0, 0.0, 0.0)
    position, velocity = universal_step(position, velocity, mu, since_perihelion)
    return math.atan2(position[1], position[0]), math.hypot(position[0], position[1])


def comet_orbital_state(orbit_type, perihelion_distance_au, eccentricity, inclination_deg,
                        arg_periapsis_deg, ascending_node_deg, primary_mass_solar,
                        mean_anomaly_rad=None, parabolic_mean_anomaly_value=None,
                        orbital_period_years=None, eccentric_anomaly_rad=None):
    """
    A star-bound comet's full orbital state -- 3D position, distance from
    its primary, and orbital speed -- at whatever point along its orbit
    `mean_anomaly_rad`/`parabolic_mean_anomaly_value` (whichever applies
    to `orbit_type`) describes.

    Dispatches on `orbit_type` ("elliptical" uses Kepler's equation and a
    finite semi-major axis derived from `perihelion_distance_au`/
    `eccentricity`; "parabolic" uses Barker's equation and an infinite
    semi-major axis) and reuses `orbits.orbital_position_au` for the actual
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
        eccentric_anomaly_rad (float, optional): The eccentric anomaly for
            `mean_anomaly_rad` when the caller solved a batch of them at
            once (`solve_eccentric_anomalies`); solved here when left out.
            Only used for `orbit_type == "elliptical"`.

    Returns:
        dict: `position_x_au`, `position_y_au`, `position_z_au`,
             `distance_au`, `orbital_speed_kms`, and the velocity relative
             to the primary, `velocity_x_au_per_year`,
             `velocity_y_au_per_year`, `velocity_z_au_per_year`
             (`state_vectors.state_from_elements` for the same orbit).

    Raises:
        ValueError: If `orbit_type` isn't `"elliptical"` or `"parabolic"`,
            or the anomaly required for that type wasn't given.
    """
    if orbit_type == "elliptical":
        if mean_anomaly_rad is None:
            raise ValueError("comet_orbital_state: mean_anomaly_rad is required for an elliptical orbit.")
        semi_major_axis_au = perihelion_distance_au / (1 - eccentricity)
        if eccentricity > KEPLER_MAX_ECCENTRICITY:
            true_anomaly_rad, distance_au = _near_parabolic_anomaly_and_distance(
                mean_anomaly_rad, eccentricity, perihelion_distance_au, primary_mass_solar)
        else:
            true_anomaly_rad, distance_au = true_anomaly_and_distance_elliptical(
                mean_anomaly_rad, eccentricity, semi_major_axis_au, eccentric_anomaly_rad
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
    eccentricity_used = 1.0 if orbit_type == "parabolic" else eccentricity
    _position, velocity = state_from_elements(
        gravitational_parameter_au3_yr2(primary_mass_solar), perihelion_distance_au, eccentricity_used,
        math.radians(inclination_deg), math.radians(ascending_node_deg), math.radians(arg_periapsis_deg),
        true_anomaly_rad)

    return {
        "position_x_au": position_x_au,
        "position_y_au": position_y_au,
        "position_z_au": position_z_au,
        "distance_au": distance_au,
        "orbital_speed_kms": orbital_speed_kms,
        "velocity_x_au_per_year": velocity[0],
        "velocity_y_au_per_year": velocity[1],
        "velocity_z_au_per_year": velocity[2],
    }
