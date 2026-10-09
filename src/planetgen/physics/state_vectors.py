# planetgen/physics/state_vectors.py

"""
State vectors and orbital elements
==================================

The two ways to describe where a body is on a two-body orbit, and the
conversion both ways between them:

- **Elements**: the orbit's shape and orientation (periapsis distance,
  eccentricity, inclination, ascending node, argument of periapsis) and
  where the body is on it (true anomaly).
- **State vector**: the body's position and velocity relative to what it
  orbits.

`state_from_elements` and `elements_from_state` work for every conic: an
ellipse (`e < 1`, closed), a parabola (`e == 1`) and a hyperbola (`e > 1`),
so one pair of functions serves a planet, a comet and a flyby alike. They
take any consistent units: positions in `L`, velocities in `L/T` and `mu`
in `L^3/T^2` (metres and seconds, or AU and years). Angles are radians.
The rotation is the one `orbits.orbital_position_au` uses, so the two agree
for the same elements.

Degenerate cases follow the usual conventions: an orbit with no tilt
(`|node| = 0`) has its ascending node at 0, a circular orbit has its
argument of periapsis at 0 and its true anomaly measured from the ascending
node (or +x with no tilt).

The orbital update (`docs/design/orbital-updates.md`, section 10) refreshes
a body's elements from its state vector with `elements_from_state`, and
`closed_orbit_points` samples the ellipse for drawing.
"""

import math

TWO_PI = 2.0 * math.pi

PARABOLIC_TOLERANCE = 1.0e-12
"""float: How near 1 an eccentricity is to count as parabolic."""

CIRCULAR_TOLERANCE = 1.0e-10
"""float: How near 0 an eccentricity is to count as circular (and how near
0 the node vector is to count as untilted)."""


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def _scale(a, k):
    return (a[0] * k, a[1] * k, a[2] * k)


def _angle_about(a, b, axis):
    """The angle from `a` to `b` turning about the unit `axis`, in [0, 2 pi)."""
    return math.atan2(_dot(_cross(a, b), axis), _dot(a, b)) % TWO_PI


def _rotate(vector, inclination, ascending_node, arg_periapsis):
    """Perifocal (x toward periapsis, z along the angular momentum) to the
    reference frame: Rz(node) Rx(inclination) Rz(arg of periapsis)."""
    x, y, z = vector
    cw, sw = math.cos(arg_periapsis), math.sin(arg_periapsis)
    x, y = x * cw - y * sw, x * sw + y * cw
    ci, si = math.cos(inclination), math.sin(inclination)
    y, z = y * ci - z * si, y * si + z * ci
    cn, sn = math.cos(ascending_node), math.sin(ascending_node)
    return (x * cn - y * sn, x * sn + y * cn, z)


def state_from_elements(mu, periapsis_distance, eccentricity, inclination, ascending_node, arg_periapsis,
                        true_anomaly):
    """
    The position and velocity, relative to the primary, of a body on the
    orbit given by the elements, at `true_anomaly`.

    Args:
        mu (float): The primary's gravitational parameter G * M.
        periapsis_distance (float): Distance at closest approach (for a
            circle, its radius). Positive.
        eccentricity (float): 0 for a circle, below 1 an ellipse, 1 a
            parabola, above 1 a hyperbola. Not negative.
        inclination, ascending_node, arg_periapsis, true_anomaly (float):
            Radians.

    Returns:
        tuple: `(position, velocity)`, each an (x, y, z) tuple.

    Raises:
        ValueError: For a non-positive `mu` or periapsis, a negative
            eccentricity, or a true anomaly beyond a hyperbola's asymptote.
    """
    if not mu > 0.0 or not math.isfinite(mu):
        raise ValueError(f"mu must be positive, got {mu!r}")
    if not periapsis_distance > 0.0 or not math.isfinite(periapsis_distance):
        raise ValueError(f"the periapsis distance must be positive, got {periapsis_distance!r}")
    if not eccentricity >= 0.0 or not math.isfinite(eccentricity):
        raise ValueError(f"the eccentricity must not be negative, got {eccentricity!r}")
    p = periapsis_distance * (1.0 + eccentricity)
    denominator = 1.0 + eccentricity * math.cos(true_anomaly)
    if denominator <= 0.0:
        raise ValueError(f"true anomaly {true_anomaly!r} is beyond the asymptote of an orbit with e = {eccentricity!r}")
    r = p / denominator
    cos_nu, sin_nu = math.cos(true_anomaly), math.sin(true_anomaly)
    speed_scale = math.sqrt(mu / p)
    position = _rotate((r * cos_nu, r * sin_nu, 0.0), inclination, ascending_node, arg_periapsis)
    velocity = _rotate((-speed_scale * sin_nu, speed_scale * (eccentricity + cos_nu), 0.0),
                       inclination, ascending_node, arg_periapsis)
    return position, velocity


def elements_from_state(position, velocity, mu):
    """
    The orbit through `position` and `velocity` (relative to the primary):
    the osculating elements.

    Returns:
        dict: `periapsis_distance`, `eccentricity`, `inclination`,
        `ascending_node`, `arg_periapsis`, `true_anomaly` (radians, angles
        in [0, 2 pi), inclination in [0, pi]), `semi_major_axis` (negative
        for a hyperbola, `inf` for a parabola), `period` (`None` unless
        elliptical), `kind` ("elliptical", "parabolic" or "hyperbolic") and
        `specific_energy`.

    Raises:
        ValueError: For a non-positive `mu`, a non-finite vector, or a
            body moving straight at or away from the primary (no orbital
            plane) or sitting on it.
    """
    if not mu > 0.0 or not math.isfinite(mu):
        raise ValueError(f"mu must be positive, got {mu!r}")
    for value in (*position, *velocity):
        if not math.isfinite(value):
            raise ValueError(f"the state vector must be finite, got {tuple(position)!r}, {tuple(velocity)!r}")
    r = _norm(position)
    h = _cross(position, velocity)
    h_norm = _norm(h)
    if r == 0.0 or h_norm <= 1.0e-12 * r * _norm(velocity):
        raise ValueError("the body has no orbital plane (it sits on its primary or moves along the line to it)")
    h_hat = _scale(h, 1.0 / h_norm)

    eccentricity_vector = tuple(c / mu - p / r for c, p in zip(_cross(velocity, h), position))
    eccentricity = _norm(eccentricity_vector)
    p = h_norm * h_norm / mu
    periapsis_distance = p / (1.0 + eccentricity)
    inclination = math.acos(max(-1.0, min(1.0, h_hat[2])))

    node_vector = (-h_hat[1], h_hat[0], 0.0)
    node_norm = _norm(node_vector)
    if node_norm < CIRCULAR_TOLERANCE:
        ascending_node = 0.0
        reference = (1.0, 0.0, 0.0)
    else:
        ascending_node = math.atan2(h_hat[0], -h_hat[1]) % TWO_PI
        reference = _scale(node_vector, 1.0 / node_norm)

    position_hat = _scale(position, 1.0 / r)
    if eccentricity < CIRCULAR_TOLERANCE:
        arg_periapsis = 0.0
        true_anomaly = _angle_about(reference, position_hat, h_hat)
    else:
        eccentricity_hat = _scale(eccentricity_vector, 1.0 / eccentricity)
        arg_periapsis = _angle_about(reference, eccentricity_hat, h_hat)
        true_anomaly = _angle_about(eccentricity_hat, position_hat, h_hat)

    energy = 0.5 * _dot(velocity, velocity) - mu / r
    if abs(eccentricity - 1.0) <= PARABOLIC_TOLERANCE:
        kind, semi_major_axis, period = "parabolic", math.inf, None
    elif eccentricity < 1.0:
        kind = "elliptical"
        semi_major_axis = periapsis_distance / (1.0 - eccentricity)
        period = TWO_PI * math.sqrt(semi_major_axis ** 3 / mu)
    else:
        kind, period = "hyperbolic", None
        semi_major_axis = periapsis_distance / (1.0 - eccentricity)
    return {
        "periapsis_distance": periapsis_distance,
        "eccentricity": eccentricity,
        "inclination": inclination,
        "ascending_node": ascending_node,
        "arg_periapsis": arg_periapsis,
        "true_anomaly": true_anomaly,
        "semi_major_axis": semi_major_axis,
        "period": period,
        "kind": kind,
        "specific_energy": energy,
    }


def equinoctial_from_state(position, velocity, mu):
    """
    The modified equinoctial elements (Walker, Ireland and Owens 1985) of
    the orbit through `position` and `velocity`: `p` (semi-latus rectum),
    `f`, `g` (the eccentricity vector in the equinoctial frame), `h`, `k`
    (the orbit pole), and `L` (true longitude, radians in [0, 2 pi)).

    The GEN.108 guard for circular and equatorial orbits: there the
    classical argument of periapsis and ascending node are undefined
    (`elements_from_state` falls back to conventions), while these
    elements stay smooth and exact through e = 0 and i = 0. They are
    singular only for an exactly retrograde equatorial orbit (i = pi),
    which `ValueError`s.

    Returns:
        dict: `p`, `f`, `g`, `h`, `k`, `L`.
    """
    if not mu > 0.0 or not math.isfinite(mu):
        raise ValueError(f"mu must be positive, got {mu!r}")
    r = _norm(position)
    h_vec = _cross(position, velocity)
    h_norm = _norm(h_vec)
    if r == 0.0 or h_norm <= 1.0e-12 * r * _norm(velocity):
        raise ValueError("the body has no orbital plane (it sits on its primary or moves along the line to it)")
    h_hat = _scale(h_vec, 1.0 / h_norm)
    if 1.0 + h_hat[2] < 1.0e-12:
        raise ValueError("modified equinoctial elements are singular for a retrograde equatorial orbit")
    k_eq = h_hat[0] / (1.0 + h_hat[2])
    h_eq = -h_hat[1] / (1.0 + h_hat[2])
    # The equinoctial frame: f_hat, g_hat in the orbit plane, w_hat = h_hat.
    s2 = 1.0 + h_eq * h_eq + k_eq * k_eq
    f_hat = ((1.0 - k_eq * k_eq + h_eq * h_eq) / s2, 2.0 * k_eq * h_eq / s2, -2.0 * k_eq / s2)
    g_hat = (2.0 * k_eq * h_eq / s2, (1.0 + k_eq * k_eq - h_eq * h_eq) / s2, 2.0 * h_eq / s2)
    eccentricity_vector = tuple(c / mu - q / r for c, q in zip(_cross(velocity, h_vec), position))
    return {
        "p": h_norm * h_norm / mu,
        "f": _dot(eccentricity_vector, f_hat),
        "g": _dot(eccentricity_vector, g_hat),
        "h": h_eq,
        "k": k_eq,
        "L": math.atan2(_dot(position, g_hat), _dot(position, f_hat)) % TWO_PI,
    }


def state_from_equinoctial(p, f, g, h, k, true_longitude, mu):
    """The position and velocity of the modified equinoctial elements
    (`equinoctial_from_state`'s inverse)."""
    if not mu > 0.0 or not p > 0.0:
        raise ValueError(f"mu and p must be positive, got {mu!r}, {p!r}")
    cos_l, sin_l = math.cos(true_longitude), math.sin(true_longitude)
    w = 1.0 + f * cos_l + g * sin_l
    if w <= 0.0:
        raise ValueError("the true longitude is beyond the orbit's asymptote")
    r = p / w
    s2 = 1.0 + h * h + k * k
    alpha2 = h * h - k * k
    root = math.sqrt(mu / p)
    position = (r / s2 * (cos_l + alpha2 * cos_l + 2.0 * h * k * sin_l),
                r / s2 * (sin_l - alpha2 * sin_l + 2.0 * h * k * cos_l),
                2.0 * r / s2 * (h * sin_l - k * cos_l))
    velocity = (-root / s2 * (sin_l + alpha2 * sin_l - 2.0 * h * k * cos_l + g - 2.0 * f * h * k + alpha2 * g),
                -root / s2 * (-cos_l + alpha2 * cos_l + 2.0 * h * k * sin_l - f + 2.0 * g * h * k + alpha2 * f),
                2.0 * root / s2 * (h * cos_l + k * sin_l + f * h + g * k))
    return position, velocity


def mean_anomaly_from_true(true_anomaly, eccentricity):
    """The mean anomaly (radians, in [0, 2 pi)) of an elliptical orbit's
    `true_anomaly`; the inverse of `kepler.true_anomaly_and_distance_elliptical`."""
    if not 0.0 <= eccentricity < 1.0:
        raise ValueError(f"only an elliptical orbit has a mean anomaly, got e = {eccentricity!r}")
    eccentric = 2.0 * math.atan2(math.sqrt(1.0 - eccentricity) * math.sin(true_anomaly / 2.0),
                                 math.sqrt(1.0 + eccentricity) * math.cos(true_anomaly / 2.0))
    return (eccentric - eccentricity * math.sin(eccentric)) % TWO_PI


def closed_orbit_points(periapsis_distance, eccentricity, inclination, ascending_node, arg_periapsis, count=128):
    """
    `count` points evenly spaced in true anomaly round a closed (elliptical)
    orbit (its shape does not depend on mu), relative to the primary, starting at periapsis -- the ellipse to
    draw. Raises `ValueError` for an orbit that does not close.
    """
    if not 0.0 <= eccentricity < 1.0:
        raise ValueError(f"only an elliptical orbit is closed, got e = {eccentricity!r}")
    if count < 3:
        raise ValueError(f"an orbit needs at least 3 points, got {count!r}")
    return [state_from_elements(1.0, periapsis_distance, eccentricity, inclination, ascending_node, arg_periapsis,
                                TWO_PI * k / count)[0] for k in range(count)]
