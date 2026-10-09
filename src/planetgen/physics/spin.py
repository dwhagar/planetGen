# planetgen/physics/spin.py

"""
Spin and axial tilt (GEN.104)
=============================

Every rotating body's spin: a unit axis (the direction of its angular
velocity, right-handed) and a rate (its rotation period), so the spin
vector is `axis * 2 pi / period`. Drawn by the cascade in "Observational
Kinetics for Rotational Vectors.md" (docs/design/orbital-updates.md
section 6):

1. A body tidally locked to what it orbits turns once an orbit, with its
   axis on the orbit normal (no tilt).
2. Cool main-sequence stars (up to 1.3 solar masses) take their period
   from gyrochronology, `P = a t^n (B - V - c)^b` (Barnes 2007, with the
   document's fit), with 10% log-normal scatter.
3. Hot stars draw an equatorial speed log-normal round 180 km/s, capped
   at 85% of breakup, `sqrt(2/3 G M / R)`.
4. Comets (and every small body over 200 m) draw a period log-normal
   (median e^2.1 h, about 8 h), never under the 2.2-hour spin barrier.
5. Black holes draw a* from Beta(1.4, 3.6) (first-generation collapse).

Tilts: stars, gas giants and untidied moons Rayleigh(15 deg) from their
reference normal; rocky planets the impact-modified
`p(e) ~ sin e (1 + gamma cos^2 e)`; small bodies the YORP pair near 10 or
170 deg. The axis is the reference normal tilted by the obliquity at a
random precession angle.

The reference normal is the body's orbit normal (a planet round its star,
a moon round its planet, a comet round its star), in the same frame as
its orbit; for a star, or anything on its own galactic orbit, it is the
galaxy's pole (0, 0, 1).
"""

import math

from scipy import special

from planetgen.physics import constants
from planetgen.util import draw

GALACTIC_POLE = (0.0, 0.0, 1.0)
"""tuple: The galaxy's rotation axis, the reference normal of stars and
of anything on its own galactic orbit."""

STELLAR_TILT_SIGMA_DEG = 15.0
"""float: Rayleigh sigma of a star's, gas giant's or free moon's obliquity."""

IMPACT_TILT_GAMMA = 3.0
"""float: gamma in the rocky-planet tilt `sin e (1 + gamma cos^2 e)`. The
document leaves it open; 3 makes a near-upright or near-flipped axis four
times as likely per solid angle as one on its side (planetGen's choice;
gamma = 0 is isotropic)."""

YORP_TILT_DEG = (10.0, 170.0)
YORP_TILT_SIGMA_DEG = 8.0
"""Small bodies' obliquity: N(10, 8) or N(170, 8), even odds (YORP)."""

GYRO_A = 0.40
GYRO_N = 0.55
GYRO_B = 0.31
GYRO_C = 0.495
GYRO_SCATTER = 0.10
"""Gyrochronology, `P[days] = a (t[Myr])^n (B - V - c)^b exp(N(0, 0.1))`
(the document's fit to Barnes 2007): the Sun (4.6 Gyr, B - V 0.65) comes
out 23 days."""

COOL_STAR_MAX_MASS_SOL = 1.3
"""float: Heaviest star with a convective envelope deep enough to brake."""

GYRO_MIN_COLOUR_EXCESS = 0.05
"""float: Gyrochronology needs `B - V - c` positive; a star bluer than
c + 0.05 (about 6200 K) spins like a hot star instead."""

HOT_STAR_SPEED_KMS = 180.0
HOT_STAR_SPEED_SIGMA = 0.45
BREAKUP_FRACTION = 0.85
"""Hot stars: `ln v ~ N(ln 180 km/s, 0.45)`, at most 85% of breakup."""

GIANT_SPEED_RANGE_KMS = (1.0, 10.0)
"""tuple: A giant's equatorial speed, log-uniform (the RGB row of the
document's table); expansion has spun it down whatever its mass."""

WHITE_DWARF_PERIOD_MEDIAN_HOURS = 24.0
WHITE_DWARF_PERIOD_SIGMA = 1.0
"""White dwarfs (not in the document): period log-normal round a day,
spanning minutes to weeks as observed (planetGen's choice)."""

SPIN_BARRIER_HOURS = 2.2
SMALL_BODY_PERIOD_MU = 2.1
SMALL_BODY_PERIOD_SIGMA = 0.65
"""Small bodies over 200 m: `P = max(2.2 h, exp(N(2.1, 0.65)) h)`."""

BLACK_HOLE_SPIN_BETA = (1.4, 3.6)
"""tuple: First-generation black-hole spin a* ~ Beta(1.4, 3.6)."""

BLACK_HOLE_MAX_SPIN = 0.998
"""float: The Thorne limit; a* never reaches 1."""


def _unit(v):
    norm = math.sqrt(sum(c * c for c in v))
    return tuple(c / norm for c in v)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def orbit_normal(inclination_deg, ascending_node_deg):
    """The unit normal of an orbit with these elements, in the frame
    `orbits.orbital_position_au` places the body in (Rz(node) Rx(incl)
    applied to +z)."""
    i, node = math.radians(inclination_deg), math.radians(ascending_node_deg)
    return (math.sin(node) * math.sin(i), -math.cos(node) * math.sin(i), math.cos(i))


def tilted_axis(normal, obliquity_deg, precession_rad):
    """`cos e n + sin e (cos phi u + sin phi v)`, with `u`, `v` completing
    an orthonormal frame round the unit `normal`."""
    n = _unit(normal)
    helper = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _unit(_cross(n, helper))
    v = _cross(n, u)
    e = math.radians(obliquity_deg)
    c, s = math.cos(e), math.sin(e)
    cp, sp = math.cos(precession_rad), math.sin(precession_rad)
    return _unit(tuple(c * n[k] + s * (cp * u[k] + sp * v[k]) for k in range(3)))


def obliquity_between(axis, normal):
    """The angle, degrees, between a spin axis and a reference normal."""
    dot = sum(a * b for a, b in zip(_unit(axis), _unit(normal)))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot))))


# --- Tilt draws ----------------------------------------------------------------------------

def rayleigh_tilt_deg(sigma_deg=STELLAR_TILT_SIGMA_DEG):
    """An obliquity from Rayleigh(sigma), folded into [0, 180]."""
    u = draw.random()
    tilt = sigma_deg * math.sqrt(-2.0 * math.log(1.0 - u))
    tilt %= 360.0
    return 360.0 - tilt if tilt > 180.0 else tilt


def impact_tilt_deg(gamma=IMPACT_TILT_GAMMA):
    """An obliquity from `p(e) ~ sin e (1 + gamma cos^2 e)` on [0, 180]:
    in `x = cos e` that is `1 + gamma x^2` on [-1, 1], sampled by inverting
    its CDF `(x + 1) + gamma (x^3 + 1) / 3` (monotone, so bisection)."""
    total = 2.0 + 2.0 * gamma / 3.0
    target = draw.random() * total
    low, high = -1.0, 1.0
    for _ in range(60):
        mid = 0.5 * (low + high)
        if (mid + 1.0) + gamma * (mid ** 3 + 1.0) / 3.0 < target:
            low = mid
        else:
            high = mid
    return math.degrees(math.acos(max(-1.0, min(1.0, 0.5 * (low + high)))))


def yorp_tilt_deg():
    """A small body's obliquity: near 10 or near 170 degrees (YORP)."""
    centre = YORP_TILT_DEG[0] if draw.random() < 0.5 else YORP_TILT_DEG[1]
    return min(180.0, max(0.0, draw.gauss(centre, YORP_TILT_SIGMA_DEG)))


def draw_axis(normal, tilt_deg):
    """`(axis, tilt)`: the spin axis tilted `tilt_deg` from `normal` at a
    random precession angle."""
    return tilted_axis(normal, tilt_deg, draw.uniform(0.0, 2.0 * math.pi)), tilt_deg


# --- Rates ---------------------------------------------------------------------------------

def b_minus_v(temperature_k):
    """A star's B - V colour from its temperature, inverting Ballesteros
    (2012), `T = 4600 (1 / (0.92 x + 1.7) + 1 / (0.92 x + 0.62))`
    (monotone falling in x, so bisection over [-0.4, 2.5])."""
    low, high = -0.4, 2.5
    for _ in range(60):
        mid = 0.5 * (low + high)
        t = 4600.0 * (1.0 / (0.92 * mid + 1.7) + 1.0 / (0.92 * mid + 0.62))
        if t > temperature_k:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def breakup_speed_kms(mass_kg, radius_km):
    """The Roche breakup equatorial speed, `sqrt(2/3 G M / R)`, km/s."""
    return math.sqrt(2.0 / 3.0 * constants.G * mass_kg / (radius_km * 1000.0)) / 1000.0


def period_hours_from_speed(radius_km, speed_kms):
    return 2.0 * math.pi * radius_km / speed_kms / 3600.0


GIANT_CLASSES = frozenset({"0", "I", "Ia", "Iab", "Ib", "II", "III", "IV"})
"""frozenset: Yerkes classes that have left the main sequence (hypergiants
to subgiants)."""


def star_rotation_period_hours(mass_kg, radius_km, temperature_k, age_gy, yerkes_class):
    """
    A star's rotation period, hours, by the cascade: white dwarfs
    (`yerkes_class` "D...") log-normal round a day; giants and supergiants
    (`GIANT_CLASSES`) a slow log-uniform 1 to 10 km/s; cool dwarfs gyrochronology;
    the rest a hot star's log-normal speed under breakup. Never faster
    than 85% of breakup.
    """
    mass_sol = mass_kg / constants.SOLAR_MASS_TO_KG
    cls = (yerkes_class or "").strip()
    cap = BREAKUP_FRACTION * breakup_speed_kms(mass_kg, radius_km)
    if cls.startswith("D"):
        period = WHITE_DWARF_PERIOD_MEDIAN_HOURS * math.exp(draw.gauss(0.0, WHITE_DWARF_PERIOD_SIGMA))
    elif cls in GIANT_CLASSES:
        low, high = (math.log(s) for s in GIANT_SPEED_RANGE_KMS)
        period = period_hours_from_speed(radius_km, math.exp(draw.uniform(low, high)))
    else:
        colour = b_minus_v(temperature_k) - GYRO_C
        if mass_sol <= COOL_STAR_MAX_MASS_SOL and colour > GYRO_MIN_COLOUR_EXCESS:
            age_myr = max(age_gy * 1000.0, 1.0)
            days = GYRO_A * age_myr ** GYRO_N * colour ** GYRO_B * math.exp(draw.gauss(0.0, GYRO_SCATTER))
            period = days * 24.0
        else:
            speed = math.exp(draw.gauss(math.log(HOT_STAR_SPEED_KMS), HOT_STAR_SPEED_SIGMA))
            period = period_hours_from_speed(radius_km, speed)
    return max(period, period_hours_from_speed(radius_km, cap))


def small_body_period_hours():
    """A comet's or asteroid's period, hours, never under the spin barrier."""
    return max(SPIN_BARRIER_HOURS, math.exp(draw.gauss(SMALL_BODY_PERIOD_MU, SMALL_BODY_PERIOD_SIGMA)))


def black_hole_spin():
    """a* ~ Beta(1.4, 3.6) (inverse CDF), below the Thorne limit."""
    return min(BLACK_HOLE_MAX_SPIN, float(special.betaincinv(*BLACK_HOLE_SPIN_BETA, draw.random())))


def black_hole_horizon_period_s(mass_kg, spin):
    """The rotation period of a Kerr black hole's horizon, seconds:
    `2 pi / Omega_H` with `Omega_H = a* c^3 / (2 G M (1 + sqrt(1 - a*^2)))`;
    infinite for a* = 0."""
    if spin <= 0.0:
        return math.inf
    omega = spin * constants.SPEED_OF_LIGHT_M_S ** 3 / (
        2.0 * constants.G * mass_kg * (1.0 + math.sqrt(1.0 - spin * spin)))
    return 2.0 * math.pi / omega


SPIN_FIELDS = ("spin_axis_x", "spin_axis_y", "spin_axis_z", "axial_tilt_deg")
"""tuple: The attributes (and columns) every spinning object gets."""


def set_spin(body, normal, tilt_deg):
    """Sets `body`'s spin axis (tilted `tilt_deg` from `normal` at a random
    precession angle) and `axial_tilt_deg`."""
    axis, tilt = draw_axis(normal, tilt_deg)
    body.spin_axis_x, body.spin_axis_y, body.spin_axis_z = axis
    body.axial_tilt_deg = tilt


def spin_values(body):
    """`body`'s `SPIN_FIELDS` values, `None` for one it hasn't got (an
    object restored from data saved before GEN.104)."""
    return tuple(getattr(body, name, None) for name in SPIN_FIELDS)
