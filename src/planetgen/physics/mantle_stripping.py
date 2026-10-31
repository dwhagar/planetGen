# planetgen/physics/mantle_stripping.py

"""
Collisional mantle stripping: how a giant impact leaves a planet richer in iron (GEN.1040; oddity W10 of docs/design/systems/weird-worlds.md)
=========================================================================================================================================

Mercury is about 70 percent iron by mass, Kepler-107 c about 70 percent (Bonomo et al. 2019), where an Earth-like planet is a third. A head-on,
fast impact between two differentiated planets blows off mostly the outer rock, because the shock is weaker at depth and the iron is the
heaviest and the most tightly bound. This module is that one step as pure functions, from the bodies and the speed to the mass and core mass
fraction of the largest remnant. It reads no database, draws no random number and moves nothing; GEN.1047 calls it for the planets the oddity
table picks (row W10, `physics/data/oddities.json`), with the impacts the planet architecture recorded.

**The law** (Marcus, Stewart, Sasselov and Hernquist 2010, ApJL 712, L73, read; arXiv:0907.0234). With `M_tot` the mass of both bodies, `mu`
their reduced mass and `V` the impact speed:

- `Q_R = 0.5 mu V^2 / M_tot`, the specific energy of the impact, and `Q*_RD = 1e-4 R_C1^1.2 V^0.8` (erg/g; cm and cm/s) the catastrophic
  disruption threshold of Stewart and Leinhardt (2009) in the gravity regime, where `R_C1` is the radius of `M_tot` at a density of 1 g/cm3
  (equation 4). Marcus found their simulated threshold below it by less than a factor of 2; `Q_STAR_SCALE` is that factor, taken from the
  one published giant-impact run for which the masses, speed and result are all printed (Bonomo 2019, Fig. 3: two 10.5 Earth mass bodies,
  62.5 km/s, a remnant of 8.2 Earth masses and 64 percent iron); 0.6 gives 8.28 Earth masses and 64.2 percent (`MT4`).
- The largest remnant is `M_lr / M_tot = 1 - 0.5 x`, `x = Q_R / Q*_RD` (equation 2), for `x` from 0 to 2; at `x = 1` half the mass is lost.
- The iron fraction of that remnant is `M_Fe / M_lr = 0.33 + 0.25 x^1.65` (equation 3), for bodies that began 33 percent iron.

**What is added here.** Equation 3 is for one starting composition. For another, the same fit is read as how much iron and how much rock a body of
that fit keeps: from `M_lr` and equation 3, a body of 33 percent iron keeps the share `r_Fe = (0.33 + 0.25 x^1.65)(1 - x/2) / 0.33` of its iron
and `r_rock = (0.67 - 0.25 x^1.65)(1 - x/2) / 0.67` of its rock, and a body that began with the core fraction `f` keeps `r_Fe` of its iron and
`r_rock` of its rock (`MT5`); it is the 33 percent fit exactly where Marcus fitted it, and is never worse than the conservation of iron.

**What is not modelled.** Only head-on and near head-on impacts strip. From 30 to 45 degrees the impact is a hit-and-run (the target is left
nearly whole, Marcus 2010 section 3.2), so `strip_probability` gives the share of random impacts that are inside `IMPACT_ANGLE_DEG`
and the caller applies it to each impact it counts. The remnant is not cooled or relaxed, a mass of less than a tenth of the target's is not
tracked as debris, and nothing is said about the composition of the lost mass.

Use:

    result = impact(mass_kg, core_fraction, projectile_mass_kg, speed_ms)           # one impact
    result = strip(mass_kg, core_fraction, [(projectile_kg, speed_ms), ...])        # a series, each on the last remnant
    speed = speed_for_core_fraction(mass_kg, projectile_kg, core_fraction, 0.70)    # what it takes to reach Mercury's iron
"""

import dataclasses
import math

from planetgen.physics import constants
from planetgen.util import detmath

Q_STAR_SCALE = 0.6
"""float: The factor on the Stewart and Leinhardt threshold for planet-sized bodies (MT4)."""

X_MAX = 2.0
"""float: The largest `x = Q_R / Q*_RD` equation 2 holds for; at 2 nothing is left (Marcus 2010 section 3.1)."""

MAX_CORE_MASS_FRACTION = 0.85
"""float: The core fraction no remnant is allowed to pass (MT6): the most iron-rich rocky planet read is about 0.8 (Plotnykov and Valencia 2020)."""

CORE_FRACTION_CEILING = math.nextafter(MAX_CORE_MASS_FRACTION, 0.0)
"""float: What a remnant past `MAX_CORE_MASS_FRACTION` is held to: the largest float under it, so the strict limit of the W10 row (`lt 0.85`) holds."""

IMPACT_ANGLE_DEG = 30.0
"""float: The impact angle, from head-on, beyond which a collision is a hit-and-run and strips nothing (MT7)."""

EARTH_LIKE_CORE_FRACTION = 0.33
"""float: The core fraction of the bodies Marcus simulated (67 percent forsterite, 33 percent iron)."""

_IRON_FIT = (0.33, 0.25, 1.65)
_REFERENCE_DENSITY = 1000.0
_ACCRETION, _EROSION, _DISRUPTION, _DISPERSAL = "accretion", "erosion", "disruption", "dispersal"
REGIMES = (_ACCRETION, _EROSION, _DISRUPTION, _DISPERSAL)
"""tuple: What an impact does to the target: it grows (the remnant is heavier than the target), it loses mass, it loses half or more of the total,
or the remnant is nothing (`x` at 2 or more)."""


@dataclasses.dataclass(frozen=True)
class Remnant:
    """The largest remnant of one impact (or of a series: `steps` then holds each)."""
    mass_kg: float
    core_mass_fraction: float
    mass_lost_kg: float
    energy_ratio: float
    regime: str
    clamped: bool
    steps: tuple = ()


def reduced_energy(mass_kg, projectile_kg, speed_ms):
    """`Q_R = 0.5 mu V^2 / M_tot` in J/kg: the kinetic energy of the reduced mass over the total mass."""
    total = mass_kg + projectile_kg
    return 0.5 * (mass_kg * projectile_kg / total) * speed_ms ** 2 / total


def disruption_threshold(total_kg, speed_ms):
    """`Q*_RD` in J/kg (equation 4 with `Q_STAR_SCALE`): 1e-4 R^1.2 V^0.8 erg/g, which in SI is 1e-4 R^1.2 V^0.8 J/kg, `R` the radius in metres
    of `total_kg` at 1000 kg/m3."""
    radius = (3.0 * total_kg / (4.0 * math.pi * _REFERENCE_DENSITY)) ** (1.0 / 3.0)
    return Q_STAR_SCALE * 1e-4 * radius ** 1.2 * speed_ms ** 0.8


def energy_ratio(mass_kg, projectile_kg, speed_ms):
    """`x = Q_R / Q*_RD` of an impact of `projectile_kg` on `mass_kg` at `speed_ms`."""
    return reduced_energy(mass_kg, projectile_kg, speed_ms) / disruption_threshold(mass_kg + projectile_kg, speed_ms)


def speed_for_energy_ratio(mass_kg, projectile_kg, x):
    """The impact speed that gives the energy ratio `x`: `Q_R` goes as `V^2` and `Q*_RD` as `V^0.8`, so `x` as `V^1.2`."""
    unit = energy_ratio(mass_kg, projectile_kg, 1.0)
    return (x / unit) ** (1.0 / 1.2)


def _fit_iron(x):
    """Equation 3: the iron fraction of the remnant of two 33 percent iron bodies."""
    a, b, c = _IRON_FIT
    return a + b * x ** c


def retention(x):
    """`(r_Fe, r_rock)`: the share of its iron and of its rock a body keeps at the energy ratio `x` (see the module docstring). Both are held to
    0 to 1, and `x` to 0 to `X_MAX`."""
    x = min(max(x, 0.0), X_MAX)
    left = 1.0 - 0.5 * x
    iron = _fit_iron(x) * left / EARTH_LIKE_CORE_FRACTION
    rock = (1.0 - _fit_iron(x)) * left / (1.0 - EARTH_LIKE_CORE_FRACTION)
    return min(max(iron, 0.0), 1.0), min(max(rock, 0.0), 1.0)


def core_fraction_after(core_fraction, x):
    """The core fraction of the remnant of a body that began with `core_fraction` when the energy ratio is `x`. Equation 3 where `core_fraction` is
    0.33; 1 once the rock is gone. The remnant's mass is equation 2's whatever its make-up, so the iron in it is held to all the iron the body had
    (`core_fraction / (1 - x/2)` of the remnant), which binds a body poor in iron."""
    iron, rock = retention(x)
    kept_iron = core_fraction * iron
    kept_rock = (1.0 - core_fraction) * rock
    if kept_iron + kept_rock <= 0.0:
        return core_fraction
    fraction = kept_iron / (kept_iron + kept_rock)
    return min(fraction, core_fraction / (1.0 - 0.5 * min(x, X_MAX - 1e-9)))


def regime(x, mass_kg, projectile_kg):
    """What the impact does (`REGIMES`): the target grows while `1 - x/2` of the total is above its own mass, and 'dispersal' once `x` is at 2."""
    if x >= X_MAX:
        return _DISPERSAL
    if x >= 1.0:
        return _DISRUPTION
    if 1.0 - 0.5 * x > mass_kg / (mass_kg + projectile_kg):
        return _ACCRETION
    return _EROSION


def impact(mass_kg, core_fraction, projectile_kg, speed_ms, projectile_core_fraction=None):
    """
    The largest remnant of `projectile_kg` striking `mass_kg` head-on at `speed_ms` (the speed at contact: the speed at infinity and the
    escape speed of the pair added in quadrature, `contact_speed_ms`).

    The pair is taken as one body of the mass-weighted core fraction (the projectile has `projectile_core_fraction`, default the target's).
    The remnant's core fraction is held just under `MAX_CORE_MASS_FRACTION` (`CORE_FRACTION_CEILING`); `clamped` says it was. At `x` of 2 or more the remnant has no mass.

    Raises:
        ValueError: For a mass, a core fraction outside 0 to 1, or a speed that is not positive.
    """
    if mass_kg <= 0.0 or projectile_kg <= 0.0 or speed_ms <= 0.0:
        raise ValueError("the masses and the speed must be positive")
    projectile_core = core_fraction if projectile_core_fraction is None else projectile_core_fraction
    if not (0.0 <= core_fraction <= 1.0 and 0.0 <= projectile_core <= 1.0):
        raise ValueError("a core fraction is a number from 0 to 1")
    total = mass_kg + projectile_kg
    combined = (core_fraction * mass_kg + projectile_core * projectile_kg) / total
    x = energy_ratio(mass_kg, projectile_kg, speed_ms)
    kind = regime(x, mass_kg, projectile_kg)
    if kind == _DISPERSAL:
        return Remnant(0.0, combined, total, x, kind, False)
    remnant = total * (1.0 - 0.5 * min(x, X_MAX))
    fraction = core_fraction_after(combined, x)
    clamped = fraction > MAX_CORE_MASS_FRACTION
    return Remnant(remnant, min(fraction, CORE_FRACTION_CEILING), mass_kg - remnant, x, kind, clamped)


def strip(mass_kg, core_fraction, impacts, projectile_core_fraction=None):
    """
    The remnant of a series of impacts, each on the remnant of the one before: `impacts` is a sequence of `(projectile_kg, speed_ms)`. An empty
    series gives the body unchanged. A series stops at a dispersal. `mass_lost_kg` is the whole body's loss (the mass before less the mass after,
    which is negative when the body grew) and `steps` holds each impact's `Remnant`.
    """
    mass, core = mass_kg, core_fraction
    steps = []
    for projectile_kg, speed_ms in impacts:
        step = impact(mass, core, projectile_kg, speed_ms, projectile_core_fraction)
        steps.append(step)
        mass, core = step.mass_kg, step.core_mass_fraction
        if step.regime == _DISPERSAL:
            break
    last = steps[-1] if steps else None
    return Remnant(mass, core, mass_kg - mass, last.energy_ratio if last else 0.0, last.regime if last else _ACCRETION,
                   any(step.clamped for step in steps), tuple(steps))


def energy_ratio_for_core_fraction(core_fraction, target_fraction):
    """The energy ratio `x` at which a body that began with `core_fraction` is left with `target_fraction` (found by bisection: the core fraction
    rises with `x`). 0 if it already has that much; `X_MAX` if no `x` below it gets there."""
    if target_fraction <= core_fraction:
        return 0.0
    if core_fraction_after(core_fraction, X_MAX * (1.0 - 1e-9)) < target_fraction:
        return X_MAX
    low, high = 0.0, X_MAX
    for _ in range(80):
        middle = 0.5 * (low + high)
        if core_fraction_after(core_fraction, middle) < target_fraction:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def speed_for_core_fraction(mass_kg, projectile_kg, core_fraction, target_fraction):
    """The impact speed in m/s that leaves `target_fraction` of iron in the remnant of `projectile_kg` striking `mass_kg`."""
    return speed_for_energy_ratio(mass_kg, projectile_kg, energy_ratio_for_core_fraction(core_fraction, target_fraction))


def escape_speed_ms(mass_kg, projectile_kg, radius_m, projectile_radius_m):
    """The mutual escape speed at contact, `sqrt(2 G (M + m) / (R + r))`."""
    return math.sqrt(2.0 * constants.G * (mass_kg + projectile_kg) / (radius_m + projectile_radius_m))


def contact_speed_ms(mass_kg, projectile_kg, radius_m, projectile_radius_m, speed_at_infinity_ms):
    """The speed at contact: the speed at infinity and the escape speed in quadrature."""
    return detmath.hypot(speed_at_infinity_ms, escape_speed_ms(mass_kg, projectile_kg, radius_m, projectile_radius_m))


def strip_probability(angle_deg=None):
    """The share of random impacts within `angle_deg` of head-on (default `IMPACT_ANGLE_DEG`): the impact parameter is `sin(angle)` of the
    contact distance and impacts are uniform in area, so the share is `sin^2(angle)`."""
    return math.sin(math.radians(IMPACT_ANGLE_DEG if angle_deg is None else angle_deg)) ** 2
