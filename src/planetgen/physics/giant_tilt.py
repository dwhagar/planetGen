# planetgen/physics/giant_tilt.py

"""
The tilt of an ice giant from a giant impact, and the planes that follow it (GEN.1053, G1 of docs/design/systems/weird-worlds.md, registry 6.78)
===================================================================================================================================================

Uranus lies on its side (98 degrees) with its rings and its five large moons in the equator. The explanation of [R2300, R3149, R4118] is a giant impact: the impactor adds its orbital angular momentum to the spin of the planet, and the
planet's axis turns to the sum. This module is that sum and the tail of the draw it makes.

**The angular momentum an impact brings.** Two bodies in contact at the sum of their radii, `R_c`, at the angle `theta` from head-on, arrive at
`v_c = sqrt(v_inf^2 + 2 G (M + m) / R_c)` and carry `L = mu v_c R_c sin(theta)`, with `mu = M m / (M + m)`. [R2300] run an impactor of 1, 2 and 3
Earth masses on a proto-Uranus of `14.536 - m` at `v_inf` of 5 km/s, with L from 1 to 10e36 kg m^2/s, and find that 2 and 3 Earth masses spin the planet at least as fast as
Uranus (17.24 h) from L of 2e36, and 1 Earth mass never does. A 2 Earth-mass impactor on a 12.5 Earth-mass target grazing gives 7e36 here, inside their range (a 1 Earth-mass one gives 3.7e36, and in their runs does not spin the planet up: the module's contact point is cruder than their SPH bodies).

**The turn of the axis.** The planet spins with `L_0 = I omega_0` before and has `L_0 + f L` after (`f`, `RETAINED_SHARE`, is the share the planet keeps: debris and an
ejected atmosphere carry the rest). With `rho = f L / L_0` and `alpha` the angle between the two vectors the new axis is at

    cos(tilt) = (1 + rho cos(alpha)) / sqrt(1 + rho^2 + 2 rho cos(alpha))

from the old one. An impact that brings less than the planet's own spin (`rho < 1`) turns it by at most `asin(rho)` (Neptune's 28 degrees takes `rho` of 0.47 or more);
one that brings more can turn it to any angle, up to 180. The direction of `L` is isotropic in the draw (`cos(alpha)` uniform from -1 to 1) and the impact parameter is
uniform in area (`sin(theta) = sqrt(u)`), the geometry of a random hit. The rate is what the impactors of the recorded mergers give, not a rate of its own; the share of
ice giants of 5 to 30 Earth masses that have such a merger is the stand-in `MERGER_SHARE`.

**The planes.** The rings and the moons form from the debris the impact puts around the planet, which orbits in the new equator: the plane of both is the equator,
so their tilt to the orbit is the planet's own (R3 and G1 of the design). [R4118] (quoted by [R3149]) make an equatorial debris disk of
1 to 3 percent of the planet's mass; [R2300] find the debris clouds "quite spherical" (axis ratios 0.7 to 1), so the module takes the equator by
construction and does not model the disk. [R3149] add that the prograde moons need an obliquity of 4 to 8 degrees or more before the last impact, so Uranus took at least two impacts;
the module is one impact, and the earlier tilt is the `initial_tilt_deg` that the caller gives (the draw of `spin.rayleigh_tilt_deg`).

Use:

    ratio = angular_momentum_ratio(mass_earth, radius_earth, impactor_earth, spin_period_h, u_parameter)
    tilt = tilt_deg(ratio, cos_alpha)
    period_h = final_period_h(mass_earth, radius_earth, ratio, cos_alpha)
"""

import math

from planetgen.physics import constants
from planetgen.util import draw

MOMENT_OF_INERTIA_FACTOR = 0.225
"""float: `I / (M R^2)` of an ice giant: Uranus 0.225, Neptune 0.23 (recalled from the [R1890] interior models, row GT2)."""

RETAINED_SHARE = 0.9
"""float: The share of the impact's angular momentum the planet keeps (the rest leaves with the debris and the ejected atmosphere; [R2300] find at least
90 percent of the atmosphere stays and the debris is a few percent of the mass, row GT2)."""

DEFAULT_IMPACT_SPEED_KM_S = 5.0
"""float: The speed at infinity, km/s ([R2300: Sec. 2.3], tested between 1 and 9: "does not significantly affect the results", row GT1)."""

PRE_IMPACT_SPIN_BREAKUP = 0.2
"""float: The spin of the planet before the impact over its breakup spin: 0.2 is a day of 15 hours for an ice giant (Uranus has 0.17 of breakup today, 17.24 h, and Neptune 0.18).
[R2300] start from none ("the impacts we consider are violent enough to dominate"); a spin only makes a larger impactor needed (row GT2)."""

MERGER_SHARE = 0.55
"""float: The share of ice giants of 5 to 30 Earth masses that have a recorded giant merger of the outer chain in their history. Stand-in: tuned so that the share of
all such giants tilted over 60 degrees is 0.2 (0.36 of the mergers do it for the impactors of `IMPACTOR_SHARE`) (the Solar System has one of five: Uranus; Neptune's 28 degrees is the
Rayleigh body of the tilt, row GT3, source tier T4, tweak level L3)."""

IMPACTOR_SHARE = (0.02, 0.25)
"""tuple: The impactor's mass over the planet's, log-uniform: 0.3 to 3.6 Earth masses on Uranus's 14.5 ([R2300]: 1, 2 and 3 Earth masses; [R3149]:
giant impacts "rather common" in the growth of the ice giants, no distribution; row GT3)."""

MASS_RANGE_EARTH = (5.0, 30.0)
"""tuple: The planet masses the module draws for, Earth masses (the design table, row GT3)."""

LARGE_TILT_DEG = 60.0
"""float: The tilt over which the design calls a tilt large (G1)."""

_G = constants.G
_R_EARTH_M = constants.EARTH_RADIUS_KM * 1.0e3
_M_EARTH = constants.EARTH_MASS_TO_KG


def contact_radius_m(radius_earth, mass_earth, impactor_earth):
    """The sum of the radii at contact, m: the target's and an impactor of the same density."""
    return radius_earth * _R_EARTH_M * (1.0 + (impactor_earth / mass_earth) ** (1.0 / 3.0))


def impact_angular_momentum(mass_earth, radius_earth, impactor_earth, sin_theta, v_inf_km_s=DEFAULT_IMPACT_SPEED_KM_S):
    """`L = mu v_c R_c sin(theta)`, kg m^2/s: the angular momentum of an impact at contact, `sin_theta` from 0 (head-on) to 1 (grazing)."""
    m, n = mass_earth * _M_EARTH, impactor_earth * _M_EARTH
    contact = contact_radius_m(radius_earth, mass_earth, impactor_earth)
    speed = math.sqrt((v_inf_km_s * 1.0e3) ** 2 + 2.0 * _G * (m + n) / contact)
    return m * n / (m + n) * speed * contact * sin_theta


def spin_angular_momentum(mass_earth, radius_earth, spin_fraction=PRE_IMPACT_SPIN_BREAKUP):
    """`L_0 = I omega_0`, kg m^2/s, for a spin of `spin_fraction` of the breakup rate `sqrt(G M / R^3)`."""
    mass, radius = mass_earth * _M_EARTH, radius_earth * _R_EARTH_M
    return MOMENT_OF_INERTIA_FACTOR * mass * radius ** 2 * spin_fraction * math.sqrt(_G * mass / radius ** 3)


def angular_momentum_ratio(mass_earth, radius_earth, impactor_earth, sin_theta, v_inf_km_s=DEFAULT_IMPACT_SPEED_KM_S,
                           spin_fraction=PRE_IMPACT_SPIN_BREAKUP):
    """`rho = f L / L_0`: the angular momentum the planet keeps from the impact over what it spun with."""
    return (RETAINED_SHARE * impact_angular_momentum(mass_earth, radius_earth, impactor_earth, sin_theta, v_inf_km_s)
            / spin_angular_momentum(mass_earth, radius_earth, spin_fraction))


def tilt_deg(ratio, cos_alpha):
    """The angle the axis turns through, degrees, when the planet gains `ratio` of its spin in a direction at `cos_alpha` to its axis."""
    along = 1.0 + ratio * cos_alpha
    length = math.sqrt(max(1.0 + ratio * ratio + 2.0 * ratio * cos_alpha, 1.0e-30))
    return math.degrees(math.acos(max(-1.0, min(1.0, along / length))))


def probability_of_tilt_over(ratio, threshold_deg=LARGE_TILT_DEG, steps=2000):
    """The chance, for a direction isotropic on the sphere, of a tilt over `threshold_deg` when the planet gains `ratio` of its spin (midpoint rule on `cos_alpha`)."""
    count = sum(1 for i in range(steps) if tilt_deg(ratio, -1.0 + 2.0 * (i + 0.5) / steps) > threshold_deg)
    return count / steps


def impactor_mass_earth(mass_earth, u):
    """An impactor mass, Earth masses, log-uniform over `IMPACTOR_SHARE` of the planet's mass for `u` in [0, 1]."""
    low, high = IMPACTOR_SHARE
    return mass_earth * low * (high / low) ** u


def draw_tilt_deg(mass_earth, radius_earth, u_mass, u_parameter, u_direction):
    """
    The tilt a giant merger gives, degrees, from three unit draws: the impactor's mass (`impactor_mass_earth`), the impact parameter (`sin(theta) = sqrt(u)`, uniform in
    area) and the direction of the angular momentum (`cos(alpha) = 2 u - 1`).
    """
    ratio = angular_momentum_ratio(mass_earth, radius_earth, impactor_mass_earth(mass_earth, u_mass), math.sqrt(u_parameter))
    return tilt_deg(ratio, 2.0 * u_direction - 1.0)


def share_tilted_over(mass_earth, radius_earth, threshold_deg=LARGE_TILT_DEG, grid=24):
    """The share of mergers of a planet that tilt it over `threshold_deg`: the mean of `probability_of_tilt_over` over the impactor masses and the impact parameters (a grid)."""
    total = 0.0
    for i in range(grid):
        impactor = impactor_mass_earth(mass_earth, (i + 0.5) / grid)
        for j in range(grid):
            sin_theta = math.sqrt((j + 0.5) / grid)
            ratio = angular_momentum_ratio(mass_earth, radius_earth, impactor, sin_theta)
            total += probability_of_tilt_over(ratio, threshold_deg, steps=200)
    return total / (grid * grid)


def large_tilt_chance(mass_earth, radius_earth, has_merger=True):
    """The chance an ice giant is tilted over 60 degrees: 0 with no recorded merger, else `MERGER_SHARE` times `share_tilted_over`."""
    if not has_merger:
        return 0.0
    return MERGER_SHARE * share_tilted_over(mass_earth, radius_earth)


def final_period_h(mass_earth, radius_earth, ratio, cos_alpha, spin_fraction=PRE_IMPACT_SPIN_BREAKUP):
    """The spin period after the impact, hours: the length of the vector sum `L_0 (1 + rho)` over `I`, held at the breakup rate."""
    mass, radius = mass_earth * _M_EARTH, radius_earth * _R_EARTH_M
    breakup = math.sqrt(_G * mass / radius ** 3)
    omega = min(spin_fraction * breakup * math.sqrt(1.0 + ratio * ratio + 2.0 * ratio * cos_alpha), breakup)
    return 2.0 * math.pi / omega / 3600.0




def combined_tilt_deg(initial_tilt_deg, turn_deg, cos_azimuth):
    """
    The obliquity after the axis turns through `turn_deg` from where it stood at `initial_tilt_deg` to the orbit normal, degrees: the turn is about an axis at
    an azimuth round the old one with cosine `cos_azimuth`, so `cos(e) = cos(e0) cos(t) + sin(e0) sin(t) cos(phi)`. From an upright planet the obliquity is the turn.
    """
    e0, t = math.radians(initial_tilt_deg), math.radians(turn_deg)
    cos_e = math.cos(e0) * math.cos(t) + math.sin(e0) * math.sin(t) * cos_azimuth
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_e))))


def tilt_after_merger_deg(mass_kg, radius_km, initial_tilt_deg):
    """
    The obliquity of a giant planet, degrees, given the tilt it had before any giant impact (GEN.1165): a giant of `MASS_RANGE_EARTH` has a recorded merger with the chance
    `MERGER_SHARE` (the share of the ice giants whose history holds one) and then takes the turn of `draw_tilt_deg` (impactor mass, impact parameter and direction drawn) on
    top of the earlier tilt; any other giant keeps it. A giant over the range is not turned because the impact brings a small share of its spin (`angular_momentum_ratio`
    stays under 1), so it is not drawn at all. Draws nothing outside the range.
    """
    mass_earth = mass_kg / _M_EARTH
    if not MASS_RANGE_EARTH[0] <= mass_earth <= MASS_RANGE_EARTH[1] or draw.random() >= MERGER_SHARE:
        return initial_tilt_deg
    turn = draw_tilt_deg(mass_earth, radius_km / constants.EARTH_RADIUS_KM, draw.random(), draw.random(), draw.random())
    return combined_tilt_deg(initial_tilt_deg, turn, math.cos(draw.uniform(0.0, 2.0 * math.pi)))
