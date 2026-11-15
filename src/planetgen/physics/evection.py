# planetgen/physics/evection.py

"""
The evection resonance: a moon that drifts out to where its orbit turns once a year (GEN.1750)
================================================================================================

The planet's equatorial bulge turns a moon's orbit (its apsides) at `nu = (3/2) J2 n_m (R / a)^2`, which falls steeply with distance (`a^(-7/2)`). Where `nu` equals the planet's
own mean motion round its star, `n_p`, the star's pull and the bulge's lock together: the evection resonance ([R5289]; [R0958]; [R5290]). Its distance follows from the planet's spin (`J2` goes as the spin squared, `kozai_bodies.bulge_precession_rad_s`), so it sits farther out for a fast planet and moves in
as the planet slows.

A moon whose tide-driven drift (or the planet's slowing) carries the resonance across it is either **caught** or **passes**. The resonance's coupling is `eps = (15/4) n_p^2 / n_m` and
the sweep `|d nu / dt|` carries the moon through it; capture is sure when `rho = |d nu / dt| / eps^2` is small. The [R5289] rule (capture likely under 1 km a year of
`da/dt`, unlikely over 10 km a year, at 4.6 Earth radii) comes out as `rho` 0.54 and 5.4 for the Earth's mass and a year's orbit, with no constant fitted: `RHO_SURE` and
`RHO_PASS` below are those two figures, and the moon is caught below their geometric midpoint (`caught`).

**What a caught moon does** (`held_eccentricity`; the hold itself is in `moon_orbits.evolve`). It stays at the resonance and follows it inward as the planet slows (`a` goes as the
spin to the power 4/7), the planet's angular momentum leaving by the resonance to the star's orbit, and its eccentricity rises to the ceiling [R5289] reached (0.5), cut short
where the periapsis would meet the planet. The model's tide history is circular; this is the one place it holds an orbit.

**What it does not do.** The exit from the resonance (dissipation faster than libration; the resonant despinning of [R0958] of the planet) is not modelled: a caught moon stays caught
for the rest of the age. A retrograde moon, an inward-drifting one that the resonance does not sweep outward, and the variation resonance are left out. Stand-ins are rows EV1 to EV3 of
section 6.302 of the gap plan, tagged in the code (`[Stand-in 6.302 EV1, tier T7, tweak L0.]` and so on).
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, kozai_bodies

RHO_SURE = 0.5
"""float: `rho` under which capture is likely: the 1 km a year of [R5289] at 4.6 Earth radii gives `rho` 0.54 ([R5291], p. 21, second-hand; the coupling and
the sweep are [R5290], equations 14 and 33 to 35). [Stand-in 6.302 EV1, tier T7, tweak L0.]"""
RHO_PASS = 5.0
"""float: `rho` over which the moon passes through: their 10 km a year gives 5.4. [Stand-in 6.302 EV1, tier T7, tweak L0.]"""
CAUGHT_BELOW = math.sqrt(RHO_SURE * RHO_PASS)
"""float: The `rho` that divides caught from passing: the geometric midpoint of the two (a moon is caught at or below it). No draw, so no new stream. [Stand-in 6.302 EV1, tier T7, tweak L0.]"""
ECCENTRICITY_CEILING = 0.5
"""float: The eccentricity a moon caught in the evection resonance reaches ([R5289] reached 0.5; the limit is unsettled, [R5291], p. 27 to 29, second-hand).
[Stand-in 6.302 EV3, tier T6, tweak L1.]"""
PERIAPSIS_MARGIN = 2.0
"""float: A caught moon's periapsis stays outside this many planet radii (the eccentricity is `min(ceiling, 1 - margin R / a)`), where a moon is lost to the planet ([R5290] eq. 41 is
1 R). [Stand-in 6.302 EV3, tier T6, tweak L1.]"""
SPIN_EXPONENT = 4.0 / 7.0
"""float: A moon held at the resonance has `a` going as `omega^(4/7)` (from `nu` going as `omega^2 a^(-7/2)`); it stays held for the rest of the age, with its tide still braking the planet and
the planet's angular momentum leaving to the star's orbit. [Stand-in 6.302 EV2, tier T6, tweak L1.]"""


@dataclass(frozen=True)
class Resonance:
    """
    The evection resonance of one planet.

    Attributes:
        giant (bool): Whether the planet is a gas giant (its `J2` per rotation parameter).
        star_rate (float): The planet's mean motion round its star, rad/s.
    """

    giant: bool
    star_rate: float


def distance_m(resonance, spin, mass_kg, radius_m):
    """
    The semi-major axis (m) at which a moon's bulge precession equals the planet's mean motion, for the planet's `spin` (rad/s): `a^(7/2) = (3/2) f omega^2 R^5 / (sqrt(G M) n_p)`.
    Zero for a planet that is not spinning.
    """
    if spin == 0.0 or resonance.star_rate <= 0.0:
        return 0.0
    f = kozai_bodies.J2_PER_ROTATION_PARAMETER_GIANT if resonance.giant else kozai_bodies.J2_PER_ROTATION_PARAMETER_ROCKY
    value = 1.5 * f * spin ** 2 * radius_m ** 5 / (math.sqrt(constants.G * mass_kg) * resonance.star_rate)
    return value ** (2.0 / 7.0)


def resonance_precession(resonance, spin, mass_kg, radius_m, axis_m):
    """The bulge's precession (rad/s) of an orbit of `axis_m` for the planet's `spin`: `kozai_bodies.bulge_precession_rad_s`."""
    return kozai_bodies.bulge_precession_rad_s(resonance.giant, spin, mass_kg, radius_m, axis_m)


def sweep_ratio(resonance, mass_kg, axis_m, nu_rate):
    """`rho = |d nu / dt| / eps^2` for a moon at `axis_m` that the resonance crosses at the rate `nu_rate` (rad/s per s)."""
    moon_rate = math.sqrt(constants.G * mass_kg / axis_m ** 3)
    coupling = 3.75 * resonance.star_rate ** 2 / moon_rate
    return abs(nu_rate) / coupling ** 2


def caught(rho):
    """Whether a moon crossed at the sweep ratio `rho` is held by the resonance."""
    return rho <= CAUGHT_BELOW


def nu_rate(spin, spin_rate, axis_m, axis_rate_m_s, nu):
    """`d nu / dt` from the spin's and the axis's rates: `nu` goes as `omega^2 a^(-7/2)`, so `nu (2 d ln omega / dt - 7/2 d ln a / dt)`."""
    return nu * (2.0 * spin_rate / spin - 3.5 * axis_rate_m_s / axis_m)


def held_eccentricity(axis_m, radius_m):
    """[Stand-in 6.302 EV3, tier T6, tweak L1.] The eccentricity of a moon held at the resonance at `axis_m`: the ceiling, or less where its periapsis would come within `PERIAPSIS_MARGIN` radii of the planet's center."""
    return max(0.0, min(ECCENTRICITY_CEILING, 1.0 - PERIAPSIS_MARGIN * radius_m / axis_m))
