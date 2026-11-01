# planetgen/physics/rogue_moons.py

"""
Whether a moon is carried out with the planet that is thrown out, and what it does there (GEN.1112, RP4)
=========================================================================================================

A rocky rogue "that kept a moon" drew its tidal heat from a log-uniform range (`ROGUE_TIDAL_FLUX_RANGE_W_M2`), and which rogues kept a moon was not asked (`rogue-planet-review.md` row 7). Two
things decide it: whether the moon is still bound when the planet leaves, and how hard it then flexes the planet.

**Survival.** A planet thrown out by a scatterer passes it at the closest approach `r_min` of the series of encounters that unbind it. A moon whose orbital time is shorter than the encounter
time follows its planet through ([R2031], eq. 4):

    a_moon^3 / (G m_planet)  <  r_min^3 / (G (m_planet + m_scatterer))    that is    a_moon < r_min (m_planet / (m_planet + m_scatterer))^(1/3).

Their N-body runs (circular restricted and finite-mass pairs of planets, point masses, a star of one solar mass) give the closest approach in units of the pair's mutual Hill radius `R_H`: for a
planet of 10 Jupiter masses ejecting one of 1 Jupiter mass round a solar-mass star (`m1 / M` about 0.01, `m1 / m2` 10), the median is 0.089 `R_H`, 48 percent of ejections stay beyond 0.1 `R_H` and 82 percent beyond
0.01 `R_H` (their conclusions; the regular moons of the Solar System orbit at about 0.01 `R_H`). `closest_approach_survival` is those two points joined on a logarithmic axis; `closest_approach_hill`
inverts it for a draw. The ejection record of RP1 (GEN.1119) gives the real `r_min` and the criterion above is then exact (`moon_survives`). Closer encounters happen for small
inclinations and weaker scatterers (their sections 4.2.1 to 4.2.4); those runs are not in the fit, which is the weakest part of this module.

**The moon's orbit afterwards.** [R1032] threw Earth-Moon pairs out with a Jupiter: in 4.6 percent of 2,700 runs the pair left bound, the median eccentricity
of the survivors was 0.13 (from 0.05 before) at a median axis of 30 Earth radii. A moon on an eccentric orbit flexes the planet and the flexing is heat; the circularisation term of their
eq. 1 (the constant eccentricity tide, [R3464] in the form of `tidal_heating.peale_power_w`) is `10.5 (k2/Q) G m_moon^2 n R_planet^5 e^2 / a^6` for the mean motion `n` of the pair. It
gives about 5e-3 W/m2 for their median pair with their rheology (Q of 12) and 1e-3 with the rock of `tidal_heating` (Q of 100), both inside the old default's range of 1e-4 to 1e-2 and near Earth's own tide. The larger early term, the planet's spin being braked by the moon (their eq. 2, 100 times
the Earth's radiogenic heat for a few hundred Myr), is not in this module: it needs the planet's spin and is part of the spin-down work (GEN.842). The flux is held constant over the age, as `tidal_heating` does.

Pure functions of the record's numbers; the random numbers come from the caller. Nothing calls them yet (GEN.1111 creates the rogue rows). Stand-ins are in
`docs/design/research/gap-filling-plan.md` section 6.112.
"""

import math

from planetgen.physics import constants, tidal_heating
from planetgen.util.checks import finite_domain

SURVIVAL_AT_0P01_HILL = 0.82
"""float: The share of ejections whose closest approach stayed beyond 0.01 mutual Hill radii ([R2031], fiducial runs). Source tier T2, tweak level L1 (their number, applied to every mass ratio)."""
SURVIVAL_AT_0P1_HILL = 0.48
"""float: And beyond 0.1 mutual Hill radii (the same). Source tier T2, tweak level L1."""
EJECTED_MOON_ECCENTRICITY = 0.13
"""float: The median eccentricity of a moon carried out with its planet ([R1032] lunar-mass moon at 30 Earth radii, a Jupiter as scatterer). Source tier T3, tweak level L1 (one setup, applied to every moon)."""


class RogueMoonError(ValueError):
    """An input the model cannot use."""


_SLOPE_PER_DECADE = SURVIVAL_AT_0P01_HILL - SURVIVAL_AT_0P1_HILL


def mutual_hill_radius_m(inner_axis_m, outer_axis_m, inner_mass_kg, outer_mass_kg, host_mass_kg):
    """The mutual Hill radius of two planets, m: `(a1 + a2) / 2 * ((m1 + m2) / 3 M)^(1/3)`."""
    if min(inner_axis_m, outer_axis_m, inner_mass_kg, outer_mass_kg, host_mass_kg) <= 0.0:
        raise RogueMoonError("the axes and masses must be positive")
    return 0.5 * (inner_axis_m + outer_axis_m) * ((inner_mass_kg + outer_mass_kg) / (3.0 * host_mass_kg)) ** (1.0 / 3.0)


def closest_approach_survival(r_min_over_hill):
    """The share of ejections whose closest approach stayed beyond `r_min_over_hill` mutual Hill radii: 0.82 at 0.01 and 0.48 at 0.1, linear in the logarithm, held to 0 to 1."""
    if not r_min_over_hill > 0.0:
        raise RogueMoonError("the distance must be positive")
    share = SURVIVAL_AT_0P01_HILL - _SLOPE_PER_DECADE * (math.log10(r_min_over_hill) + 2.0)
    return min(1.0, max(0.0, share))


@finite_domain()
def closest_approach_hill(u):
    """A closest approach drawn from `closest_approach_survival`, in mutual Hill radii: the distance beyond which a share `1 - u` of the ejections stayed. `u` uniform, 0 to 1."""
    if not 0.0 <= u <= 1.0:
        raise RogueMoonError("u must be from 0 to 1")
    survival = 1.0 - u
    return 10.0 ** ((SURVIVAL_AT_0P01_HILL - survival) / _SLOPE_PER_DECADE - 2.0)


def survival_distance_m(r_min_m, planet_mass_kg, scatterer_mass_kg):
    """The widest moon orbit that follows the planet through a closest approach of `r_min_m`, m ([R2031] eq. 4)."""
    if min(r_min_m, planet_mass_kg, scatterer_mass_kg) <= 0.0:
        raise RogueMoonError("the distance and masses must be positive")
    return r_min_m * (planet_mass_kg / (planet_mass_kg + scatterer_mass_kg)) ** (1.0 / 3.0)


def moon_survives(moon_axis_m, r_min_m, planet_mass_kg, scatterer_mass_kg):
    """Whether a moon of orbital radius `moon_axis_m` stays with its planet through a closest approach of `r_min_m`."""
    if not moon_axis_m > 0.0:
        raise RogueMoonError("the moon's orbit must be positive")
    return moon_axis_m < survival_distance_m(r_min_m, planet_mass_kg, scatterer_mass_kg)


def survival_probability(moon_axis_m, planet_mass_kg, scatterer_mass_kg, hill_radius_m):
    """The chance a moon of orbital radius `moon_axis_m` stays, when the closest approach is not recorded: the share of ejections whose `r_min` is wide enough, from `closest_approach_survival`."""
    if not moon_axis_m > 0.0 or not hill_radius_m > 0.0:
        raise RogueMoonError("the moon's orbit and the Hill radius must be positive")
    needed = moon_axis_m / (planet_mass_kg / (planet_mass_kg + scatterer_mass_kg)) ** (1.0 / 3.0)
    return closest_approach_survival(needed / hill_radius_m)


@finite_domain()
def circularisation_flux_w_m2(planet_mass_kg, planet_radius_m, moon_mass_kg, moon_axis_m, eccentricity=EJECTED_MOON_ECCENTRICITY, k2_over_q=None):
    """
    The heat a moon on an eccentric orbit raises in its planet, spread over the planet's surface, W/m2: `10.5 (k2/Q) G m_moon^2 n R^5 e^2 / a^6` ([R1032] eq. 1).

    Args:
        planet_mass_kg, planet_radius_m (float): The planet's.
        moon_mass_kg, moon_axis_m (float): The moon's, and its semi-major axis.
        eccentricity (float): The moon's orbital eccentricity after the ejection.
        k2_over_q (float, optional): The planet's; for a rocky planet `tidal_heating.rigid_k2` over `tidal_heating.ROCK_TIDAL_Q` when not given.
    """
    if min(planet_mass_kg, planet_radius_m, moon_mass_kg, moon_axis_m) <= 0.0 or not 0.0 <= eccentricity < 1.0:
        raise RogueMoonError("the masses, radius and axis must be positive and the eccentricity from 0 to under 1")
    if k2_over_q is None:
        k2_over_q = tidal_heating.rigid_k2(planet_mass_kg, planet_radius_m) / tidal_heating.ROCK_TIDAL_Q
    n = math.sqrt(constants.G * (planet_mass_kg + moon_mass_kg) / moon_axis_m ** 3)
    power = 10.5 * k2_over_q * constants.G * moon_mass_kg ** 2 * n * planet_radius_m ** 5 * eccentricity ** 2 / moon_axis_m ** 6
    return power / (4.0 * math.pi * planet_radius_m ** 2)
