# planetgen/physics/rogue_kick.py

"""
How fast an ejected planet leaves, and how fast it moves through the Galaxy (GEN.1112 RP6, GEN.1113)
====================================================================================================

The rogue planets had the speed of a star on a generic galactic orbit and no ejection kick (`rogue-planet-review.md` row 14). A planet thrown out of its system by a
scattering planet leaves the star with a speed at infinity that depends on the scatterer, not on the star or on the planet it throws. This module gives that speed from the
scatterer's mass and orbit, the survival of the ejected planet's moons, and the velocity the rogue has now relative to the local standard of rest, using the drift
(`rogue_drift`) for the years between. Nothing calls it yet (GEN.1111 creates the rogue rows, and GEN.1119 gives each ejection its scatterer and orbit). The reading behind
every number is in `docs/design/research/rogue-ejection-kick.md`.

**The speed** (GEN.1113 replaced the earlier energy-step guess). [R2031] ([R2031: Eq. 19], 24 and 26) find from N-body runs of scattering planets that
the speeds at infinity peak near

    v_c = 3 km/s (m_p / M_Jup)^(1/2) (a / 1 AU)^(-1/2) m_p / (m_p + m_ej)

for a scatterer of mass `m_p` on an orbit `a` ejecting a body of mass `m_ej`, independent of the star's mass for stars 100 to 1000 times the scatterer's. Kramer and Profumo
2026 ([R2470: Sec. 4.1] to 4.2) find a mean of 3.4 to 3.6 km/s for a 10 Jupiter-mass scatterer at 5.2 AU with no trend over four decades in the ejected mass, and 1.15
km/s for 1 Jupiter mass at 5.2 AU: the mean is 0.87 of `v_c` at both. The speeds spread from zero to about 2.5 `v_c` (Huang and Lai Sec. 4.3), so `v_inf = v_c X` with `X`
log-normal of mean `MEAN_OVER_SCALE` and 99th percentile 2.5 (`SPEED_SIGMA_LN`), drawn by the inverse of the normal distribution from a uniform number `u`. [R3937] assume a kick falling with the ejected mass as `M^-1/2` instead (a stated prescription, not a result), and the two N-body papers contradict it;
their medians of 0.7 to 10.7 km/s bracket the speeds of giants scattering between 1 and 30 AU here. The tail that eccentric giants and three-planet systems add (Kramer and
Profumo: up to 7 km/s mean, 80 km/s at the first outbound step for a giant of eccentricity 0.9; Huang and Lai: 5 to 8 `v_c` for three planets) is **not** modelled yet.
A scatterer whose Safronov number (`ejects`) is under 1 only stirs.

**Moons.** [R2031] (Sec. 1 and Eq. 4) judge a moon by how close the planets came: it stays bound if its orbital radius is under `MOON_LOSS_FRACTION` of the Hill
sphere that the closest approach `r_min` leaves it, `r_moon < f r_min (m_ej / m_12)^(1/3)`. Their runs give `r_min` in units of the mutual Hill radius with a median of
0.089 for a 10 Jupiter-mass scatterer ejecting a 1 Jupiter-mass planet, falling as the scatterer's mass to the 1.13 (0.034 at 5, 0.0066 at 1 Jupiter masses at a fixed ratio
of 10), and a spread that puts 48 percent of runs above 0.1 and 82 percent above 0.01 (a log-normal of ln-width `CLOSEST_APPROACH_SIGMA_LN`). `closest_approach_m` draws it and `moon_survival_probability` is
the share of ejections that keep a moon at a given radius.

**The direction** is isotropic in the star's frame (`isotropic_direction`): Boss answered (GEN.1113, 5:29 PM EDT) that it is isotropic, added to the planet's actual motion and the
system's. Neither N-body paper reports a direction statistic.

**The velocity through the Galaxy.** The rogue's velocity is its host star's plus the drift's relative velocity (`rogue_drift.drift`): the speed at infinity is the speed relative
to the star at ejection, and the Galaxy's tide turns it over the epicycle after. `rogue_velocity` adds the host's velocity relative to the local standard of rest (the host's
own, from the population's velocity dispersion, drawn by the caller) to that relative velocity.

Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.87.
"""

import math
from statistics import NormalDist

from planetgen.physics import constants, rogue_drift, unit_registry as units
from planetgen.util.checks import finite_domain

_kg_to_mjup = units.converter("kg", "Mjup")

SCATTER_SPEED_KM_S = 3.0
"""float: The scale of the ejection speed for a scatterer of one Jupiter mass at 1 AU, km/s: `v_c = SCATTER_SPEED_KM_S (m_p / M_Jup)^(1/2) (a / 1 AU)^(-1/2) m_p / (m_p + m_ej)`
([R2031: Eq. 26], read first-hand). Source tier T1, tweak level L0."""
MEAN_OVER_SCALE = 0.87
"""float: The mean ejection speed in units of `v_c`: [R2470] give 1.15 km/s for 1 Jupiter mass and 3.4 to 3.6 for 10 at 5.2 AU, where `v_c` is 1.32 and 4.16 (0.87 and 0.84). Source tier T1, tweak level L1."""
SPEED_SIGMA_LN = 0.51
"""float: The ln-width of the ejection speed in units of `v_c`, set so the 99th percentile is 2.5 `v_c`, the top of the range of [R2031] (Sec. 4.3, a figure reading, not a table). Source tier T6, tweak level L2."""
MOON_LOSS_FRACTION = 0.5
"""float: `f` of Huang and Lai's Eq. 4: a moon stays bound if its orbital radius is under `f` times the Hill sphere that the closest approach leaves it ("a high survival probability requires f ~ 1/2"). Source tier T1, tweak level L1."""
CLOSEST_APPROACH_MEDIAN_HILL = 0.089
"""float: The median closest approach of the ejected planet to the scatterer before ejection, in mutual Hill radii, for a scatterer of `CLOSEST_APPROACH_REFERENCE_JUPITERS` ejecting a tenth of its mass (Huang and Lai Fig. 2). Source tier T1, tweak level L2."""
CLOSEST_APPROACH_REFERENCE_JUPITERS = 10.0
"""float: The scatterer's mass, Jupiter masses, at which the median closest approach is `CLOSEST_APPROACH_MEDIAN_HILL`."""
CLOSEST_APPROACH_MASS_INDEX = 1.13
"""float: The median closest approach goes as the scatterer's mass to this power (0.089, 0.034, 0.0066 at 10, 5, 1 Jupiter masses at a fixed ratio of 10; Huang and Lai Fig. 2). Source tier T3, tweak level L3."""
CLOSEST_APPROACH_SIGMA_LN = 2.35
"""float: The ln-width of the closest approach: 48 percent of runs above 0.1 and 82 percent above 0.01 mutual Hill radii at a median of 0.089 give 2.3 to 2.4 (Huang and Lai Sec. 7). Source tier T1, tweak level L1."""
EJECTION_SAFRONOV_NUMBER = 1.0
"""float: The Safronov number above which a scatterer ejects what it scatters and below which it only stirs ([R1417, R3867] use the same cut). Source tier T6, tweak level L0."""
_KM = 1.0e3


class RogueKickError(ValueError):
    """An input the model cannot use."""


def escape_speed_km_s(host_mass_kg, distance_au):
    """The escape speed from a star at a distance, km/s."""
    if not host_mass_kg > 0.0 or not distance_au > 0.0:
        raise RogueKickError("the mass and the distance must be positive")
    return math.sqrt(2.0 * constants.G * host_mass_kg / (distance_au * constants.AU_M)) / _KM


def circular_speed_km_s(host_mass_kg, distance_au):
    """The circular orbital speed round a star at a distance, km/s."""
    return escape_speed_km_s(host_mass_kg, distance_au) / math.sqrt(2.0)


def planet_escape_speed_km_s(planet_mass_kg, planet_radius_m):
    """The escape speed from the surface of the scattering planet, km/s."""
    if not planet_mass_kg > 0.0 or not planet_radius_m > 0.0:
        raise RogueKickError("the planet's mass and radius must be positive")
    return math.sqrt(2.0 * constants.G * planet_mass_kg / planet_radius_m) / _KM


def safronov_number(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au):
    """`Theta = v_esc,p^2 / (2 v_orb^2) = (M_p / M_*)(a / R_p)`: above 1 a scatterer ejects, below 1 it stirs."""
    return (planet_escape_speed_km_s(planet_mass_kg, planet_radius_m) / circular_speed_km_s(host_mass_kg, distance_au)) ** 2 / 2.0


def ejects(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au):
    """Whether the scatterer ejects what it scatters (its Safronov number is above `EJECTION_SAFRONOV_NUMBER`)."""
    return safronov_number(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au) > EJECTION_SAFRONOV_NUMBER


def scatter_speed_km_s(planet_mass_kg, distance_au, ejected_mass_kg=0.0):
    """`v_c`, the scale of the speeds a scatterer gives what it ejects, km/s: `SCATTER_SPEED_KM_S (m_p / M_Jup)^(1/2) (a / 1 AU)^(-1/2) m_p / (m_p + m_ej)`."""
    if not planet_mass_kg > 0.0 or not distance_au > 0.0 or ejected_mass_kg < 0.0:
        raise RogueKickError("the scatterer's mass and distance must be positive and the ejected mass not negative")
    return (SCATTER_SPEED_KM_S * math.sqrt(_kg_to_mjup(planet_mass_kg) / distance_au)
            * planet_mass_kg / (planet_mass_kg + ejected_mass_kg))


@finite_domain()
def ejection_speed_km_s(u, planet_mass_kg, distance_au, ejected_mass_kg=0.0):
    """
    The speed at infinity relative to the star of a body ejected by a scatterer, km/s.

    Args:
        u (float): A uniform random number, 0 to 1 (excluding 1).
        planet_mass_kg (float): The scatterer's mass.
        distance_au (float): Where the scattering happens: the scatterer's distance from the star.
        ejected_mass_kg (float): The mass of the body thrown out, which lowers the speed a little when it is a sizeable part of the scatterer's.
    """
    if not 0.0 <= u < 1.0:
        raise RogueKickError("u must be at least 0 and under 1")
    scale = scatter_speed_km_s(planet_mass_kg, distance_au, ejected_mass_kg)
    z = NormalDist().inv_cdf(max(u, 1.0e-12))
    mu = math.log(MEAN_OVER_SCALE) - 0.5 * SPEED_SIGMA_LN ** 2
    return scale * math.exp(mu + SPEED_SIGMA_LN * z)


def mutual_hill_radius_m(host_mass_kg, distance_au, planet_mass_kg, ejected_mass_kg):
    """The mutual Hill radius of the scatterer and the ejected planet at the scatterer's distance, metres."""
    if not host_mass_kg > 0.0 or not distance_au > 0.0 or planet_mass_kg <= 0.0 or ejected_mass_kg <= 0.0:
        raise RogueKickError("the masses and the distance must be positive")
    return distance_au * constants.AU_M * ((planet_mass_kg + ejected_mass_kg) / (3.0 * host_mass_kg)) ** (1.0 / 3.0)


def _closest_approach_median_hill(planet_mass_kg):
    return CLOSEST_APPROACH_MEDIAN_HILL * (_kg_to_mjup(planet_mass_kg) / CLOSEST_APPROACH_REFERENCE_JUPITERS) ** CLOSEST_APPROACH_MASS_INDEX


@finite_domain()
def closest_approach_m(u, host_mass_kg, distance_au, planet_mass_kg, ejected_mass_kg):
    """The closest approach of the ejected planet to the scatterer before ejection, metres, from a uniform number `u` (0 to 1, excluding both ends are clamped)."""
    if not 0.0 <= u < 1.0:
        raise RogueKickError("u must be at least 0 and under 1")
    z = NormalDist().inv_cdf(max(u, 1.0e-12))
    hill = mutual_hill_radius_m(host_mass_kg, distance_au, planet_mass_kg, ejected_mass_kg)
    return hill * _closest_approach_median_hill(planet_mass_kg) * math.exp(CLOSEST_APPROACH_SIGMA_LN * z)


def moon_survives(moon_radius_m, closest_approach, planet_mass_kg, ejected_mass_kg):
    """Whether a moon at `moon_radius_m` round the ejected planet stays bound through a closest approach (Huang and Lai Eq. 4)."""
    if not moon_radius_m > 0.0 or not closest_approach > 0.0 or planet_mass_kg <= 0.0 or ejected_mass_kg <= 0.0:
        raise RogueKickError("the radii and masses must be positive")
    return moon_radius_m < MOON_LOSS_FRACTION * closest_approach * (ejected_mass_kg / (ejected_mass_kg + planet_mass_kg)) ** (1.0 / 3.0)


def moon_survival_probability(moon_radius_m, host_mass_kg, distance_au, planet_mass_kg, ejected_mass_kg):
    """The share of ejections that keep a moon at `moon_radius_m`: the chance the closest approach is wide enough, from the log-normal of `closest_approach_m`."""
    if not moon_radius_m > 0.0:
        raise RogueKickError("the moon's radius must be positive")
    hill = mutual_hill_radius_m(host_mass_kg, distance_au, planet_mass_kg, ejected_mass_kg)
    needed = moon_radius_m / (MOON_LOSS_FRACTION * (ejected_mass_kg / (ejected_mass_kg + planet_mass_kg)) ** (1.0 / 3.0))
    z = math.log(needed / (hill * _closest_approach_median_hill(planet_mass_kg))) / CLOSEST_APPROACH_SIGMA_LN
    return 1.0 - NormalDist().cdf(z)


def isotropic_direction(u_cos, u_phi):
    """A direction uniform on the sphere `(outward, along the rotation, up)` from two uniform numbers, 0 to 1."""
    if not (0.0 <= u_cos <= 1.0 and 0.0 <= u_phi <= 1.0):
        raise RogueKickError("the random numbers must be from 0 to 1")
    z = 2.0 * u_cos - 1.0
    r = math.sqrt(max(0.0, 1.0 - z * z))
    phi = 2.0 * math.pi * u_phi
    return (r * math.cos(phi), r * math.sin(phi), z)


def kick(speed_km_s, direction):
    """The kick `(u, v, w)` of a speed in a direction, in `rogue_drift`'s frame."""
    return rogue_drift.kick_from_direction(speed_km_s, direction)


def rogue_velocity(host_velocity_km_s, kick_km_s, time_myr, galaxy=rogue_drift.SOLAR_NEIGHBOURHOOD):
    """
    The rogue's velocity relative to the local standard of rest, `time_myr` after the kick: the host star's velocity plus the velocity relative to the star that the drift gives.

    Args:
        host_velocity_km_s (sequence): The host star's `(U, V, W)` relative to the local standard of rest now, km/s: outward, along the rotation, up.
        kick_km_s (sequence): The kick at ejection.
        time_myr (float): The time since ejection.
    """
    if len(host_velocity_km_s) != 3 or any(not math.isfinite(c) for c in host_velocity_km_s):
        raise RogueKickError("a host velocity is three finite components")
    here = rogue_drift.drift(kick_km_s, time_myr, galaxy)
    return (host_velocity_km_s[0] + here.vx_km_s, host_velocity_km_s[1] + here.vy_km_s, host_velocity_km_s[2] + here.vz_km_s)


def relative_speed_km_s(velocity, reference=(0.0, 0.0, 0.0)):
    """The speed of a body with a velocity relative to an observer with the reference velocity (both relative to the same standard of rest), km/s."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(velocity, reference)))
