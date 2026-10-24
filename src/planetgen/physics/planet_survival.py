# planetgen/physics/planet_survival.py

"""
Which planet orbits outlive a place's passing stars (GEN.325; docs/design/galaxy/realistic-galactic-core.md, section 7)
===============================================================================================================

A planet of semi-major axis `a` around a star in a field of other systems is stripped when one passes close enough. The
number of such passes in a time `t` is the cluster note's survival rule (docs/design/galaxy/globular-clusters.md, section 7):

    N(a) = n pi a^2 v t (1 + 2 G M / (a v^2))

with `n` the number of systems per volume, `v` the relative speed of the passers-by (`sqrt(6)` times the field's
one-dimensional dispersion `sigma`) and `M` the mass of the host, the first term the geometric cross section and the
second its gravitational focusing. The planet survives with probability `exp(-N)`, so the largest orbit that keeps its
planets is where `N = 1`. `N` is a quadratic in `a`, `c a^2 + b a` with `c = n pi v t` and `b = 2 pi n G M t / v`, so
that orbit is its positive root, in closed form.

`planet_survival_limit_au` takes any density and dispersion, so the same function serves a cluster at a King profile's
`n(r)` and `sigma(r)`, the nuclear disc, and the solar neighbourhood; the caller supplies the place.
"""

import math

from planetgen import tuning
from planetgen.physics import constants

VELOCITY_FACTOR = math.sqrt(6.0)
"""float: The relative speed of a passing system in units of the field's one-dimensional dispersion (the cluster note)."""


def _terms(systems_per_pc3, dispersion_kms, age_gy, host_mass_solar):
    """`(c, b)` of `N(a) = c a^2 + b a` with `a` in metres, from the place's numbers in SI; `(0, 0)` where nothing passes."""
    if systems_per_pc3 <= 0.0 or age_gy <= 0.0 or dispersion_kms <= 0.0:
        return 0.0, 0.0
    density = systems_per_pc3 / constants.PARSEC_M ** 3
    speed = VELOCITY_FACTOR * dispersion_kms * 1000.0
    seconds = age_gy * 1.0e9 * constants.SECONDS_PER_YEAR
    cross = math.pi * density * speed * seconds
    focus = 2.0 * math.pi * density * constants.G * host_mass_solar * constants.SOLAR_MASS_TO_KG * seconds / speed
    return cross, focus


def expected_passes(systems_per_pc3, dispersion_kms, semi_major_axis_au, age_gy=None, host_mass_solar=None):
    """
    The expected number of ionising passes `N(a)` of an orbit over `age_gy` (the survival rule above).

    Args:
        systems_per_pc3 (float): Systems per cubic parsec at the place.
        dispersion_kms (float): The field's one-dimensional velocity dispersion there, km/s.
        semi_major_axis_au (float): The orbit's semi-major axis, AU.
        age_gy (float, optional): The time the planet has been exposed, Gy; `tuning.PLANET_SURVIVAL_AGE_GY` when omitted.
        host_mass_solar (float, optional): The host's mass, solar masses; `tuning.PLANET_SURVIVAL_HOST_MASS_SOLAR` when omitted.

    Returns:
        float: Passes; the survival probability of the orbit is `exp(-passes)`.
    """
    age_gy = tuning.PLANET_SURVIVAL_AGE_GY if age_gy is None else age_gy
    host = tuning.PLANET_SURVIVAL_HOST_MASS_SOLAR if host_mass_solar is None else host_mass_solar
    cross, focus = _terms(systems_per_pc3, dispersion_kms, age_gy, host)
    metres = semi_major_axis_au * constants.AU_M
    return cross * metres ** 2 + focus * metres


def planet_survival_limit_au(systems_per_pc3, dispersion_kms, age_gy=None, host_mass_solar=None):
    """
    The largest semi-major axis whose expected number of ionising passes is one: the orbit beyond which a planet has
    been stripped, as the generator's cut (GEN.325).

    Args:
        systems_per_pc3 (float): Systems per cubic parsec at the place.
        dispersion_kms (float): The field's one-dimensional velocity dispersion there, km/s.
        age_gy (float, optional): The time the planet has been exposed, Gy; `tuning.PLANET_SURVIVAL_AGE_GY` when omitted.
        host_mass_solar (float, optional): The host's mass, solar masses; `tuning.PLANET_SURVIVAL_HOST_MASS_SOLAR` when omitted.

    Returns:
        float: AU; infinite where nothing passes (no density, no dispersion, or no age).
    """
    age_gy = tuning.PLANET_SURVIVAL_AGE_GY if age_gy is None else age_gy
    host = tuning.PLANET_SURVIVAL_HOST_MASS_SOLAR if host_mass_solar is None else host_mass_solar
    cross, focus = _terms(systems_per_pc3, dispersion_kms, age_gy, host)
    if cross <= 0.0:
        return math.inf
    metres = 2.0 / (focus + math.sqrt(focus * focus + 4.0 * cross))
    return metres / constants.AU_M
