# planetgen/physics/cluster_unbinding.py

"""
Planets a birth cluster strips from their stars (GEN.1112, RP7; GEN.291)
========================================================================

Most stars are born in clusters, and a planet on a wide orbit can be pulled off its star by a passing one. Nothing in the generator did this (`rogue-planet-review.md` row 16): the rogues
were drawn from a rate, and the planets of a cluster star kept every orbit they were given. This module gives the share of planets of a given orbit that a cluster has unbound by an age.

**The rate.** Spurzem, Giersz, Heggie and Lin 2009 (ApJ 697:458, arXiv:astro-ph/0612757) integrate the cross section for an impulsive encounter to change a planet's binding energy by its own
size, which unbinds it, and get (their eqs. 36 and 37, equal-mass stars):

    sigma = (40 pi / 3) G m a / V^2        and        tau = 1 / (n sigma V) = 3 V / (40 pi G rho a)

for the planet's semi-major axis `a`, the stellar mass density `rho = n m` and the typical relative speed `V` of an encounter. The cross section grows with the orbit (the gravitational focusing of a
slow encounter does the work) and the time falls with the density: 22 Myr for a planet at 1,000 AU in a cluster of 100 solar masses per cubic parsec and 2 km/s, 0.2 Myr at ten thousand. The share
still bound after a time is `exp(-t / tau)` for a density held that long; a cluster that expands and thins is a series of such pieces (`survival_fraction`).

**The relative speed** of two stars drawn from a Maxwellian of one-axis dispersion `sigma_1D` is `sqrt(6) sigma_1D` rms (`relative_speed_km_s`).

**The density** a planet sees is the cluster's: the half-mass density `3 M / (8 pi r_h^3)` (`half_mass_density`), from the mass and half-mass radius `open_cluster_evolution` already gives a generated cluster at an age.

**Where it stops holding.** The cross section is the impulsive, non-tidal one: it holds for planets that are soft against the cluster (orbital speed under `V`, `soft`). A harder orbit is less
easily unbound than the formula says, by a factor the paper's N-body runs do not give; the module returns the formula for every orbit and `soft` says which are inside its derivation, so the
number for a hard planet is an upper limit. Planets and rogues of the cluster's own free-floating population, and exchanges that put a planet on another star, are not counted.

Pure functions; nothing calls them yet (GEN.1111 creates the rogue rows). Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.113.
"""

import math

from planetgen.physics import constants, unit_registry as units
from planetgen.util.checks import finite_domain

_msun_to_kg = units.converter("Msun", "kg")

_KM = 1.0e3
_MYR_S = 1.0e6 * constants.SECONDS_PER_YEAR
_G_PC_KM2_S2_PER_MSUN = _msun_to_kg(constants.G) / constants.PARSEC_M / _KM ** 2
"""float: G in pc (km/s)^2 per solar mass, 4.30e-3."""
_AU_PC = constants.AU_M / constants.PARSEC_M
_PC_PER_KM_S_MYR = constants.PARSEC_M / _KM / _MYR_S
"""float: One pc/(km/s) in Myr, 0.978."""


class ClusterUnbindingError(ValueError):
    """An input the model cannot use."""


def relative_speed_km_s(dispersion_1d_km_s):
    """The rms relative speed of two stars of a Maxwellian of one-axis dispersion `dispersion_1d_km_s`, km/s: `sqrt(6)` times it."""
    if not dispersion_1d_km_s > 0.0:
        raise ClusterUnbindingError("the dispersion must be positive")
    return math.sqrt(6.0) * dispersion_1d_km_s


def half_mass_density(cluster_mass_msun, half_mass_radius_pc):
    """The stellar mass density inside the half-mass radius, solar masses per cubic parsec: `3 M / (8 pi r_h^3)`."""
    if not cluster_mass_msun > 0.0 or not half_mass_radius_pc > 0.0:
        raise ClusterUnbindingError("the mass and radius must be positive")
    return 3.0 * cluster_mass_msun / (8.0 * math.pi * half_mass_radius_pc ** 3)


def soft(semi_major_axis_au, host_mass_msun, speed_km_s):
    """Whether the planet is soft against the encounters (its circular speed is under `speed_km_s`): the regime the cross section was derived for."""
    if not semi_major_axis_au > 0.0 or not host_mass_msun > 0.0 or not speed_km_s > 0.0:
        raise ClusterUnbindingError("the axis, mass and speed must be positive")
    return math.sqrt(_G_PC_KM2_S2_PER_MSUN * host_mass_msun / (semi_major_axis_au * _AU_PC)) < speed_km_s


@finite_domain()
def cross_section_au2(semi_major_axis_au, star_mass_msun, speed_km_s):
    """The cross section to unbind a planet, AU^2: `(40 pi / 3) G m a / V^2` (Spurzem et al. 2009 eq. 36)."""
    if not semi_major_axis_au > 0.0 or not star_mass_msun > 0.0 or not speed_km_s > 0.0:
        raise ClusterUnbindingError("the axis, mass and speed must be positive")
    sigma_pc2 = 40.0 * math.pi / 3.0 * _G_PC_KM2_S2_PER_MSUN * star_mass_msun * (semi_major_axis_au * _AU_PC) / speed_km_s ** 2
    return sigma_pc2 / _AU_PC ** 2


@finite_domain()
def unbinding_time_myr(semi_major_axis_au, density_msun_pc3, speed_km_s):
    """The time for a cluster to unbind a planet, Myr: `3 V / (40 pi G rho a)` (Spurzem et al. 2009 eq. 37)."""
    if not semi_major_axis_au > 0.0 or not density_msun_pc3 > 0.0 or not speed_km_s > 0.0:
        raise ClusterUnbindingError("the axis, density and speed must be positive")
    time = 3.0 * speed_km_s / (40.0 * math.pi * _G_PC_KM2_S2_PER_MSUN * density_msun_pc3 * semi_major_axis_au * _AU_PC)
    return time * _PC_PER_KM_S_MYR


@finite_domain()
def survival_fraction(semi_major_axis_au, history):
    """
    The share of planets of an orbit still bound after a cluster history.

    Args:
        semi_major_axis_au (float): The planet's semi-major axis.
        history (sequence): `(duration_myr, density_msun_pc3, speed_km_s)` pieces, in order; the density and speed are held through each. A piece of zero density (the cluster gone) unbinds nothing.
    """
    rate = 0.0
    for duration_myr, density, speed in history:
        if duration_myr < 0.0:
            raise ClusterUnbindingError("a piece of the history cannot last a negative time")
        if density > 0.0 and duration_myr > 0.0:
            rate += duration_myr / unbinding_time_myr(semi_major_axis_au, density, speed)
    return math.exp(-rate)


def unbound_fraction(semi_major_axis_au, history):
    """The share of planets of an orbit the cluster has unbound, `1 - survival_fraction`."""
    return 1.0 - survival_fraction(semi_major_axis_au, history)
