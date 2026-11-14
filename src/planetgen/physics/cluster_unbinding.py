# planetgen/physics/cluster_unbinding.py

"""
Planets a birth cluster strips from their stars (GEN.1112, RP7; GEN.291)
========================================================================

Most stars are born in clusters, and a planet on a wide orbit can be pulled off its star by a passing one. Nothing in the generator did this (`rogue-planet-review.md` row 16): the rogues
were drawn from a rate, and the planets of a cluster star kept every orbit they were given. This module gives the share of planets of a given orbit that a cluster has unbound by an age.

**The rate.** [R4184] (ApJ 697:458, [R4184]) integrate the cross section for an impulsive encounter to change a planet's binding energy by its own
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

**The generator reads it** (GEN.1288). A star that belongs to a generated or catalogue cluster (`stars.cluster_id`) is built with the cluster's history (`history_of`): from the end of its disk (`DISK_LIFETIME_MYR`) to the cluster's age, in
`HISTORY_PIECES` equal pieces, each with the density and the relative speed of the cluster then. The cluster's stored mass and half-mass radius are the ones it has now; its radius grew from `birth radius` by the mid-point of
`tuning.OPEN_CLUSTER_EXPANSION_FACTOR_RANGE` over `tuning.OPEN_CLUSTER_EXPANSION_AGE_GY` (`open_cluster_evolution.half_mass_radii_pc`), so the birth radius is worked back from it, and the mass is held at the present value (the
cluster lost some, so the early density is a lower bound). The speed is the virial one of a Plummer sphere: `sigma_1D^2 = G M / (6 r_vir)` with `r_vir = VIRIAL_OVER_HALF_MASS r_h`. `unbinding_epoch_myr` takes a planet's uniform draw
and says when, if at all, the cluster pulled it off: the draw is turned into a cumulative hazard and walked through the pieces, counting only the pieces where the planet is soft. The planet leaves with the speed that just unbinds it
(`unbinding_kick_km_s`).

Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.113.
"""

import math

from planetgen import tuning
from planetgen.physics import constants, unit_registry as units
from planetgen.util.checks import finite_domain

_msun_to_kg = units.converter("Msun", "kg")

_KM = 1.0e3
_MYR_S = constants.SECONDS_PER_MYR
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
    """The cross section to unbind a planet, AU^2: `(40 pi / 3) G m a / V^2` ([R4184] eq. 36)."""
    if not semi_major_axis_au > 0.0 or not star_mass_msun > 0.0 or not speed_km_s > 0.0:
        raise ClusterUnbindingError("the axis, mass and speed must be positive")
    sigma_pc2 = 40.0 * math.pi / 3.0 * _G_PC_KM2_S2_PER_MSUN * star_mass_msun * (semi_major_axis_au * _AU_PC) / speed_km_s ** 2
    return sigma_pc2 / _AU_PC ** 2


@finite_domain()
def unbinding_time_myr(semi_major_axis_au, density_msun_pc3, speed_km_s):
    """The time for a cluster to unbind a planet, Myr: `3 V / (40 pi G rho a)` ([R4184] eq. 37)."""
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


DISK_LIFETIME_MYR = 3.0
"""float: When the planets exist: the end of the disk (`scatter.DISK_LIFETIME_YR`, 3 Myr); the cluster's encounters before then meet a disk, not planets."""

HISTORY_PIECES = 8
"""int: How many equal pieces the cluster's life is cut into, each held at its mid-point's density and speed (stand-in CU4)."""

VIRIAL_OVER_HALF_MASS = 1.30
"""float: A Plummer sphere's virial radius over its half-mass radius, `16 / (3 pi) / 1.305` (stand-in CU5)."""

UNBINDING_KICK_OVER_CIRCULAR = math.sqrt(2.0) - 1.0
"""float: The speed a planet unbound by an encounter leaves its star with, over its circular speed there: the impulse that just lifts it to the escape speed `sqrt(2)` times the circular one (stand-in CU6)."""


def virial_dispersion_km_s(cluster_mass_msun, half_mass_radius_pc):
    """The one-axis dispersion of a Plummer cluster in virial equilibrium, km/s: `sqrt(G M / (6 r_vir))`, `r_vir = VIRIAL_OVER_HALF_MASS r_h`."""
    if not cluster_mass_msun > 0.0 or not half_mass_radius_pc > 0.0:
        raise ClusterUnbindingError("the mass and radius must be positive")
    return math.sqrt(_G_PC_KM2_S2_PER_MSUN * cluster_mass_msun / (6.0 * VIRIAL_OVER_HALF_MASS * half_mass_radius_pc))


def history_of(cluster_mass_msun, half_mass_radius_pc, age_gy):
    """
    The history of a cluster as `survival_fraction` takes it, for a planet that formed when the disk ended: `HISTORY_PIECES` equal pieces from `DISK_LIFETIME_MYR` to `age_gy`, each `(duration_myr, half-mass density, relative speed)`.

    Args:
        cluster_mass_msun (float): The mass the cluster has now.
        half_mass_radius_pc (float): Its half-mass radius now.
        age_gy (float): Its age.

    Returns:
        tuple: The pieces; empty for a cluster no older than the disk.
    """
    if not cluster_mass_msun > 0.0 or not half_mass_radius_pc > 0.0:
        raise ClusterUnbindingError("the mass and radius must be positive")
    start, end = DISK_LIFETIME_MYR, age_gy * 1000.0
    if end <= start:
        return ()
    low, high = tuning.OPEN_CLUSTER_EXPANSION_FACTOR_RANGE
    factor = 0.5 * (low + high)
    span_myr = tuning.OPEN_CLUSTER_EXPANSION_AGE_GY * 1000.0
    born_pc = half_mass_radius_pc / (1.0 + (factor - 1.0) * min(end / span_myr, 1.0))
    width = (end - start) / HISTORY_PIECES
    pieces = []
    for k in range(HISTORY_PIECES):
        middle = start + (k + 0.5) * width
        radius = born_pc * (1.0 + (factor - 1.0) * min(middle / span_myr, 1.0))
        speed = relative_speed_km_s(virial_dispersion_km_s(cluster_mass_msun, radius))
        pieces.append((width, half_mass_density(cluster_mass_msun, radius), speed))
    return tuple(pieces)


def unbinding_epoch_myr(semi_major_axis_au, host_mass_msun, history, u):
    """
    When the cluster pulls a planet off its star, Myr after the cluster formed, or `None` if it keeps it.

    The planet's uniform draw `u` is its place in the order of fates: the cumulative hazard it takes to unbind it is `-ln(1 - u)`, and the history is walked piece by piece (hazard `duration / unbinding_time`) until that much has
    been taken. A piece in which the planet is not soft (`soft`: the formula's own validity) adds none. Over many planets the share unbound is `unbound_fraction` over the soft pieces, and the times follow the density of the history.
    """
    if not 0.0 <= u < 1.0:
        raise ClusterUnbindingError("the draw must be in [0, 1)")
    need = -math.log1p(-u)
    now = DISK_LIFETIME_MYR
    for duration_myr, density, speed in history:
        if density > 0.0 and duration_myr > 0.0 and soft(semi_major_axis_au, host_mass_msun, speed):
            rate = 1.0 / unbinding_time_myr(semi_major_axis_au, density, speed)
            if need <= rate * duration_myr:
                return now + need / rate
            need -= rate * duration_myr
        now += duration_myr
    return None


def unbinding_kick_km_s(circular_speed_km_s):
    """The speed at which a planet unbound by the cluster leaves its star, km/s: `UNBINDING_KICK_OVER_CIRCULAR` times its circular speed."""
    return UNBINDING_KICK_OVER_CIRCULAR * circular_speed_km_s
