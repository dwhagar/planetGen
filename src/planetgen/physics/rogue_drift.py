# planetgen/physics/rogue_drift.py

"""
Where an ejected planet is now: the closed-form drift in the Galaxy (GEN.1114)
===============================================================================

A planet thrown out of its system is on very nearly the star's galactic orbit, a few km/s different, and slowly goes its own way. Over a Gyr that is kiloparsecs along the
orbit, so a rogue planet is nowhere near the star it left (`docs/design/systems/rogue-planet-drift.md`, sections 1 and 3). This module is the pure calculation the virtual-donor
placement (option D) needs: where the rogue is relative to the star it was ejected from, and how fast it moves relative to it, a time after the kick. Nothing calls it yet
(GEN.1111 and GEN.1115 place rogues with it). The kick is an input, so the open questions about its direction and about storage do not hold the module.

**The equations.** In the frame that turns with the star's circular orbit at `Omega = A - B` (Oort's constants), `x` outward, `y` along the rotation and `z` up, a body near the
star obeys Hill's equations ([R0398, R0399], *Galactic Dynamics*, section 3.2; the same as the tidal-stream and epicycle derivations):

    x'' = 2 Omega y' + 4 Omega A x,    y'' = -2 Omega x',    z'' = -nu^2 z,    kappa^2 = -4 B (A - B).

A kick `(u, v, w)` relative to the star at time 0 gives, exactly,

    x = (2 Omega v / kappa^2)(1 - cos kappa t) + (u / kappa) sin kappa t
    y = -(4 Omega A / kappa^2) v t + (4 Omega^2 v / kappa^3) sin kappa t - (2 Omega u / kappa^2)(1 - cos kappa t)
    z = (w / nu) sin nu t

and the velocity relative to the star `x' = (2 Omega v / kappa) sin kappa t + u cos kappa t`, `y' = v - 2 Omega x`, `z' = w cos nu t`. The radial and vertical kicks only make
bounded oscillations (amplitudes `u / kappa`, `w / nu`); the along-orbit kick `v` changes the guiding radius by `2 Omega v / kappa^2` and so the orbital period, and the rogue
slides along the orbit at `4 Omega A / kappa^2` (1.285 for the Sun's neighbourhood) times `v`: 3.9 pc/Myr for 3 km/s.

**What it is good for.** The solution is linear, valid while the separation is a small part of the galactocentric radius (`LINEAR_LIMIT_PC`, 1 to 2 kpc; `Drift.linear` says). Beyond
that the along-orbit offset is an arc: `ring_drift` turns `y` into an azimuth `y / R` wrapped to one turn (the guiding radius shift, an angle, is all that matters when the azimuth
is fully mixed), with `x` and `z` from the closed form, both bounded, so no integrator is needed even at 10 Gyr. A bar, spiral arms or a halo would change the large
excursions; the orbit integrator behind `docs/design/galaxy/galactic-potential.md` replaces this with the same interface when it exists.

**What it does not do.** No secular heating of the floaters by encounters ([R2575], `rogue-planet-drift.md` section 2), no dependence of the Oort constants and the
vertical frequency on galactocentric radius (the caller passes them for the star's position; `Galaxy.from_frequencies` takes the potential's own `Omega`, `kappa` and `nu`),
no star-planet gravity after the kick (the planet is gone). Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.89.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants
from planetgen.util.checks import finite_domain

OORT_A_KM_S_KPC = 15.3
"""float: Oort's constant A, km/s/kpc ([R0497], quoted in `rogue-planet-drift.md` section 3.1). Source tier T1, tweak level L0."""
OORT_B_KM_S_KPC = -11.9
"""float: Oort's constant B, km/s/kpc [R0497]. Source tier T1, tweak level L0."""
VERTICAL_FREQUENCY_KM_S_KPC = 76.0
"""float: The vertical oscillation frequency `nu = sqrt(4 pi G rho)` for a midplane density of 0.1 solar masses per cubic parsec, km/s/kpc. Source tier T6, tweak level L1."""
LINEAR_LIMIT_PC = 1500.0
"""float: The separation up to which the linear equations are taken as good (the design note's "1 to 2 kpc"). Source tier T6, tweak level L1."""

_KM_S_TO_PC_MYR = 1.0e3 * constants.SECONDS_PER_YEAR * 1.0e6 / constants.PARSEC_M
_MYR_TO_KPC_PER_KM_S = _KM_S_TO_PC_MYR / 1.0e3


class RogueDriftError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class Galaxy:
    """
    The epicyclic constants at a star's place in the Galaxy, all in km/s/kpc.

    Attributes:
        omega, kappa, nu (float): The circular, radial and vertical frequencies.
        oort_a (float): `A = Omega - kappa^2 / (4 Omega)`, the shear.
    """

    omega: float
    kappa: float
    nu: float
    oort_a: float

    @classmethod
    def from_oort(cls, a=OORT_A_KM_S_KPC, b=OORT_B_KM_S_KPC, nu=VERTICAL_FREQUENCY_KM_S_KPC):
        """From Oort's constants: `Omega = A - B`, `kappa^2 = -4 B (A - B)`."""
        omega = a - b
        if not omega > 0.0 or not -4.0 * b * omega > 0.0 or not nu > 0.0:
            raise RogueDriftError("Oort's constants must give a rotating, stable orbit (A - B > 0, B < 0) and the vertical frequency must be positive")
        return cls(omega, math.sqrt(-4.0 * b * omega), nu, a)

    @classmethod
    def from_frequencies(cls, omega, kappa, nu):
        """From the potential's own frequencies (`galactic-potential.md` section 7.1): `A = Omega - kappa^2 / (4 Omega)`."""
        if not (omega > 0.0 and kappa > 0.0 and nu > 0.0):
            raise RogueDriftError("the frequencies must be positive")
        return cls(omega, kappa, nu, omega - kappa * kappa / (4.0 * omega))

    @property
    def slide_ratio(self):
        """How fast a rogue slides along the orbit per unit of along-orbit kick, `4 Omega A / kappa^2` (1.285 near the Sun)."""
        return 4.0 * self.omega * self.oort_a / self.kappa ** 2

    @property
    def guiding_radius_shift_per_kick_kpc_per_km_s(self):
        """The change of the guiding radius per unit of along-orbit kick, `2 Omega / kappa^2`, kpc per (km/s)."""
        return 2.0 * self.omega / self.kappa ** 2


SOLAR_NEIGHBOURHOOD = Galaxy.from_oort()
"""Galaxy: The Sun's neighbourhood from the Oort constants above."""


@dataclass(frozen=True)
class Drift:
    """
    The rogue relative to the star it was ejected from, a time after the kick.

    Attributes:
        x_pc, y_pc, z_pc (float): Outward, along the rotation and up, pc.
        vx_km_s, vy_km_s, vz_km_s (float): The velocity relative to the star.
        separation_pc (float): The distance from the star.
        linear (bool): The separation is inside `LINEAR_LIMIT_PC`.
    """

    x_pc: float
    y_pc: float
    z_pc: float
    vx_km_s: float
    vy_km_s: float
    vz_km_s: float
    separation_pc: float
    linear: bool


@dataclass(frozen=True)
class RingDrift:
    """
    The same drift with the along-orbit offset as an angle round the galactic centre.

    Attributes:
        radial_pc, vertical_pc (float): The offsets outward and up (bounded).
        azimuth_rad (float): The angle along the rotation, wrapped to `(-pi, pi]`.
        turns (int): The whole turns the rogue is ahead (positive along the rotation) of the star, before the wrap.
        guiding_radius_shift_pc (float): How far the rogue's guiding radius moved from the star's.
    """

    radial_pc: float
    vertical_pc: float
    azimuth_rad: float
    turns: int
    guiding_radius_shift_pc: float


def _check(time_myr, kick):
    if not time_myr >= 0.0:
        raise RogueDriftError("the time since the ejection must not be negative")
    if len(kick) != 3 or any(not math.isfinite(c) for c in kick):
        raise RogueDriftError("a kick is three finite velocity components (outward, along the rotation, up), km/s")


@finite_domain()
def drift(kick_km_s, time_myr, galaxy=SOLAR_NEIGHBOURHOOD):
    """
    The position and velocity of a rogue relative to its star `time_myr` after a kick.

    Args:
        kick_km_s (sequence): The kick `(u, v, w)` relative to the star, km/s: outward, along the galactic rotation, up.
        time_myr (float): The time since the kick, Myr.
        galaxy (Galaxy): The epicyclic constants.

    Returns:
        Drift: The offset (pc) and relative velocity (km/s).
    """
    _check(time_myr, kick_km_s)
    u, v, w = kick_km_s
    t = time_myr * _MYR_TO_KPC_PER_KM_S                          # kpc per (km/s), the unit of the frequencies' inverse
    om, k, nu, a = galaxy.omega, galaxy.kappa, galaxy.nu, galaxy.oort_a
    kt, nt = k * t, nu * t
    x_kpc = 2.0 * om * v / k ** 2 * (1.0 - math.cos(kt)) + u / k * math.sin(kt)
    y_kpc = (-4.0 * om * a / k ** 2 * v * t + 4.0 * om ** 2 * v / k ** 3 * math.sin(kt) - 2.0 * om * u / k ** 2 * (1.0 - math.cos(kt)))
    z_kpc = w / nu * math.sin(nt)
    vx = 2.0 * om * v / k * math.sin(kt) + u * math.cos(kt)
    vy = v - 2.0 * om * x_kpc
    vz = w * math.cos(nt)
    x, y, z = 1.0e3 * x_kpc, 1.0e3 * y_kpc, 1.0e3 * z_kpc
    separation = math.sqrt(x * x + y * y + z * z)
    return Drift(x, y, z, vx, vy, vz, separation, separation <= LINEAR_LIMIT_PC)


@finite_domain()
def ring_drift(kick_km_s, time_myr, galactocentric_radius_pc, galaxy=SOLAR_NEIGHBOURHOOD):
    """
    The drift with the along-orbit offset as an angle `y / R` round the galactic centre, wrapped to one turn: good for a rogue that has gone far along the orbit, where the
    linear offset in pc is an arc.

    Args:
        kick_km_s, time_myr, galaxy: As for `drift`.
        galactocentric_radius_pc (float): The star's distance from the galactic centre.
    """
    if not galactocentric_radius_pc > 0.0:
        raise RogueDriftError("the galactocentric radius must be positive")
    here = drift(kick_km_s, time_myr, galaxy)
    angle = here.y_pc / galactocentric_radius_pc
    turns = math.floor((angle + math.pi) / (2.0 * math.pi))
    wrapped = angle - 2.0 * math.pi * turns
    shift = 1.0e3 * kick_km_s[1] * galaxy.guiding_radius_shift_per_kick_kpc_per_km_s
    return RingDrift(here.x_pc, here.z_pc, wrapped, int(turns), shift)


def donor_offset_pc(kick_km_s, time_myr, galaxy=SOLAR_NEIGHBOURHOOD):
    """
    Where the donor star is, relative to the rogue, now: minus the rogue's offset from it. The closed form does not depend on where the star is, so a sector that holds a
    rogue puts the donor at its position plus this (option D, `rogue-planet-drift.md` section 5).
    """
    here = drift(kick_km_s, time_myr, galaxy)
    return (-here.x_pc, -here.y_pc, -here.z_pc)


def along_orbit_speed_pc_myr(v_km_s, galaxy=SOLAR_NEIGHBOURHOOD):
    """The speed the rogue slides along the orbit at for an along-orbit kick `v_km_s`, pc/Myr (negative: it falls behind for a forward kick)."""
    return -galaxy.slide_ratio * v_km_s * _KM_S_TO_PC_MYR


def kick_from_direction(speed_km_s, direction):
    """The kick `(u, v, w)` of a speed in a direction (any non-zero vector `(outward, along the rotation, up)`)."""
    if not speed_km_s >= 0.0:
        raise RogueDriftError("the speed must not be negative")
    norm = math.sqrt(sum(c * c for c in direction))
    if len(direction) != 3 or not norm > 0.0:
        raise RogueDriftError("a direction is a non-zero vector of three components")
    return tuple(speed_km_s * c / norm for c in direction)
