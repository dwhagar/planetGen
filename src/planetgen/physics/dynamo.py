"""
The field of a planetary dynamo from the power that drives it and the spin (GEN.966, the lock-free half of GEN.736).

A pure module: a core's convecting power, its shell and the spin of the body in, a dipole moment and a field out. It reads nothing from a body
and draws nothing. The power comes from the core's own thermal history (the heat across its boundary, `interior_heat`/`core_heat`, GEN.699 and
GEN.838) or, for a giant, from its luminosity (`giant_power_w`); the spin is the body's (`spin`, GEN.731). The values stored on a body are set by
GEN.967, which replaces the design defaults of `magnetism` with these.

1. **Strength.** The buoyancy-flux scaling of [R3355], which [R0784] extend: the
   dipole's equatorial field at the top of a dipolar dynamo's region does not depend on the spin or the conductivity, only on the buoyancy flux
   `F = P / (rho V)` (W/kg, the power per mass of the region), its thickness `D` and the density `rho`:
   `B = gamma sqrt(rho mu0 / 2) (F D)^(1/3)` with `gamma` about 0.2, and the moment is `4 pi r_o^3 B / mu0` (the design note section 3.1). The same
   as `B = c sqrt(mu0) rho^(1/6) (P D / V)^(1/3)`.
2. **Dipole share.** A dynamo turning slowly for its power is multipolar: the local Rossby number
   `Ro_l = 0.09 (Ra* / Ra*_Earth)^0.41`, with `Ra* = F / (Omega^3 D^2)` the modified flux Rayleigh number of the region, is the dipole's gate,
   and the share of the dipole is `f_dip = 0.15 + 0.85 / (1 + (Ro_l / 0.12)^4)`, the rule GEN.86 used before this module (`magnetism` now stores the moment with it applied). The scaling above
   holds for Earth and every dynamo that is dipolar as Earth is; a slower one has its moment times `f_dip` over Earth's. With the power held at its
   Earth value `Ro_l` is `0.09 (P_rot / 24 h)^1.23`: the exponent `1.23` is `3 * 0.41`.
3. **Lifetime.** The dynamo runs while the power is above zero (`lifetime`): the core's own history says when it started and when it stopped.
4. **Giants.** The metallic region of a giant is driven by a share of the heat the planet radiates (`GIANT_POWER_SHARE`, stand-in DY2).

Every number that is not a published value is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.30).
"""

import math
from dataclasses import dataclass

VACUUM_PERMEABILITY = 4e-7 * math.pi
"""float: mu0, T m/A."""

DIPOLE_COEFFICIENT = 0.2
"""float: `gamma`, the coefficient of the dipole's field in the buoyancy-flux scaling ([R3355]: 145 numerical dynamos, 0.2, which
puts the solar-system moments within about 15 percent; second-hand)."""
EARTH_OUTER_CORE_M = 3.480e6
EARTH_INNER_CORE_M = 1.221e6
EARTH_OUTER_CORE_DENSITY_KG_M3 = 1.1e4
"""float: PREM's outer core: its outer radius, the inner core's radius and the outer core's mean density."""
EARTH_SPIN_RAD_S = 7.2921e-5
"""float: Earth's spin rate, rad/s (the sidereal day)."""
EARTH_FLUX_DEPTH_M2_S2 = 1.4e-6
"""float: Earth's `F D` for the scaling, m2/s2, the value that gives its 7.8e22 A m2 (design note section 3.1, a check): `F` is about 6.2e-13 W/kg, which
is `EARTH_DYNAMO_POWER_W` of 1.15 TW in the outer core. STAND-IN DY3: the published range of the power is 0.2 to 3 TW."""
EARTH_ROSSBY = 0.09
"""float: The local Rossby number of Earth's dynamo [R3355], the unit of `local_rossby`."""
DIPOLE_ROSSBY = 0.12
"""float: The local Rossby number above which a dynamo is multipolar ([R0784]: 0.12; the band is 0.1 to 0.12)."""
FLUX_RAYLEIGH_EXPONENT = 0.41
"""float: The exponent of the flux Rayleigh number in the Rossby number ([R0784]: 0.41)."""
MIN_DIPOLE_FRACTION = 0.15
"""float: STAND-IN DY4: the share of the field a multipolar dynamo keeps in its dipole (the rule of `magnetism` since GEN.86; the simulations give 0.05 to 0.2)."""

GIANT_POWER_SHARE = 0.1
"""float: STAND-IN DY2: the share of a giant's radiated heat that drives its dynamo. With Jupiter's 3.35e17 W it gives 1.6e27 A m2, the 1.5e27 observed; the
same share gives Saturn six times its 4.6e25 (a stable helium-rain layer weakens Saturn's, which this does not model)."""


def region_volume_m3(outer_m, inner_m=0.0):
    """The volume of the shell of the dynamo region, m3."""
    return 4.0 / 3.0 * math.pi * (outer_m ** 3 - inner_m ** 3)


def buoyancy_flux(power_w, outer_m, inner_m, density_kg_m3):
    """`F = P / (rho V)`, the power per mass of the region, W/kg (m2/s3)."""
    return power_w / (density_kg_m3 * region_volume_m3(outer_m, inner_m))


def flux_rayleigh(power_w, outer_m, inner_m, density_kg_m3, spin_rad_s):
    """The modified flux Rayleigh number `Ra* = F / (Omega^3 D^2)` of the region, dimensionless (0 for no power)."""
    thickness = outer_m - inner_m
    if power_w <= 0.0 or thickness <= 0.0 or spin_rad_s <= 0.0 or density_kg_m3 <= 0.0:
        return 0.0
    return buoyancy_flux(power_w, outer_m, inner_m, density_kg_m3) / (spin_rad_s ** 3 * thickness ** 2)


EARTH_BUOYANCY_FLUX_M2_S3 = EARTH_FLUX_DEPTH_M2_S2 / (EARTH_OUTER_CORE_M - EARTH_INNER_CORE_M)
EARTH_DYNAMO_POWER_W = (EARTH_BUOYANCY_FLUX_M2_S3 * EARTH_OUTER_CORE_DENSITY_KG_M3 * region_volume_m3(EARTH_OUTER_CORE_M, EARTH_INNER_CORE_M))
"""float: The power Earth's dynamo is driven with in this scaling, W: about 1.15 TW."""
_EARTH_RAYLEIGH = EARTH_BUOYANCY_FLUX_M2_S3 / (EARTH_SPIN_RAD_S ** 3 * (EARTH_OUTER_CORE_M - EARTH_INNER_CORE_M) ** 2)


def local_rossby(power_w, outer_m, inner_m, density_kg_m3, spin_rad_s):
    """The local Rossby number of the dynamo region, `0.09 (Ra* / Ra*_Earth)^0.41` (infinite for no spin, 0 for no power)."""
    if spin_rad_s <= 0.0:
        return math.inf
    return EARTH_ROSSBY * (flux_rayleigh(power_w, outer_m, inner_m, density_kg_m3, spin_rad_s) / _EARTH_RAYLEIGH) ** FLUX_RAYLEIGH_EXPONENT


def dipole_fraction(rossby):
    """`f_dip`: the share of the field at the top of the region that is dipole: 1 for a dynamo that turns fast for its power (a small Rossby
    number), falling through `DIPOLE_ROSSBY` to `MIN_DIPOLE_FRACTION` for one that is multipolar."""
    if rossby <= 0.0:
        return 1.0
    return MIN_DIPOLE_FRACTION + (1.0 - MIN_DIPOLE_FRACTION) / (1.0 + (rossby / DIPOLE_ROSSBY) ** 4)


def dipolar_field_t(power_w, outer_m, inner_m, density_kg_m3):
    """The dipole's equatorial field at the top of a dipolar dynamo's region, `gamma sqrt(rho mu0 / 2) (F D)^(1/3)`, T."""
    flux_depth = buoyancy_flux(power_w, outer_m, inner_m, density_kg_m3) * (outer_m - inner_m)
    return DIPOLE_COEFFICIENT * math.sqrt(density_kg_m3 * VACUUM_PERMEABILITY / 2.0) * flux_depth ** (1.0 / 3.0)


@dataclass(frozen=True)
class Field:
    """
    A dynamo's field.

    Attributes:
        moment_a_m2 (float): The dipole moment, A m2 (0 while the dynamo does not run).
        rossby (float): The local Rossby number of the region.
        dipole_fraction (float): `f_dip`, the dipole's share of the field at the top of the region (the rule of `dipole_fraction`).
        dipolar_field_t (float): The scaling's dipole field at the top of the region for a dipolar dynamo, T.
        outer_m (float): The region's outer radius, m.
    """

    moment_a_m2: float
    rossby: float
    dipole_fraction: float
    dipolar_field_t: float
    outer_m: float

    @property
    def top_equatorial_t(self):
        """The dipole's equatorial field at the top of the region, T."""
        return VACUUM_PERMEABILITY * self.moment_a_m2 / (4.0 * math.pi * self.outer_m ** 3)

    def surface_equatorial_t(self, radius_m):
        """The dipole's equatorial field at the surface of a body of `radius_m`, T."""
        return VACUUM_PERMEABILITY * self.moment_a_m2 / (4.0 * math.pi * radius_m ** 3)


def field(power_w, outer_m, inner_m, density_kg_m3, spin_rad_s, attenuation=1.0):
    """
    The `Field` of a dynamo region.

    Args:
        power_w: The power that drives the convection of the region, W (0 or less: the dynamo does not run).
        outer_m: The region's outer radius (the core's, or a giant's metallic layer's), m.
        inner_m: Its inner radius (the inner core's, 0 for none), m.
        density_kg_m3: Its mean density.
        spin_rad_s: The body's spin rate, rad/s.
        attenuation: The share of the field that gets through a layer on top of the region that does not convect (Mercury's), `1` for none.
            Nothing here says when a layer is stable: the caller does.
    """
    thickness = outer_m - inner_m
    if power_w <= 0.0 or thickness <= 0.0 or outer_m <= 0.0 or density_kg_m3 <= 0.0:
        return Field(0.0, 0.0 if power_w <= 0.0 else math.inf, 1.0, 0.0, max(outer_m, 0.0))
    rossby = local_rossby(power_w, outer_m, inner_m, density_kg_m3, spin_rad_s)
    share = dipole_fraction(rossby)
    dipolar = dipolar_field_t(power_w, outer_m, inner_m, density_kg_m3)
    weaker = min(1.0, share / dipole_fraction(EARTH_ROSSBY))
    moment = 4.0 * math.pi * outer_m ** 3 * dipolar * weaker * max(min(attenuation, 1.0), 0.0) / VACUUM_PERMEABILITY
    return Field(moment, rossby, share, dipolar, outer_m)


def giant_power_w(luminosity_w):
    """The power that drives a giant's dynamo from the heat it radiates, W (`GIANT_POWER_SHARE`)."""
    return GIANT_POWER_SHARE * max(luminosity_w, 0.0)


@dataclass(frozen=True)
class Lifetime:
    """
    When a dynamo ran.

    Attributes:
        started_gy (float or None): The first age it ran at; `None` if it never did.
        ended_gy (float or None): The age it stopped, taken as the first age given after the last one it ran at; `None` while it still runs or never did.
        running (bool): Whether it runs at the last age given.
    """

    started_gy: float | None
    ended_gy: float | None
    running: bool


def lifetime(ages_gy, powers_w):
    """
    When the dynamo ran, from the power of the core at ascending ages (`ages_gy`, Gy, with `powers_w` the power at each): the first age with power
    above zero is the start; if a later age has none, the dynamo ended at the first such age after the last that ran.
    """
    ran = [power > 0.0 for power in powers_w]
    if not any(ran):
        return Lifetime(None, None, False)
    first = ran.index(True)
    last = len(ran) - 1 - ran[::-1].index(True)
    if last == len(ran) - 1:
        return Lifetime(ages_gy[first], None, True)
    return Lifetime(ages_gy[first], ages_gy[last + 1], False)
