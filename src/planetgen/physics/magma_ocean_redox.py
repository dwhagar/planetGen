# planetgen/physics/magma_ocean_redox.py

"""
The oxygen fugacity a magma ocean leaves at its surface, from the pressure at its base (GEN.1005, part 1; docs/design/composition/magma-ocean-redox-model.md).

Rows a1 and a2 of the atmosphere possibility table (`atmosphere_recipes`) drew the redox of a rocky body's outgassing mantle from a line through three Solar System
bodies (the Moon, Mars and the Earth). This module carries the physics those three points sample, in terms of the mass of the body alone:

1. **Where the ocean ends.** A magma ocean on a body of `mass_earth` Earth masses reaches the pressure `base_pressure_gpa` at its base: a power law of the mass, capped at the
   140 GPa where the first-principles calculations of Deng et al. 2020 stop. [Stand-ins MR1, MR2.]
2. **What the pressure does to iron.** Ferrous iron in the melt disproportionates, 3 FeO = Fe + Fe2O3, and the product with the smaller molar volume is Fe2O3 (the volume difference
   stays positive to 140 GPa, Deng et al. 2020), so the higher the pressure at the base the more of the iron is ferric when the melt equilibrates with metal. Deng et al. 2020
   write that the ferric share drops by a factor of two when the base moves from 55 to 25 GPa, a power of the pressure of 0.88; Hirschmann 2022 gives a slope of about 0.5 from its own
   ocean-base conditions, and the model takes the middle, 0.75, and carries the bracket 0.5 to 1.0 into the width (`FE3_PRESSURE_EXPONENT`, `FE3_EXPONENT_BRACKET`). The ferric share
   rides up with the melt and is the share at the surface, where the ocean meets the gas. [Stand-ins MR3, MR4.]
3. **What the ferric share says about the oxygen.** 4 FeO + O2 = 2 Fe2O3 puts the oxygen fugacity as the fourth power of the ratio Fe3+/Fe2+, so dIW (log10 of the fugacity above
   the iron-wuestite buffer, the scale of Hirschmann 2021) is `4 log10(Fe3+/Fe2+) + 5.62`. The 5.62 is read off the pairs of Hirschmann 2022 for the Earth's ocean: Fe3+/FeT 0.034 gives
   dIW -0.2 and 0.103 gives +1.8, the 4 log10 relation reproduces both to 0.06. The buffer's own temperature term is not in it (the surface temperature is one value). [Stand-in MR5.]
4. **The spread.** The ferric share at the reference pressure is a range of the literature (0.034 to 0.103) and the power is a bracket (0.5 to 1.0); both are uniform, their log
   spreads add in quadrature, and `diw_sigma` carries that as the standard deviation of the draw. [Stand-in MR6.]

`later_mantle_diw` is the later upper mantle's value: the magma-ocean surface of the body (at most that of one Earth mass: nothing in the sources lets the later mantle rise
beyond) plus an oxidation increment that grows with the mass. [Stand-in MR7.]

Tests against measurement (never inputs, `tests/test_magma_ocean_redox.py`): the magma-ocean surface of the Moon, Mars and the Earth from Hirschmann 2022 (-2.2 to -1.5, -1.4 to
+0.2 and -1.2 to +1.8) and Deng et al. 2020 (the base pressures 5, 14 and 55 GPa; the surface values IW-2, IW-0.3 and IW+2), each inside the model's one-sigma band.

Limits. The ocean-base pressure is a power of the mass, not an interior model: a body whose ocean is shallower or deeper than the Solar System's three rocky ones gets the wrong
pressure. The ferric share is that of the equilibrium with metal at the base and is carried to the surface unchanged (no crystallisation, no loss of H2). The graphite-flotation
alternative (Keppler and Golabek 2019, down to IW-5) is another mechanism and is not drawn. Above 5.4 Earth masses the base pressure is held at 140 GPa, so the surface redox stops
rising there: no source reaches further.
"""

import math

from planetgen.util import detmath
from planetgen.util.checks import finite_domain

BASE_PRESSURE_1_EARTH_GPA = 55.0
"""float: The pressure at the base of a magma ocean on a body of one Earth mass (Deng et al. 2020: 55 GPa; Hirschmann 2022 puts it at 28 to 53). [Stand-in MR1, tier T7, tweak L2.]"""

BASE_PRESSURE_EXPONENT = 0.55
"""float: The base pressure goes as the mass to this power: through the three base pressures of Deng et al. 2020 (5, 14 and 55 GPa at 0.0123, 0.107 and 1 Earth masses) it gives 4.9
and 16.1 GPa at the first two. The scaling is inferred from three bodies and not from an interior model. [Stand-in MR2, tier T7, tweak L2.]"""

CALCULATION_LIMIT_GPA = 140.0
"""float: The pressure up to which the first-principles equation of state is computed (Deng et al. 2020: up to 140 GPa); the base pressure is held here above 5.4 Earth masses.
[Stand-in MR1, tier T7, tweak L0.]"""

FE3_FRACTION_AT_1_EARTH_BASE = 0.068
"""float: Fe3+/FeT at the surface of an ocean whose base is at `BASE_PRESSURE_1_EARTH_GPA`: the middle of the 0.034 to 0.103 that Hirschmann 2022 gives for the Earth. [Stand-in MR3, tier T3, tweak L1.]"""

FE3_FRACTION_RANGE = (0.034, 0.103)
"""tuple: The range of `FE3_FRACTION_AT_1_EARTH_BASE` (Hirschmann 2022, section 6.3), uniform: the source of the width."""

FE3_PRESSURE_EXPONENT = 0.75
"""float: Fe3+/Fe2+ goes as the base pressure to this power: the middle of 0.5 (Hirschmann 2022's slope) and 1.0, which brackets Deng's 0.88 (a factor of two from 55 to 25 GPa). [Stand-in MR4, tier T3, tweak L1.]"""

FE3_EXPONENT_BRACKET = (0.5, 1.0)
"""tuple: The bracket of `FE3_PRESSURE_EXPONENT`, uniform: the second source of the width."""

DIW_PER_LOG_FE3_RATIO = 4.0
"""float: dIW changes by four per decade of Fe3+/Fe2+ (4 FeO + O2 = 2 Fe2O3: fO2 goes as the fourth power of the ratio). [Stand-in MR5, tier T5, tweak L0.]"""

DIW_AT_UNIT_FE3_RATIO = 5.62
"""float: dIW at Fe3+/Fe2+ = 1: read off the pairs of Hirschmann 2022 for a magma-ocean surface (0.034 gives -0.2, 0.103 gives +1.8). [Stand-in MR5, tier T3, tweak L1.]"""

OXIDATION_INCREMENT = (1.2, 2.4)
"""tuple: How far the later upper mantle sits above the magma-ocean surface of its body, log10 units: 1.2 for a body of 0.1 Earth masses or less, 2.4 at one Earth mass and above, a
straight line in the log of the mass between (the sources give a gap of 0 to 3 that grows with mass, Hirschmann 2022 section 4 and Deng et al. 2020). [Stand-in MR7, tier T7, tweak L2.]"""

LATER_MANTLE_SIGMA = 1.0
"""float: The width of the later mantle's value (the half-range of the published later-mantle values is about 1). [Stand-in MR7, tier T7, tweak L2.]"""


@finite_domain()
def base_pressure_gpa(mass_earth):
    """
    The pressure at the base of the magma ocean of a rocky body.

    Args:
        mass_earth (float): The body's mass in Earth masses, above zero.

    Returns:
        float: GPa, at most `CALCULATION_LIMIT_GPA`.
    """
    if mass_earth <= 0.0:
        raise ValueError(f"a body has a mass above zero, not {mass_earth}")
    return min(BASE_PRESSURE_1_EARTH_GPA * mass_earth ** BASE_PRESSURE_EXPONENT, CALCULATION_LIMIT_GPA)


@finite_domain()
def diw_from_fe3_fraction(fe3_fraction):
    """dIW of a melt whose ferric share of the iron is `fe3_fraction` (between 0 and 1, exclusive)."""
    if not 0.0 < fe3_fraction < 1.0:
        raise ValueError(f"the ferric share of the iron is between 0 and 1, not {fe3_fraction}")
    return DIW_PER_LOG_FE3_RATIO * math.log10(fe3_fraction / (1.0 - fe3_fraction)) + DIW_AT_UNIT_FE3_RATIO


@finite_domain()
def fe3_fraction(mass_earth, exponent=FE3_PRESSURE_EXPONENT, at_reference=FE3_FRACTION_AT_1_EARTH_BASE):
    """Fe3+/FeT at the surface of the magma ocean of a body (the ratio Fe3+/Fe2+ goes as the base pressure to `exponent`)."""
    ratio_reference = at_reference / (1.0 - at_reference)
    ratio = ratio_reference * (base_pressure_gpa(mass_earth) / BASE_PRESSURE_1_EARTH_GPA) ** exponent
    return ratio / (1.0 + ratio)


@finite_domain()
def surface_diw(mass_earth):
    """The mean dIW at the magma-ocean surface of a rocky body of `mass_earth` Earth masses."""
    return diw_from_fe3_fraction(fe3_fraction(mass_earth))


def _log_spread(low, high):
    """The standard deviation of a uniform quantity between `low` and `high`."""
    return (high - low) / math.sqrt(12.0)


@finite_domain()
def diw_sigma(mass_earth):
    """The standard deviation of the surface dIW: the log spreads of the reference ferric share and of the pressure power, in quadrature, times four."""
    low, high = FE3_FRACTION_RANGE
    share = _log_spread(math.log10(low / (1.0 - low)), math.log10(high / (1.0 - high)))
    power = _log_spread(*FE3_EXPONENT_BRACKET) * math.log10(base_pressure_gpa(mass_earth) / BASE_PRESSURE_1_EARTH_GPA)
    return DIW_PER_LOG_FE3_RATIO * detmath.hypot(share, power)


@finite_domain()
def later_mantle_diw(mass_earth):
    """
    The mean dIW of the later upper mantle: the magma-ocean surface of the body (of at most one Earth mass) plus `OXIDATION_INCREMENT`.

    Args:
        mass_earth (float): The body's mass in Earth masses, above zero.

    Returns:
        float: dIW.
    """
    low, high = OXIDATION_INCREMENT
    t = min(max(math.log10(min(mass_earth, 1.0)) + 1.0, 0.0), 1.0)
    return surface_diw(min(mass_earth, 1.0)) + low + (high - low) * t
