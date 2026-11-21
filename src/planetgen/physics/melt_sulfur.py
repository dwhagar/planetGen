# planetgen/physics/melt_sulfur.py

"""
The sulfur a lava world's melt holds (GEN.1814; stand-ins section 6.319)
=======================================================================

`sulfur_balance` gives the sulfur two solid sinks in the crust, one calcium or iron per sulfur. A body whose surface is molten has a third
and larger one: the melt itself. Silicate melt dissolves sulfur as sulfide up to the sulfur content at sulfide saturation (SCSS), which rises
with temperature and falls with pressure. [R5340] fit 42 values from the literature and their own diamond-anvil runs (1 atm to 40 GPa, 1,573 to
4,605 K, melts normalised to 8.1 wt percent FeO; no FeO term is needed over 2 to 8.1 percent):

    ln S[ppm] = A + B / T + C P[GPa] / T          A 11.79, B -7073, C -92        (R2 0.86; each fitted with an error, in the registry)

This module turns it into a capacity in mol of sulfur per square metre of the surface, like the crust capacities of `sulfur_balance`:

* **The melt column.** The molten layer reaches down to the pressure at which the solidus has risen to the surface temperature: the solidus the tidal heating already uses
  (`tidal_heating.solidus_k`, the Simon-Glatzel curve through 1,393 K at the surface, 2,200 K at 10 GPa and 4,180 K at 135 GPa, [R1375], registry TH13), inverted, so with the surface
  temperature `T` the melt reaches `a ((T / T0)^(1 / b) - 1)` GPa, at most the pressure at the base of the body's magma ocean (`magma_ocean_redox.base_pressure_gpa`: a layer cannot be
  deeper than the ocean). The mass of the column is that pressure over the gravity. A surface under `MELT_MIN_K` has no melt and no capacity.
* **The temperature.** The melt is taken at the surface temperature all the way down (a lower bound: the melt gets hotter with depth, and SCSS rises with temperature), held at the fit's top
  `FIT_MAX_K` for a body hotter than that. The pressure term runs on through the column to its base, past the fit's 40 GPa on a deep ocean.
* **The capacity** is the column's mean SCSS (the integral of the law over the pressure, `scss_mean_ppm`) times its mass, over the molar mass of sulfur.

What it does not do: SCSS is a limit when a sulfide stands in equilibrium with the melt; how much sulfur dissolves from a gas depends on the oxygen and sulfur fugacities of both, so the
capacity is an upper bound; sulfate (the form sulfur takes at high oxygen fugacity) is more soluble and is not added; the melt keeps what it takes (no burden, no lifetime: it is molten); a
surface that freezes later gives the sulfur back to nothing here (the air history of a world is not re-solved for it).

Every number that is not read from a source is a design constant with its tier and tweak level in `docs/design/research/gap-filling-plan.md` (section 6.319).
"""

import math

from planetgen.physics import constants, elements, magma_ocean_redox, mixture_lava, tidal_heating

SCSS_A = 11.79
"""float: The constant of the SCSS law, ln of ppm S ([R5340: Eq. 2], fitted with an error of 0.62). Confident, tier 2: read from the journal's open page; its 2024 corrigendum could not be read."""

SCSS_B_K = -7073.0
"""float: The temperature term of the SCSS law, K (fitted with an error of 993)."""

SCSS_C_K_PER_GPA = -92.0
"""float: The pressure term of the SCSS law, K per GPa: `C P / T` (fitted with an error of 36)."""

FIT_MAX_K = 4605.0
"""float: The hottest run of the fit, K ([R5340], supplementary section S.5); a hotter melt is held at it, so the capacity of a body above it is a floor."""

S_MOLAR_MASS_KG = elements.molar_mass("S")
"""float: Kilograms per mole of sulfur."""


def scss_ppm(temperature_k, pressure_gpa):
    """The sulfur content at sulfide saturation of a basalt melt at `temperature_k` and `pressure_gpa`, ppm by mass (`SCSS_A` ... `SCSS_C_K_PER_GPA`); the temperature is held at `FIT_MAX_K`."""
    t = min(temperature_k, FIT_MAX_K)
    return math.exp(SCSS_A + SCSS_B_K / t + SCSS_C_K_PER_GPA * pressure_gpa / t)


def scss_mean_ppm(temperature_k, base_gpa):
    """The mean of `scss_ppm` over the melt column, which is uniform in pressure from the surface to `base_gpa` (the mass of a column is its pressure over the gravity)."""
    if base_gpa <= 0.0:
        return scss_ppm(temperature_k, 0.0)
    t = min(temperature_k, FIT_MAX_K)
    k = -SCSS_C_K_PER_GPA / t
    return math.exp(SCSS_A + SCSS_B_K / t) * (1.0 - math.exp(-k * base_gpa)) / (k * base_gpa)


def solidus_depth_gpa(temperature_k):
    """The pressure, GPa, at which the solidus of `tidal_heating.solidus_k` has risen to `temperature_k`: `a ((T / T0)^(1 / b) - 1)`; 0 at or under the surface solidus."""
    if temperature_k <= tidal_heating.SOLIDUS_SURFACE_K:
        return 0.0
    return tidal_heating.SOLIDUS_PRESSURE_SCALE_GPA * ((temperature_k / tidal_heating.SOLIDUS_SURFACE_K) ** (1.0 / tidal_heating.SOLIDUS_EXPONENT) - 1.0)


def melt_base_gpa(temperature_k, mass_kg):
    """
    The pressure at the base of a body's molten layer, GPa: where the solidus (`solidus_depth_gpa`) has risen to the surface temperature, at most the pressure at the base of the
    body's magma ocean (`magma_ocean_redox.base_pressure_gpa` of its mass). 0 for a surface under `MELT_MIN_K`.
    """
    if temperature_k is None or temperature_k <= mixture_lava.MELT_MIN_K or mass_kg <= 0.0:
        return 0.0
    ocean = magma_ocean_redox.base_pressure_gpa(mass_kg / constants.EARTH_MASS_TO_KG)
    return min(solidus_depth_gpa(temperature_k), ocean)


def capacity_mol_m2(temperature_k, gravity_m_s2, radius_m):
    """
    The sulfur a body's molten layer can hold, mol per square metre of its surface; 0 when the surface is not molten.

    Args:
        temperature_k (float or None): The surface temperature.
        gravity_m_s2 (float): The surface gravity.
        radius_m (float): The radius (the mass of the body is `g R^2 / G`).
    """
    if not gravity_m_s2 or gravity_m_s2 <= 0.0 or not radius_m or radius_m <= 0.0:
        return 0.0
    mass_kg = gravity_m_s2 * radius_m ** 2 / constants.G
    base = melt_base_gpa(temperature_k, mass_kg)
    if base <= 0.0:
        return 0.0
    column_kg_m2 = base * 1.0e9 / gravity_m_s2
    return column_kg_m2 * scss_mean_ppm(temperature_k, base) * 1.0e-6 / S_MOLAR_MASS_KG
