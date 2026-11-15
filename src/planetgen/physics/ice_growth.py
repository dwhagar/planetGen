# planetgen/physics/ice_growth.py

"""
How thick the ice and the snow on it are (GEN.1728, docs/design/habitability/ice-and-snow-growth.md)
====================================================================================================

`physics/ice_light.py` (GEN.1726) says how much light a lid of a given thickness lets through, and takes the thickness from its caller. This module is where the caller
gets it: the ice a world's ocean or lake wears and the snow lying on it, from the three things the generator already knows about a world, its mean surface temperature,
the swing of its year and the starlight its surface absorbs (and from the heat the water gives up, which is the world's own).

The ice. It grows from below, as the cold above carries heat out through it, and the thicker it is the slower it grows: Stefan's law, `H = sqrt(2 k u / (rho L))` for `u` the
freezing degree days (the temperature shortfall under freezing, added up over the freezing season) and `k`, `rho`, `L` the ice's conductivity, density and latent heat. Snow on
the ice insulates it, and the air takes its share of the drop, so the thickness is the form `H = sqrt((k he)^2 + 2 k u / (rho L)) - k he` with `he = h_s / k_s + 1 / h_air`
([R5277], equations 18 and 21). The degree days come from the year's cycle as a sine curve about the mean, so a world needs only its mean and its swing. A lid that never thaws
(its warmest month under freezing) has no summer to start over and settles where the cold the ice conducts down balances the heat the water gives up from below
(`perennial_thickness_m`, the balance of [R5279]); the thickness it can never pass is the conduction limit of the ice shell of an icy body (`hydrosphere.ice_shell_thickness_km`).

The snow. The snow that falls in a year is the share of the starlight the surface absorbs that goes into the water cycle (the Earth's own share, the same energy
limit the weathering law of [R1660] uses), turned to a depth of snow of the density measured on Arctic ice, and only the share that stays on the ice counts (the rest blows
into leads and sublimes): the Arctic's 17 cm of water falling in the cold season leaves 11 cm of water standing in May, 34 cm of snow of 300 kg per cubic metre [R5278].
A world with a short freezing season, a faint star or a dry cold air gets little. The snow of the freezing season is the snow on the ice.

What it leaves out. The melt of the summer (a world whose summer thaws part of the lid keeps the seasonal thickness, which is the thicker of the two cases, and so errs toward
less light under the ice); snow-ice, the slush that forms where the weight of the snow pushes the ice under the water; ridging and drift; the salt of the ice (the melting point is
the caller's `freezing_k`, taken from the world's own water). Every constant is a stand-in in section 6.298 of `docs/design/research/gap-filling-plan.md`.

Use:

    lid = lid_of(mean_k=258.0, amplitude_k=18.0, absorbed_flux_w_m2=150.0, ocean_flux_w_m2=2.0)
    lid.ice_m, lid.snow_m, lid.perennial                  # a first-year lid: about 1 m of ice, 0.1 m of snow, False
    ice_light.transmission(lid.ice_m, lid.snow_m)         # the light under it (GEN.1726)
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, continental_weathering, hydrosphere

FREEZING_K = constants.WATER_MELTING_K
"""float: The freezing point of fresh water, K: the default `freezing_k` (a salty ocean passes its own, which is lower)."""

SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.25

ICE_DENSITY_KG_M3 = hydrosphere.ICE_DENSITY_KG_M3
"""float: The density of the ice, reused from the hydrosphere ([R5277] uses 917)."""

LATENT_HEAT_J_KG = 3.352e5
"""float: The latent heat of fusion of water, joules per kg: 3.88 W d per kg, which is 3.352e5 J per kg [R5277]. [IG1, T1, L0.]"""

AIR_COUPLING_W_M2_K = 17.5
"""float: `h_air`, the heat the ice or snow surface gives to the air per metre squared and kelvin, W per m2 and K: [R5277] gives 10 to 25, the middle. [IG2, T1, L0.]"""
AIR_COUPLING_RANGE_W_M2_K = (10.0, 25.0)
"""tuple: The range of `h_air`, from the same source. [IG2, T1, L0.]"""

SNOW_DENSITY_KG_M3 = 300.0
"""float: The density of snow lying on sea ice, kg per cubic metre: 300 averaged over the season, 250 in September to 320 in May [R5278]. [IG3, T1, L0.]"""
SNOW_DENSITY_RANGE_KG_M3 = (250.0, 320.0)
"""tuple: The range of the snow's density, the Arctic season's. [IG3, T1, L0.]"""

SNOW_CONDUCTIVITY_COEFFICIENT = 2.22362
SNOW_CONDUCTIVITY_EXPONENT = 1.885
"""floats: The snow's conductivity `k_s = 2.22362 rho^1.885` in W per m and K, for `rho` in Mg per cubic metre ([R5277]: 0.23 at 300 kg per m3). [IG4, T1, L0.]"""

OCEAN_FLUX_W_M2 = 2.0
"""float: The heat the water gives up to the bottom of the ice in the Arctic, W per m2: the ocean must supply 1 to 2 kcal per cm2 and year, about 1.3 to 2.7, and the standard case is 2.0 [R5279]. The
default `ocean_flux_w_m2` of a world with no better figure; an icy body passes its own interior heat. [IG5, T1, L0.]"""

ABLATION_M_YR = 0.4
"""float: The ice a year wears off the top of a lid that does not thaw, metres of ice per year: 40 cm in the Arctic's standard case ([R5279]: surface ablation 40 cm, bottom accretion 45 cm). The ablation of a world
with another air is not known; this is the one measured case. [IG8, T7, L1.]"""

PRECIPITATION_ENERGY_SHARE = continental_weathering.PRECIPITATION_REF_M_YR * continental_weathering.latent_j_m3(continental_weathering.REFERENCE_K) / continental_weathering.SECONDS_PER_YEAR / continental_weathering.EARTH_ABSORBED_FLUX_W_M2
"""float: The share of the absorbed starlight that goes into the water cycle: the Earth's, 0.99 m of rain a year at 288 K for 238 W per m2 absorbed, which is the same energy limit the weathering law uses [R1660]; about 0.33. [IG6, T7, L1.]"""

SNOW_RETENTION = 11.0 / 17.0
"""float: The share of the snow that falls on the ice that is still on it at the end of the season: 11 cm of water standing in May [R5278] of 17 cm that fell (5 cm of water a month in September and October,
1 cm a month from November to March, 2 cm in April; the rates [R5278] quotes); the rest blows into leads and sublimes. About 0.65. [IG7, T7, L1.]"""

WATER_DENSITY_KG_M3 = 1000.0


class IceGrowthError(ValueError):
    """An input that is not a temperature, a flux or a depth."""


@dataclass(frozen=True)
class Lid:
    """The lid of a world's water: the ice thickness, the snow depth on it, both in metres, and whether the lid is perennial (never thaws)."""

    ice_m: float
    snow_m: float
    perennial: bool


def _check_positive(value, name):
    if not (value > 0.0 and math.isfinite(value)):
        raise IceGrowthError("%s must be a finite number above 0, not %r" % (name, value))


def _check_nonneg(value, name):
    if not (value >= 0.0 and math.isfinite(value)):
        raise IceGrowthError("%s must be a finite number of 0 or more, not %r" % (name, value))


def snow_conductivity_w_m_k(density_kg_m3=SNOW_DENSITY_KG_M3):
    """The snow's thermal conductivity, W per m and K ([R5277], IG4)."""
    _check_positive(density_kg_m3, "the snow density")
    return SNOW_CONDUCTIVITY_COEFFICIENT * (density_kg_m3 / 1000.0) ** SNOW_CONDUCTIVITY_EXPONENT


def ice_conductivity_w_m_k(temperature_k):
    """The ice's thermal conductivity, W per m and K: `k = A / T` of the hydrosphere's ice shell (2.07 at the melting point, as [R5277])."""
    _check_positive(temperature_k, "the temperature")
    return hydrosphere.ICE_CONDUCTIVITY_A_W_M / temperature_k


def frozen_share(mean_k, amplitude_k, freezing_k=FREEZING_K):
    """The share of the year the surface is under freezing, for a sine year of mean `mean_k` and half-range `amplitude_k`: 1 when the warmest month is under it, 0 when it never freezes."""
    _check_nonneg(amplitude_k, "the amplitude of the year")
    if amplitude_k == 0.0:
        return 1.0 if mean_k < freezing_k else 0.0
    x = (freezing_k - mean_k) / amplitude_k
    if x >= 1.0:
        return 1.0
    if x <= -1.0:
        return 0.0
    return 1.0 - math.acos(x) / math.pi


def freezing_degree_days(mean_k, amplitude_k, freezing_k=FREEZING_K, days_per_year=DAYS_PER_YEAR):
    """The freezing degree days of one year, K days: the shortfall under `freezing_k`, added up over the time the surface is under it, for a sine year of mean `mean_k` and half-range `amplitude_k`."""
    _check_nonneg(amplitude_k, "the amplitude of the year")
    if amplitude_k == 0.0:
        return max(freezing_k - mean_k, 0.0) * days_per_year
    x = (freezing_k - mean_k) / amplitude_k
    if x >= 1.0:
        return (freezing_k - mean_k) * days_per_year
    if x <= -1.0:
        return 0.0
    phi = math.acos(x)
    return amplitude_k * days_per_year / (2.0 * math.pi) * (x * (2.0 * math.pi - 2.0 * phi) + 2.0 * math.sin(phi))


def winter_growth_m(degree_days_k_day, snow_m=0.0, ice_temperature_k=FREEZING_K - 10.0, snow_density_kg_m3=SNOW_DENSITY_KG_M3, air_coupling_w_m2_k=AIR_COUPLING_W_M2_K):
    """
    The ice a freezing season grows, metres: `sqrt((k he)^2 + 2 k u / (rho L)) - k he` ([R5277] eq. 21), `u` the freezing degree days (converted to K s), `he = h_s / k_s + 1 / h_air` the
    resistance of the snow and the air, `k` the ice's conductivity at `ice_temperature_k`. With no snow `he` is the air's alone.
    """
    _check_nonneg(degree_days_k_day, "the freezing degree days")
    _check_nonneg(snow_m, "the snow depth")
    _check_positive(air_coupling_w_m2_k, "the air coupling")
    k = ice_conductivity_w_m_k(ice_temperature_k)
    he = snow_m / snow_conductivity_w_m_k(snow_density_kg_m3) + 1.0 / air_coupling_w_m2_k
    u = degree_days_k_day * SECONDS_PER_DAY
    return math.sqrt((k * he) ** 2 + 2.0 * k * u / (ICE_DENSITY_KG_M3 * LATENT_HEAT_J_KG)) - k * he


def conduction_limit_m(cold_k, flux_w_m2, freezing_k=FREEZING_K):
    """The thickness the ice cannot pass, metres: where the cold of `cold_k` carried through the ice balances the heat `flux_w_m2` the water gives up (the ice shell of the hydrosphere)."""
    _check_positive(flux_w_m2, "the heat flux")
    return hydrosphere.ice_shell_thickness_km(flux_w_m2, cold_k, freezing_k) * 1000.0


def perennial_thickness_m(mean_k, ocean_flux_w_m2=OCEAN_FLUX_W_M2, snow_m=0.0, ablation_m_yr=ABLATION_M_YR, freezing_k=FREEZING_K, snow_density_kg_m3=SNOW_DENSITY_KG_M3):
    """
    The thickness a lid settles at when it never thaws, metres: the ice grows from below as fast as the surface wears it away and the water warms it. The heat the ice conducts is
    `k (Tf - T) / H` and the part that freezes new ice is what is left after the water's flux and the ablation `rho L a`, so
    `H = k ((Tf - T) / (F + rho L a) - h_s / k_s)` ([R5279]'s balance in its annual mean form). With no ablation it is the linear form of the conduction limit. 0 when the mean is at or over freezing.
    """
    _check_positive(ocean_flux_w_m2, "the heat flux")
    _check_nonneg(ablation_m_yr, "the ablation")
    if mean_k >= freezing_k:
        return 0.0
    k = ice_conductivity_w_m_k(0.5 * (freezing_k + mean_k))
    wear = ICE_DENSITY_KG_M3 * LATENT_HEAT_J_KG * ablation_m_yr / (DAYS_PER_YEAR * SECONDS_PER_DAY)
    return max(k * ((freezing_k - mean_k) / (ocean_flux_w_m2 + wear) - snow_m / snow_conductivity_w_m_k(snow_density_kg_m3)), 0.0)


def snowfall_m_yr(absorbed_flux_w_m2, mean_k):
    """
    The snow that falls in a year, metres of water: the water cycle's share (`PRECIPITATION_ENERGY_SHARE`) of the absorbed starlight, over the latent heat of the water at the mean
    temperature (the weathering law's own `latent_j_m3`), so a faint star or a very cold air gets little.
    """
    _check_nonneg(absorbed_flux_w_m2, "the absorbed starlight")
    if mean_k <= continental_weathering.LATENT_OFFSET_K:
        return 0.0
    return PRECIPITATION_ENERGY_SHARE * absorbed_flux_w_m2 / continental_weathering.latent_j_m3(mean_k) * continental_weathering.SECONDS_PER_YEAR


def snow_depth_m(absorbed_flux_w_m2, mean_k, amplitude_k, freezing_k=FREEZING_K, snow_density_kg_m3=SNOW_DENSITY_KG_M3):
    """
    The depth of the snow on the ice at the end of the freezing season, metres: the water that falls in the freezing share of the year, the share that stays (`SNOW_RETENTION`), turned to
    snow of `snow_density_kg_m3`. 0 for a year that never freezes.
    """
    _check_positive(snow_density_kg_m3, "the snow density")
    water_m = snowfall_m_yr(absorbed_flux_w_m2, mean_k) * frozen_share(mean_k, amplitude_k, freezing_k) * SNOW_RETENTION
    return water_m * WATER_DENSITY_KG_M3 / snow_density_kg_m3


def lid_of(mean_k, amplitude_k, absorbed_flux_w_m2, ocean_flux_w_m2=OCEAN_FLUX_W_M2, freezing_k=FREEZING_K, snow_density_kg_m3=SNOW_DENSITY_KG_M3, ablation_m_yr=ABLATION_M_YR):
    """
    The ice and the snow on a world's water: the snow of the freezing season, and the ice that season grows (half the snow insulates it, since it falls through the season), held under the
    conduction limit of the coldest month. A year whose warmest month is under freezing never thaws and settles at the perennial thickness (the thicker of the two, no summer wears it).
    A year that never freezes has no lid.
    """
    _check_nonneg(amplitude_k, "the amplitude of the year")
    if frozen_share(mean_k, amplitude_k, freezing_k) <= 0.0:
        return Lid(0.0, 0.0, False)
    snow = snow_depth_m(absorbed_flux_w_m2, mean_k, amplitude_k, freezing_k, snow_density_kg_m3)
    cold_k = max(mean_k - amplitude_k, 1.0)
    degree_days = freezing_degree_days(mean_k, amplitude_k, freezing_k)
    seasonal = winter_growth_m(degree_days, 0.5 * snow, 0.5 * (freezing_k + cold_k), snow_density_kg_m3)
    seasonal = min(seasonal, conduction_limit_m(cold_k, ocean_flux_w_m2, freezing_k))
    if mean_k + amplitude_k < freezing_k:
        return Lid(max(seasonal, perennial_thickness_m(mean_k, ocean_flux_w_m2, snow, ablation_m_yr, freezing_k, snow_density_kg_m3)), snow, True)
    return Lid(seasonal, snow, False)
