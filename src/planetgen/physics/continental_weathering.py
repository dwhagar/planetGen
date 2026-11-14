# planetgen/physics/continental_weathering.py

"""
The carbon balance of a world with land and an ocean: continental silicate weathering with its published limits, and the seafloor sink beside it (GEN.1713)
===========================================================================================================================================================

Volcanoes give CO2; two sinks take it back; the air holds the pressure at which they balance (`surface_pressure.balance_root_bar`). This
module is the balance of a body with land and an ocean that is not sealed (`waterworld_carbon.sink_kind`).

* *Land.* The weathering law of Graham and Pierrehumbert 2020 [R1660: Eqs. 24, 40 to 44, Table 1], the form of Maher and Chamberlain 2014 that
  adds to the kinetic rate of the Walker law (Walker, Hays and Kasting 1981; Lehmer et al. 2020 [R2628]) the three limits that a rate has
  on real ground. Per unit of land

      w = alpha / { [k_ref exp((T - T_ref) / T_e) (p / p_ref)^beta]^-1 + m A t_s + alpha [q Lambda p^n]^-1 }.

  The first term is the kinetics of dissolution (the Walker law's temperature and pressure dependence); `m A t_s` is the supply limit (the
  fresh rock the soil can offer; `w` cannot pass `alpha / (m A t_s)`); the last term is the thermodynamic limit (the water that runs off
  holds at most `Lambda p^n` of dissolved load before clay precipitates). The runoff is `q = Gamma P`, and the precipitation `P` rises 3
  percent per K but is held under the energy limit `S_abs / L(T)`: the starlight the surface absorbs cannot evaporate more water than its
  latent heat allows (Eqs. 41 to 44). On a hot world that limit, not the Clausius-Clapeyron rise, sets the rain, and so the weathering.
* *Seafloor.* The dissolution of basalt in the pore water of the ocean floor (`waterworld_carbon.seafloor_rate_over_modern`, [R2489]),
  at the modern share of the volcanoes' CO2 (`waterworld_carbon.SEAFLOOR_SHARE_OF_OUTGASSING`, 0.45 / 6.6). Its limit is the supply of
  fresh crust, which is the spreading term `(q / q_E)^beta` already in that law (Chambers 2020 [R5254]: "limited by the supply of fresh
  oceanic crust, and proportional to the seafloor spreading rate"; Hoening, Tosi and Spohn 2019 give no further cap).
* *The balance.* In units of the modern outgassing, with `s` the seafloor share,

      (1 - s) (land / land_E) w_rel(p, T) + s r_sea(p, T, q) = (q / q_E)^m,

  where `w_rel` is `w` over its own value at the Earth's reference (289 K, 288 microbar, 238 W per square metre absorbed), so that the Earth is
  the balance at its own pressure whatever the drawn constants are (the Walker law's identity, which Lehmer et al. keep).

*Normalisation to the Earth.* GP20's Table 1 does not reproduce the modern weathering it is meant to describe: Eq. 40 with Table 1 at the
modern reference gives 0.0184 mol per square metre of land per year, 0.0053 per square metre of the Earth's surface (land share gamma 0.3), where their reference flux
V_ref over the Earth's area is 0.0147: a factor of 0.36. The generator therefore uses the law's shape (its temperature, pressure, runoff and limits) and
sets its level by the Earth identity above (stand-in AM110). The thermodynamic limit `Lambda`, which GP20 find dominates, is drawn over the
range they tested (`LAMBDA_SCALE`, a tenth to ten times its Table 1 value).

Every number that is not read from a source is a design constant with its source tier and tweak level in
`docs/design/research/gap-filling-plan.md` (GEN.725 rule).
"""

import math

from planetgen.physics import surface_pressure, waterworld_carbon
from planetgen.physics.mixture_base import bracket

# --- Graham and Pierrehumbert 2020, Table 1 (read) -----------------------------------------------------------------------------------

THERMODYNAMIC_CONSTANT_MOL_M3 = 1.4
"""float: `Lambda`, the dissolved load the runoff holds at one bar of CO2, mol per cubic metre (1.4e-3 mol per litre, Table 1; READ, tier 1)."""

THERMODYNAMIC_EXPONENT = 0.316
"""float: `n` of the equilibrium concentration `C_eq = Lambda p^n` (Table 1; READ, tier 1)."""

ALPHA = 3.39e5
"""float: `alpha` of the thermodynamic limit, which includes the clay-precipitation factor (Table 1; READ, tier 1)."""

SOIL_AGE_YR = 1.0e5
"""float: `t_s`, the age of the soil of the supply limit, years (Table 1; READ, tier 1)."""

MOLAR_MASS_KG_MOL = 0.27
"""float: `m`, the molar mass of the dissolving rock, kg per mol (Table 1; READ, tier 1)."""

SPECIFIC_AREA_M2_KG = 100.0
"""float: `A`, the specific surface area of the weathering grains, square metres per kg (Table 1; READ, tier 1)."""

KINETIC_RATE_REF = 8.7e-6
"""float: `k_eff,ref`, the kinetic rate at the reference temperature and pressure, mol per square metre per year (Table 1; READ, tier 1)."""

KINETIC_PCO2_EXPONENT = 0.2
"""float: `beta`, the exponent of the CO2 pressure in the kinetic rate (Table 1; READ, tier 1; the paper finds it and `T_e` hardly matter)."""

KINETIC_T_E_K = 11.1
"""float: `T_e`, the e-folding temperature of the kinetic rate, K (Table 1; READ, tier 1)."""

RUNOFF_SHARE = 0.2
"""float: `Gamma`, the share of the precipitation that runs off (Table 1; READ, tier 1)."""

PRECIPITATION_PER_K = 0.03
"""float: `epsilon`, the rise of precipitation per K (Table 1; READ, tier 1)."""

PRECIPITATION_REF_M_YR = 0.99
"""float: `p_ref`, the precipitation at the reference temperature, metres per year (Table 1; READ, tier 1)."""

REFERENCE_K = 288.0
"""float: `T_ref` of the paper's Table 1, K."""

REFERENCE_PCO2_BAR = 280.0e-6
"""float: `pCO2_ref` of the paper's Table 1, bar."""

LATENT_COEFFICIENT_J_M3 = 1.918e9
"""float: The latent heat of the energy limit, `L(T) = 1.918e9 (T / (T - 33.91))^2` joules per cubic metre of liquid water (Eq. 44; READ, tier 1)."""

LATENT_OFFSET_K = 33.91
"""float: The offset of the same fit, K."""

SECONDS_PER_YEAR = 3.15576e7
"""float: The Julian year, to turn the energy limit's metres per second into metres per year."""

EARTH_ABSORBED_FLUX_W_M2 = 238.0
"""float: The starlight the Earth's surface absorbs, `S (1 - a) / 4` for 1360 W per square metre and a Bond albedo of 0.3: the reference of the Earth identity."""

# --- the drawn constants --------------------------------------------------------------------------------------------------------------

LAMBDA_SCALE = (0.1, 10.0, "log")
"""tuple: The factor of `Lambda` (the thermodynamic limit): the range GP20 tested, a tenth to ten times their value; they find it the one constant that matters
(stand-in AM109, tier T1, tweak level L0)."""


def draw_cycle(constants):
    """`(lambda_scale, m)`: the factor of the thermodynamic limit and the exponent of the heat flow in the outgassing, drawn and recorded in `constants`."""
    return (bracket(constants, "weathering_lambda_scale", LAMBDA_SCALE), bracket(constants, "outgassing_exponent", surface_pressure.OUTGASSING_EXPONENT))


# --- the land law -----------------------------------------------------------------------------------------------------------------

def runoff_m_yr(temperature_k, absorbed_flux_w_m2):
    """
    The runoff, metres per year: `Gamma` times the precipitation, which rises `PRECIPITATION_PER_K` per K from `PRECIPITATION_REF_M_YR` but is held under the
    energy limit `S_abs / L(T)` (GP20 Eqs. 41 to 44). Nil where the linear rise has run out (under about 255 K) or the temperature is under the offset of the fit.
    """
    if temperature_k <= LATENT_OFFSET_K:
        return 0.0
    precipitation = PRECIPITATION_REF_M_YR * (1.0 + PRECIPITATION_PER_K * (temperature_k - REFERENCE_K))
    latent = LATENT_COEFFICIENT_J_M3 * (temperature_k / (temperature_k - LATENT_OFFSET_K)) ** 2
    energy_limit = max(absorbed_flux_w_m2, 0.0) / latent * SECONDS_PER_YEAR
    return RUNOFF_SHARE * max(min(precipitation, energy_limit), 0.0)


def rate_mol_m2_yr(temperature_k, pco2_bar, absorbed_flux_w_m2, lambda_scale=1.0):
    """The silicate weathering per square metre of land and year, mol (GP20 Eq. 40, the law of the module's docstring)."""
    runoff = runoff_m_yr(temperature_k, absorbed_flux_w_m2)
    pco2 = max(pco2_bar, surface_pressure.MIN_PCO2_BAR)
    if runoff <= 0.0:
        return 0.0
    kinetic = KINETIC_RATE_REF * math.exp((temperature_k - REFERENCE_K) / KINETIC_T_E_K) * (pco2 / REFERENCE_PCO2_BAR) ** KINETIC_PCO2_EXPONENT
    supply = MOLAR_MASS_KG_MOL * SPECIFIC_AREA_M2_KG * SOIL_AGE_YR
    thermodynamic = ALPHA / (runoff * lambda_scale * THERMODYNAMIC_CONSTANT_MOL_M3 * pco2 ** THERMODYNAMIC_EXPONENT)
    return ALPHA / (1.0 / kinetic + supply + thermodynamic)


def rate_over_modern(temperature_k, pco2_bar, absorbed_flux_w_m2, lambda_scale=1.0):
    """`rate_mol_m2_yr` over its value at the Earth's reference (`TS_MODERN_K`, `CO2_MODERN_BAR`, `EARTH_ABSORBED_FLUX_W_M2`) with the same `lambda_scale`: the Earth identity."""
    modern = rate_mol_m2_yr(surface_pressure.TS_MODERN_K, surface_pressure.CO2_MODERN_BAR, EARTH_ABSORBED_FLUX_W_M2, lambda_scale)
    return rate_mol_m2_yr(temperature_k, pco2_bar, absorbed_flux_w_m2, lambda_scale) / modern


# --- the balance ------------------------------------------------------------------------------------------------------------------

def pco2_bar(supply_bar, heat_flux_w_m2, land_fraction, temperature_k, absorbed_flux_w_m2, cycle, seafloor_law=None):
    """
    The CO2 pressure at which the land and the seafloor together take as much as the volcanoes give (the stable root).

    `(1 - s) (land / land_E) w_rel + s r_sea = (q / q_E)^m`; both sinks rise with the pressure, so the root is unique. The land share is the caller's (the
    ice-free land, `surface_air._weathering_land`).

    Args:
        supply_bar (float): The CO2 the volcanoes have delivered, as a pressure: all of it if the sinks cannot keep up.
        heat_flux_w_m2 (float): The body's heat flow `q`.
        land_fraction (float): The share of the surface that weathers (above sea level and not under ice).
        temperature_k (float): The surface temperature.
        absorbed_flux_w_m2 (float): The starlight the surface absorbs at that temperature (`Setting.absorbed_flux_at`), the energy limit's supply.
        cycle (tuple): `draw_cycle`.
        seafloor_law (SeafloorLaw): `waterworld_carbon.draw_seafloor_law`, or `None` for a body with no ocean floor (no seafloor sink).

    Returns:
        float: The CO2 pressure, bar: at most `supply_bar`.
    """
    lambda_scale, exponent = cycle
    outgassing = max(heat_flux_w_m2 / surface_pressure.HEAT_FLUX_EARTH_W_M2, 0.0) ** exponent
    if outgassing <= 0.0:
        return surface_pressure.MIN_PCO2_BAR if supply_bar > surface_pressure.MIN_PCO2_BAR else max(supply_bar, 0.0)
    seafloor_share = waterworld_carbon.SEAFLOOR_SHARE_OF_OUTGASSING if seafloor_law is not None else 0.0
    land = max(land_fraction, 0.0) / surface_pressure.LAND_FRACTION_EARTH

    def excess(p_bar):
        taken = (1.0 - seafloor_share) * land * rate_over_modern(temperature_k, p_bar, absorbed_flux_w_m2, lambda_scale)
        if seafloor_law is not None:
            taken += seafloor_share * waterworld_carbon.seafloor_rate_over_modern(p_bar, temperature_k, heat_flux_w_m2, seafloor_law)
        return taken - outgassing

    return surface_pressure.balance_root_bar(supply_bar, excess)
