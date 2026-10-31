# planetgen/physics/surface_pressure.py

"""
How much gas a rocky body has delivered to its surface and kept there (GEN.726; docs/design/composition/surface-pressure-supply-law.md)
=============================================================================================================================

The pressure of a hydrostatic air is its mass over the area, `P = g M_air / (4 pi R^2)`. This module is the mass: the element
budget of the mantle (`mixture_rocky.rocky_inventory`, per 100 g of mantle), the share of it that the body's age has brought to the
surface, and the sinks that take gas back.

* *Degassing.* The share of an element at the surface is `eps_e * (1 - exp(-Theta))` with `Theta = E_release / (M_mantle e_proc)`:
  the mantle processed through melting follows the heat the interior has carried out (`interior_heat.Heat.released_j`), the picture
  of Tosi et al. 2017 (arXiv:1707.06051). `eps_e` is the share the melt cannot hold dissolved (1 for CO2, N2 and argon; water is
  drawn). `PROCESSING_ENERGY_J_PER_KG` is the one calibrated number: it makes an Earth twin of the interior model give the surface
  carbon of the Earth (1e22 mol), and Venus (about the same carbon in its air) is the check.
* *Carbonate-silicate cycle.* With liquid water at the surface, silicate weathering and volcanic outgassing balance
  (Walker, Hays and Kasting 1981); the continental law of Lehmer, Catling and Krissansen-Totton 2020 (PMC7708846) with their ranges
  gives the CO2 pressure of the balance (`carbonate_pco2_bar`). The temperature depends on the pressure through the greenhouse, so
  the caller passes it as a function.
* *Ocean.* Water vapour is capped by the saturation pressure and sulfur gases are taken up by an ocean; both are the callers' (they
  have the temperature).

* *Impacts.* The late accretion delivers and erodes gas before the degassing has done much (`impact_history`, GEN.739); `impact_inventory`
  turns its balance into the same element budget and the caller adds it (`surface_air` adds the carbon and the nitrogen, GEN.819).

Escape is not subtracted: the retention rows (`d1` to `d5`) and the Jeans filter of `mixture_escape` own it.

Every number that is not read from a source is a design constant with its source tier and tweak level in
`docs/design/research/gap-filling-plan.md` (GEN.725 rule).
"""

import math

from planetgen.physics import impact_history
from planetgen.physics.mixture_base import bracket

# --- degassing --------------------------------------------------------------------------------------------------------------

EARTH_DEGASSED_SHARE = 0.36
"""float: The share of the Ortenzi CO2 inventory (log-middle 0.032 weight percent of a 4.0e24 kg mantle) that the surface carbon of the
Earth, 1e22 mol or 4.6e20 kg of CO2, is. An Earth twin of the interior model is made to reach it (`THETA_EARTH`)."""

THETA_EARTH = -math.log(1.0 - EARTH_DEGASSED_SHARE)
"""float: `Theta` of the Earth twin."""

PROCESSING_ENERGY_J_PER_KG = 7.81e6
"""float: Heat released per kilogram of mantle for one e-folding of the undegassed share: the Earth twin of the interior model (mid rheology,
4.568 Gy) releases 3.49e6 J per kilogram of its mantle, which is `THETA_EARTH` of this (calibrated, tier 5, tweak level 1; it was 9.17e6 against 4.09e6
before the mantle's ceiling became its liquidus, GEN.735, which lowers the heat the young Earth twin released)."""

DEGAS_EFFICIENCY = {"C": 1.0, "N": 1.0, "Ar": 1.0, "S": 0.5}
"""dict: The share of an element the melt cannot hold dissolved. CO2 and N2 are far less soluble than water in basalt (Tosi et al.
2017); sulfur stays partly in the melt as sulfide (a stand-in, tier 6)."""

DEGAS_EFFICIENCY_H = (0.05, 0.5, "log")
"""tuple: The same for hydrogen, drawn: water is highly soluble and the magma ocean retains much of it (Tosi et al. 2017; Bower et al.
2022; Marty 2012: up to 90 percent of the water in the mantle). Stand-in, tier 6 (a marginal educated guess)."""

# --- the carbonate-silicate cycle -------------------------------------------------------------------------------------------

WEATHERING_ALPHA = (0.1, 0.5, "uniform")
"""tuple: The exponent of the CO2 pressure in the silicate weathering rate (Lehmer et al. 2020 Table 1, geologic constraints)."""

WEATHERING_T_E_K = (10.0, 40.0, "uniform")
"""tuple: The e-folding temperature of the weathering rate, K (Lehmer et al. 2020 Table 1, laboratory and Phanerozoic constraints)."""

OUTGASSING_EXPONENT = (1.0, 2.0, "uniform")
"""tuple: The exponent of the heat flow in the volcanic outgassing flux (Lehmer et al. 2020 Table 1)."""

LAND_FRACTION_EARTH = 0.29
"""float: The land share of the Earth's surface, the reference of the weathering law."""

CO2_MODERN_BAR = 288e-6
"""float: The preindustrial CO2 pressure of the Earth (Lehmer et al. 2020)."""

TS_MODERN_K = 289.0
"""float: The preindustrial surface temperature of the Earth (Lehmer et al. 2020)."""

HEAT_FLUX_EARTH_W_M2 = 0.0614
"""float: The heat flow of the Earth twin in the interior model (a mid draw of the rheology), the reference of the outgassing flux; it
is under the Earth's measured 0.09 because the stagnant lid loses heat less freely than a mobile one (`interior_heat`)."""

MIN_PCO2_BAR = 1e-12
"""float: The floor of the balance's CO2 pressure (a world whose weathering takes all it is given)."""

LIQUID_WATER_K = (273.15, 647.0)
"""tuple: The surface temperatures between which an ocean can exist (freezing point and critical point of water)."""


def degassed_share(released_j, mantle_kg):
    """`1 - exp(-Theta)`, the share of an easily degassed element at the surface after the interior has released `released_j` J."""
    if mantle_kg <= 0.0 or released_j <= 0.0:
        return 0.0
    return -math.expm1(-released_j / (mantle_kg * PROCESSING_ENERGY_J_PER_KG))


def draw_efficiencies(constants):
    """`{element: eps}` with the hydrogen share drawn (and recorded in `constants`)."""
    return {**DEGAS_EFFICIENCY, "H": bracket(constants, "degas_efficiency_h", DEGAS_EFFICIENCY_H)}


def impact_inventory(balance, mantle_kg):
    """
    The element budget the late accretion left at the surface (GEN.739): the air and the part held in the ground, ocean or ice of
    `impact_history.veneer_balance`, in moles of atoms per 100 g of mantle like the budget the degassed share of `degassed_share` of the mantle's inventory gives, to be added to it. Gas the impactors
    brought is already at the surface, so no degassed share applies; the sinks that follow do.

    Args:
        balance (dict): The result of `impact_history.veneer_balance`.
        mantle_kg (float): The mantle's mass.

    Returns:
        dict: `{"H", "C", "N"}` moles of atoms per 100 g of mantle.
    """
    if mantle_kg <= 0.0:
        return {"H": 0.0, "C": 0.0, "N": 0.0}
    held = {species: balance["air_kg"][species] + balance["kept_kg"][species] for species in impact_history.SPECIES}
    moles = {species: kg * 1000.0 / impact_history.MOLAR_MASS[species] for species, kg in held.items()}
    per_100_g = 0.1 / mantle_kg
    return {"H": 2.0 * moles["H2O"] * per_100_g, "C": moles["CO2"] * per_100_g, "N": 2.0 * moles["N2"] * per_100_g}


def pressure_pa(air_mass_kg, gravity_m_s2, radius_m):
    """The surface pressure of an air of `air_mass_kg`: `g M / (4 pi R^2)`."""
    return gravity_m_s2 * air_mass_kg / (4.0 * math.pi * radius_m ** 2)


def air_mass_kg(pressure_pa_, gravity_m_s2, radius_m):
    """The mass of the air that makes `pressure_pa_` (the inverse of `pressure_pa`)."""
    return pressure_pa_ * 4.0 * math.pi * radius_m ** 2 / gravity_m_s2


def draw_cycle(constants):
    """`(alpha, t_e_k, m)` of the weathering law, drawn and recorded in `constants`."""
    return (bracket(constants, "weathering_alpha", WEATHERING_ALPHA), bracket(constants, "weathering_t_e_k", WEATHERING_T_E_K),
            bracket(constants, "outgassing_exponent", OUTGASSING_EXPONENT))


def carbonate_pco2_bar(supply_bar, heat_flux_w_m2, land_fraction, temperature_of, alpha, t_e_k, exponent):
    """
    The CO2 pressure at which silicate weathering takes as much as the volcanoes give (the stable root).

    The weathering rate is `rho f_land (p / p_mod)^alpha exp((Ts - Ts_mod) / T_e)` (Lehmer et al. 2020 Eq. 4) and the outgassing
    `F_mod (q / q_E)^m`; both scale with the area, so it cancels, and the Earth is the steady state at `p_mod` and `Ts_mod` by
    construction. The surface temperature `Ts(p)` rises with the pressure, so the weathering rate does too and the root is unique.

    Args:
        supply_bar (float): The CO2 the volcanoes have delivered, as a pressure: all of it if weathering cannot keep up.
        heat_flux_w_m2 (float): The body's heat flow `q`.
        land_fraction (float): The share of the surface above sea level.
        temperature_of (callable): `Ts(p_bar)`, the surface temperature with that much CO2 in the air.
        alpha, t_e_k, exponent (float): The law's draws (`draw_cycle`).

    Returns:
        float: The CO2 pressure, bar: at most `supply_bar` (a world with little CO2 has it all taken to the floor of weathering
        only when the rate there exceeds the outgassing).
    """
    if supply_bar <= MIN_PCO2_BAR:
        return max(supply_bar, 0.0)
    outgassing = max(heat_flux_w_m2 / HEAT_FLUX_EARTH_W_M2, 0.0) ** exponent
    land = max(land_fraction / LAND_FRACTION_EARTH, 1e-6)

    def excess(p_bar):
        weathering = land * (p_bar / CO2_MODERN_BAR) ** alpha * math.exp((temperature_of(p_bar) - TS_MODERN_K) / t_e_k)
        return weathering - outgassing

    if outgassing <= 0.0:
        return MIN_PCO2_BAR
    if excess(supply_bar) <= 0.0:
        return supply_bar
    low, high = math.log(MIN_PCO2_BAR), math.log(supply_bar)
    if excess(MIN_PCO2_BAR) >= 0.0:
        return MIN_PCO2_BAR
    for _ in range(80):
        middle = 0.5 * (low + high)
        if excess(math.exp(middle)) > 0.0:
            high = middle
        else:
            low = middle
        if high - low < 1e-9:
            break
    return math.exp(0.5 * (low + high))
