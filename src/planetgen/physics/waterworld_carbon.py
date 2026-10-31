# planetgen/physics/waterworld_carbon.py

"""
The carbon sink of a world with no land: ocean carbonate chemistry under a high-pressure-ice seal, and seafloor weathering (GEN.814)
=====================================================================================================================================

Silicate weathering (`surface_pressure.carbonate_pco2_bar`) needs land. A water-rich world has none: the ocean is deeper than the
tallest mountain. What takes its carbon dioxide back out of the air depends on how deep the ocean is.

* **A sealed world** (seafloor pressure at or above `SEAL_PRESSURE_PA`, 1 to 3 GPa). The weight of the water stops the mantle melting
  (Kite et al. 2009; Kite and Ford 2018 appendix A: the solidus rises with pressure, the melt column shortens, the mantle under the lid
  is colder) and high-pressure ice (ice VI and VII) sits on the seafloor, so the ocean and the mantle no longer exchange carbon. About
  `SEAL_TIME_GY` after the last giant impact the atmosphere, the ocean and the shallow crust are a closed system for carbon (Kite and
  Ford 2018 section 2 item 4). Its pCO2 is set by ocean chemistry alone: the carbon inventory is shared between the air, the dissolved
  carbon dioxide, bicarbonate and carbonate, and calcite, and the cations the water leached from the crust early on (sodium and calcium,
  `leached_cations_mol_m2`) fix the alkalinity. Calcium that pairs one carbon with one calcium takes the carbon out; sodium raises the pH
  and so the share of the dissolved carbon that is bicarbonate (`ocean_pco2`). The pCO2 follows the temperature (the carbon dioxide
  leaves a warmer ocean), so a caller solves it with its greenhouse for the temperature, as Kite and Ford do.
* **An unsealed world with no land** (seafloor pressure under the seal, no land above the sea). Volcanoes work and seafloor basalt
  dissolves and precipitates carbonate in the pore space of the crust. The law is that of Krissansen-Totton and Catling (2017):
  the dissolution rate is Arrhenius in the pore-water temperature, a power of the hydrogen ion molality (so of the pCO2) and a power of
  the spreading, which here is the heat flow like the outgassing. Its modern size is a small share of the modern outgassing
  (`seafloor_pco2_bar`), so the balance is set by the temperature far more than by the pressure (the pCO2 enters to the power 0.2, the Arrhenius law with 75 kJ/mol): at the modern heat flow a 300 K world balances at about half a bar, a 288 K one needs hundreds of bar.

What this does not do (stated): the ocean is an ideal solution (no activity coefficients: seawater's first dissociation constant is
about three times the fresh water's, so a real salty ocean holds more bicarbonate at a given pCO2), the constants of the carbonate system
are the fresh-water fits of Plummer and Busenberg (1982) held to 273 to 373 K (warmer oceans use their 373 K value); only calcite
precipitates (no magnesium or iron carbonate, which Kite and Ford find silicates take instead; no dolomite, no hydrous carbonates,
no clathrates, no brucite or portlandite); the cations are a fixed leach of a basalt crust of a chosen thickness, not a reaction-transport
model; the temperature is the surface ocean's (the gradient with depth is ignored); the seal is a pressure threshold, not a solidus
calculation, and the ocean is taken hydrostatic at the surface gravity; the carbon given to the closed system is whatever the caller
degassed by `SEAL_TIME_GY`; no carbon comes back from the mantle afterwards; the pCO2 of a 10 to 50 percent water world is only as good as
that inventory (Kite and Ford's carbon is 1e-5 to 1e-3 of the planet mass); the seafloor balance has one pore-space law and not the
whole of Krissansen-Totton and Catling's inverse model (no Michaelis-Menten variant, no pelagic fraction, no carbonate dissolution
feedback). The nitrogen of the air is not here: the caller adds it.

Use:

    sealed = is_sealed(seafloor_pressure_pa(ocean_column_kg_m2, gravity), seal_pa)
    cations = leached_cations_mol_m2(altered_km)                             # (sodium, calcium), mol per m2
    carbon = ocean_pco2(carbon_mol_m2, ts_k, ocean_column_kg_m2, gravity, sodium_mol_m2=cations[0], calcium_mol_m2=cations[1])
    carbon.pco2_bar, carbon.ph, carbon.precipitated_mol_m2

Every number that is not read from a source is a design constant with its source tier and tweak level in
`docs/design/research/gap-filling-plan.md` (GEN.725 rule).
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, elements
from planetgen.physics.mixture_base import bracket
from planetgen.physics.sulfur_balance import CRUST_CALCIUM_WT
from planetgen.physics.surface_pressure import (CO2_MODERN_BAR, HEAT_FLUX_EARTH_W_M2, MIN_PCO2_BAR, TS_MODERN_K)

# --- the seal ---------------------------------------------------------------------------------------------------------------

SEAL_PRESSURE_PA = (1.0e9, 3.0e9, "uniform")
"""tuple: The seafloor pressure from which the mantle no longer melts enough to exchange carbon with the ocean, Pa: Kite and Ford 2018 impose
an abrupt 1 GPa cutoff and note it is smaller for stagnant-lid convection (about 1 to 3 GPa; Kite et al. 2009; Noack et al. 2016 put it
nearer 3 GPa). A bracket of the sources that differ (tier 3, tweak level 0)."""

SEAL_TIME_GY = 0.1
"""float: How long after the last giant impact a sealed world still exchanges carbon with its mantle, Gy: Kite and Ford 2018 start the closed
system at 0.1 Gy ("shuts down within O(1e8) yr"; READ, tier 1, tweak level 0). The caller takes the carbon degassed up to this age."""

ROCK_DENSITY_KG_M3 = 3000.0
"""float: The density of the crust, kg/m3 (Kite and Ford 2018 below their Eq. 4; READ, tier 1)."""

ALTERED_CRUST_KM = (0.1, 50.0, "log")
"""tuple: The thickness of crust the early ocean leaches completely, km of rock: Kite and Ford 2018 section 3.1 (2c) put the reasonable range at
0.1 km (Earth-like alteration of the top of the extrusives) to 50 km (pervasive alteration of a thick stagnant-lid extrusive pile); this
already takes the intrusive to extrusive ratio into account. A bracket of their argument (tier 3, tweak level 0)."""

CRUST_SODIUM_WT = 0.028
"""float: The sodium mass fraction of the crust that is leached (Gale et al. 2013 through Kite and Ford 2018 Eq. 4, Na0; READ, tier 1)."""


@dataclass(frozen=True)
class WaterWorld:
    """
    The draws of a world with no land.

    Attributes:
        seal_pa (float): The seafloor pressure of the seal, Pa.
        altered_km (float): The leached crust thickness, km.
    """

    seal_pa: float
    altered_km: float


def draw_waterworld(constants_):
    """The seal pressure and the leached crust thickness, drawn and recorded in `constants_` (`waterworld_seal_pa`, `waterworld_altered_km`)."""
    return WaterWorld(seal_pa=bracket(constants_, "waterworld_seal_pa", SEAL_PRESSURE_PA),
                      altered_km=bracket(constants_, "waterworld_altered_km", ALTERED_CRUST_KM))


def seafloor_pressure_pa(ocean_column_kg_m2, gravity_m_s2):
    """The pressure at the seafloor, `g m` of an ocean of `m` kg per square metre (hydrostatic, the surface gravity throughout)."""
    return max(ocean_column_kg_m2, 0.0) * gravity_m_s2


def is_sealed(seafloor_pa, seal_pa):
    """Whether the weight of the ocean has stopped the carbon exchange with the mantle."""
    return seafloor_pa >= seal_pa


OCEAN_CHEMISTRY, CONTINENTAL_WEATHERING, SEAFLOOR_WEATHERING = "ocean-chemistry", "continental-weathering", "seafloor-weathering"
"""str: The three carbon sinks a world with an ocean can have."""


def sink_kind(sealed, land_fraction):
    """
    Which sink sets the carbon dioxide of a world with an ocean: `OCEAN_CHEMISTRY` if the seafloor is sealed (whatever land there is, the
    mantle no longer feeds the air), else `CONTINENTAL_WEATHERING` (`surface_pressure.carbonate_pco2_bar`) if there is land, else
    `SEAFLOOR_WEATHERING` (`seafloor_pco2_bar`).
    """
    if sealed:
        return OCEAN_CHEMISTRY
    return CONTINENTAL_WEATHERING if land_fraction > 0.0 else SEAFLOOR_WEATHERING


def leached_cations_mol_m2(altered_km):
    """
    The cations an early ocean leaches from `altered_km` of crust, mol per square metre: `(sodium, calcium)`.

    Kite and Ford 2018 Eq. 4: `[Na] = Na0 / M_Na * z rho / (z_oc rho_w)`, here per square metre (the ocean's share is the division by its
    mass in `ocean_pco2`). Calcium is the crust's calcium (`sulfur_balance.CRUST_CALCIUM_WT`, one value for both sinks).
    """
    mass = max(altered_km, 0.0) * 1.0e3 * ROCK_DENSITY_KG_M3
    return (mass * CRUST_SODIUM_WT / elements.molar_mass("Na"), mass * CRUST_CALCIUM_WT / elements.molar_mass("Ca"))


# --- the carbonate system ---------------------------------------------------------------------------------------------------

TEMPERATURE_RANGE_K = (273.15, 373.15)
"""tuple: The temperatures the carbonate constants are used at (Plummer and Busenberg 1982 give 0 to 90 C; hotter oceans use the 373 K values)."""

BAR_PA = 1.0e5
"""float: One bar, Pa."""


CARBONATE_FITS = {
    "k0": (108.3865, 0.01985076, -6919.53, -40.45154, 669365.0),
    "k1": (-356.3094, -0.06091964, 21834.37, 126.8339, -1684915.0),
    "k2": (-107.8871, -0.03252849, 5151.79, 38.92561, -563713.9),
    "ksp": (-171.9065, -0.077993, 2839.319, 71.595, 0.0),
    "kw": (6.0875, -0.01706, -4470.99, 0.0, 0.0),
}
"""dict: `log10 K = a + b T + c / T + d log10 T + e / T^2` for the Henry constant of carbon dioxide (mol/(kg bar)), the first and second
dissociation constants of carbonic acid, the solubility product of calcite and the ion product of water (mol and kg): the fresh-water fits
of Plummer and Busenberg (1982) for the carbonate system and Harned and Hamer (1933) for water, valid 0 to 90 C. They reproduce at 298.15 K
the known 0.0340, 4.45e-7 (pK 6.35), 4.69e-11 (pK 10.33), 3.3e-9 (pK 8.48) and 1.0e-14 (READ from memory, not re-checked against the papers;
tier 3, tweak level 0)."""


def _log10_constants(temperature_k):
    """`(log K0, log K1, log K2, log Ksp(calcite), log Kw)` at `temperature_k`, clamped to `TEMPERATURE_RANGE_K` (`CARBONATE_FITS`)."""
    t = min(max(temperature_k, TEMPERATURE_RANGE_K[0]), TEMPERATURE_RANGE_K[1])
    lt = math.log10(t)
    return tuple(a + b * t + c / t + d * lt + e / t ** 2 for a, b, c, d, e in
                 (CARBONATE_FITS[name] for name in ("k0", "k1", "k2", "ksp", "kw")))


@dataclass(frozen=True)
class OceanCarbon:
    """
    The carbon of a closed atmosphere, ocean and crust.

    Attributes:
        pco2_bar (float): The carbon dioxide partial pressure in the air, bar.
        ph (float): The pH of the ocean (ideal solution).
        atmosphere_mol_m2 (float): Carbon in the air, mol per square metre.
        dissolved_mol_m2 (float): Dissolved carbon dioxide, bicarbonate and carbonate, mol per square metre.
        precipitated_mol_m2 (float): Calcite, mol per square metre.
        saturated (bool): Whether calcite stands in equilibrium with the ocean (otherwise the calcium is all dissolved).
    """

    pco2_bar: float
    ph: float
    atmosphere_mol_m2: float
    dissolved_mol_m2: float
    precipitated_mol_m2: float
    saturated: bool


def _chemistry(pco2_bar, constants_, sodium, calcium):
    """The ocean at a given pCO2: `(H, co2, hco3, co3, calcium dissolved)`, mol per kg of water, by bisection on the charge balance."""
    log_k0, log_k1, log_k2, log_ksp, log_kw = constants_
    co2 = 10.0 ** log_k0 * pco2_bar
    k1, k2, ksp, kw = 10.0 ** log_k1, 10.0 ** log_k2, 10.0 ** log_ksp, 10.0 ** log_kw

    def state(ph):
        h = 10.0 ** -ph
        hco3 = k1 * co2 / h
        co3 = k2 * hco3 / h
        dissolved = calcium if calcium * co3 <= ksp else ksp / co3
        return h, hco3, co3, dissolved

    def residual(ph):
        h, hco3, co3, dissolved = state(ph)
        return sodium + 2.0 * dissolved + h - kw / h - hco3 - 2.0 * co3

    low, high = -2.0, 16.0
    for _ in range(100):
        middle = 0.5 * (low + high)
        if residual(middle) > 0.0:
            low = middle
        else:
            high = middle
    ph = 0.5 * (low + high)
    h, hco3, co3, dissolved = state(ph)
    return h, co2, hco3, co3, dissolved


def ocean_pco2(carbon_mol_m2, temperature_k, ocean_column_kg_m2, gravity_m_s2, *, sodium_mol_m2, calcium_mol_m2,
               air_molar_mass_kg_mol=0.04401):
    """
    The carbon dioxide pressure of a closed atmosphere, ocean and crust that holds `carbon_mol_m2` of carbon in all.

    Args:
        carbon_mol_m2 (float): The carbon of the closed system, mol per square metre of the surface (air, dissolved and calcite together).
        temperature_k (float): The surface ocean's temperature.
        ocean_column_kg_m2 (float): The ocean's mass per square metre.
        gravity_m_s2 (float): The surface gravity (the air's weight per mole follows from it).
        sodium_mol_m2, calcium_mol_m2 (float): The leached cations (`leached_cations_mol_m2`).
        air_molar_mass_kg_mol (float): The mean molar mass of the air, kg/mol, which sets the moles of carbon dioxide per pascal
            (`p / (g M)`); 44 g/mol for a carbon dioxide air as Kite and Ford take it.

    Returns:
        OceanCarbon: The partition; the three carbon reservoirs add up to `carbon_mol_m2`.

    Raises:
        ValueError: The ocean or the gravity is not positive.
    """
    if ocean_column_kg_m2 <= 0.0 or gravity_m_s2 <= 0.0 or air_molar_mass_kg_mol <= 0.0:
        raise ValueError("the ocean column, the gravity and the air's molar mass must be positive")
    carbon = max(carbon_mol_m2, 0.0)
    ocean = ocean_column_kg_m2
    sodium, calcium = max(sodium_mol_m2, 0.0) / ocean, max(calcium_mol_m2, 0.0) / ocean
    constants_ = _log10_constants(temperature_k)
    per_pascal = 1.0 / (gravity_m_s2 * air_molar_mass_kg_mol)            # mol of carbon dioxide above one m2 per pascal of its pressure

    def totals(pco2_bar):
        h, co2, hco3, co3, dissolved = _chemistry(pco2_bar, constants_, sodium, calcium)
        atmosphere = pco2_bar * BAR_PA * per_pascal
        in_water = (co2 + hco3 + co3) * ocean
        calcite = (calcium - dissolved) * ocean
        return atmosphere, in_water, calcite, h, dissolved

    def inventory(pco2_bar):
        atmosphere, in_water, calcite, _, _ = totals(pco2_bar)
        return atmosphere + in_water + calcite

    if carbon <= 0.0 or inventory(MIN_PCO2_BAR) >= carbon:
        pco2 = MIN_PCO2_BAR if carbon > 0.0 else 0.0
    else:
        low, high = math.log(MIN_PCO2_BAR), math.log(1.0e4)
        while inventory(math.exp(high)) < carbon and high < 50.0:
            high += 5.0
        for _ in range(200):
            middle = 0.5 * (low + high)
            if inventory(math.exp(middle)) < carbon:
                low = middle
            else:
                high = middle
            if high - low < 1e-12:
                break
        pco2 = math.exp(0.5 * (low + high))
    atmosphere, in_water, calcite, h, dissolved = totals(pco2) if pco2 > 0.0 else (0.0, 0.0, 0.0, 1e-7, calcium)
    total = atmosphere + in_water + calcite
    scale = carbon / total if total > 0.0 else 0.0                       # the bisection's last digits, so the reservoirs add up exactly
    return OceanCarbon(pco2_bar=pco2, ph=-math.log10(h), atmosphere_mol_m2=atmosphere * scale, dissolved_mol_m2=in_water * scale,
                       precipitated_mol_m2=calcite * scale, saturated=calcite > 0.0)


# --- seafloor weathering ----------------------------------------------------------------------------------------------------

SEAFLOOR_SHARE_OF_OUTGASSING = 0.45 / 6.6
"""float: The modern seafloor carbonate precipitation over the modern outgassing, `0.45 / 6.6` Tmol C per year (Krissansen-Totton and Catling
2017 Table 2 and their outgassing posterior; READ, tier 1, tweak level 0). It is the size of the modern sink against the volcanoes."""

SEAFLOOR_ACTIVATION_KJ_MOL = (40.0, 110.0, "uniform")
"""tuple: The activation energy of basalt dissolution, kJ/mol: the prior of Krissansen-Totton and Catling 2017 (posterior 75, 53 to 97; READ, tier 1)."""

SEAFLOOR_PH_EXPONENT = (0.0, 0.5, "uniform")
"""tuple: The exponent of the hydrogen ion molality, `gamma`: prior of Krissansen-Totton and Catling 2017 (posterior 0.27, 0.11 to 0.43; READ, tier 1)."""

SEAFLOOR_SPREADING_EXPONENT = (0.0, 1.0, "uniform")
"""tuple: The exponent of the spreading (here the heat flow) in the dissolution rate, `beta`: prior of Krissansen-Totton and Catling 2017
(posterior 0.47, 0.15 to 0.82; READ, tier 1)."""

PH_PER_PCO2_EXPONENT = 1.34
"""float: The hydrogen ion molality of seawater rises as the pCO2 to the 1/1.34 (the Brady and Gislason 1997 approximation of Krissansen-Totton and
Catling 2017 Eq. 1, as restated: `beta_S = gamma / 1.34`; second-hand, tier 2)."""

DEEP_OCEAN_GRADIENT = 1.02
"""float: The deep-ocean temperature rises 1.02 K per K of the surface's (Krissansen-Totton and Catling 2017 Eq. 12, best fit; READ, tier 1)."""

DEEP_OCEAN_REFERENCE_K = 274.037
"""float: The deep ocean is 274.037 K when the surface is 285 K (the same equation; the line passes through it)."""

PORE_OVER_DEEP_K = 9.0
"""float: The pore water of the basalt is 9 K warmer than the deep ocean (Krissansen-Totton and Catling 2017; READ, tier 1)."""

REFERENCE_SURFACE_K = 285.0
"""float: The surface temperature the deep-ocean relation is anchored at, K."""


@dataclass(frozen=True)
class SeafloorLaw:
    """
    The draws of the seafloor weathering law.

    Attributes:
        activation_j_mol (float): `E_bas`, J/mol.
        ph_exponent (float): `gamma`.
        spreading_exponent (float): `beta`.
    """

    activation_j_mol: float
    ph_exponent: float
    spreading_exponent: float


def draw_seafloor_law(constants_):
    """The seafloor law's three draws, recorded in `constants_` (`seafloor_activation_kj_mol`, `seafloor_ph_exponent`, `seafloor_spreading_exponent`)."""
    return SeafloorLaw(activation_j_mol=1.0e3 * bracket(constants_, "seafloor_activation_kj_mol", SEAFLOOR_ACTIVATION_KJ_MOL),
                       ph_exponent=bracket(constants_, "seafloor_ph_exponent", SEAFLOOR_PH_EXPONENT),
                       spreading_exponent=bracket(constants_, "seafloor_spreading_exponent", SEAFLOOR_SPREADING_EXPONENT))


def pore_temperature_k(surface_temperature_k):
    """The temperature of the basalt's pore water under a surface ocean at `surface_temperature_k`, K."""
    deep = DEEP_OCEAN_GRADIENT * surface_temperature_k + (DEEP_OCEAN_REFERENCE_K - DEEP_OCEAN_GRADIENT * REFERENCE_SURFACE_K)
    return deep + PORE_OVER_DEEP_K


def seafloor_rate_over_modern(pco2_bar, surface_temperature_k, heat_flux_w_m2, law):
    """
    The seafloor carbonate precipitation per unit area over the Earth's modern one.

    `(pCO2 / pCO2_modern)^(gamma / 1.34) exp[E / R (1 / T_pore_modern - 1 / T_pore)] (q / q_modern)^beta`.
    """
    modern_pore = pore_temperature_k(TS_MODERN_K)
    pore = pore_temperature_k(surface_temperature_k)
    if pore <= 0.0:  # a surface this cold (a far, starved body) has no pore water to react: the extrapolated gradient has run out
        return 0.0
    return ((max(pco2_bar, MIN_PCO2_BAR) / CO2_MODERN_BAR) ** (law.ph_exponent / PH_PER_PCO2_EXPONENT)
            * math.exp(law.activation_j_mol / constants.R * (1.0 / modern_pore - 1.0 / pore))
            * (max(heat_flux_w_m2, 0.0) / HEAT_FLUX_EARTH_W_M2) ** law.spreading_exponent)


def seafloor_pco2_bar(supply_bar, heat_flux_w_m2, temperature_of, law, outgassing_exponent):
    """
    The CO2 pressure at which seafloor weathering takes as much as the volcanoes give (the stable root), for a world with no land and a
    seafloor under the seal pressure.

    The modern seafloor precipitation is `SEAFLOOR_SHARE_OF_OUTGASSING` of the modern outgassing, and the outgassing scales as
    `(q / q_modern)^m` (`surface_pressure.OUTGASSING_EXPONENT`), so the balance is
    `SEAFLOOR_SHARE_OF_OUTGASSING * seafloor_rate_over_modern = (q / q_modern)^m`.

    Args:
        supply_bar (float): The CO2 the volcanoes have delivered, as a pressure: all of it if the seafloor cannot keep up.
        heat_flux_w_m2 (float): The body's heat flow.
        temperature_of (callable): `Ts(p_bar)`, the surface temperature with that much CO2 in the air.
        law (SeafloorLaw): `draw_seafloor_law`.
        outgassing_exponent (float): `m` of the outgassing law.

    Returns:
        float: The CO2 pressure, bar, at most `supply_bar`.
    """
    if supply_bar <= MIN_PCO2_BAR:
        return max(supply_bar, 0.0)
    outgassing = max(heat_flux_w_m2 / HEAT_FLUX_EARTH_W_M2, 0.0) ** outgassing_exponent
    if outgassing <= 0.0:
        return MIN_PCO2_BAR

    def excess(p_bar):
        return SEAFLOOR_SHARE_OF_OUTGASSING * seafloor_rate_over_modern(p_bar, temperature_of(p_bar), heat_flux_w_m2, law) - outgassing

    if excess(supply_bar) <= 0.0:
        return supply_bar
    if excess(MIN_PCO2_BAR) >= 0.0:
        return MIN_PCO2_BAR
    low, high = math.log(MIN_PCO2_BAR), math.log(supply_bar)
    for _ in range(80):
        middle = 0.5 * (low + high)
        if excess(math.exp(middle)) > 0.0:
            high = middle
        else:
            low = middle
        if high - low < 1e-9:
            break
    return math.exp(0.5 * (low + high))
