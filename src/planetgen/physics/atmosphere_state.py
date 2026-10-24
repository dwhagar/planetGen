# planetgen/physics/atmosphere_state.py

"""
The air's density, scale height and greenhouse warming from its mixture (GEN.423,
docs/design/composition/atmosphere-composition-parts/03-infrared-greenhouse-thermal.md sections 1, 5.4 and 6.2)
======================================================================================================================

Until now the drawn `atm_molar_density`, the atmosphere density and the greenhouse multiplier of `physics/planets.py` were
three separate random draws that never saw the gases. Here every one of them is a function of the mixture (the dict of
`physics/mixture.py`: gas-table formula to partial pressure, Pa), the surface temperature and the gravity:

* mean molar mass `M = sum x_i M_i` (`mixture.mean_molar_mass_kg_mol`);
* surface density `rho = P M / (R T)` (ideal gas);
* scale height `H = R T / (M g)`, so the surface pressure is `rho g H` as the barometric formula says;
* column mass `P / g`, the mass above one square metre;
* the infrared optical depth `tau0` of the whole column, the sum over the gases of `kappa_i m_i`, with `m_i = x_i M_i P / (M g)`
  the mass of gas `i` above one square metre (a well-mixed hydrostatic column) and `kappa_i` its grey mass absorption
  coefficient, m2/kg. Doubling a gas's partial pressure doubles its share of `tau0`; a heavier or weaker world (more `g`) has
  less column per pascal.
* the warming `G = (T_s / T_e)^4 - 1`, the variable the repo stores, from `physics/greenhouse.py` (Robinson and Catling) with
  that `tau0`, the heat-capacity ratio of the mixture, and the lapse-rate scaling.

Use:

    state = atmosphere_state.air_state(mixture, temperature_k, gravity_m_s2)
    state.density_kg_m3, state.scale_height_m, state.mean_molar_mass_kg_mol, state.column_mass_kg_m2
    warm = atmosphere_state.greenhouse_warming(mixture, gravity_m_s2, absorbed_flux_w_m2)
    warm.optical_depth, warm.surface_temperature_k, warm.factor

Limits (stated, lane E's GEN.703 replaces the opacity part): `kappa_i` are stand-ins read from the research's statements of
how thick a pure gas is (`INFRARED_MASS_ABSORPTION_M2_KG`); there is no collision-induced absorption, so a dense N2 or H2 air
shows none of its pressure-induced warming; no condensation (the stored water vapour is taken as it is); a single grey band
with the pressure exponent `n` of 1 (opacity linear in the column, the line-dominated regime).
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, gases, greenhouse, mixture, species_thermo

INFRARED_MASS_ABSORPTION_M2_KG = {
    "H2O": 1.0e-2,
    "CO2": 3.0e-4,
    "CO": 1.0e-5,
}
"""
dict: The grey infrared mass absorption coefficient of a gas, m2 per kg of that gas. STAND-IN (labelled, replaced by GEN.703's
band opacities): read from Koll 2019 (arXiv:1907.13145, text of his Fig. 4, part 03 section 6.2), who finds a pure H2O air
optically thick above about 1e-2 bar (kappa about 1e-2 m2/kg at g near 10), a pure CO2 air of order unity at 1 bar (taken as
3 optical depths at 1 bar, kappa 3e-4), and CO thin even at 1 bar (taken as 1e-5). The values are not tuned to Earth or Venus.
"""
GENERIC_BAND_M2_KG = 1.0e-3
"""float: STAND-IN, a gas the table flags as an infrared absorber (`ir == "ir"`) with no entry above: between CO2 and H2O."""
WEAK_BAND_M2_KG = 1.0e-5
"""float: STAND-IN, a gas the table flags `"weak"`: as thin as CO."""
LAPSE_SCALING = 0.8
"""float: STAND-IN, the lapse-rate scaling alpha of Robinson and Catling (0.6 to 0.9: Earth 0.6, Venus 0.8, Titan 0.77, Jupiter 0.85)."""
PRESSURE_EXPONENT = 1.0
"""float: The exponent n of the closure; 1 because the opacity here is linear in the column (the line-dominated regime)."""


@dataclass(frozen=True)
class AirState:
    """
    The air of one body at one temperature and gravity.

    Attributes:
        total_pressure_pa (float): Sum of the partial pressures.
        mean_molar_mass_kg_mol (float): `sum x_i M_i`.
        density_kg_m3 (float): `P M / (R T)` at the surface.
        scale_height_m (float): `R T / (M g)`.
        column_mass_kg_m2 (float): `P / g`.
    """

    total_pressure_pa: float
    mean_molar_mass_kg_mol: float
    density_kg_m3: float
    scale_height_m: float
    column_mass_kg_m2: float


@dataclass(frozen=True)
class Warming:
    """
    The greenhouse warming of one air.

    Attributes:
        optical_depth (float): `tau0`, the grey infrared optical depth of the column.
        effective_temperature_k (float): `T_e`, the airless grey temperature of the absorbed flux.
        surface_temperature_k (float): `T_s` of the closure.
        factor (float): `G = (T_s / T_e)^4 - 1`, the repo's greenhouse variable.
        gamma (float): The heat-capacity ratio of the mixture the closure used.
        radiative_convective_boundary (float or None): `tau_rc` (`None` for a purely radiative atmosphere).
        by_gas (dict): Each absorbing gas's formula to its share of `tau0`.
    """

    optical_depth: float
    effective_temperature_k: float
    surface_temperature_k: float
    factor: float
    gamma: float
    radiative_convective_boundary: float | None
    by_gas: dict


def air_state(partial_pressures_pa, temperature_k, gravity_m_s2):
    """
    The mean molar mass, density, scale height and column mass of a mixture.

    Raises:
        ValueError: The mixture is empty, or the temperature or gravity is not positive.
    """
    _check_positive(temperature_k=temperature_k, gravity_m_s2=gravity_m_s2)
    total = mixture.total_pressure_pa(partial_pressures_pa)
    molar_mass = mixture.mean_molar_mass_kg_mol(partial_pressures_pa)
    return AirState(total_pressure_pa=total, mean_molar_mass_kg_mol=molar_mass,
                    density_kg_m3=total * molar_mass / (constants.R * temperature_k),
                    scale_height_m=constants.R * temperature_k / (molar_mass * gravity_m_s2),
                    column_mass_kg_m2=total / gravity_m_s2)


def gas_column_mass_kg_m2(partial_pressures_pa, gravity_m_s2):
    """Each gas's mass above one square metre of a well-mixed hydrostatic column, `x_i M_i P / (M g)`, kg/m2 (they add up to `P / g`)."""
    _check_positive(gravity_m_s2=gravity_m_s2)
    total = mixture.total_pressure_pa(partial_pressures_pa)
    molar_mass = mixture.mean_molar_mass_kg_mol(partial_pressures_pa)
    return {formula: (pressure / total) * gases.gas(formula).molar_mass * total / (molar_mass * gravity_m_s2)
            for formula, pressure in partial_pressures_pa.items()}


def mass_absorption_m2_kg(formula):
    """The grey infrared mass absorption coefficient of a gas, m2/kg: its entry in `INFRARED_MASS_ABSORPTION_M2_KG`, else the
    table's flag (`"ir"` generic, `"weak"` thin, `"cia"` and `"none"` zero: no line absorption)."""
    if formula in INFRARED_MASS_ABSORPTION_M2_KG:
        return INFRARED_MASS_ABSORPTION_M2_KG[formula]
    flag = gases.gas(formula).ir
    if flag == "ir":
        return GENERIC_BAND_M2_KG
    if flag == "weak":
        return WEAK_BAND_M2_KG
    return 0.0


def optical_depth_by_gas(partial_pressures_pa, gravity_m_s2):
    """Each gas's share of the surface infrared optical depth, `kappa_i m_i`; gases with no line absorption are left out."""
    shares = {}
    for formula, column in gas_column_mass_kg_m2(partial_pressures_pa, gravity_m_s2).items():
        kappa = mass_absorption_m2_kg(formula)
        if kappa > 0.0:
            shares[formula] = kappa * column
    return shares


def infrared_optical_depth(partial_pressures_pa, gravity_m_s2):
    """`tau0`, the grey infrared optical depth of the whole column: the sum of `optical_depth_by_gas`."""
    return sum(optical_depth_by_gas(partial_pressures_pa, gravity_m_s2).values())


def heat_capacity_ratio(partial_pressures_pa, temperature_k):
    """
    `gamma = cp / (cp - R)` of the mixture at `temperature_k` (mole-fraction mean of the NASA fits' `cp / R`; the temperature is held
    inside the fits' range, so the cold end of a fit stands in below it).
    """
    cp_over_r = 0.0
    for formula, fraction in mixture.mole_fractions(partial_pressures_pa).items():
        thermo = species_thermo.species(formula)
        low, high = thermo.range_k
        cp_over_r += fraction * thermo.cp_over_r(min(max(temperature_k, low), high))
    return cp_over_r / (cp_over_r - 1.0)


def greenhouse_warming(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2):
    """
    The warming of an air under an absorbed flux (the globally averaged W/m2 of `greenhouse.absorbed_flux_w_m2`).

    The heat-capacity ratio is taken at the effective temperature (the airless grey temperature), the lowest the surface can be;
    an air with no line absorber gives `tau0` zero and `G` zero.

    Raises:
        ValueError: The mixture is empty, or the flux or gravity is not positive.
    """
    _check_positive(absorbed_flux_w_m2=absorbed_flux_w_m2, gravity_m_s2=gravity_m_s2)
    shares = optical_depth_by_gas(partial_pressures_pa, gravity_m_s2)
    tau0 = sum(shares.values())
    effective = greenhouse.effective_temperature_k(absorbed_flux_w_m2)
    gamma = heat_capacity_ratio(partial_pressures_pa, effective)
    surface, boundary = greenhouse.surface_temperature_k(tau0, absorbed_flux_w_m2, gamma, LAPSE_SCALING, PRESSURE_EXPONENT)
    return Warming(optical_depth=tau0, effective_temperature_k=effective, surface_temperature_k=surface,
                   factor=greenhouse.greenhouse_factor(surface, effective), gamma=gamma, radiative_convective_boundary=boundary,
                   by_gas=shares)


def _check_positive(**values):
    for name, value in values.items():
        if not (value > 0.0 and math.isfinite(value)):
            raise ValueError(f"{name} must be a positive number, not {value!r}")
