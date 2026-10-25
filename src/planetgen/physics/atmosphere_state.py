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
* the infrared optical depth `tau0` of the whole column, built on a wavenumber grid from the gases' bands, their far wings, the water
  vapour continuum and the collision-induced pairs (`physics/infrared_opacity.py`, GEN.703), at the surface temperature;
* the warming `G = (T_s / T_e)^4 - 1`, the variable the repo stores, from `physics/greenhouse.py` (Robinson and Catling) with
  that `tau0`, the heat-capacity ratio of the mixture, and the lapse-rate scaling.

Use:

    state = atmosphere_state.air_state(mixture, temperature_k, gravity_m_s2)
    state.density_kg_m3, state.scale_height_m, state.mean_molar_mass_kg_mol, state.column_mass_kg_m2
    warm = atmosphere_state.greenhouse_warming(mixture, gravity_m_s2, absorbed_flux_w_m2)
    warm.optical_depth, warm.surface_temperature_k, warm.factor

Limits (stated): see `physics/infrared_opacity.py`: the band and wing opacities are stand-ins until real line data are in; the
stored water vapour is taken as it is (no condensation); the closure's pressure exponent is 2 for every air.
"""

import math
from dataclasses import dataclass

from planetgen.physics import (
    constants,
    gases,
    greenhouse,
    infrared_opacity,
    mixture,
    species_thermo,
)


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
        by_gas (dict): Each absorbing gas's formula to the `tau0` of the air with only that gas's bands and wings kept.
        by_part (dict): The `tau0` of the air with only bands, far wings, the water continuum or the collision-induced pairs kept.
        lapse_scaling (float): The `alpha` of the closure (`infrared_opacity.lapse_scaling`).
    """

    optical_depth: float
    effective_temperature_k: float
    surface_temperature_k: float
    factor: float
    gamma: float
    radiative_convective_boundary: float | None
    by_gas: dict
    by_part: dict
    lapse_scaling: float


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


def heat_capacity_ratio(partial_pressures_pa, temperature_k):
    """
    `gamma = cp / (cp - R)` of the mixture at `temperature_k` (mole-fraction mean of the NASA fits' `cp / R`; the temperature is held
    inside the fits' range, so the cold end of a fit stands in below it).
    """
    cp_over_r = 0.0
    for formula, fraction in mixture.mole_fractions(partial_pressures_pa).items():
        if fraction == 0.0:
            continue
        thermo = species_thermo.gas(formula)
        low, high = thermo.range_k
        cp_over_r += fraction * thermo.cp_over_r(min(max(temperature_k, low), high))
    return cp_over_r / (cp_over_r - 1.0)


def greenhouse_warming(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2):
    """
    The warming of an air under an absorbed flux (the globally averaged W/m2 of `greenhouse.absorbed_flux_w_m2`).

    `tau0` comes from `infrared_opacity.optical_depth` at the surface temperature, which the closure itself gives, so the two are
    solved together: steps of the fixed point from the effective temperature, each a contraction of about a tenth (`tau0` depends on
    the temperature only through the Planck weights and the continuum), so it settles in four to six steps; damped steps follow if
    it does not. The heat-capacity ratio is taken at the effective temperature (the airless grey temperature, the lowest the surface
    can be); an air with no absorber gives `tau0` zero and `G` zero.

    Raises:
        ValueError: The mixture is empty, or the flux or gravity is not positive.
    """
    _check_positive(absorbed_flux_w_m2=absorbed_flux_w_m2, gravity_m_s2=gravity_m_s2)
    effective = greenhouse.effective_temperature_k(absorbed_flux_w_m2)
    gamma = heat_capacity_ratio(partial_pressures_pa, effective)
    temperature = effective
    for step in range(infrared_opacity.ITERATIONS):
        opacity = infrared_opacity.optical_depth(partial_pressures_pa, temperature, gravity_m_s2)
        lapse = infrared_opacity.lapse_scaling(opacity)
        surface, boundary = greenhouse.surface_temperature_k(opacity.tau0, absorbed_flux_w_m2, gamma, lapse,
                                                             infrared_opacity.PRESSURE_EXPONENT)
        settled = abs(surface - temperature) < infrared_opacity.TOLERANCE_K
        if settled:
            break
        temperature = surface if step < infrared_opacity.UNDAMPED_STEPS else 0.5 * (temperature + surface)
    return Warming(optical_depth=opacity.tau0, effective_temperature_k=effective, surface_temperature_k=surface,
                   factor=greenhouse.greenhouse_factor(surface, effective), gamma=gamma, radiative_convective_boundary=boundary,
                   by_gas=opacity.by_gas, by_part=opacity.by_part, lapse_scaling=lapse)


def _check_positive(**values):
    for name, value in values.items():
        if not (value > 0.0 and math.isfinite(value)):
            raise ValueError(f"{name} must be a positive number, not {value!r}")
