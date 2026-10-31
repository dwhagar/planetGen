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
* the surface temperature of the air under an absorbed flux and the warming `G = (T_s / T_e)^4 - 1`, the variable the repo stores,
  read from the table fitted to the radiative-convective model of the ExoMol and HITRAN data (`physics/greenhouse_table.py`, GEN.850).

Use:

    state = atmosphere_state.air_state(mixture, temperature_k, gravity_m_s2)
    state.density_kg_m3, state.scale_height_m, state.mean_molar_mass_kg_mol, state.column_mass_kg_m2
    warm = atmosphere_state.greenhouse_warming(mixture, gravity_m_s2, absorbed_flux_w_m2)
    warm.surface_temperature_k, warm.factor

Limits (stated): see `physics/greenhouse_table.py` and `docs/design/composition/greenhouse-table.md`: the table is a fit to a clear-sky
one-dimensional model, held at its edge outside the airs it was fitted on.
"""

import math
from dataclasses import dataclass

from planetgen.physics import (
    constants,
    greenhouse_table,
    mixture,
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


def greenhouse_warming(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share=0.0, aloft_pressure_pa=None):
    """
    The warming of an air under an absorbed flux (the globally averaged W/m2 of `greenhouse.absorbed_flux_w_m2`), read from the table
    fitted to the radiative-convective model (`greenhouse_table.warming`, GEN.850). `aloft_share` of the flux is taken by a haze at
    `aloft_pressure_pa` instead of the ground (GEN.851).

    Raises:
        ValueError: The mixture is empty, or the flux or gravity is not positive.
    """
    return greenhouse_table.warming(partial_pressures_pa, gravity_m_s2, absorbed_flux_w_m2, aloft_share, aloft_pressure_pa)


def _check_positive(**values):
    for name, value in values.items():
        if not (value > 0.0 and math.isfinite(value)):
            raise ValueError(f"{name} must be a positive number, not {value!r}")
