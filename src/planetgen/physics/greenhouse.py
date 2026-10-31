# planetgen/physics/greenhouse.py

"""
The greenhouse variable and the flux behind it (GEN.443, docs/design/composition/radiation-model.md section 4)
=========================================================================================================

The repo's greenhouse variable is `G = (T_s / T_e)^4 - 1`, with `T_e` the effective temperature of the absorbed flux. The surface temperature of an air
comes from `greenhouse_table` (GEN.850), the fit to the radiative-convective model of the ExoMol and HITRAN data; this module holds the three
definitions around it:

* `effective_temperature_k(absorbed_flux)` is `T_e`;
* `absorbed_flux_w_m2(stellar_flux, bond_albedo, redistribution)` is the globally averaged flux a planet absorbs, `(1 - A) S f`, `f` 1/4 for a planet
  that spreads the heat (a tidally locked world with a hot dayside uses a larger `f`);
* `greenhouse_factor(T_s, T_e)` is `G`.
"""

from planetgen.physics import constants


def effective_temperature_k(absorbed_flux_w_m2):
    """`T_e`, K, the temperature that radiates the absorbed flux: `(F / sigma)^(1/4)`."""
    return (absorbed_flux_w_m2 / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25


def absorbed_flux_w_m2(stellar_flux_w_m2, bond_albedo, redistribution=0.25):
    """The flux a planet absorbs averaged over its surface, W/m2: `(1 - A) S f`, `f` 1/4 for a planet that spreads the
    heat (a tidally locked world with a hot dayside uses a larger `f`)."""
    return (1.0 - bond_albedo) * stellar_flux_w_m2 * redistribution


def greenhouse_factor(surface_temperature_k, effective_k):
    """`G = (T_s / T_e)^4 - 1`, the variable the repo stores."""
    return (surface_temperature_k / effective_k) ** 4 - 1.0
