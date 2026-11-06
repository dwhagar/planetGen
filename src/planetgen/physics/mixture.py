# planetgen/physics/mixture.py

"""
An atmosphere's mixture, read through one accessor (GEN.411, docs/design/composition/atmosphere-colour.md section 7, the
contracts of docs/design/composition/composition-spine-contracts.md)
=============================================================================================================

The mixture is a plain dict of gas-table formula (`physics/gases.py`: `"CO2"`, `"N2"`, `"He"`, `"NH3"`) to partial
pressure in Pa. Every lane that reads the air reads it here, so the storage can change underneath: since DB.71 the child
tables `planet_gases` and `moon_gases` hold it, any gas of the table, the balance gas included.

Use:

    mix = mixture.partial_pressures_pa(planet)       # {"N2": 78084.0, "O2": 20946.0, ...} Pa, zeros dropped
    mean_molar_mass = mixture.mean_molar_mass_kg_mol(mix)
    fractions = mixture.mole_fractions(mix)
"""

import re

import numpy as np

from planetgen.physics import atmosphere, constants, gases, species_thermo


def partial_pressures_pa(body):
    """
    `body`'s partial pressures.

    Returns:
        dict: Gas-table formula to Pa, only the gases with a pressure above zero. Empty for an airless body or one
        saved before GEN.85.
    """
    return {formula: pa for formula, pa in atmosphere.gas_partials_pa(body).items() if pa > 0.0}


def total_pressure_pa(mixture):
    """The sum of the partial pressures, Pa."""
    return sum(mixture.values())


def mole_fractions(mixture):
    """
    Each gas's share of the molecules.

    Raises:
        ValueError: The mixture is empty.
    """
    total = total_pressure_pa(mixture)
    if total <= 0.0:
        raise ValueError("an empty mixture has no mole fractions")
    return {formula: pressure / total for formula, pressure in mixture.items()}


def mean_molar_mass_kg_mol(mixture):
    """
    The mean molar mass, the sum of x_i M_i with M_i from the gas table, kg/mol.

    Raises:
        ValueError: The mixture is empty.
        planetgen.physics.gases.UnknownGas: A key is not a row of the gas table.
    """
    return sum(fraction * gases.gas(formula).molar_mass for formula, fraction in mole_fractions(mixture).items())


def _atoms(formula):
    """The number of atoms in a molecule written without brackets (`"HC3N"` is 5)."""
    return sum(int(count or 1) for _symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula))


def cp_over_r_without_fit(formula, temperature_k):
    """
    The heat capacity at constant pressure over R of a gas that has no NASA Glenn fit (a trace gas such as CH3Cl; stand-in row IR9, GEN.852), at `temperature_k` (a number or an array):

    1. the heat capacities the gas table gives at 100, 300 and 1000 K (`gases.Gas.cp_j_mol_k`), interpolated linearly in temperature and held at the end values beyond them (a polyatomic gas
       gains heat capacity above 1000 K, so this is low there);
    2. when the table has none, a rigid-rotor count of the atoms: 5/2 for one, 7/2 for two, and 4 + (3N - 6) / 2 for a molecule of N atoms, half its vibrations awake.
    """
    known = sorted((float(t), v / constants.R) for t, v in gases.gas(formula).cp_j_mol_k.items() if v is not None)
    if known:
        return np.interp(temperature_k, [t for t, _ in known], [v for _, v in known])
    atoms = _atoms(formula)
    return np.full_like(np.asarray(temperature_k, dtype=float), {1: 2.5, 2: 3.5}.get(atoms, 4.0 + 0.5 * (3 * atoms - 6)))


def cp_over_r(formula, temperature_k):
    """The heat capacity over R of a gas at `temperature_k`: its NASA Glenn fit, held at the fit's ends outside it, or `cp_over_r_without_fit` for a gas the fits do not hold."""
    try:
        thermo = species_thermo.gas(formula)
    except KeyError:
        return cp_over_r_without_fit(formula, temperature_k)
    low, high = thermo.range_k
    return thermo.cp_over_r(min(max(temperature_k, low), high))


def heat_capacity_ratio(mixture, temperature_k):
    """
    `gamma = cp / (cp - R)` of the mixture at `temperature_k` (mole-fraction mean of the NASA fits' `cp / R`; the temperature is held
    inside the fits' range, so the cold end of a fit stands in below it; a gas with no fit takes `cp_over_r_without_fit`).
    """
    cp_over_r_sum = 0.0
    for formula, fraction in mole_fractions(mixture).items():
        if fraction == 0.0:
            continue
        cp_over_r_sum += fraction * float(cp_over_r(formula, temperature_k))
    return cp_over_r_sum / (cp_over_r_sum - 1.0)
