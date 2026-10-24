# planetgen/physics/mixture.py

"""
An atmosphere's mixture, read through one accessor (GEN.411, docs/design/composition/atmosphere-colour.md section 7, the
contracts of docs/design/composition/composition-spine-contracts.md)
=============================================================================================================

The mixture is a plain dict of gas-table formula (`physics/gases.py`: `"CO2"`, `"N2"`, `"He"`, `"NH3"`) to partial
pressure in Pa. Every lane that reads the air reads it here, so the storage can change underneath: today the ten
`p_<gas>_kpa` columns (`physics/atmosphere.SPECIES`) hold it, and GEN.411 and DB.71 move it to a child table that holds
any gas of the table, the balance gas included. A caller never reads the columns.

Use:

    mix = mixture.partial_pressures_pa(planet)       # {"N2": 78084.0, "O2": 20946.0, ...} Pa, zeros dropped
    mean_molar_mass = mixture.mean_molar_mass_kg_mol(mix)
    fractions = mixture.mole_fractions(mix)
"""

from planetgen.physics import atmosphere, gases

KPA_TO_PA = 1000.0

STORED_FORMULAS = {key: gases.stored(key).formula for key in atmosphere.SPECIES}
"""dict: Each stored column species (`"o2"`) and its formula in the gas table (`"O2"`)."""


def partial_pressures_pa(body):
    """
    `body`'s partial pressures.

    Returns:
        dict: Gas-table formula to Pa, only the gases with a pressure above zero. Empty for an airless body or one
        saved before GEN.85.
    """
    mixture = {}
    for key, kpa in atmosphere.partial_pressures_kpa(body).items():
        if kpa > 0.0:
            mixture[STORED_FORMULAS[key]] = kpa * KPA_TO_PA
    return mixture


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
