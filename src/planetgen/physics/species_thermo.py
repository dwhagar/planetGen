# planetgen/physics/species_thermo.py

"""
Thermodynamics of the gases and condensed phases of the atmosphere chemistry (GEN.451, GEN.454, GEN.499)
=======================================================================================================

The one source of free energies for every equilibrium the atmosphere chain computes: the outgassing speciation over a
mantle of a given oxygen fugacity (`physics/outgassing.py`), the giant-planet and brown-dwarf chemistry (GEN.454), the
lava-world vapour (GEN.499) and the buffer reactions of `physics/oxygen_fugacity.py`. It holds the NASA Glenn 9-coefficient
fits ([R2968], NASA/TP-2002-211556) of 87 species: ideal gases and a few condensed phases, each with
its elements, its molar mass and the temperature intervals of its fit.

The data are `physics/data/nasa_glenn_species.json`, read from the NASA CEA database `thermo.inp` by
`scripts/build_species_thermo.py` and checked in (the raw database is not in the repository). The loader refuses a species
without elements, a molar mass more than 0.01 g/mol from the sum over its elements (`physics.elements`), a fit that does
not start where the last ended, and a name listed twice.

    H/RT = -a1/T^2 + a2 ln(T)/T + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
    S/R  = -a1/(2 T^2) - a2/T + a3 ln T + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2
    G/RT = H/RT - S/R                                  (standard state: ideal gas at 1 bar)

Use:

    g_over_rt("CO2", 1200.0)                                   # standard Gibbs energy over RT
    log10_k({"CO2": 1, "CO": -1, "O2": -0.5}, 1200.0)          # log10 of the equilibrium constant, bar to the net power
    species("COS").elements                                    # {"C": 1.0, "O": 1.0, "S": 1.0}
    species("COS").formula                                     # "OCS", the formula of the gas table (`physics/gases.py`)

A temperature outside every fit of a species raises `ValueError` (the fits start at 200 or 300 K and stop at 6000 or 20000 K,
condensed phases hold over a few hundred kelvin only): a caller that needs a colder body does not use a Gibbs
equilibrium (the frost atmospheres use the saturation curves of `physics/saturation.py`).
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen.physics import elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "nasa_glenn_species.json")
"""str: The data file."""

FORMULA_OF_NAME = {"C2H2,acetylene": "C2H2", "COS": "OCS", "AL": "Al", "ALO": "AlO", "AL2O": "Al2O", "ALO2": "AlO2", "AL2O2": "Al2O2"}
"""dict: The names in the NASA database that the gas table (`physics/gases.py`) writes differently."""

MOLAR_MASS_TOLERANCE_G_MOL = 0.01
"""float: How far the database's molar mass may be from the sum of the element table's atomic weights."""

GAS_CONSTANT = 8.314462618
"""float: J/(mol K), exact in the SI since 2019."""


@dataclass(frozen=True)
class Fit:
    """One temperature interval of a species' fit: the seven coefficients `a` and the two integration constants `b`."""

    t_min_k: float
    t_max_k: float
    a: tuple
    b: tuple


@dataclass(frozen=True)
class Species:
    """
    One gas or condensed phase.

    Attributes:
        name (str): Its name in the NASA database (`"NaOH"`, `"H2O(L)"`).
        formula (str): Its formula as the gas table writes it (`"OCS"` where the database says `"COS"`); the name otherwise.
        elements (dict): Atoms per formula unit, symbol -> count.
        condensed (bool): True for a liquid or a solid.
        molar_mass (float): kg/mol, from the element table.
        formation_enthalpy_j_mol (float): At 298.15 K, from the database.
        fits (tuple): The `Fit` intervals, in rising temperature.
    """

    name: str
    formula: str
    elements: dict
    condensed: bool
    molar_mass: float
    formation_enthalpy_j_mol: float
    fits: tuple

    @property
    def range_k(self):
        """(lowest, highest) temperature the fits cover, K."""
        return self.fits[0].t_min_k, self.fits[-1].t_max_k

    def fit(self, temperature_k):
        """The `Fit` that covers `temperature_k`; `ValueError` outside every interval."""
        for fit in self.fits:
            if fit.t_min_k <= temperature_k <= fit.t_max_k:
                return fit
        low, high = self.range_k
        raise ValueError(f"{self.name}: no thermodynamic fit at {temperature_k} K (the fits cover {low:g} to {high:g} K)")

    def h_over_rt(self, temperature_k):
        """The molar enthalpy over R T."""
        t = temperature_k
        a1, a2, a3, a4, a5, a6, a7 = self.fit(t).a
        b1 = self.fit(t).b[0]
        return -a1 / t ** 2 + a2 * math.log(t) / t + a3 + a4 * t / 2 + a5 * t ** 2 / 3 + a6 * t ** 3 / 4 + a7 * t ** 4 / 5 + b1 / t

    def s_over_r(self, temperature_k):
        """The molar entropy over R, at the 1 bar standard state."""
        t = temperature_k
        a1, a2, a3, a4, a5, a6, a7 = self.fit(t).a
        b2 = self.fit(t).b[1]
        return -a1 / (2 * t ** 2) - a2 / t + a3 * math.log(t) + a4 * t + a5 * t ** 2 / 2 + a6 * t ** 3 / 3 + a7 * t ** 4 / 4 + b2

    def cp_over_r(self, temperature_k):
        """The heat capacity over R."""
        t = temperature_k
        a1, a2, a3, a4, a5, a6, a7 = self.fit(t).a
        return a1 / t ** 2 + a2 / t + a3 + a4 * t + a5 * t ** 2 + a6 * t ** 3 + a7 * t ** 4

    def g_over_rt(self, temperature_k):
        """The standard molar Gibbs energy over R T (the chemical potential of the pure phase at 1 bar)."""
        return self.h_over_rt(temperature_k) - self.s_over_r(temperature_k)

    def supercooled_g_over_rt(self, temperature_k, below_k):
        """
        The Gibbs energy over R T of a liquid below its melting point, from its own lowest fit extended down to `below_k` K
        (a liquid's heat capacity changes little over several hundred kelvin: the database's liquid oxides are constant-Cp fits).

        Raises:
            ValueError: `temperature_k` is below `below_k`, or above the fits.
        """
        first = self.fits[0]
        if temperature_k >= first.t_min_k:
            return self.g_over_rt(temperature_k)
        if temperature_k < below_k:
            raise ValueError(f"{self.name}: the supercooled liquid is only taken down to {below_k:g} K, not {temperature_k:g} K")
        t = temperature_k
        a1, a2, a3, a4, a5, a6, a7 = first.a
        b1, b2 = first.b
        h = -a1 / t ** 2 + a2 * math.log(t) / t + a3 + a4 * t / 2 + a5 * t ** 2 / 3 + a6 * t ** 3 / 4 + a7 * t ** 4 / 5 + b1 / t
        s = -a1 / (2 * t ** 2) - a2 / t + a3 * math.log(t) + a4 * t + a5 * t ** 2 / 2 + a6 * t ** 3 / 3 + a7 * t ** 4 / 4 + b2
        return h - s


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for key in ("title", "source", "form", "species"):
        if not document.get(key):
            raise DataError(f"{DATA_FILE}: missing '{key}'")
    table = {}
    for name, row in document["species"].items():
        if name in table:
            raise DataError(f"{name} is listed twice")
        atoms = row.get("elements")
        if not atoms:
            raise DataError(f"{name}: no elements")
        mass_g_mol = sum(count * elements.element(symbol).atomic_weight for symbol, count in atoms.items())
        if abs(mass_g_mol - row["molar_mass_g_mol"]) > MOLAR_MASS_TOLERANCE_G_MOL:
            raise DataError(f"{name}: molar mass {row['molar_mass_g_mol']} g/mol is not the sum of its elements, {mass_g_mol:.4f}")
        fits = tuple(Fit(f["t_min_k"], f["t_max_k"], tuple(f["a"]), tuple(f["b"])) for f in row["intervals"])
        if not fits:
            raise DataError(f"{name}: no fit")
        for earlier, later in zip(fits, fits[1:]):
            if earlier.t_max_k != later.t_min_k:
                raise DataError(f"{name}: the fits do not join at {earlier.t_max_k} K")
        table[name] = Species(name, FORMULA_OF_NAME.get(name, name), dict(atoms), bool(row["condensed"]), mass_g_mol * elements.KG_PER_G,
                              row["formation_enthalpy_j_mol"], fits)
    return table


SPECIES = _load()
"""dict: Every species by its name in the NASA database."""

GASES = tuple(name for name, row in SPECIES.items() if not row.condensed)
"""tuple: The names of the gases."""

CONDENSED = tuple(name for name, row in SPECIES.items() if row.condensed)
"""tuple: The names of the condensed phases."""


def species(name):
    """The `Species` called `name`; `KeyError` for a name the file does not hold."""
    return SPECIES[name]


BY_FORMULA = {row.formula: row for row in SPECIES.values() if not row.condensed}
"""dict: Every gas by its formula as the gas table writes it (`"OCS"` for the database's `"COS"`)."""


def gas(formula):
    """The gas whose formula (as the gas table writes it) is `formula`; `KeyError` for a gas the file does not hold."""
    return BY_FORMULA[formula]


def g_over_rt(name, temperature_k):
    """The standard Gibbs energy over R T of species `name` at `temperature_k`."""
    return SPECIES[name].g_over_rt(temperature_k)


def log10_k(stoichiometry, temperature_k):
    """
    log10 of the equilibrium constant, with the 1 bar standard state, of a reaction.

    Args:
        stoichiometry (dict): Species name -> coefficient, products positive and reactants negative
            (`{"CO2": 1, "CO": -1, "O2": -0.5}` for CO + 1/2 O2 = CO2).
        temperature_k (float): Kelvin.
    """
    return -math.fsum(n * g_over_rt(name, temperature_k) for name, n in stoichiometry.items()) / math.log(10.0)
