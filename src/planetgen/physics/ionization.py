# planetgen/physics/ionization.py

"""
The ionization of an air (GEN.1580)
===================================

The hottest surface a body's air can have. A greenhouse that runs away heats the surface until the air itself stops being a gas of molecules: the
molecules break into atoms (well under 10,000 K at the pressures of an air) and the atoms lose their electrons, a plasma, which radiates as no
molecular air does. No radiative model of this repository reaches there (the greenhouse table fits to 4,000 K), and no generated body is meant
to: the ceiling is the temperature the air's main gas is **half ionized** at the air's pressure, so a body that arrives at it is a mistake in the
model to explain, not a state to keep (`surface_air.Air.at_ceiling`, the warning `planets` logs, the test that fails on it).

The law is the Saha equation [R5229] for the first ionization of the atoms of the main gas, a plasma of one electron per ion,
`x^2 / (1 - x^2) = K / P` for the ionized share `x`, with `K = 2 (g_ion / g_atom) (m_e k T / (2 pi hbar^2))^(3/2) k T exp(-chi / k T)`: the
electron's translation, the Boltzmann factor of the first ionization energy `chi` ([R5228]; `data/ionization_energies.json`), and `P` the air's
total pressure (the atoms' share of it is below the whole, which moves `T` by a few percent, because `chi / k T` is large).
The statistical weights are taken as `g_ion / g_atom = 1` (a stand-in; the ground-state ratio is 1/2 to 6, which moves `T` by under 2 percent).

The main gas is the air's largest partial pressure, and its atom of the lowest first ionization energy is the one that ionizes first (carbon for
carbon dioxide, hydrogen for water vapor and hydrogen, nitrogen for nitrogen, helium for helium).
"""

import functools
import json
import math
from pathlib import Path

from planetgen.physics import constants, elements

HALF_IONIZED = 0.5
"""float: The ionized share of the main gas's atoms that defines the ceiling (a half-ionized gas is a plasma on any definition)."""

HOTTEST_BRACKET_K = 1.0e7
"""float: The hottest temperature the half-ionization search looks at."""

COOLEST_BRACKET_K = 500.0
"""float: The coolest temperature the search looks at."""

STATISTICAL_WEIGHT_RATIO = 1.0
"""float: `g_ion / g_atom` of the Saha equation. DESIGN stand-in, tier educated guess from physics, tweak L0 (the true ratio is 1/2 to 6 and moves the ceiling under 2 percent)."""

ISOTOPE_HYDROGEN = {"D": "H", "T": "H"}
"""dict: The hydrogen isotopes a gas formula writes with their own letter (`HD`, `HDO`); the electrons are hydrogen's."""

_DATA = Path(__file__).with_name("data") / "ionization_energies.json"


@functools.lru_cache(maxsize=None)
def _energies_ev():
    """The first ionization energies, eV, by element symbol."""
    with _DATA.open(encoding="utf-8") as handle:
        data = json.load(handle)
    for key in ("retrieved", "sources", "elements"):
        if not data.get(key):
            raise ValueError(f"{_DATA.name} has no {key!r}")
    return {symbol: float(item["ev"]) for symbol, item in data["elements"].items()}


def first_ionization_energy_ev(symbol):
    """The energy that takes the first electron off the atom of `symbol`, eV ([R5228]); deuterium and tritium are hydrogen.

    Raises:
        KeyError: The element is not one an air can hold (not in `data/ionization_energies.json`).
    """
    return _energies_ev()[ISOTOPE_HYDROGEN.get(symbol, symbol)]


def ionized_share(temperature_k, pressure_pa, energy_ev):
    """The share of the atoms ionized once at `temperature_k` and `pressure_pa` for an atom of first ionization energy `energy_ev` (Saha, [R5229]), in [0, 1]."""
    kt = constants.BOLTZMANN * temperature_k
    exponent = -energy_ev * constants.ELECTRON_VOLT_J / kt
    if exponent < -700.0:
        return 0.0
    translation = (constants.ELECTRON_MASS_KG * kt / (2.0 * math.pi * constants.REDUCED_PLANCK ** 2)) ** 1.5
    k_saha = 2.0 * STATISTICAL_WEIGHT_RATIO * translation * kt * math.exp(exponent)
    return math.sqrt(k_saha / (k_saha + pressure_pa))


def half_ionization_temperature_k(pressure_pa, energy_ev):
    """The temperature at which the atoms of first ionization energy `energy_ev` are half ionized at `pressure_pa`, K (bisection; the share rises with `T`).

    Raises:
        ValueError: The pressure is not positive, or the atoms are not half ionized below `HOTTEST_BRACKET_K`.
    """
    if not (pressure_pa > 0.0 and math.isfinite(pressure_pa)):
        raise ValueError(f"pressure_pa must be a positive number, not {pressure_pa!r}")
    low, high = COOLEST_BRACKET_K, HOTTEST_BRACKET_K
    if ionized_share(high, pressure_pa, energy_ev) < HALF_IONIZED:
        raise ValueError(f"the atoms are not half ionized below {HOTTEST_BRACKET_K:g} K at {pressure_pa:g} Pa")
    for _ in range(100):
        middle = math.sqrt(low * high)
        if ionized_share(middle, pressure_pa, energy_ev) < HALF_IONIZED:
            low = middle
        else:
            high = middle
        if high - low < 1.0e-6 * low:
            break
    return math.sqrt(low * high)


def main_gas_energy_ev(partials_pa):
    """The first ionization energy of the air's main gas, eV: its atom of the lowest one (the largest partial pressure names the gas)."""
    formula = max(partials_pa, key=lambda name: partials_pa[name])
    atoms = {symbol for symbol in _atoms(formula)}
    return min(first_ionization_energy_ev(symbol) for symbol in atoms)


def _atoms(formula):
    """The element symbols of a gas formula, the hydrogen isotopes (`HD`, `HDO`) as hydrogen."""
    try:
        return elements.parse_formula(formula)
    except ValueError:
        return elements.parse_formula(formula.replace("D", "H").replace("T", "H"))


def ceiling_of(partials_pa):
    """`ceiling_k` of an air, or infinity for a body with no air (the bound a stored temperature is checked against)."""
    return ceiling_k(partials_pa) if partials_pa else math.inf


def ceiling_k(partials_pa, pressure_pa=None):
    """The hottest surface the air of `partials_pa` (formula -> partial pressure, Pa) may have, K: the half-ionization temperature of its main gas at
    its total pressure (`pressure_pa`, or the sum of the partials when not given).

    Raises:
        ValueError: The air is empty.
    """
    if not partials_pa:
        raise ValueError("an air with no gas has no ionization ceiling")
    total = math.fsum(partials_pa.values()) if pressure_pa is None else pressure_pa
    return half_ionization_temperature_k(total, main_gas_energy_ev(partials_pa))
