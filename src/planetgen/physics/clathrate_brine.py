# planetgen/physics/clathrate_brine.py

"""
Clathrate hydrates and eutectic brines (GEN.471, docs/design/body-composition-parts/C-volatile-and-giants.md and G4-ices-brines-licences.md)
===========================================================================================================================================

Two things keep an icy crust or an ocean moon's water in a form other than ice and vapour:

* **Clathrate hydrates** hold CO, CO2, CH4, H2S and N2 in cages of water ice (and NH3 as a monohydrate) and buffer an atmosphere's
  supply of those gases. A hydrate is stable above its dissociation pressure `ln P(Pa) = A0 + A1 / T` (Marboeuf et al. 2014 Table 5,
  from Hersant 2004), with a hydrate number of 6 for all. The curve is the equilibrium of
  vapour, ice and hydrate: it holds below the freezing point of water, and for CO2 below its quadruple point Q1 (271.60 K, 1.044 MPa,
  Semenov et al. 2023), above which liquid water takes the place of ice. The second quadruple point (283.13 K, 4.494 MPa) is the
  limit of the hydrate itself and is outside what this module covers (`OutOfRange`).
* **Eutectic brines** keep water liquid below 273.15 K: a salt lowers the freezing point of the solution down to its eutectic
  temperature at the eutectic composition. The perchlorates of NaClO4, Mg(ClO4)2 and Ca(ClO4)2 reach 236, 206 and 198 K. Authors differ by 1 to 7 K, so a
  draw (`draw_eutectic_k`) scatters about the stated value by the stated spread.

Every constant carries its source key and source state in `clathrate_brine.json`; the linear liquidus between pure water and the eutectic is a
parametrisation of this module (the data hold no more than the eutectic point). Nothing here is drawn except `draw_eutectic_k`, which reads only `util.draw`.

Use in the generator (GEN.762): `physics/atmosphere.py` lets an icy body's hydrate skin hold the guest gases above their
dissociation pressure (`hydrate_hold_pa`), and `physics/hydrosphere.py` melts the base of an ice lid over a brine at the drawn
eutectic temperature (`draw_brine`).

Use:

    dissociation_pressure_pa("CH4", 100.0)       # Pa
    hydrate_stable("N2", 90.0, 1e5)              # bool
    freezing_point_k("Mg(ClO4)2", 20.0)          # K, at 20 weight percent
    liquid_range_k("NaClO4")                     # (eutectic K, 273.15 K)
"""

import json
import math
import os

from planetgen.physics import elements, gases
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "clathrate_brine.json")
"""str: The data file."""


class OutOfRange(ValueError):
    """A temperature above the curve's equilibrium with ice (the hydrate meets liquid water there)."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    sources = data.get("sources") or {}
    for key in ("clathrates", "brines", "hydrate_number", "ice_point_k"):
        if not data.get(key):
            raise DataError(f"clathrate_brine.json states no {key}")
    for section in ("clathrates", "brines"):
        for name, row in data[section].items():
            if row.get("source") not in sources:
                raise DataError(f"{section} {name!r} names the source {row.get('source')!r}, which the file does not state")
            if row.get("state") == "partly verified" and not row.get("limit"):
                raise DataError(f"{section} {name!r} is partly verified and names no limit")
    return data


_DATA = _load()
SOURCES = _DATA["sources"]
"""dict: Source key -> citation."""
HYDRATE_NUMBER = _DATA["hydrate_number"]
"""int: Water molecules per guest molecule (Marboeuf et al. 2014)."""
ICE_POINT_K = _DATA["ice_point_k"]
"""float: The freezing point of pure water, K."""
CLATHRATES = tuple(_DATA["clathrates"])
"""tuple: The guest species."""
BRINES = tuple(_DATA["brines"])
"""tuple: The salts."""


def _clathrate(species):
    try:
        return _DATA["clathrates"][species]
    except KeyError:
        raise KeyError(f"no clathrate data for {species!r}; the species are {', '.join(CLATHRATES)}") from None


def _brine(salt):
    try:
        return _DATA["brines"][salt]
    except KeyError:
        raise KeyError(f"no brine data for {salt!r}; the salts are {', '.join(BRINES)}") from None


def structure(species):
    """The hydrate structure of `species`: `I`, `II` or `monohydrate`."""
    return _clathrate(species)["structure"]


def upper_temperature_k(species):
    """The highest temperature of the ice-hydrate-vapour curve: the freezing point of water, or the quadruple point Q1 for CO2."""
    return _clathrate(species).get("upper_k", ICE_POINT_K)


def dissociation_pressure_pa(species, temperature_k):
    """
    The pressure below which the hydrate of `species` dissociates at `temperature_k`: `exp(A0 + A1 / T)`.

    Raises:
        ValueError: A temperature that is not positive.
        OutOfRange: Above `upper_temperature_k(species)`, where liquid water changes the equilibrium.
    """
    if not temperature_k > 0.0:
        raise ValueError("the temperature must be positive")
    row = _clathrate(species)
    if temperature_k > upper_temperature_k(species) + 1e-9:
        raise OutOfRange(f"{species} hydrate: {temperature_k} K is above {upper_temperature_k(species)} K, where liquid water replaces the ice")
    return math.exp(row["a0"] + row["a1"] / temperature_k)


def hydrate_stable(species, temperature_k, pressure_pa):
    """True where the hydrate of `species` is stable: the partial pressure of the guest gas is at least the dissociation pressure."""
    return pressure_pa >= dissociation_pressure_pa(species, temperature_k)


def dissociation_temperature_k(species, pressure_pa):
    """The temperature above which the hydrate of `species` dissociates at this guest partial pressure, K (the inverse of the curve)."""
    if not pressure_pa > 0.0:
        raise ValueError("the pressure must be positive")
    row = _clathrate(species)
    temperature = row["a1"] / (math.log(pressure_pa) - row["a0"])
    if temperature > upper_temperature_k(species) + 1e-9:
        raise OutOfRange(f"{species} hydrate at {pressure_pa} Pa dissociates above {upper_temperature_k(species)} K, where liquid water replaces the ice")
    return temperature


def eutectic_k(salt):
    """The eutectic temperature of the salt-water system, K (the stated value)."""
    return _brine(salt)["eutectic_k"]


def eutectic_wt_percent(salt):
    """The eutectic composition in weight percent, or `None` where the source was not reached."""
    return _brine(salt)["wt_percent"]


def draw_eutectic_k(salt):
    """The eutectic temperature of one body's brine: a normal draw about the stated value with the stated spread (authors differ by 1 to 7 K)."""
    row = _brine(salt)
    return draw.gauss(row["eutectic_k"], row["sd_k"])


def freezing_point_k(salt, wt_percent):
    """
    The freezing point of a solution of `wt_percent` of the salt, K: linear from `ICE_POINT_K` at no salt to the eutectic at the eutectic
    composition (a parametrisation: the data hold only the end points), and the eutectic beyond it.

    Raises:
        ValueError: A percentage outside 0 to 100, or a salt whose eutectic composition was not reached.
    """
    if not 0.0 <= wt_percent <= 100.0:
        raise ValueError("the composition must be 0 to 100 weight percent")
    row = _brine(salt)
    if row["wt_percent"] is None:
        raise ValueError(f"the eutectic composition of {salt} was not reached by its source; the freezing point curve is not defined")
    share = min(wt_percent / row["wt_percent"], 1.0)
    return ICE_POINT_K + (row["eutectic_k"] - ICE_POINT_K) * share


def liquid_range_k(salt):
    """The temperatures at which some solution of the salt is liquid: `(eutectic K, 273.15 K)`."""
    return (eutectic_k(salt), ICE_POINT_K)


HYDRATE_SKIN_RANGE_M = (1.0, 1000.0)
"""tuple: The depth of ice, m, whose water cages a body's atmosphere can reach (drawn log-uniform per body). A planetGen stand-in:
the cages form where gas meets ice, and how deep that goes (kinetics, regolith, cracks) has no measured figure; the bracket runs from
a frost skin to a regolith."""
ICE_DENSITY_KG_M3 = 917.0
"""float: Water ice Ih, kg/m3 (the density `water_layer.ICE_LID_DENSITY_KG_M3` uses)."""


def hydrate_hold_pa(partials_pa, temperature_k, gravity_m_s2, skin_m, ice_density_kg_m3=ICE_DENSITY_KG_M3):
    """
    The partial pressure, Pa, of each guest gas that the water cages of a `skin_m` layer of ice hold at `temperature_k`.

    `partials_pa` maps a stored species (`"co2"`) to its partial pressure. A guest above its dissociation pressure forms hydrate until
    its pressure falls to that curve, so the demand of each is the excess over the curve; the cages (one guest per `HYDRATE_NUMBER`
    waters) are shared, so when the demand of all guests together is more than the skin holds, every guest gives up the same share of
    its excess. The mole count of a partial pressure is `p / (M g)` with `M` the mean molar mass of the listed gases (a gas not listed,
    such as a helium balance, is not counted). Only the stored species that form hydrates (`CLATHRATES`) are held, and a guest above its
    curve's upper temperature (`upper_temperature_k`) is not.

    Returns:
        dict: Species -> the pressure held, Pa (only the species held, and only above zero).
    """
    total = sum(partials_pa.values())
    if not (total > 0.0 and gravity_m_s2 > 0.0 and temperature_k > 0.0 and skin_m > 0.0):
        return {}
    mean_molar_mass = sum(pressure * gases.stored(species).molar_mass for species, pressure in partials_pa.items()) / total
    excess = {}
    for species, pressure in partials_pa.items():
        formula = gases.stored(species).formula
        if pressure <= 0.0 or formula not in _DATA["clathrates"]:
            continue
        try:
            floor = dissociation_pressure_pa(formula, temperature_k)
        except OutOfRange:
            continue
        if pressure > floor:
            excess[species] = pressure - floor
    if not excess:
        return {}
    demand_mol_m2 = sum(excess.values()) / (mean_molar_mass * gravity_m_s2)
    cages_mol_m2 = ice_density_kg_m3 * skin_m / (HYDRATE_NUMBER * gases.gas("H2O").molar_mass)
    share = min(1.0, cages_mol_m2 / demand_mol_m2)
    return {species: value * share for species, value in excess.items()}


def draw_brine():
    """One body's brine: a salt (equal weights, a planetGen stand-in) and the eutectic temperature of its saturated solution, K
    (`draw_eutectic_k`). Returns `(salt, eutectic K)`."""
    salt = draw.choice(BRINES)
    return salt, draw_eutectic_k(salt)


def water_activity(salt):
    """The water activity of the salt's eutectic brine, or `None` where the source gave none."""
    return _brine(salt)["water_activity"]
