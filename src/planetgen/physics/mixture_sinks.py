# planetgen/physics/mixture_sinks.py

"""
What an ocean takes out of a drawn mixture (GEN.757; docs/design/composition/mixture-draw.md)
================================================================================

The mixture draw (`physics/mixture_draw.py`) says what the air is made of before the surface touches it. A surface ocean then
limits two things the mixture alone does not know:

* **The water vapour.** Over liquid water or ice the vapour cannot stand above its saturation pressure at the surface temperature
  (`physics/saturation.py`): the rest condenses into the ocean. The draw's water share is capped at that pressure; above the critical
  point there is no cap.
* **The sulfur gases.** Sulfur dioxide and hydrogen sulfide dissolve in the ocean in proportion to their partial pressure (Henry's
  law, `physics/gases.py`, mol/(kg bar)). An ocean of `m` kg per square metre over an air whose column holds `p / (g M)` moles of the
  gas per square metre keeps `1 / (1 + kH m g M / P0)` of it in the air, with `P0` one bar in pascal.

What stays out of here: the outgassed inventory and the carbonate and ocean sinks of the carbon, water and nitrogen the mantle gave
(GEN.726, `physics/surface_pressure.py`; in the rocky draw the inventory's weathering sink, `mixture_rocky.INVENTORY`), the freezing
out of what the cold surface condenses (the pressure sum of GEN.411), and the clouds. These functions only cap what the drawn mixture can
show.

Use:

    sunk = apply_ocean_sinks(mixture.mole_fractions, pressure_pa, ts_k, gravity_m_s2, ocean_present=True, ocean_column_kg_m2=2.7e6)
    sunk.mole_fractions, sunk.pressure_pa           # the air after, renormalised; the pressure is the sum of what is left
    sunk.condensed_water_pa, sunk.dissolved_pa      # what left the air, by partial pressure before the sink

Limits (stated): the Henry law is the physical solubility only (`Henry.dissociates` is True for both sulfur gases: the dissolved sulfur
dioxide goes on to bisulfite and sulfite, so real oceans take more); it is given for 273.15 to 323.15 K and the nearest end stands in
outside (a warmer ocean takes less, so the uptake of a hot ocean is an upper bound); the ocean is taken as well mixed and at the surface
temperature; the saturation cap is the full saturation (a relative humidity of one: a bound, not Earth's 0.6 to 0.8).
"""

import math
from dataclasses import dataclass

from planetgen.physics import gases, mixture, saturation

SULFUR_GASES = ("SO2", "H2S")
"""tuple: The gas-table species the ocean's solubility takes sulfur out as."""

BAR_PA = 1.0e5
"""float: One bar, Pa (the unit of the Henry constant's pressure)."""


@dataclass(frozen=True)
class OceanSinks:
    """
    A mixture after an ocean.

    Attributes:
        mole_fractions (dict): Gas-table formula -> mole fraction of what is left, summing to 1 (empty when nothing is left).
        pressure_pa (float): The sum of the partial pressures left.
        condensed_water_pa (float): The water vapour pressure that left the air to the ocean, Pa.
        dissolved_pa (dict): Each sulfur gas that dissolved: its partial pressure taken, Pa (before the renormalisation).
    """

    mole_fractions: dict
    pressure_pa: float
    condensed_water_pa: float
    dissolved_pa: dict


def water_vapour_cap_pa(ts_k, ocean_present=True):
    """
    The most water vapour the air over an ocean (or ice) at `ts_k` can hold, Pa: the saturation pressure of water; `inf` with no ocean
    to hold the balance and above the critical point (about 647 K), where there is no liquid to limit it.
    """
    if not ocean_present:
        return math.inf
    return saturation.pressure_pa("H2O", ts_k)


def sulfur_remaining_fraction(formula, ts_k, gravity_m_s2, ocean_column_kg_m2):
    """
    The share of a sulfur gas that stays in the air over an ocean of `ocean_column_kg_m2` (kg of water per square metre of the whole
    surface) at `ts_k`: `1 / (1 + kH m g M / P0)`. 1 for no ocean.

    Raises:
        ValueError: `formula` is not a sulfur gas of `SULFUR_GASES`, or the gravity is not positive.
    """
    if formula not in SULFUR_GASES:
        raise ValueError(f"{formula!r} is not one of the sulfur gases {SULFUR_GASES}")
    if not gravity_m_s2 > 0.0:
        raise ValueError(f"the gravity must be positive, not {gravity_m_s2!r}")
    if not ocean_column_kg_m2 or ocean_column_kg_m2 <= 0.0:
        return 1.0
    gas = gases.gas(formula)
    low, high = gases.HENRY_RANGE_K
    kh = gas.henry.kh(min(max(ts_k, low), high))
    return 1.0 / (1.0 + kh * ocean_column_kg_m2 * gravity_m_s2 * gas.molar_mass / BAR_PA)


def apply_ocean_sinks(mole_fractions, pressure_pa, ts_k, gravity_m_s2, *, ocean_present, ocean_column_kg_m2=0.0):
    """
    `mole_fractions` (a mixture at `pressure_pa`) after the ocean's water cap and sulfur uptake (see the module text).

    Raises:
        ValueError: The pressure, the temperature or the gravity is not positive.
    """
    if not (pressure_pa > 0.0 and ts_k > 0.0 and gravity_m_s2 > 0.0):
        raise ValueError("the pressure, the temperature and the gravity must be positive")
    partial = {formula: fraction * pressure_pa for formula, fraction in mole_fractions.items() if fraction > 0.0}
    condensed, dissolved = 0.0, {}
    if ocean_present:
        water = partial.get("H2O", 0.0)
        cap = water_vapour_cap_pa(ts_k, True)
        if water > cap:
            condensed = water - cap
            partial["H2O"] = cap
        for formula in SULFUR_GASES:
            if formula in partial:
                kept = sulfur_remaining_fraction(formula, ts_k, gravity_m_s2, ocean_column_kg_m2)
                dissolved[formula] = partial[formula] * (1.0 - kept)
                partial[formula] *= kept
    partial = {formula: value for formula, value in partial.items() if value > 0.0}
    if not partial:
        return OceanSinks({}, 0.0, condensed, dissolved)
    total = mixture.total_pressure_pa(partial)
    return OceanSinks({formula: value / total for formula, value in partial.items()}, total, condensed, dissolved)
