# planetgen/physics/mixture_sinks.py

"""
What an ocean takes out of a drawn mixture (GEN.757; docs/design/composition/mixture-draw.md)
================================================================================

The mixture draw (`physics/mixture_draw.py`) says what the air is made of before the surface touches it. A surface ocean then
takes the sulfur gases: sulfur dioxide and hydrogen sulfide dissolve in the ocean in proportion to their partial pressure (Henry's
law, `physics/gases.py`, mol/(kg bar)). An ocean of `m` kg per square metre over an air whose column holds `p / (g M)` moles of the
gas per square metre keeps `1 / (1 + kH m g M / P0)` of it in the air, with `P0` one bar in pascal (`sulfur_remaining_fraction`,
read by `surface_air` once the sulfur sinks of the crust have taken theirs).

What stays out of here: the water vapour cap over an ocean (the saturation of `physics/saturation.py`, applied by the air solve), the
outgassed inventory and the carbonate and ocean sinks of the carbon, water and nitrogen the mantle gave (GEN.726,
`physics/surface_pressure.py`; in the rocky draw the inventory's weathering sink, `mixture_rocky.INVENTORY`), the freezing out of what
the cold surface condenses (the pressure sum of GEN.411), and the clouds.

Limits (stated): the Henry law is the physical solubility only (`Henry.dissociates` is True for both sulfur gases: the dissolved sulfur
dioxide goes on to bisulfite and sulfite, so real oceans take more); it is given for 273.15 to 323.15 K and the nearest end stands in
outside (a warmer ocean takes less, so the uptake of a hot ocean is an upper bound); the ocean is taken as well mixed and at the
surface temperature.
"""

from planetgen.physics import gases

SULFUR_GASES = ("SO2", "H2S")
"""tuple: The gas-table species the ocean's solubility takes sulfur out as."""

BAR_PA = 1.0e5
"""float: One bar, Pa (the unit of the Henry constant's pressure)."""


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


