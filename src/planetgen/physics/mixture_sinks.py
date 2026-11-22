# planetgen/physics/mixture_sinks.py

"""
What an ocean takes out of a drawn mixture (GEN.757; docs/design/composition/mixture-draw.md)
================================================================================

The mixture draw (`physics/mixture_draw.py`) says what the air is made of before the surface touches it. A surface ocean then
takes the sulfur gases (GEN.757) and the nitrogen (GEN.1834): sulfur dioxide, hydrogen sulfide and nitrogen dissolve in the ocean in proportion
to their partial pressure (Henry's law, `physics/gases.py`, mol/(kg bar)). An ocean of `m` kg per square metre over an air whose column holds `p / (g M)` moles of the
gas per square metre keeps `1 / (1 + kH m g M / P0)` of it in the air, with `P0` one bar in pascal (`remaining_fraction`,
read by `surface_air`; for sulfur once the sulfur sinks of the crust have taken theirs). The share does not depend on the pressure of the gas, so one
law serves a thin air and a thick one. An Earth ocean keeps 0.6 percent of its air's nitrogen dissolved; an ocean 10 times as deep, 6 percent.

What stays out of here: the water vapor cap over an ocean (the saturation of `physics/saturation.py`, applied by the air solve), the
outgassed inventory and the carbonate and ocean sinks of the carbon and water, and the nitrogen the crust, the sediment, the life and the subduction
hold (they are inside the nitrogen's degassing efficiency NB3, GEN.1833), (GEN.726,
`physics/surface_pressure.py`; in the rocky draw the inventory's weathering sink, `mixture_rocky.INVENTORY`), the freezing out of what
the cold surface condenses (the pressure sum of GEN.411), and the clouds.

Limits (stated): the Henry law is the physical solubility only (`Henry.dissociates` is True for both sulfur gases: the dissolved sulfur
dioxide goes on to bisulfite and sulfite, so real oceans take more; nitrogen does not dissociate, but the salt of the sea takes about 7 percent of
what fresh water holds, which is left out); it is given for 273.15 to 323.15 K and the nearest end stands in
outside for the sulfur gases (a warmer ocean takes less, so the uptake of a hot ocean is an upper bound), while nitrogen takes nothing outside it
(`RANGE_BOUND_GASES`); the ocean is taken as well mixed and at the
surface temperature.
"""

from planetgen.physics import constants, gases

SULFUR_GASES = ("SO2", "H2S")
"""tuple: The gas-table species the ocean's solubility takes sulfur out as."""

DISSOLVED_GASES = (*SULFUR_GASES, "N2")
"""tuple: The gas-table species an ocean takes out of the air by their Henry law (the sulfur gases, GEN.757, and nitrogen, GEN.1834)."""

RANGE_BOUND_GASES = ("N2",)
"""tuple: The dissolved gases whose uptake stops outside the Henry law's range (`gases.HENRY_RANGE_K`) instead of taking the nearest end (GEN.1834). Nitrogen is the
gas of every air: a sea under 273.15 K is ice-covered and cut off from the air, a sea over 323.15 K is past the data, and a law carried there over a column of a
hundred kilometers of water would take the whole air of a frozen water world (on a sample of 101 systems the 1 to 5 MPa airs of frozen water worlds fell to almost nothing)."""

BAR_PA = constants.BAR_PA
"""float: One bar, Pa (the unit of the Henry constant's pressure)."""


def remaining_fraction(formula, ts_k, gravity_m_s2, ocean_column_kg_m2):
    """
    The share of a gas of `DISSOLVED_GASES` that stays in the air over an ocean of `ocean_column_kg_m2` (kg of water per square metre of the whole
    surface) at `ts_k`: `1 / (1 + kH m g M / P0)`. 1 for no ocean.

    Raises:
        ValueError: `formula` is not one of `DISSOLVED_GASES`, or the gravity is not positive.

    A gas of `RANGE_BOUND_GASES` (nitrogen) stays whole in the air outside the Henry law's temperature range (the other gases take the nearest end).
    """
    if formula not in DISSOLVED_GASES:
        raise ValueError(f"{formula!r} is not one of the dissolved gases {DISSOLVED_GASES}")
    if not gravity_m_s2 > 0.0:
        raise ValueError(f"the gravity must be positive, not {gravity_m_s2!r}")
    if not ocean_column_kg_m2 or ocean_column_kg_m2 <= 0.0:
        return 1.0
    gas = gases.gas(formula)
    low, high = gases.HENRY_RANGE_K
    if formula in RANGE_BOUND_GASES and not low <= ts_k <= high:
        return 1.0
    kh = gas.henry.kh(min(max(ts_k, low), high))
    return 1.0 / (1.0 + kh * ocean_column_kg_m2 * gravity_m_s2 * gas.molar_mass / BAR_PA)
