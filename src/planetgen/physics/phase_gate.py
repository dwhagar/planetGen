# planetgen/physics/phase_gate.py

"""
The phase gate of a class family: is the family's liquid a liquid at the body's own pressure and temperature (GEN.1221)
======================================================================================================================

Boss, 10:55 AM EDT 2026-10-31: a family class applies only where its solvent is liquid at the body's stored pressure and at a temperature inside the liquid
window, computed at that pressure (design: reports/class-families-and-positions-design-2026-10-31.md section 2). The windows are those of
`surface_liquids` (the saturation curves of `saturation`, the melting points of the table):

| family | liquid | window at the pressure |
|---|---|---|
| `water` | water | freezing to the boiling point at the pressure (`hydrosphere.water_boiling_k`, the same curve) |
| `co2` | carbon dioxide | the triple point (216.6 K, from 518 kPa) to the critical point (304.1 K); above the critical pressure (7.38 MPa) a dense supercritical fluid counted as the family up to `CO2_FLUID_MAX_K` |
| `methane_ethane` | the Titan sea | from the lower freezing point of the two to the boiling point of ethane at the pressure |
| `sulfuric_acid` | sulfuric acid | its freezing point (283.5 K) to its boiling point at the pressure |
| `ammonia_water` | the ammonia-water eutectic liquid | its eutectic (176 K) to its boiling point at the pressure |

`window_k` gives `(low, high)` or `None` where the liquid exists at no temperature at that pressure (carbon dioxide under 518 kPa); `holds` asks about a body's temperature; `third` says where in the window
the temperature stands (`cold`, `ideal`, `hot`: the lower, middle and upper third, the positions of the class code digits 2, 0 and 1).

Use:

    phase_gate.window_k("co2", 1.0e6)            # (216.59, 233.03)
    phase_gate.holds("water", 288.0, 101325.0)   # True
    phase_gate.third("water", 288.0, 101325.0)   # "ideal"
"""

from planetgen.physics import hydrosphere, surface_liquids

FAMILIES = ("water", "co2", "methane_ethane", "sulfuric_acid", "ammonia_water")
"""tuple: The families whose liquid has a window (the oxygen-water family and the anoxic water family both use `water`)."""

CO2_FLUID_MAX_K = 400.0
"""float: STAND-IN PG1 (L0, the top of the dense carbon dioxide family's window in `life_classes.json`, 217 to 400 K): the supercritical fluid above the critical point still counts as the family up to this temperature."""

_LIQUIDS = {"water": ("water",), "co2": ("carbon_dioxide",), "methane_ethane": ("methane", "ethane"), "sulfuric_acid": ("sulfuric_acid",),
            "ammonia_water": ("ammonia_water",)}


def window_k(family, pressure_pa):
    """`(low, high)` K of the liquid of `family` at `pressure_pa`, or `None` where it is a liquid at no temperature at that pressure."""
    if family not in _LIQUIDS:
        raise ValueError(f"{family!r} is not a family with a phase window: {FAMILIES}")
    if not pressure_pa > 0.0:
        return None
    if family == "water":
        low, high = surface_liquids.melting_k("water"), hydrosphere.water_boiling_k(pressure_pa)
        return (low, high) if high > low else None
    windows = [window for window in (surface_liquids.window_k(liquid, pressure_pa) for liquid in _LIQUIDS[family]) if window is not None]
    if not windows:
        return None
    low, high = min(window[0] for window in windows), max(window[1] for window in windows)
    if family == "co2" and surface_liquids.boiling_k("carbon_dioxide", pressure_pa) >= surface_liquids.liquid("carbon_dioxide").critical_k:
        high = max(high, CO2_FLUID_MAX_K)       # at or above the critical pressure the liquid has become a dense fluid, which the family counts
    return low, high


def holds(family, temperature_k, pressure_pa):
    """Whether the liquid of `family` is stable at `temperature_k` and `pressure_pa`."""
    window = window_k(family, pressure_pa)
    return window is not None and window[0] <= temperature_k <= window[1]


def third(family, temperature_k, pressure_pa):
    """`cold`, `ideal` or `hot`: the lower, middle or upper third of the window at the pressure that `temperature_k` stands in; `None` outside the window."""
    window = window_k(family, pressure_pa)
    if window is None or not window[0] <= temperature_k <= window[1]:
        return None
    width = (window[1] - window[0]) / 3.0
    if temperature_k < window[0] + width:
        return "cold"
    return "ideal" if temperature_k < window[0] + 2.0 * width else "hot"
