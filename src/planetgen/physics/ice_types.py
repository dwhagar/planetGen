# planetgen/physics/ice_types.py

"""
The kind of ice a body is made of, and how strong that ice is next to ice Ih (GEN.1177, registry 6.129)
=========================================================================================================

The roundness engine (`roundness.py`) has one ice: water ice Ih at its creep strength (about 1 MPa near 158 K). The strength research
(`reports/planetary-material-strength-2026-10-31.md` sections 2b, 9, 11c, 14) gives the other ices of an icy body as ratios to that ice: the
high-pressure phases of water (ratios by temperature from the creep laws of [R1190], GEN.1183: II 1.7 to 6.5, III about 0.5, V about 8, VI 1.3 to 5.5), methane clathrate hydrate (20 to 30 times), solid
carbon dioxide (0.1 to 0.5 MPa) and ammonia dihydrate (weaker than water ice near its melting point, stronger below 135 K). The numbers are
in `data/roundness_strengths.json` under `ice_types`; this module reads them and names a body's ice type from its components. Nothing else reads
this module yet, so no generated world changes; the engine will take `strength_ratio(...)` as a factor on its ice strength when the types are wired.

* `strength_ratio(ice_type, temperature_k)`: the factor on the Ih strength at the same shape-set temperature (a type with an absolute
  strength, solid CO2, has its own `strength_pa` and no ratio).
* `strength_pa(ice_type, ice_ih_strength_pa, temperature_k)`: the strength of that ice, Pa, for a given Ih strength.
* `water_phase(pressure_pa, temperature_k)`: the phase of water ice at a central pressure.
* `ice_type_of(components, temperature_k, central_pressure_pa)`: the type of a body from the ice-group shares of its components.

Use:

    kind = ice_types.ice_type_of(makeup.components, temperature_k=120.0)          # 'ice_ih', 'ch4_clathrate', 'co2_ice', ...
    sigma = ice_types.strength_pa(kind, roundness.ICE_STRENGTH_PA, 120.0)
"""

import json
import math
import os

from planetgen.physics import small_body_composition

DATA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "roundness_strengths.json")
"""str: The file the ice types are read from (the same file as the roundness engine's strengths)."""

WATER_PHASES = ("ice_ih", "ice_ii", "ice_iii", "ice_v", "ice_vi")
"""tuple: The phases of water ice that `water_phase` names (ice VII and ice X take the ice VI value)."""

_WATER_ROWS = frozenset({"water_ice", "ice_ii", "ice_iii", "ice_v", "ice_vi", "ice_vii", "ice_x"})
_PHASE_OF_ROW = {"water_ice": "ice_ih", "ice_ii": "ice_ii", "ice_iii": "ice_iii", "ice_v": "ice_v", "ice_vi": "ice_vi", "ice_vii": "ice_vi", "ice_x": "ice_vi"}


class IceTypeError(ValueError):
    """An ice-type table that is refused."""


def _load(path=DATA_FILE):
    with open(path, encoding="utf-8") as handle:
        document = json.load(handle)
    section = document.get("ice_types")
    if not isinstance(section, dict):
        raise IceTypeError("the file has no 'ice_types' section")
    types = section.get("types") or {}
    for name in WATER_PHASES + ("ch4_clathrate", "co2_ice", "nh3_dihydrate"):
        if name not in types:
            raise IceTypeError(f"ice type {name!r} is missing")
    for name, entry in types.items():
        given = [key for key in ("strength_ratio", "strength_pa", "ratio_by_temperature_k") if key in entry]
        if len(given) != 1:
            raise IceTypeError(f"{name} needs exactly one of strength_ratio, strength_pa, ratio_by_temperature_k")
        if "strength_ratio" in entry and not entry["strength_ratio"] > 0.0:
            raise IceTypeError(f"{name} strength_ratio must be positive")
        if "strength_pa" in entry and not entry["strength_pa"] > 0.0:
            raise IceTypeError(f"{name} strength_pa must be positive")
        if "ratio_by_temperature_k" in entry:
            points = entry["ratio_by_temperature_k"]
            if len(points) < 2 or any(b[0] <= a[0] for a, b in zip(points, points[1:])) or any(r <= 0.0 for _, r in points):
                raise IceTypeError(f"{name} ratio_by_temperature_k needs rising temperatures and positive ratios")
    return types, section.get("phases") or {}, section.get("composition") or {}, float(document["materials"]["ice"]["strength_pa"])


TYPES, PHASES, COMPOSITION, ICE_IH_REFERENCE_PA = _load()
"""The loaded `ice_types.types`, `.phases` and `.composition` of `roundness_strengths.json`, and the strength of `materials.ice` there, Pa (the default ice of `roundness.py`, the unit of `strength_ratio` for an absolute type)."""

ICE_TYPES = tuple(TYPES)
"""tuple: The names of the ice types."""


def _ratio_at(points, temperature_k):
    """The ratio of a `[temperature, ratio]` list, log-linear in the ratio, held at the end values."""
    if temperature_k is None:
        temperature_k = points[len(points) // 2][0]
    if temperature_k <= points[0][0]:
        return points[0][1]
    if temperature_k >= points[-1][0]:
        return points[-1][1]
    for (t0, r0), (t1, r1) in zip(points, points[1:]):
        if temperature_k <= t1:
            share = (temperature_k - t0) / (t1 - t0)
            return math.exp(math.log(r0) + share * (math.log(r1) - math.log(r0)))
    return points[-1][1]


def strength_ratio(ice_type, temperature_k=None):
    """
    The strength of `ice_type` over the strength of ice Ih at the same shape-set temperature. A type with a temperature table (ammonia dihydrate)
    is read at `temperature_k` (the middle point of its table when none is given); solid CO2, which has an absolute strength, is expressed
    against `ICE_IH_REFERENCE_PA`, the default ice strength of the engine.
    """
    entry = _entry(ice_type)
    if "strength_ratio" in entry:
        return float(entry["strength_ratio"])
    if "ratio_by_temperature_k" in entry:
        return _ratio_at(entry["ratio_by_temperature_k"], temperature_k)
    return entry["strength_pa"] / ICE_IH_REFERENCE_PA


def _entry(ice_type):
    if ice_type not in TYPES:
        raise IceTypeError(f"unknown ice type {ice_type!r}; the types are {sorted(TYPES)}")
    return TYPES[ice_type]


def strength_pa(ice_type, ice_ih_strength_pa, temperature_k=None):
    """The strength of `ice_type`, Pa, for the strength of ice Ih `ice_ih_strength_pa` at the same temperature; solid CO2 keeps its own strength."""
    entry = _entry(ice_type)
    if "strength_pa" in entry:
        return float(entry["strength_pa"])
    return ice_ih_strength_pa * strength_ratio(ice_type, temperature_k)


def water_phase(pressure_pa, temperature_k=None):
    """
    The phase of water ice at a central pressure (and temperature for the narrow wedges of III and V): ice Ih under 0.21 GPa, then II (the cold
    phase) to 0.63 GPa, except III at 240 to 256 K from 0.21 to 0.35 GPa and V at 250 K and over from 0.35 GPa, then VI to 2.2 GPa and VI
    again beyond it (the value of VII and X, strength research section 14). Boundaries recalled from the water phase diagram (row IT8).
    """
    phases = PHASES
    if pressure_pa is None or pressure_pa < phases["ih_below_pa"]:
        return "ice_ih"
    if pressure_pa >= phases["ii_below_pa"]:
        return "ice_vi"
    if temperature_k is not None:
        low_k, high_k = phases["iii_temperature_k"]
        p_low, p_high = phases["iii_pressure_pa"]
        if low_k <= temperature_k <= high_k and p_low <= pressure_pa < p_high:
            return "ice_iii"
        v_low, v_high = phases["v_pressure_pa"]
        if temperature_k >= phases["v_above_k"] and v_low <= pressure_pa < v_high:
            return "ice_v"
    return "ice_ii"


def _ice_group_masses(components):
    """`(water, others)`: the mass of every water ice row by phase name, and of every other ice row (ammonia, methane, carbon dioxide...) by row id."""
    water, others = {}, {}
    for row, mass in components.items():
        if mass <= 0.0:
            continue
        group = small_body_composition.group_of(row)
        if row in _WATER_ROWS:
            phase = _PHASE_OF_ROW[row]
            water[phase] = water.get(phase, 0.0) + mass
        elif group == "ice":
            others[row] = others.get(row, 0.0) + mass
    return water, others


def ice_type_of(components, temperature_k=None, central_pressure_pa=None):
    """
    The ice type of a body from the mass fractions of its components (`SmallBodyMakeup.components`, a row id of the shared materials table to a
    mass fraction; a makeup is also accepted). `None` for a body with no ice. The ice group is the water rows (the phase named by the row, or by
    `water_phase` for plain water ice at `central_pressure_pa`) and every other ice row. A body is named:

    * `ch4_clathrate` when the methane that water can cage (`ch4_mass_in_hydrate` of the hydrate by mass) is `clathrate_min_share_of_ice` of the
      ice or more (hydrate matters from 20 to 30 percent by volume, strength research section 11c);
    * `nh3_dihydrate` when ammonia is `ammonia_min_share_of_ice` of the ice or more and water is there to make the dihydrate;
    * otherwise the ice that has the most mass: a water phase, or the row id of another ice (`co2_ice`); a row with no strength of its own
      (nitrogen, carbon monoxide, ethane...) is named by its row id and `strength_ratio` refuses it, so the caller decides.
    """
    components = getattr(components, "components", components)
    water, others = _ice_group_masses(components)
    water_total, others_total = sum(water.values()), sum(others.values())
    total = water_total + others_total
    if total <= 0.0:
        return None
    caged = min(others.get("ch4_ice", 0.0) / COMPOSITION["ch4_mass_in_hydrate"], water_total) if water_total > 0.0 else 0.0
    if caged / total >= COMPOSITION["clathrate_min_share_of_ice"]:
        return "ch4_clathrate"
    if water_total > 0.0 and others.get("nh3_ice", 0.0) / total >= COMPOSITION["ammonia_min_share_of_ice"]:
        return "nh3_dihydrate"
    if others and max(others.values()) > water_total:
        return max(others, key=others.get)
    if not water:
        return None
    if set(water) == {"ice_ih"}:
        return water_phase(central_pressure_pa, temperature_k)
    return max(water, key=water.get)
