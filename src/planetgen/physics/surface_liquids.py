# planetgen/physics/surface_liquids.py

"""
The liquids a surface can hold (GEN.1035; the draft is reports/design-drafts/ocean-tides-and-liquids-2026-10-30.md)
==================================================================================================================

A pure module: a temperature, a pressure and what a world holds in, and which liquids are stable there, which of them a world has, and a liquid's density
and viscosity out. It reads nothing from a body and draws nothing; the stored `surface_liquid`, its coverage and the description text are the wiring that
follows (GEN.1035, lock work), and the ocean tide (GEN.1034) reads the density from here.

`physics/data/surface_liquids.json` is the table. A row is a name, a kind (cryogenic, hydrocarbon, molecular, aqueous, brine, acid, sulfur, salt, metal,
sulfide, silicate, carbonate), the melting point, the boiling point at one bar and its heat of vaporization, the critical temperature, a density with its
thermal expansion, a viscosity curve, a surface tension, and the tier (T1 to T9) and tweak level (L0 to L4) of its numbers with their source. A row for
mercury, a molten salt or a new lava is a line in the file. Mercury, the alkali metals and the other very exotic rows are here so that the table is
complete; nothing in the generator puts them on a surface until a rule for their occurrence is built.

1. **Window.** A liquid is stable between its melting point and the temperature at which it boils at the surface pressure. The boiling temperature of a row
   that names a `saturation` species is the root of the saturation curve of `physics/saturation.py` (the Wagner fits to NIST, 0.3 percent), and the liquid
   is absent below its triple pressure (carbon dioxide is a liquid only above 5.2 bar); the other rows follow Clausius-Clapeyron from their boiling point at
   one bar, `ln(P / 1 bar) = -(dH / R)(1 / T - 1 / T_b)`. A row with a critical temperature is capped at it. The silicate and carbonate rows are liquid
   above their liquidus; between the solidus and the liquidus a lava is a slush, which this module does not model.
2. **Supply.** A liquid is on a world only if something supplies it (`candidates`): a volatile whose partial pressure in the air reaches its saturation pressure
   at the surface temperature condenses (`condensing`), ranked by how far past saturation it is; otherwise the amount of its ice or ore the body holds (the
   `reservoir`) ranks it. The rule is physics for the first and a stand-in for the second (OT3 in the registry).
3. **Seas of several liquids.** `ideal_sea` gives the mole fractions of a sea in equilibrium with an air by Raoult's law, `x_i = p_i / P_sat,i(T)`, with an
   involatile solvent (ethane, on Titan) taking the rest: Titan's lakes are methane, ethane and nitrogen (Cordier et al. 2009, recalled). Ideal mixing
   overstates the nitrogen the cold lakes hold, so a permanent gas keeps `SOLUBILITY_FACTOR` of its ideal share; the result is a stand-in (OT3).
4. **Density and viscosity.** The density is linear from the reference point in the expansion coefficient, held to 30 percent of the reference at least (a
   liquid near its critical point is lighter than that line says); the viscosity is interpolated linearly in its logarithm between the row's points and
   held at the ends. A mixture's density is the mole-weighted molar volume (`sea_density_kg_m3`).

Every number that is not a verified measurement is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.69, rows LQ1 to LQ10).
"""

import json
import math
import os
from dataclasses import dataclass
from typing import Optional

from planetgen.physics import elements, gases, saturation

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "surface_liquids.json")

GAS_CONSTANT_J_MOL_K = 8.314
"""float: R, J/(mol K)."""
ONE_BAR_PA = 1.0e5
SOLUBILITY_FACTOR = {"N2": 0.3}
"""dict: STAND-IN OT3: the share of Raoult's ideal mole fraction a permanent gas keeps in a cold hydrocarbon sea (ideal mixing overstates the nitrogen Titan's lakes
hold, the 0.1 to 0.2 of Cordier et al. 2009, recalled)."""
MIN_DENSITY_SHARE = 0.3
"""float: STAND-IN: the density line is held at or above this share of the reference density (a liquid near its critical point)."""
KINDS = ("cryogenic", "hydrocarbon", "molecular", "aqueous", "brine", "acid", "sulfur", "salt", "metal", "sulfide", "silicate", "carbonate")
TIERS = tuple(f"T{n}" for n in range(1, 10))
TWEAKS = tuple(f"L{n}" for n in range(5))


class SurfaceLiquidError(ValueError):
    """A table that is refused, or a liquid that is not in it."""


@dataclass(frozen=True)
class Liquid:
    """
    A row of the table.

    Attributes:
        id (str): The row's key.
        name (str): What it is called.
        kind (str): One of `KINDS`.
        saturation (str or None): The formula of the saturation curve in `physics/saturation.py` that gives its window, else `None`.
        melting_k, boiling_k (float): The melting point and the boiling point at one bar, K.
        vaporization_kj_mol (float): The heat of vaporization, kJ/mol, for the Clausius-Clapeyron curve of a row with no saturation curve.
        critical_k (float or None): The critical temperature, K.
        density_at_k, density_kg_m3, expansion_per_k (float): The density at a temperature and its fall per kelvin (relative).
        viscosity_mpa_s (tuple): `(K, mPa s)` points, ascending in temperature or not; interpolated in the logarithm.
        surface_tension_mn_m (float): At the reference temperature.
        tier, tweak (str): The source tier and tweak level of the row's numbers.
        source, note (str): Where they came from.
    """

    id: str
    name: str
    kind: str
    saturation: Optional[str]
    melting_k: float
    boiling_k: float
    vaporization_kj_mol: float
    critical_k: Optional[float]
    density_at_k: float
    density_kg_m3: float
    expansion_per_k: float
    viscosity_mpa_s: tuple
    surface_tension_mn_m: float
    tier: str
    tweak: str
    source: str
    note: str


def _row(entry):
    for key in ("id", "name", "kind", "melting_k", "boiling_k", "vaporization_kj_mol", "density_kg_m3", "expansion_per_k", "viscosity_mpa_s",
                "surface_tension_mn_m", "tier", "tweak", "source"):
        if key not in entry:
            raise SurfaceLiquidError(f"{entry.get('id', entry)!r}: the row has no {key!r}")
    name = entry["id"]
    if entry["kind"] not in KINDS or entry["tier"] not in TIERS or entry["tweak"] not in TWEAKS or not entry["source"]:
        raise SurfaceLiquidError(f"{name}: a kind from {KINDS}, a tier T1..T9, a tweak L0..L4 and a source are required, got {entry!r}")
    density = entry["density_kg_m3"]
    points = tuple((float(t), float(v)) for t, v in entry["viscosity_mpa_s"])
    if not points or any(v <= 0.0 for _, v in points) or not density["value"] > 0.0 or (entry.get("saturation") is None and not entry["boiling_k"] > entry["melting_k"]):
        raise SurfaceLiquidError(f"{name}: positive density and viscosity and a boiling point above the melting point (for a row with no saturation curve) are required")
    formula = entry.get("saturation")
    if formula is not None and not saturation.has(formula):
        raise SurfaceLiquidError(f"{name}: {formula!r} has no saturation curve")
    return Liquid(name, entry["name"], entry["kind"], formula, float(entry["melting_k"]), float(entry["boiling_k"]), float(entry["vaporization_kj_mol"]),
                  entry.get("critical_k"), float(density["at_k"]), float(density["value"]), float(entry["expansion_per_k"]), points,
                  float(entry["surface_tension_mn_m"]), entry["tier"], entry["tweak"], entry["source"], entry.get("note", ""))


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    rows = tuple(_row(entry) for entry in document["liquids"])
    if len({row.id for row in rows}) != len(rows):
        raise SurfaceLiquidError("a liquid id is used twice")
    return document["version"], rows


VERSION, LIQUIDS = _load()
_BY_ID = {row.id: row for row in LIQUIDS}


def liquid(liquid_id):
    """The row `liquid_id`; raises `SurfaceLiquidError`."""
    try:
        return _BY_ID[liquid_id]
    except KeyError:
        raise SurfaceLiquidError(f"no liquid {liquid_id!r} in physics/data/surface_liquids.json") from None


def by_kind(kind):
    """The rows of one kind, in file order."""
    return tuple(row for row in LIQUIDS if row.kind == kind)


# --- the window ---------------------------------------------------------------------------------------------------------

def boiling_k(liquid_id, pressure_pa):
    """
    The temperature at which `liquid_id` boils at `pressure_pa`, K, capped at its critical temperature; `None` when the liquid cannot exist at that pressure
    (a row with a saturation curve and a pressure under its triple pressure: it sublimes instead).
    """
    row = liquid(liquid_id)
    if not pressure_pa > 0.0:
        raise ValueError("a pressure must be positive")
    if row.saturation is not None:
        triple = saturation.triple_point_k(row.saturation)
        if saturation.pressure_pa(row.saturation, triple + 1e-6) > pressure_pa:      # under the triple pressure the solid goes straight to vapour
            return None
        boil = saturation.condensation_temperature_k(row.saturation, pressure_pa)
        if boil is None:
            boil = row.critical_k
    else:
        inverse = 1.0 / row.boiling_k - GAS_CONSTANT_J_MOL_K * math.log(pressure_pa / ONE_BAR_PA) / (row.vaporization_kj_mol * 1000.0)
        boil = 1.0 / inverse if inverse > 0.0 else math.inf
    return min(boil, row.critical_k) if row.critical_k else boil


def melting_k(liquid_id):
    """The melting point (the liquidus of a lava), K; the table's value, not a function of pressure to first order."""
    row = liquid(liquid_id)
    if row.saturation is not None:
        return saturation.triple_point_k(row.saturation)
    return row.melting_k


def window_k(liquid_id, pressure_pa):
    """`(melting, boiling)` of `liquid_id` at `pressure_pa`, K, or `None` when the liquid is not stable at any temperature at that pressure."""
    boil = boiling_k(liquid_id, pressure_pa)
    low = melting_k(liquid_id)
    if boil is None or not boil > low:
        return None
    return low, boil


def is_stable(liquid_id, temperature_k, pressure_pa):
    """Whether `liquid_id` is a liquid at `temperature_k` and `pressure_pa`."""
    window = window_k(liquid_id, pressure_pa)
    return window is not None and window[0] <= temperature_k <= window[1]


def stable_liquids(temperature_k, pressure_pa, kinds=None):
    """The ids of the liquids stable at the temperature and pressure (those of `kinds`, if given), in file order."""
    return tuple(row.id for row in LIQUIDS if (kinds is None or row.kind in kinds) and is_stable(row.id, temperature_k, pressure_pa))


# --- density and viscosity -------------------------------------------------------------------------------------------------

def density_kg_m3(liquid_id, temperature_k):
    """The density of a liquid, kg/m3: the reference value, falling by `expansion_per_k` per kelvin, held to `MIN_DENSITY_SHARE` of it."""
    row = liquid(liquid_id)
    scale = 1.0 - row.expansion_per_k * (temperature_k - row.density_at_k)
    return row.density_kg_m3 * max(scale, MIN_DENSITY_SHARE)


def viscosity_pa_s(liquid_id, temperature_k):
    """The dynamic viscosity, Pa s: linear in the logarithm between the row's points, held at the end ones."""
    points = sorted(liquid(liquid_id).viscosity_mpa_s)
    if temperature_k <= points[0][0]:
        value = points[0][1]
    elif temperature_k >= points[-1][0]:
        value = points[-1][1]
    else:
        value = points[-1][1]
        for (t0, v0), (t1, v1) in zip(points, points[1:]):
            if temperature_k <= t1:
                share = (temperature_k - t0) / (t1 - t0)
                value = math.exp(math.log(v0) + share * (math.log(v1) - math.log(v0)))
                break
    return value * 1e-3


# --- what a world has ------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Candidate:
    """A liquid a world could hold: its id, `basis` (`"condensing"` for a volatile at saturation in the air, `"reservoir"` for one supplied by the ice or
    ore the body holds) and `score` (the partial pressure over the saturation pressure, or the reservoir's amount)."""

    liquid_id: str
    basis: str
    score: float


def condensing(temperature_k, partial_pressures_pa):
    """
    The volatiles whose partial pressure in the air has reached the saturation pressure of their liquid at `temperature_k`, as `{formula: ratio}` (partial over
    saturation, 1 or more): those condense into a liquid. A species under its triple point or over its critical temperature, or with no liquid row, is left out.
    """
    found = {}
    for formula, pressure in partial_pressures_pa.items():
        row = next((r for r in LIQUIDS if r.saturation == formula), None)
        if row is None or not pressure > 0.0:
            continue
        low = saturation.triple_point_k(formula)
        if not low < temperature_k <= (row.critical_k or math.inf):
            continue
        saturated = saturation.pressure_pa(formula, temperature_k)
        if math.isfinite(saturated) and pressure >= saturated:
            found[formula] = pressure / saturated
    return found


def candidates(temperature_k, pressure_pa, partial_pressures_pa=None, reservoir=None):
    """
    The liquids a world at `temperature_k` and `pressure_pa` could hold, best first: those stable there (`is_stable`) that a volatile at saturation in the air
    supplies (`condensing`, ranked by the ratio) or that the body's reservoir supplies (`reservoir`, `{liquid id: amount}`, ranked by amount, after the
    condensing ones). The amounts are the caller's (an ice fraction of the body, a share of its ore).
    """
    out = []
    for formula, ratio in condensing(temperature_k, partial_pressures_pa or {}).items():
        row = next(r for r in LIQUIDS if r.saturation == formula)
        if is_stable(row.id, temperature_k, pressure_pa):
            out.append(Candidate(row.id, "condensing", ratio))
    out.sort(key=lambda c: -c.score)
    supplied = [Candidate(liquid_id, "reservoir", amount) for liquid_id, amount in (reservoir or {}).items()
                if amount > 0.0 and is_stable(liquid_id, temperature_k, pressure_pa) and all(c.liquid_id != liquid_id for c in out)]
    supplied.sort(key=lambda c: -c.score)
    return tuple(out + supplied)


def dominant(temperature_k, pressure_pa, partial_pressures_pa=None, reservoir=None):
    """The best of `candidates`, or `None` when the world holds no surface liquid."""
    found = candidates(temperature_k, pressure_pa, partial_pressures_pa, reservoir)
    return found[0] if found else None


def ideal_sea(temperature_k, partial_pressures_pa, solvent="C2H6"):
    """
    The mole fractions of a sea in equilibrium with an air, by Raoult's law: every volatile `x = p / P_sat(T)`, and the involatile `solvent` (ethane, the product
    of methane's photolysis) taking the rest. `{}` when the volatiles alone fill the sea less than it holds (`sum x < 1`) and there is no solvent to add,
    or when they overfill it (`sum x > 1`: the air would be condensing more than a sea can hold, and the volatiles are scaled to sum to 1).
    """
    fractions = {}
    for formula, pressure in partial_pressures_pa.items():
        if formula == solvent or not pressure > 0.0 or not saturation.has(formula):
            continue
        low = saturation.triple_point_k(formula)
        row = next((r for r in LIQUIDS if r.saturation == formula), None)
        if row is None or not low < temperature_k <= (row.critical_k or math.inf):
            continue
        saturated = saturation.pressure_pa(formula, temperature_k)
        if math.isfinite(saturated) and saturated > 0.0:
            fractions[formula] = min(SOLUBILITY_FACTOR.get(formula, 1.0) * pressure / saturated, 1.0)
    total = sum(fractions.values())
    if total > 1.0:
        return {formula: x / total for formula, x in fractions.items()}
    if solvent is None or not any(r.saturation == solvent for r in LIQUIDS):
        return {formula: x / total for formula, x in fractions.items()} if total > 0.0 else {}
    rest = 1.0 - total
    if rest > 1e-9:
        fractions[solvent] = rest
    return fractions


def sea_density_kg_m3(fractions, temperature_k):
    """The density of a mixed sea, kg/m3: the mole-weighted molar volume of its pure liquids (ideal mixing; no excess volume)."""
    if not fractions:
        raise ValueError("a sea needs at least one liquid")
    mass = volume = 0.0
    for formula, x in fractions.items():
        row = next((r for r in LIQUIDS if r.saturation == formula), None)
        if row is None:
            raise SurfaceLiquidError(f"{formula!r} has no liquid row")
        molar_mass = gases.gas(formula).molar_mass
        mass += x * molar_mass
        volume += x * molar_mass / density_kg_m3(row.id, temperature_k)
    return mass / volume
