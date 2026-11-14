# planetgen/physics/remnant_composition.py

"""
What a supernova remnant is made of: the ejecta of its progenitor and the interstellar gas it has swept up (GEN.1271)
===================================================================================================================

A remnant of classes R to W holds two things, reported apart and as a total (docs/design/composition/ and the GEN.1268
research note, section 6): the **ejecta** of the star that exploded, and the **swept-up** interstellar gas inside the
blast wave. This module works both out and returns the values; nothing is derived at read time (Boss, 8:49 PM EST on
2026-11-04: the composition is stored, GEN.1276), so the result carries the **progenitor's initial mass** too.

* **Swept-up mass.** `swept_up_mass_sol`: 4/3 pi r^3 n (1.4 m_H) from the stored radius and density. At the age and density
  of Cassiopeia A (350 years, about 3 particles per cm^3) this gives about 7 solar masses against the 7.9 the XMM-Newton budget
  finds [R5097] (tests/test_remnant_composition.py).
* **Type Ia (class W).** The ejecta are the thermonuclear model W7: 1.37 solar masses of the elements of [R5096: Table 3],
  mostly iron-group and intermediate-mass, no hydrogen and no helium. No core is left. The progenitor is the white dwarf's
  main-sequence star, drawn from the IMF over `type_ia_progenitor_range_sol` (a stand-in, SNC5).
* **Core-collapse (R, S, T, U, V).** The progenitor's initial mass is drawn from the IMF over the windows of its stored core
  (`remnant_progenitors.kind_windows`: a neutron star's or a black hole's; any core-collapse mass when none is stored). Its
  mass at death is its initial mass less its wind (`wind_loss_sol`) when it keeps its hydrogen envelope, or the mass of its helium
  core (`helium_core_mass_sol`) when it is **stripped**, which every progenitor of 40 solar masses or more is [R4251] and a third of
  the lighter ones are (a binary partner takes the envelope, SNC4). The ejecta are that mass less the compact core
  (`remnant_mass.draw_mass_sol`). The elements the explosion made, [R5096: Table 3]'s mean over 10 to 50 solar masses, are
  scaled with the progenitor's mass so that the mean over the IMF in that window is the table's (`yield_scale`): oxygen,
  neon, magnesium as M^2.5, the other light metals as M^1.5, the iron group not at all. The rest of the inner ejecta is
  helium and, for an unstripped star, the envelope keeps the progenitor's birth hydrogen, helium and metals
  (`stellar_abundances.composition` at its [Fe/H]).
* **The swept-up gas** has the composition of the place (`ambient`). Until the shared nebula model (GEN.1269) gives it, the
  caller passes the gas's mass fractions, or the proto-solar mix scaled to `fe_h` is used.

Every draw is keyed on the object (`key`, its ID), in streams of their own, so a remnant reads the same every time and no
existing stream moves. The constants are in `physics/data/remnant_yields.json` with their source and, for the stand-ins,
their registry rows (SNC1, SNC4 to SNC7, gap-filling-plan.md section 6.177); the presupernova table is [R4251: Table 2].
"""

import functools
import json
import math
import os
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import elements, imf, remnant_mass, remnant_progenitors, stellar_abundances
from planetgen.physics import unit_registry as units
from planetgen.util import draw

_g_to_msun = units.converter("g", "Msun")

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "remnant_yields.json")

TYPE_IA = "Type Ia"
CORE_COLLAPSE = "core-collapse"

MEAN_PARTICLE_MASS_OVER_HYDROGEN = 1.4
"""float: The mass of the swept-up gas per hydrogen nucleus, in hydrogen masses (the helium and the metals add 40 percent)."""

_PC_PER_LY = 1.0 / 3.261563777
_CM_PER_PC = 3.0856775814913673e18
_HYDROGEN_MASS_G = 1.6735575e-24


class RemnantCompositionError(ValueError):
    """A remnant's kind or core is not one this model knows."""


@functools.lru_cache(maxsize=1)
def _data():
    with open(DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def _stand_ins():
    return _data()["stand_ins"]


def table_elements():
    """`{symbol: {"type_ii_10_50": solar masses, "w7": solar masses}}`: the sourced table [R5096: Table 3] by element."""
    return _data()["elements_solar_masses"]


# --- Swept-up gas ------------------------------------------------------------------------------------------------------

def swept_up_mass_sol(radius_ly, density_cm3):
    """The interstellar mass, solar masses, inside a remnant of radius `radius_ly` in gas of `density_cm3` hydrogen nuclei per cm^3:
    the sphere's volume times the density times 1.4 hydrogen masses."""
    radius_cm = radius_ly * _PC_PER_LY * _CM_PER_PC
    volume_cm3 = 4.0 / 3.0 * math.pi * radius_cm ** 3
    grams = volume_cm3 * density_cm3 * MEAN_PARTICLE_MASS_OVER_HYDROGEN * _HYDROGEN_MASS_G
    return _g_to_msun(grams)


# --- The progenitor ----------------------------------------------------------------------------------------------------

def _presupernova_rows():
    return _data()["presupernova_table"]["rows"]


def _along_table(initial_mass_sol, column):
    """The table's column (1 final mass, 2 helium core) at an initial mass: linear between its rows, held at the heaviest row, and under
    the lightest row the power law of the first two rows."""
    rows = _presupernova_rows()
    if initial_mass_sol >= rows[-1][0]:
        return rows[-1][column]
    if initial_mass_sol < rows[0][0]:
        exponent = math.log(rows[1][column] / rows[0][column]) / math.log(rows[1][0] / rows[0][0])
        return rows[0][column] * (initial_mass_sol / rows[0][0]) ** exponent
    for (m0, *v0), (m1, *v1) in zip(rows, rows[1:]):
        if m0 <= initial_mass_sol <= m1:
            return v0[column - 1] + (v1[column - 1] - v0[column - 1]) * (initial_mass_sol - m0) / (m1 - m0)
    raise AssertionError("the rows cover every mass")


def helium_core_mass_sol(initial_mass_sol):
    """The mass, solar masses, of a core-collapse progenitor that has lost its hydrogen envelope: the helium core of [R4251: Table 2], linear
    between its rows (10, 12, 15, 20, 25, 30 solar masses, then the bare stars of 60 and more; the 30 to 60 span is interpolated)."""
    return _along_table(initial_mass_sol, 2)


def wind_loss_sol(initial_mass_sol):
    """The mass, solar masses, a progenitor loses in its wind before it explodes: its initial mass less the final mass of [R4251: Table 2],
    linear between the rows; under 10 solar masses the share the 10 solar mass model loses."""
    rows = _presupernova_rows()
    if initial_mass_sol < rows[0][0]:
        return initial_mass_sol * (1.0 - rows[0][1] / rows[0][0])
    return initial_mass_sol - _along_table(initial_mass_sol, 1)


def is_stripped(initial_mass_sol, u_binary):
    """Whether the progenitor explodes without its hydrogen envelope: always above `stripped_above_mass_sol`, else when the uniform
    draw `u_binary` falls under `binary_stripped_share` (a companion took the envelope)."""
    s = _stand_ins()
    return initial_mass_sol >= s["stripped_above_mass_sol"] or u_binary < s["binary_stripped_share"]


def progenitor_mass_sol(key, progenitor_type, compact_kind, failed_share=None):
    """
    The initial mass, solar masses, of a remnant's progenitor, drawn from the IMF in a stream of the object's own (`key`).
    A Type Ia remnant's is the white dwarf's main-sequence star (`type_ia_progenitor_range_sol`); a core-collapse one's lies
    in the windows of its stored core (`compact_kind`: `"neutron_star"`, `"black_hole"` or `None` for any core-collapse mass).
    """
    rng = draw.Stream(f"remnant-composition:{key}:progenitor")
    if progenitor_type == TYPE_IA:
        low, high = _stand_ins()["type_ia_progenitor_range_sol"]
        return imf.sample(low, high, rng)
    if progenitor_type != CORE_COLLAPSE:
        raise RemnantCompositionError(f"a progenitor is {TYPE_IA!r} or {CORE_COLLAPSE!r}, not {progenitor_type!r}")
    kinds = {"neutron_star": remnant_progenitors.NEUTRON_STAR, "black_hole": remnant_progenitors.BLACK_HOLE}
    if compact_kind is None:
        return imf.sample(*tuning.CORE_COLLAPSE_MASS_RANGE_SOL, rng)
    if compact_kind not in kinds:
        raise RemnantCompositionError(f"a core is 'neutron_star', 'black_hole' or None, not {compact_kind!r}")
    mass = remnant_progenitors.sample_mass_sol(kinds[compact_kind], rng=rng, failed_share=failed_share)
    return mass if mass is not None else imf.sample(*tuning.CORE_COLLAPSE_MASS_RANGE_SOL, rng)


# --- What the explosion made -------------------------------------------------------------------------------------------

def _exponent_of(symbol):
    s = _stand_ins()
    return s["yield_exponent"][s["yield_group"].get(symbol, "iron_group")]


@functools.lru_cache(maxsize=None)
def _mean_power(exponent):
    """The IMF-weighted mean of M^exponent over the table's mass window."""
    low, high = _stand_ins()["table_mass_window_sol"]
    return imf.integral(exponent, low, high) / imf.integral(0.0, low, high)


def yield_scale(symbol, initial_mass_sol):
    """What a progenitor of `initial_mass_sol` makes of `symbol` in units of the table's mean (the IMF-weighted mean of this
    over the table's mass window is 1)."""
    exponent = _exponent_of(symbol)
    return initial_mass_sol ** exponent / _mean_power(exponent)


def processed_masses(initial_mass_sol):
    """`{symbol: solar masses}` of the elements a core-collapse explosion of a progenitor of `initial_mass_sol` made: the table's
    mean over 10 to 50 solar masses times `yield_scale`."""
    return {symbol: row["type_ii_10_50"] * yield_scale(symbol, initial_mass_sol) for symbol, row in table_elements().items()}


def type_ia_masses():
    """`{symbol: solar masses}` of a Type Ia explosion's ejecta: model W7 of [R5096: Table 3], 1.37 solar masses in all."""
    masses = {"H": 0.0, "He": 0.0}
    masses.update({symbol: row["w7"] for symbol, row in table_elements().items()})
    return masses


# --- The whole remnant -------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class RemnantComposition:
    """
    A remnant's composition.

    Attributes:
        progenitor_mass_sol (float): The progenitor's initial mass (stored with the remnant, GEN.1276).
        stripped (bool): Whether the progenitor exploded without its hydrogen envelope (always true for Type Ia, which has none).
        compact_mass_sol (float): The mass of the core left behind (0 for Type Ia).
        ejecta_mass_sol (float): The mass thrown out.
        ejecta (dict): Element (`"H"`, `"He"`, the table's symbols, `"other metals"`) to solar masses; sums to `ejecta_mass_sol`.
        swept_mass_sol (float): The interstellar mass inside the blast wave.
        swept (dict): Element to solar masses of the swept-up gas; sums to `swept_mass_sol`.
    """

    progenitor_mass_sol: float
    stripped: bool
    compact_mass_sol: float
    ejecta_mass_sol: float
    ejecta: dict
    swept_mass_sol: float
    swept: dict

    @property
    def total_mass_sol(self):
        return self.ejecta_mass_sol + self.swept_mass_sol

    @property
    def ejecta_share(self):
        """The ejecta's share of the mass inside the remnant's radius."""
        total = self.total_mass_sol
        return self.ejecta_mass_sol / total if total > 0 else 0.0

    def fractions(self):
        """`(ejecta, swept, total)`: element mass fractions of each part and of the two together (each sums to 1)."""
        def of(masses, mass):
            return {symbol: value / mass for symbol, value in masses.items()} if mass > 0 else {}
        total = {}
        for part in (self.ejecta, self.swept):
            for symbol, value in part.items():
                total[symbol] = total.get(symbol, 0.0) + value
        return of(self.ejecta, self.ejecta_mass_sol), of(self.swept, self.swept_mass_sol), of(total, self.total_mass_sol)


def ambient_masses(swept_mass_sol, fe_h=0.0, ambient=None):
    """
    The swept-up gas split by element: `ambient` (element to mass fraction, summing to 1; the nebula model's gas, GEN.1269) if
    given, else the birth composition of gas of metallicity `fe_h` (`stellar_abundances.composition`), with the metals the table
    does not list in `"other metals"`.
    """
    fractions = dict(ambient) if ambient is not None else dict(stellar_abundances.composition(fe_h, {}).elements)
    listed = set(table_elements()) | {"H", "He"}
    out = {"H": 0.0, "He": 0.0}
    other = 0.0
    for symbol, value in fractions.items():
        if symbol in listed:
            out[symbol] = out.get(symbol, 0.0) + value
        else:
            other += value
    if other > 0:
        out["other metals"] = other
    scale = swept_mass_sol / sum(out.values())
    return {symbol: value * scale for symbol, value in out.items()}


def _core_collapse(initial_mass_sol, compact_kind, stripped, fe_h, rng):
    """`(compact_mass_sol, ejecta masses)` of a core-collapse explosion."""
    s = _stand_ins()
    final = helium_core_mass_sol(initial_mass_sol) if stripped else initial_mass_sol - wind_loss_sol(initial_mass_sol)
    final = min(final, initial_mass_sol)
    compact = remnant_mass.draw_mass_sol("black-hole" if compact_kind == "black_hole" else "neutron-star", rng=rng)
    compact = min(compact, (1.0 - s["least_ejecta_share"]) * final)
    ejecta_total = final - compact
    core = min(helium_core_mass_sol(initial_mass_sol), final)
    inner = max(core - compact, 0.0)                 # the burned shells, thrown out
    envelope = max(final - core, 0.0)                # the hydrogen envelope, for an unstripped star
    processed = processed_masses(initial_mass_sol)
    cap = s["processed_share_cap"] * max(inner, 1e-9)
    processed_total = sum(processed.values())
    if processed_total > cap:
        factor = cap / processed_total
        processed = {symbol: value * factor for symbol, value in processed.items()}
        processed_total = cap
    birth = stellar_abundances.composition(fe_h, {})
    listed = set(table_elements())
    ejecta = {"H": 0.0, "He": inner - processed_total}
    for symbol, value in processed.items():
        ejecta[symbol] = ejecta.get(symbol, 0.0) + value
    for symbol, fraction in birth.elements.items():
        mass = fraction * envelope
        if symbol in listed or symbol in ("H", "He"):
            ejecta[symbol] = ejecta.get(symbol, 0.0) + mass
        else:
            ejecta["other metals"] = ejecta.get("other metals", 0.0) + mass
    # The inner burned shells also carry the star's own unburned metals: left out, they are inside `processed` already.
    scale = ejecta_total / sum(ejecta.values())
    return compact, {symbol: value * scale for symbol, value in ejecta.items()}


def remnant_composition(key, progenitor_type, compact_kind, radius_ly, density_cm3, fe_h=0.0, ambient=None, failed_share=None):
    """
    The composition of a supernova remnant: ejecta, swept-up gas and the progenitor's mass.

    Args:
        key (str): The object's ID (every draw is keyed on it, in streams of this module's own).
        progenitor_type (str): `"Type Ia"` or `"core-collapse"`.
        compact_kind (str): The stored core: `"neutron_star"`, `"black_hole"` or `None` (none stored, or a Type Ia).
        radius_ly (float): The remnant's stored radius, light-years.
        density_cm3 (float): The ambient density, hydrogen nuclei per cm^3.
        fe_h (float): The [Fe/H] of the progenitor's birth gas and, when `ambient` is not given, of the swept-up gas.
        ambient (dict): The swept-up gas's element mass fractions (GEN.1269's), or `None` for the birth composition at `fe_h`.
        failed_share (float): The galaxy's failed-explosion share, for a black hole's progenitor windows (the reference when `None`).

    Returns:
        RemnantComposition
    """
    mass = progenitor_mass_sol(key, progenitor_type, compact_kind, failed_share)
    swept_mass = swept_up_mass_sol(radius_ly, density_cm3)
    swept = ambient_masses(swept_mass, fe_h, ambient)
    if progenitor_type == TYPE_IA:
        ejecta = type_ia_masses()
        return RemnantComposition(mass, True, 0.0, sum(ejecta.values()), ejecta, swept_mass, swept)
    rng = draw.Stream(f"remnant-composition:{key}:explosion")
    stripped = is_stripped(mass, draw.Stream(f"remnant-composition:{key}:binary").random())
    compact, ejecta = _core_collapse(mass, compact_kind, stripped, fe_h, rng)
    return RemnantComposition(mass, stripped, compact, sum(ejecta.values()), ejecta, swept_mass, swept)
