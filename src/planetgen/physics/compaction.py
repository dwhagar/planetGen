# planetgen/physics/compaction.py

"""
The mixing rule and the compaction of small bodies (GEN.442, docs/design/small-bodies-composition.md section 4)
=================================================================================================================

A small body is a heap of grains with voids between them. Its bulk density is the zero-pressure mixing rule of its grains
with the void share taken out:

    rho_grain = 1 / sum(w_i / rho_i),      rho = rho_grain (1 - porosity)

with `w_i` the mass fractions (adding to 1) and `rho_i` the grain densities of the shared materials table (GEN.427; a
caller may state a density the table does not hold, for example the silicate dust of a comet).

How much of the void a body keeps depends on its size. Bierson and Nimmo (2019, Icarus 326, 10, Table 1 and Eqs. 1 and 2)
give the brittle law, with the pressure `P` in megapascals: ice grains close their pores as `max(0.60 exp(-0.1 P), 0.20)` and
silicate grains as `min(0.6 P^-0.11, 0.6)`. The paper's own text says the law is a lower bound on the density above about 200
to 500 km radius, where viscous compaction decides, and the law alone leaves the density 1.3 to 1.8 times low up to 3e23 kg
(design note section 4.2). Boss's decision S9 (card tapped 8:13 PM EDT on 2026-10-21): the porosity is multiplied by a weight
that falls from 1 at `TAPER_START_KG` (1e19) to 0 at `TAPER_END_KG` (1e21) in log mass, so the density is continuous in mass
and a body above 1e21 kg has no macroporosity; the endpoints are flagged stand-ins until viscous compaction replaces them.

The reference integration (a uniform sphere in hydrostatic balance, grains compressed by third-order Birch-Murnaghan
equations of state, the porosity laws mixed by grain volume) is `scripts/build_compaction_table.py`; it writes
`physics/data/compaction_table.json`, the radius and bulk density of an ice-and-rock body by mass and rock fraction, which
`compacted_density_kg_m3` interpolates (the fast path of design note section 4.3; a test compares it with the integration).
The table covers 1e12 to 1e21 kg and rock fractions 0.30 to 0.95. Above 1e21 kg (`OutOfRange`) the structure solver of GEN.429
takes over with the same grains.

Use:

    grain_density_kg_m3({"water ice": 0.3, "forsterite": 0.7})          # kg/m3 of the grains; times (1 - porosity) for the heap
    compacted_density_kg_m3(radius_km=400.0, rock_fraction=0.7)        # kg/m3, the tapered brittle law
    compacted_radius_km(mass_kg=6.3e20, rock_fraction=0.7)
"""

import bisect
import json
import math
import os

from planetgen.physics import elements, materials
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "compaction_table.json")
"""str: The data file `scripts/build_compaction_table.py` writes."""

ICE_GRAIN_KG_M3 = 920.0
SILICATE_GRAIN_KG_M3 = 3500.0
"""float: The grain densities of Bierson and Nimmo (2019) Table 1."""

TAPER_START_KG = 1e19
TAPER_END_KG = 1e21
"""float: Where the porosity starts to close and where it is closed (stand-ins, decision S9)."""

SOURCE_STATE = "stand-in"
"""str: The taper endpoints are a stand-in until viscous compaction (Bierson and Nimmo Eq. 3 with its sign checked) replaces them."""


class OutOfRange(ValueError):
    """A body outside the table: heavier than `TAPER_END_KG`, where the structure solver (GEN.429) decides."""






def taper_weight(mass_kg):
    """1 at or below `TAPER_START_KG`, 0 at or above `TAPER_END_KG`, linear in log10 of the mass between."""
    if mass_kg <= 0.0:
        raise ValueError("the mass must be positive")
    span = math.log10(TAPER_END_KG) - math.log10(TAPER_START_KG)
    return min(max((math.log10(TAPER_END_KG) - math.log10(mass_kg)) / span, 0.0), 1.0)


def grain_density_kg_m3(fractions, densities=None):
    """
    The zero-pressure density of the grains, 1 / sum(w_i / rho_i).

    Args:
        fractions (dict): Component name (a row of the materials table) -> mass fraction; scaled to add to 1.
        densities (dict or None): Component name -> kg/m3 for a component to take from here and not from the table.

    Raises:
        ValueError: No fraction, a negative one, or a component with no grain density (a group row has none).
        KeyError: A component that is neither in `densities` nor a row of the table.
    """
    if not fractions:
        raise ValueError("no components")
    densities = densities or {}
    total = math.fsum(fractions.values())
    if not total > 0.0 or any(value < 0.0 for value in fractions.values()):
        raise ValueError(f"the mass fractions {fractions!r} are negative or all zero")
    volume = 0.0
    for name, fraction in fractions.items():
        if fraction == 0.0:
            continue
        if name in densities:
            rho = float(densities[name])
        else:
            row = materials.material(name)
            try:
                rho = row.prop("rho_grain_g_cm3").value * 1000.0
            except KeyError:
                raise ValueError(f"the material {name!r} has no grain density in the materials table; state one") from None
            if rho is None:
                raise ValueError(f"the material {name!r} has no grain density in the materials table; state one")
        volume += fraction / total / rho
    return 1.0 / volume




def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    for key in ("title", "method", "sources", "taper_start_kg", "taper_end_kg", "log10_mass_kg", "columns"):
        if not data.get(key):
            raise DataError(f"compaction_table.json states no {key}")
    if data["taper_start_kg"] != TAPER_START_KG or data["taper_end_kg"] != TAPER_END_KG:
        raise DataError("compaction_table.json was built with other taper endpoints than the module states")
    log_mass = data["log10_mass_kg"]
    if any(b <= a for a, b in zip(log_mass, log_mass[1:])):
        raise DataError("the masses of compaction_table.json do not rise")
    columns = []
    for column in data["columns"]:
        radius, density = column["radius_km"], column["density_kg_m3"]
        if not (len(radius) == len(density) == len(log_mass)):
            raise DataError(f"rock fraction {column.get('rock_fraction')}: the lengths do not match the masses")
        if any(b <= a for a, b in zip(radius, radius[1:])):
            raise DataError(f"rock fraction {column['rock_fraction']}: the radii do not rise with the mass")
        columns.append((float(column["rock_fraction"]), tuple(math.log(r) for r in radius), tuple(math.log(d) for d in density)))
    if len(columns) < 2 or any(b[0] <= a[0] for a, b in zip(columns, columns[1:])):
        raise DataError("the rock fractions of compaction_table.json are fewer than two or do not rise")
    return data, tuple(log_mass), tuple(columns)


_DATA, _LOG_MASS, _COLUMNS = _load()
SOURCES = _DATA["sources"]
"""dict: Source key -> citation, as the data file states them."""
ROCK_FRACTION_RANGE = (_COLUMNS[0][0], _COLUMNS[-1][0])
"""tuple: The rock fractions the table covers."""
MASS_RANGE_KG = (10.0 ** _LOG_MASS[0], 10.0 ** _LOG_MASS[-1])
"""tuple: The masses the table covers."""


def _check_rock(rock_fraction):
    low, high = ROCK_FRACTION_RANGE
    if not low - 1e-12 <= rock_fraction <= high + 1e-12:
        raise ValueError(f"the rock fraction {rock_fraction!r} is outside the table's {low} to {high}")


def _blend(rock_fraction):
    """The two columns around `rock_fraction` and the weight of the upper one."""
    rocks = [c[0] for c in _COLUMNS]
    i = min(max(bisect.bisect_right(rocks, rock_fraction) - 1, 0), len(_COLUMNS) - 2)
    low, high = _COLUMNS[i], _COLUMNS[i + 1]
    return low, high, (rock_fraction - low[0]) / (high[0] - low[0])


def _interp(xs, ys, x):
    i = min(max(bisect.bisect_right(xs, x) - 1, 0), len(xs) - 2)
    return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i])


def max_radius_km(rock_fraction):
    """The radius of a body of `TAPER_END_KG` at this rock fraction: the largest the table takes."""
    _check_rock(rock_fraction)
    low, high, t = _blend(rock_fraction)
    return math.exp(low[1][-1] * (1.0 - t) + high[1][-1] * t)


def compacted_density_kg_m3(radius_km, rock_fraction):
    """
    The bulk density of an ice-and-rock body of this radius, from the tapered brittle law of Bierson and Nimmo (the table of
    `scripts/build_compaction_table.py`). Below the table's smallest radius (about 0.6 km) the density is that of zero
    pressure, which is where the table starts.

    Raises:
        OutOfRange: The radius is beyond that of a 1e21 kg body of this mix, where the structure solver decides.
        ValueError: A rock fraction outside the table, or a radius that is not positive.
    """
    if not radius_km > 0.0:
        raise ValueError("the radius must be positive")
    _check_rock(rock_fraction)
    low, high, t = _blend(rock_fraction)
    log_r = math.log(radius_km)
    if log_r > low[1][-1] * (1.0 - t) + high[1][-1] * t + 1e-12:
        raise OutOfRange(f"a radius of {radius_km} km is heavier than {TAPER_END_KG:.0e} kg at rock fraction {rock_fraction}: "
                         f"the structure solver decides")
    values = []
    for column in (low, high):
        x = min(max(log_r, column[1][0]), column[1][-1])
        values.append(_interp(column[1], column[2], x))
    return math.exp(values[0] * (1.0 - t) + values[1] * t)


def compacted_radius_km(mass_kg, rock_fraction):
    """The radius of an ice-and-rock body of this mass (see `compacted_density_kg_m3`).

    Raises:
        OutOfRange: The mass is outside the table's `MASS_RANGE_KG`.
    """
    _check_rock(rock_fraction)
    if not MASS_RANGE_KG[0] <= mass_kg <= MASS_RANGE_KG[1] * (1.0 + 1e-12):
        raise OutOfRange(f"a mass of {mass_kg} kg is outside the table's {MASS_RANGE_KG[0]:.0e} to {MASS_RANGE_KG[1]:.0e} kg")
    low, high, t = _blend(rock_fraction)
    log_m = math.log10(mass_kg)
    values = [_interp(_LOG_MASS, column[1], log_m) for column in (low, high)]
    return math.exp(values[0] * (1.0 - t) + values[1] * t)
