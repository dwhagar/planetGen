# planetgen/physics/brown_dwarf_types.py

"""
Effective temperature against spectral type of brown dwarfs (GEN.517, docs/design/composition/brown-dwarf-types.md)
=====================================================================================================

The empirical map of [R2374: Table 13] (525 L, T and Y dwarfs within 20 pc), in the data file
`brown_dwarf_teff_map.json`: three polynomials in the spectral type `s` (L0 = 0, T0 = 10, Y0 = 20) with their scatter in Kelvin. The
type is a function of the effective temperature with that scatter, not a weighted draw; the soft windows of rows g1 to g5 of the
atmosphere recipes (`scripts/atmosphere_recipes_data.py`) are the same map seen from the temperature, and their edges are the
temperatures of the type boundaries below.

Use:

    teff_of_type(5.0)                    # 1613.5 K, an L5
    boundary_teff("t-y")                 # 460 K, the temperature of a Y0
    type_probabilities(700.0)            # {"T": ..., "late T": ..., "Y": ...} at a temperature, with the map's scatter
    number_of_teff(1000.0)               # 12.1, the mean relation read backwards (None above L0)
    drawn_number(1000.0, rng)            # the type of one dwarf at that temperature, the map's scatter drawn (GEN.760)

The generator reads `drawn_number` for the type of every brown dwarf (`brown_dwarf.spectral_type`, GEN.760); the soft
class windows of the atmosphere recipes read the same scatter through `type_probabilities`.
"""

import json
import math
import os
import statistics
from dataclasses import dataclass

from planetgen.physics import elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "brown_dwarf_teff_map.json")
"""str: The data file."""

_PHI = statistics.NormalDist().cdf


@dataclass(frozen=True)
class Piece:
    """One polynomial of the map: `Teff = sum(c_k s^k)` for `s_min <= s <= s_max`, with scatter `rms_k` Kelvin."""

    name: str
    s_min: float
    s_max: float
    coefficients: tuple
    rms_k: float

    def teff(self, s):
        return sum(c * s ** k for k, c in enumerate(self.coefficients))


@dataclass(frozen=True)
class Map:
    pieces: tuple
    boundaries: dict
    source: str


def load(path=DATA_FILE):
    """Reads and checks the map at `path`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the brown dwarf map {path}: {error}") from error
    pieces = []
    for raw in document.get("pieces") or []:
        coefficients = tuple(raw.get("coefficients") or ())
        if not coefficients or not all(isinstance(c, (int, float)) and math.isfinite(c) for c in coefficients):
            raise DataError(f"piece {raw.get('name')!r} has no finite coefficients")
        if not 0 < raw.get("rms_k", 0) < 1000 or not raw["s_min"] < raw["s_max"]:
            raise DataError(f"piece {raw.get('name')!r} has no scatter or no range")
        pieces.append(Piece(raw["name"], float(raw["s_min"]), float(raw["s_max"]), coefficients, float(raw["rms_k"])))
    if not pieces or any(a.s_max != b.s_min for a, b in zip(pieces, pieces[1:])):
        raise DataError("the pieces of the map are missing or not contiguous")
    boundaries = {b["id"]: b for b in document.get("boundaries") or []}
    if not boundaries:
        raise DataError("the map has no boundaries")
    for b in boundaries.values():
        if not pieces[0].s_min <= b["s"] <= pieces[-1].s_max:
            raise DataError(f"boundary {b['id']} is outside the map")
    return Map(tuple(pieces), boundaries, document.get("source", ""))


MAP = load()
"""Map: The empirical map."""


def piece_of(s):
    """The piece of the map that holds the spectral type `s`."""
    if not MAP.pieces[0].s_min <= s <= MAP.pieces[-1].s_max:
        raise DataError(f"spectral type {s} is outside the map ({MAP.pieces[0].s_min} to {MAP.pieces[-1].s_max})")
    for piece in MAP.pieces:
        if s <= piece.s_max:
            return piece
    return MAP.pieces[-1]


def teff_of_type(s):
    """The effective temperature in Kelvin of spectral type `s` (L0 = 0, T0 = 10, Y0 = 20)."""
    return piece_of(s).teff(s)


def scatter_k(s):
    """The scatter in Kelvin of the effective temperature at spectral type `s`."""
    return piece_of(s).rms_k


def boundary_teff(boundary_id):
    """The temperature of a type boundary (`m-l`, `l-t`, `t-late-t`, `t-y`)."""
    try:
        return teff_of_type(MAP.boundaries[boundary_id]["s"])
    except KeyError:
        raise DataError(f"boundary {boundary_id!r} is none of {', '.join(MAP.boundaries)}") from None


def boundary_scatter_k(boundary_id):
    """The scatter in Kelvin at a type boundary."""
    return scatter_k(MAP.boundaries[boundary_id]["s"])


def type_probabilities(teff_k):
    """
    `{class: probability}` of a brown dwarf at `teff_k`, the type being the map's with Gaussian scatter at the boundary: late M, L, T,
    late T and Y. The classes form a chain of boundaries (`m-l`, `l-t`, `t-late-t`, `t-y`); the probability of being cooler than a
    boundary is `Phi((boundary - teff) / scatter)`.
    """
    def cooler(boundary_id):
        """The chance the dwarf is cooler than the boundary, so of the class on the cool side of it or a cooler one."""
        return _PHI((boundary_teff(boundary_id) - teff_k) / boundary_scatter_k(boundary_id))
    m_l, l_t, t_late, t_y = (cooler(b) for b in ("m-l", "l-t", "t-late-t", "t-y"))
    return {"late M": 1.0 - m_l, "L": m_l - l_t, "T": l_t - t_late, "late T": t_late - t_y, "Y": t_y}


def number_of_teff(teff_k):
    """
    The spectral type number (L0 = 0, T0 = 10, Y0 = 20) whose mean temperature on the map is `teff_k`: the map read
    backwards. None when the temperature is that of L0 or hotter (a late-M dwarf, named by the stellar scale); colder than
    the last piece's end it is that end.
    """
    if teff_k >= teff_of_type(MAP.pieces[0].s_min):
        return None
    low, high = MAP.pieces[0].s_min, MAP.pieces[-1].s_max
    if teff_k <= teff_of_type(high):
        return high
    for _ in range(60):
        middle = 0.5 * (low + high)
        if teff_of_type(middle) > teff_k:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def drawn_number(teff_k, rng):
    """
    The type of one brown dwarf at `teff_k`, as `(number, temperature_k)`: the temperature it reads on the map is its own
    plus a Gaussian offset of the scatter of the piece its mean type is on ([R2374]: the dwarfs of one type
    spread by that many Kelvin), and the number is the map read backwards at that temperature. `number` is None when that
    temperature is the one of a late-M dwarf; `temperature_k` is the temperature the type was read at. The L/T transition,
    where the map is nearly flat, spreads a type over several subtypes, as the census does.
    """
    mean = number_of_teff(teff_k)
    if mean is None:
        scatter = MAP.pieces[0].rms_k
    else:
        scatter = scatter_k(mean)
    shifted = teff_k + rng.gauss(0.0, scatter)
    return number_of_teff(shifted), shifted
