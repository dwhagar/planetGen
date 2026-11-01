# planetgen/datasets/mamajek_sequence.py

"""
The mean dwarf spectral-type sequence of Pecaut and Mamajek (GEN.943): the effective temperature and the luminosity of a main-sequence
star of each spectral type, O3V to M9V.

The file is `EEM_dwarf_UBVIJHK_colors_Teff.txt`, version 2022.04.16, from Eric Mamajek's page at the University of Rochester
(https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt); much of it is Table 5 of [R3469], which is the paper to cite. Each row of its first table is a type (`G2V`, `K5V`, `M0.5V`) with the
effective temperature, the colours, and for most types the bolometric correction, luminosity and radius.

The generator reads the whole-number types only (O3 to M9, 67 rows): `spectral_letter_and_subclass` gives the type whose temperature is
nearest a star's, `temperature_k` the temperature of a type, and `dwarf_log_luminosity` the luminosity of a dwarf of a temperature (the
reference a subdwarf is under-luminous against). The sequence is of solar-metallicity dwarfs near the zero-age main sequence's middle
age; the hot end, O2 to O0, is not in it and is a stand-in of `tuning.HOT_O_SUBCLASS_TEMPERATURES_K`.

The readers raise `fetch.DataMissing` when the file has not been downloaded.
"""

import bisect
import functools
import math
import re
from dataclasses import dataclass

from planetgen import tuning
from planetgen.datasets import fetch

FILE = "EEM_dwarf_UBVIJHK_colors_Teff.txt"

DATASET = fetch.Dataset(
    key="mamajek_dwarf_sequence",
    title="Mamajek mean dwarf stellar colour and effective temperature sequence (2022.04.16)",
    source="[R3469], version 2022.04.16 of E. Mamajek's table",
    pieces=(fetch.Piece(FILE, "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt",
                        "1de2edeec17bb3346e0e4e70b999de5ee29947df38474e64cddb7cfacc164b7f", 55680),),
    replaces="the subclasses linear in temperature inside a letter, with a K/M boundary at 3,700 K",
    fill="PS12",
)

LETTERS = "OBAFGKM"
_ROW = re.compile(r"^([OBAFGKM])(\d)V\s+(\d+)\s")


@dataclass(frozen=True)
class SpectralType:
    """One whole-number type of the sequence: letter, subclass (0 hottest), effective temperature in kelvin and the log10 of the
    luminosity in solar units (None for a hot-end stand-in)."""

    letter: str
    subclass: int
    teff_k: float
    log_luminosity: float = None


def parse(text):
    """The whole-number dwarf types O3V to M9V of the first table of `text`, hottest first."""
    rows = []
    for line in text.splitlines():
        if line.startswith("#SpT") and rows:
            break  # the second table is a different one
        match = _ROW.match(line)
        if match:
            tokens = line.split()
            log_l = float(tokens[4]) if len(tokens) > 4 and tokens[4] not in ("...", ".....") else None
            rows.append(SpectralType(match.group(1), int(match.group(2)), float(match.group(3)), log_l))
    if len(rows) != 67:
        raise ValueError(f"the Mamajek table holds {len(rows)} whole-number types, not the 67 (O3 to M9) expected")
    return rows


@functools.lru_cache(maxsize=None)
def _load():
    """Every type, hottest first: the stand-in O0 to O2 rows (`tuning.HOT_O_SUBCLASS_TEMPERATURES_K`) and the table's O3V to M9V."""
    with open(fetch.locate(DATASET).require().path(FILE), encoding="utf-8") as handle:
        table = parse(handle.read())
    hot = [SpectralType("O", sub, kelvin) for sub, kelvin in sorted(tuning.HOT_O_SUBCLASS_TEMPERATURES_K.items())]
    return tuple(hot + table)


@functools.lru_cache(maxsize=None)
def _edges():
    """The temperatures where the nearest type changes (midpoints of neighbours), coolest first, with the types coolest first."""
    ordered = tuple(reversed(_load()))
    return tuple((a.teff_k + b.teff_k) / 2.0 for a, b in zip(ordered, ordered[1:])), ordered


def spectral_letter_and_subclass(temperature_k):
    """The letter and subclass of the type whose temperature is nearest `temperature_k` (O0 above the hottest, M9 below the coolest)."""
    edges, ordered = _edges()
    found = ordered[bisect.bisect_right(edges, temperature_k)]
    return found.letter, found.subclass


def temperature_k(letter, subclass):
    """The temperature of the type `letter` `subclass`."""
    for found in _load():
        if found.letter == letter and found.subclass == subclass:
            return found.teff_k
    raise KeyError(f"{letter}{subclass} is not a type of the sequence")


def letter_range_k(letter):
    """`(low, high)`: the temperatures between which the nearest type is of `letter` (the hottest O and the coolest M reach as
    far past their last type as the step to its neighbour)."""
    edges, ordered = _edges()
    inside = [i for i, found in enumerate(ordered) if found.letter == letter]
    low_i, high_i = min(inside), max(inside)
    low = edges[low_i - 1] if low_i > 0 else ordered[0].teff_k - (edges[0] - ordered[0].teff_k)
    high = edges[high_i] if high_i < len(edges) else ordered[-1].teff_k + (ordered[-1].teff_k - edges[-1])
    return low, high


@functools.lru_cache(maxsize=None)
def _luminosity_curve():
    rows = sorted((found.teff_k, found.log_luminosity) for found in _load() if found.log_luminosity is not None)
    return tuple(t for t, _ in rows), tuple(l for _, l in rows)


def dwarf_log_luminosity(teff_k):
    """The log10 luminosity (solar) of a solar-metallicity dwarf of `teff_k`, interpolated in log temperature along the sequence and
    held at its ends."""
    temperatures, logs = _luminosity_curve()
    if teff_k <= temperatures[0]:
        return logs[0]
    if teff_k >= temperatures[-1]:
        return logs[-1]
    i = bisect.bisect_right(temperatures, teff_k)
    t0, t1 = math.log(temperatures[i - 1]), math.log(temperatures[i])
    return logs[i - 1] + (logs[i] - logs[i - 1]) * (math.log(teff_k) - t0) / (t1 - t0)
