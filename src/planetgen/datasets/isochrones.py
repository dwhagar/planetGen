# planetgen/datasets/isochrones.py

"""
The main sequence of the MIST isochrone grid, against [Fe/H] (GEN.548).

Source: MESA Isochrones and Stellar Tracks (MIST) v1.2, the "basic" isochrones with rotation v/vcrit = 0.4, solar-scaled
([a/Fe] = 0), 15 values of [Fe/H] from -4.00 to +0.50 and 107 ages from 10^5 to 10^10.3 years (Choi et al. 2016, ApJ 823,
102, doi:10.3847/0004-637X/823/2/102; Dotter 2016, ApJS 222, 8, doi:10.3847/0067-0049/222/1/8; mist.science). The
archive is 220 MB and holds every phase; setup downloads it once, checks its sha256, keeps only the main-sequence rows
(phase 0) of each isochrone (initial mass, log L, log Teff, log R: about 340,000 rows, 4 MB) in `main_sequence.npz`
and removes the archive. The model choices (v/vcrit = 0.4, solar-scaled) are MIST's recommended grid and a default
of this build; the grid is the data, not a result tuned to any known star.

The readers interpolate: linearly in [Fe/H] between the two grid values around it (clamped at the ends of the grid),
linearly in log age between the two isochrone ages around the age, and linearly in initial mass along the isochrone.
They answer `None` when the star has left the main sequence at that age or is beyond the grid's ages, and raise
`fetch.StandIn` when the file has not been downloaded.
"""

import functools
import io
import re
import tarfile
from dataclasses import dataclass

import numpy as np

from planetgen.datasets import fetch

ARCHIVE = "MIST_v1.2_vvcrit0.4_basic_isos.txz"
DERIVED = "main_sequence.npz"

DATASET = fetch.Dataset(
    key="isochrones",
    title="MIST v1.2 isochrones, main sequence (v/vcrit 0.4, solar-scaled)",
    source="Choi et al. 2016 (ApJ 823, 102) and Dotter 2016 (ApJS 222, 8), MIST v1.2, mist.science",
    licence="open; cite Choi et al. 2016 and Dotter 2016",
    pieces=(fetch.Piece(ARCHIVE, "https://mist.science/data/tarballs_v1.2/" + ARCHIVE,
                        "bb3f42743f75676fa8bf6ca63c4cb2723a0839b77204efebcae39fa4610a92c8", 221507784, transient=True),),
    standin="the Baraffe et al. 2015 track at solar metallicity and the fitted relations of the stellar-composition note",
    fill="SC1",
    derived=(DERIVED,),
    derive=lambda folder: derive(folder),
    derive_version=1,
)

_COLUMNS = {"log_age": 1, "initial_mass": 2, "log_l": 7, "log_teff": 10, "log_r": 11, "phase": 24}
_ABUNDANCES = re.compile(rb"^#\s+([0-9.]+)\s+([0-9.Ee+-]+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)\s+([0-9.]+)\s*$")


def parse_isochrone_file(lines):
    """Reads one MIST `.iso` file (an iterable of byte lines): returns `(y_init, z_init, fe_h, rows)` where rows is a
    float32 array of the main-sequence rows with columns log age, initial mass, log L, log Teff, log R."""
    y_init = z_init = fe_h = None
    rows = []
    previous_was_header = False
    for line in lines:
        text = line.strip()
        if not text:
            continue
        if text.startswith(b"#"):
            if previous_was_header and y_init is None:
                match = _ABUNDANCES.match(text)
                if match:
                    y_init, z_init, fe_h = float(match.group(1)), float(match.group(2)), float(match.group(3))
            previous_was_header = b"Yinit" in text
            continue
        previous_was_header = False
        parts = text.split()
        if len(parts) <= _COLUMNS["phase"]:
            raise ValueError(f"not a MIST isochrone row: {len(parts)} columns, expected 25")
        if float(parts[_COLUMNS["phase"]]) != 0.0:
            continue
        rows.append([float(parts[_COLUMNS[name]]) for name in ("log_age", "initial_mass", "log_l", "log_teff", "log_r")])
    if y_init is None or not rows:
        raise ValueError("not a MIST isochrone file: no abundance line or no main-sequence rows")
    return y_init, z_init, fe_h, np.asarray(rows, dtype=np.float32)


def derive(folder):
    """The compact main-sequence file from the downloaded archive: `{DERIVED: npz bytes}`."""
    found = []
    with tarfile.open(f"{folder}/{ARCHIVE}") as archive:
        for member in archive:
            if member.isfile() and member.name.endswith(".iso"):
                found.append(parse_isochrone_file(archive.extractfile(member)))
    if not found:
        raise ValueError("the archive holds no .iso file")
    found.sort(key=lambda entry: entry[2])
    arrays = {"fe_h": np.array([entry[2] for entry in found]), "y_init": np.array([entry[0] for entry in found]),
              "z_init": np.array([entry[1] for entry in found])}
    for i, entry in enumerate(found):
        arrays[f"rows_{i}"] = entry[3]
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    return {DERIVED: buffer.getvalue()}


@dataclass(frozen=True)
class _Isochrone:
    """One [Fe/H]: its main-sequence blocks, one per age."""

    fe_h: float
    y_init: float
    z_init: float
    log_ages: np.ndarray
    blocks: tuple

    @staticmethod
    def _on_block(block, mass):
        masses = block[:, 0]
        if mass < masses[0] or mass > masses[-1]:
            return None
        return np.array([np.interp(mass, masses, block[:, column]) for column in (1, 2, 3)])

    def point(self, mass, log_age):
        """`(log L, log Teff, log R)` on this metallicity's isochrone at `log_age`, or None if off the main sequence."""
        ages = self.log_ages
        if log_age < ages[0] or log_age > ages[-1]:
            return None
        upper = int(np.searchsorted(ages, log_age))
        if ages[upper] == log_age:
            found = self._on_block(self.blocks[upper], mass)
            return None if found is None else tuple(found)
        younger, older = self._on_block(self.blocks[upper - 1], mass), self._on_block(self.blocks[upper], mass)
        if younger is None or older is None:
            return None
        weight = (log_age - ages[upper - 1]) / (ages[upper] - ages[upper - 1])
        return tuple(younger * (1 - weight) + older * weight)

    def lifetime_log_age(self, mass):
        """log10 of the age (years) at which a star of `mass` leaves the main sequence, or None outside the grid."""
        ends = np.array([block[-1, 0] for block in self.blocks])
        order = np.argsort(ends)
        if mass < ends[order][0] or mass > ends[order][-1]:
            return None
        return float(np.interp(mass, ends[order], self.log_ages[order]))


@functools.lru_cache(maxsize=4)
def _load(path):
    with np.load(path) as data:
        grid = []
        for i, fe_h in enumerate(data["fe_h"]):
            rows = data[f"rows_{i}"]
            ages = np.unique(rows[:, 0])
            blocks = tuple(rows[rows[:, 0] == age][:, 1:] for age in ages)
            grid.append(_Isochrone(float(fe_h), float(data["y_init"][i]), float(data["z_init"][i]), ages.astype(float), blocks))
    return tuple(grid)


def _grid():
    located = fetch.locate(DATASET).require()
    return _load(located.path(DERIVED))


def metallicities():
    """The [Fe/H] values of the grid, ascending."""
    return tuple(isochrone.fe_h for isochrone in _grid())


def _bracket(grid, fe_h):
    values = [isochrone.fe_h for isochrone in grid]
    fe_h = min(max(fe_h, values[0]), values[-1])
    upper = min(int(np.searchsorted(values, fe_h)), len(values) - 1)
    lower = upper if values[upper] == fe_h else upper - 1
    if lower == upper:
        return grid[upper], grid[upper], 0.0
    weight = (fe_h - values[lower]) / (values[upper] - values[lower])
    return grid[lower], grid[upper], weight


def main_sequence_point(mass_sol, age_gy, fe_h):
    """`(log10 L/Lsun, log10 Teff/K, log10 R/Rsun)` of a main-sequence star of initial mass `mass_sol` solar masses
    and `age_gy` Gy at `fe_h`, or None when it has left the main sequence, is not yet on it, or its age is outside
    the grid (10^5 to about 10^10.3 years)."""
    log_age = float(np.log10(age_gy * 1e9))
    low, high, weight = _bracket(_grid(), fe_h)
    first, second = low.point(mass_sol, log_age), high.point(mass_sol, log_age)
    if first is None or second is None:
        return None
    return tuple(float(a * (1 - weight) + b * weight) for a, b in zip(first, second))


def main_sequence_lifetime_gy(mass_sol, fe_h):
    """The main-sequence lifetime in Gy of a star of initial mass `mass_sol` at `fe_h`, or None when the grid does
    not reach it (below about 0.8 solar masses the lifetime is longer than the oldest isochrone)."""
    low, high, weight = _bracket(_grid(), fe_h)
    first, second = low.lifetime_log_age(mass_sol), high.lifetime_log_age(mass_sol)
    if first is None or second is None:
        return None
    return float(10 ** (first * (1 - weight) + second * weight) / 1e9)


def abundances(fe_h):
    """`(Y, Z)` the grid starts with at `fe_h`, interpolated between the two grid values around it."""
    low, high, weight = _bracket(_grid(), fe_h)
    return (low.y_init * (1 - weight) + high.y_init * weight, low.z_init * (1 - weight) + high.z_init * weight)
