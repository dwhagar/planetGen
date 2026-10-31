# planetgen/datasets/isochrones.py

"""
The main sequence of the MIST isochrone grid, against [Fe/H] (GEN.548), and the tracks after it (GEN.935).

Source: MESA Isochrones and Stellar Tracks (MIST) v1.2, the "basic" isochrones with rotation v/vcrit = 0.4, solar-scaled
([a/Fe] = 0), 15 values of [Fe/H] from -4.00 to +0.50 and 107 ages from 10^5 to 10^10.3 years (Choi et al. 2016, ApJ 823,
102, doi:10.3847/0004-637X/823/2/102; Dotter 2016, ApJS 222, 8, doi:10.3847/0067-0049/222/1/8; mist.science). The
archive is 220 MB and holds every phase; setup downloads it once, checks its sha256, keeps the main-sequence rows
(phase 0) of each isochrone (initial mass, log L, log Teff, log R and the present mass: about 340,000 rows, 5 MB) in
`main_sequence.npz` and, in `post_main_sequence.npz`, the tracks from the end of the main sequence on that
`evolved_tracks` builds from the rest of the rows (GEN.935), and removes the archive. The model choices (v/vcrit = 0.4, solar-scaled) are MIST's recommended grid and a default
of this build; the grid is the data, not a result tuned to any known star.

The readers interpolate: linearly in [Fe/H] between the two grid values around it (clamped at the ends of the grid),
linearly in log age between the two isochrone ages around the age, and linearly in initial mass along the isochrone.
They answer `None` when the star has left the main sequence at that age or is beyond the grid's ages, and raise
`fetch.DataMissing` when the file has not been downloaded.
"""

import bisect
import functools
import io
import itertools
import math
import os
import re
import tarfile
import time
from dataclasses import dataclass

import numpy as np

from planetgen.datasets import fetch

ARCHIVE = "MIST_v1.2_vvcrit0.4_basic_isos.txz"
DERIVED = "main_sequence.npz"
POST_DERIVED = "post_main_sequence.npz"

DATASET = fetch.Dataset(
    key="isochrones",
    title="MIST v1.2 isochrones, main sequence (v/vcrit 0.4, solar-scaled)",
    source="Choi et al. 2016 (ApJ 823, 102) and Dotter 2016 (ApJS 222, 8), MIST v1.2, mist.science",
    pieces=(fetch.Piece(ARCHIVE, "https://mist.science/data/tarballs_v1.2/" + ARCHIVE,
                        "bb3f42743f75676fa8bf6ca63c4cb2723a0839b77204efebcae39fa4610a92c8", 221507784, transient=True),),
    replaces="the Baraffe et al. 2015 track at solar metallicity and the fitted relations of the stellar-composition note",
    fill="SC1",
    fallback="no MIST point, so no luminosity, temperature, radius or lifetime for a star (an open item on the star)",
    derived=(DERIVED, POST_DERIVED),
    derive=lambda folder: derive(folder),
    derive_version=3,
)

_COLUMNS = {"log_age": 1, "initial_mass": 2, "star_mass": 3, "log_l": 7, "log_teff": 10, "log_r": 11, "h1": 13, "c12": 16,
            "phase": 24}
_ABUNDANCES = re.compile(rb"^#\s+([0-9.]+)\s+([0-9.Ee+-]+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)\s+([0-9.]+)\s*$")
FIRST_POST_EEP = 454
"""int: The first equivalent evolutionary point after the main sequence (MIST's terminal-age main sequence): a row from it on
belongs to `post_main_sequence.npz`."""


def parse_isochrone_file(lines):
    """Reads one MIST `.iso` file (an iterable of byte lines): returns `(y_init, z_init, fe_h, rows, post_rows, pre_edge)` where
    `rows` is a float32 array of the main-sequence rows with columns log age, initial mass, log L, log Teff, log R and
    present mass, and `post_rows` a float64 array of the rows from `FIRST_POST_EEP` on with the columns of
    `evolved_tracks.ROW_COLUMNS` (equivalent evolutionary point, log age, initial mass, present mass, log L, log Teff,
    log R, surface hydrogen and carbon fractions, phase), and `pre_edge` the heaviest star still contracting at each age (GEN.947: the
    rows of phase -1), a float64 array of `(log age, initial mass)` from the youngest age to the oldest."""
    y_init = z_init = fe_h = None
    rows, post = [], []
    contracting = {}
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
        phase = float(parts[_COLUMNS["phase"]])
        if phase == -1.0:
            age, mass = float(parts[_COLUMNS["log_age"]]), float(parts[_COLUMNS["initial_mass"]])
            contracting[age] = max(contracting.get(age, 0.0), mass)
        elif phase == 0.0:
            rows.append([float(parts[_COLUMNS[name]])
                         for name in ("log_age", "initial_mass", "log_l", "log_teff", "log_r", "star_mass")])
        elif float(parts[0]) >= FIRST_POST_EEP:
            post.append([float(parts[0])] + [float(parts[_COLUMNS[name]]) for name in (
                "log_age", "initial_mass", "star_mass", "log_l", "log_teff", "log_r", "h1", "c12", "phase")])
    if y_init is None or not rows:
        raise ValueError("not a MIST isochrone file: no abundance line or no main-sequence rows")
    pre_edge = np.asarray(sorted(contracting.items()), dtype=np.float64).reshape(-1, 2)
    return (y_init, z_init, fe_h, np.asarray(rows, dtype=np.float32), np.asarray(post, dtype=np.float64).reshape(-1, 10),
            pre_edge)


def derive(folder):
    """The compact files from the downloaded archive: `{DERIVED: npz bytes, POST_DERIVED: npz bytes}`."""
    from planetgen.datasets import evolved_tracks
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
        arrays[f"pre_edge_{i}"] = entry[5]
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    post = {"fe_h": arrays["fe_h"]}
    for i, entry in enumerate(found):
        post[f"tracks_{i}"], post[f"ms_age_gy_{i}"] = evolved_tracks.build(entry[4])
        post[f"ms_age_gy_{i}"] = post[f"ms_age_gy_{i}"] / 1.0e9
    post_buffer = io.BytesIO()
    np.savez_compressed(post_buffer, **post)
    return {DERIVED: buffer.getvalue(), POST_DERIVED: post_buffer.getvalue()}


@dataclass(frozen=True)
class _Isochrone:
    """One [Fe/H]: its main-sequence blocks, one per age."""

    fe_h: float
    y_init: float
    z_init: float
    log_ages: np.ndarray
    blocks: tuple
    pre_edge: np.ndarray = None
    """`(log age, heaviest initial mass still contracting)` rows, youngest first (GEN.947); empty when the file holds none."""

    @staticmethod
    def _on_block(block, mass):
        masses = block[:, 0]
        if mass < masses[0] or mass > masses[-1]:
            return None
        upper = int(masses.searchsorted(mass))
        # The block is float32 and the mass a float64: a mass a hair past the last (or before the first) row compares
        # equal to it in float32 and the search puts it just outside, so it takes that row.
        if upper >= len(masses):
            return block[-1, 1:5]
        if upper == 0 or masses[upper] == mass:
            return block[upper, 1:5]
        weight = (mass - masses[upper - 1]) / (masses[upper] - masses[upper - 1])
        return block[upper - 1, 1:5] * (1.0 - weight) + block[upper, 1:5] * weight

    def point(self, mass, log_age):
        """`(log L, log Teff, log R, present mass)` on this metallicity's isochrone at `log_age`, or None if off the main
        sequence."""
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

    def span(self, mass):
        """`(first, last)` log10 age (years) of the isochrones this mass has a main-sequence row on: from the age it
        reaches the main sequence to the age it leaves it (the oldest isochrone for a mass that does not leave in
        the grid); None when no isochrone holds it (lighter than the grid's lowest mass or heavier than its highest)."""
        low, high = self._mass_limits()
        held = np.nonzero((low <= mass) & (mass <= high))[0]
        if not len(held):
            return None
        return float(self.log_ages[held[0]]), float(self.log_ages[held[-1]])

    def clamped_point(self, mass, log_age):
        """
        `(log L, log Teff, log R, present mass)` of `mass` at `log_age`, held to the mass's own span (GEN.652): an age before the first
        isochrone it is on takes that row, and one past the last takes that row carried on in age at the slope of the
        two last rows, up to the age the mass leaves the main sequence between isochrones (`lifetime_log_age`), then
        held. The isochrones come 0.05 dex of age apart, so the last row alone would step with the mass; carried on,
        the end of the main sequence is as smooth in the mass as the rows are. None for a mass the grid does not hold.
        """
        span = self.span(mass)
        if span is None:
            return None
        first, last = span
        if log_age <= last:
            return self.point(mass, max(log_age, first))
        end = self.point(mass, last)
        earlier = self.point(mass, last - _END_SLOPE_STEP_DEX) if last - _END_SLOPE_STEP_DEX >= first else None
        leaves = self.lifetime_log_age(mass)
        if earlier is None or leaves is None or leaves <= last:
            return end
        reach = min(log_age, leaves) - last
        return tuple(float(e + (e - o) / _END_SLOPE_STEP_DEX * reach) for e, o in zip(end, earlier))

    def _mass_limits(self):
        cached = _LIMITS.get(id(self))
        if cached is None:
            cached = (np.array([block[0, 0] for block in self.blocks]), np.array([block[-1, 0] for block in self.blocks]))
            _LIMITS[id(self)] = cached
        return cached

    def lifetime_log_age(self, mass):
        """log10 of the age (years) at which a star of `mass` leaves the main sequence, or None outside the grid."""
        ends = np.array([block[-1, 0] for block in self.blocks])
        order = np.argsort(ends)
        if mass < ends[order][0] or mass > ends[order][-1]:
            return None
        return float(np.interp(mass, ends[order], self.log_ages[order]))


_END_SLOPE_STEP_DEX = 0.05
"""float: How far back in log age the second row is read for the slope the end of the main sequence is carried on at
(the isochrones' own spacing)."""

_LIMITS = {}
"""dict: Each isochrone's lightest and heaviest main-sequence mass at each age (`_Isochrone._mass_limits`), by `id`."""


@functools.lru_cache(maxsize=4)
def _load(path):
    with np.load(path) as data:
        grid = []
        for i, fe_h in enumerate(data["fe_h"]):
            rows = data[f"rows_{i}"]
            ages = np.unique(rows[:, 0])
            blocks = tuple(rows[rows[:, 0] == age][:, 1:] for age in ages)
            grid.append(_Isochrone(float(fe_h), float(data["y_init"][i]), float(data["z_init"][i]), ages.astype(float), blocks,
                                   data[f"pre_edge_{i}"]))
    return tuple(grid)


_SEEN = {}
"""dict: The grid last read (`_grid`): `{"stamp", "grid", "dataset", "env", "checked", "token"}`."""

RECHECK_SECONDS = 2.0
"""float: How long `_grid` trusts the file it read before it looks at the file's stamp again (a reader in a loop, GEN.652:
every main-sequence star's state, asks millions of times, and the files change only when a download replaces them)."""

_TOKENS = itertools.count(1)


def forget():
    """Drops the grid read so far; the next reader checks the dataset and reads it again."""
    _SEEN.clear()


def grid_token():
    """A number that changes whenever `_grid` reads a different grid, for a caller that remembers what it worked out
    from the grid (`physics.main_sequence`)."""
    seen = _SEEN
    if seen and seen["dataset"] is DATASET and time.monotonic() - seen["checked"] < RECHECK_SECONDS:
        return seen["token"]
    _grid()
    return _SEEN["token"]


def _grid():
    """The grid. The dataset is checked against the catalogue (`fetch.locate`: sizes and hashes) when the derived file is
    first read or has changed on disk; a reader in a loop only compares the dataset it read from and, every
    `RECHECK_SECONDS`, looks at the cache folder and stats the file (a test that moves the folder calls `forget`)."""
    seen = _SEEN
    now = time.monotonic()
    if seen and seen["dataset"] is DATASET:
        if now - seen["checked"] < RECHECK_SECONDS:
            return seen["grid"]
        try:
            stat = os.stat(seen["path"])
        except OSError:
            stat = None
        if (stat is not None and seen["env"] == os.environ.get(fetch.CACHE_DIR_ENV)
                and (stat.st_mtime_ns, stat.st_size) == seen["stamp"]):
            seen["checked"] = now
            return seen["grid"]
    located = fetch.locate(DATASET).require()
    path = located.path(DERIVED)
    grid = _load(path)
    stat = os.stat(path)
    _SEEN.clear()
    _EDGES.clear()
    _LIMITS.clear()
    _SEEN.update(grid=grid, dataset=DATASET, env=os.environ.get(fetch.CACHE_DIR_ENV), path=path, checked=now,
                 stamp=(stat.st_mtime_ns, stat.st_size), token=next(_TOKENS))
    return grid


def metallicities():
    """The [Fe/H] values of the grid, ascending."""
    return tuple(isochrone.fe_h for isochrone in _grid())


def _bracket(grid, fe_h):
    values = [isochrone.fe_h for isochrone in grid]
    fe_h = min(max(fe_h, values[0]), values[-1])
    upper = min(bisect.bisect_left(values, fe_h), len(values) - 1)
    lower = upper if values[upper] == fe_h else upper - 1
    if lower == upper:
        return grid[upper], grid[upper], 0.0
    weight = (fe_h - values[lower]) / (values[upper] - values[lower])
    return grid[lower], grid[upper], weight


EDGE_MASS_RANGE_SOL = (0.1, 250.0)
"""tuple: The masses (solar) `edge_point` tabulates its main-sequence arrival and departure rows over."""

EDGE_MASS_STEPS = 3000
"""int: How many log-spaced masses `edge_point`'s tables hold across `EDGE_MASS_RANGE_SOL` (0.3% of mass apart)."""

_EDGE_LOW_LOG = math.log10(EDGE_MASS_RANGE_SOL[0])
_EDGE_STEPS_PER_DEX = EDGE_MASS_STEPS / (math.log10(EDGE_MASS_RANGE_SOL[1]) - _EDGE_LOW_LOG)

_EDGES = {}
"""dict: Each isochrone's arrival and departure rows over the mass grid (`_edge_table`), by `id`."""


def _edge_table(isochrone):
    """`(arrival, departure)`, each a list over the mass grid of `(log L, log Teff, log R, present mass)` or None where the grid holds
    no main-sequence row of that mass, of one isochrone (`_Isochrone.clamped_point` at its first and last age)."""
    cached = _EDGES.get(id(isochrone))
    if cached is None:
        low, high = (np.log10(m) for m in EDGE_MASS_RANGE_SOL)
        arrival, departure = [], []
        for log_mass in np.linspace(low, high, EDGE_MASS_STEPS + 1):
            mass = float(10.0 ** log_mass)
            span = isochrone.span(mass)
            arrival.append(None if span is None else tuple(float(v) for v in isochrone.clamped_point(mass, span[0])))
            departure.append(None if span is None else tuple(float(v) for v in isochrone.clamped_point(mass, 99.0)))
        cached = _EDGES[id(isochrone)] = (arrival, departure)
    return cached


def edge_point(mass_sol, fe_h, departure):
    """
    `(log10 L/Lsun, log10 Teff/K, log10 R/Rsun, present mass in Msun)` of a star of initial mass `mass_sol` and `fe_h` as it arrives on the
    main sequence (`departure` False) or leaves it (True): `main_sequence_point(..., clamp=True)` at its first and last
    age, read from a table of 3000 masses, interpolated in log mass and in [Fe/H] (GEN.652: a table build and every
    star's record read these a million times; the rows themselves cost a few hundred microseconds each). None for a mass
    off the table or between two the grid holds no row for.
    """
    low, high = EDGE_MASS_RANGE_SOL
    if not low <= mass_sol <= high:
        return None
    position = (math.log10(mass_sol) - _EDGE_LOW_LOG) * _EDGE_STEPS_PER_DEX
    index = min(int(position), EDGE_MASS_STEPS - 1)
    along = position - index
    below, above, weight = _bracket(_grid(), fe_h)
    found = []
    for isochrone in (below, above):
        rows = _edge_table(isochrone)[1 if departure else 0]
        near, far = rows[index], rows[index + 1]
        if near is None or far is None:
            return None
        found.append(tuple(n + (f - n) * along for n, f in zip(near, far)))
    first, second = found
    return tuple(a + (b - a) * weight for a, b in zip(first, second))


def _arrival_log_age(isochrone, log_mass):
    """The log age at which `isochrone`'s heaviest contracting star has fallen to the mass (log10), interpolated, or None outside
    the masses the rows hold."""
    edge = isochrone.pre_edge
    if edge is None or len(edge) < 2:
        return None
    ages, log_masses = edge[:, 0], np.log10(edge[:, 1])
    if log_mass > log_masses.max() or log_mass < log_masses.min():
        return None
    # The heaviest contracting mass falls with age; np.interp wants rising x, so the masses are read from the oldest age down.
    order = np.argsort(log_masses, kind="stable")
    return float(np.interp(log_mass, log_masses[order], ages[order]))


def pre_main_sequence_arrival_yr(mass_sol, fe_h):
    """
    The age in years at which a star of initial mass `mass_sol` and birth `fe_h` has finished contracting onto the main sequence
    in MIST (GEN.947): the age at which the heaviest star still in phase -1 of the isochrones (equivalent evolutionary points
    before 202) has fallen to that mass, interpolated between the isochrone ages (0.05 dex apart) and between the two nearest
    [Fe/H]. None when the mass is outside the rows (heavier than the contracting stars of the youngest isochrone, or lighter than
    those of the oldest). Raises `fetch.DataMissing` when the isochrones were not downloaded.
    """
    log_mass = math.log10(mass_sol)
    below, above, weight = _bracket(_grid(), fe_h)
    found = [_arrival_log_age(isochrone, log_mass) for isochrone in (below, above)]
    if found[0] is None or found[1] is None:
        return None
    return 10.0 ** (found[0] + (found[1] - found[0]) * weight)


def lowest_arrival_mass_sol(fe_h):
    """The lightest initial mass, solar masses, for which `pre_main_sequence_arrival_yr` has an age at `fe_h`: the lightest star still contracting in the
    oldest isochrone that has one, the larger of the two grid metallicities' (GEN.1173), or None when the rows have no contracting star."""
    below, above, _ = _bracket(_grid(), fe_h)
    lowest = []
    for isochrone in (below, above):
        edge = isochrone.pre_edge
        if edge is None or len(edge) < 2:
            return None
        lowest.append(float(np.log10(edge[:, 1]).min()))
    return float(10.0 ** max(lowest)) * (1.0 + 1e-9)


def main_sequence_point(mass_sol, age_gy, fe_h, clamp=False):
    """`(log10 L/Lsun, log10 Teff/K, log10 R/Rsun, present mass in Msun)` of a main-sequence star of initial mass `mass_sol` solar masses
    and `age_gy` Gy at `fe_h`, or None when it has left the main sequence, is not yet on it, or its age is outside
    the grid (10^5 to about 10^10.3 years). With `clamp` the age is held to the mass's own span on each isochrone, the
    first row it has on the main sequence to the last (GEN.652): a star that the generator's phase times keep on the
    main sequence a little past or before MIST's own edge of it takes the row at that edge; None then only for a
    mass the grid does not hold."""
    log_age = float(np.log10(max(age_gy, 1e-300) * 1e9))
    low, high, weight = _bracket(_grid(), fe_h)
    if clamp:
        first, second = (isochrone.clamped_point(mass_sol, log_age) for isochrone in (low, high))
    else:
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


_POST = {}
"""dict: The tracks last read (`post_main_sequence_grid`): `{"token", "grid"}`."""


def post_main_sequence_grid():
    """The `evolved_tracks.Grid` of the post-main-sequence tracks (GEN.935), read from the derived file once for the grid
    `_grid` has read: it follows the main sequence's checks of the dataset and the cache folder, and is read again when
    that grid is. Raises `fetch.DataMissing` when the isochrones were not downloaded."""
    from planetgen.datasets import evolved_tracks
    _grid()
    token = _SEEN["token"]
    if _POST.get("token") != token:
        path = fetch.locate(DATASET).require().path(POST_DERIVED)
        with np.load(path) as data:
            fe_h = data["fe_h"]
            grid = evolved_tracks.Grid(fe_h, [data[f"tracks_{i}"] for i in range(len(fe_h))],
                                       [data[f"ms_age_gy_{i}"] for i in range(len(fe_h))])
        _POST.clear()
        _POST.update(token=token, grid=grid)
    return _POST["grid"]


def post_main_sequence_track(mass_sol, fe_h):
    """The `evolved_tracks.Track` of a star of initial mass `mass_sol` and birth `fe_h` from the end of its main sequence on
    (the four nearest tracks blended, GEN.935), or None when no track of the grid reaches the thermal pulses."""
    return post_main_sequence_grid().track(mass_sol, fe_h)



def post_main_sequence_light(mass_sol, fe_h):
    """`(frac, log_l)` of the track of `post_main_sequence_track`, for the readers that need only those two columns."""
    return post_main_sequence_grid().light(mass_sol, fe_h)


def helium_star_relation():
    """`(intercept, slope, sigma)` of the luminosity of a Wolf-Rayet star against its present mass in MIST
    (`evolved_tracks.Grid.helium_star_relation`, GEN.939); None when the tracks hold none. Raises `fetch.DataMissing` when the
    isochrones were not downloaded."""
    return post_main_sequence_grid().helium_star_relation()
