# planetgen/datasets/baraffe_isochrones.py

"""
The [R0254] (BHAC15) isochrones for low-mass stars and brown dwarfs (GEN.936, GEN.948): at each of 30 ages from
0.5 Myr to 10 Gyr, the effective temperature, luminosity, surface gravity and radius of stars from 0.01 to 1.4 solar masses at
solar metallicity [R0254].

The file is `BHAC15_iso.2mass`, from Isabelle Baraffe's page at ENS Lyon
(https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_iso.2mass). It is one block per age, headed `! t (Gyr) = ...`,
with a row of mass, effective temperature, log luminosity, log g, radius, lithium ratio and three 2MASS absolute magnitudes. The
reader keeps the first five columns. The lightest masses are missing from the older blocks (the 10 Gyr block starts at 0.07
solar masses); `interpolate` answers None where the grid does not reach.

The pre-main-sequence times are three constants fitted to this file by hand. The brown dwarf's temperature and radius at an age
(`physics/brown_dwarf.py`, GEN.948) and the main sequence under the MIST grid's 0.1 solar masses (`physics/main_sequence.py`,
GEN.936) are read from it.

The reader raises `fetch.DataMissing` when the file has not been downloaded.
"""

import functools
import math
import re
from dataclasses import dataclass

import numpy as np

from planetgen import tuning
from planetgen.datasets import fetch

FILE = "BHAC15_iso.2mass"

DATASET = fetch.Dataset(
    key="baraffe_isochrones",
    title="[R0254] (BHAC15) isochrones of low-mass stars and brown dwarfs",
    source="[R0254]",
    pieces=(fetch.Piece(FILE, "https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_iso.2mass",
                        "e86198b15e0f28437f7ee83bcebe93cbbf2b02ee2ab6097295eee3155a28bcf2", 64951),),
    replaces="the nine brown dwarf radius points and the power-law extension of the main sequence under 0.1 solar masses (the three "
             "pre-main-sequence constants of tuning.py stay as fitted to it)",
    fill="SC6",
    fallback="the fitted pre-main-sequence constants of tuning.py; a brown dwarf's radius and the main sequence under 0.1 solar masses are left out",
)

_AGE = re.compile(r"!\s+t \(Gyr\)\s*=\s*([0-9.]+)")


@dataclass(frozen=True)
class Isochrone:
    """One age of the grid, ascending in mass."""

    age_gyr: float
    mass_solar: np.ndarray
    teff_k: np.ndarray
    log_luminosity_solar: np.ndarray
    log_g: np.ndarray
    radius_solar: np.ndarray


def parse(text):
    """The isochrones of a BHAC15 file, youngest first: a block per `! t (Gyr) =` line, rows of nine numbers."""
    blocks, rows = [], None
    for line in text.splitlines():
        age = _AGE.match(line)
        if age:
            rows = []
            blocks.append((float(age.group(1)), rows))
        elif rows is not None and line.strip() and not line.startswith(("!", "-")):
            numbers = line.split()
            if len(numbers) != 9:
                raise ValueError(f"a BHAC15 row has nine columns, not {len(numbers)}: {line!r}")
            rows.append([float(value) for value in numbers[:5]])
    if not blocks or not all(rows for _, rows in blocks):
        raise ValueError("a BHAC15 file with no age blocks, or an age block with no rows")
    isochrones = []
    for age, rows in blocks:
        table = np.asarray(rows)
        if np.any(np.diff(table[:, 0]) <= 0):
            raise ValueError(f"the masses of the {age} Gyr block do not rise")
        isochrones.append(Isochrone(age, table[:, 0], table[:, 1], table[:, 2], table[:, 3], table[:, 4]))
    if any(a.age_gyr >= b.age_gyr for a, b in zip(isochrones, isochrones[1:])):
        raise ValueError("the ages of a BHAC15 file do not rise")
    return tuple(isochrones)


@functools.lru_cache(maxsize=2)
def _load(path):
    with open(path, encoding="utf-8") as handle:
        return parse(handle.read())


def isochrones():
    """The 30 isochrones, youngest first."""
    return _load(fetch.locate(DATASET).require().path(FILE))


def interpolate(mass_solar, age_gyr, grid=None):
    """
    The state of a star of `mass_solar` at `age_gyr`: linear in log mass at each of the two bracketing ages, then linear in log
    age between them.

    Returns:
        dict: `teff_k`, `log_luminosity_solar`, `log_g` and `radius_solar`; None when the grid does not reach the mass at either
            bracketing age, or the age is outside 0.5 Myr to 10 Gyr.
    """
    grid = isochrones() if grid is None else grid
    ages = [iso.age_gyr for iso in grid]
    if not ages[0] <= age_gyr <= ages[-1] or mass_solar <= 0:
        return None
    upper = min(int(np.searchsorted(ages, age_gyr)), len(grid) - 1)
    lower = max(upper - 1, 0) if ages[upper] > age_gyr else upper
    weight = 0.0 if lower == upper else math.log(age_gyr / ages[lower]) / math.log(ages[upper] / ages[lower])
    states = []
    for iso in (grid[lower], grid[upper]):
        if not iso.mass_solar[0] <= mass_solar <= iso.mass_solar[-1]:
            return None
        log_mass = np.log(iso.mass_solar)
        states.append([float(np.interp(math.log(mass_solar), log_mass, column))
                       for column in (iso.teff_k, iso.log_luminosity_solar, iso.log_g, iso.radius_solar)])
    first, second = states
    mixed = [a + weight * (b - a) for a, b in zip(first, second)]
    return dict(zip(("teff_k", "log_luminosity_solar", "log_g", "radius_solar"), mixed))


def cooling_track(mass_sol, age_gyr):
    """
    `(teff_k, radius_solar)` of a body of `mass_sol` solar masses at `age_gyr`, for the objects that only cool (GEN.948): `interpolate`
    where the table holds the mass, the first age's values before it. The table loses a mass at the age it has cooled below the
    models' range (the lightest after 0.01 Gyr, 0.05 solar masses after 0.625 Gyr): past that age a mass the table still holds
    at `tuning.BROWN_DWARF_TRACK_END_GY` holds its last values (the objects near the hydrogen-burning limit settle on a plateau),
    and one that ended sooner cools on at [R0606]'s `tuning.BROWN_DWARF_TEMPERATURE_AGE_EXPONENT` with its last radius.

    Raises:
        ValueError: The mass is outside the table (0.01 to 1.4 solar masses).
    """
    grid = isochrones()
    if not grid[0].mass_solar[0] <= mass_sol <= grid[0].mass_solar[-1]:
        raise ValueError(f"no cooling track for {mass_sol} solar masses (the models run {grid[0].mass_solar[0]} to {grid[0].mass_solar[-1]})")
    last_age = grid[max(index for index, iso in enumerate(grid) if iso.mass_solar[0] <= mass_sol)].age_gyr
    state = interpolate(mass_sol, min(max(age_gyr, grid[0].age_gyr), last_age), grid)
    temperature = state["teff_k"]
    if age_gyr > last_age and last_age < tuning.BROWN_DWARF_TRACK_END_GY:
        temperature *= (age_gyr / last_age) ** tuning.BROWN_DWARF_TEMPERATURE_AGE_EXPONENT
    return temperature, state["radius_solar"]
