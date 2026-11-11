# planetgen/datasets/methane_visible.py

"""
Methane absorption in the visible and near infrared at 1 nm resolution, 400 to 1050 nm (GEN.949, AD12).

Source: Karkoschka (1998, Icarus 133, 134, [R2259]), the volume GBAT_0001 of the NASA Planetary Data System Atmospheres Node ("Spectrophotometry of the Jovian
Planets and Titan at 300- to 1050-nm Wavelength: the Methane Spectrum"; the data are those of Karkoschka 1994 and 1998, 1 nm resolution at 0.4 nm
sampling from ESO spectrophotometry of 6 to 10 July 1995). One file is read, `data/1995low.tab`: the methane absorption coefficient in 1/(km-amagat) at each
wavelength, one effective spectrum for the cold giants (the paper puts its temperature near 100 K). The volume states no licence or use restriction.

The gas table holds methane as 35 coarse bins whose band means are right on average but flatten the structure of the bands: the strong ones (727, 790, 840,
890 nm) are a hundred times deeper than their cell mean and nearly black, the windows between them nearly clear, and the colour of a giant is the light
that gets through the windows. `physics/uv_cross_sections.methane_optical_depth` reads this table: the depth of a cell is the optical depth of the mean
transmission of its samples, `-ln(mean(exp(-N sigma)))` for a vertical column of `N` molecules per cm2, not the mean of the depth.

The file is small and kept as it is (no derived file). The gas table's bins are the labelled fallback for a run that disregards missing data.
"""

import functools

import numpy as np

from planetgen.datasets import fetch

TABLE = "1995low.tab"
SOURCE = "Karkoschka (1998), NASA PDS volume GBAT_0001, data/1995low.tab"

LOSCHMIDT_CM3 = 2.6867773e19
"""float: The molecules in one cm3 at one amagat (STP: 273.15 K, 1 atm), the unit of the coefficients' column (1/(km-amagat))."""

CM_PER_KM = 1.0e5
"""float: Centimeters in a kilometer."""

RANGE_NM = (400.0, 1050.0)
"""tuple: The cells the table is used for, nm: under 400 nm the file holds no methane absorption (the gas table's bins stand), over 1050 nm it ends."""

DATASET = fetch.Dataset(
    key="methane_visible",
    title="Methane absorption coefficients, 300 to 1050 nm at 1 nm resolution (Karkoschka 1998, NASA PDS GBAT_0001)",
    source="[R2259] Karkoschka (1998), Icarus 133, 134; NASA PDS Atmospheres Node volume GBAT_0001 (no licence stated)",
    pieces=(fetch.Piece(TABLE, "https://atmos.nmsu.edu/PDS/data/gbat_0001/data/1995low.tab",
                        "d8071e071c2c0a60052ec6c164e151bd42847acf5d9bfbe1e9b2f76bf49a788b", 101250),),
    replaces="the 35 coarse band means of the gas table for methane in the visible, which flatten the strong bands and the windows of a giant's colour",
    fill="AD12",
    fallback="the 35 band means of the gas table (physics/data/gases.json) for methane",
)


@functools.lru_cache(maxsize=2)
def _load(path):
    """`(air wavelength nm, cross section cm2 per molecule)` of every sample of the table, sorted by wavelength."""
    table = np.loadtxt(path)
    if table.ndim != 2 or table.shape[1] < 3 or np.any(np.diff(table[:, 1]) <= 0) or np.any(table[:, 2] < 0):
        raise ValueError(f"{path}: expected rising air wavelengths in column 2 and non-negative coefficients in column 3")
    return table[:, 1].copy(), table[:, 2] / (CM_PER_KM * LOSCHMIDT_CM3)


def samples(located):
    """The table as `(wavelength nm, cross section cm2 per molecule)` arrays from a `fetch.Located`."""
    return _load(located.path(TABLE))


def covered(edges):
    """Which cells of `edges` the table serves: the whole cell inside `RANGE_NM`."""
    edges = np.asarray(edges, dtype=float)
    return (edges[:-1] >= RANGE_NM[0]) & (edges[1:] <= RANGE_NM[1])


def cell_optical_depth(located, edges, number_column_cm2):
    """
    The methane optical depth of each cell of `edges` for a vertical column of `number_column_cm2` molecules per cm2: `-ln` of the mean transmission
    of the samples inside the cell (a table 1 nm wide steps over cells of 10 nm). NaN for a cell the table does not serve (`covered`) or that holds no sample.
    """
    edges = np.asarray(edges, dtype=float)
    if not number_column_cm2 >= 0.0:
        raise ValueError("the column must be zero or more")
    nm, sigma = samples(located)
    depth = np.full(len(edges) - 1, np.nan)
    first = np.searchsorted(nm, edges[:-1], side="left")
    last = np.searchsorted(nm, edges[1:], side="left")
    for cell in np.flatnonzero(covered(edges)):
        if last[cell] <= first[cell]:
            continue
        tau = number_column_cm2 * sigma[first[cell]:last[cell]]
        floor = tau.min()
        depth[cell] = floor - np.log(np.mean(np.exp(floor - tau)))
    return depth
