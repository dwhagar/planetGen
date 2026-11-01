# planetgen/physics/uv_cross_sections.py

"""
Absorption cross sections on the wavelength grid of the radiation model (GEN.443; docs/design/radiation-model.md section 3).

The grid (`standard_edges`) has 1 nm cells from 100 to 400 nm, 10 nm cells from 400 to 1100 nm and 100 nm cells from 1100 to
4000 nm. A gas's cross section in a cell is read from the gas table's 35 coarse bins (`physics/gases.py`,
`Gas.uv`, part 02 section 2.3: [R1873, R2040] band means, 10 to 1100 nm), and is zero where the table has no data
(beyond 1100 nm the near infrared bands of H2O, CO2 and CH4 are not in the table: the near infrared band is Rayleigh and
aerosol only, a labelled stand-in until the correlated-k grids of R2 are built).

**Ozone** is the exception: a band mean over 20 nm is wrong where the Huggins band falls by three decades across the UV-B
(300 to 320 nm), and the UV index is made there. So `ozone_cross_section_cm2` reads the high-resolution cross sections
of the Brion, Daumont, Malicet and Jenouvrier group from the MPI-Mainz UV/VIS Spectral Atlas (`datasets/ozone_uv.py`,
218, 243 and 295 K between 195 and 345 nm at 0.01 nm, and 295 K between 195 and 830 nm at 1 nm; Zenodo record 6951,
CC BY 4.0, [R2302, R2303]) when the setup step has downloaded them, interpolated linearly in temperature between
the three, and the gas table's bins when it has not and the run disregards missing data.
"""

import functools

import numpy as np

from planetgen.datasets import fetch, ozone_uv
from planetgen.physics import gases

EDGES_NM = np.concatenate([np.arange(100.0, 400.0, 1.0), np.arange(400.0, 1100.0, 10.0), np.arange(1100.0, 4001.0, 100.0)])
"""ndarray: The cell edges, nm."""

DOBSON_CM2 = 2.687e16
"""float: Molecules per square centimetre in one Dobson unit."""

DEFAULT_OZONE_TEMPERATURE_K = 220.0
"""float: The temperature of the ozone layer of a body that has none worked out, K (Earth's stratosphere at the layer, 210 to 230 K)."""


def standard_edges():
    """The cell edges of the radiation model's wavelength grid, nm (`EDGES_NM`)."""
    return EDGES_NM


def centres(edges):
    """The cell centres of a grid given by its `edges`, nm."""
    edges = np.asarray(edges, dtype=float)
    return 0.5 * (edges[:-1] + edges[1:])


@functools.lru_cache(maxsize=None)
def _bin_values(formula):
    entry = gases.gas(formula).uv
    if entry is None:
        return None
    return np.array([0.0 if value is None else value for value in entry.cross_section_cm2])


_BIN_EDGES = None


def _bin_edges():
    global _BIN_EDGES
    if _BIN_EDGES is None:
        _BIN_EDGES = np.array([low for low, _high in gases.UV_BINS_NM] + [gases.UV_BINS_NM[-1][1]])
    return _BIN_EDGES


def band_mean_cross_section_cm2(formula, edges):
    """
    The cross section of `formula` in each cell of `edges` from the gas table's bins (zero where the table has no data), cm2 per
    molecule. A gas the table has no ultraviolet entry for (He, Ne, N2 in the ultraviolet bins read as zero by the research) is zero.
    """
    values = _bin_values(formula)
    middle = centres(edges)
    if values is None:
        return np.zeros_like(middle)
    bins = _bin_edges()
    index = np.searchsorted(bins, middle, side="right") - 1
    inside = (index >= 0) & (index < len(values))
    return np.where(inside, values[np.clip(index, 0, len(values) - 1)], 0.0)


def cross_section_cm2(formula, edges):
    """The cross section of `formula` (not ozone) on the cells of `edges`, cm2 per molecule (`band_mean_cross_section_cm2`)."""
    return band_mean_cross_section_cm2(formula, edges)


def ozone_cross_section_cm2(edges, temperature_k=None):
    """
    The ozone cross section on the cells of `edges` at `temperature_k` (K; `DEFAULT_OZONE_TEMPERATURE_K` when `None`), cm2 per molecule.

    The downloaded high-resolution tables (`datasets/ozone_uv.py`) are used: cell means, linear in temperature between 218, 243
    and 295 K, held outside that range, the gas table's bins under 195 nm. A run that disregards missing data
    takes the gas table's bins everywhere when the download is absent (`fetch.stand_ins_used` lists it); any other run raises
    `fetch.DataMissing`.
    """
    bins = band_mean_cross_section_cm2("O3", edges)
    located = fetch.locate(ozone_uv.DATASET)
    if located.use_stand_in("ozone cross sections from the 20 nm bins of the gas table"):
        return bins
    return ozone_uv.cross_section_cm2(located, edges, DEFAULT_OZONE_TEMPERATURE_K if temperature_k is None else temperature_k, bins)
