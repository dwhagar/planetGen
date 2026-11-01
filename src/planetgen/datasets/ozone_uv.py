# planetgen/datasets/ozone_uv.py

"""
Ozone absorption cross sections at high resolution, 195 to 830 nm, at three temperatures (GEN.443, AD9).

Source: the MPI-Mainz UV/VIS Spectral Atlas of Gaseous Molecules, the frozen snapshot of 2013-06-24 ([R2303], Earth Syst. Sci. Data 5, 365, [R2303]; Zenodo record 6951,
[R2302], CC BY 4.0; 126 MB). Nine of its 5,756 files are read, four of the Brion, Daumont, Malicet and
Jenouvrier group [R0543, R1000, R2832] and five of [R0599]:

- 218 K, 243 K and 295 K from 195 to 345 nm at 0.01 nm (the Hartley and Huggins bands; the temperature changes the UV-B
  cross section by a factor of three at 330 nm), and
- 295 K from 195 to 830 nm at 1 nm (the Chappuis band, 410 to 760 nm, peak 5.1e-21 cm2 near 602 nm), and
- 220, 240, 260, 280 and 298 K from 407 to 762 nm ([R0599]; GEN.748): the Chappuis band changes with the
  temperature by up to a third at its blue edge, and these five tables give the change as a ratio to the 298 K one.

The archive is a transient piece: setup downloads it once, checks its sha256, writes `ozone_uv.npz` with the four tables and
removes the archive. The gas table's 35 coarse bins (`physics/gases.py`) are the labelled fallback: a 20 nm band mean is
wrong where the Huggins band falls by three decades across the UV-B (300 to 320 nm), so the UV index it gives is off by a
large factor, and a run that disregards missing data says so.

The reader returns cell means on the radiation model's wavelength grid (`physics/uv_cross_sections.py`).
"""

import functools
import io
import zipfile

import numpy as np

from planetgen.datasets import fetch

ARCHIVE = "uv-vis-spectral-atlas-2013-06-24.zip"
DERIVED = "ozone_uv.npz"
SOURCE = "Brion-team 0.01 nm (MPI-Mainz UV/VIS Spectral Atlas)"

MEMBERS = {
    "218": "cross_sections/Ozone/O3_Brion-team(1992-1995)_218K_195.00-345.00nm.txt",
    "243": "cross_sections/Ozone/O3_Brion-team(1992-1995)_243K_195.00-345.00nm.txt",
    "295": "cross_sections/Ozone/O3_Brion-team(1992-1995)_295K_195.00-345.00nm.txt",
    "vis": "cross_sections/Ozone/O3_Brion-team(1992-1998)_295K_195-830nm(1nm).txt",
    "bt220": "cross_sections/Ozone/O3_BurkholderTalukdar(1994)_220K_407-762nm.txt",
    "bt240": "cross_sections/Ozone/O3_BurkholderTalukdar(1994)_240K_407-762nm.txt",
    "bt260": "cross_sections/Ozone/O3_BurkholderTalukdar(1994)_260K_407-762nm.txt",
    "bt280": "cross_sections/Ozone/O3_BurkholderTalukdar(1994)_280K_407-762nm.txt",
    "bt298": "cross_sections/Ozone/O3_BurkholderTalukdar(1994)_298K_407-762nm.txt",
}
"""dict: The archive members read, by the key they are kept under."""

TEMPERATURES_K = (218.0, 243.0, 295.0)
"""tuple: The temperatures of the three ultraviolet tables."""
CHAPPUIS_TEMPERATURES_K = (220.0, 240.0, 260.0, 280.0, 298.0)
"""tuple: The temperatures of the five Chappuis tables of [R0599], keys `bt220` to `bt298`."""
CHAPPUIS_RANGE_NM = (407.0, 762.0)
"""tuple: What those tables cover; the ratio to the 298 K table is 1 outside it."""

ULTRAVIOLET_RANGE_NM = (195.0, 345.0)
"""tuple: What the three temperature tables cover."""
VISIBLE_RANGE_NM = (345.0, 830.0)
"""tuple: What the 295 K table adds above them."""


def derive(folder):
    """Writes `ozone_uv.npz` from the verified archive in `folder`: `nm_<key>` and `sigma_<key>` for each of `MEMBERS`."""
    arrays = {}
    with zipfile.ZipFile(f"{folder}/{ARCHIVE}") as archive:
        for key, member in MEMBERS.items():
            rows = [line.split() for line in archive.read(member).decode("latin-1").splitlines() if line.strip()]
            table = np.array([[float(a), float(b)] for a, b in rows])
            if table.shape[1] != 2 or np.any(np.diff(table[:, 0]) <= 0) or np.any(table[:, 1] < 0):
                raise ValueError(f"{member}: expected rising wavelengths and non-negative cross sections")
            arrays[f"nm_{key}"], arrays[f"sigma_{key}"] = table[:, 0], table[:, 1]
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    return {DERIVED: buffer.getvalue()}


DATASET = fetch.Dataset(
    key="ozone_uv",
    title="Ozone absorption cross sections, 195 to 830 nm at 218, 243 and 295 K, and the Chappuis band at 220 to 298 K (Brion, Daumont, Malicet; Burkholder and Talukdar; MPI-Mainz atlas)",
    source="[R2303] (ESSD 5, 365), MPI-Mainz UV/VIS Spectral Atlas snapshot 2013-06-24, Zenodo [R2302], CC BY 4.0",
    pieces=(fetch.Piece(ARCHIVE, "https://zenodo.org/records/6951/files/" + ARCHIVE + "?download=1",
                        "bea06f6f390b774b5f0aeffda1b912a6446ba287578d66486414c728c148e6c6", 126277969, transient=True),),
    replaces="the 20 nm band means of the gas table for ozone, which understate the transmission in the UV-B and misplace the UV index",
    fill="AD9",
    fallback="the 35 band means of the gas table (physics/data/gases.json) for ozone",
    derived=(DERIVED,),
    derive=lambda folder: derive(folder),
    derive_version=2,
)


@functools.lru_cache(maxsize=2)
def _load(path):
    with np.load(path) as stored:
        return {name: np.array(stored[name]) for name in stored.files}


def tables(located):
    """The four tables as a dict of arrays (`nm_218`, `sigma_218`, ..., `nm_vis`, `sigma_vis`) from a `fetch.Located`."""
    return _load(located.path(DERIVED))


def _cell_means(nm, sigma, edges):
    """The mean of the piecewise-linear table over each cell of `edges`; NaN for a cell the table does not cover."""
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (sigma[1:] + sigma[:-1]) * np.diff(nm))])
    low, high = edges[:-1], edges[1:]
    inside = (low >= nm[0]) & (high <= nm[-1])
    area = np.interp(high, nm, cumulative) - np.interp(low, nm, cumulative)
    return np.where(inside, area / (high - low), np.nan)


def chappuis_ratio(located, edges, temperature_k):
    """
    The cross section of the Chappuis band at `temperature_k` over the 298 K one, on the cells of `edges` (GEN.748): the ratio of the cell
    means of the tables of [R0599], linear in temperature between 220, 240, 260, 280 and 298 K, held outside that range,
    and 1 where the tables do not reach (under 407 nm, over 762 nm) or hold no ozone absorption.
    """
    edges = np.asarray(edges, dtype=float)
    data = tables(located)
    temperatures = CHAPPUIS_TEMPERATURES_K
    t = min(max(temperature_k, temperatures[0]), temperatures[-1])
    upper = next(i for i, value in enumerate(temperatures) if value >= t)
    lower = max(upper - 1, 0)
    weight = 0.0 if upper == lower else (t - temperatures[lower]) / (temperatures[upper] - temperatures[lower])
    reference = _cell_means(data["nm_bt298"], data["sigma_bt298"], edges)
    mix = np.zeros_like(reference)
    for index, share in ((lower, 1.0 - weight), (upper, weight)):
        key = "bt%d" % temperatures[index]
        mix = mix + share * _cell_means(data["nm_" + key], data["sigma_" + key], edges)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(np.isfinite(reference) & (reference > 0.0), mix / reference, 1.0)
    return np.where(np.isfinite(ratio), ratio, 1.0)


def cross_section_cm2(located, edges, temperature_k, fallback):
    """
    The ozone cross section on the cells of `edges` at `temperature_k`, cm2 per molecule.

    Args:
        located (fetch.Located): The downloaded dataset.
        edges (ndarray): The cell edges, nm.
        temperature_k (float): The temperature; linear between 218, 243 and 295 K, held outside (the Chappuis band, 220 to 298 K).
        fallback (ndarray): The gas-table value of each cell, used where the tables do not reach (under 195 nm, over 830 nm).
    """
    edges = np.asarray(edges, dtype=float)
    data = tables(located)
    result = np.array(fallback, dtype=float)
    middle = 0.5 * (edges[:-1] + edges[1:])
    t_low, t_mid, t_high = TEMPERATURES_K
    t = min(max(temperature_k, t_low), t_high)
    if t <= t_mid:
        first, second, weight = "218", "243", (t - t_low) / (t_mid - t_low)
    else:
        first, second, weight = "243", "295", (t - t_mid) / (t_high - t_mid)
    ultraviolet = (1.0 - weight) * _cell_means(data[f"nm_{first}"], data[f"sigma_{first}"], edges) \
        + weight * _cell_means(data[f"nm_{second}"], data[f"sigma_{second}"], edges)
    visible = _cell_means(data["nm_vis"], data["sigma_vis"], edges) * chappuis_ratio(located, edges, temperature_k)
    in_ultraviolet = (middle >= ULTRAVIOLET_RANGE_NM[0]) & (middle < ULTRAVIOLET_RANGE_NM[1]) & np.isfinite(ultraviolet)
    in_visible = (middle >= VISIBLE_RANGE_NM[0]) & (middle < VISIBLE_RANGE_NM[1]) & np.isfinite(visible)
    result = np.where(in_ultraviolet, ultraviolet, result)
    return np.where(in_visible, visible, result)
