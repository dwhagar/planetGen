# planetgen/datasets/ir_ktables.py

"""
Correlated-k absorption tables of the infrared-active gases (GEN.715, AD10): the opacity of water vapour, carbon dioxide and methane
from line lists, at every pressure and temperature the greenhouse model needs.

Source: the ExoMolOP opacity tables of the ExoMol project [R0795], the
ARCiS-format correlated-k tables at resolving power 1000 from 0.3 to 50 micron, 22 pressures (1e-5 to 100 bar), 27 temperatures
(100 to 3400 K) and 20 Gauss-Legendre g points per wavelength:

- H2O from the POKAZATEL line list [R3571],
- CO2 from the UCL-4000 line list [R4854],
- CH4 from the YT34to10 line list [R4855].

Each file is 370 to 465 MB. The three are transient pieces: setup downloads them once, checks their sha256, reduces them to the
129 bands of `BAND_EDGES_CM1` (20 cm-1 wide from 10 to 5 micron, then 50, 100 and 250 cm-1 wide) and 10 g points (six under
`G_SPLIT`, four over it, so the strongest 10 percent of the lines are resolved) on the temperatures to `T_MAX_K`, writes the
2 MB each into `ir_ktables.npz`, and removes the archives. The raw tables are never kept.

The reduction is the band mean of the k-distribution: all the (wavenumber, g) values of a band are sorted by k, the cumulative
opacity is cut at the new g edges, and each interval's mean is the new k (this keeps the band's total absorption).

Limits (stated, and in the stand-in registry of docs/design/research/gap-filling-plan.md section 6): the tables cut each line at 25 cm-1 from its
centre, so the far wings of a dense carbon dioxide air (Venus) are missing (`co2_far_wing` adds them in `line_opacity`); the
broadening is the tables' own (air-like for H2O and CO2, so an air of hydrogen and helium is under-broadened); the temperature grid
holds 100 K as its lowest point and 100 bar as its highest pressure, beyond which a lookup holds the edge value.

The gas table's band centres and `line_opacity`'s plateaus are the labelled fallback for a run that disregards missing data
(`fetch.Located.use_stand_in`).
"""

import functools
import gzip
import io

import numpy as np

from planetgen.datasets import fetch

BLOCK_BYTES = 2880
"""int: The size of a FITS block."""

DERIVED = "ir_ktables.npz"
SOURCE = "ExoMolOP R1000 correlated-k tables [R0795]"

GASES = {
    "H2O": ("ir_ktable_h2o.fits.gz", "https://exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__R1000_0.3-50mu.ktable.ARCiS.fits.gz",
            "1918fd739ada3f25a930d2794157df2668def25e330543d3a717b2afec8b469a", 465284416),
    "CO2": ("ir_ktable_co2.fits.gz", "https://exomol.com/db/CO2/12C-16O2/UCL-4000/12C-16O2__UCL-4000.R1000_0.3-50mu.ktable.ARCiS.fits.gz",
            "2c4134cf478818590ca544f18d3906ac931d00e79aa1b639e4b93d90788b3682", 421823286),
    "CH4": ("ir_ktable_ch4.fits.gz", "https://exomol.com/db/CH4/12C-1H4/YT34to10/12C-1H4__YT34to10.R1000_0.3-50mu.ktable.ARCiS.fits.gz",
            "7d8fd39d9247c6710053cd38b18083f7114fdd90c180ee9b70fd15e55d520b20", 370707002),
}
"""dict: Each gas's formula to (file name in the cache, URL, sha256, size in bytes)."""

BAND_EDGES_CM1 = np.concatenate([np.arange(200.0, 1000.0, 20.0), np.arange(1000.0, 2500.0, 50.0), np.arange(2500.0, 6000.0, 100.0),
                                 np.arange(6000.0, 12001.0, 250.0)])
"""ndarray: The 130 edges, cm-1, of the 129 bands the tables are reduced to (200 cm-1 = 50 micron, 12000 cm-1 = 0.83 micron)."""

G_POINTS = (6, 4)
"""tuple: Gauss-Legendre points under and over `G_SPLIT`."""
G_SPLIT = 0.9
"""float: The cumulative share where the g grid changes from its coarse part to its fine part (the strongest lines)."""

T_MAX_K = 2000.0
"""float: The hottest temperature kept (the 100 K steps up to it; the tables go to 3400 K)."""


def g_grid():
    """The reduced g grid: `(g, weights, edges)`: mid-points of the g intervals, their weights (they add up to 1) and the 11 edges."""
    low, high = G_POINTS
    x, w = np.polynomial.legendre.leggauss(low)
    lower_g, lower_w = (x + 1.0) / 2.0 * G_SPLIT, w / 2.0 * G_SPLIT
    x, w = np.polynomial.legendre.leggauss(high)
    upper_g, upper_w = G_SPLIT + (x + 1.0) / 2.0 * (1.0 - G_SPLIT), w / 2.0 * (1.0 - G_SPLIT)
    g = np.concatenate([lower_g, upper_g])
    weights = np.concatenate([lower_w, upper_w])
    return g, weights, np.concatenate([[0.0], np.cumsum(weights)])


def _read_header(stream):
    """The cards of the next FITS header in `stream` as a dict (numbers parsed), or `None` at the end of the file."""
    cards = {}
    while True:
        block = stream.read(BLOCK_BYTES)
        if len(block) < BLOCK_BYTES:
            return None
        for start in range(0, BLOCK_BYTES, 80):
            card = block[start:start + 80].decode("ascii", "replace")
            if card.startswith("END"):
                return cards
            if card[8:10] == "= ":
                value = card[10:].split("/")[0].strip()
                try:
                    cards[card[:8].strip()] = int(value)
                except ValueError:
                    try:
                        cards[card[:8].strip()] = float(value)
                    except ValueError:
                        cards[card[:8].strip()] = value.strip("'").strip()


def _shape(header):
    return [header[f"NAXIS{i + 1}"] for i in range(header["NAXIS"])][::-1]


def _skip_padding(stream, nbytes):
    stream.read((-nbytes) % BLOCK_BYTES)


def read_ktable(stream):
    """
    Reads an ARCiS-format k-table FITS stream (gzip already open): the primary array of shape (pressure, temperature, g, wavelength)
    as float32 and the three extension axes (temperature K, pressure bar, wavelength cm).

    Raises:
        ValueError: The file is not the layout described above.
    """
    header = _read_header(stream)
    if header is None or header.get("NAXIS") != 4 or header.get("BITPIX") != -64:
        raise ValueError("the k-table is not a four-dimensional array of 64-bit floats")
    n_p, n_t, n_g, n_l = _shape(header)
    if (n_p, n_t, n_g, n_l) != (header.get("NP"), header.get("NT"), header.get("NG"), header.get("NLAM")):
        raise ValueError("the k-table's axes do not match its NP, NT, NG and NLAM cards")
    data = np.empty((n_p, n_t, n_g, n_l), dtype=np.float32)
    row_bytes = n_t * n_g * n_l * 8
    for ip in range(n_p):
        raw = stream.read(row_bytes)
        if len(raw) != row_bytes:
            raise ValueError("the k-table ends before its data does")
        data[ip] = np.frombuffer(raw, dtype=">f8").reshape(n_t, n_g, n_l)
    _skip_padding(stream, n_p * row_bytes)
    axes = []
    for expected in (n_t, n_p, n_l):
        extension = _read_header(stream)
        if extension is None or _shape(extension) != [expected]:
            raise ValueError("the k-table's axis extensions are not (temperature, pressure, wavelength)")
        raw = stream.read(expected * 8)
        if len(raw) != expected * 8:
            raise ValueError("the k-table ends inside an axis")
        axes.append(np.frombuffer(raw, dtype=">f8").astype(float))
        _skip_padding(stream, expected * 8)
    temperature, pressure, wavelength = axes
    if np.any(np.diff(temperature) <= 0) or np.any(np.diff(pressure) <= 0) or np.any(np.diff(wavelength) <= 0):
        raise ValueError("the k-table's axes must rise")
    return data, temperature, pressure, wavelength


def reduce_distribution(k, weights, edges):
    """
    The band k-distribution on the g `edges`: `k` is (..., n) opacities with (n,) `weights` (adding to 1); they are sorted ascending
    along the last axis, the cumulative integral of k over the weight is cut at `edges`, and each interval's mean is returned (..., len(edges)-1).
    """
    order = np.argsort(k, axis=-1)
    k_sorted = np.take_along_axis(k, order, axis=-1)
    w_sorted = weights[order]
    cumulative_w = np.concatenate([np.zeros(w_sorted.shape[:-1] + (1,)), np.cumsum(w_sorted, axis=-1)], axis=-1)
    cumulative_kw = np.concatenate([np.zeros(w_sorted.shape[:-1] + (1,)), np.cumsum(k_sorted * w_sorted, axis=-1)], axis=-1)
    size = cumulative_w.shape[-1]
    at_edges = np.empty(k.shape[:-1] + (len(edges),))
    for j, edge in enumerate(edges):
        index = np.clip((cumulative_w < edge - 1.0e-12).sum(axis=-1), 1, size - 1)
        w0 = np.take_along_axis(cumulative_w, (index - 1)[..., None], -1)[..., 0]
        w1 = np.take_along_axis(cumulative_w, index[..., None], -1)[..., 0]
        k0 = np.take_along_axis(cumulative_kw, (index - 1)[..., None], -1)[..., 0]
        k1 = np.take_along_axis(cumulative_kw, index[..., None], -1)[..., 0]
        share = np.where(w1 > w0, (edge - w0) / np.maximum(w1 - w0, 1.0e-300), 0.0)
        at_edges[..., j] = k0 + share * (k1 - k0)
    return np.diff(at_edges, axis=-1) / np.diff(edges)


def reduce_table(data, temperature, pressure, wavelength, edges_cm1=None):
    """
    Reduces a table read by `read_ktable` to the bands of `BAND_EDGES_CM1` and the g grid of `g_grid`.

    Returns:
        tuple: `(k, temperature_k, pressure_bar)` with `k` of shape (band, temperature, pressure, g) in cm2 per molecule, float32,
        the temperatures kept up to `T_MAX_K`.
    """
    edges_cm1 = BAND_EDGES_CM1 if edges_cm1 is None else edges_cm1
    keep = temperature <= T_MAX_K
    data = data[:, keep]
    temperature = temperature[keep]
    wavenumber = 1.0 / wavelength
    step = np.abs(np.gradient(wavenumber))
    _nodes, node_weights = np.polynomial.legendre.leggauss(data.shape[2])
    node_weights = node_weights / 2.0
    _g, _w, g_edges = g_grid()
    result = np.full((len(edges_cm1) - 1, len(temperature), len(pressure), len(g_edges) - 1), np.nan, dtype=np.float32)
    for band in range(len(edges_cm1) - 1):
        inside = (wavenumber >= edges_cm1[band]) & (wavenumber < edges_cm1[band + 1])
        if not inside.any():
            raise ValueError(f"the k-table has no wavelength in the band {edges_cm1[band]} to {edges_cm1[band + 1]} cm-1")
        sub = data[:, :, :, inside].astype(np.float64)                      # (P, T, G, n)
        k = sub.transpose(0, 1, 3, 2).reshape(len(pressure), len(temperature), -1)
        weights = np.tile(node_weights, int(inside.sum())) * np.repeat(step[inside], len(node_weights))
        reduced = reduce_distribution(k, weights / weights.sum(), g_edges)  # (P, T, g)
        result[band] = reduced.transpose(1, 0, 2)
    if np.isnan(result).any():
        raise ValueError("the reduced k-table has gaps")
    return result, temperature, pressure


def derive(folder):
    """Writes `ir_ktables.npz` from the verified archives in `folder`: `k_<gas>` (band, temperature, pressure, g), shared `temperature_k`,
    `pressure_bar`, `g`, `weights` and `edges_cm1`."""
    g, weights, _edges = g_grid()
    arrays = {"g": g, "weights": weights, "edges_cm1": BAND_EDGES_CM1}
    for gas, (name, _url, _sha, _size) in GASES.items():
        with gzip.open(f"{folder}/{name}", "rb") as stream:
            data, temperature, pressure, wavelength = read_ktable(stream)
        k, temperature_k, pressure_bar = reduce_table(data, temperature, pressure, wavelength)
        if "temperature_k" in arrays and (not np.array_equal(arrays["temperature_k"], temperature_k)
                                          or not np.array_equal(arrays["pressure_bar"], pressure_bar)):
            raise ValueError(f"{gas}: the table's temperature and pressure axes differ from the others'")
        arrays["temperature_k"], arrays["pressure_bar"] = temperature_k, pressure_bar
        arrays[f"k_{gas}"] = k
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    return {DERIVED: buffer.getvalue()}


DATASET = fetch.Dataset(
    key="ir_ktables",
    title="Infrared correlated-k tables of H2O, CO2 and CH4 (ExoMolOP, R1000, line lists POKAZATEL, UCL-4000 and YT34to10)",
    source="[R0795], ExoMolOP R1000 ARCiS k-tables from exomol.com",
    pieces=tuple(fetch.Piece(name, url, sha, size, transient=True) for name, url, sha, size in GASES.values()),
    replaces="the band-strength and wing stand-ins of the grey infrared opacity: plateau opacities per band centre of the gas table",
    fill="AD10",
    fallback="none at generation: the fitted greenhouse table (physics/greenhouse_table.py) carries the answer; the tables are read to fit and check it",
    derived=(DERIVED,),
    derive=lambda folder: derive(folder),
    derive_version=1,
)


@functools.lru_cache(maxsize=2)
def load(path):
    with np.load(path) as stored:
        return {name: np.array(stored[name]) for name in stored.files}


def tables(located):
    """The reduced tables as a dict of arrays (`k_H2O`, `k_CO2`, `k_CH4`, `temperature_k`, `pressure_bar`, `g`, `weights`, `edges_cm1`)
    from a `fetch.Located`."""
    return load(located.path(DERIVED))

