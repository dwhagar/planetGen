# tests/ir_datasets_support.py

"""Fixtures for the infrared datasets of GEN.715 (`datasets/ir_ktables.py`, `datasets/hitran_cia.py`), built offline: a small k-table
FITS file for each of the three gases with a made-up but monotone opacity (a line band per gas that grows with the pressure and the g
point), and a CIA file for each pair with a made-up coefficient. They go through the real download step with a fake opener (the sha256 of
the fixture bytes replaces the catalogue's), so no test reaches the network and no downloaded data is committed."""

import dataclasses
import gzip
import hashlib
import io

import numpy as np

from planetgen.datasets import fetch, hitran_cia, ir_ktables
from planetgen.physics import line_opacity

BLOCK = 2880
PRESSURES_BAR = np.logspace(-2.0, 2.0, 5)
TEMPERATURES_K = np.array([100.0, 300.0, 1000.0, 2500.0])
"""ndarray: The fixture's axes; 2500 K is above `ir_ktables.T_MAX_K`, so the reducer drops it."""
LINE_CENTRE_CM1 = {"H2O": 1600.0, "CO2": 667.0, "CH4": 1306.0}
LINE_WIDTH_CM1 = {"H2O": 300.0, "CO2": 60.0, "CH4": 100.0}
LINE_STRENGTH = {"H2O": 3.0e-19, "CO2": 1.0e-18, "CH4": 5.0e-19}
"""dicts: Where each fixture gas absorbs, how wide the band is, and its peak k (cm2 per molecule at one bar, 300 K, the strongest g point)."""
FLOOR = 1.0e-30
"""float: The k away from the band, cm2 per molecule."""


def _card(key, value):
    if isinstance(value, str):
        text = f"{key:<8}= '{value}'"
    elif isinstance(value, bool):
        text = f"{key:<8}= {'T' if value else 'F':>20}"
    else:
        text = f"{key:<8}= {value:>20}"
    return text.ljust(80).encode("ascii")


def _header(cards):
    raw = b"".join(_card(key, value) for key, value in cards) + b"END".ljust(80)
    return raw + b" " * ((-len(raw)) % BLOCK)


def _padded(data):
    return data + b"\0" * ((-len(data)) % BLOCK)


def wavelengths_cm():
    """The fixture's wavelength axis, cm, ascending, fine enough to put several samples in every band of the reduced grid."""
    return 1.0 / np.linspace(12300.0, 150.0, 1500)


def line_k(gas, pressure_bar, temperature_k, g_nodes, wavelength_cm):
    """The fixture's k (cm2 per molecule) for axes broadcast to (pressure, temperature, g, wavelength)."""
    nu = (1.0 / wavelength_cm)[None, None, None, :]
    shape = np.exp(-(((nu - LINE_CENTRE_CM1[gas]) / LINE_WIDTH_CM1[gas]) ** 2))
    strength = LINE_STRENGTH[gas] * (pressure_bar[:, None, None, None] ** 0.5) * (300.0 / temperature_k[None, :, None, None]) ** 0.5
    return FLOOR + strength * shape * 10.0 ** (-6.0 * (1.0 - g_nodes[None, None, :, None]))


def ktable_gz(gas):
    """An ARCiS-layout k-table FITS file, gzipped, for `gas`."""
    wavelength = wavelengths_cm()
    nodes, _weights = np.polynomial.legendre.leggauss(20)
    g_nodes = (nodes + 1.0) / 2.0
    k = line_k(gas, PRESSURES_BAR, TEMPERATURES_K, g_nodes, wavelength)
    n_p, n_t, n_g, n_l = k.shape
    primary = _header([("SIMPLE", True), ("BITPIX", -64), ("NAXIS", 4), ("NAXIS1", n_l), ("NAXIS2", n_g), ("NAXIS3", n_t), ("NAXIS4", n_p),
                       ("EXTEND", True), ("NP", n_p), ("NT", n_t), ("NG", n_g), ("NLAM", n_l)])
    body = _padded(k.astype(">f8").tobytes())
    out = primary + body
    for axis in (TEMPERATURES_K, PRESSURES_BAR, wavelength):
        out += _header([("XTENSION", "IMAGE"), ("BITPIX", -64), ("NAXIS", 1), ("NAXIS1", len(axis)), ("PCOUNT", 0), ("GCOUNT", 1)])
        out += _padded(axis.astype(">f8").tobytes())
    return gzip.compress(out, 1)


def cia_text(pair):
    """A CIA file for `pair` with two temperatures over 10 to 4500 cm-1; k falls with the wavenumber and the temperature."""
    first, second = pair.split("-")
    lines = []
    wavenumber = np.arange(10.0, 4500.0, 10.0)
    for temperature in (200.0, 300.0):
        k = 1.0e-45 * np.exp(-wavenumber / 1500.0) * (300.0 / temperature)
        lines.append(f"{first + '-' + second:>20}{wavenumber[0]:10.3f}{wavenumber[-1]:10.3f}{len(wavenumber):7d}{temperature:7.1f}{k.max():10.3E}{0.1:7.3f}")
        lines += [f" {nu:10.4f}  {value:10.3E}" for nu, value in zip(wavenumber, k)]
    return ("\n".join(lines) + "\n").encode("ascii")


def payloads():
    """`{piece name: bytes}` of every file of both datasets."""
    files = {name: ktable_gz(gas) for gas, (name, _url, _sha, _size) in ir_ktables.GASES.items()}
    files.update({name: cia_text(pair) for pair, (name, _sha, _size) in hitran_cia.FILES.items()})
    return files


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def install(monkeypatch, cache_folder, files=None):
    """
    Points the cache folder at `cache_folder`, replaces both datasets by copies whose pieces are the fixtures, fetches them and forgets
    everything the readers remembered. Returns the two fixture datasets.
    """
    files = payloads() if files is None else files
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(cache_folder))
    fetch._VERIFIED.clear()
    for loader in (ir_ktables.load, hitran_cia.load, line_opacity._build_tables, line_opacity._build_pairs):
        loader.cache_clear()
    installed = []
    for module in (ir_ktables, hitran_cia):
        dataset = module.DATASET
        pieces = tuple(dataclasses.replace(piece, url="https://example.org/" + piece.name, size=len(files[piece.name]),
                                           sha256=hashlib.sha256(files[piece.name]).hexdigest()) for piece in dataset.pieces)
        fixture = dataclasses.replace(dataset, pieces=pieces)
        monkeypatch.setattr(module, "DATASET", fixture)
        served = {piece.url: files[piece.name] for piece in pieces}
        result = fetch.fetch(fixture, opener=lambda request, timeout=None, served=served: _Response(served[request.full_url]))
        assert result.status == "fetched", result.reasons
        installed.append(fixture)
    return installed


def forget(monkeypatch, cache_folder):
    """Points the cache folder at an empty `cache_folder` (so both datasets are missing) and forgets what the readers remembered."""
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(cache_folder))
    fetch._VERIFIED.clear()
    for loader in (ir_ktables.load, hitran_cia.load, line_opacity._build_tables, line_opacity._build_pairs):
        loader.cache_clear()

