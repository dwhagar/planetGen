# tests/stellar_spectra_fixtures.py

"""Fixture files for the stellar spectra dataset of GEN.970 (`datasets/stellar_spectra.py`), built offline: an ATLAS9 flux file and PHOENIX
spectra whose models are blackbodies, so that a band's share has an exact answer."""

import re

import numpy as np

from planetgen.datasets import stellar_spectra
from tests import datasets_support as support

STELLAR_SPECTRA_KURUCZ_MODELS = ((3500.0, 4.0), (5000.0, 4.0), (5000.0, 4.5), (6000.0, 4.5), (10000.0, 4.0))
"""tuple: The (T_eff, log g) of the fixture's ATLAS9 models (the 3500 K dwarfs are the PHOENIX fixture's)."""


def _planck_per_cm(wavelength_cm, temperature):
    """The surface flux of a blackbody, erg/s/cm^2/cm, at the wavelengths."""
    h, c, k = 6.62607015e-27, 2.99792458e10, 1.380649e-16
    return np.pi * 2.0 * h * c ** 2 / wavelength_cm ** 5 / np.expm1(h * c / (wavelength_cm * k * temperature))


def _kurucz_text():
    """An ATLAS9 flux file whose models are blackbodies: the 22 program lines, 1221 wavelengths (8 to a line), then a title and
    153 lines each of H_nu and its continuum for every model."""
    wavelength_nm = np.round(np.geomspace(9.09, 160000.0, 1221), 2)
    lines = ["NEWODF GRID  [ 0.0]  VTURB 2.0 KM/S"] + ["C     program text"] * 21
    wavelength_lines = [wavelength_nm[i:i + 8] for i in range(0, 1221, 8)]
    lines += ["".join(f"{value:10.2f}" for value in row) for row in wavelength_lines]
    for teff, logg in STELLAR_SPECTRA_KURUCZ_MODELS:
        # F_nu = pi B_nu, H_nu = F_nu / 4pi = B_nu / 4;  B_nu = B_lambda lambda^2 / c
        wavelength_cm = wavelength_nm * 1e-7
        hnu = np.maximum(_planck_per_cm(wavelength_cm, teff) / np.pi * wavelength_cm ** 2 / 2.99792458e10 / 4.0, 1e-99)
        lines.append(f"TEFF {teff:7.0f}.  GRAVITY {logg:.5f} LTE  fixture")
        for values in (hnu, hnu):
            lines += ["".join(f"{value:10.4E}" for value in values[i:i + 8]) for i in range(0, 1221, 8)]
    return ("\n".join(lines) + "\n").encode("ascii")


def _fits_vector(values, bits, **cards):
    """A FITS file of one primary vector, big-endian floats."""
    header = [("SIMPLE", "T"), ("BITPIX", str(bits)), ("NAXIS", "1"), ("NAXIS1", str(len(values)))] + [(k, repr(v) if isinstance(v, str) else str(v)) for k, v in cards.items()]
    text = "".join(f"{key:<8}= {value:>20}".ljust(80) for key, value in header) + "END".ljust(80)
    text = text.ljust(-(-len(text) // 2880) * 2880)
    data = np.asarray(values, dtype=f">f{bits // -8}").tobytes()
    return text.encode("ascii") + data + b"\0" * (-len(data) % 2880)


def _stellar_spectra_payloads():
    wavelength_a = np.geomspace(500.0, 55000.0, 4000)
    payloads = {stellar_spectra.KURUCZ: _kurucz_text(), stellar_spectra.PHOENIX_WAVE: _fits_vector(wavelength_a, -64)}
    for piece in stellar_spectra.DATASET.pieces:
        match = re.fullmatch(r"lte(\d{5})-(\d\.\d\d)-0\.0\.PHOENIX.*", piece.name)
        if match:
            teff, logg = float(match.group(1)), float(match.group(2))
            flux = _planck_per_cm(wavelength_a * 1e-8, teff).astype(np.float32)
            payloads[piece.name] = _fits_vector(flux, -32, PHXTEFF=teff, PHXLOGG=logg)
    return payloads


def install_stellar_spectra(monkeypatch):
    return support.install(monkeypatch, stellar_spectra, stellar_spectra.DATASET, _stellar_spectra_payloads())
