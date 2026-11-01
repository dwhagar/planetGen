# planetgen/datasets/stellar_spectra.py

"""
Model-atmosphere spectra of stars by effective temperature and gravity (GEN.970, parent GEN.745).

Two grids, both solar metallicity, read as downloaded and reduced at setup:

- **[R0681] ATLAS9 fluxes** (the `fp00k2odfnew` grid, [M/H] = 0, microturbulence 2 km/s, mixing length
  1.25, new opacity distribution functions; in IAU Symp. 210, poster A20, [R0681]), from Fiorella Castelli's
  page at INAF Trieste: 476 models from 3500 to 50000 K and log g 0 to 5, 1221 wavelengths from 9 nm to 160 micron.
- **[R2058] PHOENIX-ACES-AGSS-COND-2011** (A&A 553, A6, [R2058]), high-resolution spectra from the
  Goettingen server, for the stars ATLAS9 stops short of: 2300 to 3500 K at 200 K steps and log g 4.5, 5.0 and 5.5 (the
  M dwarfs), 500 A to 5.5 micron.

The pieces are transient: setup downloads them (130 MB), checks their sha256, reduces each model to the **share of its
bolometric power in each of 1500 bins** equally spaced in the logarithm of the wavelength from 9 nm to 100 micron
(`stellar_spectra.npz`, a few MB), and removes them. A model's share is its flux in the bin over sigma T^4; the part of the bin
range a grid does not reach (below 9 nm or 50 nm, above 5.5 micron for PHOENIX, 160 micron for ATLAS9) takes the blackbody of the
model's temperature, and the shares are renormalised to sum to one (the correction is stored as `closure` and the reduction
refuses a model it moves by more than 10 percent).

What it answers: the share of a star's bolometric power in a wavelength band, for a star of any effective temperature and
gravity inside the grids (`band_fractions`), and the whole share vector (`spectrum`). Between models the shares are
interpolated linearly in log10 T_eff and in log g (their logarithms, so that an ultraviolet tail that falls by a decade between two
models is a straight step). Outside the grids the lookup answers `None`, and the caller takes the blackbody.

Limits, each a stand-in of the registry (docs/design/research/gap-filling-plan.md section 6.33): the models hold no chromosphere
or corona, so the ultraviolet of a cool or active star is the photosphere's alone (`physics.stellar_library` adds the measured
emission of MUSCLES); the PHOENIX grid is dwarfs only (log g 4.5 to 5.5), so a giant cooler than 3500 K reads the log g 4.5 model; the grids are solar metallicity.
"""

import functools
import io
import math
import re

import numpy as np

from planetgen.datasets import fetch
from planetgen.physics import stellar_spectrum

DERIVED = "stellar_spectra.npz"
KURUCZ = "fp00k2odfnew.pck"
PHOENIX_WAVE = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"

EDGES_NM = np.geomspace(9.0, 1.0e5, 1501)
"""ndarray: The 1501 bin edges, nm, equally spaced in log wavelength."""

SIGMA_CGS = 5.670374419e-5
C_CGS = 2.99792458e10
CLOSURE_LIMIT = 0.10
"""float: A model whose integrated flux is more than this from sigma T^4 is refused by the reduction."""

_KURUCZ_URL = "https://wwwuser.oats.inaf.it/fiorella.castelli/grids/gridp00k2odfnew/" + KURUCZ
_PHOENIX_URL = "https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/"
_PHOENIX_TEFF = (2300, 2500, 2700, 2900, 3100, 3300, 3500)
_PHOENIX_LOGG = (4.5, 5.0, 5.5)


def _phoenix_name(teff, logg):
    return f"lte{teff:05d}-{logg:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"


_PHOENIX_SHA = {
    "WAVE": ("9befb317a6d37bf8cf031e9ef1ba8a942f2f13786c953348f437aaffadf60f73", 12556800),
    (2300, 4.5): "ba86fee94209dccedf32c25f7076586b89118d4948b87c64491cc5810b128b6b",
    (2300, 5.0): "13055a4de101ce814b5b82ca7c33fa54aed21257daf9d7e548d410c1edb5d8aa",
    (2300, 5.5): "5d7f0018261f9601687453e14bec4329c36849f5c8955cbc740758966ab48fa3",
    (2500, 4.5): "2789c43690111255feb7b224052c77bb413e4189a503d6344e889343e78883bf",
    (2500, 5.0): "342e38064fddbdd9ea6782fb1a467634d35d099e9e5978a0cb5e36cf4ea23675",
    (2500, 5.5): "411f10314c760e46b3c4d730b808a1ca47d94ccc4b1a45f0a8e9311f8a599837",
    (2700, 4.5): "466686b6b45ea0200a6356109f4c9e4a63537b090f55e639f1c7aa14bccd1237",
    (2700, 5.0): "3c7e3759e834331d0ff482bf130a5c6651fd88d6fcd8fc24ff1dcc5ef4ab3382",
    (2700, 5.5): "4c3078d777c11fc64ca715b501f1dba38ce3d23b4945ad994036ac7a66296707",
    (2900, 4.5): "30cd42862558bfdc935a6e5d9cc9db3ab5c010d72b7c35a782b5476c4ed4731e",
    (2900, 5.0): "0fe6ea037adafae3f48a024579a97a26b1823b81fdd373b024eccbf327a06cfc",
    (2900, 5.5): "09bd2295bdc4952e2500b201e86337b3ec4a8a0d36a64e1c38ff3148a73e7c64",
    (3100, 4.5): "7c9b450e5ee4eb6f947b3a1227355612254325c5d6f6702f9080b802297bbbfc",
    (3100, 5.0): "f776354b070a47422f4db50f178fca51a008f8eb91c08a1a193d09f25776197e",
    (3100, 5.5): "1b2e307ac8abaf041f0b143f507965fb798417bf7b1af825371dd8842bd6c299",
    (3300, 4.5): "ac43a2768adba14337414c03c461873abe6d4fb59d091ed870360b8a5b377378",
    (3300, 5.0): "50b8b3cd200b4f34fa6fba09b5075dcb696e8a4e44f2f5e88923dad1beefcfd6",
    (3300, 5.5): "3ed65e4e3d83ced566d84ffbc7bce726f967e87958467ef8809f03c23cad3749",
    (3500, 4.5): "96a687f9a4dbd6f0b657392cf5e56842999322c2f5592041527c41cd4f75a234",
    (3500, 5.0): "c272c0b4478f9bf53d6b821dedd5d30c4bc3e40eaadf40c18c32de345cb5b1ca",
    (3500, 5.5): "f621bf6bf3d4cc3c92a1d0177e96a253144e220c59ac4c3952a7b881d8dd0f68",
}
_PHOENIX_SPECTRUM_BYTES = 6287040

PIECES = (
    (fetch.Piece(KURUCZ, _KURUCZ_URL, "dc935dc09f83064a09f282a024053ac7815ba47c6600161e4e1ed96871ec7eea", 11821049, transient=True),
     fetch.Piece(PHOENIX_WAVE, _PHOENIX_URL + PHOENIX_WAVE, *_PHOENIX_SHA["WAVE"], transient=True))
    + tuple(fetch.Piece(_phoenix_name(t, g), f"{_PHOENIX_URL}PHOENIX-ACES-AGSS-COND-2011/Z-0.0/{_phoenix_name(t, g)}",
                        _PHOENIX_SHA[(t, g)], _PHOENIX_SPECTRUM_BYTES, transient=True)
            for t in _PHOENIX_TEFF for g in _PHOENIX_LOGG))


# --- reading the sources ---------------------------------------------------------------------------------------------------

def read_kurucz(text):
    """The models of an ATLAS9 `fp*.pck` flux file: `(wavelength_nm, [(teff, logg, hnu), ...])`, `hnu` in erg/cm^2/s/Hz/sr."""
    lines = text.splitlines()
    skip = 22
    wave_lines = math.ceil(1221 / 8)
    wave = np.array([float(line[i:i + 10]) for line in lines[skip:skip + wave_lines] for i in range(0, len(line.rstrip()), 10)])
    if wave.size != 1221:
        raise ValueError(f"the ATLAS9 wavelength table holds {wave.size} values, not 1221")
    flux_lines = math.ceil(1221 / 8)
    models, at = [], skip + wave_lines
    while at < len(lines) and lines[at].strip():
        match = re.search(r"TEFF\s+([\d.]+)\s+GRAVITY\s+([\d.]+)", lines[at])
        if not match:
            raise ValueError(f"an ATLAS9 model title is expected at line {at + 1}, found {lines[at]!r}")
        block = lines[at + 1:at + 1 + flux_lines]
        hnu = np.array([float(line[i:i + 10]) for line in block for i in range(0, len(line.rstrip()), 10)])
        if hnu.size != 1221:
            raise ValueError(f"model {match.group(0)}: {hnu.size} fluxes, not 1221")
        models.append((float(match.group(1)), float(match.group(2)), hnu))
        at += 1 + 2 * flux_lines
    return wave, models


def read_fits_vector(data):
    """The one-dimensional primary array of a FITS file (bytes) and its header cards as a dict."""
    cards, at = {}, 0
    while True:
        block = data[at:at + 2880]
        at += 2880
        for i in range(0, 2880, 80):
            card = block[i:i + 80].decode("ascii")
            key = card[:8].strip()
            if key == "END":
                count, bits = int(cards["NAXIS1"]), int(cards["BITPIX"])
                return np.frombuffer(data, dtype=f">f{abs(bits) // 8}", count=count, offset=at).astype(float), cards
            if card[8:10] == "= ":
                cards[key] = card[10:].split("/")[0].strip()
        if not block:
            raise ValueError("a FITS header with no END card")


def _bin_powers(wavelength_cm, flux_per_cm):
    """The power, erg/s/cm^2, of a spectrum in each of `EDGES_NM`'s bins: the trapezoid cumulative integral, read at the edges."""
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (flux_per_cm[1:] + flux_per_cm[:-1]) * np.diff(wavelength_cm))])
    return np.diff(np.interp(EDGES_NM * 1e-7, wavelength_cm, cumulative, left=0.0, right=cumulative[-1]))


def _shares(wavelength_nm, flux_per_cm, teff):
    """`(shares, closure)` of one model: the share of sigma T^4 in each bin, the blackbody filling the bins the model does not
    reach, normalised to one; and the integrated flux over sigma T^4 before that."""
    wavelength_cm = np.asarray(wavelength_nm) * 1e-7
    power = _bin_powers(wavelength_cm, np.asarray(flux_per_cm))
    reached = (EDGES_NM[1:] > wavelength_nm[0]) & (EDGES_NM[:-1] < wavelength_nm[-1])
    covered = power.sum() / (SIGMA_CGS * teff ** 4)
    shares = power / (SIGMA_CGS * teff ** 4)
    outside = stellar_spectrum.blackbody_fraction(EDGES_NM[:-1], EDGES_NM[1:], teff)
    shares = np.where(reached, shares, outside)
    closure = float(shares.sum())
    return shares / closure, covered


def derive(folder):
    """Writes `stellar_spectra.npz` from the verified pieces in `folder`: `edges_nm`, `teff_k`, `logg`, `source` (0 ATLAS9, 1 PHOENIX),
    `closure` and `shares` (model, bin)."""
    with open(f"{folder}/{KURUCZ}", encoding="ascii") as handle:
        wave_nm, kurucz = read_kurucz(handle.read())
    rows = []
    for teff, logg, hnu in kurucz:
        if teff == _PHOENIX_TEFF[-1] and logg in _PHOENIX_LOGG:
            continue                     # the 3500 K dwarfs are PHOENIX's: one model for each temperature and gravity
        flux = 4.0 * math.pi * hnu * C_CGS / (wave_nm * 1e-7) ** 2
        rows.append((teff, logg, 0) + _shares(wave_nm, flux, teff))
    with open(f"{folder}/{PHOENIX_WAVE}", "rb") as handle:
        wave_a, _ = read_fits_vector(handle.read())
    for teff in _PHOENIX_TEFF:
        for logg in _PHOENIX_LOGG:
            with open(f"{folder}/{_phoenix_name(teff, logg)}", "rb") as handle:
                flux, cards = read_fits_vector(handle.read())
            if float(cards["PHXTEFF"]) != teff or float(cards["PHXLOGG"]) != logg:
                raise ValueError(f"{_phoenix_name(teff, logg)} is the model of {cards['PHXTEFF']} K, log g {cards['PHXLOGG']}")
            rows.append((float(teff), logg, 1) + _shares(wave_a * 0.1, flux, teff))
    for teff, logg, _source, _shares_row, covered in rows:
        if abs(covered - 1.0) > CLOSURE_LIMIT:
            raise ValueError(f"the model of {teff:.0f} K, log g {logg}: its flux is {covered:.3f} of sigma T^4")
    order = sorted(range(len(rows)), key=lambda i: (rows[i][0], rows[i][1]))
    buffer = io.BytesIO()
    np.savez_compressed(
        buffer, edges_nm=EDGES_NM, teff_k=np.array([rows[i][0] for i in order]), logg=np.array([rows[i][1] for i in order]),
        source=np.array([rows[i][2] for i in order], dtype=np.int8), closure=np.array([rows[i][4] for i in order]),
        shares=np.array([rows[i][3] for i in order], dtype=np.float32))
    return {DERIVED: buffer.getvalue()}


DATASET = fetch.Dataset(
    key="stellar_spectra",
    title="Model-atmosphere spectra: Castelli and Kurucz ATLAS9 fluxes (3500 to 50000 K) and Husser et al. PHOENIX-ACES (2300 to 3500 K)",
    source="[R0681] (IAU Symp. 210, A20), [R2058]",
    pieces=PIECES,
    replaces="the blackbody of the stellar spectrum with a seeded UV-C share (R17) and the BT-Settl chromaticity offsets",
    fill="SL1",
    fallback="the blackbody of the effective temperature with the seeded UV shares of physics/stellar_spectrum.py",
    derived=(DERIVED,),
    derive=lambda folder: derive(folder),
    derive_version=1,
)


# --- the lookup --------------------------------------------------------------------------------------------------------------

@functools.lru_cache(maxsize=2)
def load(path):
    """The reduced library as a dict of arrays (see `derive`)."""
    with np.load(path) as stored:
        return {name: np.array(stored[name]) for name in stored.files}


def library(located=None):
    """The reduced library, from a `fetch.Located` or, when none is given, from the cache folder (raises `fetch.DataMissing`)."""
    located = located or fetch.locate(DATASET).require()
    return load(located.path(DERIVED))


_FLOOR = 1e-30


def _row(rows, logg, lib):
    """The log-shares of the models of one effective temperature interpolated to `logg`, clamped to the gravities it has."""
    gravities = lib["logg"][rows]
    shares = np.log(np.maximum(lib["shares"][rows].astype(float), _FLOOR))
    g = min(max(logg, gravities[0]), gravities[-1])
    upper = int(np.searchsorted(gravities, g))
    if upper == 0 or gravities[upper] == g:
        return shares[upper]
    t = (g - gravities[upper - 1]) / (gravities[upper] - gravities[upper - 1])
    return shares[upper - 1] * (1.0 - t) + shares[upper] * t


def spectrum(teff_k, logg, lib=None):
    """The share of a star's bolometric power in each of `EDGES_NM`'s bins, interpolated from the models, or None when the
    effective temperature is outside 2300 to 50000 K. A gravity outside a temperature's models takes the nearest."""
    lib = lib or library()
    temperatures = np.unique(lib["teff_k"])
    if not temperatures[0] <= teff_k <= temperatures[-1]:
        return None
    # at 3500 K the rows are ATLAS9's below log g 4.5 and PHOENIX's from there up; below 3500 K only PHOENIX's dwarfs are held
    upper = int(np.searchsorted(temperatures, teff_k))
    if upper == 0 or temperatures[upper] == teff_k:
        rows = np.flatnonzero(lib["teff_k"] == temperatures[upper])
        return _normalised(np.exp(_row(rows, logg, lib)))
    low, high = temperatures[upper - 1], temperatures[upper]
    t = (math.log10(teff_k) - math.log10(low)) / (math.log10(high) - math.log10(low))
    mix = (1.0 - t) * _row(np.flatnonzero(lib["teff_k"] == low), logg, lib) + t * _row(np.flatnonzero(lib["teff_k"] == high), logg, lib)
    return _normalised(np.exp(mix))


def _normalised(shares):
    return shares / shares.sum()


def band_fractions(teff_k, logg, bands_nm, lib=None):
    """The share of a star's bolometric power in each `(low, high)` band, nm (`high` may be infinity), or None outside the
    library's temperatures. Inside a bin the power is spread evenly in the logarithm of the wavelength."""
    shares = spectrum(teff_k, logg, lib)
    if shares is None:
        return None
    cumulative = np.concatenate([[0.0], np.cumsum(shares)])
    log_edges = np.log(EDGES_NM)

    def below(wavelength):
        if wavelength >= EDGES_NM[-1]:
            return 1.0
        if wavelength <= EDGES_NM[0]:
            return 0.0
        return float(np.interp(math.log(wavelength), log_edges, cumulative))

    return tuple(below(high) - below(low) for low, high in bands_nm)
