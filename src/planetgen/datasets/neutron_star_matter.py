# planetgen/datasets/neutron_star_matter.py

"""
Neutron-star matter: the equation of state and the outer-crust composition (GEN.548).

Two sets, both on Zenodo under CC BY 4.0:

- Twenty-one unified cold-matter equations of state, from the compilation of [R2834] (A&A 689,
  A242, [R4970]; Zenodo [R2834]): relativistic mean-field cores constrained by
  nuclear properties and chiral effective field theory, spanning a range of the symmetry-energy slope, with the BSk22
  outer crust and a compressible-liquid-drop inner crust. Each file lists baryon density (fm^-3), energy density and
  pressure (MeV fm^-3); a draw of the equation of state gives the spread of the maximum mass and radius.
- The ground-state outer-crust composition for twenty nuclear mass models, from [R4028] (A&A 710, A322,
  [R4971]; Zenodo [R4028]): for each model the sequence of nuclei with the
  baryon density and pressure at which one gives way to the next, ending at neutron drip. Only the unmagnetised
  (`B0`) files are read.

The readers raise `fetch.DataMissing` when the files have not been downloaded.
"""

import csv
import functools
import io
import re
import zipfile
from dataclasses import dataclass

import numpy as np

from planetgen.datasets import fetch

COMPOSITIONS = "compositions.zip"
EOS_COUNT = 21

_EOS_SHA = {
    1: ("4ead1aca30093a9a081b996d2f2c7bd117af4a540934cf9f7b47f6047973fb2d", 10432),
    2: ("cd156b18f711428bef861f89353ad1ef3b34c2b01faf1b5b53913ed82e2bd880", 10257),
    3: ("3c79ecefdd7fb3c5b7851adca02d0bcc038378e1fed531d398731d60daf94e5e", 10298),
    4: ("ed684d0ca92c6a2c6620a5523d9a857a375ac657f34af94bcec7cc3be9f71976", 10380),
    5: ("61ce386895ee892e48aa5bb02d58ac4605b953ba97033170491a2a70e713ea86", 10312),
    6: ("7cb7b1afe60ba0808f52746afd02d1d2ec2d681bf6d6b88e05526d1f45ffcb4f", 10416),
    7: ("0633432bf664254a2bfacd04bad198a3d5515a8e93344bbfe0144006d6b74c9d", 10396),
    8: ("f393d474a56068c43ae20c6561640e3b257edf7b71c96fb79865415f40d99d6b", 10414),
    9: ("ed8331968c856ff42913d3b203afcba193898d2e77b07f3d3eea6bf7ccad23bc", 10457),
    10: ("66cd99e3faae9f7b1486321ba1c8e259740d6f9750ba6c75bfbf4dc8a0235b5e", 10426),
    11: ("226f9ecf402afa29da8366df1b7a60d4c14c28675a03d72d7c48c249b0a7446e", 10327),
    12: ("fba0e6cb98f3d4809a23e92f5baa923654890414f1395ea6f25fe0776b139c06", 10415),
    13: ("87fde1cc032c9434381f151ece7e2bd9108d17639d17b318565a7e8c7390a4d0", 10337),
    14: ("266253f48b655bcecd2fa5efefb1c66d5617ad1693e51443831fd90e808314ca", 10436),
    15: ("dd5d26d6e9a6bd692ad3126d31f1be61d10e6bd6bee62556157d6e47451e3447", 10436),
    16: ("ac7d85a97d58c91bdd94e65064fac506053af1a6c83db7149212188c12e32182", 10415),
    17: ("83dff0536d23e7cfc309bed64c5f0649baddf2c0922e9f0daeb8d5b1c949163c", 10334),
    18: ("1af53f722ad5cf4beca89edb98e54eadcf19da3d37ab988da47cf1eeedf876d0", 10382),
    19: ("48572f99a5f9212f58867e888597656c90ce067d86b38feb3882ba0ba691305d", 10395),
    20: ("ba8a13174d0e0bf99928c64bb84fcc4b4816c0b65b3489ab49ce3800817b0193", 10420),
    21: ("dddeb91f7780dfeaa9e455ad934b7ecf8d4c9619f9f5feb6f4f18fc1ce7b6eca", 10480),
}

DATASET = fetch.Dataset(
    key="neutron_star_matter",
    title="Neutron-star equations of state [R2834] and outer-crust compositions [R4028]",
    source="[R2834, R4028], Zenodo",
    pieces=(fetch.Piece(COMPOSITIONS, "https://zenodo.org/records/20085135/files/compositions.zip?download=1",
                        "b98680f3b1d55786aea7db077e8be9e00f1949925d5b8858104eea4c9839a39c", 1830018),)
    + tuple(fetch.Piece(f"EOS{n}.txt", f"https://zenodo.org/records/10950421/files/EOS{n}.txt?download=1", sha, size)
            for n, (sha, size) in _EOS_SHA.items()),
    replaces="the nickel-to-krypton crust table and the 10 km radius floor of the stellar-composition note",
    fill="SC3",
    fallback="the one outer-crust sequence of [R0722: Table 2], in place of the 20 nuclear mass models",
)

_AMU_G = 1.66053906660e-24
"""float: The atomic mass unit in grams [R4361], the baryon mass for the density conversion."""

_FM3_PER_CM3 = 1e-39


@dataclass(frozen=True)
class EquationOfState:
    """One tabulated equation of state, ascending in density."""

    index: int
    baryon_density_fm3: np.ndarray
    energy_density_mev_fm3: np.ndarray
    pressure_mev_fm3: np.ndarray


@dataclass(frozen=True)
class CrustLayer:
    """One layer of the outer crust: its nucleus, and the baryon density and pressure at its bottom."""

    z: int
    a: int
    upper_density_fm3: float
    upper_pressure_mev_fm3: float


def parse_equation_of_state(index, text):
    """Reads one equation-of-state file: a `#` header line, then rows of density, energy density and pressure."""
    rows = [[float(token) for token in line.split()] for line in text.splitlines() if line.strip() and not line.startswith("#")]
    table = np.asarray(rows)
    if table.ndim != 2 or table.shape[1] != 3 or np.any(np.diff(table[:, 0]) <= 0):
        raise ValueError(f"EOS{index}: expected three columns with rising density")
    return EquationOfState(index, table[:, 0], table[:, 1], table[:, 2])


def parse_crust(text):
    """Reads one `composition_B0_<model>.csv`: the layers from the surface down to neutron drip."""
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows:
        raise ValueError("an outer-crust composition file with no rows")
    layers = []
    for row in rows:
        layers.append(CrustLayer(int(row["Z1"]), int(row["A1"]), float(row["n1max[fm-3]"]),
                                 float(row["P12[MeV fm-3]"])))
    return tuple(layers)


@functools.lru_cache(maxsize=4)
def _load_compositions(path):
    models = {}
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            match = re.search(r"([^/]+)/composition_B0_\1\.csv$", name)
            if match:
                models[match.group(1)] = parse_crust(archive.read(name).decode("utf-8"))
    if not models:
        raise ValueError("the archive holds no unmagnetised composition file")
    return models


@functools.lru_cache(maxsize=4)
def _load_eos(folder):
    tables = []
    for n in range(1, EOS_COUNT + 1):
        with open(f"{folder}/EOS{n}.txt", encoding="utf-8") as handle:
            tables.append(parse_equation_of_state(n, handle.read()))
    return tuple(tables)


def equations_of_state():
    """The twenty-one equations of state, by their number in the compilation (1 to 21)."""
    return _load_eos(fetch.locate(DATASET).require().folder)


def crust_models():
    """The nuclear mass models whose outer crust is tabulated, sorted by name."""
    return tuple(sorted(_load_compositions(fetch.locate(DATASET).require().path(COMPOSITIONS))))


def outer_crust(model):
    """The layers of `model`'s outer crust, surface first; the last layer ends at neutron drip."""
    return _load_compositions(fetch.locate(DATASET).require().path(COMPOSITIONS))[model]


def baryon_density_fm3(mass_density_g_cm3):
    """The baryon number density in fm^-3 of matter at `mass_density_g_cm3`."""
    return mass_density_g_cm3 / _AMU_G * _FM3_PER_CM3


def mass_density_g_cm3(baryon_density_fm3):
    """The mass density in g/cm^3 of matter at the baryon number density `baryon_density_fm3` (the inverse of
    `baryon_density_fm3`)."""
    return baryon_density_fm3 / _FM3_PER_CM3 * _AMU_G


def crust_nucleus(model, mass_density_g_cm3):
    """`(Z, A)` of the nucleus that fills the outer crust of `model` at `mass_density_g_cm3`, or None at and beyond
    neutron drip, where free neutrons appear and the nuclei are no longer the whole story."""
    density = baryon_density_fm3(mass_density_g_cm3)
    layers = outer_crust(model)
    for layer in layers:
        if density <= layer.upper_density_fm3:
            return layer.z, layer.a
    return None
