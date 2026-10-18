# planetgen/datasets/white_dwarf_models.py

"""
White-dwarf cooling sequences with their layers (GEN.548).

Source: Bédard, Bergeron, Brassard and Fontaine 2020, ApJ 901, 93 (doi:10.3847/1538-4357/abafbe, arXiv:2008.07469),
23 masses from 0.20 to 1.30 solar masses for a thick (q_He = 1e-2, q_H = 1e-4) and a thin (q_He = 1e-2, q_H = 1e-10)
hydrogen layer, on Antoine Bédard's page at the University of Warwick (`sequences.tar`, 2.5 MB, 46 text files). Each
model lists the effective temperature, log g, radius, cooling age, photon luminosity, central conditions, the
crystallised fraction, the neutrino luminosity and the log of the total mass fractions of hydrogen, helium, carbon and
oxygen. The core of these models is homogeneous and equimassic carbon and oxygen: the C:O split of 1:1 is the
sequences' own assumption, not a measurement, and the readers say so. The sequences are carbon-oxygen white dwarfs;
a helium-core or an oxygen-neon-core white dwarf is not in them.

The readers raise `fetch.StandIn` when the file has not been downloaded, and answer `None` outside the masses
(0.20 to 1.30) of the grid. Temperatures outside a sequence's range are clamped to its ends.
"""

import functools
import re
import tarfile
from dataclasses import dataclass

import numpy as np

from planetgen.datasets import fetch
from planetgen.physics import constants

ARCHIVE = "sequences.tar"

DATASET = fetch.Dataset(
    key="white_dwarf_models",
    title="Bédard et al. 2020 white-dwarf cooling sequences",
    source="Bédard, Bergeron, Brassard and Fontaine 2020, ApJ 901, 93 (doi:10.3847/1538-4357/abafbe)",
    licence="open; cite Bédard et al. 2020",
    pieces=(fetch.Piece(ARCHIVE, "https://warwick.ac.uk/fac/sci/physics/research/astro/people/antoinebedard/models/"
                        + ARCHIVE, "5c99114fcba0087395fcde2779ad7bdc9b7ace224abb71322d45165facae45ac", 2519040),),
    standin="the Chandrasekhar-type law of the stellar-composition note, an even carbon-oxygen split and a fixed layer split",
    fill="SC2",
)

ENVELOPES = ("thick", "thin")
"""tuple: The two hydrogen-layer cases of the sequences."""

_NAME = re.compile(r"seq_(\d{3})_(thick|thin)\.txt$")
_NUMBER = re.compile(r"[-+]?\d+\.?\d*(?:[Ee][-+]?\d+)?")
_HEADER_LINES = 5
_VALUES_PER_MODEL = 16


@dataclass(frozen=True)
class Sequence:
    """One cooling sequence: arrays over its models, hottest first."""

    mass_sol: float
    envelope: str
    teff_k: np.ndarray
    log_g: np.ndarray
    radius_cm: np.ndarray
    age_yr: np.ndarray
    log_mass_fractions: dict

    def radius_at(self, teff_k):
        """The radius in cm at `teff_k`, clamped to the sequence's temperature range."""
        return float(np.interp(teff_k, self.teff_k[::-1], self.radius_cm[::-1]))


def parse_sequence(text, mass_sol, envelope):
    """Reads one sequence file (`text`): three lines per model: six, five and five numbers."""
    values = []
    for line in text.splitlines()[_HEADER_LINES:]:
        values.extend(float(token) for token in _NUMBER.findall(line))
    if not values or len(values) % _VALUES_PER_MODEL:
        raise ValueError(f"seq {mass_sol} {envelope}: {len(values)} numbers is not a whole number of models")
    models = np.asarray(values).reshape(-1, _VALUES_PER_MODEL)
    fractions = {name: models[:, column] for name, column in (("hydrogen", 12), ("helium", 13), ("carbon", 14), ("oxygen", 15))}
    return Sequence(mass_sol, envelope, models[:, 1], models[:, 2], models[:, 3], models[:, 4], fractions)


@functools.lru_cache(maxsize=4)
def _load(path):
    sequences = {}
    with tarfile.open(path) as archive:
        for member in archive:
            match = _NAME.search(member.name)
            if match and member.isfile():
                mass = int(match.group(1)) / 100.0
                text = archive.extractfile(member).read().decode("ascii")
                sequences[(mass, match.group(2))] = parse_sequence(text, mass, match.group(2))
    if not sequences:
        raise ValueError("the archive holds no sequence files")
    return sequences


def _sequences():
    return _load(fetch.locate(DATASET).require().path(ARCHIVE))


def masses_sol():
    """The masses of the grid in solar masses, ascending."""
    return tuple(sorted({mass for mass, _ in _sequences()}))


def sequence(mass_sol, envelope="thick"):
    """The sequence of exactly that grid mass (see `masses_sol`)."""
    return _sequences()[(round(mass_sol, 2), envelope)]


def _bracket(mass_sol):
    masses = masses_sol()
    if mass_sol < masses[0] or mass_sol > masses[-1]:
        return None
    upper = int(np.searchsorted(masses, mass_sol - 1e-12))
    if abs(masses[upper] - mass_sol) < 1e-9:
        return masses[upper], masses[upper], 0.0
    return masses[upper - 1], masses[upper], (mass_sol - masses[upper - 1]) / (masses[upper] - masses[upper - 1])


def radius_solar(mass_sol, teff_k, envelope="thick"):
    """The radius in solar radii of a carbon-oxygen white dwarf of `mass_sol` at `teff_k`, interpolated linearly in mass
    between the two grid masses around it, or None outside 0.20 to 1.30 solar masses."""
    bracket = _bracket(mass_sol)
    if bracket is None:
        return None
    low, high, weight = bracket
    radius_cm = (sequence(low, envelope).radius_at(teff_k) * (1 - weight) + sequence(high, envelope).radius_at(teff_k) * weight)
    return radius_cm / (constants.SOLAR_RADIUS_M * 100.0)


def layer_mass_fractions(mass_sol, envelope="thick"):
    """The total mass fractions `{hydrogen, helium, carbon, oxygen}` of the hottest model of the nearest grid mass: the
    hydrogen and helium layers are the case's (1e-4 or 1e-10, and 1e-2) and the core is split evenly between carbon and
    oxygen by the sequences' own assumption. None outside the grid."""
    bracket = _bracket(mass_sol)
    if bracket is None:
        return None
    low, high, weight = bracket
    nearest = high if weight > 0.5 else low
    first = sequence(nearest, envelope)
    return {name: float(10.0 ** values[0]) for name, values in first.log_mass_fractions.items()}
