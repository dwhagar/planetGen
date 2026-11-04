# planetgen/datasets/white_dwarf_models.py

"""
White-dwarf cooling sequences with their layers (GEN.548).

Source: [R0308],
23 masses from 0.20 to 1.30 solar masses for a thick (q_He = 1e-2, q_H = 1e-4) and a thin (q_He = 1e-2, q_H = 1e-10)
hydrogen layer, on Antoine Bédard's page at the University of Warwick (`sequences.tar`, 2.5 MB, 46 text files). Each
model lists the effective temperature, log g, radius, cooling age, photon luminosity, central conditions, the
crystallised fraction, the neutrino luminosity and the log of the total mass fractions of hydrogen, helium, carbon and
oxygen. The core of these models is homogeneous and equimassic carbon and oxygen: the C:O split of 1:1 is the
sequences' own assumption, not a measurement, and the readers say so. The sequences are carbon-oxygen white dwarfs;
a helium-core or an oxygen-neon-core white dwarf is not in them.

The readers raise `fetch.DataMissing` when the file has not been downloaded, and answer `None` outside the masses
(0.20 to 1.30) of the grid. Temperatures outside a sequence's range are clamped to its ends.

GEN.742: `cooling_state` answers the other way round, for a white dwarf of a given mass and cooling age: its effective
temperature, photon luminosity, radius and the crystallised share of its core, interpolated in cooling age on the two
grid masses around it and then in mass (logarithms of the temperature and luminosity, so that a decade of cooling is a
straight step). The cooling age of the sequences is counted from a temperature near 100,000 K; the few tens of thousands of
years from the end of the asymptotic branch to there are not counted. A white dwarf older than a sequence's last model
(the sequences stop near 1,500 K and a few millionths of a solar luminosity) is held at that last model and the answer says so
(`extrapolated`): below that the models do not go, and the star is a black dwarf for every use the generator has. Above 1.30 solar masses (the oxygen-neon white dwarfs, not in the sequences) the 1.30 sequence stands in and the
answer says so.
"""

import functools
import re
import tarfile
from dataclasses import dataclass

import math

import numpy as np

from planetgen.datasets import fetch
from planetgen.physics import constants

ARCHIVE = "sequences.tar"

DATASET = fetch.Dataset(
    key="white_dwarf_models",
    title="[R0308] white-dwarf cooling sequences",
    source="[R0308]",
    pieces=(fetch.Piece(ARCHIVE, "https://warwick.ac.uk/fac/sci/physics/research/astro/people/antoinebedard/models/"
                        + ARCHIVE, "5c99114fcba0087395fcde2779ad7bdc9b7ace224abb71322d45165facae45ac", 2519040),),
    replaces="the Chandrasekhar-type law of the stellar-composition note, an even carbon-oxygen split and a fixed layer split",
    fill="SC2",
    fallback="the radius of the zero-temperature Chandrasekhar table (Sirius B check, no temperature or layers), no layer fractions, "
             "and the Mestel power law of the cooling (the stand-in that came before the sequences: luminosity of the cooling age, one radius per mass)",
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
    luminosity_erg_s: np.ndarray = None
    crystallised_fraction: np.ndarray = None

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
    return Sequence(mass_sol, envelope, models[:, 1], models[:, 2], models[:, 3], models[:, 4], fractions,
                    models[:, 5], models[:, 9])


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


CRYSTAL_MASS_SOL = 1.30
"""float: The heaviest sequence; a heavier white dwarf is read from it (an oxygen-neon core, which the sequences do not hold)."""


@dataclass(frozen=True)
class CoolingState:
    """A white dwarf at one cooling age: `teff_k`, `luminosity_sol` (photon luminosity), `radius_solar`, `crystallised_fraction` of the
    core's mass (0 to 1), `extrapolated` (older than the sequences' last model) and `mass_clamped` (heavier than the grid)."""

    teff_k: float
    luminosity_sol: float
    radius_solar: float
    crystallised_fraction: float
    extrapolated: bool
    mass_clamped: bool


def _at_age(sequence, age_yr):
    """`(log teff, log luminosity in erg/s, radius in cm, crystallised share, extrapolated)` of one sequence at `age_yr`."""
    ages = sequence.age_yr
    log_teff = np.log(sequence.teff_k)
    log_lum = np.log(sequence.luminosity_erg_s)
    if age_yr <= ages[-1]:
        return (float(np.interp(age_yr, ages, log_teff)), float(np.interp(age_yr, ages, log_lum)),
                float(np.interp(age_yr, ages, sequence.radius_cm)), float(np.interp(age_yr, ages, sequence.crystallised_fraction)),
                False)
    # Past the last model the sequence has nothing to say: the white dwarf is held at the last model (under 1,500 K and about a
    # millionth of a solar luminosity, a black dwarf for every use the generator has) and the answer says it was held.
    return (float(log_teff[-1]), float(log_lum[-1]), float(sequence.radius_cm[-1]), float(sequence.crystallised_fraction[-1]), True)


def cooling_state(mass_sol, cooling_age_gy, envelope="thick"):
    """
    The state of a carbon-oxygen white dwarf of `mass_sol` solar masses `cooling_age_gy` after its cooling age began, from the sequences.

    Raises:
        fetch.DataMissing: The sequences were not downloaded.
        ValueError: A mass or age that is not positive.
    """
    if not (mass_sol > 0.0 and cooling_age_gy > 0.0):
        raise ValueError("a white dwarf has a positive mass and cooling age")
    clamped = mass_sol > CRYSTAL_MASS_SOL
    bracket = _bracket(min(max(mass_sol, masses_sol()[0]), CRYSTAL_MASS_SOL))
    low, high, weight = bracket
    age_yr = max(cooling_age_gy * 1.0e9, sequence(low, envelope).age_yr[1])
    first = _at_age(sequence(low, envelope), age_yr)
    second = first if weight == 0.0 else _at_age(sequence(high, envelope), age_yr)
    mixed = [first[i] * (1 - weight) + second[i] * weight for i in range(4)]
    solar_erg_s = constants.SOLAR_LUMINOSITY * 1.0e7
    return CoolingState(teff_k=math.exp(mixed[0]), luminosity_sol=math.exp(mixed[1]) / solar_erg_s,
                        radius_solar=mixed[2] / (constants.SOLAR_RADIUS_M * 100.0), crystallised_fraction=min(max(mixed[3], 0.0), 1.0),
                        extrapolated=bool(first[4] or second[4]), mass_clamped=clamped)
