# planetgen/physics/planet_occurrence.py

"""
Observed planet occurrence (GEN.728, GEN.828, GEN.829, GEN.729, TEST.395, docs/design/composition/planet-architecture.md)
========================================================================================================================

The numbers the planets of a generated sky are checked against: how many stars of each kind have planets, how many planets
those have, how rare hot Jupiters and cold giants are, how far apart neighbours sit, and how the Kepler survey saw the
multiplicity. `physics/data/planet_occurrence.json` holds them, each with the paper it comes from, the label of how far it was
checked (READ in the source, or 2ND-HAND in a paper that quotes it) and the registry tier.

A `Band` is one observed figure with its 1-sigma interval. `compare` puts a simulated figure against it: the distance in sigma
where sigma combines the published uncertainty and the sampling error of the simulation, so a figure from 100 systems is not
held to the precision of 30,000 stars. A published figure with no printed uncertainty is widened by the file's
`default_relative_width`, so it still has an interval.

The module draws nothing and writes nothing. Use:

    band = BANDS["inner_fraction.G"]
    result = compare(0.28, band, sampling_sigma=0.04)
    result.ok, result.z
    host_group(0.62)                                # "K"
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import detmath

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "planet_occurrence.json")
"""str: The data file."""

LABELS = ("READ", "2ND-HAND")
"""tuple: How far a figure was checked: read in the source, or read in a paper that quotes it (named in `citing`)."""

SIGMAS = 2.0
"""float: `compare` passes a figure within this many sigma."""


@dataclass(frozen=True)
class Band:
    """One observed figure: `value` with its 1-sigma `minus` and `plus` (from the paper, or the default width when it printed none)."""

    id: str
    host: str
    measure: str
    unit: str
    value: float
    minus: float
    plus: float
    printed: bool
    definition: str
    source: str
    label: str
    tier: str

    def interval(self, sigmas=1.0):
        """The `(low, high)` interval `sigmas` sigma either side of the value, floored at zero."""
        return max(0.0, self.value - sigmas * self.minus), self.value + sigmas * self.plus


@dataclass(frozen=True)
class Comparison:
    """A simulated figure against a band: `z` its distance in combined sigma (signed, positive above), `ok` within `SIGMAS`."""

    band: Band
    simulated: float
    z: float
    ok: bool


@dataclass(frozen=True)
class Multiplicity:
    """The Kepler transit multiplicity: stars with 1, 2, ... transiting planets (`counts`, key 1 upward), and the survey behind it."""

    stars_surveyed: int
    systems: int
    counts: dict
    source: str
    label: str

    def shares(self):
        """The share of the planet-bearing systems with each number of transiting planets."""
        return {k: v / self.systems for k, v in self.counts.items()}


def _finite(owner, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise DataError(f"{owner}: {value!r} is not a finite number")
    return float(value)


def load(path=DATA_FILE):
    """Reads and checks the file at `path`: `(hosts, bands, multiplicity)`, `bands` a dict by id."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the planet occurrence file {path}: {error}") from error
    width = _finite("default_relative_width", document.get("default_relative_width"))
    if not 0.0 < width < 1.0:
        raise DataError(f"default_relative_width {width} is not between 0 and 1")
    hosts = {}
    for name, entry in (document.get("hosts") or {}).items():
        low, high = _finite(f"host {name} low", entry.get("low")), _finite(f"host {name} high", entry.get("high"))
        if not 0.0 < low < high:
            raise DataError(f"host {name}: mass range {low} to {high} is not increasing")
        hosts[name] = (low, high)
    bands = {}
    for raw in document.get("bands") or []:
        band_id = raw.get("id")
        if not band_id or band_id in bands:
            raise DataError(f"band id {band_id!r} is missing or repeated")
        value = _finite(band_id, raw.get("value"))
        printed = raw.get("minus") is not None or raw.get("plus") is not None
        minus = _finite(band_id, raw["minus"]) if raw.get("minus") is not None else None
        plus = _finite(band_id, raw["plus"]) if raw.get("plus") is not None else None
        if minus is None:
            minus = plus if plus is not None else value * width
        if plus is None:
            plus = minus
        if raw.get("host") not in (*hosts, "all") or raw.get("label") not in LABELS or min(minus, plus) < 0.0:
            raise DataError(f"band {band_id}: host, label or uncertainty is not valid")
        if raw["label"] == "2ND-HAND" and not raw.get("citing"):
            raise DataError(f"band {band_id}: a 2ND-HAND figure names the paper that quotes it (citing)")
        for name in ("measure", "unit", "definition", "source", "tier"):
            if not raw.get(name):
                raise DataError(f"band {band_id}: {name} is missing")
        bands[band_id] = Band(band_id, raw["host"], raw["measure"], raw["unit"], value, minus, plus, printed, raw["definition"],
                              raw["source"], raw["label"], raw["tier"])
    raw = document.get("transit_multiplicity") or {}
    counts = {int(k): int(v) for k, v in (raw.get("counts_by_transiting_planets") or {}).items()}
    if not counts or sum(counts.values()) != raw.get("systems_with_a_transiting_planet"):
        raise DataError("the transit multiplicity counts do not add to the systems surveyed")
    multiplicity = Multiplicity(int(raw["stars_surveyed"]), int(raw["systems_with_a_transiting_planet"]), counts, raw["source"], raw["label"])
    return hosts, bands, multiplicity


HOSTS, BANDS, TRANSIT_MULTIPLICITY = load()
"""The host groups (name -> `(low, high)` solar masses), the bands by id and the transit multiplicity."""


def host_group(mass_solar):
    """The name of the host group (`HOSTS`) a star of `mass_solar` solar masses belongs to; `None` outside them."""
    for name, (low, high) in HOSTS.items():
        if low <= mass_solar < high:
            return name
    return None


def bands_of(measure, host=None):
    """The bands that state `measure`, for `host` (and the ones that hold for every host) when it is given."""
    return [b for b in BANDS.values() if b.measure == measure and (host is None or b.host in (host, "all"))]


def compare(simulated, band, sampling_sigma=0.0):
    """
    `simulated` against `band`: `z` is the difference over the root sum of squares of the band's uncertainty (the side facing the
    simulated figure) and `sampling_sigma`, the sampling error of the simulation.
    """
    published = band.minus if simulated < band.value else band.plus
    sigma = detmath.hypot(published, sampling_sigma)
    z = (simulated - band.value) / sigma if sigma > 0.0 else (0.0 if simulated == band.value else math.copysign(math.inf, simulated - band.value))
    return Comparison(band, simulated, z, abs(z) <= SIGMAS)


def share_sigma(share, n):
    """The sampling error of a share `share` of `n` systems (binomial)."""
    return math.sqrt(max(share * (1.0 - share), 1.0 / (4.0 * max(n, 1))) / max(n, 1))


def multiplicity_chi2(simulated_counts, multiplicity=TRANSIT_MULTIPLICITY):
    """
    The chi-square of a simulated transit multiplicity against the observed one, both normalised to the same number of
    planet-bearing systems: `(chi2, degrees_of_freedom)`. `simulated_counts` maps the number of transiting planets (1 upward) to
    systems; bins with fewer than 5 expected systems are pooled with the last one so no cell is empty.
    """
    total = sum(simulated_counts.values())
    if total <= 0:
        raise ValueError("no simulated system has a transiting planet")
    keys = sorted(multiplicity.counts)
    expected = {k: multiplicity.counts[k] / multiplicity.systems * total for k in keys}
    observed = {k: simulated_counts.get(k, 0) + 0.0 for k in keys}
    observed[keys[-1]] += sum(v for k, v in simulated_counts.items() if k > keys[-1])
    while len(keys) > 2 and expected[keys[-1]] < 5.0:
        last = keys.pop()
        expected[keys[-1]] += expected.pop(last)
        observed[keys[-1]] += observed.pop(last)
    chi2 = sum((observed[k] - expected[k]) ** 2 / expected[k] for k in keys)
    return chi2, len(keys) - 1
