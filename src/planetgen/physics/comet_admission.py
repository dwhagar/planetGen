# planetgen/physics/comet_admission.py

"""
The comet admission model (GEN.507, docs/design/composition/tunable-gaps.md sections 3.5 and 4, docs/design/composition/comet-admission.md)
===============================================================================================================

Which comets have a published coma composition and at what heliocentric distance it was measured, so that the comet volatile ratios
(GEN.495, `comet_population`) and their heliocentric-distance term are compared with the comets that were measured and not with every
comet drawn. The compositions of the database [R3782] are of comets near the Sun and bright.

Two things, both fitted by `scripts/fit_comet_admission.py` to the database joined to the JPL Small-Body Database orbits:

1. **Admission.** The share of known comets of each group (`jupiter_family`, `halley_type`, `long_period`; a parabolic comet is long
   period) and perihelion bin (0.3 to 1, 1 to 2, 2 to 8 AU) that have a measured composition: 93 of 2,422. `admission_probability` is
   that share (the Jeffreys estimate `(k + 1/2) / (n + 1)`, so an empty cell is not zero).
2. **Distance.** `ln Rh = intercept + slope ln q + scatter eps` for two kinds of species, parent molecules (`parent`) and radicals and
   dust (`radical`); a comet is measured at its heliocentric distance `Rh`, never inside its perihelion, so the draw is truncated at
   `Rh >= q` (one random number by the inverse CDF).

The item a target's sampler returns is a mapping with `period_class` (or None for a parabolic comet) and `perihelion_distance_au`.
`observed_sampler` adds `rh_au` (one draw) and the selection model `comet_admission` returns the admission probability.

What it conditions on: a known comet. A comet nobody found is not in the fit, and the admission shares are those of found comets. The
nucleus size is not an input (the orbit table has no nucleus size for most comets); the activity of a drawn comet is not an input.
Validity: perihelion distance 0.3 to 8 AU.

Use:

    comet = {"period_class": "jupiter_family", "perihelion_distance_au": 1.3}
    admission_probability(comet)                    # the chance it has a measured composition
    measurement_distance_au(comet, "parent")        # where it would be measured (needs a bound draw stream)
"""

import json
import math
import os
import statistics
from dataclasses import dataclass

from planetgen.generation import tuning_harness
from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "comet_admission.json")
"""str: The data file."""

MODEL_ID = "comet-admission"
"""str: The id of the selection model in `tuning_harness.SELECTION_MODEL_ITEMS`."""

GROUPS = ("jupiter_family", "halley_type", "long_period")
KINDS = ("parent", "radical")

_PHI = statistics.NormalDist().cdf
_PHI_INV = statistics.NormalDist().inv_cdf


@dataclass(frozen=True)
class Cell:
    """The comets of one group and perihelion bin: how many are known and how many have a measured composition."""

    group: str
    q_min_au: float
    q_max_au: float
    comets: int
    measured: int

    @property
    def probability(self):
        """The Jeffreys estimate of the share measured."""
        return (self.measured + 0.5) / (self.comets + 1.0)


@dataclass(frozen=True)
class Distance:
    """`ln Rh = intercept + slope ln q + scatter eps`, from `rows` measurements of `comets` comets."""

    intercept: float
    slope: float
    scatter: float
    rows: int
    comets: int


@dataclass(frozen=True)
class Model:
    """The fitted model."""

    cells: tuple
    distance: dict
    known_comets: int
    measured_comets: int
    source: str
    retrieved: str


def _finite(owner, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise DataError(f"{owner}: {value!r} is not a finite number")
    return float(value)


def load(path=DATA_FILE):
    """Reads and checks the fitted model at `path`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the comet admission model {path}: {error}") from error
    cells = []
    for raw in document.get("admission") or []:
        comets, measured = raw.get("comets"), raw.get("measured")
        if raw.get("group") not in GROUPS or not isinstance(comets, int) or not isinstance(measured, int) or not 0 <= measured <= comets:
            raise DataError(f"admission cell {raw!r} is not a group and a count")
        low, high = _finite("q_min_au", raw.get("q_min_au")), _finite("q_max_au", raw.get("q_max_au"))
        if not 0 < low < high:
            raise DataError(f"admission cell {raw!r} has no perihelion range")
        cells.append(Cell(raw["group"], low, high, comets, measured))
    for group in GROUPS:
        ranges = sorted((c.q_min_au, c.q_max_au) for c in cells if c.group == group)
        if not ranges or any(a[1] != b[0] for a, b in zip(ranges, ranges[1:])):
            raise DataError(f"the perihelion bins of {group} are missing or not contiguous")
    distance = {}
    for kind in KINDS:
        raw = (document.get("distance") or {}).get(kind)
        if raw is None:
            raise DataError(f"the distance model has no {kind}")
        scatter = _finite(f"{kind} scatter", raw.get("scatter"))
        if scatter <= 0:
            raise DataError(f"{kind}: the scatter {scatter} is not positive")
        distance[kind] = Distance(_finite(f"{kind} intercept", raw.get("intercept")), _finite(f"{kind} slope", raw.get("slope")), scatter,
                                  int(raw.get("rows", 0)), int(raw.get("comets", 0)))
    return Model(cells=tuple(cells), distance=distance, known_comets=int(document.get("known_comets", 0)),
                 measured_comets=int(document.get("measured_comets", 0)), source=document.get("source", ""),
                 retrieved=document.get("retrieved", ""))


MODEL = load()
"""Model: The fitted model."""


def validity():
    """The range of the fit."""
    low = min(c.q_min_au for c in MODEL.cells)
    high = max(c.q_max_au for c in MODEL.cells)
    return {"perihelion_distance_au": (low, high)}


def group_of(item):
    """The group of a comet item: its `period_class`, or `long_period` for a parabolic comet (`None`)."""
    group = item.get("period_class")
    if group is None:
        return "long_period"
    if group not in GROUPS:
        raise DataError(f"period class {group!r} is none of {', '.join(GROUPS)}")
    return group


def cell_of(item):
    """The `Cell` of a comet item; a perihelion outside the fit's range takes the nearest bin."""
    group, q = group_of(item), item["perihelion_distance_au"]
    cells = sorted((c for c in MODEL.cells if c.group == group), key=lambda c: c.q_min_au)
    for cell in cells:
        if q < cell.q_max_au:
            return cell
    return cells[-1]


def admission_probability(item):
    """The chance a comet `item` has a published composition."""
    return cell_of(item).probability


def measurement_distance_au(item, kind="parent"):
    """
    The heliocentric distance a comet would be measured at, in AU: one number from the bound draw stream, by the inverse CDF of a
    normal truncated so that the distance is not under the perihelion distance.
    """
    try:
        line = MODEL.distance[kind]
    except KeyError:
        raise DataError(f"species kind {kind!r} is none of {', '.join(KINDS)}") from None
    q = item["perihelion_distance_au"]
    mean = line.intercept + line.slope * math.log(q)
    low = _PHI((math.log(q) - mean) / line.scatter)
    u = min(max(low + draw.random() * (1.0 - low), 1e-12), 1.0 - 1e-12)
    return math.exp(mean + line.scatter * _PHI_INV(u))


def observed_sampler(sampler, kind="parent"):
    """A sampler for `tuning_harness.simulate`: each item of `sampler(i)` gets `rh_au`, the distance it would be measured at."""
    def sample(i):
        item = sampler(i)
        return None if item is None else dict(item, rh_au=measurement_distance_au(item, kind))
    return sample


def comet_admission(item):
    """The selection model: the probability the comet has a published composition."""
    return admission_probability(item)


def expected_measured(items):
    """The number of `items` a database would hold: the sum of their admission probabilities."""
    return sum(admission_probability(item) for item in items)


tuning_harness.register_selection_model(MODEL_ID, comet_admission)
