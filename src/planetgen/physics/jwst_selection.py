# planetgen/physics/jwst_selection.py

"""
The JWST atmosphere target-selection model (GEN.506, docs/design/composition/tunable-gaps.md sections 3.3, 3.4 and 4, docs/design/composition/jwst-selection.md)
==========================================================================================================================================

Which planets JWST has observed for their atmospheres, so that a tuned share of clear, retained or hazy atmospheres is compared with the
observed list and not with the population. JWST chooses for signal: a nearby, hot, large-scale-height planet around a bright host.

The model is a logistic regression of the label "has a published JWST transmission or emission spectrum" on the transmission
spectroscopy metric `TSM` ([R2315], equation 1, scale factors of its Table 1), the emission metric `ESM` (equation 4) and an
offset for each radius class (under 1.5, 1.5 to 4, 4 to 10 and over 10 Earth radii), fitted to the NASA Exoplanet Archive's 4,618
transiting planets with the needed quantities, 78 of which have such a spectrum (`scripts/fit_jwst_selection.py`). The offsets make the
model reproduce the number observed in each class, which is how the programs were chosen (rocky planets, sub-Neptunes, giants).

The metrics need the host's apparent J and K magnitudes. A generated system has none, so `with_magnitudes` draws them from what the
archive's hosts show: J a truncated normal, J-K a line in log10 of the host's effective temperature with a scatter. A target's sampler
returns a mapping with `radius_earth`, `mass_earth`, `teq_k`, `host_radius_solar` and `host_teff_k`; `observed_sampler` adds the
magnitudes (three draws) and the selection model `jwst_target` returns the probability the planet is in the observed list.

What it conditions on: a known transiting planet with a mass. Whether it was found and whether a JWST program is its own accident are in
the fit (the magnitudes of the hosts are those of planets that were found), not separated. Validity: transiting planets, J between the
archive's extremes, the four radius classes.

Use:

    item = observed_sampler(sampler)(i)                  # adds j_mag and k_mag
    probability(item)                                    # the chance it has a JWST spectrum
    tsm(item), esm(item)                                 # the two metrics
"""

import json
import math
import os
import statistics
from dataclasses import dataclass

from planetgen.generation import tuning_harness
from planetgen.physics import constants as physical_constants
from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "jwst_selection.json")
"""str: The data file."""

MODEL_ID = "jwst-target-selection"
"""str: The id of the selection model in `tuning_harness.SELECTION_MODEL_ITEMS`."""

TSM_SCALE = ((1.5, 0.190), (2.75, 1.26), (4.0, 1.28), (float("inf"), 1.15))
"""tuple: `(upper radius in Earth radii, scale factor)` of [R2315: Table 1]."""

ESM_PREFACTOR = 4.29e6
DAYSIDE_OVER_EQUILIBRIUM = 1.10
RADIUS_RATIO = 6371.0 / (physical_constants.SOLAR_RADIUS_M / 1000)
"""float: Earth radii per solar radius."""

_PHI_INV = statistics.NormalDist().inv_cdf
_PHI = statistics.NormalDist().cdf


@dataclass(frozen=True)
class Model:
    """
    The fitted model.

    Attributes:
        planets, observed (int): The fitted sample and those with a spectrum.
        intercept, tsm_slope, esm_slope (float): The logistic coefficients (per dex of the metric).
        class_edges (tuple): The upper radius of each class in Earth radii (the last is infinite).
        class_offsets (tuple): The offset of each class (the first is 0).
        mean_j, sd_j, min_j, max_j (float): The host J magnitude of the sample.
        colour_intercept, colour_slope, colour_scatter (float): J-K = intercept + slope log10(Teff) + scatter eps.
        calibration (tuple): Per class `{planets, observed, predicted}` of the fit.
        auc (float): The area under the ROC curve of the fit.
        source, retrieved (str): Where the sample came from and when.
    """

    planets: int
    observed: int
    intercept: float
    tsm_slope: float
    esm_slope: float
    class_edges: tuple
    class_offsets: tuple
    mean_j: float
    sd_j: float
    min_j: float
    max_j: float
    colour_intercept: float
    colour_slope: float
    colour_scatter: float
    calibration: tuple
    auc: float
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
        raise DataError(f"cannot read the JWST selection model {path}: {error}") from error
    coefficients = document.get("coefficients") or {}
    classes = document.get("radius_classes_earth") or []
    offsets = coefficients.get("class_offsets") or []
    if len(classes) < 2 or len(offsets) != len(classes):
        raise DataError("the radius classes and their offsets do not match")
    edges = []
    previous = 0.0
    for low, high in classes:
        if _finite("class edge", low) != previous or (high is not None and _finite("class edge", high) <= previous):
            raise DataError("the radius classes are not contiguous and increasing")
        edges.append(math.inf if high is None else float(high))
        previous = edges[-1]
    magnitudes = document.get("host_magnitudes") or {}
    sd_j = _finite("sd_j", magnitudes.get("sd_j"))
    scatter = _finite("colour_scatter", magnitudes.get("colour_scatter"))
    min_j, max_j = _finite("min_j", magnitudes.get("min_j")), _finite("max_j", magnitudes.get("max_j"))
    if sd_j <= 0 or scatter < 0 or not min_j < max_j:
        raise DataError("the host magnitude model is out of range")
    planets, observed = document.get("planets"), document.get("observed")
    if not isinstance(planets, int) or not isinstance(observed, int) or not 0 < observed < planets:
        raise DataError(f"the counts {observed!r} of {planets!r} are not a sample")
    calibration = tuple({k: v for k, v in c.items() if k in ("planets", "observed", "predicted")} for c in document.get("calibration") or [])
    if len(calibration) != len(classes):
        raise DataError("the calibration does not cover every radius class")
    return Model(planets=planets, observed=observed, intercept=_finite("intercept", coefficients.get("intercept")),
                 tsm_slope=_finite("log10_tsm", coefficients.get("log10_tsm")), esm_slope=_finite("log10_esm", coefficients.get("log10_esm")),
                 class_edges=tuple(edges), class_offsets=tuple(_finite("offset", v) for v in offsets),
                 mean_j=_finite("mean_j", magnitudes.get("mean_j")), sd_j=sd_j, min_j=min_j, max_j=max_j,
                 colour_intercept=_finite("colour_intercept", magnitudes.get("colour_intercept")),
                 colour_slope=_finite("colour_slope", magnitudes.get("colour_slope_per_log10_teff")), colour_scatter=scatter,
                 calibration=calibration, auc=_finite("auc", document.get("auc_logistic")), source=document.get("source", ""),
                 retrieved=document.get("retrieved", ""))


MODEL = load()
"""Model: The fitted model."""


def _planck_7_5(temperature):
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    return 1.0 / math.expm1(h * c / (7.5e-6 * k * temperature))


def radius_class(radius_earth):
    """The index of the radius class of a planet."""
    for index, edge in enumerate(MODEL.class_edges):
        if radius_earth < edge:
            return index
    return len(MODEL.class_edges) - 1


def tsm(item):
    """The transmission spectroscopy metric of [R2315] (a mapping with `radius_earth`, `mass_earth`, `teq_k`,
    `host_radius_solar` and `j_mag`)."""
    scale = next(s for limit, s in TSM_SCALE if item["radius_earth"] < limit)
    return (scale * item["radius_earth"] ** 3 * item["teq_k"] / (item["mass_earth"] * item["host_radius_solar"] ** 2)
            * 10.0 ** (-item["j_mag"] / 5.0))


def esm(item):
    """The emission spectroscopy metric of [R2315] (`radius_earth`, `teq_k`, `host_radius_solar`, `host_teff_k`, `k_mag`)."""
    return (ESM_PREFACTOR * _planck_7_5(DAYSIDE_OVER_EQUILIBRIUM * item["teq_k"]) / _planck_7_5(item["host_teff_k"])
            * (item["radius_earth"] / item["host_radius_solar"] * RADIUS_RATIO) ** 2 * 10.0 ** (-item["k_mag"] / 5.0))


def probability(item):
    """The chance a planet `item` (with `j_mag` and `k_mag`) has a published JWST spectrum."""
    eta = (MODEL.intercept + MODEL.tsm_slope * math.log10(tsm(item)) + MODEL.esm_slope * math.log10(esm(item))
           + MODEL.class_offsets[radius_class(item["radius_earth"])])
    return 1.0 / (1.0 + math.exp(-eta))


def with_magnitudes(item):
    """
    The item with a host J and K magnitude from what the archive's hosts show: three numbers from the bound draw stream (J by the
    inverse CDF of a normal truncated to the archive's range, then a normal for the J-K scatter).
    """
    low = _PHI((MODEL.min_j - MODEL.mean_j) / MODEL.sd_j)
    high = _PHI((MODEL.max_j - MODEL.mean_j) / MODEL.sd_j)
    u = min(max(low + draw.random() * (high - low), 1e-12), 1.0 - 1e-12)
    j_mag = MODEL.mean_j + MODEL.sd_j * _PHI_INV(u)
    colour = (MODEL.colour_intercept + MODEL.colour_slope * math.log10(item["host_teff_k"])
              + MODEL.colour_scatter * draw.gauss(0.0, 1.0))
    return dict(item, j_mag=j_mag, k_mag=j_mag - colour)


def observed_sampler(sampler):
    """A sampler for `tuning_harness.simulate`: each item of `sampler(i)` gets host magnitudes (`None` stays `None`)."""
    def sample(i):
        item = sampler(i)
        return None if item is None else with_magnitudes(item)
    return sample


def jwst_target(item):
    """The selection model: the probability the planet is in the list of planets with a JWST spectrum."""
    return probability(item)


def expected_observed(items):
    """The number of `items` (each with magnitudes) a JWST list would hold: the sum of their probabilities."""
    return sum(probability(item) for item in items)


tuning_harness.register_selection_model(MODEL_ID, jwst_target)
