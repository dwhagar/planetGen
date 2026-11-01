# planetgen/physics/wd_pollution.py

"""
The white-dwarf pollution detection model (GEN.508, docs/design/composition/tunable-gaps.md sections 3.2 and 4, docs/design/composition/wd-pollution-detection.md)
===========================================================================================================================

What the 24 polluted white dwarfs with an oxygen fugacity ([R1152: Table 1], [R1153: Table 4]) can and cannot show.
The measured quantity is dIW = 2 log10(x_FeO / x_Fe,metal) of the parent body, with x_Fe,metal = 0.85 and x_FeO the mole fraction of
FeO in the rock, found by Monte Carlo draws of the measured abundances. Those papers' own rule decides what is reported:

* 50 percent or more of the draws give a negative x_FeO: no median, an upper limit only (`limit`);
* 16.5 to 50 percent: a median and an upper bound, no lower bound (`upper`);
* under 16.5 percent: a median with both bounds (`both`).

So a body with little FeO relative to the error on the abundances is reported with a large error, or not at all, and the reported
medians of the reduced end are scattered by the measurement and not only by the bodies. A tuned width is compared with that sample
and not with the bodies. The model:

1. the abundance error on x_FeO is absolute: `x_meas = x_true + sigma z0` with `z0` standard normal, and `sigma` is log-normal across
   stars. Its two parameters are fitted by maximum likelihood to the sample: the 16 stars with both bounds give their `sigma` directly
   (half the spread of x_FeO between the 1-sigma bounds), the 6 stars with a median and an upper bound give `sigma >= x/z` (the
   16.5 percent rule), and the 2 stars with a limit only give no more than the largest of those (a convention, stated in the doc);
2. the reported value is `2 log10(x_meas / 0.85)`, with the class the rule above gives from `x_meas / sigma`;
3. the model's selection (`reported_median`) keeps a star whose class is `both` or `upper` and drops a `limit`.

Validity: effective temperatures 10,120 to 21,200 K, hydrogen and helium atmospheres (DA, DB, DZ), parent bodies with oxygen measured.
Not modelled: the settling of the elements out of the atmosphere ([R1152] section 3.2, Figs. 8 to 11), which can make a body in
a declining phase look more reduced; a body observed in its build-up or steady-state phase is the model's case. A tuned width is
an upper bound where the sample may hold stars observed after accretion has waned.

Use:

    stars = SAMPLE.stars                                  # the 24 stars
    item = observe(-1.25)                                 # one body measured: an Observation (needs a bound draw stream)
    expected_sample(-1.25, 0.2)                           # the sample a population N(-1.25, 0.2) would show, no random numbers
    tuning_harness.SELECTION_MODELS["wd-pollution-detection"]
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

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "wd_pollution_sample.json")
"""str: The data file."""

MODEL_ID = "wd-pollution-detection"
"""str: The id of the selection model in `tuning_harness.SELECTION_MODEL_ITEMS`."""

RECOVERIES = ("both", "upper", "limit")
"""tuple: What a star yields, from the most to the least: both bounds, a median and an upper bound, an upper limit only."""

_PHI_INV = statistics.NormalDist().inv_cdf


def _phi(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


@dataclass(frozen=True)
class Star:
    """
    One polluted white dwarf.

    Attributes:
        name (str): As the paper has it.
        type (str): The spectral type (DAZ, DBAZ, DZBA, ...).
        teff_k (float): The effective temperature.
        median_diw (float or None): The median dIW; None for an upper limit only.
        plus (float): The upper 1-sigma bound above the median (the limit itself for a star with no median).
        minus (float or None): The lower bound below the median; None when it cannot be recovered.
        source (str): Key into `Sample.sources`.
        table (str): The table of the paper.
    """

    name: str
    type: str
    teff_k: float
    median_diw: float | None
    plus: float
    minus: float | None
    source: str
    table: str

    @property
    def recovery(self):
        """`both`, `upper` or `limit` (`RECOVERIES`)."""
        if self.median_diw is None:
            return "limit"
        return "both" if self.minus is not None else "upper"


@dataclass(frozen=True)
class Sample:
    """
    The sample.

    Attributes:
        stars (tuple): The `Star`s in the papers' order.
        x_metal (float): The mole fraction of Fe in the metal.
        negative_upper, negative_limit (float): The negative-draw fractions from which a star yields an upper bound only and a limit only.
        sources (dict): Key -> `{"citation", "label"}`.
    """

    stars: tuple
    x_metal: float
    negative_upper: float
    negative_limit: float
    sources: dict

    @property
    def values(self):
        """The reported medians."""
        return tuple(s.median_diw for s in self.stars if s.median_diw is not None)

    @property
    def counts(self):
        """`{recovery: stars}`."""
        return {r: sum(1 for s in self.stars if s.recovery == r) for r in RECOVERIES}

    @property
    def teff_range(self):
        temps = [s.teff_k for s in self.stars]
        return min(temps), max(temps)


def _finite(owner, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise DataError(f"{owner}: {value!r} is not a finite number")
    return float(value)


def load(path=DATA_FILE):
    """Reads and checks the sample at `path`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the white-dwarf pollution sample {path}: {error}") from error
    sources = document.get("sources") or {}
    for key, source in sources.items():
        if not source.get("citation") or not source.get("label"):
            raise DataError(f"source {key} has no citation or label")
    x_metal = _finite("x_metal", document.get("x_metal"))
    negative_upper = _finite("negative_fraction_median_and_upper", document.get("negative_fraction_median_and_upper"))
    negative_limit = _finite("negative_fraction_limit_only", document.get("negative_fraction_limit_only"))
    if not 0 < x_metal <= 1 or not 0 < negative_upper < negative_limit < 1:
        raise DataError("x_metal or the negative-draw fractions are out of range")
    stars, seen = [], set()
    for raw in document.get("stars") or []:
        name = raw.get("name")
        if not name or name in seen:
            raise DataError(f"star {name!r} is missing or listed twice")
        seen.add(name)
        if raw.get("source") not in sources:
            raise DataError(f"{name}: source {raw.get('source')!r} is not in the sources")
        median = raw.get("median_diw")
        minus = raw.get("minus")
        plus = _finite(f"{name} plus", raw.get("plus"))
        if plus <= 0:
            raise DataError(f"{name}: the upper bound {plus} is not positive")
        if median is None and minus is not None:
            raise DataError(f"{name}: an upper limit has no lower bound")
        if minus is not None and _finite(f"{name} minus", minus) <= 0:
            raise DataError(f"{name}: the lower bound {minus} is not positive")
        stars.append(Star(name=name, type=raw.get("type", ""), teff_k=_finite(f"{name} teff_k", raw.get("teff_k")),
                          median_diw=None if median is None else _finite(f"{name} median", median), plus=plus,
                          minus=None if minus is None else float(minus), source=raw["source"], table=str(raw.get("table", ""))))
    if not stars:
        raise DataError("the sample has no stars")
    return Sample(stars=tuple(stars), x_metal=x_metal, negative_upper=negative_upper, negative_limit=negative_limit, sources=sources)


SAMPLE = load()
"""Sample: The 24 polluted white dwarfs of [R1152] and 2023."""


def x_of(diw, x_metal=None):
    """The FeO mole fraction of the rock for a dIW."""
    return (SAMPLE.x_metal if x_metal is None else x_metal) * 10.0 ** (diw / 2.0)


def diw_of(x, x_metal=None):
    """The dIW of a positive FeO mole fraction."""
    return 2.0 * math.log10(x / (SAMPLE.x_metal if x_metal is None else x_metal))


# --- the abundance error ------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ErrorModel:
    """
    The log-normal distribution of the absolute error `sigma` on x_FeO across stars.

    Attributes:
        mu, sigma (float): The mean and the standard deviation of ln(sigma_x).
        exact (int): Stars whose error is read from both bounds.
        censored (int): Stars known only to have an error at least some size.
        flagged (bool): Whether the censored stars are used (False: the optimistic bracket from the 16 stars with both bounds alone).
    """

    mu: float
    sigma: float
    exact: int
    censored: int
    flagged: bool


def _censor_scale(sample):
    return -_PHI_INV(sample.negative_upper)


def error_data(sample=SAMPLE):
    """`(exact, censored)`: the sigma_x read from both bounds, and the lower limits on sigma_x of the stars without a lower bound.

    A star with a median and an upper bound has at least 16.5 percent negative draws, so `sigma_x >= x / z` with `z` the normal
    quantile of that fraction. A star with a limit only has no median to give a bound: it takes the largest bound of the others
    (a convention: the data say only that its error is large).
    """
    z = _censor_scale(sample)
    exact, censored = [], []
    for star in sample.stars:
        if star.recovery == "both":
            exact.append((x_of(star.median_diw + star.plus, sample.x_metal) - x_of(star.median_diw - star.minus, sample.x_metal)) / 2.0)
        elif star.recovery == "upper":
            censored.append(x_of(star.median_diw, sample.x_metal) / z)
    if censored:
        censored += [max(censored)] * sum(1 for s in sample.stars if s.recovery == "limit")
    return exact, censored


def _log_likelihood(mu, sigma, exact, censored):
    total = 0.0
    for value in exact:
        z = (math.log(value) - mu) / sigma
        total += -math.log(value * sigma * math.sqrt(2.0 * math.pi)) - 0.5 * z * z
    for bound in censored:
        total += math.log(max(0.5 * math.erfc((math.log(bound) - mu) / sigma / math.sqrt(2.0)), 1e-300))
    return total


def fit_error_model(sample=SAMPLE, flagged=True):
    """
    The maximum-likelihood `ErrorModel` of the sample, by a deterministic grid that is refined three times.

    Args:
        sample (Sample): The stars.
        flagged (bool): Use the stars without a lower bound as lower limits on the error. False gives the stars with both bounds
            alone: the error is then underestimated, because those are the stars whose error was small (a bracket, not a fit).
    """
    exact, censored = error_data(sample)
    if not flagged:
        censored = []
    if len(exact) < 3:
        raise DataError("the error model needs at least 3 stars with both bounds")
    lo_mu, hi_mu, lo_s, hi_s = -6.0, 0.0, 0.05, 3.0
    best = None
    for _ in range(4):
        grid = [(lo_mu + (hi_mu - lo_mu) * i / 20.0, lo_s + (hi_s - lo_s) * j / 20.0) for i in range(21) for j in range(21)]
        score, mu, sigma = max((_log_likelihood(m, s, exact, censored), m, s) for m, s in grid)
        best = (mu, sigma)
        step_mu, step_s = (hi_mu - lo_mu) / 10.0, (hi_s - lo_s) / 10.0
        lo_mu, hi_mu = mu - step_mu, mu + step_mu
        lo_s, hi_s = max(sigma - step_s, 0.02), sigma + step_s
    return ErrorModel(mu=best[0], sigma=best[1], exact=len(exact), censored=len(censored), flagged=flagged)


ERROR_MODEL = fit_error_model()
"""ErrorModel: The fit to the whole sample (the default)."""

BOUNDED_ERROR_MODEL = fit_error_model(flagged=False)
"""ErrorModel: The fit to the 16 stars with both bounds alone; the optimistic end of the bracket."""

VALIDITY = {"teff_k": SAMPLE.teff_range}
"""dict: The range of effective temperature of the sample."""


# --- observing a body ---------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Observation:
    """
    A body as the papers' method would report it.

    Attributes:
        diw_true (float): The body's dIW.
        diw (float or None): The reported median; None when only an upper limit is reported.
        recovery (str): One of `RECOVERIES`.
        negative_fraction (float): The fraction of the Monte Carlo draws with a negative FeO mole fraction.
    """

    diw_true: float
    diw: float | None
    recovery: str
    negative_fraction: float

    def as_item(self):
        """The mapping a selection model and a target read (`diw` is the target's field)."""
        return {"diw_true": self.diw_true, "diw": self.diw, "recovery": self.recovery, "negative_fraction": self.negative_fraction}


def classify(x_meas, sigma_x, sample=SAMPLE):
    """`(recovery, negative fraction)` of a measurement `x_meas` with error `sigma_x`."""
    negative = _phi(-x_meas / sigma_x)
    if negative >= sample.negative_limit:
        return "limit", negative
    return ("upper" if negative >= sample.negative_upper else "both"), negative


def observe(diw_true, model=None):
    """
    One body measured: two numbers from the bound draw stream, the error of the star and the realisation of the abundances.

    Args:
        diw_true (float): The body's dIW.
        model (ErrorModel, optional): `ERROR_MODEL` by default.
    """
    model = model or ERROR_MODEL
    sigma_x = math.exp(model.mu + model.sigma * draw.gauss(0.0, 1.0))
    x_meas = x_of(diw_true) + sigma_x * draw.gauss(0.0, 1.0)
    recovery, negative = classify(x_meas, sigma_x)
    diw = None if recovery == "limit" else diw_of(x_meas)
    return Observation(diw_true=float(diw_true), diw=diw, recovery=recovery, negative_fraction=negative)


def observed_sampler(sampler, model=None):
    """A sampler for `tuning_harness.simulate`: each draw of `sampler(i)` (a dIW) becomes `observe(...).as_item()`; give the target
    `field: "diw"` and `selection_model: "wd-pollution-detection"`."""
    return lambda i: observe(sampler(i), model).as_item()


def reported_median(item):
    """The selection model: 1 for a body whose median is reported (`both` or `upper`), 0 for an upper limit only."""
    recovery = item["recovery"]
    if recovery not in RECOVERIES:
        raise DataError(f"recovery {recovery!r} is none of {', '.join(RECOVERIES)}")
    return 0.0 if recovery == "limit" else 1.0


tuning_harness.register_selection_model(MODEL_ID, reported_median)


# --- the sample a population would show, with no random numbers -----------------------------------------------------------


@dataclass(frozen=True)
class ExpectedSample:
    """
    What a population would show.

    Attributes:
        fractions (dict): The share of each `RECOVERIES`.
        reported_mean, reported_sd (float): The mean and the standard deviation of the reported medians.
    """

    fractions: dict
    reported_mean: float
    reported_sd: float


def _grid(mean, sd, count):
    """`count` equal-probability points of a normal (the midpoint of each quantile slice)."""
    if sd == 0:
        return [mean]
    return [mean + sd * _PHI_INV((k + 0.5) / count) for k in range(count)]


def class_probabilities(diw_true, model=None, nodes=200):
    """`{recovery: probability}` of a body with dIW `diw_true`, by quadrature over the error of the star."""
    model = model or ERROR_MODEL
    x = x_of(diw_true)
    z_upper = -_PHI_INV(SAMPLE.negative_upper)
    both = upper_or_both = 0.0
    for sigma_x in (math.exp(model.mu + model.sigma * z) for z in _grid(0.0, 1.0, nodes)):
        both += _phi(x / sigma_x - z_upper)
        upper_or_both += _phi(x / sigma_x)
    both /= nodes
    upper_or_both /= nodes
    return {"both": both, "upper": upper_or_both - both, "limit": 1.0 - upper_or_both}


def expected_sample(mean, sd, model=None, bodies=25, errors=40, realisations=120):
    """
    The sample a population of bodies with dIW normal(mean, sd) shows: the class shares, and the mean and the standard
    deviation of the reported medians. Quadrature over the body, the error of the star and the abundance realisation, so the
    result has no random numbers and no seed.
    """
    model = model or ERROR_MODEL
    counts = {r: 0.0 for r in RECOVERIES}
    total = kept = first = second = 0.0
    sigmas = [math.exp(model.mu + model.sigma * z) for z in _grid(0.0, 1.0, errors)]
    z0s = _grid(0.0, 1.0, realisations)
    for diw_true in _grid(mean, sd, bodies):
        x = x_of(diw_true)
        for sigma_x in sigmas:
            for z0 in z0s:
                x_meas = x + sigma_x * z0
                recovery, _negative = classify(x_meas, sigma_x)
                counts[recovery] += 1
                total += 1
                if recovery != "limit":
                    value = diw_of(x_meas)
                    kept += 1
                    first += value
                    second += value * value
    mean_reported = first / kept
    variance = max(second / kept - mean_reported ** 2, 0.0)
    return ExpectedSample(fractions={r: counts[r] / total for r in RECOVERIES}, reported_mean=mean_reported, reported_sd=math.sqrt(variance))
