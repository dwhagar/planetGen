# planetgen/physics/kepler_recovery.py

"""
The Kepler transit-plus-recovery selection model (GEN.504, docs/design/composition/tunable-gaps.md sections 3.1 and 4, docs/design/composition/kepler-recovery.md)
==========================================================================================================================================

Whether the Kepler pipeline would have found a transiting planet, so that a tuned feature of the Kepler planets (the radius valley, the
period-radius plane, the metallicity null of the California-Kepler Survey) is compared with the planets Kepler could see and not with
the planets that exist. Three factors multiply:

1. **Geometry.** A circular orbit transits when the impact parameter `b` is under 1: `R*/a` of the orientations.
2. **Signal.** The multiple event statistic of the transits, `MES = g(b) k^2 sqrt(N) / CDPP(T)`, where `k = Rp/R*`, `g(b)` the effective
   depth over `k^2` (limb darkening and aperture correction, `c0 + c2 b^2 + c4 b^4`), `N` the number of valid transits (`factor` times the
   data span over the period), `T` the transit duration and `CDPP(T)` the host's combined differential photometric precision at that
   duration (log-linear between the tabulated 1.5, 3, 6, 12 and 15 hours). A host with no measured CDPP gets one from the host's Kepler
   magnitude and radius: `ln CDPP = a0 + a1 Kp + a2 Kp^2 + a3 ln R* + scatter z` with `z` one standard normal for all durations.
3. **Recovery.** The pipeline found `c0 GammaCDF(MES; a, b) S((N - n0) / wn) S((ln Tmax - ln T) / wd)` of the injected planets: the form of
   [R0594], with a three-transit gate and a long-duration cut.

All of it is fitted by `scripts/fit_kepler_recovery.py` to the 146,257 Kepler DR25 on-target injections (INJ1) joined to the DR25 stellar
table. It conditions on a Kepler target star, circular orbits and the pipeline's detection, not on vetting (the Robovetter removes more),
and not on a planet's mass.

The item a target's sampler returns is a mapping with `period_d`, `radius_earth`, `host_radius_solar` and `host_mass_solar`; optionally
`kepmag` and `cdpp_z` (else `observed_sampler` draws them, two numbers from the bound stream) and `impact_parameter` (else the geometry
and the impact parameter are averaged over). Validity: periods 0.5 to 500 days, radius ratios 0.005 to 0.07, the Kepler magnitude of the
targets.

Use:

    item = observed_sampler(sampler)(i)                 # adds kepmag and cdpp_z
    detection_probability(item)                         # the chance Kepler would have found the planet
    expected_mes(item, 0.3)                             # its signal at impact parameter 0.3
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen.generation import tuning_harness
from planetgen.physics import constants as physical_constants
from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "kepler_recovery.json")
"""str: The data file."""

MODEL_ID = "kepler-recovery"
"""str: The id of the selection model in `tuning_harness.SELECTION_MODEL_ITEMS`."""

G_SI = physical_constants.G
M_SUN_KG = physical_constants.SOLAR_MASS_TO_KG
R_SUN_M = physical_constants.SOLAR_RADIUS_M
R_EARTH_M = physical_constants.EARTH_RADIUS_KM * 1000.0
DAY_S = 86400.0
IMPACT_NODES = 10
"""int: The midpoints of `b` in [0, 1) the quadrature averages over."""


@dataclass(frozen=True)
class Recovery:
    """`P = c0 GammaCDF(MES; a, b) S((N - n0) / wn) S((ln Tmax - ln T) / wd)`."""

    c0: float
    a: float
    b: float
    n0: float
    wn: float
    ln_t_max: float
    wd: float


@dataclass(frozen=True)
class NoiseLine:
    """`ln CDPP(duration) = intercept + slope_kepmag Kp + slope_kepmag2 Kp^2 + slope_ln_radius ln R* + scatter z`."""

    duration_h: float
    intercept: float
    slope_kepmag: float
    slope_kepmag2: float
    slope_ln_radius: float
    scatter: float


@dataclass(frozen=True)
class Model:
    """The fitted model."""

    span_days: float
    transit_count_factor: float
    limb: tuple
    recovery: Recovery
    noise: tuple
    kepmag_quantiles: tuple
    period_range_d: tuple
    radius_ratio_range: tuple
    injections: int
    recovered: int
    predicted_recovered: float
    generated_host: dict
    calibration: tuple
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
        raise DataError(f"cannot read the Kepler recovery model {path}: {error}") from error
    raw = document.get("recovery") or {}
    recovery = Recovery(*(_finite(f"recovery {name}", raw.get(name)) for name in ("c0", "a", "b", "n0", "wn", "ln_t_max", "wd")))
    if not 0 < recovery.c0 <= 1 or min(recovery.a, recovery.b, recovery.wn, recovery.wd) <= 0:
        raise DataError(f"the recovery curve {raw!r} is out of range")
    limb = tuple(_finite(f"limb {name}", (document.get("limb") or {}).get(name)) for name in ("c0", "c2", "c4"))
    if limb[0] <= 0:
        raise DataError(f"limb c0 {limb[0]} is not positive")
    lines = []
    for entry in (document.get("noise") or {}).get("durations") or []:
        line = NoiseLine(*(_finite(f"noise {name}", entry.get(name)) for name in
                           ("duration_h", "intercept", "slope_kepmag", "slope_kepmag2", "slope_ln_radius", "scatter")))
        if line.scatter <= 0 or (lines and line.duration_h <= lines[-1].duration_h):
            raise DataError(f"noise line {entry!r} has no scatter or is not in order of duration")
        lines.append(line)
    if len(lines) < 2:
        raise DataError("the noise model needs at least two durations")
    quantiles = tuple(_finite("kepmag quantile", q) for q in (document.get("noise") or {}).get("kepmag_quantiles") or [])
    if len(quantiles) < 3 or any(b < a for a, b in zip(quantiles, quantiles[1:])):
        raise DataError("the Kepler magnitude quantiles are missing or not sorted")
    span = _finite("span_days", document.get("span_days"))
    factor = _finite("transit_count_factor", document.get("transit_count_factor"))
    if span <= 0 or factor <= 0:
        raise DataError("the data span and the transit count factor must be positive")
    validity = document.get("validity") or {}
    period_range = tuple(_finite("period range", v) for v in validity.get("period_d") or ())
    ratio_range = tuple(_finite("radius ratio range", v) for v in validity.get("radius_ratio") or ())
    if len(period_range) != 2 or len(ratio_range) != 2 or not 0 < period_range[0] < period_range[1] or not 0 < ratio_range[0] < ratio_range[1]:
        raise DataError("the validity ranges are missing")
    counts = document.get("injections") or {}
    return Model(span_days=span, transit_count_factor=factor, limb=limb, recovery=recovery, noise=tuple(lines),
                 kepmag_quantiles=quantiles, period_range_d=period_range, radius_ratio_range=ratio_range,
                 injections=int(counts.get("used", 0)), recovered=int(counts.get("recovered", 0)),
                 predicted_recovered=_finite("predicted_recovered", document.get("predicted_recovered")),
                 generated_host=dict(document.get("generated_host") or {}), calibration=tuple(document.get("calibration") or ()),
                 source=document.get("source", ""), retrieved=document.get("retrieved", ""))


MODEL = load()
"""Model: The fitted model."""


def validity():
    """The range of the fit."""
    return {"period_d": MODEL.period_range_d, "radius_ratio": MODEL.radius_ratio_range,
            "kepmag": (MODEL.kepmag_quantiles[0], MODEL.kepmag_quantiles[-1])}


def _sigmoid(z):
    return 0.5 * (1.0 + math.tanh(0.5 * z))


def gamma_cdf(x, shape):
    """The regularized lower incomplete gamma function P(shape, x) (the CDF of a unit-scale gamma variable)."""
    if x <= 0.0:
        return 0.0
    log_prefactor = shape * math.log(x) - x - math.lgamma(shape)
    if x < shape + 1.0:
        term = total = 1.0 / shape
        n = shape
        for _ in range(1000):
            n += 1.0
            term *= x / n
            total += term
            if abs(term) < abs(total) * 1e-14:
                break
        return min(1.0, total * math.exp(log_prefactor))
    tiny = 1e-300
    b = x + 1.0 - shape
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - shape)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return max(0.0, 1.0 - math.exp(log_prefactor) * h)


def radius_ratio(item):
    """`k = Rp / R*` of a planet item."""
    return item["radius_earth"] * R_EARTH_M / (item["host_radius_solar"] * R_SUN_M)


def scaled_semi_major_axis(item):
    """`a / R*` from Kepler's third law."""
    period_s = item["period_d"] * DAY_S
    a_m = (G_SI * item["host_mass_solar"] * M_SUN_KG * period_s * period_s / (4.0 * math.pi ** 2)) ** (1.0 / 3.0)
    return a_m / (item["host_radius_solar"] * R_SUN_M)


def transit_duration_h(item, impact_parameter):
    """The first-to-fourth-contact duration in hours of a circular orbit; 0 when the planet misses the star."""
    k, a_r = radius_ratio(item), scaled_semi_major_axis(item)
    chord = (1.0 + k) ** 2 - impact_parameter ** 2
    if chord <= 0.0 or a_r <= 1.0 + k:
        return 0.0
    return item["period_d"] / math.pi * math.asin(min(1.0, math.sqrt(chord) / a_r)) * 24.0


def transit_count(item):
    """The valid transits in the data span (gaps included)."""
    return MODEL.transit_count_factor * MODEL.span_days / item["period_d"]


def depth_factor(impact_parameter):
    """The effective transit depth over `k^2` at an impact parameter."""
    c0, c2, c4 = MODEL.limb
    return c0 + c2 * impact_parameter ** 2 + c4 * impact_parameter ** 4


def cdpp(item, duration_h):
    """
    The host's CDPP at `duration_h` as a fraction of the flux: from the item's `kepmag` and `cdpp_z` (0 when absent) and the host
    radius, log-linear in duration between the tabulated durations and held at the first and last outside them.
    """
    try:
        kp = item["kepmag"]
    except KeyError:
        raise DataError("the item has no kepmag; use observed_sampler or give one") from None
    z = item.get("cdpp_z", 0.0)
    ln_radius = math.log(item["host_radius_solar"])
    values = [line.intercept + line.slope_kepmag * kp + line.slope_kepmag2 * kp * kp + line.slope_ln_radius * ln_radius
              + line.scatter * z for line in MODEL.noise]
    t = min(max(math.log(duration_h), math.log(MODEL.noise[0].duration_h)), math.log(MODEL.noise[-1].duration_h))
    for i in range(len(MODEL.noise) - 1):
        low, high = math.log(MODEL.noise[i].duration_h), math.log(MODEL.noise[i + 1].duration_h)
        if t <= high:
            w = (t - low) / (high - low)
            return math.exp(values[i] * (1.0 - w) + values[i + 1] * w) * 1e-6
    return math.exp(values[-1]) * 1e-6


def expected_mes(item, impact_parameter):
    """The multiple event statistic of the planet's transits at an impact parameter (0 when it does not transit)."""
    duration = transit_duration_h(item, impact_parameter)
    if duration <= 0.0:
        return 0.0
    k = radius_ratio(item)
    return depth_factor(impact_parameter) * k * k * math.sqrt(transit_count(item)) / cdpp(item, duration)


def recovery_probability(mes, count, duration_h):
    """The chance the pipeline recovers a transit signal of `mes` with `count` valid transits and duration `duration_h`."""
    r = MODEL.recovery
    if mes <= 0.0 or duration_h <= 0.0:
        return 0.0
    return (r.c0 * gamma_cdf(mes / r.b, r.a) * _sigmoid((count - r.n0) / r.wn)
            * _sigmoid((r.ln_t_max - math.log(duration_h)) / r.wd))


def recovery_at(item, impact_parameter):
    """The chance Kepler recovers the planet if it transits at `impact_parameter` (under 1; a grazing transit is not in the fit)."""
    if not 0.0 <= impact_parameter < 1.0:
        return 0.0
    duration = transit_duration_h(item, impact_parameter)
    return recovery_probability(expected_mes(item, impact_parameter), transit_count(item), duration)


def detection_probability(item):
    """
    The chance Kepler would have found a planet `item`. With an `impact_parameter` the planet's orientation is given: it transits only
    if that is under 1, and the chance is its recovery. Without one the chance is the geometric `R*/a` times the recovery averaged over
    the impact parameter.
    """
    if "impact_parameter" in item:
        return recovery_at(item, item["impact_parameter"])
    a_r = scaled_semi_major_axis(item)
    if a_r <= 1.0:
        return 0.0
    mean = sum(recovery_at(item, (i + 0.5) / IMPACT_NODES) for i in range(IMPACT_NODES)) / IMPACT_NODES
    return min(1.0, 1.0 / a_r) * mean


def with_noise(item):
    """
    The item with a Kepler magnitude and a noise draw where it has none: two numbers from the bound draw stream (the magnitude by the
    inverse CDF of the targets' quantiles, then a standard normal).
    """
    out = dict(item)
    if "kepmag" not in out:
        quantiles = MODEL.kepmag_quantiles
        position = draw.random() * (len(quantiles) - 1)
        low = min(int(position), len(quantiles) - 2)
        out["kepmag"] = quantiles[low] + (position - low) * (quantiles[low + 1] - quantiles[low])
    if "cdpp_z" not in out:
        out["cdpp_z"] = draw.gauss(0.0, 1.0)
    return out


def observed_sampler(sampler):
    """A sampler for `tuning_harness.simulate`: each item of `sampler(i)` gets a Kepler magnitude and a noise draw (`None` stays `None`)."""
    def sample(i):
        item = sampler(i)
        return None if item is None else with_noise(item)
    return sample


def kepler_recovery(item):
    """The selection model: the probability Kepler would have found the planet."""
    return detection_probability(item)


def expected_found(items):
    """The number of `items` (each with `kepmag`) Kepler would have found: the sum of their probabilities."""
    return sum(detection_probability(item) for item in items)


tuning_harness.register_selection_model(MODEL_ID, kepler_recovery)
