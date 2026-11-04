# planetgen/physics/mass_radius_sample.py

"""
The mass-radius sample selection model (GEN.505, docs/design/composition/tunable-gaps.md sections 3.1 and 4, docs/design/composition/mass-radius-sample.md)
=================================================================================================================================

Which planets enter the sample of precisely measured planets, mass better than 25 percent and radius better than 8 percent
[R2774], so that a density distribution (the water-rich share of B2, the iron branch of B4) is compared with planets that
could be measured and not with the whole drawn population.

The model is an error model. A transiting planet with a published mass has two measurement errors that decide the cut: the error of the
radial-velocity semi-amplitude K, and the error of the transit depth. Both are fitted to the NASA Exoplanet Archive's transiting planets
with a measured mass and radius (`scripts/fit_mass_radius_sample.py`): ln sigma_K is linear in ln K and ln sigma_depth is linear in
ln depth, each with a log-normal scatter and the two correlated. With `K = 0.08946 M P^(-1/3) M*^(-2/3)` (M in Earth masses, P in
years, M* in solar masses, circular edge-on orbit) and `depth = (R/R*)^2`, the planet passes when `sigma_M / M = sigma_K / K` is under
0.25 and `sigma_R / R = sigma_depth / (2 depth)` is under 0.08.

What it conditions on: the planet is a transiting planet that has a published mass and radius. Whether it transits, is detected and
is chosen for follow-up is not here (GEN.504 is the Kepler detection; there is no TESS model). The fitted dependence on K and on the
depth is weak (slopes 0.89 and 0.84 for all small planets), so the cut tilts a distribution mildly toward large K and deep transits;
for a host cooler than 4000 K the sub-sample `m-dwarf` is the fit.

The item a target's sampler returns is a mapping with `mass_earth`, `radius_earth`, `period_d`, `host_mass_solar` and `host_radius_solar`.
`observe(item)` draws the two errors (two numbers, bound stream) and returns the item with `passes`, `mass_error_share`,
`radius_error_share` and the measured mass, radius and density. `pass_probability(item)` is the same thing without random numbers.
The selection model `precise-sample` keeps an item whose `passes` is true; it is registered as `mass-radius-sample`.

Use:

    item = {"mass_earth": 5.0, "radius_earth": 1.9, "period_d": 4.0, "host_mass_solar": 0.3, "host_radius_solar": 0.3}
    pass_probability(item, "m-dwarf")                    # the chance it is in the sample
    observed_sampler(sampler, "m-dwarf")                 # a sampler for tuning_harness.simulate
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

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "mass_radius_sample.json")
"""str: The data file."""

MODEL_ID = "mass-radius-sample"
"""str: The id of the selection model in `tuning_harness.SELECTION_MODEL_ITEMS`."""

R_EARTH_KM = 6371.0
R_SUN_KM = physical_constants.SOLAR_RADIUS_M / 1000
R_EARTH_CM = 6.371e8
M_EARTH_G = 5.972e27

_PHI_INV = statistics.NormalDist().inv_cdf


def _phi(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


@dataclass(frozen=True)
class Line:
    """`ln sigma = intercept + slope ln x + scatter eps`."""

    intercept: float
    slope: float
    scatter: float

    def mean(self, ln_x):
        return self.intercept + self.slope * ln_x


@dataclass(frozen=True)
class Subsample:
    """
    One fitted sub-sample of the archive.

    Attributes:
        name, description (str): Its name and what it holds.
        n, passing (int): Planets, and those that pass the cut.
        sigma_k, sigma_depth (Line): The two error lines.
        correlation (float): The correlation of the two scatters.
        mean_ln_k, sd_ln_k, mean_ln_depth, sd_ln_depth, population_correlation (float): The Gaussian summary of the planets themselves.
    """

    name: str
    description: str
    n: int
    passing: int
    sigma_k: Line
    sigma_depth: Line
    correlation: float
    mean_ln_k: float
    sd_ln_k: float
    mean_ln_depth: float
    sd_ln_depth: float
    population_correlation: float

    @property
    def pass_fraction(self):
        return self.passing / self.n


@dataclass(frozen=True)
class Model:
    """The fitted model: the cut, the sub-samples and where they came from."""

    mass_cut: float
    radius_cut: float
    k_per_earth_mass_ms: float
    subsamples: dict
    source: str
    retrieved: str


def _finite(owner, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise DataError(f"{owner}: {value!r} is not a finite number")
    return float(value)


def _line(owner, raw):
    scatter = _finite(f"{owner} scatter", raw.get("scatter"))
    if scatter <= 0:
        raise DataError(f"{owner}: the scatter {scatter} is not positive")
    return Line(_finite(f"{owner} intercept", raw.get("intercept")), _finite(f"{owner} slope", raw.get("slope")), scatter)


def load(path=DATA_FILE):
    """Reads and checks the fitted model at `path`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the mass-radius sample model {path}: {error}") from error
    cut = document.get("cut") or {}
    mass_cut, radius_cut = _finite("mass cut", cut.get("mass_error_share")), _finite("radius cut", cut.get("radius_error_share"))
    if not 0 < mass_cut < 1 or not 0 < radius_cut < 1:
        raise DataError("the error shares of the cut are out of range")
    subsamples = {}
    for name, raw in (document.get("subsamples") or {}).items():
        correlation = _finite(f"{name} correlation", raw.get("residual_correlation"))
        population = raw.get("population") or {}
        if not -1 < correlation < 1 or not -1 < _finite(f"{name} population correlation", population.get("correlation")) < 1:
            raise DataError(f"{name}: a correlation is outside (-1, 1)")
        n, passing = raw.get("n"), raw.get("passing")
        if not isinstance(n, int) or not isinstance(passing, int) or not 0 <= passing <= n or n < 10:
            raise DataError(f"{name}: the counts {passing!r} of {n!r} are not a sample")
        subsamples[name] = Subsample(
            name=name, description=raw.get("description", ""), n=n, passing=passing,
            sigma_k=_line(f"{name} sigma_k", raw.get("ln_sigma_k") or {}),
            sigma_depth=_line(f"{name} sigma_depth", raw.get("ln_sigma_depth") or {}), correlation=correlation,
            mean_ln_k=_finite(f"{name} mean_ln_k", population.get("mean_ln_k")), sd_ln_k=_finite(f"{name} sd_ln_k", population.get("sd_ln_k")),
            mean_ln_depth=_finite(f"{name} mean_ln_depth", population.get("mean_ln_depth")),
            sd_ln_depth=_finite(f"{name} sd_ln_depth", population.get("sd_ln_depth")),
            population_correlation=_finite(f"{name} population correlation", population.get("correlation")))
    if not subsamples:
        raise DataError("the model has no sub-samples")
    return Model(mass_cut=mass_cut, radius_cut=radius_cut, k_per_earth_mass_ms=_finite("k_per_earth_mass_ms", document.get("k_per_earth_mass_ms")),
                 subsamples=subsamples, source=document.get("source", ""), retrieved=document.get("retrieved", ""))


MODEL = load()
"""Model: The fitted model."""


def _subsample(name):
    try:
        return MODEL.subsamples[name]
    except KeyError:
        raise DataError(f"sub-sample {name!r} is none of {', '.join(MODEL.subsamples)}") from None


def rv_amplitude_ms(item):
    """K in m/s of a circular, edge-on planet (a mapping with `mass_earth`, `period_d`, `host_mass_solar`)."""
    return (MODEL.k_per_earth_mass_ms * item["mass_earth"] * (item["period_d"] / 365.25) ** (-1.0 / 3.0)
            * item["host_mass_solar"] ** (-2.0 / 3.0))


def transit_depth_ppm(item):
    """The transit depth (R/R*)^2 in ppm (a mapping with `radius_earth`, `host_radius_solar`)."""
    return (item["radius_earth"] * R_EARTH_KM / (item["host_radius_solar"] * R_SUN_KM)) ** 2 * 1e6


def _thresholds(item, sub):
    """`(mass threshold, radius threshold)` in units of the scatter: the errors pass when their normal variate is below these."""
    ln_k, ln_d = math.log(rv_amplitude_ms(item)), math.log(transit_depth_ppm(item))
    c_k = (math.log(MODEL.mass_cut * rv_amplitude_ms(item)) - sub.sigma_k.mean(ln_k)) / sub.sigma_k.scatter
    c_d = (math.log(2.0 * MODEL.radius_cut * transit_depth_ppm(item)) - sub.sigma_depth.mean(ln_d)) / sub.sigma_depth.scatter
    return c_k, c_d


def both_below(c1, c2, rho, nodes=200):
    """P(Z1 < c1 and Z2 < c2) for standard normals of correlation `rho`, by quadrature over Z1."""
    if abs(rho) >= 1:
        raise DataError("the correlation must be inside (-1, 1)")
    root = math.sqrt(1.0 - rho * rho)
    total = 0.0
    low = _phi(c1)
    if low <= 0:
        return 0.0
    for k in range(nodes):
        z1 = _PHI_INV(low * (k + 0.5) / nodes)
        total += _phi((c2 - rho * z1) / root)
    return low * total / nodes


def pass_probability(item, subsample="small"):
    """The chance that a planet `item` passes the cut (mass error under 25 percent and radius error under 8 percent)."""
    sub = _subsample(subsample)
    c_k, c_d = _thresholds(item, sub)
    return both_below(c_k, c_d, sub.correlation)


@dataclass(frozen=True)
class Measurement:
    """
    A planet as the sample would report it.

    Attributes:
        mass_error_share, radius_error_share (float): sigma_M/M and sigma_R/R the measurement would have.
        passes (bool): Both are under the cut.
        mass_earth, radius_earth, density_g_cm3 (float): The measured values (the true ones with a normal error of those shares).
    """

    mass_error_share: float
    radius_error_share: float
    passes: bool
    mass_earth: float
    radius_earth: float
    density_g_cm3: float


def observe(item, subsample="small"):
    """
    The item with the measurement the sample would make: four numbers from the bound draw stream (the two correlated error
    variates, then the errors of the mass and of the radius).
    """
    sub = _subsample(subsample)
    z1, z2, zm, zr = (draw.gauss(0.0, 1.0) for _ in range(4))
    z_depth = sub.correlation * z1 + math.sqrt(1.0 - sub.correlation ** 2) * z2
    k, depth = rv_amplitude_ms(item), transit_depth_ppm(item)
    sigma_k = math.exp(sub.sigma_k.mean(math.log(k)) + sub.sigma_k.scatter * z1)
    sigma_depth = math.exp(sub.sigma_depth.mean(math.log(depth)) + sub.sigma_depth.scatter * z_depth)
    mass_share, radius_share = sigma_k / k, sigma_depth / (2.0 * depth)
    mass = max(item["mass_earth"] * (1.0 + mass_share * zm), 1e-6 * item["mass_earth"])
    radius = max(item["radius_earth"] * (1.0 + radius_share * zr), 1e-6 * item["radius_earth"])
    density = mass * M_EARTH_G / (4.0 / 3.0 * math.pi * (radius * R_EARTH_CM) ** 3)
    return Measurement(mass_error_share=mass_share, radius_error_share=radius_share,
                       passes=mass_share < MODEL.mass_cut and radius_share < MODEL.radius_cut,
                       mass_earth=mass, radius_earth=radius, density_g_cm3=density)


def observed_sampler(sampler, subsample="small"):
    """A sampler for `tuning_harness.simulate`: each item of `sampler(i)` becomes itself plus the measurement (`passes`,
    `mass_error_share`, `radius_error_share`, `measured_mass_earth`, `measured_radius_earth`, `measured_density_g_cm3`)."""
    def sample(i):
        item = sampler(i)
        if item is None:
            return None
        measured = observe(item, subsample)
        return dict(item, passes=measured.passes, mass_error_share=measured.mass_error_share,
                    radius_error_share=measured.radius_error_share, measured_mass_earth=measured.mass_earth,
                    measured_radius_earth=measured.radius_earth, measured_density_g_cm3=measured.density_g_cm3)
    return sample


def precise_sample(item):
    """The selection model: 1 for an item whose measurement passes the cut, 0 otherwise."""
    return 1.0 if item["passes"] else 0.0


def validity(subsample="small"):
    """The range of the fit: transiting planets under 4 Earth radii (and, for `m-dwarf`, hosts under 4000 K)."""
    spec = {"radius_earth": (None, 4.0)}
    if subsample == "m-dwarf":
        spec["host_teff_k"] = (None, 4000.0)
    return spec


tuning_harness.register_selection_model(MODEL_ID, precise_sample)


def population_pass_fraction(subsample="small", nodes=40):
    """The pass probability averaged over the Gaussian summary of the sub-sample's planets (ln K and ln depth): what the model
    predicts for the fraction of the archive's planets that pass; compare with `Subsample.pass_fraction`."""
    sub = _subsample(subsample)
    rho = sub.population_correlation
    root = math.sqrt(1.0 - rho * rho)
    total = 0.0
    for i in range(nodes):
        z1 = _PHI_INV((i + 0.5) / nodes)
        for j in range(nodes):
            z2 = _PHI_INV((j + 0.5) / nodes)
            ln_k = sub.mean_ln_k + sub.sd_ln_k * z1
            ln_d = sub.mean_ln_depth + sub.sd_ln_depth * (rho * z1 + root * z2)
            c_k = (math.log(MODEL.mass_cut) + ln_k - sub.sigma_k.mean(ln_k)) / sub.sigma_k.scatter
            c_d = (math.log(2.0 * MODEL.radius_cut) + ln_d - sub.sigma_depth.mean(ln_d)) / sub.sigma_depth.scatter
            total += both_below(c_k, c_d, sub.correlation, nodes=60)
    return total / (nodes * nodes)
