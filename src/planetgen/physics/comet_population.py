# planetgen/physics/comet_population.py

"""
The comet volatile population (GEN.495, docs/design/composition/body-composition-parts/G1-comets.md sections 2.3, 3.1 and 3.2)
==================================================================================================================

What a comet's coma holds, as a ratio to water in percent of the number of molecules, drawn from the population of measured
comets and not from one comet. The species are log-normal with the geometric mean and the scatter in dex of the population
([R3782] database, coma ratios at a heliocentric distance of at most 2 AU), truncated at the observed extremes;
the hyper-volatile block keeps the observed trade-off of CO and CO2 and the carbon budget:

1. S = (CO + CO2)/H2O, log-normal, median 18 percent (Harrington [R3535]), 0.2 dex: the scatter among single comets
   computed from the [R3782] database (GEN.629), not the uncertainty of the median;
2. f = CO/(CO + CO2), logit-normal with the median of the dynamical class (0.08 for Jupiter-family and Halley-type comets, 0.5
   for long-period ones) and a logit scatter of 1.5; CO = S f and CO2 = S (1 - f) are redrawn together while either lies outside
   its observed bounds (the retry is bounded and a failure is an error, never a clipped value);
3. every other species an independent truncated log-normal. The class matters only for f and for the carbon-chain-depleted flag
   (the class explains 28 percent of the variance of CO and under 10 percent for every other species);
4. the carbon-chain-depleted flag is a draw of the possibility-table loader (`CARBON_CHAIN_TABLE`, gated rows per
   dynamical class: 0.37 Jupiter-family, 0.185 long-period, 0.25 Halley-type).

An interstellar comet draws its CO and CO2 from the wider ranges one measured object forces (CO/H2O up to 2, CO2/H2O to 10; both
log-uniform) and the other species from the population. The one comet whose bulk ratios were measured in situ is a test target in the
test, never a centre or a preset. Nothing here reads a Solar System comet: the numbers are the population's.

Every row carries its source, label and source state (`SPECIES`, with the partly verified ones naming their limit). The draw reads only `draw` (the one random wrapper), in a fixed order, so a comet's
volatiles are the same at any worker count. The comet record does not store them yet (GEN.442, GEN.455, GEN.457 do), so nothing
generated changes and there is no epoch of its own.

Use:

    draw_volatiles("jupiter_family")                 # a Volatiles: percent of H2O by species, the flag
    draw_volatiles("long_period", interstellar=True)
    SPECIES["CH3OH"].geo_mean_percent                # the population geometric mean
"""

import json
import math
import os
import statistics
from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation import possibility_table as pt
from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "comet_population.json")
"""str: The data file."""

PERIOD_CLASSES = tuple(tuning.COMET_PERIOD_CLASSES)
"""tuple: The dynamical classes of a bound comet (`tuning.COMET_PERIOD_CLASSES`)."""

MAX_ATTEMPTS = 1000
"""int: How many times the hyper-volatile block may be redrawn before the draw is an error."""


class PopulationError(RuntimeError):
    """The hyper-volatile block could not be drawn inside the observed bounds."""


@dataclass(frozen=True)
class Species:
    """
    One species of the coma.

    Attributes:
        name (str): The formula.
        geo_mean_percent (float): The population geometric mean, percent of H2O.
        sigma_dex (float): The scatter of log10 of the ratio.
        lo_percent, hi_percent (float): The observed extremes the draw is truncated at.
        n (int): How many comets the statistics rest on.
        source (str): Key into `SOURCES`.
        label (str): READ, COMPUTED, ...
        state (str): The source state (`possibility_table.STATES`).
        limit (str): For a partly verified row, its limit.
    """

    name: str
    geo_mean_percent: float
    sigma_dex: float
    lo_percent: float
    hi_percent: float
    n: int
    source: str
    label: str
    state: str
    limit: str

    def draw(self):
        """One truncated log-normal draw in percent, with one random number (inverse transform)."""
        mu = math.log10(self.geo_mean_percent)
        low = _phi((math.log10(self.lo_percent) - mu) / self.sigma_dex)
        high = _phi((math.log10(self.hi_percent) - mu) / self.sigma_dex)
        q = min(max(low + draw.random() * (high - low), 1e-16), 1.0 - 1e-16)
        x = mu + self.sigma_dex * statistics.NormalDist().inv_cdf(q)
        return min(max(10.0 ** x, self.lo_percent), self.hi_percent)

    def cdf(self, percent):
        """The probability the truncated log-normal lies at or below `percent`."""
        mu = math.log10(self.geo_mean_percent)
        low = _phi((math.log10(self.lo_percent) - mu) / self.sigma_dex)
        high = _phi((math.log10(self.hi_percent) - mu) / self.sigma_dex)
        if percent <= self.lo_percent:
            return 0.0
        if percent >= self.hi_percent:
            return 1.0
        return (_phi((math.log10(percent) - mu) / self.sigma_dex) - low) / (high - low)


@dataclass(frozen=True)
class Volatiles:
    """
    A comet's drawn volatiles.

    Attributes:
        percent (dict): Percent of H2O by species, CO and CO2 first, then the others in the file's order.
        carbon_chain_depleted (bool): The flag.
        depleted_probability (float): The probability the flag was drawn with.
    """

    percent: dict
    carbon_chain_depleted: bool
    depleted_probability: float


def _phi(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def _positive(owner, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise DataError(f"{owner}: {value!r} is not a positive number")
    return float(value)


def _provenance(owner, item, sources):
    for key in ("source", "label", "state"):
        if not item.get(key):
            raise DataError(f"{owner} has no {key}")
    if item["source"] not in sources:
        raise DataError(f"{owner}: the source {item['source']!r} is not in the file's sources")
    if item["state"] not in (pt.VERIFIED, pt.PARTLY_VERIFIED):
        raise DataError(f"{owner}: the state {item['state']!r} is neither verified nor partly verified")
    if item["state"] == pt.PARTLY_VERIFIED and not item.get("limit"):
        raise DataError(f"{owner} is partly verified and names no limit")


def _species(name, item, sources):
    _provenance(name, item, sources)
    lo, hi = _positive(name, item["lo_percent"]), _positive(name, item["hi_percent"])
    mean, sigma = _positive(name, item["geo_mean_percent"]), _positive(name, item["sigma_dex"])
    if not lo < mean < hi:
        raise DataError(f"{name}: the geometric mean {mean} is not inside its bounds {lo} to {hi}")
    return Species(name, mean, sigma, lo, hi, int(item["n"]), item["source"], item["label"], item["state"],
                   item.get("limit", ""))


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        return _parse(json.load(handle))


def _parse(data):
    for key in ("title", "retrieved", "part", "sources", "total_carbon_oxide", "co_fraction", "carbon_oxides", "species",
                "carbon_chain_depleted", "interstellar"):
        if not data.get(key):
            raise DataError(f"comet_population.json states no {key}")
    sources = data["sources"]
    total = data["total_carbon_oxide"]
    _provenance("total_carbon_oxide", total, sources)
    _positive("total_carbon_oxide", total["median_percent"])
    _positive("total_carbon_oxide", total["sigma_dex"])
    fraction = data["co_fraction"]
    _provenance("co_fraction", fraction, sources)
    for period in PERIOD_CLASSES:
        median = fraction["median"].get(period)
        if not (isinstance(median, (int, float)) and 0.0 < median < 1.0):
            raise DataError(f"co_fraction: the median of {period} is not a fraction between 0 and 1")
    _positive("co_fraction", fraction["logit_sigma"])
    flag = data["carbon_chain_depleted"]
    _provenance("carbon_chain_depleted", flag, sources)
    for period in PERIOD_CLASSES:
        value = flag["probability"].get(period)
        if not (isinstance(value, (int, float)) and 0.0 < value < 1.0):
            raise DataError(f"carbon_chain_depleted: the probability for {period} is not between 0 and 1")
    if set(data["interstellar"]) != {"co_over_h2o", "co2_over_h2o"}:
        raise DataError("interstellar must hold co_over_h2o and co2_over_h2o")
    for name, item in data["interstellar"].items():
        _provenance(f"interstellar {name}", item, sources)
        if not 0 < _positive(name, item["lo"]) < _positive(name, item["hi"]):
            raise DataError(f"interstellar {name}: the range is not 0 < low < high")
    oxides = {name: _species(name, item, sources) for name, item in data["carbon_oxides"].items()}
    if set(oxides) != {"CO", "CO2"}:
        raise DataError("carbon_oxides must hold CO and CO2")
    species = {s["name"]: _species(s["name"], s, sources) for s in data["species"]}
    if len(species) != len(data["species"]):
        raise DataError("a species appears twice")
    if set(species) & {"CO", "CO2"}:
        raise DataError("CO and CO2 are drawn as the hyper-volatile block, not as independent species")
    return data, oxides, species


_DATA, CARBON_OXIDES, SPECIES = _load()
SOURCES = _DATA["sources"]
"""dict: Source key -> citation."""

TOTAL_CARBON_OXIDE = _DATA["total_carbon_oxide"]
CO_FRACTION = _DATA["co_fraction"]
INTERSTELLAR = _DATA["interstellar"]
DEPLETED = _DATA["carbon_chain_depleted"]


def _carbon_chain_rows():
    rows = []
    for period, probability in DEPLETED["probability"].items():
        for depleted, weight in ((True, probability), (False, 1.0 - probability)):
            kind = "carbon-chain-depleted" if depleted else "carbon-chain-normal"
            rows.append(pt.Row(f"{kind}-{period}", ("comet",), weight, SOURCES[DEPLETED["source"]], state=DEPLETED["state"],
                               conditions=(pt.Condition("period_class", equals=period),), limit=DEPLETED["limit"],
                               notes=f"{kind} for a {period} comet"))
    return rows


CARBON_CHAIN_TABLE = pt.REGISTRY.register(pt.Table(
    "comet-carbon-chain", 1, _carbon_chain_rows(),
    description="Whether a comet is depleted in carbon-chain species (C2 and C3), by dynamical class (`period_class`)."))
"""Table: The carbon-chain-depleted flag, one pair of gated rows per dynamical class."""


def _gauss():
    return draw.gauss(0.0, 1.0)


def _carbon_oxides(period_class, interstellar, ice_ratios=None):
    if interstellar:
        co = INTERSTELLAR["co_over_h2o"]
        co2 = INTERSTELLAR["co2_over_h2o"]
        return (100.0 * math.exp(draw.uniform(math.log(co["lo"]), math.log(co["hi"]))),
                100.0 * math.exp(draw.uniform(math.log(co2["lo"]), math.log(co2["hi"]))))
    co_bounds, co2_bounds = CARBON_OXIDES["CO"], CARBON_OXIDES["CO2"]
    if ice_ratios is None:
        total_median = TOTAL_CARBON_OXIDE["median_percent"]
        median = CO_FRACTION["median"][period_class]
        co_range = (co_bounds.lo_percent, co_bounds.hi_percent)
        co2_range = (co2_bounds.lo_percent, co2_bounds.hi_percent)
    else:
        # the median comes from the disc zone the comet formed in; the observed extremes keep their distance from the population's
        # geometric mean (the bounds move with the median, so the rejection loop is as likely to end as with the population)
        co2_zone, co_zone = (100.0 * ratio for ratio in ice_ratios)
        if not (co2_zone > 0.0 and co_zone > 0.0):
            return 0.0, 0.0     # a star with no carbon freezes no carbon oxide
        total_median = co_zone + co2_zone
        median = co_zone / total_median
        co_range = tuple(co_zone * edge / co_bounds.geo_mean_percent for edge in (co_bounds.lo_percent, co_bounds.hi_percent))
        co2_range = tuple(co2_zone * edge / co2_bounds.geo_mean_percent for edge in (co2_bounds.lo_percent, co2_bounds.hi_percent))
    log_median = math.log10(total_median)
    logit = math.log(median / (1.0 - median))
    for _ in range(MAX_ATTEMPTS):
        total = 10.0 ** (log_median + TOTAL_CARBON_OXIDE["sigma_dex"] * _gauss())
        share = 1.0 / (1.0 + math.exp(-(logit + CO_FRACTION["logit_sigma"] * _gauss())))
        co, co2 = total * share, total * (1.0 - share)
        if co_range[0] <= co <= co_range[1] and co2_range[0] <= co2 <= co2_range[1]:
            return co, co2
        draw.redrawn("comet.composition")
    raise PopulationError(f"no CO and CO2 inside the observed bounds in {MAX_ATTEMPTS} draws for {period_class}")


def median_volatiles(period_class):
    """
    The volatiles at the medians of the population, with no draw: percent of H2O by species, CO and CO2 from the median total and the
    median CO share of the class, the others at their geometric means. A stand-in for a body that needs the proportions of its
    non-water ices and has no draw of its own to spend.

    Raises:
        ValueError: An unknown class.
    """
    if period_class not in PERIOD_CLASSES:
        raise ValueError(f"period class {period_class!r} is none of {', '.join(PERIOD_CLASSES)}")
    total = TOTAL_CARBON_OXIDE["median_percent"]
    share = CO_FRACTION["median"][period_class]
    percent = {"CO": total * share, "CO2": total * (1.0 - share)}
    for name, species in SPECIES.items():
        percent[name] = species.geo_mean_percent
    return percent


def draw_volatiles(period_class, interstellar=False, ice_ratios=None):
    """
    The volatiles of one comet.

    Args:
        period_class (str): `jupiter_family`, `halley_type` or `long_period` (an interstellar comet takes the class of its
            carbon-chain flag from the same table; the class does not move its CO and CO2, which are drawn from the ranges
            the measured interstellar object needs).
        interstellar (bool): True for an interstellar comet.
        ice_ratios (tuple or None): `(CO2/H2O, CO/H2O)` of the ice of the disc zone the comet formed in (`comet_ice_budget.ice_ratios`, GEN.595):
            the median CO and CO2 come from it (and not from the Solar System's measured median of the class) and the scatter stays the
            population's. `None` keeps the population's median.

    Returns:
        Volatiles: percent of H2O by species and the carbon-chain-depleted flag.

    Raises:
        ValueError: an unknown class. PopulationError: the hyper-volatile block found no draw inside the bounds.
    """
    if period_class not in PERIOD_CLASSES:
        raise ValueError(f"period class {period_class!r} is none of {', '.join(PERIOD_CLASSES)}")
    co, co2 = _carbon_oxides(period_class, interstellar, ice_ratios)
    percent = {"CO": co, "CO2": co2}
    for name, species in SPECIES.items():
        percent[name] = species.draw()
    pick = CARBON_CHAIN_TABLE.draw("comet", {"period_class": period_class})
    depleted = pick.name.startswith("carbon-chain-depleted")
    return Volatiles(percent, depleted, pick.probability if depleted else 1.0 - pick.probability)
