# planetgen/physics/saturation.py

"""
The saturation pressure of the gases and ices (GEN.470): `saturation(species, temperature_k)`, the pressure at which the
solid and its vapour are in balance, with the standing of the answer.

The curves are the Fray and Schmitt 2009 polynomials `ln P[bar] = A0 + sum Ai / T^i` (i = 1 to 6) as quoted in the
supplementary table of Lisse et al. (arXiv:2009.02277), and for water ice the Feistel and Wagner 2007 form that IAPWS
R14-08 adopts. `data/saturation_pressure.json` holds the coefficients and says where each came from. This is the one place
a vapour-pressure curve lives; the callers (GEN.452, the cold-end work of GEN.497) read it.

* **A curve is only as good as its window, and the windows were not read.** The paper is behind a paywall and the quoted
  table has no lower limits. Where a species has more than one fit (alpha, beta and gamma O2, alpha and beta N2, CO and
  CH3OH, two Ne, two CO2) the next takes over at the temperature where the two cross, or where they come closest when they
  do not (`from_k` in the file). The switch points land within 1 K of the known solid transitions (N2 35.6 K, CO 61.55 K,
  CH3OH 157.4 K, O2 23.9 K and 43.8 K) but they are derived, not read, so the choice is a stand-in. The two CO2 fits do not
  cross and are 1.1 percent apart at their switch; CO2-1 serves the cold because CO2-2 puts 1e-9 bar at 64.7 K, 21 K colder
  than the latent heat of sublimation allows. H2S-2 and the Ne-1 and CH3OH-1 fits above their switch are not used.
* **Below 50 K the answer is `cold-unverified`.** The polynomials are fits to laboratory data that stop short of the
  coldest ices and the cold sublimation is wrong by model error, not by a sign (docs/design/atmosphere-composition.md,
  finding 11). The laboratory-corrected cold end and the Span and Wagner CO2 curve are GEN.497, which extends this module.
* **Below its `floor_k` a polynomial turns downwards** and its value means nothing (`below-fit`); the pressure returned
  is the polynomial's at the floor, which is below 1e-30 bar for every species and cannot overflow.
* **Above the triple point the liquid curve takes over** (`phase` is `liquid`): the Wagner fits to the NIST saturation
  tables of docs/design/atmosphere-composition-parts/01-gas-table-thermodynamics.md section 6.2 (hydrogen cyanide by the
  Perry correlation), good to 0.3 percent of the table. The solid and liquid curves do not meet exactly at the triple
  point: the step is the "Err at Tt" of that part (SO2 +3.8 percent, benzene +7 percent, N2 -2.8 percent), accepted as
  the published spread, so a pressure inverted across the triple point can have two answers within it. Acetylene has no
  liquid curve (`above-triple`, the value is the polynomial carried on). Above the critical temperature there is no
  saturation: the status is `above-critical` and the pressure is infinite.
* **Known departures at the triple point** from the gas table's pressure: SO2 3.8 percent over, benzene 7 percent over,
  acetylene 3.6 percent under; the other species are within 3 percent (src/tests/test_saturation.py).
* **Hydrogen's solid** is Clausius-Clapeyron from its triple point with a latent heat of sublimation of 1.02 kJ/mol (the
  gas table's 0.90 of vapourisation plus 0.12 of fusion, a stand-in, unchecked against Woolley 1948). Helium and the species
  without a curve raise `UnknownSpecies`.

The functions are pure.
"""

import json
import math
import os
from dataclasses import dataclass
from typing import Optional, Tuple

from planetgen.physics import elements, gases

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "saturation_pressure.json")

BAR_TO_PA = 1.0e5
"""float: Pascals in a bar; the Fray and Schmitt polynomials give `ln P` in bar."""

GAS_CONSTANT_J_MOL_K = 8.314

OK = "ok"
COLD_UNVERIFIED = "cold-unverified"
BELOW_FIT = "below-fit"
ABOVE_TRIPLE = "above-triple"
ABOVE_CRITICAL = "above-critical"
STATUSES = (OK, COLD_UNVERIFIED, BELOW_FIT, ABOVE_TRIPLE, ABOVE_CRITICAL)

SOLID = "solid"
LIQUID = "liquid"


class UnknownSpecies(KeyError):
    """A species with no saturation curve."""


@dataclass(frozen=True)
class Curve:
    """
    One polynomial.

    Attributes:
        curve_id (str): The paper's name for the fit (`"N2-2"`).
        from_k (float): The temperature from which this fit is the one used (0 for the first).
        floor_k (float or None): Below this the polynomial turns downwards and means nothing; `None` when it never does.
        phase (str): The solid phase the fit describes (`"alpha"`, `"beta"`) or `""`.
        coefficients (tuple): A0 to A6.
    """

    curve_id: str
    from_k: float
    floor_k: Optional[float]
    phase: str
    coefficients: Tuple[float, ...]


@dataclass(frozen=True)
class Liquid:
    """
    The liquid branch, from the triple point to `t_max_k`.

    Attributes:
        form (str): `"wagner"` (`P = Pc exp((a x + b x^1.5 + c x^2.5 + d x^5) / Tr)`, `x = 1 - T/Tc`, `Tr = T/Tc`) or
            `"dippr101"` (`ln P = A + B/T + C ln T + D T^E`).
        t_max_k (float): The critical temperature, or the end of the correlation.
        pc_pa (float or None): The critical pressure (Wagner).
        coefficients (tuple): `a, b, c, d` (Wagner) or `A, B, C, D, E` (DIPPR 101).
    """

    form: str
    t_max_k: float
    pc_pa: Optional[float]
    coefficients: Tuple[float, ...]


@dataclass(frozen=True)
class Species:
    """
    One species.

    Attributes:
        formula, name (str): As the gas table writes them.
        form (str): The solid's: `"fray-schmitt"`, `"feistel-wagner"` or `"clausius-clapeyron"`.
        selection (str): How a second fit was placed (`"single"`, `"crossing"`, `"closest-approach"`).
        curves (tuple): The `Curve`s, coldest first (empty for water).
        water_coefficients (tuple): e0 to e6 of the water form (empty otherwise).
        sublimation_kj_mol (float or None): The latent heat of the Clausius-Clapeyron solid.
        liquid (Liquid or None): The liquid branch; `None` where there is none.
        note (str): What is known to be wrong or unread about the row.
    """

    formula: str
    name: str
    form: str
    selection: str
    curves: Tuple[Curve, ...]
    water_coefficients: Tuple[float, ...]
    sublimation_kj_mol: Optional[float]
    liquid: Optional[Liquid]
    note: str


@dataclass(frozen=True)
class Saturation:
    """
    The answer: `pressure_pa` (infinite above the critical point), how far to trust it (`status`, one of `STATUSES`), the
    fit it came from (`curve_id`) and the `phase` of the condensate (`SOLID` or `LIQUID`).
    """

    pressure_pa: float
    status: str
    curve_id: str
    phase: str


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    rows = {}
    for formula, row in data["species"].items():
        curves = tuple(
            Curve(item["id"], float(item["from_k"]), item["floor_k"], item["phase"], tuple(float(a) for a in item["coefficients"]))
            for item in row.get("curves", ()))
        if row["form"] == "fray-schmitt":
            if not curves or curves[0].from_k != 0.0 or any(b.from_k <= a.from_k for a, b in zip(curves, curves[1:])):
                raise ValueError(f"saturation_pressure.json: {formula} curves must start at 0 K and switch upwards")
            if any(len(curve.coefficients) != 7 for curve in curves):
                raise ValueError(f"saturation_pressure.json: {formula} needs the seven coefficients A0 to A6")
        elif row["form"] not in ("feistel-wagner", "clausius-clapeyron"):
            raise ValueError(f"saturation_pressure.json: {formula} has an unknown form {row['form']!r}")
        liquid = None
        if row.get("liquid"):
            item = row["liquid"]
            if item["form"] not in ("wagner", "dippr101"):
                raise ValueError(f"saturation_pressure.json: {formula} has an unknown liquid form {item['form']!r}")
            liquid = Liquid(item["form"], float(item["t_max_k"]), item.get("pc_pa"), tuple(float(a) for a in item["coefficients"]))
        rows[formula] = Species(formula, row["name"], row["form"], row["selection"], curves,
                                tuple(float(e) for e in row.get("coefficients", ())), row.get("sublimation_kj_mol"), liquid,
                                row.get("note", ""))
    return data["cold_unverified_below_k"], rows


COLD_UNVERIFIED_BELOW_K, _TABLE = _load()
"""float: Below this temperature the answer is `cold-unverified` (model error of the fits in the cold, finding 11)."""

WATER_TRIPLE_POINT_K = 273.16
WATER_TRIPLE_POINT_PA = 611.657


def species(formula):
    """The row of `formula`; raises `UnknownSpecies`."""
    try:
        return _TABLE[formula]
    except KeyError:
        raise UnknownSpecies(f"{formula!r} has no saturation curve (physics/data/saturation_pressure.json)") from None


def formulas():
    """The formulas of the table, in file order."""
    return tuple(_TABLE)


def has(formula):
    """Whether `formula` has a saturation curve."""
    return formula in _TABLE


def curve_at(formula, temperature_k):
    """The `Curve` that serves `formula` at `temperature_k` (the last one whose `from_k` is not above it)."""
    row = species(formula)
    if row.form != "fray-schmitt":
        raise ValueError(f"{formula} is not a Fray and Schmitt species")
    chosen = row.curves[0]
    for curve in row.curves[1:]:
        if temperature_k >= curve.from_k:
            chosen = curve
    return chosen


def _fray_schmitt_pa(coefficients, temperature_k):
    ln_bar = coefficients[0]
    for power in range(1, 7):
        ln_bar += coefficients[power] / temperature_k ** power
    return math.exp(ln_bar) * BAR_TO_PA


def _water_pa(coefficients, temperature_k):
    theta = temperature_k / WATER_TRIPLE_POINT_K
    h = 0.0
    for power, e in enumerate(coefficients):
        h += e * theta ** power
    return WATER_TRIPLE_POINT_PA * math.exp(1.5 * math.log(theta) + (1.0 - 1.0 / theta) * h)


def _clausius_clapeyron_pa(formula, sublimation_kj_mol, temperature_k):
    row = gases.gas(formula)
    return row.triple_point_pa * math.exp(
        -sublimation_kj_mol * 1000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / temperature_k - 1.0 / row.triple_point_k))


def _liquid_pa(liquid, temperature_k):
    if liquid.form == "wagner":
        a, b, c, d = liquid.coefficients
        reduced = temperature_k / liquid.t_max_k
        x = 1.0 - reduced
        return liquid.pc_pa * math.exp((a * x + b * x ** 1.5 + c * x ** 2.5 + d * x ** 5) / reduced)
    a, b, c, d, e = liquid.coefficients
    return math.exp(a + b / temperature_k + c * math.log(temperature_k) + d * temperature_k ** e)


def triple_point_k(formula):
    """The triple-point temperature of `formula` (water's from this table, the others from the gas table); `None` if unknown."""
    if species(formula).form == "feistel-wagner":
        return WATER_TRIPLE_POINT_K
    return gases.gas(formula).triple_point_k


def saturation(formula, temperature_k):
    """
    The saturation pressure of `formula` at `temperature_k`, Pa, with its status and phase.

    Below the triple point the solid's curve, above it the liquid's, up to the critical temperature. The status is the
    first that applies of `above-critical` (no saturation: the pressure is infinite), `above-triple` (over the triple
    point of a species with no liquid curve: the solid's polynomial carried on), `below-fit` (under the fit's `floor_k`,
    where the value is the one at the floor), `cold-unverified` (under `COLD_UNVERIFIED_BELOW_K`), else `ok`.

    Raises:
        UnknownSpecies: There is no curve for `formula`.
        ValueError: The temperature is not above zero.
    """
    row = species(formula)
    if not temperature_k > 0:
        raise ValueError("a temperature must be positive")
    triple = triple_point_k(formula)
    liquid = row.liquid
    if liquid is not None and triple is not None and temperature_k > triple:
        if temperature_k > liquid.t_max_k:
            return Saturation(math.inf, ABOVE_CRITICAL, f"{formula}-{liquid.form}", LIQUID)
        return Saturation(_liquid_pa(liquid, temperature_k), OK, f"{formula}-{liquid.form}", LIQUID)
    floor = None
    if row.form == "feistel-wagner":
        pressure = _water_pa(row.water_coefficients, temperature_k)
        curve_id = "H2O-FW"
    elif row.form == "clausius-clapeyron":
        pressure = _clausius_clapeyron_pa(formula, row.sublimation_kj_mol, temperature_k)
        curve_id = f"{formula}-CC"
    else:
        curve = curve_at(formula, temperature_k)
        floor = curve.floor_k
        pressure = _fray_schmitt_pa(curve.coefficients, temperature_k if floor is None else max(temperature_k, floor))
        curve_id = curve.curve_id
    if triple is not None and temperature_k > triple:
        status = ABOVE_TRIPLE
    elif floor is not None and temperature_k < floor:
        status = BELOW_FIT
    elif temperature_k < COLD_UNVERIFIED_BELOW_K:
        status = COLD_UNVERIFIED
    else:
        status = OK
    return Saturation(pressure, status, curve_id, SOLID)


def pressure_pa(formula, temperature_k):
    """`saturation(...).pressure_pa`, for a caller that has no use for the status."""
    return saturation(formula, temperature_k).pressure_pa
