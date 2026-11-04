# planetgen/physics/saturation.py

"""
The saturation pressure of the gases and ices (GEN.470): `saturation(species, temperature_k)`, the pressure at which the
solid and its vapour are in balance, with the standing of the answer.

The curves are the [R1445] polynomials `ln P[bar] = A0 + sum Ai / T^i` (i = 1 to 6) as quoted in the
supplementary table of Lisse et al. [R2711], for water ice the [R1345] form that IAPWS
R14-08 adopts, and the laboratory-corrected cold ends and the carbon dioxide equation described below.
`data/saturation_pressure.json` holds the coefficients and says where each came from. This is the one place a
vapour-pressure curve lives; the callers (GEN.452) read it.

* **A curve is only as good as its window, and the windows were not read.** The paper is behind a paywall and the quoted
  table has no lower limits. Where a species has more than one fit (alpha, beta and gamma O2, alpha and beta N2, CO and
  CH3OH, two Ne) the next takes over at the temperature where the two cross (`from_k` in the file). The switch points land within 1 K of the known solid transitions (N2 35.6 K, CO 61.55 K,
  CH3OH 157.4 K, O2 23.9 K and 43.8 K) but they are derived, not read, so the choice is a stand-in. H2S-2 and the Ne-1 and
  CH3OH-1 fits above their switch are not used.
* **Below 50 K the answer is `cold-unverified`** where only a Fray and Schmitt polynomial serves. Those fits are fits to
  laboratory data that stop short of the coldest ices, and the cold sublimation is wrong by model error, not by a sign
  (docs/design/composition/atmosphere-composition.md, finding 11).
* **The laboratory-corrected cold end (GEN.497)** replaces the polynomial for nitrogen (alpha, 20 to 35 K), carbon
  monoxide (25 to 54 K), methane (30 to 51 K), oxygen (beta, 30 to 43 K) and ammonia (80 to 162 K) with the quartz crystal
  microbalance fits of [R0418, R1708], `ln P[bar] = A - B/T + C ln T + D2 T + D3 T^2 + D4 T^3`,
  up to the temperature at which they meet the Fray and Schmitt polynomial (nitrogen: the top of
  its window, where the Fray and Schmitt alpha fit is about 2.1 times higher, a step left as published). Inside its window
  the answer is `ok`. Below the window the curve is carried on by Clausius-Clapeyron with the latent heat the fit has at its
  edge (`below-fit`); that stand-in is what Boss's BC12 decision asks to be flagged. The Fray and Schmitt cold end is up to
  15 times too high for carbon monoxide near 25 K. `tuning.COLD_END_LABORATORY_CORRECTION` is the switch for the 2026
  Grundy preprint, which multiplies all five by 1.19 +/- 0.03 (the extra gauge-pressure factor of its Table 3): on, the
  answer is the 2026 set (A and B of the preprint, the other constants from the 2024 papers or fitted to them, max 0.6
  percent); off, the 2024 set (for nitrogen, carbon monoxide and methane the same curve divided by 1.19, because their 2024
  constants C to D4 are in no source read). Both sets are in the file.
* **Carbon dioxide is [R4166] equation 3.12** (`span-wagner`): `ln(P/Pt) = (Tt/T)(a1 th + a2 th^1.9 + a3
  th^2.9)`, `th = 1 - T/Tt`, good to 0.1 percent at 150 to 200 K (842.65 Pa at 150 K, 103971 Pa at 195 K, 101325 Pa at
  194.6855 K). The equation was fitted to data above 154 K and matched to the laboratory data above 85 K; below 85 K the
  value is the equation's, flagged `below-fit`. It replaces the two Fray and Schmitt fits, which are 7 percent high at 125 K
  and 41 percent at 100 K.
* **Below its `floor_k` a polynomial turns downwards** and its value means nothing (`below-fit`); the pressure returned
  is the polynomial's at the floor, which is below 1e-30 bar for every species and cannot overflow.
* **Above the triple point the liquid curve takes over** (`phase` is `liquid`): the Wagner fits to the NIST saturation
  tables of docs/design/composition/atmosphere-composition-parts/01-gas-table-thermodynamics.md section 6.2 (hydrogen cyanide by the
  Perry correlation), good to 0.3 percent of the table. The solid and liquid curves do not meet exactly at the triple
  point: the step is the "Err at Tt" of that part (SO2 +3.8 percent, benzene +7 percent, N2 -2.8 percent), accepted as
  the published spread, so a pressure inverted across the triple point can have two answers within it. Acetylene has no
  liquid curve (`above-triple`, the value is the polynomial carried on). Above the critical temperature there is no
  saturation: the status is `above-critical` and the pressure is infinite.
* **Known departures at the triple point** from the gas table's pressure: SO2 3.8 percent over, benzene 7 percent over,
  acetylene 3.6 percent under; the other species are within 3 percent (src/tests/test_saturation.py).
* **Hydrogen's solid** is Clausius-Clapeyron from its triple point with a latent heat of sublimation of 1.02 kJ/mol (the
  gas table's 0.90 of vapourisation plus 0.12 of fusion, a stand-in, unchecked against any measurement). Helium and the species
  without a curve raise `UnknownSpecies`.

The functions are pure.
"""

import json
import math
import os
from dataclasses import dataclass
from typing import Optional, Tuple

from planetgen import tuning
from planetgen.physics import elements, gases

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "saturation_pressure.json")

BAR_TO_PA = 1.0e5
"""float: Pascals in a bar; the Fray and Schmitt polynomials give `ln P` in bar."""

GAS_CONSTANT_J_MOL_K = 8.314

PRESSURE_FLOOR_PA = 1.0e-300
"""float: The least pressure the cold-end and carbon dioxide curves return, so a caller's logarithm of a very cold ice is
finite where the exponential would underflow to zero."""

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
class ColdEnd:
    """
    The laboratory-corrected cold end of a Fray and Schmitt species (GEN.497): `ln P[bar] = A - B/T + C ln T + D2 T +
    D3 T^2 + D4 T^3`.

    Attributes:
        curve_id (str): `"N2-QCM"` and so on.
        phase (str): The solid phase the fit describes.
        from_k (float): The bottom of the fit's window; below it the curve is carried on by Clausius-Clapeyron.
        corrected, uncorrected (tuple): A, B, C, D2, D3, D4 of the 2026 and the 2024 set (`tuning.COLD_END_LABORATORY_CORRECTION`).
        hand_over_corrected_k, hand_over_uncorrected_k (float): Where the set meets the Fray and Schmitt polynomial and gives way to it.
    """

    curve_id: str
    phase: str
    from_k: float
    corrected: Tuple[float, ...]
    uncorrected: Tuple[float, ...]
    hand_over_corrected_k: float
    hand_over_uncorrected_k: float

    def active(self):
        """`(coefficients, hand_over_k)` of the set the switch picks."""
        if tuning.COLD_END_LABORATORY_CORRECTION:
            return self.corrected, self.hand_over_corrected_k
        return self.uncorrected, self.hand_over_uncorrected_k


@dataclass(frozen=True)
class SpanWagner:
    """
    [R4166] equation 3.12: `ln(P/Pt) = (Tt/T)(a1 th + a2 th^1.9 + a3 th^2.9)`, `th = 1 - T/Tt`.

    Attributes:
        t_t_k, p_t_pa (float): The triple point.
        coefficients (tuple): a1, a2, a3.
        floor_k (float): Below this the equation is an extrapolation of the data it was matched to (`below-fit`).
    """

    t_t_k: float
    p_t_pa: float
    coefficients: Tuple[float, ...]
    floor_k: float


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
        form (str): The solid's: `"fray-schmitt"`, `"feistel-wagner"`, `"span-wagner"`, `"clausius-clapeyron"` or `"edmister"` (the critical-point
            correlation of a species with no triple point in the gas table, GEN.524).
        selection (str): How a second fit was placed (`"single"` or `"crossing"`).
        curves (tuple): The `Curve`s, coldest first (empty for water and carbon dioxide).
        cold_end (ColdEnd or None): The laboratory-corrected cold end that serves below the polynomial.
        span_wagner (SpanWagner or None): The carbon dioxide equation.
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
    cold_end: Optional[ColdEnd]
    span_wagner: Optional[SpanWagner]
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
        elif row["form"] not in ("feistel-wagner", "clausius-clapeyron", "span-wagner", "edmister"):
            raise ValueError(f"saturation_pressure.json: {formula} has an unknown form {row['form']!r}")
        cold_end = None
        if row.get("cold_end"):
            item = row["cold_end"]
            if row["form"] != "fray-schmitt" or any(len(item[key]) != 6 for key in ("corrected", "uncorrected")):
                raise ValueError(f"saturation_pressure.json: {formula} cold_end needs a Fray and Schmitt species and six constants")
            cold_end = ColdEnd(item["id"], item["phase"], float(item["from_k"]), tuple(float(a) for a in item["corrected"]),
                               tuple(float(a) for a in item["uncorrected"]), float(item["hand_over_k"]["corrected"]),
                               float(item["hand_over_k"]["uncorrected"]))
        span_wagner = None
        if row["form"] == "span-wagner":
            item = row["span_wagner"]
            span_wagner = SpanWagner(float(item["t_t_k"]), float(item["p_t_pa"]), tuple(float(a) for a in item["coefficients"]),
                                     float(item["floor_k"]))
        liquid = None
        if row.get("liquid"):
            item = row["liquid"]
            if item["form"] not in ("wagner", "dippr101"):
                raise ValueError(f"saturation_pressure.json: {formula} has an unknown liquid form {item['form']!r}")
            liquid = Liquid(item["form"], float(item["t_max_k"]), item.get("pc_pa"), tuple(float(a) for a in item["coefficients"]))
        rows[formula] = Species(formula, row["name"], row["form"], row["selection"], curves, cold_end, span_wagner,
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


def _lobo_ferreira_bar(coefficients, temperature_k):
    a, b, c, d2, d3, d4 = coefficients
    return math.exp(a - b / temperature_k + c * math.log(temperature_k) + d2 * temperature_k + d3 * temperature_k ** 2
                    + d4 * temperature_k ** 3)


def _lobo_ferreira_slope(coefficients, temperature_k):
    """d ln P / dT of the Lobo and Ferreira form."""
    _a, b, c, d2, d3, d4 = coefficients
    return b / temperature_k ** 2 + c / temperature_k + d2 + 2.0 * d3 * temperature_k + 3.0 * d4 * temperature_k ** 2


def _cold_end_pa(cold_end, coefficients, temperature_k):
    """The cold end's pressure, and whether `temperature_k` is under its window (then Clausius-Clapeyron from the edge)."""
    if temperature_k >= cold_end.from_k:
        return _lobo_ferreira_bar(coefficients, temperature_k) * BAR_TO_PA, False
    edge = cold_end.from_k
    latent_over_r = _lobo_ferreira_slope(coefficients, edge) * edge ** 2
    pressure = _lobo_ferreira_bar(coefficients, edge) * BAR_TO_PA * math.exp(-latent_over_r * (1.0 / temperature_k - 1.0 / edge))
    return max(pressure, PRESSURE_FLOOR_PA), True


def _span_wagner_pa(equation, temperature_k):
    theta = 1.0 - temperature_k / equation.t_t_k
    a1, a2, a3 = equation.coefficients
    pressure = equation.p_t_pa * math.exp(equation.t_t_k / temperature_k * (a1 * theta + a2 * theta ** 1.9 + a3 * theta ** 2.9))
    return max(pressure, PRESSURE_FLOOR_PA)


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


EDMISTER_COEFFICIENT = 5.373
"""float: the Edmister correlation: `ln(P/Pc) = 5.373 (1 + omega) (1 - Tc/T)`, the vapour pressure from the critical point and the acentric factor."""


def _edmister_pa(formula, temperature_k):
    row = gases.gas(formula)
    return row.critical_pressure_pa * math.exp(
        EDMISTER_COEFFICIENT * (1.0 + row.acentric_factor) * (1.0 - row.critical_temperature_k / temperature_k))


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
    if row.cold_end is not None:
        coefficients, hand_over = row.cold_end.active()
        if temperature_k < hand_over:
            pressure, under = _cold_end_pa(row.cold_end, coefficients, temperature_k)
            return Saturation(pressure, BELOW_FIT if under else OK,
                              f"{row.cold_end.curve_id}-CC" if under else row.cold_end.curve_id, SOLID)
    if row.form == "span-wagner":
        pressure = _span_wagner_pa(row.span_wagner, temperature_k)
        if triple is not None and temperature_k > triple:
            return Saturation(pressure, ABOVE_TRIPLE, "CO2-SW", SOLID)
        return Saturation(pressure, BELOW_FIT if temperature_k < row.span_wagner.floor_k else OK, "CO2-SW", SOLID)
    if row.form == "feistel-wagner":
        pressure = _water_pa(row.water_coefficients, temperature_k)
        curve_id = "H2O-FW"
    elif row.form == "clausius-clapeyron":
        pressure = _clausius_clapeyron_pa(formula, row.sublimation_kj_mol, temperature_k)
        curve_id = f"{formula}-CC"
    elif row.form == "edmister":
        pressure = _edmister_pa(formula, temperature_k)
        curve_id = f"{formula}-Edmister"
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


def condensation_temperature_k(formula, pressure_pa):
    """
    The temperature at which `formula` saturates at `pressure_pa` (a gas of mole fraction `x` in air at `P` has the partial
    pressure `x P`): the root of `pressure_pa(formula, T)`, found by bisection. `None` when the pressure is over the
    critical pressure (the gas does not condense into a liquid; hydrogen cyanide, with no critical pressure stored, is
    capped at its liquid curve's end).

    A pressure inside the step between the solid and liquid curves at the triple point (`saturation`'s docstring) has
    two roots; the one the bisection lands on is returned.
    """
    if not pressure_pa > 0:
        raise ValueError("a pressure must be positive")
    row = species(formula)
    if row.liquid is not None:
        top, top_pa = row.liquid.t_max_k, _liquid_pa(row.liquid, row.liquid.t_max_k)
    else:
        top = triple_point_k(formula)
        top_pa = saturation(formula, top).pressure_pa
    if pressure_pa > top_pa:
        return None
    low, high = 0.5, top
    for _ in range(200):
        middle = 0.5 * (low + high)
        if saturation(formula, middle).pressure_pa < pressure_pa:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)
