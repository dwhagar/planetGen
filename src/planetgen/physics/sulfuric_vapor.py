# planetgen/physics/sulfuric_vapor.py

"""
The vapor pressure of sulfuric acid and water over aqueous sulfuric acid, and the window in which the liquid condenses (GEN.1424)
================================================================================================================================

`sulfuric_cloud` and `acid_medium` take the acid vapor as made in the air and stand the cloud base at a fixed temperature, because the
curve that says where sulfuric acid and water condense together was not reached. This module is that curve. It does not change a stored
value; `GEN.1425` (the cloud step of the Venus chemistry set) reads it.

* **The solution.** Zeleznik's Gibbs-energy model of aqueous sulfuric acid ([R4877], 200 to 350 K, 0.1 MPa) gives the chemical
  potentials of acid and water against the acid mole fraction `x` and the temperature (`ln_acid_activity`, `ln_water_activity`, equations 12
  and 13 with the parameters of its Table 6, checked against its Table 7). The activity of a component is its chemical potential against the
  pure liquid at the same temperature, so the vapor pressure over the solution is the pure liquid's times the activity.
* **Water.** The pure liquid curve is the liquid curve of `saturation` at and above the triple point and [R4289] equation 1 (supercooled
  water, 185 to 260 K, carried to the triple point) below it. A water pressure in the air gives the water activity `p / p0`, and the model
  gives the acid mole fraction at that activity (`composition`): the weight percent of an aerosol in equilibrium with the air.
* **Acid.** The pure liquid curve is the Ayers form with the [R2516] low-temperature correction as three open papers
  print it [R0976, R1516, R5168]; the primary paper was not opened, so the tier is T3. The refit constant of Noppel and others
  2002 is the other end of the bracket (1.28 times higher at every temperature), reported by `condensation_temperature_k(..., alternative=True)`.
* **The window.** An air condenses the liquid where the acid vapor pressure it holds is above the pressure over the solution that its water
  fixes (`condenses`); the warmest temperature at which it does is `condensation_temperature_k`, and the cloud lies where the air is cooler
  than that. Above 350 K the solution model is extrapolated and the status says so (`EXTRAPOLATED`).

The functions are pure. The numbers are in `physics/data/sulfuric_vapor.json` with their tier and tweak level.
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen.physics import constants, elements, saturation

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "sulfuric_vapor.json")

OK = "ok"
EXTRAPOLATED = "extrapolated"
"""str: The status of a temperature over the 350 K top of the solution model (or under its 200 K bottom): the value is the model carried on."""
WATER_SATURATED = "water-saturated"
"""str: The air holds as much water as the pure liquid can, so the solution is pure water."""

ATMOSPHERE_PA = constants.STANDARD_ATMOSPHERE_PA
"""float: Pascals in a standard atmosphere; the pure-acid curve is written in atmospheres."""

BISECTIONS = 80
MOLE_FRACTION_FLOOR = 1.0e-15
"""float: The least acid mole fraction the logarithm is taken at (the acid potential is minus infinity at none)."""


class SulfuricVaporError(ValueError):
    """A rule file that is refused."""


def _check(name, entry):
    if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) \
            or not entry.get("basis"):
        raise SulfuricVaporError(f"{name}: every entry needs a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        _check(name, entry)
        if "value" not in entry:
            raise SulfuricVaporError(f"{name}: no value")
    _check("zeleznik_1991", document["zeleznik_1991"])
    for name, coefficients in document["zeleznik_1991"]["parameters"].items():
        if len(coefficients) != 5:
            raise SulfuricVaporError(f"{name}: five coefficients (1, T, T^2, 1/T, ln T) are needed")
    return document


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]
VERSION = DOCUMENT["version"]
_ZELEZNIK = DOCUMENT["zeleznik_1991"]
_PARAMETERS = _ZELEZNIK["parameters"]
MODEL_FROM_K = _ZELEZNIK["valid_from_k"]
MODEL_TO_K = _ZELEZNIK["valid_to_k"]


def number(name):
    """The value of the rule file's number `name`."""
    return NUMBERS[name]["value"]


def _polynomial(coefficients, temperature_k):
    return coefficients[0] + coefficients[1] * temperature_k + coefficients[2] * temperature_k ** 2 + coefficients[3] / temperature_k \
        + coefficients[4] * math.log(temperature_k)


def _terms(temperature_k):
    """The (j, k, i, mu, eps) terms of Zeleznik's equation 12, 0-based indices (acid is 0, water is 1)."""
    mu = {name: _polynomial(_PARAMETERS[name], temperature_k) for name in _PARAMETERS if name.startswith("mu_")}
    eps = {name: _polynomial(_PARAMETERS[name], temperature_k) for name in _PARAMETERS if name.startswith("eps_")}
    return (
        (0, 0, 0, mu["mu_111"], eps["eps_111"]),
        (0, 1, 0, mu["mu_121"], eps["eps_121"]),
        (1, 0, 0, mu["mu_121"], eps["eps_211"]),
        (1, 1, 0, mu["mu_221"], eps["eps_221"]),
        (0, 1, 1, mu["mu_122"], eps["eps_122"]),
        (1, 0, 1, mu["mu_122"], eps["eps_212"]),
    )


def _clamped(acid_fraction):
    """The acid mole fraction held off 0 and 1 by `MOLE_FRACTION_FLOOR` (the logarithms are taken of both components)."""
    return min(max(acid_fraction, MOLE_FRACTION_FLOOR), 1.0 - MOLE_FRACTION_FLOOR)


def relative_gibbs(acid_fraction, temperature_k):
    """`(G/RT, d(G/RT)/dx)`, the relative molar Gibbs energy of the solution (Zeleznik equation 12, sign restored) and its slope in the acid mole
    fraction `x`, against the reference of the paper (not the pure liquid at the temperature; see `ln_acid_activity`)."""
    x = _clamped(acid_fraction)
    fractions = (x, 1.0 - x)
    slopes = (1.0, -1.0)
    value = slope = 0.0
    for j, k, i, mu, eps in _terms(temperature_k):
        phi, dphi = (1.0, 0.0) if i == 0 else (x * (1.0 - x), 1.0 - 2.0 * x)
        xj, xk = fractions[j], fractions[k]
        product = xj * xk
        dproduct = slopes[j] * xk + xj * slopes[k]
        inner = mu + eps * math.log(xj)
        dinner = eps * slopes[j] / xj
        value += phi * inner * product
        slope += dphi * inner * product + phi * (dinner * product + inner * dproduct)
    return -value, -slope


def relative_potentials(acid_fraction, temperature_k):
    """`(mu_acid/RT, mu_water/RT)`, the relative chemical potentials of Zeleznik's Table 7 (against the paper's reference, not the pure liquids)."""
    g, dg = relative_gibbs(acid_fraction, temperature_k)
    x = _clamped(acid_fraction)
    return g + (1.0 - x) * dg, g - x * dg


def ln_acid_activity(acid_fraction, temperature_k):
    """ln of the activity of sulfuric acid in the solution, against pure liquid acid at the same temperature."""
    return relative_potentials(acid_fraction, temperature_k)[0] - relative_potentials(1.0, temperature_k)[0]


def ln_water_activity(acid_fraction, temperature_k):
    """ln of the activity of water in the solution, against pure liquid water (supercooled below the freezing point) at the same temperature."""
    return relative_potentials(acid_fraction, temperature_k)[1] - relative_potentials(0.0, temperature_k)[1]


def weight_percent(acid_fraction):
    """The weight percent of sulfuric acid of a solution of acid mole fraction `acid_fraction`."""
    acid = acid_fraction * number("acid_molar_mass_kg_mol")
    water = (1.0 - acid_fraction) * number("water_molar_mass_kg_mol")
    return 100.0 * acid / (acid + water)


def mole_fraction(weight_percent_acid):
    """The acid mole fraction of a solution of `weight_percent_acid` percent sulfuric acid."""
    acid = weight_percent_acid / number("acid_molar_mass_kg_mol")
    water = (100.0 - weight_percent_acid) / number("water_molar_mass_kg_mol")
    return acid / (acid + water)


def _status(temperature_k):
    return OK if MODEL_FROM_K <= temperature_k <= MODEL_TO_K else EXTRAPOLATED


def water_pressure_pa(temperature_k):
    """The vapor pressure of pure liquid water, Pa: [R4289] equation 1 (supercooled) below the triple point, the liquid curve of `saturation` above."""
    if not temperature_k > 0:
        raise ValueError("a temperature must be positive")
    if temperature_k >= saturation.WATER_TRIPLE_POINT_K:
        return saturation.pressure_pa("H2O", temperature_k)
    ln_mbar = number("water_curve_a") - number("water_curve_b") / temperature_k - number("water_curve_c") / temperature_k ** 2 \
        + number("water_curve_d") / temperature_k ** 3
    return 100.0 * math.exp(ln_mbar)


def pure_acid_pressure_pa(temperature_k, alternative=False):
    """The vapor pressure of pure liquid sulfuric acid, Pa (the Ayers curve with the Kulmala and Laaksonen correction). `alternative` takes the refit
    constant of Noppel and others 2002, the other end of the bracket."""
    if not temperature_k > 0:
        raise ValueError("a temperature must be positive")
    exponent = number("pure_acid_exponent_k")
    if alternative:
        reference, constant = number("pure_acid_alternative_reference_k"), number("pure_acid_alternative_constant")
    else:
        reference = number("pure_acid_reference_k")
        constant = number("pure_acid_constant") - exponent / reference
    correction = exponent * (1.0 / reference - 1.0 / temperature_k + number("pure_acid_curvature")
                             / (number("pure_acid_critical_k") - reference)
                             * (1.0 + math.log(reference / temperature_k) - reference / temperature_k))
    return ATMOSPHERE_PA * math.exp(constant + correction)


@dataclass(frozen=True)
class Solution:
    """The solution in equilibrium with an air's water: its water activity, acid mole fraction, weight percent and the status of the model there."""

    water_activity: float
    acid_fraction: float
    weight_percent: float
    status: str


def composition(temperature_k, water_pressure_in_air_pa):
    """
    The solution in equilibrium with an air holding `water_pressure_in_air_pa` of water at `temperature_k`.

    The water activity is the air's water pressure over the pure liquid's. At or above 1 the air is saturated with water and the solution is pure water
    (`WATER_SATURATED`, 0 percent acid). An air with no water gives the pure acid. The acid mole fraction is where Zeleznik's water activity equals
    that (it falls the whole way from 1 at none to 0 at pure acid, so there is one answer).
    """
    activity = water_pressure_in_air_pa / water_pressure_pa(temperature_k)
    if activity >= 1.0:
        return Solution(1.0, 0.0, 0.0, WATER_SATURATED)
    if activity <= 0.0:
        return Solution(0.0, 1.0, 100.0, _status(temperature_k))
    target = math.log(activity)
    low, high = 0.0, 1.0 - MOLE_FRACTION_FLOOR
    for _ in range(BISECTIONS):
        middle = 0.5 * (low + high)
        if ln_water_activity(middle, temperature_k) > target:
            low = middle
        else:
            high = middle
    x = 0.5 * (low + high)
    return Solution(activity, x, weight_percent(x), _status(temperature_k))


def acid_pressure_over_solution_pa(temperature_k, water_pressure_in_air_pa, alternative=False):
    """The sulfuric acid vapor pressure over the solution that `water_pressure_in_air_pa` of water fixes at `temperature_k`, Pa: the pure liquid's times the
    acid activity. A pure-water solution holds none."""
    solution = composition(temperature_k, water_pressure_in_air_pa)
    if solution.status == WATER_SATURATED:
        return 0.0
    return pure_acid_pressure_pa(temperature_k, alternative) * math.exp(ln_acid_activity(solution.acid_fraction, temperature_k))


def condenses(temperature_k, acid_pressure_in_air_pa, water_pressure_in_air_pa, alternative=False):
    """Whether sulfuric acid and water condense together at `temperature_k`: the air holds more acid vapor than the solution its water fixes."""
    return acid_pressure_in_air_pa > acid_pressure_over_solution_pa(temperature_k, water_pressure_in_air_pa, alternative)


def condensation_temperature_k(acid_pressure_in_air_pa, water_pressure_in_air_pa, alternative=False, low_k=MODEL_FROM_K - 15.0, high_k=500.0):
    """
    The warmest temperature at which an air holding the two pressures condenses the liquid, K, or `None` when it does not condense anywhere from
    `low_k` to `high_k` (no acid vapor, or too little for even the coldest). The two pressures are those of the air (the partial pressures stay as the
    air is cooled, so the caller passes mixing ratio times the pressure of the level). Above 350 K the result is the model carried on; check
    `status_at` for that.
    """
    if acid_pressure_in_air_pa <= 0.0:
        return None
    if not condenses(low_k, acid_pressure_in_air_pa, water_pressure_in_air_pa, alternative):
        return None
    if condenses(high_k, acid_pressure_in_air_pa, water_pressure_in_air_pa, alternative):
        return high_k
    low, high = low_k, high_k
    for _ in range(BISECTIONS):
        middle = 0.5 * (low + high)
        if condenses(middle, acid_pressure_in_air_pa, water_pressure_in_air_pa, alternative):
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def status_at(temperature_k):
    """`OK` inside the 200 to 350 K of Zeleznik's model, `EXTRAPOLATED` outside it."""
    return _status(temperature_k)


def window(levels, acid_ratio, water_ratio, alternative=False):
    """
    The levels of an air column where the liquid condenses: `levels` is a sequence of `(pressure_pa, temperature_k)`, `acid_ratio` and `water_ratio` the mole
    mixing ratios of sulfuric acid vapor and water (taken as the same at every level, above a cloud they would be depleted; the caller says what its
    column holds). Returns a list of booleans, one per level, and is how the cloud step reads the cloud base and top.
    """
    return [condenses(temperature, acid_ratio * pressure, water_ratio * pressure, alternative) for pressure, temperature in levels]
