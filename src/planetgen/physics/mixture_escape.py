# planetgen/physics/mixture_escape.py

"""
Thermal escape of the light gases from a small body's mixture (GEN.693)
=======================================================================

A small warm world cannot keep hydrogen or helium: the gas outgassed at the surface leaves from the exobase faster than the
volcanoes supply it. The filter works on the mole fractions `mixture_draw.draw_mixture` takes from the families that outgas
(`rocky`, `venus`, `mars`, `abiotic-oxygen`); the giants, the envelopes, the lava vapour (supplied by a
magma ocean faster than it leaves), the ice-supported air and the Titan-like air are not filtered.

The Jeans parameter of a species at the exobase is `lambda = m v_esc^2 / (2 k T_exo)`, the square of the escape speed over the most
probable speed. A species is kept over the age of the body when its escape speed is at least six times its rms speed (the rule of thumb
of Jeans escape for 4.5 Gyr, Catling and Kasting 2017, chapter 5), which is `lambda = 54`; the age moves the threshold by the log of
the time (the Jeans flux falls as `exp(-lambda)`, so a body twice as old needs `ln 2` more). The share of a species that stays is a
logistic in `lambda` of width `ESCAPE_WIDTH`, so that the filter is smooth rather than a cut and no species flips from all to none
between two neighbouring bodies; a species with less than `LOST_BELOW` of it left is gone. What stays is renormalised; if less than `LOST_TOTAL` of the gas stays the body is airless.

The exobase temperature is the stand-in `EXOBASE_RATIO` times the equilibrium temperature, scaled by the star's XUV to the power
`XUV_EXPONENT`: a gas that radiates well in the infrared (carbon dioxide) keeps the thermosphere cold (Mars and Venus, 1.2 to 1.5 times
the equilibrium temperature), an air of nitrogen and oxygen does not (Earth, 4 times). Pressure is not changed: the total comes from
the pressure stage, which already carries the loss of the atmosphere (the cosmic shoreline of rows d).
"""

import dataclasses
import math

from planetgen.physics import gases
from planetgen.physics.constants import BOLTZMANN
from planetgen.physics.mixture_base import MIN_FRACTION, MixtureError

AVOGADRO = 6.02214076e23

FAMILIES = ("rocky", "venus", "mars", "abiotic-oxygen")
"""tuple: The families that are filtered. Not `titan` (its hydrogen trace is the steady state of photochemical production and escape,
measured by Niemann et al. 2010, which a loss-only filter would remove) and not `frost` (an air held by the vapour pressure of its ice,
which the ice resupplies faster than it leaves)."""

EXOBASE_RATIO = {"cool": 1.3, "hot": 4.0}
"""dict: The exobase temperature over the equilibrium temperature: `cool` when carbon dioxide is half the gas or more (Mars about 1.2,
Venus about 1.4), `hot` otherwise (Earth about 4 at the mean of the solar cycle); stand-ins fitted to the Solar System, replaced when the
thermosphere model of GEN.475's neighbours exists (NO FIXED CONTENT: they are ratios, not the planets)."""

XUV_EXPONENT = 0.25
"""float: The exobase temperature goes as the XUV flux (in Earth's) to this power, clamped to 0.5 .. 3.0 times; a stand-in."""

LAMBDA_KEPT = 54.0
"""float: `lambda` at which half of a species is kept for 4.5 Gyr: the escape speed is six times the rms speed, `6**2 * 3/2`."""

ESCAPE_WIDTH = 4.0
"""float: The width of the logistic in `lambda` (a factor `e` in the share kept for 4 in `lambda`); design."""

LOST_BELOW = 1e-3
"""float: A species of which less than this share stays is gone (not left as a trace that renormalising would promote to the main gas)."""

LOST_TOTAL = 1e-2
"""float: A body whose gas has less than this share left after the loss is airless: what stays of a mostly escaped air is not an air."""

REFERENCE_AGE_GY = 4.5


def exobase_temperature_k(subject, mole_fractions):
    """The exobase temperature of the body, K (see the module text)."""
    teq = subject.get("teq_k") or subject.get("ts_k")
    if not teq:
        raise MixtureError("the escape filter needs the equilibrium temperature (teq_k) or the surface temperature")
    co2 = mole_fractions.get("CO2", 0.0)
    ratio = EXOBASE_RATIO["cool" if co2 >= 0.5 else "hot"]
    xuv = max(subject.get("xuv_index") or 1.0, 1e-3)
    return teq * ratio * min(3.0, max(0.5, xuv ** XUV_EXPONENT))


def jeans_parameter(formula, v_esc_km_s, exobase_k):
    """`lambda = m v_esc^2 / (2 k T)` for the species `formula` (gas-table name)."""
    mass = gases.gas(formula).molar_mass / AVOGADRO
    return mass * (v_esc_km_s * 1e3) ** 2 / (2.0 * BOLTZMANN * exobase_k)


def kept_share(parameter, age_gy):
    """The share of a species with Jeans parameter `parameter` that a body of age `age_gy` keeps, 0..1."""
    threshold = LAMBDA_KEPT + math.log(max(age_gy, 1e-3) / REFERENCE_AGE_GY)
    x = (parameter - threshold) / ESCAPE_WIDTH
    return 1.0 / (1.0 + math.exp(-x)) if x > -700.0 else 0.0


def apply(subject, result):
    """
    `result` (a `MixtureDraw`) without the gases the body loses.

    Args:
        subject (Subject): With `v_esc_kms`, `age_gy` and `teq_k` (else `ts_k`); `xuv_index` is optional (1).
        result (MixtureDraw): The family's mixture.

    Returns:
        MixtureDraw: The mixture renormalised, with the share each species kept in `constants["escape_kept"]`; an empty mixture
        (airless) if nothing stays.
    """
    if result.family not in FAMILIES or not result.mole_fractions:
        return result
    v_esc = subject.get("v_esc_kms")
    age = subject.get("age_gy")
    if v_esc is None or age is None:
        raise MixtureError("the escape filter needs the quantities v_esc_kms and age_gy")
    exobase = exobase_temperature_k(subject, result.mole_fractions)
    kept = {formula: kept_share(jeans_parameter(formula, v_esc, exobase), age) for formula in result.mole_fractions}
    left = {formula: fraction * kept[formula] for formula, fraction in result.mole_fractions.items() if kept[formula] >= LOST_BELOW}
    left = {formula: value for formula, value in left.items() if value > MIN_FRACTION}
    total = sum(left.values())
    if total < LOST_TOTAL:
        total = 0.0
    constants = {**result.constants, "exobase_k": exobase, "escape_kept": {f: round(k, 6) for f, k in kept.items() if k < 0.999999}}
    if total <= 0.0:
        return dataclasses.replace(result, mole_fractions={}, family="airless", constants=constants)
    return dataclasses.replace(result, mole_fractions={formula: value / total for formula, value in left.items()}, constants=constants)
