# planetgen/physics/co2_far_wing.py

"""
The far wings of carbon dioxide (GEN.1019, docs/design/composition/line-opacity.md; replaces the flat wing of stand-in IR3)
===========================================================================================================================

The ExoMolOP k-tables cut every line at 25 cm-1 from its centre, so a dense carbon dioxide air (Venus: 92 bar) is too clear between its
bands. This module is the opacity of what the cutoff leaves out: for every band `j` of strength `S_j`, the Lorentz tail of its lines beyond
the cutoff, times the sub-Lorentzian factor `chi` of Perrin and Hartmann 1989 (JQSRT 42, 311), out to a reach.

    kappa_i = sum_j S_j  < (1/pi) [x_c gamma_c chi_cc(d) + (1 - x_c) gamma_f chi_cn(d)] / (gamma^2 + d^2) >    for cutoff <= d <= reach

with `d` the distance between a point of band `i` and a point of band `j` (averaged over both bands, a band's strength spread evenly over its width), `gamma = x_c gamma_c + (1 - x_c) gamma_f` the Lorentz half width at the pressure and
temperature, and `chi` the four-piece exponential of Simonetti and others 2021 (their equations 12 to 15):
`chi = 1` to `sigma1`, then `exp(-b1 (d - sigma1))`, `exp(-b1 (sigma2 - sigma1) - b2 (d - sigma2))`, and beyond `sigma3` the same with `b3`;
`b_i(T) = alpha_i + beta_i exp(-epsilon_i T)`, held at 190 to 800 K, the temperatures it was measured at.

The band strengths `S_j` are not stored: the caller passes them (`line_opacity` reads them off the carbon dioxide k-table at its lowest
pressure, where the lines are narrow and the band mean times the band width is the sum of the line strengths). Every other number is in
`physics/data/co2_far_wing.json` with its source tier and tweak level; the registry rows are FW1 to FW5 of the gap plan.

Limits (what the wing is not): the factor is applied to every band, not only the 4.3 micron one it was fitted to; the tables' own line
shape inside 25 cm-1 is not made sub-Lorentzian (a few percent of a band's absorption at the pressures that matter); the width is the
carbon dioxide and nitrogen one for every gas; the temperature exponent and the widths are recalled, not read.

Units: strengths in cm per molecule (cm-1 times cm2), result in cm2 per molecule of carbon dioxide; wavenumbers in cm-1; pressure in Pa.
"""

import json
import os

import numpy as np

from planetgen.physics import elements
from planetgen.util import detmath

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "co2_far_wing.json")

ATMOSPHERE_PA = 101325.0
"""float: One standard atmosphere, Pa, the pressure the half widths are stated at."""
REFERENCE_K = 296.0
"""float: The temperature the half widths are stated at (HITRAN's)."""
REFERENCE_PRESSURE_BAR = 1.0e-5
"""float: The pressure of the k-table the band strengths are read at, bar: its lowest, where the lines are narrow and the band mean times the
width is the sum of the line strengths."""
GRID_STEP_CM1 = 0.5
"""float: The step of the grid the tail is integrated on, cm-1."""


class FarWingError(ValueError):
    """A rule file that is refused."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) \
                or not entry.get("basis") or "value" not in entry:
            raise FarWingError(f"{name}: every number needs a value, a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")
    return document


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]
VERSION = DOCUMENT["version"]


def number(name):
    """The value of the rule file's number `name`."""
    return NUMBERS[name]["value"]


def chi_parameters(pair):
    """The Perrin and Hartmann parameters of `pair` (`"co2_co2"` or `"co2_n2"`) as arrays: `(sigma, alpha, beta, epsilon)`, three of each."""
    entry = number(f"chi_{pair}")
    return tuple(np.array(entry[key], dtype=float) for key in ("sigma_cm1", "alpha", "beta", "epsilon"))


def chi(distance_cm1, temperature_k, pair):
    """
    The sub-Lorentzian factor at `distance_cm1` from a line centre for `pair`, at `temperature_k` (held at the measured range, `b_i` at least 0
    so the factor never exceeds 1). `distance_cm1` and `temperature_k` broadcast.
    """
    sigma, alpha, beta, epsilon = chi_parameters(pair)
    temperature = np.clip(np.asarray(temperature_k, dtype=float), number("temperature_low_k"), number("temperature_high_k"))
    b = [np.maximum(alpha[i] + beta[i] * np.exp(-epsilon[i] * temperature), 0.0) for i in range(3)]
    distance = np.asarray(distance_cm1, dtype=float)
    first = np.clip(distance - sigma[0], 0.0, sigma[1] - sigma[0])
    second = np.clip(distance - sigma[1], 0.0, sigma[2] - sigma[1])
    third = np.maximum(distance - sigma[2], 0.0)
    return np.exp(-(b[0] * first + b[1] * second + b[2] * third))


def half_width_cm1(temperature_k, pressure_pa, co2_fraction):
    """The Lorentz half width, cm-1, of a carbon dioxide line in air with this share of carbon dioxide."""
    scale = (np.asarray(pressure_pa, dtype=float) / ATMOSPHERE_PA) * (REFERENCE_K / np.asarray(temperature_k, dtype=float)) ** number("gamma_exponent")
    return scale * (co2_fraction * number("gamma_self_cm1_atm") + (1.0 - co2_fraction) * number("gamma_foreign_cm1_atm"))


def band_strengths(k_mean_cm2, edges_cm1):
    """
    The line strength in each band, cm per molecule: the band-mean cross section times the band width.

    Args:
        k_mean_cm2 (ndarray): The mean cross section of each band (..., bands), cm2 per molecule, at low pressure.
        edges_cm1 (ndarray): The band edges, cm-1 (bands + 1).
    """
    return np.asarray(k_mean_cm2, dtype=float) * np.diff(edges_cm1)


def _second_integral(distance, tail, step):
    """`F2(u)` on the grid `distance`: the integral from 0 of `F1(v) dv`, `F1` the integral from 0 of the tail `f` (trapezoids)."""
    first = np.concatenate([[0.0], np.cumsum(0.5 * (tail[1:] + tail[:-1]) * step)])
    return np.concatenate([[0.0], np.cumsum(0.5 * (first[1:] + first[:-1]) * step)])


def _pair_matrix(edges_cm1, tail_of_distance):
    """
    The weight of every band pair, shape (band i, band j): the strength of band `j` spread evenly over its width and the wing averaged evenly
    over band `i`, so the weight is the double integral of the tail `f(d)` over the two intervals, `F2(b - c) - F2(a - c) - F2(b - d) + F2(a - d)`
    (`F2` the second integral of `f`, even) over both widths. It does not depend on how the bands are sampled.

    Args:
        edges_cm1 (ndarray): The band edges.
        tail_of_distance (callable): `f(distance_cm1)`, the line-shape tail per cm-1 on an array of distances (0 inside the cutoff and past the reach).
    """
    edges = np.asarray(edges_cm1, dtype=float)
    width = np.diff(edges)
    step = GRID_STEP_CM1
    distance = np.arange(0.0, edges[-1] - edges[0] + 2.0 * step, step)
    f2 = _second_integral(distance, tail_of_distance(distance), step)
    lower, upper = edges[:-1], edges[1:]

    def at(u):
        return np.interp(np.abs(u), distance, f2)

    pair = at(upper[:, None] - lower[None, :]) - at(lower[:, None] - lower[None, :]) \
        - at(upper[:, None] - upper[None, :]) + at(lower[:, None] - upper[None, :])
    return pair / (width[:, None] * width[None, :])


def _window(distance):
    return (distance >= number("cutoff_cm1")) & (distance <= number("reach_cm1"))


def cross_section_cm2(strengths, edges_cm1, temperature_k, pressure_pa, co2_fraction):
    """
    The far-wing cross section of carbon dioxide, cm2 per molecule, in each band, for a stack of layers.

    Args:
        strengths (ndarray): Each layer's band strengths (layers, bands), cm per molecule (`band_strengths`).
        edges_cm1 (ndarray): The band edges, cm-1.
        temperature_k, pressure_pa, co2_fraction (ndarray): One value per layer: temperature, total pressure, carbon dioxide mole fraction.

    Returns:
        ndarray: (layers, bands).
    """
    strengths = np.atleast_2d(np.asarray(strengths, dtype=float))
    temperature = np.atleast_1d(np.asarray(temperature_k, dtype=float))
    pressure = np.atleast_1d(np.asarray(pressure_pa, dtype=float))
    fraction = np.atleast_1d(np.asarray(co2_fraction, dtype=float))
    result = np.zeros(strengths.shape)
    for layer in range(len(temperature)):
        if fraction[layer] <= 0.0 or pressure[layer] <= 0.0:
            continue
        x, t = fraction[layer], temperature[layer]
        gamma = half_width_cm1(t, pressure[layer], x)
        scale = (pressure[layer] / ATMOSPHERE_PA) * (REFERENCE_K / t) ** number("gamma_exponent")

        def tail(distance):
            weight = x * number("gamma_self_cm1_atm") * chi(distance, t, "co2_co2") \
                + (1.0 - x) * number("gamma_foreign_cm1_atm") * chi(distance, t, "co2_n2")
            return np.where(_window(distance), scale * weight / (np.pi * (gamma ** 2 + distance ** 2)), 0.0)

        result[layer] = detmath.matmul(_pair_matrix(edges_cm1, tail), strengths[layer][:, None])[:, 0]
    return np.maximum(result, 0.0)


def minor_gas_cross_section_cm2(gas, strengths, edges_cm1, temperature_k, pressure_pa, co2_fraction):
    """
    The far-wing cross section, cm2 per molecule, of a minor absorber (`"H2O"`, `"SO2"`, `"CO"`) broadened by a carbon dioxide air, in each
    band: its Lorentz tail beyond the cutoff with the half width `minor_gamma_cm1_atm` and the carbon dioxide self factor (or none, when
    `minor_uses_chi` is 0), weighted by the carbon dioxide share of the air so that an air without carbon dioxide adds nothing (its water
    far wing is the continuum of `line_opacity`). Arguments as `cross_section_cm2`.
    """
    strengths = np.atleast_2d(np.asarray(strengths, dtype=float))
    temperature = np.atleast_1d(np.asarray(temperature_k, dtype=float))
    pressure = np.atleast_1d(np.asarray(pressure_pa, dtype=float))
    fraction = np.atleast_1d(np.asarray(co2_fraction, dtype=float))
    gamma0 = number("minor_gamma_cm1_atm")[gas]
    result = np.zeros(strengths.shape)
    for layer in range(len(temperature)):
        if fraction[layer] <= 0.0 or pressure[layer] <= 0.0:
            continue
        t = temperature[layer]
        scale = (pressure[layer] / ATMOSPHERE_PA) * (REFERENCE_K / t) ** number("gamma_exponent")
        gamma = gamma0 * scale

        def tail(distance):
            shape = chi(distance, t, "co2_co2") if number("minor_uses_chi") else 1.0
            return np.where(_window(distance), fraction[layer] * gamma * shape / (np.pi * (gamma ** 2 + distance ** 2)), 0.0)

        result[layer] = detmath.matmul(_pair_matrix(edges_cm1, tail), strengths[layer][:, None])[:, 0]
    return np.maximum(result, 0.0)
