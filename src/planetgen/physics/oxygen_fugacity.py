# planetgen/physics/oxygen_fugacity.py

"""
One oxygen-fugacity scale for the generator (GEN.498).

A rock's oxygen fugacity is quoted in log units above a buffer, and the buffers differ: [R4154] put
Earth's magma ocean at IW+0.5, and the same gas ratios are IW+0.0 on the iron-wustite calibration of Hirschmann
(2021). Mixing values from different papers without converting shifts the gas ratios by a factor of three.

* **The reference is the [R1951] (Geochim. Cosmochim. Acta 313, 74-84, [R1951])
  iron-wustite buffer** (`iw_log_fo2`, 1000 to 3000 K, 1e-4 to 100 GPa, root mean square 0.00647 log units, largest
  mismatch 0.0283). `delta_iw` everywhere in the generator, the stored `mantle_delta_iw` included, is
  `log10 fO2 - iw_log_fo2(T, P)` on it. Every coefficient and its sign, and the iron-phase boundary, were checked
  against the paper's Table 1. A value quoted on the [R4154] scale is 0.5 log units above the same state here
  at 2173 K and 1 bar (derived from the gas ratios of that paper; consistent with the JANAF offset Hirschmann
  finds in Sossi's iron-wustite); convert it before storing it.
* **What is verified** is in `SOURCES`: Hirschmann's Table 1, the NASA Glenn coefficients and Sossi's equation 4
  against the papers' text.
* **Gas ratios** follow from the oxygen fugacity and the reactions `CO + 1/2 O2 = CO2` and `H2 + 1/2 O2 = H2O`
  (ideal gas, 1 bar standard state, the NASA Glenn fits of `physics/species_thermo.py`). At 2173 K the CO2/CO ratio
  of 1 is at IW+1.37 and the H2O/H2 ratio of 1 at IW-0.06 (docs/design/composition/body-composition-research.md section 15 and
  the G3 research). The mixture draw takes its gas ratios from the equilibrium of `species_thermo`, so no function here
  turns a delta-IW into one (the buffer conversions of [R1464] and the gas-ratio helpers were removed by TEST.524:
  nothing called them).
* **Melt iron** (`fe3_fraction`): [R4154] equation 4,
  `log10(X_FeO1.5 / X_FeO) = 0.252 dIW - 1.530`, peridotite liquid, 1 bar, 2173 K, dIW on the Sossi scale, measured
  from IW-1.5 to IW+6.4.

The functions are pure and take temperatures in kelvin and pressures in gigapascals (`P_GPa`).
"""

import math

IW_CALIBRATED_K = (1000.0, 3000.0)
IW_CALIBRATED_GPA = (1e-4, 100.0)
"""tuple: The range the [R1951] fit is calibrated over; `iw_log_fo2` extrapolates smoothly outside it."""

SOSSI_FE_SLOPE = 0.252
SOSSI_FE_INTERCEPT = -1.530
"""float: [R4154] equation 4: log10(FeO1.5/FeO) = SLOPE dIW + INTERCEPT."""

# [R1951: Table 1]: m = m0 + m1 P + m2 P^2 + m3 P^3 + m4 sqrt(P), P in GPa, for the parameters a b c d
# (face-centred cubic and body-centred cubic iron) and e f g h (hcp iron).
_FCC = {
    "a": (6.844864, 1.175691e-1, 1.143873e-3, 0.0, 0.0),
    "b": (5.791364e-4, -2.891434e-4, -2.737171e-7, 0.0, 0.0),
    "c": (-7.971469e-5, 3.198005e-5, 0.0, 1.059554e-10, 2.014461e-7),
    "d": (-2.769002e4, 5.285977e2, -2.919275, 0.0, 0.0),
}
_HCP = {
    "a": (8.463095, -3.000307e-3, 7.213445e-5, 0.0, 0.0),
    "b": (1.148738e-3, -9.352312e-5, 5.161592e-7, 0.0, 0.0),
    "c": (-7.448624e-4, -6.329325e-6, 0.0, -1.407339e-10, 1.830014e-4),
    "d": (-2.782082e4, 5.285977e2, -8.473231e-1, 0.0, 0.0),
}
_HCP_BOUNDARY = (-18.64, 0.04359, -5.069e-6)

IW_FIT_RMS = 0.00647
IW_FIT_MAX_MISMATCH = 0.0283
"""float: [R1951: Table 1] footnote: the fit's root mean square and largest mismatch against the thermodynamic calculation, log units."""

SOURCES = (
    ("hirschmann2021 IW fit: coefficients, signs, phase boundary, range", "[R1951: Table 1]", "VERIFIED", "the paper's text"),
    ("NASA Glenn fits of CO, CO2, H2, H2O, O2", "[R2968], NASA/TP-2002-211556, Appendix D", "VERIFIED", "the report's text, every interval and coefficient"),
    ("[R4154] equation 4: melt Fe3+/Fe2+ against delta-IW at 2173 K, 1 bar", "[R4154], Sci. Adv. 6, eabd1387", "VERIFIED", "the paper's text"),
)
"""tuple: One row per constant set: what it is, where it comes from, how far it is checked, and against what."""


def _poly(coefficients, pressure_gpa):
    m0, m1, m2, m3, m4 = coefficients
    return m0 + m1 * pressure_gpa + m2 * pressure_gpa ** 2 + m3 * pressure_gpa ** 3 + m4 * math.sqrt(pressure_gpa)


def iw_log_fo2(temperature_k, pressure_gpa=1e-4):
    """
    log10 of the oxygen fugacity (bar) of the iron-wustite buffer, [R1951]: the reference scale.

    Args:
        temperature_k (float): Kelvin, calibrated 1000 to 3000.
        pressure_gpa (float): Gigapascals, calibrated 1e-4 (1 bar) to 100.
    """
    if not temperature_k > 0 or not pressure_gpa >= 0:
        raise ValueError("a temperature must be positive and a pressure not negative")
    boundary = _HCP_BOUNDARY[0] + _HCP_BOUNDARY[1] * temperature_k + _HCP_BOUNDARY[2] * temperature_k ** 2
    p = _HCP if pressure_gpa > boundary else _FCC
    return (
        _poly(p["a"], pressure_gpa)
        + _poly(p["b"], pressure_gpa) * temperature_k
        + _poly(p["c"], pressure_gpa) * temperature_k * math.log(temperature_k)
        + _poly(p["d"], pressure_gpa) / temperature_k
    )


def log_fo2(delta_iw, temperature_k, pressure_gpa=1e-4):
    """log10 fO2 (bar) at `delta_iw` log units above the reference iron-wustite buffer."""
    return iw_log_fo2(temperature_k, pressure_gpa) + delta_iw


def delta_iw(log_fo2_bar, temperature_k, pressure_gpa=1e-4):
    """The log units a log10 fO2 (bar) lies above the reference iron-wustite buffer."""
    return log_fo2_bar - iw_log_fo2(temperature_k, pressure_gpa)


def fe3_fraction(sossi_delta_iw):
    """Fe3+/Fe(total) of a peridotite liquid at 1 bar and 2173 K at `sossi_delta_iw` ([R4154] equation 4)."""
    ratio = 10.0 ** (SOSSI_FE_SLOPE * sossi_delta_iw + SOSSI_FE_INTERCEPT)
    return ratio / (1.0 + ratio)


