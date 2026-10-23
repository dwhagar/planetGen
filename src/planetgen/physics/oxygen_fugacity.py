# planetgen/physics/oxygen_fugacity.py

"""
One oxygen-fugacity scale for the generator (GEN.498).

A rock's oxygen fugacity is quoted in log units above a buffer, and the buffers differ: Sossi et al. (2020) put
Earth's magma ocean at IW+0.5, and the same gas ratios are IW+0.0 on the iron-wustite calibration of Hirschmann
(2021); the Frost (1991) formula is 0.003 to 0.024 log units from it at 1 bar; QFM stands 3.1 to 3.6 log units above
IW. Mixing values from different papers without converting shifts the gas ratios by a factor of three.

* **The reference is the Hirschmann (2021, Geochim. Cosmochim. Acta 313, 74-84, doi 10.1016/j.gca.2021.08.039)
  iron-wustite buffer** (`iw_log_fo2`, 1000 to 3000 K, 1e-4 to 100 GPa, root mean square 0.00647 log units, largest
  mismatch 0.0283). `delta_iw` everywhere in the generator, the stored `mantle_delta_iw` included, is
  `log10 fO2 - iw_log_fo2(T, P)` on it. Every coefficient and its sign, and the iron-phase boundary, were checked
  against the paper's Table 1.
* **Other buffers** are converted by `buffer_offset` (Frost 1991 fits of FMQ and NNO, whose own iron-wustite is
  `frost_log_fo2("IW", ...)`) and other authors' scales by `convert_delta_iw`; the `sossi2020` scale is the
  reference minus a constant 0.5 at 2173 K, 1 bar, where the gas ratios of that paper were computed. Sossi's IW rests
  on the NIST-JANAF free energy of FeO, which Hirschmann shows puts IW 0.2 to 1.1 log units too reducing (0.4 to 0.7
  at 2000 K): the offset is consistent with that, and its constant is derived from the gas ratios, not read.
* **What is verified** is in `SOURCES`: Hirschmann's Table 1, the NASA Glenn coefficients and Sossi's equation 4
  against the papers' text; the Frost IW against Hirschmann (0.003 to 0.024 log units at 1 bar from 1500 to 2200 K).
  **Stand-ins:** the Frost FMQ and NNO constants (the Frost 1991 text is not available; the form
  `log fO2 = A/T + B + C (P - 1)/T` is confirmed by Anenburg 2026, equation 4, which also finds the FMQ fit offset from
  the thermodynamic value with a wrong pressure slope, and the NNO one good at 1 bar). The constants are the
  high-temperature set (above 846 K) of the PetroBuffer code, checked only by the gap they leave between NNO and FMQ
  at 1 GPa (0.15 to 0.37 over 1000 to 1873 K against Anenburg's "only 0.3"). Above about 1 GPa the FMQ and NNO
  offsets are not reliable; the Hirschmann IW is.
* **Gas ratios** (`co2_to_co`, `h2o_to_h2`) follow from the oxygen fugacity and the reactions
  `CO + 1/2 O2 = CO2` and `H2 + 1/2 O2 = H2O` (ideal gas, 1 bar standard state, the NASA Glenn fits of
  `data/nasa_glenn_redox_gases.json`): `log10(f_CO2/f_CO) = log K1 + 0.5 log10 fO2`. At 2173 K the CO2/CO ratio of 1
  is at IW+1.37 and the H2O/H2 ratio of 1 at IW-0.06 (docs/design/body-composition-research.md section 15 and the
  G3 research).
* **Melt iron** (`fe3_fraction`, `delta_iw_from_fe3_fraction`): Sossi et al. (2020) equation 4,
  `log10(X_FeO1.5 / X_FeO) = 0.252 dIW - 1.530`, peridotite liquid, 1 bar, 2173 K, dIW on the Sossi scale, measured
  from IW-1.5 to IW+6.4.

The functions are pure and take temperatures in kelvin and pressures in gigapascals (`P_GPa`) except the Frost form,
which is in bar as published.
"""

import json
import math
import os

from planetgen.physics import elements

REFERENCE_SCALE = "hirschmann2021"
"""str: The scale every stored delta-IW is on."""

IW_CALIBRATED_K = (1000.0, 3000.0)
IW_CALIBRATED_GPA = (1e-4, 100.0)
"""tuple: The range the Hirschmann (2021) fit is calibrated over; `iw_log_fo2` extrapolates smoothly outside it."""

SOSSI_OFFSET_2173K = 0.5
"""float: The Sossi (2020) scale minus the reference at 2173 K and 1 bar, from the gas ratios (CO2/CO 0.21 and H2O/H2 1.05 at IW+0.5 there equal the reference's at IW+0.0)."""

SOSSI_T_K = 2173.0
SOSSI_FE_SLOPE = 0.252
SOSSI_FE_INTERCEPT = -1.530
"""float: Sossi et al. (2020) equation 4: log10(FeO1.5/FeO) = SLOPE dIW + INTERCEPT."""

# Hirschmann (2021) Table 1: m = m0 + m1 P + m2 P^2 + m3 P^3 + m4 sqrt(P), P in GPa, for the parameters a b c d
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

FROST_BUFFERS = {
    "IW": (-27489.0, 6.702, 0.055),
    "FMQ": (-25096.3, 8.735, 0.110),
    "NNO": (-24930.0, 9.360, 0.046),
}
"""dict: Frost (1991, Rev. Mineral. 25) buffers, `log10 fO2 = A/T + B + C (P - 1)/T` with P in bar. FMQ is the quartz-fayalite-magnetite buffer above 846 K (the same as QFM; Frost gives another set below). Constants second hand, from PetroBuffer 0.0.1; IW is cross-checked against Hirschmann, FMQ and NNO are stand-ins (`SOURCES`)."""

IW_FIT_RMS = 0.00647
IW_FIT_MAX_MISMATCH = 0.0283
"""float: Hirschmann (2021) Table 1 footnote: the fit's root mean square and largest mismatch against the thermodynamic calculation, log units."""

SOURCES = (
    ("hirschmann2021 IW fit: coefficients, signs, phase boundary, range", "Hirschmann 2021 Table 1", "VERIFIED", "the paper's text"),
    ("NASA Glenn fits of CO, CO2, H2, H2O, O2", "McBride, Zehe and Gordon 2002, NASA/TP-2002-211556, Appendix D", "VERIFIED", "the report's text, every interval and coefficient"),
    ("Sossi 2020 equation 4: melt Fe3+/Fe2+ against delta-IW at 2173 K, 1 bar", "Sossi et al. 2020, Sci. Adv. 6, eabd1387", "VERIFIED", "the paper's text"),
    ("Sossi 2020 gas ratios at his IW+0.5 (H2O/H2 1.05, CO2/CO 0.21)", "Sossi et al. 2020", "VERIFIED", "the paper's text; the NASA fits at reference IW+0.0 give 4 percent and 2 percent from them"),
    ("Sossi scale 0.5 below the reference", "derived from the gas ratios above", "DERIVED", "consistent with the JANAF offset of Hirschmann 2021 (0.4 to 0.7 at 2000 K); Sossi's own IW formula is in his Supplementary Materials, not read"),
    ("Frost 1991 IW constants", "Frost 1991, Rev. Mineral. 25 (text not available); second hand from PetroBuffer", "CROSS-CHECKED", "agrees with the Hirschmann IW to 0.003 to 0.024 log units at 1 bar, 1500 to 2200 K"),
    ("Frost 1991 FMQ and NNO constants", "Frost 1991 (text not available); second hand from PetroBuffer", "STAND-IN", "the form is Anenburg 2026 equation 4; NNO minus FMQ at 1 GPa is 0.15 to 0.37 against Anenburg's 0.3; the absolute values are not checked"),
)
"""tuple: One row per constant set: what it is, where it comes from, how far it is checked, and against what."""

SCALES = (REFERENCE_SCALE, "frost1991", "sossi2020")
"""tuple: The scales `convert_delta_iw` knows."""

GAS_REACTIONS = {
    "co2_co": {"CO2": 1.0, "CO": -1.0, "O2": -0.5},
    "h2o_h2": {"H2O": 1.0, "H2": -1.0, "O2": -0.5},
}
"""dict: The reactions, products positive."""

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "nasa_glenn_redox_gases.json")


def _poly(coefficients, pressure_gpa):
    m0, m1, m2, m3, m4 = coefficients
    return m0 + m1 * pressure_gpa + m2 * pressure_gpa ** 2 + m3 * pressure_gpa ** 3 + m4 * math.sqrt(pressure_gpa)


def iw_log_fo2(temperature_k, pressure_gpa=1e-4):
    """
    log10 of the oxygen fugacity (bar) of the iron-wustite buffer, Hirschmann (2021): the reference scale.

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


def frost_log_fo2(buffer, temperature_k, pressure_bar=1.0):
    """log10 of the oxygen fugacity (bar) of `buffer` ("IW", "FMQ" or "NNO") by Frost (1991)."""
    a, b, c = FROST_BUFFERS[buffer]
    return a / temperature_k + b + c * (pressure_bar - 1.0) / temperature_k


def log_fo2(delta_iw, temperature_k, pressure_gpa=1e-4):
    """log10 fO2 (bar) at `delta_iw` log units above the reference iron-wustite buffer."""
    return iw_log_fo2(temperature_k, pressure_gpa) + delta_iw


def delta_iw(log_fo2_bar, temperature_k, pressure_gpa=1e-4):
    """The log units a log10 fO2 (bar) lies above the reference iron-wustite buffer."""
    return log_fo2_bar - iw_log_fo2(temperature_k, pressure_gpa)


def buffer_offset(buffer, temperature_k, pressure_gpa=1e-4):
    """
    How far a Frost (1991) `buffer` ("FMQ", "NNO", "IW") lies above the reference iron-wustite buffer, log units:
    3.1 to 3.6 for FMQ at 1 bar between 2200 and 1500 K.
    """
    return frost_log_fo2(buffer, temperature_k, pressure_gpa * 1e4) - iw_log_fo2(temperature_k, pressure_gpa)


def convert_delta_iw(value, from_scale, to_scale, temperature_k=SOSSI_T_K, pressure_gpa=1e-4):
    """
    A delta-IW from one scale to another at a temperature and pressure.

    `frost1991` is log units above the Frost iron-wustite fit, `sossi2020` above the Sossi et al. IW (the reference
    minus `SOSSI_OFFSET_2173K` where it was calibrated: the offset is applied at every temperature, as the paper's
    own scale is not published as a function of temperature).
    """
    for scale in (from_scale, to_scale):
        if scale not in SCALES:
            raise ValueError(f"unknown scale {scale!r}; the scales are {', '.join(SCALES)}")
    return value + _scale_origin(from_scale, temperature_k, pressure_gpa) - _scale_origin(to_scale, temperature_k, pressure_gpa)


def _scale_origin(scale, temperature_k, pressure_gpa):
    """log10 fO2 of a scale's zero, relative to the reference zero (so a value on the scale is `origin + value` above the reference)."""
    if scale == REFERENCE_SCALE:
        return 0.0
    if scale == "frost1991":
        return frost_log_fo2("IW", temperature_k, pressure_gpa * 1e4) - iw_log_fo2(temperature_k, pressure_gpa)
    return -SOSSI_OFFSET_2173K  # sossi2020: its zero lies 0.5 below the reference's, so IW+0.5 there is IW+0.0 here


def _load_gases():
    with open(DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)["gases"]


_GASES = _load_gases()


def _g_over_rt(gas, temperature_k):
    for fit in _GASES[gas]:
        if fit["t_min_k"] <= temperature_k <= fit["t_max_k"]:
            break
    else:
        raise ValueError(f"{gas}: no thermodynamic fit at {temperature_k} K (200 to 20000 K; H2O to 6000 K)")
    a1, a2, a3, a4, a5, a6, a7 = fit["a"]
    b1, b2 = fit["b"]
    t = temperature_k
    h_rt = -a1 / t ** 2 + a2 * math.log(t) / t + a3 + a4 * t / 2 + a5 * t ** 2 / 3 + a6 * t ** 3 / 4 + a7 * t ** 4 / 5 + b1 / t
    s_r = -a1 / (2 * t ** 2) - a2 / t + a3 * math.log(t) + a4 * t + a5 * t ** 2 / 2 + a6 * t ** 3 / 3 + a7 * t ** 4 / 4 + b2
    return h_rt - s_r


def log_k(reaction, temperature_k):
    """log10 of the equilibrium constant (bar to the net power) of a `GAS_REACTIONS` key at `temperature_k`."""
    stoichiometry = GAS_REACTIONS[reaction]
    return -math.fsum(n * _g_over_rt(gas, temperature_k) for gas, n in stoichiometry.items()) / math.log(10.0)


def co2_to_co(delta_iw_value, temperature_k, pressure_gpa=1e-4):
    """f_CO2 / f_CO in gas in equilibrium with a melt at `delta_iw_value` on the reference scale (ideal gas)."""
    return 10.0 ** (log_k("co2_co", temperature_k) + 0.5 * log_fo2(delta_iw_value, temperature_k, pressure_gpa))


def h2o_to_h2(delta_iw_value, temperature_k, pressure_gpa=1e-4):
    """f_H2O / f_H2 in gas in equilibrium with a melt at `delta_iw_value` on the reference scale (ideal gas)."""
    return 10.0 ** (log_k("h2o_h2", temperature_k) + 0.5 * log_fo2(delta_iw_value, temperature_k, pressure_gpa))


def delta_iw_where_equal(reaction, temperature_k, pressure_gpa=1e-4):
    """The delta-IW at which the two gases of `reaction` ("co2_co", "h2o_h2") have equal fugacity."""
    return -2.0 * log_k(reaction, temperature_k) - iw_log_fo2(temperature_k, pressure_gpa)


def fe3_fraction(sossi_delta_iw):
    """Fe3+/Fe(total) of a peridotite liquid at 1 bar and 2173 K at `sossi_delta_iw` (Sossi et al. 2020 equation 4)."""
    ratio = 10.0 ** (SOSSI_FE_SLOPE * sossi_delta_iw + SOSSI_FE_INTERCEPT)
    return ratio / (1.0 + ratio)


def delta_iw_from_fe3_fraction(fraction):
    """The inverse of `fe3_fraction`: delta-IW on the Sossi scale from Fe3+/Fe(total), 0 < fraction < 1."""
    if not 0.0 < fraction < 1.0:
        raise ValueError("Fe3+/Fe(total) must lie between 0 and 1")
    return (math.log10(fraction / (1.0 - fraction)) - SOSSI_FE_INTERCEPT) / SOSSI_FE_SLOPE
