# planetgen/physics/stability.py

"""
Hierarchical stability (GEN.261): the tests a nested system of stars and
compact objects has to pass, as pure functions of masses (solar masses),
axes (AU), eccentricities and angles (radians, except where a name says
`_deg`).

- **Triple stability** -- one function, `is_triple_stable`, with a switch:
  Mardling and Aarseth 2001 (`"ma01"`, for planets and wherever the mass
  ratios leave the fitted range of the newer fit), Eggleton and Kiseleva
  1995 (`"ek95"`, a cross-check) and Vynatheya et al. 2022, Equation 4
  (`"eq4"`, for star-only triples inside its fitted range). `"auto"` picks
  Equation 4 inside its range, Mardling and Aarseth outside it where their
  `q_out <= 5` holds, and Eggleton and Kiseleva beyond that (GEN.1028). The
  forms and their unit-test values are in
  docs/design/systems/exotic-environments-planets-and-compact-binaries.md 2.3 and
  docs/design/systems/multistar-and-compact-systems.md 5.1.
- **Kozai-Lidov screen** -- `kozai_verdict`: the quadrupole maximum
  eccentricity, the Kozai time against the age, the octupole widening for
  eccentric outer orbits, and the relativistic quench (a cycle slower than the
  inner orbit's general-relativity precession never happens).
- **Peters merger time** -- `peters_merger_time_yr`: when a compact pair
  merges by gravitational waves.

Notation: a triple is an inner pair of masses `m1 >= m2` (`q_in = m2 / m1`)
on an orbit `(a_in, e_in)` and an outer body `m3` on `(a_out, e_out)`, with
`q_out = m3 / (m1 + m2)` and mutual inclination `i`.
"""

import math

from planetgen.physics import constants, unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

_m_to_au = units.converter("m", "au")
_au_to_m = units.converter("au", "m")
_s_to_yr = units.converter("s", "yr")
_SOLAR_GM_OVER_C2_AU = _m_to_au(_msun_to_kg(constants.G) / constants.SPEED_OF_LIGHT_M_S ** 2)
"""float: Half the Schwarzschild radius of the Sun, `G M / c^2`, in AU (9.87e-9)."""

_PETERS_COEFFICIENT_S = (5.0 / 256.0) * constants.SPEED_OF_LIGHT_M_S ** 5 / constants.G ** 3
"""float: The constant of the circular merger time `t = C a^4 / (m1 m2 (m1 + m2))` in SI units."""

METHODS = ("auto", "ma01", "ek95", "eq4")
"""tuple: The criteria `is_triple_stable` accepts."""

MA01_Q_OUT_MAX = 5.0
"""float: The largest `q_out` Mardling and Aarseth 2001 give their criterion for (Equation 90: "holds for q_out <= 5"; beyond it the failure is
exchange, not escape). `required_axis_ratio` refuses `"ma01"` past it and `"auto"` takes `"ek95"` there (GEN.1028, registry row GM12)."""

KOZAI_WINDOW_COS2 = 3.0 / 5.0
"""float: Kozai-Lidov needs `cos^2 i < 3/5`: a mutual inclination of 39.2 to 140.8 degrees."""


def orbital_period_yr(a_au, total_mass_solar):
    """Kepler's third law: the period in years of an orbit of axis `a_au` about `total_mass_solar`."""
    return math.sqrt(a_au ** 3 / total_mass_solar)


def axis_au(period_yr, total_mass_solar):
    """The semi-major axis in AU of an orbit of `period_yr` about `total_mass_solar`."""
    return (total_mass_solar * period_yr ** 2) ** (1.0 / 3.0)


def mutual_inclination(inc1, node1, inc2, node2):
    """The angle (radians) between two orbital planes given each one's inclination and ascending node."""
    cos_i = math.cos(inc1) * math.cos(inc2) + math.sin(inc1) * math.sin(inc2) * math.cos(node1 - node2)
    return math.acos(max(-1.0, min(1.0, cos_i)))


def _q_out(m1, m2, m3):
    return m3 / (m1 + m2)


def mardling_aarseth_pericentre_ratio(q_out, e_out, i):
    """Mardling and Aarseth 2001 (Equation 2 of Vynatheya et al.): the smallest `a_out (1 - e_out) / a_in` that is stable."""
    return 2.8 * ((1.0 + q_out) * (1.0 + e_out) / math.sqrt(1.0 - e_out)) ** 0.4 * (1.0 - 0.3 * i / math.pi)


def mardling_aarseth_axis_ratio(q_out, e_out, i):
    """The smallest stable `a_out / a_in` by Mardling and Aarseth."""
    return mardling_aarseth_pericentre_ratio(q_out, e_out, i) / (1.0 - e_out)


def eggleton_kiseleva_ycrit(q_out, q_in):
    """Eggleton and Kiseleva 1995: the critical `Y = a_out (1 - e_out) / (a_in (1 + e_in))`, with `q = q_out`."""
    q13 = q_out ** (1.0 / 3.0)
    return 1.0 + 3.7 / q13 - 2.2 / (1.0 + q13) + 1.4 / q_in ** (1.0 / 3.0) * (q13 - 1.0) / (q13 + 1.0)


def effective_inner_eccentricity(e_in, i):
    """Equation 3 of Vynatheya et al.: `max(e_in, e_in_avg)` with the Kozai-Lidov time-averaged eccentricity.
    Outside the Kozai window (`cos^2 i > 3/5`) the maximum is taken as 0, so the result is `e_in`."""
    window = 1.0 - 5.0 / 3.0 * math.cos(i) ** 2
    e_max = math.sqrt(window) if window > 0.0 else 0.0
    return max(e_in, 0.5 * e_max ** 2)


def eq4_ycrit(q_out, e_in, e_out, i):
    """Equation 4 of Vynatheya et al. 2022: the critical `Y~` (with the effective inner eccentricity)."""
    e_eff = effective_inner_eccentricity(e_in, i)
    base = 2.4 * ((1.0 + q_out) / ((1.0 + e_eff) * math.sqrt(1.0 - e_out))) ** 0.4
    return base * (((1.0 - 0.2 * e_eff + e_out) / 8.0) * (math.cos(i) - 1.0) + 1.0)


def eq4_axis_ratio(q_out, e_in, e_out, i):
    """The smallest stable `a_out / a_in` by Equation 4."""
    e_eff = effective_inner_eccentricity(e_in, i)
    return eq4_ycrit(q_out, e_in, e_out, i) * (1.0 + e_eff) / (1.0 - e_out)


def eq4_in_range(q_in, q_out, alpha):
    """Whether the masses and axis ratio are inside the range Equation 4 was fitted on."""
    return 1e-2 <= q_in <= 1.0 and 1e-2 <= q_out <= 1e2 and 1e-4 < alpha < 1.0


def ma01_in_range(q_out):
    """Whether Mardling and Aarseth's criterion is given for this outer mass ratio (`MA01_Q_OUT_MAX`)."""
    return q_out <= MA01_Q_OUT_MAX


def criterion(m1, m2, m3, alpha=None):
    """
    The criterion `"auto"` takes for a triple: `"eq4"` where `eq4_in_range` holds (`alpha`, `a_in / a_out`, defaults to a mid-range value
    when not given), else `"ma01"` where `ma01_in_range` holds (a star-sized inner pair with a tiny or small outer body), else `"ek95"`
    (a planet pair about a star, a star pair about a black hole: the outer body outweighs the pair past what either fit covers).
    """
    heavy, light = max(m1, m2), min(m1, m2)
    q_out = _q_out(m1, m2, m3)
    if eq4_in_range(light / heavy, q_out, 0.1 if alpha is None else alpha):
        return "eq4"
    return "ma01" if ma01_in_range(q_out) else "ek95"


def required_axis_ratio(m1, m2, m3, e_in, e_out, i, method="auto", alpha=None):
    """
    The smallest stable `a_out / a_in` for the triple by `method` (`METHODS`).

    `"auto"` is the one `criterion` picks. `"ma01"` asked for by name is refused past `MA01_Q_OUT_MAX` (GEN.1028).
    """
    if method not in METHODS:
        raise ValueError(f"unknown stability method {method!r}; choose from {METHODS}")
    heavy, light = max(m1, m2), min(m1, m2)
    q_in = light / heavy
    q_out = _q_out(m1, m2, m3)
    if method == "auto":
        method = criterion(m1, m2, m3, alpha)
    if method == "ma01":
        if not ma01_in_range(q_out):
            raise ValueError(f"Mardling and Aarseth 2001 holds for q_out <= {MA01_Q_OUT_MAX:g}, got {q_out:.4g}; use method='auto'")
        return mardling_aarseth_axis_ratio(q_out, e_out, i)
    if method == "eq4":
        return eq4_axis_ratio(q_out, e_in, e_out, i)
    return eggleton_kiseleva_ycrit(q_out, q_in) * (1.0 + e_in) / (1.0 - e_out)


def is_triple_stable(m1, m2, m3, a_in, a_out, e_in, e_out, i, method="auto", safety=1.0):
    """Whether the triple is stable: `a_out / a_in` exceeds `required_axis_ratio` times `safety`."""
    return a_out / a_in > safety * required_axis_ratio(m1, m2, m3, e_in, e_out, i, method, alpha=a_in / a_out)


def max_inner_axis(m1, m2, m3, a_out, e_in, e_out, i, method="auto", safety=1.0):
    """The largest inner axis (AU) the triple allows: `a_out` over the required ratio, over `safety`.
    The ratio depends on `alpha = a_in / a_out` only through the choice of fit, so the answer is
    refined once with the alpha it implies."""
    ratio = required_axis_ratio(m1, m2, m3, e_in, e_out, i, method)
    a_in = a_out / (ratio * safety)
    return a_out / (required_axis_ratio(m1, m2, m3, e_in, e_out, i, method, alpha=a_in / a_out) * safety)


def kozai_timescale_yr(p_in_yr, p_out_yr, total_mass, m3, e_out):
    """The Kozai-Lidov time `(2 / 3 pi) (P_out^2 / P_in) (M / m3) (1 - e_out^2)^1.5`, years."""
    return (2.0 / (3.0 * math.pi)) * (p_out_yr ** 2 / p_in_yr) * (total_mass / m3) * (1.0 - e_out ** 2) ** 1.5


def gr_precession_period_yr(a_in, e_in, m_in, p_in_yr):
    """The period of the inner orbit's general-relativity apsidal precession, years:
    `P_in c^2 a (1 - e^2) / (3 G M)`."""
    return p_in_yr * a_in * (1.0 - e_in ** 2) / (3.0 * _SOLAR_GM_OVER_C2_AU * m_in)


def kozai_max_eccentricity(e_in, i):
    """The quadrupole-order maximum inner eccentricity, or `e_in` outside the Kozai window."""
    x = (5.0 / 3.0) * (1.0 - e_in ** 2) * math.cos(i) ** 2
    return math.sqrt(1.0 - x) if x < 1.0 else e_in


def kozai_verdict(m1, m2, m3, a_in, a_out, e_in, e_out, i, age_yr, extent_sum_au):
    """
    The Kozai-Lidov screen for a triple, with the relativistic quench.

    Returns `"safe"`, `"kl_active"` (cycles run but the pair keeps clear of
    its stars' extent: the pair is kept and marked) or `"unsafe"` (the
    cycle brings the pericentre inside four times `extent_sum_au`, a
    merger). A cycle needs all of: a mutual inclination in the Kozai window
    (39.2 to 140.8 degrees), a Kozai time shorter than the age, and a Kozai
    time shorter than the inner orbit's relativistic precession period. The
    octupole widening for an eccentric outer orbit lifts the maximum to
    0.98 for inclinations of 25 to 90 degrees (and the mirror).
    """
    m_in = m1 + m2
    x = (5.0 / 3.0) * (1.0 - e_in ** 2) * math.cos(i) ** 2
    mutual = min(i, math.pi - i)
    eps = (a_in / a_out) * e_out / (1.0 - e_out ** 2)
    octupole = eps > 0.01 and math.radians(25.0) < mutual < math.radians(90.0)
    if x >= 1.0 and not octupole:
        return "safe"
    e_max = kozai_max_eccentricity(e_in, i)
    if octupole:
        e_max = max(e_max, 0.98)
    p_in = orbital_period_yr(a_in, m_in)
    p_out = orbital_period_yr(a_out, m_in + m3)
    t_kl = kozai_timescale_yr(p_in, p_out, m_in + m3, m3, e_out)
    if t_kl >= age_yr or t_kl >= gr_precession_period_yr(a_in, e_in, m_in, p_in):
        return "safe"
    return "unsafe" if a_in * (1.0 - e_max) < 4.0 * extent_sum_au else "kl_active"


def peters_merger_time_yr(m1, m2, a_au, e=0.0):
    """
    The gravitational-wave merger time of a pair, in years (Peters 1964):
    `(5/256) c^5 a^4 / (G^3 m1 m2 (m1 + m2))` for a circular orbit, times
    `(1 - e^2)^3.5` for an eccentric one (the approximate form, good to a
    factor of a few at high eccentricity).
    """
    seconds = _PETERS_COEFFICIENT_S * _au_to_m(a_au) ** 4 / (
        (_msun_to_kg(m1)) * (_msun_to_kg(m2))
        * (_msun_to_kg(m1 + m2)))
    return _s_to_yr(seconds * (1.0 - e ** 2) ** 3.5)


def survives_gravitational_waves(m1, m2, a_au, e, age_yr):
    """Whether a compact pair has not yet merged: its Peters time exceeds the system's age."""
    return peters_merger_time_yr(m1, m2, a_au, e) > age_yr
