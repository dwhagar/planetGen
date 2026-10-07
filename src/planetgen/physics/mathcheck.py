# planetgen/physics/mathcheck.py

"""
Math Check: the Gate in Front of Everything Else
==================================================

A fixed list of checks that the generator's math gives the right answers,
run before anything trusts it (TEST.63): pytest runs it first and stops
the suite if it fails (`tests/test_math_check.py`, `conftest.py`), and the
website runs it once at startup and warns admins if it fails
(`planetgen/web/__init__.py`).

Each `Check` names the function it calls, the value it expects, the
tolerance it allows and where that expected value comes from (a textbook
figure, a paper, or an exact identity). There are three groups:

- **reference** (TEST.64): known answers from real astronomy -- the Sun,
  Earth's and Jupiter's orbits, Earth's Hill sphere, the habitable zone
  and snow line, white dwarf sizes, the Sun's Schwarzschild radius and
  galactic orbit, Holman & Wiegert's stability limits, and the Kepler and
  Barker equations against known solutions.
- **invariant** (TEST.65): things true for any input -- unit conversions
  round-trip, the constants agree with each other, luminosity rises and
  lifetime falls with mass, a Kepler orbit conserves energy and angular
  momentum, the sector grid tiles each ring exactly and finds every cell
  again from its own center, density is 1 at its calibration point, and
  nothing returns NaN or infinity over a fixed sweep.
- **distribution** (TEST.66): a few thousand seeded draws from each
  sampler (initial mass function, star ages, Poisson sector counts, the
  bounded bell, the planet class table) land on their intended shares by
  a chi-square test, so a broken sampler fails even when every single
  value looks fine.

Pure: no database, network or files; every random draw is seeded, and the
global `random` module's state is put back afterwards. The whole run takes
well under 5 seconds. `python -m planetgen.physics.mathcheck` prints the
report and exits 1 if anything failed.
"""

import logging
import math
import random
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Callable

from planetgen.galaxy import galactic_orbit
from planetgen.physics import formation
from planetgen.physics import orbits
from planetgen.physics import units
from planetgen.util import random as sampling
from planetgen.galaxy import density as galaxyDensity, geometry, sector as spaceSector
from planetgen.physics import constants as pc, kepler, stellar_evolution
from planetgen import tuning
from planetgen.util import log
from planetgen.physics import planets

SEED = 20261001
"""int: The seed every distribution check starts from (each check adds its
own offset), so a run is the same every time."""

CHI_SQUARE_P_MIN = 1e-4
"""float: A distribution check fails when its chi-square p-value drops
below this. The draws are seeded, so a pass or fail is the same every run;
this only decides how far off a sampler may be before it counts as
broken (a fair sampler lands below 1e-4 one seed in ten thousand)."""


@dataclass(frozen=True)
class Check:
    """
    One check.

    Attributes:
        name (str): Short unique name (`snake_case`), shown in reports.
        group (str): "reference", "invariant" or "distribution".
        function (str): The function or functions under test.
        compute (callable): No-argument callable returning the actual value.
        expected (float): What `compute` should return.
        tolerance (float): How far off it may be.
        source (str): Where `expected` comes from.
        mode (str): How `actual` is compared with `expected`:
            "rel" -- `|actual - expected| <= tolerance * |expected|`;
            "abs" -- `|actual - expected| <= tolerance`;
            "max" -- `actual <= expected` (an error bound; `tolerance` unused);
            "min" -- `actual >= expected` (a p-value floor; `tolerance` unused).
        unit (str): Unit of `actual`/`expected`, for the report.
    """
    name: str
    group: str
    function: str
    compute: Callable[[], float]
    expected: float
    tolerance: float
    source: str
    mode: str = "rel"
    unit: str = ""


@dataclass
class Result:
    """The outcome of one `Check`."""
    check: Check
    passed: bool
    actual: float = math.nan
    error: str = ""
    seconds: float = 0.0
    detail: str = field(default="")

    @property
    def name(self):
        return self.check.name

    def describe(self):
        """One line: what was expected, what came back, and why that is a
        pass or a failure."""
        c = self.check
        unit = f" {c.unit}" if c.unit else ""
        if self.error:
            outcome = f"raised {self.error}"
        else:
            outcome = f"got {self.actual:.6g}{unit}"
        if c.mode == "rel":
            want = f"{c.expected:.6g}{unit} within {c.tolerance:.3g} (relative)"
        elif c.mode == "abs":
            want = f"{c.expected:.6g}{unit} within {c.tolerance:.3g}{unit}"
        elif c.mode == "max":
            want = f"at most {c.expected:.3g}{unit}"
        else:
            want = f"at least {c.expected:.3g}{unit}"
        status = "ok  " if self.passed else "FAIL"
        return f"{status} {c.name}: {c.function} -- expected {want}, {outcome} [{c.source}]"


def _compare(check, actual):
    if not isinstance(actual, (int, float)) or math.isnan(actual):
        return False
    if check.mode == "rel":
        return abs(actual - check.expected) <= check.tolerance * abs(check.expected)
    if check.mode == "abs":
        return abs(actual - check.expected) <= check.tolerance
    if check.mode == "max":
        return actual <= check.expected
    if check.mode == "min":
        return actual >= check.expected
    raise ValueError(f"unknown comparison mode {check.mode!r}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rel_err(a, b):
    """Relative difference of `a` from `b` (absolute when `b` is 0)."""
    return abs(a - b) / abs(b) if b else abs(a - b)


def _sweep_max(pairs):
    """The largest relative error over `(actual, expected)` pairs."""
    return max(_rel_err(a, b) for a, b in pairs)


def _log_sweep(low, high, count):
    """`count` points spaced evenly in log between `low` and `high`."""
    ratio = (high / low) ** (1 / (count - 1))
    return [low * ratio ** i for i in range(count)]


def _regularized_gamma_q(a, x):
    """Q(a, x) = 1 - P(a, x), the upper regularized incomplete gamma
    function (Numerical Recipes 6.2: series below a + 1, continued
    fraction above)."""
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:
        term = total = 1.0 / a
        ap = a
        for _ in range(1000):
            ap += 1
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return max(0.0, 1.0 - total * math.exp(-x + a * math.log(x) - gln))
    tiny = 1e-300
    b = x + 1 - a
    c = 1 / tiny
    d = 1 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def chi_square_p_value(observed, expected_shares):
    """
    Pearson's chi-square goodness-of-fit p-value for counts `observed`
    against `expected_shares` (any positive weights, normalized here).
    Bins expecting fewer than 5 draws are merged into their neighbor first,
    so the test stays valid for rare categories.

    Returns:
        float: The p-value (1.0 for a perfect fit, near 0 for a bad one).
    """
    n = sum(observed)
    total_share = sum(expected_shares)
    bins = []
    obs_acc = exp_acc = 0.0
    for obs, share in zip(observed, expected_shares):
        obs_acc += obs
        exp_acc += n * share / total_share
        if exp_acc >= 5:
            bins.append((obs_acc, exp_acc))
            obs_acc = exp_acc = 0.0
    if exp_acc > 0 or obs_acc > 0:
        if bins:
            last_obs, last_exp = bins.pop()
            bins.append((last_obs + obs_acc, last_exp + exp_acc))
        else:
            bins.append((obs_acc, exp_acc))
    if len(bins) < 2:
        return 1.0
    statistic = sum((o - e) ** 2 / e for o, e in bins)
    return _regularized_gamma_q((len(bins) - 1) / 2, statistic / 2)


def _bin_counts(values, edges):
    """How many of `values` fall in each `[edges[i], edges[i+1])` (the
    last bin closed)."""
    counts = [0] * (len(edges) - 1)
    for v in values:
        for i in range(len(edges) - 1):
            if edges[i] <= v < edges[i + 1] or (i == len(edges) - 2 and v == edges[-1]):
                counts[i] += 1
                break
    return counts


@contextmanager
def _seeded_global_random(seed):
    """
    Seeds the global `random` module for samplers that draw from it
    (`sampling.sample_bounded_bell`, `planets._choose_weighted_planet_class`),
    then puts its state back. The debug log's per-draw tracing and DEBUG
    lines are switched off meanwhile, so a check run never floods the log
    with thousands of rolls.
    """
    state = random.getstate()
    tracing = bool(log._original_random_functions)
    if tracing:
        log._trace_random(False)
    previous_disable = logging.root.manager.disable
    logging.disable(logging.DEBUG)
    try:
        random.seed(seed)
        yield
    finally:
        logging.disable(previous_disable)
        if tracing:
            log._trace_random(True)
        random.setstate(state)


# ---------------------------------------------------------------------------
# TEST.64: reference values
# ---------------------------------------------------------------------------

def _sun_teff_from_si():
    """Stefan-Boltzmann with SI constants: T = (L / (4 pi R^2 sigma))^(1/4)."""
    return (pc.SOLAR_LUMINOSITY / (4 * math.pi * pc.SOLAR_RADIUS_M ** 2
                                   * pc.STEFAN_BOLTZMANN_CONSTANT)) ** 0.25


def _schwarzschild_radius_km(mass_solar):
    """r_s = 2 G M / c^2, the formula compactRemnant.BlackHole uses."""
    return 2 * pc.G * mass_solar * pc.SOLAR_MASS_TO_KG / pc.SPEED_OF_LIGHT_M_S ** 2 / 1000


def _reference_checks():
    sun_kg = pc.SOLAR_MASS_TO_KG
    period = planets.calculate_orbital_period_years
    gc_ly = pc.GALACTIC_CENTER_DISTANCE_LY
    return [
        Check("sun_luminosity", "reference", "stellar_evolution.main_sequence_luminosity_sol(1)",
              lambda: stellar_evolution.main_sequence_luminosity_sol(1.0), 1.0, 1e-9,
              "definition: 1 M_sun on the main sequence gives 1 L_sun", unit="L_sun"),
        Check("sun_main_sequence_lifetime", "reference", "stellar_evolution.main_sequence_lifetime_gy(1)",
              lambda: stellar_evolution.main_sequence_lifetime_gy(1.0), 10.0, 0.1,
              "the Sun's main-sequence lifetime, about 10 Gy (Sackmann et al. 1993)", unit="Gy"),
        Check("sun_effective_temperature", "reference", "stellar_evolution.effective_temperature_k(1, 1)",
              lambda: stellar_evolution.effective_temperature_k(1.0, 1.0), 5772.0, 1e-3,
              "IAU 2015 Resolution B3: nominal solar T_eff 5,772 K", unit="K"),
        Check("sun_effective_temperature_si", "reference",
              "Stefan-Boltzmann from SOLAR_LUMINOSITY, SOLAR_RADIUS_M, STEFAN_BOLTZMANN_CONSTANT",
              _sun_teff_from_si, 5772.0, 0.005,
              "IAU 2015 Resolution B3: nominal solar T_eff 5,772 K", unit="K"),
        Check("earth_orbital_period", "reference", "planets.calculate_orbital_period_years(1 AU, 1 M_sun)",
              lambda: period(1.0, sun_kg), 1.0, 1e-6,
              "Kepler's third law: 1 AU around 1 M_sun is 1 year", unit="yr"),
        Check("jupiter_orbital_period", "reference", "planets.calculate_orbital_period_years(5.2026 AU, 1 M_sun)",
              lambda: period(5.2026, sun_kg), 11.862, 0.002,
              "Jupiter: a = 5.2026 AU, sidereal period 11.862 yr (NASA planetary fact sheet)", unit="yr"),
        Check("earth_orbital_speed_vis_viva", "reference", "kepler.vis_viva_speed_kms(1, 1, 1)",
              lambda: kepler.vis_viva_speed_kms(1.0, 1.0, 1.0), 29.78, 0.001,
              "Earth's mean orbital speed, 29.78 km/s (NASA planetary fact sheet)", unit="km/s"),
        Check("earth_orbital_speed_circular", "reference", "orbits.circular_orbital_speed_kms(1, 1)",
              lambda: orbits.circular_orbital_speed_kms(1.0, 1.0), 29.78, 0.001,
              "Earth's mean orbital speed, 29.78 km/s (NASA planetary fact sheet)", unit="km/s"),
        Check("earth_hill_sphere", "reference", "orbits.calculate_hill_sphere(1 AU, M_earth, M_sun)",
              lambda: orbits.calculate_hill_sphere(pc.AU_M, pc.EARTH_MASS_TO_KG, sun_kg) / 1e9, 1.5, 0.01,
              "Earth's Hill sphere, about 1.5 million km (Murray & Dermott, Solar System Dynamics)",
              unit="million km"),
        Check("habitable_zone_inner_1_lsun", "reference", "orbits.calculate_habitable_zone(L_sun)[0]",
              lambda: orbits.calculate_habitable_zone(pc.SOLAR_LUMINOSITY)[0], 0.95, 0.01,
              "inner habitable zone edge for the Sun, about 0.95 AU (Kasting et al. 1993)", unit="AU"),
        Check("habitable_zone_outer_1_lsun", "reference", "orbits.calculate_habitable_zone(L_sun)[1]",
              lambda: orbits.calculate_habitable_zone(pc.SOLAR_LUMINOSITY)[1], 1.37, 0.01,
              "outer habitable zone edge for the Sun, about 1.37 AU (Kasting et al. 1993)", unit="AU"),
        Check("snow_line_1_lsun", "reference", "formation.snow_line_au(L_sun)",
              lambda: formation.snow_line_au(pc.SOLAR_LUMINOSITY), 2.7, 1e-6,
              "snow line at 2.7 AU for 1 L_sun (Hayashi 1981)", unit="AU"),
        Check("white_dwarf_0_6_msun_radius", "reference", "stellar_evolution.white_dwarf_radius_km(0.6)",
              lambda: stellar_evolution.white_dwarf_radius_km(0.6), pc.EARTH_RADIUS_KM, 0.15,
              "a typical 0.6 M_sun white dwarf is about Earth-sized (Shapiro & Teukolsky 1983)", unit="km"),
        Check("sirius_b_radius", "reference", "stellar_evolution.white_dwarf_radius_km(1.02)",
              lambda: stellar_evolution.white_dwarf_radius_km(1.02), 5840.0, 0.05,
              "Sirius B: 1.02 M_sun, 0.0084 R_sun = 5,840 km (Bond et al. 2017)", unit="km"),
        Check("sun_schwarzschild_radius", "reference", "2 G M_sun / c^2 (compactRemnant.BlackHole)",
              lambda: _schwarzschild_radius_km(1.0), 2.953, 0.002,
              "Schwarzschild radius of 1 M_sun, 2.953 km (2GM/c^2 with GM_sun = 1.3271e20 m^3/s^2)",
              unit="km"),
        Check("sun_galactic_radius", "reference", "physical_constants.GALACTIC_CENTER_DISTANCE_LY",
              lambda: units.ly_to_pc(gc_ly) / 1000, 8.2, 0.05,
              "the Sun's distance from the galactic center, 8.2 kpc (GRAVITY Collaboration 2019)",
              unit="kpc"),
        Check("sun_galactic_orbital_speed", "reference", "galactic_orbit.calculate_galactic_orbit(Sun)[0]",
              lambda: galactic_orbit.calculate_galactic_orbit(gc_ly)[0], 225.0, 0.1,
              "the Sun's circular speed around the galaxy, about 220-230 km/s (IAU 1985; Reid et al. 2014)",
              unit="km/s"),
        Check("sun_galactic_year", "reference", "galactic_orbit.calculate_galactic_orbit(Sun)[1]",
              lambda: galactic_orbit.calculate_galactic_orbit(gc_ly)[1] * 1000, 230.0, 0.1,
              "the galactic year, about 230 million years", unit="My"),
        Check("holman_wiegert_s_type_equal_circular", "reference",
              "orbits.holman_wiegert_critical_semimajor_axis(1, 0.5, 0)",
              lambda: orbits.holman_wiegert_critical_semimajor_axis(1.0, 0.5, 0.0), 0.274, 1e-9,
              "Holman & Wiegert 1999, AJ 117:621, eq. 1 at mu = 0.5, e = 0", unit="a_bin"),
        Check("holman_wiegert_s_type_eccentric", "reference",
              "orbits.holman_wiegert_critical_semimajor_axis(1, 0.3, 0.5)",
              lambda: orbits.holman_wiegert_critical_semimajor_axis(1.0, 0.3, 0.5), 0.14505, 1e-9,
              "Holman & Wiegert 1999, AJ 117:621, eq. 1 at mu = 0.3, e = 0.5", unit="a_bin"),
        Check("holman_wiegert_p_type_equal_circular", "reference",
              "orbits.holman_wiegert_circumbinary_a_crit_au(1, 0.5, 0)",
              lambda: orbits.holman_wiegert_circumbinary_a_crit_au(1.0, 0.5, 0.0), 2.3875, 1e-9,
              "Holman & Wiegert 1999, AJ 117:621, eq. 3 at mu = 0.5, e = 0", unit="a_bin"),
        Check("holman_wiegert_p_type_eccentric", "reference",
              "orbits.holman_wiegert_circumbinary_a_crit_au(1, 0.3, 0.3)",
              lambda: orbits.holman_wiegert_circumbinary_a_crit_au(1.0, 0.3, 0.3), 3.361141, 1e-9,
              "Holman & Wiegert 1999, AJ 117:621, eq. 3 at mu = 0.3, e = 0.3", unit="a_bin"),
        Check("kepler_equation_meeus_30a", "reference", "kepler.solve_eccentric_anomaly(5 deg, 0.1)",
              lambda: math.degrees(kepler.solve_eccentric_anomaly(math.radians(5.0), 0.1)),
              5.554589, 1e-6, "Meeus, Astronomical Algorithms, example 30.a: E = 5.554589 deg",
              mode="abs", unit="deg"),
        Check("kepler_equation_high_eccentricity", "reference",
              "kepler.solve_eccentric_anomaly(M(E = 1 rad), 0.99)",
              lambda: kepler.solve_eccentric_anomaly(1.0 - 0.99 * math.sin(1.0), 0.99), 1.0, 1e-9,
              "exact: E = 1 rad gives M = 1 - 0.99 sin 1 by Kepler's equation", mode="abs", unit="rad"),
        Check("barker_equation_d_1", "reference", "kepler.solve_barker_equation(4/3)",
              lambda: kepler.solve_barker_equation(4 / 3), 1.0, 1e-12,
              "exact: D = 1 solves D^3 + 3D = 3 * 4/3", mode="abs"),
        Check("barker_equation_d_2", "reference", "kepler.solve_barker_equation(14/3)",
              lambda: kepler.solve_barker_equation(14 / 3), 2.0, 1e-12,
              "exact: D = 2 solves D^3 + 3D = 3 * 14/3", mode="abs"),
        Check("barker_equation_d_minus_3", "reference", "kepler.solve_barker_equation(-12)",
              lambda: kepler.solve_barker_equation(-12.0), -3.0, 1e-12,
              "exact: D = -3 solves D^3 + 3D = 3 * -12", mode="abs"),
    ]


# ---------------------------------------------------------------------------
# TEST.65: identities and invariants
# ---------------------------------------------------------------------------

def _unit_round_trips():
    pairs = []
    for x in _log_sweep(1e-9, 1e12, 200):
        pairs.append((units.ly_to_pc(units.pc_to_ly(x)), x))
        pairs.append((units.pc_to_ly(units.ly_to_pc(x)), x))
        pairs.append((units.au_to_ly(units.ly_to_au(x)), x))
        pairs.append((units.ly_to_au(units.au_to_ly(x)), x))
        pairs.append((units.mpc_to_pc(units.pc_to_mpc(x)), x))
        pairs.append((units.pc_to_mpc(units.mpc_to_pc(x)), x))
        pairs.append((units.milliparsecs_to_ly(units.ly_to_milliparsecs(x)), x))
        pairs.append((units.ly_to_milliparsecs(units.milliparsecs_to_ly(x)), x))
        pairs.append((x * pc.AU_TO_KM / pc.AU_TO_KM, x))
    return _sweep_max(pairs)


def _unit_chain_consistency():
    """pc -> ly, pc -> mpc -> ly and pc -> AU -> ly must all agree."""
    pairs = []
    for x in _log_sweep(1e-6, 1e6, 100):
        direct = units.pc_to_ly(x)
        pairs.append((units.milliparsecs_to_ly(units.pc_to_mpc(x)), direct))
        pairs.append((units.au_to_ly(x * pc.AU_PER_PARSEC), direct))
    return _sweep_max(pairs)


def _speed_of_light_consistency():
    """LIGHTYEAR_M, SPEED_OF_LIGHT_M_S and SPEED_OF_LIGHT_KMS describe one c."""
    c_from_ly = pc.LIGHTYEAR_M / (365.25 * 86400)
    return max(_rel_err(pc.SPEED_OF_LIGHT_M_S, c_from_ly),
               _rel_err(pc.SPEED_OF_LIGHT_KMS * 1000, c_from_ly))


def _kepler_units_vs_si():
    """keplerMotion's mu = 4 pi^2 AU^3/yr^2 per M_sun against G * M_sun in SI."""
    gm_au3_yr2 = pc.G * pc.SOLAR_MASS_TO_KG * pc.SECONDS_PER_YEAR ** 2 / pc.AU_M ** 3
    return _rel_err(gm_au3_yr2, kepler.gravitational_parameter_au3_yr2(1.0))


def _mass_monotonic_violations():
    """How many steps along a mass sweep break the expected direction:
    luminosity and radius up, lifetime down."""
    masses = _log_sweep(0.08, 150.0, 400)
    bad = 0
    for low, high in zip(masses, masses[1:]):
        if not stellar_evolution.main_sequence_luminosity_sol(high) > stellar_evolution.main_sequence_luminosity_sol(low):
            bad += 1
        if not stellar_evolution.main_sequence_lifetime_gy(high) < stellar_evolution.main_sequence_lifetime_gy(low):
            bad += 1
        if not stellar_evolution.main_sequence_radius_sol(high) > stellar_evolution.main_sequence_radius_sol(low):
            bad += 1
    return bad


def _kepler_orbit_conservation():
    """
    Steps a body around Kepler orbits with `true_anomaly_and_distance_elliptical`
    and differentiates its position numerically: the specific orbital
    energy must stay -mu/2a and the specific angular momentum
    sqrt(mu a (1 - e^2)). Returns the largest relative error.
    """
    mu = kepler.gravitational_parameter_au3_yr2(1.0)
    worst = 0.0
    for a, e in ((1.0, 0.0), (1.0, 0.3), (5.2, 0.7), (30.0, 0.95)):
        n = kepler.mean_motion_per_year(a, 1.0)
        energy = -mu / (2 * a)
        momentum = math.sqrt(mu * a * (1 - e * e))
        dt = 1e-6 / n

        def position(t):
            nu, r = kepler.true_anomaly_and_distance_elliptical(n * t, e, a)
            return r * math.cos(nu), r * math.sin(nu)

        for step in range(24):
            t = (step + 0.37) / 24 * 2 * math.pi / n
            (x0, y0), (x1, y1) = position(t - dt), position(t + dt)
            x, y = position(t)
            vx, vy = (x1 - x0) / (2 * dt), (y1 - y0) / (2 * dt)
            r = math.hypot(x, y)
            worst = max(worst,
                        _rel_err((vx * vx + vy * vy) / 2 - mu / r, energy),
                        _rel_err(x * vy - y * vx, momentum))
    return worst


def _kepler_equation_residual():
    """Max |E - e sin E - M| over a sweep of M and e."""
    worst = 0.0
    for e in (0.0, 0.1, 0.5, 0.8, 0.9, 0.99, 0.999):
        for i in range(73):
            m = i * 2 * math.pi / 72
            ecc = kepler.solve_eccentric_anomaly(m, e)
            residual = (ecc - e * math.sin(ecc) - m) % (2 * math.pi)
            worst = max(worst, min(residual, 2 * math.pi - residual))
    return worst


def _ring_volume_error():
    """Sum of a ring's cell volumes against the annulus they tile."""
    edge = tuning.DEFAULT_SECTOR_EDGE_PC
    worst = 0.0
    for ring in list(range(0, 200)) + list(range(200, 4000, 97)):
        cell = geometry.SectorCell.for_ring(ring, edge)
        total = cell.volume * geometry.ring_sector_count(ring)
        annulus = math.pi * ((ring + 1) ** 2 - ring ** 2) * edge ** 2 * edge
        worst = max(worst, _rel_err(total, annulus))
    return worst


def _sector_address_round_trip_failures():
    """Cells whose own center is not found again by `sector_address_at`,
    or whose slots don't cover `ring_sector_count` exactly."""
    edge = tuning.DEFAULT_SECTOR_EDGE_PC
    failures = 0
    for ring in list(range(0, 60)) + list(range(60, 4000, 211)):
        n = geometry.ring_sector_count(ring)
        slots = range(n) if n <= 400 else list(range(0, n, max(1, n // 97))) + [n - 1]
        for layer in (-7, -1, 0, 1, 7):
            for slot in slots:
                center = geometry.sector_position_pc(ring, layer, slot, edge)
                if geometry.sector_address_at(center, edge) != (ring, layer, slot):
                    failures += 1
        try:
            geometry.sector_position_pc(ring, 0, n, edge)
            failures += 1  # slot n must be out of range
        except ValueError:
            pass
    return failures


def _ring_count_rule_violations():
    """Rings whose slot count is not a multiple of the master wedge count,
    or whose centerline arc strays outside 0.94-1.065 edges."""
    bad = 0
    for ring in range(0, 4000):
        n = geometry.ring_sector_count(ring)
        if n % geometry.ring_master_count(ring):
            bad += 1
        arc = 2 * math.pi * (ring + 0.5) / n
        if not 0.94 <= arc <= 1.065:
            bad += 1
    return bad


def _default_galaxy_shape():
    """`generate.py plan`'s default galaxy shape."""
    return galaxyDensity.build_galaxy_shape(
        disk_scale_length_pc=2800.0, disk_scale_height_pc=350.0, bulge_scale_radius_pc=200.0,
        bulge_amplitude=1.0, arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)


def _density_at_calibration():
    shape = _default_galaxy_shape()
    radius = 2.82 * shape.disk_scale_length_pc
    theta_arm = (shape.spiral_reference_angle_rad
                 + math.log(radius / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad))
    theta = theta_arm + math.pi / shape.arm_count
    return galaxyDensity.relative_density((radius * math.cos(theta), radius * math.sin(theta), 0.0), shape)


def _population_split_error():
    """`population_densities` must sum to `relative_density` everywhere."""
    shape = _default_galaxy_shape()
    worst = 0.0
    for r in (0.0, 50.0, 1000.0, 8000.0, 15000.0):
        for z in (0.0, 30.0, 400.0, 3000.0):
            for theta in (0.0, 1.0, 2.5, 4.0):
                p = (r * math.cos(theta), r * math.sin(theta), z)
                total = galaxyDensity.relative_density(p, shape)
                split = sum(galaxyDensity.population_densities(p, shape).values())
                worst = max(worst, abs(split - total) / max(total, 1e-300))
    return worst


def _non_finite_outputs():
    """How many outputs over a fixed input sweep are NaN or infinite."""
    shape = _default_galaxy_shape()
    outputs = []
    for m in _log_sweep(0.08, 150.0, 60):
        lum = stellar_evolution.main_sequence_luminosity_sol(m)
        rad = stellar_evolution.main_sequence_radius_sol(m)
        outputs += [lum, rad, stellar_evolution.main_sequence_lifetime_gy(m),
                    stellar_evolution.effective_temperature_k(lum, rad),
                    *orbits.calculate_habitable_zone(lum * pc.SOLAR_LUMINOSITY),
                    formation.snow_line_au(lum * pc.SOLAR_LUMINOSITY),
                    stellar_evolution.white_dwarf_radius_km(min(m, 1.4)),
                    _schwarzschild_radius_km(m)]
    for a in _log_sweep(0.01, 1e5, 40):
        outputs += [planets.calculate_orbital_period_years(a, pc.SOLAR_MASS_TO_KG),
                    kepler.vis_viva_speed_kms(a, a, 1.0),
                    orbits.calculate_hill_sphere(a * pc.AU_M, pc.EARTH_MASS_TO_KG, pc.SOLAR_MASS_TO_KG)]
        for e in (0.0, 0.5, 0.99):
            outputs += list(kepler.true_anomaly_and_distance_elliptical(a, e, a))
    for ly in (1.0, 100.0, 25800.0, 60000.0, 1e6):
        outputs += list(galactic_orbit.calculate_galactic_orbit(ly))
    for mp in (-1e9, -10.0, -1e-9, 0.0, 1e-9, 10.0, 1e9):
        outputs += list(kepler.true_anomaly_and_distance_parabolic(mp, 1.0))
    for r in (0.0, 1.0, 1e3, 1e4, 3e4):
        for z in (0.0, 100.0, 1e4):
            outputs.append(galaxyDensity.relative_density((r, 0.0, z), shape))
    return sum(1 for v in outputs if not math.isfinite(v))


def _invariant_checks():
    exact = "exact identity"
    return [
        Check("unit_round_trips", "invariant",
              "units.pc_to_ly/ly_to_pc, ly_to_au/au_to_ly, pc_to_mpc/mpc_to_pc, ly_to_milliparsecs/milliparsecs_to_ly",
              _unit_round_trips, 1e-12, 0.0, f"{exact}: converting there and back changes nothing",
              mode="max", unit="relative error"),
        Check("unit_chains_agree", "invariant", "pc -> ly directly, through mpc and through AU",
              _unit_chain_consistency, 1e-12, 0.0, f"{exact}: every route between two units agrees",
              mode="max", unit="relative error"),
        Check("parsec_in_au", "invariant", "physical_constants.AU_PER_PARSEC",
              lambda: pc.AU_PER_PARSEC, 648000 / math.pi, 1e-12,
              "IAU 2015 Resolution B2: 1 pc = 648,000/pi AU", unit="AU"),
        Check("parsec_in_ly", "invariant", "units.pc_to_ly(1)",
              lambda: units.pc_to_ly(1.0), 3.261564, 1e-6,
              "1 pc = 3.261564 ly (IAU 2012 AU, Julian-year light-year)", unit="ly"),
        Check("speed_of_light_consistent", "invariant",
              "physical_constants.SPEED_OF_LIGHT_M_S, SPEED_OF_LIGHT_KMS, LIGHTYEAR_M",
              _speed_of_light_consistency, 1e-12, 0.0,
              "exact: c = 299,792,458 m/s (SI) and 1 ly = c x 365.25 days",
              mode="max", unit="relative error"),
        Check("kepler_units_match_si", "invariant", "kepler.gravitational_parameter_au3_yr2(1) vs G * M_sun",
              _kepler_units_vs_si, 1e-3, 0.0,
              "Kepler's third law in AU, yr and M_sun (mu = 4 pi^2) against G and M_sun in SI",
              mode="max", unit="relative error"),
        Check("mass_sequence_monotonic", "invariant",
              "stellar_evolution.main_sequence_luminosity_sol/radius_sol/lifetime_gy",
              _mass_monotonic_violations, 0, 0.0,
              "on the main sequence, heavier stars are brighter, larger and shorter-lived",
              mode="max", unit="violations"),
        Check("kepler_orbit_conserves_energy_and_momentum", "invariant",
              "kepler.true_anomaly_and_distance_elliptical along an orbit",
              _kepler_orbit_conservation, 1e-6, 0.0,
              "two-body problem: specific energy -mu/2a and angular momentum sqrt(mu a (1-e^2)) are constant",
              mode="max", unit="relative error"),
        Check("kepler_equation_residual", "invariant", "kepler.solve_eccentric_anomaly",
              _kepler_equation_residual, 1e-9, 0.0, f"{exact}: the solution satisfies M = E - e sin E",
              mode="max", unit="rad"),
        Check("ring_cells_fill_annulus", "invariant", "geometry.SectorCell.volume x ring_sector_count",
              _ring_volume_error, 1e-12, 0.0, f"{exact}: a ring's cells tile its annulus with no gap or overlap",
              mode="max", unit="relative error"),
        Check("sector_address_round_trip", "invariant",
              "geometry.sector_address_at(sector_position_pc(...))",
              _sector_address_round_trip_failures, 0, 0.0,
              f"{exact}: every cell's own center lies in that cell, and slot N is out of range",
              mode="max", unit="failures"),
        Check("ring_sector_count_rule", "invariant", "geometry.ring_sector_count",
              _ring_count_rule_violations, 0, 0.0,
              "design (docs/design/galaxy-coordinate-system.md): a multiple of the master wedges, "
              "arc about one edge", mode="max", unit="violations"),
        Check("density_is_1_at_calibration", "invariant", "galaxyDensity.relative_density at the calibration point",
              _density_at_calibration, 1.0, 1e-12,
              "definition (galaxyDensity.build_galaxy_shape): 1.0 at the inter-arm calibration point"),
        Check("population_densities_sum", "invariant", "galaxyDensity.population_densities",
              _population_split_error, 1e-12, 0.0, f"{exact}: the populations add up to relative_density",
              mode="max", unit="relative error"),
        Check("no_nan_or_infinity", "invariant", "stellar, orbit, galactic and density functions over a sweep",
              _non_finite_outputs, 0, 0.0, "every physical quantity here is finite for real inputs",
              mode="max", unit="non-finite values"),
    ]


# ---------------------------------------------------------------------------
# TEST.66: distributions
# ---------------------------------------------------------------------------

KROUPA_BREAKS_SOL = (0.08, 0.5, 150.0)
"""tuple: Kroupa (2001) IMF segment edges, M_sun -- written out here, not
read from `program_constants`, so a changed constant fails the check."""

KROUPA_SLOPES = (1.3, 2.3)
"""tuple: Kroupa (2001) slope alpha of each segment (dN/dM ~ M^-alpha)."""


def _kroupa_cdf(m):
    """The fraction of Kroupa IMF stars below `m` M_sun."""
    def integral(a, b, alpha):
        return (b ** (1 - alpha) - a ** (1 - alpha)) / (1 - alpha)
    (lo, mid, hi), (a1, a2) = KROUPA_BREAKS_SOL, KROUPA_SLOPES
    k2 = mid ** (a2 - a1)  # continuity at the break
    total = integral(lo, mid, a1) + k2 * integral(mid, hi, a2)
    if m <= mid:
        return integral(lo, m, a1) / total
    return (integral(lo, mid, a1) + k2 * integral(mid, m, a2)) / total


def _imf_p_value():
    rng = random.Random(SEED + 1)
    draws = [stellar_evolution.sample_imf_mass_sol(rng=rng) for _ in range(4000)]
    edges = [0.08, 0.15, 0.3, 0.5, 0.8, 1.5, 3.0, 8.0, 150.0]
    shares = [_kroupa_cdf(b) - _kroupa_cdf(a) for a, b in zip(edges, edges[1:])]
    return chi_square_p_value(_bin_counts(draws, edges), shares)


def _star_age_p_value():
    """Disk ages are uniform over 0-10 Gy (constant star formation)."""
    rng = random.Random(SEED + 2)
    draws = [stellar_evolution.sample_star_age_gy(rng=rng) for _ in range(3000)]
    edges = [float(i) for i in range(11)]
    return chi_square_p_value(_bin_counts(draws, edges), [1.0] * 10)


def _bulge_age_p_value():
    """Bulge ages are uniform over 8-12 Gy."""
    rng = random.Random(SEED + 3)
    draws = [stellar_evolution.sample_star_age_gy(rng=rng, population="bulge") for _ in range(2000)]
    edges = [8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0]
    return chi_square_p_value(_bin_counts(draws, edges), [1.0] * 8)


def _poisson_p_value():
    """Knuth's sampler at a typical sector mean against the Poisson pmf."""
    rng = random.Random(SEED + 4)
    mean = 6.5
    draws = [spaceSector._sample_poisson_count(mean, rng=rng) for _ in range(4000)]
    top = 20
    counts = [0] * (top + 1)
    for d in draws:
        counts[min(d, top)] += 1
    pmf = [math.exp(-mean) * mean ** k / math.factorial(k) for k in range(top)]
    pmf.append(1.0 - sum(pmf))
    return chi_square_p_value(counts, pmf)


def _poisson_large_mean_z():
    """Above `_POISSON_NORMAL_APPROX_MEAN` the sampler switches to a normal
    approximation: its sample mean must sit within a few standard errors."""
    rng = random.Random(SEED + 5)
    mean, n = 800.0, 2000
    draws = [spaceSector._sample_poisson_count(mean, rng=rng) for _ in range(n)]
    return abs(sum(draws) / n - mean) / math.sqrt(mean / n)


def _truncated_normal_cdf(x, mean, sd, lo, hi):
    def phi(v):
        return 0.5 * (1 + math.erf((v - mean) / (sd * math.sqrt(2))))
    return (phi(x) - phi(lo)) / (phi(hi) - phi(lo))


def _bounded_bell_p_value():
    """`sample_bounded_bell(0, 1, 0.27)` is a normal(0.27, 0.09) cut to [0, 1]."""
    with _seeded_global_random(SEED + 6):
        draws = [sampling.sample_bounded_bell(0.0, 1.0, 0.27) for _ in range(4000)]
    mean, sd = 0.27, 0.27 / 3
    edges = [i / 20 for i in range(21)]
    shares = [_truncated_normal_cdf(b, mean, sd, 0.0, 1.0) - _truncated_normal_cdf(a, mean, sd, 0.0, 1.0)
              for a, b in zip(edges, edges[1:])]
    return chi_square_p_value(_bin_counts(draws, edges), shares)


def _planet_class_p_value():
    """Unrestricted planet class draws follow `PLANET_CLASS_PROBABILITIES`."""
    table = tuning.PLANET_CLASS_PROBABILITIES
    classes = sorted(table, key=lambda c: -table[c])
    with _seeded_global_random(SEED + 7):
        draws = [planets._choose_weighted_planet_class(classes) for _ in range(5000)]
    counts = [draws.count(c) for c in classes]
    return chi_square_p_value(counts, [table[c] for c in classes])


def _distribution_checks():
    chi = f"chi-square goodness of fit, p >= {CHI_SQUARE_P_MIN:g}"
    return [
        Check("imf_matches_kroupa", "distribution", "stellar_evolution.sample_imf_mass_sol",
              _imf_p_value, CHI_SQUARE_P_MIN, 0.0,
              f"Kroupa 2001, MNRAS 322:231 (alpha 1.3 below 0.5 M_sun, 2.3 above); {chi}", mode="min",
              unit="p"),
        Check("star_ages_uniform", "distribution", "stellar_evolution.sample_star_age_gy",
              _star_age_p_value, CHI_SQUARE_P_MIN, 0.0,
              f"constant disk star formation over 0-10 Gy (tuning.STAR_FORMATION_AGE_RANGE_GY); {chi}",
              mode="min", unit="p"),
        Check("bulge_ages_uniform", "distribution", "stellar_evolution.sample_star_age_gy(population='bulge')",
              _bulge_age_p_value, CHI_SQUARE_P_MIN, 0.0,
              f"bulge stars 8-12 Gy old (tuning.STELLAR_POPULATION_AGE_RANGES_GY); {chi}",
              mode="min", unit="p"),
        Check("sector_counts_poisson", "distribution", "spaceSector._sample_poisson_count(6.5)",
              _poisson_p_value, CHI_SQUARE_P_MIN, 0.0,
              f"Poisson pmf with mean 6.5; {chi}", mode="min", unit="p"),
        Check("sector_counts_large_mean", "distribution", "spaceSector._sample_poisson_count(800)",
              _poisson_large_mean_z, 4.0, 0.0,
              "Poisson mean 800 (normal approximation): sample mean within 4 standard errors",
              mode="max", unit="standard errors"),
        Check("bounded_bell_shape", "distribution", "sampling.sample_bounded_bell(0, 1, 0.27)",
              _bounded_bell_p_value, CHI_SQUARE_P_MIN, 0.0,
              f"normal(0.27, 0.09) truncated to [0, 1]; {chi}", mode="min", unit="p"),
        Check("planet_class_shares", "distribution", "planets._choose_weighted_planet_class",
              _planet_class_p_value, CHI_SQUARE_P_MIN, 0.0,
              f"tuning.PLANET_CLASS_PROBABILITIES; {chi}", mode="min", unit="p"),
    ]


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def all_checks():
    """Every check, in the order they run: reference, invariant, distribution."""
    return _reference_checks() + _invariant_checks() + _distribution_checks()


def run_check(check):
    """Runs one check; an exception is a failure, never a crash."""
    start = time.perf_counter()
    try:
        actual = check.compute()
    except Exception as exc:  # noqa: BLE001 -- any error is a failed check
        return Result(check, False, error=f"{type(exc).__name__}: {exc}",
                      seconds=time.perf_counter() - start)
    actual = float(actual)
    return Result(check, _compare(check, actual), actual=actual, seconds=time.perf_counter() - start)


def run_all(checks=None):
    """Runs `checks` (default: `all_checks()`) and returns their `Result`s."""
    return [run_check(check) for check in (all_checks() if checks is None else checks)]


def failures(results):
    """The failed `Result`s."""
    return [r for r in results if not r.passed]


def format_report(results, verbose=False):
    """
    A plain-text report: one line per failure (or per check, `verbose`),
    then a summary line.
    """
    failed = failures(results)
    lines = [r.describe() for r in results if verbose or not r.passed]
    seconds = sum(r.seconds for r in results)
    if failed:
        lines.append(f"Math check FAILED: {len(failed)} of {len(results)} checks failed "
                     f"({', '.join(r.name for r in failed)}) in {seconds:.2f} s.")
    else:
        lines.append(f"Math check passed: all {len(results)} checks in {seconds:.2f} s.")
    return "\n".join(lines)


_STARTUP_RESULTS = None


def startup_failures():
    """
    The failed checks of this process's one startup run, running it the
    first time it's asked for (the website calls this once while starting,
    then reads the cached answer for its admin warning).
    """
    global _STARTUP_RESULTS
    if _STARTUP_RESULTS is None:
        _STARTUP_RESULTS = run_all()
        failed = failures(_STARTUP_RESULTS)
        if failed:
            log.error(format_report(_STARTUP_RESULTS))
        else:
            log.debug(format_report(_STARTUP_RESULTS))
    return failures(_STARTUP_RESULTS)


def main(argv=None):
    """`python -m planetgen.physics.mathcheck [-v]`: prints the report, exits 1
    if any check failed."""
    argv = sys.argv[1:] if argv is None else argv
    results = run_all()
    print(format_report(results, verbose="-v" in argv or "--verbose" in argv))
    return 1 if failures(results) else 0


if __name__ == "__main__":
    sys.exit(main())
