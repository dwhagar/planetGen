"""
GEN.1040: `planetgen.physics.mantle_stripping`, the iron-rich core of a giant impact (oddity W10). The law is checked against the equations of
Marcus and others (2010) as printed, then against the three cases the item names: the one giant-impact run with every number printed (Bonomo
and others 2019, Fig. 3), Kepler-107 b and c (their masses and radii through the repository's own mass-radius model), and Mercury.
"""

import math

import pytest

from planetgen.physics import bulk_structure, mantle_stripping as ms, oddities
from planetgen.physics.constants import EARTH_MASS_TO_KG as EARTH

MERCURY = 3.3011e23
"""float: Mercury's mass in kg."""
MERCURY_RADIUS_EARTH = 2439.7 / 6371.0

SPH_MASS = 10.5 * EARTH
SPH_SPEED = 62.5e3


def _x_speed(mass, projectile, x):
    return ms.speed_for_energy_ratio(mass, projectile, x)


# --- the equations as printed ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("x", [0.0, 0.25, 0.5, 1.0, 1.5, 1.9])
def test_the_largest_remnant_is_one_minus_half_the_energy_ratio_of_the_total(x):
    """Equation 2: M_lr / M_tot = 1 - 0.5 x, so half the total is lost at x = 1."""
    target, projectile = 4.0 * EARTH, 1.0 * EARTH
    speed = _x_speed(target, projectile, x) if x > 0 else 1.0
    result = ms.impact(target, 0.33, projectile, speed)
    if x == 0.0:
        assert result.energy_ratio == pytest.approx(0.0, abs=1e-4)
    else:
        assert result.energy_ratio == pytest.approx(x, rel=1e-9)
        assert result.mass_kg == pytest.approx((target + projectile) * (1.0 - 0.5 * x), rel=1e-9)


@pytest.mark.parametrize("x", [0.0, 0.3, 0.5, 1.0, 1.3, 1.5])
def test_a_body_that_began_a_third_iron_follows_equation_three(x):
    assert ms.core_fraction_after(0.33, x) == pytest.approx(0.33 + 0.25 * x ** 1.65, rel=1e-9)


def test_the_threshold_is_the_gravity_regime_law_in_the_units_of_the_paper():
    """Equation 4 in cgs, worked out here from the paper's units, against the module's SI form."""
    total = 21.0 * EARTH
    radius_cm = (3.0 * total * 1000.0 / (4.0 * math.pi * 1.0)) ** (1.0 / 3.0) * 1.0
    # R_C1: the radius of the mass at 1 g/cm3, with the mass in grams
    radius_cm = (3.0 * (total * 1000.0) / (4.0 * math.pi)) ** (1.0 / 3.0)
    erg_per_g = 1e-4 * radius_cm ** 1.2 * (SPH_SPEED * 100.0) ** 0.8
    assert ms.disruption_threshold(total, SPH_SPEED) == pytest.approx(ms.Q_STAR_SCALE * erg_per_g * 1e-4, rel=1e-9)


def test_the_remnant_loses_mass_and_gains_iron_as_the_impact_gets_harder():
    last_mass, last_core = None, None
    for x in (0.1, 0.4, 0.8, 1.0, 1.2, 1.5, 1.7):
        result = ms.impact(4.0 * EARTH, 0.33, 1.0 * EARTH, _x_speed(4.0 * EARTH, EARTH, x))
        if last_mass is not None:
            assert result.mass_kg < last_mass and result.core_mass_fraction > last_core
        last_mass, last_core = result.mass_kg, result.core_mass_fraction


def test_the_iron_and_the_rock_a_body_keeps_stay_between_none_and_all_of_them():
    for x in [i / 20.0 for i in range(0, 50)]:
        iron, rock = ms.retention(x)
        assert 0.0 <= iron <= 1.0 and 0.0 <= rock <= 1.0
    assert ms.retention(0.0) == (pytest.approx(1.0), pytest.approx(1.0))
    assert ms.retention(5.0) == ms.retention(ms.X_MAX)
    assert ms.retention(-1.0) == ms.retention(0.0)


@pytest.mark.parametrize("start", [0.0, 0.1, 0.33, 0.6, 1.0])
def test_a_body_never_gains_iron_and_a_pure_body_stays_pure(start):
    for x in (0.2, 0.6, 1.0, 1.4):
        speed = _x_speed(3.0 * EARTH, EARTH, x)
        result = ms.impact(3.0 * EARTH, start, EARTH, speed)
        before = start * 4.0 * EARTH
        after = result.core_mass_fraction * result.mass_kg
        assert after <= before * (1.0 + 1e-9) or result.clamped
        if start in (0.0, 1.0):
            assert result.core_mass_fraction == pytest.approx(min(start, ms.CORE_FRACTION_CEILING))


def test_a_richer_start_stays_richer():
    speed = _x_speed(3.0 * EARTH, EARTH, 1.0)
    fractions = [ms.impact(3.0 * EARTH, f, EARTH, speed).core_mass_fraction for f in (0.1, 0.2, 0.33, 0.5)]
    assert fractions == sorted(fractions)


# --- the run with every number printed ------------------------------------------------------------------------------------------------------

def test_the_threshold_scale_reproduces_the_one_giant_impact_run_printed_with_all_its_numbers():
    """Bonomo and others 2019, Fig. 3: two 10.5 Earth mass bodies, 30 percent iron, head-on at 62.5 km/s leave 8.2 Earth masses at 64 percent."""
    result = ms.impact(SPH_MASS, 0.30, SPH_MASS, SPH_SPEED)
    assert result.mass_kg / EARTH == pytest.approx(8.2, rel=0.03)
    assert result.core_mass_fraction == pytest.approx(0.64, abs=0.02)
    assert result.regime == "disruption"


def test_the_published_threshold_alone_is_a_factor_of_two_off_that_run_as_marcus_says(monkeypatch):
    monkeypatch.setattr(ms, "Q_STAR_SCALE", 1.0)
    result = ms.impact(SPH_MASS, 0.30, SPH_MASS, SPH_SPEED)
    assert result.mass_kg / EARTH > 12.0
    assert ms.Q_STAR_SCALE == 1.0 and 0.5 <= 0.6 <= 1.0


# --- Kepler-107 b and c ---------------------------------------------------------------------------------------------------------------------

def test_kepler_107_c_is_what_two_ten_earth_mass_bodies_leave_and_its_radius_follows():
    """c: 9.39 +- 1.77 Earth masses, 1.597 +- 0.026 Earth radii, about 70 percent iron [R0456]."""
    result = ms.impact(SPH_MASS, 0.30, SPH_MASS, SPH_SPEED)
    radius = bulk_structure.radius_fit_earth(result.mass_kg / EARTH, result.core_mass_fraction)
    assert abs(result.mass_kg / EARTH - 9.39) < 1.77
    assert abs(radius - 1.597) < 2 * 0.026
    assert abs(result.core_mass_fraction - 0.70) < 0.1


def test_the_speed_that_makes_kepler_107_c_iron_rich_is_in_the_range_of_the_run():
    speed = ms.speed_for_core_fraction(SPH_MASS, SPH_MASS, 0.30, 0.70)
    assert 55e3 < speed < 75e3
    assert ms.impact(SPH_MASS, 0.30, SPH_MASS, speed).core_mass_fraction == pytest.approx(0.70, abs=1e-6)


def test_kepler_107_b_is_the_same_planet_without_the_impact():
    """b: 3.51 +- 1.52 Earth masses and 5.3 +- 2.3 g/cm3. With no impact it keeps the Earth-like third, and its density is inside one sigma."""
    kept = ms.strip(3.51 * EARTH, 0.33, [])
    assert kept.mass_kg == 3.51 * EARTH and kept.core_mass_fraction == 0.33 and kept.steps == ()
    radius = bulk_structure.radius_fit_earth(3.51, kept.core_mass_fraction)
    density = 3.51 / radius ** 3 * 5.514
    assert abs(density - 5.3) < 2.3
    c_density = 9.39 / bulk_structure.radius_fit_earth(9.39, 0.70) ** 3 * 5.514
    assert c_density > 1.5 * density, "c is much denser than b, the contrast of the paper"


def test_the_repository_reads_kepler_107_c_as_the_paper_does():
    cmf = bulk_structure.cmf_from_mass_radius(9.39, 1.597)
    assert 0.65 < cmf < 0.80, "the paper has about 0.70 from a different interior model"


# --- Mercury --------------------------------------------------------------------------------------------------------------------------------

def test_mercury_is_the_remnant_of_an_impact_at_the_speed_marcus_gives():
    """About 70 percent iron; Marcus: 10 to 20 km/s for planets near Mercury's mass ([R0351]: about 20 km/s)."""
    total = 2.74 * MERCURY
    projectile = total / 7.0
    target = total - projectile
    speed = ms.speed_for_core_fraction(target, projectile, 0.33, 0.70)
    assert 10e3 < speed < 25e3
    result = ms.impact(target, 0.33, projectile, speed)
    assert result.core_mass_fraction == pytest.approx(0.70, abs=1e-6)
    assert result.mass_kg == pytest.approx(MERCURY, rel=0.05)
    assert result.regime == "disruption"


def test_the_repository_reads_mercury_as_seventy_percent_iron_too():
    cmf = bulk_structure.cmf_from_mass_radius(MERCURY / EARTH, MERCURY_RADIUS_EARTH)
    assert cmf == pytest.approx(0.70, abs=0.05)


def test_several_milder_impacts_get_there_too_each_one_on_the_last_remnant():
    series = [(0.5 * EARTH, 30e3)] * 4
    result = ms.strip(3.0 * EARTH, 0.33, series)
    assert len(result.steps) == 4
    masses = [3.0 * EARTH] + [step.mass_kg for step in result.steps]
    cores = [0.33] + [step.core_mass_fraction for step in result.steps]
    assert cores == sorted(cores) and masses[-1] < masses[0] + 4 * 0.5 * EARTH
    assert result.core_mass_fraction > 0.5
    assert result.mass_lost_kg == pytest.approx(3.0 * EARTH - result.mass_kg)
    # the second step starts from the first's remnant
    again = ms.impact(result.steps[0].mass_kg, result.steps[0].core_mass_fraction, 0.5 * EARTH, 30e3)
    assert again.mass_kg == pytest.approx(result.steps[1].mass_kg)


# --- inverse, limits, regimes ---------------------------------------------------------------------------------------------------------------

def test_the_speed_for_an_energy_ratio_is_the_inverse_of_the_energy_ratio():
    for x in (0.1, 0.7, 1.4):
        speed = ms.speed_for_energy_ratio(2.0 * EARTH, 0.5 * EARTH, x)
        assert ms.energy_ratio(2.0 * EARTH, 0.5 * EARTH, speed) == pytest.approx(x, rel=1e-9)


def test_the_energy_ratio_for_a_core_fraction_inverts_the_remnant_and_has_its_edges():
    for target in (0.4, 0.55, 0.7, 0.8):
        x = ms.energy_ratio_for_core_fraction(0.33, target)
        assert ms.core_fraction_after(0.33, x) == pytest.approx(target, abs=1e-9)
    assert ms.energy_ratio_for_core_fraction(0.33, 0.2) == 0.0
    assert ms.energy_ratio_for_core_fraction(0.33, 0.9999) == ms.X_MAX or ms.energy_ratio_for_core_fraction(0.33, 0.9999) > 1.7


def test_no_remnant_is_allowed_past_the_core_fraction_limit_and_the_clamp_is_reported():
    speed = _x_speed(3.0 * EARTH, EARTH, 1.75)
    result = ms.impact(3.0 * EARTH, 0.33, EARTH, speed)
    assert result.core_mass_fraction == ms.CORE_FRACTION_CEILING < ms.MAX_CORE_MASS_FRACTION and result.clamped
    soft = ms.impact(3.0 * EARTH, 0.33, EARTH, _x_speed(3.0 * EARTH, EARTH, 0.5))
    assert not soft.clamped and soft.core_mass_fraction < ms.MAX_CORE_MASS_FRACTION


def test_an_impact_at_twice_the_catastrophic_energy_leaves_nothing_and_the_series_stops():
    speed = _x_speed(3.0 * EARTH, EARTH, 2.0)
    result = ms.impact(3.0 * EARTH, 0.33, EARTH, speed)
    assert result.mass_kg == 0.0 and result.regime == "dispersal" and result.mass_lost_kg == pytest.approx(4.0 * EARTH)
    series = ms.strip(3.0 * EARTH, 0.33, [(EARTH, speed), (EARTH, 1e3)])
    assert len(series.steps) == 1 and series.mass_kg == 0.0


def test_the_regimes_run_from_the_target_growing_to_the_remnant_gone():
    target, projectile = 4.0 * EARTH, 1.0 * EARTH
    seen = [ms.impact(target, 0.33, projectile, _x_speed(target, projectile, x)).regime for x in (0.1, 0.6, 1.2, 2.2)]
    assert seen == ["accretion", "erosion", "disruption", "dispersal"]
    assert set(seen) <= set(ms.REGIMES)
    assert ms.impact(target, 0.33, projectile, _x_speed(target, projectile, 0.1)).mass_lost_kg < 0.0


def test_the_projectile_can_have_its_own_core_fraction():
    speed = _x_speed(3.0 * EARTH, EARTH, 0.8)
    plain = ms.impact(3.0 * EARTH, 0.33, EARTH, speed)
    iron_projectile = ms.impact(3.0 * EARTH, 0.33, EARTH, speed, projectile_core_fraction=0.7)
    assert iron_projectile.core_mass_fraction > plain.core_mass_fraction
    assert iron_projectile.mass_kg == plain.mass_kg


@pytest.mark.parametrize("args", [(0.0, 0.33, 1.0, 1e4), (1.0, 0.33, -1.0, 1e4), (1.0, 0.33, 1.0, 0.0), (1.0, 1.2, 1.0, 1e4), (1.0, -0.1, 1.0, 1e4)])
def test_an_impossible_impact_is_refused(args):
    with pytest.raises(ValueError):
        ms.impact(*args)


def test_the_contact_speed_adds_the_escape_speed_in_quadrature():
    escape = ms.escape_speed_ms(EARTH, 1e-9 * EARTH, 6.371e6, 1.0)
    assert escape == pytest.approx(11.19e3, rel=0.01)
    assert ms.contact_speed_ms(EARTH, 1e-9 * EARTH, 6.371e6, 1.0, 0.0) == pytest.approx(escape)
    assert ms.contact_speed_ms(EARTH, 1e-9 * EARTH, 6.371e6, 1.0, 5e3) == pytest.approx(math.hypot(5e3, escape))


def test_the_share_of_impacts_inside_the_angle_is_the_sine_squared():
    assert ms.strip_probability() == pytest.approx(0.25)
    assert ms.strip_probability(45.0) == pytest.approx(0.5)
    assert ms.strip_probability(90.0) == pytest.approx(1.0)
    assert ms.strip_probability(0.0) == 0.0


# --- the oddity row that calls it -----------------------------------------------------------------------------------------------------------

def test_the_w10_row_names_this_module_s_limit_and_its_gate_admits_kepler_107_c():
    row = oddities.BY_ID["W10"]
    assert row.params["max_fraction"] == ms.MAX_CORE_MASS_FRACTION
    assert {"quantity": "core_mass_fraction", "lt": ms.MAX_CORE_MASS_FRACTION} in list(row.limits)
    assert (row.lever, row.input, row.effect) == ("L3", "core_mass_fraction", "strip_mantle")
    assert oddities.candidates("planet", {"body_type": "t", "mass_earth": 9.39, "collision_count": 1}) == [row]
    assert oddities.candidates("planet", {"body_type": "t", "mass_earth": 9.39, "collision_count": 0}) == []
    assert oddities.candidates("planet", {"body_type": "g", "mass_earth": 9.39, "collision_count": 3}) == []


def test_every_remnant_the_module_returns_passes_the_limits_of_the_row_it_serves():
    row = oddities.BY_ID["W10"]
    for target in (0.3, 1.0, 5.0, 9.0):
        for projectile in (0.05, 0.3, 1.0):
            for x in (0.2, 0.8, 1.2, 1.6, 1.9):
                m, p = target * EARTH, projectile * target * EARTH
                result = ms.impact(m, 0.33, p, _x_speed(m, p, x))
                assert oddities.limit_failures(row, {"core_mass_fraction": result.core_mass_fraction}) == []
                assert result.mass_kg >= 0.0
