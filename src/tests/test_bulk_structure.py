"""
Radius, density and gravity from the bulk fractions (GEN.429, `planetgen.physics.bulk_structure`).

The check numbers are part B's (docs/design/body-composition-parts/B-rocky.md section 3.9) and part C's (C-volatile-and-giants.md
section 2.4); they are test targets, never constants of the generator.
"""

import math

import pytest

from planetgen.physics import bulk_composition as bc
from planetgen.physics import bulk_structure as bs


def test_the_solver_reproduces_part_bs_earth_sized_check():
    result = bs.layered_structure(1.0, 0.3, 0.7, 0.0, x_fe=0.11, core_scale=0.9)
    assert result["radius"] == pytest.approx(0.9996, abs=2e-4)
    assert result["central_pressure"] == pytest.approx(3.7049e11, rel=1e-3)
    assert result["cmb_pressure"] == pytest.approx(1.4507e11, rel=1e-3)
    assert result["gravity"] == pytest.approx(9.8278, abs=2e-3)
    assert result["escape"] == pytest.approx(11188.3, abs=2.0)
    assert result["moment_of_inertia"] == pytest.approx(0.331, abs=1e-3)


def test_the_solver_limits_are_the_coreless_and_the_pure_iron_planet():
    assert 1.07 <= bs.layered_structure(1.0, 0.0, 1.0, 0.0, 0.0)["radius"] <= 1.085
    assert 0.76 <= bs.layered_structure(1.0, 1.0, 0.0, 0.0, 0.0)["radius"] <= 0.78


def test_the_solver_water_layer_matches_part_cs_integration():
    # C 2.4: pure water at 1 Earth mass 1.401 ([R4882]: 1.41); 50 percent water on CMF 0.325 rock 1.2301
    assert bs.layered_structure(1.0, 0.0, 0.0, 1.0)["radius"] == pytest.approx(1.401, rel=0.01)
    assert bs.layered_structure(1.0, 0.1625, 0.3375, 0.5)["radius"] == pytest.approx(1.2301, rel=0.01)


def test_the_radius_grows_with_mass_and_falls_with_the_core():
    radii = [bs.radius_fit_earth(mass, 0.3) for mass in (0.1, 0.5, 1.0, 5.0, 20.0)]
    assert radii == sorted(radii)
    cores = [bs.radius_fit_earth(1.0, cmf) for cmf in (0.0, 0.2, 0.4, 0.8)]
    assert cores == sorted(cores, reverse=True)
    assert bs.radius_fit_earth(1.0, 0.3, x_fe=0.2) < bs.radius_fit_earth(1.0, 0.3, x_fe=0.0)


@pytest.mark.parametrize("mass,cmf,x_fe,scale", [(0.3, 0.1, 0.05, 0.85), (1.0, 0.3, 0.11, 0.9), (4.0, 0.5, 0.2, 1.0),
                                                 (15.0, 0.25, 0.0, 0.8)])
def test_the_fast_fit_agrees_with_the_solver(mass, cmf, x_fe, scale):
    solver = bs.layered_structure(mass, cmf, 1.0 - cmf, 0.0, x_fe, scale)["radius"]
    assert bs.radius_fit_earth(mass, cmf, x_fe, scale) == pytest.approx(solver, rel=0.005)


def test_the_inverse_of_the_fast_fit_round_trips():
    radius = bs.radius_fit_earth(3.0, 0.37, 0.15, 0.85)
    assert bs.cmf_from_mass_radius(3.0, radius, 0.15, 0.85) == pytest.approx(0.37, abs=1e-3)
    # M 5, R 2.3 is too big for a dry rock
    assert bs.cmf_from_mass_radius(5.0, 2.3) is None


def test_the_water_closed_form_reproduces_part_cs_limits():
    assert bs.water_world_radius_earth(1.0, 0.325, 0.0) == pytest.approx(0.993, abs=0.002)
    assert bs.water_world_radius_earth(1.0, 0.0, 1.0) == pytest.approx(1.378, abs=0.005)
    # against the integrated reference within the fit's stated 3.6 percent
    for mass, cmf, water in ((0.5, 0.325, 0.25), (2.0, 0.325, 0.1), (8.0, 0.1, 0.5)):
        solver = bs.layered_structure(mass, cmf * (1 - water), (1 - cmf) * (1 - water), water)["radius"]
        assert bs.water_world_radius_earth(mass, cmf, water) == pytest.approx(solver, rel=0.036)


def test_the_envelope_form_is_near_lopez_and_fortneys_table():
    # 1 Gyr, enhanced opacity, flux 10 (the paper's table, CSV of GEN.425); the fit's own accuracy is rms 6.3 percent
    table = [(10.0, 2.4, 0.01, 2.17), (10.0, 2.4, 0.05, 3.47), (10.0, 8.5, 0.01, 2.37), (10.0, 8.5, 0.05, 3.31),
             (10.0, 8.5, 0.10, 4.27)]
    for flux, mass, envelope, radius in table:
        assert bs.envelope_radius_earth(mass, envelope, flux, 1.0) == pytest.approx(radius, rel=0.12)
    assert bs.core_radius_earth(16.0) == pytest.approx(2.0)


def test_the_one_bar_level_is_below_the_20_millibar_level():
    top = bs.envelope_radius_earth(5.0, 0.05, 1.0, 4.5)
    bottom = bs.envelope_radius_earth(5.0, 0.05, 1.0, 4.5, temperature_k=278.0, gravity_m_s2=15.0)
    assert 0.0 < top - bottom < 0.2


def test_the_mixing_rule_gives_a_heap_of_grains():
    rock = bs.mixing_density_kg_m3(0.0, 1.0, 0.0, 0.0)
    ice = bs.mixing_density_kg_m3(0.0, 0.0, 1.0, 0.0)
    assert rock == pytest.approx(3189.0, rel=0.01) and ice == pytest.approx(917.0, rel=0.01)
    assert bs.mixing_density_kg_m3(0.0, 0.7, 0.3, 0.0) == pytest.approx(1.0 / (0.7 / rock + 0.3 / ice))
    assert bs.mixing_density_kg_m3(0.0, 0.7, 0.3, 0.0, porosity=0.4) == pytest.approx(
        0.6 * bs.mixing_density_kg_m3(0.0, 0.7, 0.3, 0.0))


def _bulk(**shares):
    values = dict(metal=0.0, rock=0.0, water=0.0, ice=0.0, envelope=0.0)
    values.update(shares)
    return bc.BulkFractions(**values)


def test_structure_takes_the_closed_form_in_range_and_the_solver_outside():
    earth = _bulk(metal=0.3, rock=0.7, rock_feo=0.08)
    dry = bs.structure(earth, bs.M_EARTH)
    assert dry.method == "closed-dry"
    assert dry.radius_km == pytest.approx(6371.0, rel=0.03)
    assert dry.gravity_g == pytest.approx(1.0, rel=0.07)
    assert dry.density_g_cm3 == pytest.approx(5.5, rel=0.1)
    solver = bs.structure(earth, bs.M_EARTH, prefer_solver=True)
    assert solver.method == "solver" and solver.moment_of_inertia == pytest.approx(0.33, abs=0.02)
    assert solver.radius_km == pytest.approx(dry.radius_km, rel=0.01)
    # a body of 0.05 Earth masses is below part B's 0.1: the small-body fit
    small = bs.structure(earth, 0.05 * bs.M_EARTH)
    assert small.method == "closed-small" and 2000.0 < small.radius_km < 3200.0

    ocean = bs.structure(_bulk(metal=0.15, rock=0.35, water=0.5), 2.0 * bs.M_EARTH)
    assert ocean.method == "closed-water" and ocean.density_g_cm3 < dry.density_g_cm3

    sub_neptune = bs.structure(_bulk(metal=0.15, rock=0.8, envelope=0.05), 5.0 * bs.M_EARTH, flux_earth=10.0, age_gy=1.0)
    assert sub_neptune.method == "closed-envelope" and sub_neptune.radius_km > 2.5 * 6371.0 * 0.8


def test_a_heap_of_grains_and_a_giant():
    small = bs.structure(_bulk(rock=0.7, water=0.3, porosity=0.3), 1.0e19)
    assert small.method == "mixing"
    assert small.density_g_cm3 == pytest.approx(0.7 * bs.mixing_density_kg_m3(0.0, 0.7, 0.3, 0.0) / 1000.0, rel=1e-9)
    assert small.gravity_g < 0.01
    with pytest.raises(bs.OutsideRange, match="giant"):
        bs.structure(_bulk(rock=0.1, envelope=0.9), 100.0 * bs.M_EARTH)
    with pytest.raises(bs.OutsideRange):
        bs.layered_structure(1e6, 0.3, 0.7, 0.0)


def test_the_feo_conversion_matches_part_bs_example():
    # B 3.9 item 2: 8 wt percent FeO in the mantle is x_Fe 0.1073 (that example has Si/(Fe+Mg) of the stoichiometry,
    # here the single-silicate-per-cation closure)
    assert 0.08 < bs.x_fe_from_feo_mass_fraction(0.08) < 0.14
    assert bs.x_fe_from_feo_mass_fraction(0.0) == pytest.approx(0.0, abs=1e-9)


def test_the_water_isotherm_is_continuous_and_rising():
    pressures = [10.0 ** (5.0 + 0.1 * i) for i in range(0, 100)]
    densities = [bs.water_density(p) for p in pressures]
    assert densities == sorted(densities)
    assert densities[0] == pytest.approx(996.56)
    assert bs.water_density(1.0e9) == pytest.approx(1351.1, rel=0.01)
    assert not any(math.isnan(d) for d in densities)


@pytest.mark.parametrize("mass,cmf,water", [(0.002, 0.3, 0.0), (0.025, 0.2, 0.5), (0.1, 0.4, 0.0), (0.2, 0.1, 0.8), (0.0015, 0.0, 1.0)])
def test_the_small_body_fit_agrees_with_the_solver(mass, cmf, water):
    dry = 1.0 - water
    solver = bs.layered_structure(mass, cmf * dry, (1.0 - cmf) * dry, water, 0.11, 0.9)["radius"]
    assert bs.small_body_radius_earth(mass, cmf, water) == pytest.approx(solver, rel=0.035)


def test_the_solver_branch_is_still_there_when_asked_for():
    earth = _bulk(metal=0.3, rock=0.7)
    assert bs.structure(earth, 0.05 * bs.M_EARTH, prefer_solver=True).method == "solver"
