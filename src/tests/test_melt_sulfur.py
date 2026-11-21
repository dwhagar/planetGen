"""
The sulfur a lava world's melt holds (GEN.1814, `planetgen.physics.melt_sulfur`): the SCSS law of [R5340], the molten column and its place in the sulfur balance.
"""

import math

import pytest

from planetgen.physics import constants, magma_ocean_redox, melt_sulfur as ms, mixture_lava, sulfur_balance as sb, tidal_heating

LIFETIMES = sb.Lifetimes(wet_s=sb.DAY_S, dry_s=2.0e6 * sb.YEAR_S)
AGE_S = 4.5e9 * sb.YEAR_S


def test_the_law_reproduces_the_measured_runs_of_its_own_paper_within_a_factor_of_one_and_a_half():
    # [R5340: Table 1] (diamond anvil runs, thermal pressures): 4300 K and 43 GPa gave 10,837 ppm, 3927 K and 44 GPa 11,806 ppm, 4605 K and 53 GPa 6,979 ppm.
    for temperature, pressure, measured in ((4300.0, 43.0, 10837.0), (3927.0, 44.0, 11806.0), (4605.0, 53.0, 6979.0)):
        assert ms.scss_ppm(temperature, pressure) == pytest.approx(measured, rel=0.5)         # the fit is to 42 values with R2 0.86: a run is not on the line


def test_the_content_rises_with_temperature_and_falls_with_pressure():
    assert ms.scss_ppm(3000.0, 10.0) > ms.scss_ppm(2000.0, 10.0)
    assert ms.scss_ppm(3000.0, 10.0) < ms.scss_ppm(3000.0, 0.0)


def test_a_body_hotter_than_the_fit_is_held_at_its_hottest_run():
    assert ms.scss_ppm(9000.0, 20.0) == ms.scss_ppm(ms.FIT_MAX_K, 20.0)


def test_the_mean_over_the_column_lies_between_its_top_and_its_base():
    top, base = ms.scss_ppm(3000.0, 0.0), ms.scss_ppm(3000.0, 30.0)
    assert base < ms.scss_mean_ppm(3000.0, 30.0) < top
    assert ms.scss_mean_ppm(3000.0, 0.0) == top
    # the closed form is the mean of the law over equal steps of pressure
    steps = [ms.scss_ppm(3000.0, 30.0 * (i + 0.5) / 2000) for i in range(2000)]
    assert ms.scss_mean_ppm(3000.0, 30.0) == pytest.approx(sum(steps) / len(steps), rel=1e-5)


def test_a_surface_under_the_melting_point_has_no_melt():
    assert ms.melt_base_gpa(mixture_lava.MELT_MIN_K, 6.0e24) == 0.0
    assert ms.melt_base_gpa(None, 6.0e24) == 0.0
    assert ms.capacity_mol_m2(1000.0, 9.8, 6.4e6) == 0.0
    assert ms.capacity_mol_m2(2000.0, 0.0, 6.4e6) == 0.0


def test_the_melt_reaches_where_the_solidus_meets_the_surface_temperature_and_no_deeper_than_the_ocean():
    mass = 1.0 * constants.EARTH_MASS_TO_KG
    shallow = ms.melt_base_gpa(tidal_heating.solidus_k(3.0), mass)
    assert shallow == pytest.approx(3.0)            # the solidus is 1,750 K at 3 GPa: a surface at that temperature melts to 3 GPa
    ocean = magma_ocean_redox.base_pressure_gpa(1.0)
    assert ms.melt_base_gpa(6000.0, mass) == pytest.approx(ocean)


def test_a_deeper_melt_holds_more_sulfur():
    earth_mass = constants.EARTH_MASS_TO_KG
    radius = 6.371e6
    gravity = constants.G * earth_mass / radius ** 2
    assert 0.0 < ms.capacity_mol_m2(1800.0, gravity, radius) < ms.capacity_mol_m2(2800.0, gravity, radius)


def test_the_capacity_is_the_column_times_its_mean_content_over_the_molar_mass():
    mass, radius = 1.824e24, 4.697e6
    gravity = constants.G * mass / radius ** 2
    base = ms.melt_base_gpa(2860.0, mass)
    expected = base * 1.0e9 / gravity * ms.scss_mean_ppm(2860.0, base) * 1.0e-6 / 0.03206
    assert ms.capacity_mol_m2(2860.0, gravity, radius) == pytest.approx(expected, rel=1e-3)
    assert expected == pytest.approx(1.18e9, rel=0.02)         # measured 2026-11-20 (reports/gen1698-h2s-so2-sinks-2026-11-19) for the 2,860 K body of seed 137


def test_the_melt_takes_what_the_crust_could_not_and_the_four_reservoirs_sum_to_the_sulfur():
    crust = sb.crust_column_kg_m2(3.49e6 * sb.EARTH_MANTLE_KG, 6.371e6)
    degassed = 1.0e9
    plain = sb.balance(degassed, age_s=AGE_S, ocean_present=False, oxidized_share=0.5, crust_kg_m2=crust, lifetimes=LIFETIMES)
    melted = sb.balance(degassed, age_s=AGE_S, ocean_present=False, oxidized_share=0.5, crust_kg_m2=crust, lifetimes=LIFETIMES, melt_capacity_mol_m2=3.0e8)
    assert plain.melt_mol_m2 == 0.0
    assert melted.melt_mol_m2 == pytest.approx(3.0e8)
    assert melted.air_mol_m2 == pytest.approx(plain.air_mol_m2 - 3.0e8)
    assert melted.sulfate_mol_m2 == plain.sulfate_mol_m2 and melted.sulfide_mol_m2 == plain.sulfide_mol_m2       # the crust sinks come first
    assert melted.air_mol_m2 + melted.sulfate_mol_m2 + melted.sulfide_mol_m2 + melted.melt_mol_m2 == pytest.approx(degassed)


def test_the_melt_takes_all_the_excess_when_its_capacity_is_larger_and_nothing_when_there_is_none():
    crust = sb.crust_column_kg_m2(3.49e6 * sb.EARTH_MANTLE_KG, 6.371e6)
    degassed = 1.0e9
    big = sb.balance(degassed, age_s=AGE_S, ocean_present=False, oxidized_share=0.5, crust_kg_m2=crust, lifetimes=LIFETIMES, melt_capacity_mol_m2=1.0e12)
    assert big.air_mol_m2 < 1.0e-3 * degassed         # what is left is the crust sinks' burden, not the excess
    assert big.melt_mol_m2 < degassed
    none = sb.balance(1.0e3, age_s=AGE_S, ocean_present=False, oxidized_share=0.5, crust_kg_m2=crust, lifetimes=LIFETIMES, melt_capacity_mol_m2=1.0e12)
    assert none.melt_mol_m2 == 0.0       # the crust holds all of it: the melt is for the excess only


def test_a_hot_body_in_the_solver_keeps_less_sulfur_gas_with_the_melt_than_without(monkeypatch):
    """The lava world of seed 137 (2,860 K, 5.25e8 mol/m2 degassed, crust capacity 2.3e7): H2S and SO2 stood at 30,000 kPa before the melt."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import atmosphere
    from planetgen.util import draw

    def sulfur_gas_kpa():
        draw.set_run_seed(137)
        config = SystemConfig()
        config.PLANETS = True
        system = StarSystem(system_config=config)
        body = next(b for p in system.all_planets if p.body_type != "a" for b in [p, *(p.moons or ())] if b.name == "Iltauris I")
        return atmosphere.partial_kpa(body, "h2s") + atmosphere.partial_kpa(body, "so2")

    with_melt = sulfur_gas_kpa()
    monkeypatch.setattr(ms, "capacity_mol_m2", lambda *a, **k: 0.0)
    without = sulfur_gas_kpa()
    assert without > 1000.0
    assert with_melt < 0.95 * without       # the melt takes about half of the sulfur and the body cools from 2,860 to about 1,900 K (the sulfur gas was its greenhouse), so the gas falls by a tenth
