"""
Where the degassed sulfur goes (GEN.815, `planetgen.physics.sulfur_balance`): the sulfate and sulfide balance. The Earth is a test
target, never an input of the code.
"""

import math

import pytest

from planetgen.physics import sulfur_balance as sb
from planetgen.physics import surface_pressure as sp

AGE_S = 4.5e9 * sb.YEAR_S
EARTH_RADIUS_M = 6.371e6
EARTH_RELEASED_J = 3.49e6 * sb.EARTH_MANTLE_KG
"""float: The heat the Earth twin of the interior model has released: 3.49e6 J per kilogram of its mantle (`surface_pressure`)."""
LIFETIMES = sb.Lifetimes(wet_s=1.0 * sb.DAY_S, dry_s=2.0e6 * sb.YEAR_S)


def _crust():
    return sb.crust_column_kg_m2(EARTH_RELEASED_J, EARTH_RADIUS_M)


def _balance(degassed=1.0e7, ocean=True, oxidized=0.9, crust=None, age=AGE_S):
    return sb.balance(degassed, age_s=age, ocean_present=ocean, oxidized_share=oxidized,
                      crust_kg_m2=_crust() if crust is None else crust, lifetimes=LIFETIMES)


def test_the_earth_twin_builds_the_crust_of_the_earth():
    area = 4.0 * math.pi * EARTH_RADIUS_M ** 2
    assert sb.crust_column_kg_m2(EARTH_RELEASED_J, EARTH_RADIUS_M) * area == pytest.approx(sb.EARTH_CRUST_KG, rel=0.01)
    assert sb.crust_column_kg_m2(0.0, EARTH_RADIUS_M) == 0.0


def test_a_hotter_interior_has_built_more_crust():
    assert sb.crust_column_kg_m2(2 * EARTH_RELEASED_J, EARTH_RADIUS_M) == pytest.approx(2 * _crust())


@pytest.mark.parametrize("degassed", [0.0, 1.0e3, 1.0e7, 1.0e9, 1.0e11])
@pytest.mark.parametrize("ocean", [True, False])
@pytest.mark.parametrize("oxidized", [0.0, 0.3, 1.0])
def test_the_sulfur_is_conserved(degassed, ocean, oxidized):
    result = _balance(degassed, ocean, oxidized)
    assert result.air_mol_m2 + result.sulfate_mol_m2 + result.sulfide_mol_m2 == pytest.approx(degassed, rel=1e-12, abs=1e-9)
    assert min(result.air_mol_m2, result.sulfate_mol_m2, result.sulfide_mol_m2) >= 0.0


def test_a_wet_world_with_a_crust_keeps_almost_no_sulfur_in_the_air():
    result = _balance(1.0e7, True)
    assert result.air_mol_m2 / 1.0e7 < 1.0e-9
    assert result.lifetime_s == LIFETIMES.wet_s


def test_a_dry_world_keeps_more_than_a_wet_one_for_the_longer_lifetime():
    wet, dry = _balance(1.0e7, True), _balance(1.0e7, False)
    assert dry.air_mol_m2 > wet.air_mol_m2
    assert dry.air_mol_m2 / 1.0e7 == pytest.approx(LIFETIMES.dry_s / AGE_S, rel=1e-6)


def test_a_body_younger_than_the_lifetime_keeps_all_its_sulfur_in_the_air():
    young = _balance(1.0e7, False, age=1.0e5 * sb.YEAR_S)
    assert young.air_mol_m2 == pytest.approx(1.0e7)
    assert _balance(1.0e7, True, age=0.0).air_mol_m2 == pytest.approx(1.0e7)


def test_without_a_crust_nothing_is_bound():
    result = _balance(1.0e7, True, crust=0.0)
    assert result.air_mol_m2 == pytest.approx(1.0e7)
    assert result.sulfate_mol_m2 == 0.0 and result.sulfide_mol_m2 == 0.0


def test_sulfur_beyond_the_capacity_of_the_sinks_stays_in_the_air_whole():
    capacity_sulfate, capacity_sulfide = sb.sink_capacities_mol_m2(_crust())
    degassed = 10.0 * (capacity_sulfate + capacity_sulfide)
    result = _balance(degassed, True, 0.5)
    assert result.air_mol_m2 >= degassed - capacity_sulfate - capacity_sulfide
    assert result.sulfate_mol_m2 == pytest.approx(capacity_sulfate, rel=1e-6)
    assert result.sulfide_mol_m2 == pytest.approx(capacity_sulfide, rel=1e-6)


def test_an_oxidized_gas_fills_the_sulfate_sink_and_a_reduced_one_the_sulfide_sink():
    oxidized, reduced = _balance(1.0e7, True, 1.0), _balance(1.0e7, True, 0.0)
    assert oxidized.sulfide_mol_m2 == 0.0 and oxidized.sulfate_mol_m2 > 0.99e7
    assert reduced.sulfate_mol_m2 == 0.0 and reduced.sulfide_mol_m2 > 0.99e7


def test_a_sink_that_is_full_does_not_take_the_other_sinks_share():
    capacity_sulfate, _ = sb.sink_capacities_mol_m2(_crust())
    result = _balance(10.0 * capacity_sulfate, True, 1.0)
    assert result.sulfide_mol_m2 == 0.0
    assert result.air_mol_m2 > 8.0 * capacity_sulfate


def test_the_capacity_is_one_calcium_and_one_iron_per_sulfur():
    sulfate, sulfide = sb.sink_capacities_mol_m2(1000.0)
    assert sulfate == pytest.approx(1000.0 * sb.CRUST_CALCIUM_WT / 0.040078, rel=1e-3)
    assert sulfide == pytest.approx(1000.0 * sb.CRUST_IRON_WT / 0.055845, rel=1e-3)


def test_the_oxidized_share_reads_the_gas():
    assert sb.oxidized_share_of({"SO2": 0.002, "H2S": 0.002}) == pytest.approx(0.5)
    assert sb.oxidized_share_of({"SO2": 1.0e-3, "S2": 1.0e-3}) == pytest.approx(1.0 / 3.0)
    assert sb.oxidized_share_of({"SO3": 1.0, "CO2": 5.0}) == 1.0
    assert sb.oxidized_share_of({"H2S": 1.0e-6, "OCS": 1.0e-6, "N2": 1.0}) == 0.0
    assert sb.oxidized_share_of({"N2": 1.0}) == 1.0


def test_the_lifetime_draws_stay_in_their_brackets_and_are_recorded():
    recorded = {}
    life = sb.draw_lifetimes(recorded)
    assert sb.WET_LIFETIME_DAYS[0] * sb.DAY_S <= life.wet_s <= sb.WET_LIFETIME_DAYS[1] * sb.DAY_S
    assert sb.DRY_LIFETIME_YEARS[0] * sb.YEAR_S <= life.dry_s <= sb.DRY_LIFETIME_YEARS[1] * sb.YEAR_S
    assert recorded == pytest.approx({"sulfur_lifetime_wet_days": life.wet_s / sb.DAY_S, "sulfur_lifetime_dry_years": life.dry_s / sb.YEAR_S}, rel=1e-12)


def test_the_crust_of_the_earth_twin_binds_the_sulfur_of_a_bulk_silicate_earth_degassing():
    """The mantle's 250 ppm of sulfur with the degassed share of `surface_pressure` fits in the sinks."""
    mantle_sulfur_mol_m2 = 250.0e-6 * sb.EARTH_MANTLE_KG / 0.032065 / (4.0 * math.pi * EARTH_RADIUS_M ** 2)
    degassed = mantle_sulfur_mol_m2 * sp.degassed_share(EARTH_RELEASED_J, sb.EARTH_MANTLE_KG) * sp.DEGAS_EFFICIENCY["S"]
    capacity_sulfate, capacity_sulfide = sb.sink_capacities_mol_m2(_crust())
    assert degassed < capacity_sulfate + capacity_sulfide


def test_the_dry_lifetime_range_holds_the_two_published_figures():
    """GEN.1562 N1: 1.9 Myr [R1343] and 20 Myr, quoted by [R0157], at Venus's surface are the bracket starts within a tenth of the first and ends at the second."""
    low, high = sb.DRY_LIFETIME_YEARS[:2]
    assert low <= 1.9e6 * 1.1 and high >= 2.0e7


def test_a_cooler_dry_surface_sulfates_slower_by_the_activation_energy_of_the_rate_law():
    """GEN.1562 N2: tau(T) = tau(735 K) exp(15248 (1/T - 1/735)): 5e5 times at 450 K, the same at 735 K, a hotter surface shorter."""
    tau = 2.0e6 * sb.YEAR_S
    assert sb.sulfation_lifetime_s(tau, 735.0) == pytest.approx(tau)
    assert sb.sulfation_lifetime_s(tau, 450.0) / tau == pytest.approx(math.exp(15248.0 * (1 / 450.0 - 1 / 735.0)), rel=1e-12)
    assert 4.0e5 < sb.sulfation_lifetime_s(tau, 450.0) / tau < 6.0e5
    assert sb.sulfation_lifetime_s(tau, 900.0) < tau
    assert sb.sulfation_lifetime_s(tau, None) == tau
    assert math.isfinite(sb.sulfation_lifetime_s(tau, 3.0)), "the exponent is held, no overflow at a few kelvin"


def test_a_cool_dry_body_keeps_its_oxidized_sulfur_in_the_air():
    hot = sb.balance(1.0e7, age_s=AGE_S, ocean_present=False, oxidized_share=1.0, crust_kg_m2=_crust(), lifetimes=LIFETIMES, temperature_k=735.0)
    cool = sb.balance(1.0e7, age_s=AGE_S, ocean_present=False, oxidized_share=1.0, crust_kg_m2=_crust(), lifetimes=LIFETIMES, temperature_k=320.0)
    assert cool.air_mol_m2 > 100.0 * hot.air_mol_m2
    assert cool.air_mol_m2 + cool.sulfate_mol_m2 + cool.sulfide_mol_m2 == pytest.approx(1.0e7)


def test_water_removes_at_most_one_mole_of_sulfur_per_mole_of_water_and_the_rest_meets_the_dry_lifetime():
    """GEN.1562 N6: a thin layer of water (no ocean) acid-izes only the sulfur it has water for; at 2.75 m and over the wet lifetime holds for all."""
    crust = _crust()
    kw = dict(age_s=AGE_S, ocean_present=False, oxidized_share=1.0, crust_kg_m2=crust, lifetimes=LIFETIMES, temperature_k=320.0)
    bound = _balance(1.0e7, False, oxidized=1.0).sulfate_mol_m2 + _balance(1.0e7, False, oxidized=1.0).air_mol_m2
    dry = sb.balance(1.0e7, **kw)
    thin = sb.balance(1.0e7, water_layer_m=0.05, **kw)                     # 0.05 m of water: 2.8e3 mol/m2 of water
    wet = sb.balance(1.0e7, water_layer_m=sb.WATER_LAYER_FULL_M, **kw)
    assert wet.air_mol_m2 < thin.air_mol_m2 < dry.air_mol_m2
    assert thin.sulfate_mol_m2 + thin.air_mol_m2 == pytest.approx(dry.sulfate_mol_m2 + dry.air_mol_m2)
    assert dry.air_mol_m2 - thin.air_mol_m2 <= 0.05 * 1000.0 / sb.WATER_MOLAR_MASS_KG * 1.0001, "no more sulfur is taken out than there is water"
    assert bound > 0.0


def test_a_water_layer_changes_nothing_for_an_ocean_world_and_none_keeps_the_old_rule():
    kw = dict(age_s=AGE_S, oxidized_share=0.9, crust_kg_m2=_crust(), lifetimes=LIFETIMES)
    assert sb.balance(1.0e7, ocean_present=True, water_layer_m=0.001, **kw).air_mol_m2 == pytest.approx(sb.balance(1.0e7, ocean_present=True, **kw).air_mol_m2)
    assert sb.balance(1.0e7, ocean_present=False, **kw).air_mol_m2 == pytest.approx(_balance(1.0e7, False).air_mol_m2)


def test_the_redox_shift_is_applied_to_the_magma_ocean_and_the_later_mantle():
    """GEN.1562 N4: `mixture_rocky.draw_redox` returns the recipe's first value plus `REDOX_SHIFT_DEX`, and the building blocks' and the white dwarf debris's values are not moved."""
    from types import SimpleNamespace
    from planetgen.physics import mixture_rocky
    recipe = lambda **values: SimpleNamespace(values=values)           # noqa: E731
    subject = {"age_gy": 4.5, "long_lived_geology": False}
    assert mixture_rocky.draw_redox(subject, recipe(a1=0.5), {}) == pytest.approx(0.5 + mixture_rocky.REDOX_SHIFT_DEX)
    assert mixture_rocky.draw_redox(subject, recipe(a2=3.0), {}) == pytest.approx(3.0 + mixture_rocky.REDOX_SHIFT_DEX)
    assert mixture_rocky.draw_redox(subject, recipe(a3=-1.0), {}) == -1.0
    assert mixture_rocky.draw_redox(subject, recipe(a4=-1.2), {}) == -1.2
