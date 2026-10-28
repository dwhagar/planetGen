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
EARTH_RELEASED_J = 4.09e6 * sb.EARTH_MANTLE_KG
"""float: The heat the Earth twin of the interior model has released: 4.09e6 J per kilogram of its mantle (`surface_pressure`)."""
LIFETIMES = sb.Lifetimes(wet_s=1.0 * sb.DAY_S, dry_s=2.0e6 * sb.YEAR_S)


def _crust():
    return sb.crust_column_kg_m2(EARTH_RELEASED_J, EARTH_RADIUS_M)


def _balance(degassed=1.0e7, ocean=True, oxidised=0.9, crust=None, age=AGE_S):
    return sb.balance(degassed, age_s=age, ocean_present=ocean, oxidised_share=oxidised,
                      crust_kg_m2=_crust() if crust is None else crust, lifetimes=LIFETIMES)


def test_the_earth_twin_builds_the_crust_of_the_earth():
    area = 4.0 * math.pi * EARTH_RADIUS_M ** 2
    assert sb.crust_column_kg_m2(EARTH_RELEASED_J, EARTH_RADIUS_M) * area == pytest.approx(sb.EARTH_CRUST_KG, rel=0.01)
    assert sb.crust_column_kg_m2(0.0, EARTH_RADIUS_M) == 0.0


def test_a_hotter_interior_has_built_more_crust():
    assert sb.crust_column_kg_m2(2 * EARTH_RELEASED_J, EARTH_RADIUS_M) == pytest.approx(2 * _crust())


@pytest.mark.parametrize("degassed", [0.0, 1.0e3, 1.0e7, 1.0e9, 1.0e11])
@pytest.mark.parametrize("ocean", [True, False])
@pytest.mark.parametrize("oxidised", [0.0, 0.3, 1.0])
def test_the_sulfur_is_conserved(degassed, ocean, oxidised):
    result = _balance(degassed, ocean, oxidised)
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


def test_an_oxidised_gas_fills_the_sulfate_sink_and_a_reduced_one_the_sulfide_sink():
    oxidised, reduced = _balance(1.0e7, True, 1.0), _balance(1.0e7, True, 0.0)
    assert oxidised.sulfide_mol_m2 == 0.0 and oxidised.sulfate_mol_m2 > 0.99e7
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


def test_the_oxidised_share_reads_the_gas():
    assert sb.oxidised_share_of({"SO2": 0.002, "H2S": 0.002}) == pytest.approx(0.5)
    assert sb.oxidised_share_of({"SO2": 1.0e-3, "S2": 1.0e-3}) == pytest.approx(1.0 / 3.0)
    assert sb.oxidised_share_of({"SO3": 1.0, "CO2": 5.0}) == 1.0
    assert sb.oxidised_share_of({"H2S": 1.0e-6, "OCS": 1.0e-6, "N2": 1.0}) == 0.0
    assert sb.oxidised_share_of({"N2": 1.0}) == 1.0


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
