# tests/test_comet_ice_budget.py

"""GEN.595: the ices of a comet from the zone of the disc it formed in (`physics/comet_ice_budget.py`; Oberg, Murray-Clay and Bergin 2011)."""

import math

import pytest

from planetgen.physics import comet_ice_budget as cib

SOLAR = cib.SOLAR_LUMINOSITY_W


def test_the_snow_lines_of_a_sun_like_star_are_those_oberg_et_al_give():
    """Their text puts the CO2 and CO lines at 10 and 40 AU for the average disc."""
    assert cib.ice_line_au("H2O") == pytest.approx(1.9, abs=0.05)
    assert cib.ice_line_au("CO2") == pytest.approx(10.3, abs=0.2)
    assert cib.ice_line_au("CO") == pytest.approx(41.0, abs=0.5)
    assert cib.disk_temperature_k(cib.ice_line_au("CO2")) == pytest.approx(cib.ICE_LINES_K["CO2"])


def test_a_more_luminous_star_pushes_every_line_out_by_the_same_factor():
    factor = (4.0) ** (cib.LUMINOSITY_EXPONENT / cib.DISK_INDEX)
    for ice in cib.ICE_LINES_K:
        assert cib.ice_line_au(ice, 4.0 * SOLAR) == pytest.approx(factor * cib.ice_line_au(ice))
    assert factor > 1.0


def test_the_nebula_holds_a_third_as_much_co2_and_one_and_two_thirds_as_much_co_as_water():
    co2, co = cib.nebular_ratios()
    assert co2 == pytest.approx(0.3 / 0.9) and co == pytest.approx(1.5 / 0.9)


def test_carbon_scales_the_carbon_ices_and_oxygen_the_water():
    co2, co = cib.nebular_ratios()
    assert cib.nebular_ratios(carbon=2.0) == pytest.approx((2.0 * co2, 2.0 * co))
    rock = cib.condensate_shares.ROCK_OXYGEN_SHARE
    assert cib.nebular_ratios(oxygen=1.0 + (1.0 - rock) * 0.5)[0] == pytest.approx(co2 / 1.5)
    with pytest.raises(ValueError):
        cib.nebular_ratios(oxygen=rock)            # all the oxygen in rock: no water to be a ratio of
    with pytest.raises(ValueError):
        cib.nebular_ratios(carbon=-1.0)


def test_each_zone_keeps_what_its_temperature_lets_it():
    co2, co = cib.nebular_ratios()
    trapped = cib.TRAPPED_FRACTION
    assert cib.zone_of(1.0) == "dry" and cib.ice_ratios(1.0) is None
    assert cib.zone_of(5.0) == "water" and cib.ice_ratios(5.0) == pytest.approx((trapped * co2, trapped * co))
    assert cib.zone_of(20.0) == "carbon dioxide" and cib.ice_ratios(20.0) == pytest.approx((co2, trapped * co))
    assert cib.zone_of(60.0) == "carbon monoxide" and cib.ice_ratios(60.0) == pytest.approx((co2, co))
    assert cib.ice_ratios(20.0, trapped=0.0) == pytest.approx((co2, 0.0))
    with pytest.raises(ValueError):
        cib.ice_ratios(20.0, trapped=1.5)


def test_the_ratios_step_up_outward_at_each_line():
    previous = (0.0, 0.0)
    for distance in (3.0, 15.0, 60.0):
        now = cib.ice_ratios(distance)
        assert now[0] >= previous[0] and now[1] >= previous[1]
        previous = now


def test_equal_masses_form_in_equal_intervals_of_ln_r():
    shares = cib.zone_shares(5.0, 50.0)
    assert sum(shares.values()) == pytest.approx(1.0) and shares["dry"] == 0.0
    assert shares["water"] == pytest.approx(math.log(cib.ice_line_au("CO2") / 5.0) / math.log(10.0))
    assert shares["carbon monoxide"] == pytest.approx(math.log(50.0 / cib.ice_line_au("CO")) / math.log(10.0))
    wide = cib.zone_shares(0.5, 500.0)
    assert sum(wide.values()) == pytest.approx(1.0) and all(share > 0.0 for share in wide.values())
    with pytest.raises(ValueError):
        cib.zone_shares(10.0, 5.0)


def test_borisov_is_a_comet_of_the_carbon_monoxide_zone():
    """2I/Borisov had CO/H2O of 130 to 155 percent (Bodewits et al. 2020): the zone beyond the CO line holds 167 percent for a Sun-like budget, to within 30 percent."""
    _, co = cib.ice_ratios(100.0)
    for measured in (1.30, 1.55):
        assert abs(math.log(co / measured)) < math.log(1.3)


def test_the_median_comet_of_the_solar_system_lies_between_the_zones_not_in_one():
    """CO2/H2O 10.9 and CO/H2O 2.5 percent are under the CO2 zone's (0.5 and 0.7 dex) and over the warm zone's (CO2 8 times under; CO 5 times over): the model's absolute level for the
    comets of the Solar System is the stand-in that stays, and the order (more CO2 than CO) is the CO2 zone's, not the water zone's."""
    median_co2, median_co = 0.109, 0.025
    co2_zone, warm_zone = cib.ice_ratios(20.0), cib.ice_ratios(5.0)
    assert 0.4 < math.log10(co2_zone[0] / median_co2) < 0.6 and 0.6 < math.log10(co2_zone[1] / median_co) < 0.8
    assert warm_zone[0] < median_co2 and warm_zone[1] > median_co
    assert median_co2 > median_co and co2_zone[0] > co2_zone[1] and warm_zone[0] < warm_zone[1]


def test_the_co2_of_3i_atlas_is_beyond_any_zone_of_this_budget():
    """CO2/H2O 7.6 (Cordiner et al. 2025) is over twenty times the most carbon and oxygen of a solar star give."""
    assert 7.6 / cib.nebular_ratios()[0] > 20.0
    assert 7.6 / cib.nebular_ratios(carbon=2.0)[0] > 10.0


def test_a_comet_forms_between_the_snow_line_and_the_edge_of_the_disc_in_equal_masses_per_interval_of_ln_r():
    """GEN.595 wiring: the uniform draw is a position in ln r, so u = 0 is the inner edge (the water zone) and u near 1 the outer (the CO zone)."""
    inner = cib.formation_ratios(SOLAR, 2.7, 50.0, 1.0, 1.0, 0.0)
    outer = cib.formation_ratios(SOLAR, 2.7, 50.0, 1.0, 1.0, 0.999)
    assert inner == cib.ice_ratios(2.7, SOLAR)
    assert outer == cib.ice_ratios(50.0 ** 0.999 * 2.7 ** 0.001, SOLAR) and outer[1] > inner[1]
    shares = cib.zone_shares(2.7, 50.0, SOLAR)
    hits = [cib.zone_of(2.7 * (50.0 / 2.7) ** ((i + 0.5) / 4000), SOLAR) for i in range(4000)]
    for zone in ("water", "carbon dioxide", "carbon monoxide"):
        assert hits.count(zone) / 4000 == pytest.approx(shares[zone], abs=0.01)


def test_a_disc_that_ends_inside_the_water_line_still_gives_the_ice_at_its_inner_edge():
    assert cib.formation_ratios(SOLAR, 1.0, 1.5, 1.0, 1.0, 0.5) is not None
    assert cib.formation_ratios(SOLAR, 1.0, 0.5, 1.0, 1.0, 0.5) is not None


def test_a_formation_draw_must_be_a_fraction():
    with pytest.raises(ValueError):
        cib.formation_ratios(SOLAR, 2.7, 50.0, 1.0, 1.0, 1.0)
