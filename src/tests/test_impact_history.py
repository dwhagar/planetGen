"""
Late accretion, impact delivery and impact erosion of an air (GEN.739, `planetgen.physics.impact_history`,
docs/design/composition/impact-erosion-and-late-delivery.md).
"""

import math

import pytest

from planetgen.physics import impact_history as ih, surface_pressure as sp
from planetgen.util import draw

EARTH = dict(mass_kg=5.97e24, radius_m=6.371e6, temperature_k=288.0, orbital_speed_ms=29780.0)
VENUS = dict(mass_kg=4.87e24, radius_m=6.052e6, temperature_k=2500.0, orbital_speed_ms=35020.0)
MARS = dict(mass_kg=6.39e23, radius_m=3.39e6, temperature_k=223.0, orbital_speed_ms=24070.0)
M_EARTH = 5.97e24

# the surface limits of [R3872] section 2.3 (pascals)
EARTH_CAPS = {"H2O": 1.7e-2 * 1e5, "CO2": 1.0 * 1e5}
MARS_CAPS = {"H2O": 3.9e-6 * 1e5, "CO2": 3.0 * 1e5}

EARTH_N2_AIR_KG = 3.9e18          # [R2923] (read, surface-pressure-supply-law.md)
EARTH_OCEAN_KG = 1.4e21
VENUS_N2_AIR_KG = 1.08e19         # 4.77e20 kg of air, 3.5 percent N2 by volume
MARS_N2_AIR_KG = 4.4e14           # 2.5e16 kg of air, 2.7 percent N2 by volume


def _balance(planet, scale, caps=None, fraction=None):
    mean = ih.late_accreted_percent_mean(planet["mass_kg"] / M_EARTH) / 100.0
    volatiles = {s: share * scale for s, share in ih.VENEER_VOLATILES.items()}
    return ih.veneer_balance(accreted_fraction=fraction or mean, volatiles=volatiles, caps_pa=caps, **planet)


def test_the_mean_late_accreted_mass_follows_the_four_means_of_brasser_et_al_2016_table_1_inside_their_scatter():
    table = {"Mercury": (0.055, 2.7, 4.5), "Venus": (0.815, 0.8, 0.5), "Earth": (1.0, 0.7, 0.5), "Mars": (0.107, 1.4, 2.7)}
    for name, (mass, percent, sd) in table.items():
        assert abs(ih.late_accreted_percent_mean(mass) - percent) < sd, name
    assert ih.late_accreted_percent_mean(1.0) == pytest.approx(0.7)
    assert 0.5 <= ih.late_accreted_percent_mean(1.0) <= 0.8, "[R4598]: the highly siderophile elements ask 0.5 to 0.8 percent"


def test_the_drawn_fraction_is_recorded_scattered_and_capped():
    constants = {}
    with draw.bound(7):
        fractions = [ih.late_accreted_fraction(1.0, constants) for _ in range(300)]
    assert constants["late_accreted_percent"] <= ih.LATE_ACCRETED_PERCENT_MAX
    assert max(fractions) <= ih.LATE_ACCRETED_PERCENT_MAX / 100.0 and min(fractions) > 0.0
    assert math.exp(sum(math.log(f) for f in fractions) / len(fractions)) == pytest.approx(0.007, rel=0.2)


def test_the_impact_speed_is_the_published_one_for_the_inner_planets():
    """[R0518: Table 1]: Venus 24.3, Earth 21.1, Mars 13.4 km/s (0.6 of the orbital speed for the approach)."""
    for planet, published in ((VENUS, 24.3), (EARTH, 21.1), (MARS, 13.4)):
        speed = ih.impact_speed_ms(planet["mass_kg"], planet["radius_m"], planet["orbital_speed_ms"]) / 1e3
        assert speed == pytest.approx(published, rel=0.16)


def test_a_thin_air_loses_less_per_impactor_mass_and_a_big_impactor_less_than_a_small_one():
    thick = ih.erosion_per_impactor_mass(1.2, 8500.0, 21e3, 11.2e3)
    thin = ih.erosion_per_impactor_mass(0.12, 8500.0, 21e3, 11.2e3)
    assert thick / 100.0 < thin < thick / 3.0, "at least linear in the density (a thick air also slows the small impactors)"
    assert ih.erosion_per_impactor_mass(0.0, 8500.0, 21e3, 11.2e3) == 0.0
    small = ih.erosion_efficiency(1.2, 8500.0, 1e3, 21e3, 11.2e3)
    large = ih.erosion_efficiency(1.2, 8500.0, 1e6, 21e3, 11.2e3)
    assert small > large > 0.0
    # Svetsov's leading term for a 1000 km body hitting the present Earth, 3 rho0 / (4 rho_imp) * 2H / D, times the angle factor and the share
    # of the plume that outruns the escape speed (x = v_esc / (5.1 v), 5.1 = sqrt(4 gamma / (gamma - 1)))
    plume = ih._plume_fraction(11.2e3 / (21e3 * math.sqrt(26.0)))
    assert large == pytest.approx(3 * 1.2 / (4 * 3320.0) * 2 * 8500.0 / 1e6 * ih.ANGLE_FACTOR * plume, rel=0.1)


def test_a_slow_impact_that_cannot_launch_the_plume_past_escape_erodes_nothing():
    assert ih.erosion_efficiency(1.2, 8500.0, 1e5, 1e3, 11.2e3) == pytest.approx(0.0, abs=1e-12)


def test_the_balance_conserves_mass_and_keeps_every_partial_pressure_under_its_cap():
    result = _balance(EARTH, 0.3, EARTH_CAPS)
    area = 4 * math.pi * EARTH["radius_m"] ** 2
    gravity = ih.G * EARTH["mass_kg"] / EARTH["radius_m"] ** 2
    for species in ih.SPECIES:
        total = result["air_kg"][species] + result["kept_kg"][species] + result["lost_kg"][species]
        assert total == pytest.approx(result["delivered_kg"][species], rel=1e-9), species
        cap = EARTH_CAPS.get(species)
        if cap is not None:
            assert result["air_kg"][species] * gravity / area <= cap * (1 + 1e-9)


def test_with_no_impactor_volatiles_nothing_arrives_and_nothing_is_eroded_from_an_empty_air():
    result = ih.veneer_balance(accreted_fraction=0.007, volatiles={}, **EARTH)
    assert all(v == 0.0 for kind in ("air_kg", "kept_kg", "lost_kg") for v in result[kind].values())


def test_an_existing_air_is_eroded_toward_what_the_impactors_replace():
    start = {"N2": 5e19}
    result = ih.veneer_balance(accreted_fraction=0.007, volatiles={"N2": 0.0002}, air_kg=start, **EARTH)
    assert result["air_kg"]["N2"] < start["N2"] / 5.0
    assert result["lost_kg"]["N2"] > 0.8 * start["N2"]


def test_the_earth_twin_keeps_an_n2_air_that_the_observed_one_lies_inside_the_bracket_of_the_veneer_scale():
    """[R3872] and this module: the veneer scale 0.1 (enstatite-like) to 1 (carbonaceous) brackets the Earth's N2."""
    low, high = (_balance(EARTH, scale, EARTH_CAPS)["air_kg"]["N2"] for scale in (0.1, 1.0))
    assert low < EARTH_N2_AIR_KG < high
    assert high / low > 5.0


def test_the_earth_twin_keeps_a_share_of_the_ocean_in_the_water_the_impactors_brought():
    low, high = (_balance(EARTH, scale, EARTH_CAPS)["kept_kg"]["H2O"] for scale in (0.1, 1.0))
    assert low < EARTH_OCEAN_KG < high, "the veneer alone is a few tenths to a few times the ocean (an enstatite-like one is the lower part)"


def test_venus_with_no_ocean_and_no_carbonate_sink_holds_the_whole_delivery_that_the_impacts_did_not_erode():
    venus = _balance(VENUS, 1.0)
    assert all(v == 0.0 for v in venus["kept_kg"].values())
    assert venus["air_kg"]["N2"] > 0.2 * VENUS_N2_AIR_KG, "within a factor of five of the N2 the planet has"
    assert venus["air_kg"]["N2"] < 5.0 * VENUS_N2_AIR_KG


def test_mars_would_hold_far_more_n2_than_it_has_so_the_thin_air_is_the_escape_of_the_retention_rows_not_the_impacts():
    """The check [R3872] report: impacts alone leave Mars with a larger share of N2 than Venus or the Earth."""
    mars = _balance(MARS, 0.1, MARS_CAPS)["air_kg"]["N2"]
    assert mars > 100.0 * MARS_N2_AIR_KG
    assert mars / MARS["mass_kg"] > _balance(EARTH, 0.1, EARTH_CAPS)["air_kg"]["N2"] / EARTH["mass_kg"]


def test_the_impactor_vapour_that_escapes_is_a_tenth_to_a_fifth_for_the_inner_planets():
    for planet in (EARTH, VENUS, MARS):
        zeta = _balance(planet, 0.1, fraction=1e-4)["escape_share"]
        assert 0.05 < zeta < 0.25


def test_the_supply_law_takes_the_late_accretion_as_atoms_per_100_g_of_the_mantle():
    result = _balance(EARTH, 0.1, EARTH_CAPS)
    mantle_kg = 4.0e24
    atoms = sp.impact_inventory(result, mantle_kg)
    held_n2 = result["air_kg"]["N2"] + result["kept_kg"]["N2"]
    assert atoms["N"] == pytest.approx(2.0 * held_n2 * 1000.0 / 28.013 * 0.1 / mantle_kg)
    held_h2o = result["air_kg"]["H2O"] + result["kept_kg"]["H2O"]
    assert atoms["H"] == pytest.approx(2.0 * held_h2o * 1000.0 / 18.015 * 0.1 / mantle_kg)
    assert sp.impact_inventory(result, 0.0) == {"H": 0.0, "C": 0.0, "N": 0.0}


def test_the_veneer_scale_is_drawn_inside_its_bracket_and_recorded():
    constants = {}
    with draw.bound(3):
        volatiles = ih.veneer_volatiles(constants)
    scale = constants["veneer_volatile_scale"]
    assert ih.VENEER_VOLATILE_SCALE[0] <= scale <= ih.VENEER_VOLATILE_SCALE[1]
    assert volatiles["N2"] == pytest.approx(ih.VENEER_VOLATILES["N2"] * scale)
