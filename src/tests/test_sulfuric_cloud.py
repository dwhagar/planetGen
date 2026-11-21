# tests/test_sulfuric_cloud.py

"""GEN.1019: the sulfuric acid cloud of an air (`planetgen.physics.sulfuric_cloud`): the supply gates, the base, the limit and the draw."""

import json

import pytest

from planetgen.physics import acid_medium, air_column, aerosols
from planetgen.physics import sulfuric_cloud as sc

VENUS_PRESSURE_PA = 92.0e5
VENUS_FLUX_W_M2 = 2600.0
"""float: The sunlight at the top of the air, for the cold top (a test air, not stored anywhere)."""


def column(so2=150.0e-6, h2o=30.0e-6, pressure=VENUS_PRESSURE_PA, gravity=8.87):
    parts = {"CO2": 0.965 * pressure, "N2": 0.035 * pressure}
    if so2:
        parts["SO2"] = so2 * pressure
    if h2o:
        parts["H2O"] = h2o * pressure
    return air_column.from_partials(parts, gravity)


@pytest.fixture(scope="module")
def skin():
    return sc.skin_k(VENUS_FLUX_W_M2)


def test_a_thick_air_with_sulfur_dioxide_and_a_trace_of_water_gets_a_cloud_of_the_published_depth(skin):
    for tau_draw in (0.0, 0.5, 1.0):
        tau, cover = sc.draw(column(), 735.0, tau_draw, 0.5, skin)
        assert sc.number("optical_depth_low") <= tau <= sc.number("optical_depth_high") and cover == sc.number("cover")


def test_an_air_with_no_sulfur_gets_no_acid_cloud(skin):
    assert sc.draw(column(so2=0.0), 735.0, 0.5, 0.5, skin) == (0.0, 0.0)
    assert sc.acid_column_limit_kg_m2(column(so2=0.0), 735.0, skin) == 0.0


def test_an_air_without_water_or_without_an_oxidiser_gets_none(skin):
    assert sc.draw(column(h2o=0.0), 735.0, 0.5, 0.5, skin) == (0.0, 0.0)
    no_oxidant = air_column.from_partials({"N2": 1.0e5, "SO2": 1.0e-3 * 1.0e5, "H2O": 1.0e-3 * 1.0e5}, 9.81)
    assert not sc.supplied(no_oxidant) and sc.draw(no_oxidant, 300.0, 0.5, 0.5, skin) == (0.0, 0.0)


def test_the_supply_gates_are_those_of_the_acid_medium():
    assert sc.supplied(column(so2=2.0 * acid_medium.number("so2_min_share")))
    assert not sc.supplied(column(so2=0.5 * acid_medium.number("so2_min_share")))


def test_the_base_is_the_warmest_level_where_the_acid_condenses_and_below_the_top(skin):
    found = sc.window(column(), 735.0, skin)
    assert found is not None and 345.0 <= found.base_temperature_k <= 372.0, "the published base is 345 to 367 K"
    assert sc.base_pressure_pa(column(), 735.0, skin) == found.base_pressure_pa
    assert 0.0 < found.top_pressure_pa < found.base_pressure_pa < VENUS_PRESSURE_PA
    pressure, temperature = aerosols.profile(column(), 735.0, skin_k=skin)
    assert temperature[pressure >= found.base_pressure_pa][-1] == pytest.approx(found.base_temperature_k, abs=1.0)


def test_a_wetter_air_condenses_a_weaker_acid_at_a_warmer_base(skin):
    wet = sc.window(column(h2o=3.0e-4), 735.0, skin)
    dry = sc.window(column(h2o=3.0e-6), 735.0, skin)
    assert wet is not None and dry is not None and wet.base_temperature_k > dry.base_temperature_k
    assert wet.weight_percent < dry.weight_percent


def test_the_infrared_cloud_carries_the_base_and_strength_of_its_window(skin):
    found = sc.infrared_cloud(column(), 735.0, 0.5, skin)
    window = sc.window(column(), 735.0, skin)
    assert found[2] == window.base_temperature_k and found[3] == window.weight_percent and found[0] > 0.0
    assert sc.infrared_cloud(column(so2=0.0), 735.0, 0.5, skin) is None


def test_an_air_that_stays_warm_has_no_base():
    """With no cold top (skin held at the ground) nothing is cooler than the base temperature."""
    assert sc.base_pressure_pa(column(), 735.0, 735.0) is None


def test_a_cool_ground_puts_the_base_at_the_ground():
    cool = column(pressure=1.0e5, gravity=9.81, so2=1.0e-4, h2o=1.0e-3)
    assert sc.base_pressure_pa(cool, 300.0) == pytest.approx(cool.pressure_pa)


def test_a_little_sulfur_limits_the_depth(skin):
    """The depth cannot exceed what the sulfur above the base makes: a thin air with a few parts per million gives a thin cloud."""
    thin = column(so2=2.0e-6, h2o=1.0e-3, pressure=1.0e4, gravity=9.81)
    tau, cover = sc.draw(thin, 300.0, 1.0, 0.0)
    limit = sc.optical_depth(sc.acid_column_limit_kg_m2(thin, 300.0), sc.number("radius_low_um"))
    assert 0.0 < tau == pytest.approx(limit) and tau < sc.number("optical_depth_low")
    assert sc.draw(column(so2=2.0e-6), 735.0, 0.0, 0.0, skin)[0] == sc.number("optical_depth_low")     # a thick air has sulfur to spare


def test_the_limit_is_the_sulfur_above_the_base_made_into_acid(skin):
    air = column()
    base = sc.base_pressure_pa(air, 735.0, skin)
    expected = air.fraction("SO2") * base / 8.87 * sc.ACID_MOLAR_MASS_KG_MOL / air.mean_molar_mass_kg_mol
    assert sc.acid_column_limit_kg_m2(air, 735.0, skin) == pytest.approx(expected)


def test_the_depth_of_a_mass_goes_inversely_with_the_droplet_radius():
    assert sc.optical_depth(0.03, 1.0) == pytest.approx(2.0 * sc.optical_depth(0.03, 2.0))
    assert sc.optical_depth(0.03, 1.0) == pytest.approx(3.0 * 2.0 * 0.03 / (4.0 * 1800.0 * 1.0e-6))


def test_the_published_venus_cloud_is_seven_to_a_few_hundred_grams_per_square_metre():
    """Not a target: 20 to 30 deep in droplets of 0.3 to 3.85 micron is 0.007 to 0.3 kg/m2, the order of the model loadings of Dai and others 2022."""
    low = 20.0 * 4.0 * 1800.0 * 0.3e-6 / (3.0 * 2.0)
    high = 30.0 * 4.0 * 1800.0 * 3.85e-6 / (3.0 * 2.0)
    assert 0.005 < low < 0.01 and 0.1 < high < 0.5


def test_every_number_has_a_tier_and_a_tweak_and_a_bad_file_is_refused(tmp_path, monkeypatch):
    for name, entry in sc.NUMBERS.items():
        assert entry["tier"].startswith("T") and entry["tweak"].startswith("L") and entry["basis"], name
    with open(sc.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    del document["numbers"]["cover"]["tweak"]
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    monkeypatch.setattr(sc, "DATA_FILE", str(path))
    with pytest.raises(sc.SulfuricCloudError):
        sc._load()


COLD_AIR = {"H2O": 0.0912, "CO2": 3459.2, "SO2": 0.5654, "H2S": 184.9, "N2": 5949.6, "Ar": 308.8, "OCS": 0.1072, "CH4": 53.12}
"""dict: A cold sulfurous air of the 5,000-system census (a class E world, the ground at 198 K, 9.9 kPa): its acid condenses at every level it can, and every one of them is under the freezing point of its droplets (GEN.1644)."""


def test_a_cloud_stands_only_where_its_droplets_are_above_their_freezing_point():
    """GEN.1644: the acid cloud is a liquid cloud, so a level whose droplets are at or below the freezing point of their strength is not in it; the cold air that the 185 K floor gave a cloud (the census body had a 23.5 deep one) has none, the Venus air keeps its cloud with a base above its freezing point."""
    import math

    from planetgen.physics import acid_freezing, sulfur_cloud_set, sulfuric_vapor
    cold = air_column.from_partials(COLD_AIR, 1.4559 * 9.80665)
    temperature = 197.95
    skin = aerosols.SKIN_FACTOR * temperature / (1.0 + 0.3549) ** 0.25
    pressure, level_k = aerosols.profile(cold, temperature, skin_k=skin)
    acid, water = cold.fraction("SO2") * sulfur_cloud_set.number("acid_vapor_share_of_so2"), cold.fraction("H2O")
    condensing = []
    for p, t in zip(pressure, level_k):
        water_pa = water * float(p)
        saturation = aerosols.saturation_pa("H2O", float(t))
        if math.isfinite(saturation):
            water_pa = min(water_pa, saturation)
        if sulfuric_vapor.condenses(float(t), acid * float(p), water_pa):
            condensing.append(acid_freezing.is_frozen(float(t), sulfuric_vapor.composition(float(t), water_pa).acid_fraction))
    assert condensing and all(condensing), "the acid condenses somewhere, and only as a solid"
    assert sc.window(cold, temperature, skin) is None
    assert sc.draw(cold, temperature, 0.5, 0.5, skin) == (0.0, 0.0)
    venus = sc.window(column(), 735.0, sc.skin_k(VENUS_FLUX_W_M2))
    assert venus is not None
    assert venus.base_temperature_k > acid_freezing.freezing_point_of_weight_percent(venus.weight_percent)
