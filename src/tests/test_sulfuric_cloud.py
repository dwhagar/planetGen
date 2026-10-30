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


def test_the_base_is_cooler_than_the_base_temperature_and_below_the_top(skin):
    base = sc.base_pressure_pa(column(), 735.0, skin)
    assert 1.0e4 < base < VENUS_PRESSURE_PA
    pressure, temperature = aerosols.profile(column(), 735.0, skin_k=skin)
    assert temperature[pressure >= base][-1] <= sc.number("base_temperature_k") + 1.0


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
