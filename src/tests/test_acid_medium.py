"""
The acid medium reader (GEN.888): where a body holds liquid concentrated sulfuric acid.
"""

import copy
import json

import pytest

from planetgen.physics import acid_medium, aerosols, atmosphere

GAS = atmosphere.GAS_FIELD


def _body(surface_k, gravity_g, gases, flux=None):
    return {"surface_temperature_k": surface_k, "gravity_g": gravity_g, "energy_flux_w_m2": flux, GAS: dict(gases)}


def _venus(**share):
    total = 9.2e6
    mix = {"CO2": 0.965 * total, "N2": 0.035 * total, "SO2": 1.5e-4 * total, "H2O": 3e-5 * total}
    mix.update({formula: fraction * total for formula, fraction in share.items()})
    return _body(737.0, 0.904, mix, flux=2600.0)


def test_venus_like_air_holds_a_cloud_layer():
    found = acid_medium.medium(_venus())
    assert found.level == acid_medium.CLOUD and found.layer == "cloud"
    assert 283.0 <= found.temperature_k <= 373.0
    assert 30000.0 <= found.pressure_pa <= 300000.0
    assert found.weight_percent >= 81.0


def test_a_wet_venus_has_only_a_dilute_acid_so_none():
    assert acid_medium.level(_venus(H2O=3e-3)) == acid_medium.NONE


def test_earth_like_air_has_none():
    earth = _body(288.0, 1.0, {"N2": 78000.0, "O2": 21000.0, "H2O": 1500.0, "Ar": 900.0, "CO2": 40.0}, flux=1361.0)
    assert acid_medium.level(earth) == acid_medium.NONE


def test_a_warm_dry_sulfurous_ground_holds_a_pool():
    dry = _body(330.0, 1.0, {"CO2": 80000.0, "N2": 20000.0, "SO2": 1000.0, "H2O": 2.0})
    found = acid_medium.medium(dry)
    assert found.level == acid_medium.POOL and found.layer == "pool"
    assert found.temperature_k == 330.0 and found.weight_percent >= 81.0


@pytest.mark.parametrize("without", ["SO2", "H2O"])
def test_a_missing_supply_makes_none(without):
    body = _venus()
    del body[GAS][without]
    assert acid_medium.level(body) == acid_medium.NONE


def test_no_oxidiser_makes_none():
    body = _body(330.0, 1.0, {"N2": 100000.0, "SO2": 1000.0, "H2O": 2.0})
    assert acid_medium.level(body) == acid_medium.NONE


def test_an_airless_body_has_none():
    assert acid_medium.level(_body(300.0, 1.0, {})) == acid_medium.NONE
    assert acid_medium.level({"surface_temperature_k": None, "gravity_g": 1.0, GAS: {"CO2": 1.0}}) == acid_medium.NONE


def test_a_hot_dry_ground_is_outside_the_liquid_window():
    hot = _body(500.0, 1.0, {"CO2": 80000.0, "SO2": 1000.0, "H2O": 2.0})
    assert acid_medium.medium(hot).layer != "pool"


def test_strength_falls_as_the_air_wets():
    strengths = [acid_medium.weight_percent(activity) for activity in (1e-6, 1e-4, 1e-3, 1e-2, 1.0)]
    assert strengths == sorted(strengths, reverse=True)
    assert acid_medium.weight_percent(0.0) == 100.0 and acid_medium.weight_percent(1.0) >= 0.0
    assert acid_medium.weight_percent(acid_medium.number("reference_water_activity")) == pytest.approx(acid_medium.number("reference_wt_percent"))


def test_the_skin_argument_lowers_the_cold_top_only():
    from planetgen.physics import air_column, constants
    column = air_column.from_partials({"CO2": 9.2e6}, 0.904 * constants.EARTH_GRAVITY)
    _, default = aerosols.profile(column, 737.0)
    _, cold = aerosols.profile(column, 737.0, skin_k=200.0)
    assert default.min() == pytest.approx(737.0 * aerosols.SKIN_FACTOR) and cold.min() == pytest.approx(200.0)
    assert default[0] == cold[0] == 737.0


def test_every_number_has_a_tier_a_tweak_and_a_basis():
    with open(acid_medium.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    assert document["numbers"]
    for name, entry in document["numbers"].items():
        assert entry["tier"] in {f"T{n}" for n in range(1, 10)} and entry["tweak"] in {f"L{n}" for n in range(5)} and entry["basis"], name


def test_a_bad_entry_is_refused(tmp_path, monkeypatch):
    document = copy.deepcopy(acid_medium.DOCUMENT)
    document["numbers"]["so2_min_share"]["tier"] = "T10"
    path = tmp_path / "acid_medium.json"
    path.write_text(json.dumps(document))
    monkeypatch.setattr(acid_medium, "DATA_FILE", str(path))
    with pytest.raises(acid_medium.AcidMediumError):
        acid_medium._load()
