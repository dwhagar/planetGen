"""GEN.471: clathrate dissociation curves and eutectic brines (`planetgen/physics/clathrate_brine.py`)."""

import json
import math

import pytest

from planetgen.physics import clathrate_brine as cb
from planetgen.physics.elements import DataError
from planetgen.util import draw


def test_the_curves_are_the_marboeuf_table():
    assert set(cb.CLATHRATES) == {"CO", "CO2", "CH4", "H2S", "N2", "NH3"}
    assert cb.dissociation_pressure_pa("CH4", 100.0) == pytest.approx(math.exp(22.6864 - 2176.2 / 100.0))
    assert cb.dissociation_pressure_pa("N2", 80.0) == pytest.approx(math.exp(22.719 - 1696.05 / 80.0))
    assert cb.HYDRATE_NUMBER == 6
    assert [cb.structure(s) for s in ("CO", "CO2", "CH4", "H2S", "N2", "NH3")] == ["I", "I", "I", "I", "II", "monohydrate"]


def test_the_co2_curve_meets_the_measured_quadruple_point():
    # [R4015]: Q1 is 271.60 K at 1.044 MPa; the curve is the fit that ends there.
    assert cb.dissociation_pressure_pa("CO2", 271.60) == pytest.approx(1.044e6, rel=0.03)
    assert cb.upper_temperature_k("CO2") == 271.60


def test_the_pressure_rises_with_temperature_and_the_inverse_agrees():
    for species in cb.CLATHRATES:
        low, high = 90.0, 150.0
        assert cb.dissociation_pressure_pa(species, high) > cb.dissociation_pressure_pa(species, low)
        pressure = cb.dissociation_pressure_pa(species, 120.0)
        assert cb.dissociation_temperature_k(species, pressure) == pytest.approx(120.0)
        assert cb.hydrate_stable(species, 120.0, 1.01 * pressure) and not cb.hydrate_stable(species, 120.0, 0.99 * pressure)


def test_a_curve_stops_where_liquid_water_takes_over():
    with pytest.raises(cb.OutOfRange):
        cb.dissociation_pressure_pa("CO2", 275.0)
    with pytest.raises(cb.OutOfRange):
        cb.dissociation_pressure_pa("CH4", 280.0)
    with pytest.raises(cb.OutOfRange):
        cb.dissociation_temperature_k("CH4", 1e9)
    with pytest.raises(ValueError):
        cb.dissociation_pressure_pa("CH4", 0.0)
    with pytest.raises(KeyError):
        cb.dissociation_pressure_pa("He", 100.0)


def test_a_brine_stays_liquid_below_the_freezing_point_of_water():
    for salt in cb.BRINES:
        low, high = cb.liquid_range_k(salt)
        assert low < high == cb.ICE_POINT_K and low < 273.15 - 30.0
    assert cb.eutectic_k("NaClO4") == 236.0 and cb.eutectic_k("Mg(ClO4)2") == 206.0 and cb.eutectic_k("Ca(ClO4)2") == 198.0
    assert cb.eutectic_wt_percent("Mg(ClO4)2") == 44.0 and cb.eutectic_wt_percent("Ca(ClO4)2") is None


def test_the_freezing_point_falls_to_the_eutectic_and_stays_there():
    values = [cb.freezing_point_k("Mg(ClO4)2", w) for w in (0.0, 10.0, 22.0, 44.0, 60.0)]
    assert values[0] == pytest.approx(273.15) and values[3] == pytest.approx(206.0) and values[4] == pytest.approx(206.0)
    assert values == sorted(values, reverse=True)
    with pytest.raises(ValueError):
        cb.freezing_point_k("Ca(ClO4)2", 10.0)
    with pytest.raises(ValueError):
        cb.freezing_point_k("NaClO4", 120.0)


def test_the_eutectic_draw_is_normal_about_the_stated_value_and_deterministic():
    with draw.bound(1):
        values = [cb.draw_eutectic_k("NaClO4") for _ in range(400)]
    assert sum(values) / len(values) == pytest.approx(236.0, abs=0.2)
    assert max(values) - min(values) > 2.0
    with draw.bound(2):
        a = [cb.draw_eutectic_k("Mg(ClO4)2") for _ in range(5)]
    with draw.bound(2):
        b = [cb.draw_eutectic_k("Mg(ClO4)2") for _ in range(5)]
    assert a == b


def test_the_data_file_names_a_source_for_every_row_and_a_limit_for_every_partly_verified_one(monkeypatch, tmp_path):
    with open(cb.DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    assert all(row["source"] in data["sources"] for section in ("clathrates", "brines") for row in data[section].values())
    broken = json.loads(json.dumps(data))
    broken["brines"]["NaClO4"]["source"] = "nowhere"
    path = tmp_path / "clathrate_brine.json"
    path.write_text(json.dumps(broken), encoding="utf-8")
    monkeypatch.setattr(cb, "DATA_FILE", str(path))
    with pytest.raises(DataError):
        cb._load()
    broken = json.loads(json.dumps(data))
    del broken["brines"]["Ca(ClO4)2"]["limit"]
    path.write_text(json.dumps(broken), encoding="utf-8")
    with pytest.raises(DataError):
        cb._load()
