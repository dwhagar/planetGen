# tests/test_sedimentary_deposits.py

"""GEN.1378: the sedimentary and biogenic rows of the deposit table (`physics/sedimentary_deposits.py`)."""

import pytest

from planetgen.physics import ore_deposits as od, phosphorus_budget as pb, sedimentary_deposits as sd
from planetgen.util import draw


def _state(**over):
    """An Earth-like finished body at 4.5 Gyr: a neutral surface ocean over 71 percent, life from 1.5 Gyr, complex cells at 3.0, animals at 3.9, flora scoring 5."""
    base = dict(age_gy=4.5, hydrosphere="surface ocean", ocean_class="neutral", ocean_fraction=0.71, land_fraction=0.29, temperature_k=288.0, co2_pa=40.0, feo=0.081,
                crust_p2o5_wt=0.14, phosphate_umol_l=2.3, photosynthesis_gy=1.5, complex_cells_gy=3.0, multicellular_gy=3.9, flora=5, salts={"halite": 0.7, "gypsum": 0.3})
    base.update(over)
    return sd.State(**base)


def _rows(state, seed="a"):
    return {kind.id: sd.assess(kind, state, draw.Stream(f"{seed}:{kind.id}")) for kind in sd.TYPES}


def _present_somewhere(state, kind, tries=40):
    return [row for row in (sd.assess(sd.TYPE_BY_ID[kind], state, draw.Stream(f"{i}:{kind}")) for i in range(tries)) if row["status"] == od.PRESENT]


# --- the gates ------------------------------------------------------------------------------------------------------------
def test_an_earth_like_body_admits_all_seven_types():
    for kind in sd.TYPES:
        assert _present_somewhere(_state(), kind.id), kind.id


def test_a_lifeless_body_has_only_the_deposits_that_need_no_life():
    state = _state(photosynthesis_gy=None, complex_cells_gy=None, multicellular_gy=None, flora=0)
    for kind in (sd.BIF, sd.PHOSPHORITE, sd.CHERT, sd.COAL, sd.KEROGEN):
        assert sd.TYPE_BY_ID[kind].setting(state) is None, kind
    assert sd.TYPE_BY_ID[sd.EVAPORITE].setting(state) is not None and sd.TYPE_BY_ID[sd.LIMESTONE].setting(state) is not None


def test_no_iron_formation_before_the_oxygen_rise_or_in_the_wrong_ocean():
    assert sd._bif(_state(photosynthesis_gy=4.6)) is None                          # the first photosynthesis comes after the body's age
    for bad in ("acid sulfate", "soda", "chloride brine", "ice-sealed"):
        assert sd._bif(_state(ocean_class=bad)) is None, bad
    window = sd._bif(_state())
    assert window.start_gy == 1.5 and window.end_gy == pytest.approx(1.5 + sd.BIF_SPAN_GY)
    assert sd._bif(_state(photosynthesis_gy=4.0)).end_gy == 4.5                      # the span is cut at the body's age


def test_an_iron_poor_rock_makes_a_lean_formation_that_is_not_stored():
    assert sd.assess(sd.TYPE_BY_ID[sd.BIF], _state(feo=0.02), draw.Stream("lean"))["status"] == od.ABSENT_BY_STATE
    assert all(0.25 <= row["grade"] <= 0.35 for row in _present_somewhere(_state(), sd.BIF))


def test_phosphorite_needs_life_an_open_ocean_and_phosphate():
    assert sd._phosphorite(_state(phosphate_umol_l=0.5)) is None
    assert sd._phosphorite(_state(hydrosphere="ice-covered ocean")) is None
    assert sd._phosphorite(_state(photosynthesis_gy=None)) is None
    assert sd._phosphorite(_state()).start_gy == pytest.approx(1.5 + sd.PHOSPHORITE_LAG_GY)
    richer = sd._phosphorite(_state(phosphate_umol_l=9.2)).opportunity
    assert richer == pytest.approx(4 * sd._phosphorite(_state()).opportunity)
    poor_crust = sd.assess(sd.TYPE_BY_ID[sd.PHOSPHORITE], _state(crust_p2o5_wt=0.02), draw.Stream("p"))
    assert poor_crust["status"] == od.ABSENT_BY_STATE                                   # 200 x 0.02 percent is 4 percent, under the 18 percent floor of the ore


def test_an_evaporite_takes_its_key_and_grade_from_the_salt_sequence():
    setting = sd._evaporite(_state())
    assert setting.key == "salt" and setting.grade == 0.7
    assert sd._evaporite(_state(temperature_k=250.0)) is None and sd._evaporite(_state(land_fraction=0.0)) is None and sd._evaporite(_state(salts={})) is None
    row = _present_somewhere(_state(), sd.EVAPORITE)[0]
    assert row["key_element"] == "salt" and row["grade"] == 0.7


def test_coal_needs_land_flora_and_multicellular_life():
    assert sd._coal(_state(flora=2)) is None and sd._coal(_state(multicellular_gy=None)) is None and sd._coal(_state(land_fraction=0.0)) is None
    assert sd._coal(_state()).start_gy == 3.9


def test_the_biogenic_types_follow_their_milestones():
    assert sd._chert(_state(complex_cells_gy=None)) is None and sd._chert(_state()).start_gy == 3.0
    assert sd._kerogen(_state()).start_gy == pytest.approx(1.5 + sd.KEROGEN_LAG_GY)
    assert sd._limestone(_state(co2_pa=0.0)) is None
    assert sd._limestone(_state()).opportunity == pytest.approx(sd.LIMESTONE_BIOTIC_BOOST)
    assert sd._limestone(_state(complex_cells_gy=None)).opportunity == pytest.approx(1.0)


def test_a_row_has_the_stored_fields_and_a_window_inside_the_age():
    for kind, row in _rows(_state()).items():
        assert set(row) == set(od.ROW_FIELDS)
        if row["status"] == od.PRESENT:
            low, high = sd.TYPE_BY_ID[kind].grade_range
            assert low <= row["grade"] <= high and row["endowment_kg"] > 0.0
            assert 0.0 <= row["age_window_gyr_start"] < row["age_window_gyr_end"] <= 4.5
        else:
            assert row["endowment_kg"] is None and row["grade"] is None


def test_the_same_stream_gives_the_same_row_and_a_type_does_not_move_another():
    assert _rows(_state(), "x") == _rows(_state(), "x")
    first = sd.assess(sd.TYPE_BY_ID[sd.COAL], _state(), draw.Stream("same"))
    sd.assess(sd.TYPE_BY_ID[sd.BIF], _state(), draw.Stream("other"))
    assert sd.assess(sd.TYPE_BY_ID[sd.COAL], _state(), draw.Stream("same")) == first


def test_the_ids_fit_the_columns_and_do_not_collide_with_the_magmatic_ones():
    ids = [kind.id for kind in sd.TYPES]
    assert len(set(ids)) == len(ids) == len(sd.SEDIMENTARY_IDS) and not set(ids) & set(od.TYPE_BY_ID)
    assert all(len(name) <= 24 for name in ids) and all(len(kind.key_element) <= 8 for kind in sd.TYPES) and len("salt") <= 8


# --- a generated body ------------------------------------------------------------------------------------------------------
def test_generate_replaces_its_rows_and_keeps_the_magmatic_ones(monkeypatch):
    from types import SimpleNamespace
    body = SimpleNamespace(**{od.ATTRIBUTE: [{"deposit_type": od.CHROMITITE, "status": od.PRESENT}, {"deposit_type": sd.BIF, "status": "stale"}], "mass": 1.0, "radius": 1.0, "distance": 1.0,
                              "star": SimpleNamespace(age=4.5)})
    monkeypatch.setattr(sd, "state_of", lambda b: _state())
    sd.generate(body)
    sd.generate(body)
    rows = getattr(body, od.ATTRIBUTE)
    assert [row["deposit_type"] for row in rows] == [od.CHROMITITE] + [kind.id for kind in sd.TYPES]
    assert next(row for row in rows if row["deposit_type"] == sd.BIF)["status"] != "stale"


def test_a_body_the_ore_pass_did_not_assess_is_left_alone(monkeypatch):
    from types import SimpleNamespace
    body = SimpleNamespace(**{od.ATTRIBUTE: None})
    sd.generate(body)
    assert getattr(body, od.ATTRIBUTE) is None
