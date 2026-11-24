"""GEN.1227: the coma class, the dominant volatiles and the carbon-chain flag of an interstellar comet (`physics/comet_coma.py`), read off the volatiles its nucleus was drawn with and stored
beside it (`interstellar_comets.coma_class`, `coma_volatiles`, `carbon_chain_depleted`)."""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import InterstellarComet
from planetgen.physics import comet_coma, small_body_composition
from planetgen.util import draw

pytestmark = pytest.mark.db


class _Makeup:
    """The one thing `coma_of` reads of a makeup: its `detail`."""

    def __init__(self, volatiles, depleted=False):
        self.detail = {"volatiles": volatiles, "carbon_chain_depleted": depleted}


def test_the_more_abundant_carbon_oxide_names_the_class():
    assert comet_coma.coma_of(_Makeup({"CO": 3.0, "CO2": 9.0, "CH4": 0.2})).coma_class == comet_coma.CO2_DOMINATED
    assert comet_coma.coma_of(_Makeup({"CO": 12.0, "CO2": 9.0, "CH4": 0.2})).coma_class == comet_coma.CO_DOMINATED
    assert comet_coma.coma_of(_Makeup({"CO": 5.0, "CO2": 5.0})).coma_class == comet_coma.CO2_DOMINATED      # equal numbers: the first class


def test_the_dominant_volatiles_are_the_three_most_abundant_after_water_most_abundant_first():
    coma = comet_coma.coma_of(_Makeup({"CO": 3.0, "CO2": 9.0, "CH4": 0.2, "CH3OH": 1.5, "NH3": 0.7}))
    assert coma.volatiles == ("CO2", "CO", "CH3OH")
    assert coma.text() == "CO2, CO, CH3OH"
    assert len(comet_coma.coma_of(_Makeup({"CO": 3.0, "CO2": 9.0})).volatiles) == 2        # fewer volatiles than the count: all of them


def test_a_tie_between_two_volatiles_is_broken_by_formula_so_the_text_is_one_text():
    coma = comet_coma.coma_of(_Makeup({"CO": 2.0, "CO2": 2.0, "CH4": 2.0, "NH3": 2.0}))
    assert coma.volatiles == ("CH4", "CO", "CO2")


def test_the_carbon_chain_flag_is_the_drawn_one():
    assert comet_coma.coma_of(_Makeup({"CO": 1.0, "CO2": 2.0}, depleted=True)).carbon_chain_depleted is True
    assert comet_coma.coma_of(_Makeup({"CO": 1.0, "CO2": 2.0}, depleted=False)).carbon_chain_depleted is False


def test_a_makeup_with_no_drawn_volatiles_is_refused():
    for detail in ({}, {"volatiles": {}}, {"volatiles": {"CO": 1.0}}, {"volatiles": {"CO2": 1.0}}):
        makeup = _Makeup({})
        makeup.detail = detail
        with pytest.raises(comet_coma.ComaError):
            comet_coma.coma_of(makeup)


def test_reading_the_coma_draws_nothing():
    """The same seed gives the same comet whether or not the coma is read, and the coma is a function of the makeup alone."""
    with draw.bound(5):
        makeup = small_body_composition.draw_comet(interstellar=True)
        after = draw.random()
    with draw.bound(5):
        small_body_composition.draw_comet(interstellar=True)
        again = draw.random()
    assert after == again
    assert comet_coma.coma_of(makeup) == comet_coma.coma_of(makeup)


def test_a_drawn_interstellar_comet_carries_the_coma_of_its_own_volatiles():
    for seed in range(40):
        with draw.bound(seed):
            comet = InterstellarComet(SystemConfig())
        percent = comet.makeup.detail["volatiles"]
        expected = comet_coma.CO2_DOMINATED if percent["CO2"] >= percent["CO"] else comet_coma.CO_DOMINATED
        assert comet.coma_class == expected
        assert comet.coma_volatiles.split(", ")[0] == max(percent, key=lambda formula: (percent[formula], formula))
        assert 1 <= len(comet.coma_volatiles.split(", ")) <= comet_coma.DOMINANT_COUNT
        assert comet.carbon_chain_depleted == comet.makeup.detail["carbon_chain_depleted"]


def test_both_classes_and_both_flags_come_up_in_a_few_dozen_comets():
    seen_class, seen_flag = set(), set()
    for seed in range(120):
        with draw.bound(seed):
            comet = InterstellarComet(SystemConfig())
        seen_class.add(comet.coma_class)
        seen_flag.add(comet.carbon_chain_depleted)
    assert seen_class == set(comet_coma.COMA_CLASSES)
    assert seen_flag == {True, False}


def test_a_comet_survives_a_round_trip_through_its_dict():
    with draw.bound(3):
        comet = InterstellarComet(SystemConfig())
    back = InterstellarComet.from_dict(comet.to_dict(), SystemConfig())
    assert (back.coma_class, back.coma_volatiles, back.carbon_chain_depleted) == (comet.coma_class, comet.coma_volatiles, comet.carbon_chain_depleted)


def test_a_comet_read_back_from_a_record_with_no_coma_has_none():
    with draw.bound(3):
        comet = InterstellarComet(SystemConfig())
    data = comet.to_dict()
    for key in ("coma_class", "coma_volatiles", "carbon_chain_depleted"):
        data.pop(key, None)
    back = InterstellarComet.from_dict(data, SystemConfig())
    assert back.coma_class is None and back.coma_volatiles is None and back.carbon_chain_depleted is None


def test_the_row_keeps_the_coma_and_the_sheet_shows_it(mysql_config):
    from planetgen.describe import sheet
    with draw.bound(11):
        comet = InterstellarComet(SystemConfig())
    conn = store.get_connection(mysql_config)
    try:
        comet_id = store.insert_interstellar_comet(conn, comet)
        conn.commit()
        row = conn.execute("SELECT coma_class, coma_volatiles, carbon_chain_depleted FROM interstellar_comets WHERE id = ?", (comet_id,)).fetchone()
    finally:
        conn.close()
    assert row["coma_class"] == comet.coma_class and row["coma_volatiles"] == comet.coma_volatiles
    assert row["carbon_chain_depleted"] == int(comet.carbon_chain_depleted)
    columns = [spec.key for spec in sheet.specs_for("interstellar_comet")]
    assert {"coma_class", "coma_volatiles", "carbon_chain_depleted"} <= set(columns)
