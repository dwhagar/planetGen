# tests/test_carbide_deposits.py

"""GEN.1423: the moissanite and cohenite rows of the deposit table for a carbon-rich body (`physics/carbide_deposits.py`)."""

from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_draw, carbide_deposits as cd, constants, ore_deposits as od, sedimentary_deposits as sd
from planetgen.util import draw

BY_ID = {kind.id: kind for kind in cd.TYPES}


def _state(**over):
    base = dict(age_gy=4.5, carbon_rock=True, radius_km=constants.EARTH_RADIUS_KM, crust_km=35.0, metal_mass_fraction=0.325)
    base.update(over)
    return cd.State(**base)


def _present(state, kind, tries=40):
    return [row for row in (sd.assess(BY_ID[kind], state, draw.Stream(f"{i}:{kind}")) for i in range(tries)) if row["status"] == od.PRESENT]


def test_a_carbon_rich_body_admits_both_minerals():
    for kind in cd.TYPES:
        assert _present(_state(), kind.id), kind.id


def test_a_silicate_body_has_neither_mineral():
    for kind in cd.TYPES:
        assert kind.setting(_state(carbon_rock=False)) is None
        assert sd.assess(kind, _state(carbon_rock=False), draw.Stream("x"))["status"] == od.ABSENT_BY_STATE


def test_cohenite_needs_metal_and_moissanite_needs_a_crust():
    assert cd._cohenite(_state(metal_mass_fraction=0.0)) is None and cd._moissanite(_state(metal_mass_fraction=0.0)) is not None
    assert cd._moissanite(_state(crust_km=0.0)) is None and cd._cohenite(_state(crust_km=0.0)) is not None


def test_the_moissanite_grade_is_the_sic_share_of_the_rock_the_bulk_draw_makes():
    low, high = bulk_draw.RANGES["graphite_over_rock"]
    assert cd.SIC_GRADE == pytest.approx((1.0 / (1.0 + high), 1.0 / (1.0 + low)))
    rows = _present(_state(), cd.MOISSANITE)
    assert all(cd.SIC_GRADE[0] <= row["grade"] <= cd.SIC_GRADE[1] and row["key_element"] == "SiC" for row in rows)


def test_the_cohenite_grade_stops_at_the_carbon_saturated_share_of_the_metal():
    assert cd.COHENITE_GRADE[1] == pytest.approx(4.3 / 6.69, rel=1e-3)
    assert 12.011 / 179.55 * 100.0 == pytest.approx(cd.COHENITE_C_WT, abs=0.01)         # Fe3C is 6.69 percent carbon
    rows = _present(_state(), cd.COHENITE)
    assert all(cd.COHENITE_GRADE[0] <= row["grade"] <= cd.COHENITE_GRADE[1] and row["key_element"] == "Fe3C" for row in rows)


def test_a_larger_body_and_a_thicker_crust_have_more_moissanite_and_more_metal_more_cohenite():
    assert cd._moissanite(_state(radius_km=2 * constants.EARTH_RADIUS_KM)).opportunity == pytest.approx(4 * cd._moissanite(_state()).opportunity)
    assert cd._moissanite(_state(crust_km=70.0)).opportunity == pytest.approx(2 * cd._moissanite(_state()).opportunity)
    assert cd._cohenite(_state(metal_mass_fraction=0.65)).opportunity == pytest.approx(2 * cd._cohenite(_state()).opportunity)


def test_the_windows_are_inside_the_age():
    for kind in cd.TYPES:
        setting = kind.setting(_state())
        assert 0.0 <= setting.start_gy < setting.end_gy <= 4.5


def test_the_ids_fit_the_columns_and_do_not_collide_with_other_slices():
    from planetgen.physics import hydrothermal_deposits as hd
    ids = [kind.id for kind in cd.TYPES]
    assert len(set(ids)) == len(ids) == len(cd.CARBIDE_IDS)
    assert not set(ids) & (set(od.TYPE_BY_ID) | set(sd.SEDIMENTARY_IDS) | set(hd.HYDROTHERMAL_IDS))
    assert all(len(name) <= 24 for name in ids) and all(len(kind.key_element) <= 8 for kind in cd.TYPES)


def test_generate_replaces_its_rows_and_keeps_the_others(monkeypatch):
    body = SimpleNamespace(**{od.ATTRIBUTE: [{"deposit_type": od.CHROMITITE, "status": od.PRESENT}, {"deposit_type": cd.COHENITE, "status": "stale"}],
                              "mass": 1.0, "radius": 1.0, "distance": 1.0, "star": SimpleNamespace(age=4.5)})
    monkeypatch.setattr(cd, "state_of", lambda b: _state())
    cd.generate(body)
    cd.generate(body)
    rows = getattr(body, od.ATTRIBUTE)
    assert [row["deposit_type"] for row in rows] == [od.CHROMITITE] + [kind.id for kind in cd.TYPES]
    assert next(row for row in rows if row["deposit_type"] == cd.COHENITE)["status"] != "stale"


def test_a_body_the_ore_pass_did_not_assess_is_left_alone():
    body = SimpleNamespace(**{od.ATTRIBUTE: None})
    cd.generate(body)
    assert getattr(body, od.ATTRIBUTE) is None


# --- generated bodies ------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def sun():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.star import Star
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("carbide-deposits:star"):
        star = Star(cfg)
    star.age = 4.5
    return star


def _planets(star, count=8, carbon=False):
    from planetgen.generation.config import SystemConfig
    from tests import class_support
    found = []
    for seed in range(count):
        planet = class_support.draw_planet(star, 1.6, seed, SystemConfig(), mass=1.0 * 5.97e24, zone_override="c", moon_count=0)
        planet.carbon_rock = carbon
        cd.generate(planet)
        found.append(planet)
    return found


def test_carbon_rich_bodies_get_the_minerals_and_oxygen_rich_ones_never_do(sun):
    present = {kind: 0 for kind in cd.CARBIDE_IDS}
    for planet in _planets(sun, carbon=True):
        rows = {row["deposit_type"]: row for row in planet.deposit_types}
        assert set(cd.CARBIDE_IDS) <= set(rows)
        for kind in cd.CARBIDE_IDS:
            present[kind] += rows[kind]["status"] == od.PRESENT
    assert all(count > 0 for count in present.values()), present
    for planet in _planets(sun, carbon=False):
        rows = {row["deposit_type"]: row for row in planet.deposit_types}
        assert all(rows[kind]["status"] == od.ABSENT_BY_STATE for kind in cd.CARBIDE_IDS)


def test_a_generated_body_made_again_has_the_same_rows_and_the_run_stream_does_not_move(sun):
    planet = _planets(sun, 1, carbon=True)[0]
    before = [dict(row) for row in planet.deposit_types]
    with draw.bound("carbide-deposits:outer"):
        cd.generate(planet)
        mark = draw.random()
    with draw.bound("carbide-deposits:outer"):
        assert mark == draw.random()
    assert before == planet.deposit_types
