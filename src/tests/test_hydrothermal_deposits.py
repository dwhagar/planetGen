# tests/test_hydrothermal_deposits.py

"""GEN.1381: the hydrothermal and uranium rows of the deposit table (`physics/hydrothermal_deposits.py`)."""

from types import SimpleNamespace

import pytest

from planetgen.physics import constants, hydrothermal_deposits as hd, interior_heat, ore_deposits as od, sedimentary_deposits as sd
from planetgen.util import draw

BY_ID = {kind.id: kind for kind in hd.TYPES}


def _state(**over):
    """An Earth-like finished body at 4.5 Gyr: mobile lid, oxidised mantle, 71 percent ocean, 35 km crust, heat flow of the Earth's, oxygen of the air, life."""
    base = dict(age_gy=4.5, lid=interior_heat.MOBILE, redox="oxidized", wet=True, ocean=True, ocean_fraction=0.71, land_fraction=0.29, radius_km=constants.EARTH_RADIUS_KM,
                heat_flux_w_m2=0.087, crust_km=35.0, sulfur_index=0.8, felsic_share=0.3, magma_ratio=1.0, o2_pa=21278.0, reductant=True,
                crust_ppm=dict(hd.CRUST_PPM), u_ppm=0.0218)
    base.update(over)
    return hd.State(**base)


def _present_somewhere(state, kind, tries=40):
    return [row for row in (sd.assess(BY_ID[kind], state, draw.Stream(f"{i}:{kind}")) for i in range(tries)) if row["status"] == od.PRESENT]


# --- the gates ------------------------------------------------------------------------------------------------------------
def test_an_earth_like_body_admits_all_six_types():
    for kind in hd.TYPES:
        assert _present_somewhere(_state(), kind.id), kind.id


def test_a_stagnant_or_dry_body_has_no_hydrothermal_row():
    for state in (_state(lid=interior_heat.STAGNANT), _state(wet=False, ocean=False)):
        for kind in (hd.PORPHYRY, hd.VMS, hd.OROGENIC, hd.EPITHERMAL):
            assert BY_ID[kind].setting(state) is None, kind
    assert BY_ID[hd.SEDEX].setting(_state(wet=False)) is None and BY_ID[hd.SEDEX].setting(_state(lid=interior_heat.STAGNANT)) is not None


def test_porphyry_needs_an_oxidised_mantle_and_a_thick_crust():
    assert hd._porphyry(_state(redox="reduced")) is None and hd._porphyry(_state(redox="intermediate")) is None
    assert hd._porphyry(_state(crust_km=hd.PORPHYRY_MIN_CRUST_KM - 1.0)) is None


def test_vms_needs_a_liquid_ocean_and_heat():
    assert hd._vms(_state(ocean=False)) is None
    assert hd._vms(_state(heat_flux_w_m2=0.01)) is None


def test_sedex_needs_sulfur_a_basin_and_old_sediment():
    assert hd._sedex(_state(sulfur_index=0.1)) is None and hd._sedex(_state(land_fraction=0.0)) is None and hd._sedex(_state(age_gy=0.3)) is None


def test_orogenic_gold_needs_continental_crust():
    assert hd._orogenic(_state(felsic_share=0.0)) is None and hd._orogenic(_state(crust_km=10.0)) is None


def test_epithermal_gold_is_kept_only_in_the_youngest_tenth_of_the_age():
    setting = hd._epithermal(_state())
    assert setting.start_gy == pytest.approx(0.9 * 4.5) and setting.end_gy == pytest.approx(4.5)
    assert hd._epithermal(_state(magma_ratio=0.1)) is None


def test_roll_front_uranium_needs_oxygen_land_granite_and_a_reductant():
    for bad in (dict(o2_pa=0.0), dict(land_fraction=0.0), dict(felsic_share=0.0), dict(reductant=False), dict(wet=False)):
        assert hd._roll_front(_state(**bad)) is None, bad
    assert hd._roll_front(_state(o2_pa=2000.0)).opportunity < hd._roll_front(_state()).opportunity


def test_a_copper_poor_star_makes_a_lean_porphyry_that_is_not_stored():
    lean = _state(crust_ppm={**hd.CRUST_PPM, "Cu": 5.0e-6})
    assert sd.assess(BY_ID[hd.PORPHYRY], lean, draw.Stream("lean"))["status"] == od.ABSENT_BY_STATE
    assert all(0.003 <= row["grade"] <= 0.010 for row in _present_somewhere(_state(), hd.PORPHYRY))


def test_the_uranium_grade_is_held_to_the_rock_s_uranium():
    poor = _state(u_ppm=0.001)
    assert sd.assess(BY_ID[hd.ROLL_FRONT], poor, draw.Stream("poor"))["status"] == od.ABSENT_BY_STATE
    rows = _present_somewhere(_state(), hd.ROLL_FRONT)
    assert rows and all(row["key_element"] == "U3O8" for row in rows)


def test_a_larger_hotter_body_has_more_deposits_to_draw_from():
    big = hd._porphyry(_state(radius_km=2 * constants.EARTH_RADIUS_KM)).opportunity
    assert big == pytest.approx(4 * hd._porphyry(_state()).opportunity)
    assert hd._porphyry(_state(heat_flux_w_m2=10.0)).opportunity == pytest.approx(hd.HEAT_RATIO_CAP * hd._porphyry(_state(heat_flux_w_m2=0.087)).opportunity)


def test_the_ids_fit_the_columns_and_do_not_collide_with_other_slices():
    ids = [kind.id for kind in hd.TYPES]
    assert len(set(ids)) == len(ids) == len(hd.HYDROTHERMAL_IDS) and not set(ids) & set(od.TYPE_BY_ID) and not set(ids) & set(sd.SEDIMENTARY_IDS)
    assert all(len(name) <= 24 for name in ids) and all(len(kind.key_element) <= 8 for kind in hd.TYPES)


# --- a generated body ------------------------------------------------------------------------------------------------------
def test_generate_replaces_its_rows_and_keeps_the_others(monkeypatch):
    body = SimpleNamespace(**{od.ATTRIBUTE: [{"deposit_type": od.CHROMITITE, "status": od.PRESENT}, {"deposit_type": sd.BIF, "status": od.PRESENT},
                                              {"deposit_type": hd.VMS, "status": "stale"}], "mass": 1.0, "radius": 1.0, "distance": 1.0, "star": SimpleNamespace(age=4.5)})
    monkeypatch.setattr(hd, "state_of", lambda b: _state())
    hd.generate(body)
    hd.generate(body)
    rows = getattr(body, od.ATTRIBUTE)
    assert [row["deposit_type"] for row in rows] == [od.CHROMITITE, sd.BIF] + [kind.id for kind in hd.TYPES]
    assert next(row for row in rows if row["deposit_type"] == hd.VMS)["status"] != "stale"


def test_a_body_the_ore_pass_did_not_assess_is_left_alone():
    body = SimpleNamespace(**{od.ATTRIBUTE: None})
    hd.generate(body)
    assert getattr(body, od.ATTRIBUTE) is None


def test_the_star_ratio_follows_the_iron_for_a_metal_the_set_does_not_hold():
    star = SimpleNamespace(metallicity_fe_h=None, age=4.5)
    assert od.over_iron(star, "Au") == 1.0 and od.over_iron(star, "Cu") == 1.0
    assert od.over_iron(SimpleNamespace(metallicity_fe_h=0.0, age=4.5), "Au") == 1.0


# --- generated bodies ------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def sun():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.star import Star
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("hydrothermal-deposits:star"):
        star = Star(cfg)
    star.age = 4.5
    return star


def _planets(star, count=8):
    from planetgen.generation.config import SystemConfig
    from tests import class_support
    found = []
    for seed in range(count):
        planet = class_support.draw_planet(star, 1.6, seed, SystemConfig(), mass=1.0 * 5.97e24, zone_override="c", moon_count=0)
        hd.generate(planet)
        found.append(planet)
    return found


def test_generated_bodies_get_six_rows_and_only_wet_mobile_ones_have_the_arc_types(sun):
    arc = (hd.PORPHYRY, hd.VMS, hd.OROGENIC, hd.EPITHERMAL)
    seen_present = 0
    for planet in _planets(sun):
        rows = {row["deposit_type"]: row for row in planet.deposit_types}
        assert set(hd.HYDROTHERMAL_IDS) <= set(rows)
        for kind in arc:
            if rows[kind]["status"] == od.PRESENT:
                seen_present += 1
                assert interior_heat.lid_of(planet) == interior_heat.MOBILE and planet.hydrosphere in od.LIQUID_WATER_STATES, kind
    assert seen_present > 0


def test_a_generated_body_made_again_has_the_same_rows_and_the_run_stream_does_not_move(sun):
    planet = _planets(sun, 1)[0]
    before = [dict(row) for row in planet.deposit_types]
    with draw.bound("hydrothermal-deposits:outer"):
        hd.generate(planet)
        mark = draw.random()
    with draw.bound("hydrothermal-deposits:outer"):
        assert mark == draw.random()
    assert before == planet.deposit_types
