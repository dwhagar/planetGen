# tests/test_life_air.py

"""GEN.891, GEN.857, GEN.858, GEN.859: `physics.life_air`, the gases of a life class's ledger, and their wiring in `generation.life`."""

import copy
import json

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, biogenic_oxygen as bo, life_air, life_classes

from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy


def _class(number, variant=None):
    return next(row for row in life_classes.CLASSES if row.number == number and (variant is None or row.variant == variant))


def _worlds(seeds, star="G2V"):
    for seed in seeds:
        config = SystemConfig()
        config.STAR_TYPE = star
        config.PLANETS = True
        config.HABITABLE_WORLD = True
        with deterministic_entropy(seed):
            system = StarSystem(system_config=config)
        for planet in system.planets:
            if planet.body_type == "a":
                continue
            for body in (planet, *planet.moons):
                if body.body_type == "t":
                    yield body


@pytest.fixture(scope="module")
def airy_world():
    """A generated rocky world with an oxygen-free, carbon dioxide and nitrogen air of its own making stored on it."""
    for body in _worlds(range(0, 12)):
        if (body.atmospheric_pressure or 0) > 20000.0 and body.surface_temperature and body.greenhouse_factor is not None:
            body = copy.deepcopy(body)
            bo.store_air(body, {"N2": 70000.0, "CO2": 30000.0})
            return body
    pytest.skip("no world with air in these seeds")


@pytest.fixture
def world(airy_world):
    body = copy.deepcopy(airy_world)
    body.life_age_photosynthesis_gyr = 0.5
    body.star.age = 3.0
    return body


def _methane(body):
    return atmosphere.partial_kpa(body, "ch4")


def test_the_ledger_file_carries_a_tier_a_tweak_and_a_basis_on_every_number_and_class():
    assert life_air.VERSION == 1
    assert set(life_air.CLASSES) == {str(n) for n in range(1, 10)}
    for entry in (*life_air.NUMBERS.values(), *life_air.CLASSES.values()):
        assert entry["tier"] and entry["tweak"] and entry["basis"]
    assert {row["status"] for row in life_air.CLASSES.values()} <= set(life_air.STATUSES)


def test_a_ledger_number_without_a_basis_is_refused(tmp_path, monkeypatch):
    with open(life_air.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    del document["numbers"]["methane_lifetime_yr"]["basis"]
    path = tmp_path / "life_ledger.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    monkeypatch.setattr(life_air, "DATA_FILE", str(path))
    with pytest.raises(life_air.LifeLedgerError):
        life_air._load()


def test_the_published_earth_methane_flux_sits_inside_the_archean_range():
    assert 9.0 <= life_air.number("earth_donor_methane_tmol_yr") <= 67.0


def test_a_pending_class_makes_no_methane_and_a_body_without_a_radius_none(world, monkeypatch):
    monkeypatch.setitem(life_air.CLASSES["2"], "status", "pending")
    assert life_air.methane_flux_mol_yr(world, _class(2), {"donor_supply": 1.0}) == 0.0
    monkeypatch.setitem(life_air.CLASSES["2"], "status", "built")
    bare = copy.copy(world)
    bare.radius = None
    assert life_air.methane_flux_mol_yr(bare, _class(2)) == 0.0


def test_the_methane_of_a_supply_class_scales_with_the_area_and_the_supply(world):
    area = bo.area_m2(world) / bo.EARTH_AREA_M2
    flux = life_air.methane_flux_mol_yr(world, _class(2), {"donor_supply": 0.5})
    assert flux == pytest.approx(0.5 * area * life_air.number("earth_donor_methane_tmol_yr") * 1e12)
    assert life_air.methane_flux_mol_yr(world, _class(2), {"donor_supply": 1.0}) == pytest.approx(2.0 * flux)
    hydrogen = life_air.methane_flux_mol_yr(world, _class(3, "a"), {"h2_supply": 1.0})
    assert hydrogen == pytest.approx(area * life_air.number("earth_h2_tmol_yr") * 1e12 / life_air.number("h2_per_methane"))


def test_the_methane_the_air_holds_is_the_flux_times_its_life_within_the_share_cap():
    moles = {"N2": 1.0e9, "CO2": 5.0e8}
    assert life_air.methane_target_mol(1.0, moles) == pytest.approx(life_air.number("methane_lifetime_yr"))
    assert life_air.methane_target_mol(1e12, moles) == pytest.approx(life_air.number("methane_max_share") * 1.5e9)
    carbon_poor = {"N2": 1.0e12, "CO2": 1.0e6}
    assert life_air.methane_target_mol(1e12, carbon_poor) == pytest.approx((1.0 - life_air.number("carbon_left_share")) * 1.0e6)


def _force_flux(monkeypatch, mol_yr):
    monkeypatch.setattr(life_air, "methane_flux_mol_yr", lambda body, life_class, supplies=None: mol_yr)


def test_producers_put_methane_in_the_oxygen_free_air_and_take_its_carbon_from_the_carbon_dioxide(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    before = bo.air_moles(world)
    assert life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    after = bo.air_moles(world)
    gained = after.get("CH4", 0.0) - before.get("CH4", 0.0)
    assert gained > 0.0
    assert before["CO2"] - after["CO2"] == pytest.approx(gained, rel=1e-6)
    assert after["CO2"] >= before["CO2"] * life_air.number("carbon_left_share") * (1 - 1e-6)
    assert after["N2"] == pytest.approx(before["N2"]) if "N2" in before else True
    partials = atmosphere.gas_partials_pa(world)
    assert sum(partials.values()) == pytest.approx(world.atmospheric_pressure, rel=2e-3)


def test_the_methane_is_held_within_the_share_cap(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e15)
    life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    moles = bo.air_moles(world)
    assert moles["CH4"] / sum(moles.values()) <= life_air.number("methane_max_share") * (1 + 1e-6)
    assert moles["CO2"] > 0.0


def test_applying_twice_changes_nothing_more(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    assert life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    once = atmosphere.gas_partials_pa(world)
    assert not life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    assert atmosphere.gas_partials_pa(world) == once


def test_no_methane_before_the_producers_and_none_for_a_lifeless_star_age(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    world.star.age = 0.2
    assert not life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    assert _methane(world) == pytest.approx(0.0, abs=1e-9)


def test_an_oxygenic_class_gives_no_methane_once_its_oxygen_has_risen(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    oxygenic = _class(1, "a")
    assert life_air.row_of(oxygenic)["methane"]["until"] == "oxygen"
    assert not life_air.apply_gases(world, [life_air.Living(oxygenic)], True)
    assert life_air.apply_gases(world, [life_air.Living(oxygenic)], False)


def test_no_methane_into_an_air_that_holds_free_oxygen(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    moles = bo.air_moles(world)
    moles["O2"] = 0.2 * sum(moles.values())
    bo.store_air(world, bo.partials_pa(world, moles))
    assert not life_air.apply_gases(world, [life_air.Living(_class(2))], False)


def test_hydrogen_eaters_draw_hydrogen_down_to_the_cap(world, monkeypatch):
    _force_flux(monkeypatch, 0.0)
    moles = bo.air_moles(world)
    moles["H2"] = 0.3 * sum(moles.values())
    bo.store_air(world, bo.partials_pa(world, moles))
    assert "H2" in life_air.row_of(_class(3, "a"))["draws_down"]
    assert life_air.apply_gases(world, [life_air.Living(_class(3, "a"))], False)
    after = bo.air_moles(world)
    assert after["H2"] / sum(after.values()) == pytest.approx(life_air.number("hydrogen_max_share"), rel=1e-3)


def _add_gas(world, formula, share):
    moles = bo.air_moles(world)
    moles[formula] = share / (1.0 - share) * sum(moles.values())
    bo.store_air(world, bo.partials_pa(world, moles))


def test_the_sulfur_eaters_draw_hydrogen_sulfide_down_to_the_cap_and_make_no_gas_of_it(world, monkeypatch):
    _force_flux(monkeypatch, 0.0)
    _add_gas(world, "H2S", 0.01)
    before = bo.air_moles(world)
    assert life_air.row_of(_class(4, "a"))["status"] == "built" and not life_air.row_of(_class(4, "a")).get("products")
    assert life_air.apply_gases(world, [life_air.Living(_class(4, "a"))], False)
    after = bo.air_moles(world)
    assert after["H2S"] / sum(after.values()) == pytest.approx(life_air.number("h2s_max_share"), rel=1e-3)
    assert after.get("SO2", 0.0) == pytest.approx(before.get("SO2", 0.0))


def test_the_acid_eaters_turn_the_hydrogen_sulfide_they_take_into_the_same_moles_of_sulfur_dioxide(world, monkeypatch):
    _force_flux(monkeypatch, 0.0)
    _add_gas(world, "H2S", 0.01)
    before = bo.air_moles(world)
    for acid_class in (_class(5, "a"), _class(9)):
        body = copy.deepcopy(world)
        assert life_air.apply_gases(body, [life_air.Living(acid_class)], False)
        after = bo.air_moles(body)
        taken = before["H2S"] - after.get("H2S", 0.0)
        assert taken > 0.0
        assert after.get("SO2", 0.0) - before.get("SO2", 0.0) == pytest.approx(taken * life_air.number("sulfur_dioxide_per_h2s"), rel=1e-6)
        assert 0.0 < after["H2S"] / sum(after.values()) <= life_air.number("h2s_max_share") * (1.0 + 1e-6)


def test_the_cold_hydrogenators_turn_two_thirds_of_a_mole_of_methane_out_of_each_mole_of_hydrogen(world, monkeypatch):
    _force_flux(monkeypatch, 0.0)
    _add_gas(world, "H2", 0.3)
    before = bo.air_moles(world)
    assert life_air.apply_gases(world, [life_air.Living(_class(7))], False)
    after = bo.air_moles(world)
    taken = before["H2"] - after["H2"]
    assert taken > 0.0 and 0.0 < after["H2"] / sum(after.values()) <= life_air.number("hydrogen_max_share") * (1.0 + 1e-6)
    assert after.get("CH4", 0.0) - before.get("CH4", 0.0) == pytest.approx(taken * life_air.number("methane_per_h2_hydrogenated"), rel=1e-6)


def test_a_ledger_product_of_a_gas_the_class_does_not_draw_down_is_refused(tmp_path, monkeypatch):
    with open(life_air.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    document["classes"]["4"]["products"] = {"H2": {"gas": "CH4", "per": "methane_per_h2_hydrogenated"}}
    path = tmp_path / "life_ledger.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    monkeypatch.setattr(life_air, "DATA_FILE", str(path))
    with pytest.raises(life_air.LifeLedgerError):
        life_air._load()


def test_a_giant_or_a_body_without_air_is_left_alone(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    giant = copy.copy(world)
    giant.body_type = "g"
    assert not life_air.apply_gases(giant, [life_air.Living(_class(2))], False)
    airless = copy.copy(world)
    monkeypatch.setattr(bo, "air_moles", lambda body: {})
    assert not life_air.apply_gases(airless, [life_air.Living(_class(2))], False)


def test_the_generated_worlds_keep_one_air_and_methane_beside_life_leaves_no_free_oxygen():
    for body in _worlds(range(0, 10)):
        if not body.atmospheric_pressure:
            continue
        partials = atmosphere.gas_partials_pa(body)
        assert sum(partials.values()) == pytest.approx(body.atmospheric_pressure, rel=2e-3), body.name
        if body.life_class and atmosphere.partial_kpa(body, "ch4") > 0.05:
            assert atmosphere.partial_kpa(body, "o2") < 0.5, body.name


def test_a_methane_under_a_pascal_is_not_worth_redoing_the_air(world, monkeypatch):
    """The carbon dioxide of a carbon-poor air bounds the methane; a gain under `biogenic_oxygen.MINIMUM_PA` leaves the air as it was."""
    moles = bo.air_moles(world)
    moles["CO2"] = 1.0e-9 * sum(moles.values())
    bo.store_air(world, bo.partials_pa(world, moles))
    before = atmosphere.gas_partials_pa(world)
    _force_flux(monkeypatch, 1.0e12)
    assert not life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    assert atmosphere.gas_partials_pa(world) == before


def test_a_world_scored_again_after_life_changed_its_air_keeps_its_flux_swing(world):
    """The swing the insolation gave a body is an input of its score; `habitability_world.generate` keeps it (the rebuild after the gases
    cleared it and left the body unscored by the swing, which the rescoring command then disagreed with)."""
    from planetgen.physics import habitability_world
    world.flux_variation = 0.37
    habitability_world.generate(world, 0.0)
    assert world.flux_variation == 0.37


def test_the_stored_class_is_the_class_of_the_air_life_left():
    """GEN.980: the gases of the ledger change the air after the class was derived; the class is derived again, so a stored class is the one
    its stored conditions give (these seeds held worlds of a hydrogen-eating class whose class stayed that of the air before life)."""
    from planetgen.physics import planet_classes
    seen = 0
    for seed in (1, 8, 9, 14):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.all_planets:
            if planet.body_type == "a":
                continue
            for body in (planet, *planet.moons):
                if body.body_type != "g" and body.planet_class is not None:
                    seen += 1
                    assert planet_classes.classify(body) == body.planet_class, (seed, body.name)
    assert seen > 0


# --- UX.237: the air through the age of the star ------------------------------------------------------------------

STAGES = {"abiogenesis": 0.2, "photosynthesis": 0.5, "complex_cells": 1.5, "multicellularity": 2.5, "technological_civilization": None}


def _steps(history, formula):
    return next(gas["steps"] for gas in history["gases"] if gas["formula"] == formula)


def test_the_air_of_an_oxygenic_world_is_methane_before_the_rise_and_oxygen_after_it(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    monkeypatch.setattr(bo, "profile_of", lambda body, *args: bo.Profile(1.0, 0.0, 0.3e9, 1.0, 1.0, 0.3))
    assert bo.apply(world, 2.0e19)
    history = life_air.air_history(world, [life_air.Living(_class(1), 1.0, STAGES['photosynthesis'])], STAGES)
    assert history["producers_gyr"] == 0.5 and history["rose_gyr"] == pytest.approx(0.8) and history["age_gyr"] == 3.0
    oxygen, methane, carbon = (_steps(history, formula) for formula in ("O2", "CH4", "CO2"))
    assert [age for age, _ in oxygen] == [0.0, 0.5, pytest.approx(0.8)]
    assert oxygen[0][1] == 0.0 and oxygen[1][1] == 0.0
    assert oxygen[2][1] == pytest.approx(atmosphere.partial_kpa(world, "o2"), rel=1e-3)
    assert methane[0][1] == 0.0 and methane[1][1] > 0.0 and methane[2][1] == pytest.approx(0.0, abs=1e-6)
    assert carbon[0][1] > carbon[1][1], "the methane's carbon was carbon dioxide before the producers"
    assert [stage["key"] for stage in history["stages"]] == ["abiogenesis", "photosynthesis", "complex_cells", "multicellularity"]
    assert history["speculative"] is None


def test_a_stage_the_star_is_too_young_for_is_not_listed(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    assert life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    history = life_air.air_history(world, [life_air.Living(_class(2), 1.0, STAGES['photosynthesis'])], {**STAGES, "multicellularity": 9.0})
    assert "multicellularity" not in [stage["key"] for stage in history["stages"]]


def test_a_class_that_keeps_its_methane_has_the_same_air_from_the_producers_on(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    assert life_air.apply_gases(world, [life_air.Living(_class(2))], False)
    history = life_air.air_history(world, [life_air.Living(_class(2), 1.0, STAGES['photosynthesis'])], STAGES)
    assert history["rose_gyr"] is None and {gas["formula"] for gas in history["gases"]} == {"CO2", "CH4"}
    methane = _steps(history, "CH4")
    assert [age for age, _ in methane] == [0.0, 0.5]
    assert methane[0][1] == 0.0 and methane[1][1] == pytest.approx(_methane(world), rel=1e-3)
    ledger = life_air.row_of(_class(2))
    assert history["speculative"] == ({"tier": ledger["tier"], "tweak": ledger["tweak"]} if ledger["tweak"] in ("L3", "L4") else None)


@pytest.mark.parametrize("change", ["giant", "no_producers", "producers_after_now", "no_change"])
def test_a_world_with_nothing_to_chart_has_no_history(world, change):
    ages = dict(STAGES)
    if change == "giant":
        world.body_type = "g"
    elif change == "no_producers":
        ages["photosynthesis"] = None
    elif change == "producers_after_now":
        ages["photosynthesis"] = 5.0
    assert life_air.air_history(world, [life_air.Living(_class(2), 1.0, ages['photosynthesis'])], ages) is None


# ---- GEN.1213: the air of the whole mix of living classes

def _per_class_flux(monkeypatch, fluxes):
    monkeypatch.setattr(life_air, "methane_flux_mol_yr", lambda body, life_class, supplies=None: fluxes.get(life_class.number, 0.0))


def test_the_methane_of_two_classes_is_their_fluxes_times_their_shares(world, monkeypatch):
    _per_class_flux(monkeypatch, {2: 1.0e12, 8: 3.0e11})
    before = bo.air_moles(world)
    assert life_air.apply_gases(world, [life_air.Living(_class(2), 0.6), life_air.Living(_class(8), 0.4)], False)
    after = bo.air_moles(world)
    assert after["CH4"] == pytest.approx(life_air.methane_target_mol(0.6e12 + 0.12e12, before), rel=1e-6)


def test_a_class_whose_producers_are_not_yet_begun_adds_nothing(world, monkeypatch):
    _per_class_flux(monkeypatch, {2: 1.0e12, 8: 1.0e12})
    late = life_air.Living(_class(8), 0.5, producers_gyr=world.star.age + 1.0)
    early = life_air.Living(_class(2), 0.5, producers_gyr=0.5)
    assert life_air.apply_gases(world, [early, late], False)
    alone = copy.deepcopy(world)
    bo.store_air(alone, {"N2": 70000.0, "CO2": 30000.0})
    assert life_air.apply_gases(alone, [life_air.Living(_class(2), 0.5, producers_gyr=0.5)], False)
    assert _methane(world) == pytest.approx(_methane(alone), rel=1e-9)


def test_one_class_with_the_whole_share_is_the_old_single_class_air(world, monkeypatch):
    _per_class_flux(monkeypatch, {2: 1.0e12})
    split = copy.deepcopy(world)
    assert life_air.apply_gases(world, [life_air.Living(_class(2), 1.0)], False)
    assert life_air.apply_gases(split, [life_air.Living(_class(2), 0.25), life_air.Living(_class(2), 0.75)], False)
    assert _methane(world) == pytest.approx(_methane(split), rel=1e-9)


def _mock_profile(monkeypatch, source=2.0, reductant=0.5, store=2.0e9):
    monkeypatch.setattr(bo, "profile_of", lambda body, share=1.0, extra_store_mol=0.0: bo.Profile(source * share, reductant, store + extra_store_mol, 1.0, 1.0, 1.0))


def test_the_oxygen_rise_of_one_class_is_its_start_plus_the_delay(world, monkeypatch):
    _mock_profile(monkeypatch)
    assert bo.rise_age_gyr(world, [(0.5, 1.0)]) == pytest.approx(0.5 + 2.0e9 / 1.5 / 1e9)


def test_the_oxygen_rise_of_two_classes_follows_the_summed_source(world, monkeypatch):
    _mock_profile(monkeypatch)
    # 0.5 of the source from 0 Gyr gives a surplus of 0.5 a year for a Gyr; the second half from 1 Gyr makes it 1.5 a year
    rise = bo.rise_age_gyr(world, [(0.0, 0.5), (1.0, 0.5)])
    assert rise == pytest.approx(1.0 + (2.0e9 - 0.5e9) / 1.5 / 1e9)
    assert rise > bo.rise_age_gyr(world, [(0.0, 1.0)]), "a class that carries only half of the production fills the store later"


def test_a_source_that_the_reductants_eat_never_makes_a_rise(world, monkeypatch):
    _mock_profile(monkeypatch, source=1.0, reductant=2.0)
    assert bo.rise_age_gyr(world, [(0.5, 1.0)]) is None
    assert bo.rise_age_gyr(world, []) is None


def test_the_methane_of_the_classes_before_the_rise_delays_it(world, monkeypatch):
    _mock_profile(monkeypatch, source=2.0e13, reductant=0.5e13, store=2.0e22)
    _per_class_flux(monkeypatch, {1: 5.0e12})
    world.star.age = 9.0
    oxygenic = life_air.Living(_class(1, "a"), 1.0, producers_gyr=0.5)
    rise, target = life_air.oxygen_state(world, [oxygenic])
    assert rise is not None and target == 1.0
    standing = life_air.standing_methane_mol(world, [oxygenic], rise)
    assert standing > 0.0
    assert rise == pytest.approx(0.5 + (2.0e22 + 2.0 * standing) / 1.5e13 / 1e9, rel=1e-9)
    _per_class_flux(monkeypatch, {1: 0.0})
    assert life_air.oxygen_state(world, [oxygenic])[0] < rise


def test_no_oxygen_when_the_star_is_younger_than_the_rise(world, monkeypatch):
    _mock_profile(monkeypatch)
    world.star.age = 1.0
    assert life_air.oxygen_state(world, [life_air.Living(_class(1, "a"), 1.0, producers_gyr=0.5)]) == (None, 0.0)


def test_the_history_of_two_classes_steps_at_each_start_and_at_the_rise(world, monkeypatch):
    _per_class_flux(monkeypatch, {2: 1.0e12, 8: 5.0e12})
    living = [life_air.Living(_class(2), 0.5, producers_gyr=0.5), life_air.Living(_class(8), 0.5, producers_gyr=1.5)]
    assert life_air.apply_gases(world, living, False)
    history = life_air.air_history(world, living, STAGES)
    methane = _steps(history, "CH4")
    assert [age for age, _ in methane] == [0.0, 0.5, 1.5]
    assert methane[0][1] == 0.0 and 0.0 < methane[1][1] < methane[2][1]
    assert methane[2][1] == pytest.approx(_methane(world), rel=1e-3)
    assert [(entry["life_class"], entry["share"]) for entry in history["classes"]] == [(2, 0.5), (8, 0.5)]
    assert history["producers_gyr"] == 0.5
