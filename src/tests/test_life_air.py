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
    assert {row["status"] for row in life_air.CLASSES.values()} == set(life_air.STATUSES)


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


def test_a_pending_class_makes_no_methane_and_a_body_without_a_radius_none(world):
    pending = next(n for n, row in life_air.CLASSES.items() if row["status"] == "pending")
    assert life_air.methane_flux_mol_yr(world, _class(int(pending))) == 0.0
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
    assert life_air.apply_gases(world, _class(2), False)
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
    life_air.apply_gases(world, _class(2), False)
    moles = bo.air_moles(world)
    assert moles["CH4"] / sum(moles.values()) <= life_air.number("methane_max_share") * (1 + 1e-6)
    assert moles["CO2"] > 0.0


def test_applying_twice_changes_nothing_more(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    assert life_air.apply_gases(world, _class(2), False)
    once = atmosphere.gas_partials_pa(world)
    assert not life_air.apply_gases(world, _class(2), False)
    assert atmosphere.gas_partials_pa(world) == once


def test_no_methane_before_the_producers_and_none_for_a_lifeless_star_age(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    world.star.age = 0.2
    assert not life_air.apply_gases(world, _class(2), False)
    assert _methane(world) == pytest.approx(0.0, abs=1e-9)


def test_an_oxygenic_class_gives_no_methane_once_its_oxygen_has_risen(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    oxygenic = _class(1, "a")
    assert life_air.row_of(oxygenic)["methane"]["until"] == "oxygen"
    assert not life_air.apply_gases(world, oxygenic, True)
    assert life_air.apply_gases(world, oxygenic, False)


def test_no_methane_into_an_air_that_holds_free_oxygen(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    moles = bo.air_moles(world)
    moles["O2"] = 0.2 * sum(moles.values())
    bo.store_air(world, bo.partials_pa(world, moles))
    assert not life_air.apply_gases(world, _class(2), False)


def test_hydrogen_eaters_draw_hydrogen_down_to_the_cap(world, monkeypatch):
    _force_flux(monkeypatch, 0.0)
    moles = bo.air_moles(world)
    moles["H2"] = 0.3 * sum(moles.values())
    bo.store_air(world, bo.partials_pa(world, moles))
    assert "H2" in life_air.row_of(_class(3, "a"))["draws_down"]
    assert life_air.apply_gases(world, _class(3, "a"), False)
    after = bo.air_moles(world)
    assert after["H2"] / sum(after.values()) == pytest.approx(life_air.number("hydrogen_max_share"), rel=1e-3)


def test_a_giant_or_a_body_without_air_is_left_alone(world, monkeypatch):
    _force_flux(monkeypatch, 1.0e12)
    giant = copy.copy(world)
    giant.body_type = "g"
    assert not life_air.apply_gases(giant, _class(2), False)
    airless = copy.copy(world)
    monkeypatch.setattr(bo, "air_moles", lambda body: {})
    assert not life_air.apply_gases(airless, _class(2), False)


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
    assert not life_air.apply_gases(world, _class(2), False)
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
