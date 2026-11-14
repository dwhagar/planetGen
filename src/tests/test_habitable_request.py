# tests/test_habitable_request.py

"""GEN.911, GEN.916: a request for a habitable world asks for conditions (temperate, air, a surface ocean), never for a class letter, and
aims its orbit by the temperature the draw solved."""

import types

import pytest

from planetgen import tuning
from planetgen.generation import habitable_request as hr
from planetgen.generation import system as systemData
from planetgen.generation.config import SystemConfig
from planetgen.physics import habitability
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _body(temperature=290.0, pressure=101325.0, hydrosphere="surface ocean", planet_class="C", moons=(), body_type="t", ocean_fraction=0.71):
    return types.SimpleNamespace(surface_temperature=temperature, atmospheric_pressure=pressure, hydrosphere=hydrosphere,
                                 planet_class=planet_class, moons=list(moons), body_type=body_type, ocean_fraction=ocean_fraction)


def test_a_temperate_world_with_air_and_a_sea_meets_the_request_whatever_its_class_letter():
    for letter in ("A", "C", "M", "Z"):
        assert hr.meets(_body(planet_class=letter))


@pytest.mark.parametrize("change", [
    {"temperature": 250.0}, {"temperature": 340.0}, {"temperature": None},
    {"pressure": 5000.0}, {"pressure": 40000.0}, {"pressure": 300000.0}, {"pressure": None},
    {"ocean_fraction": 0.1}, {"ocean_fraction": 0.97}, {"ocean_fraction": 1.0}, {"ocean_fraction": 0.0}, {"ocean_fraction": None}, {"hydrosphere": "ice-covered ocean"}, {"hydrosphere": "hycean"},
    {"hydrosphere": None}, {"body_type": "g"},
])
def test_a_world_missing_any_condition_does_not_meet_it_even_in_a_class_that_can_carry_life(change):
    assert not hr.meets(_body(planet_class="M", **change))


def test_the_edges_of_the_temperature_pressure_and_water_cover_ranges_count():
    low, high = tuning.REQUEST_TEMPERATURE_K
    assert hr.meets(_body(temperature=low)) and hr.meets(_body(temperature=high))
    low_pa, high_pa = tuning.REQUEST_PRESSURE_PA
    assert hr.meets(_body(pressure=low_pa)) and hr.meets(_body(pressure=high_pa))
    assert not hr.meets(_body(pressure=low_pa - 1.0)) and not hr.meets(_body(pressure=high_pa + 1.0))
    low_cover, high_cover = tuning.REQUEST_WATER_COVER
    assert hr.meets(_body(ocean_fraction=low_cover)) and hr.meets(_body(ocean_fraction=high_cover))
    assert not hr.meets(_body(ocean_fraction=low_cover - 0.001)) and not hr.meets(_body(ocean_fraction=high_cover + 0.001))
    assert tuning.REQUEST_PRESSURE_PA == (50000.0, 250000.0) and tuning.REQUEST_WATER_COVER == (0.20, 0.95)


def test_a_moon_counts_and_a_belt_is_skipped():
    moon = _body()
    planet = _body(temperature=200.0, moons=[moon])
    belt = types.SimpleNamespace(body_type="a", moons=[])
    assert not hr.meets(planet)
    assert hr.count([planet, belt, _body()]) == 2


def test_the_class_letter_does_not_decide_the_count():
    assert hr.count([_body(planet_class="M")]) == 1
    assert hr.count([_body(temperature=200.0, planet_class="M")]) == 0


@pytest.fixture(scope="module")
def system():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = False
    draw.set_run_seed(5)
    return StarSystem(system_config=cfg)


class _Draw:
    """A planet drawn at the distance it was given that meets the request with chance `share` (the mass the slot gave stays)."""
    share = 0.1
    seen = []

    def __init__(self, config, star, zone, distance, **kwargs):
        self.body_type = "t"
        self.moons = []
        self.surface_temperature = 290.0
        self.atmospheric_pressure = 101325.0
        self.ocean_fraction = 0.71
        self.hydrosphere = "surface ocean" if draw.random() < type(self).share else "dry"
        type(self).seen.append((distance, kwargs.get("mass")))

    def add_moons(self, count):
        self.moons_added = count


def _drawn_in_ecosphere(system, monkeypatch, share, runs):
    monkeypatch.setattr(systemData, "Planet", type("Draw", (_Draw,), {"share": share, "seen": []}))
    planets = []
    for _ in range(runs):
        system._request_draws_left = 100                # the budget of one system, per planet here
        planets.append(system._planet_in_ecosphere(system.star, system.star.habitable_zone, 1.0, mass=5.0e24))
    return planets, systemData.Planet.seen


def test_a_request_raises_the_odds_of_a_habitable_planet_in_the_ecosphere_and_moves_nothing(system, monkeypatch):
    draw.set_run_seed(11)
    planets, seen = _drawn_in_ecosphere(system, monkeypatch, 0.1, 400)
    share = sum(1 for p in planets if hr.meets(p)) / len(planets)
    odds = share / (1.0 - share)
    assert 0.1 / 0.9 * tuning.REQUEST_ODDS * 0.6 < odds < 0.1 / 0.9 * tuning.REQUEST_ODDS * 1.5     # the odds of the chance, times the multiplier
    assert share > 0.35 and share < 1.0                  # raised, not forced: a run is left without one
    assert {d for d, _ in seen} == {1.0} and {m for _, m in seen} == {5.0e24}


def test_a_request_that_no_draw_can_meet_keeps_the_last_planet_and_places_nothing(system, monkeypatch):
    draw.set_run_seed(12)
    planets, seen = _drawn_in_ecosphere(system, monkeypatch, 0.0, 3)
    assert not any(hr.meets(p) for p in planets)
    assert len(seen) >= 3 and {d for d, _ in seen} == {1.0}
    assert not system._request_unreachable                # a habitable request is no placement; it just comes out unmet


def test_a_planet_with_no_surface_is_kept_as_drawn(system, monkeypatch):
    drawn = []

    class Giant:
        body_type = "g"
        moons = []

        def __init__(self, config, star, zone, distance, **kwargs):
            drawn.append(distance)

        def add_moons(self, count):
            pass

    monkeypatch.setattr(systemData, "Planet", Giant)
    system._planet_in_ecosphere(system.star, system.star.habitable_zone, 1.0)
    assert drawn == [1.0]


def test_a_request_for_a_letter_alone_does_not_move_the_orbit(system, monkeypatch):
    seen = []

    class Recorder:
        planet_class = "C"
        moons = []
        body_type = "t"
        surface_temperature = 100.0

        def __init__(self, config, star, zone, distance, **kwargs):
            seen.append(distance)

        def add_moons(self, count):
            pass

    monkeypatch.setattr(systemData, "Planet", Recorder)
    system.build_planet(system.star, system.star.habitable_zone, 1.0, planet_class="D")
    assert set(seen) == {1.0}


def test_a_requested_habitable_world_that_does_not_come_out_is_said_to_be_unmet_and_not_made():
    answers = set()
    for seed in range(3):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "K5V"
        cfg.BINARY_SYSTEM = False
        cfg.HABITABLE_WORLD = True
        draw.set_run_seed(seed)
        built = StarSystem(system_config=cfg)
        assert "a habitable world" not in built.unmet_requirements           # an odds-only request is never a failed requirement
        met = hr.count(built.planets) > 0
        answers.add(met)
        assert built.habitable_request == ("aim" if met else "unmet")
    assert answers                                                          # some draws may come out unmet; none is forced


def test_the_answer_is_aim_exactly_when_a_body_meets_the_conditions_of_the_request():
    """GEN.1248: the answer came from the PHI-6 index while the request steers by conditions, so a system with a world said
    `unmet` and one without said `fallback`. It is now `aim` when `count` finds a body and `unmet` otherwise; a system with no
    request has no answer, and a moon of a giant counts as a body."""
    for seed in range(12):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V" if seed % 2 else "K5V"
        cfg.BINARY_SYSTEM = False
        cfg.HABITABLE_WORLD = True
        draw.set_run_seed(seed)
        built = StarSystem(system_config=cfg)
        found = hr.count(built.all_planets) > 0
        assert built.habitable_request == ("aim" if found else "unmet"), (seed, found)
    plain = SystemConfig()
    plain.STAR_TYPE = "K5V"
    plain.BINARY_SYSTEM = False
    draw.set_run_seed(0)
    assert StarSystem(system_config=plain).habitable_request is None
    moon = types.SimpleNamespace(body_type="t", surface_temperature=288.0, atmospheric_pressure=1.0e5, hydrosphere="surface ocean", ocean_fraction=0.7, moons=[])
    giant = types.SimpleNamespace(body_type="g", moons=[moon])
    holder = types.SimpleNamespace(system_config=types.SimpleNamespace(HABITABLE_WORLD=True), all_planets=[giant])
    assert StarSystem.habitable_request_status(holder) == "aim"
    holder.all_planets = [types.SimpleNamespace(body_type="g", moons=[])]
    assert StarSystem.habitable_request_status(holder) == "unmet"


def test_a_barred_habitable_world_leaves_no_habitable_body_anywhere_in_the_system():
    """GEN.1275, GEN.1387: `HABITABLE_WORLD` False barred only an ecosphere planet that met the request, so 0.35 of such systems still held a habitable
    body (a moon, an outer planet, a body whose air changed in the life pass). It now bars every planet and moon rated habitable, remade from the same
    conditions when it comes out so and dropped when none of the attempts leaves it. The rating comes from the conditions, never from the letter."""
    for seed in range(8):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V" if seed % 2 else "K5V"
        cfg.BINARY_SYSTEM = False
        cfg.HABITABLE_WORLD = False
        draw.set_run_seed(seed)
        built = StarSystem(system_config=cfg)
        assert built.hab_count == 0, seed
        for body in built.all_planets:
            assert not habitability.is_habitable(getattr(body, "phi4", None)), (seed, body)
            for moon in getattr(body, "moons", []):
                assert not habitability.is_habitable(getattr(moon, "phi4", None)), (seed, moon.planet_class)


def test_a_system_that_bars_habitable_worlds_keeps_out_a_body_that_meets_the_request_or_is_rated_habitable():
    """GEN.1275 on the bar of GEN.1250, GEN.1387: one test for both (`StarSystem._barred`); the class letter decides neither, and life on a world a human
    cannot live on is not barred."""
    meets = dict(body_type="t", surface_temperature=288.0, atmospheric_pressure=1.0e5, hydrosphere="surface ocean", ocean_fraction=0.7)
    cold = dict(body_type="t", surface_temperature=120.0, atmospheric_pressure=1.0e5, hydrosphere="ice")
    assert StarSystem._barred(types.SimpleNamespace(planet_class="D", **meets))
    assert StarSystem._barred(types.SimpleNamespace(planet_class="D", phi4=0.9, **cold))
    assert not StarSystem._barred(types.SimpleNamespace(planet_class="E", phi4=0.2, life_class="3", **cold))
    assert not StarSystem._barred(types.SimpleNamespace(planet_class="D", **cold))
