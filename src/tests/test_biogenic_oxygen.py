# tests/test_biogenic_oxygen.py

"""GEN.855: `physics.biogenic_oxygen`, the oxygen oxygenic life puts in the air, and its wiring in `generation.life`."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import evolution, life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, biogenic_oxygen as bo, gases, life_air, life_classes

from tests.fuzz_support import deterministic_entropy


def _worlds(seeds):
    for seed in seeds:
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        config.PLANETS = True
        config.HABITABLE_WORLD = True
        with deterministic_entropy(seed):
            system = StarSystem(system_config=config)
        for planet in system.planets:
            if planet.body_type != "g" and planet.body_type != "a":
                yield planet


OXYGENIC_WORLDS_WANTED = 1
"""int: How many worlds whose oxygenic life built oxygen into the air the sample must hold for the tests that read them. TEST.664: it was 3 worlds of oxygenic life that could, and a sample of three took 600 s on
4 vCPUs (a G2V system with a request for a habitable world costs 8 to 18 s and about one in 25 holds such a world); one that holds the oxygen it built shows the rise and the air it left."""

FIRST_SEEDS = 4
"""int: The seeds always drawn; more follow until `OXYGENIC_WORLDS_WANTED` worlds hold oxygen their life built."""

LAST_SEED = 80
"""int: The seed after which the sample stops growing (the tests then fail on the missing oxygenic worlds)."""

DECISIONS = {}
"""dict: `{id(world): (rise_gyr, took_oxygen)}` as the generator decided it: the age the oxygen of the world's classes rises (`life_air.oxygen_state`, `None` when it never does) and whether the air took
oxygen at that call (`generation.life._build_oxygen`). The last call for a world is the one its stored air follows."""


def _makes_oxygen(world):
    """Whether the world has oxygenic life whose oxygen the generator put in the air (`DECISIONS`): a request for a habitable world only raises the odds (GEN.1130), so a sample of them holds worlds that cannot, and the tests
    read the ones that can."""
    return bool(world.life_class_id and evolution.releases_oxygen(world.life_class_id)) and DECISIONS.get(id(world), (None, False))[1]


@pytest.fixture(scope="module")
def worlds():
    patch = pytest.MonkeyPatch()
    real = life._build_oxygen

    def recording(planet, living):
        rise, _target = life_air.oxygen_state(planet, living)
        took = real(planet, living)
        DECISIONS[id(planet)] = (rise, took)
        return took

    patch.setattr(life, "_build_oxygen", recording)
    try:
        held = list(_worlds(range(0, FIRST_SEEDS)))
        seed = FIRST_SEEDS
        while seed < LAST_SEED and sum(1 for w in held if _makes_oxygen(w)) < OXYGENIC_WORLDS_WANTED:
            held += list(_worlds([seed]))
            seed += 1
    finally:
        patch.undo()
    return held


def test_earth_numbers_balance_at_the_oxygen_of_earths_air():
    assert bo.WEATHERING_PER_YR * bo.EARTH_O2_MOL == pytest.approx(bo.EARTH_BURIAL_MOL_YR - bo.EARTH_REDUCTANT_MOL_YR)
    # GEN.1463: the store is the surplus summed over the 0.714 Gyr after the biomarkers (2.7 Ga) with the reductant at 2.2 times its present value falling to 1.9 times (the model's heat
    # history of an Earth-mass body, whose heat now is 0.76 of the 44 TW Earth gives off), so the delay at today's reductant is 0.3 to 1 Gyr of the published bracket and shorter than it was when R was held at today's value
    now_delay = bo.EARTH_STORE_MOL / (bo.EARTH_BURIAL_MOL_YR - bo.EARTH_REDUCTANT_MOL_YR) / 1e9
    then_delay = bo.EARTH_STORE_MOL / (bo.EARTH_BURIAL_MOL_YR - 2.2 * 0.76 * bo.EARTH_REDUCTANT_MOL_YR) / 1e9
    assert 0.3 < now_delay < 1.0
    assert then_delay == pytest.approx(0.714, abs=0.12)


def test_a_more_reduced_mantle_puts_out_more_reductants_inside_the_range():
    assert bo.redox_factor(None) == 1.0
    assert bo.redox_factor(bo.EARTH_MANTLE_DELTA_IW) == pytest.approx(1.0)
    assert bo.redox_factor(bo.EARTH_MANTLE_DELTA_IW - 1.0) == pytest.approx(10 ** bo.REDOX_DECADE_PER_UNIT)
    low, high = bo.REDOX_FACTOR_RANGE
    assert bo.redox_factor(-50.0) == high and bo.redox_factor(50.0) == low


def test_the_sea_makes_nothing_outside_its_liquid_range_and_most_near_the_reference():
    cold, _, _, hot = bo.OCEAN_SPAN_K
    assert bo.temperature_factor(None) == 0.0
    assert bo.temperature_factor(cold) == 0.0 and bo.temperature_factor(hot) == 0.0
    assert bo.temperature_factor(bo.REFERENCE_TEMPERATURE_K) == pytest.approx(1.0)
    assert 0.0 < bo.temperature_factor(295.0)


def test_burning_takes_the_oxygen_the_reduced_gases_want_and_leaves_their_oxides():
    moles = {"N2": 10.0, "CH4": 1.0, "H2": 2.0, "CO": 1.0, "O2": 0.5}
    burnt, taken = bo.burn(moles)
    assert taken == pytest.approx(1.0 * 2.0 + 2.0 * 0.5 + 1.0 * 0.5)
    assert set(burnt) == {"N2", "CO2", "O2"} and burnt["CO2"] == pytest.approx(2.0) and burnt["O2"] == 0.5


class _Body:
    """A bare body for the fire limit: Earth's radius and gravity unless given."""

    def __init__(self, gravity=1.0, radius=6371.0):
        self.gravity, self.radius = gravity, radius


def _fire_limit_state(body, moles):
    limit = bo._fire_limit_mol(body, moles)
    air = bo.partials_pa(body, {**moles, "O2": limit})
    return limit, air["O2"], air["O2"] / sum(air.values())


def test_the_fire_limit_is_a_pressure_of_oxygen_on_a_thick_air():
    from planetgen.physics import mixture_oxygen
    body = _Body()
    column = 1.0e6 / (9.80665 * 0.028) * bo.area_m2(body)       # about 10 bar of nitrogen over the surface, mol
    _, pascal, share = _fire_limit_state(body, {"N2": column})
    assert pascal == pytest.approx(mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA, rel=1e-6)
    assert share < 0.05                                                 # a thick air holds 35 kPa as a small share


def test_a_thin_air_may_hold_a_larger_share_than_35_percent():
    from planetgen.physics import mixture_oxygen
    body = _Body()
    for kilopascal, ceiling_first in ((1.0, True), (14.0, False), (50.0, False)):
        column = kilopascal * 1.0e3 / (9.80665 * 0.028) * bo.area_m2(body)         # nitrogen over the surface, mol
        _, pascal, share = _fire_limit_state(body, {"N2": column})
        if ceiling_first:                                  # 1 kPa of nitrogen: the 95 percent ceiling is reached before the pressure cap
            assert share == pytest.approx(mixture_oxygen.ABIOTIC_O2_MAX, rel=1e-6)
            assert pascal < mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA
        else:
            assert pascal == pytest.approx(mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA, rel=1e-6)
            assert 0.35 < share < mixture_oxygen.ABIOTIC_O2_MAX
    assert _fire_limit_state(body, {"N2": 14.0e3 / (9.80665 * 0.028) * bo.area_m2(body)})[2] > _fire_limit_state(
        body, {"N2": 50.0e3 / (9.80665 * 0.028) * bo.area_m2(body)})[2]            # the thinner air holds the larger share


def _air_after_apply(monkeypatch, partials, cap, target_mol):
    """The air `apply` stores on a bare body holding `partials` (Pa) whose retention regime keeps at most `cap` Pa (`None` for no cap)."""
    stored = {}
    monkeypatch.setattr(bo, "store_air", lambda body, air: stored.update(air))
    body = _Body()
    setattr(body, atmosphere.GAS_FIELD, dict(partials))
    body._retention_cap_pa = cap
    return body, stored, bo.apply(body, target_mol)


def test_the_oxygen_goes_in_under_the_pressure_the_retention_regime_keeps(monkeypatch):
    # GEN.1340: a moderate-regime air (10 kPa cap) of nitrogen took 22 to 35 kPa of oxygen on top, 23 to 45 kPa of air that was 75 to 95 percent oxygen
    target = 40.0e3 / (9.80665 * 0.032) * bo.area_m2(_Body())            # the mol of 40 kPa of oxygen, far more than the air has room for
    _body, free, took_free = _air_after_apply(monkeypatch, {"N2": 10.0e3}, None, target)
    assert took_free and sum(free.values()) > 20.0e3                        # no cap known: the old behaviour, the oxygen piles on
    _body, held, took = _air_after_apply(monkeypatch, {"N2": 10.0e3}, 10.0e3, target)
    assert took and sum(held.values()) == pytest.approx(10.0e3, rel=1e-9)   # the regime's cap, the oxygen and nitrogen shares kept
    assert held["O2"] / held["N2"] == pytest.approx(free["O2"] / free["N2"], rel=1e-9)
    assert held["O2"] < 10.0e3


def test_an_air_under_its_cap_is_not_touched_by_it(monkeypatch):
    target = 5.0e3 / (9.80665 * 0.032) * bo.area_m2(_Body())
    _body, free, _ = _air_after_apply(monkeypatch, {"N2": 80.0e3}, None, target)
    _body, held, _ = _air_after_apply(monkeypatch, {"N2": 80.0e3}, 1.0e6, target)
    assert held == pytest.approx(free)


def test_the_fire_limit_of_an_air_with_nothing_else_is_zero():
    assert bo._fire_limit_mol(_Body(), {"O2": 5.0}) == 0.0


def test_the_oxygen_rises_after_the_delay_and_not_before(monkeypatch):
    monkeypatch.setattr(bo, "profile_of", lambda body, share=1.0, extra_store_mol=0.0: bo.Profile(2.0 * share, 1.0, 4.0e9 + extra_store_mol, 3.0, 2.0e4, 4.0))
    assert bo.rise_age_gyr(None, [(0.5, 1.0)]) == pytest.approx(4.5)
    assert bo.rise_age_gyr(None, []) is None


def test_a_world_without_a_sunlit_sea_makes_no_oxygen(worlds):
    seen = 0
    for world in worlds:
        if getattr(world, "hydrosphere", None) not in bo.LIVING_HYDROSPHERES:
            assert bo.productivity(world, 5.1e14) == 0.0
            seen += 1
    assert seen > 3


def test_only_oxygenic_life_holds_oxygen_it_made_and_only_in_a_sunlit_sea(worlds):
    seen = 0
    for world in worlds:
        if not world.life_class_id:
            continue
        o2_kpa = atmosphere.partial_kpa(world, "o2")
        if o2_kpa > 0.2 and not evolution.releases_oxygen(world.life_class_id):
            # an abiotic oxygen air is allowed, but an oxygenic world is the only one that grows its own
            assert world.life_age_complex_cells_gyr is not None
        if evolution.releases_oxygen(world.life_class_id):
            assert world.hydrosphere in bo.LIVING_HYDROSPHERES, (world.name, world.hydrosphere)
            seen += 1
    assert seen >= 1


def test_the_stored_air_stays_one_air_after_the_oxygen_goes_in(worlds):
    for world in worlds:
        if world.life_class_id and evolution.releases_oxygen(world.life_class_id) and atmosphere.partial_kpa(world, "o2") > 1.0:
            partials = atmosphere.gas_partials_pa(world)
            assert sum(partials.values()) == pytest.approx(world.atmospheric_pressure, rel=2e-3), world.name
            assert not (set(partials) & set(bo.REDUCED_GASES)) or partials["O2"] < 1e-3 * sum(partials.values())
            assert all(gases.gas(f) for f in partials)
            return
    pytest.skip("no oxygen made in these seeds")


def test_an_oxygenic_class_eats_light():
    oxygenic = [row for row in life_classes.CLASSES if row.releases_oxygen]
    assert oxygenic
    for row in oxygenic:
        assert row.light and row.metabolism == "oxygenic", row.id


def test_an_oxygenic_world_past_its_rise_holds_oxygen_and_one_before_it_holds_none(worlds):
    """GEN.873, GEN.1213: the air follows the age the oxygen of all the world's classes rises (`life_air.oxygen_state`): a world whose star is older than the rise holds the oxygen, one whose
    oxygen never rises, or rises after its star's age, holds what its air had before (none made by life). TEST.664: the test compared the star's age less the photosynthesis age with the delay of
    the headline class alone (`biogenic_oxygen.profile_of(world)`), which GEN.1213 replaced: the rise also counts the methane the earlier classes made (it takes oxygen to burn) and the producers' shares,
    so four of twelve sampled worlds disagreed with it in both directions while each one's air agreed with the rise the generator decided."""
    past = before = 0
    for world in worlds:
        if not (world.life_class_id and evolution.releases_oxygen(world.life_class_id)) or id(world) not in DECISIONS:
            continue
        rise, took = DECISIONS[id(world)]
        o2 = atmosphere.partial_kpa(world, "o2")
        if rise is not None and rise <= world.star.age:
            assert took and o2 > 0.5, (world.name, o2, rise, world.star.age)
            past += 1
        else:
            assert not took and o2 < 0.2, (world.name, o2, rise, world.star.age)
            before += 1
    assert past + before >= 1


class _Sink:
    """A carbon sink that holds more CO2 the colder the surface (the weathering rate falls with the temperature): 0.001 bar at 300 K, doubling every 10 K colder."""

    supply_bar = 1.0
    solved_bar = 0.001
    solved_temperature_k = 300.0

    @staticmethod
    def kept_at(temperature_k):
        return 0.001 * 2.0 ** ((300.0 - temperature_k) / 10.0)


def _fake_air(monkeypatch, body, partials, greenhouse_k_per_kpa=None):
    """Patches the air reads and the store of `bo` to a body that keeps its partials in `body.air` and has a temperature that falls by `greenhouse_k_per_kpa` per kPa of CO2 less than 1 kPa."""
    monkeypatch.setattr(atmosphere, "gas_partials_pa", lambda b: dict(b.air))

    def store(b, air):
        b.air = dict(air)
        if greenhouse_k_per_kpa is not None:
            b.surface_temperature = 280.0 + greenhouse_k_per_kpa * air.get("CO2", 0.0) / 1000.0
    monkeypatch.setattr(bo, "store_air", store)


def test_a_world_with_no_carbon_sink_is_left_as_it_is(monkeypatch):
    body = SimpleNamespace(air={"N2": 80000.0, "O2": 20000.0, "CO2": 100.0}, surface_temperature=285.0, _co2_sink=None)
    _fake_air(monkeypatch, body, body.air)
    assert bo.rebalance_co2(body) is False and body.air["CO2"] == 100.0


def test_the_co2_follows_the_sink_to_the_temperature_the_oxygen_rise_left(monkeypatch):
    """The air was solved at 300 K with 100 Pa of CO2; the rise left it at 280 K, where the sink holds four times as much. With no greenhouse from the CO2 the
    fixed point is the sink's answer at 280 K: 100 + (400 - 100) Pa."""
    body = SimpleNamespace(air={"N2": 80000.0, "O2": 20000.0, "CO2": 100.0}, surface_temperature=280.0, _co2_sink=_Sink())
    _fake_air(monkeypatch, body, body.air)
    assert bo.rebalance_co2(body) is True
    assert body.air["CO2"] == pytest.approx(400.0, rel=0.03)
    assert body.air["N2"] == 80000.0 and body.air["O2"] == 20000.0


def test_the_co2_and_the_temperature_settle_against_each_other(monkeypatch):
    """With a greenhouse (20 K per kPa of CO2) the warmer air holds less: the pair settles where the CO2 is what the sink holds at the temperature that CO2 gives."""
    body = SimpleNamespace(air={"N2": 80000.0, "O2": 20000.0, "CO2": 100.0}, surface_temperature=280.0, _co2_sink=_Sink())
    _fake_air(monkeypatch, body, body.air, greenhouse_k_per_kpa=20.0)
    bo.rebalance_co2(body)
    held = _Sink.kept_at(body.surface_temperature) * 1.0e5 - _Sink.solved_bar * 1.0e5 + 100.0
    assert body.air["CO2"] == pytest.approx(held, rel=0.1)
    assert body.air["CO2"] > 100.0


def test_a_warmer_air_drives_the_co2_down_to_nothing_not_below(monkeypatch):
    body = SimpleNamespace(air={"N2": 80000.0, "O2": 20000.0, "CO2": 100.0}, surface_temperature=340.0, _co2_sink=_Sink())
    _fake_air(monkeypatch, body, body.air)
    bo.rebalance_co2(body)
    assert 0.0 <= body.air.get("CO2", 0.0) < 100.0


def test_the_sink_of_a_generated_world_gives_the_pressure_it_solved_at_the_temperature_it_solved_at(worlds):
    sinks = [planet._co2_sink for planet in worlds if getattr(planet, "_co2_sink", None) is not None]
    assert sinks, "no generated world kept a carbon sink"
    for sink in sinks:
        assert sink.kept_at(sink.solved_temperature_k) == pytest.approx(sink.solved_bar, rel=1e-9)


def test_the_carbon_dioxide_takes_only_the_room_the_retention_cap_leaves(monkeypatch):
    """A 10 kPa cap on an air of 6.5 kPa of oxygen and nitrogen: the sink asks for far more CO2, but the oxygen is not squeezed out; the CO2 stops at the cap."""
    body = SimpleNamespace(air={"N2": 1000.0, "O2": 6500.0, "CO2": 100.0}, surface_temperature=250.0, _co2_sink=_Sink(), _retention_cap_pa=10000.0)
    _fake_air(monkeypatch, body, body.air)
    bo.rebalance_co2(body)
    assert body.air["O2"] == 6500.0 and body.air["N2"] == 1000.0
    assert sum(body.air.values()) <= 10000.0 + 1e-9
    assert body.air["CO2"] > 100.0
