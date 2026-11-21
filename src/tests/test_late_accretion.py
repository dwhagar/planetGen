"""GEN.1493 slices 1 and 2: the air reads the system's late accretion (`generation/late_accretion.py`, `planets.remake_air`), from the belts it has, the stretches it lost and what the giants did to them."""

import types

import pytest

from planetgen.generation import late_accretion
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import impact_history, late_supply, planets as planet_physics, scatter
from planetgen.util import draw
from planetgen.generation.build_record import BuildRecord
from tests.fuzz_support import deterministic_entropy

EARTH_KG = 5.9722e24


def _earth_like(seed):
    """A rocky planet of a G star at 1 AU on its own, with the late accretion the first air drew (`None` for one that kept an envelope: since GEN.1710 a rocky planet the escape law
    leaves an envelope has no ground, its air is solved at the 1 bar level like a giant's, and the late accretion is not drawn for it)."""
    seen = []
    real = impact_history.late_accreted_fraction

    def record(mass_earth, constants):
        fraction = real(mass_earth, constants)
        seen.append(fraction * 100.0)
        return fraction

    impact_history.late_accreted_fraction = record
    try:
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        with deterministic_entropy(seed):
            star = Star(cfg)
            planet = Planet(cfg, star, star.habitable_zone, 1.0, mass=EARTH_KG, moon_count=0)
    finally:
        impact_history.late_accreted_fraction = real
    return planet, (seen[-1] if seen else None)


def _air(planet):
    return (planet.atmospheric_pressure, planet.surface_temperature, dict(planet.atmosphere) if isinstance(planet.atmosphere, dict) else planet.atmosphere)


def test_making_the_air_again_with_the_late_accretion_it_drew_gives_the_planet_it_was():
    """The air solves twice (a first pass, then the one that is kept), each with its own draw of the late accretion; the replay holds both to the one percent, so a
    planet that drew nearly the same twice comes back as it was and one that drew differently comes back within a hair of it."""
    checked = 0
    for seed in range(1, 11):       # GEN.1710: about a half of the Earth-mass planets of a G star at 1 AU keep an envelope and so have no late accretion to replay
        planet, percent = _earth_like(seed)
        if percent is None or planet.body_type != "t" or planet.build.air_replay is None or planet.build.air_replay[3] is None:
            continue       # the last of these: a planet that came out with no air has none to make again (`remake_air` leaves it, GEN.1595)
        before = _air(planet)
        assert planet_physics.remake_air(planet, percent, 0.0)
        after = _air(planet)
        assert after[0] == pytest.approx(before[0], rel=0.01)
        assert after[1] == pytest.approx(before[1], rel=0.01)
        checked += 1
    assert checked >= 3


def test_making_the_air_again_twice_with_one_percent_is_the_same_air():
    planet, percent = _earth_like(4)
    state = planet.build.air_replay
    assert planet_physics.remake_air(planet, percent, 0.0)
    first = _air(planet)
    planet.build.air_replay = state
    assert planet_physics.remake_air(planet, percent, 0.0)
    assert _air(planet) == first


def test_a_different_late_accretion_changes_the_air_and_not_the_run_stream():
    # TEST.714: a seed whose planet keeps a temperate air (one that came out a runaway steam world of 30 MPa holds that air whatever it accreted late, and which seeds do moves with every
    # change to the star draws, GEN.1623), found by search
    from tests.fuzz_support import first_seed_with
    _seed, (planet, percent) = first_seed_with(_earth_like, lambda made: made[1] is not None and made[0].body_type == "t" and made[0].surface_temperature < 400.0, seeds=range(1, 20),
                                               what="an Earth-like planet of a G star at 1 AU with a temperate air")
    before = _air(planet)
    draw.set_run_seed(77)
    state = draw.getstate()
    assert planet_physics.remake_air(planet, percent * 3.0, 0.0)
    assert draw.getstate() == state                      # no draw of the run's stream moved
    assert _air(planet) != before


def test_a_moon_or_a_planet_with_no_replay_is_left_alone():
    planet, percent = _earth_like(2)
    planet.build.air_replay = None
    assert not planet_physics.remake_air(planet, percent)
    moon = types.SimpleNamespace(is_moon=True, build=BuildRecord(air_replay=(None, None, None)))
    assert not planet_physics.remake_air(moon, percent)


def test_the_floor_and_the_ceiling_hold():
    planet = types.SimpleNamespace(mass=EARTH_KG, radius=6371.0, distance=1.0)
    alone = [types.SimpleNamespace(body_type="t", is_moon=False, mass=EARTH_KG, radius=6371.0, distance=1.0)]
    mean = impact_history.late_accreted_percent_mean(1.0)
    assert late_accretion.percent_for(planet, alone, 1.0) == pytest.approx(late_accretion.FLOOR_SHARE * mean)        # no belt: the floor
    huge = types.SimpleNamespace(body_type="a", is_moon=False, primordial_mass_earth=1.0e4, lower_limit=2.0, upper_limit=3.0)
    assert late_accretion.percent_for(planet, [*alone, huge], 1.0) == impact_history.LATE_ACCRETED_PERCENT_MAX      # a vast belt: the ceiling


def test_a_giant_beside_a_belt_leaves_the_planet_inside_less_than_the_same_system_without_it():
    planet = types.SimpleNamespace(mass=EARTH_KG, radius=6371.0, distance=1.0)
    inner = types.SimpleNamespace(body_type="t", is_moon=False, mass=EARTH_KG, radius=6371.0, distance=1.0)
    belt = types.SimpleNamespace(body_type="a", is_moon=False, primordial_mass_earth=10.0, lower_limit=2.1, upper_limit=3.3)
    giant = types.SimpleNamespace(body_type="g", is_moon=False, mass=1.9e27, radius=69911.0, distance=5.2)
    assert late_accretion.percent_for(planet, [inner, belt, giant], 1.0) < 0.3 * late_accretion.percent_for(planet, [inner, belt], 1.0)


def test_a_whole_system_has_its_rocky_planets_set_by_the_belts_and_leaves_the_run_stream_to_the_rest():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = False
    found = 0
    for seed in range(1, 12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=cfg)
        for planet in system.all_planets:
            if planet.body_type == "t" and not planet.is_moon and planet.build.late_accreted_percent is not None:
                found += 1
                assert 0.0 < planet.build.late_accreted_percent <= impact_history.LATE_ACCRETED_PERCENT_MAX
                assert planet.build.air_replay is None
    assert found > 5


def _drawn_percent(settled):
    """The late accretion percent the air was made with: the first draw of the stream where its late accretion began."""
    from planetgen.physics import surface_air

    setting, late_state = settled[0], settled[1]
    stream = draw.Stream(0)
    stream.setstate(late_state)
    design = {}
    with draw.bound(stream):
        surface_air._late_accretion(setting, design)
    return design["late_accreted_percent"]


def test_the_gate_leaves_the_air_when_the_late_accretion_barely_moves_what_the_volcanoes_gave_and_makes_it_when_it_moves_it_much():
    planet, percent = _earth_like(2)
    before = _air(planet)
    state = planet.build.air_replay
    drawn = _drawn_percent(state[3])    # the percent the air was made with: it moves nothing
    assert not planet_physics.remake_air(planet, drawn, 1.0e-9)
    assert _air(planet) == before and planet.build.late_accreted_percent == drawn and planet.build.air_replay is None
    planet.build.air_replay = state
    assert planet_physics.remake_air(planet, percent * 100.0, 1.0e-9)    # a gate every change passes
    assert _air(planet) != before


def test_a_planet_with_no_air_is_not_made_again():
    planet = types.SimpleNamespace(is_moon=False, build=BuildRecord(air_replay=(None, None, None, None)))      # the air found nothing to hold, so it kept no setting
    assert not planet_physics.remake_air(planet, 1.0, 0.0)
    assert planet.build.air_replay is None and planet.build.late_accreted_percent == 1.0


def _stub(body_type="t", deferred=True, **kwargs):
    star = kwargs.pop("star")
    fields = dict(body_type=body_type, is_moon=False, mass=EARTH_KG, radius=6371.0, distance=1.0, star=star, build=BuildRecord(deferred_climate=(0, None) if deferred else None))
    fields.update(kwargs)
    return types.SimpleNamespace(**fields)


def test_the_first_air_reads_the_late_accretion_the_system_gives_and_only_for_a_planet_whose_air_is_still_to_be_made(monkeypatch):
    """PERF.271: `give_first_air` sets the number `apply` would take, on a rocky planet built bare; a planet already made, a giant and a belt are left."""
    monkeypatch.setattr(late_accretion.planet_type, "has_ground", lambda body: body.body_type == "t")
    star = types.SimpleNamespace(mass=1.989e30, starless=False)
    bare = _stub(star=star)
    made = _stub(star=star, deferred=False)
    giant = _stub("g", star=star, mass=1.9e27, radius=69911.0)
    belt = types.SimpleNamespace(body_type="a", is_moon=False, primordial_mass_earth=10.0, lower_limit=2.1, upper_limit=3.3)
    monkeypatch.setattr(late_accretion, "_same_star", lambda belt_body, planet: True)
    planets = [bare, made, giant, belt]
    late_accretion.give_first_air(planets)
    assert bare.build.late_accreted_percent == late_accretion.percent_in(bare, planets) > 0.0
    for left in (made, giant):
        assert left.build.late_accreted_percent is None
    assert not hasattr(belt, "build")      # a belt keeps no build record


def test_a_whole_system_makes_few_airs_again_now_the_first_air_reads_the_late_accretion(monkeypatch):
    """PERF.271: the second pass finds the planets done unless the belts it promotes moved the carbon and the nitrogen past the gate."""
    remade = []
    real = planet_physics.remake_air

    def count(planet, percent, gate=planet_physics.LATE_ACCRETION_GATE):
        made = real(planet, percent, gate)
        remade.append(made)
        return made

    monkeypatch.setattr(planet_physics, "remake_air", count)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "K5V"
    cfg.PLANETS = True
    cfg.MOONS = False
    for seed in range(1, 6):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            StarSystem(system_config=cfg)
    assert remade and sum(remade) < 0.5 * len(remade)


def _origin(primordial_earth, stirred=False, inner_au=30.0, outer_au=50.0):
    return types.SimpleNamespace(primordial_mass_earth=primordial_earth, stirred=stirred, inner_au=inner_au, outer_au=outer_au)


def _fates(**shares):
    fates = dict.fromkeys(scatter.POPULATIONS, 0.0)
    fates.update(shares)
    fates["cold_classical"] = 1.0 - sum(v for k, v in fates.items() if k != "cold_classical")
    return fates


def test_a_stretch_the_giants_cleared_of_everything_is_a_reservoir_the_planets_inside_it_still_read():
    planet = types.SimpleNamespace(mass=EARTH_KG, radius=6371.0, distance=1.0)
    alone = [types.SimpleNamespace(body_type="t", is_moon=False, mass=EARTH_KG, radius=6371.0, distance=1.0)]
    mean = impact_history.late_accreted_percent_mean(1.0)
    gone = late_accretion.early_of_origin(_origin(20.0), 30.0, 50.0, _fates(scattered_disk=0.6, centaurs=0.3, ejected=0.1))
    assert gone.kind == "belt" and gone.mass_kg == pytest.approx(20.0 * EARTH_KG) and gone.inward_share == pytest.approx(late_supply.CROSSING_GAIN * 0.9)
    assert late_accretion.percent_for(planet, alone, 1.0) == pytest.approx(late_accretion.FLOOR_SHARE * mean)        # no belt left: the floor
    assert late_accretion.percent_for(planet, alone, 1.0, [gone]) > late_accretion.percent_for(planet, alone, 1.0)     # the lost stretch feeds it


def test_a_belt_is_read_with_what_the_giants_did_to_it():
    quiet = types.SimpleNamespace(body_type="a", is_moon=False, primordial_mass_earth=5.0, lower_limit=2.1, upper_limit=3.3)
    swept = types.SimpleNamespace(body_type="a", is_moon=False, primordial_mass_earth=5.0, lower_limit=2.1, upper_limit=3.3, scatter_fates=_fates(scattered_disk=0.5, centaurs=0.2), stirred=True)
    assert late_accretion.early_of_belt(quiet).inward_share == pytest.approx(late_supply.INWARD_SHARE)
    assert late_accretion.early_of_belt(swept).inward_share == pytest.approx(late_supply.inward_share_of(swept.scatter_fates, True))
    assert late_accretion.early_of_belt(swept).inward_share > late_accretion.early_of_belt(quiet).inward_share


def test_lost_reservoirs_are_kept_by_star_and_read_only_by_the_planets_of_that_star():
    star, other = types.SimpleNamespace(mass=1.989e30), types.SimpleNamespace(mass=1.989e30)
    planet = types.SimpleNamespace(mass=EARTH_KG, radius=6371.0, distance=1.0, star=star, body_type="t", is_moon=False)
    gone = late_accretion.early_of_origin(_origin(20.0), 30.0, 50.0, _fates(scattered_disk=0.6, centaurs=0.3))
    with_it = late_accretion.percent_in(planet, [planet], {id(star): [gone]})
    assert with_it > late_accretion.percent_in(planet, [planet], {id(other): [gone]}) == late_accretion.percent_in(planet, [planet], None)


def test_build_hands_on_the_stretch_a_giant_cleared_of_everything_and_leaves_out_one_nothing_was_thrown_from(monkeypatch):
    from planetgen.generation import belt_birth
    from planetgen.physics import belt_origin
    from tests.test_systems import make_config
    with draw.bound(5):
        system = StarSystem(system_config=make_config("G2V", PLANETS=False))
    swept = belt_origin.BeltOrigin("outer-disk", (), 30.0, 50.0, 2.7, 20.0, 0.0, 18.0, 2.0, False, (), "")
    dry = belt_origin.BeltOrigin("outer-disk", (), 60.0, 80.0, 2.7, 20.0, 0.0, 18.0, 2.0, False, (), "")
    monkeypatch.setattr(belt_origin, "build", lambda *args, **kwargs: [swept, dry])
    fates = {key: 0.0 for key in scatter.POPULATIONS}
    fates.update(scattered_disk=0.6, centaurs=0.3, ejected=0.1)
    monkeypatch.setattr(belt_birth, "_fates", lambda star, mass, planets, origin, edge, instability=0.0: ((fates, "", 30.0, None) if origin is swept else (None, None, None, None)))
    lost = []
    belts = belt_birth.build(system.system_config, system.star, [], 0.1, 500.0, 500.0, lost=lost)
    assert belts == [] and len(lost) == 1
    assert lost[0].mass_kg == pytest.approx(20.0 * EARTH_KG) and (lost[0].lower_au, lost[0].upper_au) == (30.0, 50.0)
    assert lost[0].inward_share == pytest.approx(late_supply.CROSSING_GAIN * 0.9)
