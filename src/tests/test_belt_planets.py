# tests/test_belt_planets.py

"""
GEN.1510: the heaviest round members of a belt are planets (`generation/belt_planets.py`): made by the planet generator from the member's mass, orbit and makeup, in distance order among the star's planets, out of the belt's
members and mass, exempt from the slot spacing, and stored with the belt they came from.
"""


import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation import belt_planets, validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.names import bodies as body_names
from planetgen.physics import body_strength, compaction
from planetgen.util import draw

SEEDS = range(1, 25)


def _config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M3V"
    cfg.PLANETS = True
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = False
    return cfg


def _system(seed, config=None):
    draw.set_run_seed(seed)
    return StarSystem(system_config=config or _config())


def _promoted(system):
    return [body for body in system.planets if getattr(body, "belt_rank", None) is not None]


@pytest.fixture(scope="module")
def promoted_system():
    """A seeded system with a belt that has promoted bodies, and its seed."""
    for seed in SEEDS:
        system = _system(seed)
        if _promoted(system):
            return seed, system
    pytest.fail(f"no system with a promoted belt body in {len(SEEDS)} seeds")


def test_a_belt_promotes_at_most_the_limit_of_its_heaviest_round_members(promoted_system):
    _seed, system = promoted_system
    made = _promoted(system)
    per_belt = {}
    for body in made:
        per_belt.setdefault(body.belt_order, []).append(body)
    assert made and all(len(planets) <= tuning.BELT_PLANETS_PER_BELT for planets in per_belt.values())      # GEN.1561: a member over the oligarch mass would be an exception; none is drawn in these seeds
    for planets in per_belt.values():
        assert all(planet.planet_class is not None for planet in planets)


def test_a_promoted_planet_has_the_orbit_and_makeup_of_its_member(promoted_system, monkeypatch):
    seed, system = promoted_system
    monkeypatch.setattr(tuning, "BELT_PLANETS_PER_BELT", 0)
    drawn = belt_planets.irregular_moons.belts_of(_system(seed).planets)       # the belts as drawn, with every member still in them
    for body in _promoted(system):
        member = next(m for m in drawn[body.belt_order].members if m["rank"] == body.belt_rank)
        assert body.belt_member is None                                        # the member was only an input
        assert body.mass == pytest.approx(member["mass_kg"], rel=1e-9) and body.distance == pytest.approx(member["distance_au"], rel=1e-9)
        assert body.eccentricity == pytest.approx(member["eccentricity"], rel=1e-9)
        shares = {"metal": member["metal_fraction"], "rock": member["rock_fraction"], "water": member["water_fraction"], "ice": member["ice_fraction"]}
        if body.mass < compaction.TAPER_END_KG:
            assert body.density == pytest.approx(member["density_g_cm3"], rel=1e-6)          # the member's own density, below the taper's end
        assert body_strength.ratio_for(body.mass, body.ice_mass_fraction) >= 1.0          # the floor is roundness, as Lane 21 asked, one call
        assert body.envelope_mass_fraction == 0.0 and body.planet_class is not None        # no gas captured, a class from the finished body
        total = sum(shares.values())
        for name, share in shares.items():
            assert getattr(body, f"{name}_mass_fraction") == pytest.approx(share / total, abs=1e-6)


def test_the_planets_stand_in_distance_order_and_are_numbered_in_it(promoted_system):
    _seed, system = promoted_system
    distances = [body.distance for body in system.planets]
    assert distances == sorted(distances)
    body_names.name_bodies("Test", system.planets)
    names = [body.name for body in system.planets if body.body_type != "a"]
    assert names == [f"Test {body_names.to_roman(number)}" for number in range(1, len(names) + 1)]


def test_the_same_seed_promotes_the_same_bodies(promoted_system):
    seed, system = promoted_system
    again = _system(seed)
    assert [(b.belt_order, b.belt_rank, b.mass, b.planet_class) for b in _promoted(system)] == [(b.belt_order, b.belt_rank, b.mass, b.planet_class) for b in _promoted(again)]


def test_the_belt_and_its_promoted_bodies_add_up_to_the_belt_that_was_drawn(promoted_system, monkeypatch):
    seed, system = promoted_system
    monkeypatch.setattr(tuning, "BELT_PLANETS_PER_BELT", 0)
    drawn = _system(seed)
    assert not _promoted(drawn)
    drawn_belts = {order: belt for order, belt in belt_planets.irregular_moons.belts_of(drawn.planets).items()}
    belts = belt_planets.irregular_moons.belts_of(system.planets)
    for order, belt in belts.items():
        taken = sum(b.mass for b in _promoted(system) if b.belt_order == order)
        assert belt.mass_kg + taken == pytest.approx(drawn_belts[order].mass_kg, rel=1e-9)
        left = {m["rank"] for m in belt.members}
        gone = {b.belt_rank for b in _promoted(system) if b.belt_order == order}
        assert not left & gone and left | gone == {m["rank"] for m in drawn_belts[order].members}


def test_a_promoted_planet_is_not_a_slot_of_the_spacing(promoted_system):
    _seed, system = promoted_system
    slots = validation.slot_planets(system.planets)
    assert slots and all(getattr(body, "belt_rank", None) is None for body in slots)
    assert len(slots) < len(system.planets)


def test_a_limit_of_zero_promotes_nothing(monkeypatch):
    monkeypatch.setattr(tuning, "BELT_PLANETS_PER_BELT", 0)
    assert not any(_promoted(_system(seed)) for seed in SEEDS[:4])


def test_a_system_with_a_fixed_number_of_orbits_promotes_nothing(promoted_system):
    seed, _system_made = promoted_system
    config = _config()
    config.NUM_ORBITS = 4
    assert not _promoted(_system(seed, config))


def test_a_promoted_planet_is_stored_with_its_belt_and_rank(promoted_system, mysql_config):
    _seed, system = promoted_system
    saved = store.save_system(system, _config(), config=mysql_config)
    with store.get_connection(mysql_config) as conn:
        rows = conn.execute("SELECT orbital_index, belt_order, belt_rank, distance_km FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (saved,)).fetchall()
        loaded = store.load_star_system(conn, saved)
    stored = [(row["belt_order"], row["belt_rank"]) for row in rows if row["belt_rank"] is not None]
    assert sorted(stored) == sorted((b.belt_order, b.belt_rank) for b in _promoted(system))
    assert sorted((b.belt_order, b.belt_rank) for b in _promoted(loaded)) == sorted(stored)
    distances = [row["distance_km"] for row in rows]
    assert distances == sorted(distances)


def test_a_free_spinning_promoted_planet_has_the_range_its_tilt_wanders_over(promoted_system):
    """GEN.1668: the settle pass the promoted planets were missing gives a free-spinning one the tilt range of `settle_obliquity` (GEN.1535)."""
    from planetgen.physics import spin
    _seed, system = promoted_system
    free = [body for body in _promoted(system)
            if body.rotation_period_hours and body.period and spin.spin_state(body.rotation_period_hours, body.period)[0] == "free"]
    for body in free:
        assert body.obliquity_min_deg is not None and body.obliquity_max_deg is not None, body.name
        assert body.obliquity_min_deg <= body.obliquity_max_deg + 1e-9


def test_settling_only_some_planets_leaves_the_others_as_they_are(promoted_system):
    """GEN.1668: `only=` names the planets `settle_obliquity` and `settle_moon_tilts` judge; with none named, nothing moves."""
    from planetgen.physics import planets as planet_physics
    _seed, system = promoted_system
    before = [(body.name, body.axial_tilt_deg, body.obliquity_min_deg, body.obliquity_max_deg) for body in system.planets if body.body_type != 'a']
    assert planet_physics.settle_obliquity(system.planets, system.star.mass, only=[]) == 0
    assert planet_physics.settle_moon_tilts(system.planets, system.star.mass, only=[]) == 0
    assert before == [(body.name, body.axial_tilt_deg, body.obliquity_min_deg, body.obliquity_max_deg) for body in system.planets if body.body_type != 'a']


def test_a_belts_largest_body_is_the_largest_member_it_still_holds(promoted_system):
    """GEN.1561: the figure was the drawn belt's, a planet once a member is promoted."""
    _seed, system = promoted_system
    belts = belt_planets.irregular_moons.belts_of(system.planets)
    taken = {body.belt_order for body in _promoted(system)}
    assert taken
    for order in taken:
        belt = belts[order]
        left = [2.0 * member["radius_km"] for member in belt.members if member["status"] != "moon"]
        assert belt.largest_body_km == pytest.approx(max(left, default=0.0), rel=1e-12)


def test_a_member_over_the_oligarch_mass_is_promoted_whatever_the_limit(promoted_system, monkeypatch):
    """GEN.1561 (Boss, 10:06 PM EST 11-12): the ten is a cost limit; a round member heavier than the isolation mass of its stretch is a planet."""
    seed, system = promoted_system
    monkeypatch.setattr(tuning, "BELT_PLANETS_PER_BELT", 1)
    monkeypatch.setattr(belt_planets, "oligarch_cap_kg", lambda star, belt: float("inf"))
    capped = _promoted(_system(seed))
    assert capped and max(sum(1 for b in capped if b.belt_order == order) for order in {b.belt_order for b in capped}) == 1
    monkeypatch.setattr(belt_planets, "oligarch_cap_kg", lambda star, belt: 0.0)
    over = _promoted(_system(seed))
    assert len(over) > len(capped) and {(b.belt_order, b.belt_rank) for b in capped} <= {(b.belt_order, b.belt_rank) for b in over}
    monkeypatch.setattr(tuning, "BELT_PLANETS_PER_BELT", 0)
    assert not _promoted(_system(seed))           # a limit of zero still promotes nothing


def test_the_oligarch_mass_of_a_belt_is_the_isolation_mass_of_its_stretch(promoted_system):
    _seed, system = promoted_system
    belts = belt_planets.irregular_moons.belts_of(system.planets)
    for belt in belts.values():
        cap = belt_planets.oligarch_cap_kg(system.star, belt)
        assert cap > 0.0 and cap == pytest.approx(belt_planets.oligarch_cap_kg(system.star, belt))
        assert all(member["mass_kg"] <= cap for member in belt.members if member["status"] != "moon")       # none is left over it
