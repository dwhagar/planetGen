# tests/test_stable_body_ids.py -- DB.75: a body keeps its ID when the system's content is drawn again.

import functools

from planetgen.db import store
from planetgen.galaxy import object_uid
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _system(seed):
    draw.set_run_seed(seed)
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.PLANETS = True
    config.MAX_PLANETS = True
    config.BINARY_SYSTEM = False
    return StarSystem(system_config=config), config


@functools.lru_cache(maxsize=None)
def _systems_by_planet_count(seeds=range(1, 80), enough=4):
    """Generated systems, keyed by how many orbits (planets and belts) they hold; it stops at `enough` different counts."""
    found = {}
    for seed in seeds:
        system, config = _system(seed)
        count = len(system.planets)
        if count >= 2 and count not in found:
            found[count] = (system, config)
            if len(found) >= enough:
                break
    return found


def _numbers(conn, system_id):
    """`{slot: body number}` of the bodies of a system."""
    out = {}
    for slot, uid in store._slot_uids(conn, system_id).items():
        out[slot] = object_uid.body_of(object_uid.from_bytes(uid))
    return out


def test_a_slot_names_a_body_by_its_kind_and_place_not_by_its_row(mysql_config):
    by_count = _systems_by_planet_count()
    system, config = next(iter(by_count.values()))
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        slots = store._body_slots(conn, system_id)
        kinds = {slot[0] for slot in slots.values()}
        assert "star" in kinds and "planet" in kinds or "belt" in kinds
        assert len(set(slots.values())) == len(slots), "two bodies share a slot"
    finally:
        conn.close()


def test_regenerating_a_system_keeps_the_id_of_every_body_that_keeps_its_place(mysql_config):
    by_count = _systems_by_planet_count()
    counts = sorted(by_count)
    assert len(counts) >= 3, f"the seeds gave systems with only {counts} orbits"
    first, second, third = (by_count[counts[i]] for i in (1, 0, 2))   # a few orbits, fewer, more
    system_id = store.save_system(first[0], first[1], config=mysql_config)
    conn = store.open_write(mysql_config)
    try:
        before = _numbers(conn, system_id)
        ever = dict(before)       # {slot: the number it has had}
        for content in (second, third):
            with conn:
                store.replace_star_system_content(conn, system_id, content[0], content[1])
            after = _numbers(conn, system_id)
            assert len(set(after.values())) == len(after), "two bodies took the same number"
            for slot in set(before) & set(after):
                assert after[slot] == before[slot], f"{slot} changed its ID"
            for slot in set(after) - set(before):
                if slot in ever:      # a slot an earlier draw left empty (DB.75): the object it was
                    assert after[slot] == ever[slot], f"{slot} did not get its number back"
                else:
                    assert after[slot] not in ever.values(), f"{slot} took a number an earlier body had held"
            ever.update(after)
            before = after
    finally:
        conn.close()


def test_a_fresh_system_draws_its_ids_in_the_order_it_always_did(mysql_config):
    system, config = _system(3)
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        numbers = sorted(_numbers(conn, system_id).values())
        assert numbers == list(range(1, len(numbers) + 1))
    finally:
        conn.close()


def test_a_slot_left_empty_by_one_draw_gets_its_number_back_in_a_later_one(mysql_config):
    by_count = _systems_by_planet_count()
    counts = sorted(by_count)
    more, fewer = by_count[counts[-1]], by_count[counts[0]]
    system_id = store.save_system(more[0], more[1], config=mysql_config)
    conn = store.open_write(mysql_config)
    try:
        original = _numbers(conn, system_id)
        for content in (fewer, more):
            with conn:
                store.replace_star_system_content(conn, system_id, content[0], content[1])
            if content is fewer:
                assert len(_numbers(conn, system_id)) <= len(original)
        again = _numbers(conn, system_id)
        assert again == original, "the same content drawn again after a smaller draw changed the numbers"
        remembered = {row["slot_key"]: row["body_number"] for row in conn.execute(
            "SELECT slot_key, body_number FROM body_numbers WHERE star_system_id = ?", (system_id,)).fetchall()}
        # the map holds every slot any draw filled, the smaller system's own slots too
        assert {key: remembered.get(key) for key in map(store._slot_key, original)} == {store._slot_key(slot): number for slot, number in original.items()}
    finally:
        conn.close()


def test_the_numbers_go_with_the_system_and_are_written_by_a_draw_not_a_first_save(mysql_config):
    system, config = _system(3)
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.open_write(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM body_numbers WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == 0
        with conn:
            store.replace_star_system_content(conn, system_id, system, config)
        assert conn.execute("SELECT COUNT(*) AS n FROM body_numbers WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == len(_numbers(conn, system_id))
        with conn:
            conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
        assert conn.execute("SELECT COUNT(*) AS n FROM body_numbers WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_slot_is_written_as_the_path_of_its_host_chain():
    assert store._slot_key(("star", "primary", 0)) == "star:primary:0"
    assert store._slot_key(("planet", "primary", "orbit:3")) == "planet:primary:orbit:3"
    assert store._slot_key(("planet", "primary", 3)) == "planet:primary:3"
    assert store._slot_key(("planet", "primary", "belt:1:7")) == "planet:primary:belt:1:7"
    assert store._slot_key(("belt", "primary", "rank:0")) == "belt:primary:rank:0"
    assert store._slot_key(("moon", ("planet", "system", "orbit:1"), "born:2")) == "moon:planet:system:orbit:1:born:2"


SEARCH = range(1, 120)
"""range: The seeds searched for a system with the bodies a test needs; a search that finds none fails with the count tried."""


@functools.lru_cache(maxsize=None)
def _seed_with_three_planned_planets():
    for seed in SEARCH:
        system, _config = _system(seed)
        slots = sorted(planet.plan_slot for planet in system.planets if planet.body_type != "a" and planet.plan_slot is not None)
        if len(slots) >= 3:
            return seed
    raise AssertionError(f"none of the {len(SEARCH)} systems tried held three planned planets")


@functools.lru_cache(maxsize=None)
def _seed_with_two_born_moons():
    for seed in SEARCH:
        system, _config = _system(seed)
        for planet in system.planets:
            born = planet.build.born_moons if planet.body_type != "a" else None
            kept = [moon for moon in (born or ()) if any(moon is have for have in planet.moons)]
            if planet.body_type != "a" and planet.plan_slot is not None and len(kept) >= 2:
                return seed
    raise AssertionError(f"none of the {len(SEARCH)} systems tried held a planet with two born moons")


def test_a_draw_that_trims_a_planet_leaves_the_id_of_every_planet_it_keeps(mysql_config):
    """DB.75b: the draw keeps one planet fewer, from the middle of the plan, and the planets after it keep their IDs (an index in the list would have moved them up one place)."""
    seed = _seed_with_three_planned_planets()
    system, config = _system(seed)
    system_id = store.save_system(system, config, config=mysql_config)
    again, _config = _system(seed)
    planned = [planet for planet in again.planets if planet.body_type != "a" and planet.plan_slot is not None]
    gone = planned[1]
    again.planets.remove(gone)
    conn = store.open_write(mysql_config)
    try:
        before = _numbers(conn, system_id)
        with conn:
            store.replace_star_system_content(conn, system_id, again, config)
        after = _numbers(conn, system_id)
        dropped = {slot for slot in before if slot[0] == "planet" and slot[2] == f"orbit:{gone.plan_slot}"}
        assert len(dropped) == 1, "the trimmed planet was not found under its slot"
        kept = {slot for slot in before if slot[0] in ("planet", "belt") and slot not in dropped}
        assert kept, "no planet was left to keep"
        for slot in kept:
            assert after.get(slot) == before[slot], f"{slot} did not keep its ID when {gone.plan_slot} was trimmed"
        assert not (dropped & set(after)), "the trimmed planet's slot is still held"
    finally:
        conn.close()


def test_a_draw_that_trims_a_moon_leaves_the_id_of_every_moon_it_keeps(mysql_config):
    """DB.75b: the draw keeps a planet's first moon out, and the planet's other moons keep their IDs (their place in the born list does not move, their index does)."""
    seed = _seed_with_two_born_moons()
    system, config = _system(seed)
    system_id = store.save_system(system, config, config=mysql_config)
    again, _config = _system(seed)
    planet = next(p for p in again.planets if p.body_type != "a" and p.plan_slot is not None
                  and len([m for m in (p.build.born_moons or ()) if any(m is have for have in p.moons)]) >= 2)
    planet.moons.remove(planet.moons[0])
    conn = store.open_write(mysql_config)
    try:
        before = _numbers(conn, system_id)
        with conn:
            store.replace_star_system_content(conn, system_id, again, config)
        after = _numbers(conn, system_id)
        hosts = {slot for slot in before if slot[0] == "planet" and slot[2] == f"orbit:{planet.plan_slot}"}
        assert len(hosts) == 1, "the planet was not found under its slot"
        host = hosts.pop()
        moons = {slot for slot in before if slot[0] == "moon" and slot[1] == host}
        assert len(moons) >= 2, f"the planet's moons are not stored under {host}"
        survivors = {slot for slot in moons if slot in after}
        assert survivors and len(survivors) == len(moons) - 1, "a draw with one moon fewer must keep all the others under their own slots"
        for slot in survivors:
            assert after[slot] == before[slot], f"{slot} did not keep its ID"
    finally:
        conn.close()


def test_the_history_of_a_drawn_system_goes_with_its_bodies(mysql_config):
    """DB.75b: a draw moves `body_history` and `lost_moons` to the kept system, so the slots a second draw reads are the first draw's."""
    first = _system(3)
    second = _system(4)
    system_id = store.save_system(first[0], first[1], config=mysql_config)
    conn = store.open_write(mysql_config)
    try:
        with conn:
            store.replace_star_system_content(conn, system_id, second[0], second[1])
        rows = conn.execute("SELECT body_type, body_id FROM body_history WHERE star_system_id = ?", (system_id,)).fetchall()
        planets = {row["id"] for row in conn.execute("SELECT id FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()}
        assert rows, "the drawn system left no body_history"
        assert {row["body_id"] for row in rows if row["body_type"] == store.effect_ledger.BODY_CODES["planets"]} <= planets
        assert conn.execute("SELECT COUNT(*) AS n FROM body_history WHERE star_system_id NOT IN (SELECT id FROM star_systems)").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_stored_planet_gets_its_slot_in_the_plan_back_so_an_edit_writes_it_again(mysql_config):
    seed = _seed_with_three_planned_planets()
    system, config = _system(seed)
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    made = sorted(planet.plan_slot for planet in system.planets if planet.body_type != "a" and planet.plan_slot is not None)
    got = sorted(planet.plan_slot for planet in loaded.planets if planet.body_type != "a" and planet.plan_slot is not None)
    assert made and got == made
