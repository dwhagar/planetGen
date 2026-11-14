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
    assert store._slot_key(("planet", "primary", 3)) == "planet:primary:3"
    assert store._slot_key(("moon", ("planet", "system", 1), 2)) == "moon:planet:system:1:2"
