# tests/test_orbital_events.py

"""
GEN.197: the orbital event writer. A pair's window is asked of `closest_approach`, the hit is resolved by `collision_outcome`, and `record_event` writes the
`orbital_events` row and its bodies; the ID rules keep the heavier body's ID in a merger, an ejected planet's ID, and never hand out a destroyed body's again.
"""

import math

import pytest

from planetgen.admin import activity_log
from planetgen.db import orbital_events as events
from planetgen.db import store
from planetgen.galaxy import object_uid as uid
from planetgen.physics import closest_approach as ca
from planetgen.physics import collision_outcome as outcomes
from planetgen.physics import constants

EARTH_KG, EARTH_M = constants.EARTH_MASS_TO_KG, 6.371e6


def _uid(serial, sector=(1, 2, 3)):
    return uid.to_bytes(uid.pack(*sector, uid.SERIAL_RUNTIME, serial, 0))


def _party(serial, mass=EARTH_KG, radius=EARTH_M, position=(0.0, 0.0, 0.0), velocity=(0.0, 0.0, 0.0), kind="T", name=""):
    key = events.key_of(_uid(serial))
    return events.Party(outcomes.Body(key, kind, mass, radius, position, velocity), _uid(serial), name or f"Body {serial}", "terrestrial", None)


def _segment(position, velocity, epoch=0.0):
    return ca.Segment(position, velocity, epoch)


def test_a_head_on_pair_is_a_hit_with_the_miss_distance_and_the_separation():
    x = _segment((0.0, 0.0, 0.0), (1e4, 0.0, 0.0))
    y = _segment((2.0e8, 3.0e6, 0.0), (0.0, 0.0, 0.0))
    contact = events.find_contact(x, y, 2 * EARTH_M, 0.0, 1.0e5)
    assert contact is not None and not contact.already_overlapping and not contact.window_capped
    assert contact.impact_parameter_m == pytest.approx(3.0e6) and contact.separation_start_m == pytest.approx(math.hypot(2.0e8, 3.0e6))
    assert 1.8e4 < contact.t_hit_s < 1.9e4 and contact.window_start_s == 0.0 and contact.window_end_s == 1.0e5


def test_a_pair_that_passes_wide_or_whose_window_has_not_begun_is_no_hit():
    x = _segment((0.0, 0.0, 0.0), (1e4, 0.0, 0.0))
    assert events.find_contact(x, _segment((2.0e8, 5.0e7, 0.0), (0.0, 0.0, 0.0)), 2 * EARTH_M, 0.0, 1.0e5) is None
    assert events.find_contact(x, _segment((2.0e8, 0.0, 0.0), (0.0, 0.0, 0.0), epoch=2.0e5), 2 * EARTH_M, 0.0, 1.0e5) is None


def test_a_pair_already_touching_at_the_start_of_its_window_is_a_hit_at_once():
    x = _segment((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    contact = events.find_contact(x, _segment((1.0e7, 0.0, 0.0), (0.0, 0.0, 0.0)), 2 * EARTH_M, 5.0, 9.0)
    assert contact.already_overlapping and contact.t_hit_s == 5.0


def test_a_window_longer_than_a_sector_crossing_is_cut_and_says_so():
    x = _segment((0.0, 0.0, 0.0), (1e4, 0.0, 0.0))
    y = _segment((2.0e8, 0.0, 0.0), (0.0, 0.0, 0.0))
    contact = events.find_contact(x, y, 2 * EARTH_M, 0.0, 1.0e5, omega_s=0.0, edge_m=1.0e8)
    assert contact.window_capped


def _slow_earth_pair():
    a = _party(1, position=(0.0, 0.0, 0.0), velocity=(0.0, 0.0, 0.0))
    b = _party(2, mass=0.5 * EARTH_KG, radius=0.8 * EARTH_M, position=(1.0e7, 0.0, 0.0), velocity=(-2.0e3, 0.0, 0.0))
    contact = events.Contact(1000.0, 1.0e6, 1.0e7, 0.0, 5000.0, False, False)
    return a, b, contact


def test_two_planets_that_merge_make_a_merger_event_that_conserves_mass_and_momentum():
    a, b, contact = _slow_earth_pair()
    outcome = events.resolve_contact(a, b, contact, always_destroy_terrestrial=False)
    assert outcome.code == "merger"
    event, bodies = events.event_row(outcome, [a, b], contact, 6000.0, product_ids=events.product_uids(outcome, [a, b]))
    assert event["kind"] == "merger" and event["outcome_code"] == "merger" and event["seed_string"] == outcome.seed
    assert event["mass_error"] < 1e-12 and event["momentum_error"] < 1e-12 and not event["numeric_warning"]
    assert event["impact_parameter_km"] == pytest.approx(1.0e3) and event["window_s"] == 5000.0
    assert event["merged_mass_kg"] == pytest.approx(1.5 * EARTH_KG)
    participants = [row for row in bodies if row["role"] == "participant"]
    products = [row for row in bodies if row["role"] == "product"]
    assert {row["fate"] for row in participants} == {"survived", "absorbed"} and len(products) == 1
    heavier = max((a, b), key=lambda p: p.body.mass_kg)
    assert products[0]["uid"] == heavier.uid and products[0]["fate"] == "survived"       # a merger keeps the heavier body's ID
    assert all(row["fate"] == "survived" for row in participants if row["uid"] == heavier.uid)


def test_a_disruption_destroys_both_and_creates_a_field_with_a_new_id():
    a = _party(1, position=(0.0, 0.0, 0.0), velocity=(3.0e4, 0.0, 0.0))
    b = _party(2, position=(1.0e7, 0.0, 0.0), velocity=(-3.0e4, 0.0, 0.0))
    contact = events.Contact(10.0, 0.0, 1.0e7, 0.0, 100.0, False, False)
    outcome = events.resolve_contact(a, b, contact, always_destroy_terrestrial=True)
    assert outcome.debris is not None and not outcome.survivors
    with pytest.raises(ValueError):
        events.product_uids(outcome, [a, b])
    new = _uid(99)
    ids = events.product_uids(outcome, [a, b], [new])
    assert ids == {"field": new} and new not in (a.uid, b.uid)
    event, bodies = events.event_row(outcome, [a, b], contact, 200.0, product_ids=ids)
    assert event["kind"] == "collision" and {row["fate"] for row in bodies if row["role"] == "participant"} == {"destroyed"}
    field = [row for row in bodies if row["role"] == "product"]
    assert len(field) == 1 and field[0]["fate"] == "created" and field[0]["uid"] == new and field[0]["mass_kg"] == pytest.approx(2 * EARTH_KG)
    assert event["mass_error"] < 1e-12


def test_an_ejected_planet_keeps_its_id_and_a_party_must_be_named_by_its_id():
    a, b, contact = _slow_earth_pair()
    assert events.ejected_uid(a.uid) == a.uid
    wrong = events.Party(outcomes.Body("not-the-id", "T", EARTH_KG, EARTH_M), a.uid)
    with pytest.raises(ValueError):
        events.resolve_contact(wrong, b, contact)


def test_the_activity_log_gets_one_line_for_a_recorded_event(monkeypatch):
    lines = []
    monkeypatch.setattr(activity_log, "event", lambda category, action, **fields: lines.append((category, action, fields)))

    class Recorder:
        def __init__(self):
            self.sql = []

        def execute(self, sql, params=()):
            self.sql.append((sql, params))
            return type("Cursor", (), {"lastrowid": 7})()

    a, b, contact = _slow_earth_pair()
    outcome = events.resolve_contact(a, b, contact, always_destroy_terrestrial=False)
    conn = Recorder()
    assert events.record_event(conn, outcome, [a, b], contact, 6000.0, run_id=3) == 7
    assert conn.sql[0][0].startswith("INSERT INTO orbital_events") and sum(1 for s, _ in conn.sql if "orbital_event_bodies" in s) == 3
    assert lines == [("GEN", "orbital_event", {"id": 7, "kind": "merger", "code": "merger", "regime": events.regime_of(outcome.numbers),
                                               "bodies": f"{events.key_of(a.uid)},{events.key_of(b.uid)}"})]


@pytest.mark.db
def test_a_recorded_event_reads_back_and_a_destroyed_bodys_id_is_never_handed_out_again(mysql_config):
    a = _party(1, position=(0.0, 0.0, 0.0), velocity=(3.0e4, 0.0, 0.0))
    b = _party(2, position=(1.0e7, 0.0, 0.0), velocity=(-3.0e4, 0.0, 0.0))
    contact = events.Contact(10.0, 0.0, 1.0e7, 0.0, 100.0, False, False)
    outcome = events.resolve_contact(a, b, contact, always_destroy_terrestrial=True)
    conn = store.get_connection(mysql_config)
    try:
        first = store._runtime_uids(conn, (1, 2, 3), 1)[0]
        event_id = events.record_event(conn, outcome, [a, b], contact, 200.0, new_field_uid=first, run_id=None, sectors_dirtied=2)
        conn.commit()
        row = conn.execute("SELECT kind, outcome_code, regime, seed_string, mass_error, sectors_dirtied, impact_parameter_km, flags FROM orbital_events WHERE id = ?",
                           (event_id,)).fetchone()
        assert row["kind"] == "collision" and row["outcome_code"] == "disruption" and row["seed_string"] == outcome.seed and row["sectors_dirtied"] == 2
        assert row["mass_error"] < 1e-12 and row["impact_parameter_km"] == 0.0
        rows = conn.execute("SELECT role, ordinal, uid, fate, mass_kg FROM orbital_event_bodies WHERE event_id = ? ORDER BY role DESC, ordinal", (event_id,)).fetchall()
        assert [(r["role"], r["fate"]) for r in rows] == [("product", "created"), ("participant", "destroyed"), ("participant", "destroyed")]
        assert bytes(rows[0]["uid"]) == first
        later = store._runtime_uids(conn, (1, 2, 3), 3)
        assert first not in later and len(set(later)) == 3                                   # the counter only moves up
        conn.execute("DELETE FROM orbital_events WHERE id = ?", (event_id,))
        conn.commit()
    finally:
        conn.close()
