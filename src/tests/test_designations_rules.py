# tests/test_designations_rules.py

"""GEN.1138 (pure part): the rules of designation history (`planetgen/galaxy/designations.py`; docs/design/architecture/designation-history.md
section 3): a retired ID is never reused, a retired name is not reserved, a change of system ID moves every body, an old ID resolves to the
current object or says it is gone, and a member's designation follows its parent."""

import pytest

from planetgen.galaxy import designations as d
from planetgen.galaxy import object_uid as ou

SYSTEM = ou.pack(1200, 3, 77, ou.SERIAL_GENERATED, 5, 0)
BODIES = [ou.with_body(SYSTEM, n) for n in (1, 2, 3, 40)]
NEW_SYSTEM = ou.pack(1200, 3, 77, ou.SERIAL_RUNTIME, 9, 0)


def retired(designation, kind=d.KIND_UID, epoch=1, gone=None, row=7, object_kind="star_systems", reason="rename"):
    return d.Retired(kind, designation, object_kind, row, None, epoch, reason, gone)


def test_a_new_id_that_is_live_or_was_ever_an_id_is_refused_in_words():
    assert d.check_new_uid(NEW_SYSTEM, {SYSTEM}, {BODIES[0]}) == NEW_SYSTEM
    with pytest.raises(d.DesignationError, match="already an object's ID"):
        d.check_new_uid(SYSTEM, {SYSTEM}, set())
    with pytest.raises(d.DesignationError, match="reserved for ever"):
        d.check_new_uid(BODIES[0], {SYSTEM}, {BODIES[0]})


def test_the_reasons_are_a_closed_list_without_a_move():
    assert "moved" not in d.REASONS and "collision bump" in d.REASONS
    assert d.check_reason("rename") == "rename"
    with pytest.raises(d.DesignationError):
        d.check_reason("moved")


def test_a_system_renumbering_gives_every_body_its_number_under_the_new_id():
    moved = d.cascade(SYSTEM, NEW_SYSTEM, BODIES)
    assert list(moved) == BODIES
    for old, new in moved.items():
        assert ou.body_of(new) == ou.body_of(old) and ou.same_system(new, NEW_SYSTEM) and not ou.same_system(new, SYSTEM)
    assert len(set(moved.values())) == len(BODIES)


@pytest.mark.parametrize("bodies,message", [
    ([SYSTEM], "system's own ID"),
    ([ou.with_body(ou.pack(1, 0, 0, ou.SERIAL_GENERATED, 1, 0), 1)], "not a body of"),
])
def test_a_cascade_refuses_a_body_that_is_not_of_the_system(bodies, message):
    with pytest.raises(d.DesignationError, match=message):
        d.cascade(SYSTEM, NEW_SYSTEM, bodies)


def test_a_cascade_refuses_an_unchanged_system_id_and_a_body_id_given_as_the_system():
    with pytest.raises(d.DesignationError):
        d.cascade(SYSTEM, SYSTEM, BODIES)
    with pytest.raises(d.DesignationError):
        d.cascade(BODIES[0], NEW_SYSTEM, [])


def test_a_live_id_names_its_holder_a_retired_one_redirects_a_deleted_ones_is_gone_and_a_stranger_is_unknown():
    live = lambda uid: ("star_systems", 7) if uid == NEW_SYSTEM else None
    assert d.resolve_uid(NEW_SYSTEM, live, []) == d.Resolution(d.CURRENT, "star_systems", 7)
    redirect = d.resolve_uid(SYSTEM, live, [retired(SYSTEM)])
    assert (redirect.status, redirect.object_row_id) == (d.REDIRECT, 7)
    gone = d.resolve_uid(SYSTEM, live, [retired(SYSTEM, gone=3)])
    assert gone.status == d.GONE and gone.object_row_id == 7
    assert d.resolve_uid(BODIES[0], live, [retired(SYSTEM)]).status == d.UNKNOWN


def test_a_name_is_held_by_its_live_holder_first_then_the_latest_past_holder_that_still_exists():
    rows = [retired("Voranthis", d.KIND_NAME, epoch=1, row=1), retired("Voranthis", d.KIND_NAME, epoch=5, row=2),
            retired("Voranthis", d.KIND_NAME, epoch=9, row=3, gone=10)]
    assert d.resolve_name("Voranthis", [("star_systems", 8)], rows).object_row_id == 8
    past = d.resolve_name("Voranthis", [], rows)
    assert (past.status, past.object_row_id) == (d.REDIRECT, 2) and [r.object_row_id for r in past.formerly] == [3, 2, 1]
    only_gone = d.resolve_name("Voranthis", [], [rows[2]])
    assert only_gone.status == d.GONE and only_gone.object_row_id == 3
    assert d.resolve_name("Nobody", [], rows).status == d.UNKNOWN


def test_a_belt_members_old_designation_follows_the_parent_the_rank_kept():
    belt_old = ou.with_body(SYSTEM, 4)
    belt_new = ou.with_body(NEW_SYSTEM, 4)
    live = lambda uid: ("belts", 11) if uid == belt_new else None
    parent = lambda uid: d.resolve_uid(uid, live, [retired(belt_old, object_kind="belts", row=11)])
    current_of = lambda kind, row: belt_new
    old_text = ou.designation(belt_old, 0x42)
    found, text = d.resolve_member(old_text, parent, current_of)
    assert found.status == d.REDIRECT and text == ou.designation(belt_new, 0x42)
    found, text = d.resolve_member(ou.designation(belt_new, 0x42), parent, current_of)
    assert found.status == d.CURRENT and text == ou.designation(belt_new, 0x42)
    stranger = ou.designation(ou.with_body(ou.pack(5, 0, 0, ou.SERIAL_GENERATED, 1, 0), 2), 7)
    assert d.resolve_member(stranger, parent, current_of) == (d.Resolution(d.UNKNOWN), None)
    with pytest.raises(ValueError):
        d.resolve_member("not a member", parent, current_of)


def test_the_formerly_lines_are_newest_first_with_ids_printed_and_names_as_text():
    rows = [retired(SYSTEM, epoch=2, reason="regeneration"), retired("Old Name", d.KIND_NAME, epoch=4, reason="collision bump")]
    assert d.former_lines(rows) == [(d.KIND_NAME, "Old Name", "collision bump", 4), (d.KIND_UID, ou.format_id(SYSTEM), "regeneration", 2)]
    assert d.former_lines([]) == []


@pytest.mark.parametrize("found,method,status", [
    (d.Resolution(d.REDIRECT, "star_systems", 7), "GET", 301),
    (d.Resolution(d.REDIRECT, "star_systems", 7), "head", 301),
    (d.Resolution(d.REDIRECT, "star_systems", 7), "PATCH", 409),
    (d.Resolution(d.REDIRECT, "star_systems", 7), "DELETE", 409),
    (d.Resolution(d.GONE, "star_systems", 7), "GET", 410),
    (d.Resolution(d.GONE, "star_systems", 7), "PATCH", 410),
    (d.Resolution(d.CURRENT, "star_systems", 7), "PATCH", None),
    (d.Resolution(d.UNKNOWN), "GET", None),
])
def test_a_read_of_a_retired_designation_redirects_a_write_is_refused_and_a_deleted_one_is_gone(found, method, status):
    assert d.status_of(method, found) == status
