# tests/test_member_designation.py

"""GEN.1065: the designation of a belt member, `<belt ID>.<rank hex>` (`galaxy/object_uid.designation`, docs/design/systems/belt-evolution.md section 9)."""

import pytest

from planetgen.galaxy import object_uid as ids
from planetgen.names import naming_key

BELT = ids.pack(250, 7, 3, ids.SERIAL_GENERATED, 5, body=3)
SYSTEM = ids.pack(250, 7, 3, ids.SERIAL_GENERATED, 5, body=0)


def test_a_designation_is_the_belts_printed_id_a_point_and_the_rank_in_hex():
    assert ids.designation(BELT, 0x42) == f"{ids.format_id(BELT)}.0042"
    assert ids.designation(BELT, 1) == f"{ids.format_id(BELT)}.0001"
    assert ids.designation(BELT, 0xFFFF) == f"{ids.format_id(BELT)}.FFFF"
    assert ids.designation(BELT, 0x10000) == f"{ids.format_id(BELT)}.10000"
    assert ids.designation(BELT, ids.MAX_MEMBER_RANK) == f"{ids.format_id(BELT)}.FFFFFFFF"


def test_it_is_a_pure_function_of_the_two_numbers():
    first = [ids.designation(BELT, rank) for rank in (1, 2, 300, 100000)]
    assert first == [ids.designation(BELT, rank) for rank in (1, 2, 300, 100000)]
    assert len(set(first)) == 4
    assert ids.designation(BELT, 5) != ids.designation(ids.with_body(BELT, 4), 5)


def test_it_reads_no_stream_and_no_naming_key():
    from planetgen.util import draw
    with draw.bound("designation-test"):
        before = draw.random()
    with draw.bound("designation-test"):
        ids.designation(BELT, 12)
        assert draw.random() == before


def test_parse_designation_round_trips_it_through_parse_id():
    for belt in (BELT, ids.with_body(BELT, 4095), ids.pack(0, 0, 0, ids.SERIAL_RUNTIME, 0, body=1)):
        for rank in (1, 2, 0xFF, 0xFFFF, 0x10000, 299999, ids.MAX_MEMBER_RANK):
            text = ids.designation(belt, rank)
            assert ids.parse_designation(text) == (belt, rank)
            assert ids.parse_id(text.split(".")[0]) == belt
            assert ids.designation(*ids.parse_designation(text)) == text


def test_parse_designation_ignores_case_and_refuses_a_second_text_for_one_member():
    text = ids.designation(BELT, 0xAB)
    assert ids.parse_designation(text.lower()) == (BELT, 0xAB)
    for bad in (text.replace(".00AB", ".AB"), text.replace(".00AB", ".000AB"), text.replace(".00AB", ".00ab ").strip() + "x", text + ".1",
                ids.format_id(BELT), text.replace(".", ","), text.replace(".00AB", "."), text.replace(".00AB", ".0000"), ids.format_id(BELT) + ".GGGG"):
        with pytest.raises(ValueError):
            ids.parse_designation(bad)


def test_parse_id_still_refuses_a_designation():
    with pytest.raises(ValueError):
        ids.parse_id(ids.designation(BELT, 3))


def test_a_system_id_and_an_out_of_range_rank_are_refused():
    with pytest.raises(ValueError, match="belt"):
        ids.designation(SYSTEM, 1)
    for rank in (0, -1, ids.MAX_MEMBER_RANK + 1, 1.5, True, "3"):
        with pytest.raises(ValueError):
            ids.designation(BELT, rank)
    with pytest.raises(ValueError):
        ids.parse_designation(f"{ids.format_id(SYSTEM)}.0001")


def test_a_rank_above_the_individual_count_addresses_a_body_drawn_on_demand():
    assert not ids.is_on_demand(300, 300) and ids.is_on_demand(301, 300)
    assert not ids.is_on_demand(1, 1) and ids.is_on_demand(1, 0)
    # the designation of an on-demand body is the same kind of text, so it can be addressed without being stored
    assert ids.parse_designation(ids.designation(BELT, 301)) == (BELT, 301)


def test_belt_members_are_left_out_of_the_naming_key():
    """No salad name, no codec name: the key renames only the kinds it lists, and a belt member is none of them."""
    assert not [kind for kind in naming_key.KINDS if "belt" in kind or "member" in kind]


# --- captured bodies (TEST.495) -------------------------------------------------------------------------------------------------

PLANET = ids.with_body(BELT, 1)
OTHER_SYSTEM_PLANET = ids.pack(250, 7, 3, ids.SERIAL_GENERATED, 6, body=1)


def test_a_captured_belt_member_keeps_its_belt_designation_and_gains_a_parent_link():
    captured = ids.captured_member(BELT, 0x42, PLANET)
    assert captured.designation == ids.designation(BELT, 0x42)           # capture renumbers nothing
    assert captured.captured_by == ids.format_id(PLANET)
    assert ids.parse_captured(captured) == (BELT, 0x42, PLANET)


def test_an_outside_reservoir_body_prints_the_planets_id_and_a_rank():
    captured = ids.outside_reservoir_member(PLANET, 7)
    assert captured.designation == f"{ids.format_id(PLANET)}.0007" and captured.captured_by == ids.format_id(PLANET)
    parent, rank, planet = ids.parse_captured(captured)
    assert (parent, rank, planet) == (PLANET, 7, PLANET)
    assert ids.body_of(ids.parse_id(captured.designation.split(".")[0])) == ids.body_of(PLANET)      # the 12-bit body field is untouched


def test_both_are_pure_and_round_trip_through_parse_id():
    for belt, planet in ((BELT, PLANET), (ids.with_body(BELT, 4095), ids.with_body(BELT, 4094))):
        for rank in (1, 300, 0x10000, ids.MAX_MEMBER_RANK):
            a, b = ids.captured_member(belt, rank, planet), ids.captured_member(belt, rank, planet)
            assert a == b and ids.parse_captured(a) == (belt, rank, planet)
            c = ids.outside_reservoir_member(planet, rank)
            assert c == ids.outside_reservoir_member(planet, rank) and ids.parse_captured(c) == (planet, rank, planet)
            assert ids.parse_id(a.captured_by) == planet


def test_a_capture_across_systems_or_by_the_belt_itself_or_by_a_system_id_is_refused():
    with pytest.raises(ValueError, match="not in the system"):
        ids.captured_member(BELT, 1, OTHER_SYSTEM_PLANET)
    with pytest.raises(ValueError, match="own members"):
        ids.captured_member(BELT, 1, BELT)
    with pytest.raises(ValueError, match="system's own ID"):
        ids.captured_member(BELT, 1, SYSTEM)
    with pytest.raises(ValueError, match="system's own ID"):
        ids.outside_reservoir_member(SYSTEM, 1)
    with pytest.raises(ValueError):
        ids.outside_reservoir_member(PLANET, 0)


def test_a_captured_link_to_another_system_does_not_parse():
    bad = ids.Captured(ids.designation(BELT, 3), ids.format_id(OTHER_SYSTEM_PLANET))
    with pytest.raises(ValueError, match="not in the system"):
        ids.parse_captured(bad)


def test_same_system_ignores_only_the_body_number():
    assert ids.same_system(BELT, PLANET) and ids.same_system(BELT, SYSTEM)
    assert not ids.same_system(BELT, OTHER_SYSTEM_PLANET)
