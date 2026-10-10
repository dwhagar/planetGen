# tests/test_object_uid_layout.py

"""GEN.170: the 80-bit object ID (`galaxy/object_uid.py`): fixed length, a
round trip through every form, and a layout that grows with the galaxy."""

import pytest

from planetgen.galaxy import object_uid as oid

EDGES = [
    (0, 0, 0, oid.SERIAL_GENERATED, 0, 0),
    (3762, -1020, 23039, oid.SERIAL_GENERATED, 67_108_863, 4095),
    (1, 1020, 5, oid.SERIAL_RUNTIME, 12345, 1),
    (4095, 2047, 65535, oid.SERIAL_FIELD, (1 << 26) - 1, 4095),
    (4095, -2048, 0, oid.SERIAL_RUNTIME, 0, 0),
]


def test_the_default_layout_is_80_bits_and_ten_bytes():
    layout = oid.DEFAULT_LAYOUT
    assert (layout.sector_bits, layout.total_bits, layout.length) == (40, 80, 10)
    assert layout.sector_bits // 4 == 10 and layout.serial_bits // 4 == 7 and layout.body_bits // 4 == 3


@pytest.mark.parametrize("fields", EDGES)
def test_an_id_round_trips_through_every_form(fields):
    value = oid.pack(*fields)
    assert tuple(oid.unpack(value)) == fields
    text = oid.format_id(value)
    assert len(text) == 22 and text.count("-") == 2
    assert oid.parse_id(text) == value == oid.parse_id(text.replace("-", "").lower())
    raw = oid.to_bytes(value)
    assert len(raw) == 10 and oid.from_bytes(raw) == value


def test_every_id_has_the_same_length():
    values = [oid.pack(*fields) for fields in EDGES]
    assert {len(oid.format_id(value).replace("-", "")) for value in values} == {20}
    assert {len(oid.to_bytes(value)) for value in values} == {10}


def test_the_printed_form_groups_sector_serial_and_body():
    value = oid.pack(0xFE8 >> 0, -2048 + 0x100, 0x0A2B, oid.SERIAL_GENERATED, 5, 0x1A)
    sector, serial, body = oid.format_id(value).split("-")
    assert (len(sector), len(serial), len(body)) == (10, 7, 3)
    assert serial == "0000005" and body == "01A"
    assert oid.format_id(oid.pack(0, 0, 0, oid.SERIAL_RUNTIME, 3)).split("-")[1] == "4000003"


def test_the_ids_of_one_sector_sort_together():
    a = oid.pack(5, 0, 7, oid.SERIAL_GENERATED, 1, 0)
    b = oid.pack(5, 0, 7, oid.SERIAL_GENERATED, 1, 9)
    c = oid.pack(5, 0, 8, oid.SERIAL_GENERATED, 0, 0)
    assert a < b < c


@pytest.mark.parametrize("kwargs", [
    dict(ring=4096, layer=0, slot=0, serial_kind=0, serial=0, body=0),
    dict(ring=0, layer=2048, slot=0, serial_kind=0, serial=0, body=0),
    dict(ring=0, layer=-2049, slot=0, serial_kind=0, serial=0, body=0),
    dict(ring=0, layer=0, slot=65536, serial_kind=0, serial=0, body=0),
    dict(ring=0, layer=0, slot=0, serial_kind=0, serial=1 << 26, body=0),
    dict(ring=0, layer=0, slot=0, serial_kind=0, serial=0, body=4096),
    dict(ring=-1, layer=0, slot=0, serial_kind=0, serial=0, body=0),
    dict(ring=0, layer=0, slot=0, serial_kind=3, serial=0, body=0),
])
def test_a_field_that_does_not_fit_is_refused(kwargs):
    with pytest.raises(ValueError):
        oid.pack(**kwargs)


@pytest.mark.parametrize("text", ["", "xyz", "FE81000A2B-0000005", "FE81000A2B-000005-0000", "FE81000A2B00000050",
                                  "FE81000A2B-0000005-000-1", "FE81000A2B-0000005-00G"])
def test_text_that_is_not_an_id_is_refused(text):
    with pytest.raises(ValueError):
        oid.parse_id(text)


def test_bytes_of_the_wrong_length_are_refused():
    with pytest.raises(ValueError):
        oid.from_bytes(bytes(12))


def test_a_serial_kind_that_is_not_known_is_refused_on_unpack():
    with pytest.raises(ValueError):
        oid.unpack(3 << (oid.DEFAULT_LAYOUT.serial_count_bits + oid.BODY_BITS))


def test_the_milky_way_fits_the_default_layout():
    # 3,763 rings, 2,041 layers, 23,040 slots in the outer ring (docs/design/object-id-options.md).
    assert oid.fit_layout(3762, -1020, 1020, 23040) == oid.DEFAULT_LAYOUT


@pytest.mark.parametrize("bounds,total", [
    ((4096, -10, 10, 100), 96),       # one ring too many
    ((100, -2050, 2050, 100), 96),    # layers past 11 bits either side
    ((100, -10, 10, 70000), 96),      # a ring with more than 65,536 slots
    ((10 ** 5, -5000, 5000, 10 ** 6), 96),
    ((10 ** 6, -5000, 5000, 10 ** 7), 128),
    ((2 ** 30, -(2 ** 18), 2 ** 18, 2 ** 24), 128),
])
def test_a_galaxy_that_does_not_fit_gets_a_longer_id(bounds, total):
    layout = oid.fit_layout(*bounds)
    assert layout.total_bits == total and layout.total_bits % 8 == 0 and layout.sector_bits % 4 == 0
    max_ring, min_layer, max_layer, max_slots = bounds
    for fields in [(max_ring, min_layer, max_slots - 1), (0, max_layer, 0)]:
        value = oid.pack(*fields, oid.SERIAL_GENERATED, 7, 3, layout=layout)
        assert tuple(oid.unpack(value, layout))[:3] == fields
        assert len(oid.to_bytes(value, layout)) == total // 8
        assert oid.parse_id(oid.format_id(value, layout), layout) == value


def test_a_galaxy_too_big_for_128_bits_is_refused():
    with pytest.raises(ValueError):
        oid.fit_layout(2 ** 60, -(2 ** 40), 2 ** 40, 2 ** 40)
