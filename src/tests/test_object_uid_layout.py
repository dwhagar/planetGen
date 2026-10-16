# tests/test_object_uid_layout.py

"""GEN.170, GEN.230: the 96-bit object ID (`galaxy/object_uid.py`): fixed length, a
round trip through every form, and a layout that holds every cell of the grid's reach."""

import pytest

from planetgen.galaxy import object_uid as oid

EDGES = [
    (0, 0, 0, oid.SERIAL_GENERATED, 0, 0),
    (3762, -1020, 23039, oid.SERIAL_GENERATED, 67_108_863, 4095),
    (1, 1020, 5, oid.SERIAL_RUNTIME, 12345, 1),
    (65535, 102942, 411599, oid.SERIAL_FIELD, (1 << 26) - 1, 4095),
    (65535, -102942, 0, oid.SERIAL_RUNTIME, 0, 0),
    ((1 << 19) - 1, (1 << 17) - 1, (1 << 19) - 1, oid.SERIAL_GENERATED, 7, 1),
    ((1 << 19) - 1, -(1 << 17), 0, oid.SERIAL_GENERATED, 7, 1),
]


def test_the_layout_is_96_bits_and_twelve_bytes():
    layout = oid.LAYOUT
    assert (layout.shell_bits, layout.band_bits, layout.wedge_bits) == (19, 18, 19)
    assert (layout.sector_bits, layout.total_bits, layout.length) == (56, 96, 12)
    assert layout.sector_bits // 4 == 14 and layout.serial_bits // 4 == 7 and layout.body_bits // 4 == 3


@pytest.mark.parametrize("fields", EDGES)
def test_an_id_round_trips_through_every_form(fields):
    value = oid.pack(*fields)
    assert tuple(oid.unpack(value)) == fields
    text = oid.format_id(value)
    assert len(text) == 26 and text.count("-") == 2
    assert oid.parse_id(text) == value == oid.parse_id(text.replace("-", "").lower())
    raw = oid.to_bytes(value)
    assert len(raw) == 12 and oid.from_bytes(raw) == value


def test_every_id_has_the_same_length():
    values = [oid.pack(*fields) for fields in EDGES]
    assert {len(oid.format_id(value).replace("-", "")) for value in values} == {24}
    assert {len(oid.to_bytes(value)) for value in values} == {12}


def test_the_printed_form_groups_sector_serial_and_body():
    value = oid.pack(0xFE8 >> 0, -131072 + 0x100, 0x0A2B, oid.SERIAL_GENERATED, 5, 0x1A)
    sector, serial, body = oid.format_id(value).split("-")
    assert (len(sector), len(serial), len(body)) == (14, 7, 3)
    assert serial == "0000005" and body == "01A"
    assert oid.format_id(oid.pack(0, 0, 0, oid.SERIAL_RUNTIME, 3)).split("-")[1] == "4000003"


def test_the_ids_of_one_sector_sort_together():
    a = oid.pack(5, 0, 7, oid.SERIAL_GENERATED, 1, 0)
    b = oid.pack(5, 0, 7, oid.SERIAL_GENERATED, 1, 9)
    c = oid.pack(5, 0, 8, oid.SERIAL_GENERATED, 0, 0)
    assert a < b < c


@pytest.mark.parametrize("kwargs", [
    dict(shell=1 << 19, band=0, wedge=0, serial_kind=0, serial=0, body=0),
    dict(shell=0, band=1 << 17, wedge=0, serial_kind=0, serial=0, body=0),
    dict(shell=0, band=-(1 << 17) - 1, wedge=0, serial_kind=0, serial=0, body=0),
    dict(shell=0, band=0, wedge=1 << 19, serial_kind=0, serial=0, body=0),
    dict(shell=0, band=0, wedge=0, serial_kind=0, serial=1 << 26, body=0),
    dict(shell=0, band=0, wedge=0, serial_kind=0, serial=0, body=4096),
    dict(shell=-1, band=0, wedge=0, serial_kind=0, serial=0, body=0),
    dict(shell=0, band=0, wedge=0, serial_kind=3, serial=0, body=0),
])
def test_a_field_that_does_not_fit_is_refused(kwargs):
    with pytest.raises(ValueError):
        oid.pack(**kwargs)


@pytest.mark.parametrize("text", ["", "xyz", "FE8100000A2B-0000005", "FE8100000A2B-000005-0000", "FE8100000A2B00000050",
                                  "FE8100000A2B-0000005-000-1", "FE8100000A2B-0000005-00G",
                                  "FE81000A2B-0000005-000"])
def test_text_that_is_not_an_id_is_refused(text):
    with pytest.raises(ValueError):
        oid.parse_id(text)


def test_bytes_of_the_wrong_length_are_refused():
    with pytest.raises(ValueError):
        oid.from_bytes(bytes(10))


def test_a_serial_kind_that_is_not_known_is_refused_on_unpack():
    with pytest.raises(ValueError):
        oid.unpack(3 << (oid.LAYOUT.serial_count_bits + oid.BODY_BITS))


def test_the_milky_way_and_the_whole_reach_fit_the_layout():
    # Shell 3,764 is where the default outline ends, bands -1602..1602 at the outline; shell 65,535 is the 262 kpc
    # reach, whose polar cap is band 102,942 and whose equatorial band holds 411,600 wedges.
    from planetgen.galaxy import skeleton, spherical
    spherical.check_reach(skeleton.DEFAULT_REACH_SHELL)
    assert spherical.band_range(65535)[1] == 102942 and spherical.wedge_count(65535, 0) == 411600
    value = oid.pack(65535, -102942, 411599, oid.SERIAL_GENERATED, 7, 3)
    assert tuple(oid.unpack(value))[:3] == (65535, -102942, 411599)


def test_the_last_addressable_shell_is_the_last_one_the_band_field_holds():
    from planetgen.galaxy import spherical
    spherical.check_reach(spherical.MAX_ADDRESSABLE_SHELL)
    with pytest.raises(ValueError):
        spherical.check_reach(spherical.MAX_ADDRESSABLE_SHELL + 1)
    with pytest.raises(ValueError):
        spherical.check_reach(-1)
