# tests/test_named_regions.py

"""
GEN.825: the name, the anchor and the framing sphere of a region (`galaxy/named_regions.py`): a 9-block of the block ladder
named "Word Word Region" from its 81-block and its place in it, by the codec kind `region`, with no database read.
"""

import math
import random
import zlib

import pytest

from planetgen import tuning
from planetgen.galaxy import drill, named_regions, spherical
from planetgen.names import naming_key

KEY = "DEADBEEF"
OTHER_KEY = "0BADF00D"
EDGE = tuning.DEFAULT_SECTOR_EDGE_PC


def _cell(shell, band, wedge):
    """The 9-block holding the sector `(shell, band, wedge)`."""
    return drill.drill_chain_of(shell, band, wedge)[drill.DRILL_LEVELS.index(9)]


CELLS = [_cell(*spherical.sector_address_at(point, EDGE)) for point in (
    (8000.0, 15.0, 3.0), (-3000.0, 4000.0, 80.0), (20.0, -10.0, 5.0), (0.0, 0.0, 3000.0), (500.0, 500.0, -500.0),
    (12000.0, -9000.0, 40.0), (1.0, 1.0, 1.0), (0.0, 0.0, -1800.0))]


def test_the_codec_names_the_region_kind_and_the_object_kinds_are_as_before():
    assert "region" in naming_key.KINDS
    assert naming_key.codec_domain(KEY, "region") == f"{KEY}:region"


@pytest.mark.parametrize("cell", CELLS, ids=drill.format_drill_key)
def test_a_name_is_the_same_wherever_it_is_asked_and_changes_with_the_key(cell):
    name = named_regions.region_name(cell, KEY)
    assert name.endswith(" Region") and len(name.split()) == 3
    assert named_regions.region_name(tuple(cell), KEY) == name
    random.random()
    assert named_regions.region_name(cell, KEY) == name
    assert named_regions.region_name(cell, OTHER_KEY) != name


def test_golden_names():
    """A name changes only when the codec, the hash or the order of the cells does: that is an epoch of the names."""
    names = {drill.format_drill_key(cell): named_regions.region_name(cell, KEY) for cell in CELLS[:4]}
    assert names == GOLDEN


GOLDEN = {
    "9.222.0.0": "Kanbee Get Region",
    "9.138.2.342": "Tamun Nonb Region",
    "9.0.0.2": "Tuhen Nanr Region",
    "9.83.112.0": "Mooda Nanl Region",
}


def test_a_region_is_a_nine_block_only():
    with pytest.raises(ValueError):
        named_regions.region_name(drill.DrillBlock(27, 5, 0, 0), KEY)
    with pytest.raises(ValueError):
        named_regions.region_name(drill.DrillBlock(3, 5, 0, 0), KEY)


@pytest.mark.parametrize("cell", CELLS[:5], ids=drill.format_drill_key)
def test_the_regions_of_one_family_share_the_first_word_and_never_a_name(cell):
    family = named_regions.family_of(cell)
    cells = named_regions._family_cells(family)
    assert 0 < len(cells) <= 16 ** named_regions.ORDINAL_DIGITS
    assert len(set(cells)) == len(cells)
    assert all(named_regions.family_of(other) == family for other in cells)
    names = [named_regions.region_name(other, KEY) for other in cells]
    assert len(set(names)) == len(names)
    assert len({name.split()[0] for name in names}) == 1
    assert names[0].split()[0] + " Reach" == named_regions.family_name(family, KEY)
    assert [named_regions.ordinal_in_family(other) for other in cells] == list(range(len(cells)))


def test_the_sectors_of_a_cell_all_name_that_cell_and_its_family():
    for cell in CELLS[:3]:
        family = named_regions.family_of(cell)
        addresses = list(drill.block_addresses(cell))
        assert addresses
        for address in addresses[::7]:
            chain = drill.drill_chain_of(*address)
            assert cell in chain and family in chain
        assert cell in named_regions._family_cells(family)


def _sample_families(count=40):
    rng = random.Random(1)
    families = set()
    while len(families) < count:
        radius = rng.choice([30.0, 300.0, 3000.0, 9000.0, 15000.0])
        theta = math.acos(rng.uniform(-1, 1))
        phi = rng.uniform(0, 2 * math.pi)
        point = (radius * math.sin(theta) * math.cos(phi), radius * math.sin(theta) * math.sin(phi), radius * math.cos(theta))
        families.add(named_regions.family_of(_cell(*spherical.sector_address_at(point, EDGE))))
    return sorted(families)


def test_a_family_word_is_never_a_neighbours():
    """The names of neighbouring 81-blocks start with different words (a smaller-keyed neighbour's wins a tie)."""
    families = _sample_families()
    for family in families:
        neighbours = named_regions.neighbour_families(family)
        assert 8 <= len(neighbours) <= 26 and family not in neighbours
        mine = named_regions.family_digits(family, KEY)
        for other in neighbours:
            assert named_regions.family_digits(other, KEY) != mine, (family, other)


def test_a_tie_between_neighbours_is_broken_for_the_larger_key(monkeypatch):
    """With digits forced to collide often, the block with the larger key moves off every smaller neighbour's digits."""
    def few(key, text, digits):
        return f"{zlib.crc32(text.encode()) % 40:0{digits}X}"
    monkeypatch.setattr(named_regions, "_hex", few)
    named_regions._family_digits.cache_clear()
    salted = 0
    for family in _sample_families(20):
        mine = named_regions.family_digits(family, KEY)
        for other in named_regions.neighbour_families(family):
            if tuple(other) < tuple(family):
                assert mine != few(KEY, named_regions._key_text(other), named_regions.FAMILY_DIGITS)
        salted += mine != few(KEY, named_regions._key_text(family), named_regions.FAMILY_DIGITS)
    named_regions._family_digits.cache_clear()
    assert salted


def test_the_name_reads_no_database():
    for value in vars(named_regions).values():
        assert not getattr(value, "__module__", "").startswith("planetgen.db"), value
    source = open(named_regions.__file__).read()
    assert "planetgen.db" not in source and "import sqlite3" not in source and "mysql" not in source.lower()


@pytest.mark.parametrize("cell", CELLS, ids=drill.format_drill_key)
def test_the_anchor_is_inside_its_cell_and_the_framing_sphere_stays_there(cell):
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(cell, EDGE)
    x, y, z = named_regions.anchor_pc(cell, EDGE)
    r = math.sqrt(x * x + y * y + z * z)
    assert inner - 1e-9 <= r <= outer + 1e-9
    radius = named_regions.framing_radius_pc(cell, EDGE)
    assert 0 < radius <= named_regions.REGION_RADIUS_PC
    rng = random.Random(2)
    for _ in range(200):
        d = [rng.gauss(0, 1) for _ in range(3)]
        norm = math.sqrt(sum(c * c for c in d))
        px, py, pz = (a + 0.999 * radius * c / norm for a, c in zip((x, y, z), d))
        pr = math.sqrt(px * px + py * py + pz * pz)
        theta = math.atan2(math.hypot(px, py), pz)
        phi = math.atan2(py, px) % (2 * math.pi)
        assert inner - 1e-6 <= pr <= outer + 1e-6
        assert north - 1e-9 <= theta <= south + 1e-9
        if end - start < 2 * math.pi - 1e-9:
            assert start - 1e-9 <= phi <= end + 1e-9 or (end > 2 * math.pi and phi + 2 * math.pi <= end + 1e-9)


def test_a_standard_cell_is_framed_by_the_standard_fill_sphere_and_a_squeezed_one_by_what_it_can_hold():
    assert named_regions.framing_radius_pc(CELLS[0], EDGE) == named_regions.REGION_RADIUS_PC
    rows = drill.drill_band_count(9, 7000)
    squeezed = [drill.DrillBlock(9, 7000, band, 0) for band in range(-(rows // 2), rows // 2 + 1)
                if named_regions.inscribed_radius_pc(drill.DrillBlock(9, 7000, band, 0), EDGE) < named_regions.REGION_RADIUS_PC]
    assert squeezed
    for cell in squeezed:
        radius = named_regions.framing_radius_pc(cell, EDGE)
        assert 7.0 < radius == named_regions.inscribed_radius_pc(cell, EDGE) < named_regions.REGION_RADIUS_PC
