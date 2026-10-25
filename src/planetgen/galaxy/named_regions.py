# planetgen/galaxy/named_regions.py

"""
Named galaxy regions: the name, the anchor and the framing sphere of a 9-block (GEN.825;
/mnt/project-files/design-drafts/region-names-from-id-prefix-2026-10-25.md section 0.2).

A region is a 9-block cell of the Galaxy Map's block ladder (`drill`) that holds at least one generated sector; whether it
holds one is the database's to say (`region_stats`, DB.107), and nothing here reads a database. This module is the geometry
and the name, both pure functions of the block and the galaxy's naming key:

* **The anchor** of a cell is the centre of its bounds (the middle radius, colatitude and longitude): where the standard
  fill of the cell (a sphere of `REGION_RADIUS_PC`, 109 sectors) is centred.
* **The framing radius** is that standard radius, or the cell's inscribed radius where the grid squeezes the cell below it
  (next to the core and the poles), so the sphere never leaves its cell.
* **The name** is "Word Word Region": the codec kind `region` on seven hex digits, the first four from a hash of the
  cell's 81-block (the family: every cell of one 81-block starts with the same word) and the last three the cell's ordinal
  among the 9-blocks of its 81-block (so two cells of one family never share a name). The 81-block alone is "Word Reach".
  The codec spells the first four digits as the first word and the last three as the second (`names/gated_phoneme_codec.py`).
  A family word that equals a neighbouring family's is salted (`family_digits`), so the 81-blocks around one never share it.

Two blocks 324 pc apart sharing a family word is harmless; two neighbours doing so is what the salt prevents.
"""

import functools
import hashlib
import math

from planetgen import tuning
from planetgen.galaxy import drill, spherical
from planetgen.names import naming_key
from planetgen.names.gated_phoneme_codec import GatedPhonemeCodec

REGION_RADIUS_PC = tuning.DEFAULT_GENERATE_RADIUS_PC
"""float: The radius of a region's standard sphere, parsecs: the sphere a standard fill of one anchor fills (about 109 sectors)
and the smallest sphere a region is drawn with. The cell it sits in is 36 pc across at the default sector edge, so the sphere
never leaves it; where the grid squeezes the cell the radius is the cell's inscribed radius (`framing_radius_pc`)."""

REGION_LEVEL = 9
"""int: The block size of a region (sector shells deep): about 780 sectors, 36 pc across at the default edge."""

FAMILY_LEVEL = 81
"""int: The block size whose hash gives the first word of the name of every region inside it."""

FAMILY_DIGITS = 4
"""int: Hex digits of the family hash: the codec spells them as the first word."""

ORDINAL_DIGITS = 3
"""int: Hex digits of a cell's ordinal in its family: the codec spells them as the second word."""

MAX_SALT = 4095
"""int: The salts tried before the rule gives up (it needs one per neighbour at most, 26)."""

_CODEC = GatedPhonemeCodec()

REGION_KIND = "region"
"""str: The codec kind of region names (`naming_key.KINDS`)."""


def _check_region(block):
    block = drill.DrillBlock(*block)
    if block.m != REGION_LEVEL:
        raise ValueError(f"a region is a {REGION_LEVEL}-block, got size {block.m}")
    return block


def family_of(block):
    """The 81-block that holds the 9-block `block`."""
    block = _check_region(block)
    return drill.drill_parent(drill.drill_parent(block))


@functools.lru_cache(maxsize=4096)
def _family_cells(family):
    """The 9-blocks of an 81-block in the ladder's order: shells outward, then rows north to south, then wedges."""
    cells = []
    for _shell, mids in drill.drill_children(family):
        for mid in mids:
            for _child_shell, children in drill.drill_children(mid):
                cells.extend(children)
    return tuple(cells)


def ordinal_in_family(block):
    """The place of the 9-block `block` among the 9-blocks of its 81-block, from 0, in the ladder's order."""
    block = _check_region(block)
    return _family_cells(family_of(block)).index(block)


def _key_text(block):
    return drill.format_drill_key(block)


def _hex(key, text, digits):
    return hashlib.sha256(f"{naming_key.parse_key(key)}:{REGION_KIND}:{text}".encode()).hexdigest()[:digits].upper()


def anchor_pc(block, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """
    The anchor of the 9-block `block`, `(x, y, z)` parsecs in the galaxy frame: the centre of its bounds, the middle radius
    (the core cell's origin side is the origin), colatitude and longitude (`drill.drill_block_bounds`).
    """
    block = _check_region(block)
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(block, edge_pc)
    radius, colatitude, longitude = (inner + outer) / 2, (north + south) / 2, (start + end) / 2
    plane = radius * math.sin(colatitude)
    return (plane * math.cos(longitude), plane * math.sin(longitude), radius * math.cos(colatitude))


def inscribed_radius_pc(block, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """
    How far the anchor of the 9-block `block` is from the nearest face of the cell, parsecs: the smaller of half the shell
    depth, the distance to the colatitude cones and the distance to the longitude half-planes (infinite for a cell that goes
    all the way round). The grid squeezes a cell near the core and the poles, where a row has few wedges.
    """
    block = _check_region(block)
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(block, edge_pc)
    radius, colatitude = (inner + outer) / 2, (north + south) / 2
    reach = [(outer - inner) / 2, radius * math.sin((south - north) / 2)]
    arc = end - start
    if arc < 2 * math.pi * (1 - 1e-9):
        reach.append(radius * math.sin(colatitude) * math.sin(min(arc / 2, math.pi / 2)))
    return min(reach)


def framing_radius_pc(block, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The radius of the smallest sphere a region is drawn with and the standard fill of its cell fills: `REGION_RADIUS_PC`,
    or the cell's inscribed radius where that is smaller. No cell is left without one: a cell with a generated sector is a
    region."""
    return min(REGION_RADIUS_PC, inscribed_radius_pc(block, edge_pc))


@functools.lru_cache(maxsize=8192)
def _displaced_families(family, edge_pc):
    """The 81-blocks around `family`: those that hold the points one block's width from its anchor along each of the 26
    directions of the cube (diagonals included), the block's own excluded. A cell the grid squeezes may miss a neighbour; the
    block it finds are the ones beside it in the ladder."""
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(family, edge_pc)
    radius, colatitude, longitude = (inner + outer) / 2, (north + south) / 2, (start + end) / 2
    plane = radius * math.sin(colatitude)
    x, y, z = plane * math.cos(longitude), plane * math.sin(longitude), radius * math.cos(colatitude)
    width = FAMILY_LEVEL * edge_pc
    found = set()
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == dy == dz == 0:
                    continue
                point = (x + dx * width, y + dy * width, z + dz * width)
                if math.sqrt(sum(c * c for c in point)) > spherical.MAX_ADDRESSABLE_SHELL * edge_pc:
                    continue
                shell, band, wedge = spherical.sector_address_at(point, edge_pc)
                found.add(drill.drill_chain_of(shell, band, wedge)[drill.DRILL_LEVELS.index(FAMILY_LEVEL)])
    found.discard(family)
    return found


def neighbour_families(family, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The 81-blocks around the 81-block `family`, sorted (see `_displaced_families`)."""
    return sorted(_displaced_families(drill.DrillBlock(*family), float(edge_pc)))


def family_digits(family, key, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """
    The four hex digits that spell the first word of the names of the regions in the 81-block `family`: the hash of the
    block under `key`, salted (1, 2, ...) while it equals the unsalted digits of a neighbouring 81-block with a smaller key,
    so of two neighbours that would share a word the one with the larger key changes. Pure geometry: it reads no database
    and does not depend on what has been generated.
    """
    return _family_digits(drill.DrillBlock(*family), naming_key.parse_key(key), edge_pc)


@functools.lru_cache(maxsize=8192)
def _family_digits(family, key, edge_pc):
    if family.m != FAMILY_LEVEL:
        raise ValueError(f"a family is an {FAMILY_LEVEL}-block, got size {family.m}")
    taken = {_hex(key, _key_text(other), FAMILY_DIGITS) for other in neighbour_families(family, edge_pc)
             if tuple(other) < tuple(family)}
    for salt in range(MAX_SALT + 1):
        digits = _hex(key, _key_text(family) + (f":{salt}" if salt else ""), FAMILY_DIGITS)
        if digits not in taken:
            return digits
    raise RuntimeError(f"no family digits free of the neighbours of {_key_text(family)} in {MAX_SALT} salts")


def _word(digits, key):
    return _CODEC.encode(digits, domain=naming_key.codec_domain(key, REGION_KIND)).title()


def region_name(block, key, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """
    The name of the 9-block `block` under naming key `key`: "Fardoo Beez Region". The same for the same block and key
    wherever it is asked; the regions of one 81-block start with the same word and never share a name.
    """
    block = _check_region(block)
    ordinal = ordinal_in_family(block)
    digits = family_digits(family_of(block), key, edge_pc) + f"{ordinal:0{ORDINAL_DIGITS}X}"
    return f"{_word(digits, key)} Region"


def family_name(family, key, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The name of the 81-block `family`: the first word of its regions, "Fardoo Reach"."""
    digits = family_digits(family, key, edge_pc) + "0" * ORDINAL_DIGITS
    return f"{_word(digits, key).split(' ')[0]} Reach"
