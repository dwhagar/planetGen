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
* **The extent** (GEN.827) is the sphere a region is drawn and framed with: the bounding sphere of the largest face-connected
  clump of its generated sectors (`region_extent`), never smaller than the framing sphere. Every other clump of the cell
  still counts toward the share charted; it is not inside the sphere (no line joins two unrelated places).
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
from collections import namedtuple

import numpy as np

from planetgen import tuning
from planetgen.galaxy import drill, spherical
from planetgen.names import naming_key
from planetgen.names.gated_phoneme_codec import GatedPhonemeCodec
from planetgen.util import detmath

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


def region_of(address):
    """The region (9-block) that holds the sector `(shell, band, wedge)`."""
    shell, band, wedge = (int(value) for value in address)
    return drill.drill_chain_of(shell, band, wedge)[drill.DRILL_LEVELS.index(REGION_LEVEL)]


def anchor_address(block, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """
    The sector `(shell, band, wedge)` that holds the anchor of the 9-block `block`: where the standard fill of the cell (the
    sphere of `REGION_RADIUS_PC` about the anchor) is centred (GEN.826). The anchor is the middle of the cell's bounds, so the
    sector there is in the cell; where the grid's rounding puts it in a neighbour, the sector of the cell nearest the anchor.
    """
    block = _check_region(block)
    anchor = anchor_pc(block, edge_pc)
    address = spherical.sector_address_at(anchor, edge_pc)
    if region_of(address) == block:
        return tuple(address)
    return min(drill.block_addresses(block),
               key=lambda found: (detmath.dist(spherical.sector_position_pc(*found, edge_pc), anchor), found))


def snap_to_anchor(address, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The sector the standard fill of the region holding `address` is centred on (`anchor_address` of `region_of`)."""
    return anchor_address(region_of(address), edge_pc)


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


SPHERE_ITERATIONS = 400
"""int: Steps of the bounding-sphere search (`bounding_sphere`): the sphere's radius is within 1 / sqrt(steps) = 5 percent of
the smallest one that holds the points."""

Extent = namedtuple("Extent", ["centre_pc", "radius_pc", "clumps", "sectors", "sectors_in_sphere"])
"""The extent of a region (`region_extent`): the sphere's centre `(x, y, z)` and radius in parsecs, how many face-connected
clumps the cell's generated sectors form, how many sectors there are, and how many of them the sphere holds."""


def bounding_sphere(centres_pc, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """
    A sphere that holds the sectors whose centres are `centres_pc` (an N by 3 array, parsecs): the smallest sphere through the
    centres (found by the Badoiu and Clarkson steps, within 5 percent of the smallest) with the circumradius of a sector,
    `sqrt(3) / 2` edges, added so the whole of each sector is inside. Deterministic.

    Returns:
        tuple: `(centre, radius)` as `(x, y, z)` and parsecs.
    """
    points = np.asarray(centres_pc, dtype=float).reshape(-1, 3)
    if not len(points):
        raise ValueError("a bounding sphere needs at least one sector")
    centre = points[0].copy()
    for step in range(1, SPHERE_ITERATIONS):
        farthest = points[int(np.argmax(((points - centre) ** 2).sum(axis=1)))]
        centre += (farthest - centre) / (step + 1)
    radius = float(np.sqrt(((points - centre) ** 2).sum(axis=1).max())) + math.sqrt(3) / 2 * edge_pc
    return (float(centre[0]), float(centre[1]), float(centre[2])), radius


def clumps_of(addresses):
    """The face-connected groups of the sector addresses `(shell, band, wedge)` (neighbours by `spherical.neighbor_addresses`),
    largest first, ties by the smallest address."""
    left = {tuple(address) for address in addresses}
    groups = []
    while left:
        start = min(left)
        left.remove(start)
        group, stack = [start], [start]
        while stack:
            for other in spherical.neighbor_addresses(*stack.pop()):
                if other in left:
                    left.remove(other)
                    group.append(other)
                    stack.append(other)
        groups.append(sorted(group))
    return sorted(groups, key=lambda group: (-len(group), group[0]))


def region_extent(block, addresses, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """
    The extent of the 9-block `block` holding the generated sectors `addresses`: the bounding sphere of the largest
    face-connected clump (ties by the smallest address), at least `framing_radius_pc` in radius. A sphere smaller than that is
    the framing radius about the clump's sphere centre, so a lone sector gets a sensible frame; a chain of sectors is one
    clump, and the clumps of a cell that were made apart are not joined by a sphere that covers the gap between them.

    Raises:
        ValueError: With no sector.
    """
    block = _check_region(block)
    groups = clumps_of(addresses)
    if not groups:
        raise ValueError("a region needs at least one generated sector")
    centres = np.array([spherical.sector_position_pc(*address, edge_pc) for address in groups[0]])
    centre, radius = bounding_sphere(centres, edge_pc)
    radius = max(radius, framing_radius_pc(block, edge_pc))
    held = int((((np.array([spherical.sector_position_pc(*a, edge_pc) for g in groups for a in g]) - centre) ** 2)
                .sum(axis=1) <= radius ** 2).sum())
    return Extent(centre, radius, len(groups), sum(len(group) for group in groups), held)
