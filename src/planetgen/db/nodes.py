# planetgen/db/nodes.py

"""
The sector-block nodes of the single-space Galaxy Map (MAP.260, docs/design-drafts single-space-map-zoom section 2.3).

A node is a block of the ladder (`galaxy.drill`: a sector, or a block of 3 to 177,147 shells) and what the map needs
to draw it before it opens it into its children:

* `aggregate`: what is stored about the block as a whole, which needs no object read: the sums of its generated
  sectors (`stats`, `galaxy_block_stats`) and the light of the pre-placed objects nothing was built from (`lights`,
  `galaxy_block_lights`, MAP.259). Its nebula cover is not in it: the cover changes without a sector changing, so it is
  read live beside the node (`query.galaxy_nebula_cover`, as a stage's is), never from the node cache;
* `stars`, `generated`, `points` and `clouds`: the brightest objects of each kind inside the block, `NODE_*` of
  each (a sector lists every object it holds up to a larger budget), in the same ranking the cube tiles use
  (`query.galaxy_bright_stars_in_box`, `galaxy_generated_stars_in_box`, `galaxy_point_phenomena_in_box`,
  `galaxy_scattered_points_in_box`, `galaxy_clouds_in_box`: one key for every block, so a child lists every object
  of its parent that lies inside it, MAP.154), with the unlisted rest summed in the aggregate's lights;
* `placed`: the generated sectors of a level-3 block or of a sector, for the cell under the pointer and the charted
  tint.

An object is a member of the block it lies in: a pre-placed star and a generated star by the address of its sector,
a point by the sector its position falls in (`spherical.sector_address_at`). The tile readers take a box, so a block
is read by its bounding box (`block_box`) and the objects outside the block itself are dropped (`address_in_block`);
each reader is asked for `POOL` times the budget first so the cut still has `NODE_*` of them to keep.

Every part of a node depends only on its key and the database's contents (`query.galaxy_content_stamp`), so the
node cache keys each node by key and stamp as the tiles were.
"""

import math

from planetgen.db import block_lights, query
from planetgen.galaxy import spherical
from planetgen.galaxy.drill import (
    DRILL_LEVELS, DRILL_STEP, DrillBlock, drill_block_bounds, drill_parent, format_drill_key, parse_node_key)
from planetgen.physics.units import ly_to_pc

NODE_STARS = 96
"""int: Most pre-placed bright stars a block node lists (design section 2.3: K about 96)."""

NODE_GENERATED_STARS = 96
"""int: Most stars of generated systems a block node lists."""

NODE_POINTS = 60
"""int: Most point phenomena (placed and scattered together) a block node lists."""

NODE_CLOUDS = 120
"""int: Most nebulae and supernova remnants a block node lists."""

NODE_SECTOR_STARS = 300
"""int: Most pre-placed stars a sector node lists."""

NODE_SECTOR_GENERATED_STARS = 1200
"""int: Most generated stars a sector node lists (a sector holds a few hundred; only one deep in the bulge more)."""

NODE_SECTOR_POINTS = 200
"""int: Most point phenomena a sector node lists."""

NODE_PLACED_MAX = 250
"""int: Most generated sectors the `placed` list of a level-3 block lists (it holds at most a few dozen)."""

POOL = 4
"""int: How many times a node's budget each reader is asked for, because those outside the block itself are dropped."""

PLACED_LEVELS = (DRILL_STEP, 1)
"""tuple: The levels that list their generated sectors."""

COARSE_ABOVE = 81
"""int: A block bigger than this lists only the biggest classes of scattered points and no placed point phenomena, as
tiles coarser than 256 pc did (`query.POINT_PHENOMENON_MIN_LEVEL`)."""

MAX_NODES_PER_REQUEST = 128
"""int: Most nodes one request may ask for."""

def block_box(block, edge_pc):
    """
    The axis-aligned box `(lo, hi)` (parsecs, galaxy frame) that holds a block: the extremes of
    `radius * (sin(colatitude) * cos(longitude), sin(colatitude) * sin(longitude), cos(colatitude))` over the block's
    bounds (`drill.drill_block_bounds`), widened by a thousandth of a parsec so a point on the edge is inside.
    """
    (inner, outer), (north, south), (start, end) = drill_block_bounds(tuple(block), edge_pc)
    cos_north, cos_south = math.cos(north), math.cos(south)
    z_low = min(inner * cos_south, outer * cos_south)
    z_high = max(inner * cos_north, outer * cos_north)
    sin_low = min(math.sin(north), math.sin(south))
    sin_high = 1.0 if north <= math.pi / 2 <= south else max(math.sin(north), math.sin(south))
    rho_low, rho_high = inner * sin_low, outer * sin_high
    angles = [start, end]
    turn = math.ceil(start / (math.pi / 2)) * (math.pi / 2)
    while turn <= end:
        angles.append(turn)
        turn += math.pi / 2
    xs = [rho * math.cos(angle) for angle in angles for rho in (rho_low, rho_high)]
    ys = [rho * math.sin(angle) for angle in angles for rho in (rho_low, rho_high)]
    margin = 1e-3
    return ((min(xs) - margin, min(ys) - margin, z_low - margin), (max(xs) + margin, max(ys) + margin, z_high + margin))


def ancestor_at(address, m, memo=None):
    """The block of size `m` that holds the sector `address` (`(shell, band, wedge)`), by walking `drill_parent`;
    `memo` (a dict) keeps the walks of a request."""
    if memo is not None and (address, m) in memo:
        return memo[(address, m)]
    block = DrillBlock(1, *address)
    while block.m < m:
        block = drill_parent(block)
        if memo is not None:
            memo[(address, block.m)] = block
    return block


def address_in_block(address, block, memo=None):
    """Whether the sector `address` lies in `block` (a sector is in itself)."""
    if address is None or any(value is None for value in address):
        return False
    return ancestor_at(tuple(address), block.m, memo) == DrillBlock(*block)


def _star_address(star):
    return (star["shell_index"], star["band_index"], star["wedge_index"])


def _point_address(point, edge_pc):
    return spherical.sector_address_at((point["x"], point["y"], point["z"]), edge_pc)


def _stats(conn, block):
    """`{"generated", "stats"}` of the block's generated sectors (the stored row of a bigger block, the sector's own
    figures for a sector), or `None` when it holds none."""
    if block.m == 1:
        where, params = query._sectors_where([block])
        row = conn.execute(f"SELECT {query._STAGE_LOOK_SUMS} FROM ({query._stage_sectors_with_stats(where)}) placed",
                           params).fetchone()
        if not row or not row["n"]:
            return None
        counts = {}
        query._add_stage_sums(counts, None, row)
        return {"generated": counts[None][0], "stats": query._stage_stats(counts[None])}
    return query.galaxy_block_stats(conn, block.m, [(block.shell, block.band, block.wedge)]).get(
        (block.shell, block.band, block.wedge))


def _lights(conn, block):
    """The stored light of the block's pre-placed objects, as the centre of light and the mean temperature that the
    sums give (`block_lights`), or `None` when the block holds none or its aggregate is not kept (a sector's)."""
    if block.m == 1:
        return None
    row = conn.execute(
        f"SELECT {', '.join(block_lights.SUMS)} FROM galaxy_block_lights WHERE m = ? AND shell_index = ?"
        " AND band_index = ? AND wedge_index = ?", (block.m, block.shell, block.band, block.wedge)).fetchone()
    if row is None:
        return None
    lights = {"stars": int(row["star_count"]), "points": int(row["point_count"])}
    for kind in ("star", "point"):
        total = float(row[f"{kind}_luminosity_sol"])
        lights[f"{kind}_luminosity_sol"] = total
        if total > 0:
            lights[f"{kind}_centre"] = [round(float(row[f"{kind}_lp{axis}"]) / total, 3) for axis in "xyz"]
    if float(row["star_luminosity_sol"]) > 0:
        lights["star_temperature_k"] = round(float(row["star_lt"]) / float(row["star_luminosity_sol"]))
    return lights


def node_aggregate(conn, block):
    """What is stored about a block as a whole: `generated` and `stats` of its generated sectors (`_stats`) and `lights`
    (`_lights`), each only when there is something."""
    aggregate = {}
    stats = _stats(conn, block)
    if stats is not None:
        aggregate.update(stats)
    lights = _lights(conn, block)
    if lights is not None:
        aggregate["lights"] = lights
    return aggregate


def _budgets(block):
    """`(stars, generated, points)`: how many of each kind a node of `block` lists."""
    if block.m == 1:
        return NODE_SECTOR_STARS, NODE_SECTOR_GENERATED_STARS, NODE_SECTOR_POINTS
    return NODE_STARS, NODE_GENERATED_STARS, NODE_POINTS


def _outer_shell(conn):
    """The outermost shell anything is placed in: the outline's (`galaxy_shell`), or before a plan was made the
    furthest shell of a pre-placed star or a sector; `None` when there is nothing at all."""
    outer = conn.execute("SELECT MAX(shell_index) AS shell FROM galaxy_shell").fetchone()["shell"]
    if outer is not None:
        return int(outer)
    found = [row["shell"] for row in (
        conn.execute("SELECT MAX(shell_index) AS shell FROM bright_stars").fetchone(),
        conn.execute("SELECT MAX(shell_index) AS shell FROM sectors").fetchone()) if row["shell"] is not None]
    return max(int(shell) for shell in found) if found else 0


def galaxy_nodes(conn, node_keys):
    """
    The contents of each requested node: the interactive Galaxy Map's data source once its cut is made of nodes
    (MAP.261). Every part of the result depends only on its key and the database's contents, so a caller can cache each
    node by key and stamp.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        node_keys (list[str]): `m.shell.band.wedge` keys (`drill.parse_node_key`), at most `MAX_NODES_PER_REQUEST`.

    Returns:
        dict: `nodes` (`{key: {"aggregate", "stars", "generated", "points", "clouds", "placed"}}`, `placed` only for
            a level-3 block and a sector; `incomplete` names the pieces a statement limit cut, as a tile's does),
            `edge_pc` and `has_shape`.

    Raises:
        ValueError: On a malformed key or too many keys.
    """
    parsed = [(key, parse_node_key(key)) for key in dict.fromkeys(node_keys)]
    if len(parsed) > MAX_NODES_PER_REQUEST:
        raise ValueError(f"at most {MAX_NODES_PER_REQUEST} nodes per request, got {len(parsed)}")
    skeleton = query.get_galaxy_shape(conn)
    edge_pc = skeleton.edge_pc if skeleton is not None else ly_to_pc(query.DEFAULT_SECTOR_EDGE_LY)
    margins = query._cloud_margins_pc(conn) if parsed else {}
    outer = _outer_shell(conn) if parsed else None
    reach_pc = None if outer is None else (int(outer) + 1) * edge_pc
    sample = []
    memo = {}

    def brightest():
        if not sample:
            sample.append(query.galaxy_brightest_stars(conn))
        return sample[0]

    nodes = {}
    for key, block in parsed:
        lo, hi = block_box(block, edge_pc)
        if reach_pc is not None:     # nothing is placed past the outline: a root block reads the galaxy, not its own reach
            lo = tuple(max(v, -reach_pc) for v in lo)
            hi = tuple(min(v, reach_pc) for v in hi)
        incomplete = []

        def piece(name, compute, default):
            return query._tile_piece(conn, incomplete, name, compute, default)

        def inside(address):
            return address_in_block(address, block, memo)

        aggregate = piece("aggregate", lambda: node_aggregate(conn, block), {})
        stars_limit, generated_limit, points_limit = _budgets(block)
        stars = piece("stars", lambda: query._top_by_key([
            star for star in query.galaxy_bright_stars_in_box(
                conn, lo, hi, edge_pc, limit=stars_limit * POOL, brightest=brightest)
            if inside(_star_address(star))], stars_limit), [])
        generated = []
        generated_sectors = int(aggregate.get("generated") or 0)
        if generated_sectors:
            allowance = None if block.m == 1 else max(1, int(math.ceil(
                query.SECTOR_ALLOWANCE_FUDGE * generated_limit / float(generated_sectors))))
            generated = piece("generated", lambda: [
                star for star in query.galaxy_generated_stars_in_box(
                    conn, lo, hi, 0.0, generated_sectors, limit=generated_limit * POOL, per_sector=allowance)
                if inside(_star_address(star))][:generated_limit], [])
        coarse = block.m > COARSE_ABOVE
        points = [] if coarse else piece("points", lambda: [
            point for point in query.galaxy_point_phenomena_in_box(conn, lo, hi, limit=points_limit * POOL)
            if inside(_point_address(point, edge_pc))], [])
        scattered = piece("scattered", lambda: [
            point for point in query.galaxy_scattered_points_in_box(
                conn, lo, hi, edge_pc, points_limit * POOL, coarse)
            if inside(_point_address(point, edge_pc))], [])
        points = (points + scattered)[:points_limit]
        clouds = piece("clouds", lambda: query.galaxy_clouds_in_box(
            conn, lo, hi, max_clouds=NODE_CLOUDS, margins=margins), [])
        node = {"aggregate": aggregate, "stars": stars, "generated": generated, "points": points, "clouds": clouds}
        if block.m in PLACED_LEVELS:
            node["placed"] = piece("placed", lambda: [
                sector for sector in query.galaxy_sectors_in_box(conn, lo, hi, limit=NODE_PLACED_MAX * POOL)
                if inside((sector["shell_index"], sector["band_index"], sector["wedge_index"]))][:NODE_PLACED_MAX], [])
        if incomplete:
            node["incomplete"] = incomplete
        nodes[format_drill_key(block)] = node
        if key != format_drill_key(block):
            nodes[key] = nodes[format_drill_key(block)]
    return {"nodes": nodes, "edge_pc": edge_pc, "has_shape": skeleton is not None}


def populated_node_keys(conn, floor_m):
    """
    The keys of the nodes a visit is likely to open, for the warm-up (MAP.260): every block of level `floor_m` and
    above that holds a generated sector (`galaxy_block_stats`) or a pre-placed object nothing was built from
    (`galaxy_block_lights`), the biggest blocks first and, within a level, in key order.

    Args:
        conn: A connection to the database.
        floor_m (int): The finest level listed (a level of `DRILL_LEVELS`, at least `DRILL_STEP`; the nodes of
            sectors are asked for as they are opened, not ahead).

    Returns:
        list[str]: `m.shell.band.wedge` keys, each once.
    """
    rows = conn.execute(
        "SELECT m, shell_index, band_index, wedge_index FROM galaxy_block_stats WHERE m >= ? UNION "
        "SELECT m, shell_index, band_index, wedge_index FROM galaxy_block_lights WHERE m >= ? "
        "ORDER BY m DESC, shell_index, band_index, wedge_index", (floor_m, floor_m)).fetchall()
    return [format_drill_key(DrillBlock(r["m"], r["shell_index"], r["band_index"], r["wedge_index"])) for r in rows]
