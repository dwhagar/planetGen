"""
planetgen/galaxy/region_run.py

The plan of a region run (GEN.1128; docs/design/systems/region-object-store.md sections 3, 4, 8.5 and 9): a sphere of a given radius from a
start point becomes whole sectors (`region.Region`, a sector is in when its centre is within the radius plus `BALL_EDGE_ALLOWANCE` edges, so a
sphere that touches a sector holds it), cut into units of at most `BLOCK_SECTORS` consecutive wedges of one (shell, band) arc, nearest unit first.

Nothing here reads a stream or the database, and no draw may read any of it: the tile is the sector and every draw is keyed by the sector's own
address (GEN.39), so the same seed gives the same objects whatever the start point, the radius, the cut or the order the regions are run in. A
plan only says which sectors to build and in what order.

`next_units` is the barrier of section 8.5: every unit of stages 0 to 2 first, then (once GEN.1120 turns it on) the units of stage 3, each reading
the final stored state. `skip` takes the sectors already charted away, which is how a stopped run resumes and how a second, overlapping region
costs only its new sectors.
"""
import math
from collections import namedtuple

from planetgen.galaxy import region as region_module
from planetgen.util import detmath

BLOCK_SECTORS = 27
"""int: Sectors in a unit, at most: a block of 27 (about 1,700 pc3 at a 4 pc edge) is the first guess of the design (section 4), a unit being the
retry and resume grain, to be sized from the measured seconds per sector (the PERF items own that number)."""

STAGE_GENERATE = "generate"
STAGE_TOPUP = "topup"
STAGE_DONE = "done"

Unit = namedtuple("Unit", ["shell", "band", "first", "count", "nearest_pc"])
"""One unit of work: `count` consecutive wedges from `first` (wrapping past the band's last wedge) of `band` of `shell`, and the distance of its
nearest sector from the start point. A unit's key is `(shell, band, first, count)`."""

Plan = namedtuple("Plan", ["center_pc", "radius_pc", "grow_shells", "reach_pc", "edge_pc", "block_sectors", "units", "sectors", "skipped"])
"""A region run: the start, the requested radius, the shells added (`grow_shells`), the reach the sectors were taken to, the units in the order to
run them, how many sectors they hold and how many the region held that `skip` took out."""


def key(unit):
    """A unit's identity, the same whatever the run it is in."""
    return unit.shell, unit.band, unit.first, unit.count


def reach_pc(radius_pc, edge_pc, grow_shells=0):
    """How far from the start a sector's centre may be: the radius, `grow_shells` shells of one edge more (section 9: the charted sphere grows
    by whole shells) and the sector allowance."""
    if radius_pc < 0:
        raise ValueError(f"radius_pc must be >= 0, got {radius_pc}")
    if grow_shells < 0 or int(grow_shells) != grow_shells:
        raise ValueError(f"grow_shells must be a whole number >= 0, got {grow_shells}")
    return radius_pc + (grow_shells + region_module.BALL_EDGE_ALLOWANCE) * edge_pc


def plan(center_pc, radius_pc, edge_pc, block_sectors=BLOCK_SECTORS, grow_shells=0, skip=(), bounds=None):
    """
    The units of a region run.

    Args:
        center_pc (tuple): `(x, y, z)` parsecs, the start point.
        radius_pc (float): `>= 0`, the sphere's radius before the sector allowance and the grown shells.
        edge_pc (float): The sector edge, `> 0`.
        block_sectors (int): `>= 1`, the most sectors in a unit.
        grow_shells (int): Whole shells of one edge added to the radius.
        skip (container): `(shell, band, wedge)` addresses already charted; they are left out and split the arcs they sit in.
        bounds (skeleton.GalaxyBounds, optional): Leave out the cells beyond the galaxy.

    Returns:
        Plan: `units` nearest first, ties by shell, band and first wedge.

    Raises:
        ValueError: On a negative radius or shell count, a non-positive edge or a block under one.
    """
    if block_sectors < 1 or int(block_sectors) != block_sectors:
        raise ValueError(f"block_sectors must be a whole number >= 1, got {block_sectors}")
    reach = reach_pc(radius_pc, edge_pc, grow_shells)
    ball = region_module.Region(center_pc, reach, edge_pc)
    units = []
    held = 0
    skipped = 0
    for arc in ball.band_arcs():
        if bounds is not None and not bounds.qualifies(arc.shell, arc.band):
            continue
        wedges = ball.wedge_count(arc.shell, arc.band)
        run = []
        for step in range(arc.count):
            wedge = (arc.first + step) % wedges
            held += 1
            if (arc.shell, arc.band, wedge) in skip:
                skipped += 1
                _close(units, ball, arc, run)
                run = []
                continue
            run.append(wedge)
            if len(run) == block_sectors:
                _close(units, ball, arc, run)
                run = []
        _close(units, ball, arc, run)
    units.sort(key=lambda unit: (unit.nearest_pc, unit.shell, unit.band, unit.first))
    return Plan(tuple(float(v) for v in center_pc), float(radius_pc), int(grow_shells), reach, float(edge_pc), int(block_sectors), tuple(units),
                held - skipped, skipped)


def _close(units, ball, arc, run):
    if not run:
        return
    nearest = min(detmath.dist(ball._position(arc.shell, arc.band, wedge), ball.center_pc) for wedge in run)
    units.append(Unit(arc.shell, arc.band, run[0], len(run), nearest))


def unit_items(unit, wedge_count):
    """A unit's sectors as the flat `shell, band, wedge` triples a sector stage takes (`run_galaxy.expand_generate_sectors`), from its first wedge on."""
    wedges = wedge_count(unit.shell, unit.band)
    items = []
    for step in range(unit.count):
        items.extend((unit.shell, unit.band, (unit.first + step) % wedges))
    return items


def next_units(region_plan, done_generate=(), done_topup=(), topup=False):
    """
    The barrier of design section 8.5: `(stage, units)`, the units still to run in the stage the run is in.

    Stages 0 to 2 (`STAGE_GENERATE`) are all the plan's units, in the plan's order. The top-up (`STAGE_TOPUP`, the same cut, only when `topup`)
    begins only when every generate unit is done, so a run stopped at any point resumes in the generate stage until the whole set is charted,
    and no top-up unit ever reads a half-built region. When nothing is left the stage is `STAGE_DONE`.

    Args:
        region_plan (Plan): The run.
        done_generate (container): The `key`s of the units whose stages 0 to 2 are charted.
        done_topup (container): The `key`s of the units topped up.
        topup (bool): Whether the run has a stage 3.
    """
    pending = [unit for unit in region_plan.units if key(unit) not in done_generate]
    if pending:
        return STAGE_GENERATE, pending
    if topup:
        pending = [unit for unit in region_plan.units if key(unit) not in done_topup]
        if pending:
            return STAGE_TOPUP, pending
    return STAGE_DONE, []


def cost(region_plan, seconds_per_sector=None):
    """
    What the run costs, as a dict: `sectors`, `units`, `largest_unit`, `shells` (shells the units touch), `volume_sectors` (the sphere's volume over
    one sector's, the design's scale: roughly 8,200 sectors at 50 pc and a 4 pc edge before the allowance), and, only when `seconds_per_sector` is
    given (the PERF items own that number), `seconds` run in one lane. A grown shell's share is `shell_sectors` of the last `grow_shells` shells.
    """
    shells = sorted({unit.shell for unit in region_plan.units})
    out = {
        "sectors": region_plan.sectors, "units": len(region_plan.units),
        "largest_unit": max((unit.count for unit in region_plan.units), default=0),
        "shells": len(shells),
        "volume_sectors": 4.0 / 3.0 * math.pi * (region_plan.radius_pc + region_plan.grow_shells * region_plan.edge_pc) ** 3 / region_plan.edge_pc ** 3,
        "skipped": region_plan.skipped,
    }
    if seconds_per_sector is not None:
        out["seconds"] = region_plan.sectors * seconds_per_sector
    return out
