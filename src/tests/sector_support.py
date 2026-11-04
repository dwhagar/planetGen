"""Makes the sector of one grid address for a test (PERF.149 rule 2).

`run_galaxy.ensure_sector_generated`, the "generate on a visit" entry point, was removed: no page or API route generates a
sector while it answers (an admin makes sectors by a run, `regenerate_sector_cell` or the Generate page). Many tests still
need one sector at a named address with the whole save path behind it (profile density, block totals, the bright-star
backfill, the settle), so `make_sector` is that path for a test, built from the same pieces the runs use."""

import logging

from planetgen.generation import run_galaxy, run_sector
from planetgen.db import store
from planetgen.queue import redisqueue

log = logging.getLogger(__name__)


def settle_sectors(config, sector_ids):
    """
    Saves the sector paths of `sector_ids` and the sectors around them (GEN.126), here and now, as a run of `settle-sector-paths` through the
    wrapper (`run_galaxy.settle_named_sectors`). With no Redis server for that run the paths are not saved and the next settle covers them.
    """
    try:
        run_galaxy.settle_named_sectors(sector_ids, config or store.DEFAULT_MYSQL_CONFIG)
    except redisqueue.Unavailable as exc:
        log.error(f"The sector paths of {len(sector_ids)} sector(s) were not saved: {exc}")


def backfill_bright_stars(config, center_pc):
    """The bright-star backfill (GEN.23) around one center: `run_galaxy.backfill_bright_stars_around` for a list of one."""
    return run_galaxy.backfill_bright_stars_around(config, [center_pc])


def make_sector(shell_index, band_index, wedge_index, config=None, backfill=True, settle=True, link=True, num_systems=None):
    """
    Generates and saves the sector at an address from its profile, as a run would, and returns what was made.

    Args:
        backfill (bool): Run the bright-star backfill around the new sector.
        settle (bool): Settle the new sector's paths (`settle_sectors` above).
        link (bool): Link the sector to its neighbours as it is saved.
        num_systems (int): Generate this many systems instead of the profile's count (TEST.403: a test that checks what is
            done around the sector, not what is in it, asks for one; a dense cell draws thousands).

    Returns:
        dict: `created`, `qualifies`, `in_reach`, `sector_id`, `sector_name` and, when created, `summary`. A cell that
            already holds a sector, or lies past the reach of the grid, makes nothing and says so.

    Raises:
        RuntimeError: The galaxy's skeleton has never been built.
    """
    address = (shell_index, band_index, wedge_index)
    early, args, position_pc, qualifies, edge_pc = run_galaxy._cell_to_generate(address, config, link)
    if early is not None:
        return early
    if num_systems is not None:
        args.num_systems = num_systems
        args.density = None   # a density resolves to the profile's own count and overrides this one
    sector_id, sector_name, sector = run_galaxy.generate_and_save_sector_at(args, address, position_pc, edge_pc)
    if backfill:
        backfill_bright_stars(config, position_pc)
    run_galaxy._refresh_block_totals(config, address)
    if settle:
        settle_sectors(config, [sector_id])
    return {"created": True, "qualifies": qualifies, "in_reach": True, "sector_id": sector_id, "sector_name": sector_name,
            "summary": run_sector.sector_generation_summary_lines(sector, args)}
