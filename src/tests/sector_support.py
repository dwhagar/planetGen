"""Makes the sector of one grid address for a test (PERF.149 rule 2).

`run_galaxy.ensure_sector_generated`, the "generate on a visit" entry point, was removed: no page or API route generates a
sector while it answers (an admin makes sectors by a run, `regenerate_sector_cell` or the Generate page). Many tests still
need one sector at a named address with the whole save path behind it (profile density, block totals, the bright-star
backfill, the queued settle), so `make_sector` is that path for a test, built from the same pieces the runs use."""

from planetgen.generation import run_galaxy, run_sector


def make_sector(shell_index, band_index, wedge_index, config=None, backfill=True, settle=True, link=True):
    """
    Generates and saves the sector at an address from its profile, as a run would, and returns what was made.

    Args:
        backfill (bool): Run the bright-star backfill around the new sector.
        settle (bool): Queue the settle job (`run_galaxy.queue_settle`) for the new sector.
        link (bool): Link the sector to its neighbours as it is saved.

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
    sector_id, sector_name, sector = run_galaxy.generate_and_save_sector_at(args, address, position_pc, edge_pc)
    if backfill:
        run_galaxy.backfill_bright_stars(config, position_pc)
    run_galaxy._refresh_block_totals(config, address)
    if settle:
        run_galaxy.queue_settle(config, [sector_id])
    return {"created": True, "qualifies": qualifies, "in_reach": True, "sector_id": sector_id, "sector_name": sector_name,
            "summary": run_sector.sector_generation_summary_lines(sector, args)}
