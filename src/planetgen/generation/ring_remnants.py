"""
The neutron stars and black holes below the phenomenon scatter's mass cut that the sectors of a backfill ring hold before they are filled (GEN.1682).

The mass backfill (GEN.187, `docs/design/generation/mass-backfill.md`) takes the sectors around a run's own to the stars born at 1, 2, 5 and 8 solar masses and up, but the phenomenon scatter places only the neutron
stars and black holes at or above its cut (`tuning.PHENOMENON_MIN_MASS_SOLAR`, 14 solar masses), and the share below the cut is drawn when a sector is filled (`phenomenon_scatter.below_cut_draws`, GEN.168). A ring
sector was therefore full of stars and held none of its neutron stars and 60 percent of its black holes until somebody filled it. `rows_for` draws them, with the fill's own function on the sector's own stream
(`phenomenon_scatter.ring_remnant_rows`), as scatter rows in the sector's cell: the fill builds them like any scattered row, and does not draw them again (`run_galaxy.fill_context` reads
`store.sector_holds_below_cut_remnants`). It runs where the scatter's rows are written and nothing clears them afterwards: the `write-special-phenomena` unit of a galaxy run (the clear of the scatter comes before), and the
backfill's own chunk when the galaxy already has its scatter (an on-demand visit or a run that adds sectors).
"""

from planetgen.db import store
from planetgen.galaxy import resolution as galaxyResolution, seed as galaxySeed
from planetgen.galaxy.profile import sector_profile
from planetgen.generation import phenomenon_scatter
from planetgen.generation.phenomenon_cuts import as_cuts


def rows_for(conn, skeleton, limits, seed, min_mass_solar, addresses, held=None, on_sector=None):
    """
    The scatter rows (`PHENOMENON_SCATTER_COLUMNS` order) of the below-cut neutron stars and black holes of the sectors `addresses`: per sector the draws of `phenomenon_scatter.ring_remnant_rows`, at the density the
    scatter used (`sector_profile(...).expected_systems`, as the fill's `BelowCut` has it). A sector the fill would give no remnant of its own is left out: a cell above the ceiling (`resolution.capped`: its remnants
    are dark and part of what it leaves unresolved), and one that already holds some (a run again, or a sector whose chunk was drawn before). `held` is that set of cells when the caller has it from one read
    (`store.addresses_holding_below_cut_remnants`, for the unit that works while the table has no index); without it each cell is asked of the table.

    Args:
        limits (ScatterLimits): The galaxy scatter's limits by place (`run_plan._scatter_limits_and_seed`): the mass limit a cell's resolution reads.
        seed (int): The phenomenon scatter's seed (`store.phenomenon_scatter_settings`).
        min_mass_solar (float, `phenomenon_cuts.Cuts` or None): The scatter's cut, by place for `Cuts` (GEN.377: each sector takes its own band's); `None` draws nothing.
        held (set or None): The cells that already hold below-cut remnants.
        on_sector (callable or None): Called once per address visited (a progress tick).

    Returns:
        list: The rows, in address order.
    """
    rows = []
    if min_mass_solar is None:
        return rows
    cuts = as_cuts(min_mass_solar)
    for address in addresses:
        address = tuple(address)
        if on_sector is not None:
            on_sector()
        cut = cuts.cut_at(address[0], address[1])      # GEN.377: the cut of the sector's own band
        if (address in held) if held is not None else store.sector_holds_below_cut_remnants(conn, address, cut):
            continue
        profile = sector_profile(address, skeleton)
        with galaxySeed.seeded(skeleton.galaxy_seed, "sector", address, disk_scale_length_pc=skeleton.shape.disk_scale_length_pc, shape=skeleton.shape):
            resolution = galaxyResolution.cell_resolution(profile, skeleton.cell_star_ceiling, skeleton.fill_luminosity_cap_sol, limits.mass_limit_at(address[0], address[1]))
        if resolution.capped:
            continue
        rows.extend(phenomenon_scatter.ring_remnant_rows(address, profile.position_pc, profile.shape, profile.expected_systems, cut, seed, skeleton.edge_pc))
    return rows
