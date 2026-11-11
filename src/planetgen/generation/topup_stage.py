# generation/topup_stage.py

"""
Stage 3 of a sector's generation, run on the stored sector (GEN.1120; docs/design/systems/region-object-store.md section 8.5): the rogue top-up.

When a region run has finished stages 0 to 2 for every sector, each charted sector is compared with the mass it was calculated to hold (its expected systems
times the mean mass of a system there), the part of that mass nothing in the sector holds is its shortfall, and `physics/rogue_topup.plan` makes it up with
rogue planets, no more than the Rogues level allows. Each rogue is built (`rogue_topup_build`), saved with no donor and an apparent origin
(`store.save_topup_rogues`), and the cell's ledger row notes the top-up (`chart_ledger.mark_topped_up`), so a cell keeps one. What the level leaves unmatched
is returned for the light asteroid fields (GEN.1126).
"""

import collections
import math

from planetgen.db import chart_ledger, sector_mass, store
from planetgen.galaxy import dials, resident_mass
from planetgen.galaxy import seed as galaxy_seed
from planetgen.galaxy.profile import sector_profile
from planetgen.generation.config import SystemConfig
from planetgen.generation import field_topup_build, rogue_topup_build
from planetgen.physics import rogue_topup, topup_field
from planetgen.util import draw

MAX_HOSTS = 8
"""int: The most of a sector's systems loaded as donors of the top-up's mass. STAND-IN (GEN.725): a rogue's donor only supplies the mass, the snow line and the
age of the star; a handful is as good as all of them and a sector can hold hundreds."""

TopUpResult = collections.namedtuple("TopUpResult", ("written", "remainder_kg", "shortfall_kg", "stopped_by_level", "skipped", "fields", "unmatched_kg"),
                                     defaults=(0, 0.0))
"""What a sector's top-up did: the rogues written, the mass the level left over for the fields, the shortfall it started from, whether the level ended the
drawing, the reason nothing was done (`None` when it ran), the light asteroid fields written for the left-over mass and the mass under the fields' cut."""


AUTO_LEVEL = "few"
"""str: What the Auto level of the Rogues dial means (GEN.1506; Boss, 5:26 PM EST on 2026-11-09, after the PERF.136 measure: at Normal the rogues add about 67 percent to
generating a sector's systems, at Few about 43 percent). The word the Generate page, the console help and the dial's own text give for Auto (`galaxy/dials.py`)."""


def level_word(args_or_value):
    """The tuning key of the Rogues level (`Few`) from the lever's word (`few`), or from an `args` namespace; Auto (`None`) is `AUTO_LEVEL`."""
    value = getattr(args_or_value, "rogues", args_or_value)
    return AUTO_LEVEL.capitalize() if value is None else str(value).capitalize()


def galaxy_level(conn):
    """The Rogues level the stored plan was told (`galaxy_shape_pins`), or `AUTO_LEVEL` when it was left on Auto."""
    return dials.choice_pinned(store.get_galaxy_shape_pins(conn), "rogues") or AUTO_LEVEL


def topup_sector(conn, address, profile, level, system_config, edge_pc, stream_key=None):
    """
    Runs the top-up of one charted sector.

    Args:
        conn (Connection): Not mid-`batched`.
        address (tuple): The `(shell, band, wedge)` of the sector.
        profile (SectorProfile): The sector's place, for its expected mass.
        level (str): The Rogues level, `none`, `few`, `normal` or `many` (any case).
        system_config (SystemConfig): The shared settings.
        edge_pc (float): The grid's sector edge, pc.
        stream_key (str, optional): The key the draws are seeded from; the sector's address by default, so the same state gives the same rogues.

    Returns:
        TopUpResult: Nothing is written (and `skipped` says why) for a cell that is not charted, one already topped up, or a level of none.
    """
    word = level_word(level)
    ledger = chart_ledger.row(conn, address)
    if ledger is None:
        return TopUpResult(0, 0.0, 0.0, False, "not charted")
    if ledger["topped_up_at"] is not None:
        return TopUpResult(0, 0.0, 0.0, False, "already topped up")
    sector_id = store.kept_sector_id_at(conn, *address)
    held = sector_mass.sector_mass_kg(conn, sector_id)
    if held is None:
        return TopUpResult(0, 0.0, 0.0, False, "a body has no stored mass")
    short = resident_mass.shortfall_kg(resident_mass.remaining_budget(resident_mass.expected_mass_kg(profile)), held)
    stars = conn.execute("SELECT COUNT(*) AS n FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id WHERE ss.sector_id = ?",
                         (sector_id,)).fetchone()["n"]
    system_ids = [r["id"] for r in conn.execute("SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id LIMIT ?", (sector_id, MAX_HOSTS)).fetchall()]
    hosts = [store.load_star_system(conn, system_id).star for system_id in system_ids]
    key = stream_key or "topup:" + ".".join(str(part) for part in address)
    radius_pc = math.sqrt(sum(c * c for c in profile.position_pc))
    plan = rogue_topup.plan(
        draw.Stream(key + ":plan"), word, stars, short, [host.mass for host in hosts], rogues_present=store.count_rogues_in_sector(conn, sector_id),
        sector_edge_pc=edge_pc, galactocentric_radius_pc=radius_pc)
    made = rogue_topup_build.build(plan.rogues, hosts, system_config, key, address, edge_pc)
    written = store.save_topup_rogues(conn, sector_id, made)
    # What the level left unmatched becomes the lighter fields (Boss, 5:41 PM EDT), placed and given an apparent history the same way.
    light = topup_field.plan(draw.Stream(key + ":fields"), plan.remainder_kg, [host.mass for host in hosts] if plan.remainder_kg else [], sector_edge_pc=edge_pc,
                             galactocentric_radius_pc=radius_pc)
    field_count = store.save_topup_fields(conn, sector_id, field_topup_build.build(light.fields, system_config, key, address, edge_pc))
    chart_ledger.mark_topped_up(conn, address)
    sector_mass.refresh_sectors(conn, [tuple(address)])
    return TopUpResult(written, plan.remainder_kg, short, plan.stopped_by_level, None, field_count, light.unmatched_kg)


def topup_region(config, addresses, edge_pc, level=None, system_config=None, report=None):
    """
    Stage 3 over the charted sectors of a region run, one sector after the next (the barrier of `region_run.next_units`: every sector has finished stages 0
    to 2 before the first top-up reads the stored state).

    Args:
        config (MySQLConfig): The database.
        addresses (iterable): The `(shell, band, wedge)` of each sector of the region.
        edge_pc (float): The grid's sector edge, pc.
        level (str, optional): The Rogues level of the run (`--rogues`); the galaxy's own (`galaxy_level`) when omitted.
        system_config (SystemConfig, optional): The shared settings for the bodies; the defaults when omitted.
        report (callable, optional): Called with `(address, TopUpResult)` for each sector.

    Returns:
        dict: `{"sectors": n topped up, "rogues": n written, "fields": n light fields written, "remainder_kg": mass under the fields' cut}`.
    """
    conn = store.get_connection(config)
    totals = {"sectors": 0, "rogues": 0, "fields": 0, "remainder_kg": 0.0}
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            return totals
        word = level or galaxy_level(conn)
        for address in addresses:
            address = tuple(address)
            profile = sector_profile(address, skeleton)
            with galaxy_seed.seeded(skeleton.galaxy_seed, "sector", address, disk_scale_length_pc=skeleton.shape.disk_scale_length_pc, shape=skeleton.shape):
                result = topup_sector(conn, address, profile, word, system_config or SystemConfig(), edge_pc)
            conn.commit()
            if result.skipped is None:
                totals["sectors"] += 1
                totals["rogues"] += result.written
                totals["fields"] += result.fields
                totals["remainder_kg"] += result.unmatched_kg
            if report is not None:
                report(address, result)
    finally:
        conn.close()
    return totals
