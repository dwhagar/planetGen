# tests/test_star_scatter_passes.py

"""
GEN.185: the star scatter's passes. Pass 2 places every star born at or
above the mass limit, pass 3 marks the sectors holding one at least as
bright as the luminosity floor, pass 4 places the lighter stars at least
that bright and skips the marked sectors; a sector's own draw is lighter
and dimmer than both.
"""

import random

import pytest

from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import run_plan
from planetgen.generation import star_population as sp
from planetgen.generation.config import SystemConfig
from planetgen.physics import constants
from planetgen.physics.stellar_evolution import sample_living_star
from planetgen.galaxy.geometry import sector_position_pc
from tests.test_bright_star_scatter import (
    EDGE_PC, E_VALUE, EXTENTS, SHAPE, THRESHOLD, _plan_args, _row_dict, _seed_galaxy,
)

MASS_LIMIT = 3.0
LIGHT_FLOOR = 500.0
"""The luminosity the tests that need stars lighter than the limit use: at the real floor (2,500 and up, GEN.184) a
star under 3 solar masses is rare, so the toy galaxy would hold none."""


@pytest.mark.parametrize("population", [None, "young", "intermediate"])
def test_the_mass_ranges_split_the_luminosity_fraction_exactly(population):
    whole = sp.bright_star_fraction(THRESHOLD, population)
    light = sp.bright_star_fraction(THRESHOLD, population, (None, MASS_LIMIT))
    heavy = sp.bright_star_fraction(THRESHOLD, population, (MASS_LIMIT, None))
    assert light + heavy == pytest.approx(whole, rel=1e-9)


def test_the_placed_share_is_the_mass_pass_plus_the_lighter_bright_stars():
    for population in sp.POPULATIONS:
        massive = sp.massive_star_fraction(8.0, population)
        placed = sp.placed_star_fraction(3000.0, 8.0, population)
        assert placed == pytest.approx(massive + sp.bright_star_fraction(3000.0, population, (None, 8.0)))
        # Stars of 8 to 9 solar masses are dimmer than 3000 Lsun on the main sequence, so the mass pass adds some.
        assert placed >= sp.bright_star_fraction(3000.0, population) - 1e-18
    assert sp.placed_star_fraction(3000.0, None, "young") == sp.bright_star_fraction(3000.0, "young")


def test_an_empty_mass_range_is_refused():
    with pytest.raises(ValueError):
        sp.bright_star_fraction(THRESHOLD, "young", (8.0, 8.0))


def test_sampled_stars_stay_inside_their_mass_range():
    rng = random.Random(3)
    heavy = sp.sample_bright_stars(300, sp.MASS_PASS_MIN_LUMINOSITY_SOL, "young", rng, mass_range=(MASS_LIMIT, None))
    assert all(star["initial_mass_sol"] >= MASS_LIMIT * 0.999999 for star in heavy)
    light = sp.sample_bright_stars(300, LIGHT_FLOOR, "intermediate", rng, mass_range=(None, MASS_LIMIT))
    assert all(star["initial_mass_sol"] < MASS_LIMIT for star in light)
    assert all(star["luminosity_w"] >= LIGHT_FLOOR * constants.SOLAR_LUMINOSITY * 0.999999 for star in light)


def test_a_sectors_own_stars_are_lighter_than_the_mass_limit():
    rng = random.Random(4)
    for _ in range(300):
        mass, _age, _state = sample_living_star(rng=rng, max_mass_sol=MASS_LIMIT)
        assert mass < MASS_LIMIT


def test_fill_context_draws_lighter_stars_and_counts_the_mass_pass_in_its_share():
    center = sector_position_pc(2, 0, 3, EDGE_PC)
    plain = brightStars.FillContext(center, SHAPE, min_luminosity_sol=THRESHOLD)
    limited = brightStars.FillContext(center, SHAPE, min_luminosity_sol=THRESHOLD, star_mass_limit_sol=MASS_LIMIT)
    assert limited.bright_share() >= plain.bright_share() > 0.0
    cfg = limited.apply(SystemConfig(), random.Random(1))
    assert cfg.MAX_STAR_MASS_SOL == MASS_LIMIT and cfg.MAX_STAR_LUMINOSITY_SOL == THRESHOLD
    assert plain.apply(SystemConfig(), random.Random(1)).MAX_STAR_MASS_SOL is None


def test_the_backfill_draws_only_stars_lighter_than_the_mass_limit():
    addresses = [(2, 0, slot) for slot in range(15)] + [(3, 1, slot) for slot in range(21)]
    rows = list(brightStars.backfill_cells(SHAPE, addresses, EDGE_PC, E_VALUE, 100.0, THRESHOLD, 2,
                                           mass_range=(None, MASS_LIMIT)))
    assert rows and all(_row_dict(row)["initial_mass_sol"] < MASS_LIMIT for row in rows)


def _scattered(mysql_config, mass_limit=MASS_LIMIT):
    _seed_galaxy(mysql_config)
    args = _plan_args(mysql_config, "--bright-stars-only", "--workers", "1")
    args.phenomenon_min_mass = mass_limit
    args.bright_star_min_luminosity = LIGHT_FLOOR
    summary = run_plan.scatter_bright_stars(args)
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT * FROM bright_stars").fetchall()
        limit = store.bright_star_mass_limit(conn)
    finally:
        conn.close()
    return summary, rows, limit


def test_the_scatter_runs_the_mass_pass_then_a_luminosity_pass_that_skips_marked_sectors(mysql_config):
    summary, rows, limit = _scattered(mysql_config)
    assert limit == MASS_LIMIT and summary["total"] == len(rows) > 20
    floor_w = LIGHT_FLOOR * constants.SOLAR_LUMINOSITY
    heavy = [row for row in rows if row["initial_mass_sol"] >= MASS_LIMIT]
    light = [row for row in rows if row["initial_mass_sol"] < MASS_LIMIT]
    assert heavy and light
    assert all(row["luminosity_w"] >= floor_w * 0.999999 for row in light)
    marked = {(row["ring_index"], row["layer_index"], row["ring_slot_index"])
              for row in heavy if row["luminosity_w"] >= floor_w}
    assert marked
    # Pass 4 left every marked sector alone: it already holds a star that bright.
    assert not [row for row in light if (row["ring_index"], row["layer_index"], row["ring_slot_index"]) in marked]


def test_a_band_below_the_floor_draws_only_lighter_stars(mysql_config, monkeypatch):
    _scattered(mysql_config)
    args = _plan_args(mysql_config, "--bright-stars-down-to", "200", "--workers", "1")
    run_plan.add_bright_star_band(args)
    conn = store.get_connection(mysql_config)
    try:
        band = conn.execute("SELECT initial_mass_sol FROM bright_stars WHERE luminosity_w < ?",
                            (LIGHT_FLOOR * constants.SOLAR_LUMINOSITY,)).fetchall()
    finally:
        conn.close()
    assert band and all(row["initial_mass_sol"] < MASS_LIMIT for row in band)
