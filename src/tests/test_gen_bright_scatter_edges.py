# tests/test_gen_bright_scatter_edges.py

"""
TEST.24 Bright-star scatter edge cases and TEST.25 Interrupted bright-star
scatter (docs/TODO.md).

TEST.24: a shell whose stars all fail to place (a point that always rounds
out of its cell), a sparse shell far outside the disk (which still has its
chance since GEN.78), the core shell, empty extents, and a threshold below
every white dwarf.

TEST.25: a scatter worker that fails after some of its 10,000-row commits
(the threshold and seed are written only at the end) leaves a partial
`bright_stars` table; a re-plan replaces it with exactly one scatter's
stars, and a later fill never builds a leftover star into a system (nor
does the backfill after a run, GEN.30, leave any beside its own stars).
"""

import itertools
import math
import random

import pytest

from planetgen.generation import run_galaxy
from planetgen.generation import run_plan
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.galaxy import seed as galaxySeed
from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import sector_position_pc
from planetgen import tuning

from tests import first_fill, worker_patches
from tests.sector_support import backfill_bright_stars
from tests.test_bright_star_scatter import (
    light_steps,  # noqa: F401  (autouse: quick backfill steps)
    E_VALUE, EDGE_PC, EXTENTS, SHAPE, THRESHOLD, _scatter_args, _seed_galaxy,
)

EMPTY_SHELL = 300
"""A shell of the toy galaxy far outside its disk, where only the stellar halo
is left (GEN.234, which replaced GEN.78's flat floor)."""

HEAVIEST_SHELL = 8
"""The toy outline's shell with the most expected work, so the first task of a plan scatter (a task is a run of
consecutive shells, heaviest first, `run_plan.static_shell_chunks`; each shell of the toy outline is its own task but the
core and shell 1)."""

FAILING_SHELL = 6
"""The third heaviest: the two heavier shells are drawn before it."""


class _CountingRandom(random.Random):
    """A seeded `random.Random` that counts its `random()` draws."""

    def __init__(self, seed):
        super().__init__(seed)
        self.draws = 0

    def random(self):
        self.draws += 1
        return super().random()


# --- TEST.24: shells and extents ---------------------------------------------

def test_a_shell_whose_stars_all_fail_to_place_yields_nothing(monkeypatch):
    monkeypatch.setattr(brightStars, "_stored_in_cell", lambda *args, **kwargs: False)
    reports = []
    rows = list(brightStars.scatter_shell(SHAPE, 4, 3, EDGE_PC, E_VALUE, THRESHOLD, 9,
                                          on_progress=lambda done, estimate: reports.append((done, estimate))))
    assert rows == []
    # Every candidate was tried, and nothing was "done".
    assert reports and all(done == 0 for done, _estimate in reports)
    assert reports[-1][1] == pytest.approx(0.0, abs=1e-6)


def test_a_shell_far_outside_the_disk_keeps_a_small_chance_of_a_bright_star():
    # GEN.78: no cell is zeroed for being sparse; out there the stars are the halo's.
    top = spherical.band_range(EMPTY_SHELL)[1]
    densities = []
    for band in range(top - 3, top + 1):
        plane, z = brightStars.band_center_pc(EMPTY_SHELL, band, EDGE_PC)
        densities.append(brightStars._densities_at([(plane, 0.0)], z, SHAPE)[0])
    assert all(item["halo"] > 0.0 for item in densities)
    assert all(sum(item.values()) == pytest.approx(item["halo"], rel=1e-3) for item in densities)
    fractions = brightStars.band_fractions(THRESHOLD)
    expected = brightStars.shell_expected_stars(SHAPE, EMPTY_SHELL, top, EDGE_PC, E_VALUE, fractions)
    assert expected > 0.0
    # A scatter whose outline is only that shell finishes it.
    done = []
    list(brightStars.scatter(SHAPE, [(EMPTY_SHELL, top)], EDGE_PC, E_VALUE, THRESHOLD, 3,
                             on_shell=lambda *report: done.append(report)))
    assert done == [(1, 1)]


def test_the_core_shell_places_exactly_the_core_stars_of_a_wider_outline():
    only = [row for seed in range(21, 101) for row in brightStars.scatter_shell(SHAPE, 0, 0, EDGE_PC, E_VALUE,
                                                                                  THRESHOLD, seed)]
    assert only, "the bulge's core should hold a bright star at some of these seeds"
    assert {row[:3] for row in only} == {(0, 0, 0)}

    fractions = brightStars.band_fractions(THRESHOLD)
    weight, expected = brightStars.shell_weight(SHAPE, 0, 0, EDGE_PC, E_VALUE, fractions)
    assert expected > 0.0
    assert weight == pytest.approx(expected + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS)


def test_empty_extents_scatter_nothing():
    done = []
    assert list(brightStars.scatter(SHAPE, [], EDGE_PC, E_VALUE, THRESHOLD, 1,
                                    on_shell=lambda *report: done.append(report))) == []
    assert done == []


def test_a_plan_with_no_shells_records_an_empty_scatter(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=0,
                          expected_system_count_at_density_1=E_VALUE, config=mysql_config, cell_star_ceiling=1e12)
    store.replace_galaxy_shells([], config=mysql_config)
    summary = run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    assert summary["total"] == 0
    assert set(summary["counts"]) == set(brightStars.POPULATIONS)
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
        assert store.bright_star_scatter_settings(conn)[0] == THRESHOLD
    finally:
        conn.close()


# --- TEST.24: a threshold below every white dwarf ----------------------------

BELOW_WHITE_DWARFS = (tuning.WD_LUMINOSITY_RANGE_SOL[1] * 0.5,
                      tuning.WD_LUMINOSITY_RANGE_SOL[0] * 0.5)


@pytest.mark.parametrize("threshold", BELOW_WHITE_DWARFS)
def test_a_threshold_below_every_white_dwarf_is_refused_before_any_star(threshold):
    with pytest.raises(ValueError, match="white dwarf"):
        next(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, threshold, 1))
    with pytest.raises(ValueError, match="white dwarf"):
        brightStars.band_fractions(threshold)
    with pytest.raises(ValueError, match="white dwarf"):
        next(brightStars.backfill_cells(SHAPE, [(2, 0, 0)], EDGE_PC, E_VALUE, threshold, THRESHOLD, 1))


# --- TEST.25: an interrupted scatter -------------------------------------------

def _scatter_seed(mysql_config, address="scatter"):
    """The seed a scatter (or band, `band/<floor>-<ceiling>`) of this
    galaxy draws from: derived from the galaxy's seed (GEN.39)."""
    conn = store.get_connection(mysql_config)
    try:
        return galaxySeed.short_seed(store.get_galaxy_seed(conn), "bright-stars", address)
    finally:
        conn.close()


def _count(mysql_config, sql="SELECT COUNT(*) AS n FROM bright_stars", params=()):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(sql, params).fetchone()["n"]
    finally:
        conn.close()


def _settings(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.bright_star_scatter_settings(conn)
    finally:
        conn.close()


def _stored_rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute(f"SELECT {', '.join(store.BRIGHT_STAR_COLUMNS)} FROM bright_stars").fetchall()
    finally:
        conn.close()
    return [tuple(row[column] for column in store.BRIGHT_STAR_COLUMNS) for row in rows]


def _comparable(rows):
    """Rows as compared after a round trip through MySQL (floats as stored)."""
    return sorted(tuple(float(value) if isinstance(value, float) else value for value in row) for row in rows)


def _stored_shape(mysql_config):
    """The galaxy's shape as the scatter reads it: from the database, with the masses its seed draws (GEN.469)."""
    conn = store.get_connection(mysql_config)
    try:
        return store.get_galaxy_shape(conn).shape
    finally:
        conn.close()


def _expected_scatter(mysql_config, seed):
    """What the scatter's passes give (GEN.185): the mass pass, then the lighter stars at least `THRESHOLD`
    bright, skipping the sectors the mass pass marked."""
    from planetgen.generation import star_population
    from planetgen.physics import constants
    limit = tuning.PHENOMENON_MIN_MASS_SOLAR
    mass_seed = _scatter_seed(mysql_config, "mass-scatter")
    shape = _stored_shape(mysql_config)
    heavy = list(brightStars.scatter(shape, EXTENTS, EDGE_PC, E_VALUE, star_population.MASS_PASS_MIN_LUMINOSITY_SOL,
                                     mass_seed, mass_range=(limit, None)))
    marked = {row[:3] for row in heavy if row[12] >= THRESHOLD * constants.SOLAR_LUMINOSITY}
    light = list(brightStars.scatter(shape, EXTENTS, EDGE_PC, E_VALUE, THRESHOLD, seed, skip_addresses=marked,
                                     mass_range=(None, limit)))
    return heavy + light


ROWS_BEFORE_FAILURE = 25_000
"""Rows the failing shell yields: two 10,000-row commits, then 5,000 lost."""


def _yielded_rows(scatter_shell):
    """What `heaviest_shell_fails_after_commits` yields before it fails, from
    the real `scatter_shell`."""
    template = list(scatter_shell(SHAPE, HEAVIEST_SHELL, 3, EDGE_PC, E_VALUE, THRESHOLD, 5))
    assert template
    return list(itertools.islice(itertools.cycle(template), ROWS_BEFORE_FAILURE))


def heaviest_shell_fails_after_commits(real):
    """A `brightStars.scatter_shell` whose heaviest shell yields
    `ROWS_BEFORE_FAILURE` real rows over and over, then fails
    (`worker_patches.patch_everywhere` builds it in the workers too)."""

    def failing_shell(shape, shell_index, *args, **kwargs):
        if shell_index != HEAVIEST_SHELL:
            yield from real(shape, shell_index, *args, **kwargs)
            return
        yield from _yielded_rows(real)
        raise RuntimeError("scatter worker died")

    return failing_shell


def shell_fails(real, shell):
    """A `brightStars.scatter_shell` that fails on `shell` before drawing."""

    def failing_shell(shape, shell_index, *args, **kwargs):
        if shell_index == shell:
            raise RuntimeError("scatter worker died")
        yield from real(shape, shell_index, *args, **kwargs)

    return failing_shell


def _interrupt_after_commits(mysql_config, monkeypatch):
    """Runs a plan scatter whose first task (the heaviest shell) yields
    `ROWS_BEFORE_FAILURE` real rows over and over, then fails. Returns the
    rows it yielded."""
    yielded = _yielded_rows(brightStars.scatter_shell)
    undo = worker_patches.patch_everywhere(monkeypatch, brightStars, "scatter_shell",
                                           "tests.test_gen_bright_scatter_edges:heaviest_shell_fails_after_commits")
    with pytest.raises(RuntimeError, match="scatter worker died"):
        run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    undo()
    return yielded


def _shell_rows(rows, shell_index):
    return [row for row in rows if row[0] == shell_index]


def test_a_scatter_that_fails_after_two_commits_keeps_them_and_no_seed(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    yielded = _interrupt_after_commits(mysql_config, monkeypatch)
    stored = _stored_rows(mysql_config)
    # The heaviest shell keeps its two commits and loses the rest. With one
    # worker the run stops there; with more, the shells other workers were
    # drawing at the time finish and keep their stars.
    assert _comparable(_shell_rows(stored, HEAVIEST_SHELL)) == _comparable(yielded[:20_000])
    if worker_patches.workers() == 1:
        assert len(stored) == 20_000
    # The threshold and seed are written only once every shell is in.
    assert _settings(mysql_config) is None


def test_a_re_plan_after_an_interrupted_scatter_holds_exactly_one_scatter(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    _interrupt_after_commits(mysql_config, monkeypatch)
    summary = run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))

    seed = _scatter_seed(mysql_config)
    expected = _expected_scatter(mysql_config, seed)
    stored = _stored_rows(mysql_config)
    assert summary["total"] == len(expected) == len(stored)
    assert _comparable(stored) == _comparable(expected)
    # Every star once: no two rows share a position.
    assert len({row[3:6] for row in stored}) == len(stored)
    assert _settings(mysql_config) == (THRESHOLD, seed)


def _leftover_cells(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT shell_index, band_index, wedge_index, COUNT(*) AS n, MAX(id) AS top"
                            " FROM bright_stars GROUP BY shell_index, band_index, wedge_index").fetchall()
    finally:
        conn.close()
    return {(row["shell_index"], row["band_index"], row["wedge_index"]): row["n"] for row in rows}


def _reached(mysql_config):
    """The sectors a backfill took to their own level or mass (GEN.44, GEN.187)."""
    conn = store.get_connection(mysql_config)
    try:
        masses = conn.execute("SELECT shell_index, band_index, wedge_index FROM sector_stats"
                              " WHERE bright_mass_sol IS NOT NULL").fetchall()
        return set(store.sector_bright_level_keys(conn)) | {
            (row["shell_index"], row["band_index"], row["wedge_index"]) for row in masses}
    finally:
        conn.close()


def _backfill_steps(mysql_config, address):
    """The sectors a backfill around `address` reaches (GEN.187), with their masses."""
    conn = store.get_connection(mysql_config)
    try:
        return run_galaxy.backfill_step_targets({address}, store.get_galaxy_bounds(conn))
    finally:
        conn.close()


def _nearest(cells, address):
    return min(cells, key=lambda cell: (math.dist(sector_position_pc(*cell, EDGE_PC),
                                                  sector_position_pc(*address, EDGE_PC)), cell))


def _is_built_star(entry):
    """Whether a sector's entry is the system built around a stored bright star (the members of an emission region are
    preplaced too, GEN.321, and carry no row id)."""
    return getattr(entry, "bright_star_id", None) is not None


def _fill(mysql_config, address):
    """Fills one sector, as a galaxy run does for each of its sectors
    (no backfill: since GEN.30 that runs once after the run's sectors)."""
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    return run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)


def _backfill_around(mysql_config, address):
    """The backfill a `galaxy --wedge` run of `address` ends with
    (`backfill_after_run`, GEN.30): around that sector, the default mass steps."""
    return backfill_bright_stars(mysql_config, sector_position_pc(*address, EDGE_PC))


def _unbuilt_leftover_cells(mysql_config, last_leftover):
    conn = store.get_connection(mysql_config)
    try:
        return {
            (row["shell_index"], row["band_index"], row["wedge_index"])
            for row in conn.execute("SELECT shell_index, band_index, wedge_index FROM bright_stars"
                                    " WHERE id <= ? AND star_system_id IS NULL", (last_leftover,)).fetchall()
        }
    finally:
        conn.close()


def _cells_with_stars_after(mysql_config, last_id):
    """The cells holding a bright star with an id above `last_id` (stars drawn after that point)."""
    conn = store.get_connection(mysql_config)
    try:
        return {
            (row["shell_index"], row["band_index"], row["wedge_index"])
            for row in conn.execute("SELECT DISTINCT shell_index, band_index, wedge_index FROM bright_stars"
                                    " WHERE id > ?", (last_id,)).fetchall()
        }
    finally:
        conn.close()


def _built_in(mysql_config, address):
    """`(ids of the stars built in this cell, its unbuilt star count)`."""
    conn = store.get_connection(mysql_config)
    try:
        built = [row["id"] for row in conn.execute(
            "SELECT b.id FROM bright_stars b JOIN star_systems s ON s.id = b.star_system_id"
            " WHERE b.shell_index = ? AND b.band_index = ? AND b.wedge_index = ?", address).fetchall()]
        unlinked = conn.execute(
            "SELECT COUNT(*) AS n FROM bright_stars WHERE shell_index = ? AND band_index = ?"
            " AND wedge_index = ? AND star_system_id IS NULL", address).fetchone()["n"]
    finally:
        conn.close()
    return built, unlinked


def test_a_fill_after_a_shell_failed_mid_scatter_builds_no_leftover_star(mysql_config, monkeypatch):
    # The two heaviest shells are written and committed, the third fails: no
    # threshold is recorded. A run's fill then builds none of its cell's
    # leftovers (no level to fill down to), and the backfill after it
    # (GEN.30) draws its sectors from their step masses with no ceiling.
    # Every other cell it reaches must lose its leftovers then (GEN.44: no
    # level anywhere, yet stars), or a later fill there would build them
    # as well as the backfill's stars (every bright star twice).
    _seed_galaxy(mysql_config)
    first_fill.without_the_deepest_floor(monkeypatch)   # this is about the scatter's leftovers, not the first fill's draw
    undo = worker_patches.patch_everywhere(monkeypatch, brightStars, "scatter_shell",
                                           "tests.test_gen_bright_scatter_edges:shell_fails", shell=FAILING_SHELL)
    with pytest.raises(RuntimeError, match="scatter worker died"):
        run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    undo()
    assert _settings(mysql_config) is None

    leftovers = _leftover_cells(mysql_config)
    # With more than one worker, shells handed out beside the failing one
    # may finish too; the failing shell itself holds nothing.
    shells = {address[0] for address in leftovers}
    assert leftovers and {8, 7} <= shells and FAILING_SHELL not in shells
    if worker_patches.workers() == 1:
        assert shells == {8, 7}
    last_leftover = _count(mysql_config, "SELECT MAX(id) AS n FROM bright_stars")
    address = max((address for address in leftovers if address[1] == 0), key=lambda a: (leftovers[a], a))
    ring_cells = [cell for cell in _backfill_steps(mysql_config, address) if cell in leftovers]

    _sector_id, _name, sector = _fill(mysql_config, address)
    assert [entry for entry in sector.entries if _is_built_star(entry)] == []
    assert _count(mysql_config, "SELECT COUNT(*) AS n FROM bright_stars WHERE star_system_id IS NOT NULL") == 0
    assert _backfill_around(mysql_config, address)["stars"] > 0
    # The backfill draws each sector's stars from its mass, and a sector it reaches can draw none (about 1 run
    # in 20 for the nearest one, TEST.132): the cell filled next is the nearest that was given stars.
    given = _cells_with_stars_after(mysql_config, last_leftover)
    neighbor = _nearest([cell for cell in ring_cells if cell in given], address)
    reached = _reached(mysql_config)
    assert neighbor in reached
    # The filled cell's own leftovers stay unbuilt (a filled sector never
    # gets stars); the cells the backfill reached hold none.
    assert _unbuilt_leftover_cells(mysql_config, last_leftover) & (reached | {address}) == {address}

    _sector_id, _name, sector = _fill(mysql_config, neighbor)
    preplaced = [entry for entry in sector.entries if _is_built_star(entry)]
    built, unlinked = _built_in(mysql_config, neighbor)
    assert unlinked == 0
    assert preplaced and len(built) == len(preplaced)
    # Every star built here came from the backfill, none from the
    # unfinished scatter.
    assert all(star_id > last_leftover for star_id in built)


def test_a_fill_after_failed_commits_clears_the_reached_sectors_leftovers(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    first_fill.without_the_deepest_floor(monkeypatch)   # this is about the scatter's leftovers, not the first fill's draw
    _interrupt_after_commits(mysql_config, monkeypatch)
    leftovers = _leftover_cells(mysql_config)
    # The heaviest shell's equator cell with the fewest leftovers, so the fill stays small
    # even if it did build them.
    address = min((cell for cell in leftovers if cell[1] == 0), key=lambda cell: (leftovers[cell], cell))
    neighbor = _nearest([cell for cell in _backfill_steps(mysql_config, address) if cell in leftovers], address)
    last_leftover = _count(mysql_config, "SELECT MAX(id) AS n FROM bright_stars")
    # With more workers, the shells drawn beside the heaviest left theirs too.
    assert last_leftover == 20_000 if worker_patches.workers() == 1 else last_leftover >= 20_000
    _sector_id, _name, sector = _fill(mysql_config, address)
    assert [entry for entry in sector.entries if _is_built_star(entry)] == []
    _backfill_around(mysql_config, address)

    assert _count(mysql_config, "SELECT COUNT(*) AS n FROM bright_stars WHERE star_system_id IS NOT NULL") == 0
    after = _leftover_cells(mysql_config)
    backfilled = _reached(mysql_config) | {address}
    stale_cells = _unbuilt_leftover_cells(mysql_config, last_leftover)
    assert neighbor in backfilled
    assert stale_cells, "leftovers outside the backfilled sectors should still be there"
    # Only the filled cell keeps its leftovers among the backfilled
    # sectors: it never gets stars again, so they are never built.
    assert stale_cells & backfilled == {address}
    # Leftovers outside the backfilled sectors are untouched until a
    # backfill reaches them (or a re-plan clears them).
    assert {cell for cell in leftovers if cell not in backfilled} <= set(after)


BAND_FLOOR = 300.0


def test_re_running_an_interrupted_band_holds_the_band_once(mysql_config, monkeypatch):
    # GEN.32: the re-run starts the band over, so the shells the first
    # run finished don't hold it twice.
    _seed_galaxy(mysql_config)
    seed = _scatter_seed(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "20"))
    first = _count(mysql_config)
    first_rows = _stored_rows(mysql_config)
    undo = worker_patches.patch_everywhere(monkeypatch, brightStars, "scatter_shell",
                                           "tests.test_gen_bright_scatter_edges:shell_fails", shell=FAILING_SHELL)
    with pytest.raises(RuntimeError, match="scatter worker died"):
        run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", str(BAND_FLOOR)))
    undo()
    assert _settings(mysql_config) == (THRESHOLD, seed)
    assert _count(mysql_config) > first

    band = run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", str(BAND_FLOOR)))
    assert band is not None
    assert _settings(mysql_config) == (BAND_FLOOR, seed)
    band_seed = _scatter_seed(mysql_config, f"band/{BAND_FLOOR:g}-{THRESHOLD:g}")
    expected = list(brightStars.scatter(_stored_shape(mysql_config), EXTENTS, EDGE_PC, E_VALUE, BAND_FLOOR, band_seed,
                                        max_luminosity_sol=THRESHOLD,
                                        mass_range=(None, 20.0)))
    stored = _stored_rows(mysql_config)
    assert len(stored) == first + len(expected)
    # The first scatter's stars and the band's, each once.
    assert _comparable(stored) == _comparable(first_rows + expected)
