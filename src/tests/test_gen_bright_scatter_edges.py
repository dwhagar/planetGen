# tests/test_gen_bright_scatter_edges.py

"""
TEST.24 Bright-star scatter edge cases and TEST.25 Interrupted bright-star
scatter (docs/TODO.md).

TEST.24: `brightStars._place_one` running out of redraws (a slot that
never qualifies, a point that always rounds out of its cell), a zero-weight
bin picked by float rounding, a layer where nothing qualifies,
`outer_ring=0`, empty extents, and a threshold below every white dwarf.

TEST.25: a scatter worker that fails after some of its 10,000-row commits
(the threshold and seed are written only at the end) leaves a partial
`bright_stars` table; a re-plan replaces it with exactly one scatter's
stars, and a later fill never builds a leftover star into a system (nor
does the backfill after a run, GEN.30, leave any beside its own stars).
"""

import itertools
import random

import pytest

import generate
from stellarObjects import _db, brightStars, galaxySeed
from stellarObjects.galaxyDrill import DrillBlock, drill_parent
from stellarObjects.galaxyGeometry import ring_sector_count, sector_position_pc
from stellarObjects import program_constants

from tests import worker_patches
from tests.test_bright_star_scatter import (
    E_VALUE, EDGE_PC, EXTENTS, SHAPE, THRESHOLD, _plan_args, _seed_galaxy,
)

EMPTY_LAYER = 40
"""A layer of the toy galaxy far above its disk: no slot qualifies."""


class _CountingRandom(random.Random):
    """A seeded `random.Random` that counts its `random()` draws."""

    def __init__(self, seed):
        super().__init__(seed)
        self.draws = 0

    def random(self):
        self.draws += 1
        return super().random()


class _ScriptedRandom(random.Random):
    """Returns `script`'s values first, then a seeded stream."""

    def __new__(cls, script, seed=1):
        # Before Python 3.11, Random.__new__ seeds itself with its first
        # argument, which can't be a list (TEST.76).
        return super().__new__(cls, seed)

    def __init__(self, script, seed=1):
        super().__init__(seed)
        self.script = list(script)

    def random(self):
        if self.script:
            return self.script.pop(0)
        return super().random()


# --- TEST.24: _place_one ----------------------------------------------------

def test_place_one_gives_up_after_its_redraws_when_no_slot_qualifies():
    ring_index = 3
    slots = ring_sector_count(ring_index)
    rng = _CountingRandom(4)
    spot = brightStars._place_one(rng, [1.0] * min(slots, brightStars.ANGLE_BINS), ring_index, EMPTY_LAYER,
                                  slots, SHAPE, E_VALUE, EDGE_PC)
    assert spot is None
    # Two draws (bin, angle) per try, and no point drawn in a cell.
    assert rng.draws == 2 * brightStars.SLOT_REDRAWS


def test_place_one_gives_up_when_every_point_rounds_out_of_its_cell(monkeypatch):
    calls = []

    def never_inside(*args):
        calls.append(args)
        return None

    monkeypatch.setattr(brightStars, "_point_in_cell", never_inside)
    slots = ring_sector_count(1)
    spot = brightStars._place_one(random.Random(2), [1.0] * slots, 1, 0, slots, SHAPE, E_VALUE, EDGE_PC)
    assert spot is None
    assert len(calls) == brightStars.SLOT_REDRAWS


def test_a_layer_whose_stars_all_fail_to_place_yields_nothing(monkeypatch):
    monkeypatch.setattr(brightStars, "_place_one", lambda *args, **kwargs: None)
    reports = []
    rows = list(brightStars.scatter_layer(SHAPE, 0, 3, EDGE_PC, E_VALUE, THRESHOLD, 9,
                                          on_progress=lambda done, estimate: reports.append((done, estimate))))
    assert rows == []
    # Every star still counted down its estimate, and nothing was "done".
    assert reports and all(done == 0 for done, _estimate in reports)
    assert reports[-1][1] == pytest.approx(0.0, abs=1e-6)


def test_a_zero_weight_bin_is_never_picked_by_float_rounding():
    # 0.1 + 0.2 + 0.3 rounds up, so the top of the range (`random()` just
    # under 1) is left over after every positive bin, and the old walk fell
    # through to the last bin, whose weight is zero. Ring 1's slots all
    # qualify, so only the weights keep a star out of bins 3..8.
    ring_index, layer_index = 1, 0
    slots = ring_sector_count(ring_index)
    weights = [0.1, 0.2, 0.3] + [0.0] * (slots - 3)
    assert all(generate.brightStars._qualifies(sector_position_pc(ring_index, layer_index, slot, EDGE_PC),
                                               SHAPE, E_VALUE) for slot in range(slots))
    rng = _ScriptedRandom([1.0 - 2.0 ** -53, 0.5])
    spot = brightStars._place_one(rng, weights, ring_index, layer_index, slots, SHAPE, E_VALUE, EDGE_PC)
    assert spot is not None
    assert spot[0] in (0, 1, 2)
    # A walk that ends on the last bin by rounding lands in the last bin
    # with weight, at its far edge.
    assert spot[0] == 2


def test_leading_and_trailing_zero_weight_bins_get_no_stars():
    ring_index, layer_index = 1, 0
    slots = ring_sector_count(ring_index)
    weights = [0.0, 0.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0]
    assert len(weights) == slots
    rng = random.Random(17)
    picked = set()
    for _ in range(300):
        spot = brightStars._place_one(rng, weights, ring_index, layer_index, slots, SHAPE, E_VALUE, EDGE_PC)
        assert spot is not None
        picked.add(spot[0])
    assert picked == {2, 3}


# --- TEST.24: layers and extents ---------------------------------------------

def test_a_layer_where_nothing_qualifies_draws_and_reports_nothing():
    slots, bins = brightStars._ring_bins(5, EMPTY_LAYER, SHAPE, E_VALUE, EDGE_PC)
    assert bins and all(densities is None for densities in bins)
    reports = []
    rows = list(brightStars.scatter_layer(SHAPE, EMPTY_LAYER, 8, EDGE_PC, E_VALUE, THRESHOLD, 3,
                                          on_progress=lambda *report: reports.append(report)))
    assert rows == [] and reports == []
    fractions = brightStars.band_fractions(THRESHOLD)
    assert brightStars.layer_expected_stars(SHAPE, EMPTY_LAYER, 8, EDGE_PC, E_VALUE, fractions) == 0.0
    # A scatter whose outline is only that layer still finishes it.
    done = []
    assert list(brightStars.scatter(SHAPE, [(EMPTY_LAYER, 8)], EDGE_PC, E_VALUE, THRESHOLD, 3,
                                    on_layer=lambda *report: done.append(report))) == []
    assert done == [(1, 1)]


def test_outer_ring_zero_places_exactly_the_ring_zero_stars_of_a_wider_layer():
    seed = 21
    only = list(brightStars.scatter_layer(SHAPE, 0, 0, EDGE_PC, E_VALUE, THRESHOLD, seed))
    wide = list(brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, seed))
    assert only, "the bulge's ring 0 should hold a bright star at this seed"
    assert {row[0] for row in only} == {0}
    assert {row[2] for row in only} <= set(range(ring_sector_count(0)))
    # Ring 0 is placed first from the layer's own stream, so it lands in
    # the same spots (the stars drawn into them come later in the stream).
    spots = lambda rows: sorted((row[6], row[:6]) for row in rows if row[0] == 0)  # noqa: E731
    assert spots(only) == spots(wide)

    fractions = brightStars.band_fractions(THRESHOLD)
    weight, expected = brightStars.layer_weight(SHAPE, 0, 0, EDGE_PC, E_VALUE, fractions)
    assert expected > 0.0
    assert weight == pytest.approx(expected + brightStars.RING_WEIGHT_STARS)
    assert generate._layer_slots(0) == ring_sector_count(0)


def test_empty_extents_scatter_nothing():
    done = []
    assert list(brightStars.scatter(SHAPE, [], EDGE_PC, E_VALUE, THRESHOLD, 1,
                                    on_layer=lambda *report: done.append(report))) == []
    assert done == []


def test_a_plan_with_no_layers_records_an_empty_scatter(mysql_config):
    _db.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=0,
                          expected_system_count_at_density_1=E_VALUE, config=mysql_config)
    _db.replace_galaxy_layers([], config=mysql_config)
    summary = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    assert summary["total"] == 0
    assert set(summary["counts"]) == set(brightStars.POPULATIONS)
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
        assert _db.bright_star_scatter_settings(conn)[0] == THRESHOLD
    finally:
        conn.close()


# --- TEST.24: a threshold below every white dwarf ----------------------------

BELOW_WHITE_DWARFS = (program_constants.WD_LUMINOSITY_RANGE_SOL[1] * 0.5,
                      program_constants.WD_LUMINOSITY_RANGE_SOL[0] * 0.5)


@pytest.mark.parametrize("threshold", BELOW_WHITE_DWARFS)
def test_a_threshold_below_every_white_dwarf_is_refused_before_any_star(threshold):
    with pytest.raises(ValueError, match="white dwarf"):
        next(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, threshold, 1))
    with pytest.raises(ValueError, match="white dwarf"):
        brightStars.band_fractions(threshold)
    with pytest.raises(ValueError, match="white dwarf"):
        next(brightStars.backfill_cells(SHAPE, [(2, 0, 0)], EDGE_PC, E_VALUE, threshold, THRESHOLD,
                                        random.Random(1)))


def test_a_backfill_below_every_white_dwarf_writes_nothing(mysql_config):
    _seed_galaxy(mysql_config)
    with pytest.raises(ValueError, match="white dwarf"):
        generate.backfill_bright_stars(mysql_config, sector_position_pc(4, 0, 5, EDGE_PC), radius_ly=20.0,
                                       min_luminosity_sol=BELOW_WHITE_DWARFS[0])
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
        # No block is left holding a lock row (level NULL) or a level.
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_star_blocks").fetchone()["n"] == 0
    finally:
        conn.close()


# --- TEST.25: an interrupted scatter -------------------------------------------

def _scatter_seed(mysql_config, address="scatter"):
    """The seed a scatter (or band, `band/<floor>-<ceiling>`) of this
    galaxy draws from: derived from the galaxy's seed (GEN.39)."""
    conn = _db.get_connection(mysql_config)
    try:
        return galaxySeed.short_seed(_db.get_galaxy_seed(conn), "bright-stars", address)
    finally:
        conn.close()


def _count(mysql_config, sql="SELECT COUNT(*) AS n FROM bright_stars", params=()):
    conn = _db.get_connection(mysql_config)
    try:
        return conn.execute(sql, params).fetchone()["n"]
    finally:
        conn.close()


def _settings(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        return _db.bright_star_scatter_settings(conn)
    finally:
        conn.close()


def _stored_rows(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute(f"SELECT {', '.join(_db.BRIGHT_STAR_COLUMNS)} FROM bright_stars").fetchall()
    finally:
        conn.close()
    return [tuple(row[column] for column in _db.BRIGHT_STAR_COLUMNS) for row in rows]


def _comparable(rows):
    """Rows as compared after a round trip through MySQL (floats as stored)."""
    return sorted(tuple(float(value) if isinstance(value, float) else value for value in row) for row in rows)


ROWS_BEFORE_FAILURE = 25_000
"""Rows the failing layer yields: two 10,000-row commits, then 5,000 lost."""


def _yielded_rows(scatter_layer):
    """What `layer_zero_fails_after_commits` yields before it fails, from
    the real `scatter_layer`."""
    template = list(scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, 5))
    assert template
    return list(itertools.islice(itertools.cycle(template), ROWS_BEFORE_FAILURE))


def layer_zero_fails_after_commits(real):
    """A `brightStars.scatter_layer` whose layer 0 yields
    `ROWS_BEFORE_FAILURE` real layer-0 rows over and over, then fails
    (`worker_patches.patch_everywhere` builds it in the workers too)."""

    def failing_layer(shape, layer_index, *args, **kwargs):
        if layer_index != 0:
            yield from real(shape, layer_index, *args, **kwargs)
            return
        yield from _yielded_rows(real)
        raise RuntimeError("scatter worker died")

    return failing_layer


def layer_fails(real, layer):
    """A `brightStars.scatter_layer` that fails on `layer` before drawing."""

    def failing_layer(shape, layer_index, *args, **kwargs):
        if layer_index == layer:
            raise RuntimeError("scatter worker died")
        yield from real(shape, layer_index, *args, **kwargs)

    return failing_layer


def _interrupt_after_commits(mysql_config, monkeypatch):
    """Runs a plan scatter whose first layer (layer 0, the densest, goes
    first) yields `ROWS_BEFORE_FAILURE` real layer-0 rows over and over,
    then fails. Returns the rows it yielded."""
    yielded = _yielded_rows(brightStars.scatter_layer)
    undo = worker_patches.patch_everywhere(monkeypatch, brightStars, "scatter_layer",
                                           "tests.test_gen_bright_scatter_edges:layer_zero_fails_after_commits")
    with pytest.raises(RuntimeError, match="scatter worker died"):
        generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    undo()
    return yielded


def _layer_rows(rows, layer_index):
    return [row for row in rows if row[1] == layer_index]


def test_a_scatter_that_fails_after_two_commits_keeps_them_and_no_seed(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    yielded = _interrupt_after_commits(mysql_config, monkeypatch)
    stored = _stored_rows(mysql_config)
    # Layer 0 keeps its two commits and loses the rest. With one worker
    # the run stops there; with more, the layers other workers were
    # drawing at the time finish and keep their stars.
    assert _comparable(_layer_rows(stored, 0)) == _comparable(yielded[:20_000])
    if worker_patches.workers() == 1:
        assert len(stored) == 20_000
    # The threshold and seed are written only once every layer is in.
    assert _settings(mysql_config) is None


def test_a_re_plan_after_an_interrupted_scatter_holds_exactly_one_scatter(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    _interrupt_after_commits(mysql_config, monkeypatch)
    summary = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))

    seed = _scatter_seed(mysql_config)
    expected = list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, THRESHOLD, seed))
    stored = _stored_rows(mysql_config)
    assert summary["total"] == len(expected) == len(stored)
    assert _comparable(stored) == _comparable(expected)
    # Every star once: no two rows share a position.
    assert len({row[3:6] for row in stored}) == len(stored)
    assert _settings(mysql_config) == (THRESHOLD, seed)


def _leftover_cells(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT ring_index, layer_index, ring_slot_index, COUNT(*) AS n, MAX(id) AS top"
                            " FROM bright_stars GROUP BY ring_index, layer_index, ring_slot_index").fetchall()
    finally:
        conn.close()
    return {(row["ring_index"], row["layer_index"], row["ring_slot_index"]): row["n"] for row in rows}


def _block_of(address):
    ring_index, layer_index, slot = address
    return drill_parent(DrillBlock(1, ring_index, slot, layer_index))


def _fill(mysql_config, address):
    """Fills one sector, as a galaxy run does for each of its sectors
    (no backfill: since GEN.30 that runs once after the run's sectors)."""
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 0
    return generate.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)


def _backfill_around(mysql_config, address):
    """The backfill a `galaxy --slot` run of `address` ends with
    (`backfill_after_run`, GEN.30): around that sector, default tiers."""
    return generate.backfill_bright_stars(mysql_config, sector_position_pc(*address, EDGE_PC))


def _unbuilt_leftover_cells(mysql_config, last_leftover):
    conn = _db.get_connection(mysql_config)
    try:
        return {
            (row["ring_index"], row["layer_index"], row["ring_slot_index"])
            for row in conn.execute("SELECT ring_index, layer_index, ring_slot_index FROM bright_stars"
                                    " WHERE id <= ? AND star_system_id IS NULL", (last_leftover,)).fetchall()
        }
    finally:
        conn.close()


def _built_in(mysql_config, address):
    """`(ids of the stars built in this cell, its unbuilt star count)`."""
    conn = _db.get_connection(mysql_config)
    try:
        built = [row["id"] for row in conn.execute(
            "SELECT b.id FROM bright_stars b JOIN star_systems s ON s.id = b.star_system_id"
            " WHERE b.ring_index = ? AND b.layer_index = ? AND b.ring_slot_index = ?", address).fetchall()]
        unlinked = conn.execute(
            "SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
            " AND ring_slot_index = ? AND star_system_id IS NULL", address).fetchone()["n"]
    finally:
        conn.close()
    return built, unlinked


def test_a_fill_after_a_layer_failed_mid_scatter_builds_no_leftover_star(mysql_config, monkeypatch):
    # Layers 0 and -1 are written and committed, layer 1 fails: no
    # threshold is recorded. A run's fill then builds none of its cell's
    # leftovers (no level to fill down to), and the backfill after it
    # (GEN.30) draws its blocks from their tier floors with no ceiling.
    # Its blocks' other cells must lose their leftovers then, or a later
    # fill there would build them as well as the backfill's stars (every
    # bright star twice).
    _seed_galaxy(mysql_config)
    undo = worker_patches.patch_everywhere(monkeypatch, brightStars, "scatter_layer",
                                           "tests.test_gen_bright_scatter_edges:layer_fails", layer=1)
    with pytest.raises(RuntimeError, match="scatter worker died"):
        generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    undo()
    assert _settings(mysql_config) is None

    leftovers = _leftover_cells(mysql_config)
    # With more than one worker, layers handed out beside layer 1 may
    # finish too; layer 1 itself holds nothing.
    layers = {address[1] for address in leftovers}
    assert leftovers and {0, -1} <= layers and 1 not in layers
    if worker_patches.workers() == 1:
        assert layers == {0, -1}
    last_leftover = _count(mysql_config, "SELECT MAX(id) AS n FROM bright_stars")
    address = max((address for address in leftovers if address[1] == 0), key=lambda a: (leftovers[a], a))
    block_cells = set(generate._block_addresses(_block_of(address)))
    neighbor = max((cell for cell in block_cells if cell != address and leftovers.get(cell)),
                   key=lambda cell: (leftovers[cell], cell))

    _sector_id, _name, sector = _fill(mysql_config, address)
    assert [entry for entry in sector.entries if entry.preplaced] == []
    assert _count(mysql_config, "SELECT COUNT(*) AS n FROM bright_stars WHERE star_system_id IS NOT NULL") == 0
    assert _backfill_around(mysql_config, address)["stars"] > 0
    # The filled cell's own leftovers stay unbuilt (a filled sector never
    # gets stars); the block's other cells hold none.
    assert _unbuilt_leftover_cells(mysql_config, last_leftover) & block_cells == {address}

    _sector_id, _name, sector = _fill(mysql_config, neighbor)
    preplaced = [entry for entry in sector.entries if entry.preplaced]
    built, unlinked = _built_in(mysql_config, neighbor)
    assert unlinked == 0
    assert preplaced and len(built) == len(preplaced)
    # Every star built here came from the backfill, none from the
    # unfinished scatter.
    assert all(star_id > last_leftover for star_id in built)


def test_a_fill_after_failed_commits_clears_its_blocks_leftovers(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    _interrupt_after_commits(mysql_config, monkeypatch)
    leftovers = _leftover_cells(mysql_config)
    # The layer-0 cell with the fewest leftovers, so the fill stays small
    # even if it did build them.
    address = min((cell for cell in leftovers if cell[1] == 0), key=lambda cell: (leftovers[cell], cell))
    block = _block_of(address)
    block_cells = set(generate._block_addresses(block))
    last_leftover = _count(mysql_config, "SELECT MAX(id) AS n FROM bright_stars")
    # With more workers, the layers drawn beside layer 0 left theirs too.
    assert last_leftover == 20_000 if worker_patches.workers() == 1 else last_leftover >= 20_000
    _sector_id, _name, sector = _fill(mysql_config, address)
    assert [entry for entry in sector.entries if entry.preplaced] == []
    _backfill_around(mysql_config, address)

    assert _count(mysql_config, "SELECT COUNT(*) AS n FROM bright_stars WHERE star_system_id IS NOT NULL") == 0
    after = _leftover_cells(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        backfilled = {cell for key in _db.bright_star_block_keys(conn)
                      for cell in generate._block_addresses(DrillBlock(3, *key))}
    finally:
        conn.close()
    stale_cells = _unbuilt_leftover_cells(mysql_config, last_leftover)
    assert block_cells <= backfilled
    assert stale_cells, "leftovers outside the backfilled blocks should still be there"
    # Only the filled cell keeps its leftovers inside the backfilled
    # blocks: it never gets stars again, so they are never built.
    assert stale_cells & backfilled == {address}
    # Leftovers outside the backfilled blocks are untouched until their
    # own block is reached (or a re-plan clears them).
    assert {cell for cell in leftovers if cell not in backfilled} <= set(after)


BAND_FLOOR = 300.0


def test_re_running_an_interrupted_band_holds_the_band_once(mysql_config, monkeypatch):
    # GEN.32: the re-run starts the band over, so the layers the first
    # run finished don't hold it twice.
    _seed_galaxy(mysql_config)
    seed = _scatter_seed(mysql_config)
    generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    first = _count(mysql_config)
    first_rows = _stored_rows(mysql_config)
    undo = worker_patches.patch_everywhere(monkeypatch, brightStars, "scatter_layer",
                                           "tests.test_gen_bright_scatter_edges:layer_fails", layer=1)
    with pytest.raises(RuntimeError, match="scatter worker died"):
        generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", str(BAND_FLOOR)))
    undo()
    assert _settings(mysql_config) == (THRESHOLD, seed)
    assert _count(mysql_config) > first

    band = generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", str(BAND_FLOOR)))
    assert band is not None
    assert _settings(mysql_config) == (BAND_FLOOR, seed)
    band_seed = _scatter_seed(mysql_config, f"band/{BAND_FLOOR:g}-{THRESHOLD:g}")
    expected = list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, BAND_FLOOR, band_seed,
                                        max_luminosity_sol=THRESHOLD))
    stored = _stored_rows(mysql_config)
    assert len(stored) == first + len(expected)
    # The first scatter's stars and the band's, each once.
    assert _comparable(stored) == _comparable(first_rows + expected)
