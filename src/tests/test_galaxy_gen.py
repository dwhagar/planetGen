# tests/test_galaxy_gen.py

"""
End-to-end tests for `planetgen`'s `galaxy` subcommand -- the
galaxy-scale sector generation logic (`docs/design/galaxy/spherical-sector-coordinates.md`,
`docs/TODO.md` Phase 4). `test_galaxy_spherical.py` already covers
`galaxy/spherical.py`'s pure functions (shell/band/wedge math, the
enumeration primitive) in isolation; this file is the missing piece --
actually running the `galaxy` subcommand's modes (`--shell`,
`--shell --band`, `--shell --band --wedge`, the spans, `--center-sector`,
random start) against a real, throwaway MySQL database (see `conftest.py`'s
`mysql_config` fixture) and inspecting the resulting `sectors` rows.

The flags themselves (what is refused, and with what words) are pinned in
`test_gen_cli_errors.py`; the galaxy's address maths in
`test_galaxy_spherical.py`.

Runs the real CLI entry point (`main()`, via `sys.argv`) rather than
calling internals directly -- this repo's existing convention for testing
a CLI end-to-end (as opposed to `test_sector_gen.py`, which tests
`generate_sector_name` as a plain function because it has no database/CLI
side effects to speak of).

`--num-systems 1` keeps every generated sector's own system-generation
cost to a minimum -- these tests are about sector *placement*/persistence
(shell/band/wedge addresses, positions, occupied-sector skipping), not
system content, which is already covered extensively elsewhere
(`test_systems.py`, `test_db_persistence.py`).

Every test here takes the `mysql_config` fixture -- they're skipped, not
failed, when no MySQL test server is configured/reachable.
"""

import argparse
import math
import os
import sys
from unittest import mock

import pymysql
import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import run_galaxy
from planetgen.generation import run_plan
from planetgen.generation import run_sector
from planetgen.db import query
from planetgen.db import store as _db
from planetgen.db import store
from planetgen import tuning
from planetgen.galaxy import profile as galaxy_profile, resolution as galaxy_resolution, spherical, structure_draw
from planetgen.galaxy.density import FAMILY_COLUMNS, predicted_star_count, relative_density
from planetgen.galaxy.frame import galaxy_to_local_pc
from planetgen.galaxy.spherical import (
    SectorCell, enumerate_sectors_within_radius, sector_position_pc, wedge_count,
)
from planetgen.physics.units import ly_to_pc, mpc_to_pc

from tests.toy_shape import toy_shape
from tests import worker_patches
from tests.bughunt_support import forced_system_config
from tests.sector_support import make_sector

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)

# A small, fast-to-evaluate toy shape -- same parameters test_galaxy_skeleton.py
# uses, deep enough in its own bulge at the core that every cell there
# qualifies under any reasonable threshold, with a real, findable edge
# a few dozen shells out (not real Milky-Way scale, which would make
# these tests slow for no benefit -- the generation logic
# doesn't care about scale).
# A plan of a toy galaxy holds it in relative units, density 1 at its reference radius (GEN.404).
pytestmark = pytest.mark.usefixtures("unit_density_plans")


@pytest.fixture(autouse=True)
def planetless_systems():
    """These tests are about which sectors a run makes and where, not what is in them. A toy sector holds 100 to 350 systems, and a
    system with its planets, moons and solved air takes seconds, so a sector took minutes and the file half an hour; with no planets it
    takes seconds."""
    with forced_system_config(PLANETS=False):
        yield

# A stored shape carries no thick disk, bar or halo values: the galaxy seed draws them on load (GEN.532, GEN.533). The
# toy shape is therefore built with the structure its fixed seed draws, and `_seed_skeleton` stores that seed.
_SKELETON_SEED = bytes(range(16))
_SKELETON_SHAPE = toy_shape(
    disk_scale_length_pc=40.0,
    disk_scale_height_pc=12.0,
    bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15),
    arm_amplitude=0.4,
    structure=structure_draw.draw_structure(_SKELETON_SEED),
)



def poisson_always_one(_real):
    """`_sample_poisson_count` that draws one of anything with a positive
    mean (`worker_patches.patch_everywhere` factory)."""

    def sample(mean, rng=None):
        return 1 if mean > 0 else 0

    return sample


def capture_generate_sector(_real, path):
    """A `generate_sector` that records each sector's `(density,
    num_systems)` in `path` and makes an empty sector."""

    def fake(args, galactic_center_dist_ly=None, cell=None, fill=None, cloud_field=None):
        worker_patches.append_json(path, [args.density, args.num_systems])
        from planetgen.galaxy.sector import SpaceSector
        return "Fake Sector", SpaceSector(name="Fake Sector")

    return fake


def _shells(outer_shell, high_band):
    """Outline rows `(shell, highest band)` for shells `0 .. outer_shell`,
    each reaching bands `-high_band .. high_band` (as far as the shell has
    them: the core holds only band 0, shell 1 bands -1 to 1)."""
    return [(shell, min(high_band, spherical.band_range(shell)[1])) for shell in range(outer_shell + 1)]


WIDE_GALAXY_HIGH_BAND = 50
WIDE_GALAXY_OUTER_SHELL = 999


def _plan_wide_galaxy(mysql_config, galaxy_seed=None):
    """Seeds a skeleton whose outline (bands +/-50, out to shell 999) holds
    every address these tests generate with an explicit --num-systems --
    generation now refuses any address outside the stored outline."""
    _seed_skeleton(mysql_config, shells=_shells(WIDE_GALAXY_OUTER_SHELL, WIDE_GALAXY_HIGH_BAND),
                   galaxy_seed=galaxy_seed)


def _seed_skeleton(mysql_config, shape=_SKELETON_SHAPE, outer_shell_index=999, e_value=1.0, shells=(),
                   galaxy_seed=None, reach_shell_index=65535):
    """
    Directly writes a `galaxy_shape`/`galaxy_shell` skeleton, without
    running `planetgen plan`'s own scan -- these tests exercise
    the generation logic against a skeleton it can
    already read, not `build_shell_extents`' search (covered separately in
    `test_galaxy_skeleton.py`).

    Args:
        mysql_config (MySQLConfig): The fixture's throwaway database.
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        outer_shell_index (int): See `_db.save_galaxy_shape`.
        e_value (float): `expected_system_count_at_density_1`.
        shells (iterable): `(shell_index, high_band_index)` pairs -- see
            `_db.replace_galaxy_shells`.
        galaxy_seed (bytes, optional): The galaxy's seed; the stored one when there
            is one, else the fixed `_SKELETON_SEED`, so the loaded shape matches
            `_SKELETON_SHAPE`.
        reach_shell_index (int): The outermost shell the grid has addresses for.
    """
    if galaxy_seed is None:
        conn = _db.get_connection(mysql_config)
        try:
            galaxy_seed = _db.get_galaxy_seed(conn) or _SKELETON_SEED
        finally:
            conn.close()
    _db.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=outer_shell_index,
                           expected_system_count_at_density_1=e_value, config=mysql_config,
                           galaxy_seed=galaxy_seed, reach_shell_index=reach_shell_index)
    _db.replace_galaxy_shells(list(shells), config=mysql_config)


def _mysql_argv(mysql_config):
    """The `--mysql-*` argv fragment pointing a CLI invocation at the
    fixture's throwaway database."""
    return [
        "--mysql-host", mysql_config.host,
        "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user,
        "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database,
    ]


def _run_cli(argv, planets=None):
    """Runs `generate.main()` with the `sectors` subcommand (PERF.123: `galaxy` does every stage after the sectors, which
    these tests do not want) and the given argv (excluding argv[0]), restoring `sys.argv` afterward. `planets=False`
    makes every system planet-less (cheap), as `-planets` did before galaxy
    runs lost their forcing options (GEN.51)."""
    old_argv = sys.argv
    forced = {} if planets is None else {"PLANETS": planets}
    try:
        sys.argv = ["planetgen", "sectors"] + argv
        with forced_system_config(**forced):
            generate_cli.main()
    finally:
        sys.argv = old_argv


def _run_sector_gen_cli(argv):
    """Runs `generate.main()` with the `sector` subcommand the same way --
    used to seed an "unplaced" (no galaxy position) sector for
    `test_center_sector_without_galaxy_position_is_rejected`."""
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "sector"] + argv
        generate_cli.main()
    finally:
        sys.argv = old_argv


def _all_sectors(mysql_config):
    """Returns every `sectors` row. `get_connection`'s own `_ensure_schema`
    creates the table if a rejected run never got as far as calling it
    (`run_addresses`'s `LARGE_RUN_WARNING_THRESHOLD` guard raises
    before ever calling `_db.get_connection`), so this never errors on a
    fresh, empty database -- it just returns an empty list."""
    conn = _db.get_connection(mysql_config)
    try:
        return conn.execute(
            "SELECT id, name, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc, "
            "shell_index, band_index, wedge_index FROM sectors"
        ).fetchall()
    finally:
        conn.close()


def _address(row):
    """A `sectors` row's `(shell, band, wedge)` address."""
    return row["shell_index"], row["band_index"], row["wedge_index"]


def _band_cells(shell, band):
    """Every address of one band of a shell."""
    return {(shell, band, wedge) for wedge in range(wedge_count(shell, band))}


def test_galaxy_density_shape_returns_none_before_a_skeleton_is_built(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        assert query.galaxy_density_shape(conn) is None
    finally:
        conn.close()


def test_galaxy_density_shape_serializes_the_stored_skeleton(mysql_config):
    # The Galaxy Map (html/galaxy.py, via GET /api/galaxy/shape) reads this
    # exact dict to shade its "expected density" cloud from the real model
    # instead of the illustrative fallback gradient -- confirms the shape
    # this function hands back round-trips every field `GalaxyShape`/
    # `_db.get_galaxy_shape` produce, not just that it returns *something*.
    _seed_skeleton(mysql_config, outer_shell_index=42, e_value=2.5)

    conn = _db.get_connection(mysql_config)
    try:
        shape_dict = query.galaxy_density_shape(conn)
    finally:
        conn.close()

    assert shape_dict is not None
    for field in _SKELETON_SHAPE._fields:
        if field in ("components", "masses", "structure", "system_mass_scale"):  # objects for the density model, not JSON for the map
            assert field not in shape_dict
            continue
        if field == "spheroid_orientation_rad" or field in FAMILY_COLUMNS:  # DB.66 columns the map does not read (a name the terms share keeps theirs)
            continue
        assert shape_dict[field] == pytest.approx(getattr(_SKELETON_SHAPE, field))
    assert shape_dict["edge_pc"] == pytest.approx(EDGE_PC)
    assert shape_dict["outer_shell_index"] == 42
    assert shape_dict["expected_system_count_at_density_1"] == pytest.approx(2.5)


def test_band_mode_generates_every_wedge_and_is_idempotent(mysql_config):
    _plan_wide_galaxy(mysql_config)
    # Shell 1's equator band holds exactly 12 wedges -- small enough to run
    # under LARGE_RUN_WARNING_THRESHOLD with no --limit/--yes needed.
    assert wedge_count(1, 0) == 12
    _run_cli(["--shell", "1", "--band", "0", "--num-systems", "1"] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    assert len(sectors) == 12
    assert {_address(row) for row in sectors} == _band_cells(1, 0)
    for row in sectors:
        assert row["center_x_pc"] is not None
        assert row["galactic_radius_pc"] == pytest.approx(
            (row["center_x_pc"] ** 2 + row["center_y_pc"] ** 2 + row["center_z_pc"] ** 2) ** 0.5
        )

    # Re-running the same band must skip every already-occupied wedge
    # (get_occupied_addresses) rather than generating duplicates.
    _run_cli(["--shell", "1", "--band", "0", "--num-systems", "1"] + _mysql_argv(mysql_config))
    sectors_after_rerun = _all_sectors(mysql_config)
    assert len(sectors_after_rerun) == 12
    assert {row["id"] for row in sectors_after_rerun} == {row["id"] for row in sectors}


def test_band_mode_limit_generates_only_first_n_missing_wedges(mysql_config):
    _plan_wide_galaxy(mysql_config)
    _run_cli(["--shell", "1", "--band", "0", "--limit", "1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    sectors = _all_sectors(mysql_config)
    assert len(sectors) == 1
    assert sectors[0]["shell_index"] == 1

    # A second --limit 1 run must generate the *next* missing wedge, not
    # regenerate the one that already exists.
    _run_cli(["--shell", "1", "--band", "0", "--limit", "1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    sectors_after_second_run = _all_sectors(mysql_config)
    assert len(sectors_after_second_run) == 2
    addresses = {_address(row) for row in sectors_after_second_run}
    assert len(addresses) == 2


def test_band_mode_honors_the_band_and_places_the_sectors_at_its_height(mysql_config):
    _plan_wide_galaxy(mysql_config)
    _run_cli(["--shell", "2", "--band", "-2", "--num-systems", "1"] + _mysql_argv(mysql_config))
    sectors = _all_sectors(mysql_config)
    assert {_address(row) for row in sectors} == _band_cells(2, -2)
    for row in sectors:
        center = sector_position_pc(*_address(row), EDGE_PC)
        assert (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) == pytest.approx(center)
        assert row["center_z_pc"] < 0


def test_local_neighborhood_mode_generates_expected_new_sectors_and_skips_occupied(mysql_config):
    _plan_wide_galaxy(mysql_config)
    # Seed with just one already galaxy-placed sector, shell 1 wedge 0 of the
    # equator -- --center-sector needs an already galaxy-placed sector to
    # search a neighborhood around (a 'planetgen sector'-standalone sector
    # won't do, see test_center_sector_without_galaxy_position_is_rejected
    # below).
    _run_cli(["--shell", "1", "--band", "0", "--limit", "1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    seeded = _all_sectors(mysql_config)
    assert len(seeded) == 1
    center_row = seeded[0]
    center_id = center_row["id"]
    center = (center_row["center_x_pc"], center_row["center_y_pc"], center_row["center_z_pc"])
    seed_address = _address(center_row)

    radius_pc = 8.0
    expected_addresses = {
        (shell_index, band_index, wedge_index)
        for shell_index, band_index, wedge_index, _x, _y, _z, _dist in
        enumerate_sectors_within_radius(center, radius_pc, EDGE_PC)
    }
    # Sanity check on the test's own setup: this radius must reach past the
    # seed's shell and off the equator, and cover more than just the seed
    # sector itself -- otherwise this wouldn't actually exercise cross-shell
    # neighborhood generation at all.
    assert any(shell_index != 1 for shell_index, _band, _wedge in expected_addresses)
    assert any(band_index != 0 for _shell, band_index, _wedge in expected_addresses)
    assert len(expected_addresses) > 1
    assert seed_address in expected_addresses  # the seed sector itself, at distance 0

    _run_cli([
        "--center-sector", str(center_id), "--radius-pc", str(radius_pc),
        "--num-systems", "1",
    ] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    actual_addresses = {_address(row) for row in sectors}
    assert actual_addresses == expected_addresses
    # No duplicates: exactly one row per address, including the seed's own.
    assert len(sectors) == len(actual_addresses)

    # Every newly placed sector's stored position matches the geometry it
    # was supposed to be placed at.
    by_address = {_address(row): row for row in sectors}
    for shell_index, band_index, wedge_index, x, y, z, _dist in enumerate_sectors_within_radius(
            center, radius_pc, EDGE_PC):
        row = by_address[(shell_index, band_index, wedge_index)]
        assert row["center_x_pc"] == pytest.approx(x)
        assert row["center_y_pc"] == pytest.approx(y)
        assert row["center_z_pc"] == pytest.approx(z)

    # Re-running the identical neighborhood must not create duplicates or
    # regenerate anything -- every candidate sector is now already occupied.
    _run_cli([
        "--center-sector", str(center_id), "--radius-pc", str(radius_pc),
        "--num-systems", "1",
    ] + _mysql_argv(mysql_config))
    sectors_after_rerun = _all_sectors(mysql_config)
    assert {row["id"] for row in sectors_after_rerun} == {row["id"] for row in sectors}


def test_the_neighborhood_is_a_ball_of_sectors_not_a_cylinder(mysql_config):
    """ADM.30 asked for a cylinder of X layers; the grid is spherical now, so a neighborhood is every sector the
    ball of `--radius-pc` touches, however far up or down that reaches."""
    _plan_wide_galaxy(mysql_config)
    _run_cli(["--shell", "30", "--band", "0", "--wedge", "0"] + _mysql_argv(mysql_config), planets=False)
    center_row = _all_sectors(mysql_config)[0]
    center = (center_row["center_x_pc"], center_row["center_y_pc"], center_row["center_z_pc"])
    _run_cli(["--center-sector", str(center_row["id"]), "--radius-pc", "6", "--num-systems", "1"]
             + _mysql_argv(mysql_config))
    addresses = {_address(row) for row in _all_sectors(mysql_config)}
    assert {a[1] for a in addresses} >= {-1, 0, 1}
    for address in addresses:
        assert math.dist(sector_position_pc(*address, EDGE_PC), center) <= 6.0 + EDGE_PC * math.sqrt(3)


def test_center_sector_without_galaxy_position_is_rejected(mysql_config):
    """--center-sector must refuse a sector that was never placed in a
    galaxy (e.g. one generated by 'planetgen sector's own standalone CLI) --
    planetgen's run_local_neighborhood docstring calls this out
    explicitly, since there's no center point to search a neighborhood
    around."""
    _run_sector_gen_cli(["--num-systems", "1"] + _mysql_argv(mysql_config))

    unplaced = _all_sectors(mysql_config)
    assert len(unplaced) == 1
    assert unplaced[0]["shell_index"] is None
    assert unplaced[0]["center_x_pc"] is None

    with pytest.raises(SystemExit):
        _run_cli([
            "--center-sector", str(unplaced[0]["id"]), "--radius-pc", "5.0",
            "--num-systems", "1",
        ] + _mysql_argv(mysql_config))

    # Rejected before anything was generated -- still just the one
    # unplaced sector.
    assert len(_all_sectors(mysql_config)) == 1


def test_band_mode_rejects_a_large_band_without_limit_or_yes_under_strict(mysql_config, monkeypatch):
    """A band of shell 400 holds more than the (lowered) threshold -- under
    --strict a bare --shell 400 --band 0 is refused rather than attempting to
    generate all of them."""
    monkeypatch.setattr(run_galaxy, "LARGE_RUN_WARNING_THRESHOLD", 100)
    assert wedge_count(400, 0) > run_galaxy.LARGE_RUN_WARNING_THRESHOLD
    _plan_wide_galaxy(mysql_config)

    with pytest.raises(SystemExit):
        _run_cli(["--shell", "400", "--band", "0", "--num-systems", "1", "--strict"] + _mysql_argv(mysql_config))

    assert len(_all_sectors(mysql_config)) == 0

    # --limit bypasses the guard even without --yes.
    _run_cli(["--shell", "400", "--band", "0", "--limit", "2", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == 2


@pytest.mark.parametrize("argv", [
    ["--shell", "1", "--band", "0"],
    ["--shell", "2"],
    ["--shells", "1:2"],
])
def test_a_large_run_without_limit_or_yes_warns_and_goes_ahead(mysql_config, monkeypatch, capsys, argv):
    # GEN.81: a threshold of 3 makes these "large"; they warn and run whole.
    monkeypatch.setattr(run_galaxy, "LARGE_RUN_WARNING_THRESHOLD", 3)
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 2)])
    _run_cli(argv + ["--num-systems", "1"] + _mysql_argv(mysql_config), planets=False)
    assert len(_all_sectors(mysql_config)) > 3
    assert "a very large run" in " ".join(capsys.readouterr().out.split())


def test_shell_mode_generates_every_band_of_a_shell_the_outline_holds(mysql_config):
    # Shell 3 reaches only bands -1 to 1 here; the grid has bands -4 to 4.
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 2), (3, 1)])
    _run_cli(["--shell", "3", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert {_address(row) for row in _all_sectors(mysql_config)} == (
        _band_cells(3, -1) | _band_cells(3, 0) | _band_cells(3, 1))

    # Idempotent.
    _run_cli(["--shell", "3", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == sum(wedge_count(3, band) for band in (-1, 0, 1))


def test_shell_mode_with_no_band_of_the_outline_in_it_generates_nothing(mysql_config, capsys):
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1)])
    _run_cli(["--shell", "5", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []
    assert "nothing to generate" in capsys.readouterr().out


def test_span_mode_fills_a_range_of_shells_and_bands(mysql_config):
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 2)])
    _run_cli(["--shells", "1:2", "--bands", "0:1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert {_address(row) for row in _all_sectors(mysql_config)} == (
        _band_cells(1, 0) | _band_cells(1, 1) | _band_cells(2, 0) | _band_cells(2, 1))


def _level3_block_holding(shell, band, wedge):
    """The level-3 block of the ladder that holds sector `(shell, band, wedge)`, as the key the command line takes."""
    from planetgen.galaxy import drill
    block = drill.drill_chain_of(shell, band, wedge)[-2]
    assert block.m == 3
    return block, drill.format_drill_key(block)


def test_block_mode_fills_the_sectors_the_outline_holds_inside_one_block(mysql_config):
    """GEN.310: `--block` generates the sectors of one block of the ladder that the stored outline holds, nothing else, and
    a second run finds them all there."""
    from planetgen.galaxy import drill
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 2), (3, 1), (4, 1), (5, 1)])
    block, key = _level3_block_holding(4, 0, 0)
    wanted = {address for address in drill.block_addresses(block) if address[1] in (-1, 0, 1) and address[0] <= 5}
    assert wanted and any(shell == 4 for shell, _band, _wedge in wanted)
    _run_cli(["--block", key, "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert {_address(row) for row in _all_sectors(mysql_config)} == wanted

    first = len(_all_sectors(mysql_config))
    _run_cli(["--block", key, "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == first


def test_block_mode_generates_the_inner_shells_first_and_stops_at_the_limit(mysql_config):
    from planetgen.galaxy import drill
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 2), (3, 1), (4, 1), (5, 1)])
    block, key = _level3_block_holding(4, 0, 0)
    ordered = [address for address in drill.block_addresses(block) if abs(address[1]) <= 1 and address[0] <= 5]
    assert [a[0] for a in ordered] == sorted(a[0] for a in ordered)
    _run_cli(["--block", key, "--limit", "3", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert {_address(row) for row in _all_sectors(mysql_config)} == set(ordered[:3])


def test_block_mode_with_a_block_the_outline_does_not_reach_generates_nothing(mysql_config, capsys):
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1)])
    _block, key = _level3_block_holding(40, 0, 0)
    _run_cli(["--block", key, "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []
    assert "holds no sector inside the galaxy's outline" in " ".join(capsys.readouterr().out.split())


def test_block_mode_needs_limit_or_yes_when_large_under_strict(mysql_config, monkeypatch):
    monkeypatch.setattr(run_galaxy, "LARGE_RUN_WARNING_THRESHOLD", 3)
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 2), (3, 1), (4, 1), (5, 1)])
    _block, key = _level3_block_holding(4, 0, 0)
    with pytest.raises(SystemExit):
        _run_cli(["--block", key, "--num-systems", "1", "--strict"] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []


def test_span_mode_fills_a_wrapping_wedge_arc(mysql_config):
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 2), (3, 1)])
    last = wedge_count(3, 0) - 1
    _run_cli(["--shells", "3", "--bands", "0", "--wedges", f"{last}:1", "--num-systems", "1"]
             + _mysql_argv(mysql_config))
    assert {_address(row) for row in _all_sectors(mysql_config)} == {(3, 0, last), (3, 0, 0), (3, 0, 1)}


def test_shell_mode_needs_limit_or_yes_when_large_under_strict(mysql_config):
    _plan_wide_galaxy(mysql_config)
    assert spherical.cell_count(30) > run_galaxy.LARGE_RUN_WARNING_THRESHOLD
    with pytest.raises(SystemExit):
        _run_cli(["--shell", "30", "--num-systems", "1", "--strict"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == 0
    _run_cli(["--shell", "30", "--limit", "2", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == 2


def test_single_address_with_radius_also_generates_its_neighborhood(mysql_config):
    # Shell 30 rather than the toy shape's dense bulge (shell 2 built a
    # neighborhood of thousands of systems, over a minute) -- the same
    # neighborhood logic, at a density that builds in a couple of seconds.
    _plan_wide_galaxy(mysql_config)
    _run_cli(["--shell", "30", "--band", "0", "--wedge", "4", "--radius-pc", "5"] + _mysql_argv(mysql_config))
    addresses = {_address(row) for row in _all_sectors(mysql_config)}
    assert (30, 0, 4) in addresses
    assert len(addresses) > 1
    center = sector_position_pc(30, 0, 4, EDGE_PC)
    assert addresses <= {(s, b, w) for s, b, w, *_ in enumerate_sectors_within_radius(center, 5.0, EDGE_PC)}


def test_a_block_is_refused_until_the_ladder_is_rebuilt():
    # MAP.205: the 243/81/27/9/3/1 ladder is rebuilt on the spherical grid after the switch.
    with pytest.raises(SystemExit):
        _run_cli(["--block", "3.2.1.0"] + _DUMMY_MYSQL_ARGV)


# ---------------------------------------------------------------------------
# --shell I --band J --wedge K -- single-address mode, the direct path from a
# designation/address copied out of the interactive 3D Galaxy Map
# (planetgen/web/maps/galaxymap3d.py) into this script. Thin wrapper around
# the sector made at one address (covered in isolation below) plus its
# own argparse validation (pinned in test_gen_cli_errors.py, which needs no
# real database at all).
# ---------------------------------------------------------------------------

_DUMMY_MYSQL_ARGV = [
    "--mysql-host", "localhost", "--mysql-port", "3306",
    "--mysql-user", "x", "--mysql-password", "x", "--mysql-database", "x",
]


def test_wedge_mode_generates_exactly_the_one_requested_sector(mysql_config):
    _seed_skeleton(mysql_config, shells=_shells(3, 3))

    _run_cli(["--shell", "2", "--band", "2", "--wedge", "1"] + _mysql_argv(mysql_config), planets=False)

    sectors = _all_sectors(mysql_config)
    assert len(sectors) == 1
    assert _address(sectors[0]) == (2, 2, 1)
    assert sectors[0]["center_x_pc"] is not None

    # Re-running the identical address must reuse it, not duplicate it.
    _run_cli(["--shell", "2", "--band", "2", "--wedge", "1"] + _mysql_argv(mysql_config), planets=False)
    sectors_after_rerun = _all_sectors(mysql_config)
    assert len(sectors_after_rerun) == 1
    assert sectors_after_rerun[0]["id"] == sectors[0]["id"]


def test_wedge_mode_rejects_an_out_of_range_wedge(mysql_config):
    _seed_skeleton(mysql_config, shells=_shells(2, 2))

    with pytest.raises(SystemExit):
        _run_cli(["--shell", "1", "--band", "0", "--wedge", str(wedge_count(1, 0))] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []


def test_wedge_mode_rejects_an_address_when_there_is_no_outline(mysql_config):
    # No stored shell at all -- the plan is empty, so there is nothing to follow.
    _seed_skeleton(mysql_config, shells=[])

    with pytest.raises(SystemExit):
        _run_cli(["--shell", "5000", "--band", "0", "--wedge", "0"] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []


# ---------------------------------------------------------------------------
# make_sector (tests/sector_support.py) -- the sector of one address, built on top of the 'planetgen plan' skeleton (galaxy_shape/galaxy_shell),
# rather than an explicit --shell/--center-sector batch.
# ---------------------------------------------------------------------------

def test_a_sector_made_raises_without_a_skeleton(mysql_config):
    with pytest.raises(RuntimeError):
        make_sector(0, 0, 0, config=mysql_config)


def test_a_sector_made_creates_then_reuses_the_same_sector(mysql_config):
    # The core sits deep in this toy shape's own bulge -- it clears even a
    # demanding threshold.
    _seed_skeleton(mysql_config, shells=_shells(2, 1))

    first = make_sector(0, 0, 0, config=mysql_config)
    assert first["created"] is True, first
    assert first["qualifies"] is True, first
    assert first["sector_id"] is not None, first
    assert first["sector_name"], first

    second = make_sector(0, 0, 0, config=mysql_config)
    assert second["created"] is False, second
    assert second["qualifies"] is True
    assert second["sector_id"] == first["sector_id"], (first, second)

    assert len(_all_sectors(mysql_config)) == 1, _all_sectors(mysql_config)


def test_a_sector_generated_on_a_visit_is_in_the_block_totals(mysql_config):
    """DB.83: the blocks above a sector made on the spot hold it at once, the same totals `planetgen
    fill-block-totals` (`block_stats.rebuild`) would compute."""
    from planetgen.db import block_stats
    columns = "m, shell_index, band_index, wedge_index, sector_count, systems, expected_systems, stars, habitable_worlds"
    _seed_skeleton(mysql_config, shells=_shells(2, 1))
    assert make_sector(0, 0, 0, config=mysql_config)["created"] is True
    conn = _db.get_connection(mysql_config)
    try:
        kept = [tuple(r.values()) for r in conn.execute(f"SELECT {columns} FROM galaxy_block_stats ORDER BY 1, 2, 3, 4").fetchall()]
        block_stats.rebuild(conn)
        rebuilt = [tuple(r.values()) for r in conn.execute(f"SELECT {columns} FROM galaxy_block_stats ORDER BY 1, 2, 3, 4").fetchall()]
    finally:
        conn.close()
    assert kept and kept == rebuilt, (kept, rebuilt)


def test_a_sector_made_refuses_only_beyond_the_reach(mysql_config):
    # GEN.235: the reach is the only refusal. A shell past it has no address, and neither has a band past what
    # its shell holds; the skeleton's reach is exact, so this is a certain "no", no density check needed.
    _seed_skeleton(mysql_config, shells=_shells(2, 1), reach_shell_index=100)

    no_content = {"created": False, "qualifies": False, "in_reach": False, "sector_id": None, "sector_name": None}
    assert make_sector(5000, 0, 0, config=mysql_config) == no_content
    assert make_sector(101, 0, 0, config=mysql_config) == no_content
    assert make_sector(2, 9, 0, config=mysql_config) == no_content
    assert _all_sectors(mysql_config) == []


def test_a_sector_made_saves_a_halo_sector_from_its_profile_even_when_empty(mysql_config, monkeypatch):
    """GEN.235: a cell the outline does not hold but the reach does is generated from its profile when asked, and
    saved and marked generated though its draw comes out empty; the second ask returns the same sector."""
    _seed_skeleton(mysql_config, shells=_shells(2, 1), reach_shell_index=300)
    address = (200, 0, 3)     # 800 pc out in the toy galaxy (a scale length is 40 pc): the halo
    densities = []
    real = run_galaxy.generate_and_save_sector_at

    def spy(args, *rest, **kwargs):
        densities.append(args.density)
        return real(args, *rest, **kwargs)

    monkeypatch.setattr(run_galaxy, "generate_and_save_sector_at", spy)
    first = make_sector(*address, config=mysql_config, backfill=False, settle=False)
    assert first["created"] is True and first["qualifies"] is False and first["in_reach"] is True, first
    assert first["sector_id"] is not None
    conn = _db.get_connection(mysql_config)
    try:
        skeleton = _db.get_galaxy_shape(conn)
    finally:
        conn.close()
    expected = galaxy_profile.sector_profile(address, skeleton)
    assert expected.region == "halo" and 0.0 < expected.density < 1e-3
    assert densities == [expected.density]
    assert {_address(row) for row in _all_sectors(mysql_config)} == {address}

    again = make_sector(*address, config=mysql_config, backfill=False, settle=False)
    assert again["created"] is False and again["sector_id"] == first["sector_id"]


def test_a_sector_made_generates_sparse_wedges_within_a_stored_band(mysql_config):
    """
    GEN.76: inside the galaxy's stored outline every address generates,
    however sparse. The predicted density only sets the draw's expected
    count; a wedge far below one expected star still gets a sector, and
    it is marked generated (`sector_stats.bright_level_sol` = 0) like
    any other.
    """
    shell_index = 5
    n_k = wedge_count(shell_index, 0)

    densities = sorted(
        ((relative_density(sector_position_pc(shell_index, 0, i, EDGE_PC), _SKELETON_SHAPE), i)
         for i in range(n_k)),
        reverse=True,
    )
    densest_wedge = densities[0][1]
    sparsest_wedge = densities[-1][1]
    threshold_rho = (densities[0][0] + densities[-1][0]) / 2.0
    assert densities[0][0] >= threshold_rho > densities[-1][0], (
        "test setup needs a real spread of densities within this band"
    )

    # The sparsest wedge predicts under one star at this E value -- the
    # old rule skipped it.
    _seed_skeleton(mysql_config, e_value=1.0 / threshold_rho, shells=_shells(shell_index, 5))

    dense_result = make_sector(shell_index, 0, densest_wedge, config=mysql_config)
    assert dense_result["qualifies"] is True
    assert dense_result["created"] is True

    sparse_result = make_sector(shell_index, 0, sparsest_wedge, config=mysql_config)
    assert sparse_result["qualifies"] is True
    assert sparse_result["created"] is True
    assert sparse_result["sector_id"] is not None

    assert len(_all_sectors(mysql_config)) == 2

    conn = _db.get_connection(mysql_config)
    try:
        levels = _db.sector_bright_levels(conn, [(shell_index, 0, densest_wedge), (shell_index, 0, sparsest_wedge)])
    finally:
        conn.close()
    assert levels.get((shell_index, 0, sparsest_wedge)) == 0
    assert levels.get((shell_index, 0, densest_wedge)) == 0


def test_a_sector_made_generates_a_near_empty_density_sector(mysql_config):
    """
    GEN.76: a sector far out in the disk, where the density predicts
    almost nothing, still runs its draw, is saved, and is marked
    generated -- never refused.
    """
    address = (900, 0, 0)
    assert relative_density(sector_position_pc(*address, EDGE_PC), _SKELETON_SHAPE) < 0.01, (
        "test setup needs a nearly empty position"
    )
    _seed_skeleton(mysql_config, shells=_shells(999, 0))

    result = make_sector(*address, config=mysql_config)
    assert result["created"] is True and result["qualifies"] is True

    conn = _db.get_connection(mysql_config)
    try:
        levels = _db.sector_bright_levels(conn, [address])
    finally:
        conn.close()
    assert levels.get(address) == 0


def test_a_neighborhood_from_a_sparse_centre_at_the_edge_fills_every_sector_in_range(mysql_config):
    """
    GEN.77: a neighborhood run whose centre sector predicts almost no
    stars (here at the galaxy's outer edge) still generates every sector
    in range that the outline holds, density-driven, none skipped.
    """
    centre = (995, 0, 0)
    position = sector_position_pc(*centre, EDGE_PC)
    assert relative_density(position, _SKELETON_SHAPE) < 0.01, "test setup needs a nearly empty centre"
    _seed_skeleton(mysql_config, shells=_shells(999, 1))
    _run_cli(["--shell", "995", "--band", "0", "--wedge", "0", "--radius-pc", "6"] + _mysql_argv(mysql_config),
             planets=False)
    generated = {_address(row) for row in _all_sectors(mysql_config)}
    in_range = {(s, b, w) for s, b, w, *_ in enumerate_sectors_within_radius(position, 6.0, EDGE_PC)
                if s <= 999 and -1 <= b <= 1}
    assert centre in generated
    assert len(in_range) > 5
    assert generated == in_range


def test_a_sector_made_passes_relative_density_as_the_density_multiplier(mysql_config, monkeypatch):
    """
    A qualifying sector's actual system count should be driven by its own
    `relative_density` (from the stored skeleton), not a uniform default
    -- confirmed by intercepting `run_sector.generate_sector` and checking
    what `args.density` it was actually called with, rather than relying
    on the Poisson-sampled system count to differ across a single run.
    """
    _seed_skeleton(mysql_config, shells=_shells(2, 1))

    captured = {}

    def _fake_generate_sector(args, galactic_center_dist_ly=None, cell=None, fill=None, cloud_field=None):
        captured["density"] = args.density
        captured["num_systems"] = args.num_systems
        from planetgen.galaxy.sector import SpaceSector
        return "Fake Sector", SpaceSector(name="Fake Sector")

    monkeypatch.setattr(run_sector, "generate_sector", _fake_generate_sector)

    position = sector_position_pc(0, 0, 0, EDGE_PC)
    expected_density = relative_density(position, _SKELETON_SHAPE)

    result = make_sector(0, 0, 0, config=mysql_config)
    assert result["created"] is True
    assert captured["density"] == pytest.approx(expected_density)
    assert captured["num_systems"] is None


def test_band_batch_uses_skeleton_density_when_neither_flag_given(mysql_config, monkeypatch, tmp_path):
    """
    'planetgen galaxy --shell N --band J' with neither --density nor
    --num-systems should drive each sector's own system count from the galaxy
    skeleton's real position-based relative_density (_BatchDensity), the same
    way a sector made at one address does -- not silently fall back to a flat count shared by the whole
    band (validate_shared_generation_args's old behavior for every
    subcommand, still correct for the position-less 'sector' subcommand,
    which keeps defaulting to 10 -- galaxy mode alone leaves both None so
    this per-sector path can take over).
    """
    n_cells = wedge_count(1, 0)
    _seed_skeleton(mysql_config, shells=_shells(2, 1))

    captured_file = str(tmp_path / "captured.jsonl")
    worker_patches.patch_everywhere(monkeypatch, run_sector, "generate_sector",
                                    "tests.test_galaxy_gen:capture_generate_sector", path=captured_file)

    _run_cli(["--shell", "1", "--band", "0"] + _mysql_argv(mysql_config))

    captured = [tuple(item) for item in worker_patches.read_json_lines(captured_file)]
    assert len(captured) == n_cells
    # Workers can finish the sectors in any order: compare the sets.
    expected = sorted(relative_density(sector_position_pc(1, 0, wedge, EDGE_PC), _SKELETON_SHAPE)
                      for wedge in range(n_cells))
    assert sorted(density for density, _num in captured) == pytest.approx(expected)
    assert all(num_systems is None for _density, num_systems in captured)


def test_band_batch_explicit_num_systems_still_overrides_skeleton_density(mysql_config, monkeypatch, tmp_path):
    """An explicit --num-systems is still a uniform, intentional override
    for the whole batch -- _BatchDensity must not second-guess it."""
    n_cells = wedge_count(1, 0)
    _seed_skeleton(mysql_config, shells=_shells(2, 1))

    captured_file = str(tmp_path / "captured.jsonl")
    worker_patches.patch_everywhere(monkeypatch, run_sector, "generate_sector",
                                    "tests.test_galaxy_gen:capture_generate_sector", path=captured_file)

    _run_cli(["--shell", "1", "--band", "0", "--num-systems", "4"] + _mysql_argv(mysql_config))

    captured = [tuple(item) for item in worker_patches.read_json_lines(captured_file)]
    assert len(captured) == n_cells
    assert all(density is None and num_systems == 4 for density, num_systems in captured)


def test_sectors_table_rejects_duplicate_address(mysql_config):
    """The schema-level guarantee make_sector's own race
    recovery depends on: two sectors can never share a (shell_index,
    band_index, wedge_index) address."""
    conn = _db.get_connection(mysql_config)
    try:
        from planetgen.galaxy.sector import SpaceSector
        galaxy_position = {
            "center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 3.0,
            "galactic_radius_pc": 3.74, "shell_index": 1, "band_index": 0, "wedge_index": 0,
        }
        _db.insert_sector(conn, SpaceSector(name="First"), galaxy_position=galaxy_position)
        with pytest.raises(pymysql.err.IntegrityError):
            _db.insert_sector(conn, SpaceSector(name="Second"), galaxy_position=galaxy_position)
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Random-start mode (no --shell/--center-sector given) -- the galaxy subcommand's
# zero-argument default: pick a random, not-yet-occupied address, generate
# it, then generate every sector within --radius-pc (default 12 pc) of it.
# ---------------------------------------------------------------------------

def _plan_only_the_core_and_shell_one_equator(mysql_config):
    """Seeds a skeleton whose whole outline is the core and the equator band
    of shell 1, so every random-start draw (which only ever picks inside the
    outline, at or below `--max-shell 1`) lands on one of those 13 cells."""
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 0)])


_SMALL_OUTLINE = _band_cells(0, 0) | _band_cells(1, 0)
"""set: The 13 cells `_plan_only_the_core_and_shell_one_equator` plans."""


def test_process_args_with_no_arguments_resolves_to_random_start_mode():
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "galaxy"]
        args = generate_cli.process_args()
    finally:
        sys.argv = old_argv

    assert args.shell is None
    assert args.center_sector is None
    assert args.radius_pc is None
    assert args.max_shell is None


def test_process_args_accepts_radius_pc_and_max_shell_alone():
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "galaxy", "--radius-pc", "10", "--max-shell", "5"]
        args = generate_cli.process_args()
    finally:
        sys.argv = old_argv

    assert args.radius_pc == 10
    assert args.max_shell == 5


def test_process_args_accepts_min_start_density_alone():
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "galaxy", "--min-start-density", "1.5"]
        args = generate_cli.process_args()
    finally:
        sys.argv = old_argv

    assert args.min_start_density == 1.5


def test_random_start_mode_generates_a_seed_sector_and_its_neighborhood(mysql_config, monkeypatch):
    # --max-shell 1 puts the randomly chosen seed in the core or shell 1,
    # inside the outline -- a small --radius-pc keeps generation fast.
    _plan_wide_galaxy(mysql_config)
    radius_pc = 8.0
    _run_cli(
        ["--max-shell", "1", "--radius-pc", str(radius_pc), "--num-systems", "1"] + _mysql_argv(mysql_config)
    )

    sectors = _all_sectors(mysql_config)
    assert len(sectors) >= 1
    # The seed is the first row inserted (generate_and_save_sector_at runs
    # before run_local_neighborhood's own loop) -- every other row's
    # address must fall within radius_pc of its real galaxy-frame center.
    seed = min(sectors, key=lambda row: row["id"])
    assert seed["shell_index"] <= 1

    seed_center = (seed["center_x_pc"], seed["center_y_pc"], seed["center_z_pc"])
    expected_addresses = {
        (shell_index, band_index, wedge_index)
        for shell_index, band_index, wedge_index, _x, _y, _z, _dist in
        enumerate_sectors_within_radius(seed_center, radius_pc, EDGE_PC)
        # The wide outline holds every band of the shells this reaches.
        if abs(band_index) <= WIDE_GALAXY_HIGH_BAND
    }
    actual_addresses = {_address(row) for row in sectors}
    assert actual_addresses == expected_addresses
    assert len(sectors) == len(actual_addresses)


def test_random_neighborhoods_are_separate_and_whole(mysql_config, capsys):
    """GEN.97: `--neighborhoods 3` makes three starts, each at least twice the radius from the others."""
    import math
    _plan_wide_galaxy(mysql_config)
    radius_pc = 4.0
    _run_cli(["--neighborhoods", "3", "--max-shell", "60", "--radius-pc", str(radius_pc), "--num-systems", "1"]
             + _mysql_argv(mysql_config))
    output = capsys.readouterr()
    starts = [m for m in (output.out + output.err).splitlines() if m.startswith("Saved random starting sector")]
    assert len(starts) == 3
    sectors = _all_sectors(mysql_config)
    seeds = sorted(sectors, key=lambda row: row["id"])
    centers = [(r["center_x_pc"], r["center_y_pc"], r["center_z_pc"]) for r in sectors]
    # Every sector belongs to a neighborhood of one of the three starts.
    start_ids = []
    for message in starts:
        start_ids.append(int(message.split("sector_id=")[1].split(")")[0]))
    by_id = {r["id"]: (r["center_x_pc"], r["center_y_pc"], r["center_z_pc"]) for r in sectors}
    start_centers = [by_id[i] for i in start_ids]
    for a in range(3):
        for b in range(a + 1, 3):
            assert math.dist(start_centers[a], start_centers[b]) >= 2 * radius_pc
    assert all(min(math.dist(c, s) for s in start_centers) <= radius_pc + 1e-6 for c in centers)
    assert len(seeds) == len(sectors)


def test_avoid_filled_space_starts_in_empty_space_and_says_how_many_qualified(mysql_config, capsys):
    """GEN.186: with the flag, no start's neighborhood touches a sector that was already generated."""
    import math
    _plan_wide_galaxy(mysql_config)
    radius_pc = 4.0
    base = ["--max-shell", "60", "--radius-pc", str(radius_pc), "--num-systems", "1"] + _mysql_argv(mysql_config)
    _run_cli(["--neighborhoods", "2"] + base)
    before = _all_sectors(mysql_config)
    capsys.readouterr()
    _run_cli(["--avoid-filled-space"] + base)
    output = capsys.readouterr()
    text = output.out + output.err
    assert "Keeping away from filled space:" in text and "qualified" in text
    new = [r for r in _all_sectors(mysql_config) if r["id"] not in {b["id"] for b in before}]
    assert new
    old_centres = [(r["center_x_pc"], r["center_y_pc"], r["center_z_pc"]) for r in before]
    # The new sectors are a whole neighborhood of their own; none of them touches an old sector.
    for r in new:
        centre = (r["center_x_pc"], r["center_y_pc"], r["center_z_pc"])
        assert all(math.dist(centre, old) > 0 for old in old_centres)


def _occupy_small_outline(mysql_config, addresses):
    for address in sorted(addresses):
        run_galaxy.generate_and_save_sector_at(
            _fake_args_for_direct_generation(mysql_config), address,
            sector_position_pc(*address, EDGE_PC), EDGE_PC,
        )


def test_random_start_mode_retries_until_an_unoccupied_address_is_found(mysql_config, monkeypatch):
    # Occupy all but one of the outline's 13 cells directly -- the random
    # pick must keep retrying (not immediately fail) until it lands on the
    # one remaining unoccupied cell.
    _plan_only_the_core_and_shell_one_equator(mysql_config)
    free = (1, 0, 7)
    _occupy_small_outline(mysql_config, _SMALL_OUTLINE - {free})
    assert len(_all_sectors(mysql_config)) == 12

    _run_cli(["--max-shell", "1", "--radius-pc", "0.001", "--num-systems", "1"] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    addresses = {_address(row) for row in sectors}
    assert free in addresses
    assert len(sectors) == 13


def test_random_start_mode_gives_up_after_max_attempts_when_fully_occupied(mysql_config, monkeypatch):
    # Every cell of the outline already occupied -- random-start mode must
    # eventually give up with a clear error rather than looping forever or
    # crashing obscurely. A tiny attempts cap keeps this test fast.
    _plan_only_the_core_and_shell_one_equator(mysql_config)
    _occupy_small_outline(mysql_config, _SMALL_OUTLINE)

    monkeypatch.setattr(tuning, "RANDOM_START_MAX_PLACEMENT_ATTEMPTS", 5)

    with pytest.raises(SystemExit):
        _run_cli(["--max-shell", "1", "--radius-pc", "1.0", "--num-systems", "1"] + _mysql_argv(mysql_config))


def test_random_start_mode_respects_min_start_density(mysql_config, monkeypatch):
    """
    --min-start-density must reject an otherwise-qualifying seed address
    whose own real relative_density falls short of it, retrying until one
    that clears it is found. The threshold sits between the densest and
    sparsest cell of the outline (the arms make them differ).
    """
    _plan_only_the_core_and_shell_one_equator(mysql_config)

    densities = {
        address: relative_density(sector_position_pc(*address, EDGE_PC), _SKELETON_SHAPE)
        for address in _SMALL_OUTLINE
    }
    threshold = (max(densities.values()) + min(densities.values())) / 2
    assert min(densities.values()) < threshold < max(densities.values()), (
        "test setup needs the outline's cells to straddle the threshold"
    )

    _run_cli(
        ["--max-shell", "1", "--radius-pc", "0.001", "--min-start-density", str(threshold)]
        + _mysql_argv(mysql_config)
    )

    sectors = _all_sectors(mysql_config)
    seed = min(sectors, key=lambda row: row["id"])
    assert densities[_address(seed)] >= threshold


def test_random_start_mode_gives_up_when_min_start_density_unattainable(mysql_config, monkeypatch, capsys):
    """A --min-start-density no address in range can ever satisfy must
    still give up cleanly (not loop forever), with a message that calls
    out the density constraint specifically. That message goes through
    `log.error` (logged to the console, same as the fully-occupied case
    above) rather than `SystemExit`'s own args, so it's asserted against
    captured stdout instead of `pytest.raises(..., match=...)`."""
    _seed_skeleton(mysql_config, shells=_shells(1, 1))

    monkeypatch.setattr(tuning, "RANDOM_START_MAX_PLACEMENT_ATTEMPTS", 5)

    with pytest.raises(SystemExit):
        _run_cli(
            ["--max-shell", "1", "--radius-pc", "0.001", "--min-start-density", "1000000"]
            + _mysql_argv(mysql_config)
        )

    assert "min-start-density" in capsys.readouterr().out
    assert _all_sectors(mysql_config) == []


def _fake_args_for_direct_generation(mysql_config):
    """A minimal args namespace shaped like generate_cli.process_args()'s own
    output, for directly calling generate_and_save_sector_at in a test
    without going through the CLI -- mirrors run_galaxy._default_generation_args
    but pointed at the fixture's own throwaway database, with num_systems
    pinned to 1 (this file's own "keep it fast" convention -- see the
    module docstring)."""
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.density = None
    args.num_systems = 1
    return args


# ---------------------------------------------------------------------------
# End-to-end against a *real* 'planetgen plan' skeleton -- every test
# above either hand-seeds a skeleton via _seed_skeleton (bypassing
# build_shell_extents' own scan entirely) or passes an explicit
# --num-systems/--density that makes _BatchDensity.resolve a no-op (see
# its own docstring), so none of them exercise the actual "plan, then
# galaxy with neither flag given" workflow an operator runs -- exactly the
# combination that once shipped a regression saving real (0-system,
# 0-phenomenon) sectors. Since GEN.76 every cell inside the outline is
# generated however sparse, and generate_sector's fallback keeps each one
# from being empty. These tests close that gap: a real skeleton, real CLI
# generation, no bypass.
# ---------------------------------------------------------------------------

REFERENCE_SHAPE_ARGV = [
    "--disk-scale-length", "2600pc", "--flatness", "8.667",
    "--bulge-size", "0.6077", "--bulge-share", "0.595", "--arm-count", "2", "--pitch-angle-deg", "15",
]
"""list: The Milky Way's published shape as `planetgen plan` options. A plan draws its shape from the seed (GEN.401);
a test whose counts were written for one galaxy's density pins this one."""

_PLAN_SHAPE_ARGV = [
    "--disk-scale-length", "40pc", "--flatness", "3.333",
    "--bulge-size", "0.25", "--bulge-share", "0.0245",
    "--arm-count", "2", "--pitch-angle-deg", "15", "--arm-contrast", "0.4",
]
"""list: The same shape `_SKELETON_SHAPE` above is built from, expressed as
`planetgen plan` CLI args instead -- a small, fast-to-scan toy galaxy
(its real, found outer edge lands a few dozen shells out, not the ~3,800 a
real Milky-Way-scale build reaches) that still has genuine bulge/disk/arm
structure, not a degenerate single-shell case."""


def _build_real_skeleton(mysql_config, extra_argv=()):
    """
    Runs `planetgen plan`'s own skeleton command (`run_plan._skeleton_on_queue`) -- the real scan
    (`galaxySkeleton.build_shell_extents`), not the hand-seeded
    `_seed_skeleton` shortcut above -- against the fixture's throwaway
    database, and returns its summary dict (`outer_shell_index`, etc.).
    Calls `build_skeleton` directly rather than going through `main()`/
    `sys.argv` purely to get that return value back without parsing
    stdout; it still does the same real persisting
    (`_db.replace_galaxy_shells`/`_db.save_galaxy_shape`) `run_plan`
    itself does.
    """
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args(_PLAN_SHAPE_ARGV + list(extra_argv) + _mysql_argv(mysql_config))
    generate_cli.validate_plan_args(args, parser)
    # The plan draws the galaxy's naming key into the control database (GEN.70); the fixture's database holds the
    # control schema too, so the draw works and no "control database can't be opened" line is logged (TEST.191).
    # A toy of 40 pc scale length forms half a million young open clusters of 50 solar masses and up (GEN.244: the count
    # follows the star formation rate density), which every plan test would draw, store and load; only the heavy end is kept.
    with mock.patch.dict(os.environ, {_db.CONTROL_DB_ENV_VAR: mysql_config.database}), \
            mock.patch.object(tuning, "OPEN_CLUSTER_MASS_RANGE_SOL", (2.0e4, 1.0e6)):
        _db.get_control_connection(mysql_config, ensure_schema=True).close()
        return run_plan._skeleton_on_queue(args)


def test_a_real_plan_draws_the_galaxys_naming_key_into_the_control_database(mysql_config, capsys):
    """TEST.191: the plan tests have a control database to draw into, so the naming path runs for real."""
    _build_real_skeleton(mysql_config)
    assert "No naming key drawn" not in capsys.readouterr().out
    conn = _db.get_control_connection(mysql_config)
    try:
        rows = conn.execute("SELECT naming_key FROM galaxy_naming WHERE database_name = ?",
                            (mysql_config.database,)).fetchall()
    finally:
        conn.close()
    assert len(rows) == 1 and rows[0]["naming_key"]


def _sector_system_counts(mysql_config, sector_ids):
    """`{sector_id: system_count}` for every id in `sector_ids` (0 for one
    with no systems at all)."""
    if not sector_ids:
        return {}
    conn = _db.get_connection(mysql_config)
    try:
        placeholders = ", ".join(["?"] * len(sector_ids))
        rows = conn.execute(
            f"SELECT sector_id, COUNT(*) AS cnt FROM star_systems "
            f"WHERE sector_id IN ({placeholders}) GROUP BY sector_id",
            tuple(sector_ids),
        ).fetchall()
    finally:
        conn.close()
    counts = {sector_id: 0 for sector_id in sector_ids}
    for row in rows:
        counts[row["sector_id"]] = row["cnt"]
    return counts


DENSE_CORE_SEED = "96901074FC015B42F2E80F1434609687"
"""str: A galaxy seed whose toy galaxy holds cells at 1,000 times the local density within 55 shells. A galaxy drawn at
random is too light for one in three draws (its densest cell stays under that), and the run then finds no start
(TEST.225)."""


@pytest.mark.parametrize("start, plan, radius_ly", [([], [], 25.0), (["--min-start-density", "1000"], ["--seed", DENSE_CORE_SEED], 18.0)],
                         ids=["random", "dense_core"])
def test_random_start_neighborhood_matches_the_real_skeleton_plan(mysql_config, monkeypatch, start, plan, radius_ly):
    """
    The standard workflow this project's own docs describe -- 'plan' once,
    then 'galaxy' with no flags: pick a random location, generate every
    not-yet-generated sector out to `--radius-pc` (the real default is
    `tuning.DEFAULT_GENERATE_RADIUS_PC`, 12 pc;
    trimmed to 25 ly here so this test runs in a reasonable time) -- run
    for real, against a real skeleton, with neither `--density` nor
    `--num-systems` given so every sector's own system count is driven
    end-to-end by its real galaxy-frame position. planet-less systems are forced (`planets=False`) so
    each system skips its own planet/moon tree -- system *count* per
    sector (what this test actually checks) doesn't depend on that, and
    skipping it keeps this test's runtime independent of how dense a
    shell the random draw happens to land in.

    Confirms two things no test above does:
    1. Every candidate cell the run actually saved, and every one it
       didn't, agrees with an independent recomputation -- done here,
       against the real stored skeleton -- of whether that exact cell
       is inside its shell's stored extent. Inside, every cell is saved,
       however sparse (GEN.76); outside, none is.
    2. The aggregate system count actually generated across the
       neighborhood's denser sectors (one or more predicted stars) is in
       the right statistical ballpark of what the plan's own density
       predicted at those same positions -- not a flat,
       position-independent count (`_default_generation_args`'s own
       num_systems=10 default, in particular, would badly fail this).
       A cell the plan predicts above the cell star ceiling is not filled
       with every system (GEN.281): it is taken down to its luminosity cap
       and no further, and the scatter holds the rest, which this run does not
       draw, so such a cell is left out of the sum and checked against the
       ceiling instead. A start in the dense core (`--min-start-density`)
       lands in such cells (TEST.225).
       Sparser sectors are left out of the sum: one with nothing drawn
       at all gets one fallback system, which would swamp their tiny
       predictions.
    """
    radius_pc = ly_to_pc(radius_ly)
    summary = _build_real_skeleton(mysql_config, plan)
    assert summary["outer_shell_index"] > 0, "the toy shape's own skeleton should find real content"

    _run_cli([
        "--max-shell", "55", "--radius-pc", str(radius_pc),
    ] + start + _mysql_argv(mysql_config), planets=False)

    sectors = _all_sectors(mysql_config)
    assert sectors, "random-start mode should have generated at least the seed sector"
    seed = min(sectors, key=lambda row: row["id"])
    center = (seed["center_x_pc"], seed["center_y_pc"], seed["center_z_pc"])

    candidates = list(enumerate_sectors_within_radius(center, radius_pc, EDGE_PC))
    assert len(candidates) > 1, "the 25 ly neighborhood should reach beyond just the seed sector itself"

    by_address = {_address(row): row for row in sectors}

    conn = _db.get_connection(mysql_config)
    try:
        skeleton = _db.get_galaxy_shape(conn)
        high_bands = {shell: high for shell, high, _cells in _db.get_galaxy_shells(conn)}

        def inside_outline(shell_index, band_index):
            return shell_index in high_bands and abs(band_index) <= high_bands[shell_index]

        expected_total = 0.0
        generated_sector_ids = []
        capped_sector_ids = []
        for shell_index, band_index, wedge_index, x, y, z, _dist in candidates:
            address = (shell_index, band_index, wedge_index)
            saved = address in by_address
            inside = inside_outline(shell_index, band_index)
            assert saved == inside, (
                f"address {address} was {'saved' if saved else 'skipped'} by 'planetgen galaxy', but "
                f"it is {'inside' if inside else 'outside'} the real stored skeleton's outline -- "
                f"generation has drifted from the plan it's supposed to follow."
            )
            predicted = predicted_star_count(
                (x, y, z), skeleton.shape, skeleton.expected_system_count_at_density_1,
            )
            if inside and predicted >= 1.0:
                if galaxy_resolution.cell_resolution(
                        galaxy_profile.sector_profile(address, skeleton), skeleton.cell_star_ceiling,
                        skeleton.fill_luminosity_cap_sol).capped:
                    capped_sector_ids.append(by_address[address]["id"])
                    continue
                generated_sector_ids.append(by_address[address]["id"])
                expected_total += predicted
        ceiling = skeleton.cell_star_ceiling
    finally:
        conn.close()

    # GEN.281: a cell over the ceiling is taken down to its cap and holds no more than the ceiling.
    for sector_id, count in _sector_system_counts(mysql_config, capped_sector_ids).items():
        assert count <= ceiling, f"sector {sector_id} holds {count} systems, over the cell star ceiling {ceiling}"
    counts = _sector_system_counts(mysql_config, generated_sector_ids)
    actual_total = sum(counts.values())

    if expected_total >= 5.0:
        # The sum of independent Poisson draws is itself ~Poisson(expected_total)
        # -- a generous 3x-either-way band comfortably absorbs real
        # sampling noise while still catching a real "density isn't
        # driving this any more" regression (a flat count, or the wrong
        # multiplier, would miss by far more than 3x here).
        assert expected_total / 3.0 <= actual_total <= expected_total * 3.0, (
            f"generated {actual_total} systems across this neighborhood, but the real skeleton "
            f"predicted {expected_total:.1f} -- generation doesn't look density-driven any more."
        )


def test_a_shell_beyond_the_real_skeletons_outer_edge_generates_from_its_profiles(mysql_config):
    """
    Deterministic companion to the neighborhood test above (that one's
    outcome depends on where the random draw lands; this one doesn't): a
    shell and band named well beyond the real skeleton's own discovered outer edge
    are inside the grid's reach, so (GEN.235) their sectors are generated from their
    profiles, with or without an explicit --num-systems, and saved though they come
    out empty.
    """
    summary = _build_real_skeleton(mysql_config)
    beyond_edge_shell = summary["outer_shell_index"] + 10

    for extra, limit in (([], "3"), (["--num-systems", "1"], "5")):
        _run_cli([
            "--shell", str(beyond_edge_shell), "--band", "0", "--limit", limit,
        ] + extra + _mysql_argv(mysql_config), planets=False)

    sectors = _all_sectors(mysql_config)
    assert len(sectors) == 3 + 5       # each run takes its limit of the not-yet-generated sectors
    assert {row["shell_index"] for row in sectors} == {beyond_edge_shell}


# ---------------------------------------------------------------------------
# Guaranteed non-empty: a qualifying sector's own system count and each
# phenomenon type's own count are independent Poisson draws, so all of them
# landing on zero simultaneously is a real, expected outcome at low means
# (e.g. ~13% at mean 2, and only gets more likely approaching the mean-~1
# qualification threshold itself) -- generate_sector forces exactly one
# system onto an otherwise-completely-empty, --density-driven sector rather
# than leave it with nothing in it at all. No database needed -- these call
# generate_sector directly, not through the CLI/mysql_config fixture.
# ---------------------------------------------------------------------------

def test_generate_sector_is_never_left_with_nothing_in_it():
    """
    Runs generate_sector directly, many times, at a mean of exactly ~1
    system/sector (right at the qualification threshold, where an empty
    Poisson draw is common -- P(0) = 1/e = 36.8%) -- across 200 draws, an
    unpatched version would produce a genuinely empty sector far more
    often than this test could plausibly miss by chance alone.
    """
    from planetgen.galaxy.sector import SpaceSector

    e_value = SpaceSector(name="calibration").expected_system_count()
    args = run_galaxy._default_generation_args()
    args.density = 1.0 / e_value
    args.num_systems = None

    for _ in range(200):
        _sector_name, sector = run_sector.generate_sector(args)
        assert sector.entries or sector.phenomena, (
            "generate_sector produced a sector with nothing in it at all"
        )


def test_generate_sector_does_not_force_content_onto_an_explicit_zero():
    """The guaranteed-non-empty fallback only overrides a --density-driven
    zero -- an explicit --num-systems 0 is a deliberate request it must
    never second-guess."""
    args = run_galaxy._default_generation_args()
    args.density = None
    args.num_systems = 0

    _sector_name, sector = run_sector.generate_sector(args)
    assert sector.entries == []
    assert sector.phenomena == []


# ---------------------------------------------------------------------------
# Sector-view bounds: every star system AND every placed exotic phenomenon,
# across a real, multi-shell test galaxy, must land inside the actual
# spherical cell its own sector occupies in the galaxy frame -- its shell's
# radial bounds, its band's colatitude cone and its wedge's longitude --
# not just a same-named number compared to a constant.
#
# Star systems store `position_x/y/z_mpc` in the sector's local frame
# (`galaxy.frame.sector_orientation`: +X radial, +Y tangential, +Z
# north), which is the frame `SectorCell.contains` reads. Phenomena already
# carry an absolute `center_x/y/z_pc` (`_db.insert_sector` computes it via
# `_galaxy_placement_from_sector_offset`), turned into the local frame the
# same way before the check.
# ---------------------------------------------------------------------------

_CHECKABLE_PHENOMENON_TABLES = (
    ("black_holes", "black hole"),
    ("neutron_stars", "neutron star"),
    ("nebulae", "nebula"),
    ("asteroid_fields", "asteroid field"),
    ("supernova_remnants", "supernova remnant"),
    ("rogue_planets", "rogue planet"),
    ("interstellar_comets", "interstellar comet"),
)
"""tuple: `(table_name, label)` for every phenomenon type that gets its own
real `center_x/y/z_pc` -- the only ones a coordinate-bounds check is even
possible for."""

_BOUNDS_SHRINK = 1.0 - 1e-7
"""float: A sector-local point is pulled this much toward the sector center
before the cell-membership check below, purely to absorb rounding (positions
are stored to a thousandth of a parsec), not to paper over a real violation."""

_BOUNDS_TEST_BANDS = [(1, 0), (2, 0), (3, 0), (4, 0), (2, 1), (3, -1)]
"""list: `(shell, band)` batches the bounds test generates -- the first four
shells' equator bands plus two off-equator ones."""


def _sector_cell_contexts(mysql_config):
    """
    For every galaxy-placed sector: its absolute center, its cell and its
    `(shell, band, wedge)` address.
    """
    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT id, shell_index, band_index, wedge_index, "
            "center_x_pc, center_y_pc, center_z_pc FROM sectors "
            "WHERE center_x_pc IS NOT NULL"
        ).fetchall()
    finally:
        conn.close()

    contexts = {}
    for row in rows:
        center_pc = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        address = _address(row)
        contexts[row["id"]] = {
            "center_pc": center_pc,
            "cell": SectorCell.for_address(*address, EDGE_PC),
            "address": address,
        }
    return contexts


def _bounds_violation(sector_ctx, local_position_pc, label):
    """`None` if `local_position_pc` (the sector's local frame, parsecs)
    falls inside `sector_ctx`'s real cell; otherwise a description of which
    cell it was meant to be in."""
    shell, band, wedge = sector_ctx["address"]
    pulled = tuple(value * _BOUNDS_SHRINK for value in local_position_pc)
    if sector_ctx["cell"].contains(pulled):
        return None
    return f"{label} (shell {shell}, band {band}, wedge {wedge}): {tuple(local_position_pc)} pc is outside its cell"


def test_stars_and_phenomena_fit_within_their_sectors_real_cells(mysql_config, monkeypatch):
    """
    Generates a real test galaxy (the first shells' equator bands plus two
    off-equator ones, see `_BOUNDS_TEST_BANDS`) through the real
    `planetgen galaxy --shell --band` pipeline, then checks every star system's
    and every checkable phenomenon's position against the real spherical cell
    its own sector occupies.

    Real astrophysical rates make a black hole or nebula genuinely rare
    per sector, so `_sample_poisson_count` is monkeypatched to always
    return 1 (for every positive mean), forcing every one of the seven
    phenomenon types to appear in every sector -- deterministic, full
    coverage of every checkable type's own placement path.

    `planets=False` keeps each system's own generation cheap (position, not
    planet/moon content, is what this test cares about).
    """
    worker_patches.patch_everywhere(monkeypatch, run_sector, "_sample_poisson_count",
                                    "tests.test_galaxy_gen:poisson_always_one")
    _plan_wide_galaxy(mysql_config)

    for shell_index, band_index in _BOUNDS_TEST_BANDS:
        _run_cli(
            ["--shell", str(shell_index), "--band", str(band_index), "--num-systems", "3"]
            + _mysql_argv(mysql_config), planets=False
        )

    sectors = _all_sectors(mysql_config)
    expected_sector_count = sum(wedge_count(shell, band) for shell, band in _BOUNDS_TEST_BANDS)
    assert len(sectors) == expected_sector_count

    contexts = _sector_cell_contexts(mysql_config)
    assert len(contexts) == expected_sector_count

    conn = _db.get_connection(mysql_config)
    try:
        star_rows = conn.execute(
            "SELECT sector_id, id AS system_id, position_x_mpc, position_y_mpc, position_z_mpc "
            "FROM star_systems WHERE position_x_mpc IS NOT NULL"
        ).fetchall()

        phenomenon_rows = []
        # A supernova remnant's core is filed under its remnant's sector but
        # has drifted off by its birth kick, often out of that cell.
        core_filter = {
            "black_holes": " AND id NOT IN (SELECT compact_remnant_black_hole_id FROM supernova_remnants"
                           " WHERE compact_remnant_black_hole_id IS NOT NULL)",
            "neutron_stars": " AND id NOT IN (SELECT compact_remnant_neutron_star_id FROM supernova_remnants"
                             " WHERE compact_remnant_neutron_star_id IS NOT NULL)",
            # A field cloud (the dark family) is filed under the sector holding its centroid, and a cloud
            # tens of light years across has its centre in a neighbouring cell whenever it sits near a face.
            "nebulae": " AND nebula_type <> 'dark'",
        }
        for table, label in _CHECKABLE_PHENOMENON_TABLES:
            rows = conn.execute(
                f"SELECT sector_id, id, center_x_pc, center_y_pc, center_z_pc FROM {table} "
                f"WHERE center_x_pc IS NOT NULL{core_filter.get(table, '')}"
            ).fetchall()
            phenomenon_rows.extend((label, row) for row in rows)
    finally:
        conn.close()

    assert star_rows, "expected at least some star systems"
    seen_labels = {label for label, _row in phenomenon_rows}
    # Isolated asteroid fields have a zero research rate (they disperse), so
    # sector generation never makes one; placement is still checked for the rest.
    assert seen_labels == {label for _table, label in _CHECKABLE_PHENOMENON_TABLES} - {"asteroid field"}

    violations = []

    for row in star_rows:
        ctx = contexts[row["sector_id"]]
        local = (mpc_to_pc(row["position_x_mpc"]), mpc_to_pc(row["position_y_mpc"]),
                 mpc_to_pc(row["position_z_mpc"]))
        violation = _bounds_violation(ctx, local, f"star system {row['system_id']}")
        if violation:
            violations.append(violation)

    for label, row in phenomenon_rows:
        ctx = contexts[row["sector_id"]]
        absolute_pc = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        violation = _bounds_violation(ctx, galaxy_to_local_pc(ctx["center_pc"], absolute_pc),
                                      f"{label} {row['id']}")
        if violation:
            violations.append(violation)

    total_checked = len(star_rows) + len(phenomenon_rows)
    assert not violations, (
        f"{len(violations)} of {total_checked} generated stars/phenomena fell outside their own "
        f"sector's real cell:\n" + "\n".join(violations[:50])
    )


def test_sector_map_draws_no_point_object_outside_its_own_sector(mysql_config, monkeypatch):
    """
    MAP.45: what the Sector Map draws for a sector (`query.phenomena_near_sector`,
    rendered by `planetgen/web/maps/starmap.py`) keeps every point-like object -- a
    rogue planet, comet, black hole or neutron star -- inside that sector's
    own cell. A neighbor's rogue planets used to be pulled in by a sphere
    sized for the old cube and drawn outside the wireframe. Only a cloud
    may come from a neighbor, and it is marked `neighbor` on the map.

    Shell 1's equator band and the first three wedges of shell 2's sit side
    by side, so each sector has generated neighbors whose rogue planets would
    show up if the search still took them.
    """
    from planetgen.web.maps.starmap import map_scene_data

    _plan_wide_galaxy(mysql_config)
    _run_cli(["--shell", "1", "--band", "0", "--num-systems", "10"] + _mysql_argv(mysql_config), planets=False)
    _run_cli(["--shell", "2", "--band", "0", "--limit", "3", "--num-systems", "10"] + _mysql_argv(mysql_config),
             planets=False)

    contexts = _sector_cell_contexts(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        drawn = {sector_id: query.phenomena_near_sector(conn, sector_id) for sector_id in contexts}
        edges = {row["id"]: row["edge_mpc"] for row in conn.execute("SELECT id, edge_mpc FROM sectors").fetchall()}
    finally:
        conn.close()

    violations = []
    rogue_planets = 0
    for sector_id, phenomena in drawn.items():
        ctx = contexts[sector_id]
        for entry in phenomena:
            if not entry["home"]:
                assert entry["radius_ly"] > 0, f"{entry['type']} {entry['id']} is a neighbor's point object"
                continue
            rogue_planets += entry["type"] == "rogue_planet"
            if entry["type"] in ("black_hole", "neutron_star"):
                continue  # a supernova core drifts off by its kick; the nucleus sits on the axis
            # The offsets are along the galaxy's axes; the cell's are the sector's own.
            point = tuple(center + ly_to_pc(entry[f"offset_{axis}_ly"])
                          for center, axis in zip(ctx["center_pc"], "xyz"))
            local = galaxy_to_local_pc(ctx["center_pc"], point)
            violation = _bounds_violation(ctx, local, f"{entry['type']} {entry['id']}")
            if violation:
                violations.append(violation)

        scene = map_scene_data(lambda name, **params: "#", edges[sector_id], ctx["address"], ctx["center_pc"],
                               [], phenomena=phenomena)
        for entry, phenomenon in zip(scene["clouds"], phenomena):
            assert entry.get("neighbor", False) == (not phenomenon["home"])

    assert rogue_planets, "expected the sectors to generate some rogue planets"
    assert not violations, "\n".join(violations[:50])


def test_sector_cell_used_for_generation_matches_the_grid():
    """`SectorCell.for_address` in light-years (what generation samples in)
    is the same cell as the grid's parsec bounds, scaled."""
    for shell, band, wedge in ((0, 0, 0), (3, 1, 2), (3, -4, 0), (40, 0, 5), (40, 12, 30)):
        cell = SectorCell.for_address(shell, band, wedge, tuning.DEFAULT_SECTOR_EDGE_LY)
        r_lo, r_hi = spherical.shell_bounds_pc(shell, EDGE_PC)
        assert ly_to_pc(cell.r_inner) == pytest.approx(r_lo)
        assert ly_to_pc(cell.r_outer) == pytest.approx(r_hi)
        north, south = spherical.band_colatitude_bounds(shell, band)
        assert (cell.colat_north, cell.colat_south) == pytest.approx((north, south))
        start, end = spherical.wedge_longitude_bounds(shell, band, wedge)
        assert (cell.lon_start, cell.lon_end) == pytest.approx((start, end))
        assert tuple(ly_to_pc(value) for value in cell.center) == pytest.approx(
            sector_position_pc(shell, band, wedge, EDGE_PC))


# ---------------------------------------------------------------------------
# Validation before generation: every galaxy path checks an address against
# the stored outline first, and generates nothing outside it.
# ---------------------------------------------------------------------------

def test_galaxy_refuses_to_run_before_a_plan(mysql_config):
    with pytest.raises(SystemExit):
        _run_cli(["--shell", "1", "--band", "0", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []


_SMALL_PLAN = [(0, 0), (1, 1), (2, 1)]
"""list: The outline the checks below plan: the core, shell 1 and shell 2, the last two reaching bands -1 to 1."""


@pytest.mark.parametrize("argv", [
    ["--shell", "6", "--band", "0"],               # past the reach (shell 5)
    ["--shell", "2", "--band", "9"],               # a band shell 2 does not hold (-3 to 3)
    ["--shell", "6", "--band", "0", "--wedge", "0"],
    ["--shell", "2", "--band", "9", "--wedge", "0"],
])
@pytest.mark.parametrize("strict", [[], ["--strict"]])
def test_addresses_beyond_the_reach_are_refused_before_generating(mysql_config, argv, strict):
    # GEN.235: the reach is the only refusal, with or without --strict.
    _seed_skeleton(mysql_config, shells=_SMALL_PLAN, reach_shell_index=5)
    with pytest.raises(SystemExit):
        _run_cli(argv + ["--num-systems", "1"] + strict + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []


@pytest.mark.parametrize("argv, expected", [
    (["--shell", "3", "--band", "0", "--wedge", "0"], {(3, 0, 0)}),
    (["--shell", "2", "--band", "2", "--wedge", "0"], {(2, 2, 0)}),
    (["--shell", "2", "--band", "2"], _band_cells(2, 2)),
])
def test_an_address_named_outside_the_outline_is_generated_from_its_profile(mysql_config, capsys, argv, expected):
    # GEN.235: an address the user named is generated anywhere in the reach, with no warning.
    _seed_skeleton(mysql_config, shells=_SMALL_PLAN)
    _run_cli(argv + ["--strict"] + _mysql_argv(mysql_config), planets=False)
    assert {_address(row) for row in _all_sectors(mysql_config)} == expected
    assert "WARNING:" not in capsys.readouterr().out


def test_neighborhood_is_trimmed_to_the_outline_even_with_explicit_num_systems(mysql_config):
    # A galaxy of four shells and the equator band: a 12 pc sphere around
    # shell 2 reaches well past it on every side, and only the cells inside
    # get generated.
    outline = [(0, 0), (1, 0), (2, 0), (3, 0)]
    _seed_skeleton(mysql_config, shells=outline)
    run_galaxy.generate_and_save_sector_at(
        _fake_args_for_direct_generation(mysql_config), (2, 0, 0), sector_position_pc(2, 0, 0, EDGE_PC), EDGE_PC,
    )
    center_id = _all_sectors(mysql_config)[0]["id"]
    center = sector_position_pc(2, 0, 0, EDGE_PC)

    _run_cli(["--center-sector", str(center_id), "--radius-pc", "12", "--num-systems", "1"]
             + _mysql_argv(mysql_config))

    limits = dict(outline)
    sphere = {c[:3] for c in enumerate_sectors_within_radius(center, 12.0, EDGE_PC)}
    inside = {a for a in sphere if a[0] in limits and abs(a[1]) <= limits[a[0]]}
    assert inside != sphere, "test setup: the sphere must reach past the outline"
    assert {_address(row) for row in _all_sectors(mysql_config)} == inside


def test_neighborhood_around_a_sector_outside_the_outline_takes_the_reach_and_says_how_big_it_is(mysql_config, capsys):
    # GEN.235: the pull from a cell inside the outline stays in it (the test above); the pull from a halo sector
    # takes every cell of the reach within the radius, each from its own profile, and reports the size first.
    _seed_skeleton(mysql_config, shells=_SMALL_PLAN)
    far = (40, 0, 0)
    center = sector_position_pc(*far, EDGE_PC)
    run_galaxy.generate_and_save_sector_at(
        _fake_args_for_direct_generation(mysql_config), far, center, EDGE_PC,
    )
    center_id = _all_sectors(mysql_config)[0]["id"]

    _run_cli(["--center-sector", str(center_id), "--radius-pc", "12", "--num-systems", "1", "--strict"]
             + _mysql_argv(mysql_config))
    sphere = {c[:3] for c in enumerate_sectors_within_radius(center, 12.0, EDGE_PC)}
    assert {_address(row) for row in _all_sectors(mysql_config)} == sphere
    assert len(sphere) > 20
    out = capsys.readouterr().out
    assert f"{len(sphere):,} sector(s) lie within 12 pc of" in out and "outside the galaxy's outline" in out
    assert "systems in all" in out
    # The API job's version takes the same cells (nothing left to add, the sphere is generated).
    counts = run_galaxy.generate_sector_neighborhood(center_id, radius_ly=12.0 / ly_to_pc(1.0), config=mysql_config)
    assert counts["outside_galaxy"] == 0


def test_random_start_and_its_neighborhood_stay_inside_a_real_plan(mysql_config):
    summary = _build_real_skeleton(mysql_config)
    assert summary["shell_count"] > 0
    _run_cli(["--radius-pc", "10", "--num-systems", "1"] + _mysql_argv(mysql_config), planets=False)

    conn = _db.get_connection(mysql_config)
    try:
        bounds = _db.get_galaxy_bounds(conn)
    finally:
        conn.close()
    sectors = _all_sectors(mysql_config)
    assert sectors
    for row in sectors:
        assert bounds.qualifies(row["shell_index"], row["band_index"]), _address(row)


def test_wedge_mode_reports_the_sector_it_made_with_its_stats(mysql_config, capsys):
    """ADM.47: the one sector named outright is reported like every sector of a batch."""
    _seed_skeleton(mysql_config, shells=_shells(3, 3))
    _run_cli(["--shell", "2", "--band", "2", "--wedge", "1"] + _mysql_argv(mysql_config))
    out = capsys.readouterr().out
    assert "Saved sector" in out
    assert "Systems:" in out and "Totals:" in out and "Star density: actual" in out


def test_a_batch_says_how_many_sectors_it_runs_and_reports_each_with_its_stats(mysql_config, capsys):
    _seed_skeleton(mysql_config, shells=_shells(3, 3))
    _run_cli(["--shell", "2", "--band", "2", "--limit", "2"] + _mysql_argv(mysql_config))
    out = capsys.readouterr().out
    assert "Generating 2 sector(s) with" in out
    assert out.count("Saved sector") == 2 and out.count("Totals:") == 2 and out.count("Star density: actual") == 2


def test_the_freshness_state_does_not_build_the_galaxy_shape(mysql_config, monkeypatch):
    """A page view's freshness check (`galaxy_content_state`) read the whole shape, whose component field works out a mean
    system mass for every stored cluster: minutes on a big galaxy, past the web daemon's 60-second limit, for every page.
    It reads what the shape is made from instead, and still moves when the stored shape does."""
    _seed_skeleton(mysql_config, outer_shell_index=42, e_value=2.5)

    def not_here(*args, **kwargs):
        raise AssertionError("the freshness state built the galaxy shape")

    monkeypatch.setattr(query, "galaxy_density_shape", not_here)
    conn = _db.get_connection(mysql_config)
    try:
        first = query.galaxy_content_state(conn)
        assert query.galaxy_content_state(conn) == first
        conn.execute("UPDATE galaxy_shape SET outer_shell_index = 43 WHERE id = 1")
        conn.commit()
        changed = query.galaxy_content_state(conn)
    finally:
        conn.close()
    assert changed["base"] != first["base"] and changed["tile_base"] != first["tile_base"]
