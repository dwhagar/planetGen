# tests/test_galaxy_gen.py

"""
End-to-end tests for `generate.py`'s `galaxy` subcommand -- the
galaxy-scale sector generation logic (`docs/design/galaxy-coordinate-system.md`
section 8, `docs/TODO.md` Phase 4). `test_galaxy_geometry.py` already
covers `galaxyGeometry`'s pure functions (shell math, the enumeration
primitive) in isolation; this file is the missing piece -- actually
running the `galaxy` subcommand's two CLI modes (`--shell`,
`--center-sector`) against a real, throwaway MySQL database (see
`conftest.py`'s `mysql_config` fixture) and inspecting the resulting
`sectors` rows.

Runs the real CLI entry point (`main()`, via `sys.argv`) rather than
calling internals directly -- this repo's existing convention for testing
a CLI end-to-end (as opposed to `test_sector_gen.py`, which tests
`generate_sector_name` as a plain function because it has no database/CLI
side effects to speak of).

`--num-systems 1` keeps every generated sector's own system-generation
cost to a minimum -- these tests are about sector *placement*/persistence
(shell/slot addresses, positions, occupied-slot skipping), not system
content, which is already covered extensively elsewhere
(`test_systems.py`, `test_db_persistence.py`).

Every test here takes the `mysql_config` fixture -- they're skipped, not
failed, when no MySQL test server is configured/reachable.
"""

import argparse
import math
import sys

import pymysql
import pytest

import generate as galaxyGen
import generate as sectorGen
import queryDb
from stellarObjects import _db, program_constants
from stellarObjects.galaxyDensity import build_galaxy_shape, predicted_star_count, relative_density
from stellarObjects.galaxyGeometry import (
    enumerate_sectors_within_radius, sector_position_pc, shell_sector_count,
)
from stellarObjects.sectorGeometry import cube_orientation
from stellarObjects.utils import ly_to_pc, mpc_to_pc

EDGE_PC = ly_to_pc(program_constants.DEFAULT_SECTOR_EDGE_LY)

# A small, fast-to-evaluate toy shape -- same parameters test_galaxy_skeleton.py
# uses, deep enough in its own bulge at shell 0 that every slot there
# qualifies under any reasonable threshold, with a real, findable edge
# a few hundred shells out (not real Milky-Way scale, which would make
# these tests slow for no benefit -- ensure_sector_generated's own logic
# doesn't care about scale).
_SKELETON_SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0,
    disk_scale_height_pc=12.0,
    bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15),
    arm_amplitude=0.4,
)


def _seed_skeleton(mysql_config, shape=_SKELETON_SHAPE, outer_shell_index=999, e_value=1.0, bands=()):
    """
    Directly writes a `galaxy_shape`/`galaxy_shell_band` skeleton, without
    running `generate.py plan`'s own scan -- these tests exercise
    `ensure_sector_generated`'s own logic against a skeleton it can
    already read, not `find_shell_bands`'s search (covered separately in
    `test_galaxy_skeleton.py`).

    Args:
        mysql_config (MySQLConfig): The fixture's throwaway database.
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        outer_shell_index (int): See `_db.save_galaxy_shape`.
        e_value (float): `expected_system_count_at_density_1`.
        bands (iterable): `(shell_index, band_index, slot_index_min,
            slot_index_max)` tuples -- see `_db.replace_galaxy_shell_bands`.
    """
    _db.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=outer_shell_index,
                           expected_system_count_at_density_1=e_value, config=mysql_config)
    _db.replace_galaxy_shell_bands(list(bands), config=mysql_config)


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


def _run_cli(argv):
    """Runs `generate.main()` with the `galaxy` subcommand and the given
    argv (excluding argv[0]), restoring `sys.argv` afterward."""
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy"] + argv
        galaxyGen.main()
    finally:
        sys.argv = old_argv


def _run_sector_gen_cli(argv):
    """Runs `generate.main()` with the `sector` subcommand the same way --
    used to seed an "unplaced" (no galaxy position) sector for
    `test_center_sector_without_galaxy_position_is_rejected`."""
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "sector"] + argv
        sectorGen.main()
    finally:
        sys.argv = old_argv


def _all_sectors(mysql_config):
    """Returns every `sectors` row. `get_connection`'s own `_ensure_schema`
    creates the table if a rejected run never got as far as calling it
    (`run_shell_batch`'s `LARGE_SHELL_WARNING_THRESHOLD` guard raises
    before ever calling `_db.get_connection`), so this never errors on a
    fresh, empty database -- it just returns an empty list."""
    conn = _db.get_connection(mysql_config)
    try:
        return conn.execute(
            "SELECT id, name, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc, "
            "shell_index, shell_slot_index FROM sectors"
        ).fetchall()
    finally:
        conn.close()


def test_galaxy_density_shape_returns_none_before_a_skeleton_is_built(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        assert queryDb.galaxy_density_shape(conn) is None
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
        shape_dict = queryDb.galaxy_density_shape(conn)
    finally:
        conn.close()

    assert shape_dict is not None
    for field in _SKELETON_SHAPE._fields:
        assert shape_dict[field] == pytest.approx(getattr(_SKELETON_SHAPE, field))
    assert shape_dict["edge_pc"] == pytest.approx(EDGE_PC)
    assert shape_dict["outer_shell_index"] == 42
    assert shape_dict["expected_system_count_at_density_1"] == pytest.approx(2.5)


def test_shell_batch_mode_generates_every_slot_and_is_idempotent(mysql_config):
    # Shell 0 holds exactly 3 slots -- small enough to run under
    # LARGE_SHELL_WARNING_THRESHOLD with no --limit/--yes needed.
    assert shell_sector_count(0) == 3
    _run_cli(["--shell", "0", "--num-systems", "1"] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    assert len(sectors) == 3
    addresses = {(row["shell_index"], row["shell_slot_index"]) for row in sectors}
    assert addresses == {(0, 0), (0, 1), (0, 2)}
    for row in sectors:
        assert row["center_x_pc"] is not None
        assert row["galactic_radius_pc"] == pytest.approx(
            (row["center_x_pc"] ** 2 + row["center_y_pc"] ** 2 + row["center_z_pc"] ** 2) ** 0.5
        )

    # Re-running the same shell must skip every already-occupied slot
    # (get_occupied_shell_slots) rather than generating duplicates.
    _run_cli(["--shell", "0", "--num-systems", "1"] + _mysql_argv(mysql_config))
    sectors_after_rerun = _all_sectors(mysql_config)
    assert len(sectors_after_rerun) == 3
    assert {row["id"] for row in sectors_after_rerun} == {row["id"] for row in sectors}


def test_shell_batch_mode_limit_generates_only_first_n_missing_slots(mysql_config):
    _run_cli(["--shell", "0", "--limit", "1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    sectors = _all_sectors(mysql_config)
    assert len(sectors) == 1
    assert sectors[0]["shell_index"] == 0

    # A second --limit 1 run must generate the *next* missing slot, not
    # regenerate the one that already exists.
    _run_cli(["--shell", "0", "--limit", "1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    sectors_after_second_run = _all_sectors(mysql_config)
    assert len(sectors_after_second_run) == 2
    addresses = {(row["shell_index"], row["shell_slot_index"]) for row in sectors_after_second_run}
    assert len(addresses) == 2


def test_local_neighborhood_mode_generates_expected_new_slots_and_skips_occupied(mysql_config):
    # Seed with just one already galaxy-placed sector, shell 0 slot 0 --
    # --center-sector needs an already galaxy-placed sector to search a
    # neighborhood around (a 'generate.py sector'-standalone sector won't do, see
    # test_center_sector_without_galaxy_position_is_rejected below).
    _run_cli(["--shell", "0", "--limit", "1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    seeded = _all_sectors(mysql_config)
    assert len(seeded) == 1
    center_row = seeded[0]
    center_id = center_row["id"]
    center = (center_row["center_x_pc"], center_row["center_y_pc"], center_row["center_z_pc"])

    radius_pc = 8.0
    expected_addresses = {
        (shell_index, slot_index)
        for shell_index, slot_index, _x, _y, _z, _dist in
        enumerate_sectors_within_radius(center, radius_pc, EDGE_PC)
    }
    # Sanity check on the test's own setup: this radius must reach past
    # shell 0 (3 slots total) into at least shell 1, and cover more than
    # just the seed sector itself -- otherwise this wouldn't actually
    # exercise cross-shell neighborhood generation at all.
    assert any(shell_index != 0 for shell_index, _slot in expected_addresses)
    assert len(expected_addresses) > 1
    assert (0, 0) in expected_addresses  # the seed sector itself, at distance 0

    _run_cli([
        "--center-sector", str(center_id), "--radius-pc", str(radius_pc),
        "--num-systems", "1",
    ] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    actual_addresses = {(row["shell_index"], row["shell_slot_index"]) for row in sectors}
    assert actual_addresses == expected_addresses
    # No duplicates: exactly one row per address, including the seed's own.
    assert len(sectors) == len(actual_addresses)

    # Every newly placed sector's stored position matches the geometry it
    # was supposed to be placed at.
    by_address = {(row["shell_index"], row["shell_slot_index"]): row for row in sectors}
    for shell_index, slot_index, x, y, z, _dist in enumerate_sectors_within_radius(center, radius_pc, EDGE_PC):
        row = by_address[(shell_index, slot_index)]
        assert row["center_x_pc"] == pytest.approx(x)
        assert row["center_y_pc"] == pytest.approx(y)
        assert row["center_z_pc"] == pytest.approx(z)

    # Re-running the identical neighborhood must not create duplicates or
    # regenerate anything -- every candidate slot is now already occupied.
    _run_cli([
        "--center-sector", str(center_id), "--radius-pc", str(radius_pc),
        "--num-systems", "1",
    ] + _mysql_argv(mysql_config))
    sectors_after_rerun = _all_sectors(mysql_config)
    assert {row["id"] for row in sectors_after_rerun} == {row["id"] for row in sectors}


def test_center_sector_without_galaxy_position_is_rejected(mysql_config):
    """--center-sector must refuse a sector that was never placed in a
    galaxy (e.g. one generated by 'generate.py sector's own standalone CLI) --
    generate.py's run_local_neighborhood docstring calls this out
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


def test_shell_batch_mode_rejects_large_shell_without_limit_or_yes(mysql_config):
    """Shell 20 holds 3,631 slots (comfortably above
    LARGE_SHELL_WARNING_THRESHOLD's 2,000) -- a bare --shell 20 must be
    refused rather than silently attempting to generate all of them."""
    assert shell_sector_count(20) > galaxyGen.LARGE_SHELL_WARNING_THRESHOLD

    with pytest.raises(SystemExit):
        _run_cli(["--shell", "20", "--num-systems", "1"] + _mysql_argv(mysql_config))

    assert len(_all_sectors(mysql_config)) == 0

    # --limit bypasses the guard even without --yes.
    _run_cli(["--shell", "20", "--limit", "2", "--num-systems", "1"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == 2


# ---------------------------------------------------------------------------
# --shell K --slot N -- single-address mode, the direct path from a
# designation/address copied out of the interactive 3D Galaxy Map
# (html/lib/galaxymap3d.py) into this script. Thin wrapper around
# ensure_sector_generated (already covered in isolation below) plus its
# own argparse validation -- the validation tests need no real database at
# all (parser.error fires during argument parsing, before any connection
# attempt), so they use a dummy --mysql-* argv fragment and run
# unconditionally, unlike every mysql_config-fixture test in this file.
# ---------------------------------------------------------------------------

_DUMMY_MYSQL_ARGV = [
    "--mysql-host", "localhost", "--mysql-port", "3306",
    "--mysql-user", "x", "--mysql-password", "x", "--mysql-database", "x",
]


def test_slot_requires_shell():
    with pytest.raises(SystemExit):
        _run_cli(["--slot", "3"] + _DUMMY_MYSQL_ARGV)


@pytest.mark.parametrize("extra_arg", [
    ["--limit", "1"],
    ["--yes"],
    ["--radius-pc", "5.0"],
    ["--density", "2.0"],
    ["--num-systems", "1"],
])
def test_slot_rejects_flags_that_only_apply_to_other_modes(extra_arg):
    with pytest.raises(SystemExit):
        _run_cli(["--shell", "0", "--slot", "0"] + extra_arg + _DUMMY_MYSQL_ARGV)


def test_slot_rejects_negative_index():
    with pytest.raises(SystemExit):
        _run_cli(["--shell", "0", "--slot", "-1"] + _DUMMY_MYSQL_ARGV)


def test_slot_mode_generates_exactly_the_one_requested_slot(mysql_config):
    shell_index = 0
    n_0 = shell_sector_count(shell_index)
    _seed_skeleton(mysql_config, bands=[(shell_index, 0, 0, n_0 - 1)])

    _run_cli(["--shell", str(shell_index), "--slot", "1"] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    assert len(sectors) == 1
    assert (sectors[0]["shell_index"], sectors[0]["shell_slot_index"]) == (shell_index, 1)
    assert sectors[0]["center_x_pc"] is not None

    # Re-running the identical address must reuse it, not duplicate it.
    _run_cli(["--shell", str(shell_index), "--slot", "1"] + _mysql_argv(mysql_config))
    sectors_after_rerun = _all_sectors(mysql_config)
    assert len(sectors_after_rerun) == 1
    assert sectors_after_rerun[0]["id"] == sectors[0]["id"]


def test_slot_mode_rejects_out_of_range_slot(mysql_config):
    shell_index = 0
    n_0 = shell_sector_count(shell_index)
    _seed_skeleton(mysql_config, bands=[(shell_index, 0, 0, n_0 - 1)])

    with pytest.raises(SystemExit):
        _run_cli(["--shell", str(shell_index), "--slot", str(n_0)] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []


def test_slot_mode_rejects_a_non_qualifying_slot(mysql_config):
    # No stored band at all -- ensure_sector_generated's own "certain no"
    # path (see test_ensure_sector_generated_reports_no_content_outside_every_stored_band).
    _seed_skeleton(mysql_config, bands=[])

    with pytest.raises(SystemExit):
        _run_cli(["--shell", "5000", "--slot", "0"] + _mysql_argv(mysql_config))
    assert _all_sectors(mysql_config) == []


# ---------------------------------------------------------------------------
# ensure_sector_generated -- the visit-triggered lazy-generation entry point
# built on top of the 'generate.py plan' skeleton (galaxy_shape/galaxy_shell_band),
# rather than an explicit --shell/--center-sector batch.
# ---------------------------------------------------------------------------

def test_ensure_sector_generated_raises_without_a_skeleton(mysql_config):
    with pytest.raises(RuntimeError):
        galaxyGen.ensure_sector_generated(0, 0, config=mysql_config)


def test_ensure_sector_generated_creates_then_reuses_the_same_sector(mysql_config):
    n_0 = shell_sector_count(0)
    # Shell 0 sits deep in this toy shape's own bulge -- every slot there
    # clears even a demanding threshold.
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])

    first = galaxyGen.ensure_sector_generated(0, 0, config=mysql_config)
    assert first["created"] is True
    assert first["qualifies"] is True
    assert first["sector_id"] is not None
    assert first["sector_name"]

    second = galaxyGen.ensure_sector_generated(0, 0, config=mysql_config)
    assert second["created"] is False
    assert second["qualifies"] is True
    assert second["sector_id"] == first["sector_id"]

    assert len(_all_sectors(mysql_config)) == 1


def test_ensure_sector_generated_reports_no_content_outside_every_stored_band(mysql_config):
    # A shell with NO stored band at all (e.g. beyond the galaxy's outer
    # edge, or simply never populated) -- find_shell_bands's own bound is
    # exact, so this must be a certain "no", no live density check needed.
    _seed_skeleton(mysql_config, bands=[])

    result = galaxyGen.ensure_sector_generated(5000, 0, config=mysql_config)
    assert result == {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}
    assert _all_sectors(mysql_config) == []


def test_ensure_sector_generated_checks_exact_density_within_a_stored_band(mysql_config):
    """
    A stored candidate band is a safe *superset*, not an exact membership
    list (see `galaxySkeleton`'s own module docstring) -- a slot inside
    the band can still fail the real, exact check. This seeds a band
    deliberately wider than reality (the whole shell) and confirms
    ensure_sector_generated still tells a genuinely dense slot apart from
    a genuinely sparse one within it, rather than trusting the band alone.
    """
    shell_index = 5
    n_k = shell_sector_count(shell_index)

    # Ground truth, computed directly (not via find_shell_bands) --
    # sort every slot in this shell by its own exact relative_density.
    densities = sorted(
        ((relative_density(sector_position_pc(shell_index, i, EDGE_PC), _SKELETON_SHAPE), i)
         for i in range(n_k)),
        reverse=True,
    )
    densest_slot = densities[0][1]
    sparsest_slot = densities[-1][1]
    threshold_rho = (densities[0][0] + densities[-1][0]) / 2.0
    assert densities[0][0] >= threshold_rho > densities[-1][0], (
        "test setup needs a real spread of densities within this shell"
    )

    # A deliberately over-wide band -- the whole shell -- with a
    # threshold picked so only some of it genuinely qualifies.
    _seed_skeleton(mysql_config, e_value=1.0 / threshold_rho, bands=[(shell_index, 0, 0, n_k - 1)])

    dense_result = galaxyGen.ensure_sector_generated(shell_index, densest_slot, config=mysql_config)
    assert dense_result["qualifies"] is True
    assert dense_result["created"] is True

    sparse_result = galaxyGen.ensure_sector_generated(shell_index, sparsest_slot, config=mysql_config)
    assert sparse_result == {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}

    assert len(_all_sectors(mysql_config)) == 1


def test_ensure_sector_generated_passes_relative_density_as_the_density_multiplier(mysql_config, monkeypatch):
    """
    A qualifying sector's actual system count should be driven by its own
    `relative_density` (from the stored skeleton), not a uniform default
    -- confirmed by intercepting `sectorGen.generate_sector` and checking
    what `args.density` it was actually called with, rather than relying
    on the Poisson-sampled system count to differ across a single run.
    """
    n_0 = shell_sector_count(0)
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])

    captured = {}

    def _fake_generate_sector(args, galactic_center_dist_ly=None):
        captured["density"] = args.density
        captured["num_systems"] = args.num_systems
        from stellarObjects.spaceSector import SpaceSector
        return "Fake Sector", SpaceSector(name="Fake Sector")

    monkeypatch.setattr(sectorGen, "generate_sector", _fake_generate_sector)

    position = sector_position_pc(0, 0, EDGE_PC)
    expected_density = relative_density(position, _SKELETON_SHAPE)

    result = galaxyGen.ensure_sector_generated(0, 0, config=mysql_config)
    assert result["created"] is True
    assert captured["density"] == pytest.approx(expected_density)
    assert captured["num_systems"] is None


def test_shell_batch_uses_skeleton_density_when_neither_flag_given(mysql_config, monkeypatch):
    """
    'generate.py galaxy --shell N' with neither --density nor --num-systems
    should drive each sector's own system count from the galaxy skeleton's
    real position-based relative_density (_BatchDensity), the same way
    ensure_sector_generated already does for a single lazily-generated
    sector -- not silently fall back to a flat count shared by the whole
    shell (validate_shared_generation_args's old behavior for every
    subcommand, still correct for the position-less 'sector' subcommand,
    which keeps defaulting to 10 -- galaxy mode alone leaves both None so
    this per-sector path can take over).
    """
    n_0 = shell_sector_count(0)
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])

    captured = []

    def _fake_generate_sector(args, galactic_center_dist_ly=None):
        captured.append((args.density, args.num_systems))
        from stellarObjects.spaceSector import SpaceSector
        return "Fake Sector", SpaceSector(name="Fake Sector")

    monkeypatch.setattr(sectorGen, "generate_sector", _fake_generate_sector)

    _run_cli(["--shell", "0"] + _mysql_argv(mysql_config))

    assert len(captured) == n_0
    for slot_index, (density, num_systems) in enumerate(captured):
        expected_density = relative_density(sector_position_pc(0, slot_index, EDGE_PC), _SKELETON_SHAPE)
        assert density == pytest.approx(expected_density)
        assert num_systems is None


def test_shell_batch_explicit_num_systems_still_overrides_skeleton_density(mysql_config, monkeypatch):
    """An explicit --num-systems is still a uniform, intentional override
    for the whole batch -- _BatchDensity must not second-guess it."""
    n_0 = shell_sector_count(0)
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])

    captured = []

    def _fake_generate_sector(args, galactic_center_dist_ly=None):
        captured.append((args.density, args.num_systems))
        from stellarObjects.spaceSector import SpaceSector
        return "Fake Sector", SpaceSector(name="Fake Sector")

    monkeypatch.setattr(sectorGen, "generate_sector", _fake_generate_sector)

    _run_cli(["--shell", "0", "--num-systems", "4"] + _mysql_argv(mysql_config))

    assert len(captured) == n_0
    assert all(density is None and num_systems == 4 for density, num_systems in captured)


def test_ensure_sector_generated_recovers_from_a_concurrent_insert_race(mysql_config, monkeypatch):
    """
    Simulates the race `sectors`'s `UNIQUE (shell_index, shell_slot_index)`
    constraint (schema.sql's "v8" note) exists to catch: another caller's
    `INSERT` lands between this call's own "not yet generated" check and
    its own `INSERT`. Forces `generate_and_save_sector_at` to raise
    `pymysql.err.IntegrityError` (what a real UNIQUE-constraint violation
    raises) after a sector already exists at the address, and confirms
    `ensure_sector_generated` recovers by returning that sector rather
    than propagating the error.
    """
    n_0 = shell_sector_count(0)
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])

    # A sector generated at a different address, standing in for "the
    # concurrent winner's row" the mocked recovery re-check below returns
    # regardless of which address it's actually asked about.
    winner = galaxyGen.ensure_sector_generated(0, 1, config=mysql_config)
    assert winner["created"] is True

    calls = {"n": 0}

    def _fake_get_sector_id_at(conn, shell_index, shell_slot_index):
        calls["n"] += 1
        if calls["n"] == 1:
            return None  # ensure_sector_generated's own initial "not yet generated" check
        return winner["sector_id"]  # its post-IntegrityError recovery re-check

    monkeypatch.setattr(galaxyGen._db, "get_sector_id_at", _fake_get_sector_id_at)

    def _fake_generate_and_save_sector_at(*_args, **_kwargs):
        raise pymysql.err.IntegrityError(1062, "Duplicate entry for key 'shell_index'")

    monkeypatch.setattr(galaxyGen, "generate_and_save_sector_at", _fake_generate_and_save_sector_at)

    result = galaxyGen.ensure_sector_generated(0, 0, config=mysql_config)
    assert result == {
        "created": False, "qualifies": True,
        "sector_id": winner["sector_id"], "sector_name": None,
    }


def test_sectors_table_rejects_duplicate_shell_slot_address(mysql_config):
    """The schema-level guarantee ensure_sector_generated's own race
    recovery depends on: two sectors can never share a (shell_index,
    shell_slot_index) address."""
    conn = _db.get_connection(mysql_config)
    try:
        from stellarObjects.spaceSector import SpaceSector
        galaxy_position = {
            "center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 3.0,
            "galactic_radius_pc": 3.74, "shell_index": 0, "shell_slot_index": 0,
            "vertices_pc": {"inner": [], "outer": []},
        }
        _db.insert_sector(conn, SpaceSector(name="First"), galaxy_position=galaxy_position)
        with pytest.raises(pymysql.err.IntegrityError):
            _db.insert_sector(conn, SpaceSector(name="Second"), galaxy_position=galaxy_position)
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Random-start mode (no --shell/--center-sector given) -- the galaxy subcommand's
# zero-argument default: pick a random, not-yet-occupied address, generate
# it, then generate every sector within --radius-pc (default 100 ly) of it.
# ---------------------------------------------------------------------------

def test_pick_random_shell_index_stays_within_bounds():
    # No DB needed -- a pure function. Run many draws to catch an
    # off-by-one at either boundary, not just the common case.
    for _ in range(500):
        shell_index = galaxyGen._pick_random_shell_index(0, EDGE_PC)
        assert shell_index == 0

    max_shell_index = 10
    seen = set()
    # Volume weighting means shell 0's own share of the sphere is tiny --
    # measured empirically at ~0.00086 for this max_shell_index -- so a
    # small draw count can genuinely (if rarely) never land there by
    # chance; that's expected behavior, not a bug. 100,000 draws pushes
    # the chance of missing it below 1e-40, so this stays a reliable check
    # of the RNG's own reachability (an off-by-one that made shell 0 or
    # max_shell_index literally unreachable) rather than a coin flip.
    for _ in range(100_000):
        shell_index = galaxyGen._pick_random_shell_index(max_shell_index, EDGE_PC)
        assert 0 <= shell_index <= max_shell_index
        seen.add(shell_index)
    assert seen == set(range(max_shell_index + 1))


def test_pick_random_shell_index_is_volume_weighted_toward_outer_shells():
    # A uniformly-by-volume pick within a sphere of max radius R has a
    # median radius at R * 0.5^(1/3) (~0.7937 R) -- the closed-form median
    # of the r^3-uniform CDF _pick_random_shell_index's own r = R*u^(1/3)
    # sampling implements. With max_shell_index=4253 (this module's own
    # real-galaxy default), the median shell index should land close to
    # that same fraction of the range, and overwhelmingly more draws
    # should fall in the outer half than the inner half -- the opposite of
    # what a naive uniform-by-index pick would give.
    max_shell_index = 4253
    samples = sorted(galaxyGen._pick_random_shell_index(max_shell_index, EDGE_PC) for _ in range(3000))
    median = samples[len(samples) // 2]
    expected_median = max_shell_index * 0.5 ** (1 / 3)
    assert median == pytest.approx(expected_median, rel=0.05)

    outer_half_fraction = sum(1 for s in samples if s > max_shell_index / 2) / len(samples)
    assert outer_half_fraction > 0.8


def test_process_args_with_no_arguments_resolves_to_random_start_mode():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy"]
        args = galaxyGen.process_args()
    finally:
        sys.argv = old_argv

    assert args.shell is None
    assert args.center_sector is None
    assert args.radius_pc is None
    assert args.max_shell is None


def test_process_args_radius_pc_rejected_with_shell():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--shell", "0", "--radius-pc", "10"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()
    finally:
        sys.argv = old_argv


def test_process_args_max_shell_rejected_with_shell_or_center_sector():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--shell", "0", "--max-shell", "10"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()

        sys.argv = ["generate.py", "galaxy", "--center-sector", "1", "--radius-pc", "5", "--max-shell", "10"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()
    finally:
        sys.argv = old_argv


def test_process_args_max_shell_must_be_non_negative():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--max-shell", "-1"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()
    finally:
        sys.argv = old_argv


def test_process_args_accepts_radius_pc_and_max_shell_alone():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--radius-pc", "10", "--max-shell", "5"]
        args = galaxyGen.process_args()
    finally:
        sys.argv = old_argv

    assert args.radius_pc == 10
    assert args.max_shell == 5


def test_process_args_min_start_density_rejected_with_shell_or_center_sector():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--shell", "0", "--min-start-density", "1.5"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()

        sys.argv = [
            "generate.py", "galaxy", "--center-sector", "1", "--radius-pc", "5",
            "--min-start-density", "1.5",
        ]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()
    finally:
        sys.argv = old_argv


def test_process_args_min_start_density_must_be_positive():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--min-start-density", "0"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()

        sys.argv = ["generate.py", "galaxy", "--min-start-density", "-1"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()
    finally:
        sys.argv = old_argv


def test_process_args_min_start_density_rejected_with_density_or_num_systems():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--min-start-density", "1.5", "--num-systems", "5"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()

        sys.argv = ["generate.py", "galaxy", "--min-start-density", "1.5", "--density", "2.0"]
        with pytest.raises(SystemExit):
            galaxyGen.process_args()
    finally:
        sys.argv = old_argv


def test_process_args_accepts_min_start_density_alone():
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "galaxy", "--min-start-density", "1.5"]
        args = galaxyGen.process_args()
    finally:
        sys.argv = old_argv

    assert args.min_start_density == 1.5


def test_random_start_mode_generates_a_seed_sector_and_its_neighborhood(mysql_config):
    # --max-shell 0 pins the randomly chosen seed to one of shell 0's 3
    # slots (deterministic scope, matching this file's other small-radius
    # neighborhood tests) -- a small --radius-pc keeps generation fast.
    radius_pc = 8.0
    _run_cli(
        ["--max-shell", "0", "--radius-pc", str(radius_pc), "--num-systems", "1"] + _mysql_argv(mysql_config)
    )

    sectors = _all_sectors(mysql_config)
    assert len(sectors) >= 1
    # The seed is the first row inserted (generate_and_save_sector_at runs
    # before run_local_neighborhood's own loop) -- every other row's
    # address must fall within radius_pc of its real galaxy-frame center.
    seed = min(sectors, key=lambda row: row["id"])
    assert seed["shell_index"] == 0

    seed_center = (seed["center_x_pc"], seed["center_y_pc"], seed["center_z_pc"])
    expected_addresses = {
        (shell_index, slot_index)
        for shell_index, slot_index, _x, _y, _z, _dist in
        enumerate_sectors_within_radius(seed_center, radius_pc, EDGE_PC)
    }
    actual_addresses = {(row["shell_index"], row["shell_slot_index"]) for row in sectors}
    assert actual_addresses == expected_addresses
    assert len(sectors) == len(actual_addresses)


def test_random_start_mode_retries_until_an_unoccupied_address_is_found(mysql_config):
    # Occupy 2 of shell 0's 3 slots directly -- with --max-shell 0, the
    # random pick must keep retrying (not immediately fail) until it lands
    # on the one remaining unoccupied slot.
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            for slot_index in (0, 1):
                position_pc = sector_position_pc(0, slot_index, EDGE_PC)
                galaxyGen.generate_and_save_sector_at(
                    _fake_args_for_direct_generation(mysql_config), 0, slot_index, position_pc, EDGE_PC,
                )
    finally:
        conn.close()
    assert len(_all_sectors(mysql_config)) == 2

    _run_cli(["--max-shell", "0", "--radius-pc", "0.001", "--num-systems", "1"] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    addresses = {(row["shell_index"], row["shell_slot_index"]) for row in sectors}
    # The only address that could have been newly generated is slot 2 --
    # the sole remaining unoccupied one in shell 0.
    assert (0, 2) in addresses
    assert len(sectors) == 3


def test_random_start_mode_gives_up_after_max_attempts_when_fully_occupied(mysql_config, monkeypatch):
    # Every shell-0 slot already occupied, and --max-shell 0 forces every
    # draw into shell 0 -- random-start mode must eventually give up with
    # a clear error rather than looping forever or crashing obscurely.
    # A tiny attempts cap keeps this test fast.
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            for slot_index in range(shell_sector_count(0)):
                position_pc = sector_position_pc(0, slot_index, EDGE_PC)
                galaxyGen.generate_and_save_sector_at(
                    _fake_args_for_direct_generation(mysql_config), 0, slot_index, position_pc, EDGE_PC,
                )
    finally:
        conn.close()

    monkeypatch.setattr(program_constants, "RANDOM_START_MAX_PLACEMENT_ATTEMPTS", 5)

    with pytest.raises(SystemExit):
        _run_cli(["--max-shell", "0", "--radius-pc", "1.0", "--num-systems", "1"] + _mysql_argv(mysql_config))


def test_random_start_mode_respects_min_start_density(mysql_config):
    """
    --min-start-density must reject an otherwise-qualifying seed address
    whose own real relative_density falls short of it, retrying until one
    that clears it is found -- shell 0's own 3 slots split cleanly under
    _SKELETON_SHAPE (~83.6 for slots 0/1, ~63.1 for slot 2; verified
    directly, not assumed), so a threshold of 75.0 must always land the
    seed on slot 0 or 1, never slot 2.
    """
    n_0 = shell_sector_count(0)
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])

    slot_2_density = relative_density(sector_position_pc(0, 2, EDGE_PC), _SKELETON_SHAPE)
    slot_01_density = relative_density(sector_position_pc(0, 0, EDGE_PC), _SKELETON_SHAPE)
    threshold = 75.0
    assert slot_2_density < threshold < slot_01_density, (
        "test setup needs shell 0's slots to straddle the threshold"
    )

    _run_cli(
        ["--max-shell", "0", "--radius-pc", "0.001", "--min-start-density", str(threshold)]
        + _mysql_argv(mysql_config)
    )

    sectors = _all_sectors(mysql_config)
    seed = min(sectors, key=lambda row: row["id"])
    assert seed["shell_index"] == 0
    assert seed["shell_slot_index"] in (0, 1)
    seed_position = (seed["center_x_pc"], seed["center_y_pc"], seed["center_z_pc"])
    assert relative_density(seed_position, _SKELETON_SHAPE) >= threshold


def test_random_start_mode_gives_up_when_min_start_density_unattainable(mysql_config, monkeypatch, capsys):
    """A --min-start-density no address in range can ever satisfy must
    still give up cleanly (not loop forever), with a message that calls
    out the density constraint specifically. That message goes through
    `log.error` (logged to the console, same as the fully-occupied case
    above) rather than `SystemExit`'s own args, so it's asserted against
    captured stdout instead of `pytest.raises(..., match=...)`."""
    n_0 = shell_sector_count(0)
    _seed_skeleton(mysql_config, bands=[(0, 0, 0, n_0 - 1)])

    monkeypatch.setattr(program_constants, "RANDOM_START_MAX_PLACEMENT_ATTEMPTS", 5)

    with pytest.raises(SystemExit):
        _run_cli(
            ["--max-shell", "0", "--radius-pc", "0.001", "--min-start-density", "1000000"]
            + _mysql_argv(mysql_config)
        )

    assert "min-start-density" in capsys.readouterr().out
    assert _all_sectors(mysql_config) == []


def _fake_args_for_direct_generation(mysql_config):
    """A minimal args namespace shaped like galaxyGen.process_args()'s own
    output, for directly calling generate_and_save_sector_at in a test
    without going through the CLI -- mirrors galaxyGen._default_generation_args
    but pointed at the fixture's own throwaway database, with num_systems
    pinned to 1 (this file's own "keep it fast" convention -- see the
    module docstring)."""
    args = galaxyGen._default_generation_args(config=mysql_config)
    args.density = None
    args.num_systems = 1
    return args


# ---------------------------------------------------------------------------
# End-to-end against a *real* 'generate.py plan' skeleton -- every test
# above either hand-seeds a skeleton via _seed_skeleton (bypassing
# find_shell_bands's own scan entirely) or passes an explicit
# --num-systems/--density that makes _BatchDensity.resolve a no-op (see
# its own docstring), so none of them exercise the actual "plan, then
# galaxy with neither flag given" workflow an operator runs -- exactly the
# combination that shipped a regression where --shell/--center-sector/
# random-start mode saved a real (0-system, 0-phenomenon) sector for every
# not-yet-occupied slot regardless of how far below the qualification
# threshold its own position's density fell, instead of skipping it the
# way ensure_sector_generated always did. These tests close that gap: a
# real skeleton, real CLI generation, no bypass.
# ---------------------------------------------------------------------------

_PLAN_SHAPE_ARGV = [
    "--disk-scale-length-pc", "40", "--disk-scale-height-pc", "12",
    "--bulge-scale-radius-pc", "10", "--bulge-amplitude", "2.0",
    "--arm-count", "2", "--pitch-angle-deg", "15", "--arm-amplitude", "0.4",
]
"""list: The same shape `_SKELETON_SHAPE` above is built from, expressed as
`generate.py plan` CLI args instead -- a small, fast-to-scan toy galaxy
(its real, found outer edge lands a few dozen shells out, not the ~4,100 a
real Milky-Way-scale build reaches) that still has genuine bulge/disk/arm
structure, not a degenerate single-shell case."""


def _build_real_skeleton(mysql_config, extra_argv=()):
    """
    Runs `generate.py plan`'s own `build_skeleton()` -- the real scan
    (`galaxySkeleton.find_shell_bands`), not the hand-seeded
    `_seed_skeleton` shortcut above -- against the fixture's throwaway
    database, and returns its summary dict (`outer_shell_index`, etc.).
    Calls `build_skeleton` directly rather than going through `main()`/
    `sys.argv` purely to get that return value back without parsing
    stdout; it still does the same real persisting
    (`_db.replace_galaxy_shell_bands`/`_db.save_galaxy_shape`) `run_plan`
    itself does.
    """
    parser = argparse.ArgumentParser(prefix_chars='-+')
    galaxyGen.add_plan_arguments(parser)
    args = parser.parse_args(_PLAN_SHAPE_ARGV + list(extra_argv) + _mysql_argv(mysql_config))
    galaxyGen.validate_plan_args(args, parser)
    return galaxyGen.build_skeleton(args)


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


def test_random_start_neighborhood_matches_the_real_skeleton_plan(mysql_config):
    """
    The standard workflow this project's own docs describe -- 'plan' once,
    then 'galaxy' with no flags: pick a random location, generate every
    not-yet-generated sector out to `--radius-pc` (the real default is
    `program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY`, 100 ly;
    trimmed to 25 ly here so this test runs in a reasonable time) -- run
    for real, against a real skeleton, with neither `--density` nor
    `--num-systems` given so every sector's own system count is driven
    end-to-end by its real galaxy-frame position. `-planets` is forced so
    each system skips its own planet/moon tree -- system *count* per
    sector (what this test actually checks) doesn't depend on that, and
    skipping it keeps this test's runtime independent of how dense a
    shell the random draw happens to land in.

    Confirms two things no test above does:
    1. Every candidate slot the run actually saved, and every one it
       didn't, agrees with an independent recomputation -- done here,
       against the real stored skeleton -- of whether that exact position
       qualifies (`predicted_star_count >= 1.0` within a stored candidate
       band). This is the exact regression: a batch/neighborhood run
       saving a sector regardless of qualification.
    2. The aggregate system count actually generated across the
       neighborhood is in the right statistical ballpark of what the
       plan's own density predicted at those same positions -- not a flat,
       position-independent count (`_default_generation_args`'s own
       num_systems=10 default, in particular, would badly fail this).
    """
    radius_pc = ly_to_pc(25.0)
    summary = _build_real_skeleton(mysql_config)
    assert summary["outer_shell_index"] > 0, "the toy shape's own skeleton should find real content"

    _run_cli([
        "--max-shell", "55", "--radius-pc", str(radius_pc), "-planets",
    ] + _mysql_argv(mysql_config))

    sectors = _all_sectors(mysql_config)
    assert sectors, "random-start mode should have generated at least the seed sector"
    seed = min(sectors, key=lambda row: row["id"])
    center = (seed["center_x_pc"], seed["center_y_pc"], seed["center_z_pc"])

    candidates = list(enumerate_sectors_within_radius(center, radius_pc, EDGE_PC))
    assert len(candidates) > 1, "the 25 ly neighborhood should reach beyond just the seed sector itself"

    by_address = {(row["shell_index"], row["shell_slot_index"]): row for row in sectors}

    conn = _db.get_connection(mysql_config)
    try:
        skeleton = _db.get_galaxy_shape(conn)
        bands_cache = {}

        def really_qualifies(shell_index, slot_index, position_pc):
            if shell_index not in bands_cache:
                bands_cache[shell_index] = _db.get_galaxy_shell_bands(conn, shell_index)
            bands = bands_cache[shell_index]
            if not any(lo <= slot_index <= hi for lo, hi in bands):
                return False
            return predicted_star_count(
                position_pc, skeleton.shape, skeleton.expected_system_count_at_density_1,
            ) >= 1.0

        expected_total = 0.0
        generated_sector_ids = []
        for shell_index, slot_index, x, y, z, _dist in candidates:
            address = (shell_index, slot_index)
            saved = address in by_address
            qualifies = really_qualifies(shell_index, slot_index, (x, y, z))
            assert saved == qualifies, (
                f"address {address} was {'saved' if saved else 'skipped'} by 'generate.py galaxy', but "
                f"an independent recomputation against the real stored skeleton says it "
                f"{'qualifies' if qualifies else 'does not qualify'} -- generation has drifted from the "
                f"plan it's supposed to follow."
            )
            if qualifies:
                generated_sector_ids.append(by_address[address]["id"])
                expected_total += predicted_star_count(
                    (x, y, z), skeleton.shape, skeleton.expected_system_count_at_density_1,
                )
    finally:
        conn.close()

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


def test_shell_batch_generates_nothing_beyond_the_real_skeletons_outer_edge(mysql_config):
    """
    Deterministic companion to the neighborhood test above (that one's
    outcome depends on where the random draw lands; this one doesn't): a
    shell chosen well beyond the real skeleton's own discovered outer edge
    must generate exactly zero sectors, every slot skipped -- confirming
    `generate.py galaxy --shell` actually prunes on the real plan rather
    than (as the empty-sectors regression did) saving a sector for every
    slot regardless of position.
    """
    summary = _build_real_skeleton(mysql_config)
    beyond_edge_shell = summary["outer_shell_index"] + 10

    _run_cli([
        "--shell", str(beyond_edge_shell), "--limit", "25", "-planets",
    ] + _mysql_argv(mysql_config))

    assert _all_sectors(mysql_config) == []


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
    from stellarObjects.spaceSector import SpaceSector

    e_value = SpaceSector(name="calibration").expected_system_count()
    args = galaxyGen._default_generation_args()
    args.density = 1.0 / e_value
    args.num_systems = None

    for _ in range(200):
        _sector_name, sector = galaxyGen.generate_sector(args)
        assert sector.entries or sector.phenomena, (
            "generate_sector produced a sector with nothing in it at all"
        )


def test_generate_sector_does_not_force_content_onto_an_explicit_zero():
    """The guaranteed-non-empty fallback only overrides a --density-driven
    zero -- an explicit --num-systems 0 is a deliberate request it must
    never second-guess."""
    args = galaxyGen._default_generation_args()
    args.density = None
    args.num_systems = 0

    _sector_name, sector = galaxyGen.generate_sector(args)
    assert sector.entries == []
    assert sector.phenomena == []


# ---------------------------------------------------------------------------
# Sector-view bounds: every star system AND every placed exotic phenomenon,
# across a real, multi-shell test galaxy, must land within the actual
# coordinate region its own sector's cube occupies in the galaxy frame --
# not just a same-named number compared to a constant.
#
# A sector's cube is NOT axis-aligned with the galaxy's global X/Y/Z: its
# local +Z points radially outward from the galactic center and its local
# +X/+Y follow from that (`sectorGeometry.cube_orientation`), the same
# convention `html/lib/starmap.py`'s `_rotate_to_galaxy_frame` applies to
# every star dot before drawing it against the sector's own wedge outline.
# So "is this point inside the sector" only means something once it's
# re-expressed along *that* sector's own local axes -- `_bounds_violation`
# below does exactly that (re-project the point's offset from the sector's
# center onto its own local axes, then compare each component to
# +/- half the sector's edge), for both:
#   - every star system's `position_x/y/z_mpc` (already stored in that
#     local frame -- `SpaceSector.add_system` samples it directly there),
#     rotated into the galaxy frame the same way the Sector Map does, then
#     re-checked against the real cube; and
#   - every checkable phenomenon's already-absolute `center_x/y/z_pc`
#     (black holes, neutron stars, nebulae, asteroid fields -- the four
#     `sector.phenomena` types `_db.insert_sector` gives a real galaxy-frame
#     center via `_galaxy_placement_from_sector_offset`; a supernova
#     remnant/rogue planet/comet has no coordinate columns of its own at
#     all -- see those tables' own "no placement columns" schema notes --
#     so there is nothing to check for those three).
# ---------------------------------------------------------------------------

_CHECKABLE_PHENOMENON_TABLES = (
    ("black_holes", "black hole"),
    ("neutron_stars", "neutron star"),
    ("nebulae", "nebula"),
    ("asteroid_fields", "asteroid field"),
)
"""tuple: `(table_name, label)` for every phenomenon type that gets its own
real `center_x/y/z_pc` (see this section's own docstring) -- the only ones
a coordinate-bounds check is even possible for."""

_BOUNDS_TOLERANCE_PC = 1e-6
"""float: Floating-point slack for the cube-membership check below -- far
finer than any real generation-scale distance here, purely to absorb
rotation/projection round-off, not to paper over a real violation."""


def _project_onto_axes(vector, axes):
    """`vector`'s own components along each of `axes` (an orthonormal
    basis, e.g. `cube_orientation`'s return value) -- a plain change of
    basis via dot products, the inverse of how a local offset gets rotated
    into the galaxy frame in the first place."""
    return tuple(sum(v * a for v, a in zip(vector, axis)) for axis in axes)


def _sector_cube_contexts(mysql_config):
    """
    For every galaxy-placed sector, the real coordinate region its cube
    occupies: its absolute galaxy-frame center, half-edge (parsecs), and
    its own local cube axes (`cube_orientation`) -- see this section's own
    docstring for why the axes matter and aren't just the global X/Y/Z.

    Returns:
        dict: `{sector_id: {"center_pc": (x, y, z), "half_edge_pc": float,
              "axes": (local_x, local_y, local_z), "shell_index": int,
              "shell_slot_index": int}}`.
    """
    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT id, shell_index, shell_slot_index, edge_mpc, "
            "center_x_pc, center_y_pc, center_z_pc FROM sectors "
            "WHERE center_x_pc IS NOT NULL"
        ).fetchall()
    finally:
        conn.close()

    contexts = {}
    for row in rows:
        center_pc = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        contexts[row["id"]] = {
            "center_pc": center_pc,
            "half_edge_pc": mpc_to_pc(row["edge_mpc"]) / 2.0,
            "axes": cube_orientation(center_pc),
            "shell_index": row["shell_index"],
            "shell_slot_index": row["shell_slot_index"],
        }
    return contexts


def _bounds_violation(sector_ctx, absolute_position_pc, label):
    """`None` if `absolute_position_pc` (a plain galaxy-frame `(x, y, z)`,
    parsecs) falls inside `sector_ctx`'s real cube; otherwise a
    human-readable description of exactly how far outside, and on which of
    the sector's own local axes, it landed."""
    center_pc = sector_ctx["center_pc"]
    delta = tuple(absolute_position_pc[i] - center_pc[i] for i in range(3))
    local_coords = _project_onto_axes(delta, sector_ctx["axes"])
    half_edge = sector_ctx["half_edge_pc"]

    for axis_name, value in zip("xyz", local_coords):
        if not (-half_edge - _BOUNDS_TOLERANCE_PC <= value <= half_edge + _BOUNDS_TOLERANCE_PC):
            return (
                f"{label} (shell {sector_ctx['shell_index']}, slot {sector_ctx['shell_slot_index']}): "
                f"local-{axis_name}={value:.6f} pc, outside +/-{half_edge:.6f} pc"
            )
    return None


def test_stars_and_phenomena_fit_within_their_sectors_real_bounds_across_shells_0_through_4(
    mysql_config, monkeypatch,
):
    """
    Generates a real test galaxy spanning shells ("rings") 0 through 4 --
    518 sectors total (3 + 28 + 79 + 154 + 254, see `shell_sector_count`)
    -- through the real `generate.py galaxy --shell` pipeline, then checks
    every star system's and every checkable phenomenon's actual galaxy-frame
    position against the real coordinate region its own sector's cube
    occupies (see this section's own docstring -- not merely a same-named
    number compared to a constant, but the point re-expressed along that
    sector's own, generally non-axis-aligned, local cube axes).

    Real astrophysical rates make a black hole or nebula genuinely rare
    per sector (see `PHENOMENON_RATE_PER_STAR_SYSTEM`'s own docstring --
    a mean well under 1 per sector for either), so leaving their counts to
    chance would make this test's own phenomenon coverage flaky. `_sample_
    poisson_count` is monkeypatched to always return 1 (for every positive
    mean it's asked for) instead, forcing every one of the seven phenomenon
    types to appear at least once per sector -- deterministic, full
    coverage of every checkable type's own placement path, not a hope that
    a large enough run rolls at least one of each by luck.

    `-planets` keeps each system's own generation cheap (system/phenomenon
    *position*, not planet/moon content, is what this test cares about),
    so the whole 518-sector galaxy (with every sector's phenomena forced
    on top) still generates in well under a minute.
    """
    monkeypatch.setattr(galaxyGen, "_sample_poisson_count", lambda mean, rng=None: 1 if mean > 0 else 0)

    for shell_index in range(5):
        _run_cli(
            ["--shell", str(shell_index), "--num-systems", "3", "-planets"] + _mysql_argv(mysql_config)
        )

    sectors = _all_sectors(mysql_config)
    expected_sector_count = sum(shell_sector_count(k) for k in range(5))
    assert len(sectors) == expected_sector_count

    contexts = _sector_cube_contexts(mysql_config)
    assert len(contexts) == expected_sector_count

    conn = _db.get_connection(mysql_config)
    try:
        star_rows = conn.execute(
            "SELECT sector_id, id AS system_id, position_x_mpc, position_y_mpc, position_z_mpc "
            "FROM star_systems WHERE position_x_mpc IS NOT NULL"
        ).fetchall()

        phenomenon_rows = []
        for table, label in _CHECKABLE_PHENOMENON_TABLES:
            rows = conn.execute(
                f"SELECT sector_id, id, center_x_pc, center_y_pc, center_z_pc FROM {table} "
                f"WHERE center_x_pc IS NOT NULL"
            ).fetchall()
            phenomenon_rows.extend((label, row) for row in rows)
    finally:
        conn.close()

    assert star_rows, "expected at least some star systems across shells 0-4"
    # The forced-count monkeypatch above should make every checkable type
    # show up at least once -- confirms this test's own setup actually
    # exercises all four, not just whichever ones happened to place.
    seen_labels = {label for label, _row in phenomenon_rows}
    assert seen_labels == {label for _table, label in _CHECKABLE_PHENOMENON_TABLES}

    violations = []

    for row in star_rows:
        ctx = contexts[row["sector_id"]]
        local_x, local_y, local_z = ctx["axes"]
        # position_x/y/z_mpc is already this system's position along the
        # sector's own local cube axes (SpaceSector.add_system's own
        # sampling frame -- see spaceSector.py's module docstring), so its
        # absolute galaxy-frame position is the sector's own center plus
        # that local offset rotated into the galaxy frame -- the exact
        # transform html/lib/starmap.py's `_rotate_to_galaxy_frame` applies
        # at render time, reproduced here rather than assumed correct.
        lx = mpc_to_pc(row["position_x_mpc"])
        ly_ = mpc_to_pc(row["position_y_mpc"])
        lz = mpc_to_pc(row["position_z_mpc"])
        rotated = tuple(
            lx * local_x[i] + ly_ * local_y[i] + lz * local_z[i]
            for i in range(3)
        )
        absolute_pc = tuple(ctx["center_pc"][i] + rotated[i] for i in range(3))

        violation = _bounds_violation(ctx, absolute_pc, f"star system {row['system_id']}")
        if violation:
            violations.append(violation)

    for label, row in phenomenon_rows:
        ctx = contexts[row["sector_id"]]
        absolute_pc = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        violation = _bounds_violation(ctx, absolute_pc, f"{label} {row['id']}")
        if violation:
            violations.append(violation)

    total_checked = len(star_rows) + len(phenomenon_rows)
    assert not violations, (
        f"{len(violations)} of {total_checked} generated stars/phenomena fell outside their own "
        f"sector's real (rotated) cube bounds:\n" + "\n".join(violations[:50])
    )
