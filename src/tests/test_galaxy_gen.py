# tests/test_galaxy_gen.py

"""
End-to-end tests for `galaxyGen.py` -- the galaxy-scale sector generation
CLI script (`docs/design/galaxy-coordinate-system.md` section 8,
`docs/TODO.md` Phase 4). `test_galaxy_geometry.py` already covers
`galaxyGeometry`'s pure functions (shell math, the enumeration primitive)
in isolation; this file is the missing piece -- actually running
`galaxyGen.py`'s two CLI modes (`--shell`, `--center-sector`) against a
real, throwaway MySQL database (see `conftest.py`'s `mysql_config`
fixture) and inspecting the resulting `sectors` rows.

Runs the real CLI entry point (`main()`, via `sys.argv`) rather than
calling internals directly -- this repo's existing convention for testing
a CLI script end-to-end (as opposed to `test_sector_gen.py`, which tests
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

import sys

import pytest

import galaxyGen
import sectorGen
from stellarObjects import _db, program_constants
from stellarObjects.galaxyGeometry import enumerate_sectors_within_radius, shell_sector_count
from stellarObjects.utils import ly_to_pc

EDGE_PC = ly_to_pc(program_constants.DEFAULT_SECTOR_EDGE_LY)


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
    """Runs `galaxyGen.main()` with the given argv (excluding argv[0]),
    restoring `sys.argv` afterward."""
    old_argv = sys.argv
    try:
        sys.argv = ["galaxyGen.py"] + argv
        galaxyGen.main()
    finally:
        sys.argv = old_argv


def _run_sector_gen_cli(argv):
    """Runs `sectorGen.main()` the same way -- used to seed an "unplaced"
    (no galaxy position) sector for `test_center_sector_without_galaxy_position_is_rejected`."""
    old_argv = sys.argv
    try:
        sys.argv = ["sectorGen.py"] + argv
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
    # neighborhood around (a sectorGen.py-standalone sector won't do, see
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
    galaxy (e.g. one generated by sectorGen.py's own standalone CLI) --
    galaxyGen.py's run_local_neighborhood docstring calls this out
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
