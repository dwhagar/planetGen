"""
Bright-star pre-placement (`planetgen.generation.bright_stars`, `planetgen
plan`'s scatter step and the fill in `generate_sector`): every scattered
star lands inside a qualifying cell at or above the threshold, the same
seed gives the same stars, and filling a sector builds a system around
each of its stars and links it back.
"""

import argparse
import math
import random
import re

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.galaxy import structure_draw
from planetgen.generation import run_common
from planetgen.generation import run_galaxy
from planetgen.generation import run_plan
from planetgen.util import log
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation.star_population import bright_star_fraction
from planetgen.physics import constants
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.density import predicted_star_count
from planetgen.galaxy.profile import profile_at
from planetgen.galaxy.spherical import neighbor_addresses, sector_address_at, sector_position_pc
from planetgen.physics.units import ly_to_pc
from planetgen import tuning
from tests import first_fill
from tests.toy_shape import toy_shape
from tests import sector_support
from tests.sector_support import backfill_bright_stars, make_sector

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
THRESHOLD = 2500.0
LIGHT_STEP_MASSES = (5.0, 8.0, 12.0, 16.0)
STEP_MASSES = LIGHT_STEP_MASSES
"""The mass steps the tests backfill with (GEN.187). The toy galaxy holds about 2,000 systems per sector at density 1,
so the real steps (1, 2, 5 and 8 solar masses; 9% of all stars at 1) would draw hundreds of thousands of stars in a
backfill; heavier ones keep the same four-step walk quick (`tests/test_mass_backfill.py` checks the real values)."""


@pytest.fixture(autouse=True)
def light_steps(monkeypatch):
    monkeypatch.setattr(tuning, "BRIGHT_STAR_BACKFILL_RING_MASSES_SOL", LIGHT_STEP_MASSES)


# A small toy galaxy (the same one test_galaxy_gen.py seeds), with a high
# systems-per-sector value so a few shells hold enough bright stars to test.
# A seed whose drawn structure has a stellar halo of the Milky Way's order (2.7% of the stars): the halo tier of the
# scatter covers the toy galaxy's clusters and satellites only where a halo is.
SEED = bytes([39]) * 16
SHAPE = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    structure=structure_draw.draw_structure(SEED),   # a stored shape carries no structure: its seed draws it
)
E_VALUE = 2000.0
SMALL_E_VALUE = 20.0
"""A galaxy for a test that makes whole sectors and checks what is done around them, not what is in them (TEST.403): at
`E_VALUE` a sector's fill builds a system around each of the thousand pre-placed stars down to the backfill's floor,
minutes of generation a sector."""
EXTENTS = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3)]
"""`(shell, highest band)` of the toy outline: a small ball, bands -3 to 3 from shell 3 out."""


def _profile(center):
    """The toy galaxy's sector profile at `center` (what a fill reads, GEN.235)."""
    return profile_at(center, SHAPE, EDGE_PC, E_VALUE)


def _scatter(seed=7, **kwargs):
    return list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, THRESHOLD, seed, **kwargs))


def _row_dict(row):
    return dict(zip(store.BRIGHT_STAR_COLUMNS, row))


def test_scattered_stars_sit_in_their_own_cell_inside_the_outline_and_are_bright():
    rows = _scatter()
    assert len(rows) > 20
    for row in map(_row_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
        assert sector_address_at(point, EDGE_PC) == address
        assert address[0] in dict(EXTENTS) and abs(address[1]) <= dict(EXTENTS)[address[0]]
        assert row["luminosity_w"] >= THRESHOLD * constants.SOLAR_LUMINOSITY * 0.99
        assert row["population"] in brightStars.POPULATIONS


def test_the_same_seed_gives_the_same_stars():
    assert _scatter(seed=11) == _scatter(seed=11)
    assert _scatter(seed=11) != _scatter(seed=12)


def test_skipped_addresses_get_no_stars():
    rows = _scatter()
    skip = {(row[0], row[1], row[2]) for row in rows[:5]}
    assert not [row for row in _scatter(skip_addresses=skip) if (row[0], row[1], row[2]) in skip]


def test_the_count_matches_the_expected_bright_share():
    # Expected total = the sum over the outline's cells of E * cell volume * density * fraction (every band and
    # angle bin sampled); the Poisson total lands within 5 sigma of it, plus the binning's own few percent.
    rows = _scatter(seed=3)
    fractions = brightStars.band_fractions(THRESHOLD)
    expected = sum(brightStars.shell_expected_stars(SHAPE, shell, high, EDGE_PC, E_VALUE, fractions,
                                                    band_samples=2 * high + 1, angle_bins=brightStars.ANGLE_BINS)
                   for shell, high in EXTENTS)
    assert expected > 20
    assert abs(len(rows) - expected) < 5 * math.sqrt(expected) + 0.1 * expected + 5


def test_a_band_scatter_draws_only_stars_between_its_limits():
    low = 100.0
    assert _scatter(seed=5, max_luminosity_sol=THRESHOLD) == []
    rows_low = list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, low, 5, max_luminosity_sol=THRESHOLD))
    assert len(rows_low) > len(_scatter(seed=5)) > 0
    for row in map(_row_dict, rows_low):
        assert low * constants.SOLAR_LUMINOSITY * 0.99 <= row["luminosity_w"]
        assert row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY
    # The band's expected share is the difference of the two fractions.
    for population in brightStars.POPULATIONS:
        band = brightStars.bright_band_fraction(low, THRESHOLD, population)
        assert band == pytest.approx(bright_star_fraction(low, population)
                                     - bright_star_fraction(THRESHOLD, population))
        assert band > 0


def test_fill_context_caps_dim_stars_and_sets_their_population():
    fill = brightStars.FillContext(_profile(sector_position_pc(2, 0, 3, EDGE_PC)), min_luminosity_sol=THRESHOLD)
    assert 0.0 < fill.bright_share() < 0.05
    cfg = fill.apply(SystemConfig(), random.Random(1))
    assert cfg.POPULATION in brightStars.POPULATIONS
    assert cfg.MAX_STAR_LUMINOSITY_SOL == THRESHOLD

    no_scatter = brightStars.FillContext(_profile(sector_position_pc(2, 0, 3, EDGE_PC)))
    assert no_scatter.bright_share() == 0.0
    assert no_scatter.apply(SystemConfig()).MAX_STAR_LUMINOSITY_SOL is None


SCATTER_KINDS = ("--mass", "--luminosity", "--phenomena", "--nebula-field", "--band")


def _scatter_args(mysql_config, *extra):
    """The arguments of a `planetgen scatter` run against `mysql_config` at the tests' threshold; the scatter functions
    take their kind from the function called, so `--mass` stands in when `extra` names none."""
    extra = list(extra)
    if not any(kind in extra for kind in SCATTER_KINDS):
        extra.insert(0, "--mass")
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_scatter_arguments(parser)
    args = parser.parse_args([
        "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database, "--bright-star-min-luminosity", str(THRESHOLD), *extra,
    ])
    generate_cli.validate_scatter_args(args, parser)
    return args


def _clear_args(mysql_config, *flags):
    """The arguments of a `planetgen clear-scatter` run against `mysql_config`."""
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_clear_scatter_arguments(parser)
    return parser.parse_args(["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
                              "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
                              "--mysql-database", mysql_config.database, *flags])


def _seed_galaxy(mysql_config, e_value=E_VALUE):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8,
                          expected_system_count_at_density_1=e_value, config=mysql_config, cell_star_ceiling=1e12,
                          galaxy_seed=SEED)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)


def test_plan_scatter_stores_stars_and_fill_builds_their_systems(mysql_config):
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    summary = run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    assert summary["total"] > 20

    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == summary["total"]
        assert store.bright_star_scatter_settings(conn)[0] == THRESHOLD
        target = conn.execute(
            "SELECT shell_index, band_index, wedge_index, COUNT(*) AS n FROM bright_stars"
            " GROUP BY shell_index, band_index, wedge_index ORDER BY n DESC LIMIT 1").fetchone()
        address = (target["shell_index"], target["band_index"], target["wedge_index"])
        stars = store.bright_stars_for_sector(conn, *address)
    finally:
        conn.close()

    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    position = sector_position_pc(*address, EDGE_PC)
    sector_id, _name, sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)
    # The members of an emission region are preplaced too (GEN.321) but are no stored star: the built stars carry the
    # id of their row.
    preplaced = [entry for entry in sector.entries if getattr(entry, "bright_star_id", None) is not None]
    # The backfill (GEN.23) adds this sector's 100-500 L_sun stars first.
    assert len(preplaced) >= len(stars)

    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_stars_for_sector(conn, *address) == []
        rows = conn.execute(
            "SELECT b.id AS id, b.luminosity_w AS stored_w, b.metallicity_fe_h AS stored_fe_h FROM bright_stars b"
            " JOIN star_systems ss ON ss.id = b.star_system_id WHERE ss.sector_id = ?", (sector_id,)).fetchall()
        built = {row["id"]: conn.execute(
            "SELECT s.luminosity_w AS w, s.metallicity_fe_h AS fe_h FROM stars s JOIN bright_stars b"
            " ON b.star_system_id = s.star_system_id WHERE b.id = ?", (row["id"],)).fetchall() for row in rows}
    finally:
        conn.close()
    assert len(rows) == len(preplaced)
    for row in rows:
        # The row's star is a member of the system its tree built; which role it has is the tree's business (the
        # primary and secondary are the top pair of the shape, GEN.631), so it is found by what the row stored.
        (star,) = [star for star in built[row["id"]] if star["w"] == pytest.approx(row["stored_w"])]
        # GEN.158: the system built around a scattered star keeps the [Fe/H] the scatter drew it.
        assert row["stored_fe_h"] is not None and star["fe_h"] == pytest.approx(row["stored_fe_h"])


def test_the_scatter_reports_only_the_shells_that_drew_stars(mysql_config, monkeypatch):
    # GEN.79: "all 5685 layers were generated" hid that only the middle
    # ones held any stars. Shell 8 here draws none.
    _seed_galaxy(mysql_config)
    real = brightStars.scatter_shell
    monkeypatch.setattr(brightStars, "scatter_shell",
                        lambda shape, shell_index, *args, **kwargs:
                        iter(()) if shell_index == 8 else real(shape, shell_index, *args, **kwargs))
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--workers", "1"))
    (line,) = [message for message in messages
               if message.startswith("Bright stars: stars landed in ")]
    landed = int(line.split("stars landed in ")[1].split(" of ")[0])
    assert 0 < landed < len(EXTENTS) and f"of {len(EXTENTS)} shells" in line
    assert not line.endswith("to 8).")


def test_scatter_always_leaves_filled_sectors_out(mysql_config):
    # GEN.30: no refusal and no --force needed; a filled sector never gets stars.
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (0, 0, 0)
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)

    summary = run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    assert summary["total"] > 0
    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_stars_for_sector(conn, *address) == []
    finally:
        conn.close()


def test_going_down_a_band_keeps_the_old_stars_and_adds_only_the_band(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    first = run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "20"))
    # Nor the first fill's own draw down to the backfill's deepest floor (GEN.633).
    first_fill.without_the_deepest_floor(monkeypatch)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (0, 0, 0)
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        before = {row["id"]: row["luminosity_w"] for row in conn.execute("SELECT id, luminosity_w FROM bright_stars").fetchall()}
        seed = store.bright_star_scatter_settings(conn)[1]
    finally:
        conn.close()
    assert len(before) == first["total"]

    band = run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", "100"))
    assert band["total"] > first["total"]
    assert (band["from_luminosity_sol"], band["to_luminosity_sol"]) == (THRESHOLD, 100.0)

    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_star_scatter_settings(conn) == (100.0, seed)
        rows = conn.execute("SELECT id, luminosity_w, shell_index, band_index, wedge_index"
                            " FROM bright_stars").fetchall()
    finally:
        conn.close()
    after = {row["id"]: row["luminosity_w"] for row in rows}
    assert {key: after[key] for key in before} == before
    added = [row for row in rows if row["id"] not in before]
    assert len(added) == band["total"]
    for row in added:
        assert row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY
        assert row["luminosity_w"] >= 100.0 * constants.SOLAR_LUMINOSITY * 0.99
        assert (row["shell_index"], row["band_index"], row["wedge_index"]) != address

    # At or above the stored level there is nothing to add.
    assert run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", "100")) is None
    assert run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", "300")) is None


def test_going_down_a_band_needs_a_scatter_first(mysql_config):
    _seed_galaxy(mysql_config)
    assert run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", "100")) is None
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
    finally:
        conn.close()


def test_scatter_options_conflict():
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_scatter_arguments(parser)
    with pytest.raises(SystemExit):      # exactly one kind
        parser.parse_args([])
    with pytest.raises(SystemExit):
        parser.parse_args(["--mass", "--luminosity"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--band", "--mass", "--down-to", "100"])
    args = parser.parse_args(["--luminosity", "--bright-star-min-luminosity", "49"])
    with pytest.raises(SystemExit):
        generate_cli.validate_scatter_args(args, parser)
    for argv in (["--band"], ["--band", "--down-to", "49"], ["--mass", "--down-to", "100"],
                 ["--luminosity", "--phenomenon-min-mass", "12"], ["--mass", "--compact-min-mass", "2"]):
        with pytest.raises(SystemExit):
            generate_cli.validate_scatter_args(parser.parse_args(argv), parser)
    generate_cli.validate_scatter_args(parser.parse_args(["--band", "--down-to", "100"]), parser)


# --- GEN.23: the backfill around a generated sector --------------------

FLOOR = 100.0


def test_band_stars_stay_inside_their_band():
    from planetgen.generation.star_population import bright_band_fraction, sample_bright_stars
    stars = sample_bright_stars(40, FLOOR, "young", random.Random(5), max_luminosity_sol=THRESHOLD)
    for star in stars:
        luminosity = star["luminosity_w"] / constants.SOLAR_LUMINOSITY
        assert FLOOR * 0.99 <= luminosity < THRESHOLD
    assert bright_band_fraction(FLOOR, THRESHOLD, "young") == pytest.approx(
        bright_star_fraction(FLOOR, "young") - bright_star_fraction(THRESHOLD, "young"))
    assert bright_band_fraction(FLOOR, None, "young") == bright_star_fraction(FLOOR, "young")
    with pytest.raises(ValueError):
        sample_bright_stars(1, THRESHOLD, max_luminosity_sol=FLOOR)


def test_backfilled_cells_get_band_stars_in_their_own_cell():
    addresses = [(2, 0, wedge) for wedge in range(12)] + [(3, 1, wedge) for wedge in range(12)]
    rows = list(brightStars.backfill_cells(SHAPE, addresses, EDGE_PC, E_VALUE, FLOOR, THRESHOLD, 2))
    assert rows
    for row in map(_row_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
        assert address in addresses
        assert sector_address_at(point, EDGE_PC) == address
        assert FLOOR * 0.99 <= row["luminosity_w"] / constants.SOLAR_LUMINOSITY < THRESHOLD

    # The count matches the band's share of each cell's expected stars.
    expected = 0.0
    for address in addresses:
        center = sector_position_pc(*address, EDGE_PC)
        if predicted_star_count(center, SHAPE, E_VALUE) >= 1.0:
            expected += E_VALUE * sum(
                density * (bright_star_fraction(FLOOR, population)
                           - bright_star_fraction(THRESHOLD, population))
                for population, density in brightStars._densities(center, SHAPE).items())
    assert abs(len(rows) - expected) < 5 * math.sqrt(expected) + 5


def _level_rows(conn):
    """Every sector a backfill took to its own level (GEN.44)."""
    return {(row["shell_index"], row["band_index"], row["wedge_index"]): row["bright_level_sol"]
            for row in conn.execute("SELECT * FROM sector_stats WHERE bright_level_sol > 0").fetchall()}


def _mass_rows(conn):
    """Every sector a mass backfill took down (GEN.187), with its mass."""
    return {(row["shell_index"], row["band_index"], row["wedge_index"]): row["bright_mass_sol"]
            for row in conn.execute("SELECT * FROM sector_stats WHERE bright_mass_sol IS NOT NULL").fetchall()}


def _cell_stars(conn, address):
    return conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE shell_index = ? AND band_index = ?"
                        " AND wedge_index = ?", address).fetchone()["n"]


def test_the_backfill_has_its_own_bar_whose_eta_counts_down(mysql_config, monkeypatch):
    """PERF.28: the backfill adds its own bar, ticks it once per sector it
    visits, and the bar's ETA changes as they finish, down to 0."""
    monkeypatch.setattr(tuning, "PROGRESS_BAR_SECONDS", 0.0)   # the bar draws at once, whatever stats are recorded
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    progress = run_common._generation_progress()
    etas = []
    update = progress.update

    def watched_update(task_id, **kwargs):
        before = progress._tasks[task_id].completed
        update(task_id, **kwargs)
        task = progress._tasks[task_id]
        if task.completed > before:   # the closing update to "full" moves nothing
            etas.append(task.fields["rate"].eta(task.total - task.completed))

    progress.update = watched_update
    summary = run_galaxy.backfill_bright_stars_around(
        mysql_config, [sector_position_pc(4, 0, 5, EDGE_PC)], progress=progress)
    assert summary["sectors"] > 0
    (task,) = progress.tasks
    assert task.description.startswith("Bright-star backfill")
    assert task.completed == task.total >= summary["sectors"]
    assert len(etas) == task.total
    assert len({eta for eta in etas if eta is not None}) > 1
    assert etas[-1] == 0


def test_a_single_sector_run_leaves_the_backfill_to_the_end_of_the_run(mysql_config, monkeypatch):
    """PERF.28: `galaxy --wedge` generates its sector without the inline
    backfill (which stalled the sector bar at 0 of 1); the run's `backfill` stage
    does it afterwards, with its own bar."""
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    calls = []
    monkeypatch.setattr(sector_support, "backfill_bright_stars", lambda *args, **kwargs: calls.append(args))
    # shell 30, not the dense shell 4: the same two calls take two minutes there and fifteen near the centre (TEST.417)
    assert make_sector(30, 0, 5, config=mysql_config, backfill=False, num_systems=1)["created"]
    assert calls == []
    assert make_sector(30, 0, 6, config=mysql_config, num_systems=1)["created"]
    assert len(calls) == 1


def test_backfill_fills_sectors_once_and_skips_filled_sectors(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "20"))
    center = sector_position_pc(4, 0, 5, EDGE_PC)

    conn = store.get_connection(mysql_config)
    try:
        before = conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"]
        limit = store.bright_star_mass_limit(conn)
    finally:
        conn.close()
    assert limit is not None and limit > max(STEP_MASSES)

    first = backfill_bright_stars(mysql_config, center)
    assert first["sectors"] > 0 and first["stars"] > 0
    conn = store.get_connection(mysql_config)
    try:
        masses = _mass_rows(conn)
        assert len(masses) == first["sectors"] and set(masses.values()) == set(STEP_MASSES)
        # The nearest step takes the lowest mass; each sector's level is the scatter's.
        assert set(_level_rows(conn).values()) == {THRESHOLD}
        added = conn.execute(
            "SELECT luminosity_w, initial_mass_sol, shell_index, band_index, wedge_index FROM bright_stars"
            " ORDER BY id").fetchall()[before:]
        assert len(added) == first["stars"]
        for row in added:
            address = (row["shell_index"], row["band_index"], row["wedge_index"])
            assert masses[address] <= row["initial_mass_sol"] < limit
            assert row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY
        assert (4, 0, 5) not in masses  # the generated sector itself is not backfilled
        neighbor = min(address for address, mass in masses.items() if mass == STEP_MASSES[0])
        assert store.bright_star_fill_mass_limit(conn, *neighbor) == STEP_MASSES[0]
        assert store.bright_star_fill_level(conn, *neighbor) == THRESHOLD
        assert store.bright_star_fill_mass_limit(conn, 8, 0, 24) == limit  # a sector nobody reached
        # PERF.11: a row a backfill made carries the sector's expected density.
        stats = store.get_sector_stats(conn, *neighbor)
        assert stats["expected_systems"] == pytest.approx(stats["relative_density"] * E_VALUE)
        assert stats["actual_systems"] is None
    finally:
        conn.close()

    # Already that deep: nothing more is drawn.
    assert backfill_bright_stars(mysql_config, center) == {"sectors": 0, "stars": 0}

    # A sector filled before the backfill reaches it gets no new stars.
    far = (8, 0, 0)
    near_far = next(address for address in neighbor_addresses(*far) if address[0] == 7)
    conn = store.get_connection(mysql_config)
    try:
        far_center = sector_position_pc(*far, EDGE_PC)
        conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index,"
                     " center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                     ("Filled", EDGE_PC * 1000, *far, *far_center, math.hypot(far_center[0], far_center[1])))
        conn.commit()
        stars_there = _cell_stars(conn, far)
    finally:
        conn.close()
    backfill_bright_stars(mysql_config, sector_position_pc(*near_far, EDGE_PC))
    conn = store.get_connection(mysql_config)
    try:
        assert _cell_stars(conn, far) == stars_there
        assert far not in _mass_rows(conn)
    finally:
        conn.close()

    # A new plan scatter forgets every sector's level and mass.
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "20"))
    conn = store.get_connection(mysql_config)
    try:
        assert _level_rows(conn) == {} and _mass_rows(conn) == {}
    finally:
        conn.close()


def test_going_down_a_band_gives_backfilled_sectors_only_what_they_lack(mysql_config):
    # A staged scatter (PERF.5) after a mass backfill (GEN.187): a sector the
    # backfill took down to a mass already holds every star born from that
    # mass up, so a band below the scatter's floor adds there only the stars
    # born lighter than that.
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "20"))
    assert backfill_bright_stars(mysql_config, sector_position_pc(4, 0, 5, EDGE_PC))["sectors"]
    conn = store.get_connection(mysql_config)
    try:
        masses = _mass_rows(conn)
    finally:
        conn.close()

    def stars_since(last_id):
        conn = store.get_connection(mysql_config)
        try:
            rows = conn.execute("SELECT id, luminosity_w, initial_mass_sol, shell_index, band_index, wedge_index"
                                " FROM bright_stars WHERE id > ? ORDER BY id", (last_id,)).fetchall()
            top = conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM bright_stars").fetchone()["n"]
            return rows, top, _level_rows(conn), _mass_rows(conn)
        finally:
            conn.close()

    _rows, last, _levels, _masses = stars_since(0)
    run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", "400"))
    added, last, levels, after = stars_since(last)
    assert added and after == masses
    assert {levels[address] for address in masses} == {400.0}
    inside = [row for row in added if (row["shell_index"], row["band_index"], row["wedge_index"]) in masses]
    assert inside
    for row in inside:
        mass = masses[(row["shell_index"], row["band_index"], row["wedge_index"])]
        assert row["initial_mass_sol"] < mass
        assert 400.0 * constants.SOLAR_LUMINOSITY * 0.99 <= row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY

    band = run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", str(FLOOR)))
    added, last, levels, after = stars_since(last)
    assert len(added) == band["total"] and after == masses
    assert {levels[address] for address in masses} == {FLOOR}
    for row in added:
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
        if address in masses:
            assert row["initial_mass_sol"] < masses[address]
            assert row["luminosity_w"] < 400.0 * constants.SOLAR_LUMINOSITY


def _database_now(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.database_now(conn)
    finally:
        conn.close()


def _distance_ly(a, b):
    from planetgen.physics.units import pc_to_ly
    return pc_to_ly(math.dist(sector_position_pc(*a, EDGE_PC), sector_position_pc(*b, EDGE_PC)))


def test_the_backfill_waits_for_the_run_and_fills_down_to_the_step_masses(mysql_config):
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)  # no plan scatter: the backfill draws every luminosity
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    address = (3, 0, 4)
    position = sector_position_pc(*address, EDGE_PC)
    started = _database_now(mysql_config)
    _sector_id, _name, _sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        # GEN.30: nothing around it until the run is done; the sector itself went to the deepest floor at its fill
        # (GEN.633).
        assert _mass_rows(conn) == {address: STEP_MASSES[0]}
    finally:
        conn.close()

    args.shell, args.band, args.wedge = address
    summary = run_galaxy._backfill_on_queue(args, started)
    assert summary["sectors"] > 0 and summary["stars"] > 0
    conn = store.get_connection(mysql_config)
    try:
        masses = _mass_rows(conn)
        # GEN.187: each step a face further out takes the next mass up.
        assert masses == {**run_galaxy.backfill_step_targets({address}, store.get_galaxy_bounds(conn)),
                          address: STEP_MASSES[0]}
        assert set(masses.values()) == set(STEP_MASSES)
        assert _level_rows(conn) == {}  # no scatter, so no luminosity level to record
        assert store.bright_stars_for_sector(conn, *address) == []
        assert store.get_sector_stats(conn, *address)["bright_level_sol"] == 0.0
        for row in conn.execute("SELECT initial_mass_sol, shell_index, band_index, wedge_index"
                                " FROM bright_stars").fetchall():
            assert row["initial_mass_sol"] >= masses[(row["shell_index"], row["band_index"], row["wedge_index"])]
    finally:
        conn.close()

    neighbor = (3, 0, 5)
    fill = run_galaxy._fill_context(args, neighbor)
    assert fill.star_mass_limit_sol == STEP_MASSES[0] and fill.min_luminosity_sol is None


def test_backfill_does_nothing_below_the_scatters_mass_limit(mysql_config):
    # Every star from the scatter's mass limit up is placed already, so a step whose mass is not below it has
    # nothing to add.
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        store.record_bright_star_scatter(conn, THRESHOLD, 1, mass_limit_sol=STEP_MASSES[0])
        conn.commit()
    finally:
        conn.close()
    assert backfill_bright_stars(mysql_config, sector_position_pc(3, 0, 4, EDGE_PC)) == {
        "sectors": 0, "stars": 0}
    conn = store.get_connection(mysql_config)
    try:
        store.record_bright_star_scatter(conn, THRESHOLD, 1, mass_limit_sol=STEP_MASSES[1])
        conn.commit()
    finally:
        conn.close()
    result = backfill_bright_stars(mysql_config, sector_position_pc(3, 0, 4, EDGE_PC))
    assert result["sectors"] > 0
    conn = store.get_connection(mysql_config)
    try:
        assert set(_mass_rows(conn).values()) == {STEP_MASSES[0]}
    finally:
        conn.close()


def test_concurrent_backfills_draw_each_sector_once(mysql_config):
    import threading
    _seed_galaxy(mysql_config)
    center = sector_position_pc(4, 0, 5, EDGE_PC)
    results, errors = [], []

    def run():
        try:
            results.append(backfill_bright_stars(mysql_config, center))
        except Exception as exc:  # pragma: no cover - reported below
            errors.append(exc)

    threads = [threading.Thread(target=run) for _ in range(3)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == sum(
            result["stars"] for result in results)
        assert sum(result["sectors"] for result in results) == len(_mass_rows(conn))
    finally:
        conn.close()


def _stored_stars(conn):
    """Every bright star as its drawn columns, sorted (ids and times left out)."""
    columns = ", ".join(store.BRIGHT_STAR_COLUMNS)
    return sorted(tuple(row.values()) for row in conn.execute(f"SELECT {columns} FROM bright_stars").fetchall())


def test_a_sector_taken_down_in_steps_gets_the_stars_of_one_draw(mysql_config):
    # GEN.44 (Boss 03:25Z), by mass since GEN.187: a sector a run three faces away backfilled to 5 solar masses
    # and a nearer run later took to 1 draws only the lighter bands the second time, and ends up with exactly the
    # stars one backfill down to 1 gives.
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    target = (3, 0, 4)
    conn = store.get_connection(mysql_config)
    try:
        around = run_galaxy.backfill_step_targets({target}, store.get_galaxy_bounds(conn))
    finally:
        conn.close()
    deep = min(address for address, mass in around.items() if mass == STEP_MASSES[2])
    near = min(neighbor_addresses(*target))

    def target_stars():
        conn = store.get_connection(mysql_config)
        try:
            return _mass_rows(conn).get(target), [row for row in _stored_stars(conn) if row[:3] == target]
        finally:
            conn.close()

    backfill_bright_stars(mysql_config, sector_position_pc(*deep, EDGE_PC))
    held, first = target_stars()
    assert held == STEP_MASSES[2]
    second = backfill_bright_stars(mysql_config, sector_position_pc(*near, EDGE_PC))
    held, stepped = target_stars()
    assert held == STEP_MASSES[0] and second["stars"] > 0 and len(stepped) > len(first)
    assert all(row in stepped for row in first)

    conn = store.get_connection(mysql_config)
    try:
        store.clear_bright_stars(conn)
    finally:
        conn.close()
    backfill_bright_stars(mysql_config, sector_position_pc(*near, EDGE_PC))
    assert target_stars() == (STEP_MASSES[0], stepped)


def test_a_sector_holding_stars_but_no_level_is_wiped_and_drawn_again(mysql_config):
    # GEN.44: stars in a sector still at -1, with no galaxy scatter and no mass, are
    # what a failed run left; the backfill wipes them and draws it whole.
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    address = (1, 0, 2)
    center = sector_position_pc(*min(neighbor_addresses(*address)), EDGE_PC)
    backfill_bright_stars(mysql_config, center)
    conn = store.get_connection(mysql_config)
    try:
        clean = _stored_stars(conn)
        store.clear_bright_stars(conn)
        edges = brightStars.mass_band_edges(None)
        stray = list(brightStars.backfill_mass_cells(SHAPE, [address], EDGE_PC, E_VALUE, STEP_MASSES[0], None, None,
                                                     99, edges))[:1] or [
            (*address, 0, 0, 0, "young", "B", "III", 1e31, 1e6, 9000.0, 2000.0 * constants.SOLAR_LUMINOSITY,
             0.1, 0.2, 5.0, 0.15, 7)]
        store.insert_bright_stars(conn, stray)
        conn.commit()
        assert store.sector_bright_levels(conn, [address]).get(address, -1.0) == -1.0
        assert address not in _mass_rows(conn)
    finally:
        conn.close()
    backfill_bright_stars(mysql_config, center)
    conn = store.get_connection(mysql_config)
    try:
        assert _stored_stars(conn) == clean
        assert _mass_rows(conn)[address] == STEP_MASSES[0]
    finally:
        conn.close()


def test_a_filled_sector_records_its_stats_and_a_delete_puts_its_level_back(mysql_config):
    # PERF.11 and GEN.44: a fill writes the sector's expected and actual
    # density, its stars' mean temperature and luminosity, and level 0; the
    # galaxy keeps a decaying average of actual against expected.
    from planetgen.db import edits
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    address = (2, 0, 3)
    position = sector_position_pc(*address, EDGE_PC)
    backfill_bright_stars(mysql_config, sector_position_pc(*min(neighbor_addresses(*address)), EDGE_PC))
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 4
    sector_id, _name, _sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        stats = store.get_sector_stats(conn, *address)
        systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?",
                               (sector_id,)).fetchone()["n"]
        stars = conn.execute("SELECT temperature_k, luminosity_w FROM stars st JOIN star_systems ss"
                             " ON ss.id = st.star_system_id WHERE ss.sector_id = ?", (sector_id,)).fetchall()
        assert stats["bright_level_sol"] == 0.0 and stats["level_before_fill_sol"] == -1.0 and stats["bright_mass_sol"] == STEP_MASSES[0]
        assert stats["actual_systems"] == systems > 0 and stats["actual_stars"] == len(stars)
        assert stats["expected_systems"] == pytest.approx(stats["relative_density"] * SMALL_E_VALUE)
        assert stats["mean_temperature_k"] == pytest.approx(sum(row["temperature_k"] for row in stars) / len(stars))
        assert stats["mean_luminosity_sol"] == pytest.approx(
            sum(row["luminosity_w"] for row in stars) / len(stars) / constants.SOLAR_LUMINOSITY)
        assert stats["filled_at"] is not None
        # DB.14: the raw facts the map colors from, not a baked color.
        ages = conn.execute("SELECT st.age_gy FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id"
                            " WHERE ss.sector_id = ?", (sector_id,)).fetchall()
        assert stats["mean_age_gy"] == pytest.approx(sum(row["age_gy"] for row in ages) / len(ages))
        assert stats["total_luminosity_sol"] == pytest.approx(
            sum(row["luminosity_w"] for row in stars) / constants.SOLAR_LUMINOSITY)
        assert not {"fill_share", "color_r", "color_g", "color_b"} & set(stats)
        average, samples = store.galaxy_density_ratio(conn)
        assert samples == 1 and average == pytest.approx(systems / stats["expected_systems"])

        edits.delete_sector_with_contents(conn, sector_id)
        conn.commit()
        stats = store.get_sector_stats(conn, *address)
        assert stats["bright_level_sol"] == -1.0 and stats["actual_systems"] is None
        assert stats["mean_age_gy"] is None and stats["total_luminosity_sol"] is None
        assert store.galaxy_density_ratio(conn) == (pytest.approx(average), 1)
    finally:
        conn.close()


def test_the_density_ratio_decays_toward_newer_fills(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index)"
                     " VALUES ('Empty', 13046, 1, 0, 0)")
        conn.execute("UPDATE galaxy_shape SET density_ratio_avg = 1.0, density_ratio_samples = ?",
                     (store.DENSITY_RATIO_WINDOW * 5,))
        sector_id = conn.execute("SELECT id FROM sectors WHERE name = 'Empty'").fetchone()["id"]
        store.record_sector_stats(conn, sector_id, (1, 0, 0), sector_position_pc(1, 0, 0, EDGE_PC))
        average, samples = store.galaxy_density_ratio(conn)
        # No systems against some expected: one fill moves a long-running
        # average 1 / DENSITY_RATIO_WINDOW of the way to 0.
        assert samples == store.DENSITY_RATIO_WINDOW * 5 + 1
        assert average == pytest.approx(1.0 - 1.0 / store.DENSITY_RATIO_WINDOW)
    finally:
        conn.close()


class _Everywhere:
    """Bounds that contain every cell."""

    @staticmethod
    def qualifies(_shell, _band):
        return True


class _InsideShell(_Everywhere):
    """Bounds that stop at one shell, and at the equator band."""

    def __init__(self, outer):
        self.outer = outer

    def qualifies(self, shell, band):
        return shell <= self.outer and band == 0


def test_the_default_luminosity_floor_is_9000():
    assert tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL == 9000.0


def test_each_step_a_face_further_out_takes_the_next_mass():
    start = (6, 0, 9)
    targets = run_galaxy.backfill_step_targets({start}, _Everywhere())
    assert start not in targets
    # Step 1 is exactly the face neighbours (no diagonals); the others sit one face further each.
    step_one = {address for address, mass in targets.items() if mass == STEP_MASSES[0]}
    assert step_one == set(neighbor_addresses(*start))
    assert set(targets.values()) == set(STEP_MASSES)
    step_two = {address for address, mass in targets.items() if mass == STEP_MASSES[1]}
    assert step_two and step_two.isdisjoint(step_one)
    for address in step_two:
        assert set(neighbor_addresses(*address)) & step_one
        assert not set(neighbor_addresses(*address)) & {start}


def test_steps_count_from_every_sector_of_the_run_and_stop_at_the_outline():
    run = {(3, 0, 4), (3, 0, 5)}
    targets = run_galaxy.backfill_step_targets(run, _Everywhere())
    assert run.isdisjoint(targets)
    # A sector next to either of the run's sectors is in step 1.
    assert all(targets[a] == STEP_MASSES[0] for a in neighbor_addresses(3, 0, 4) if a not in run)
    assert all(targets[a] == STEP_MASSES[0] for a in neighbor_addresses(3, 0, 5) if a not in run)
    bounded = run_galaxy.backfill_step_targets(run, _InsideShell(5))
    assert bounded and all(shell <= 5 and band == 0 for shell, band, _wedge in bounded)
    assert set(run_galaxy.backfill_step_targets(run, _Everywhere(), step_masses=(3.0,)).values()) == {3.0}


def test_backfill_runs_from_the_edge_of_every_generated_sector(mysql_config):
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    near, far = (2, 0, 0), (7, 0, 20)
    started = _database_now(mysql_config)
    for address in (near, far):
        run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)

    run_galaxy._backfill_on_queue(args, started)
    conn = store.get_connection(mysql_config)
    try:
        every = _mass_rows(conn)
        assert store.bright_stars_for_sector(conn, *far) == []
        assert store.bright_stars_for_sector(conn, *near) == []
        expected = run_galaxy.backfill_step_targets({near, far}, store.get_galaxy_bounds(conn))
    finally:
        conn.close()
    # The union of both runs' steps, the nearer step winning where they overlap; each generated sector went to the
    # deepest floor at its own fill (GEN.633).
    expected = {**expected, near: STEP_MASSES[0], far: STEP_MASSES[0]}
    assert every == expected
    assert every[min(neighbor_addresses(*near))] == STEP_MASSES[0]
    assert every[min(neighbor_addresses(*far))] == STEP_MASSES[0]


@pytest.fixture
def web_progress(tmp_path, monkeypatch):
    """What a run the Generate page started writes to its progress file
    (`planetgen.queue.progress_file`), in order."""
    from planetgen.queue import progress_file

    monkeypatch.setenv(progress_file.ENV_VAR, str(tmp_path / "progress.json"))
    reports = []
    real = progress_file.report

    def record(completed, total=None, description=None, **kwargs):
        reports.append((description, completed, total))
        real(completed, total, description, **kwargs)

    monkeypatch.setattr(progress_file, "report", record)
    return reports


def test_the_backfill_publishes_progress_to_the_web_from_the_start(mysql_config, web_progress):
    """ADM.26: the Generate page's progress file shows the backfill's bar
    (unmeasured while it finds its sectors, then counting them) and ends
    at the full count."""
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    del web_progress[:]
    with run_common._generation_progress() as progress:
        summary = run_galaxy.backfill_bright_stars_around(
            mysql_config, [sector_position_pc(4, 0, 5, EDGE_PC)], progress=progress)
    assert summary["sectors"] > 0
    backfill = [report for report in web_progress if report[0].startswith("Bright-star backfill")]
    assert backfill[0] == ("Bright-star backfill (finding sectors)", 0, None)
    measured = [report for report in backfill if report[2] is not None]
    assert measured[0][1:] == (0, measured[0][2]) and measured[0][2] >= summary["sectors"]
    assert measured[-1][1] == measured[-1][2]


def test_a_backfill_with_nothing_to_draw_leaves_no_unmeasured_bar(mysql_config, web_progress):
    _seed_galaxy(mysql_config)
    with run_common._generation_progress() as progress:
        run_galaxy.backfill_bright_stars_around(
            mysql_config, [sector_position_pc(4, 0, 5, EDGE_PC)], progress=progress)
        run_galaxy.backfill_bright_stars_around(
            mysql_config, [sector_position_pc(4, 0, 5, EDGE_PC)], progress=progress)
        assert not [task for task in progress.tasks if task.total is None]


def test_topping_up_backfilled_sectors_shows_a_bar_on_the_web(mysql_config, web_progress):
    """ADM.26: the band run's second phase (backfilled sectors getting what
    they lack) reports its own progress."""
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    assert backfill_bright_stars(mysql_config, sector_position_pc(4, 0, 5, EDGE_PC))["sectors"]
    del web_progress[:]
    run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", str(FLOOR)))
    topping = [report for report in web_progress if report[0] == "Topping up backfilled sectors"]
    assert topping and topping[0][1:] == (0, topping[0][2]) and topping[0][2] > 0
    assert topping[-1][1] == topping[-1][2]


def test_a_sector_built_around_a_bright_star_gives_it_planets_and_moons(mysql_config):
    # GEN.127 (GitHub #513): generating a sector that holds a backfilled bright star builds that star's whole
    # system (planets, belts, moons), as GEN.72 made it, not just the star.
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT shell_index, band_index, wedge_index, COUNT(*) AS n FROM bright_stars"
            " GROUP BY shell_index, band_index, wedge_index ORDER BY n DESC LIMIT 1").fetchall()
    finally:
        conn.close()
    address = (rows[0]["shell_index"], rows[0]["band_index"], rows[0]["wedge_index"])
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    position = sector_position_pc(*address, EDGE_PC)
    sector_id, _name, sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)
    assert [entry for entry in sector.entries if entry.preplaced]
    conn = store.get_connection(mysql_config)
    try:
        systems = conn.execute(
            "SELECT ss.id AS id, COUNT(p.id) AS planets FROM star_systems ss JOIN bright_stars b"
            " ON b.star_system_id = ss.id LEFT JOIN planets p ON p.star_system_id = ss.id"
            " WHERE ss.sector_id = ? GROUP BY ss.id", (sector_id,)).fetchall()
    finally:
        conn.close()
    assert systems
    # Not every star keeps planets (a close binary, a giant), but a sector's worth of bright systems has some.
    assert sum(row["planets"] for row in systems) > 0, [dict(row) for row in systems]


def test_the_scatter_and_the_backfill_log_the_stars_added_to_each_shell_by_type(mysql_config, monkeypatch):
    # GEN.131: like the sector fill's summary, each shell says how many stars it was given and of what kind.
    _seed_galaxy(mysql_config)
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    summary = run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--workers", "1"))
    lines = [message for message in messages if message.startswith(("Bright stars, shell", "Massive stars, shell"))]
    assert lines and all(": added " in line and "-type" in line or "white dwarf" in line for line in lines)
    # Each line says which unit of the stage it is and how many stars the cut expected there, so the sizes can be read.
    assert all(re.search(r" \(unit [\d,]+ of [\d,]+, about [\d,]+ stars expected\): added ", line) for line in lines)
    added = sum(int(line.split(": added ")[1].split(" stars")[0].replace(",", "")) for line in lines)
    assert added == summary["total"]
    # Each line's kinds add up to its own total.
    for line in lines:
        total, kinds = line.split(": added ")[1].split(" stars: ")
        parts = kinds.rstrip(".").split(", ")
        assert sum(int(part.split(" ", 1)[0]) for part in parts) == int(total.replace(",", ""))

    messages.clear()
    center = sector_position_pc(4, 0, 5, EDGE_PC)
    result = backfill_bright_stars(mysql_config, center)
    lines = [message for message in messages if message.startswith("Bright-star backfill, shell ")]
    assert result["stars"] > 0 and lines
    assert sum(int(line.split(": added ")[1].split(" stars")[0].replace(",", "")) for line in lines) == result["stars"]


def test_a_galaxy_run_scatters_the_stars_then_the_phenomena(mysql_config, monkeypatch):
    """`galaxy` follows the sectors and the link with the star scatter and then the phenomena (GEN.185's last pass), at
    the limit it was given, all as stages of one run (PERF.172); it has no option to leave either out (PERF.123)."""
    _seed_galaxy(mysql_config)
    monkeypatch.setattr(run_common, "_edge_pc", lambda: EDGE_PC)
    monkeypatch.setattr(run_galaxy, "_plan_sectors", lambda args, edge_pc: [
        run_galaxy.SectorsFollower(args, "start", "Starting sector", edge_pc),
        run_galaxy.SectorsFollower(args, "neighborhood", "Neighborhood", edge_pc)])
    runs = []
    monkeypatch.setattr(run_galaxy.followed, "run_followed", lambda args, followers, title, seed=None: runs.append(followers))
    _parser, parsers = generate_cli.build_parser()
    args = parsers["galaxy"].parse_args([
        "--phenomenon-min-mass", "8", "--mysql-host", mysql_config.host,
        "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
        "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database])
    run_galaxy.run_galaxy(args)
    [followers] = runs
    assert [follower.key for follower in followers] == [
        "start", "neighborhood", "link", "mass", "luminosity", "phenomena", "phenomena-index", "nebula-draw",
        "block-lights", "backfill", "population", "territories", "settle", "run-blocks"]
    by_key = {follower.key: follower for follower in followers}
    # The limit is read when the stage opens (PERF.212): the followers carry what to read it from, which gives the limit it was given.
    mass, phenomena = by_key["mass"].params()["settings"], by_key["phenomena"].params()["settings"]
    assert run_plan._bright_star_settings(mass)["mass"]["mass_limit_sol"] == 8.0
    assert run_plan._phenomena_settings(phenomena)["phenomena"]["min_mass_solar"] == 8.0
