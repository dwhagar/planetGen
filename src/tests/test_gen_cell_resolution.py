"""
GEN.281: a cell that expects more systems than the cell star ceiling is filled only with the stars at or above the
fill luminosity cap; what is left is unresolved, worked out from the address. `galaxy.resolution` is the one rule the
fill, the stats, the API and the sector page read.
"""

import argparse
import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.api import routes
from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.galaxy import resolution, structure_draw
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.profile import profile_at, sector_profile
from planetgen.galaxy.spherical import cell_volume, sector_position_pc
from planetgen.generation import bright_stars as brightStars, run_galaxy
from planetgen.generation.star_population import placed_star_fraction
from planetgen.physics import constants
from planetgen.physics.units import ly_to_pc
from tests.toy_shape import toy_shape

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SEED = bytes(range(16))
SHAPE = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    structure=structure_draw.draw_structure(SEED),   # a stored shape carries no structure: its seed draws it
)
E_VALUE = 2000.0
EXTENTS = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3)]


def _fake_profile(expected, mix, volume=EDGE_PC ** 3):
    return SimpleNamespace(expected_systems=expected, mix=mix, volume_pc3=volume, edge_pc=EDGE_PC, shape=SHAPE)


def test_a_cell_at_or_under_the_ceiling_is_filled_whole():
    for expected in (10.0, 49_999.0, 50_000.0):
        found = resolution.cell_resolution(_fake_profile(expected, {"old": 1.0}), 50_000.0, 100.0)
        assert not found.capped
        assert found.resolved_systems == found.expected_systems == expected
        assert found.unresolved_systems == 0.0 and found.unresolved_mass_solar == 0.0


def test_a_cell_over_the_ceiling_resolves_the_stars_at_or_above_the_cap_and_the_rest_adds_up():
    mix = {"old": 3.0, "bulge": 1.0}
    found = resolution.cell_resolution(_fake_profile(50_001.0, mix), 50_000.0, 100.0)
    assert found.capped
    laws = brightStars.age_laws_of(SHAPE)
    share = (3.0 * placed_star_fraction(100.0, None, "old", laws) + 1.0 * placed_star_fraction(100.0, None, "bulge")) / 4.0
    assert found.resolved_systems == pytest.approx(50_001.0 * share)
    assert found.resolved_systems + found.unresolved_systems == pytest.approx(50_001.0)
    assert found.unresolved_mass_solar == pytest.approx(found.unresolved_systems * tuning.MEAN_SYSTEM_MASS_SOLAR)
    # A step and not a slope: just over the ceiling keeps about 0.2% of the systems.
    assert found.resolved_systems < 0.01 * found.expected_systems


def test_a_cell_the_cap_would_still_fill_past_the_ceiling_raises_its_cap_until_it_fits():
    """The nuclear cluster's cell expects 15 million systems; the 100 solar luminosity cap would still make about 2,700, more than a 1,000 ceiling (GEN.502: the shorter giant branch leaves fewer bright old stars)."""
    profile = _fake_profile(1.5e7, {"nuclear_cluster": 1.0})
    floor = resolution.cell_resolution(profile, 1e12, 100.0)
    assert not floor.capped
    found = resolution.cell_resolution(profile, 1_000.0, 100.0)
    assert found.capped and found.luminosity_cap_sol > 100.0
    assert found.resolved_systems <= 1_000.0
    # The lowest rung that fits: the one below it would not.
    below = resolution.cell_resolution(profile, 1e12, 100.0).resolved_systems
    one_down = found.luminosity_cap_sol / resolution.CAP_LADDER_STEP
    assert one_down < 100.0 or placed_star_fraction(one_down, None, "nuclear_cluster") * 1.5e7 > 1_000.0
    assert found.resolved_systems + found.unresolved_systems == pytest.approx(1.5e7)
    assert below > 1_000.0
    # A cell the cap already fits keeps the cap it was given.
    assert resolution.cell_resolution(_fake_profile(20_000.0, {"old": 1.0}), 5_000.0, 100.0).luminosity_cap_sol == 100.0


def test_a_higher_cap_resolves_fewer_and_the_scatters_mass_limit_adds_back_the_massive_stars():
    mix = {"young": 1.0, "old": 1.0}
    profile = _fake_profile(1e6, mix)
    caps = [resolution.cell_resolution(profile, 1e5, cap).resolved_systems for cap in (100.0, 1000.0, 9000.0)]
    assert caps[0] > caps[1] > caps[2]
    with_mass = resolution.cell_resolution(profile, 1e5, 9000.0, mass_limit_sol=14.0).resolved_systems
    assert with_mass > caps[2]


def test_a_cells_own_volume_scales_what_it_expects():
    half = _fake_profile(1000.0, {"old": 1.0}, volume=EDGE_PC ** 3 / 2)
    assert resolution.expected_in_cell(half) == pytest.approx(500.0)
    assert resolution.cell_resolution(half, 600.0, 100.0).capped is False
    assert resolution.cell_resolution(half, 400.0, 100.0).capped is True


def test_an_empty_cell_resolves_nothing_and_is_never_capped():
    found = resolution.cell_resolution(_fake_profile(0.0, {"old": 0.0}), 0.5, 100.0)
    assert not found.capped and found.expected_systems == 0.0


def _plan_args(*extra):
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args(list(extra))
    generate_cli.validate_plan_args(args, parser)
    return args


def test_the_plan_options_default_to_none_and_are_checked(capsys):
    args = _plan_args()
    assert args.cell_star_ceiling is None and args.fill_luminosity_cap is None
    args = _plan_args("--cell-star-ceiling", "20000", "--fill-luminosity-cap", "300")
    assert (args.cell_star_ceiling, args.fill_luminosity_cap) == (20000.0, 300.0)
    for argv, message in ((["--cell-star-ceiling", "0"], "--cell-star-ceiling must be at least 1"),
                          (["--fill-luminosity-cap", "49"], "--fill-luminosity-cap"),
                          (["--fill-luminosity-cap", "nan"], "must be a finite number")):
        with pytest.raises(SystemExit):
            _plan_args(*argv)
        assert message in capsys.readouterr().err



@pytest.mark.db
def test_the_settings_are_stored_with_the_galaxy_and_default_to_the_tunings(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=E_VALUE,
                            config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    assert (stored.cell_star_ceiling, stored.fill_luminosity_cap_sol) == (
        tuning.CELL_STAR_CEILING, tuning.FILL_LUMINOSITY_CAP_SOL)
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=E_VALUE,
                            config=mysql_config, cell_star_ceiling=1234.0, fill_luminosity_cap_sol=321.0)
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    assert (stored.cell_star_ceiling, stored.fill_luminosity_cap_sol) == (1234.0, 321.0)


def _seed_galaxy(mysql_config, ceiling, cap):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=E_VALUE,
                            config=mysql_config, cell_star_ceiling=ceiling, fill_luminosity_cap_sol=cap,
                            galaxy_seed=SEED)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)


CROWDED = (7, 0, 0)
"""A cell of the toy galaxy that expects about 29,000 systems (its bulge and disk are tiny, so every cell is crowded)."""


def _crowded_cell():
    return CROWDED, profile_at(sector_position_pc(*CROWDED, EDGE_PC), SHAPE, EDGE_PC, E_VALUE).expected_systems


@pytest.mark.db
def test_a_cell_above_the_ceiling_is_filled_down_to_the_cap_and_no_further(mysql_config):
    address, expected = _crowded_cell()
    cap = 500.0
    _seed_galaxy(mysql_config, ceiling=100.0, cap=cap)
    result = run_galaxy.ensure_sector_generated(*address, config=mysql_config, backfill=False)
    assert result["created"]
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        found = resolution.cell_resolution(sector_profile(address, skeleton), 100.0, cap)
        assert found.capped and found.expected_systems == pytest.approx(expected, rel=0.01)
        stars = conn.execute(
            "SELECT s.luminosity_w AS w FROM stars s JOIN star_systems ss ON ss.id = s.star_system_id"
            " WHERE ss.sector_id = ? AND s.role IN ('primary', 'single')", (result["sector_id"],)).fetchall()
        stats = conn.execute(
            "SELECT bright_level_sol, level_before_fill_sol, expected_systems, actual_systems FROM sector_stats"
            " WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", address).fetchone()
    finally:
        conn.close()
    # Only the stars at or above the cap were made: far fewer systems than the cell expects, none dimmer.
    assert 0 < len(stars) < expected / 4
    assert min(row["w"] for row in stars) >= cap * constants.SOLAR_LUMINOSITY * (1 - 1e-9)
    # The cell is taken to the cap, and what it was aimed at is what it resolves.
    assert stats["level_before_fill_sol"] == cap and stats["bright_level_sol"] == 0
    assert stats["expected_systems"] == pytest.approx(found.resolved_systems)
    assert stats["actual_systems"] == len(stars)


@pytest.mark.db
@pytest.mark.parametrize("cap", [50.0, 100.0])
def test_a_capped_cell_run_without_the_scatter_makes_the_systems_the_resolution_gives(mysql_config, cap):
    """GEN.592: a cell over the ceiling, filled with no star scatter drawn, makes the number of systems
    `cell_resolution` says it resolves (a Poisson draw about it), not none."""
    address, _expected = _crowded_cell()
    _seed_galaxy(mysql_config, ceiling=60.0, cap=cap)
    result = run_galaxy.ensure_sector_generated(*address, config=mysql_config, backfill=False)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        found = resolution.cell_resolution(sector_profile(address, skeleton), 60.0, cap)
        made = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?",
                            (result["sector_id"],)).fetchone()["n"]
    finally:
        conn.close()
    assert found.capped and found.resolved_systems >= 10.0
    assert abs(made - found.resolved_systems) <= 5.0 * math.sqrt(found.resolved_systems) + 2.0


def test_a_capped_cell_draws_no_system_of_its_own_and_an_uncapped_one_is_as_before():
    profile = sector_profile(CROWDED, SimpleNamespace(
        shape=SHAPE, edge_pc=EDGE_PC, expected_system_count_at_density_1=E_VALUE, galaxy_seed=None))
    capped = brightStars.FillContext(profile, min_luminosity_sol=500.0, resolution=resolution.cell_resolution(
        profile, 100.0, 500.0))
    whole = brightStars.FillContext(profile, min_luminosity_sol=500.0, resolution=resolution.cell_resolution(
        profile, 1e9, 500.0))
    unset = brightStars.FillContext(profile, min_luminosity_sol=500.0)
    assert capped.capped and capped.bright_share() == 1.0
    assert not whole.capped and not unset.capped
    assert whole.bright_share() == unset.bright_share() < 0.01


@pytest.mark.db
def test_the_api_reports_what_a_capped_cell_leaves_unresolved_and_nothing_for_a_whole_one(mysql_config):
    address, _expected = _crowded_cell()
    _seed_galaxy(mysql_config, ceiling=100.0, cap=500.0)
    detail = {"placed": True, "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]}
    conn = store.get_connection(mysql_config)
    try:
        body = routes.unresolved_of(conn, detail)
        assert body["ceiling_systems"] == 100.0 and body["luminosity_cap_sol"] == 500.0
        assert body["unresolved_systems"] > 0.99 * body["expected_systems"]
        assert body["unresolved_mass_solar"] == pytest.approx(body["unresolved_systems"] * tuning.MEAN_SYSTEM_MASS_SOLAR,
                                                              rel=1e-3)
        assert routes.unresolved_of(conn, {**detail, "placed": False}) is None
    finally:
        conn.close()
    _seed_galaxy(mysql_config, ceiling=1e9, cap=500.0)
    conn = store.get_connection(mysql_config)
    try:
        assert routes.unresolved_of(conn, detail) is None
    finally:
        conn.close()
