"""
GEN.303: `generation.cell_cap` in config.json is the deployment's cell star ceiling (GEN.281, GEN.364). A plan reads it when the run
starts and stores it with the galaxy; `--cell-star-ceiling` sets one galaxy's own and wins; at the default a galaxy is as before.
"""

import argparse
import json
import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.galaxy import resolution, structure_draw
from planetgen.generation import run_plan
from planetgen.physics.units import ly_to_pc
from planetgen.util import settings
from tests.toy_shape import toy_shape

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0, arm_count=2,
    pitch_angle_rad=math.radians(15), arm_amplitude=0.4, structure=structure_draw.draw_structure(bytes(range(16))),
)


@pytest.fixture
def config(tmp_path, monkeypatch):
    """A config.json of the test's own; writes `{"generation": {...}}` and drops the settings cache."""
    path = tmp_path / "config.json"
    monkeypatch.setattr(settings.logpaths, "CONFIG_PATH", str(path))

    def write(**generation):
        path.write_text(json.dumps({"generation": generation}) if generation else "{}")
        settings.reset_cache()

    write()
    yield write
    settings.reset_cache()


def _plan_args(*extra):
    parser = argparse.ArgumentParser(prefix_chars="-+")
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args(list(extra))
    generate_cli.validate_plan_args(args, parser)
    return args


def _params(*extra):
    return run_plan.skeleton_params(_plan_args(*extra), SimpleNamespace(host="h", port=1, user="u", password="p", database="d"))


def _fake_profile(expected):
    return SimpleNamespace(expected_systems=expected, mix={"old": 1.0}, volume_pc3=EDGE_PC ** 3, edge_pc=EDGE_PC, shape=SHAPE)


def test_the_option_defaults_to_the_tuning_and_is_held_between_500_and_50000(config):
    assert settings.get_settings().generation.cell_cap == tuning.CELL_STAR_CEILING == 5000
    for value in (500, 5000, 50000):
        config(cell_cap=value)
        assert settings.get_settings().generation.cell_cap == value
    for value in (499, 50001, 0, -5, "many"):
        config(cell_cap=value)
        with pytest.raises(settings.SettingsError):
            settings.get_settings()


def test_the_plan_takes_the_setting_at_the_start_of_the_run_and_the_option_on_the_command_line_wins(config):
    assert _params()["cell_star_ceiling"] == tuning.CELL_STAR_CEILING            # as before at the default
    config(cell_cap=1500)
    assert _params()["cell_star_ceiling"] == 1500.0
    assert _params("--cell-star-ceiling", "20000")["cell_star_ceiling"] == 20000.0
    config(cell_cap=40000)
    assert _params()["cell_star_ceiling"] == 40000.0


def test_a_lower_cap_stores_fewer_systems_in_a_crowded_cell(config):
    """The cell is filled down to the fill luminosity cap above the ceiling; a lower ceiling caps more cells and a capped one holds fewer systems."""
    expected = 20_000.0
    fill_cap = tuning.FILL_LUMINOSITY_CAP_SOL
    held = {}
    for cap in (50000, 5000, 500):
        config(cell_cap=cap)
        ceiling = _params()["cell_star_ceiling"]
        held[cap] = resolution.cell_resolution(_fake_profile(expected), ceiling, fill_cap).resolved_systems
    assert held[50000] == expected                              # the cell is under a high ceiling: filled whole
    assert held[500] <= held[5000] < held[50000]
    assert held[500] <= 500.0 and held[5000] <= 5000.0


def test_the_option_is_in_the_example_file_the_docs_and_the_help():
    example = settings.example_dict()["generation"]
    assert example == {"cell_cap": 5000}
    assert "generation.cell_cap" in settings.docs_table()
    parser = argparse.ArgumentParser(prefix_chars="-+")
    generate_cli.add_plan_arguments(parser)
    assert "generation.cell_cap" in parser.format_help().replace("\n", " ").replace("  ", " ")
