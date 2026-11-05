# tests/test_capture_windows.py

"""
GEN.1082: the capture windows the planet architecture records for a system (`physics/capture_windows.py`): its giant-giant instabilities and their close encounters.
"""

import pytest

from planetgen.physics import capture_windows as cw
from planetgen.physics import constants, irregular_population, planet_architecture as pa
from planetgen.util import draw


def _seed(runaway):
    return pa.Seed(5.0, 1.0e27 if runaway else 1.0e25, 1.0e25, runaway, 1, False, 5.0)


def _architecture(thrown):
    lost = pa.Lost(0.0, 0.0, 0, 0.0, 0, (), (), tuple(thrown), 2.7, constants.SOLAR_MASS_TO_KG)
    return pa.Architecture((), 1.0, 12.0, 0, lost)


def test_a_system_with_no_thrown_pair_has_no_windows():
    assert cw.from_architecture(_architecture([])) == cw.NO_WINDOWS == (0, 0.0)


def test_only_a_pair_with_a_giant_on_both_sides_is_an_instability():
    windows = cw.from_architecture(_architecture([(_seed(False), _seed(True)), (_seed(True), _seed(True)), (_seed(False), _seed(False))]))
    assert windows.instabilities == 1 and windows.encounters == cw.ENCOUNTERS_PER_INSTABILITY == irregular_population.REFERENCE_ENCOUNTERS


def test_each_instability_brings_the_encounters_of_the_reference_one():
    windows = cw.from_architecture(_architecture([(_seed(True), _seed(True))] * 3))
    assert windows.instabilities == 3 and windows.encounters == pytest.approx(3.0 * cw.ENCOUNTERS_PER_INSTABILITY)


def test_the_level_for_the_scatter_step_is_zero_with_none_and_one_from_the_first():
    assert cw.instability_level(cw.NO_WINDOWS) == 0.0
    assert cw.instability_level(cw.Windows(1, 130.0)) == 1.0
    assert cw.instability_level(cw.Windows(4, 520.0)) == 1.0


def test_some_but_not_most_architectures_hold_an_instability():
    held = 0
    for i in range(300):
        architecture = pa.build(constants.SOLAR_MASS_TO_KG, 2.7, 1.0, 0.0, 50.0, stream=draw.Stream(f"windows:{i}"))
        windows = cw.from_architecture(architecture)
        if windows.instabilities:
            held += 1
            assert sum(seed.runaway for seed in architecture.planets) >= 1 and windows.encounters > 0.0
    assert 3 <= held <= 60


def test_a_system_gives_its_star_the_windows_of_its_architecture():
    from planetgen.generation.system import StarSystem
    from tests.test_belt_birth import _config
    config = _config()
    for seed in range(1, 8):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        for star, architecture in system.architectures:
            assert system.capture_windows(star) == cw.from_architecture(architecture)
    assert system.capture_windows(object()) == cw.NO_WINDOWS


def test_an_instability_excites_a_share_of_the_cold_outer_disc():
    from types import SimpleNamespace
    from planetgen.generation import belt_birth
    giant = SimpleNamespace(body_type='g', distance=30.0, mass=300.0 * constants.EARTH_MASS_TO_KG, radius=70000.0, formation_distance_au=30.0)
    origin = SimpleNamespace(origin="outer-disc", inner_au=40.0)
    calm, _, _ = belt_birth._fates(constants.SOLAR_MASS_TO_KG, [giant], origin, 150.0, 0.0)
    stirred, _, _ = belt_birth._fates(constants.SOLAR_MASS_TO_KG, [giant], origin, 150.0, 1.0)
    assert calm["hot_classical"] == 0.0 < stirred["hot_classical"]
    assert stirred["cold_classical"] < calm["cold_classical"]
