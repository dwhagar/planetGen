# tests/test_surface_air_start_independence.py

"""
GEN.1555 (built ahead by Lane 07, rewritten for Lane 14 after TODO's decision): what the surface air solve (`surface_air.solve`) promises about its starting temperature.

GEN.1403 (water-limited weathering) showed the solver choosing the other branch of a bistable air: a G2V body at 2.44 au that sat on `MOIST_LIMIT_K` (375 K, 3.5 bar)
went to 145 K and airless, and a second (seed 3085, 4.74 au) went from 313 K to 120 K, after a small change in weathering. TODO decided that real bistability stays:
a world can truly have a snowball and a warm answer, and the generator starts from one fixed production start (`planets` from the airless temperature, `air_branches` from a cold or a warm
start), so the answer is a function of the setting. The 6 bodies of the G2V census that sit on 375 K are the clamped fixed point under their own drawn design constants; a probe that draws
other constants finds cold answers for them, which is not a second answer of the same setting. So the promise is narrower than "one answer from every start":

* `test_the_same_setting_start_and_draw_state_give_the_same_answer_twice`: the solve is a function of the setting, the start and the draw state, for every captured setting
  (including seed 3085, which is two-valued by start and stays so).
* `test_a_giant_ends_at_one_temperature_whatever_the_start`: a giant has the temperature alone to settle (its pressure is the reference level), and the iteration of
  `MAX_GIANT_STEPS` brings it to one answer; before GEN.1555 a hot giant ended at 1,381 K from a 120 K start and 1,633 K from a 600 K one.
* `test_the_setting_that_is_single_valued_today_stays_so`: seed 3032 keeps its one answer.

The settings are captured from the real generation of G2V seeds (the sweep of reports/gen1403-sweep.py: `SystemConfig` with a G2V star, planets on, comets and binaries off,
`draw.bound(seed)`), then solved from starts below and above the answer with ONE draw stream, so only the start differs. `scripts/gen1555_moist_limit_census.py` scores a change over 100
systems: the bodies sitting exactly on 375.0 K (probed under the production draw state as well as a fresh stream), the bodies near the boiling point of their own pressure, and the
share of air bodies whose answer depends on the start.
"""

import sys

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import surface_air
from planetgen.util import draw

STARTS_K = (100.0, 150.0, 200.0, 250.0, 300.0, 375.0, 450.0, 600.0, 900.0)
"""tuple: Start temperatures of the probe, from well below the cold answer to well above the hot one."""

TARGETS = ((3032, 2.439), (3085, 4.738), (3009, 0.04884), (3010, 0.16039), (3000, 1.01804))
"""tuple: `(seed, orbit in au)` of the two bodies GEN.1403 moved and of three giants (a hot one at 0.049 au, one at 0.16 au and one at 1.0 au); the settings are those of the final solve of the
planet at that distance."""
GIANT_SEEDS = (3009, 3010, 3000)
"""tuple: The seeds of `TARGETS` whose setting is a giant's."""

SAME_TEMPERATURE_K = 5.0
SAME_PRESSURE_RATIO = 1.25
"""float: Two answers are one fixed point when they are within these (a solver tolerance is 0.1 K, so this is wide on purpose)."""
AIR_PA = 1.0e4
"""float: The pressure ratio is compared only for answers with at least this much air; a trace of air under 10 kPa moves by tens of percent with the start and is not a branch."""


def _solve(setting, start):
    """The answer of `solve` from `start` with the probe stream, as `(temperature_k, pressure_pa)`, or None for no air."""
    with draw.bound("gen1555-probe"):
        air = surface_air.solve(setting, start)
    return None if air is None else (air.temperature_k, air.pressure_pa)


def _spread(setting):
    """The answers of `solve` from every start, with the same draw stream, as `(start, temperature_k, pressure_pa)`; a start with no answer is `(start, None, 0.0)`."""
    answers = []
    for start in STARTS_K:
        with draw.bound("gen1555-probe"):
            air = surface_air.solve(setting, start)
        answers.append((start, None if air is None else air.temperature_k, 0.0 if air is None else air.pressure_pa))
    return answers


def _one_fixed_point(answers):
    temperatures = [t for _, t, _ in answers if t is not None]
    pressures = [p for _, t, p in answers if t is not None and p >= AIR_PA]
    if len(temperatures) != len(answers):
        return False
    if max(temperatures) - min(temperatures) > SAME_TEMPERATURE_K:
        return False
    return not pressures or max(pressures) <= SAME_PRESSURE_RATIO * min(pressures)


@pytest.fixture(scope="module")
def settings():
    """The final `Setting` the generation solved for the planet at each target distance, by seed."""
    captured = {seed: [] for seed, _ in TARGETS}
    current = [None]
    real = surface_air.solve

    def recording(setting, start_temperature_k):
        planet = sys._getframe(1).f_locals.get("planet")
        captured[current[0]].append((getattr(planet, "distance", None), setting))
        return real(setting, start_temperature_k)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(surface_air, "solve", recording)
        for seed, _ in TARGETS:
            current[0] = seed
            config = SystemConfig()
            config.STAR_TYPE = "G2V"
            config.PLANETS = True
            config.COMETS = False
            config.BINARY_SYSTEM = False
            with draw.bound(seed):
                StarSystem(system_config=config)
    found = {}
    for seed, au in TARGETS:
        rows = [setting for distance, setting in captured[seed] if distance is not None and abs(distance - au) < max(0.02, 0.01 * au)]
        assert rows, f"seed {seed} has no solve at {au} au: the generation moved, so the seeds of GEN.1555 must be found again (scripts/gen1555_moist_limit_census.py)"
        found[seed] = rows[-1]
    return found


@pytest.mark.slow
def test_the_setting_that_is_single_valued_today_stays_so(settings):
    answers = _spread(settings[3032])
    assert _one_fixed_point(answers), answers


@pytest.mark.slow
def test_the_same_setting_start_and_draw_state_give_the_same_answer_twice(settings):
    for seed, setting in settings.items():
        for start in (150.0, 300.0, 600.0):
            first, second = _solve(setting, start), _solve(setting, start)
            assert first == second, (seed, start, first, second)


@pytest.mark.slow
def test_a_giant_ends_at_one_temperature_whatever_the_start(settings):
    for seed in GIANT_SEEDS:
        assert settings[seed].giant, seed
        answers = _spread(settings[seed])
        assert all(t is not None for _, t, _ in answers), (seed, answers)
        temperatures = [t for _, t, _ in answers]
        assert max(temperatures) - min(temperatures) <= SAME_TEMPERATURE_K, (seed, [(s, round(t, 1)) for s, t, _ in answers])
