# tests/test_surface_air_start_independence.py

"""
GEN.1555 (reproduction, built ahead by Lane 07 for Lane 14): the surface air solve (`surface_air.solve`) reaches one fixed point whatever temperature it starts from.

GEN.1403 (water-limited weathering) showed the solver choosing the other branch of a bistable air: a G2V body at 2.44 au that sat on `MOIST_LIMIT_K` (375 K, 3.5 bar)
went to 145 K and airless, and a second (seed 3085, 4.74 au) went from 313 K to 120 K, after a small change in weathering. Less weathering cannot cool a real world, so this
is the solver landing on a different fixed point, not the feedback. The fixed point has more than one answer for these settings and which one it finds depends on where it
starts (`planets` starts from the airless temperature, `air_branches` from a cold or a warm start).

The two settings are captured from the real generation of G2V seeds 3032 and 3085 (the sweep of reports/gen1403-sweep.py: `SystemConfig` with a G2V star, planets on,
comets and binaries off, `draw.bound(seed)`), then solved from starts below and above the answer with ONE draw stream, so only the start differs.

* `test_a_setting_has_one_fixed_point_whatever_it_starts_from` is the target of the fix. It is marked `xfail(strict=True)` because seed 3085 is start-dependent today
  (starts of 100, 150 and 200 K end at 117 K with 2 bar, 250 to 900 K end at 174 to 223 K with 6 to 24 bar); Lane 14 removes the marker with the fix. A fix that makes it pass
  for 3085 must keep 3032 on its single answer (`test_the_setting_that_is_single_valued_today_stays_so`, which passes today).
* `scripts/gen1555_moist_limit_census.py` scores a fix over 100 systems: the bodies sitting exactly on 375.0 K, those near the boiling point of their own pressure, and the
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

TARGETS = ((3032, 2.439), (3085, 4.738))
"""tuple: `(seed, orbit in au)` of the two bodies GEN.1403 moved; the settings are those of the final solve of the planet at that distance."""

SAME_TEMPERATURE_K = 5.0
SAME_PRESSURE_RATIO = 1.25
"""float: Two answers are one fixed point when they are within these (a solver tolerance is 0.1 K, so this is wide on purpose)."""
AIR_PA = 1.0e4
"""float: The pressure ratio is compared only for answers with at least this much air; a trace of air under 10 kPa moves by tens of percent with the start and is not a branch."""


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
        rows = [setting for distance, setting in captured[seed] if distance is not None and abs(distance - au) < 0.02]
        assert rows, f"seed {seed} has no solve at {au} au: the generation moved, so the seeds of GEN.1555 must be found again (scripts/gen1555_moist_limit_census.py)"
        found[seed] = rows[-1]
    return found


@pytest.mark.slow
def test_the_setting_that_is_single_valued_today_stays_so(settings):
    answers = _spread(settings[3032])
    assert _one_fixed_point(answers), answers


@pytest.mark.slow
@pytest.mark.xfail(strict=True, reason="GEN.1555: the air solve is bistable near the moist limit; seed 3085 ends at 117 K or 174 to 223 K by where it starts. Lane 14 removes this mark with the fix.")
def test_a_setting_has_one_fixed_point_whatever_it_starts_from(settings):
    for seed, _ in TARGETS:
        answers = _spread(settings[seed])
        assert _one_fixed_point(answers), (seed, [(s, None if t is None else round(t, 1), round(p / 1e5, 2)) for s, t, p in answers])
