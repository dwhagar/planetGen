"""
The haze in the air (GEN.851): the fixed point of the air reads the haze of the air it has (`aerosols.haze_of`) and takes the sunlight the haze
absorbs aloft, and the body stores the same haze depth. The bodies are generated, never set: a world whose air is methane-rich enough to haze is
found among the first seeds.
"""

import dataclasses
import sys

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import aerosols, constants, planet_type, surface_air as sa
from planetgen.util import draw


def _generate(seed):
    captured = {}
    real = sa.solve

    def spy(setting, start_temperature_k):
        air = real(setting, start_temperature_k)
        frame = sys._getframe(1)
        while frame is not None:
            if frame.f_code.co_filename.endswith("air_branches.py"):
                return air          # GEN.1049: the second solve of oddity W20 tries the other branch; the body stores the first air unless that branch is taken
            frame = frame.f_back
        captured["setting"], captured["start"], captured["air"] = setting, start_temperature_k, air
        return air

    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound(f"haze-in-the-air-star:{seed}"):
        star = Star(config)
    inner, outer = star.habitable_zone
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(sa, "solve", spy)
        with draw.bound(seed):
            planet = Planet(config, star, star.habitable_zone, (0.2 + (seed % 40) / 10.0) * 0.5 * (inner + outer),
                            mass=constants.EARTH_MASS_TO_KG * (0.5 + seed % 7), moon_count=0)
    return planet, captured


@pytest.fixture(scope="module")
def hazy():
    for seed in range(60):
        planet, captured = _generate(seed)
        if planet.tau_haze and planet.tau_haze > 0.0:
            return planet, captured
    pytest.fail("no body of the first 60 seeds has a haze")


def test_the_haze_the_air_was_solved_with_is_the_haze_the_body_stores(hazy):
    planet, captured = hazy
    assert captured["air"].haze.tau == pytest.approx(planet.tau_haze, rel=0.02)
    assert captured["air"].haze.absorbed_share > 0.0 and captured["air"].haze.pressure_pa < captured["air"].pressure_pa


def test_the_setting_carries_the_draws_of_the_bodys_own_stream(hazy):
    planet, captured = hazy
    assert captured["setting"].haze_draws == aerosols.draw_haze_draws(aerosols.haze_stream_of(planet))
    assert captured["setting"].star_temperature_k == planet.star.temperature


def _replays(captured):
    """The body's own air solved again from the same draws with its haze and without it (TEST.778: both repeat every draw of the body's own solve, `Air.entry_state`, so they are the same world)."""
    with draw.bound("haze-in-the-air-replay"):
        draw.setstate(captured["air"].entry_state)
        with_haze = sa.solve(captured["setting"], captured["start"])
    with draw.bound("haze-in-the-air-replay"):
        draw.setstate(captured["air"].entry_state)
        clear = sa.solve(dataclasses.replace(captured["setting"], haze_draws=None), captured["start"])
    return with_haze, clear


def _same_world(with_haze, clear):
    """Both replays settle and neither sits inside the ramp of the photochemical loss (GEN.1732): there the loss takes a share of the ammonia and methane that depends on the temperature, so an air with its haze (at 229 K, seed 106) and the same air without it (a cold branch at 134 K, no loss) are two airs with two fixed points and the haze does not decide which the solve finds. Outside the ramp the loss is whole or absent and the two replays are the same air."""
    half_ramp = sa.PHOTOCHEMICAL_RAMP_K / 2.0
    return with_haze.converged and clear.converged and all(abs(air.temperature_k - sa.PHOTOCHEMICAL_MIN_K) > half_ramp for air in (with_haze, clear))


@pytest.fixture(scope="module")
def hazy_with_ground():
    for seed in range(400):      # GEN.1646 moved the recipe draws (the first hazy body with ground was seed 98, it had been inside the first 60); GEN.1732 (the photochemical loss of ammonia and methane) moved them again: under 1 percent of the bodies of this sweep are hazy with ground (0 of the first 200 seeds, 5 of seeds 200 to 799). The first one that settles outside the loss's ramp (`_same_world`) is seed 273
        planet, captured = _generate(seed)
        # TEST.778: an air that did not settle (`Air.converged` False, a knee, GEN.1344) has no fixed point to compare: its haze switches on and off between two answers and the
        # stored air is their mean, so the replay with the haze is not "the same air". GEN.1754 moved the water of rocky bodies and with it the first hazy body with ground
        # (seed 98 to seed 15, a 257 K body whose haze switches at a temperature, GEN.1344's knee); the test speaks of airs that settle.
        if planet.tau_haze and planet.tau_haze >= 1.0 and planet_type.has_ground(planet) and captured["air"].converged and not captured["air"].knee and _same_world(*(replays := _replays(captured))):
            return planet, captured, replays
    pytest.fail("no body with ground of the first 400 seeds has a haze")


def test_taking_the_haze_out_of_the_air_warms_the_ground_and_never_cools_it(hazy_with_ground):
    """A haze that absorbs sunlight aloft shields the ground (the anti-greenhouse): the same air, solved on the same stream without the haze,
    is as warm or warmer. A body with ground: a planet under a deep envelope takes one step of the solve from its starting temperature (GEN.1334, GEN.1394) and is
    not a fixed point to compare."""
    _, captured, (with_haze, clear) = hazy_with_ground
    # TEST.778: both solves repeat every draw of the body's own solve (`Air.entry_state`, as `air_branches.other_branch` does), so they are the same world: with a stream of its own the
    # replay drew another inventory and ended with no haze at all (a replayed tau of 0.0 in every body tried), and the test compared two airs that did not differ.
    assert with_haze.converged and clear.converged
    assert with_haze.haze.tau > 0.0, "the replay holds the haze the body stores, or the comparison says nothing"
    assert clear.haze is aerosols.NO_HAZE
    assert clear.temperature_k >= with_haze.temperature_k - 0.5


def test_an_air_of_gases_with_no_rayleigh_entry_has_no_haze():
    """Sodium (class B's second balance gas) has no Rayleigh cross section: an air of it alone is an exosphere, not a column, and does not haze."""
    from planetgen.physics import aerosols

    draws = aerosols.HazeDraws(1.0e4, 0.1, 1.0e-7)
    assert aerosols.haze_of({"Na": 1.0e-3}, 9.8, 300.0, draws, 5772.0, aerosols.FarUv(1.0, 1.0)) is aerosols.NO_HAZE
