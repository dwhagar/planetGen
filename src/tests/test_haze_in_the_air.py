"""
The haze in the air (GEN.851): the fixed point of the air reads the haze of the air it has (`aerosols.haze_of`) and takes the sunlight the haze
absorbs aloft, and the body stores the same haze depth. The bodies are generated, never set: a world whose air is methane-rich enough to haze is
found among the first seeds.
"""

import dataclasses

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import aerosols, constants, surface_air as sa
from planetgen.util import draw


def _generate(seed):
    captured = {}
    real = sa.solve

    def spy(setting, start_temperature_k):
        air = real(setting, start_temperature_k)
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


def test_taking_the_haze_out_of_the_air_warms_the_ground_and_never_cools_it(hazy):
    """A haze that absorbs sunlight aloft shields the ground (the anti-greenhouse): the same air, solved on the same stream without the haze,
    is as warm or warmer."""
    _, captured = hazy
    with draw.bound("haze-in-the-air-replay"):
        with_haze = sa.solve(captured["setting"], captured["start"])
    with draw.bound("haze-in-the-air-replay"):
        clear = sa.solve(dataclasses.replace(captured["setting"], haze_draws=None), captured["start"])
    assert with_haze.converged and clear.converged
    assert clear.haze is aerosols.NO_HAZE
    assert clear.temperature_k >= with_haze.temperature_k - 0.5


def test_an_air_of_gases_with_no_rayleigh_entry_has_no_haze():
    """Sodium (class B's second balance gas) has no Rayleigh cross section: an air of it alone is an exosphere, not a column, and does not haze."""
    from planetgen.physics import aerosols

    draws = aerosols.HazeDraws(1.0e4, 0.1, 1.0e-7)
    assert aerosols.haze_of({"Na": 1.0e-3}, 9.8, 300.0, draws, 5772.0) is aerosols.NO_HAZE
