# tests/test_life_scores.py

"""
GEN.649: the four life scores (Microbial, Flora, Fauna, Sentient), each 0 (none) to 9 (abundant), from the evolutionary path
a body reached and its surface conditions (`physics/life_scores.py`), stored on every scored planet and moon.
"""

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability_world as hw
from planetgen.physics import life_scores
from planetgen.util import draw

MILESTONES = tuning.LIFE_MILESTONES


def test_no_milestone_is_no_life():
    assert life_scores.scores(1.0, 1.0, None) == (0, 0, 0, 0)


def test_a_kind_of_life_needs_its_stage():
    assert life_scores.scores(0.9, 0.9, "abiogenesis") == (9, 0, 0, 0)
    microbial, flora, fauna, sentient = life_scores.scores(0.9, 0.9, "photosynthesis")
    assert (microbial, 0 < flora <= tuning.LIFE_SCORE_SIMPLE_FLORA_MAX, fauna, sentient) == (9, True, 0, 0)
    assert life_scores.scores(0.9, 0.9, "complex_cells")[2:] == (0, 0)
    assert life_scores.scores(0.9, 0.9, "multicellularity") == (9, 9, 9, 0)
    assert life_scores.scores(0.9, 0.9, "technological_civilization", 4.5, 4.6)[3] == 9


def test_the_surface_sets_how_much():
    assert life_scores.scores(0.2, 0.2, "multicellularity") == (2, 2, 2, 0)
    assert life_scores.scores(0.9, 0.1, "multicellularity")[:3] == (9, 9, 1)  # fauna follows the complex-life likelihood
    assert life_scores.scores(0.04, 0.0, "multicellularity")[:3] == (1, 1, 0)


def test_a_likelihood_under_the_line_is_none():
    assert life_scores.scores(tuning.HABITABILITY_LIFE_SCORE_MIN / 2, 0.9, "multicellularity") == (0, 0, 0, 0)
    assert life_scores.scores(None, None, "multicellularity") == (0, 0, 0, 0)


def test_fauna_never_exceeds_the_flora_it_feeds_on():
    for bio in (0.05, 0.3, 0.6, 1.0):
        for cpx in (0.0, 0.05, 0.3, 0.6, 1.0):
            _, flora, fauna, _ = life_scores.scores(bio, min(bio, cpx), "multicellularity")
            assert fauna <= flora


def test_a_young_civilization_scores_in_proportion_and_never_below_one():
    full = life_scores.scores(0.9, 0.9, "technological_civilization", 1.0, 1.0 + tuning.LIFE_SCORE_SENTIENT_FULL_GYR)[3]
    half = life_scores.scores(0.9, 0.9, "technological_civilization", 1.0, 1.0 + tuning.LIFE_SCORE_SENTIENT_FULL_GYR / 2)[3]
    new = life_scores.scores(0.9, 0.9, "technological_civilization", 1.0, 1.0)[3]
    assert full == 9 and 1 < half < full and new == 1


def test_life_with_no_timeline_is_microbial_only():
    assert life_scores.stage_of("microbial", None, False) == "abiogenesis"
    assert life_scores.stage_of("none", None, False) is None
    assert life_scores.stage_of("microbial", "photosynthesis", True) == "photosynthesis"
    assert life_scores.scores(0.7, 0.7, life_scores.stage_of("microbial", None, False)) == (7, 0, 0, 0)


def test_every_stage_key_is_known():
    assert [life_scores.stage_index(key) for key in MILESTONES] == list(range(len(MILESTONES)))
    assert life_scores.stage_index(None) == -1


def test_generated_bodies_carry_scores_that_agree_with_their_life():
    seen = {"scored": 0, "alive": 0, "giant": 0}
    for seed in range(12):
        with draw.bound(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets + system.secondary_planets:
            for body in (planet, *getattr(planet, "moons", ())):
                if body.body_type == "g":
                    seen["giant"] += 1
                    assert all(getattr(body, name) is None for name in life_scores.FIELDS)
                    continue
                if body.body_type == "a" or body.phi4 is None:
                    continue
                seen["scored"] += 1
                scores = tuple(getattr(body, name) for name in life_scores.FIELDS)
                assert all(0 <= value <= tuning.LIFE_SCORE_MAX for value in scores)
                assert scores[2] <= scores[1]
                if body.life_chemical is None:
                    assert scores == (0, 0, 0, 0)
                else:  # a planet or moon, in any zone, with a chemistry has an evolutionary path (Boss, 2026-10-23)
                    assert body.evolutionary_pace is not None
                if scores[3]:
                    assert body.life_stage == "technological_civilization"
                if scores[2]:
                    assert body.life_stage in ("multicellularity", "technological_civilization")
                if scores[0]:
                    seen["alive"] += 1
                    assert body.phi_bio >= tuning.HABITABILITY_LIFE_SCORE_MIN
    assert seen["scored"] > 10 and seen["giant"] > 0
