# tests/test_population_tech.py

"""POP.8: the draw of a species' six tech indices from its age, and the tech level (`planetgen.population.tech`)."""

import statistics

import pytest

from planetgen import tuning
from planetgen.population import tech
from planetgen.util import draw


def _draws(age, count=600, label="t"):
    return [tech.draw_indices(age, draw.Stream(f"{label}:{age}:{i}")) for i in range(count)]


def test_the_weights_are_the_technological_assessments_and_sum_to_one():
    assert tuning.TECH_WEIGHTS == (0.25, 0.20, 0.20, 0.15, 0.10, 0.10)
    assert sum(tuning.TECH_WEIGHTS) == pytest.approx(1.0)
    assert len(tuning.TECH_DOMAINS) == len(tuning.TECH_WEIGHTS) == len(tuning.TECH_DOMAIN_OFFSETS) == 6


def test_the_assessments_worked_example_gives_4_25():
    indices = dict(zip(tech.DOMAINS, (4, 3, 6, 5, 4, 3)))
    assert tech.tech_level(indices) == 4.25
    assert tech.band_name(4.25) == "Fusion"


def test_a_species_with_the_same_stream_draws_the_same_indices():
    first = tech.draw_indices(5_000.0, draw.Stream("tech:42"))
    assert first == tech.draw_indices(5_000.0, draw.Stream("tech:42"))


def test_every_index_is_an_integer_from_0_to_7_at_any_age():
    for age in (100.0, 266.0, 3_000.0, 1e6, 7.5e9):
        for indices in _draws(age, 200):
            assert set(indices) == set(tech.DOMAINS)
            assert all(isinstance(v, int) and 0 <= v <= 7 for v in indices.values())


def test_a_draw_obeys_the_prerequisite_margins():
    for age in (150.0, 2_000.0, 1e5, 1e9):
        for indices in _draws(age, 300, "clamp"):
            # Rounding can add up to 1 over the unrounded clamp, and a lopsided species ignores it.
            if tech.is_lopsided(indices):
                continue
            assert indices["information"] <= indices["materials"] + 2
            assert indices["propulsion"] <= min(indices["energy"], indices["materials"], indices["information"]) + 2
            assert indices["defense"] <= min(indices["energy"], indices["materials"]) + 2


def test_the_median_species_at_earths_age_has_earths_tech_level():
    levels = [tech.tech_level(indices) for indices in _draws(266.0, 1000)]
    assert 2.9 <= statistics.median(levels) <= 3.3
    assert tech.band_name(statistics.median(levels)) == "Atomic"


def test_tech_grows_with_age_and_old_species_stop_short_of_the_top():
    medians = [statistics.median(tech.tech_level(i) for i in _draws(age, 300)) for age in (100.0, 1_000.0, 50_000.0, 1e8)]
    assert medians == sorted(medians) and medians[0] < 3.0 and medians[-1] > 5.7
    assert medians[-1] < 6.6, "the oldest median must stay below the top band"


def test_collapses_and_stagnation_only_ever_shorten_the_effective_age():
    ages = [tech.effective_age(1e5, draw.Stream(f"eff:{i}")) for i in range(600)]
    assert all(0 < age <= 1e5 for age in ages)
    assert sum(1 for age in ages if age < 0.99e5) / len(ages) > 0.15
    young = [tech.effective_age(100.0, draw.Stream(f"young:{i}")) for i in range(100)]
    assert all(age <= 100.0 for age in young) and sum(1 for age in young if age == 100.0) > 70


def test_bands_name_the_level_and_the_limiting_domain_is_the_lowest_foundation_domain():
    assert [tech.band_name(level) for level in (0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.0)] == [
        "Primal", "Agrarian", "Industrial", "Atomic", "Fusion", "Gravitic", "Antimatter", "Transcendent", "Transcendent"]
    jump_gates = dict(zip(tech.DOMAINS, (3, 0, 3, 3, 5, 3)))
    assert tech.limiting_domain(jump_gates) == "materials"
    assert tech.limiting_domain(dict.fromkeys(tech.DOMAINS, 3)) == "energy"
    assert tech.is_lopsided(jump_gates) and not tech.is_lopsided(dict.fromkeys(tech.DOMAINS, 3))
