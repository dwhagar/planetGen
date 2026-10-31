# tests/test_life_readings.py

"""
GEN.1210: the three words for each kind of life on a world (`planetgen.physics.life_readings`): abundance from the stored 0 to 9 score,
trend from the conditions and the timeline, and security from the time the living class's conditions last as the star brightens.
"""

import math

import pytest

from planetgen.datasets.fetch import DataMissing
from planetgen.physics import atmosphere, life_readings, main_sequence

from tests.test_life_classes import world

SUN = {"mass_sol": 1.0, "age_gyr": 4.6, "fe_h": 0.0}
"""dict: The host star of the Earth-like world, the Sun today."""

LIFE = dict(life_class=1, life_class_variant="a", life_stage="multicellularity", life_microbial=9, life_flora=6, life_fauna=4,
            life_sentient=0, life_system_age_gyr=4.6, life_age_abiogenesis_gyr=0.5, life_age_multicellularity_gyr=3.6)
"""dict: The living world of Class 1a on Earth's conditions: microbes dominant, flora frequent, fauna occasional, no one sentient."""


def living(**over):
    return world(**{**LIFE, **over})


def by_kind(values, star=SUN):
    return {reading.kind: reading for reading in life_readings.readings(values, star)}


@pytest.mark.parametrize("score, word", [(0, "Absent"), (1, "Rare"), (2, "Rare"), (3, "Occasional"), (4, "Occasional"), (5, "Frequent"),
                                         (6, "Frequent"), (7, "Abundant"), (8, "Abundant"), (9, "Dominant")])
def test_each_score_has_its_dafor_word(score, word):
    assert life_readings.abundance_word(score) == word


def test_a_score_that_is_not_stored_has_no_word():
    assert life_readings.abundance_word(None) is None


def test_a_cool_living_world_reads_as_the_example_says():
    cool = living(energy_flux_w_m2=64.0)       # 80 percent of Earth's light: far from the runaway greenhouse
    microbes = by_kind(cool)["microbial"]
    assert (microbes.abundance, microbes.trend, microbes.security) == ("Dominant", "Stable", "Secure")
    assert life_readings.sentence(microbes) == "Microbes: dominant, stable, secure."
    assert life_readings.sentence(by_kind(cool)["sentient"]) == "Sentient life: absent."
    assert microbes.years_left_gyr is None


def test_a_world_near_the_runaway_greenhouse_is_marginal_and_a_closer_one_precarious_then_dying_out():
    # The Sun brightens about 7 percent in a Gyr (`main_sequence`), and the sea boils off when the light passes the runaway limit.
    words = {}
    for share in (1.0, 1.08, 1.10, 1.125):
        read = by_kind(living(energy_flux_w_m2=80.0 * share))["microbial"]
        words[share] = (read.security, read.years_left_gyr)
    assert words[1.0][0] == "Marginal" and 1.0 < words[1.0][1] <= 2.0
    assert words[1.08][0] == "Precarious" and words[1.08][1] <= 0.5
    assert words[1.125][0] == "Dying out" and words[1.125][1] <= life_readings.DYING_GYR
    assert words[1.0][1] > words[1.08][1] >= words[1.10][1] >= words[1.125][1]


def test_a_star_that_brightens_makes_a_declining_trend_and_a_steady_one_does_not():
    hot = by_kind(living(energy_flux_w_m2=80.0 * 1.08))
    assert {read.trend for read in hot.values() if read.score} == {"Declining"}
    cool = by_kind(living(energy_flux_w_m2=64.0))
    assert {read.trend for read in cool.values() if read.score} == {"Stable"}


def test_a_kind_that_reached_its_milestone_recently_is_still_growing():
    values = living(energy_flux_w_m2=64.0, life_age_multicellularity_gyr=4.5)
    reads = by_kind(values)
    assert reads["flora"].trend == "Growing" and reads["fauna"].trend == "Growing"
    assert reads["microbial"].trend == "Stable"


def test_conditions_that_improve_make_a_growing_trend(monkeypatch):
    def likelihood(values, row, ratio):
        return 0.4 if ratio == 1.0 else 0.6

    monkeypatch.setattr(life_readings, "_likelihood", likelihood)
    assert by_kind(living())["microbial"].trend == "Growing"


def test_a_world_with_no_life_reads_absent_for_every_kind_with_no_trend_and_no_security():
    barren = world(life_class=None, life_class_variant=None, life_stage=None, life_microbial=0, life_flora=0, life_fauna=0, life_sentient=0)
    for read in life_readings.readings(barren, SUN):
        assert (read.abundance, read.trend, read.security) == ("Absent", None, None)
        assert life_readings.sentence(read).endswith(": absent.")


def test_scores_with_no_living_class_still_read_absent():
    ageless = living(life_stage=None)       # a class the conditions suit, with no life begun (GEN.1001)
    assert {read.abundance for read in life_readings.readings(ageless, SUN)} == {"Absent"}


def test_a_gas_giant_or_a_body_with_no_kind_scores_has_no_readings():
    assert life_readings.readings(world(body_type="g", life_microbial=None), SUN) == []
    assert life_readings.readings(world(life_microbial=None, life_flora=None, life_fauna=None, life_sentient=None), SUN) == []


def test_a_living_class_that_scores_nothing_now_is_dying_out():
    # An air with no carbon dioxide fails the class's gate: the life is there from before, the conditions are gone.
    values = living(**{atmosphere.GAS_FIELD: {"O2": 21278.0, "N2": 78084.0, "Ar": 934.0, "H2O": 1000.0}})
    assert by_kind(values)["microbial"].security == "Dying out"


def test_without_a_star_the_trend_is_the_timeline_and_the_security_is_todays_score(monkeypatch):
    reads = by_kind(living(life_age_multicellularity_gyr=4.5), star=None)
    assert reads["flora"].trend == "Growing" and reads["microbial"].trend == "Stable"
    assert reads["microbial"].security == "Secure" and reads["microbial"].years_left_gyr is None
    monkeypatch.setattr(life_readings, "_likelihood", lambda values, row, ratio: 0.2)    # a score of 2 of 9
    assert by_kind(living(), star=None)["microbial"].security == "Marginal"


def test_without_the_isochrones_the_star_is_not_read_and_nothing_breaks(monkeypatch):
    def missing(*args, **kwargs):
        raise DataMissing("isochrones")

    monkeypatch.setattr(main_sequence, "state", missing)
    read = by_kind(living())["microbial"]
    assert read.trend == "Stable" and read.security == "Secure" and read.years_left_gyr is None


def test_a_star_that_leaves_the_main_sequence_ends_the_time_left():
    old = {"mass_sol": 1.0, "age_gyr": 9.9, "fe_h": 0.0}       # a Sun-like star's main sequence is about 10 Gyr
    read = by_kind(living(energy_flux_w_m2=64.0), star=old)["microbial"]
    assert read.years_left_gyr is not None and read.years_left_gyr <= 0.5 and read.security in ("Precarious", "Dying out")


def test_the_readings_change_no_value_and_draw_nothing():
    values = living()
    before = {key: value for key, value in values.items()}
    first = life_readings.readings(values, SUN)
    assert values == before and life_readings.readings(values, SUN) == first


def test_star_of_takes_an_object_a_mapping_or_nothing():
    class Star:
        mass = 1.988409870698051e30
        age = 4.6
        metallicity_fe_h = 0.1

    assert life_readings.star_of(Star()) == {"mass_sol": pytest.approx(1.0), "age_gyr": 4.6, "fe_h": 0.1}
    assert life_readings.star_of(SUN) is SUN and life_readings.star_of(None) is None and life_readings.star_of({"mass_sol": 1.0}) is None


def test_the_search_steps_end_at_the_marginal_horizon():
    assert math.isclose(life_readings.STEP_GYR * round(life_readings.MARGINAL_GYR / life_readings.STEP_GYR), life_readings.MARGINAL_GYR)
