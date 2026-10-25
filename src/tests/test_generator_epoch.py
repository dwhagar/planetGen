# tests/test_generator_epoch.py

"""GEN.764: every epoch says which stages of generation it changed (`affects`), and the accessors answer what a result
made under an older epoch has to redo."""

import json

import pytest

from planetgen.util import generator_epoch


def test_every_epoch_after_the_first_names_the_stages_it_changed():
    stages = generator_epoch.stages()
    assert stages == ("plan", "scatter", "fill", "population")
    assert generator_epoch.affects(1) == ()
    for number in range(2, generator_epoch.generator_epoch() + 1):
        listed = generator_epoch.affects(number)
        assert listed and set(listed) <= set(stages), (number, listed)


def test_the_file_is_in_the_form_every_writer_leaves_it():
    text = generator_epoch.EPOCHS_FILE.read_text(encoding="utf-8")
    assert generator_epoch.format_epochs(json.loads(text)) == text


def test_a_result_one_epoch_behind_redoes_what_that_epoch_changed():
    current = generator_epoch.generator_epoch()
    assert generator_epoch.stages_between(current) == ()
    assert generator_epoch.stages_between(current - 1) == tuple(
        s for s in generator_epoch.stages() if s in generator_epoch.affects(current))


def test_the_stages_are_the_union_in_the_order_of_a_galaxy():
    assert generator_epoch.stages_between(2, 4) == ("scatter", "fill", "population")   # epochs 3 and 4
    assert generator_epoch.stages_between(58, 59) == ("population",)


def test_a_plan_epoch_between_means_a_full_regeneration():
    assert generator_epoch.needs_full_regeneration(52, 53)      # epoch 53 changed the nuclear disc
    assert not generator_epoch.needs_full_regeneration(58, 60)


def test_a_result_newer_than_the_target_or_an_unknown_epoch_is_refused():
    with pytest.raises(ValueError, match="newer"):
        generator_epoch.stages_between(9, 5)
    with pytest.raises(ValueError, match="no generator epoch"):
        generator_epoch.affects(0)


def test_short_lists_stay_on_one_line():
    text = generator_epoch.format_epochs({"stages": ["plan", "fill"], "epochs": [{"epoch": 2, "affects": ["fill"]}]})
    assert '"affects": ["fill"]' in text and '"stages": ["plan", "fill"]' in text
