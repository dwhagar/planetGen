# tests/test_rogue_auto_level.py

"""GEN.1506: the Auto level of the Rogues dial means Few; the four explicit levels are unchanged, and the Generate page and the console help say so."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.galaxy import dials
from planetgen.generation import topup_stage


def test_auto_means_few_and_the_explicit_levels_are_unchanged():
    assert topup_stage.AUTO_LEVEL == "few" and topup_stage.level_word(None) == "Few"
    assert [topup_stage.level_word(word) for word in ("none", "few", "normal", "many")] == ["None", "Few", "Normal", "Many"]
    assert topup_stage.level_word(SimpleNamespace(rogues=None)) == "Few" and topup_stage.level_word(SimpleNamespace(rogues="many")) == "Many"
    assert set(topup_stage.level_word(word) for word in dials.lever("rogues").choices) == set(tuning.ROGUE_LEVEL_CEILINGS)


def test_a_plan_left_on_auto_tops_up_at_few_and_a_pinned_plan_keeps_its_word(monkeypatch):
    monkeypatch.setattr(topup_stage.store, "get_galaxy_shape_pins", lambda conn: [])
    assert topup_stage.galaxy_level(object()) == "few"
    monkeypatch.setattr(topup_stage.dials, "choice_pinned", lambda pins, name: "many")
    assert topup_stage.galaxy_level(object()) == "many"


def test_the_dial_text_and_the_console_help_say_auto_is_few():
    lever = dials.lever("rogues")
    assert lever.auto_text.startswith("few") and lever.choices == ("none", "few", "normal", "many")
    import argparse

    from planetgen.cli import generate

    def help_of(parser):
        for action in parser._actions:
            if "--rogues" in action.option_strings:
                return action.help
            if isinstance(action, argparse._SubParsersAction):
                for sub in action.choices.values():
                    found = help_of(sub)
                    if found:
                        return found
        return None

    assert help_of(generate.build_parser()[0]).rstrip().endswith("else few.")


def test_the_ceilings_behind_the_words_are_unchanged():
    assert tuning.ROGUE_LEVEL_CEILINGS == {"None": 0.0, "Few": 2.5, "Normal": 4.3, "Many": 10.0}
