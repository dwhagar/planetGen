# tests/test_galaxy_dials.py

"""
The lever table (GEN.601; planetgen/galaxy/dials.py): the one definition of the galaxy shape inputs. The command-line
flags, the Generate form's fields, the help lines and the table in docs/cli.md are all built from it, and these tests
fail when one is missing for a row, or exists without a row.
"""

import argparse
import pathlib

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.galaxy import dials
from planetgen.web import generate_page

ROOT = pathlib.Path(__file__).resolve().parents[2]


def _plan_parser():
    _parser, parsers = generate_cli.build_parser()
    return parsers["plan"]


def _shape_group_flags(parser):
    for group in parser._action_groups:
        if group.title.startswith("galaxy shape"):
            return {option for action in group._group_actions for option in action.option_strings}
    raise AssertionError("no galaxy shape group")


# --- the table itself -------------------------------------------------------------------------------------------

def test_every_row_is_well_formed():
    names = [row.name for row in dials.LEVERS]
    flags = [row.flag for row in dials.LEVERS]
    assert len(set(names)) == len(names) and len(set(flags)) == len(flags)
    dial_keys = {key for key, _label in dials.DIALS}
    for row in dials.LEVERS:
        assert row.kind in dials.KINDS, row.name
        assert row.dial in dial_keys, row.name
        assert row.flag == "--" + row.name.replace("_", "-"), row.name
        assert row.law and row.auto_text and row.help and row.label, row.name
        assert row.source.strip() and row.status in dials.STATUSES, row.name
        if row.kind == "enum":
            assert row.choices and row.minimum is None and row.maximum is None, row.name
        else:
            assert row.minimum is not None and row.maximum is not None and row.minimum < row.maximum, row.name
            assert row.step and row.step > 0, row.name
        if row.kind == "length":
            assert row.units and any(unit.factor == 1.0 for unit in row.units), row.name


def test_an_auto_value_stands_only_for_a_row_with_a_fixed_auto():
    for row in dials.LEVERS:
        if row.auto_value is not None:
            assert row.problem(row.auto_value) is None, row.name
            assert not row.auto_text.startswith("drawn"), row.name


# --- the command line -------------------------------------------------------------------------------------------

def test_the_plan_flags_are_the_table():
    parser = _plan_parser()
    assert _shape_group_flags(parser) == {row.flag for row in dials.LEVERS}


def test_no_flag_has_a_default_other_than_auto():
    parsed = vars(_plan_parser().parse_args([]))
    for row in dials.LEVERS:
        assert parsed[row.dest] is None, row.flag


def test_every_help_line_comes_from_its_row_and_names_the_range_and_its_source():
    parser = _plan_parser()
    actions = {option: action for action in parser._actions for option in action.option_strings}
    for row in dials.LEVERS:
        assert actions[row.flag].help == row.help_text(), row.flag
        assert row.range_text() in row.help_text() and row.source in row.help_text(), row.flag
        assert row.status in row.help_text(), row.flag


def test_a_value_inside_the_range_is_taken_and_one_outside_is_refused_with_the_range_and_source(capsys):
    parser = _plan_parser()
    for row in dials.LEVERS:
        if row.kind == "enum":
            args = parser.parse_args([row.flag, row.choices[0]])
            assert getattr(args, row.dest) == row.choices[0]
            with pytest.raises(SystemExit):
                parser.parse_args([row.flag, "no-such-choice"])
            err = capsys.readouterr().err
            assert row.range_text() in err and row.source in err
            continue
        unit = "pc" if row.kind == "length" else ""
        for edge in (row.minimum, row.maximum):
            value = int(edge) if row.kind == "int" else edge
            assert getattr(parser.parse_args([row.flag, f"{value}{unit}"]), row.dest) == value
        for outside in (row.minimum - 1, row.maximum + 1):
            value = int(outside) if row.kind == "int" else outside
            with pytest.raises(SystemExit):
                parser.parse_args([f"{row.flag}={value}{unit}"])
            err = capsys.readouterr().err
            assert row.range_text() in err and row.source in err and row.status in err, row.flag


def test_a_bare_number_for_a_length_is_refused_with_the_units(capsys):
    row = dials.lever("disk_scale_length")
    with pytest.raises(SystemExit):
        _plan_parser().parse_args([row.flag, "2600"])
    assert row.unit_names() in capsys.readouterr().err


def test_a_length_is_taken_in_its_units_and_kept_in_the_base_unit():
    row = dials.lever("disk_scale_length")
    assert row.parse("2600pc") == 2600.0 and row.parse("2.6kpc") == pytest.approx(2600.0)
    assert row.parse(" 3 kpc ") == pytest.approx(3000.0)


def test_a_relative_unit_is_kept_as_given_and_gated_by_its_own_range():
    edge = dials.Lever(
        name="edge", flag="--edge", dial="edge", label="Edge", kind="length", unit="pc", minimum=1.0, maximum=262000.0,
        step=1.0, units=(dials.Unit("L", None, 3.0, 20.0), dials.Unit("pc", 1.0), dials.Unit("kpc", 1000.0),
                         dials.Unit("Mpc", 1.0e6)),
        law="where the galaxy ends", auto_text="the one-system-per-sector rule", help="Outer edge.",
        source="test row", status="placeholder")
    assert edge.parse("12L") == dials.RelativeLength(12.0, "L")
    assert edge.parse("0.03Mpc") == pytest.approx(30000.0)
    for bad in ("2L", "21L", "12", "0.3Mpc", "12 furlong"):
        with pytest.raises(argparse.ArgumentTypeError):
            edge.parse(bad)


@pytest.mark.parametrize("raw", ["nan", "inf", "-inf", "1e999", ""])
def test_a_non_finite_number_is_refused(raw):
    for row in dials.LEVERS:
        if row.kind in ("float", "int"):
            with pytest.raises(argparse.ArgumentTypeError):
                row.parse(raw)


def test_a_retired_flag_is_not_a_flag():
    parser = _plan_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["--disk-scale-length-pc", "2600"])


# --- the Generate form ------------------------------------------------------------------------------------------

def test_the_form_has_a_row_per_lever():
    assert generate_page.PLAN_FIELDS == dials.LEVERS


# --- the docs ---------------------------------------------------------------------------------------------------

def test_the_docs_carry_the_table_as_built():
    text = (ROOT / "docs" / "cli.md").read_text(encoding="utf-8")
    start = text.index(dials.DOCS_BEGIN.split(" (")[0])
    end = text.index(dials.DOCS_END) + len(dials.DOCS_END)
    assert text[start:end] == dials.docs_block(), "docs/cli.md is stale: run python3 -m planetgen.galaxy.dials --write"
    for row in dials.LEVERS:
        assert f"`{row.flag}`" in text[start:end]


def test_no_help_text_quotes_a_milky_way_number():
    """GEN.610: nothing fixed from the Milky Way; the help a user reads names no stand-in value."""
    for row in dials.LEVERS:
        for text in (row.help, row.law, row.auto_text):
            for banned in ("2.5e7", "1.5e9", "Milky Way", "Solar System"):
                assert banned not in text, (row.name, banned)


def test_the_gap_filling_plan_carries_the_range_registry_as_built():
    text = (ROOT / "docs" / "design" / "gap-filling-plan.md").read_text(encoding="utf-8")
    start = text.index(dials.REGISTRY_BEGIN.split(" (")[0])
    end = text.index(dials.REGISTRY_END) + len(dials.REGISTRY_END)
    assert text[start:end] == dials.registry_block(), "the plan is stale: run python3 -m planetgen.galaxy.dials --write"
    for row in dials.LEVERS:
        assert f"`{row.flag}`" in text[start:end]
