"""Typed Generate commands (ADM.121, web/generate_commands.py): each line fills the Add sectors form so that it builds the
command the old form built, and a line that cannot be read says why."""

import pytest

from planetgen.web import generate_commands, generate_page


@pytest.mark.parametrize("line, argv", [
    ("add sectors", []),
    ("add sectors radius 50", ["--radius-pc", "50.0"]),
    ("Add Sectors radius 50pc density 1.5", ["--radius-pc", "50.0", "--min-start-density", "1.5"]),
    ("add sectors radius 100ly", ["--radius-pc", "30.7"]),
    ("add sectors around sector 9 radius 20", ["--center-sector", "9", "--radius-pc", "20.0"]),
    ("add sectors around shell 40 band -1 wedge 7 radius 12",
     ["--shell", "40", "--band", "-1", "--wedge", "7", "--radius-pc", "12.0"]),
    ("add sectors around shell 40 wedge 7 radius 12", ["--shell", "40", "--band", "0", "--wedge", "7", "--radius-pc", "12.0"]),
    ("add sectors around -10 0 radius 20", ["--shell", "3", "--band", "0", "--wedge", "9", "--radius-pc", "20.0"]),
    ("add sectors around 10 0 4 radius 20", ["--shell", "3", "--band", "1", "--wedge", "0", "--radius-pc", "20.0"]),
    ("add sectors shell 5", ["--shell", "5"]),
    ("add sectors shell 5 band -2 limit 40", ["--shell", "5", "--band", "-2", "--limit", "40"]),
    ("sectors shells 10:20 bands -2:2", ["--shells=10:20", "--bands=-2:2"]),
    ("add sectors shells 7 bands 0 wedges 4:9 limit 3", ["--shells=7", "--bands=0", "--wedges=4:9", "--limit", "3"]),
])
def test_a_typed_line_fills_the_form_the_old_form_would_have_filled(line, argv):
    command = generate_commands.parse(line)
    assert command.task == "add-sectors" and command.review
    built, _ = generate_page.galaxy_argv(command.fields, 4.0)
    assert built == argv


@pytest.mark.parametrize("line, task, fields, review", [
    ("check", "check", {"action": "check_db"}, True),
    ("check deep", "check", {"action": "check_db", "deep": "1"}, True),
    ("scatter", "scatters", {}, False),
    ("scatters", "scatters", {}, False),
    ("new galaxy", "new-galaxy", {}, False),
    ("plan", "plan-only", {}, False),
    ("plan only", "plan-only", {}, False),
    ("wipe", "wipe", {}, False),
])
def test_the_other_commands_open_a_task(line, task, fields, review):
    command = generate_commands.parse(line)
    assert (command.task, command.fields, command.review) == (task, fields, review)


@pytest.mark.parametrize("line", ["", "hello", "show me the settings", "addendum", "checking account", "wiped"])
def test_text_that_is_not_a_command_is_left_alone(line):
    assert generate_commands.parse(line) is None


@pytest.mark.parametrize("line, message", [
    ("add sectors foo", "I don't know 'foo'"),
    ("add sectors radius", "'radius' needs a value"),
    ("add sectors radius 5 radius 6", "given twice"),
    ("add sectors radius big", "Radius needs a number"),
    ("add sectors radius nan", "finite"),
    ("add sectors sector 9", "needs a radius"),
    ("add sectors sector x radius 5", "Sector needs a number"),
    ("add sectors around 1 radius 5", "x y or x y z"),
    ("add sectors around 1 2 shell 3 radius 5", "cannot be combined"),
    ("add sectors sector 9 shell 3 radius 5", "cannot be combined"),
    ("add sectors shells 3:5 radius 5", "cannot be combined"),
    ("add sectors wedge 3 radius 2", "needs its shell"),
    ("add sectors band 1", "goes with a shell"),
    ("check carefully", "only 'deep'"),
    ("wipe everything", "opens its form"),
])
def test_a_line_that_cannot_be_read_says_why(line, message):
    with pytest.raises(generate_commands.CommandError, match=message):
        generate_commands.parse(line)
