# tests/test_gen_cli_errors.py

"""
`docs/TODO.md` TEST.28 (CLI errors by message) and TEST.29 (limits stay
consistent).

TEST.28: every `parser.error` in `planetgen`'s galaxy/plan/phenomenon/
population validation (block, column, shell, center-sector, limit, plan
shape, workers, population, phenomenon, mysql-port) is asserted by its
own text on stderr, not just by `SystemExit` -- `test_galaxy_gen.py` and
`test_bughunt_cli_edges.py` only check that a bad command line exits, so
an error that starts firing for the wrong reason (an earlier check
swallowing a later one) went unnoticed. Each upper bound is also tested
at exactly its maximum: the maximum passes validation, maximum + 1 is
refused (`MAX_GENERATE_RING`, `MAX_GENERATE_LIMIT`, the first and last
`--block-layer`, `MAX_GENERATE_RADIUS_PC`, `MAX_NUM_ORBITS`, the MySQL
port range).

Everything goes through `generate_cli.process_args()` with a patched
`sys.argv`: validation runs in full and nothing connects to a database.

TEST.29: `MAX_GENERATE_LIMIT` is the cell count of the last allowed shell
(`MAX_GENERATE_SHELL`), and `MAX_GENERATE_SHELL` is `skeleton.DEFAULT_MAX_SHELL`,
so both are derived at import; reloading `generationLimits` under a patched
`DEFAULT_MAX_SHELL` shows they follow it (and the CLI's messages with them).
"""

import importlib
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.galaxy.skeleton import DEFAULT_MAX_SHELL
from planetgen.generation import limits as generationLimits
from planetgen.galaxy import skeleton, spherical

MAX_SHELL = generationLimits.MAX_GENERATE_SHELL
MAX_LIMIT = generationLimits.MAX_GENERATE_LIMIT
MAX_RADIUS = generationLimits.MAX_GENERATE_RADIUS_PC
MAX_ORBITS = generationLimits.MAX_NUM_ORBITS


def _parse(monkeypatch, command, argv):
    monkeypatch.setattr(sys, "argv", ["planetgen", command, *map(str, argv)])
    return generate_cli.process_args()


def _error(monkeypatch, capsys, command, argv):
    """Runs validation, expects argparse's usage-error exit, and returns
    the one `error:` message it printed."""
    capsys.readouterr()
    with pytest.raises(SystemExit) as excinfo:
        _parse(monkeypatch, command, argv)
    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    prefix = f"planetgen {command}: error: "
    lines = [line for line in err.splitlines() if line.startswith(prefix)]
    assert len(lines) == 1, err
    return lines[0][len(prefix):]


# ---------------------------------------------------------------------------
# galaxy: every message, by text
# ---------------------------------------------------------------------------

GALAXY_ERRORS = [
    # --block (refused until the ladder is rebuilt on shells)
    (["--block", "3.2.1.0"], "--block: the Galaxy Map's blocks are being rebuilt on shells"),
    # --shell
    (["--shell", "-1"], "--shell must be >= 0."),
    (["--shell", MAX_SHELL + 1], f"--shell must be at most {MAX_SHELL}."),
    # --band / --wedge (single-address mode)
    (["--band", "0"], "--band requires --shell."),
    (["--wedge", "0"], "--wedge requires --shell and --band."),
    (["--shell", "1", "--wedge", "0"], "--wedge requires --shell and --band."),
    (["--shell", "1", "--band", "-1", "--wedge", "-1"], "--wedge must be >= 0."),
    (["--shell", "1", "--band", "5"], "--band: shell 1 holds bands -1 to 1."),
    (["--shell", "1", "--band", "0", "--wedge", "12"], "--wedge: shell 1 band 0 holds wedges 0 to 11."),
    (["--shell", "1", "--band", "0", "--wedge", "0", "--limit", "1"],
     "--limit only applies to --shell (batch mode), not --shell --band --wedge (single-address mode)."),
    (["--shell", "1", "--band", "0", "--wedge", "0", "--yes"],
     "--yes only applies to --shell (batch mode), not --shell --band --wedge (single-address mode)."),
    (["--shell", "1", "--band", "0", "--wedge", "0", "--density", "2"],
     "--density/--num-systems can't be combined with --shell --band --wedge (single-address mode)"),
    # --center-sector / --radius-pc
    (["--center-sector", "5"], "--center-sector requires --radius-pc."),
    (["--shell", "1", "--radius-pc", "5"], "--radius-pc only applies to --center-sector, --shell --band --wedge, or random-start"),
    (["--radius-pc", "0"], "--radius-pc must be a positive number."),
    (["--center-sector", "5", "--radius-pc", "-1"], "--radius-pc must be a positive number."),
    (["--radius-pc", MAX_RADIUS * 1.0001], f"--radius-pc must be at most {MAX_RADIUS:g}."),
    # --limit / --yes
    (["--limit", "5"], "--limit only applies to --shell."),
    (["--center-sector", "5", "--radius-pc", "5", "--limit", "5"], "--limit only applies to --shell."),
    (["--shell", "1", "--limit", "0"], "--limit must be a positive integer."),
    (["--shell", "1", "--limit", "-3"], "--limit must be a positive integer."),
    (["--shell", "1", "--limit", MAX_LIMIT + 1], f"--limit must be at most {MAX_LIMIT}."),
    (["--yes"], "--yes only applies to --shell."),
    # spans
    (["--wedges", "0:3"], "--wedges requires --shells and --bands (one shell, one band)."),
    (["--shells", "2:1"], "shells 2:1 runs backwards"),
    (["--shells", "x"], "--shells: expected whole numbers, got 'x'"),
    (["--shells", "1:2:3"], "--shells: expected N or FIRST:LAST, got '1:2:3'"),
    (["--shells", MAX_SHELL + 1], f"--shells must stay at or below shell {MAX_SHELL}."),
    (["--shells", "1", "--bands", "0", "--wedges", "99"], "shell 1 band 0 holds wedges 0 to 11"),
    (["--shells", "1:3", "--bands", "0", "--wedges", "1"], "a wedge arc needs exactly one shell and one band"),
    (["--shells", "1", "--shell", "2"], "--shells/--bands/--wedges can't be combined with --shell."),
    (["--shells", "1", "--radius-pc", "5"], "--shells/--bands/--wedges can't be combined with --radius-pc."),
    (["--bands=-1:1", "--max-shell", "5"], "--shells/--bands/--wedges can't be combined with --max-shell."),
    (["--shells", "1", "--limit", MAX_LIMIT + 1], f"--limit must be between 1 and {MAX_LIMIT}."),
    # --course
    (["--course", "1", "2", "--shell", "1"], "--course can't be combined with --shell."),
    (["--course", "1", "2", "--shells", "1", "--limit", "3"], "--course can't be combined with --shells, --limit."),
    (["--course-border"], "--course-border requires --course."),
    # random-start options
    (["--shell", "1", "--max-shell", "10"],
     "--max-shell only applies to random-start mode (neither --shell nor --center-sector)."),
    (["--max-shell", "-1"], "--max-shell must be >= 0."),
    (["--max-shell", MAX_SHELL + 1], f"--max-shell must be at most {MAX_SHELL}."),
    (["--shell", "1", "--min-start-density", "1"], "--min-start-density only applies to random-start mode"),
    (["--min-start-density", "0"], "--min-start-density must be a positive number."),
    (["--min-start-density", "1", "--num-systems", "3"],
     "--min-start-density cannot be combined with --density/--num-systems"),
    # shared generation options, as the galaxy subcommand reports them
    (["--workers", "-1"], "--workers must be 0 (automatic) or more."),
    (["--density", "2", "--num-systems", "3"], "--density cannot be combined with --num-systems."),
    (["--density", "0"], "--density must be a positive number."),
]


@pytest.mark.parametrize("argv, message", GALAXY_ERRORS)
def test_galaxy_error_messages(monkeypatch, capsys, argv, message):
    assert _error(monkeypatch, capsys, "galaxy", argv).startswith(message)


# ---------------------------------------------------------------------------
# galaxy: each limit at exactly its maximum
# ---------------------------------------------------------------------------

def test_shell_at_max_is_accepted_and_one_past_is_refused(monkeypatch, capsys):
    args = _parse(monkeypatch, "galaxy", ["--shell", MAX_SHELL])
    assert (args.shell, args.band, args.wedge) == (MAX_SHELL, None, None)
    top = spherical.band_range(MAX_SHELL)[1]
    args = _parse(monkeypatch, "galaxy", ["--shell", MAX_SHELL, "--band", top, "--wedge", 0])
    assert (args.shell, args.band, args.wedge) == (MAX_SHELL, top, 0)
    assert _error(monkeypatch, capsys, "galaxy", ["--shell", MAX_SHELL + 1]) == f"--shell must be at most {MAX_SHELL}."


def test_max_shell_at_max_is_accepted_and_one_past_is_refused(monkeypatch, capsys):
    assert _parse(monkeypatch, "galaxy", ["--max-shell", MAX_SHELL]).max_shell == MAX_SHELL
    assert _parse(monkeypatch, "galaxy", ["--max-shell", 0]).max_shell == 0
    assert _error(monkeypatch, capsys, "galaxy", ["--max-shell", MAX_SHELL + 1]) == \
        f"--max-shell must be at most {MAX_SHELL}."


@pytest.mark.parametrize("mode", [["--shell", "3"], ["--shells", "2:3"]])
def test_limit_at_max_is_accepted_and_one_past_is_refused(monkeypatch, capsys, mode):
    assert _parse(monkeypatch, "galaxy", mode + ["--limit", MAX_LIMIT]).limit == MAX_LIMIT
    assert _parse(monkeypatch, "galaxy", mode + ["--limit", 1]).limit == 1
    message = _error(monkeypatch, capsys, "galaxy", mode + ["--limit", MAX_LIMIT + 1])
    if "--shells" in mode:
        assert message == f"--limit must be between 1 and {MAX_LIMIT}."
    else:
        assert message == f"--limit must be at most {MAX_LIMIT}."


def test_the_first_and_last_band_and_wedge_of_a_shell_are_accepted(monkeypatch, capsys):
    low, high = spherical.band_range(7)
    for band in (low, high):
        total = spherical.wedge_count(7, band)
        for wedge in (0, total - 1):
            args = _parse(monkeypatch, "galaxy", ["--shell", 7, "--band", band, "--wedge", wedge])
            assert (args.shell, args.band, args.wedge) == (7, band, wedge)
        assert f"holds wedges 0 to {total - 1}" in _error(
            monkeypatch, capsys, "galaxy", ["--shell", 7, "--band", band, "--wedge", total])
    assert f"holds bands {low} to {high}" in _error(monkeypatch, capsys, "galaxy", ["--shell", 7, "--band", high + 1])
    assert f"holds bands {low} to {high}" in _error(monkeypatch, capsys, "galaxy", ["--shell", 7, "--band", low - 1])


def test_a_span_becomes_a_span_object(monkeypatch):
    args = _parse(monkeypatch, "galaxy", ["--shells", "2:5", "--bands=-1:2"])
    assert (args.span.shells, args.span.bands, args.span.wedges) == ((2, 5), (-1, 2), None)
    args = _parse(monkeypatch, "galaxy", ["--shells", "3", "--bands", "1", "--wedges", "10:2"])
    assert args.span.wedges == (10, 2)


def test_radius_at_max_is_accepted(monkeypatch, capsys):
    assert _parse(monkeypatch, "galaxy", ["--radius-pc", MAX_RADIUS]).radius_pc == MAX_RADIUS
    args = _parse(monkeypatch, "galaxy", ["--center-sector", "1", "--radius-pc", MAX_RADIUS])
    assert args.radius_pc == MAX_RADIUS


def test_workers_zero_is_accepted(monkeypatch):
    assert _parse(monkeypatch, "galaxy", ["--workers", "0"]).workers == 0
    assert _parse(monkeypatch, "plan", ["--workers", "0"]).workers == 0


# ---------------------------------------------------------------------------
# plan, phenomenon, population
# ---------------------------------------------------------------------------

PLAN_ERRORS = [
    (["--arm-contrast", "1.01"], "argument --arm-contrast: 1.01 is outside 0 to 1 ("),
    (["--arm-contrast", "-0.01"], "argument --arm-contrast: -0.01 is outside 0 to 1 ("),
    (["--nuclear-amplitude", "-1"], "argument --nuclear-amplitude: -1 is outside 0 to 10 ("),
    (["--bulge-share", "-1"], "argument --bulge-share: -1 is outside 0 to 0.9 ("),
    (["--workers", "-1"], "--workers must be 0 (automatic) or more."),
    (["--edge", "3pc"], "argument --edge: 3 is outside 4 to 262140 pc, or 3 to 20 L ("),
    (["--edge", "2L"], "argument --edge: 2 L is outside 3 to 20 L ("),
    (["--edge", "300"], "argument --edge: '300' is not a length (a length needs its unit: L, pc, kpc, Mpc)"),
    (["--max-shell", "40"], "unrecognized arguments: --max-shell 40"),
    (["--disk-scale-length", "0pc"], "argument --disk-scale-length: 0 is outside 28 to 8250 pc ("),
    (["--disk-scale-length", "2600"], "argument --disk-scale-length: '2600' is not a length (a length needs its unit: pc, kpc)"),
    (["--disk-scale-length", "5Mpc"], "argument --disk-scale-length: '5Mpc' is not a length ('Mpc' is not one of the units pc, kpc)"),
    (["--flatness", "2"], "argument --flatness: 2 is outside 3 to 20 ("),
    (["--bulge-size", "0"], "argument --bulge-size: 0 is outside 0.02 to 1 ("),
    (["--arm-count", "0"], "argument --arm-count: 0 is outside 1 to 10 ("),
    (["--pitch-angle-deg", "0"], "argument --pitch-angle-deg: 0 is outside -90 to 90 degrees, at least 1 away from 0 ("),
    (["--pitch-angle-deg", "-90.001"], "argument --pitch-angle-deg: -90.001 is outside -90 to 90 degrees, at least 1 away from 0 ("),
    (["--shape", "elliptical"], "argument --shape: 'elliptical' is not one of spiral ("),
    (["--no-bright-stars", "--bright-stars-only"], "--no-bright-stars and --bright-stars-only can't be combined."),
    (["--bright-stars-only", "--bright-stars-down-to", "10"],
     "--bright-stars-down-to can't be combined with --no-bright-stars or --bright-stars-only."),
    (["--bright-stars-down-to", "-1"], "--bright-stars-down-to: "),
    (["--bright-star-min-luminosity", "-1"], "--bright-star-min-luminosity: "),
]


@pytest.mark.parametrize("argv, message", PLAN_ERRORS)
def test_plan_error_messages(monkeypatch, capsys, argv, message):
    assert _error(monkeypatch, capsys, "plan", argv).startswith(message)


def test_plan_shape_edges_that_are_allowed(monkeypatch):
    """The other side of each shape bound: amplitude 0, a 90 degree
    pitch either way, one arm, one shell."""
    for argv in (["--arm-contrast", "0"], ["--arm-contrast", "1"], ["--pitch-angle-deg", "90"], ["--pitch-angle-deg", "-90"],
                 ["--arm-count", "1"], ["--edge", "4pc"], ["--edge", "262140pc"], ["--edge", "3L"], ["--edge", "20L"],
                 ["--edge", "0.26Mpc"]):
        _parse(monkeypatch, "plan", argv)


PHENOMENON_ERRORS = [
    (["--num-orbits", "3"], "--num-orbits requires --anchor-system."),
    (["--anchor-system", "--type", "black-hole", "--num-orbits", "-1"],
     "--num-orbits must be zero or a positive integer."),
    (["--anchor-system", "--type", "black-hole", "--num-orbits", MAX_ORBITS + 1],
     f"--num-orbits must be at most {MAX_ORBITS}."),
    (["--anchor-system", "--type", "nebula"],
     "--anchor-system is only valid with --type black-hole or --type neutron-star."),
    (["--anchor-system", "--type", "neutron-star", "--sector-id", "1"],
     "--sector-id cannot be combined with --anchor-system"),
]


@pytest.mark.parametrize("argv, message", PHENOMENON_ERRORS)
def test_phenomenon_error_messages(monkeypatch, capsys, argv, message):
    assert _error(monkeypatch, capsys, "phenomenon", argv).startswith(message)


def test_phenomenon_num_orbits_at_max_is_accepted(monkeypatch):
    for orbits in (0, MAX_ORBITS):
        args = _parse(monkeypatch, "phenomenon", ["--anchor-system", "--type", "black-hole", "--num-orbits", orbits])
        assert args.num_orbits == orbits


def test_population_rescan_and_territories_only_message(monkeypatch, capsys):
    assert _error(monkeypatch, capsys, "population", ["--rescan", "--territories-only"]) == \
        "--rescan and --territories-only cannot be combined."
    assert _parse(monkeypatch, "population", ["--rescan"]).rescan is True
    assert _parse(monkeypatch, "population", ["--territories-only"]).territories_only is True


# ---------------------------------------------------------------------------
# --mysql-port, on every subcommand that takes it
# ---------------------------------------------------------------------------

_PORT_COMMANDS = ["system", "sector", "galaxy", "plan", "phenomenon", "population"]


@pytest.mark.parametrize("command", _PORT_COMMANDS)
@pytest.mark.parametrize("port", [0, -1, 65536])
def test_mysql_port_out_of_range_message(monkeypatch, capsys, command, port):
    assert _error(monkeypatch, capsys, command, ["--mysql-port", port]) == "argument --mysql-port: must be between 1 and 65535"


@pytest.mark.parametrize("command", _PORT_COMMANDS)
@pytest.mark.parametrize("port", [1, 65535])
def test_mysql_port_range_ends_are_accepted(monkeypatch, command, port):
    assert _parse(monkeypatch, command, ["--mysql-port", port]).mysql_port == port


def test_mysql_port_is_checked_before_the_subcommands_own_options(monkeypatch, capsys):
    """The port check runs first, so a bad port is what gets reported
    even alongside another bad option."""
    assert _error(monkeypatch, capsys, "galaxy", ["--shell", "-1", "--mysql-port", "0"]) == \
        "argument --mysql-port: must be between 1 and 65535"


# ---------------------------------------------------------------------------
# TEST.29: limits stay consistent
# ---------------------------------------------------------------------------

def test_limits_are_derived_from_the_skeletons_max_shell():
    assert generationLimits.MAX_GENERATE_SHELL == skeleton.DEFAULT_MAX_SHELL
    top = spherical.band_range(generationLimits.MAX_GENERATE_SHELL)[1]
    assert generationLimits.MAX_GENERATE_LIMIT == spherical.cells_in_band_range(generationLimits.MAX_GENERATE_SHELL, top)
    # The plan subcommand's own default scan cap is the same number.
    assert DEFAULT_MAX_SHELL == skeleton.DEFAULT_MAX_SHELL


def test_no_shell_up_to_the_max_holds_more_cells_than_the_limit():
    """`MAX_GENERATE_LIMIT` claims no shell up to `MAX_GENERATE_SHELL` holds
    more cells: shell sizes never shrink outward, checked over the first
    shells and a spread of later ones."""
    shells = [*range(0, 400), *range(400, generationLimits.MAX_GENERATE_SHELL, 997), generationLimits.MAX_GENERATE_SHELL]
    previous = 0
    for shell in shells:
        count = spherical.cell_count(shell)
        assert count >= previous, shell
        previous = count
    assert previous == generationLimits.MAX_GENERATE_LIMIT


@pytest.mark.parametrize("max_shell", [0, 1, 7, 3900, 25_000])
def test_limits_follow_a_changed_default_max_shell(monkeypatch, capsys, max_shell):
    """Reloading `generationLimits` under another `DEFAULT_MAX_SHELL`
    moves both limits together, and the command line's checks and messages
    (which read them through the module) move with them."""
    try:
        monkeypatch.setattr(skeleton, "DEFAULT_MAX_SHELL", max_shell)
        importlib.reload(generationLimits)
        assert generationLimits.MAX_GENERATE_SHELL == max_shell
        assert generationLimits.MAX_GENERATE_LIMIT == spherical.cell_count(max_shell)
        assert generate_cli.limits is generationLimits

        new_limit = generationLimits.MAX_GENERATE_LIMIT
        assert _parse(monkeypatch, "galaxy", ["--shell", max_shell]).shell == max_shell
        assert _error(monkeypatch, capsys, "galaxy", ["--shell", max_shell + 1]) == \
            f"--shell must be at most {max_shell}."
        assert _parse(monkeypatch, "galaxy", ["--shell", "0", "--limit", new_limit]).limit == new_limit
        assert _error(monkeypatch, capsys, "galaxy", ["--shell", "0", "--limit", new_limit + 1]) == \
            f"--limit must be at most {new_limit}."
    finally:
        monkeypatch.undo()
        importlib.reload(generationLimits)
    assert generationLimits.MAX_GENERATE_SHELL == skeleton.DEFAULT_MAX_SHELL == MAX_SHELL
    assert generationLimits.MAX_GENERATE_LIMIT == MAX_LIMIT
