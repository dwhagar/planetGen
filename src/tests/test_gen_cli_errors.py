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

TEST.29: `MAX_GENERATE_LIMIT` is `ring_sector_count(MAX_GENERATE_RING)`
and `MAX_GENERATE_RING` is `skeleton.DEFAULT_MAX_RING`, so both are
derived at import; reloading `generationLimits` under a patched
`DEFAULT_MAX_RING` shows they follow it (and the CLI's messages with them).
"""

import importlib
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.galaxy.skeleton import DEFAULT_MAX_RING
from planetgen.generation import run_galaxy
from planetgen.generation import limits as generationLimits
from planetgen.galaxy import skeleton
from planetgen.galaxy.drill import parse_drill_key
from planetgen.galaxy.geometry import ring_sector_count

MAX_RING = generationLimits.MAX_GENERATE_RING
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
    # --block
    (["--block", "nope"], "--block: not a block key: 'nope'"),
    (["--block", "3.2.1"], "--block: not a block key: '3.2.1'"),
    (["--block", "1.0.0.0"], "--block: block size must be one of (243, 27, 3), got 1"),
    (["--block", "9.0.0.0"], "--block: block size must be one of (243, 27, 3), got 9"),
    (["--block", "3.0.99.0"], "--block: no such block: '3.0.99.0'"),
    (["--block", "3.-1.0.0"], "--block: no such block: '3.-1.0.0'"),
    (["--block", "3.2.1.0", "--layer", "0"], "--block can't be combined with --layer."),
    (["--block", "3.2.1.0", "--slot", "0", "--radius-pc", "5", "--shell"],
     "--block can't be combined with --slot, --radius-pc, --shell."),
    (["--block", "3.2.1.0", "--max-ring", "5", "--min-start-density", "1", "--column"],
     "--block can't be combined with --max-ring, --min-start-density, --column."),
    (["--block", "3.2.1.0", "--block-layer", "5"], "--block-layer must be one of the block's layers, -1 to 1."),
    (["--block", "3.2.1.0", "--limit", "0"], f"--limit must be between 1 and {MAX_LIMIT}."),
    (["--block-layer", "0"], "--block-layer requires --block."),
    (["--block-layer", "0", "--ring", "1"], "--block-layer requires --block."),
    # --ring
    (["--ring", "-1"], "--ring must be >= 0."),
    (["--ring", MAX_RING + 1], f"--ring must be at most {MAX_RING}."),
    # --column / --shell
    (["--column", "--ring", "1"], "--column requires --ring and --slot."),
    (["--column", "--slot", "1"], "--column requires --ring and --slot."),
    (["--shell"], "--shell requires --ring."),
    (["--column", "--shell", "--ring", "1", "--slot", "0"], "--column and --shell can't be combined."),
    (["--column", "--ring", "1", "--slot", "0", "--layer", "0"],
     "--column and --shell cover every layer, so they don't take --layer."),
    (["--shell", "--ring", "1", "--layer", "0"],
     "--column and --shell cover every layer, so they don't take --layer."),
    (["--shell", "--ring", "1", "--slot", "0"], "--shell covers every slot of the ring; use --column for one slot."),
    (["--column", "--ring", "1", "--slot", "-1"], "--slot must be >= 0."),
    (["--column", "--ring", "1", "--slot", "0", "--limit", "2"], "--limit and --yes don't apply to --column."),
    (["--column", "--ring", "1", "--slot", "0", "--yes"], "--limit and --yes don't apply to --column."),
    # --layer / --slot (single-address mode)
    (["--layer", "0"], "--layer requires --ring."),
    (["--slot", "0"], "--slot requires --ring."),
    (["--ring", "1", "--slot", "-1"], "--slot must be >= 0."),
    (["--ring", "1", "--slot", "0", "--limit", "1"],
     "--limit only applies to --ring (batch mode), not --ring --slot (single-address mode)."),
    (["--ring", "1", "--slot", "0", "--yes"],
     "--yes only applies to --ring (batch mode), not --ring --slot (single-address mode)."),
    (["--ring", "1", "--slot", "0", "--density", "2"],
     "--density/--num-systems can't be combined with --ring --slot (single-address mode)"),
    # --center-sector / --radius-pc
    (["--center-sector", "5"], "--center-sector requires --radius-pc."),
    (["--ring", "1", "--radius-pc", "5"], "--radius-pc only applies to --center-sector, --ring --slot, or random-start"),
    (["--radius-pc", "0"], "--radius-pc must be a positive number."),
    (["--center-sector", "5", "--radius-pc", "-1"], "--radius-pc must be a positive number."),
    (["--radius-pc", MAX_RADIUS * 1.0001], f"--radius-pc must be at most {MAX_RADIUS:g}."),
    # --limit / --yes
    (["--limit", "5"], "--limit only applies to --ring."),
    (["--center-sector", "5", "--radius-pc", "5", "--limit", "5"], "--limit only applies to --ring."),
    (["--ring", "1", "--limit", "0"], "--limit must be a positive integer."),
    (["--ring", "1", "--limit", "-3"], "--limit must be a positive integer."),
    (["--ring", "1", "--limit", MAX_LIMIT + 1], f"--limit must be at most {MAX_LIMIT}."),
    (["--yes"], "--yes only applies to --ring."),
    # random-start options
    (["--ring", "1", "--max-ring", "10"],
     "--max-ring only applies to random-start mode (neither --ring nor --center-sector)."),
    (["--max-ring", "-1"], "--max-ring must be >= 0."),
    (["--max-ring", MAX_RING + 1], f"--max-ring must be at most {MAX_RING}."),
    (["--ring", "1", "--min-start-density", "1"], "--min-start-density only applies to random-start mode"),
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


def test_single_sector_block_message_is_unreachable_through_parse_drill_key():
    """`--block takes a block, not a single sector` guards `m == 1`, but
    `parse_drill_key` already refuses size 1 (`DRILL_LEVELS[:-1]`), so a
    one-sector key gets the block-size message instead."""
    with pytest.raises(ValueError, match="block size must be one of"):
        parse_drill_key("1.0.0.0")


# ---------------------------------------------------------------------------
# galaxy: each limit at exactly its maximum
# ---------------------------------------------------------------------------

def test_ring_at_max_is_accepted_and_one_past_is_refused(monkeypatch, capsys):
    args = _parse(monkeypatch, "galaxy", ["--ring", MAX_RING])
    assert (args.ring, args.layer) == (MAX_RING, 0)
    args = _parse(monkeypatch, "galaxy", ["--ring", MAX_RING, "--slot", ring_sector_count(MAX_RING) - 1])
    assert args.ring == MAX_RING
    assert _error(monkeypatch, capsys, "galaxy", ["--ring", MAX_RING + 1]) == f"--ring must be at most {MAX_RING}."


def test_max_ring_at_max_is_accepted_and_one_past_is_refused(monkeypatch, capsys):
    assert _parse(monkeypatch, "galaxy", ["--max-ring", MAX_RING]).max_ring == MAX_RING
    assert _parse(monkeypatch, "galaxy", ["--max-ring", 0]).max_ring == 0
    assert _error(monkeypatch, capsys, "galaxy", ["--max-ring", MAX_RING + 1]) == \
        f"--max-ring must be at most {MAX_RING}."


@pytest.mark.parametrize("mode", [["--ring", "3"], ["--shell", "--ring", "3"], ["--block", "3.2.1.0"]])
def test_limit_at_max_is_accepted_and_one_past_is_refused(monkeypatch, capsys, mode):
    assert _parse(monkeypatch, "galaxy", mode + ["--limit", MAX_LIMIT]).limit == MAX_LIMIT
    assert _parse(monkeypatch, "galaxy", mode + ["--limit", 1]).limit == 1
    message = _error(monkeypatch, capsys, "galaxy", mode + ["--limit", MAX_LIMIT + 1])
    if "--block" in mode:
        assert message == f"--limit must be between 1 and {MAX_LIMIT}."
    else:
        assert message == f"--limit must be at most {MAX_LIMIT}."


@pytest.mark.parametrize("key", ["3.2.1.0", "3.5.0.2", "3.0.0.-1", "27.1.0.0", "243.0.0.0"])
def test_first_and_last_block_layer_are_accepted_and_their_neighbors_refused(monkeypatch, capsys, key):
    layers = run_galaxy.block_layers(parse_drill_key(key))
    first, last = layers[0], layers[-1]
    assert layers == list(range(first, last + 1))
    for layer in (first, last):
        args = _parse(monkeypatch, "galaxy", ["--block", key, "--block-layer", layer])
        assert args.block_layer == layer and args.block == parse_drill_key(key)
    expected = f"--block-layer must be one of the block's layers, {first} to {last}."
    for layer in (first - 1, last + 1):
        assert _error(monkeypatch, capsys, "galaxy", ["--block", key, "--block-layer", layer]) == expected


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
    (["--arm-amplitude", "1"], "--arm-amplitude must be in [0, 1)."),
    (["--arm-amplitude", "-0.01"], "--arm-amplitude must be in [0, 1)."),
    (["--workers", "-1"], "--workers must be 0 (automatic) or more."),
    (["--max-ring", "0"], "--max-ring must be a positive integer."),
    (["--disk-scale-length-pc", "0"], "--disk-scale-length-pc must be a positive number."),
    (["--disk-scale-height-pc", "-1"], "--disk-scale-height-pc must be a positive number."),
    (["--bulge-scale-radius-pc", "0"], "--bulge-scale-radius-pc must be a positive number."),
    (["--arm-count", "0"], "--arm-count must be a positive integer."),
    (["--pitch-angle-deg", "0"], "--pitch-angle-deg must be non-zero and at most 90 degrees in magnitude."),
    (["--pitch-angle-deg", "-90.001"], "--pitch-angle-deg must be non-zero and at most 90 degrees in magnitude."),
    (["--calibration-radius-pc", "0"], "--calibration-radius-pc must be a positive number."),
    (["--calibration-radius-pc", "1e7"], "these galaxy shape parameters can't be normalized ("),
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
    pitch either way, one arm, one ring."""
    for argv in (["--arm-amplitude", "0"], ["--pitch-angle-deg", "90"], ["--pitch-angle-deg", "-90"],
                 ["--arm-count", "1"], ["--max-ring", "1"], ["--max-ring", MAX_RING]):
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
    assert _error(monkeypatch, capsys, "galaxy", ["--ring", "-1", "--mysql-port", "0"]) == \
        "argument --mysql-port: must be between 1 and 65535"


# ---------------------------------------------------------------------------
# TEST.29: limits stay consistent
# ---------------------------------------------------------------------------

def test_limits_are_derived_from_the_skeletons_max_ring():
    assert generationLimits.MAX_GENERATE_RING == skeleton.DEFAULT_MAX_RING
    assert generationLimits.MAX_GENERATE_LIMIT == ring_sector_count(generationLimits.MAX_GENERATE_RING)
    # The plan subcommand's own default scan cap is the same number.
    assert DEFAULT_MAX_RING == skeleton.DEFAULT_MAX_RING


def test_no_ring_up_to_the_max_holds_more_slots_than_the_limit():
    """`MAX_GENERATE_LIMIT` claims no ring up to `MAX_GENERATE_RING` holds
    more slots: ring sizes never shrink outward, checked over every ring."""
    previous = 0
    for ring_index in range(generationLimits.MAX_GENERATE_RING + 1):
        count = ring_sector_count(ring_index)
        assert count >= previous, ring_index
        previous = count
    assert previous == generationLimits.MAX_GENERATE_LIMIT


@pytest.mark.parametrize("max_ring", [0, 1, 7, 3900, 250_000])
def test_limits_follow_a_changed_default_max_ring(monkeypatch, capsys, max_ring):
    """Reloading `generationLimits` under another `DEFAULT_MAX_RING`
    moves both limits together, and the command line's checks and messages
    (which read them through the module) move with them."""
    try:
        monkeypatch.setattr(skeleton, "DEFAULT_MAX_RING", max_ring)
        importlib.reload(generationLimits)
        assert generationLimits.MAX_GENERATE_RING == max_ring
        assert generationLimits.MAX_GENERATE_LIMIT == ring_sector_count(max_ring)
        assert generate_cli.limits is generationLimits

        new_limit = generationLimits.MAX_GENERATE_LIMIT
        assert _parse(monkeypatch, "galaxy", ["--ring", max_ring]).ring == max_ring
        assert _error(monkeypatch, capsys, "galaxy", ["--ring", max_ring + 1]) == \
            f"--ring must be at most {max_ring}."
        assert _parse(monkeypatch, "galaxy", ["--ring", "0", "--limit", new_limit]).limit == new_limit
        assert _error(monkeypatch, capsys, "galaxy", ["--ring", "0", "--limit", new_limit + 1]) == \
            f"--limit must be at most {new_limit}."
    finally:
        monkeypatch.undo()
        importlib.reload(generationLimits)
    assert generationLimits.MAX_GENERATE_RING == skeleton.DEFAULT_MAX_RING == MAX_RING
    assert generationLimits.MAX_GENERATE_LIMIT == MAX_LIMIT
