"""OPS.196: no command goes silent by being written as one loop of its own. Every registered command gives the engine `pieces` and `do` (`queue/engine.py`), which count them into the unit's beat and
check for a cancel or a pause between two. The last of them (the clears and drops of OPS.200 and the reset of OPS.206) moved with OPS.196, so the list that held the commands still to convert is gone
and a command written with a `unit` function of its own fails here. The silence budget (every command run on a small galaxy, its longest stretch with no movement held to a budget) is
`test_longest_silence.py`."""

from planetgen.queue import commands


def test_every_command_gives_the_engine_its_pieces():
    unit_of_its_own = sorted(kind for kind, spec in commands.REGISTRY.items() if spec.pieces is None)
    assert not unit_of_its_own, f"commands with a loop of their own (give them `pieces` and `do`, so the engine counts them and can stop them between two): {unit_of_its_own}"


def test_the_commands_converted_last_are_the_clears_the_drops_and_the_reset():
    assert {kind for kind, spec in commands.REGISTRY.items() if spec.pieces is not None} >= {
        "clean-nebula-field", "clear-scatter", "empty-phenomenon-scatter", "drop-phenomenon-indexes", "clear-bright-stars", "empty-nebula-field", "reset-galaxy",
        "plan-skeleton", "finish-phenomena", "refresh-block-lights", "warm-cluster-field", "warm-object-index", "warm-galaxy-map", "warm-galaxy-nodes", "warm-galaxy-deep-nodes",
        "warm-lists", "warm-sectors", "warm-systems", "warm-phenomena", "warm-search", "end-map-rebuild-maintenance", "run-population", "run-territories",
        "link-neighbors", "settle-sector-paths", "generate-unplaced-sector"}
